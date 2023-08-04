from pymol import cmd, stored
import re

cmd.set("retain_order", value="1")

def gen_cp2k_qmmm(filename, qmsele, mmsele):
    """
    Generate &QMMM section of CP2K input

    N.B.: run "set retain_order, 1" command  before loading QMMM object to PyMol
          to avoid atoms indices mismatch between PyMol and CP2K

    N.B.: The output is written to file "filename". "filename" will be overwritten!

    USAGE:
            gen_cp2k_qmmm(filename, qm_part_selection_name, the_whole_system_object)

    """
    with open(filename, 'w') as f:
        f.write("&QMMM\n")

        # QMMM CELL setup
        # n.b. 6A padding
        # n.b. cell should be cubic in case of PSOLVER=WAVELET
        ext = cmd.get_extent(qmsele)
        cell = (ext[1][0] - ext[0][0] + 6, ext[1][1] - ext[0][1] + 6, ext[1][2] - ext[0][2] + 6)
        f.write("  &CELL\n    ABC %.2f %.2f %.2f\n  &END CELL\n" % cell)

        # KINDS
        stored.elems = []
        cmd.iterate_state(1, qmsele, "stored.elems.append(elem)")
        kinds = set(stored.elems)
        for kind in kinds:
            f.write("  &QM_KIND %s\n" % kind)
            stored.ids = []
            cmd.iterate_state(1, qmsele + ' & elem ' + kind, "stored.ids.append(ID)")
            f.write("    MM_INDEX")
            for id in stored.ids:
                f.write(" %d" % id)
            f.write("\n  &END QM_KIND\n")

        # LINKS
        stored.qm_ids = []
        cmd.iterate_state(1, qmsele, "stored.qm_ids.append(ID)")
        model = cmd.get_model(mmsele)
        for bond in model.bond:
            [ia, ib] = bond.index
            a = model.atom[ia].id
            b = model.atom[ib].id
            link = 0
            if ((a in stored.qm_ids) and not (b in stored.qm_ids)):
                qmid = a
                mmid = b
                link = 1
            if ((b in stored.qm_ids) and not (a in stored.qm_ids)):
                mmid = a
                qmid = b
                link = 1
            if link:
                f.write("  &LINK\n    QM_INDEX %d\n    MM_INDEX %d\n    QM_KIND H\n  &END LINK\n" % (qmid, mmid))

        f.write("&END QMMM\n")

cmd.extend("gen_cp2k_qmmm", gen_cp2k_qmmm)

re_mmindex = re.compile("^\s*MM_INDEX\s+([\s\d]+)$")

def load_cp2k_qmmm(filename, obj_name, qmsele_name):
    """
    Load &QMMM section of CP2K input into selection

    N.B.: run "set retain_order, 1" command  before loading QMMM object to PyMol
          to avoid atoms indices mismatch between PyMol and CP2K

    USAGE:
            load_cp2k_qmmm(filename, obj_name, selection_name)
    """
    atcount = 0
    with open(filename, 'r') as inp:
        in_qmkind = 0
        for str_line in inp:
            if (in_qmkind):
                if re.match(r"\s*&END", str_line):
                    in_qmkind = 0
            if (in_qmkind):
                str_mmindex = re_mmindex.match(str_line)
                if str_mmindex:
                    for id in str_mmindex.group(1).split():
                        if atcount:
                            cmd.select(qmsele_name, obj_name + " & (" + qmsele_name + " | " + "index " + id + ")")
                            atcount += 1
                        else:
                            cmd.select(qmsele_name, obj_name + " & index " + id)
                            atcount += 1
            if re.match(r'\s*&QM_KIND', str_line):
                in_qmkind = 1
    print(atcount, "atoms selected")

cmd.extend("load_cp2k_qmmm", load_cp2k_qmmm)

