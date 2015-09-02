# PyMol script to generate/load part of the &QMMM section of CP2K input
# Copyright (C) 2015 Peter Mamonov <pmamonov@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# USAGE:
#   gen_cp2k_qmmm filename, qm_part_selection_name, the_whole_system_object
#
#   N.B.: The output is written to file "filename". "Filename" is overwritten!
#
# REFERENCES:
#   CP2K: http://cp2k.berlios.de/
#   PyMol: http://www.pymol.org/

from pymol import cmd, stored
import re

cmd.set("retain_order", value="1")

def gen_cp2k_qmmm (filename, qmsele, mmsele):
	"""
Generate &QMMM section of CP2K input

N.B.: run "set retain_order, 1" command  before loading QMMM object to PyMol
      to avoid atoms indices mismatch between PyMol and CP2K

N.B.: The output is written to file "filename". "filename" will be overwritten!

USAGE:
        gen_cp2k_qmmm filename, qm_part_selection_name, the_whole_system_object

	"""
	f = open(filename, 'w')
	print >> f, "&QMMM"

# QMMM CELL setup
# n.b. 6A padding
# n.b. cell should be cubic in case of PSOLVER=WAVELET
	ext = cmd.get_extent(qmsele)
	cell = (ext[1][0]-ext[0][0] + 6, ext[1][1]-ext[0][1] + 6, ext[1][2]-ext[0][2] + 6)
	print >> f, "  &CELL\n    ABC %.2f %.2f %.2f\n  &END CELL" % cell

# KINDS
	stored.elems=[]
	cmd.iterate_state(1, qmsele, "stored.elems.append(elem)")
	kinds = set(stored.elems)
	for kind in kinds:
		print >> f, "  &QM_KIND", kind
		stored.ids=[]
		cmd.iterate_state(1, qmsele + ' & elem '+kind, "stored.ids.append(ID)")
		print >> f, "    MM_INDEX", 
		for id in stored.ids: 
			print >> f, id,
		print >> f, "\n  &END QM_KIND"

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
		if link: print >> f, "  &LINK\n    QM_INDEX %d\n    MM_INDEX %d\n    QM_KIND H\n  &END LINK" % (qmid, mmid)
 
	print >> f, "&END QMMM"

cmd.extend("gen_cp2k_qmmm", gen_cp2k_qmmm)

re_mmindex = re.compile("^\s*MM_INDEX\s+([\s\d]+)$")

def load_cp2k_qmmm(filename, obj_name, qmsele_name):
	"""
Load &QMMM section of CP2K input into selection

N.B.: run "set retain_order, 1" command  before loading QMMM object to PyMol
      to avoid atoms indices mismatch between PyMol and CP2K

USAGE:
        load_cp2k_qmmm filename, obj_name, selection_name
	"""
	atcount=0
	inp = open(filename, 'r')
	in_qmkind=0
	for str in inp:
		if (in_qmkind):
			if (re.match("\s*&END", str)):
				in_qmkind=0
		if (in_qmkind):
			str_mmindex=re_mmindex.match(str)
			if (str_mmindex):
				for id in str_mmindex.group(1).split():
					if (atcount):
						cmd.select(qmsele_name, obj_name+" & ("+qmsele_name+" | "+"index "+id+")" )
						atcount += 1
					else:
						cmd.select(qmsele_name, obj_name+" & index "+id)
						atcount += 1
		if (re.match('\s*&QM_KIND', str)):
			in_qmkind=1
	print atcount, "atoms selected"

cmd.extend("load_cp2k_qmmm", load_cp2k_qmmm)
