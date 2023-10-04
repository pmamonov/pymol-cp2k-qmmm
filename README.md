PyMol script to generate/load part of the &QMMM section of CP2K input

Copyright (C) 2015 Peter Mamonov <pmamonov@gmail.com>

Originally written by Peter Mamonov. Adapted to Python3 and PyMOL2.5 by Miquel Canyelles.

# USAGE

```
PyMOL> import sys
PyMOL> sys.path.append("/path/to/pymol-cp2k-qmmm")
PyMOL> import pymol_cp2k_qmmm

PyMOL> set retain_order, 1
PyMOL> load_cp2k_qmmm filename, the_whole_system_object, qm_part_selection_name

PyMOL> gen_cp2k_qmmm filename, qm_part_selection_name, the_whole_system_object
```

N.B.: The output is written to file "filename". "Filename" is overwritten!

# REFERENCES

- CP2K: http://cp2k.berlios.de/
- PyMol: http://www.pymol.org/
