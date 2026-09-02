#!/publicfs10/fs10-m9/home/m9s005343/software/anaconda3/bin/python
from pymatgen.io.vasp.inputs import Kpoints
from pymatgen.core import Structure
from pymatgen.symmetry.bandstructure import HighSymmKpath

struct=Structure.from_file("POSCAR")
kpath=HighSymmKpath(struct,symprec=0.00001)
kpts=Kpoints.automatic_linemode(divisions=30,ibz=kpath)
kpts.write_file("KPOINTS")
