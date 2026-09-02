#!/publicfs10/fs10-m9/home/m9s005343/software/anaconda3/bin/python

from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.core.structure import Structure
s=Structure.from_file('POSCAR')
s1=SpacegroupAnalyzer(s,symprec=0.0001).get_primitive_standard_structure(international_monoclinic=False)
s1.to(filename='POSCAR.vasp',fmt='poscar')
