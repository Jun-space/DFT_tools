#!/home/dj/software/anaconda3_2019.07/envs/abipy/bin/python
from pymatgen.core.structure import Structure
import sys

s=Structure.from_file('%s'%sys.argv[1])
s.to(fmt='poscar',filename='POSCAR')
