#!/publicfs10/fs10-m9/home/m9s005343/software/anaconda3/bin/python

from pymatgen.core.structure import Structure
import sys

inputfilename=sys.argv[1]
outputfilename=sys.argv[2]
s=Structure.from_file(inputfilename)
s.to(fmt='poscar',filename=outputfilename)
