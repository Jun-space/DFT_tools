#!/bin/bash
/ur/bin/python2 /home/dj/Documents/abinit/scripts/cif2vasp/cif2vasp.py $1
sed -n '1,5p' POSCAR >> POSCAR2
cat POSCAR >> POSCAR2 
sed -i '7,10d' POSCAR2
mv POSCAR2 POSCAR
