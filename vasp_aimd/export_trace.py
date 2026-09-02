# -*- coding: utf-8 -*-
"""
Created on Wed Apr  9 20:10:32 2025

@author: dengjun
"""

import numpy as np
from pymatgen.io.vasp import Xdatcar
import matplotlib.pyplot as plt
from pymatgen.core.structure import Structure
import sys

xd=Xdatcar('XDATCAR')
pressure=int(sys.argv[1])
temperature=int(sys.argv[2])

H_symbols=[]
H_coords=[]
C_symbols=[]
C_coords=[]
O_symbols=[]
O_coords=[]
O1=[]
H1=[]
C1=[]
nstruct=len(xd.structures)
step=20
for num in np.arange(0,nstruct,step):
    structure=xd.structures[num]
    H_symbols.extend(structure.atomic_numbers[:96])
    H_coords.extend(structure.frac_coords[:96])
    H1.extend(structure.frac_coords[0,:])
    C_symbols.extend(structure.atomic_numbers[96:144])
    C_coords.extend(structure.frac_coords[96:144])
    C1.extend(structure.frac_coords[96,:])
    O_symbols.extend(structure.atomic_numbers[144:288])
    O_coords.extend(structure.frac_coords[144:288])
    O1.extend(structure.frac_coords[144,:])

matrix=structure.lattice.matrix

s_H=Structure(matrix,H_symbols,H_coords)
s_H.to(fmt='cif',filename='%dbar_%dK_H.cif'%(pressure,temperature))
s_C=Structure(matrix,C_symbols,C_coords)
s_C.to(fmt='cif',filename='%dbar_%dK_C.cif'%(pressure,temperature))
s_O=Structure(matrix,O_symbols,O_coords)
s_O.to(fmt='cif',filename='%dbar_%dK_O.cif'%(pressure,temperature))

# S_H1=Structure(matrix,len(H1)*[1],H1)
    
symbols=[]
symbols.extend(H_symbols)
symbols.extend(C_symbols)
symbols.extend(O_symbols)

coords=[]
coords.extend(H_coords)
coords.extend(C_coords)
coords.extend(O_coords)
s=Structure(matrix,symbols,coords)
s.to(fmt='cif',filename='%dbar_%dK.cif'%(pressure,temperature))
