# -*- coding: utf-8 -*-
"""
Created on Thu May 29 10:16:04 2025

@author: dengjun
"""

import os
from pymatgen.core.structure import Structure
import numpy as np

cp_numb=[int(i) for i in os.popen('grep \'Critical point no\' cps.out | awk \'{print $5}\'').readlines()]
types=os.popen('grep \'Type :\' cps.out | awk \'{print $3}\'').readlines()
rho=[float(i) for i in os.popen('grep \'Field value (f):\' cps.out | awk \'{print $4}\'').readlines()]
delta_rho=[[float(j) for j in i.split()] for i in os.popen('grep \'Gradient (grad f):\' cps.out | awk \'{print $4,$5,$6}\'').readlines()]
delta_rho_norm=[float(i) for i in os.popen('grep \'Gradient norm (|grad f|):\' cps.out | awk \'{print $5}\'').readlines()]
delta2_rho=[float(i) for i in os.popen('grep \'Laplacian (del2 f):\' cps.out | awk \'{print $4}\'').readlines()]
eigenvalues=[[float(j) for j in i.split()] for i in os.popen('grep \'Hessian eigenvalues:\' cps.out | awk \'{print $3,$4,$5}\'').readlines()]
coord=[[float(j) for j in i.split()] for i in os.popen('grep \'Crystallographic coordinates:\' cps.out | awk \'{print $3,$4,$5}\'').readlines()]
cps={}
datas=[]
f=open('cps.dat','a')
print('#num, type, rho, delta_rho_norm, delta2_rho, lambda1, lambda2, lambda3',file=f)
for i in cp_numb:
    cps[i]={'num':i,'type':types[i-1][:-1],'rho':rho[i-1],'delta_rho':delta_rho[i-1],'delta_rho_norm':delta_rho_norm[i-1],\
            'delta2_rho':delta2_rho[i-1],'eigenvalue':eigenvalues[i-1],'coordinate':coord[i-1]}
    print(i,types[i-1][:-1],rho[i-1],delta_rho_norm[i-1],delta2_rho[i-1],eigenvalues[i-1][0],eigenvalues[i-1][1],eigenvalues[i-1][2],file=f)
f.close()
np.save('cps.npy',np.array(cps))

s=Structure.from_file('POSCAR')
lattice=s.lattice.matrix
atoms=s.frac_coords
symbols=s.atomic_numbers

bcp=[]
for i in cp_numb:
    if cps[i]['type'] == '(3,-1)':
        bcp.append(cps[i]['coordinate'])

atoms=atoms.tolist()
atoms.extend(bcp)

symbols=list(symbols)
symbols.extend(len(bcp)*[118])
crystal=Structure(lattice,symbols,atoms)
crystal.to(fmt='poscar',filename='bcp.vasp')