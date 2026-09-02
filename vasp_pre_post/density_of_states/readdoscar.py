# -*- coding: utf-8 -*-
"""
Created on Mon Mar 13 09:37:04 2023

@author: DengJun
"""


import numpy as np


print('reading DOSCAR')
f=open('DOSCAR')
lines=f.readlines()
f.close()

[nion1, nion, pdos, ncdij] = [int(i) for i in lines[0].split()]
[vol,a,b,c,potdim] = [float(i) for i in lines[1].split()] #A^3, m,m,m,s
[Emin, Emax, nedos, efermi, t] = [float(i) for i in lines[5].split()]


##############read TDOS
tdos=[]
for i in np.arange(6,int(6+nedos)):
    tdos.append([float(j) for j in lines[i].split()])

tdos=np.array(tdos)

np.savetxt('tdos.dat',tdos,'%.5f')

header='Energy      s          py         pz         px         dxy        dyz        dz2        dxz        dx2        Total'
##############read PDOS
pdos={}
for i in np.arange(0,nion):
    print(i)
    pdos[i]=[]
    lstart=int((i+1)*(nedos+1)+6)
    lend=int((i+2)*(nedos+1)+5)
    for j in np.arange(lstart,lend):
        pdos[i].append([float(k) for k in lines[j].split()])
    pdos[i]=np.array(pdos[i])
    t=np.sum(pdos[i][:,1:],axis=1)
    pdos[i]=np.concatenate((pdos[i],t.reshape(len(t),1)),axis=1)
    np.savetxt('%d.dat'%i,pdos[i],fmt='%5.5f',header=header)

f=open('POSCAR')
lines=f.readlines()
f.close()

elements=lines[5].split()
nelements=[int(i) for i in lines[6].split()]
pdos_ele={}
[a,b]=np.shape(pdos[1])
for num,i in enumerate(nelements):
    pdos_ele[elements[num]]=np.zeros([a,b-1])
    if num==0:
        begin=0
        end=i
    else:
        begin=np.sum(nelements[:num])
        end=np.sum(nelements[:num+1])
    for j in range(begin,end):
        pdos_ele[elements[num]]=pdos_ele[elements[num]]+pdos[j][:,1:]
    pdos_ele[elements[num]]=np.concatenate((pdos[j][:,0].reshape(a,1),pdos_ele[elements[num]]),axis=1)
    np.savetxt(elements[num]+'.dat',pdos_ele[elements[num]],fmt='%10.5f',header=header)