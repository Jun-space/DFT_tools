# -*- coding: utf-8 -*-
"""
Created on Thu Apr 29 10:27:05 2021

@author: DengJun
"""


import numpy as np
import matplotlib.pyplot as plt
from numpy import linalg as LA
from matplotlib import rcParams


def read_CONTCAR(filename):
    file = open(filename,'r')
    file.readline()
    scale = float(file.readline())
    a1 = file.readline().split()
    a2 = file.readline().split()
    a3 = file.readline().split()
    lattice = [a1, a2, a3]
    lattice = np.array(lattice, dtype='float') * scale
    file.close()
    return lattice

def get_reciprocal(lattice):  # this funciton is to get reciprocal lattice from primitive lattice

    v = np.dot(lattice[0], np.cross(lattice[1], lattice[2]))
    reciprocal = []
    reciprocal.append(2 * np.pi * np.cross(lattice[1], lattice[2]) / v)
    reciprocal.append(2 * np.pi * np.cross(lattice[2], lattice[0]) / v)
    reciprocal.append(2 * np.pi * np.cross(lattice[0], lattice[1]) / v)
    reciprocal = np.array(reciprocal, dtype='float')
    return reciprocal

def point_scale(pt_list, a):
    # point_scale, including rpt and kpt, if it is rpt, put a as lattice, if it is kpt, put a as inversed lattice
    pt_scaled=[]
    for pt in pt_list:
        pt_scaled.append(np.dot(pt, a))
    pt_scaled = np.array(pt_scaled, dtype='float')
    return pt_scaled


file = open(r'wannier90_hr.dat','r')
file.readline()
num_wann = int(file.readline().split()[0])
nrpts = int(file.readline().split()[0])
weight = []
for i in range(int(np.ceil(nrpts / 15.0))):
    buffer = file.readline().split()
    weight = weight + buffer
weight = np.array(weight, dtype='int')
rpt = []
hamr = np.zeros((num_wann, num_wann, nrpts), dtype='complex')

for i in range(nrpts):
    for j in range(num_wann):
        for k in range(num_wann):
            buffer = file.readline().split()
            hamr[k, j, i] = float(buffer[5]) + 1j * float(buffer[6])
    rpt = rpt + [buffer[0:3]]

rpt = np.array(rpt, dtype='int')
    
lattice = read_CONTCAR('POSCAR')
reciprocal = get_reciprocal(lattice)
hr = {'num_wann': num_wann, 'nrpts': nrpts, 'weight': weight, 'rpt': rpt, 'hamr': hamr, 'lattice': lattice, 'reciprocal': reciprocal}
file.close()
print("Reading processe is done (hr)")

###############################method1##########################
kpt_list=np.loadtxt('wannier90_band.kpt')

kpt_list = point_scale(kpt_list, hr['reciprocal'])
rpts = point_scale(hr['rpt'], hr['lattice'])
phase = np.exp(1j * np.dot(rpts, kpt_list.T))/hr['weight'][:, None]

hamk = np.tensordot(hr['hamr'], phase, axes=1)

hamk = np.rollaxis(hamk, 2, 0)
evalus= LA.eigvalsh(hamk)


#####################method2##########################
kpt_list=np.loadtxt('wannier90_band.kpt')

kpt_list = point_scale(kpt_list, hr['reciprocal'])
rpts = point_scale(hr['rpt'], hr['lattice'])
Hmn_phase2=np.zeros((num_wann,num_wann,nrpts,len(kpt_list)),dtype='complex')
for k in range(len(kpt_list)):
    print(k)
    for m in range(num_wann):
        for n in range(num_wann):
            for l in range(len(rpts)):
                Hmn_phase2[m,n,l,k]=np.exp(1j*np.dot(kpt_list[k],rpts[l]))*hr['hamr'][m,n,l]/hr['weight'][l]

hamk2=np.sum(Hmn_phase2,axis=2)
evalus2=LA.eigvalsh(np.rollaxis(hamk2, 2, 0))


#####################method3##########################
kpt_list=np.loadtxt('wannier90_band.kpt')
wannier_center=np.loadtxt('wannier_center.dat')

kpt_list = point_scale(kpt_list, hr['reciprocal'])
rpts = point_scale(hr['rpt'], hr['lattice'])
Hmn_phase3=np.zeros((num_wann,num_wann,nrpts,len(kpt_list)),dtype='complex')
for k in range(len(kpt_list)):
    print(k)
    for m in range(num_wann):
        for n in range(num_wann):
            for l in range(len(rpts)):
                Hmn_phase3[m,n,l,k]=np.exp(1j*np.dot(kpt_list[k],rpts[l]-wannier_center[m]+wannier_center[n]))*hr['hamr'][m,n,l]/hr['weight'][l]

hamk3=np.sum(Hmn_phase3,axis=2)
evalus3=LA.eigvalsh(np.rollaxis(hamk3, 2, 0))

x=np.zeros(len(kpt_list))

for i in range(len(kpt_list)-1):
    x[i+1]=x[i]+np.linalg.norm(kpt_list[i+1]-kpt_list[i])
x=x/np.max(x)
for i in np.arange(len(evalus[0,:])):
    plt.plot(x,evalus[:,i],'--',color='r',linewidth=0.5)
plt.plot(x,evalus[:,i],'--',color='r',linewidth=0.5,label='e1')

for i in np.arange(len(evalus2[0,:])):
    plt.plot(x,evalus2[:,i],'--',color='g',linewidth=0.5)
plt.plot(x,evalus2[:,i],'--',color='g',linewidth=0.5,label='e2')

for i in np.arange(len(evalus[0,:])):
    plt.plot(x,evalus3[:,i],'--',color='b',linewidth=0.5)
plt.plot(x,evalus3[:,i],'--',color='b',linewidth=0.5,label='e3')

# wannier_center=np.loadtxt('wannier_center.dat')
# rpts = point_scale(hr['rpt'], hr['lattice'])
# Rmn=np.zeros((num_wann,num_wann,nrpts,3),dtype='complex')

# for i in range(nrpts):
#     for j in range(num_wann):
#         for k in range(num_wann):
#             Rmn[k,j,i,:]=rpts[i]+wannier_center[k]-wannier_center[j]
            
# phase2 = np.exp(1j * np.dot(Rmn, kpt_list.T))/hr['weight'][:, None]

# hamk=np.zeros((len(kpt_list[:,0]),num_wann,num_wann),dtype='complex')
# #hamk = np.tensordot(hr['hamr'], phase2, axes=1)
# for i in range(len(kpt_list[:,0])):
#     for m in range(num_wann):
#         for n in range(num_wann):
#             for k in range(len(rpts)):
#                 hamk[i,m,n]=hamk[i,m,n]+np.exp(np.dot(1j*kpt_list[i],(rpts[k]+wannier_center[m]-wannier_center[n])))*hr['hamr'][m,n,k]/hr['weight'][k]
        