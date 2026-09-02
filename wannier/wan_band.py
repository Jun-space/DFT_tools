# -*- coding: utf-8 -*-
"""
Created on Wed Nov 25 11:16:13 2020

@author: DengJun
"""


import numpy as np
import matplotlib.pyplot as plt
from numpy import linalg as LA
from matplotlib import rcParams

rcParams['mathtext.fontset'] = 'custom'
rcParams['mathtext.rm'] = 'Arial'
rcParams['mathtext.it'] = 'Arial:italic'
rcParams['mathtext.bf'] = 'Arial:bold'

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

def gen_kpt_grid(nk0,nk1,nk2):
    nkpt = nk0 * nk1 * nk2
    kx = np.linspace(0,1,nk0, endpoint=False)
    ky = np.linspace(0,1,nk1, endpoint=False)
    kz = np.linspace(0,1,nk2, endpoint=False)

    kpt_list = []
    kpt_label = []
    for i in range(nk0):
        for j in range(nk1):
            for k in range(nk2):
                kpt_list.append([kx[i], ky[j], kz[k]])
                kpt_label.append(j)
    return kpt_list


def read_hr(filename):
    file = open(filename,'r')
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
    return hr

def fourier(kpt_list, hr):

    kpt_list = point_scale(kpt_list, hr['reciprocal'])
    rpts = point_scale(hr['rpt'], hr['lattice'])
    phase = np.exp(1j * np.dot(rpts, kpt_list.T))/hr['weight'][:, None]

    hamk = np.tensordot(hr['hamr'], phase, axes=1)

    return hamk

def eigenvalue(kpt_list, hr):

    hamk = fourier(kpt_list, hr)
    hamk = np.rollaxis(hamk, 2, 0)
    evalus= LA.eigvalsh(hamk)

    return evalus


dim=[90,90,1]
efermi=  0.94118204 #fermi energy
M1=[0.0,0.5,0.0]
K1=[1/3,2/3,0.0]
kx=np.linspace(M1[0],K1[0],1000)
ky=np.linspace(M1[1],K1[1],1000)
kz=np.zeros(1000)
kpath=np.concatenate([kx.reshape(len(kx),1),ky.reshape(len(ky),1),kz.reshape(len(kz),1)],axis=1).tolist()

hr=read_hr('wannier90_hr.dat')
# kgrid=gen_kpt_grid(dim[0],dim[1],dim[2])

k_evals=eigenvalue(kpath,hr)-efermi
kpath_real=point_scale(kpath,hr['reciprocal'])
x=np.zeros(len(kpath))
for i in range(len(kpath)-1):
    x[i+1]=x[i]+np.linalg.norm(kpath_real[i+1]-kpath_real[i])
x=x/np.max(x)
fig=plt.figure()
ax=fig.add_subplot(111)
# for i in np.arange(len(k_evals[0,:])):
#     ax.plot(x,k_evals[:,i],'--',color='r',linewidth=0.5)
# ax.plot(x,k_evals[:,117],'--',color='r',linewidth=0.5,label='TB')
ax.plot(x,k_evals[:,118],'--',color='r',linewidth=0.5,label='TB')
ax.plot(x,k_evals[:,119],'--',color='r',linewidth=0.5,label='TB')

