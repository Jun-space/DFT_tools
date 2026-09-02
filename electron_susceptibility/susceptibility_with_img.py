# -*- coding: utf-8 -*-
"""
Created on Thu Dec 17 14:30:45 2020

@author: DengJun
"""


import numpy as np
import matplotlib.pyplot as plt
from numpy import linalg as LA
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import rcParams
rcParams['mathtext.fontset'] = 'custom'
rcParams['mathtext.rm'] = 'Times New Roman'
rcParams['mathtext.it'] = 'Times New Roman:italic'
rcParams['mathtext.bf'] = 'Times New Roman:bold'

def delta_function(E,efermi):
    epsilon=10**-2
    x=E-efermi
    return (1/np.pi)*epsilon/(x**2+epsilon**2)

def fermi_dirac(E,efermi,T):
    kB=8.617333262145*10**-5 # eV/K
    return 1/(np.exp((E-efermi)/(kB*T))+1)

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
    kx = np.linspace(10**-5,1-10**-5,nk0, endpoint=True)
    ky = np.linspace(10**-5,1-10**-5,nk1, endpoint=True)
    kz = np.linspace(10**-5,1-10**-5,nk2, endpoint=True)

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



efermi= 1.79852 #fermi energy
hr=read_hr('wannier90_hr.dat')
band_cross_fermi1=[60,61,62,63]
band_cross_fermi2=[60,61,62,63]

nk=[50,50,1]
#kgrid=gen_kpt_grid(50,50,32)
kgrid=gen_kpt_grid(nk[0],nk[1],nk[2])
qgrid=np.array(kgrid)

k_evals=eigenvalue(kgrid,hr)

T=300
delta=0.000001
chi_real=np.zeros(len(qgrid),dtype = 'complex')
chi_imag=np.zeros(len(qgrid),dtype = 'complex')
#chi2=np.zeros(len(qgrid),dtype = 'complex')
for i in np.arange(0,len(qgrid)):
    q_evals=eigenvalue(kgrid+qgrid[i,:],hr)
    for m in band_cross_fermi1:
        for n in band_cross_fermi2:
            chi_real[i]=chi_real[i]+np.sum((fermi_dirac(k_evals[:,m],efermi,T)-fermi_dirac(q_evals[:,n],efermi,T))/(k_evals[:,m]-q_evals[:,n]+1j*delta))
            #chi2[i]=chi2[i]+np.sum((fermi_dirac(k_evals[:,m],efermi,T))*(1-fermi_dirac(q_evals[:,n],efermi,T))/(q_evals[:,n]-k_evals[:,m]+1j*delta))
            chi_imag[i]=chi_imag[i]+np.sum(np.multiply(delta_function(k_evals[:,m],efermi),delta_function(q_evals[:,n],efermi)))
    print(i)

   
chi_real=np.real(-np.real(chi_real))/(nk[1]*nk[1])
chi_imag=np.real(np.real(chi_imag))/(nk[1]*nk[1])
#chi2=np.real(chi2)

qgrid_real=point_scale(qgrid,hr['reciprocal'])
X=qgrid_real[:,0].reshape(nk[0],nk[1])
Y=qgrid_real[:,1].reshape(nk[0],nk[1])
chi_real=chi_real.reshape(nk[0],nk[1])
chi_imag=chi_imag.reshape(nk[0],nk[1])
#chi2=chi2.reshape(nk[0],nk[1])

plt.figure()
plt.contourf(X,Y,chi_real,70)
plt.axis('equal')
plt.title('real')
plt.colorbar()
fig=plt.gcf()
cm=8.5
inch=0.3937007874
fig.set_size_inches(cm*inch,cm*inch*1)
plt.savefig('real.png',format='png',dpi=600,bbox_inches='tight')

plt.figure()
plt.contourf(X,Y,chi_imag,70)
plt.axis('equal')
plt.title('imag')
plt.colorbar()
fig=plt.gcf()
cm=8.5
inch=0.3937007874
fig.set_size_inches(cm*inch,cm*inch*1)
plt.savefig('imag.png',format='png',dpi=600,bbox_inches='tight')

fig=plt.figure()
ax = Axes3D(fig)     
ax.plot_surface(X,Y,chi_real,cmap=plt.cm.bwr)

plt.show()
        
        
np.savetxt('chi_real.dat',chi_real)
np.savetxt('chi_imag.dat',chi_imag)  
np.savetxt('X.dat',X) 
np.savetxt('Y.dat',Y)
    
    

