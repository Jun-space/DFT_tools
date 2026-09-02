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



efermi= 3.39279586 #fermi energy
hr=read_hr('wannier90_hr.dat')
band_cross_fermi=[120,121,122,123,124,125,126,127]
# band_cross_fermi2=[60,61,62,63]

nk=[50,50,1]
#kgrid=gen_kpt_grid(50,50,32)
kgrid=gen_kpt_grid(nk[0],nk[1],nk[2])
#qgrid=np.array(kgrid)
qpath=np.loadtxt('qpath.dat')

k_evals=eigenvalue(kgrid,hr)

T=300
delta=0.000001
chi=np.zeros(len(qpath),dtype = 'complex')
chi2=np.zeros(len(qpath),dtype = 'complex')
for i in np.arange(0,len(qpath)):
    q_evals=eigenvalue(kgrid+qpath[i,:],hr)
    for m in band_cross_fermi:
        for n in band_cross_fermi:
            chi[i]=chi[i]+np.sum((fermi_dirac(k_evals[:,m],efermi,T)-fermi_dirac(q_evals[:,n],efermi,T))/(k_evals[:,m]-q_evals[:,n]+1j*delta))
            #chi2[i]=chi2[i]+np.sum((fermi_dirac(k_evals[:,m],efermi,T))*(1-fermi_dirac(q_evals[:,n],efermi,T))/(q_evals[:,n]-k_evals[:,m]+1j*delta))
    print(i)      
chi=np.real(-np.real(chi))
#chi2=np.real(chi2)

qpath_real=point_scale(qpath,hr['reciprocal'])
x=np.zeros(len(qpath))
for i in range(len(qpath)-1):
    x[i+1]=x[i]+np.linalg.norm(qpath_real[i+1]-qpath_real[i])


plt.figure()

font={'family':'Times New Roman','weight':'normal','size':6}   
label=["$\Gamma$",'X','M','$\Gamma$','Z','R','A','Z']
# d=33000
# plt.xticks([])
# plt.text(0-0.008, d,label[0],font)
# plt.text(high_sym[1]-0.02, d,label[1],font)
# plt.text(high_sym[2]-0.02,d,label[2],font)
# plt.text(high_sym[3]-0.02,d,label[3],font)


# high_sym=np.array([0,x[100],x[361],x[461],x[722],x[914],x[1015],x[1276],x[1376],x[1636]])

# d=0.7
# plt.text(high_sym[0]-0.005, d,label[0],font)
# plt.text(high_sym[1]-0.02, d,label[1],font)
# plt.text(high_sym[2]-0.01, d,label[2],font) 
# plt.text(high_sym[3]-0.015, d,label[3],font)
# plt.text(high_sym[4]-0.01, d,label[4],font)
# plt.text(high_sym[5]-0.01, d,label[5],font)
# plt.text(high_sym[6]-0.01, d,label[6],font)
# plt.text(high_sym[7]-0.02, d,label[7],font)
# plt.text(high_sym[8]-0.02, d,label[7],font)
# plt.text(high_sym[9]-0.02, d,label[7],font)

#  # #plot high symmetry lines
# for i in np.arange(1,10):
#     plt.plot([high_sym[i],high_sym[i]],[0.5,3],'--',color='k',linewidth=0.5)
   

#set ylabels
plt.ylabel('$\chi$ (a.u.)',{'family':'Times New Roman','weight':'normal','size':8})  

plt.plot(x,chi/(nk[0]*nk[1]*nk[2]),linewidth=0.5)
plt.xlim(0,x[-1])
# plt.axis([0,x[-1],0.8,2.4])

ax=plt.gca()
ax.minorticks_on()
ax.yaxis.set_minor_locator(plt.MultipleLocator(50000))
ax.tick_params(axis='y',direction='in',which="major",width=0.5,length=3)
ax.tick_params(axis='y',direction='in',which="minor",width=0.5,length=1.5)
ax.tick_params(labelsize=6)
labels=ax.get_xticklabels()+ax.get_yticklabels()
[label.set_fontname('Times New Roman') for label in labels]

#set line width
ax.spines['bottom'].set_linewidth(0.5)
ax.spines['left'].set_linewidth(0.5)
ax.spines['right'].set_linewidth(0.5)
ax.spines['top'].set_linewidth(0.5)

plt.xticks([])

fig=plt.gcf()
cm=8.5
inch=0.3937007874
fig.set_size_inches(cm*inch,6*inch*1)
plt.savefig('real_path.png',format='png',dpi=600,bbox_inches='tight')

x2=np.reshape(x,[len(x),1])
chi2=np.reshape(chi/(nk[0]*nk[1]*nk[2]),[len(chi2),1])
data=np.concatenate((x2,chi2),axis=1)
np.savetxt('chi_path_50_50_1.dat',data,fmt='%0.4f')