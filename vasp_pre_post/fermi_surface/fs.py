# -*- coding: utf-8 -*-
"""
Created on Tue Feb 14 12:33:57 2023

@author: DengJun
"""


import numpy as np
from scipy.interpolate import griddata
import matplotlib.pyplot as plt
from matplotlib import rcParams
from pylab import *
from matplotlib.colors import ListedColormap,LinearSegmentedColormap

rcParams['mathtext.fontset'] = 'custom'
rcParams['mathtext.rm'] = 'Times New Roman'
rcParams['mathtext.it'] = 'Times New Roman:italic'
rcParams['mathtext.bf'] = 'Times New Roman:bold'
rcParams['pdf.fonttype'] = 42
rcParams['ps.fonttype'] = 42

wr = LinearSegmentedColormap.from_list('mycmap1', ['white', 'red'])
wb = LinearSegmentedColormap.from_list('mycmap2', ['white', 'blue'])
font={'family':'Arial','weight':'normal','size':8} 

eigval=np.loadtxt('band.dat')
klist=np.loadtxt('klist.dat')
mz=np.loadtxt('mz.dat')
efermi=np.loadtxt('Fermi.dat')

eigval=eigval-efermi
# klist=klist[::2,:]
# mz=mz[::2,:]

shape=2000
X_new=np.zeros([shape,shape])
Y_new=np.zeros([shape,shape])
for i in range(shape):
    for j in range(shape):
        X_new[i,j]=1/shape*i
        Y_new[i,j]=1/shape*j
delta=abs((eigval[:,24]-eigval[:,23])*1000)
band1_new=griddata(klist[:,0:2],eigval[:,23],(X_new,Y_new),method='cubic')
band2_new=griddata(klist[:,0:2],eigval[:,24],(X_new,Y_new),method='cubic')
m1_new=griddata(klist[:,0:2],mz[:,22],(X_new,Y_new),method='cubic')
m2_new=griddata(klist[:,0:2],mz[:,23],(X_new,Y_new),method='cubic')
delta_new=griddata(klist[:,0:2],delta,(X_new,Y_new),method='cubic')
tolerance=0.001
kF=[]
m_kF=[]
delta_kF=[]
for i in range(shape):
    for j in range(shape):
        if abs(band1_new[i,j])< tolerance:
            kF.append([X_new[i,j],Y_new[i,j]])
            m_kF.append(m1_new[i,j])
            delta_kF.append(delta_new[i,j])

for i in range(shape):
    for j in range(shape):
        if abs(band2_new[i,j])< tolerance:
            kF.append([X_new[i,j],Y_new[i,j]]) 
            m_kF.append(m2_new[i,j])
            delta_kF.append(delta_new[i,j])

kF=np.array(kF)
rec_lattice=np.array([[0.287159865,0.165791825],[0.000000000,0.331583651]])*2*np.pi
kF=np.dot(kF,rec_lattice)

up=[]
dw=[]
for num,i in enumerate(m_kF):
    if i < 0:
        dw.append([kF[num,0],kF[num,1],i,delta_kF[num]])
    else:
        up.append([kF[num,0],kF[num,1],i,delta_kF[num]])

up=np.array(up)
dw=np.array(dw)

fig=plt.figure()
ax=fig.add_subplot(111)
# plt.scatter(kF[:,0],kF[:,1],0.01,c=delta_kF,cmap='bwr')
map1=ax.scatter(up[:,0],up[:,1],s=up[:,3]/100,c=up[:,3],cmap=wr)
map2=ax.scatter(dw[:,0],dw[:,1],s=dw[:,3]/100,c=dw[:,3],cmap=wb)
ax.axis('equal')

cax1=plt.axes([0.75,0.5, 0.006, 0.15])
cax2=plt.axes([0.85,0.5, 0.006, 0.15])
cb=plt.colorbar(map1,aspect=30,pad=0.1,cax=cax1)
cb.set_label(label=r'$\Delta$ (meV)',size=6,fontdict=font)
cb.ax.tick_params(labelsize=6) 
cb.outline.set_visible(True)
cb.outline.set_linewidth(0.5)

cb=plt.colorbar(map2,aspect=30,pad=0.1,cax=cax2)
cb.set_label(label=r'$\Delta$ (meV)',size=6,fontdict=font)
cb.ax.tick_params(labelsize=6) 
cb.outline.set_visible(True)
cb.outline.set_linewidth(0.5)

ax.spines['bottom'].set_visible(False)
ax.spines['left'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)


k=np.dot(np.array([1/3,1/3]),rec_lattice)
k2=np.dot(np.array([2/3,2/3]),rec_lattice)

x0=np.array([0,0])
x1=np.array([rec_lattice[0,0],rec_lattice[0,1]])
x2=np.array([rec_lattice[1,0],rec_lattice[1,1]])
x3=np.array([rec_lattice[0,0],rec_lattice[0,1]+rec_lattice[1,1]])


ax.plot([x0[0],x1[0]],[x0[1],x1[1]],'-',c='k',linewidth=0.5)
ax.plot([x0[0],x2[0]],[x0[1],x2[1]],'-',c='k',linewidth=0.5)
ax.plot([x2[0],x3[0]],[x2[1],x3[1]],'-',c='k',linewidth=0.5)
ax.plot([x1[0],x3[0]],[x1[1],x3[1]],'-',c='k',linewidth=0.5)

# ax.scatter(k[0],k[1],5,c='k')
# ax.scatter(k2[0],k2[1],5,c='k')
ax.tick_params(labelleft=False,labelbottom=False)
ax.set_xticks([])
ax.set_yticks([])


fig=plt.gcf()
cm=15
inch=0.3937007874
fig.set_size_inches(cm*inch,cm*inch)
plt.savefig('fs.png',dpi=600,bbox_inches='tight')
