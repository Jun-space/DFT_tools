# -*- coding: utf-8 -*-
"""
Created on Wed Nov 25 11:16:13 2020

@author: DengJun
"""


import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
import re

rcParams['mathtext.fontset'] = 'custom'
rcParams['mathtext.rm'] = 'Arial'
rcParams['mathtext.it'] = 'Arial:italic'
rcParams['mathtext.bf'] = 'Arial:bold'

band_wannier=np.loadtxt(r'wannier90_band.dat')
band=np.loadtxt('../band_PBE/band.dat')
efermi=np.loadtxt(r'../band_PBE/Fermi.dat')
high_sym=np.loadtxt(r'../band_PBE/high_sym.dat')
f=open(r'../band_PBE/ticks.dat','r')
ticks=f.readlines()
ticks=[i[:-1] for i in ticks]
f.close()
band_number=int(len(band[0,:])-1)

# band2=band[90:,]
# band2[:,0]=band2[:,0]-band2[0,0]
# band=band2
# high_sym=high_sym[3:,]
# high_sym=high_sym-high_sym[0]
high_sym=high_sym/high_sym[-1]
fig=plt.figure()

# efermi = 0

band[:,0]=band[:,0]/band[-1,0]
ax=fig.add_subplot(1,1,1)
xlim=[0,1]
ylim=[-3,3]
# ax.axis([xlim[0],xlim[1],ylim[0],ylim[1]])
ax.axis([0,high_sym[-1],ylim[0],ylim[1]])
band_range=[]
for ii in np.arange(1,band_number+1):
    ax.plot(band[:,0],band[:,ii]-efermi,'-',color='b',linewidth=0.5)
    band_range.append([np.min(band[:,ii]-efermi),np.max(band[:,ii]-efermi)])
ax.plot(band[:,0],band[:,ii]-efermi,'-',color='b',linewidth=0.5,label='DFT')
font={'family':'Arial','weight':'normal','size':6}    

ax.set_ylabel('Energy (eV)',{'family':'Arial','weight':'normal','size':8})

#plot high symmetry lines
for i in np.arange(0,len(high_sym)):
    ax.plot([high_sym[i],high_sym[i]],[-8,8],color='grey',linewidth=0.5,zorder=0)

#set high symmetry label
font={'family':'Arial','weight':'normal','size':6}   
label=ticks
for i in range(len(label)):
    if label[i][0]== "\\":
        label[i]="$"+label[i]+"$"
    if '_' in label[i]:
        label[i]=re.sub(r'_(.)', r'$_\1$', label[i])
# label=['$\Gamma$','M','K','$\Gamma$']
ax.set_xticks([])
d=-0.04
# high_sym=[0,0.366024,0.577353,1]
high_sym=high_sym/high_sym[-1]
for i in np.arange(1,len(label)-1):
    ax.text(high_sym[i]-0.01,d,label[i],{'family':'Arial','weight':'normal','size':6},transform=ax.transAxes)
ax.text(high_sym[0]-0.005,d,label[0],{'family':'Arial','weight':'normal','size':6},transform=ax.transAxes)
ax.text(high_sym[-1]-0.02,d,label[-1],{'family':'Arial','weight':'normal','size':6},transform=ax.transAxes)

#set yticks
ax.minorticks_on()
ax.yaxis.set_minor_locator(plt.MultipleLocator(0.5))
ax.tick_params(axis='y',direction='in',which="major",width=0.5,length=3)
ax.tick_params(axis='y',direction='in',which="minor",width=0.5,length=1.5)
ax.tick_params(labelsize=6)
labels=ax.get_xticklabels()+ax.get_yticklabels()
[label.set_fontname('Arial') for label in labels]

#set line width
ax.spines['bottom'].set_linewidth(0.5)
ax.spines['left'].set_linewidth(0.5)
ax.spines['right'].set_linewidth(0.5)
ax.spines['top'].set_linewidth(0.5)
ax.set_yticks(np.arange(ylim[0],ylim[1]+0.1,0.5))

#draw Fermi Energy
x=np.array([0,high_sym[-1]])
y=np.array([0,0])
ax.plot(x,y,'k--',linewidth=0.5,zorder=0)



######plot wannier band############
kpoints=[]
points=len(band_wannier[:,0])
for i in range(points):
    if band_wannier[i+1,0]==0:
        kpoints=i+1
        break

band=band_wannier
points=len(band[:,0])
band_number=int(points/(kpoints))
# efermi=-1.65418943-0.03601016999999995
# efermi = 3.9248023
x1=band[:,0]/np.max(band[:,0])
y1=band[:,1]-efermi
band_energy_up=[]
for ii in range(band_number):
    plt.plot(x1[(ii*kpoints):(ii*kpoints+kpoints)],y1[(ii*kpoints):(ii*kpoints+kpoints)],'--',color='r',linewidth=0.5)
    band_energy_up.append([np.min(y1[(ii*kpoints):(ii*kpoints+kpoints)]),np.max(y1[(ii*kpoints):(ii*kpoints+kpoints)])])
# plt.plot(x1[(ii*kpoints):(ii*kpoints+kpoints)],y1[(ii*kpoints):(ii*kpoints+kpoints)],'--',color='r',linewidth=0.5,label='wannier')
# plt.plot(x1[(ii*kpoints):(ii*kpoints+kpoints)],y1[(ii*kpoints):(ii*kpoints+kpoints)],'--',color='r',linewidth=0.5,label='wannier')
ii=36
plt.plot(x1[(ii*kpoints):(ii*kpoints+kpoints)],y1[(ii*kpoints):(ii*kpoints+kpoints)],'-',color='r',linewidth=0.5,label='wannier')
print(np.max(y1[(ii*kpoints):(ii*kpoints+kpoints)]))
ax.legend(prop=font)
######end wannier band#############  


#this is what you need to specify



fig=plt.gcf()
cm=8.5
inch=0.3937007874
fig.set_size_inches(cm*inch,5*inch*1.5)
plt.savefig('wannier_compare.png',format='png',dpi=600,bbox_inches='tight')
