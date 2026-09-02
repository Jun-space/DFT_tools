# -*- coding: utf-8 -*-
"""
Created on Wed Nov 25 11:16:13 2020

@author: DengJun
"""


import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib.collections import LineCollection

rcParams['mathtext.fontset'] = 'custom'
rcParams['mathtext.rm'] = 'Arial'
rcParams['mathtext.it'] = 'Arial:italic'
rcParams['mathtext.bf'] = 'Arial:bold'

mz=np.loadtxt('mz.dat')
# mx=np.loadtxt('mx.dat')
# my=np.loadtxt('my.dat')
band=np.loadtxt('band.dat')
efermi=np.loadtxt(r'Fermi.dat')
high_sym=np.loadtxt(r'high_sym.dat')
f=open(r'ticks.dat','r')
ticks=f.readlines()
ticks=[i[:-1] for i in ticks]
f.close()
band_number=int(len(band[0,:])-1)

fig=plt.figure()

ax=fig.add_subplot(1,1,1)
xlim=[0,high_sym[-1]]
ylim=[-0.5,1]
ax.axis([xlim[0],xlim[1],ylim[0],ylim[1]])
band_range=[]
norm = plt.Normalize(-1, 1)
for ii in np.arange(1,band_number+1):
    # ax.plot(band[:,0],band[:,ii]-efermi,'-',color='k',linewidth=1,zorder=1)
    x = band[:,0]
    y = band[:,ii]-efermi
    z = mz[:,ii-1]
    points = np.array([x, y]).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)
    lc = LineCollection(segments, cmap='coolwarm', norm=norm)
    lc.set_array(z)
    lc.set_linewidth(0.5)
    line = ax.add_collection(lc)
    # q=ax.scatter(band[:,0],band[:,ii]-efermi,3,c=mz[:,ii-1],cmap='bwr',vmin=-1,vmax=1)
    # band_range.append([np.min(band[:,ii]-efermi),np.max(band[:,ii]-efermi)])

font={'family':'Arial','weight':'normal','size':6}    
cb=plt.colorbar(line,aspect=40)
cb.set_label(label=r'm$_y$',size=6,fontdict=font)
cb.ax.tick_params(labelsize=6,width=0.5) 
cb.set_ticks(np.array([-1,-0.5,0,0.5,1]))
cb.outline.set_visible(True)
cb.outline.set_linewidth(0.5)
ax.set_ylabel('Energy (eV)',{'family':'Arial','weight':'normal','size':8})

#plot high symmetry lines
for i in np.arange(1,len(high_sym)):
    ax.plot([high_sym[i],high_sym[i]],[-8,8],color='grey',linewidth=0.5,zorder=0)

#set high symmetry label
font={'family':'Arial','weight':'normal','size':6}   
label=ticks
for i in range(len(label)):
    if label[i][0]== "\\":
        label[i]="$"+label[i]+"$"
ax.set_xticks([])
d=-0.05
high_sym=high_sym/high_sym[-1] 
for i in np.arange(1,len(ticks)-1):
    ax.text(high_sym[i]-0.01,d,label[i],{'family':'Arial','weight':'normal','size':6},transform=ax.transAxes)
ax.text(high_sym[0]-0.005,d,label[0],{'family':'Arial','weight':'normal','size':6},transform=ax.transAxes)
ax.text(high_sym[-1]-0.02,d,label[-1],{'family':'Arial','weight':'normal','size':6},transform=ax.transAxes)

#set yticks
ax.minorticks_on()
ax.yaxis.set_minor_locator(plt.MultipleLocator(0.25))
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

###################################################band  no soc
# band_wannier=np.loadtxt(r'wannier90_band.dat') 
# kpoints=[]
# points=len(band_wannier[:,0])
# for i in range(points):
#     if band_wannier[i+1,0]==0:
#         kpoints=i+1
#         break
# points=len(band_wannier[:,0])
# band_number=int(points/(kpoints))
# ratio=np.max(band[:,0])/np.max(band_wannier[:,0])
# x1=band_wannier[:,0]*ratio
# y1=band_wannier[:,1]-efermi
# band_energy_up=[]
# for ii in range(band_number):
#     plt.plot(x1[(ii*kpoints):(ii*kpoints+kpoints)],y1[(ii*kpoints):(ii*kpoints+kpoints)],'--',color='k',linewidth=0.5)
#     band_energy_up.append([np.min(y1[(ii*kpoints):(ii*kpoints+kpoints)]),np.max(y1[(ii*kpoints):(ii*kpoints+kpoints)])])
# plt.plot(x1[(36*kpoints):(36*kpoints+kpoints)],y1[(36*kpoints):(36*kpoints+kpoints)],'--',color='k',linewidth=0.5,label='wannier')




fig=plt.gcf()
cm=6
inch=0.3937007874
fig.set_size_inches(cm*inch,6*inch)
plt.savefig('mz.pdf',dpi=600,bbox_inches='tight')
