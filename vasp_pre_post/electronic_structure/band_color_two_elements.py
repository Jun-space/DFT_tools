# -*- coding: utf-8 -*-
"""
Created on Wed Oct 30 22:53:42 2024

@author: dengjun
"""

import numpy as np
import matplotlib.pyplot as plt
import pymongo
from matplotlib.colors import ListedColormap,LinearSegmentedColormap
from matplotlib import rcParams
from matplotlib.gridspec import GridSpec
from scipy.interpolate import griddata
from matplotlib.collections import LineCollection

fig=plt.figure()

rcParams['mathtext.fontset'] = 'custom'
rcParams['mathtext.rm'] = 'Arial'
rcParams['mathtext.it'] = 'Arial:italic'
rcParams['mathtext.bf'] = 'Arial:bold'
rcParams['pdf.fonttype'] = 42
rcParams['ps.fonttype'] = 42

font={'family':'Arial','weight':'normal','size':8}

####################s-2179 CaAlSiF band###############
band=np.loadtxt(r'band.dat')
efermi=np.loadtxt(r'Fermi.dat')
high_sym=np.loadtxt(r'high_sym.dat')
In=np.loadtxt(r'In.dat')
S=np.loadtxt(r'S.dat')

f=open(r'ticks.dat','r')
ticks=f.readlines()
ticks=[i[:-1] for i in ticks]
f.close()
band_number=int(len(band[0,:])-1)

ax=fig.add_subplot(111)
xlim=[0,high_sym[-1]]
ylim=[-2,3]
ax.axis([xlim[0],xlim[1],ylim[0],ylim[1]])
band_range=[]
# for ii in np.arange(1,band_number+1):
#     ax.plot(band[:,0],band[:,ii]-efermi,'-',color='b',linewidth=0.5)
#     band_range.append([np.min(band[:,ii]-efermi),np.max(band[:,ii]-efermi)])

In=(In[:,11])/(In[:,11]+S[:,11])
In=In.reshape(len(band[0,:])-1,len(band[:,0])).T
# Rh=Rh[48:,:]
norm = plt.Normalize(0,1)
cmap = LinearSegmentedColormap.from_list('mycmap', ['blue', 'lightgrey', 'red'])
for i in np.arange(1,band_number+1):
    x = band[:,0]
    y = band[:,i]-efermi
    z = In[:,i-1]
    points = np.array([x, y]).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)
    lc = LineCollection(segments, cmap=cmap, norm=norm)
    lc.set_array(z)
    lc.set_linewidth(0.5)
    line = ax.add_collection(lc)

cax=plt.axes([0.95,0.4, 0.01, 0.1])
cb=plt.colorbar(line,aspect=40,pad=0.01,cax=cax)
cb.ax.tick_params(labelsize=6,width=0.5) 
# cb.set_label(label=r'$\Omega_z$',size=6,fontdict=font)
cb.outline.set_visible(False)
# cb.outline.set_linewidth(0.5)
cb.set_ticks([])

ax.text(1.06,0.34,'S',{'family':'Arial','weight':'normal','size':6},transform=ax.transAxes)
ax.text(1.06,0.515,'In',{'family':'Arial','weight':'normal','size':6},transform=ax.transAxes)


#draw Fermi Energy
x=np.array([0,high_sym[-1]])
y=np.array([0,0])
ax.plot(x,y,'k--',linewidth=0.5,zorder=0)

font={'family':'Arial','weight':'normal','size':6}    

ax.set_ylabel('Energy (eV)',{'family':'Arial','weight':'normal','size':8})

#plot high symmetry lines
for i in np.arange(1,len(high_sym)):
    ax.plot([high_sym[i],high_sym[i]],[-8,8],color='grey',linewidth=0.5,zorder=2)

#set high symmetry label
font={'family':'Arial','weight':'normal','size':6}   
label=ticks
for i in range(len(label)):
    if label[i][0]== "\\":
        label[i]="$"+label[i]+"$"
ax.set_xticks([])
d=-0.05
high_sym=high_sym/high_sym[-1] 
for i in np.arange(1,len(ticks)):
    ax.text(high_sym[i]-0.01,d,label[i],{'family':'Arial','weight':'normal','size':6},transform=ax.transAxes)
ax.text(high_sym[0]-0.005,d,label[0],{'family':'Arial','weight':'normal','size':6},transform=ax.transAxes)
# label[-3]='Z'
# ax.text(high_sym[-3]-0.02,d,label[-3],{'family':'Arial','weight':'normal','size':6},transform=ax.transAxes)

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
ax.set_yticks(np.arange(ylim[0],ylim[1]+0.1,1))


fig=plt.gcf()
cm=8.5
inch=0.3937007874
fig.set_size_inches(cm*inch,5*inch*1.5)
plt.savefig('band_color.png',format='png',dpi=600,bbox_inches='tight')