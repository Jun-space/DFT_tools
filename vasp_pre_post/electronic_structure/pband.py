# -*- coding: utf-8 -*-
"""
Created on Wed Nov 25 11:16:13 2020

@author: DengJun
"""


import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
import re
from matplotlib.gridspec import GridSpec

rcParams['mathtext.fontset'] = 'custom'
rcParams['mathtext.rm'] = 'Arial'
rcParams['mathtext.it'] = 'Arial:italic'
rcParams['mathtext.bf'] = 'Arial:bold'
rcParams['pdf.fonttype'] = 42
rcParams['ps.fonttype'] = 42
rcParams['font.family'] = 'Arial'
rcParams['font.size'] = 8


fig=plt.figure()
widths = [4,2]
gs = GridSpec(1,1)

####################################################################band
band=np.loadtxt(r'band.dat')
efermi=np.loadtxt(r'Fermi.dat')
high_sym=np.loadtxt(r'high_sym.dat')
Sr=np.loadtxt('Sr.dat')
Hg=np.loadtxt('Hg.dat')
Ti=np.loadtxt('Ti.dat')
O=np.loadtxt('O.dat')

f=open(r'ticks.dat','r')
ticks=f.readlines()
ticks=[i[:-1] for i in ticks]
f.close()
band_number=int(len(band[0,:])-1)


ax=fig.add_subplot(gs[0,0])
xlim=[0,high_sym[-1]]
ylim=[-5,5]
ax.axis([xlim[0],xlim[1],ylim[0],ylim[1]])
band_range=[]
for ii in np.arange(1,band_number+1):
    ax.plot(band[:,0],band[:,ii]-efermi,'-',color='k',linewidth=0.5)
    band_range.append([np.min(band[:,ii]-efermi),np.max(band[:,ii]-efermi)])

scale=10
ax.scatter(np.tile(band[:,0],band_number),Hg[:,1]-efermi,Hg[:,2]*scale,c='r',zorder=2,edgecolor='None',label=r'Hg-$s$')
ax.scatter(np.tile(band[:,0],band_number),O[:,1]-efermi,np.sum(O[:,3:6],axis=1)*scale,c='tab:orange',zorder=2,edgecolor='None',label=r'O-$p$',alpha=0.5)
ax.scatter(np.tile(band[:,0],band_number),Ti[:,1]-efermi,np.sum(Ti[:,6:11],axis=1)*scale,c='tab:green',zorder=2,edgecolor='None',label=r'Ti-$d$',alpha=0.5)
ax.scatter(np.tile(band[:,0],band_number),Sr[:,1]-efermi,Sr[:,2]*scale,c='tab:pink',zorder=2,edgecolor='None',label=r'Sr-$s$',alpha=0.5)
ax.scatter(np.tile(band[:,0],band_number),Hg[:,1]-efermi,np.sum(Hg[:,6:11],axis=1)*scale,c='tab:blue',zorder=2,edgecolor='None',label=r'Hg-$d$',alpha=0.5)


font={'family':'Arial','weight':'normal','size':8}    
ax.legend(prop=font,frameon=False,ncol=2,markerscale=2,loc=[0.22,0.5])
ax.set_ylabel('Energy (eV)',{'family':'Arial','weight':'normal','size':8})

#plot high symmetry lines
for i in np.arange(1,len(high_sym)):
    ax.plot([high_sym[i],high_sym[i]],[-8,8],color='grey',linewidth=0.5,zorder=0)

#draw Fermi Energy
x=np.array([0,high_sym[-1]])
y=np.array([0,0])
ax.plot(x,y,'k--',linewidth=0.5,zorder=0)
	
	
	
#set high symmetry label
font={'family':'Arial','weight':'normal','size':8}   
label=ticks
for i in range(len(label)):
    if label[i][0]== "\\":
        label[i]="$"+label[i]+"$"
    if '_' in label[i]:
        label[i]=re.sub(r'_(.)', r'$_\1$', label[i])
ax.set_xticks([])
d=-0.05
high_sym=high_sym/high_sym[-1] 
# for i in np.arange(1,len(ticks)-1):
#     ax.text(high_sym[i]-0.01,d,label[i],{'family':'Arial','weight':'normal','size':8},transform=ax.transAxes)
ax.text(high_sym[0]-0.005,d,label[0],{'family':'Arial','weight':'normal','size':8},transform=ax.transAxes)
ax.text(high_sym[1]-0.01,d,label[1],{'family':'Arial','weight':'normal','size':8},transform=ax.transAxes)
ax.text(high_sym[2]-0.01,d,label[2],{'family':'Arial','weight':'normal','size':8},transform=ax.transAxes)
ax.text(high_sym[3]-0.01,d,label[3],{'family':'Arial','weight':'normal','size':8},transform=ax.transAxes)
ax.text(high_sym[4]-0.01,d,label[4],{'family':'Arial','weight':'normal','size':8},transform=ax.transAxes)
ax.text(high_sym[5]-0.055,d,label[5],{'family':'Arial','weight':'normal','size':8},transform=ax.transAxes)
ax.text(high_sym[6]-0.02,d,label[6],{'family':'Arial','weight':'normal','size':8},transform=ax.transAxes)

#set yticks 
ax.minorticks_on()
ax.yaxis.set_minor_locator(plt.MultipleLocator(1))
ax.tick_params(axis='y',direction='in',which="major",width=0.5,length=3)
ax.tick_params(axis='y',direction='in',which="minor",width=0.5,length=1.5)
ax.tick_params(labelsize=8)
labels=ax.get_xticklabels()+ax.get_yticklabels()
[label.set_fontname('Arial') for label in labels]

#set line width
ax.spines['bottom'].set_linewidth(0.5)
ax.spines['left'].set_linewidth(0.5)
ax.spines['right'].set_linewidth(0.5)
ax.spines['top'].set_linewidth(0.5)
ax.set_yticks(np.arange(-4,4.1,2))


fig=plt.gcf()
cm=8.5
inch=0.3937007874
fig.set_size_inches(cm*inch,cm*inch*1.2)
plt.savefig('pband.png',format='png',dpi=600,bbox_inches='tight')
