# -*- coding: utf-8 -*-
"""
Created on Wed Nov 25 11:16:13 2020

@author: DengJun
"""


import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams

rcParams['mathtext.fontset'] = 'custom'
rcParams['mathtext.rm'] = 'Arial'
rcParams['mathtext.it'] = 'Arial:italic'
rcParams['mathtext.bf'] = 'Arial:bold'


band_up=np.loadtxt('band_up.dat')
band_dw=np.loadtxt('band_dw.dat')
efermi=np.loadtxt(r'Fermi.dat')
high_sym=np.loadtxt(r'high_sym.dat')
f=open(r'ticks.dat','r')
ticks=f.readlines()
ticks=[i[:-1] for i in ticks]
f.close()
band_number=int(len(band_up[0,:])-1)

fig=plt.figure()

ax=fig.add_subplot(1,1,1)
xlim=[0,high_sym[-1]]
ylim=[-2,2]
ax.axis([xlim[0],xlim[1],ylim[0],ylim[1]])
x1=band_up[:,0]
band_range=[]
for ii in np.arange(1,band_number+1):
    ax.plot(band_up[:,0],band_up[:,ii]-efermi,'-',color='r',linewidth=0.5)
    band_range.append([np.min(band_up[:,ii]-efermi),np.max(band_up[:,ii]-efermi)])

for ii in np.arange(1,band_number+1):
    ax.plot(band_dw[:,0],band_dw[:,ii]-efermi,'-',color='b',linewidth=0.5)

#draw Fermi Energy
x=np.array([0,high_sym[-1]])
y=np.array([0,0])
ax.plot(x,y,'k--',linewidth=0.5,zorder=0)    

font={'family':'Arial','weight':'normal','size':6}    

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
d=-0.04
high_sym=high_sym/high_sym[-1] 
for i in np.arange(1,len(ticks)-1):
    ax.text(high_sym[i],d,label[i],{'family':'Arial','weight':'normal','size':6},horizontalalignment='center',transform=ax.transAxes)
ax.text(high_sym[0],d,label[0],{'family':'Arial','weight':'normal','size':6},horizontalalignment='left',transform=ax.transAxes)
ax.text(high_sym[-1],d,label[-1],{'family':'Arial','weight':'normal','size':6},horizontalalignment='right',transform=ax.transAxes)

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


#this is what you need to specify

fig=plt.gcf()
cm=5
inch=0.3937007874
fig.set_size_inches(cm*inch,cm*inch*1.5)
plt.savefig('band.png',format='png',dpi=600,bbox_inches='tight')
