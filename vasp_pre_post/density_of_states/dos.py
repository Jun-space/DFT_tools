# -*- coding: utf-8 -*-
"""
Created on Sun May  5 17:17:56 2019

@author: Think
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams

rcParams['mathtext.fontset'] = 'custom'
rcParams['mathtext.rm'] = 'Arial'
rcParams['mathtext.it'] = 'Arial:italic'
rcParams['mathtext.bf'] = 'Arial:bold'

color=['gray','hotpink','tab:blue','orange','forestgreen','blueviolet']

plt.figure()
###################################CrPS4########################################
#get data
tdos=np.loadtxt('tdos.dat',comments='#')
Sr=np.loadtxt('Sr.dat',comments='#')
Hg=np.loadtxt('Hg.dat',comments='#') 
Ti=np.loadtxt('Ti.dat',comments='#')
O=np.loadtxt('O.dat',comments='#') 
efermi=np.loadtxt(r'Fermi.dat')
                
plt.subplot(1,1,1)
y0=np.zeros([len(tdos[:,0])])
#draw dos
plt.plot(tdos[:,0]-efermi,tdos[:,1],color=color[0],linewidth=0.5,label='Total')
plt.fill_between(tdos[:,0]-efermi,tdos[:,1],y0,color=color[0],alpha=0.3,edgecolor='none')


plt.plot(Sr[:,0]-efermi,Sr[:,1],color='tab:pink',linewidth=0.5,label=r'Sr-$s$')
plt.fill_between(Sr[:,0]-efermi,Sr[:,1],y0,color='tab:pink',alpha=0.3,edgecolor='none')

plt.plot(Hg[:,0]-efermi,Hg[:,1],color='r',linewidth=0.5,label=r'Hg-$s$')
plt.fill_between(Hg[:,0]-efermi,Hg[:,1],y0,color='r',alpha=0.3,edgecolor='none')

plt.plot(Ti[:,0]-efermi,np.sum(Ti[:,5:10],axis=1),color='tab:green',linewidth=0.5,label=r'Ti-$d$')
plt.fill_between(Ti[:,0]-efermi,np.sum(Ti[:,5:10],axis=1),y0,color='tab:green',alpha=0.3,edgecolor='none')

plt.plot(O[:,0]-efermi,np.sum(O[:,2:5],axis=1),color='tab:orange',linewidth=0.5,label=r'O-$p$')
plt.fill_between(O[:,0]-efermi,np.sum(O[:,2:5],axis=1),y0,color='tab:orange',alpha=0.3,edgecolor='none')

plt.plot(Hg[:,0]-efermi,np.sum(Hg[:,5:10],axis=1),color='tab:blue',linewidth=0.5,label=r'Hg-$d$')
plt.fill_between(Hg[:,0]-efermi,np.sum(Hg[:,5:10],axis=1),y0,color='tab:blue',alpha=0.3,edgecolor='none')

#draw xline
x=np.array([-25,25])
y=np.array([0,0])
plt.plot(x,y,'k',linewidth=0.5)

#draw Fermi Energy
x1=np.array([0,0])
y1=np.array([-20,100])
plt.plot(x1,y1,'k--',linewidth=0.5)
font={'family':'Arial','weight':'normal','size':8}

#control Y ticks text size


#plt.legend([total,P,Cr,Se],['total','P','Cr','S'],loc=[0.5,0.05],prop={'family':'Times New Roman','size':6},frameon=False)
#plt.legend([Cr_dxy,Cr_dyz,Cr_dxz,Cr_dz2,Cr_dx2],['Cr-d$_{xy}$','Cr-d$_{yz}$','Cr-d$_{xz}$','Cr-d$_{z^2}$','Cr-d$_{x^2}$'],loc=[0.5,0.05],prop={'family':'Times New Roman','size':6},frameon=False)
plt.legend(loc=0,prop={'family':'Arial','size':6},frameon=False,ncol=2)
plt.ylabel('DOS (states/eV)',font)
plt.xlabel('Energy (eV)',font)
ax=plt.gca()
ax.minorticks_on()

ax.yaxis.set_minor_locator(plt.MultipleLocator(10))
ax.tick_params(axis='y',direction='in',which="major",width=0.5,length=3)
ax.tick_params(axis='y',direction='in',which="minor",width=0.5,length=1.5)

ax.xaxis.set_minor_locator(plt.MultipleLocator(1))
ax.tick_params(axis='x',direction='in',which="major",width=0.5,length=3)
ax.tick_params(axis='x',direction='in',which="minor",width=0.5,length=1.5)

ax.spines['bottom'].set_linewidth(0.5)
ax.spines['left'].set_linewidth(0.5)
ax.spines['right'].set_linewidth(0.5)
ax.spines['top'].set_linewidth(0.5)

plt.tick_params(labelsize=8)
labels=ax.get_xticklabels()+ax.get_yticklabels()
[label.set_fontname('Arial') for label in labels]

# ymax=np.max(tdos[])
# plt.yticks(np.arange(-10,10.1,5))
# plt.xticks([])
ax.axis([-6,6,0,50])
ax.set_xticks(np.arange(-6,6.1,2))
# plt.text(-3.78,10,'(a) CrPS$_4$',font)





fig=plt.gcf()
cm=6
inch=0.3937007874
fig.set_size_inches(cm*inch,8.5*inch*0.6*0.8)
plt.savefig('dos.png',dpi=600,bbox_inches='tight')
