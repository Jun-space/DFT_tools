# -*- coding: utf-8 -*-
"""
Created on Tue Jul 22 14:24:11 2025

@author: dengjun
"""

import yaml
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib import rcParams

rcParams['mathtext.fontset'] = 'custom'
rcParams['mathtext.rm'] = 'Arial'
rcParams['mathtext.it'] = 'Arial:italic'
rcParams['mathtext.bf'] = 'Arial:bold'
rcParams['pdf.fonttype'] = 42
rcParams['ps.fonttype'] = 42

phonon=np.loadtxt('phonon.dat')
high_sym=np.loadtxt('high_sym.dat')

f=open('ticks.dat')
ticks=f.readlines()
labels=[i[:-1] for i in ticks]
f.close()

# np.savetxt('phonon.dat',phonon,fmt='%18.8f')
# np.savetxt('high_sym.dat',high_sym,fmt='%18.8f')

# f=open('ticks.dat','w')
# for i in labels:
#     f.write(i+'\n')
# f.close()

fig=plt.figure()
ax=fig.add_subplot(111)

for i in np.arange(1,len(phonon[0,:]),1):
    ax.plot(phonon[:,0],phonon[:,i],'b',linewidth=0.5)
    
for j in high_sym:
    ax.plot([j,j],[-10,100],'grey',linewidth=0.5)
#draw Fermi Energy
x=np.array([0,high_sym[-1]])
y=np.array([0,0])
ax.plot(x,y,'k--',linewidth=0.5)

freq_min=np.min(phonon[:,1:])
freq_max=np.max(phonon[:,1:])
ylim=[freq_min-(freq_max-freq_min)*0.02,freq_max+(freq_max-freq_min)*0.02]
ax.axis([0,high_sym[-1],ylim[0],ylim[1]])
    
ax.set_ylabel('Frequency (THz)',{'family' : 'Arial','weight' : 'normal','size'   : 6})

ax.set_xticks([])
# label=["$\Gamma$",'Y',r'F$_{\rm 0}$|$\Delta_{\rm 0}$',r'$\Gamma$','Z',r'B$_{\rm 0}$|G$_{\rm 0}$','T',r'Y|$\Gamma$','S','R','Z','T']
font = {'family' : 'Arial','weight' : 'normal','size'   : 8}

# d=-7
# ax.text(0-0.003, d,labels[0],font)
# ax.text(high_sym[1]-0.04, d,labels[1],font)
# ax.text(high_sym[2]-0,d,labels[2],font)
# ax.text(high_sym[3]-0.01,d,labels[3],font)
# ax.text(high_sym[4]-0.01,d,labels[4],font)
# ax.text(high_sym[5]-0.01,d,labels[5],font)
# ax.text(high_sym[6]-0.01,d,labels[6],font)
# ax.text(high_sym[7]-0.01,d,labels[7],font)
# ax.text(high_sym[8]-0.01,d,labels[8],font)
# ax.text(high_sym[9]-0.01,d,labels[9],font)
# # ax.text(high_sym[10]-0.01,-9,labels[10],font)
# # ax.text(high_sym[11]-0.01,-9,labels[11],font)

d=-0.04
high_sym=high_sym/high_sym[-1] 
for i in np.arange(1,len(ticks)-1):
    ax.text(high_sym[i]-0.01,d,labels[i],{'family':'Arial','weight':'normal','size':6},transform=ax.transAxes)
ax.text(high_sym[0]-0.005,d,labels[0],{'family':'Arial','weight':'normal','size':6},transform=ax.transAxes)
ax.text(high_sym[-1]-0.02,d,labels[-1],{'family':'Arial','weight':'normal','size':6},transform=ax.transAxes)



ax=plt.gca()

ax.spines['bottom'].set_linewidth(0.5)
ax.spines['left'].set_linewidth(0.5)
ax.spines['right'].set_linewidth(0.5)
ax.spines['top'].set_linewidth(0.5)

#control Y ticks text size
plt.tick_params(labelsize=6)
labels=ax.get_xticklabels()+ax.get_yticklabels()
[label.set_fontname('Arial') for label in labels]

ax.minorticks_on()
ax.yaxis.set_minor_locator(plt.MultipleLocator(5))
ax.tick_params(axis='y',direction='in',which="major",width=0.5,length=3)
ax.tick_params(axis='y',direction='in',which="minor",width=0.5,length=1.5)

# ax.set_title('FE',font)
# ax.text(-0.15,1.03,'(a)',{'family':'Arial','weight':'normal','size':8},transform=ax.transAxes)

fig=plt.gcf()
cm=10
inch=0.3937007874
fig.set_size_inches(cm*inch,8.5*inch*1.2)
plt.savefig('phonon.png',format='png',dpi=600,bbox_inches='tight')
