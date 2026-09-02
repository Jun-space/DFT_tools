# -*- coding: utf-8 -*-
"""
Created on Wed Nov 25 17:53:16 2020

@author: DengJun
"""


import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams

rcParams['mathtext.fontset'] = 'custom'
rcParams['mathtext.rm'] = 'Times New Roman'
rcParams['mathtext.it'] = 'Times New Roman:italic'
rcParams['mathtext.bf'] = 'Times New Roman:bold'


f=open('fatband.projwfc_up')
lines=f.readlines()
f.close()

[orbital,kpoints,nbands]=(int(x) for x in lines[12].split())

data=[]
for i in np.arange(0,orbital):
    tmp=[]
    for line in lines[15+i*(kpoints*nbands+1):14+(i+1)*(kpoints*nbands+1)]:
        tmp.append([float(x) for x in line.split()])
    tmp=np.array(tmp)
    tmp2=tmp[:,2]
    tmp3=np.reshape(tmp2,(kpoints,nbands))
    tmp4=tmp3.T
    tmp5=np.reshape(tmp4,(kpoints*nbands))
    data=np.concatenate((data,tmp5))
data=np.reshape(data,(orbital,nbands*kpoints)).T

Te_p=np.sum(data[:,2:8],axis=1)
Ir_d=np.sum(data[:,24:35],axis=1)


d=-4.3
fig=plt.figure()
#this is what you need to specify
band_dat=np.loadtxt('IrTe2_band.out.gnu')
efermi=9.5268

#this is what you need to specify
intersection=30#high symmetry distance
number=6#high symmetry number

kpoints=[]
points=len(band_dat[:,0])
for i in range(points):
    if band_dat[i+1,0]==0:
        kpoints=i+1
        break

band=band_dat[0:int(len(band_dat[:,1])-1),:]
points=len(band[:,0])
band_number=int(points/(kpoints+1))

ax1=fig.add_subplot(1,1,1)
ax1.axis([0,3.5280,-4,4])
#plt.axis([0,soc[kpoints-1,0],-2.5,2.5])
x1=band[:,0]
y1=band[:,1]-efermi
for ii in range(band_number):
    ax1.plot(x1[(ii*kpoints):(ii*kpoints+kpoints)],y1[(ii*kpoints):(ii*kpoints+kpoints)],'-',color='k',linewidth=0.5)
    
ax1.scatter(x1,y1,Ir_d)
    
font={'family':'Times New Roman',
      'weight':'normal',
      'size':6}    
 
high_sym=np.array([0,0.5774,0.9107,1.5773,1.9507,2.5281,2.8614,3.5280])

#plot high symmetry lines
for i in np.arange(1,6):
    ax1.plot([high_sym[i],high_sym[i]],[-4,4],color='grey',linewidth=0.5)
    
#set high symmetry label
font={'family':'Times New Roman','weight':'normal','size':6}   
label=["$\Gamma$",'M','K','$\Gamma$','A','L','H','A']
ax1.set_xticks([])
ax1.text(0-0.008, d,label[0],font)
ax1.text(high_sym[1]-0.02, d,label[1],font)
ax1.text(high_sym[2]-0.02,d,label[2],font)
ax1.text(high_sym[3]+0.006,d,label[3],font)
ax1.text(high_sym[4]-0.015,d,label[4],font)
ax1.text(high_sym[5]-0.015,d,label[5],font)
ax1.text(high_sym[6]-0.015,d,label[6],font)
ax1.text(high_sym[7]-0.02,d,label[7],font)     

#set ylabels
ax1.set_ylabel('Energy (eV)',{'family':'Times New Roman','weight':'normal','size':8})
#ax1.set_title('AFM',{'family':'Times New Roman','weight':'normal','size':8})
#set yticks
ax1.minorticks_on()
ax1.yaxis.set_minor_locator(plt.MultipleLocator(0.5))
ax1.tick_params(axis='y',direction='in',which="major",width=0.5,length=3)
ax1.tick_params(axis='y',direction='in',which="minor",width=0.5,length=1.5)
ax1.tick_params(labelsize=6)
labels=ax1.get_xticklabels()+ax1.get_yticklabels()
[label.set_fontname('Times New Roman') for label in labels]

#set line width
ax1.spines['bottom'].set_linewidth(0.5)
ax1.spines['left'].set_linewidth(0.5)
ax1.spines['right'].set_linewidth(0.5)
ax1.spines['top'].set_linewidth(0.5)



#draw Fermi Energy
x=np.array([0,high_sym[7]])
y=np.array([0,0])
ax1.plot(x,y,'k--',linewidth=0.5)

fig=plt.gcf()
cm=5
inch=0.3937007874
fig.set_size_inches(cm*inch,6*inch*1.2)
plt.savefig('Ir_d_pband.png',format='png',dpi=600,bbox_inches='tight')