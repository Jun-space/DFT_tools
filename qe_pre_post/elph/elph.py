# -*- coding: utf-8 -*-
"""
Created on Mon Feb  6 17:00:19 2023

@author: DengJun
"""


import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib import rcParams
from matplotlib.collections import LineCollection

rcParams['mathtext.fontset'] = 'custom'
rcParams['mathtext.rm'] = 'Arial'
rcParams['mathtext.it'] = 'Arial:italic'
rcParams['mathtext.bf'] = 'Arial:bold'
rcParams['pdf.fonttype'] = 42
rcParams['ps.fonttype'] = 42
font={'family':'Arial', 'weight':'normal','size':8}

fig=plt.figure()
widths = [4,2,2]
gs = GridSpec(1,3,wspace=0.2,width_ratios=widths)

phonon=np.loadtxt(r'RhSeCl_444.freq.gp')
dos=np.loadtxt(r'matdyn.dos',comments='#')        
gam=np.loadtxt(r'0.0400')
# a2F=np.loadtxt('alpha2F.dat')

[a,b]=np.shape(gam)
lambda_qv=np.zeros(np.shape(gam))
frequency=phonon[:,1:]
N_ef=7.517070
pi=3.141592653589793
hbar=6.582119569*10**(-16) #eV.s


for i in np.arange(a):
    for j in np.arange(b):
        if frequency[i,j]==0:
            lambda_qv[i,j]=0
        else:
            # lam=gamma/(pi*hbar*N_ef*omega**2*2*pi)
            lambda_qv[i,j]=gam[i,j]*10**12/(pi*hbar*N_ef*(frequency[i,j]/33.356*10**12)**2*2*pi)
            # lambda_qv[i,j]=gam[i,j]

ax=fig.add_subplot(gs[0])
high_sym=np.array([0, 1.3150, 1.4261, 2.2774, 3.2858, 3.8264, 4.6776, 4.9591, 5.8890, 6.1812, 7.3851])
for i in high_sym:
    ax.plot([i,i],[-100,1500],color='grey',linewidth=0.5)

# norm = plt.Normalize(0,np.max(lambda_qv))
# for i in np.arange(1,25):
#     ax.plot(phonon[:,0],phonon[:,i],'k',linewidth=0.5)
#     #plt.scatter(phonon[:,0],phonon[:,i],lambda_qv[:,i-1]*500000+0.0035213584217470223,marker='o',c='r')
#     line=ax.scatter(phonon[:,0],phonon[:,i],lambda_qv[:,i-1]*50,c=lambda_qv[:,i-1],marker='o',cmap='inferno',norm=norm)
# norm = plt.Normalize(0,np.max(lambda_qv))
norm = plt.Normalize(0,0.01)
for i in np.arange(1,len(phonon[0,:])):
    x = phonon[:,0]
    y = phonon[:,i]/33.356
    z = lambda_qv[:,i-1]
    points = np.array([x, y]).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)
    lc = LineCollection(segments, cmap='jet', norm=norm)
    lc.set_array(z)
    lc.set_linewidth(0.5)
    line = ax.add_collection(lc)
    
# plt.colorbar()
# cb=plt.colorbar(line,aspect=30,pad=0.05,fraction=0.02)
# cb.ax.tick_params(labelright=False)
# cb.ax.tick_params(size=0)
# cb.outline.set_visible(True)
# cb.outline.set_linewidth(0.5)
cb=plt.colorbar(line,aspect=30)
cb.set_label(label=r'$\gamma_{q\upsilon}$',size=6,fontdict=font)
cb.ax.tick_params(labelsize=6,width=0.5) 
cb.set_ticks([])
cb.outline.set_visible(True)
cb.outline.set_linewidth(0.5)    
ax.plot([0,phonon[-1,0]],[0,0],'k--',linewidth=0.5)
ax.axis([0,phonon[-1,0],0,15])

label=["$\Gamma$",'X','Y',r'$\Sigma$','$\Gamma$','Z',r'$\Sigma_1$','N','P','Y$_1$','Z']

plt.xticks([])
d=30/33.356
font2={'family':'Arial', 'weight':'normal','size':6}
ax.text(0-0.01, -d,label[0],font2)
ax.text(high_sym[1]-0.08, -d,label[1],font2)
ax.text(high_sym[2]-0.08,-d,label[2],font2)
ax.text(high_sym[3]-0.08,-d,label[3],font2)
ax.text(high_sym[4]-0.08,-d,label[4],font2)
ax.text(high_sym[5]-0.08,-d,label[5],font2)
ax.text(high_sym[6]-0.08,-d,label[6],font2)
ax.text(high_sym[7]-0.1,-d,label[7],font2)
ax.text(high_sym[8]-0.1,-d,label[8],font2)
ax.text(high_sym[9]-0.1,-d,label[9],font2)
ax.text(high_sym[10]-0.1,-d,label[10],font2)


ax.set_ylabel('Frequency (THz)',{'family':'Arial','weight':'normal','size':8}  )

ax.spines['bottom'].set_linewidth(0.5)
ax.spines['left'].set_linewidth(0.5)
ax.spines['right'].set_linewidth(0.5)
ax.spines['top'].set_linewidth(0.5)

#control Y ticks text size
ax.tick_params(labelsize=6)
labels=ax.get_xticklabels()+ax.get_yticklabels()
[label.set_fontname('Arial') for label in labels]

ax.minorticks_on()
ax.yaxis.set_minor_locator(plt.MultipleLocator(100))
ax.tick_params(axis='y',direction='in',which="major",width=0.5,length=3)
ax.tick_params(axis='y',direction='in',which="minor",width=0.5,length=1.5)


###########################dos
ax=fig.add_subplot(gs[1])
Rh=dos[:,2]
Se=dos[:,3]
Cl=dos[:,4]
# F=dos[:,8]+dos[:,9]
color=['gray','hotpink','tab:blue','orange','forestgreen','blueviolet']
y0=np.zeros([len(dos[:,0])])
ax.plot(dos[:,1],dos[:,0]/33.356,'k-',linewidth=0.5,label='total')
plt.fill_betweenx(dos[:,0]/33.356,dos[:,1],y0,color=color[0],alpha=0.3,edgecolor='None')
ax.plot(Rh,dos[:,0]/33.356,'-',c=color[1],linewidth=0.5,label='Rh')
plt.fill_betweenx(dos[:,0]/33.356,Rh,y0,color=color[1],alpha=0.3,edgecolor='None')
ax.plot(Se,dos[:,0]/33.356,'-',c=color[2],linewidth=0.5,label='Se')
plt.fill_betweenx(dos[:,0]/33.356,Se,y0,color=color[2],alpha=0.3,edgecolor='None')
ax.plot(Cl,dos[:,0]/33.356,'-',c=color[3],linewidth=0.5,label='Cl')
plt.fill_betweenx(dos[:,0]/33.356,Cl,y0,color=color[3],alpha=0.3,edgecolor='None')
# ax.plot(F,dos[:,0],'-',c=color[4],linewidth=0.5,label='F')
# plt.fill_betweenx(dos[:,0],F,y0,color=color[4],alpha=0.3)
ax.legend(prop={'family':'Arial', 'weight':'normal','size':6},frameon=False,loc=[0.1,0])
ax.plot([0,2],[0,0],'k--',linewidth=0.5)

ax.spines['bottom'].set_linewidth(0.5)
ax.spines['left'].set_linewidth(0.5)
ax.spines['right'].set_linewidth(0.5)
ax.spines['top'].set_linewidth(0.5)

#control Y ticks text size
plt.tick_params(labelsize=6)
labels=ax.get_xticklabels()+ax.get_yticklabels()
[label.set_fontname('Arial') for label in labels]

ax.minorticks_on()
ax.yaxis.set_minor_locator(plt.MultipleLocator(50))
ax.tick_params(axis='y',direction='in',which="major",width=0.5,length=3)
ax.tick_params(axis='y',direction='in',which="minor",width=0.5,length=1.5)

ax.xaxis.set_minor_locator(plt.MultipleLocator(5))
ax.tick_params(axis='x',direction='in',which="major",width=0.5,length=1.5)
ax.tick_params(axis='x',direction='in',which="minor",width=0.5,length=1.5)


ax.tick_params(labelleft=False)
ax.set_xticks([])
# ax.set_yticks([])
#plt.xlabel('PDOS',{'family' : 'Arial','weight' : 'normal','size'   : 6})
# ax.text(0.1,-20,'PDOS',{'family' : 'Arial','weight' : 'normal','size'   : 6})
ax.set_xlabel('DOS',{'family':'Arial','weight':'normal','size':6}  ,labelpad=0.6)
ax.axis([0,0.06,0,15])

###alpha2F
f=open(r'alpha2F.dat')
lines=f.readlines()
f.close()

l1=[]
l2=[]
l3=[]
l4=[]
for num,line in enumerate(lines[4:]):
    values=[float(x) for x in line[:].split()]
    if num%4 == 0:
        l1.append(values)
    elif num%4 == 1:
        l2.append(values)
    elif num%4 == 2:
        l3.append(values)
    elif num%4 == 3:
        l4.append(values)
l1=np.array(l1)
l2=np.array(l2)
l3=np.array(l3)
l4=np.array(l4)

a2F=np.concatenate((l1,l2,l3,l4),axis=1)

ax=fig.add_subplot(gs[2])
ax.plot(a2F[:,20],a2F[:,0],'b-',linewidth=0.5,label='$\lambda$($\omega$)')
ax.axis([0,1,0,15])

lam=np.zeros((len(a2F[:,0]),1))
for i in np.arange(1,len(a2F)):
    x=a2F[1:i,0]
    y=a2F[1:i,20]
    lam[i]=2*np.trapz(y/x,x)

ax.plot(lam,a2F[:,0],'--',color='purple',linewidth=0.5,label=r'$\alpha^2$F($\omega$)')
ax.legend(prop=font,frameon=False,loc=[0.07,0.0],markerscale=0.6,labelspacing=0.2)

ax.spines['bottom'].set_linewidth(0.5)
ax.spines['left'].set_linewidth(0.5)
ax.spines['right'].set_linewidth(0.5)
ax.spines['top'].set_linewidth(0.5)

#control Y ticks text size
ax.tick_params(labelsize=6)
labels=ax.get_xticklabels()+ax.get_yticklabels()
[label.set_fontname('Arial') for label in labels]

ax.minorticks_on()
ax.yaxis.set_minor_locator(plt.MultipleLocator(50))
ax.tick_params(axis='y',direction='in',which="major",width=0.5,length=3)
ax.tick_params(axis='y',direction='in',which="minor",width=0.5,length=1.5)

ax.xaxis.set_minor_locator(plt.MultipleLocator(0.5))
ax.tick_params(axis='x',direction='in',which="major",width=0.5,length=1.5)
ax.tick_params(axis='x',direction='in',which="minor",width=0.5,length=1.5)
ax.set_yticks(np.arange(0,450/33.356,100))
ax.tick_params(labelleft=False)

fig=plt.gcf()
cm=15
inch=0.3937007874
fig.set_size_inches(cm*inch,8.5*inch*0.8)
plt.savefig('elph.png',format='png',dpi=600,bbox_inches='tight')