# -*- coding: utf-8 -*-
"""
Created on Wed Nov 16 09:08:47 2022

@author: DengJun
"""


import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib import rcParams

rcParams['mathtext.fontset'] = 'custom'
rcParams['mathtext.rm'] = 'Times New Roman'
rcParams['mathtext.it'] = 'Times New Roman:italic'
rcParams['mathtext.bf'] = 'Times New Roman:bold'

font={'family':'Times New Roman', 'weight':'normal','size':8}

g_005=np.loadtxt('0.0050')
g_010=np.loadtxt('0.0100')
g_015=np.loadtxt('0.0150')
g_020=np.loadtxt('0.0200')
g_025=np.loadtxt('0.0250')
g_030=np.loadtxt('0.0300')
g_035=np.loadtxt('0.0350')
g_040=np.loadtxt('0.0400')
g_045=np.loadtxt('0.0450')
g_050=np.loadtxt('0.0500')
phonon=np.loadtxt('RhSeCl_444.freq.gp')
high_sym=np.array([0, 0.5, 1, 1.7071, 1.9288, 2.4288, 2.9288, 3.6359])
for i in high_sym:
    plt.plot([i,i],[-0.01,0.01],color='grey',linewidth=0.5)

i=5
plt.plot(phonon[:,0],g_005[:,i],label='0.005')
plt.plot(phonon[:,0],g_010[:,i],label='0.010')
plt.plot(phonon[:,0],g_015[:,i],label='0.015')
plt.plot(phonon[:,0],g_020[:,i],label='0.020')
plt.plot(phonon[:,0],g_025[:,i],label='0.025')
plt.plot(phonon[:,0],g_030[:,i],label='0.030')
plt.plot(phonon[:,0],g_035[:,i],label='0.035')
plt.plot(phonon[:,0],g_040[:,i],label='0.040')
plt.plot(phonon[:,0],g_045[:,i],label='0.045')
plt.plot(phonon[:,0],g_050[:,i],label='0.050')
plt.legend()

label=["$\Gamma$",'X','S','Y','$\Gamma$','Z','U','R','T','Z']

plt.xticks([])
d=0.02
font2={'family':'Times New Roman', 'weight':'normal','size':6}
plt.text(0-0.003, -d,label[0],font2)
plt.text(high_sym[1]-0.01, -d,label[1],font2)
plt.text(high_sym[2]-0.01,-d,label[2],font2)
plt.text(high_sym[3]-0.01,-d,label[3],font2)
plt.text(high_sym[4]-0.01,-d,label[4],font2)
plt.text(high_sym[5]-0.01,-d,label[5],font2)
plt.text(high_sym[6]-0.02,-d,label[6],font2)
plt.text(high_sym[7]-0.03,-d,label[7],font2)
