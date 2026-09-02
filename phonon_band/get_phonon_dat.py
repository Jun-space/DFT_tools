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

f=open('band.yaml')
data=yaml.safe_load(f)
f.close()

labels=[]
labels.append(data['labels'][0][0])
for num in range(1,len(data['labels'])):
    if data['labels'][num][0] == data['labels'][num-1][1]:
        labels.append(data['labels'][num][0])
    else:
        labels.append(data['labels'][num-1][1]+'|'+data['labels'][num][0])
labels.append(data['labels'][-1][1])

rec_lattice=data['reciprocal_lattice']

freq=[]
distance=[]

for i in data['phonon']:
    distance.append(i['distance'])
    tmp=[]
    for j in i['band']:
        tmp.append(j['frequency'])
    tmp=np.array(tmp)
    freq.append(tmp.tolist())

freq=np.array(freq)
distance=np.array(distance)

high_sym=[]
end=0
high_sym.append(0)
for num,i in enumerate(data['segment_nqpoint']):
    end=end+i
    if end == len(distance):
        high_sym.append(distance[end-1])
    else:
        high_sym.append(distance[end])

high_sym.append(distance[-1])

high_sym=np.array(high_sym)
phonon=np.concatenate((distance.reshape(len(distance),1),freq),axis=1)


np.savetxt('phonon.dat',phonon,fmt='%18.8f')
np.savetxt('high_sym.dat',high_sym,fmt='%18.8f')

f=open('ticks.dat','w')
for i in labels:
    f.write(i+'\n')
f.close()

