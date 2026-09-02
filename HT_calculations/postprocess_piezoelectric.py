#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Nov 12 21:51:03 2022

@author: dj
"""

from pymatgen.io.vasp.outputs import Outcar
import numpy as np
import pymongo
import copy
import os

def get_piezoelectric_tensor(path):
    path=path+'/d2/piezoelectric_opp'
    # path=path.replace('jdeng/','dj/C2/')
    # path='/home/dj/C2/Documents/vasp/layer_ionic/htc7/s-5802/d2/piezoelectric_opp'
    pie_ele_tmp=os.popen('grep -A 5 \'PIEZOELECTRIC TENSOR (including local field effects)  for field in x, y, z        (C/m^2)\' %s/OUTCAR | tail -3| awk {\'print $2,$3,$4,$5,$6,$7\'}'%path).readlines()
    
    pie_ele=[]
    for num,i in enumerate(pie_ele_tmp):
        values=[float(j) for j in i.split()]
        pie_ele.append(values)
    pie_ele=np.array(pie_ele)    
    
    pie_ion_tmp=os.popen('grep -A 5 \'PIEZOELECTRIC TENSOR IONIC CONTR  for field in x, y, z        (C/m^2)\' %s/OUTCAR | tail -3| awk {\'print $2,$3,$4,$5,$6,$7\'}'%path).readlines()
    
    pie_ion=[]
    for num,i in enumerate(pie_ion_tmp):
        values=[float(j) for j in i.split()]
        pie_ion.append(values)
    pie_ion=np.array(pie_ion)
    
    #pie_tmp=(pie_ele+pie_ion)*25*10**-10*10**12 #pm/V
	pie_tmp=(pie_ele+pie_ion) #pm/V
    pie=copy.deepcopy(pie_tmp)
    pie[:,3]=pie_tmp[:,4]
    pie[:,4]=pie_tmp[:,5]
    pie[:,5]=pie_tmp[:,3]

    return pie

get_piezoelectric_tensor('.')