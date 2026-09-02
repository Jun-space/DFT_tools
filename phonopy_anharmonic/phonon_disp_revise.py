# -*- coding: utf-8 -*-
"""
Created on Thu Mar 12 09:30:13 2026

@author: dengjun
"""

import yaml
from pymatgen.io.vasp import Vasprun,Outcar
import os

n_disp=int(input('please input total number of displacements: '))

filenames=[]
for i in range(1, n_disp+1):
    if i < 1000:
        filenames.append(str(i).zfill(3))
    else:
        filenames.append(i)

converged=[]
for filename in filenames:
    try:
        vasprun=Vasprun('%s/vasprun.xml'%filename)
        ele_converge=vasprun.converged_electronic
        ionic_converge=vasprun.converged_ionic
        #print(filename,ele_converge,ionic_converge)
        if ele_converge == True:
            # print('%s/vasprun.xml '%filename)
            converged.append(int(filename))
    except:
        continue

f=open('phonopy_disp.yaml')
data=yaml.safe_load(f)
f.close()


tmp = [data['dataset']['displacements'][i-1] for i in converged]
data['dataset']['displacements'] = tmp
# tmp=[]
# for i in numbers:
#     tmp.append(data['dataset']['displacements'][i])


with open('phonopy_disp_reshape.yaml', 'w', encoding='utf-8') as f:
    yaml.dump(data, f, 
              allow_unicode=True,      # 支持中文
              default_flow_style=False, # 使用块样式（更易读）
              sort_keys=True)           # 保持字典顺序
    
    

f=open('phonopy_disp_reshape.yaml')
data2=yaml.safe_load(f)
f.close()
