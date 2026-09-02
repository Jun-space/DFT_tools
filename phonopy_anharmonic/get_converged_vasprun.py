#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  5 16:29:52 2024

@author: dj
"""

from pymatgen.io.vasp import Vasprun,Outcar
import os

filenames=[]
for i in range(1, 1001):
    if i < 1000:
        filenames.append(str(i).zfill(3))
    else:
        filenames.append(i)

for filename in filenames:
    try:
        vasprun=Vasprun('%s/vasprun.xml'%filename)
        ele_converge=vasprun.converged_electronic
        ionic_converge=vasprun.converged_ionic
        #print(filename,ele_converge,ionic_converge)
        if ele_converge == True:
            print('%s/vasprun.xml '%filename)
    except:
        continue
