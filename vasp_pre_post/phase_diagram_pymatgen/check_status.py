#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  5 16:29:52 2024

@author: dj
"""

from pymatgen.io.vasp import Vasprun,Outcar
import os

filenames=os.listdir('./structures')
for filename in filenames:
    prefix=filename.split('.')[0]
    try:
        vasprun=Vasprun('%s/scf/vasprun.xml'%prefix)
        ele_converge=vasprun.converged_electronic
        ionic_converge=vasprun.converged_ionic
        print(prefix,ele_converge,ionic_converge)
    except:
        print(prefix+' error')
