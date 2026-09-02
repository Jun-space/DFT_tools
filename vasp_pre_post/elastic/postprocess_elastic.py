#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Nov 12 15:38:51 2022

@author: dj
"""

from pymatgen.io.vasp.outputs import Outcar
import numpy as np
import pymongo
import copy

def get_elastic_tensor():
    # path=path.replace('jdeng/','dj/C2/')
    outcar=Outcar('OUTCAR')
    outcar.read_elastic_tensor()
    matrix=np.array(outcar.data['elastic_tensor'])/10  #GPa
    # matrix=np.array([[matrix[0,0],matrix[0,1],matrix[0,3]],[matrix[1,0],matrix[1,1],matrix[1,3]],[matrix[3,0],matrix[3,1],matrix[3,3]]])
    matrix2=copy.deepcopy(matrix)
    # matrix=matrix/10*25*10**-10*10**9
    matrix2[:,3]=matrix[:,4]
    matrix2[:,4]=matrix[:,5]
    matrix2[:,5]=matrix[:,3]
    matrix_tmp=copy.deepcopy(matrix2)
    matrix2[3,:]=matrix_tmp[4,:]
    matrix2[4,:]=matrix_tmp[5,:]
    matrix2[5,:]=matrix_tmp[3,:]
    return matrix2

matrix=get_elastic_tensor()
np.savetxt('elastic_tensor.dat',matrix,fmt='%12.5f')
