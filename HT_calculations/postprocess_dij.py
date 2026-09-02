#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 14 09:33:41 2022

@author: dj
"""

from pymatgen.io.vasp.outputs import Outcar
from pymatgen.core.structure import Structure
import numpy as np
import pymongo
import copy
import os

client=pymongo.MongoClient("mongodb://localhost:27017")
db=client['layers'] 
col=db['stacking_ionic']
query={'vdw':True}
datas=col.find(query)
dij_tot=[]
for data in datas:
    pie=np.array(data['d2']['piezoelectric'])
    elastic=np.array(data['d2']['elastic'])
    s=Structure.from_dict(data['d2']['structure_opt'])
    dij=np.dot(pie,np.linalg.inv(elastic))
    query={'structure_id':data['structure_id']}
    new={'$set':{'d2.dij':dij.tolist()}}
    col.update_one(query,new)
    dij_tot.append(dij)
    print(data['structure_id'],s.get_space_group_info(),abs(dij[1,1]),abs(dij[2,0]))
