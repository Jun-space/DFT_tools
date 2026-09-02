# -*- coding: utf-8 -*-
"""
Created on Sat Feb 12 21:26:15 2022

@author: DengJun
"""


from find_clusters import *
from pymatgen.core.structure import Structure
from dimension_clusters import *
import numpy as np
from output_clusters import *
from delete_am_aem import delete_am_aem
#from valence_clusters import valence_clusters
from remove_duplicate import remove_duplicate_by_Structure
from lattice_clusters import lattice_clusters
from parent_info import *
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
import os
from multiprocessing import Pool
import pymongo
import time
import json
import copy
from bson import json_util

def clusters(filename):
    s=Structure.from_file(filename)
    s=SpacegroupAnalyzer(s,symprec=0.1).get_conventional_standard_structure()
    #try:
    #    valence=[i.specie.oxi_state for i in s.sites]
    #except:
    #    s.add_oxidation_state_by_guess()
    s_delete=s.copy() 
    # s_delete=delete_am_aem(s_delete)#delete alkali and alkali earth metal in the structure
    s_delete.make_supercell([2,2,2])
    clusters_series=find_clusters(s_delete,'covalent')#find clusters in the structure and return the atomic serials
    dimensions=dimensions_by_series(s_delete,clusters_series,2,'covalent') # judge the dimension of the clusters
    structures=output_clusters(s_delete,clusters_series) # convert cluster to pymatgen structural format

    d3=[]
    d2=[]
    d1=[]
    d0=[]
    single=[]
    for i in range(len(clusters_series)):
        if dimensions[i][3] == 3:
            d3.append(i)
        elif dimensions[i][3] == 2:
            d2.append(i)
        elif dimensions[i][3] == 1:
            d1.append(i)
        elif len(clusters_series[i]) == 1:
            single.append(i)
        else:
            d0.append(i)
    delete=d3+single
    delete.sort(reverse=True)
    for i in delete:
        structures.pop(i)
    structures=remove_duplicate_by_Structure(structures) # remove duplicate clusters
    dimensions=dimensions_by_structures(structures,2,'covalent') # judge the dimension of independent clusters
    structures_old=copy.deepcopy(structures)
    for num,i in enumerate(structures):
        structures[num]=standerize_cluster(dimensions[num],structures[num])
    lattices=lattice_clusters(structures,s) # get lattice information of the clusters
    parent=parent_info(filename,s,'mp')

    clusters_info=[]
    for i in range(len(structures)):
        info={'structure':structures[i].as_dict(),'structure_old':structures_old[i].as_dict(),'parent':parent,'valence':{},'dimension':dimensions[i]}
        info=dict(info,**lattices[i])
        clusters_info.append(info)
    return clusters_info




