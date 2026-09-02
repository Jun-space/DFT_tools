# -*- coding: utf-8 -*-
"""
Created on Sun Feb 20 16:01:11 2022

@author: DengJun

return structural information for the network
"""

from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

def lattice_cluster(cluster_structure,parent):
    a=cluster_structure.lattice.a
    b=cluster_structure.lattice.b
    c=cluster_structure.lattice.c
    alpha=cluster_structure.lattice.alpha
    beta=cluster_structure.lattice.beta
    gamma=cluster_structure.lattice.gamma
    formula=cluster_structure.formula
    spg=SpacegroupAnalyzer(cluster_structure,symprec=0.1)
    crys_sys=spg.get_crystal_system()
    spg=cluster_structure.get_space_group_info()
    ntype_parent=list(set(parent.atomic_numbers))
    ntype_cluster=list(set(cluster_structure.atomic_numbers))
    if ntype_parent != ntype_cluster:
        type='ionic'
    else:
        type='vdw'
    return {'a' : a,'b' : b, 'c':c,'alpha':alpha,'beta':beta,'gamma':gamma,\
            'formula':formula,'spg':spg,'crys_sys': crys_sys,'type':type}

def lattice_clusters(cluster_structures,parent):
    info=[]
    for num,i in enumerate(cluster_structures):
        info.append(lattice_cluster(i,parent))
    return info
