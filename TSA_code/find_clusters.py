# -*- coding: utf-8 -*-
"""
Created on Sat Feb 12 21:23:49 2022

@author: DengJun
"""


from pymatgen.core.structure import Structure
import numpy as np
from radius import radius

def find_clusters(s,bond_type):
    '''
    find the clusters in the structure
    
    Parameters
    ----------
    s : TYPE: pymatgen format structure file
        DESCRIPTION: define the structure
    bond_type : TYPE: string
        DESCRIPTION: define the type of bond type, Options include 'vdw' and 'covalent'

    Returns
    -------
    TYPE: list
        DESCRIPTION: return the clusters in the structure.

    '''
    # s=s.makes
    atomic_pos=s.cart_coords
    atoms=s.atomic_numbers
    distance_matrix=s.distance_matrix
    # distance_matrix=np.zeros([len(atomic_pos[:,0]),len(atomic_pos[:,0])])
    # for i in range(len(atomic_pos[:,0])):
    #     for j in range(len(atomic_pos[:,0])):
    #         distance_matrix[i,j]=np.linalg.norm(atomic_pos[i,:]-atomic_pos[j,:])    
    
    [r_covalent,r_vdw]=radius()
    
    bonds=np.zeros([len(atomic_pos[:,0]),len(atomic_pos[:,0])])
    for i in range(len(atoms)):
        for j in range(len(atoms)):
            if i == j:
                continue
            else:
                if bond_type == 'covalent':
                    if distance_matrix[i,j] <= (r_covalent[atoms[i]]['radius']+r_covalent[atoms[j]]['radius'])*1.2:
                        bonds[i,j]=1
                if bond_type == 'vdw':
                    if distance_matrix[i,j] <= r_vdw[atoms[i]]['radius']+r_vdw[atoms[j]]['radius']-1.3:
                        bonds[i,j]=1
    
    def list_combine(l):
        if  len(l) == 1:
            return list(set(l[0]))
        else:
            t=[]
            for i in l:
                t=t+i
            return list(set(t))
    
    def find_cluster(bonds,cluster):
        tmp=[]
        for i in cluster[-1]:
            link=[]
            link.append(i)
            for j in range(len(bonds)):
                if bonds[i][j] == 1:
                    link.append(j)
            tmp=tmp+link
        new=list(set(tmp).difference(set(list_combine(cluster))))
        if new == []:
            return cluster
        cluster.append(new)
        return find_cluster(bonds,cluster)
                 
    clusters=[]
    for i in range(len(atoms)):
        cluster=[]
        if (i in list_combine(clusters)) == False:
            cluster.append([i])
            cluster=find_cluster(bonds,cluster)
            clusters.append(list_combine(cluster))    
            
    return clusters
