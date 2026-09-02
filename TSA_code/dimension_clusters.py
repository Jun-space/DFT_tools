# -*- coding: utf-8 -*-
"""
Created on Sat Feb 19 20:44:26 2022

@author: DengJun
"""



from radius import radius
from find_clusters import find_clusters
from pymatgen.core.structure import Structure
import numpy as np

def cluster_structures(s,clusters):
    cluster_structures=[]        
    for num,i in enumerate(clusters):
        s_tmp=s.copy()
        atoms=s_tmp.atomic_numbers
        re=sorted(list(set(np.arange(len(atoms))).difference(set(i))))[::-1]
        for j in re:
            s_tmp.remove(s_tmp.sites[j])
        cluster_structures.append(s_tmp)
    return cluster_structures

def dimensions_by_series(s,clusters,n,bond_type):
    
    def cluster_structures(s,clusters):
        cluster_structures=[]        
        for num,i in enumerate(clusters):
            s_tmp=s.copy()
            atoms=s_tmp.atomic_numbers
            re=sorted(list(set(np.arange(len(atoms))).difference(set(i))))[::-1]
            for j in re:
                s_tmp.remove(s_tmp.sites[j])
            cluster_structures.append(s_tmp)
        return cluster_structures
    
    def cluster_dimension(cluster_structure,n,bond_type):   
        ####judge dimension
        nnn=cluster_structure.copy()
        nnn.make_supercell([n,n,n])
        ncluster_nnn=len(find_clusters(nnn,bond_type))
        if int((n*n*n)/(ncluster_nnn)) == 1:
            d=0
        if int((n*n*n)/(ncluster_nnn)) == n:
            d=1
        if int((n*n*n)/(ncluster_nnn)) == n*n:
            d=2
        if int((n*n*n)/(ncluster_nnn)) == n*n*n:
            d=3
            
        ####periodic along x
            '''
            dx=0 means the cluster is not periodic along x direction.
            dx=1 means the cluster is periodic along x direction.
            '''
        nx=cluster_structure.copy()
        nx.make_supercell([n,1,1])
        ncluster_nx=len(find_clusters(nx,bond_type))
        if int((n)/(ncluster_nx)) == 1:
            dx=0
        if int((n)/(ncluster_nx)) == n:
            dx=1
        # print('dx: '+str(dx))   
            
        ####periodic along y
            # '''
            # dy=0 means the cluster is not periodic along y direction.
            # dy=1 means the cluster is periodic along y direction.
            # '''
        ny=cluster_structure.copy()
        ny.make_supercell([1,n,1])
        ncluster_ny=len(find_clusters(ny,bond_type))
        # print('ny: '+str(ncluster_ny))
        if int((n)/(ncluster_ny)) == 1:
            dy=0
        if int((n)/(ncluster_ny)) == n:
            dy=1
            
        ####periodic along y
            '''
            dz=0 means the cluster is not periodic along z direction.
            dz=1 means the cluster is periodic along z direction.
            '''
        nz=cluster_structure.copy()
        nz.make_supercell([1,1,n])
        ncluster_nz=len(find_clusters(nz,bond_type))
        # print('nz: '+str(ncluster_nz))
        if int((n)/(ncluster_nz)) == 1:
            dz=0
        if int((n)/(ncluster_nz)) == n:
            dz=1
        dimension=[dx,dy,dz,d]
        
        return dimension    
    
    cluster_structures=cluster_structures(s,clusters)
    d=[]
    for cluster_structure in cluster_structures:
        d.append(cluster_dimension(cluster_structure,n,bond_type))
    return d


def dimensions_by_structures(cluster_structures,n,bond_type):  
    
    def cluster_dimension(cluster_structure,n,bond_type):   
        ####judge dimension
        nnn=cluster_structure.copy()
        nnn.make_supercell([n,n,n])
        ncluster_nnn=len(find_clusters(nnn,bond_type))
        if int((n*n*n)/(ncluster_nnn)) == 1:
            d=0
        if int((n*n*n)/(ncluster_nnn)) == n:
            d=1
        if int((n*n*n)/(ncluster_nnn)) == n*n:
            d=2
        if int((n*n*n)/(ncluster_nnn)) == n*n*n:
            d=3
            
        ####periodic along x
            '''
            dx=0 means the cluster is not periodic along x direction.
            dx=1 means the cluster is periodic along x direction.
            '''
        nx=cluster_structure.copy()
        nx.make_supercell([n,1,1])
        ncluster_nx=len(find_clusters(nx,bond_type))
        if int((n)/(ncluster_nx)) == 1:
            dx=0
        if int((n)/(ncluster_nx)) == n:
            dx=1
            
        ####periodic along y
            '''
            dy=0 means the cluster is not periodic along y direction.
            dy=1 means the cluster is periodic along y direction.
            '''
        ny=cluster_structure.copy()
        ny.make_supercell([1,n,1])
        ncluster_ny=len(find_clusters(ny,bond_type))
        if int((n)/(ncluster_ny)) == 1:
            dy=0
        if int((n)/(ncluster_ny)) == n:
            dy=1
            
        ####periodic along y
            '''
            dz=0 means the cluster is not periodic along z direction.
            dz=1 means the cluster is periodic along z direction.
            '''
        nz=cluster_structure.copy()
        nz.make_supercell([1,1,n])
        ncluster_nz=len(find_clusters(nz,bond_type))
        if int((n)/(ncluster_nz)) == 1:
            dz=0
        if int((n)/(ncluster_nz)) == n:
            dz=1
        
        dimension=[dx,dy,dz,d]
        
        return dimension    
    d=[]
    for cluster_structure in cluster_structures:
        d.append(cluster_dimension(cluster_structure,n,bond_type))
    return d

# def cluster_dimension(cluster_structure,n,bond_type):   
#     ####judge dimension
#     nnn=cluster_structure.copy()
#     nnn.make_supercell([n,n,n])
#     ncluster_nnn=len(find_clusters(nnn,bond_type))
#     if int((n*n*n)/(ncluster_nnn)) == 1:
#         d=0
#     if int((n*n*n)/(ncluster_nnn)) == n:
#         d=1
#     if int((n*n*n)/(ncluster_nnn)) == n*n:
#         d=2
#     if int((n*n*n)/(ncluster_nnn)) == n*n*n:
#         d=3
        
#     ####periodic along x
#         '''
#         dx=0 means the cluster is not periodic along x direction.
#         dx=1 means the cluster is periodic along x direction.
#         '''
#     nx=cluster_structure.copy()
#     nx.make_supercell([n,1,1])
#     ncluster_nx=len(find_clusters(nx,bond_type))
#     if int((n)/(ncluster_nx)) == 1:
#         dx=0
#     if int((n)/(ncluster_nx)) == n:
#         dx=1
        
#     ####periodic along y
#         '''
#         dy=0 means the cluster is not periodic along y direction.
#         dy=1 means the cluster is periodic along y direction.
#         '''
#     ny=cluster_structure.copy()
#     ny.make_supercell([1,n,1])
#     ncluster_ny=len(find_clusters(ny,bond_type))
#     if int((n)/(ncluster_ny)) == 1:
#         dy=0
#     if int((n)/(ncluster_ny)) == n:
#         dy=1
        
#     ####periodic along y
#         '''
#         dz=0 means the cluster is not periodic along z direction.
#         dz=1 means the cluster is periodic along z direction.
#         '''
#     nz=cluster_structure.copy()
#     nz.make_supercell([1,1,n])
#     ncluster_nz=len(find_clusters(nz,bond_type))
#     if int((n)/(ncluster_nz)) == 1:
#         dz=0
#     if int((n)/(ncluster_nz)) == n:
#         dz=1
#     dimension=[dx,dy,dz,d]
    
#     return dimension  


    
    
    
    
