# -*- coding: utf-8 -*-
"""
Created on Sun Feb 20 09:28:15 2022

@author: DengJun

This code standardize the structure and normalize the aperiodic direction to 25 angstrom
"""


import numpy as np
from pymatgen.core.structure import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.core import Lattice, Structure
import copy


def output_clusters(s,clusters):
    cluster_structures=[]        
    for num,i in enumerate(clusters):
        s_tmp=s.copy()
        atoms=s_tmp.atomic_numbers
        re=sorted(list(set(np.arange(len(atoms))).difference(set(i))))[::-1]
        for j in re:
            s_tmp.remove(s_tmp.sites[j])
        s_tmp=SpacegroupAnalyzer(s_tmp).get_conventional_standard_structure()
        cluster_structures.append(s_tmp)
    return cluster_structures

def generate_poscar(cluster_structures,path):
    for num,i in enumerate(cluster_structures):
        i.to(fmt='poscar',filename=path+'/'+str(num)+'.vasp')
        
def generate_poscar_by_dict(cluster_structures,path):
    for num,i in enumerate(cluster_structures):
        s=Structure.from_dict(i['structure'])
        s.to(fmt='poscar',filename=path+'/'+str(num)+'_'+i['formula'].replace(' ','')+'.vasp')

def standerize_clusters(dimensions,cluster_structures):
    redefine_structures=[]
    for i,cluster in enumerate(cluster_structures):
        print(i)
        symbols=cluster.atomic_numbers
        a=cluster.lattice.a
        b=cluster.lattice.b
        c=cluster.lattice.c
        alpha=cluster.lattice.alpha
        beta=cluster.lattice.beta
        gamma=cluster.lattice.gamma
        matrix=cluster.lattice.matrix
        coords_frac=cluster.frac_coords
        cluster_dimension=dimensions[i]
        
        
        def is_bottom(direction,coords_frac):
            '''
            if up-bottom is larger than 0.5, the atoms are at the both sides (fractional coordination)
            '''

            bottom=np.min(coords_frac[:,direction])
            top=np.max(coords_frac[:,direction])
            if top-bottom >=0.5:
                return True
            else:
                return False
        for direction in range(3):
            if cluster_dimension[direction] == 0:
                if is_bottom(direction,coords_frac):
                    for i in range(len(coords_frac)):
                        if coords_frac[i,direction] >0.5:
                            coords_frac[i,direction]=coords_frac[i,direction]-0.5
                        else:
                            coords_frac[i,direction]=coords_frac[i,direction]+0.5
                bottom=np.min(coords_frac[:,direction])
                top=np.max(coords_frac[:,direction])
                t_avg=(bottom+top)/2
                coords_frac[:,direction]=coords_frac[:,direction]+0.5-t_avg
                coords_cart=np.dot(coords_frac,matrix)
                matrix_new=copy.deepcopy(matrix)
                matrix_new[direction,:]=25/np.linalg.norm(matrix[direction,:])*matrix[direction,:]
                delta=matrix_new[direction,:]-matrix[direction,:]
                for i in np.arange(0,len(coords_cart)):
                    coords_cart[i,:]=coords_cart[i,:]+delta
                coords_frac=np.dot(coords_cart,np.linalg.inv(matrix_new))
                bottom=np.min(coords_frac[:,direction])
                top=np.max(coords_frac[:,direction])
                t_avg=(bottom+top)/2
                coords_frac[:,direction]=coords_frac[:,direction]+0.5-t_avg
        if cluster_dimension[0] == 0:
            a=25
        if cluster_dimension[1] == 0:
            b=25
        if cluster_dimension[2] == 0:
            c=25
        lattice=Lattice.from_parameters(a=a,b=b,c=c,alpha=alpha,beta=beta,gamma=gamma)
        struct=Structure(lattice,symbols,coords_frac,coords_are_cartesian=False)
        redefine_structures.append(struct)
    return redefine_structures


def standerize_cluster(dimension,cluster_structure):
    cluster=cluster_structure
    symbols=cluster.atomic_numbers
    a=cluster.lattice.a
    b=cluster.lattice.b
    c=cluster.lattice.c
    alpha=cluster.lattice.alpha
    beta=cluster.lattice.beta
    gamma=cluster.lattice.gamma
    matrix=cluster.lattice.matrix
    coords_frac=cluster.frac_coords
    cluster_dimension=dimension
    
    
    def is_bottom(direction,coords_frac):
        '''
        if up-bottom is larger than 0.5, the atoms are at the both sides (fractional coordination)
        '''

        bottom=np.min(coords_frac[:,direction])
        top=np.max(coords_frac[:,direction])
        if top-bottom >=0.5:
            return True
        else:
            return False
    for direction in range(3):
        if cluster_dimension[direction] == 0:
            if is_bottom(direction,coords_frac):
                for i in range(len(coords_frac)):
                    if coords_frac[i,direction] >0.5:
                        coords_frac[i,direction]=coords_frac[i,direction]-0.5
                    else:
                        coords_frac[i,direction]=coords_frac[i,direction]+0.5
            bottom=np.min(coords_frac[:,direction])
            top=np.max(coords_frac[:,direction])
            t_avg=(bottom+top)/2
            coords_frac[:,direction]=coords_frac[:,direction]+0.5-t_avg
            coords_cart=np.dot(coords_frac,matrix)
            matrix_new=copy.deepcopy(matrix)
            matrix_new[direction,:]=25/np.linalg.norm(matrix[direction,:])*matrix[direction,:]
            delta=matrix_new[direction,:]-matrix[direction,:]
            for i in np.arange(0,len(coords_cart)):
                coords_cart[i,:]=coords_cart[i,:]+delta
            coords_frac=np.dot(coords_cart,np.linalg.inv(matrix_new))
            bottom=np.min(coords_frac[:,direction])
            top=np.max(coords_frac[:,direction])
            t_avg=(bottom+top)/2
            coords_frac[:,direction]=coords_frac[:,direction]+0.5-t_avg
    if cluster_dimension[0] == 0:
        a=25
    if cluster_dimension[1] == 0:
        b=25
    if cluster_dimension[2] == 0:
        c=25
    lattice=Lattice.from_parameters(a=a,b=b,c=c,alpha=alpha,beta=beta,gamma=gamma)
    struct=Structure(lattice,symbols,coords_frac,coords_are_cartesian=False)
    return struct