# -*- coding: utf-8 -*-
"""
Created on Fri Jan 27 15:49:48 2023

@author: DengJun
"""


import pymongo
import os,sys
import re
import shutil
from pymatgen.ext.matproj import MPRester
from pymatgen.io.vasp.inputs import VaspInput
from pymatgen.io.vasp.inputs import Incar, Kpoints, Poscar,Potcar
from pymatgen.io.vasp.outputs import Outcar, Vasprun
from pymatgen.io.vasp.sets import VaspInputSet,MPRelaxSet,DictSet,MPStaticSet
from pymatgen.analysis.magnetism import CollinearMagneticStructureAnalyzer
from pymatgen.core.structure import Structure
import itertools
import numpy as np
import datetime
import time
import copy


def create_mag_structures(s_dict,mag_ele):
    s = Structure.from_dict(s_dict) #get structure
    site_symbol=[i.name for i in s.species] #get site symbols
    mag_ele=mag_ele #get magnetic moment for each site
    
    ####get spin state for each element according to mag_ele
    spin={}
    for i in range(len(site_symbol)):
        spin[site_symbol[i]]=mag_ele[i]
    mag=[]
    site_symbol=[i.name for i in s.species]
    for i,site in enumerate(site_symbol):
        if spin[site]>=0.5:
            mag.append(i)
    mag_tmp=[str(i) for i in mag]

    ####### combination of afm structures
    """
    Note: if the number of magnetic atoms are odd, 
    we build a 2*2*2 supercell
    """
    if len(mag_tmp)%2==1:
        s.make_supercell([2,2,2])
        mag=[]
        site_symbol=[i.name for i in s.species]
        for i,site in enumerate(site_symbol):
            if spin[site]>=0.5:
                mag.append(i)
        mag_tmp=[str(i) for i in mag]
    # else:
    #     s.make_supercell([1,1,2])
    #     mag=[]
    #     site_symbol=[i.name for i in s.species]
    #     for i,site in enumerate(site_symbol):
    #         if spin[site]>=0.5:
    #             mag.append(i)
    #     mag_tmp=[str(i) for i in mag]
        
    combinations=[]
    for i in itertools.combinations(mag_tmp,int(len(mag_tmp)/2)):
        combinations.append(i)
    
    remove=[]
    for num1,i in enumerate(combinations):
        for num2,j in enumerate(combinations):
            tmp=list(i+j)
            if list(set(mag_tmp).difference(set(tmp))) == []:
                combinations.pop(num2)

###############get FM and AFM structures
    fm=[]
    for i in site_symbol:
        fm.append(spin[i])

###############AFM
    afm=[]
    for i in combinations:
        spins=[]
        tmp=np.ones(len(site_symbol))
        for k in list(i):
            tmp[int(k)]=tmp[int(k)]*-1
        for numb,j in enumerate(site_symbol):
            spins.append(spin[j])
        spins=list(np.array(spins)*tmp)
        afm.append(spins)
    
    mag_config=[fm]+afm
    
    mag_structures=[]
    for i in mag_config:
        s_tmp=copy.deepcopy(s)
        s_tmp.add_spin_by_site(i)
        mag_structures.append(s_tmp)
        
########remove duplicated magnetic structures
    def compare_two(s1,s2):
        mat=CollinearMagneticStructureAnalyzer(s1)
        return mat.matches_ordering(s2) 
    
    d=np.zeros([len(mag_structures),len(mag_structures)])
    for i in range(len(mag_structures)):
        for j in range(i,len(mag_structures)):
            if i == j:
                continue
            else:
                if compare_two(mag_structures[i],mag_structures[j]):
                    d[i,j]=1
                    d[j,i]=1
    
    def list_combine(l):
        if  len(l) == 1:
            return list(set(l[0]))
        else:
            t=[]
            for i in l:
                t=t+i
            return list(set(t))
    
    def find_match(d,match):
        tmp=[]
        for i in match[-1]:
            link=[]
            link.append(i)
            for j in range(len(d)):
                if d[i][j] == 1:
                    link.append(j)
            tmp=tmp+link
        new=list(set(tmp).difference(set(list_combine(match))))
        if new == []:
            return match
        match.append(new)
        return find_match(d,match)
        
    matches=[]
    for i in range(len(mag_structures)):
        match=[]
        if (i in list_combine(matches)) == False:
            match.append([i])
            match=find_match(d,match)
            matches.append(list_combine(match))
    keep=[]
    for i in matches:
        keep.append(i[0])
    remove=list(set(np.arange(len(mag_structures))).difference(set(keep)))
    remove.sort(reverse=True)
    for i in remove:
        mag_structures.pop(i)
    for num,i in enumerate(mag_structures):
        mag_structures[num]=i.as_dict()

    return mag_structures

s=Structure.from_file('CONTCAR')
s_dict=s.as_dict()
mag_ele=[1,0,0]
mag_struct=create_mag_structures(s_dict,mag_ele)
for num,i in enumerate(mag_struct):
    s=Structure.from_dict(i)
    s.to(fmt='poscar',filename='s'+str(num)+'.vasp')