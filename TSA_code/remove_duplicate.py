# -*- coding: utf-8 -*-
"""
Created on Wed Feb 23 19:11:11 2022

@author: DengJun

remove duplicated structures for a set of structures
"""

from pymatgen.analysis.structure_matcher import StructureMatcher
from pymatgen.core.structure import Structure
import numpy as np


def remove_duplicate_by_Structure(cluster_structures):
    def compare_two(s1,s2):
        mat=StructureMatcher()
        return mat.fit(s1,s2)
    
    d=np.zeros([len(cluster_structures),len(cluster_structures)])
    for num1,i in enumerate(cluster_structures):
        for num2,j in enumerate(cluster_structures):
            if num1 == num2:
                continue
            else:
                if compare_two(i,j):
                    d[num1,num2]=1
    
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
            return
        match.append(new)
        find_match(d,match)
        
    matches=[]
    for i in range(len(cluster_structures)):
        match=[]
        if (i in list_combine(matches)) == False:
            match.append([i])
            find_match(d,match)
            matches.append(list_combine(match))
    keep=[]
    for i in matches:
        keep.append(i[0])
    remove=list(set(np.arange(len(cluster_structures))).difference(set(keep)))
    remove.sort(reverse=True)
    for i in remove:
        cluster_structures.pop(i)
    
    return cluster_structures


def remove_duplicate_by_dict(cluster_dicts):
    def compare_two(s1,s2):
        mat=StructureMatcher()
        return mat.fit(s1,s2)
    
    d=np.zeros([len(cluster_dicts),len(cluster_dicts)])
    cluster_structures=[Structure.from_dict(i['structure']) for i in cluster_dicts]
    for i in range(len(cluster_structures)):
        for j in range(i,len(cluster_structures)):
            if i == j:
                continue
            else:
                if compare_two(cluster_structures[i],cluster_structures[j]):
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
    for i in range(len(cluster_structures)):
        match=[]
        if (i in list_combine(matches)) == False:
            match.append([i])
            match=find_match(d,match)
            matches.append(list_combine(match))
    keep=[]
    valences=[]
    parents=[]
    for i in matches:
        keep.append(i[0])
        v_tmp=[]
        p_tmp=[]
        for j in i:
            v_tmp.append(cluster_dicts[j]['valence'])
            p_tmp.append(cluster_dicts[j]['parent'])
        valences.append(v_tmp)
        parents.append(p_tmp)
    for num,i in enumerate(keep):
        cluster_dicts[i]['parent']=parents[num]
        cluster_dicts[i]['valence']=valences[num]
    remove=list(set(np.arange(len(cluster_structures))).difference(set(keep)))
    remove.sort(reverse=True)
    for i in remove:
        cluster_dicts.pop(i)
    
    return cluster_dicts
