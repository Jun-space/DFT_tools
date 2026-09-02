# -*- coding: utf-8 -*-
"""
Created on Sun Feb 20 09:49:43 2022

@author: DengJun
"""


from pymatgen.core.structure import Structure


##AM: Alkali metals, AEM: alkaline earth metals
def delete_am_aem(s):
    am=[3,11,19,37,55]  #AM: Alkali metals
    aem=[12,20,38,56]
    atoms=s.atomic_numbers
    counts=[]
    for num,i in enumerate(atoms):
        if (i in am):
            counts.append(num)
        if (i in aem):
            counts.append(num)
    counts.sort(reverse=True)
    for j in counts:
        s.remove(s.sites[j])
    return s