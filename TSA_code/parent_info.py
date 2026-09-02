# -*- coding: utf-8 -*-
"""
Created on Thu Feb 24 10:47:12 2022

@author: DengJun
"""

import re
from pymatgen.core.structure import Structure

def parent_info(filename,s,database):
    formula=s.formula
    if database == 'ICSD':
        id = re.sub("\D", "", filename)
    elif database == 'mp':
        id = filename.split('_')[0]
    else:
        id =filename.split('/')[1].split('.')[0]
    return {'database':database,'id':id,'formula':formula}
    