#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 23 19:23:08 2022

@author: dj
"""
from pymatgen.ext.matproj import MPRester
from pymatgen.io.vasp import Vasprun
from pymatgen.analysis.phase_diagram import PhaseDiagram, PDPlotter
from pymatgen.entries.compatibility import MaterialsProjectCompatibility
from pandas import DataFrame
import numpy as np
import copy
import os
import sys

def get_phase_diagram():
    entries=[]
    datas=[]
    filenames=os.listdir('./structures')
    for filename in filenames:
        prefix=filename.split('.')[0]
        try:
            vasprun = Vasprun('%s/scf/vasprun.xml'%prefix)
            entry = vasprun.get_computed_entry(inc_structure=True)
            entries.append(entry)
            datas.append(prefix)
            print(prefix)
        except:
            print(prefix,'error')

    pd = PhaseDiagram(entries)
    plotter = PDPlotter(pd, show_unstable=True)
    plotter.show()
    plotter.write_image('phase.svg')

    labels=['Fe','Si','B']
    print("{:<20}{:<20}{:<20}{:<20}{:<20}{:<20}{:<20}{:<20}{:<20}{:<20}".format('a', 'b', 'c', 'ehull','eformation','e_per_atom','spg', 'composition', 'formula', 'folder', 'pressure'))
    for num,entry in enumerate(entries):
        elements=entry.as_dict()['composition'].keys()
        if labels[0] in elements:
            a = entry.as_dict()['composition'][labels[0]]
        else:
            a = 0 
        if labels[1] in elements:
            b = entry.as_dict()['composition'][labels[1]]
        else:
            b = 0 
        if labels[2] in elements:
            c = entry.as_dict()['composition'][labels[2]]
        else:
            c = 0

        print("{:<20.3f}{:<20.3f}{:<20.3f}{:<20.8f}{:<20.8f}{:<20.8f}{:<20}{:<20}{:<20}".format(a/(a+b+c), b/(a+b+c), c/(a+b+c), pd.get_e_above_hull(entry=entry), pd.get_form_energy_per_atom(entry=entry),entry.energy_per_atom,entry.structure.get_space_group_info()[0] ,str(entry.composition), datas[num]))
        #if pd.get_e_above_hull(entry=entry) == 0:
        #    filename='%s-%s-%s.vasp'%(datas[num][0],entry.structure.get_space_group_info()[1],datas[num][1])
        #    entry.structure.to(filename=filename,fmt='poscar')

    return

if __name__ == "__main__":
    get_phase_diagram()
