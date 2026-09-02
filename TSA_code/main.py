# -*- coding: utf-8 -*-
"""
Created on Sun Feb 20 09:28:15 2022

@author: Deng Jun
email: 355458045@qq.com

This is main part of TSA code
It identifies the bonded networks for a given structure
Two atoms are considered bonded when dij <= k * (ri+rj), where ri and rj are covalent radius in Dalton Trans, 2832 (2008)
Scaling parameter "k" can be modified at line 48 of find_clusters, default value is 1.05

The return values are bonded networks found for the given structure, it contains the following information:
1. "structure": the structure of the network with dict format, which can be read by pymatgen.core.structure.Structure.
                the aperiodic direction is normalized to 25 angstrom
2. "structure_old": this is same with "structure but aperiodic direction is normalized.
                    This structure retains because sometimes the "normalization" code does not work for some structure.
3. "parent": source of the network, extracted in the filenames
4. "valence": possible valence state of the network
5. "dimension": list, four values, the first three show the periodicity of the network in x,y,z direction, respectively. "1" means it is periodic and "0" is not.
                the last one show the dimension of the network, which can be 0-3.
6. "a“，“b”,"c","alpha","beta","gamma": lattice constant of the network.
7. "spg": space group of the network.
8. "crys_sys": crystal system of the network.
9. "type": type of the network, which can be "vdw" or "ionic"



"""

from clusters import clusters
from pymatgen.core.structure import Structure

networks=clusters(filename='mp-001_CO3H2.vasp') #file can be "POSCAR" or "cif" format, but should be named with "mp-id+formula"
for num,network in enumerate(networks):
    s=Structure.from_dict(network['structure_old'])
    s.to(filename='%d.vasp'%num,fmt='poscar')
