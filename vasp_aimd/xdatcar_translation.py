# -*- coding: utf-8 -*-
"""
Created on Mon Apr 28 09:47:08 2025

@author: dengjun
"""

import numpy as np
from pymatgen.analysis.diffusion.analyzer import DiffusionAnalyzer as PymatgenDiffusionAnalyzer
from pymatgen.io.vasp import Xdatcar
from pymatgen.core.structure import Structure

xd=Xdatcar('XDATCAR')
structures=xd.structures
structures_trans=[]
for num,structure in enumerate(structures):
    matrix=structure.lattice.matrix
    symbols=structure.atomic_numbers
    cart_coords=structure.cart_coords
    center_of_mass = np.mean(cart_coords, axis=0)
    cart_coords_trans=cart_coords-center_of_mass+np.dot(np.array([0.5,0.5,0.5]),matrix)
    structure_trans=Structure(matrix,symbols,cart_coords_trans,coords_are_cartesian=True)
    structures_trans.append(structure_trans)
    print(num)

xd.structures=structures_trans
xd.write_file('XDATCAR-trans-center')