# -*- coding: utf-8 -*-
"""
Created on Tue Jul 22 17:10:18 2025

@author: dengjun
"""

from mechelastic.core import ELATE
import pyvista as pv
import numpy as np

elastic_tensor=np.loadtxt('elastic_tensor_standard.dat').tolist()
density=4.73675132485115*10**3
elate = ELATE(s = elastic_tensor, density = density)

# meshes = elate.plot_3D(elastic_calc = "YOUNG", show = True)
# fig = elate.plot_2D(elastic_calc = "YOUNG", show = True )
elate.print_properties()