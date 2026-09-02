# -*- coding: utf-8 -*-
"""
Created on Wed Apr  9 20:10:32 2025

@author: dengjun
"""

import numpy as np
from pymatgen.analysis.diffusion.analyzer import DiffusionAnalyzer as PymatgenDiffusionAnalyzer
from pymatgen.io.vasp import Xdatcar
import matplotlib.pyplot as plt
import sys
import matplotlib as mpl 
mpl.use('Agg')

temperature=int(sys.argv[1])
xd=Xdatcar('XDATCAR')
pymatgen_diff1 = PymatgenDiffusionAnalyzer.from_structures(xd.structures, temperature=temperature,specie='H',time_step=1,step_skip=1,smoothed=None)
pymatgen_diff2 = PymatgenDiffusionAnalyzer.from_structures(xd.structures, temperature=temperature,specie='C',time_step=1,step_skip=1,smoothed=None)
pymatgen_diff3 = PymatgenDiffusionAnalyzer.from_structures(xd.structures, temperature=temperature,specie='O',time_step=1,step_skip=1,smoothed=None)


diffusion_H_component=pymatgen_diff1.diffusivity_components
diffusion_C_component=pymatgen_diff2.diffusivity_components
diffusion_O_component=pymatgen_diff3.diffusivity_components

diffusion_component=np.array([diffusion_H_component,diffusion_C_component,diffusion_O_component])
header='                   Dx(cm^2/s)                     Dy(cm^2/s)                     Dz(cm^2/s)'
np.savetxt('diffusion_component.dat',diffusion_component,fmt='%30.18e',header=header)
