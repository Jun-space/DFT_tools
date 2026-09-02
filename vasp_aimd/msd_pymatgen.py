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

plt.figure()
plt.plot(pymatgen_diff1.dt / 1000, pymatgen_diff1.msd, label='H')
plt.plot(pymatgen_diff2.dt / 1000, pymatgen_diff2.msd, label='C')
plt.plot(pymatgen_diff3.dt / 1000, pymatgen_diff3.msd, label='O')
plt.ylabel('MSD/Å$^2$')
plt.xlabel(r'$\Delta t$/ps')
plt.show()
plt.legend()
plt.savefig('msd_tot.png',dpi=600,bbox_inches='tight',transparent=True)

plt.figure()
plt.plot(pymatgen_diff1.dt / 1000, pymatgen_diff1.msd_components[:,0],label='Hx')
plt.plot(pymatgen_diff2.dt / 1000, pymatgen_diff1.msd_components[:,1],label='Hy')
plt.plot(pymatgen_diff3.dt / 1000, pymatgen_diff1.msd_components[:,2],label='Hz')
plt.legend()
plt.ylabel('MSD/Å$^2$')
plt.xlabel(r'$\Delta t$/ps')
plt.legend()
plt.show()
plt.savefig('msd_H_comp.png',dpi=600,bbox_inches='tight',transparent=True)

plt.figure()
plt.plot(pymatgen_diff1.dt / 1000, np.sqrt(pymatgen_diff1.msd), label='H')
plt.plot(pymatgen_diff2.dt / 1000, np.sqrt(pymatgen_diff2.msd), label='C')
plt.plot(pymatgen_diff3.dt / 1000, np.sqrt(pymatgen_diff3.msd), label='O')
plt.ylabel('RMSD/Å')
plt.xlabel(r'$\Delta t$/ps')
plt.legend()
plt.show()
plt.savefig('rmsd_tot.png',dpi=600,bbox_inches='tight',transparent=True)


t=np.reshape(pymatgen_diff1.dt/1000,[len(pymatgen_diff1.dt),1])
H=np.reshape(pymatgen_diff1.msd,[len(pymatgen_diff1.dt),1])
C=np.reshape(pymatgen_diff2.msd,[len(pymatgen_diff2.dt),1])
O=np.reshape(pymatgen_diff3.msd,[len(pymatgen_diff3.dt),1])

H_comp=pymatgen_diff1.msd_components
C_comp=pymatgen_diff2.msd_components
O_comp=pymatgen_diff3.msd_components

H_data=np.concatenate((t,H,H_comp),axis=1)
C_data=np.concatenate((t,C,C_comp),axis=1)
O_data=np.concatenate((t,O,O_comp),axis=1)
header='             t(ps)         tot-MSD (A2)       tot-MSD-x (A2)       tot-MSD-y (A2)       tot-MSD-z (A2)'
np.savetxt('MSD_H_pymatgen.dat',H_data,fmt='%20.8f',header=header)
np.savetxt('MSD_C_pymatgen.dat',C_data,fmt='%20.8f',header=header)
np.savetxt('MSD_O_pymatgen.dat',O_data,fmt='%20.8f',header=header)

chg_conductivity_H=pymatgen_diff1.chg_conductivity
chg_diffusivity_H=pymatgen_diff1.chg_diffusivity
conductivity_H=pymatgen_diff1.conductivity
diffusivity_H=pymatgen_diff1.diffusivity

chg_conductivity_C=pymatgen_diff2.chg_conductivity
chg_diffusivity_C=pymatgen_diff2.chg_diffusivity
conductivity_C=pymatgen_diff2.conductivity
diffusivity_C=pymatgen_diff2.diffusivity

chg_conductivity_O=pymatgen_diff3.chg_conductivity
chg_diffusivity_O=pymatgen_diff3.chg_diffusivity
conductivity_O=pymatgen_diff3.conductivity
diffusivity_O=pymatgen_diff3.diffusivity

diffusion=np.array([[chg_conductivity_H,chg_diffusivity_H,conductivity_H,diffusivity_H],\
                    [chg_conductivity_C,chg_diffusivity_C,conductivity_C,diffusivity_C],\
                    [chg_conductivity_O,chg_diffusivity_O,conductivity_H,diffusivity_O]])

header='     chg_conductivity(mS/cm)        chg_diffusivity(cm^2/s)            conductivity(mS/cm)            diffusivity(cm^2/s)   order(H C O)'
np.savetxt('diffusion.dat',diffusion,fmt='%30.18e',header=header)


diffusion_H_component=pymatgen_diff1.diffusivity_components
diffusion_C_component=pymatgen_diff2.diffusivity_components
diffusion_O_component=pymatgen_diff3.diffusivity_components

diffusion_component=np.array([diffusion_H_component,diffusion_C_component,diffusion_O_component])

np.savetxt('diffusion_component.dat',diffusion_component,fmt='%30.18e',header=header)
