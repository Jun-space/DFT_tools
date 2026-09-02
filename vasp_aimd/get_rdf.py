#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  5 11:19:26 2024

@author: dj
"""

from pymatgen.core.structure import Structure
from pymatgen.io.vasp import Xdatcar
from vasppy.rdf import RadialDistributionFunction
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import rcParams
import matplotlib as mpl 
mpl.use('Agg')

rcParams['mathtext.fontset'] = 'custom'
rcParams['mathtext.rm'] = 'Arial'
rcParams['mathtext.it'] = 'Arial:italic'
rcParams['mathtext.bf'] = 'Arial:bold'
rcParams['pdf.fonttype'] = 42
rcParams['ps.fonttype'] = 42

font={'family':'Arial','weight':'normal','size':8}


print('read XDATCAR')
xd = Xdatcar('XDATCAR')
print('read XDATCAR susscessful!')
weights=np.zeros(9000).tolist()
weights[0:]=1*np.ones(len(weights[0:]))
indices=np.arange(0,len(xd.structures[0].atomic_numbers))
rdf_tot = RadialDistributionFunction(structures=xd.structures,indices_i=indices,weights=weights)
print('extract total rdf done')

sigma=0.05

rdf_CO=RadialDistributionFunction.from_species_strings(structures=xd.structures,species_i='C',species_j='O')
print('extract C-O rdf done')
rdf_CH=RadialDistributionFunction.from_species_strings(structures=xd.structures,species_i='C',species_j='H')
print('extract C-H rdf done')
rdf_OH=RadialDistributionFunction.from_species_strings(structures=xd.structures,species_i='O',species_j='H')
print('extract O-H rdf done')
rdf_HH=RadialDistributionFunction.from_species_strings(structures=xd.structures,species_i='H',species_j='H')
print('extract H-H rdf done')
rdf_CC=RadialDistributionFunction.from_species_strings(structures=xd.structures,species_i='C',species_j='C')
print('extract C-C rdf done')
rdf_OO=RadialDistributionFunction.from_species_strings(structures=xd.structures,species_i='O',species_j='O')
print('extract O-O rdf done')

x=rdf_tot.r
y=rdf_tot.smeared_rdf(sigma=sigma)
y_CO=rdf_CO.smeared_rdf(sigma=sigma)
y_CH=rdf_CH.smeared_rdf(sigma=sigma)
y_OH=rdf_OH.smeared_rdf(sigma=sigma)
y_HH=rdf_HH.smeared_rdf(sigma=sigma)
y_CC=rdf_CC.smeared_rdf(sigma=sigma)
y_OO=rdf_OO.smeared_rdf(sigma=sigma)
header='         r(A)             C-O             C-H             O-H             H-H             C-C             O-O             tot'
data=np.concatenate((x.reshape(len(x),1),y_CO.reshape(len(y_CO),1),\
                                         y_CH.reshape(len(y_CH),1),\
                                         y_OH.reshape(len(y_OH),1),\
                                         y_HH.reshape(len(y_HH),1),\
                                         y_CC.reshape(len(y_CC),1),\
                                         y_OO.reshape(len(y_OO),1),\
                                         y.reshape(len(y),1)),axis=1)
np.savetxt('rdf.dat',data,header=header,fmt='%15.8f')

fig=plt.figure()
ax=fig.add_subplot(331)
ax.plot(x, y, label='Tot',linewidth=0.5)
ax.set_xlabel(r'r ($\rm \AA$)',font)
ax.set_ylabel('g(r)',font)
ax.legend(prop=font,frameon=False)
ax.minorticks_on()
# ax.yaxis.set_minor_locator(plt.MultipleLocator(5))
ax.tick_params(axis='y',direction='in',which="major",width=0.5,length=3)
ax.tick_params(axis='y',direction='in',which="minor",width=0.5,length=1.5)

# ax.xaxis.set_minor_locator(plt.MultipleLocator(0.01))
ax.tick_params(axis='x',direction='in',which="major",width=0.5,length=3)
ax.tick_params(axis='x',direction='in',which="minor",width=0.5,length=1.5)
plt.tick_params(labelsize=6)
labels=ax.get_xticklabels()+ax.get_yticklabels()
[label.set_fontname('Arial') for label in labels]

#set line width
ax.spines['bottom'].set_linewidth(0.5)
ax.spines['left'].set_linewidth(0.5)
ax.spines['right'].set_linewidth(0.5)
ax.spines['top'].set_linewidth(0.5)


ax=fig.add_subplot(332)
ax.plot(x, y_CO, label='C-O',linewidth=0.5)
ax.set_xlabel(r'r ($\rm \AA$)',font)
ax.set_ylabel('g(r)',font)
ax.legend(prop=font,frameon=False)
ax.minorticks_on()
ax.yaxis.set_minor_locator(plt.MultipleLocator(5))
ax.tick_params(axis='y',direction='in',which="major",width=0.5,length=3)
ax.tick_params(axis='y',direction='in',which="minor",width=0.5,length=1.5)

# ax.xaxis.set_minor_locator(plt.MultipleLocator(0.01))
ax.tick_params(axis='x',direction='in',which="major",width=0.5,length=3)
ax.tick_params(axis='x',direction='in',which="minor",width=0.5,length=1.5)
plt.tick_params(labelsize=6)
labels=ax.get_xticklabels()+ax.get_yticklabels()
[label.set_fontname('Arial') for label in labels]

#set line width
ax.spines['bottom'].set_linewidth(0.5)
ax.spines['left'].set_linewidth(0.5)
ax.spines['right'].set_linewidth(0.5)
ax.spines['top'].set_linewidth(0.5)


ax=fig.add_subplot(333)
ax.plot(x, y_CH, label='C-H',linewidth=0.5)
ax.set_xlabel(r'r ($\rm \AA$)',font)
ax.set_ylabel('g(r)',font)
ax.legend(prop=font,frameon=False)
ax.minorticks_on()
ax.yaxis.set_minor_locator(plt.MultipleLocator(5))
ax.tick_params(axis='y',direction='in',which="major",width=0.5,length=3)
ax.tick_params(axis='y',direction='in',which="minor",width=0.5,length=1.5)

# ax.xaxis.set_minor_locator(plt.MultipleLocator(0.01))
ax.tick_params(axis='x',direction='in',which="major",width=0.5,length=3)
ax.tick_params(axis='x',direction='in',which="minor",width=0.5,length=1.5)
plt.tick_params(labelsize=6)
labels=ax.get_xticklabels()+ax.get_yticklabels()
[label.set_fontname('Arial') for label in labels]

#set line width
ax.spines['bottom'].set_linewidth(0.5)
ax.spines['left'].set_linewidth(0.5)
ax.spines['right'].set_linewidth(0.5)
ax.spines['top'].set_linewidth(0.5)


ax=fig.add_subplot(334)
ax.plot(x, y_OH, label='O-H',linewidth=0.5)
ax.set_xlabel(r'r ($\rm \AA$)',font)
ax.set_ylabel('g(r)',font)
ax.legend(prop=font,frameon=False)
ax.minorticks_on()
ax.yaxis.set_minor_locator(plt.MultipleLocator(5))
ax.tick_params(axis='y',direction='in',which="major",width=0.5,length=3)
ax.tick_params(axis='y',direction='in',which="minor",width=0.5,length=1.5)

# ax.xaxis.set_minor_locator(plt.MultipleLocator(0.01))
ax.tick_params(axis='x',direction='in',which="major",width=0.5,length=3)
ax.tick_params(axis='x',direction='in',which="minor",width=0.5,length=1.5)
plt.tick_params(labelsize=6)
labels=ax.get_xticklabels()+ax.get_yticklabels()
[label.set_fontname('Arial') for label in labels]

#set line width
ax.spines['bottom'].set_linewidth(0.5)
ax.spines['left'].set_linewidth(0.5)
ax.spines['right'].set_linewidth(0.5)
ax.spines['top'].set_linewidth(0.5)


ax=fig.add_subplot(335)
ax.plot(x, y_HH, label='H-H',linewidth=0.5)
ax.set_xlabel(r'r ($\rm \AA$)',font)
ax.set_ylabel('g(r)',font)
ax.legend(prop=font,frameon=False)
ax.minorticks_on()
ax.yaxis.set_minor_locator(plt.MultipleLocator(5))
ax.tick_params(axis='y',direction='in',which="major",width=0.5,length=3)
ax.tick_params(axis='y',direction='in',which="minor",width=0.5,length=1.5)

# ax.xaxis.set_minor_locator(plt.MultipleLocator(0.01))
ax.tick_params(axis='x',direction='in',which="major",width=0.5,length=3)
ax.tick_params(axis='x',direction='in',which="minor",width=0.5,length=1.5)
plt.tick_params(labelsize=6)
labels=ax.get_xticklabels()+ax.get_yticklabels()
[label.set_fontname('Arial') for label in labels]

#set line width
ax.spines['bottom'].set_linewidth(0.5)
ax.spines['left'].set_linewidth(0.5)
ax.spines['right'].set_linewidth(0.5)
ax.spines['top'].set_linewidth(0.5)


ax=fig.add_subplot(336)
ax.plot(x, y_CC, label='C-C',linewidth=0.5)
ax.set_xlabel(r'r ($\rm \AA$)',font)
ax.set_ylabel('g(r)',font)
ax.legend(prop=font,frameon=False)
ax.minorticks_on()
ax.yaxis.set_minor_locator(plt.MultipleLocator(5))
ax.tick_params(axis='y',direction='in',which="major",width=0.5,length=3)
ax.tick_params(axis='y',direction='in',which="minor",width=0.5,length=1.5)

# ax.xaxis.set_minor_locator(plt.MultipleLocator(0.01))
ax.tick_params(axis='x',direction='in',which="major",width=0.5,length=3)
ax.tick_params(axis='x',direction='in',which="minor",width=0.5,length=1.5)
plt.tick_params(labelsize=6)
labels=ax.get_xticklabels()+ax.get_yticklabels()
[label.set_fontname('Arial') for label in labels]

#set line width
ax.spines['bottom'].set_linewidth(0.5)
ax.spines['left'].set_linewidth(0.5)
ax.spines['right'].set_linewidth(0.5)
ax.spines['top'].set_linewidth(0.5)


ax=fig.add_subplot(337)
ax.plot(x, y_OO, label='O-O',linewidth=0.5)
ax.set_xlabel(r'r ($\rm \AA$)',font)
ax.set_ylabel('g(r)',font)
ax.legend(prop=font,frameon=False)
ax.minorticks_on()
ax.yaxis.set_minor_locator(plt.MultipleLocator(5))
ax.tick_params(axis='y',direction='in',which="major",width=0.5,length=3)
ax.tick_params(axis='y',direction='in',which="minor",width=0.5,length=1.5)

# ax.xaxis.set_minor_locator(plt.MultipleLocator(0.01))
ax.tick_params(axis='x',direction='in',which="major",width=0.5,length=3)
ax.tick_params(axis='x',direction='in',which="minor",width=0.5,length=1.5)
plt.tick_params(labelsize=6)
labels=ax.get_xticklabels()+ax.get_yticklabels()
[label.set_fontname('Arial') for label in labels]

#set line width
ax.spines['bottom'].set_linewidth(0.5)
ax.spines['left'].set_linewidth(0.5)
ax.spines['right'].set_linewidth(0.5)
ax.spines['top'].set_linewidth(0.5)
plt.show()


fig.subplots_adjust(wspace=0.3,hspace=0.5)

fig=plt.gcf()
cm=28
inch=0.3937007874
fig.set_size_inches(cm*inch,cm*inch*0.618)
plt.savefig('rdf.png',dpi=600,bbox_inches='tight',transparent=True)
