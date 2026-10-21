#!/home/dj/software/anaconda3_2019.07/envs/abipy/bin/python

import sys
import os
from abipy import abilab

abilab.enable_notebook(with_seaborn=True)

band = abilab.abiopen(sys.argv[1])
band.ebands.set_fermie_to_vbm()
a=band.ebands.plot(with_gaps=True,ylims=(-5,10),fontsize=12);

a.savefig(sys.argv[2],dpi=300)
