#!/home/dj/software/anaconda3_2019.07/envs/abipy/bin/python
import sys
import os
from abipy import abilab

optic=abilab.abiopen(sys.argv[1])
lin_re=optic.plot_linear_epsilon(what='re')
lin_im=optic.plot_linear_epsilon(what='im')
shg=optic.plot_shg()

lin_re.savefig(sys.argv[2]+'_re.png',dpi=300)
lin_im.savefig(sys.argv[2]+'_im.png',dpi=300)
shg.savefig(sys.argv[2]+'_shg.png',dpi=300)
