#!/home/dj/software/anaconda3_2019.07/envs/abipy/bin/python
import sys
import os
from abipy import abilab
import matplotlib as mpl 
mpl.use('Agg')

optic=abilab.abiopen(sys.argv[1])
lin_re=optic.plot_linear_epsilon(what='re')
#lin_re.remove()
lin_im=optic.plot_linear_epsilon(what='im')
#lin_im.remove()
shg=optic.plot_shg()
#shg.remove()

lin_re.savefig(sys.argv[2]+'_re.png',dpi=600)
lin_im.savefig(sys.argv[2]+'_im.png',dpi=600)
shg.savefig(sys.argv[2]+'_shg.png',dpi=600)
