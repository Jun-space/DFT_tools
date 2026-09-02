#!/publicfs10/fs10-m9/home/m9s005343/software/anaconda3/bin/python

import math
import sys
import os
import yaml
from pymatgen.io.vasp import Potcar
from pymatgen.io.vasp import Poscar

def potcar_symbols(s,sets):
    elements=s.site_symbols
    potcar_symbols = []
    settings = sets["POTCAR"]

    if isinstance(settings[elements[-1]],dict):
        for el in elements:
            potcar_symbols.append(settings[el]["symbol"] if el in settings else el)
    else:
        for el in elements:
            potcar_symbols.append(settings.get(el,el))
    return potcar_symbols

s=Poscar.from_file('POSCAR')

f=open('/publicfs10/fs10-m9/home/m9s005343/software/anaconda3/lib/python3.13/site-packages/pymatgen/io/vasp/MPRelaxSet.yaml','r')
#f=open('/public/home/gouhuiyang/software/anaconda3/lib/python3.8/site-packages/pymatgen/io/vasp/MPRelaxSet.yaml','r')
sets=yaml.safe_load(f)
f.close()

symbols=potcar_symbols(s,sets)

potcar=Potcar(symbols,functional='PBE_54')
#potcar=Potcar(symbols)
potcar.write_file('POTCAR')
