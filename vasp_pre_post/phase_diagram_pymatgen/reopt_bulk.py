#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 21 12:13:15 2022

@author: dj
"""

import os,sys
import re
import shutil
from pymatgen.ext.matproj import MPRester
from pymatgen.io.vasp.inputs import VaspInput
from pymatgen.io.vasp.inputs import Incar, Kpoints, Poscar,Potcar
from pymatgen.io.vasp.outputs import Outcar, Vasprun
from pymatgen.io.vasp.sets import VaspInputSet,MPRelaxSet,DictSet
from pymatgen.core.structure import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
import itertools
import numpy as np
import datetime
import time
from command import Command
from kpoints import kpoints,kpoints_2d

def get_path():
    return os.getcwd()

def generate_input_opt(fpath,filename,prefix):
    ###fetch task information and generate information for task
    s = Structure.from_file(fpath)
    #s = SpacegroupAnalyzer(s).get_primitive_standard_structure()
    jobname = prefix+'-opt'
    path=get_path()+'/'+prefix+'/opt'
    #os.makedirs(path)
    current_path=get_path()
    try:
        ####generate INCAR POSCAR KPOINTS POTCAR
        user_incar_settings = {"SYSTEM": '%s' % jobname, "PREC": 'Accurate', "NPAR": 4, 'ISTART': 0, "ICHARG": 2,"ALGO":'Normal',\
                "IALGO": 38, "ISMEAR": 1, "SIGMA": 0.05, "EDIFF": 1E-6,"EDIFFG": -0.01, "ENCUT": 520, "NSW": 200,\
                        "IBRION": 2, "ISIF": 3, "ISPIN": 2,"LDAU":'.FALSE.', "LREAL": '.FALSE.', "LWAVE": ".F.",\
                        "LCHARG":".False"}
        mpvis = MPRelaxSet(s, user_incar_settings=user_incar_settings, user_potcar_functional="PBE_54")
        mpvis.poscar.write_file('%s/POSCAR' % path)
        #mpvis.incar.write_file("%s/INCAR" % path)
        #mpvis.kpoints.write_file("%s/KPOINTS" % path)
        kpoints(0.02,path)
        mpvis.potcar.write_file('%s/POTCAR' % path)
        sub_script = ['#!/bin/sh', '#SBATCH -p amd9654', '#SBATCH -n 48', '#SBATCH -J ' + jobname,\
                '#SBATCH -N 1',\
                'source /public/env/compiler_intel2025',\
                'export export PATH=/public/home/gouhuiyang/software/vasp_5.4.4_standard_intel/bin:$PATH',\
                'ulimit -s unlimited',\
                'echo Job starts at `date`>log',\
                'mpirun -np $SLURM_NPROCS vasp_std >> log',\
                'echo Job ends at `date` >> log']
        #if os.path.exists('%s/runvasp.pbs' % path) == False:
        #    os.mknod('%s/runvasp.pbs' % path)
        #fp = open('%s/runvasp.pbs' % path, 'w')
        #for i in sub_script:
        #    fp.write(str(i) + '\n')
        #fp.close()
        os.chdir(path)
        tmp=os.popen('sbatch runvasp.pbs').readlines()[0].split()[1]
        os.chdir(current_path)
    except:
        print('fail '+filename)
    return

if __name__ == '__main__':
    filenames=os.listdir('./structures')
    fpath='./structures'
    for filename in filenames:
        prefix=filename.split('.')[0]
        fpath=prefix+'/opt/CONTCAR'
        generate_input_opt(fpath,filename,prefix)
        print(prefix)
