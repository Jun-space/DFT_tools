# -*- coding: utf-8 -*-
"""
Created on Tue Dec 26 09:47:04 2023

@author: DengJun
"""


from pymatgen.analysis.adsorption import *
from pymatgen.core.surface import Slab, SlabGenerator, generate_all_slabs, Structure, Lattice, ReconstructionGenerator
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.core.structure import Structure
from pymatgen.ext.matproj import MPRester
from pymatgen.io.vasp.sets import VaspInputSet,MPRelaxSet,DictSet,MPStaticSet
from pymatgen.io.vasp.inputs import Incar, Kpoints, Poscar,Potcar
from matplotlib import pyplot as plt
from pymatgen.io.vasp.inputs import Poscar
import sys,os
from kpoints import kpoints,kpoints_2d
from command import Command

def get_path():
    return os.getcwd()

def generate_input_scf(file,prefix,path):
    ###fetch task information and generate information for task
    s = Structure.from_file(file)
    jobname = prefix+'-scf'
    os.makedirs(path)
    current_path=get_path()
    try:
        ####generate INCAR POSCAR KPOINTS POTCAR
        user_incar_settings = {"SYSTEM": '%s' % jobname, "PREC": 'Accurate', "NPAR": 4, 'ISTART': 0, "ICHARG": 2,
                       "ALGO": 'Normal', "IALGO": 38, "ISMEAR": -5, "EDIFF": 1E-6, "EDIFFG": -0.005,\
                       "ENCUT": 520, "NSW": 0, "IBRION": -1, "ISIF": 3, "ISPIN": 2, "LDAU": '.FALSE.',\
                       "LREAL": '.FALSE.', "LWAVE":".F.","LCHARG":".F.",'LAECHG':'F','LVHAR':'F'}
        mpvis = MPStaticSet(s, user_incar_settings=user_incar_settings, user_potcar_functional="PBE_54", reciprocal_density=350)
        mpvis.poscar.write_file('%s/POSCAR' % path)
        mpvis.incar.write_file("%s/INCAR" % path)
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
        if os.path.exists('%s/runvasp.pbs' % path) == False:
            os.mknod('%s/runvasp.pbs' % path)
        fp = open('%s/runvasp.pbs' % path, 'w')
        for i in sub_script:
            fp.write(str(i) + '\n')
        fp.close()
        os.chdir(path)
        tmp=os.popen('sbatch runvasp.pbs').readlines()[0].split()[1]
        os.chdir(current_path)
    except:
        print('fail '+prefix)
    return

if __name__ == "__main__":
    filenames=os.listdir('./structures')
    for filename in filenames[27:]:
        prefix=filename.split('.')[0]
        fpath=prefix+'/opt/CONTCAR'
        path=prefix+'/scf'
        generate_input_scf(fpath,prefix,path)
        print(prefix)
