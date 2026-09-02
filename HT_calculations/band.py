#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 21 12:13:15 2022

@author: dj
"""
import pymongo
import os,sys
import re
import shutil
from pymatgen.ext.matproj import MPRester
from pymatgen.io.vasp.inputs import VaspInput
from pymatgen.io.vasp.inputs import Incar, Kpoints, Poscar,Potcar
from pymatgen.io.vasp.outputs import Outcar, Vasprun
from pymatgen.io.vasp.sets import VaspInputSet,MPRelaxSet,DictSet,MPStaticSet,MPNonSCFSet
from pymatgen.symmetry.bandstructure import HighSymmKpath
from pymatgen.core.structure import Structure
import itertools
import numpy as np
import datetime
import time
from command import Command


def generate_input_band(job_id):
    ###fetch task information and generate information for task
    client = pymongo.MongoClient("mongodb://localhost:27017")
    db = client['layers']
    col = db['vasptasks']
    query = {'job_id': job_id}
    data = col.find_one(query)
    path = data['run_info']['path']
    s = Structure.from_dict(data['task_info']['structure'])
    jobname = data['task_info']['structure_id'] + '_band'
    os.makedirs(path)
    if np.abs(data['task_info']['mag_opt']) <= 0.1:
        ispin = 1
    else:
        ispin = 2

    try:
        ####generate INCAR POSCAR KPOINTS POTCAR
        user_incar_settings = {"SYSTEM": jobname, "PREC": 'Accurate', 'ISTART': 1, "NPAR": 4, "ALGO": 'Normal',
                               "IALGO": 38, "EDIFF": 1E-6, "EDIFFG": -0.01, "ENCUT": 520, "NSW": 0, "IBRION": -1,
                               "ISIF": 3, "ISPIN": ispin, "LDAU": '.TRUE.', "LREAL": '.FALSE.'}
        mpvis = MPNonSCFSet(s, mode='line', kpoints_line_density=40, user_incar_settings=user_incar_settings,
                            potcar_functional="PBE", reciprocal_density=400)
        # kpath=HighSymmKpath(s)
        # kpts=Kpoints.automatic_linemode(divisions=30,ibz=kpath)
        mpvis.poscar.write_file('%s/POSCAR' % path)
        mpvis.incar.write_file("%s/INCAR" % path)
        mpvis.kpoints.write_file("%s/KPOINTS" % path)
        # kpts.write_file("%s/KPOINTS"%path)
        mpvis.potcar.write_file('%s/POTCAR' % path)
        sub_script = ['#!/bin/sh', '#SBATCH -p CD', '#SBATCH -n 36', '#SBATCH -J ' + jobname,
                      'module load vasp-5.4.4.pl2.2019u5', 'ln -s ../scf/CHGCAR ./', 'ln -s ../scf/WAVECAR ./',
                      'echo Job starts at `date`>log', 'srun vasp_std >> log', 'echo Job ends at `date` >> log']
        if os.path.exists('%s/runvasp.pbs' % path) == False:
            os.mknod('%s/runvasp.pbs' % path)
        fp = open('%s/runvasp.pbs' % path, 'w')
        for i in sub_script:
            fp.write(str(i) + '\n')
        fp.close()
        print('generated input file for job %d' % job_id)

        new = {'$set': {'run_info.slurm_job_name': jobname}}
        col.update_one(query, new)
        return
    except:
        new = {'$set': {'run_info.status': -1}}
        col.update_one(query, new)


def create_vasptasks_band(s_id):
    client=pymongo.MongoClient("mongodb://localhost:27017")
    db=client['layers']
    col=db['stacking_ionic']
    query={'structure_id' : s_id}
    data=col.find_one(query)
    start_job_id=db['vasptasks'].count_documents({})
    job_id=start_job_id
    path=data['path']+'/band'
    run_info={'slurm_job_id':None,\
              'slurm_job_name':None,\
                  'status':1,\
                      'path':path,\
                          'start_time':None,\
                              'end_time':None}
    task_info={'structure_id':data['structure_id'],\
               'task_type': 'band',\
                       'ads_id':data['ads_id'],\
                           'sub_id':data['sub_id'],\
                               'structure':data['structure_opt'],
                               'mag_opt':data['mag_opt']}
    post_info={'ele_converge':None,\
               'ionic_converge':None,\
                   'efermi':None,\
                       'etot':None,\
                           'contcar':None,\
                               'band':None,\
                                   'dos':None,\
                                       'tot_mag':None}
    vasptask={'job_id':job_id,\
                      'run_info':run_info,\
                     'task_info':task_info,\
                         'post_info':post_info}
    db['vasptasks'].insert_one(vasptask)
    return