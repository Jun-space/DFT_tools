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
from pymatgen.io.vasp.sets import VaspInputSet,MPRelaxSet,DictSet
from pymatgen.core.structure import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
import itertools
import numpy as np
import datetime
import time
from command import Command


def get_path():
    return os.getcwd()


def generate_input_opt(job_id):
    ###fetch task information and generate information for task
    client = pymongo.MongoClient("mongodb://localhost:27017")
    db = client['layers']
    col = db['vasptasks']
    query = {'job_id': job_id}
    data = col.find_one(query)
    s = Structure.from_dict(data['task_info']['structure'])
    path = data['run_info']['path']
    jobname = data['task_info']['structure_id']+'_opt'
    os.makedirs(path)

    try:
        ####generate INCAR POSCAR KPOINTS POTCAR
        user_incar_settings = {"SYSTEM": '%s' % jobname, "PREC": 'Accurate', "NPAR": 4, 'ISTART': 0, "ICHARG": 2,
                               "ALGO": 'Normal', "IALGO": 38, "ISMEAR": 0, "SIGMA": 0.05, "EDIFF": 1E-6,
                               "EDIFFG": -0.01, "ENCUT": 520, "NSW": 200, "IBRION": 2, "ISIF": 3, "ISPIN": 2,
                               "LDAU": '.FALSE.', "LREAL": '.FALSE.', "LWAVE": ".F.", "LCHARG": ".False"}
        mpvis = MPRelaxSet(s, user_incar_settings=user_incar_settings, potcar_functional="PBE_54")
        mpvis.poscar.write_file('%s/POSCAR' % path)
        mpvis.incar.write_file("%s/INCAR" % path)
        mpvis.kpoints.write_file("%s/KPOINTS" % path)
        mpvis.potcar.write_file('%s/POTCAR' % path)
        sub_script = ['#!/bin/sh', '#SBATCH -p CD', '#SBATCH -n 36', '#SBATCH -J ' + jobname,
                      'module load vasp-5.4.4.pl2.2019u5', 'echo Job starts at `date`>log', 'srun vasp_std >> log',
                      'echo Job ends at `date` >> log']
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

def create_vasptasks_opt(s_id):
    client=pymongo.MongoClient("mongodb://localhost:27017")
    db=client['layers']
    col=db['stacking_ionic']
    query={'structure_id' : s_id}
    data=col.find_one(query)
    start_job_id=db['vasptasks'].count_documents({})
    path=get_path() + '/' + str(s_id) + '/opt'
    job_id=start_job_id
    s=Structure.from_dict(data['structure_unrelax'])
    s=SpacegroupAnalyzer(s).get_primitive_standard_structure()
    s=s.as_dict()
    run_info={'slurm_job_id':None,\
                  'slurm_job_name':None,\
                      'status':1,\
                          'path':path,\
                              'start_time':None,\
                                  'end_time':None}
    task_info={'structure_id':data['structure_id'],\
                   'task_type': 'relax',\
                           'ads_id':data['ads_id'],\
                               'sub_id':data['sub_id'],\
                                   'structure':data['structure_id']}
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
    db['vasptasks'].insert(vasptask)
    return
