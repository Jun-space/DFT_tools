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
from phonopy.api_phonopy import Phonopy
from phonopy.interface.vasp import read_vasp,write_vasp,get_vasp_structure_lines
from pymatgen.ext.matproj import MPRester
from pymatgen.io.vasp.inputs import VaspInput
from pymatgen.io.vasp.inputs import Incar, Kpoints, Poscar,Potcar
from pymatgen.io.vasp.outputs import Outcar, Vasprun
from pymatgen.io.vasp.sets import VaspInputSet,MPRelaxSet,DictSet,MPStaticSet
from pymatgen.core.structure import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
import itertools
import numpy as np
import datetime
import time
from command import Command

def generate_input_ph_opt(job_id):
    ###fetch task information and generate information for task
    client = pymongo.MongoClient("mongodb://localhost:27017")
    db = client['layers']
    col = db['vasptasks']
    query = {'job_id': job_id}
    data = col.find_one(query)
    path = data['run_info']['path']
    s = Structure.from_dict(data['task_info']['structure'])
    jobname = data['task_info']['structure_id'] + '_ph'
    os.makedirs(path)
    if np.abs(data['task_info']['mag_opt']) <= 0.1:
        ispin = 1
    else:
        ispin = 2
    try:
        ####generate INCAR POSCAR KPOINTS POTCAR
        # user_incar_settings={"SYSTEM": '%s'%jobname,"PREC":'Accurate',"NPAR":4,'ISTART': 0, "ICHARG": 2,"ALGO":'Normal',"IALGO": 38, "ISMEAR" : 0, SIGMA":0.05,"EDIFF":1E-6,"EDIFFG":-0.001,"ENCUT":520,"NSW":200,"IBRION":2,"ISIF":3,"ISPIN":2,"LREAL":'.FALSE.',"LWAVE":".F.","LCHARG":".False","LDAU":'.TRUE.'}
        user_incar_settings = {"SYSTEM": '%s' % jobname, "PREC": 'Accurate', "NPAR": 4, 'ISTART': 0, "ICHARG": 2,
                               "ALGO": 'Normal', "IALGO": 38, "ISMEAR": 0, "SIGMA": 0.05, "EDIFF": 1E-8,
                               "EDIFFG": -0.0001, "ENCUT": 520, "NSW": 200, "IBRION": 2, "ISIF": 3, "ISPIN": ispin,
                               "LDAU": '.TRUE.', "LREAL": '.FALSE.', "LWAVE": ".F.", "LCHARG": ".False"}
        mpvis = MPRelaxSet(s, user_incar_settings=user_incar_settings, potcar_functional="PBE")
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

        new = {'$set': {'run_info.path': path, 'run_info.slurm_job_name': jobname, 'task_info.task_type': 'phonon-opt'}}
        col.update_one(query, new)
        return
    except:
        new = {'$set': {'run_info.status': -1}}
        col.update_one(query, new)


def generate_input_ph_scf(job_id):
    ###fetch task information and generate information for task
    client = pymongo.MongoClient("mongodb://localhost:27017")
    db = client['layers']
    col = db['vasptasks']
    query = {'job_id': job_id}
    data = col.find_one(query)
    path = data['run_info']['path']
    s = Structure.from_str(data['task_info']['structure'], fmt='poscar')
    # path=path=get_path()+'/scf'
    jobname = data['task_info']['structure_id'] + '_ph_' + str(data['task_info']['ph_config_id'])
    os.makedirs(path)
    if np.abs(data['task_info']['mag_opt']) <= 0.1:
        ispin = 1
    else:
        ispin = 2
    try:
        ####generate INCAR POSCAR KPOINTS POTCAR
        user_incar_settings = {"SYSTEM": '%s' % jobname, "PREC": 'Accurate', "NPAR": 4, 'ISTART': 0, "ICHARG": 2,
                               "ADDGRID": '.TRUE.', "ALGO": 'Normal', "IALGO": 38, "ISMEAR": 0, "SIGMA": 0.05,
                               "EDIFF": 1E-8, "EDIFFG": -0.0001, "ENCUT": 520, "NSW": 0, "IBRION": -1, "ISIF": 3,
                               "ISPIN": ispin, "LDAU": '.TRUE.', "LREAL": '.FALSE.', "LWAVE": ".F.", "LCHARG": ".False",
                               'LAECHG': 'F'}
        mpvis = MPStaticSet(s, user_incar_settings=user_incar_settings, potcar_functional="PBE", reciprocal_density=400)
        # mpvis.poscar.write_file('%s/POSCAR'%path)
        mpvis.incar.write_file("%s/INCAR" % path)
        mpvis.kpoints.write_file("%s/KPOINTS" % path)
        mpvis.potcar.write_file('%s/POTCAR' % path)
        sub_script = ['#!/bin/sh', '#SBATCH -n 224', '#SBATCH -J ' + jobname, 'module load intel_parallel',
                      'module load vasp-avx512', 'unset I_MPI_PMI_LIBRARY', 'ulimit -s unlimited',
                      'echo Job starts at `date`>log', 'mpirun -np $SLURM_NPROCS vasp_std >> log',
                      'echo Job ends at `date` >> log']
        if os.path.exists('%s/runvasp.pbs' % path) == False:
            os.mknod('%s/runvasp.pbs' % path)
        fp = open('%s/runvasp.pbs' % path, 'w')
        for i in sub_script:
            fp.write(str(i) + '\n')
        fp.close()
        # write poscar
        fp = open('%s/POSCAR' % path, 'w')
        fp.write(data['task_info']['structure'])
        print('generated input file for job %d' % job_id)

        new = {'$set': {'run_info.path': path, 'run_info.slurm_job_name': jobname, 'task_info.task_type': 'phonon'}}
        col.update_one(query, new)
        return
    except:
        new = {'$set': {'run_info.status': -1}}
        col.update_one(query, new)


def create_vasptasks_ph_opt(s_id):
    client = pymongo.MongoClient("mongodb://localhost:27017")
    db = client['layers']
    col = db['stacking_ionic']
    query = {'structure_id': s_id}
    data = col.find_one(query)
    start_job_id = db['vasptasks'].count_documents({})
    vasptask = {}
    job_id = start_job_id
    path = data['path']+'/phonon/opt'
    run_info = {'slurm_job_id': None, \
                'slurm_job_name': None, \
                'status': 1, \
                'path': path, \
                'start_time': None, \
                'end_time': None}
    task_info = {'structure_id': data['structure_id'], \
                 'task_type': 'phonon-opt', \
                 'ads_id': data['ads_id'], \
                 'sub_id': data['sub_id'], \
                 'structure': data['structure'], \
                 'mag_scf': data['mag_scf']}
    post_info = {'ele_converge': None, \
                 'ionic_converge': None, \
                 'efermi': None, \
                 'etot': None, \
                 'contcar': None, \
                 'band': None, \
                 'dos': None, \
                 'tot_mag': None}
    vasptask = {'job_id': job_id, \
                'run_info': run_info, \
                'task_info': task_info, \
                'post_info': post_info}
    db['vasptasks'].insert_one(vasptask)
    return


def create_vasptasks_phonon_displacement(s_id):
    client = pymongo.MongoClient("mongodb://localhost:27017")
    db = client['layers']
    col = db['stacking_ionic']
    query = {'structure_id': s_id}
    data = col.find_one(query)
    start_job_id = db['vasptasks'].count_documents({})
    vasptasks = []
    for num, i in enumerate(data['phonon']['structures']):
        job_id = start_job_id + num
        path=data['path']+ '/phonon/' + str(num)
        run_info = {'slurm_job_id': None, \
                    'slurm_job_name': None, \
                    'status': 1, \
                    'path': path, \
                    'start_time': None, \
                    'end_time': None}
        task_info = {'structure_id': data['structure_id'], \
                     'task_type': 'phonon', \
                     'ph_config_id': num, \
                     'ads_id': data['ads_id'], \
                     'sub_id': data['sub_id'], \
                     'structure': i, \
                     'mag_opt': data['mag_opt']}
        post_info = {'ele_converge': None, \
                     'ionic_converge': None, \
                     'efermi': None, \
                     'etot': None, \
                     'contcar': None, \
                     'band': None, \
                     'dos': None, \
                     'tot_mag': None, \
                     'mag_ele': None}
        vasptasks.append({'job_id': job_id, \
                          'run_info': run_info, \
                          'task_info': task_info, \
                          'post_info': post_info})
    db['vasptasks'].insert_many(vasptasks)
    return


def create_phonon_displacement(s_id):
    client = pymongo.MongoClient("mongodb://localhost:27017")
    db = client['layers']
    col = db['stacking_ionic']
    query = {'structure_id': s_id}
    data = col.find_one(query)
    path = data['path']
    unitcell = read_vasp(path + '/phonon/opt/CONTCAR')
    lattice = unitcell.get_cell()
    a = np.linalg.norm(lattice[0])
    b = np.linalg.norm(lattice[1])
    c = np.linalg.norm(lattice[2])
    na, nb, nc = int(np.ceil(15 / a)), int(np.ceil(15 / b)), int(np.ceil(15 / c))

    sc_matrix = [[na, 0, 0], [0, nb, 0], [0, 0, nc]]
    phonon = Phonopy(unitcell, supercell_matrix=sc_matrix)
    phonon.generate_displacements(distance=0.01)
    supercells = phonon.supercells_with_displacements

    # for i in range(len(supercells)):
    #     write_vasp('POSCAR-'+str(i),supercells[i])

    structures = []
    for i in range(len(supercells)):
        s_list = get_vasp_structure_lines(supercells[i])
        s_str = str()
        for j in s_list:
            s_str = s_str + j + '\n'
        structures.append(s_str)
    new = {'$set': {'phonon.structures': structures, 'phonon.dim': [na, nb, nc]}}
    col.update_one(query, new)
    return