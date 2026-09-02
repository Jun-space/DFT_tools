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
from pymatgen.io.vasp.sets import VaspInputSet,MPRelaxSet,DictSet,MPStaticSet
from pymatgen.analysis.magnetism import CollinearMagneticStructureAnalyzer
from pymatgen.core.structure import Structure
import itertools
import numpy as np
import datetime
import time
import copy
from command import Command

def generate_input_mag_opt(job_id):
    ###fetch task information and generate information for task
    client = pymongo.MongoClient("mongodb://localhost:27017")
    db = client['layers']
    col = db['vasptasks']
    query = {'job_id': job_id}
    data = col.find_one(query)
    s_path = data['run_info']['path']
    s = Structure.from_dict(data['task_info']['structure'])
    path = s_path
    # path=path=get_path()+'/scf'
    jobname = data['task_info']['structure_id'] + '_mag_' + str(data['task_info']['mag_config_id']) + '_opt'
    os.makedirs(path)

    try:
        ####generate INCAR POSCAR KPOINTS POTCAR
        user_incar_settings = {"SYSTEM": '%s' % jobname, "PREC": 'Accurate', "NPAR": 4, 'ISTART': 0, "ICHARG": 2,
                               "ALGO": 'Normal', "IALGO": 38, "ISMEAR": 0, "SIGMA": 0.05, "EDIFF": 1E-6,
                               "EDIFFG": -0.01, "ENCUT": 520, "NSW": 200, "IBRION": 2, "ISIF": 3, "ISPIN": 2,
                               "LDAU": '.TRUE.', "LREAL": '.FALSE.', "LWAVE": ".F.", "LCHARG": ".False", 'LAECHG': 'F'}
        mpvis = MPStaticSet(s, user_incar_settings=user_incar_settings, potcar_functional="PBE", reciprocal_density=350)
        mpvis.poscar.write_file('%s/POSCAR' % path)
        mpvis.incar.write_file("%s/INCAR" % path)
        mpvis.kpoints.write_file("%s/KPOINTS" % path)
        mpvis.potcar.write_file('%s/POTCAR' % path)
        sub_script = ['#!/bin/sh', '#SBATCH -p CD', '#SBATCH -n 144', '#SBATCH -J ' + jobname,
                      'module load vasp-5.4.4.pl2.2019u5', 'echo Job starts at `date`>log', 'srun vasp_std >> log',
                      'echo Job ends at `date` >> log']
        if os.path.exists('%s/runvasp.pbs' % path) == False:
            os.mknod('%s/runvasp.pbs' % path)
        fp = open('%s/runvasp.pbs' % path, 'w')
        for i in sub_script:
            fp.write(str(i) + '\n')
        fp.close()
        print('generated input file for job %d' % job_id)

        new = {'$set': {'run_info.path': path, 'run_info.slurm_job_name': jobname, 'task_info.task_type': 'mag-opt'}}
        col.update_one(query, new)
        return
    except:
        new = {'$set': {'run_info.status': -1}}
        col.update_one(query, new)

def create_vasptasks_mag_opt(s_id):
    client=pymongo.MongoClient("mongodb://localhost:27017")
    db=client['layers']
    col=db['stacking_ionic']
    query={'structure_id' : s_id}
    data=col.find_one(query)
    start_job_id=db['vasptasks'].count_documents({})
    vasptasks=[]
    path=data['path']
    for num,i in enumerate(data['mag_configs']):
        job_id=start_job_id+num
        path=data['path']+ '/GS/' + str(num) + '/opt'
        run_info={'slurm_job_id':None,\
                  'slurm_job_name':None,\
                      'status':1,\
                          'path':path,\
                              'start_time':None,\
                                  'end_time':None}
        task_info={'structure_id':data['structure_id'],\
                   'task_type': 'mag-opt',\
                       'mag_config_id':num,\
                           'ads_id':data['ads_id'],\
                               'sub_id':data['sub_id'],\
                                   'structure':i}
        post_info={'ele_converge':None,\
                   'ionic_converge':None,\
                       'efermi':None,\
                           'etot':None,\
                               'contcar':None,\
                                   'band':None,\
                                       'dos':None,\
                                           'tot_mag':None,\
                                               'mag_ele':None}
        vasptasks.append({'job_id':job_id,\
                          'run_info':run_info,\
                         'task_info':task_info,\
                             'post_info':post_info})
    db['vasptasks'].insert_many(vasptasks)
    return

def create_mag_structures(s_id):
    client = pymongo.MongoClient("mongodb://localhost:27017")
    db = client['layers']
    col = db['stacking_ionic']
    query = {'structure_id': s_id}
    ##fetch data
    data = col.find_one(query)
    s = Structure.from_dict(data['structure_opt'])  # get structure
    site_symbol = [i.name for i in s.species]  # get site symbols
    mag_ele = [i['tot'] for i in data['mag_ele']]  # get magnetic moment for each site

    ####get spin state for each element according to mag_ele
    spin = {}
    for i in range(len(site_symbol)):
        spin[site_symbol[i]] = mag_ele[i]
    mag = []
    site_symbol = [i.name for i in s.species]
    for i, site in enumerate(site_symbol):
        if spin[site] >= 0.5:
            mag.append(i)
    mag_tmp = [str(i) for i in mag]

    ####### combination of afm structures
    """
    Note: if the number of magnetic atoms are odd, 
    we build a 2*2*2 supercell
    """
    if len(mag_tmp) % 2 == 1:
        s.make_supercell([2, 2, 2])
        mag = []
        site_symbol = [i.name for i in s.species]
        for i, site in enumerate(site_symbol):
            if spin[site] >= 0.5:
                mag.append(i)
        mag_tmp = [str(i) for i in mag]
    else:
        s.make_supercell([1, 1, 2])
        mag = []
        site_symbol = [i.name for i in s.species]
        for i, site in enumerate(site_symbol):
            if spin[site] >= 0.5:
                mag.append(i)
        mag_tmp = [str(i) for i in mag]

    combinations = []
    for i in itertools.combinations(mag_tmp, int(len(mag_tmp) / 2)):
        combinations.append(i)

    remove = []
    for num1, i in enumerate(combinations):
        for num2, j in enumerate(combinations):
            tmp = list(i + j)
            if list(set(mag_tmp).difference(set(tmp))) == []:
                combinations.pop(num2)

    ###############get FM and AFM structures
    fm = []
    for i in site_symbol:
        fm.append(spin[i])

    ###############AFM
    afm = []
    for i in combinations:
        spins = []
        tmp = np.ones(len(site_symbol))
        for k in list(i):
            tmp[int(k)] = tmp[int(k)] * -1
        for numb, j in enumerate(site_symbol):
            spins.append(spin[j])
        spins = list(np.array(spins) * tmp)
        afm.append(spins)

    mag_config = [fm] + afm

    mag_structures = []
    for i in mag_config:
        s_tmp = copy.deepcopy(s)
        s_tmp.add_spin_by_site(i)
        mag_structures.append(s_tmp)

    ########remove duplicated magnetic structures
    def compare_two(s1, s2):
        mat = CollinearMagneticStructureAnalyzer(s1)
        return mat.matches_ordering(s2)

    d = np.zeros([len(mag_structures), len(mag_structures)])
    for i in range(len(mag_structures)):
        for j in range(i, len(mag_structures)):
            if i == j:
                continue
            else:
                if compare_two(mag_structures[i], mag_structures[j]):
                    d[i, j] = 1
                    d[j, i] = 1

    def list_combine(l):
        if len(l) == 1:
            return list(set(l[0]))
        else:
            t = []
            for i in l:
                t = t + i
            return list(set(t))

    def find_match(d, match):
        tmp = []
        for i in match[-1]:
            link = []
            link.append(i)
            for j in range(len(d)):
                if d[i][j] == 1:
                    link.append(j)
            tmp = tmp + link
        new = list(set(tmp).difference(set(list_combine(match))))
        if new == []:
            return match
        match.append(new)
        return find_match(d, match)

    matches = []
    for i in range(len(mag_structures)):
        match = []
        if (i in list_combine(matches)) == False:
            match.append([i])
            match = find_match(d, match)
            matches.append(list_combine(match))
    keep = []
    for i in matches:
        keep.append(i[0])
    remove = list(set(np.arange(len(mag_structures))).difference(set(keep)))
    remove.sort(reverse=True)
    for i in remove:
        mag_structures.pop(i)
    for num, i in enumerate(mag_structures):
        mag_structures[num] = i.as_dict()
    new = {'$set': {'mag_configs': mag_structures}}
    col.update_one(query, new)
    return