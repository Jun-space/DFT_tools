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
from pymatgen.core.structure import Structure
import itertools
import numpy as np
import datetime
import time
from command import Command
from opt import generate_input_opt
from scf import generate_input_scf
from band import generate_input_band
from magnetic_structure_opt import generate_input_mag_opt
from phonon import generate_input_ph_opt,generate_input_ph_scf


def submitTask(job_id):
    client = pymongo.MongoClient("mongodb://localhost:27017")
    db = client['layers']
    col = db['vasptasks']
    query = {'job_id': job_id}
    data = col.find_one(query)
    pwd = get_path()
    path = data['run_info']['path']
    try:
        os.chdir(path)
        slurm_job_id = int(os.popen('sbatch runvasp.pbs').readlines()[0].split()[3])
        print('submitted job %d' % job_id)
        status = 2
        new = {'$set': {'run_info.slurm_job_id': slurm_job_id, 'run_info.status': status}}
        col.update_one(query, new)
        os.chdir(pwd)
    except:
        print('submitted job %d failed' % job_id)
        new = {'$set': {'run_info.status': -2}}
        col.update_one(query, new)


def checkTaskStatus(job_id):
    '''
    Check whether the task is complete or not.
    If it is completed, update the status and post-process it.
    '''
    client = pymongo.MongoClient("mongodb://localhost:27017")
    db = client['layers']
    col = db['vasptasks']
    query = {'job_id': job_id}
    data = col.find_one(query)
    slurm_job_id = data['run_info']['slurm_job_id']
    ###check whether the job is completed
    c = Command('squeue -j %d' % slurm_job_id)
    c.run(10)
    if str(slurm_job_id) in c.output.split():
        #        print('job %d is still running with slurm job id %d'%(job_id,slurm_job_id))
        return
    else:
        time.sleep(2)
        c = Command("sacct -n -P --format=jobid,start,end -j %d" % slurm_job_id)
        c.run(10)
        tmp = c.output.split()[0].split("|")
        start_time = datetime.datetime.strptime(tmp[1], "%Y-%m-%dT%H:%M:%S").strftime("%Y/%m/%d %H:%M:%S")
        end_time = datetime.datetime.strptime(tmp[2], "%Y-%m-%dT%H:%M:%S").strftime("%Y/%m/%d %H:%M:%S")
        status = 0
        new = {'$set': {'run_info.start_time': start_time, 'run_info.end_time': end_time, 'run_info.status': status}}
        col.update_one(query, new)
        if data['task_info']['task_type'] == 'relax':
            post_process_opt(job_id)
        elif data['task_info']['task_type'] == 'scf':
            post_process_scf(job_id)
        elif data['task_info']['task_type'] == 'band':
            post_process_band(job_id)
        elif data['task_info']['task_type'] == 'mag-opt':
            post_process_opt(job_id)
        if data['task_info']['task_type'] == 'phonon-opt':
            post_process_opt(job_id)
        if data['task_info']['task_type'] == 'phonon':
            post_process_scf(job_id)

def post_process_opt(job_id):
    client=pymongo.MongoClient("mongodb://localhost:27017")
    db=client['layers']
    col=db['vasptasks']
    query={'job_id' : job_id}
    data=col.find_one(query)
    path=data['run_info']['path']
    try:
        vasprun=Vasprun('%s/vasprun.xml'%path)
        outcar=Outcar('%s/OUTCAR'%path)
        efermi=vasprun.efermi
        etot=vasprun.final_energy
        contcar=vasprun.final_structure.as_dict()
        ele_converge=vasprun.converged_electronic
        ionic_converge=vasprun.converged_ionic
        mag_tot=outcar.total_mag
        new={'$set':{'post_info.efermi':efermi,\
                     'post_info.etot':etot,\
                         'post_info.ele_converge':ele_converge,\
                             'post_info.ionic_converge':ionic_converge,\
                                 'post_info.tot_mag': mag_tot,\
                                     'post_info.contcar':contcar}}
        col.update_one(query,new)
        print('job %d is completed and parsed'%(job_id))
    except:
        print('job %d is completed but parsed with errors'%(job_id))
        new={'$set':{'run_info.status':-3}}
        col.update_one(query,new)

def post_process_scf(job_id):
    client=pymongo.MongoClient("mongodb://localhost:27017")
    db=client['layers']
    col=db['vasptasks']
    query={'job_id' : job_id}
    data=col.find_one(query)
    path=data['run_info']['path']
    try:
        vasprun=Vasprun('%s/vasprun.xml'%path)
        outcar=Outcar('%s/OUTCAR'%path)
        efermi=vasprun.efermi
        etot=vasprun.final_energy
        contcar=vasprun.final_structure.as_dict()
        ele_converge=vasprun.converged_electronic
        ionic_converge=vasprun.converged_ionic
        mag_tot=outcar.total_mag
        new={'$set':{'post_info.efermi':efermi,\
                     'post_info.etot':etot,\
                         'post_info.ele_converge':ele_converge,\
                             'post_info.ionic_converge':ionic_converge,\
                                 'post_info.tot_mag': mag_tot,\
                                     'post_info.contcar':contcar,\
                                         'post_info.mag_ele':outcar.magnetization}}
        col.update_one(query,new)
        print('job %d is completed and parsed'%(job_id))
    except:
        print('job %d is completed but parsed with errors'%(job_id))
        new={'$set':{'run_info.status':-3}}
        col.update_one(query,new)


def post_process_band(job_id):
    client=pymongo.MongoClient("mongodb://localhost:27017")
    db=client['layers']
    col=db['vasptasks']
    query={'job_id' : job_id}
    data=col.find_one(query)
    path=data['run_info']['path']
    try:
        vasprun=Vasprun('%s/vasprun.xml'%path,parse_projected_eigen=True)
        # outcar=Outcar('%s/OUTCAR'%path)
        efermi=vasprun.efermi
        # etot=vasprun.final_energy
        # contcar=vasprun.final_structure.as_dict()
        ele_converge=vasprun.converged_electronic
        ionic_converge=vasprun.converged_ionic
        bs_data=vasprun.get_band_structure(line_mode=1).as_dict()
        # mag_tot=outcar.total_mag
        new={'$set':{'post_info.efermi':efermi,\
                     'post_info.etot': None,\
                         'post_info.ele_converge':ele_converge,\
                             'post_info.ionic_converge':ionic_converge,\
                                 'post_info.tot_mag': None,\
                                     'post_info.band': path+'/vasprun.xml',\
                                     'post_info.contcar':None}}
        col.update_one(query,new)
        print('job %d is completed and parsed'%(job_id))
    except:
        print('job %d is completed but parsed with errors'%(job_id))
        new={'$set':{'run_info.status':-3}}
        col.update_one(query,new)


def getRunningTasks():
    client = pymongo.MongoClient("mongodb://localhost:27017")
    db = client['layers']
    col = db['vasptasks']
    query = {'run_info.status': 2}
    data = col.find(query, {'job_id': 1})
    running_ids = []
    for i in data:
        running_ids.append(i['job_id'])
    return running_ids


def getWaitingTasks():
    client = pymongo.MongoClient("mongodb://localhost:27017")
    db = client['layers']
    col = db['vasptasks']
    query = {'run_info.status': 1}
    data = col.find(query, {'job_id': 1})
    waiting_ids = []
    for i in data:
        waiting_ids.append(i['job_id'])
    return waiting_ids


def Run(tasklimit, timeinterval):
    count = 0
    while True:
        ####check the status of running tasks and update its status
        running_ids = getRunningTasks()
        for job_id in running_ids:
            checkTaskStatus(job_id)
        ###get waitting tasks and submit it
        running_ids = getRunningTasks()
        waiting_ids = getWaitingTasks()
        if len(running_ids) < tasklimit and len(waiting_ids) != 0:
            n = tasklimit - len(running_ids)
            if len(waiting_ids) < n:
                for i in range(len(waiting_ids)):
                    generate_input_opt(waiting_ids[i])
                    submitTask(waiting_ids[i])
            else:
                for i in range(n):
                    generate_input_opt(waiting_ids[i])
                    submitTask(waiting_ids[i])
        count = count + 1
        print("count: " + str(count) + "   " + time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(time.time())))
        time.sleep(timeinterval)

if __name__ == "__main__":
    Run(tasklimit=12, timeinterval=10)