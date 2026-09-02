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

client_name='tasks'
hostname='ln3'

def get_path():
    return os.getcwd()

def submitTask(job_id):
    client = pymongo.MongoClient("mongodb://%s:27018"%hostname)
    db = client[client_name]
    col = db['vasptasks']
    query = {'job_id': job_id}
    data = col.find_one(query)
    pwd = get_path()
    path = data['run_info']['path']
    try:
        os.chdir(path)
        slurm_job_id = int(os.popen('sbatch runvasp.pbs').readlines()[0].split()[3])
        #tmp=os.popen('bsub < runvasp.pbs').readlines()[0].split()[1]
        #slurm_job_id = int(re.split('<|>',tmp)[1])
        print('submitted job %d with slurm_job_id %d' % (job_id, slurm_job_id))
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
    client = pymongo.MongoClient("mongodb://%s:27018"%hostname)
    db = client[client_name]
    col = db['vasptasks']
    query = {'job_id': job_id}
    data = col.find_one(query)
    slurm_job_id = data['run_info']['slurm_job_id']
    ###check whether the job is completed
    c = Command('squeue -j %d' % slurm_job_id)
    c.run(30)
    if str(slurm_job_id) in c.output.split():
        #        print('job %d is still running with slurm job id %d'%(job_id,slurm_job_id))
        return
    else:
        import time
        time.sleep(30)
        #obtain time for slurm system
        c = Command("sacct -n -P --format=jobid,start,end -j %d" % slurm_job_id)
        #c = Command("bhist -l 2991 | grep Starting")
        c.run(30)
        tmp = c.output.split()[0].split("|")
        #print(tmp)
        start_time = datetime.datetime.strptime(tmp[1], "%Y-%m-%dT%H:%M:%S").strftime("%Y/%m/%d %H:%M:%S")
        end_time = datetime.datetime.strptime(tmp[2], "%Y-%m-%dT%H:%M:%S").strftime("%Y/%m/%d %H:%M:%S")
        ###obtain time for lsf system
        #tmp = os.popen('bhist -l %d | grep Starting | awk \'{print $5,$2,$3,$4}\''%slurm_job_id).readlines()[0]
        #t=tmp.split()[0].split(':')[0]+'-'+tmp.split()[1]+'-'+tmp.split()[2]+'-'+tmp.split()[3]
        #start_time=datetime.datetime.strptime(t, "%Y-%b-%d-%H:%M:%S").strftime("%Y/%m/%d %H:%M:%S")
        ###obtain lsf job start time
        #tmp = os.popen('bhist -l %d | grep Done | awk \'{print $5,$2,$3,$4}\''%slurm_job_id).readlines()[0]
        #t=tmp.split()[0].split(':')[0]+'-'+tmp.split()[1]+'-'+tmp.split()[2]+'-'+tmp.split()[3]
        #end_time=datetime.datetime.strptime(t, "%Y-%b-%d-%H:%M:%S").strftime("%Y/%m/%d %H:%M:%S")
        #time.sleep(5)
        status = 0
        new = {'$set': {'run_info.start_time': start_time, 'run_info.end_time': end_time, 'run_info.status': status}}
        col.update_one(query, new)

def getRunningTasks():
    client = pymongo.MongoClient("mongodb://%s:27018"%hostname)
    db = client[client_name]
    col = db['vasptasks']
    query = {'run_info.status': 2}
    data = col.find(query, {'job_id': 1})
    running_ids = []
    for i in data:
        running_ids.append(i['job_id'])
    return running_ids

def getWaitingTasks():
    client = pymongo.MongoClient("mongodb://%s:27018"%hostname)
    db = client[client_name]
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
                    submitTask(waiting_ids[i])
            else:
                for i in range(n):
                    submitTask(waiting_ids[i])
        count = count + 1
        print("count: " + str(count) + "   " + time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(time.time())))
        time.sleep(timeinterval)
if __name__ == "__main__":
    Run(tasklimit=45, timeinterval=60)
