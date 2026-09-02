#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 22 17:42:55 2022

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
from pymatgen.io.vasp import Vasprun
from pymatgen.analysis.phase_diagram import PhaseDiagram, PDPlotter
from pymatgen.entries.compatibility import MaterialsProjectCompatibility
import itertools
import numpy as np
import datetime
import time
import operator
from phonopy.api_phonopy import Phonopy
# from phonopy.structure.atoms import PhonopyAtoms
from pymatgen.core.structure import Structure
from phonopy.interface.vasp import read_vasp,write_vasp,get_vasp_structure_lines,parse_set_of_forces,Vasprun
from pymatgen.io.vasp.sets import MPNonSCFSet
from phonopy.cui.create_force_sets import create_FORCE_SETS
from phonopy.phonon.band_structure import get_band_qpoints_and_path_connections
from pymatgen.symmetry.bandstructure import HighSymmKpath
from pymatgen.io.vasp.inputs import Kpoints
import matplotlib.pyplot as plt
import matplotlib as mpl
from phonopy.structure.atoms import PhonopyAtoms

def postprocess_opt():
    client=pymongo.MongoClient("mongodb://localhost:27017")
    db=client['layers']
    col=db['vasptasks']
    query={'$and': [{'run_info.status': 0},{'post_info.ionic_converge': True},{'task_info.task_type': 'relax'}]}
    data=col.find(query)
    structure_ids=[]
    for i in data:
        query={'structure_id':i['task_info']['structure_id']}
        path=i['run_info']['path']
        # path=str.join('/',path.split('/')[:-1])
        new={'$set':{'structure_opt':i['post_info']['contcar'],'path':path,'mag_opt':i['post_info']['tot_mag']}}
        db['stacking_ionic'].update_one(query,new)
        print(i['task_info']['structure_id'])
        # Structure.from_dict(i['post_info']['contcar']).to(fmt='poscar',filename='structures_reopt2/%s.vasp'%i['task_info']['structure_id'])

def postprocess_scf():
    client = pymongo.MongoClient("mongodb://localhost:27017")
    db = client['layers']
    col = db['vasptasks']
    # query = {'post_info.ele_converge': True, 'task_info.task_type': 'scf'}
    query = {'$and': [{'run_info.status': 0}, {'post_info.ele_converge': True}, {'task_info.task_type': 'scf'}]}
    data = col.find(query)
    structure_ids = []
    for i in data:
        query = {'structure_id': i['task_info']['structure_id']}
        print(i['job_id'], i['task_info']['structure_id'])
        new = {'$set': {'mag_scf': i['post_info']['tot_mag'], 'efermi': i['post_info']['efermi'],'etot': i['post_info']['etot'],'mag_ele':i['post_info']['mag_ele']}}
        # new={'$set':{'mag_ele':i['post_info']['mag_ele']}}
        db['stacking_ionic'].update_one(query, new)


def get_ehull(path):
    vasprun = Vasprun(path + "/vasprun.xml")
    s = vasprun.final_structure.as_dict()
    labels = []
    for i in s['sites']:
        labels.append(i['label'])
    labels = list(set(labels))
    entry = vasprun.get_computed_entry(inc_structure=True)

    a = MPRester('ezvWtmTqlDJNfCVF')
    mp_entries = a.get_entries_in_chemsys(labels)
    compatibility = MaterialsProjectCompatibility()
    entry = compatibility.process_entry(entry)
    entries = compatibility.process_entries([entry] + mp_entries)
    pd = PhaseDiagram(entries)
    ehull = pd.get_e_above_hull(entry)
    return ehull

def pross_ehull():
    client = pymongo.MongoClient("mongodb://localhost:27017")
    db = client['layers']
    col = db['stacking_ionic']

    # datas = col.find({'band': {'$exists': True}})
    query={}
    datas=col.find(query)
    ids = []
    ehulls = []
    for num, data in enumerate(datas):
        try:
            query = {'structure_id': data['structure_id']}
            data = col.find_one(query)
            path = data['path'] + '/scf'
            # data['path']=data['path'].replace('/home/jdeng','/home/dj/C2')
            ehull = get_ehull(path)
            ehulls.append([data['structure_id'], ehull])
            new = {'$set': {'e_above_hull': ehull}}
            col.update_one(query, new)
            print(num, data['structure_id'], ehull)
        except:
            ehulls.append([data['structure_id'], 'error'])
            print(num, data['structure_id'], 'error')

def read_EIGENVAL(path):
    f = open(path + 'POSCAR', 'r')
    lines = f.readlines()
    f.close()
    elements = lines[5].split()
    nelements = [int(i) for i in lines[6].split()]
    a = [x for x in lines[2].split()]
    b = [x for x in lines[3].split()]
    c = [x for x in lines[4].split()]
    a11 = float(a[0])
    a12 = float(a[1])
    a13 = float(a[2])
    a21 = float(b[0])
    a22 = float(b[1])
    a23 = float(b[2])
    a31 = float(c[0])
    a32 = float(c[1])
    a33 = float(c[2])
    V = a11 * (a22 * a33 - a23 * a32) - a12 * (a21 * a33 - a23 * a31) + a13 * (a21 * a32 - a22 * a31)

    b11 = (a22 * a33 - a23 * a32) / V
    b12 = -(a21 * a33 - a23 * a31) / V
    b13 = (a21 * a32 - a31 * a22) / V

    b21 = (a32 * a13 - a33 * a12) / V
    b22 = -(a31 * a13 - a33 * a11) / V
    b23 = (a31 * a12 - a32 * a11) / V

    b31 = (a12 * a23 - a13 * a22) / V
    b32 = (a11 * a23 - a13 * a21) / V
    b33 = (a11 * a22 - a12 * a21) / V

    recipro_lattice = np.array([[b11, b12, b13], [b21, b22, b23], [b31, b32, b33]])

    efermi = float(
        os.popen('grep efermi %s../scf/vasprun.xml |awk \'{print $3}\'' % path).readlines()[0].split()[0])

    f = open(path + 'EIGENVAL', 'r')
    lines = f.readlines()
    f.close

    data = {'K': [], 'E_up': []}

    for line in lines[7:]:
        values = [float(x) for x in line[:].split()]
        if len(values) == 4:
            data['K'].append(values[0:3])
            E_up = []
        elif len(values) == 3:
            E_up.append(values[1])
        elif len(values) == 0:
            data['E_up'].append(E_up)
    data['E_up'].append(E_up)
    data['E_up'] = np.array(data['E_up'])
    data['K'] = np.array(data['K'])
    data['K'] = np.matmul(data['K'], recipro_lattice)

    f = open(path + 'KPOINTS', 'r')
    lines = f.readlines()
    f.close()
    k_path = []
    high_sym = []
    ticks = []
    num_high_sym = []
    for num, line in enumerate(lines[3:]):
        tmp = [i for i in line.split()]
        k_path.append(np.array([float(tmp[0]), float(tmp[1]), float(tmp[2])]))
        if len(tmp) == 5:
            ticks.append(tmp[4])
            high_sym.append(np.array([float(tmp[0]), float(tmp[1]), float(tmp[2])]))
            num_high_sym.append(num)

    band_data = []
    for i in np.arange(0, len(num_high_sym), 2):
        path = {'kpoints': data['K'][num_high_sym[i]:num_high_sym[i + 1] + 1, :],
                'E_up': data['E_up'][num_high_sym[i]:num_high_sym[i + 1] + 1, :]}
        band_data.append(path)

    for i in band_data:
        x = np.zeros(len(i['kpoints']))
        for j in range(len(i['kpoints']) - 1):
            x[j + 1] = x[j] + np.linalg.norm(i['kpoints'][j + 1] - i['kpoints'][j])
        x = x.reshape(len(x), 1)
        i['x'] = x
        i_band_up = np.concatenate((x, i['E_up']), axis=1)
        i['band_up'] = i_band_up

    band_up = []
    for num, i in enumerate(band_data):
        if num == 0:
            band_up = i['band_up']
        else:
            i['band_up'][:, 0] = i['band_up'][:, 0] + band_data[num - 1]['band_up'][-1, 0]
            band_up = np.concatenate((band_up, i['band_up']), axis=0)
    # np.savetxt('band.dat',band_up,'%9.06f')

    high_sym = []
    for i in num_high_sym[::2]:
        high_sym.append(band_up[i, 0])
    high_sym.append(band_up[-1, 0])
    # np.savetxt('high_sym.dat',high_sym,'%9.06f')

    tick_tmp = [ticks[0]]
    for i in np.arange(1, len(ticks) - 1, 2):
        if ticks[i] == ticks[i + 1]:
            tick_tmp.append(ticks[i])
        else:
            tick_tmp.append(ticks[i] + '|' + ticks[i + 1])
    tick_tmp.append(ticks[-1])

    ticks = tick_tmp
    # f=open('ticks.dat','w')
    # for i in ticks:
    #     f.write(i+'\n')
    # f.close()
    # efermi=np.array(float(os.popen('grep efermi %s../scf/vasprun.xml |awk \'{print $3}\''%path).readlines()[0])).reshape(1,1)
    data = {'eigval': band_up, 'ticks': ticks, 'high_sym': high_sym, 'elements': elements, 'nelements': nelements,
            'efermi': efermi}
    return data


def read_PROCAR(path):
    f = open(path + 'PROCAR')
    lines = f.readlines()
    f.close()
    nkpoints = int(lines[1].split()[3])
    nbands = int(lines[1].split()[7])
    nions = int(lines[1].split()[11])

    tmp = []
    for num, line in enumerate(lines):
        t = line.split()
        if t == []:
            continue
        elif t[0] == 'band':
            tmp.append(num)
        if len(tmp) == 2:
            break
    bcycle = tmp[1] - tmp[0]

    tmp = []
    for num, line in enumerate(lines):
        t = line.split()
        if t == []:
            continue
        elif t[0] == 'k-point':
            tmp.append(num)
        if len(tmp) == 2:
            break
    kcycle = tmp[1] - tmp[0]

    kpoints = []
    lines = lines[3:]
    mtot = []
    energy = []
    pro = []
    for i in range(nkpoints):
        # kpoints.append([float(i) for i in lines[i*kcycle].split()[3:6]])
        kpoints.append(
            [float(lines[i * kcycle][19:30]), float(lines[i * kcycle][30:41]), float(lines[i * kcycle][41:52])])
        klines = lines[2 + i * kcycle:(i + 1) * kcycle]
        for j in range(nbands):
            t = klines[(j * bcycle):(j + 1) * bcycle]
            energy.append(float(t[0].split()[4]))
            blines = klines[(3 + j * bcycle):(j + 1) * bcycle]
            pro.append(blines[0:nions])

    data = {'K': [], 'E': []}
    pband = {}
    for i in range(nions):
        pband['ion_' + str(i)] = []
        # print('ions: '+str(i))
        for j in range(nbands):
            for k in range(nkpoints):
                line = k * nbands + j
                pro_i = list(map(float, pro[line][i].split()))
                pband['ion_' + str(i)].append(pro_i)

    for i in range(nbands):
        for j in range(nkpoints):
            data['E'].append(energy[j * nbands + i])
    data['E'] = np.array(data['E'])
    for i in range(nions):
        pband['ion_' + str(i)] = np.array(pband['ion_' + str(i)])
        pband['ion_' + str(i)] = np.delete(pband['ion_' + str(i)], 0, axis=1)

    return pband


def read_EIGENVAL_spin(path):
    f = open(path + 'POSCAR', 'r')
    lines = f.readlines()
    f.close()
    elements = lines[5].split()
    nelements = [int(i) for i in lines[6].split()]
    a = [x for x in lines[2].split()]
    b = [x for x in lines[3].split()]
    c = [x for x in lines[4].split()]
    a11 = float(a[0])
    a12 = float(a[1])
    a13 = float(a[2])
    a21 = float(b[0])
    a22 = float(b[1])
    a23 = float(b[2])
    a31 = float(c[0])
    a32 = float(c[1])
    a33 = float(c[2])
    V = a11 * (a22 * a33 - a23 * a32) - a12 * (a21 * a33 - a23 * a31) + a13 * (a21 * a32 - a22 * a31)

    b11 = (a22 * a33 - a23 * a32) / V
    b12 = -(a21 * a33 - a23 * a31) / V
    b13 = (a21 * a32 - a31 * a22) / V

    b21 = (a32 * a13 - a33 * a12) / V
    b22 = -(a31 * a13 - a33 * a11) / V
    b23 = (a31 * a12 - a32 * a11) / V

    b31 = (a12 * a23 - a13 * a22) / V
    b32 = (a11 * a23 - a13 * a21) / V
    b33 = (a11 * a22 - a12 * a21) / V

    recipro_lattice = np.array([[b11, b12, b13], [b21, b22, b23], [b31, b32, b33]])

    efermi = float(
        os.popen('grep efermi %s../scf/vasprun.xml |awk \'{print $3}\'' % path).readlines()[0].split()[0])

    f = open(path + 'EIGENVAL', 'r')
    lines = f.readlines()
    f.close

    data = {'K': [], 'E_up': [], 'E_dw': []}

    for line in lines[7:]:
        values = [float(x) for x in line[:].split()]
        if len(values) == 4:
            data['K'].append(values[0:3])
            E_up = []
            E_dw = []
        elif len(values) == 5:
            E_up.append(values[1])
            E_dw.append(values[2])
        elif len(values) == 0:
            data['E_up'].append(E_up)
            data['E_dw'].append(E_dw)
    data['E_up'].append(E_up)
    data['E_dw'].append(E_dw)
    data['E_up'] = np.array(data['E_up'])
    data['E_dw'] = np.array(data['E_dw'])
    data['K'] = np.array(data['K'])
    data['K'] = np.matmul(data['K'], recipro_lattice)

    f = open(path + 'KPOINTS', 'r')
    lines = f.readlines()
    f.close()
    k_path = []
    high_sym = []
    ticks = []
    num_high_sym = []
    for num, line in enumerate(lines[3:]):
        tmp = [i for i in line.split()]
        k_path.append(np.array([float(tmp[0]), float(tmp[1]), float(tmp[2])]))
        if len(tmp) == 5:
            ticks.append(tmp[4])
            high_sym.append(np.array([float(tmp[0]), float(tmp[1]), float(tmp[2])]))
            num_high_sym.append(num)

    band_data = []
    for i in np.arange(0, len(num_high_sym), 2):
        path = {'kpoints': data['K'][num_high_sym[i]:num_high_sym[i + 1] + 1, :],
                'E_up': data['E_up'][num_high_sym[i]:num_high_sym[i + 1] + 1, :],
                'E_dw': data['E_dw'][num_high_sym[i]:num_high_sym[i + 1] + 1, :]}
        band_data.append(path)

    for i in band_data:
        x = np.zeros(len(i['kpoints']))
        for j in range(len(i['kpoints']) - 1):
            x[j + 1] = x[j] + np.linalg.norm(i['kpoints'][j + 1] - i['kpoints'][j])
        x = x.reshape(len(x), 1)
        i['x'] = x
        i_band_up = np.concatenate((x, i['E_up']), axis=1)
        i_band_dw = np.concatenate((x, i['E_dw']), axis=1)
        i['band_up'] = i_band_up
        i['band_dw'] = i_band_dw

    band_up = []
    band_dw = []
    for num, i in enumerate(band_data):
        if num == 0:
            # band_up.append(i['band_up'])
            # band_dw.append(i['band_dw'])
            band_up = i['band_up']
            band_dw = i['band_dw']
        else:
            i['band_up'][:, 0] = i['band_up'][:, 0] + band_data[num - 1]['band_up'][-1, 0]
            i['band_dw'][:, 0] = i['band_dw'][:, 0] + band_data[num - 1]['band_dw'][-1, 0]
            # band_up.append(i['band_up'])
            # band_dw.append(i['band_dw'])
            band_up = np.concatenate((band_up, i['band_up']), axis=0)
            band_dw = np.concatenate((band_dw, i['band_dw']), axis=0)
    # band_up=np.concatenate((band_up[:,0].reshape(len(band_up[:,0]),1),data['E_up']),axis=1)
    # band_dw=np.concatenate((band_dw[:,0].reshape(len(band_dw[:,0]),1),data['E_dw']),axis=1)
    # np.savetxt('band_up.dat',band_up,'%9.06f')
    # np.savetxt('band_dw.dat',band_dw,'%9.06f')

    high_sym = []
    for i in num_high_sym[::2]:
        high_sym.append(band_up[i, 0])
    high_sym.append(band_up[-1, 0])
    # np.savetxt('high_sym.dat',high_sym,'%9.06f')

    tick_tmp = [ticks[0]]
    for i in np.arange(1, len(ticks) - 1, 2):
        if ticks[i] == ticks[i + 1]:
            tick_tmp.append(ticks[i])
        else:
            tick_tmp.append(ticks[i] + '|' + ticks[i + 1])
    tick_tmp.append(ticks[-1])

    ticks = tick_tmp
    # f=open('ticks.dat','w')
    # for i in ticks:
    #     f.write(i+'\n')
    # f.close()
    # efermi=np.array(float(os.popen('grep efermi %s../scf/vasprun.xml |awk \'{print $3}\''%path).readlines()[0])).reshape(1,1)
    data = {'eigval_up': band_up, 'eigval_dw': band_dw, 'ticks': ticks, 'high_sym': high_sym, 'elements': elements,
            'nelements': nelements, 'efermi': efermi}
    return data


def read_PROCAR_spin(path):
    f = open(path + 'PROCAR')
    lines = f.readlines()
    f.close()
    nkpoints = int(lines[1].split()[3])
    nbands = int(lines[1].split()[7])
    nions = int(lines[1].split()[11])

    tmp = []
    for num, line in enumerate(lines):
        t = line.split()
        if t == []:
            continue
        elif t[0] == 'band':
            tmp.append(num)
        if len(tmp) == 2:
            break
    bcycle = tmp[1] - tmp[0]

    tmp = []
    for num, line in enumerate(lines):
        t = line.split()
        if t == []:
            continue
        elif t[0] == 'k-point':
            tmp.append(num)
        if len(tmp) == 2:
            break
    kcycle = tmp[1] - tmp[0]

    # bcycle=(nions+1)+4
    # kcycle=bcycle*nbands+3

    kpoints = []
    ####spin up
    lines = lines[3:]
    energy_up = []
    pro_up = []
    for i in range(nkpoints):
        kpoints.append(
            [float(lines[i * kcycle][19:30]), float(lines[i * kcycle][30:41]), float(lines[i * kcycle][41:52])])
        klines = lines[2 + i * kcycle:(i + 1) * kcycle]
        # print('kpoints (up) '+str(i+1))
        for j in range(nbands):
            t = klines[(j * bcycle):(j + 1) * bcycle]
            energy_up.append(float(t[0].split()[4]))
            blines = klines[(3 + j * bcycle):(j + 1) * bcycle]
            pro_up.append(blines[0:nions])

    lines = lines[kcycle * nkpoints + 1:]
    energy_dw = []
    pro_dw = []
    for i in range(nkpoints):
        # kpoints.append([float(i) for i in lines[i*kcycle].split()[3:6]])
        klines = lines[2 + i * kcycle:(i + 1) * kcycle]
        # print('kpoints (dw)'+str(i+1))
        for j in range(nbands):
            t = klines[(j * bcycle):(j + 1) * bcycle]
            energy_dw.append(float(t[0].split()[4]))
            blines = klines[(3 + j * bcycle):(j + 1) * bcycle]
            pro_dw.append(blines[0:nions])

    data = {'K': [], 'E_up': [], 'E_dw': []}
    pband = {}
    for i in range(nions):
        pband['ion_' + str(i)] = {}
        pband['ion_' + str(i)]['up'] = []
        pband['ion_' + str(i)]['dw'] = []
        # print('ions: '+str(i))
        for j in range(nbands):
            for k in range(nkpoints):
                line = k * nbands + j
                pro_i_up = list(map(float, pro_up[line][i].split()))
                pro_i_dw = list(map(float, pro_dw[line][i].split()))
                pband['ion_' + str(i)]['up'].append(pro_i_up)
                pband['ion_' + str(i)]['dw'].append(pro_i_dw)

    for i in range(nbands):
        for j in range(nkpoints):
            data['E_up'].append(energy_up[j * nbands + i])
            data['E_dw'].append(energy_dw[j * nbands + i])
    data['E_up'] = np.array(data['E_up'])
    data['E_dw'] = np.array(data['E_dw'])
    for i in range(nions):
        pband['ion_' + str(i)]['up'] = np.array(pband['ion_' + str(i)]['up'])
        pband['ion_' + str(i)]['dw'] = np.array(pband['ion_' + str(i)]['dw'])
        pband['ion_' + str(i)]['up'] = np.delete(pband['ion_' + str(i)]['up'], 0, axis=1)
        pband['ion_' + str(i)]['dw'] = np.delete(pband['ion_' + str(i)]['dw'], 0, axis=1)
    return pband

def postprocess_band():
    client = pymongo.MongoClient("mongodb://localhost:27017")
    db = client['layers']
    col = db['stacking_ionic']
    # query={'band' : { '$exists': True },'e_above_hull':{'$lt':0.1}}
    # query = {'flag': 'keep', "$or": [{'phonon.stability': 'high'}, {'phonon.stability': 'stable'}]}
    query={}
    datas = col.find(query)
    for num, data in enumerate(datas):
        path = data['band']['path']
        # path=path.replace('jdeng/','dj/C2/')
        bs_vasprun = Vasprun(path, parse_projected_eigen=True)
        bs_data = bs_vasprun.get_band_structure(line_mode=1)
        # plot=BSDOSPlotter(bs_projection='None',vb_energy_range=4,cb_energy_range=4,font="Arial")
        # plot.get_plot(bs=bs_data)
        structure_id = data['structure_id']
        # plt.savefig('fig_eah0.1/'+structure_id+'.png',format='png',dpi=600,bbox_inches='tight')
        bs_dict = bs_data.as_dict()
        bs_dict['vbm'].pop('projections')
        bs_dict['cbm'].pop('projections')
        bs_dict['cbm'].pop('kpoint_index')
        bs_dict['vbm'].pop('kpoint_index')
        new = {'$set': {'band.band_gap': bs_dict['band_gap'], 'band.is_metal': bs_dict['is_metal'],
                        'band.is_spin_polarized': bs_dict['is_spin_polarized'], 'band.vbm': bs_dict['vbm'],
                        'band.cbm': bs_dict['cbm']}}
        query = {'structure_id': structure_id}
        col.update_one(query, new)
        print(num, structure_id)

def postprocess_band_get_data():
    client = pymongo.MongoClient("mongodb://localhost:27017")
    db = client['layers']
    col = db['stacking_ionic']
    # query = {'flag': 'keep', "$or": [{'phonon.stability': 'high'}, {'phonon.stability': 'stable'}]}
    query = {}
    datas = col.find(query)
    for num, data in enumerate(datas):
        # structure_id='s-6'
        # query={'structure_id':structure_id}
        # data=col.find_one(query)
        path = data['path'] + '/band/'
        # path=path.replace('jdeng/','dj/C2/')
        if data['band']['is_spin_polarized'] == True:
            eigval = read_EIGENVAL_spin(path)
            # procar=read_PROCAR_spin(path)
            eigval['eigval_up'] = eigval['eigval_up'].tolist()
            eigval['eigval_dw'] = eigval['eigval_dw'].tolist()
            # for i in procar.keys():
            # procar[i]['up']=procar[i]['up'].tolist()
            # procar[i]['dw']=procar[i]['dw'].tolist()
            # new={'$set':{'band.eigval_up':eigval['eigval_up'],'band.eigval_dw':eigval['eigval_dw'],'band.elements':eigval['elements'],'band.high_sym':eigval['high_sym'],'band.nelements':eigval['nelements'],'band.ticks':eigval['ticks'],'band.procar':procar}}
            new = {'$set': {'band.eigval_up': eigval['eigval_up'], 'band.eigval_dw': eigval['eigval_dw'],
                            'band.elements': eigval['elements'], 'band.high_sym': eigval['high_sym'],
                            'band.nelements': eigval['nelements'], 'band.ticks': eigval['ticks'],
                            'band.efermi': eigval['efermi']}}
            query = {'structure_id': data['structure_id']}
            col.update_one(query, new)
        else:
            eigval = read_EIGENVAL(path)
            # procar=read_PROCAR(path)
            eigval['eigval'] = eigval['eigval'].tolist()
            # for i in procar.keys():
            #     procar[i]=procar[i].tolist()
            # new={'$set':{'band.eigval':eigval['eigval'],'band.elements':eigval['elements'],'band.high_sym':eigval['high_sym'],'band.nelements':eigval['nelements'],'band.ticks':eigval['ticks'],'band.procar':procar}}
            new = {'$set': {'band.eigval': eigval['eigval'], 'band.elements': eigval['elements'],
                            'band.high_sym': eigval['high_sym'], 'band.nelements': eigval['nelements'],
                            'band.ticks': eigval['ticks'], 'band.efermi': eigval['efermi']}}
            query = {'structure_id': data['structure_id']}
            col.update_one(query, new)
        print(num, data['structure_id'])


def get_phonon(s_id):
    client = pymongo.MongoClient("mongodb://localhost:27017")
    db = client['layers']
    col = db['stacking_ionic']
    query = {'structure_id': s_id}
    data = col.find_one(query)
    path = data['path'] + '/phonon'

    n = len(data['phonon']['structures'])
    filenames = []
    for i in range(0, n):
        filenames.append(path + '/' + str(i) + '/vasprun.xml')

    unitcell = read_vasp(path + '/opt/CONTCAR')
    lattice = unitcell.get_cell()
    a = np.linalg.norm(lattice[0])
    b = np.linalg.norm(lattice[1])
    c = np.linalg.norm(lattice[2])
    na, nb, nc = np.ceil(15 / a), np.ceil(15 / b), np.ceil(15 / c)

    sc_matrix = [[na, 0, 0], [0, nb, 0], [0, 0, nc]]
    phonon = Phonopy(unitcell, supercell_matrix=sc_matrix)
    phonon.generate_displacements(distance=0.01)
    supercells = phonon.supercells_with_displacements

    dataset = phonon.dataset

    force_sets = parse_set_of_forces(dataset['natom'], filenames)
    phonon.forces = force_sets
    dataset = phonon.dataset
    phonon.produce_force_constants()

    s = Structure.from_file(path + '/opt/CONTCAR')
    kpath = HighSymmKpath(s).kpath
    high_syms = []
    for i in kpath['path']:
        high_sym = []
        for j in i:
            high_sym.append(list(kpath['kpoints'][j]))
        high_syms.append(high_sym)
    labels = kpath['path']
    # kpath=HighSymmKpath(s).kpath
    # high_syms=[]
    # for i in kpath['path'][0]:
    #     high_syms.append(list(kpath['kpoints'][i]))
    # high_syms=[high_syms]

    qpoints, connections = get_band_qpoints_and_path_connections(high_syms, npoints=101)
    phonon.run_band_structure(paths=qpoints, with_eigenvectors=False, is_band_connection=False)
    phonon.run_mesh([30, 30, 30])
    phonon.run_total_dos()
    band = phonon.get_band_structure_dict()
    dos = phonon.get_total_dos_dict()
    phonon.plot_band_structure_and_dos().show()
    plt.savefig('./figs/' + s_id + '.png', dpi=300)

    data = {'unitcell': unitcell, 'sc_matrix': sc_matrix, 'dataset': dataset, "band": band, 'dos': dos,
            'force_sets': force_sets, 'labels': labels}
    # np.save('./data/' + s_id + '.npy', np.array(data))
    return data


def put_data_to_database(s_id):
    data = get_phonon(s_id)
    cell = data['unitcell'].get_cell().tolist()
    symbols = data['unitcell'].get_chemical_symbols()
    scaled_positions = data['unitcell'].get_scaled_positions().tolist()
    data.pop('unitcell')
    data['unitcell'] = {'cell': cell, 'symbols': symbols, 'scaled_positions': scaled_positions}

    ##############convert dosdata to list
    for i in data['dos'].keys():
        data['dos'][i] = data['dos'][i].tolist()

    ##############convert force_sets to list
    force_sets = data['force_sets']
    tmp = []
    for i in force_sets:
        tmp.append(i.tolist())
    data.pop('force_sets')
    data['force_sets'] = tmp

    ##############convert band to list
    for num, i in enumerate(data['band']['distances']):
        data['band']['distances'][num] = data['band']['distances'][num].tolist()

    for num, i in enumerate(data['band']['frequencies']):
        data['band']['frequencies'][num] = data['band']['frequencies'][num].tolist()

    for num, i in enumerate(data['band']['qpoints']):
        data['band']['qpoints'][num] = data['band']['qpoints'][num].tolist()

    ###############convert dataset to list
    for num, i in enumerate(data['dataset']['first_atoms']):
        data['dataset']['first_atoms'][num]['forces'] = data['dataset']['first_atoms'][num]['forces'].tolist()

    client = pymongo.MongoClient("mongodb://localhost:27017")
    db = client['layers']
    col = db['stacking_ionic']
    query = {'structure_id': s_id}
    new = {'$set': {'phonon.data': data}}
    col.update_one(query, new)
    # new={'$unset':{'phonon.data.force_sets': ""}}
    # col.update_one(query,new)
    return

def recover_phonon_from_db(s_id):
    client = pymongo.MongoClient("mongodb://localhost:27017")
    db = client['layers']
    col = db['stacking_ionic']
    query = {'structure_id': s_id}
    data = col.find_one(query)['phonon']['data']

    ##############convert dosdata to list
    for i in data['dos'].keys():
        data['dos'][i] = np.array(data['dos'][i])

    ##############convert force_sets to list
    force_sets = data['force_sets']
    tmp = []
    for i in force_sets:
        tmp.append(np.array(i))
    data.pop('force_sets')
    data['force_sets'] = tmp

    ##############convert band to list
    for num, i in enumerate(data['band']['distances']):
        data['band']['distances'][num] = np.array(data['band']['distances'][num])

    for num, i in enumerate(data['band']['frequencies']):
        data['band']['frequencies'][num] = np.array(data['band']['frequencies'][num])

    for num, i in enumerate(data['band']['qpoints']):
        data['band']['qpoints'][num] = np.array(data['band']['qpoints'][num])

    ###############convert dataset to list
    for num, i in enumerate(data['dataset']['first_atoms']):
        data['dataset']['first_atoms'][num]['forces'] = np.array(data['dataset']['first_atoms'][num]['forces'])

    unitcell = PhonopyAtoms(cell=data['unitcell']['cell'], symbols=data['unitcell']['symbols'],
                            scaled_positions=data['unitcell']['scaled_positions'])
    phonon = Phonopy(unitcell, supercell_matrix=data['sc_matrix'])
    phonon.dataset = data['dataset']
    phonon.forces = data['force_sets']
    phonon.produce_force_constants()
    phonon.run_band_structure(paths=data['band']['qpoints'], with_eigenvectors=False, is_band_connection=False)
    phonon.plot_band_structure().show()