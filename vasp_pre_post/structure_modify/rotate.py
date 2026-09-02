# -*- coding: utf-8 -*-
"""
Created on Thu Feb  1 15:16:58 2024

@author: DengJun
"""


from pymatgen.core.structure import Structure
import numpy as np
import pandas as pd
import pymongo

def rot_hex(m,n,s,shift):
    coordinates=s.frac_coords
    coordinates=coordinates+shift
    matrix=s.lattice.matrix
    symbols=s.atomic_numbers
    s=Structure(matrix,symbols,coordinates)
    matrix=np.array([[m,n,0],[-n,m-n,0],[0,0,1]])
    s.make_supercell(matrix)
    # s.to(fmt='poscar',filename='supercell2.vasp')
    return s
    

def rot_tet(m,n,s,shift):
    coordinates=s.frac_coords
    coordinates=coordinates+shift
    matrix=s.lattice.matrix
    symbols=s.atomic_numbers
    s=Structure(matrix,symbols,coordinates)
    matrix=np.array([[m,n,0],[-n,m,0],[0,0,1]])
    s.make_supercell(matrix)
    # s.to(fmt='poscar',filename='supercell2.vasp')
    return s

def make_supercell(ads,sub):
    delta=2.8 ###distaaance between layers
    s_ads=ads.copy()
    s_sub=sub.copy()
    
    ####get new a and matrix
    a_ads=s_ads.lattice.a
    a_sub=s_sub.lattice.a
    strain_ads=(a_ads-a_sub)/a_sub
    strain_sub=(a_sub-a_sub)/a_sub
    s_ads.apply_strain([strain_ads,strain_ads,0])
    s_sub.apply_strain([strain_sub,strain_sub,0])
    
    matrix=s_sub.lattice.matrix
    matrix.flags['WRITEABLE']=True
    matrix_sub=s_sub.lattice.matrix
    matrix_sub.flags['WRITEABLE']=True
    ####get atoms
    atoms=s_ads.atomic_numbers+s_sub.atomic_numbers
    
    ###get atomic positions
    sub_cart=s_sub.cart_coords
    ads_cart=s_ads.cart_coords
    
    matrix_sub=s_sub.lattice.matrix
    matrix_ads=s_ads.lattice.matrix
    matrix_t=np.concatenate((matrix_sub[[0,1],:],matrix_ads[2,:].reshape(1,3)))
    ads_cart=np.dot(np.dot(ads_cart,np.linalg.inv(matrix_ads)),matrix_t)
    
    sub_cart[:,2]=sub_cart[:,2]-np.min(sub_cart[:,2]) ##move atom at bottom
    ads_cart[:,2]=ads_cart[:,2]-np.min(ads_cart[:,2]) ##move atom at bottom
    
    sub_height=np.max(sub_cart[:,2])-np.min(sub_cart[:,2])
    ads_cart[:,2]=ads_cart[:,2]+sub_height+delta
    
    positions=np.vstack([ads_cart,sub_cart])
    
    #### new matrix
    matrix[2,2]=np.max(positions[:,2])+delta+15
    length=np.max(positions[:,2])
    positions[:,2]=positions[:,2]+(matrix[2,2]-length)/2
    ads_atoms=np.arange(1,len(ads_cart)+1)
    sub_atoms=np.arange(1,len(sub_cart)+1)+len(ads_cart)
    ####generate new structure
    s=Structure(matrix,atoms,positions,coords_are_cartesian=True)
    # s=SpacegroupAnalyzer(s).get_conventional_standard_structure()
    # s.to(fmt='poscar',filename=path+'/'+str(num)+'_'+s.formula.replace(' ','')+'.vasp')
    return s,ads_atoms,sub_atoms    
    
def find_pair(sub,ads,mismatch):
    shift=[0,0,0]
    data_sub=[]
    for n in range(10):
        for m in range(10):
            if not (m == 0 and n == 0):
                sub_tmp=sub.copy()
                sub_tmp=rot_hex(m,n,sub_tmp,shift)
                cosd=(m+n/2)/np.sqrt(m**2+n**2+m*n)
                degree=np.arccos(cosd)/np.pi*180
                data_sub.append({'lattice':sub_tmp.lattice.a,'degree':degree,'m':m,'n':n,'atoms':len(sub_tmp.atomic_numbers)})
    
    data_ads=[]
    for n in range(10):
        for m in range(10):
            if not (m == 0 and n == 0):
                ads_tmp=ads.copy()
                ads_tmp=rot_hex(m,n,ads_tmp,shift)
                cosd=(m+n/2)/np.sqrt(m**2+n**2+m*n)
                degree=np.arccos(cosd)/np.pi*180
                data_ads.append({'lattice':ads_tmp.lattice.a,'degree':degree,'m':m,'n':n,'atoms':len(ads_tmp.atomic_numbers)})
                
    combine=[]
    for i in data_sub:
        for j in data_ads:
            if abs((2*i['lattice'])/(i['lattice']+j['lattice'])-1) <= mismatch:
                twist=abs(i['degree']-j['degree'])
                misfit=abs((2*i['lattice'])/(i['lattice']+j['lattice'])-1)
                atoms=i['atoms']+j['atoms']
                combine.append({'twist': twist,'sub':i,'ads':j,'misfit':misfit,'atoms':atoms,'m_sub':i['m'],'n_sub':i['n'],'m_ads':j['m'],'n_ads':j['n']})
    return combine

def create_interfaces(ads_id,sub_id,mismatch,shift):
    client=pymongo.MongoClient("mongodb://172.32.12.2:27017")
    db=client['2DMSI']
    ads_query={'s_id':ads_id}
    sub_query={'s_id':sub_id}
    data_ads=db['semiconductors'].find_one(ads_query)
    data_sub=db['metals'].find_one(sub_query)
    ads=Structure.from_dict(data_ads['structure'])
    sub=Structure.from_dict(data_sub['surface']['structure'])
    # fname=data_ads['formula']+'-'+data_sub['formula']
    combinations=find_pair(sub,ads,mismatch)
    df=pd.DataFrame(combinations)
    row=df['atoms'].idxmin()
    data=df.iloc[row].to_dict()
    # df.to_excel('combine.xlsx')
    # shift=[0,0,0]
    s_ads_rot= rot_hex(data['m_ads'],data['n_ads'],ads,shift)  
    s_sub_rot= rot_hex(data['m_sub'],data['n_sub'],sub,[0,0,0])   
    [interface,ads_atoms,sub_atoms]=make_supercell(s_ads_rot,s_sub_rot)
    # interface.to(filename='structures/%s.vasp'%fname,fmt='poscar')
    return interface,data['atoms'],data['m_ads'],data['n_ads'],data['m_sub'],\
        data['n_sub'],data['misfit'],data['twist'],data_ads['formula'],data_sub['formula'],\
            ads_atoms,sub_atoms

def selective_dynamics(interface):
    
    coordinates=interface.frac_coords
    c=coordinates[:,2]
    c=list(set(c))
    c.sort()
    if interface.atomic_numbers[-1] == 42:
        fix=c[5]
    else:
        fix=c[2]
    selective=[]
    for i in coordinates:
        if i[2] <= fix:
            selective.append([False,False,False])
        else:
            selective.append([True,True,True])
    interface.add_site_property("selective_dynamics",selective)
    return interface
    
if __name__ == '__main__':
    client=pymongo.MongoClient("mongodb://172.32.12.2:27017")
    db=client['2DMSI']
    mismatch=0.05
    s_ids_metal=[]
    s_ids_semiconductor=[]
    query={'$or':[{'crys_sys':'trigonal'},{'crys_sys':'hexagonal'}]}
    entries=db['semiconductors'].find(query)
    for entry in entries:
        s_ids_semiconductor.append(entry['s_id'])
    entries=db['metals'].find()
    for entry in entries:
        s_ids_metal.append(entry['s_id'])
    
    interfaces=[]
    shifts=np.array([[0,0,0],[2/3,1/3,0],[1/3,2/3,0],[0.5,0,0]])
    col=db['interfaces']
    id_nums=[]
    for num1,sub_id in enumerate(s_ids_metal):
        for num2,ads_id in enumerate(s_ids_semiconductor):
            configs=[]
            for num3,shift in enumerate(shifts):
                [interface,natoms,m_ads,n_ads,m_sub,n_sub,misfit,twist,ads_formula,sub_formula,ads_atoms,sub_atoms]\
                    =create_interfaces(ads_id,sub_id,mismatch,shift)
                interface=selective_dynamics(interface)
                configs.append({'structure':interface.as_dict(),'shift':shift.tolist(),'config_id':num3})
            id_num=num1*len(s_ids_semiconductor)+num2
            if_id='if_%d'%id_num
            configs_tmp={'0':configs[0],'1':configs[1],'2':configs[2],'3':configs[3]}
            rootpath='/edahome/y00201/Documents/vasp/2D_materials/HT_calculations/interfaces/works/'+if_id
            data={'if_id':if_id,'natoms':natoms,'rootpath':rootpath,\
                  'twist_angle':twist,'misfit':misfit,'configs':configs_tmp,\
                  'sub':{'formula':sub_formula,'m':m_sub,'n':n_sub,'s_id':sub_id,\
                         'atoms':sub_atoms.tolist()},\
                  'ads':{'formula':ads_formula,'m':m_ads,'n':n_ads,'s_id':ads_id,\
                         'atoms':ads_atoms.tolist()}}
            col.insert_one(data)    
            # interfaces.append([len(interface.atomic_numbers),interface])
            print(sub_id,ads_id)
            
    
    
    
            

