# -*- coding: utf-8 -*-
"""
Created on Wed Feb 23 09:27:57 2022

@author: DengJun

Note: the unitcell of hexgonal is required, 
and the angle between two lattice should be 60 degree

"""


from pymatgen.core.structure import Structure
import numpy as np
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
import pandas as pd
import os

def generate_twist_structures(path,filename):
    
    formula=filename.split('_')[2].split('.')[0]
    os.mkdir(formula)
    
    s=Structure.from_file('%s/%s'%(path,filename))
    data1=[]
    for n in range(10):
        m=10
        while m>-n:
            t=np.array([[m,n,0],[-n,m+n,0],[0,0,1]])
            s_tmp=s.copy()
            s_tmp.make_supercell(t)
            cosd=(m+n/2)/np.sqrt(m**2+n**2+m*n)
            degree=np.arccos(cosd)/np.pi*180
            data1.append({'lattice':s_tmp.lattice.a,'degree':degree,'m':m,'n':n,'atoms':len(s_tmp.atomic_numbers)})
            m=m-1
    
    data2=[]
    for n in range(10):
        m=10
        while m>-n:
            t=np.array([[m,n,0],[-n,m+n,0],[0,0,1]])
            s_tmp=s.copy()
            s_tmp.make_supercell(t)
            cosd=(m+n/2)/np.sqrt(m**2+n**2+m*n)
            degree=np.arccos(cosd)/np.pi*180
            data2.append({'lattice':s_tmp.lattice.a,'degree':degree,'m':m,'n':n,'atoms':len(s_tmp.atomic_numbers)})
            m=m-1
    
    combine=[]
    for i in data1:
        for j in data2:
            if abs((2*i['lattice'])/(i['lattice']+j['lattice'])-1) <= 0.005:
                twist=abs(i['degree']-j['degree'])
                misfit=abs((2*i['lattice'])/(i['lattice']+j['lattice'])-1)
                atoms=i['atoms']+j['atoms']
                if twist > 60:
                    twist=abs(twist-120)
                combine.append({'twist': twist,'sub':i,'ads':j,'misfit':misfit,'atoms':atoms})
    
    
    
    d=np.zeros([len(combine),len(combine)])
    for num1,i in enumerate(combine):
        for num2,j in enumerate(combine):
            if num1 == num2:
                continue
            else:
                if abs(i['twist']-j['twist']) <= 10**-5:
                    d[num1,num2]=1
    
    def list_combine(l):
        if  len(l) == 1:
            return list(set(l[0]))
        else:
            t=[]
            for i in l:
                t=t+i
            return list(set(t))
    
    def find_match(d,match):
        tmp=[]
        for i in match[-1]:
            link=[]
            link.append(i)
            for j in range(len(d)):
                if d[i][j] == 1:
                    link.append(j)
            tmp=tmp+link
        new=list(set(tmp).difference(set(list_combine(match))))
        if new == []:
            return
        match.append(new)
        find_match(d,match)
        
    matches=[]
    for i in range(len(combine)):
        match=[]
        if (i in list_combine(matches)) == False:
            match.append([i])
            find_match(d,match)
            matches.append(list_combine(match))
    keep=[]
    for i in matches:
        keep.append(i[0])
    remove=list(set(np.arange(len(combine))).difference(set(keep)))
    remove.sort(reverse=True)
    for i in remove:
        combine.pop(i)
        

    for i in combine:
        m_sub=i['sub']['m']
        n_sub=i['sub']['n']
        t_sub=np.array([[m_sub,n_sub,0],[-n_sub,m_sub+n_sub,0],[0,0,1]])
        s_sub=s.copy()
        s_sub.make_supercell(t_sub)
        
        m_ads=i['ads']['m']
        n_ads=i['ads']['n']
        t_ads=np.array([[m_ads,n_ads,0],[-n_ads,m_ads+n_ads,0],[0,0,1]])
        s_ads=s.copy()
        s_ads.make_supercell(t_ads)
        
        #######create the twisted structure
        matrix=s_sub.lattice.matrix
        d_layer=5.8 #2.9 angstrom
        
        s_coord_sub=s_sub.frac_coords  #bottom layer
        s_coord_ads=s_ads.frac_coords  #upper layer
        thickness=(np.max(s_coord_sub[:,2])-np.min(s_coord_sub[:,2]))*matrix[2,2]
        d_layer=3+thickness
        s_coord_ads[:,2]=s_coord_ads[:,2]-d_layer/matrix[2,2]
        
        
        s_symbol_sub=s_sub.atomic_numbers
        s_symbol_ads=s_ads.atomic_numbers
        
        twist_coord=np.vstack([s_coord_sub,s_coord_ads])
        twist_symbol=s_symbol_sub+s_symbol_ads
        
        ####shift the layer to center
        tmin=np.min(twist_coord[:,2])
        tmax=np.max(twist_coord[:,2])
        t_avg=(tmin+tmax)/2
        twist_coord[:,2]=twist_coord[:,2]+0.5-t_avg
        ##########
        degree=i['twist']
        # print('%.3f'%degree)
        degrees = ['0.000','13.174','21.787','32.204','38.213','42.103','46.826','60.000']
        if str('%.3f'%degree) in degrees:
            twist=Structure(matrix,twist_symbol,twist_coord)
            twist=SpacegroupAnalyzer(twist).get_conventional_standard_structure()
            matrix=np.array([[1,0,0],[1,1,0],[0,0,1]])
            twist.make_supercell(matrix)
            i['atoms']=len(twist.atomic_numbers)
            twist.to(fmt='poscar',filename=r'%s\%.3f.vasp'%(formula,degree))
    
    df=pd.DataFrame(combine)
    df.to_excel('%s.xlsx'%formula)
        
if __name__ == "__main__":
    path='./structures_opt'
    filenames=os.listdir('./structures_opt')
    filenames=['24_645529_NbTe2.vasp']
    for filename in filenames:
        generate_twist_structures(path,filename)
        print(filename)
            
        
        
