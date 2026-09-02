#!/usr/bin/env python3                                                                                                                                                                                              
# -*- coding: utf-8 -*-
"""
Created on Sun Feb 22 09:35:13 2026
 
@author: dj
"""
 
f=open('q2r.out')
lines=f.readlines()
f.close()
 
weights=[]
qpoints=[]
for num,line in enumerate(lines):
    tmp=line.split()
    if 'nqs=' in tmp:
        weights.append(int(tmp[1]))
        q=lines[num+1].split()
        qpoints.append([float(q[1]),float(q[2]),float(q[3])])
 
filename='lambda.in'
f=open(filename,'w')
emax=20
smearing=0.2
f.write('%f %f 1 \n'%(emax,smearing))
f.write(str(len(weights))+'\n')
for i in range(len(weights)):
    f.write(str(qpoints[i][0])+' '+str(qpoints[i][1])+' '+str(qpoints[i][2])+' '+str(weights[i])+'\n')
 
for i in range(len(weights)):
    f.write('elph_dir/elph.inp_lambda.%d \n'%(i+1))
 
 
f.write('0.1')    
 
f.close()
