#!/home/software/anaconda3/bin/python
 
import math
import sys 
import os
import numpy as np
import re

f=open('POSCAR','r')
lines=f.readlines()
f.close()
 
a=[x for x in lines[2].split()]
b=[x for x in lines[3].split()]
c=[x for x in lines[4].split()]
a11=float(a[0])
a12=float(a[1])
a13=float(a[2])
a21=float(b[0])
a22=float(b[1])
a23=float(b[2])
a31=float(c[0])
a32=float(c[1])
a33=float(c[2])
V=a11*(a22*a33-a23*a32)-a12*(a21*a33-a23*a31)+a13*(a21*a32-a22*a31)

b11=(a22*a33-a23*a32)/V
b12=-(a21*a33-a23*a31)/V
b13=(a21*a32-a31*a22)/V

b21=(a32*a13-a33*a12)/V
b22=-(a31*a13-a33*a11)/V
b23=(a31*a12-a32*a11)/V

b31=(a12*a23-a13*a22)/V
b32=-(a11*a23-a13*a21)/V
b33=(a11*a22-a12*a21)/V

recipro_lattice=np.array([[b11,b12,b13],[b21,b22,b23],[b31,b32,b33]])


f=open('EIGENVAL','r')
lines=f.readlines()
f.close

data={'K':[],'E_up':[]}

for line in lines[7:]:
    values=[float(x) for x in line[:].split()]
    if len(values)==4:
        data['K'].append(values[0:3])
        E_up=[]
    elif len(values)==3:
        E_up.append(values[1])
    elif len(values)==0:
        data['E_up'].append(E_up)
data['E_up'].append(E_up)
data['E_up']=np.array(data['E_up'])
data['K']=np.array(data['K'])
data['K']=np.matmul(data['K'],recipro_lattice)


f=open('KPOINTS','r')
lines=f.readlines()
f.close()
k_miss=[]
k_path=[]
high_sym=[]
ticks=[]
num_high_sym=[]
for num,line in enumerate(lines[3:]):
    tmp=[i for i in line.split()]
    if tmp[3] == '0':
        k_path.append(line)
    else:
        k_miss.append(line)
    if len(tmp) == 6:
        ticks.append(tmp[5])
        high_sym.append(np.array([float(tmp[0]),float(tmp[1]),float(tmp[2])]))
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
np.savetxt('band.dat',band_up,'%9.06f')

high_sym = []
for i in num_high_sym[::2]:
    high_sym.append(band_up[i-len(k_miss), 0])
high_sym.append(band_up[-1, 0])
np.savetxt('high_sym.dat',high_sym,'%9.06f')

tick_tmp = [ticks[0]]
for i in np.arange(1, len(ticks) - 1, 2):
    if ticks[i] == ticks[i + 1]:
        tick_tmp.append(ticks[i])
    else:
        tick_tmp.append(ticks[i] + '|' + ticks[i + 1])
tick_tmp.append(ticks[-1])

ticks = tick_tmp

        
# data['K']=data['K'][len(k_miss):]
# data['E_up']=data['E_up'][len(k_miss):]

# x=np.zeros(len(data['K']))
# for i in range(len(data['K'])-1):
#     x[i+1]=x[i]+np.linalg.norm(data['K'][i+1]-data['K'][i])
# x=x.reshape(len(x),1)
# band_up=np.concatenate((x,data['E_up']),axis=1)
# np.savetxt('band.dat',band_up,'%9.06f')

# high_sym=np.matmul(np.array(high_sym),recipro_lattice)
# high_tmp=np.zeros(len(high_sym))
# for i in range(len(high_sym)-1):
#     high_tmp[i+1]=high_tmp[i]+np.linalg.norm(high_sym[i+1]-high_sym[i])
# high_sym=high_tmp
# np.savetxt('high_sym.dat',high_sym,'%9.06f')


#output Fermi level
efermi=np.array(float(os.popen('grep efermi vasprun.xml |awk \'{print $3}\'').readlines()[0])).reshape(1,1)
np.savetxt('Fermi.dat',efermi,'%9.06f')


f=open('ticks.dat','w')
for i in ticks:
    f.write(i+'\n')
f.close()
