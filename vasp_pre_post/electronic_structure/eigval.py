#!/home/software/anaconda3/bin/python
 
import math
import sys 
import os
import numpy as np
import re

f=open('POSCAR','r')
lines=f.readlines()
f.close
 
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
b32=(a11*a23-a13*a21)/V
b33=(a11*a22-a12*a21)/V

recipro_lattice=np.array([[b11,b12,b13],[b21,b22,b23],[b31,b32,b33]])


f=open('EIGENVAL','r')
lines=f.readlines()
f.close

data={'K':[],'E':[]}

for line in lines[7:]:
    values=[float(x) for x in line[:].split()]
    if len(values)==4:
        data['K'].append(values[0:3])
        E=[]
    elif len(values)==3:
        E.append(values[1])
    elif len(values)==0:
        data['E'].append(E)
data['E'].append(E)
data['E']=np.array(data['E'])
data['K']=np.array(data['K'])

data['K']=np.matmul(data['K'],recipro_lattice)
x=np.zeros(len(data['K']))
for i in range(len(data['K'])-1):
    x[i+1]=x[i]+np.linalg.norm(data['K'][i+1]-data['K'][i])
x=x.reshape(len(x),1)
band=np.concatenate((x,data['E']),axis=1)
np.savetxt('band.dat',band,'%9.06f')


f=open('KPOINTS','r')
lines=f.readlines()
f.close()
kspace=int(lines[1])
high_sym_num=(len(lines[4:])+1)/3+1


high_sym=[]
high_sym.append(float(x[0]))
for i in range(1,int(high_sym_num)):
    high_sym.append(float(x[kspace*i-1]))

high_sym=np.array(high_sym)

np.savetxt('high_sym.dat',high_sym,'%9.06f')


#output Fermi level
efermi=np.array(float(os.popen('grep efermi vasprun.xml |awk \'{print $3}\'').readlines()[0])).reshape(1,1)
np.savetxt('Fermi.dat',efermi,'%9.06f')

#output ticks
#output ticks
f=open("KPOINTS")
KPread=f.read()
# k2=f.readlines()
f.close()
tick_tmp=np.array(re.findall("!\s(.*)",KPread))
tick=tick_tmp[::2]
tick=list(str(i) for i in tick)
tick.append(str(tick_tmp[-1]))

f=open('ticks.dat','w')
for i in tick:
    f.write(i+'\n')
f.close()


