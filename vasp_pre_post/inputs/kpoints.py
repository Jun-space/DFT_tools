#!/publicfs10/fs10-m9/home/m9s005343/software/anaconda3/bin/python

import math
import sys
import os

f=open('POSCAR','r')
lines=f.readlines()
f.close

ratio=float(lines[1])
a=[x for x in lines[2].split()]
b=[x for x in lines[3].split()]
c=[x for x in lines[4].split()]
a11=float(a[0])*ratio
a12=float(a[1])*ratio
a13=float(a[2])*ratio
a21=float(b[0])*ratio
a22=float(b[1])*ratio
a23=float(b[2])*ratio
a31=float(c[0])*ratio
a32=float(c[1])*ratio
a33=float(c[2])*ratio
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

len_b1=math.sqrt(b11**2+b12**2+b13**2)
len_b2=math.sqrt(b21**2+b22**2+b23**2)
len_b3=math.sqrt(b31**2+b32**2+b33**2)

intersection=float(sys.argv[1])
k1=int(round(len_b1/intersection))
k2=int(round(len_b2/intersection))
k3=int(round(len_b3/intersection))

if(k1==0):
    k1=1
if(k2==0):
    k2=1
if(k3==0):
    k3=1

if(os.path.exists("KPOINTS")):
    os.remove("KPOINTS")
os.mknod("KPOINTS")
fp=open('KPOINTS','w')
fp.write(str('Auto mesh: kspacing ~')+ sys.argv[1]+'\n')
fp.write(str('0 \n'))
fp.write(str('Gamma \n'))
fp.write(str(k1))
fp.write(" ")
fp.write(str(k2))
fp.write(" ")
fp.write(str(k3))
fp.write('\n')
fp.write(str('0.0 0.0 0.0'))
fp.close

