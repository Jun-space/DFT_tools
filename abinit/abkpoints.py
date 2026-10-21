#!/home/dj/software/anaconda3_2019.07/envs/abipy/bin/python
import math
import sys


bohr=0.5291772108
a11=float(sys.argv[1])*bohr
a12=float(sys.argv[2])*bohr
a13=float(sys.argv[3])*bohr
a21=float(sys.argv[4])*bohr
a22=float(sys.argv[5])*bohr
a23=float(sys.argv[6])*bohr
a31=float(sys.argv[7])*bohr
a32=float(sys.argv[8])*bohr
a33=float(sys.argv[9])*bohr
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

intersection=float(sys.argv[10])
k1=int(round(len_b1/intersection))
k2=int(round(len_b2/intersection))
k3=int(round(len_b3/intersection))
if(k1==0):
    k1=1
if(k2==0):
    k2=1
if(k3==0):
    k3=1

print(k1, end=' ')
print(k2, end=' ')
print(k3, end=' ')
