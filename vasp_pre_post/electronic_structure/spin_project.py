#!/home/jdeng/software/anaconda3/bin/python
import numpy as np
import os

os.system('grep tot PROCAR | awk \'{print $11}\' > pro.dat')
#grep tot PROCAR | awk '{print $11}' > pro.dat

f=open('pro.dat')
lines=f.readlines()
f.close()


data_x=[]
data_y=[]
data_z=[]
mx=[]
my=[]
mz=[]
for i in np.arange(0,len(lines)):
    if (i+1)%5==3:
        data_x.append(float(lines[i]))
    if (i+1)%5==4:
        data_y.append(float(lines[i]))
    if (i+1)%5==0:
        data_z.append(float(lines[i]))
#nbands=88
#nkpoints=90
nbands=int(input("please input number of NBANDS: "))
nkpoints=int(input("please input number of KPOINTS:"))
data_x=np.array(data_x)
data_x=np.reshape(data_x,(nkpoints,nbands))
data_y=np.array(data_y)
data_y=np.reshape(data_y,(nkpoints,nbands))
data_z=np.array(data_z)
data_z=np.reshape(data_z,(nkpoints,nbands))
for j in np.arange(nbands):
    for i in np.arange(nkpoints):
        mx.append(data_x[i,j])
        my.append(data_y[i,j])
        mz.append(data_z[i,j])

mx=np.array(mx)
mx=np.reshape(mx,(nbands,nkpoints)).T
my=np.array(my)
my=np.reshape(my,(nbands,nkpoints)).T
mz=np.array(mz)
mz=np.reshape(mz,(nbands,nkpoints)).T
np.savetxt('mx.dat',mx,fmt='%.07f')
np.savetxt('my.dat',my,fmt='%.07f')
np.savetxt('mz.dat',mz,fmt='%.07f')
