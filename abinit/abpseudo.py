#!/home/dj/software/anaconda3_2019.07/envs/abipy/bin/python
import sys
import os

pseudo=str(sys.argv[1])
znucl=str(sys.argv[2])
z=str(sys.argv[3])
ntypat=int(sys.argv[4])
filename=sys.argv[5]
Z_val=str(sys.argv[6])
typat=str(sys.argv[7])
natom=int(sys.argv[8])
fp=open(filename,'a+')
pseudo=[x for x in pseudo.split()]
znucl=[float(x) for x in znucl.split()]
z=[float(x) for x in z.split()]
Z_val=[float(x) for x in Z_val.split()]
typat=[int(x) for x in typat.split()]
nband=0
for i in range(0,ntypat):
    for j in range(0,ntypat):
        if z[j]==znucl[i]:
            fp.write(str(pseudo[j]))
            fp.write('\n')
            break    
for i in range(0,natom):
    for j in range(0,ntypat):
        if znucl[typat[i]-1]==z[j]:
            nband=nband+Z_val[j]
fp.close()
print(int(nband/2))
