#!/edahome/y00201/software/anaconda3/bin/python
import numpy as np
import copy

print('reading PROCAR')
f=open('PROCAR')
lines=f.readlines()
f.close()
print('reading PROCAR complete')
nkpoints=int(lines[1].split()[3])
nbands=int(lines[1].split()[7])
nions=int(lines[1].split()[11])

# bcycle=(nions+1)+4
# kcycle=bcycle*nbands+3
tmp=[]
for num,line in enumerate(lines):
    t=line.split()
    if t == []:
        continue
    elif t[0] == 'band':
        tmp.append(num)
    if len(tmp) == 2:
        break
bcycle=tmp[1]-tmp[0]

tmp=[]
for num,line in enumerate(lines):
    t=line.split()
    if t == []:
        continue
    elif t[0] == 'k-point':
        tmp.append(num)
    if len(tmp) == 2:
        break
kcycle=tmp[1]-tmp[0]

kpoints=[]
lines=lines[3:]
mtot=[]
energy=[]
pro=[]
for i in range(nkpoints):
#    kpoints.append([float(i) for i in lines[i*kcycle].split()[3:6]])
    kpoints.append([float(lines[i*kcycle][19:30]),float(lines[i*kcycle][30:41]),float(lines[i*kcycle][41:52])])
    klines=lines[2+i*kcycle:(i+1)*kcycle]
    print('kpoints '+str(i+1))
    for j in range(nbands):
        t=klines[(j*bcycle):(j+1)*bcycle]
        energy.append(float(t[0].split()[4]))
        blines=klines[(3+j*bcycle):(j+1)*bcycle]
        pro.append(blines[0:nions])   
        
data={'K':[],'E':[]}
pband={}
for i in range(nions):
    pband[i]=[]
    print('ions: '+str(i))
    for j in range(nbands):
        for k in range(nkpoints):
            line=k*nbands+j
            pro_i=list(map(float,pro[line][i].split()))
            pband[i].append(pro_i)

for i in range(nbands):
    for j in range(nkpoints):
        data['E'].append(energy[j*nbands+i])
data['E']=np.array(data['E'])
for i in range(nions):
    pband[i]=np.array(pband[i])
    pband[i]=np.delete(pband[i],0,axis=1)
    # np.savetxt(str(i+1)+'.dat',pband[i],'%.3f')

f=open('POSCAR')
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

data['K']=np.matmul(kpoints,recipro_lattice)
x=np.zeros(len(data['K']))
for i in range(len(data['K'])-1):
    x[i+1]=x[i]+np.linalg.norm(data['K'][i+1]-data['K'][i])
x=x.reshape(len(x),1)
x_tmp=copy.deepcopy(x)
x=x_tmp
for i in range(nbands-1):
    x=np.concatenate((x,x_tmp),axis=0)
band=np.concatenate((x,data['E'].reshape(len(data['E']),1)),axis=1)

header='kpoints        Energy          s         py         pz         px        dxy        dyz        dz2        dxz        dx2      Total'
for i in range(nions):
    t=np.concatenate((band,pband[i]),axis=1)
    np.savetxt(str(i+1)+'.dat',t,'%.8f',header=header)

elements=lines[5].split()
nelements=[int(i) for i in lines[6].split()]
pelements={}
for num,i in enumerate(nelements):
    pelements[elements[num]]=np.zeros(np.shape(pband[0]))
    if num==0:
        begin=0
        end=i
    else:
        begin=np.sum(nelements[:num])
        end=np.sum(nelements[:num+1])
    # print(begin,end)
    for j in range(begin,end):
        pelements[elements[num]]=pelements[elements[num]]+pband[j]
    pelements[elements[num]]=np.concatenate((band,pelements[elements[num]]),axis=1)
    np.savetxt(elements[num]+'.dat',pelements[elements[num]],'%.8f',header=header)
        
