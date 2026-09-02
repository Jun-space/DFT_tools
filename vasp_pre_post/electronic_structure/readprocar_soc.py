#!/home/jdeng/software/anaconda3/bin/python
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

# bcycle=(nions+1)*4+4
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

tmp=[]
for num,line in enumerate(lines):
    t=line.split()
    if t == []:
        continue
    elif t[0] == 'tot':
        tmp.append(num)
    if len(tmp) == 2:
        break
m_cycle=tmp[1]-tmp[0]

kpoints=[]
lines=lines[3:]
mtot=[]
energy=[]
mx=[]
my=[]
mz=[]
for i in range(nkpoints):
    #kpoints.append([float(i) for i in lines[i*kcycle].split()[3:6]])
    kpoints.append([float(lines[i*kcycle][19:30]),float(lines[i*kcycle][30:41]),float(lines[i*kcycle][41:52])])
    klines=lines[2+i*kcycle:(i+1)*kcycle]
    print('kpoints '+str(i))
    for j in range(nbands):
        t=klines[(j*bcycle):(j+1)*bcycle]
        energy.append(float(t[0].split()[4]))
        blines=klines[(3+j*bcycle):(j+1)*bcycle]
        mtot.append(blines[0:nions])
        mx.append(blines[m_cycle:m_cycle+nions])
        my.append(blines[m_cycle*2:m_cycle*2+nions])
        mz.append(blines[m_cycle*3:m_cycle*3+nions])    
        
data={'K':[],'E':[]}
pband={}
for i in range(nions):
    pband[i]={}
    pband[i]['mx']=[]
    pband[i]['my']=[]
    pband[i]['mz']=[]
    pband[i]['mtot']=[]
    print('ions: '+str(i))
    for j in range(nbands):
        for k in range(nkpoints):
            line=k*nbands+j
            m_ix=list(map(float,mx[line][i].split()))
            m_iy=list(map(float,my[line][i].split()))
            m_iz=list(map(float,mz[line][i].split()))
            m_itot=list(map(float,mtot[line][i].split()))
            pband[i]['mx'].append(m_ix)
            pband[i]['my'].append(m_iy)
            pband[i]['mz'].append(m_iz)
            pband[i]['mtot'].append(m_itot)

for i in range(nbands):
    for j in range(nkpoints):
        data['E'].append(energy[j*nbands+i])
data['E']=np.array(data['E'])
for i in range(nions):
    pband[i]['mx']=np.array(pband[i]['mx'])
    pband[i]['my']=np.array(pband[i]['my'])
    pband[i]['mz']=np.array(pband[i]['mz'])
    pband[i]['mtot']=np.array(pband[i]['mtot'])
    pband[i]['mx']=np.delete(pband[i]['mx'],0,axis=1)
    pband[i]['my']=np.delete(pband[i]['my'],0,axis=1)
    pband[i]['mz']=np.delete(pband[i]['mz'],0,axis=1)
    pband[i]['mtot']=np.delete(pband[i]['mtot'],0,axis=1)
    # np.savetxt(str(i+1)+'_mx.dat',pband[i]['mx'],'%.3f')
    # np.savetxt(str(i+1)+'_my.dat',pband[i]['my'],'%.3f')
    # np.savetxt(str(i+1)+'_mz.dat',pband[i]['mz'],'%.3f')
    # np.savetxt(str(i+1)+'_mtot.dat',pband[i]['mtot'],'%.3f')

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
b32=(a11*a23-a13*a21)/V
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


elements=lines[5].split()
nelements=[int(i) for i in lines[6].split()]
pelements={}
for num,i in enumerate(nelements):
    pelements[elements[num]]=np.zeros(np.shape(pband[0]['mz']))
    if num==0:
        begin=0
        end=i
    else:
        begin=np.sum(nelements[:num])
        end=np.sum(nelements[:num+1])
    # print(begin,end)
    for j in range(begin,end):
        pelements[elements[num]]=pelements[elements[num]]+pband[j]['mtot']
    pelements[elements[num]]=np.concatenate((band,pelements[elements[num]]),axis=1)
    np.savetxt(elements[num]+'_soc.dat',pelements[elements[num]],'%.8f')
        
