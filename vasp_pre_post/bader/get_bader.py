from pymatgen.io.vasp import Outcar
from pymatgen.core.structure import Structure
import numpy as np


def get_bader_charge(path):
    ##########get atomic labels from POSCAR
    s=Structure.from_file(path+'/POSCAR').as_dict()
    outcar=Outcar(path+'/OUTCAR')
    outcar.read_pseudo_zval()
    zval_dict=outcar.zval_dict
    labels=[]
    for i in s['sites']:
        labels.append(i['label'])
    
    ###########total charge
    nele=outcar.nelect
    
    ##########get bader charges from ACF.dat
    f=open(path+'/ACF.dat')
    lines=f.readlines()
    f.close()
    bader_charge=[]
    for line in lines[2:len(lines)-4]:
        bader_charge.append(float(line.split()[4]))
    
    ##########get charge state, positive values mean positive charged, and vice versa
    charge=[]
    for num,i in enumerate(labels):
        charge.append(zval_dict[labels[num]]-bader_charge[num])
        
    #data={'charge':charge,'bader_charge':bader_charge,'potcar_zval':zval_dict,'labels':labels,'total_charge':nele}
    
    charge=np.array(charge)
    np.savetxt('bader.dat',charge,fmt='%10.8f')
    return

if __name__ == '__main__':
    get_bader_charge('./')
