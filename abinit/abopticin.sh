#!/bin/bash

prefix=`echo ${1%%.*}`

wfkfile=${prefix}_o_DS3_WFK
ddkfile_1=`ls ${prefix}_o_DS4_1WF*`
ddkfile_2=`ls ${prefix}_o_DS5_1WF*`
ddkfile_3=`ls ${prefix}_o_DS6_1WF*`

#obtain the point group 
python /home/dj/Documents/abinit/scripts/cif2pos.py $1
phonopy --symmetry > spacegroup.dat
length=`grep point_group_type spacegroup.dat| awk '{print $2}'| awk '{print length($0)}'`
point_group=`grep point_group_type spacegroup.dat| awk '{print $2}'| cut -c2-$[$length-1]`


#check file existence, if existence the input file, then delete it
if [ ! -f "${prefix}_optic.in" ]; then
    touch "${prefix}_optic.in"
else
    rm ${prefix}_optic.in
fi

if [ ! -f "${prefix}_optic.files" ]; then
    touch "${prefix}_optic.files"                                                                            
else
    rm ${prefix}_optic.files
fi


#############  Write the optical *.in  file  ##############

cat >> ${prefix}_optic.in <<EOF
&FILES
 ddkfile_1 = '${ddkfile_1}',
 ddkfile_2 = '${ddkfile_2}',
 ddkfile_3 = '${ddkfile_3}',
 wfkfile = '${wfkfile}'
/
&PARAMETERS
 broadening = 0.002,
 domega = 0.0003,
 maxomega = 0.5,
 scissor = 0.000,
 tolerance = 0.002
/
&COMPUTATIONS
 num_lin_comp = 3,
 lin_comp = 11,22,33,
EOF

python3 /home/dj/Documents/abinit/scripts/abi_optic_coefficient.py ${point_group} >> ${prefix}_optic.in
echo / >> ${prefix}_optic.in

#############  Write the optical *.files  file #############

cat >> ${prefix}_optic.files << !
${prefix}_optic.in
${prefix}_optic.out
${prefix}
!



rm spacegroup.dat POSCAR PPOSCAR BPOSCAR phonopy_symcells.yaml
