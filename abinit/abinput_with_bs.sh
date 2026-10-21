#!/bin/bash

#this shell script could automatically generate the input files 
#for optical calculation in abinit.

#genergate the structure file
abinp.py ebands "$1" > tmp.in
prefix=`echo ${1%%.*}`
inputfile="tmp.in"
output="${prefix}.in"
file="${prefix}.files"
python /home/dj/Documents/abinit/scripts/cif2pos.py $1
#cif2pos.sh $1

#generrate kshift
phonopy --symmetry --tolerance=1e-2 > spacegroup.dat
lattice=`grep space_group_type spacegroup.dat| awk '{print $2}'|cut -c2-2`
length=`grep point_group_type spacegroup.dat| awk '{print $2}'| awk '{print length($0)}'`
point_group=`grep point_group_type spacegroup.dat| awk '{print $2}'| cut -c2-$[$length-1]`

#how to generate the kpoints
kptopt=`grep "kptopt" $inputfile | awk '{print $2}'`

#electronic configuration
prtden=1
toldfe=1.0d-8
ecut=`grep "ecut" $inputfile | awk '{print $2}'`
nstep=100
#nband=`grep nband $inputfile | awk '{print $2}'`

#structure
a11=`grep -w -A 1 "rprim" $inputfile| tail -1| awk '{print $1}'`
a12=`grep -w -A 1 "rprim" $inputfile| tail -1| awk '{print $2}'`
a13=`grep -w -A 1 "rprim" $inputfile| tail -1| awk '{print $3}'`
a21=`grep -w -A 2 "rprim" $inputfile| tail -1| awk '{print $1}'`
a22=`grep -w -A 2 "rprim" $inputfile| tail -1| awk '{print $2}'`
a23=`grep -w -A 2 "rprim" $inputfile| tail -1| awk '{print $3}'`
a31=`grep -w -A 3 "rprim" $inputfile| tail -1| awk '{print $1}'`
a32=`grep -w -A 3 "rprim" $inputfile| tail -1| awk '{print $2}'`
a33=`grep -w -A 3 "rprim" $inputfile| tail -1| awk '{print $3}'`
natom=`grep "natom" $inputfile | awk '{print $2}'`
ntypat=`grep "ntypat" $inputfile | awk '{print $2}'`
xred=`grep -w -A $natom "xred" $inputfile`
acell=`grep -w "acell" $inputfile`

if [ $natom -le 3 ]; then
typat=`grep -w "typat" $inputfile|sed 's/ typat//g'`
else
line=`python3 /home/dj/Documents/abinit/scripts/abceil.py $natom`
typat=`grep -w -A $line "typat" $inputfile |sed 's/ typat//g'`
fi

if [ $ntypat -le 3 ]; then
znucl=`grep znucl $inputfile | sed 's/znucl //g'`
else
line=`python3 /home/dj/Documents/abinit/scripts/abceil.py $ntypat`
znucl=`grep -w -A $line znucl $inputfile |sed 's/ znucl//g'`
fi

#kgrid
#kpoints.py 0.03
#ngkpt=`sed -n '4,4p' KPOINTS`
ngkpt=`python3 /home/dj/Documents/abinit/scripts/abkpoints.py $a11 $a12 $a13 $a21 $a22 $a23 $a31 $a32 $a33 0.02`

#ebands
kptopt7=`grep kptopt2 $inputfile | awk '{print $2}'`
kptbounds7=`grep -w -A $[kptopt7*(-1)+1] "kptbounds2" $inputfile`



#check file existence, if existence the input file ,then delete it
if [ ! -f "$output" ]; then
    touch "$output"  
else
    rm $output
fi

if [ ! -f "$file" ]; then
    touch "$file"
else
    rm $file
fi




########################### wirte the *.files file #############################
cat >> $file << end
${prefix}.in
${prefix}.out
${prefix}_i
${prefix}_o
${prefix}
end
z=`grep "\"Z\":" $inputfile| awk '{print $3}'| sed 's/.$//'`
Z_val=`grep Z_val $inputfile |awk '{print $3}'|sed 's/,//g'`
pseudo=`grep "filepath" $inputfile | sed 's/#            "filepath": "//g'|sed 's/",//g'`
#python3 /home/dj/Documents/scripts/pseudo.py "$pseudo" "$znucl" "$z" "$ntypat" "$file" "${Z_val}" "$typat" "$natom"
nband="`python3 /home/dj/Documents/abinit/scripts/abpseudo.py "$pseudo" "$znucl" "$z" "$ntypat" "$file" "${Z_val}" "$typat" "$natom"`"


########################  write the *.in file  ################################

#write sixdatabase
cat >> $output << EOF

ndtset 7

#First dataset : SC run with kpoints in the IBZ
nband1 $nband
nstep1 $nstep
kptopt1 1
nbdbuf1 0
prtden1 1
getden1 0
getwfk1 0

#Second dataset : NSC run with large number of bands, and points in the IBZ
iscf2 -2
nband2 $[nband*4]
nstep2 $nstep
kptopt2 1
getwfk2 1
getden2 1

#Third dataset : NSC run with large number of bands, and points in the the full BZ
iscf3 -2
nband3 $[nband*4]
nstep3 100
kptopt3 3
getwfk3 2
getden3 1 

#Fourth dataset : ddk response function along axis 1
iscf4 -3
nband4 $[nband*4]
nstep4 1
nline4 0                   #It must be zero. Gives maximum number of line minimizations allowed in preconditioned conjugate gradient minimization for each band.
prtwf4 3
kptopt4 3
nqpt4 1                    #Determines whether one q point must be read
qpt4 0.0d0 0.0d0 0.0d0     #Only used if nqpt = 1.
rfdir4 1 0 0               #the directions to be considered for response function calculations 
rfelfd4 2                  #Turns on electric field response function calculations.2 → only the derivative of ground-state wavefunctions with respect to k
getwfk4 3

#Fifth dataset : ddk response function along axis 2
iscf5 -3
nband5 $[nband*4]
nstep5 1
nline5 0
prtwf5 3
kptopt5 3
nqpt5 1
qpt5 0.0d0 0.0d0 0.0d0
rfdir5  0 1 0
rfelfd5  2
getwfk5  3

#Sixth dataset : ddk response function along axis 3
iscf6 -3
nband6 $[nband*4]
nstep6 1
nline6 0
prtwf6 3
kptopt6 3
nqpt6 1 
qpt6 0.0d0 0.0d0 0.0d0
rfdir6 0 0 1
rfelfd6 2
getwfk6 3

#Seventh dataset : band structure
iscf7   -2
getden7 -6
kptopt7 $kptopt7
nband7 $[nband+15]
enunit7 1
ndivsm7 30
$kptbounds7
EOF

#determin the k shift type
type=`python3 /home/dj/Documents/abinit/scripts/abspace_group.py $lattice $point_group`
if [ $type == 1 ]; then  # FCC
nshiftk=4
shiftk1="0.5 0.5 0.5"
shiftk2="0.5 0.0 0.0"
shiftk3="0.0 0.5 0.0"
shiftk4="0.0 0.0 0.5"
cat >> $output << EOF

#kpoints
nshiftk $nshiftk
shiftk $shiftk1
       $shiftk2
       $shiftk3
       $shiftk4
EOF
elif [ $type == 2 ] ; then #BCC
nshiftk=2
shiftk1="0.25 0.25 0.25"
shiftk2="-0.25 -0.25 -0.25"
cat >> $output << EOF
nshiftk $nshiftk
shiftk $shiftk1
       $shiftk2
EOF
elif [ $type == 3 ]; then #hexagonal
cat >> $output << EOF
nshiftk 1
shiftk  0.0 0.0 0.5
EOF
elif [ $type == 4 ]; then # R lattice
cat >> $output << EOF
nshiftk 1
shiftk 0.5 0.5 0.5
EOF
else                     # others
nshiftk=1
shiftk="0.5 0.5 0.5"
cat >> $output << EOF

#kpoints
nshiftk $nshiftk
shiftk  $shiftk
chksymbreak 0
EOF
fi


#kpoints
cat >> $output << !
ngkpt $ngkpt
!

cat >> $output << EOF

#structure
natom $natom
ntypat $ntypat
typat $typat
znucl $znucl
$xred
$acell
rprim
      $a11 $a12 $a13
      $a21 $a22 $a23
      $a31 $a32 $a33

#electronic relaxation
diemac 5.0
ecut $ecut
ixc 11
#tnons
#toldfe 1.0d-8
tolwfr 1e-20
tolsym 1e-3
EOF

rm tmp.in BPOSCAR PPOSCAR spacegroup.dat POSCAR 
sed -i 's/ kptbounds2/kptbounds7/g' $output 
