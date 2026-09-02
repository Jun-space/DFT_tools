#!/bin/bash

prefix=`pwd`
potential=~/Documents/pseudo/PBE_new
if [ ! -f "POTCAR" ]; then
    touch "POTCAR"
else
    rm POTCAR
fi
for i in $@
do
cat $potential/$i/POTCAR >> $prefix/POTCAR
done
