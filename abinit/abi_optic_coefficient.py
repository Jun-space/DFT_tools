#!/home/dj/software/anaconda3_2019.07/envs/abipy/bin/python
import sys
import os

point_group=sys.argv[1]

if point_group=="1":
    print("num_nonlin_comp=18,")
    print("nonlin_comp=111,122,133,123,113,112,211,222,233,223,213,212,311,322,333,323,313,312,")

if point_group=='2':
    print("num_nonlin_comp=8,")
    print("nonlin_comp=123,112,211,222,233,213,323,312,")

if point_group=="m":
    print("num_nonlin_comp=10,")
    print("nonlin_comp=111,122,133,113,223,212,311,322,333,313,")

if point_group=="222":
    print("num_nonlin_comp=3,")
    print("nonlin_comp=123,213,312,")

if point_group=="mm2":
    print("num_nonlin_comp=5,")
    print("nonlin_comp=113,223,311,322,333,")

if point_group=="4":
    print("num_nonlin_comp=4")
    print("nonlin_comp=123,113,311,333")

if point_group=="-4":
    print("num_nonlin_comp=4")
    print("nonlin_comp=123,113,311,312,")

if point_group=="422":
    print("num_nonlin_comp=1")
    print("nonlin_comp=123,")

if point_group=="-42m":
    print("num_nonlin_comp=2,")
    print("nonlin_comp=123,112,")

if point_group=="4mm":
    print("num_nonlin_comp=3,")
    print("nonlin_comp=113,311,333,")

if point_group=="3":
    print("num_nonlin_comp=6,")
    print("nonlin_comp=111,123,113,112,311,333,")

if point_group=="32":
    print("num_nonlin_comp=2,")
    print("nonlin_comp=111,123,")

if point_group=="3m":
    print("num_nonlin_comp=4,")
    print("nonlin_comp=113,222,311,333,")

if point_group=="6":
    print("num_nonlin_comp=4,")
    print("nonlin_comp=123,113,311,333,")

if point_group=="6mm":
    print("num_nonlin_comp=3,")
    print("nonlin_comp=113,311,333,")

if point_group=="622":
    print("num_nonlin_comp=1,")
    print("nonlin_comp=123,")

if point_group=="-6":
    print("num_nonlin_comp=2,")
    print("nonlin_comp=111,112,")

if point_group=="-6m2":
    print("num_nonlin_comp=1,")
    print("nonlin_comp=112,")

if point_group=="23" or point_group=="-43m":
    print("num_nonlin_comp=1,")
    print("nonlin_comp=123,")





















