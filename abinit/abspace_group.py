#!/home/dj/software/anaconda3_2019.07/envs/abipy/bin/python
import sys
import os

lattice=sys.argv[1]
point_group=sys.argv[2]

if (lattice=="F") and ((point_group=="23") or (point_group=="m-3") or (point_group=="432") or (point_group=="-43m") or (point_group=="m-3m") ):
    print(1)
elif (lattice=="I") and ((point_group=="23") or (point_group=="m-3") or (point_group=="432") or (point_group=="-43m") or (point_group=="m-3m") ):
    print(2)
elif lattice=="R":
    print(4)
elif (point_group=="3") or (point_group=="-3") or (point_group=="32") or (point_group=="3m") or (point_group=="-3m") or (point_group=="6") or (point_group=="-6") or (point_group=="6/m") or (point_group=="622") or (point_group=="6mm") or (point_group=="-62m") or (point_group=="6/mmm"):
    print(3)
else:
    print(5)
