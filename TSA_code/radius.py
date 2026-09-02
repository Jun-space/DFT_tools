# -*- coding: utf-8 -*-
"""
Created on Sat Feb 12 10:26:53 2022

@author: DengJun

define the covalent radius for atoms
"""


def radius():
    r_covalent={1:{'label':'H','radius':0.31},2:{'label':'He','radius':0.28},3:{'label':'Li','radius':1.28},
                4:{'label':'Be','radius':0.96},5:{'label':'B','radius':0.84},6:{'label':'C','radius':0.76},
                7:{'label':'N','radius':0.71},8:{'label':'O','radius':0.66},9:{'label':'F','radius':0.57},
                10:{'label':'Ne','radius':0.58},11:{'label':'Na','radius':1.66},12:{'label':'Mg','radius':1.41},
                13:{'label':'Al','radius':1.21},14:{'label':'Si','radius':1.11},15:{'label':'P','radius':1.07},
                16:{'label':'S','radius':1.05},17:{'label':'Cl','radius':1.02},18:{'label':'Ar','radius':1.06},
                19:{'label':'K','radius':2.03},20:{'label':'Ca','radius':1.76},21:{'label':'Sc','radius':1.70},
                22:{'label':'Ti','radius':1.60},23:{'label':'V','radius':1.53},24:{'label':'Cr','radius':1.39},
                25:{'label':'Mn','radius':1.39},26:{'label':'Fe','radius':1.32},27:{'label':'Co','radius':1.26},
                28:{'label':'Ni','radius':1.24},29:{'label':'Cu','radius':1.32},30:{'label':'Zn','radius':1.22},
                31:{'label':'Ga','radius':1.22},32:{'label':'Ge','radius':1.20},33:{'label':'As','radius':1.19},
                34:{'label':'Se','radius':1.20},35:{'label':'Br','radius':1.20},36:{'label':'Kr','radius':1.16},
                37:{'label':'Rb','radius':2.20},38:{'label':'Sr','radius':1.95},39:{'label':'Y','radius':1.90},
                40:{'label':'Zr','radius':1.75},41:{'label':'Nb','radius':1.64},42:{'label':'Mo','radius':1.54},
                43:{'label':'Tc','radius':1.47},44:{'label':'Ru','radius':1.46},45:{'label':'Rh','radius':1.42},
                46:{'label':'Pd','radius':1.39},47:{'label':'Ag','radius':1.45},48:{'label':'Cd','radius':1.44},
                49:{'label':'In','radius':1.42},50:{'label':'Sn','radius':1.39},51:{'label':'Sb','radius':1.39},
                52:{'label':'Te','radius':1.38},53:{'label':'I','radius':1.39},54:{'label':'Xe','radius':1.40},
                55:{'label':'Cs','radius':2.44},56:{'label':'Ba','radius':2.15},57:{'label':'La','radius':2.07},
                58:{'label':'Ce','radius':2.04},59:{'label':'Pr','radius':2.03},60:{'label':'Nd','radius':2.01},
                61:{'label':'Pm','radius':1.99},62:{'label':'Sm','radius':1.98},63:{'label':'Eu','radius':1.98},
                64:{'label':'Gd','radius':1.96},65:{'label':'Tb','radius':1.94},66:{'label':'Dy','radius':1.92},
                67:{'label':'Ho','radius':1.92},68:{'label':'Er','radius':1.89},69:{'label':'Tm','radius':1.90},
                70:{'label':'Yb','radius':1.87},71:{'label':'Lu','radius':1.87},72:{'label':'Hf','radius':1.75},
                73:{'label':'Ta','radius':1.70},74:{'label':'W','radius':1.62},75:{'label':'Re','radius':1.51},
                76:{'label':'Os','radius':1.44},77:{'label':'Ir','radius':1.41},78:{'label':'Pt','radius':1.36},
                79:{'label':'Au','radius':1.36},80:{'label':'Hg','radius':1.32},81:{'label':'Tl','radius':1.45},
                82:{'label':'Pb','radius':1.46},83:{'label':'Bi','radius':1.48},84:{'label':'Po','radius':1.40},
                85:{'label':'At','radius':1.50},86:{'label':'Rn','radius':1.50},87:{'label':'Fr','radius':2.60},
                88:{'label':'Ra','radius':2.21},89:{'label':'Ac','radius':2.15},90:{'label':'Th','radius':2.06},
                91:{'label':'Pa','radius':2.00},92:{'label':'U','radius':1.96},93:{'label':'Np','radius':1.90},
                94:{'label':'Pu','radius':1.87},95:{'label':'Am','radius':1.80},96:{'label':'Cm','radius':1.69},
                97:{'label':'Bk','radius':None},98:{'label':'Cf','radius':None},99:{'label':'Es','radius':None},
                }
    
    
    r_vdw={1:{'label':'H','radius':1.20},2:{'label':'He','radius':1.43},3:{'label':'Li','radius':2.12},
                4:{'label':'Be','radius':1.98},5:{'label':'B','radius':1.91},6:{'label':'C','radius':1.77},
                7:{'label':'N','radius':1.66},8:{'label':'O','radius':1.50},9:{'label':'F','radius':1.46},
                10:{'label':'Ne','radius':1.58},11:{'label':'Na','radius':2.50},12:{'label':'Mg','radius':2.51},
                13:{'label':'Al','radius':2.25},14:{'label':'Si','radius':2.19},15:{'label':'P','radius':1.90},
                16:{'label':'S','radius':1.89},17:{'label':'Cl','radius':1.82},18:{'label':'Ar','radius':1.83},
                19:{'label':'K','radius':2.73},20:{'label':'Ca','radius':2.62},21:{'label':'Sc','radius':2.58},
                22:{'label':'Ti','radius':2.46},23:{'label':'V','radius':2.42},24:{'label':'Cr','radius':2.45},
                25:{'label':'Mn','radius':2.45},26:{'label':'Fe','radius':2.44},27:{'label':'Co','radius':2.40},
                28:{'label':'Ni','radius':2.40},29:{'label':'Cu','radius':2.38},30:{'label':'Zn','radius':2.39},
                31:{'label':'Ga','radius':2.32},32:{'label':'Ge','radius':2.29},33:{'label':'As','radius':1.88},
                34:{'label':'Se','radius':1.82},35:{'label':'Br','radius':1.86},36:{'label':'Kr','radius':2.25},
                37:{'label':'Rb','radius':3.21},38:{'label':'Sr','radius':2.84},39:{'label':'Y','radius':2.75},
                40:{'label':'Zr','radius':2.52},41:{'label':'Nb','radius':2.56},42:{'label':'Mo','radius':2.45},
                43:{'label':'Tc','radius':2.44},44:{'label':'Ru','radius':2.46},45:{'label':'Rh','radius':2.44},
                46:{'label':'Pd','radius':2.15},47:{'label':'Ag','radius':2.53},48:{'label':'Cd','radius':2.49},
                49:{'label':'In','radius':2.43},50:{'label':'Sn','radius':2.42},51:{'label':'Sb','radius':2.47},
                52:{'label':'Te','radius':1.99},53:{'label':'I','radius':2.04},54:{'label':'Xe','radius':2.06},
                55:{'label':'Cs','radius':3.48},56:{'label':'Ba','radius':3.03},57:{'label':'La','radius':2.98},
                58:{'label':'Ce','radius':2.88},59:{'label':'Pr','radius':2.92},60:{'label':'Nd','radius':2.95},
                61:{'label':'Pm','radius':None},62:{'label':'Sm','radius':2.90},63:{'label':'Eu','radius':2.87},
                64:{'label':'Gd','radius':2.83},65:{'label':'Tb','radius':2.79},66:{'label':'Dy','radius':2.87},
                67:{'label':'Ho','radius':2.81},68:{'label':'Er','radius':2.83},69:{'label':'Tm','radius':2.79},
                70:{'label':'Yb','radius':2.80},71:{'label':'Lu','radius':2.74},72:{'label':'Hf','radius':2.63},
                73:{'label':'Ta','radius':2.53},74:{'label':'W','radius':2.57},75:{'label':'Re','radius':2.49},
                76:{'label':'Os','radius':2.48},77:{'label':'Ir','radius':2.41},78:{'label':'Pt','radius':2.29},
                79:{'label':'Au','radius':2.32},80:{'label':'Hg','radius':2.45},81:{'label':'Tl','radius':2.47},
                82:{'label':'Pb','radius':2.60},83:{'label':'Bi','radius':2.54},84:{'label':'Po','radius':None},
                85:{'label':'At','radius':None},86:{'label':'Rn','radius':None},87:{'label':'Fr','radius':None},
                88:{'label':'Ra','radius':None},89:{'label':'Ac','radius':2.8},90:{'label':'Th','radius':2.93},
                91:{'label':'Pa','radius':2.88},92:{'label':'U','radius':2.71},93:{'label':'Np','radius':2.82},
                94:{'label':'Pu','radius':2.81},95:{'label':'Am','radius':2.83},96:{'label':'Cm','radius':3.05},
                97:{'label':'Bk','radius':3.40},98:{'label':'Cf','radius':3.05},99:{'label':'Es','radius':2.70}
                }

    return r_covalent,r_vdw