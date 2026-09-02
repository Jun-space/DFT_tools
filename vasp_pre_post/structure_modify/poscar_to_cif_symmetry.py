# -*- coding: utf-8 -*-
"""
Created on Thu Aug 14 22:26:27 2025

@author: dengjun
"""

from pymatgen.core import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.io.cif import CifWriter

# 读取结构（比如从 POSCAR）
structure = Structure.from_file("K2O-225-Fm-3m.vasp")

# 使用对称性分析器识别空间群（tol=0.01 是对原子位置容差）
sga = SpacegroupAnalyzer(structure, symprec=0.01)
symmetrized_structure = sga.get_refined_structure()

# 导出带有对称性的 CIF
writer = CifWriter(symmetrized_structure, symprec=0.01)
writer.write_file("symmetrized_structure.cif")
