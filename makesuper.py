#!/usr/bin/env python3
# -*- coding: utf-8 -*-

""" @author: yyyu200@163.com """
import numpy as np
from CELL import CELL

unitcell=CELL("sc.vasp")
P=np.mat([[2,0,0],[0,2,0],[0,0,2]],dtype=np.float64)    
supercell=CELL.cell2supercell(unitcell,P)
supercell.print_poscar("./tmp/super.vasp")

