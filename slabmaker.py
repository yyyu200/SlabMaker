#!/usr/bin/env python3
# -*- coding: utf-8 -*-

""" @author: yyyu200@163.com """

import numpy as np
from CELL import CELL

unitcell=CELL("Al2O3.vasp")

miller_index=[1,1,0]
slab=unitcell.makeslab(miller_index, length=50.0)

slab.print_poscar("slab.vasp")

