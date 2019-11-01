#!/usr/bin/env python3
# -*- coding: utf-8 -*-

""" @author: yyyu200@163.com """

import numpy as np
from CELL import CELL

unitcell=CELL("sc.vasp")

miller_index=[1,1,0]
slab=unitcell.makeslab(miller_index, layer=1)

slab.print_poscar("./test/slab.vasp")

