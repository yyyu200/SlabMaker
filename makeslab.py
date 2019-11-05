#!/usr/bin/env python3
# -*- coding: utf-8 -*-

""" @author: yyyu200@163.com """

from CELL import CELL

unitcell=CELL("sc.vasp")

miller_index=[1,1,1]
slab=unitcell.makeslab(miller_index, layer=2, method="brute-force")

slab.print_poscar("./test/slab.vasp")

