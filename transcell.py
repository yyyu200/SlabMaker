#!/usr/bin/env python3
# -*- coding: utf-8 -*-

""" @author: yyyu200@163.com """

from CELL import CELL

# test 4 kinds of unit to primitive cell transform

sc=CELL("sc.vasp")
fcc=CELL.unit2prim(sc,2)
fcc.print_poscar("fcc.vasp")

in2o3_unit=CELL("In2O3.vasp")
bcc=CELL.unit2prim(in2o3_unit,3)
bcc.print_poscar("bcc.vasp")

hexcell=CELL("Al2O3.vasp")
rhcell=CELL.unit2prim(hexcell,5)
rhcell.print_poscar("rh.vasp")

c1=CELL("mC.vasp")
prim=CELL.unit2prim(c1,13)
prim.print_poscar("prim.vasp")

