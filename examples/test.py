#!/usr/bin/env python3
# -*- coding: utf-8 -*-

""" @author: yyyu200@163.com """

from SlabMaker.build import CELL

outpath="./results/"
import os

if not os.path.exists(outpath):
    os.makedirs(outpath)

# test 4 kinds of unit to primitive cell transform
sc=CELL("SiC.vasp")
fcc=CELL.unit2prim(sc,2)
fcc.print_poscar(outpath+"SiC-fcc.vasp")

in2o3_unit=CELL("In2O3.vasp")
bcc=CELL.unit2prim(in2o3_unit,3)
bcc.print_poscar(outpath+"In2O3-bcc.vasp")

hexcell=CELL("Al2O3.vasp")
rhcell=CELL.unit2prim(hexcell,5)
rhcell.print_poscar(outpath+"Al2O3-rh.vasp")

c1=CELL("mC.vasp")
prim=CELL.unit2prim(c1,13)
prim.print_poscar(outpath+"mC-bcm.vasp")

# test build slab
unitcell=CELL("SiC.vasp")
miller_index=[1,1,1]
slab=unitcell.makeslab(miller_index, layer=2, method="brute-force")
slab.print_poscar(outpath+"slab.vasp")

# build supercell
unitcell=CELL("SiC.vasp")
supercell=CELL.cell2supercell(unitcell,[[2,0,0],[0,2,0],[0,0,2]])
supercell.print_poscar(outpath+"super.vasp")

# test QE to POSCAR format
sno=CELL('SnO2.inp',fmt='QE')
sno.print_poscar(outpath+"SnO2.vasp")

