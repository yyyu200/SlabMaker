# SlabMaker

Functions:

1. Transformation between primitive cell and unit cell.

2. Build slab models for crystal surfaces.

3. Build supercell given the transform matrix P.

4. Convert PW input file to POSCAR.

## Build primitive cell from unit cell

Input unit cell and output slab models

use VASP POSCAR format
```
Si C Unit cell: comment goes here
1.0
        4.3480000496         0.0000000000         0.0000000000
        0.0000000000         4.3480000496         0.0000000000
        0.0000000000         0.0000000000         4.3480000496
   Si    C
    4    4
Direct
     0.000000000         0.000000000         0.000000000
     0.000000000         0.500000000         0.500000000
     0.500000000         0.000000000         0.500000000
     0.500000000         0.500000000         0.000000000
     0.250000000         0.250000000         0.250000000
     0.750000000         0.750000000         0.250000000
     0.750000000         0.250000000         0.750000000
     0.250000000         0.750000000         0.750000000

```
save above as file 'sc.vasp'. Run

```python
from build import CELL

sc=CELL("examples/SiC.vasp")
fcc=CELL.unit2prim(sc,2)
fcc.print_poscar("fcc.vasp")
```

Will get the primitive cell file 'fcc.vasp',

```
system: Si C
1.0
 0.000000 2.174000 2.174000
 2.174000 0.000000 2.174000
 2.174000 2.174000 0.000000
Si C
1 1
Direct
 0.000000000000 0.000000000000 0.000000000000
 0.250000000000 0.250000000000 0.250000000000
```

## Build slab

```python
from build import CELL

unitcell=CELL("examples/SiC.vasp")
miller_index=[1,1,1]
slab=unitcell.makeslab(miller_index, layer=1, vacuum=15.0)
slab.print_poscar("./slab.vasp")
```

will get 
```
P1 =  [[ 1.  0.  1.]
 [ 0.  1.  1.]
 [-1. -1.  1.]]
P2 =  [[ 0.5 -0.5  0. ]
 [ 0.   0.5  0. ]
 [ 0.   0.   1. ]]
reduced slab cell
 [[ 3.07450032  0.          0.        ]
 [-1.53725016  2.66259538  0.        ]
 [ 0.          0.         21.90337725]]
reduced slab No. of atoms:  6
slab and vacuum length:  6.903377247450933 15.0 Ang.
inplane edge and angle:  3.074500319671605 3.074500319671605  Ang.  120.00000000000001  degree.
reduced slab cell area:  8.186150349361135  Ang^2.
```
and also the slab.vasp file:

```
system: Si C
1.0
    3.0745003197    0.0000000000    0.0000000000
   -1.5372501598    2.6625953808    0.0000000000
    0.0000000000    0.0000000000   21.9033772475
Si C
3 3
Direct
 0.000000000000 0.000000000000 0.342412949166
 0.666666666667 0.333333333333 0.457021713409
 0.333333333333 0.666666666667 0.571630477652
 0.000000000000 0.000000000000 0.428369522348
 0.666666666667 0.333333333333 0.542978286591
 0.333333333333 0.666666666667 0.657587050834
```

## Build supercell:

```python
import numpy as np
from build import CELL

unitcell=CELL("examples/SiC.vasp")
P=np.mat([[2,0,0],[0,2,0],[0,0,2]],dtype=np.float64)
supercell=CELL.cell2supercell(unitcell,P)
supercell.print_poscar("./super.vasp")
```

## Install

* Copy build.py to working directory.

* Or more elegantly, you can install SlabMaker as a python library instead of copy build.py to every working directory. Copy the SlabMaker folder to $HOME/local/lib or wherever you want, add the directory to the $PYTHONPATH by writing to the file ~/.bashrc:

```bash
export PYTHONPATH=$HOME/local/lib:$PYTHONPATH
```

reopen the terminal, or ```source ~/.bashrc```

then SlabMaker can be imported as:
```
from SlabMaker.build import CELL
```

run the test.py in examples directory,

```bash
python test.py
```

## Dependences:

* Numpy >= 1.15

# An introduction to the principles behind.
 
This code provide a single-python-file solution to build slab model for atomistic calculations of crystals.

A summary to how the things get done.

1. Find a cell from the definition of miller index. This cell is not necessarily the in-plane primitive cell.

2. Using a brute-force method to find the smallest unit in-plane. (1) For each atom in the cell, find in-plane atoms in the cell and the neighbor cells in certain range. Then we get a list of vectors in-plane. Find the pairs of in-plane vectors with several different areas. (2) Do the step (1) for all the N atoms in the cell. So we get N groups of list of pairs of in-plane vectors, each group consists of several pairs of vectors, find common ones in these pairs, i.e. common in-plane translation invariant vector pairs. At least one common pair is found, or maybe more than one pair is found. (3) Find the smallest and with angle closest to 90 degrees( and other preferential conditions) in the common pairs of vectors.

A more detailed demo (in Chinese):

https://yyyu200.github.io/DFTbook/blogs/2019/04/07/TransCell/

