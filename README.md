# SlabMaker

1. TransCell
Transformation between primitive cell and unit cell.

2. SlabMaker
Build slab models for crystal surfaces.

3. SuperCell
Build supercell given the transform matrix P.

4. PW input file convert to POSCAR format.

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
from CELL import CELL

sc=CELL("sc.vasp")
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
from CELL import CELL

unitcell=CELL("sc.vasp")

miller_index=[1,1,1]
slab=unitcell.makeslab(miller_index, layer=1, vacuum=15.0)

slab.print_poscar("./slab.vasp")
```

will get 
```
P1 =  [[ 1.  0.  1.]
 [ 0.  1.  1.]
 [-1. -1.  1.]]
inplane vectors:
u:
 [ 0.  -0.5  0. ]
v:
 [-0.5  0.   0. ]
angle:
 59.99999999999999
P2 =  [[ 0.  -0.5  0. ]
 [-0.5  0.   0. ]
 [ 0.   0.  -1. ]]
reduced slab cell
 [[ -1.53725016  -2.66259538   0.        ]
 [ -3.07450032   0.           0.        ]
 [  0.           0.         -21.90337725]]
reduced slab cell area:  8.186150349361137  Ang^2
```
and also the slab.vasp file:

```
system: Si C
1.0
   -1.5372501598   -2.6625953808    0.0000000000
   -3.0745003197    0.0000000000    0.0000000000
    0.0000000000    0.0000000000  -21.9033772475
Si C
3 3
Direct
 0.000000000000 0.000000000000 0.657587050834
 0.666666666667 0.666666666667 0.542978286591
 0.333333333333 0.333333333333 0.428369522348
 0.000000000000 0.000000000000 0.571630477652
 0.666666666667 0.666666666667 0.457021713409
 0.333333333333 0.333333333333 0.342412949166
```

## Build supercell:

```python
import numpy as np
from CELL import CELL

unitcell=CELL("sc.vasp")
P=np.mat([[2,0,0],[0,2,0],[0,0,2]],dtype=np.float64)
supercell=CELL.cell2supercell(unitcell,P)
supercell.print_poscar("./super.vasp")
```

## Dependences:

* Numpy >= 1.15

# An introduction of the principles behind.

https://yyyu200.github.io/DFTbook/blogs/2019/04/07/TransCell/
