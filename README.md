# SlabMaker

Consist of two part:
1. TransCell
Transformation between primitive cell and unit cell.

2. SlabMaker
Build slab models for crystal surfaces.

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

Dependences:

* Numpy >= 1.15

# An introduction of the principles behind.

https://yyyu200.github.io/DFTbook/blogs/2019/04/07/TransCell/
