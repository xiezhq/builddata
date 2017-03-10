#ifndef CNCTAREA2_H
#define	CNCTAREA2_H


extern void cnctarea(double errorpar, double *probe, double *x, double *y, double *z,
		      double *radx, int natms, double *accss);

/*-------------------
The comments below are added by Zhiqun Xie in 2005.
---------------------
errorpar	z-spacing factor (slice thickness), default errorpar = 0.05 angstrom
probe		probe = atomic radius + radius of solvent molecule, default = Rvdw + Rwater
x, y, z		coordinates for each atom
radx		van der Walls radius of atom
natms		number of atoms in the given protein unit
accss		accessible area for each atom

---------------------
Note: Lee & Rechards algorithm, refer to:

Lee B & Richards FM. The interpretation of protein structures: estimation of static accessibility.
J Mol Biol, 1971, 55:379-400.

readme.txt of NACCESS.
Hubbard SJ & Thornton JM. NACCESS. Department of Biochemistry and Molecular Biology,
University College London, 1993.

J Biol Chem, 2004, 279(22):23061-23072 for the explanation of errorpar.

-------------------
The comments above are added by Zhiqun Xie in 2005.
-----------------*/


extern void sortag(double *a, int n, int *tag);
extern void *ckalloc(int amount);
#include <stdio.h>
extern FILE *ckopen(char *name, char *mode);

#endif
