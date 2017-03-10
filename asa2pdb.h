#ifndef ASA2PDB_H
#define	ASA2PDB_H


#include "pdb.h"


/* ---------------------------------- */

#define ERRORPAR	0.05
/* z-spacing factor */

#define Rwat	1.4
/*
#define Rwat	1.46
*/
/* radius of water molecule. Most people use 1.4 as the vdw radius of water molecule.
Since we use the radii set from Gerstein's data set in calculation, it seems better to
define radius of water as 1.46(O2H2u) in our computation. HOWEVER, the calculated ASAs
of residues in extended tripeptide gly-X-gly are more close to the results in literatures
especially the results reported by Protein Engineering, 2002, 15(8):659-667, when we
use 1.4 in the calculation. For comparison, also refer to Protein Science, 2003, 12:1406-1417. */


/* ---------------------------------- */

#define HYBRID_MAX	50

typedef struct atomhybrid
{
	char	hybrid[10];
		/*
		name of atomic group defined by status of hybrid of molecular orbit.

		Usually we use the definition from Gerstein's report (JMB, 1999, 290:253-266).

		Nomenclature:
		The atomic groups found in proteins are given labels of the general form XnHmx, where
		X indicates the chemical nature of the non-hydrogen atoms; n, their valence; Hm, the
		number (m) of hydrogen atoms attached to the non-hydrogen atom; and x, the atom type
		assigned according to their volumes: s, b, or u (small, big, or unique, respectively).
		*/

	double	r;
		/* VDW radius of atom or atomic group.

		Usually we use the value from Gerstein's report (JMB, 1999, 290:253-266).
		In Gerstein's definition, atomic group subsume a heavy-atom and its covalently attached 
		hydrogen atoms into one moiety because the positions of hydrogen atoms in protein structures
		are generally not known.
		*/
} HYBRID;



#define	N_RESRAD	50
/* number of types of amino acid residues, nucleic acid residues and het groups 
probably existing in PDB database */

#define RAD_DEFAULT	1.8
/* VDW radius for unknown atom */

typedef struct atmrad
{
	char	name[5]; /* name of atom defined by PDB */
	int	namei; /* index of name of atom above */
	char	hybrid[10]; /* name of atom defined by status of hybrid of molecular orbit */
	double	r; /* VDW radius for atomic group */

	struct atmrad	*next;
} ATMRAD;

typedef struct
{
	char	namet[5]; /* multiple-letter name of residue */
	char	name; /* single-letter name of residue */
	int	namei; /* index of name of residue or het group */
	int	natm; /* number of atoms in residue */

	ATMRAD	*atm; /* the first atom of residue */
} RESRAD;


/* ---------------------------------- */

extern void	pdbArray(CHAIN *chain, int nchain, CHAIN *pdb, RES *residue, ATM *atom);


/* get hybrid status of atoms and the corresponding atomic radii

HYBRID	*hybrid:	all types of hybrids (the last element of hybrid array is the default value) 
return value:		number of types of hybrids of atoms
*/
extern int	hybrid2rad(char *file, HYBRID *hybrid);


/* read in from file the parameters for atomic VDW radii defined by Gerstein */
extern int	atm2hybrid(char *file, RESRAD *radi);

/* read in parameters for atomic VDW radii defined by Gerstein */
extern void	getRad(HYBRID *hybrid, int nhybrid, RESRAD *rad, int nresrad);

/* add radii of atoms to structure */
extern void	rad2pdb(CHAIN *chain, int nchain, RESRAD *radi, int nrad);

extern void	rad2pdbVer2(CHAIN *pdb, int nchain, RESRAD *rad, int nresrad);

extern void	freeRad(RESRAD *rad, int n);

/* count the number of elements in protein unit, such as residues and atoms.

chain0	full details for each atom, residue and chain in computed protein unit
nchain	number of chains
nres	number of residues

return: number of atoms
*/
extern int	nele(CHAIN *chain, int nchain, int *nres);


/* calculate the solvent accessible area for each atom in computed protein unit

chain0	full details for each atom, residue and chain in computed protein unit
nchain	number of chains in processed protein
*/
extern void	asa2pdb(CHAIN *chain, int nchain, int natm);
extern void	asa2pdbVer2(ATM *atm, int natm, double *accss);

/* All atoms are defined as sidechain atom except four atoms, N, CA, C, O. 
The only exception is GLY since it has no normal sidechain. So the CA atom
is defined as the sidechain of GLY for convenience of further calculation.
*/
extern void	asa4sidechain(CHAIN *chain, int nchain);


/* get data of standard solvent accessible surface area for normal residues */
extern int	rdStdAsa(char *file, STDASA *stdAsa);

/* calculate the relative solvent accessibility of residues in protein */
extern void	acc(CHAIN *chain, int nchain, STDASA *stdAsa, int nstdAsa);

/* calculate the relative solvent accessibility of N-terminal residue */
extern void	accN(CHAIN *chain, int nchain, STDASA *stdAsa, int nstdAsa);

/* calculate the relative solvent accessibility of C-terminal residue */
extern void	accC(CHAIN *chain, int nchain, STDASA *stdAsa, int nstdAsa);

#endif
