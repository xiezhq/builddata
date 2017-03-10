#ifndef PDB_H
#define PDB_H


#undef ATOM
#define ATOM
/*
Just in case, discard pdbATOM() function in pdb.c file
*/


#define BREAKDIST	2.5
#define BREAKDIST2	6.25
/*
BREAKDIST	Maximum allowed peptide bond length. If distance is greater, a polypeptide chain
		interruption is assumed.
BREAKDIST2	BREAKDIST * BREAKDIST
*/


/* ------------------------ */

typedef struct atom
{
	int	serial; /* atom serial number */
	char	chem[3]; /* chemical symbol, except for hydrogen atoms */
	char	name[5]; /* atom name */
	char	altLoc; /* alternate location indicator */

	double	x; /* orthogonal coordinates for X in angstroms */
	double	y;
	double	z;

	double	r; /* radius of van der Walls */
	double	accss; /* solvent accessible surface area */

	struct atom	*next; /* next atom in same residue */
} ATM;


typedef struct residue
{
	char	name; /* single-letter name of residue */
	int	seq; /* residue sequence number */
	char	iCode; /* code for insertion residues */
	int	natm; /* number of atoms in residue */

	int	pol; /* 1 for polar residue, -1 for nonpolar residue, 0 for others */
	int	philic; /* 1 for hydrophilic residue, -1 for hydrophobic residue, 0 for special residue */
	int	e; /* 1 for positively charged residue, -1 for negatively charged residue, 0 for others */

	double	asa; /* solvent accessible surface area */
	/* double	stdacc;*/ /* solvent accessible surface area in tripeptide */
	double	acc; /* relative solvent accessibility, 0.0 ~ 1.0 */

	double	asaSide; /* solvent accessible surface area of the sidechain */
	double	accSide;

	ATM	*atm; /* the first atom of the residue */

	struct residue	*next; /* next residue in same chain */
} RES;


typedef struct
{
	char	het; /* het='H' if HETATM chain (non-ATOM chain), else het='A' */
	char	id; /* chain identifier */
	int	nres; /* number of residues in chain */

	RES	*res;
} CHAIN;


typedef struct
{
	char	name; /* single-letter name of residue */
	double	asa; /* standard ASA of residue */
} STDASA;


/* ------------------------------------------------------- */


/* Prototype of functions */

/* translate three-letter name into single-letter name and return it. The unkonwn name
will be translated into 'X'. */
extern char	t2sRes(char *res);

/* translate single-letter name into three-letter name and return a pointer pointer to
an 4-char storage holding three-letter name. The unknown name will be translated into
'UNK'. */
extern char	*s2tRes(char sRes);


extern void	freeChain(CHAIN *chain, int nchain);


/* process ATOM records in PDB file and return the number of chains */
extern int	pdbATOM(char *file, CHAIN *chain);

/*
pdbHETATOM() processes HETATM and ATOM records in PDB file and return the number of chains.

Note: it returns number of chains in ATOM + HETATM records, namely,
if there are two chains with same identifier but in different records (HETATM and ATOM records),
the number of chains will be 2 not 1 although they have the same chain identifier.
*/
extern int	pdbHETATOM(char *file, CHAIN *chain);

/*
rmhetchain() removes chains appearing in HETATM records and returns number of left chains.
*/
extern int	rmhetchain(CHAIN *chain, int nchain);

/*
Remove chains containing nucleic acid residue and return number of left chains
*/
extern int	rmnuchain(CHAIN *chain, int nchain);

/*
Remove chains containing water molecules and return number of left chains
*/
extern int	rmWater(CHAIN *chain, int nchain);

/*
remove residues specified by "name" (single-letter name) from PDB structure pointed by "CHAIN *chain",
return number of removed residues
*/
extern int	rmres(CHAIN *chain, int nchain, char name);


#include "asa2pdb.h"
extern int	checkATOM(CHAIN *chain, int nchain, RESRAD *rad, int nresrad);

/* upper case to lower case letter conversion for character string:
*str	character string to be converse
n	the first n characters will be converse if n > 0, else converse all characters in *str
*/
void	str2lower(char *str, int n);

/* search for atoms with same coordinates on either X, Y or Z dimension */
extern void	coordIdentity(CHAIN *chain, int nchain);

#endif
