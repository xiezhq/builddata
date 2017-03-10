#ifndef DIFFINTERF_H
#define DIFFINTERF_H

#include "max.h"
#include "bool.h"
#include "asa2pdb.h"


/* ------------------------------------ */

#define INTERF_MAX	10
/* maximum of number of interfaces in a protein complex */


typedef struct chainid
{
	char	id;

	struct chainid	*next;
} CHAINID;


typedef struct
{
	double	surfCutoff; /* ratio of ASA/stdAsa, used to define surface and buried residues */

	double	interfCutoff;
		/* Residues are defined as interface residues when the monomer loses at least
		interfCutoff angstroms**2 of the ASA upon association. */

	int	nchain1, nchain2; /* number of chains in two interacting protein units */

	char	atm2hybridFile[FILE_NAME_LEN];
		/* mapping atom name to hybrid name */
	char	hybrid2radFile[FILE_NAME_LEN];
		/* mapping the hybrid of molecular orbit to radii of atoms */

	CHAINID	*chainId1, *chainId2; /* identifiers of two interacting protein units */
} PAR;


typedef struct
{
	char	name; /* single-letter name of residue */
	char	chainid;
	int	seq;
	char	iCode;

	int	pol; /* pol = 1 or -1 if polar or nonpolar residue, respectively, else pol = 0 */
	int	philic; /* 1 and -1 for hydrophilic and hydrophobic residues, respectively, else 0 */
	int	e; /* 1 and -1 for positively and negatively charged residues, respectively, else 0 */

	double	asa; /* solvent accessible surface area */
	double	acc; /* relative accessibility */

	int	intflg;
		/* intflg == TRUE if the residue is located at interface, else intflg == FALSE */
} SURF;

#endif
