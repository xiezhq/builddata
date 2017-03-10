#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <limits.h>
#include "bool.h"
#include "pdb.h"
#include "memory.h"
#include "max.h"
#include "asa2pdb.h"
#include "property4res.h"



/*
refer to PDB format description version 2.2 and Appendix5(formulas and molecular weights
for standard residues).

NOTE: Here are amino acid residues not amino acids.

Be careful to deal with the nucleic acid residues */


/*
#define N_STDRES	24
#define N_STDRES	36
#define N_STDRES	42
*/
/* number of standard residues + an unknown residue */



/*
static struct stdresidue
{
	char	*three;
	char	single;
	int	natm;
	int	pol;
	int	philic;
	int	e;
} stdres[] =
{
	{"ALA", 'A', 5,  -1, -1,  0},
	{"ARG", 'R', 11,  1,  1,  1},
	{"ASN", 'N', 8,   1,  1,  0},
	{"ASP", 'D', 8,   1,  1, -1},
	{"CYS", 'C', 6,   1,  0,  0},

	{"GLN", 'Q', 9,   1,  1,  0},
	{"GLU", 'E', 9,   1,  1, -1},
	{"GLY", 'G', 4,  -1,  0,  0},
	{"HIS", 'H', 10,  1,  1,  1},
	{"ILE", 'I', 8,  -1, -1,  0},

	{"LEU", 'L', 8,   1, -1,  0},
	{"LYS", 'K', 9,   1,  1,  1},
	{"MET", 'M', 8,  -1, -1,  0},
	{"PHE", 'F', 11, -1, -1,  0},
	{"PRO", 'P', 7,  -1,  0,  0},

	{"SER", 'S', 7,   1,  1,  0},
	{"THR", 'T', 7,   1,  1,  0},
	{"TRP", 'W', 14, -1, -1,  0},
	{"TYR", 'Y', 12,  1, -1,  0},
	{"VAL", 'V', 7,  -1, -1,  0},

	{"ASX", 'B', 8,   1,  1,  0},
	{"GLX", 'Z', 9,   1,  1,  0},
	{"CSS", 'c', 6,   0,  0,  0},

	{"ZN",  'a', 1,   1,  1,  2},
	{"HEM", 'b', 39,  0,  0,  0},
	{"ACE", 'd', 3,   0,  0,  0},
	{"HOH", 'e', 1,   1,  1,  0},
	{"WAT", 'e', 1,   1,  1,  0},
	{"DOD", 'f', 1,   1,  1,  0},

	{"A",   '1', 0,   0,  0,  0},
	{"+A",  '1', 0,   0,  0,  0},
	{"C",   '2', 0,   0,  0,  0},
	{"+C",  '2', 0,   0,  0,  0},
	{"G",   '3', 0,   0,  0,  0},
	{"+G",  '3', 0,   0,  0,  0},
	{"I",   '4', 0,   0,  0,  0},
	{"+I",  '4', 0,   0,  0,  0},
	{"T",   '5', 0,   0,  0,  0},
	{"+T",  '5', 0,   0,  0,  0},
	{"U",   '6', 0,   0,  0,  0},
	{"+U",  '6', 0,   0,  0,  0},
	{"UNK", 'X', 8,   0,  0,  0}
};



static struct stdNucleicAcid
{
	char	*three;
	char	single;
	int	natm;
} stdna[] =
{
	"  A", 'A', 0,
	" +A", 'A', 0,
	"  C", 'C', 0,
	" +C", 'C', 0,
	"  G", 'G', 0,
	" +G", 'G', 0,
	"  I", 'I', 0,
	" +I", 'I', 0,
	"  T", 'T', 0,
	" +T", 'T', 0,
	"  U", 'U', 0,
	" +U", 'U', 0,
	"UNK", 'X', 0
};
*/


/* upper case to lower case letter conversion for character string:
char	*str:	character string to be converse
int	n:	the first n characters will be converse if n > 0, else converse all characters in *str
*/
void	str2lower(char *str, int n)
{
	char	*strt;

	
	if(n > 0)
	{
		for(strt = str + n; str < strt; str ++)
		{
			*str = tolower(*str);
		}
	}
	else
	{
		while((*str++ = tolower(*str)) != '\0');
	}
}


char	t2sRes(char *res)
{
	extern struct stdresidue	stdres[];
	int	i;

	for(i = 0; stdres[i].single != 'X'; i ++)
	{
		if(strcmp(stdres[i].three, res) == 0)
			return stdres[i].single;
	}

	/* assign X as the single-letter name for any unknown residue */
	return 'X';
}


char	*s2tRes(char sRes)
{
	extern struct stdresidue	stdres[];
	int	i;
	char	*tRes;


	tRes = (char *)memalloc(4 * sizeof(char));

	for(i = 0; stdres[i].single != 'X'; i ++)
	{
		if(stdres[i].single == sRes)
		{
			strcpy(tRes, stdres[i].three);

			return tRes;
		}
	}

	/* assign UNK as the three-letter name for any unknown residue */
	tRes[0] = 'U';
	tRes[1] = 'N';
	tRes[2] = 'K';
	tRes[3] = '\0';

	return tRes;
}


static void	freeRes(RES *res)
{
	ATM	*atm, *nextatm;


	if(res != NULL)
	{
		for(atm = res->atm; atm != NULL; atm = nextatm)
		{
			nextatm = atm->next;
			free((void *)atm);
		}

		free((void *)res);
	}
}


void	freeChain(CHAIN *chain, int nchain)
{
	CHAIN	*chaint;
	RES	*res, *nextres;
	/*
	ATM	*atm, *nextatm;
	*/


	for(chaint = chain + nchain; chain < chaint && chain != NULL; chain ++)
	{
		for(res = chain->res; res != NULL; res = nextres)
		{
			/*
			for(atm = res->atm; atm != NULL; atm = nextatm)
			{
				nextatm = atm->next;
				free((void *)atm);
			}

			nextres = res->next;
			free((void *)res);
			*/

			nextres = res->next;
			freeRes(res);
		}
	}
}


/* calculate the distances between two groups of atoms
*/
/*
void	dist(double *x1, double *y1, double *z1, int n1, double *x2, double *y2, double *z2, int n2,
		double *dd)
{
	int	i, j, n;
	double	dx, dy, dz;


	n = 0;
	for(i = 0; i < n1; i ++)
	{
		for(j = 0; j < n2; j ++)
		{
			dx = x1 - x2;
			dy = y1 - y2;
			dz = z1 - z2;

			dd[n] = dx * dx + dy * dy + dz * dz;

			n ++;
		}
	}
}
*/


/* Check the atoms and residues in a structure

Activities: 
. use ASP and GLU instead of ASX and GLX, respectively, if ASX and GLX exist in PDB file
. check atoms of residue against the corresponding normal residue

Return value:
. number of chains left in PDB structure
*/
int	checkATOM(CHAIN *chain, int nchain, RESRAD *rad, int nresrad)
{
	CHAIN	*chaint;
	RES	*res;
	RES	*lastres;
	ATM	*atm;
	ATM	*lastatm;
	RESRAD	*radi, *radt;
	ATMRAD	*atmr;
	RES	*nextres;
	ATM	*atm2;
	double	dx, dy, dz;
	/*
	int	nucleicflg;
	*/
	int	i;


	radi = rad;
	radt = radi + nresrad;
	/*
	nucleicflg = FALSE;
	*/
	for(chaint = chain + nchain; chain < chaint && chain != NULL; chain ++)
	{
		lastres = NULL;
		for(res = chain->res; res != NULL; res = nextres)
		{
			/*
			 . check nucleic acid residues
			 . Note: single-letter name for nucleic acid residues are '1', '2', '3',
			   '4', '5' and '6'.
			*/
			/*
			switch(res->name)
			{
				case '1':
				case '2':
				case '3':
				case '4':
				case '5':
				case '6':
				{
					nucleicflg = TRUE;
					break;
				}
			}
			if(nucleicflg == TRUE)
				break;
			*/

			nextres = res->next;

			/* check C-terminal residue */
			/*
			if(res->next == NULL)
			{
				lastatm = NULL;
				for(atm = res->atm; atm != NULL; atm = atm->next)
				{
					if(memcmp("OXT", atm->name, 4) == 0)
					{
						if(lastatm == NULL)
							res->atm = atm->next;
						else	lastatm->next = atm->next;

						res->natm -= 1;

						free((void *)atm);
						printf("#Remark: remove OXT atom (%c)\n", chain->id);
						break;
					}
					lastatm = atm;
				}
			}
			*/

			/* . check the ambiguity of the terminal atoms of ASN and GLN and the ring
			     atoms of HIS
			   . set residue as ASP if ASX (ASP is more redundant than ASN in proteins)
			*/
			if(res->name == 'B')
			{
				for(atm = res->atm; atm != NULL; atm = atm->next)
				{
					if(*atm->name == 'A'
					&& *(atm->name + 1) == 'D'
					&& *(atm->name + 3) == '\0'
					&& (*(atm->name + 2) == '1' || *(atm->name + 2) == '2'))
					{
						*atm->name = 'O';

						printf("#Remark: set ASX as ASP (%c %d%c)\n",
							chain->id, res->seq, res->iCode);
					}
				}
			}
			/* set residue as GLU if GLX (GLU is more redundant than GLN in proteins) */
			else if(res->name == 'Z')
			{
				for(atm = res->atm; atm != NULL; atm = atm->next)
				{
					if(*atm->name == 'A'
					&& *(atm->name + 1)== 'E'
					&& *(atm->name + 3) == '\0'
					&& (*(atm->name + 2) == '1' || *(atm->name + 2) == '2'))
					{
						*atm->name = 'O';

						printf("#Remark: set GLX as GLU (%c %d%c)\n",
							chain->id, res->seq, res->iCode);
					}
				}
			}
			/* set first atom as oxygen and second as nitrogen if AD1 and AD2 exist in ASN */
			else if(res->name == 'N')
			{
				for(atm = res->atm; atm != NULL; atm = atm->next)
				{
					if(*atm->name == 'A'
					&& *(atm->name + 1)== 'D'
					&& *(atm->name + 3) == '\0')
					{
						if(*(atm->name + 2) == '1')
							*atm->name = 'O';
						else if(*(atm->name + 2) == '2')
							*atm->name = 'N';

						printf("#Remark: fixed ambiguous terminal atoms (%c %d%c)\n",
							chain->id, res->seq, res->iCode);
					}
				}
			}
			/* set first atom as oxygen and second as nitrogen if AE1 and AE2 exist in GLN */
			else if(res->name == 'Q')
			{
				for(atm = res->atm; atm != NULL; atm = atm->next)
				{
					if(*atm->name == 'A'
					&& *(atm->name + 1) == 'E'
					&& *(atm->name + 3) == '\0')
					{
						if(*(atm->name + 2) == '1')
							*atm->name = 'O';
						else if(*(atm->name + 2) == '2')
							*atm->name = 'N';

						printf("#Remark: fixed ambiguous terminal atoms (%c %d%c)\n",
							chain->id, res->seq, res->iCode);
					}
				}
			}
			/* set the four terminal atoms as ND1, CD2, CE1 and NE2 if AD1, AD2, AE1 and AE2
			   exist in HIS
			*/
			else if(res->name == 'H')
			{
				for(atm = res->atm; atm != NULL; atm = atm->next)
				{
					if(*atm->name != 'A' || *(atm->name + 3) != '\0')
						continue;

					if(*(atm->name + 1) == 'D')
					{
						if(*(atm->name + 2) == '1')
							*atm->name = 'N';
						else if(*(atm->name + 2) == '2')
							*atm->name = 'C';

						printf("#Remark: fixed ambiguous ring atoms of HIS (%c %d%c)\n",
							chain->id, res->seq, res->iCode);
					}
					else if(*(atm->name + 1) == 'E')
					{
						if(*(atm->name + 2) == '1')
							*atm->name = 'C';
						else if(*(atm->name + 2) == '2')
							*atm->name = 'N';

						printf("#Remark: fixed ambiguous ring atoms of HIS (%c %d%c)\n",
							chain->id, res->seq, res->iCode);
					}
				}
			}

			/*
			check atoms in each residue against the corresponding residue from
			parameter file (ProtOr.resi-defs.dat or fmr74.resi-defs.dat)
			*/
			for(rad = radi; rad < radt && rad != NULL; rad ++)
			{
				if(res->name != rad->name || res->name == 'X')
					continue;

				/*
				if(res->natm < rad->natm)
				{
					printf("#Warning: missing atoms (%c %d%c, %d < %d)\n",
						chain->id, res->seq, res->iCode, res->natm, rad->natm);
				}
				else if(res->natm > rad->natm)
				{
					printf("#Warning: found additional atoms (%c %d%c, %d > %d)\n",
						chain->id, res->seq, res->iCode, res->natm, rad->natm);
				}
				*/

				for(atm = res->atm; atm != NULL; atm = atm->next)
				{
					for(atmr = rad->atm; atmr != NULL; atmr = atmr->next)
					{
						if(strcmp(atm->name, atmr->name) == 0)
							break;
					}

					if(atmr == NULL)
					{
						printf("#Warning: non-standard atom (%s) in residue (%c %d%c)\n",
							atm->name, chain->id, res->seq, res->iCode);
					}
				}

				i = 0;
				for(atmr = rad->atm; atmr != NULL; atmr = atmr->next)
				{
					for(atm = res->atm; atm != NULL; atm = atm->next)
					{
						if(strcmp(atmr->name, atm->name) == 0)
						{
							i ++;
							break;
						}
					}
				}
				if(i < rad->natm)
				{
					/* exit if any residue has missing atoms */
					/*
					printf("#Error: %d missing atoms (%c %d%c %c)\n",
						rad->natm - i, chain->id, res->seq, res->iCode, res->name);

					exit(1);
					*/
					/* exit if any residue has missing atoms */

					if(i == res->natm)
					{
						printf("#Warning: %d missing atoms (%c %d%c %c)\n",
							rad->natm - i, chain->id, res->seq, res->iCode, res->name);
					}
					else if(i < res->natm)
					{
						printf("#Warning: %d missing atoms but with %d non-standard atoms (%c %d%c %c)\n",
							rad->natm - i, res->natm - i, chain->id, res->seq, res->iCode, res->name);
					}
				}
				else if(i == rad->natm && res->natm > rad->natm)
				{
					printf("#Warning: %d additional non-standard atoms (%c %d%c)\n",
						res->natm - i, chain->id, res->seq, res->iCode);
				}

				break;
			}
			if(rad == radt)
			{
				printf("#Warning: unknown residue (%c %d%c %c)\n",
					chain->id, res->seq, res->iCode, res->name);
			}

			/* check the chain break */
			if(nextres == NULL)
				break;

			for(atm = res->atm; atm != NULL; atm = atm->next)
			{
				if(*atm->name != 'C' || *(atm->name + 1) != '\0')
					continue;

				for(atm2 = nextres->atm; atm2 != NULL; atm2 = atm2->next)
				{
					if(*atm2->name != 'N' || *(atm2->name + 1) != '\0')
						continue;

					dx = atm->x - atm2->x;
					dy = atm->y - atm2->y;
					dz = atm->z - atm2->z;
					if((dx*dx + dy*dy + dz*dz) > BREAKDIST2)
					{
						printf("#Warning: chain (%c) break between %d%c and %d%c\n",
							chain->id, res->seq, res->iCode,
							nextres->seq, nextres->iCode);
					}

					break;
				}

				break;
			}

			lastres = res;
		}
		/* Remove chain containing at least one nucleic acid residue
		*/
		/*
		if(nucleicflg == TRUE)
		{
			printf("#Warning: remove nucleic acid chain (%c)\n", chain->id);

			freeChain(chain, 1);
			nchain --;
			chaint --;

			if(chaint == chain)
			{
				break;
			}

			chain->het = chaint->het;
			chain->id = chaint->id;
			chain->nres = chaint->nres;
			chain->res = chaint->res;

			nucleicflg = FALSE;
			chain --;
		}
		*/
	}

	return nchain;
}


/*
rmhetchain() removes chains appearing in HETATM records and returns number of left chains.
*/
int	rmhetchain(CHAIN *chain, int nchain)
{
	CHAIN	*chaint;


	for(chaint = chain + nchain; chain < chaint && chain != NULL; chain ++)
	{
		if(chain->het == 'H')
		{
			printf("#Warning: remove HETATM chain (chainId=%c)\n", chain->id);

			freeChain(chain, 1);
			nchain --;
			chaint --;

			if(chaint == chain)
			{
				break;
			}

			chain->het = chaint->het;
			chain->id = chaint->id;
			chain->nres = chaint->nres;
			chain->res = chaint->res;

			chain --;
		}
	}

	return nchain;
}


/* rmWater() removes chain(s) containning water molecules and returns number of left chains.
*/
int	rmWaterChain(CHAIN *chain, int nchain)
{
	CHAIN	*chaint;
	RES	*res;
	int	watFlg;


	watFlg = FALSE;
	for(chaint = chain + nchain; chain < chaint && chain != NULL; chain ++)
	{
		for(res = chain->res; res != NULL; res = res->next)
		{
			/* check water molecules
			Note: single-letter name for water are 'e' and 'f'.
			*/
			switch(res->name)
			{
				case 'e':
				case 'f':
				{
					watFlg = TRUE;
					break;
				}
			}
			if(watFlg == TRUE)
				break;
		}

		if(watFlg == TRUE)
		{
			printf("#Warning: remove water molecules (chainId=%c)\n", chain->id);

			freeChain(chain, 1);
			nchain --;
			chaint --;

			if(chaint == chain)
			{
				break;
			}

			chain->het = chaint->het;
			chain->id = chaint->id;
			chain->nres = chaint->nres;
			chain->res = chaint->res;

			watFlg = FALSE;
			chain --;
		}
	}

	return nchain;
}


/*
Remove chains containing nucleic acid residue

Return value: number of left chains
*/
int	rmnuchain(CHAIN *chain, int nchain)
{
	CHAIN	*chaint;
	RES	*res;
	int	nucleicflg;


	nucleicflg = FALSE;
	for(chaint = chain + nchain; chain < chaint && chain != NULL; chain ++)
	{
		for(res = chain->res; res != NULL; res = res->next)
		{
			/* check nucleic acid residues
			Note: single-letter name for nucleic acid residues are '1', '2', '3',
			'4', '5' and '6'.
			*/
			switch(res->name)
			{
				case '1':
				case '2':
				case '3':
				case '4':
				case '5':
				case '6':
				{
					nucleicflg = TRUE;
					break;
				}
			}
			if(nucleicflg == TRUE)
				break;
		}

		if(nucleicflg == TRUE)
		{
			printf("#Warning: remove nucleic acid chain (chainId=%c)\n", chain->id);

			freeChain(chain, 1);
			nchain --;
			chaint --;

			if(chaint == chain)
			{
				break;
			}

			chain->het = chaint->het;
			chain->id = chaint->id;
			chain->nres = chaint->nres;
			chain->res = chaint->res;

			nucleicflg = FALSE;
			chain --;
		}
	}

	return nchain;
}


/*
remove residues specified by "name" (single-letter name) from PDB structure pointed by "CHAIN *chain",
return number of removed residues
*/
int	rmres(CHAIN *chain, int nchain, char name)
{
	CHAIN	*chaint;
	RES	*res, *lastres, *nextres;
	int	nrm;


	nrm = 0;
	for(chaint = chain + nchain; chain < chaint && chain != NULL; chain ++)
	{
		lastres = NULL;
		for(res = chain->res; res != NULL; res = nextres)
		{
			nextres = res->next;
			if(res->name == name)
			{
				printf("#Warning: remove residue %c (%c %d%c)\n",
					name, chain->id, res->seq, res->iCode);

				if(res == chain->res)
					chain->res = nextres;
				else	lastres->next = nextres;

				chain->nres --;
				nrm ++;
				freeRes(res);
			}
			else	lastres = res;
		}
	}

	return nrm;
}


void coordIdentity(CHAIN *chain0, int nchain0)
{
CHAIN	*chaint, *chain;
CHAIN	*chaint1, *chain1;
RES	*res, *res1;
ATM	*atm, *atm1;

printf("Start searching for coordinate identity\n");
chain = chain0;
for(chaint = chain0 + nchain0; chain < chaint && chain != NULL; chain ++)
{
	for(res = chain->res; res != NULL; res = res->next)
	{
		for(atm = res->atm; atm != NULL; atm = atm->next)
		{
			chain1 = chain0;
			for(chaint1 = chain0 + nchain0; chain1 < chaint1; chain1 ++)
			{
				for(res1 = chain1->res; res1 != NULL; res1 = res1->next)
				{
					for(atm1 = res1->atm; atm1 != NULL; atm1 = atm1->next)
					{
						if(chain == chain1 && res == res1 && atm == atm1)
							continue;

						if(!(atm->x > atm1->x || atm->x < atm1->x))
						{
							printf("X-coord identity: "
								"%6d %6d %8.3f %8.3f\n",
								atm->serial, atm1->serial,
								atm->x, atm1->x);
						}
						if(!(atm->y > atm1->y || atm->y < atm1->y))
						{
							printf("Y-coord identity: "
								"%6d %6d %8.3f %8.3f\n",
								atm->serial, atm1->serial,
								atm->y, atm1->y);
						}
						if(!(atm->z > atm1->z || atm->z < atm1->z))
						{
							printf("Z-coord identity: "
								"%6d %6d %8.3f %8.3f\n",
								atm->serial, atm1->serial,
								atm->z, atm1->z);
						}
					}
				}
			}
		}
	}
}

printf("coordIdentity searching is finished!\n");
}


#ifdef ATOM
/*
Return number of chains if no error was found with it when parsing ATOM records in PDB file, else return zero. 
Process only data holded in ATOM record line and remove hydrogen atoms if any.
*/

#define LEN_ATOM	82
/* number of columns in ATOM record (80) + '\n' + '0', 82 = 80 + 1 + 1) */

int	pdbATOM(char *file, CHAIN *chain0)
{
	FILE	*fp;
	int	zipflg;
	char	*chp;
	char	str[LEN_ATOM];
	char	strTmp[10];
	int	intTmp;
	CHAIN	*chain;
	char	lastchain;
	int	lastSeq;
	char	lastiCode;
	char	lastaltLoc;
	RES	*res, *lastres;
	ATM	*atm, *lastatm;
	

	if((chp = strrchr(file, '.')) == NULL)
	{
		printf("#Unkown type of file: %s\n", file);
		return 0;
	}

	zipflg = FALSE;
	if(chp[1] == 'Z' || (chp[1] == 'g' && chp[2] == 'z'))
	{
		char	cmd[NAME_MAX] = "zcat ";
		int	len;

		zipflg = TRUE;
		len = strlen(file);
		memcpy(cmd+5, file, len+1);

		if((fp = popen(cmd, "r")) == NULL)
		{
			printf("#Open %s error!\n", file);
			return 0;
		}
	}
	else if((fp = fopen(file, "r")) == NULL)
	{
		printf("#Open %s error!\n", file);
		return 0;
	}

	/*
	for(; fgets(str, LEN_ATOM, fp) != NULL;)
	{
		if(memcmp("ATOM  ", str, 6) == 0)
			break;
	}
	*/

	/* initialization */
	chain = chain0 - 1;
	lastchain = '\0';
	lastSeq = INT_MIN; /* minimum value a 'signed int' can hold */
	lastiCode = '\0';
	lastaltLoc = ' ';

	res = NULL;
	lastres = NULL;
	lastatm = NULL;

	for(; fgets(str, LEN_ATOM, fp) != NULL;)
	{
		/* process only the first model in PDB file */
		if(memcmp("ENDMDL", str, 6) == 0)
			break;

		/* process only ATOM record */
		if(memcmp("ATOM  ", str, 6) != 0)
			continue;

		/* remove het-group and hydrogen atoms */
		if(str[12] != ' ' || str[13] == 'H')
			continue;

		if(str[21] != lastchain)
		{
			lastchain = str[21];

			lastSeq = INT_MIN;
			lastiCode = '\0';

			chain ++;
			chain->id = lastchain;
			chain->nres = 0;
			chain->res = NULL;
		}

		memcpy(strTmp, str+22, 4);
		strTmp[4] = '\0';
		if((intTmp = atoi(strTmp)) != lastSeq || str[26] != lastiCode)
		{
			/*
			lastSeq = atoi(strTmp);
			*/
			lastSeq = intTmp;
			lastiCode = str[26];

			chain->nres ++;

			res = (RES *)memalloc(sizeof(RES));

			sscanf(str+17, "%s", strTmp);
			res->name = t2sRes(strTmp);
			res->seq = lastSeq;
			res->iCode = lastiCode;
			res->natm = 0;
			res->atm = NULL;
			res->next = NULL;

			if(chain->res == NULL)
				chain->res = res;
			else	lastres->next = res;

			lastres = res;
		}

		if(str[16] != ' ')
		{
			if(lastaltLoc == ' ')
				lastaltLoc = str[16];
			else if(str[16] != lastaltLoc)
				continue;
		}

		/*
		if(lastaltLoc != ' ' && str[16] != lastaltLoc)
			continue;

		if(str[16] != ' ' && lastaltLoc == ' ')
		{
			lastaltLoc = str[16];
		}
		*/

		res->natm ++;

		atm = (ATM *)memalloc(sizeof(ATM));

		/*
		memcpy(atm->name, str+12, 4);
		*/
		/*
		memcpy(atm->name, str+13, 3);
		*/
		/*
		*(atm->name + 4) = '\0';
		*/
		memcpy(strTmp, str+12, 4);
		strTmp[4] = '\0';
		sscanf(strTmp, "%s", atm->name);

		/* check ambiguity of terminal atoms of ASX/GLX and ring atoms of HIS */
		/*
		if(atm->name == 'A' && *(atm->name + 3) == '\0')
		{
			if((*(atm->name + 1) == 'D' || *(atm->name + 1) == 'E'))
			&& (*(atm->name + 2) == '1' || *(atm->name + 2) == '2'))
				printf("#Found ambiguity of atoms: %s", str);
		}
		*/

		atm->altLoc = str[16];

		memcpy(strTmp, str+30, 8);
		atm->x = atof(strTmp);
		memcpy(strTmp, str+38, 8);
		atm->y = atof(strTmp);
		memcpy(strTmp, str+46, 8);
		atm->z = atof(strTmp);

		atm->next = NULL;

		if(res->atm == NULL)
			res->atm = atm;
		else	lastatm->next = atm;

		lastatm = atm;
	}

	if(zipflg == TRUE)
		pclose(fp);
	else	fclose(fp);


	return chain - chain0 + 1;
}
#endif



/*
Return number of peptide chains if no error was found when parsing ATOM and HETATM records
in *file (PDB file), else return zero.

pdbHETATOM() processes only data holded in ATOM and HETATM record lines and removes
hydrogen atoms if any.
*/

#define LEN_ATOM	82
/* number of columns in ATOM record (80) + '\n' + '0', 82 = 80 + 1 + 1) */

int	pdbHETATOM(char *file, CHAIN *chain0)
{
	FILE	*fp;
	int	zipflg;
	char	*chp;
	char	str[LEN_ATOM];
	char	strTmp[10];
	int	intTmp;
	CHAIN	*chain;
	char	lasthet;
	char	lastchain;
	int	lastSeq;
	char	lastiCode;
	char	lastaltLoc;
	RES	*res, *lastres;
	ATM	*atm, *lastatm;
	

	if((chp = strrchr(file, '.')) == NULL)
	{
		printf("#Unkown type of file: %s\n", file);
		return 0;
	}

	zipflg = FALSE;
	if(chp[1] == 'Z' || (chp[1] == 'g' && chp[2] == 'z'))
	{
		char	cmd[NAME_MAX] = "zcat ";
		int	len;

		zipflg = TRUE;
		len = strlen(file);
		memcpy(cmd+5, file, len+1);

		if((fp = popen(cmd, "r")) == NULL)
		{
			printf("#Open %s error!\n", file);
			return 0;
		}
	}
	else if((fp = fopen(file, "r")) == NULL)
	{
		printf("#Open %s error!\n", file);
		return 0;
	}

	/*
	for(; fgets(str, LEN_ATOM, fp) != NULL;)
	{
		if(memcmp("ATOM  ", str, 6) == 0)
			break;
	}
	*/

	/* initialization */
	chain = chain0 - 1;
	lasthet = '\0';
	lastchain = '\0';
	lastSeq = INT_MIN; /* minimum value a 'signed int' can hold */
	lastiCode = '\0';
	lastaltLoc = ' ';

	res = NULL;
	lastres = NULL;
	lastatm = NULL;

	for(; fgets(str, LEN_ATOM, fp) != NULL;)
	{
		/* process only the first model in PDB file */
		if(memcmp("ENDMDL", str, 6) == 0)
			break;

		/* process only ATOM and HETATM records */
		if(!(memcmp("ATOM  ", str, 6) == 0 || memcmp("HETATM", str, 6) == 0))
			continue;

		/* skip hydrogen atoms */
		if(str[13] == 'H')
		{
			printf("#Warning: hydrogen atom - %s", str);
			continue;
		}

		if(str[21] != lastchain || str[0] != lasthet)
		{
			lastchain = str[21];
			lasthet = str[0];

			lastSeq = INT_MIN;
			lastiCode = '\0';

			chain ++;
			chain->het = lasthet;
			chain->id = lastchain;
			chain->nres = 0;
			chain->res = NULL;
		}

		memcpy(strTmp, str+22, 4);
		strTmp[4] = '\0';
		if((intTmp = atoi(strTmp)) != lastSeq || str[26] != lastiCode)
		{
			/*
			lastSeq = atoi(strTmp);
			*/
			lastSeq = intTmp;
			lastiCode = str[26];

			chain->nres ++;

			res = (RES *)memalloc(sizeof(RES));

			sscanf(str+17, "%s", strTmp);
			res->name = t2sRes(strTmp);
			res->seq = lastSeq;
			res->iCode = lastiCode;
			res->natm = 0;
			res->atm = NULL;
			res->next = NULL;

			if(chain->res == NULL)
				chain->res = res;
			else	lastres->next = res;

			lastres = res;
		}

		if(str[16] != ' ')
		{
			if(lastaltLoc == ' ')
				lastaltLoc = str[16];
			else if(str[16] != lastaltLoc)
				continue;
		}

		/*
		if(lastaltLoc != ' ' && str[16] != lastaltLoc)
			continue;

		if(str[16] != ' ' && lastaltLoc == ' ')
		{
			lastaltLoc = str[16];
		}
		*/


		res->natm ++;

		atm = (ATM *)memalloc(sizeof(ATM));

		atm->serial = atoi(str+6);

		memcpy(strTmp, str+12, 2);
		strTmp[2] = '\0';
		sscanf(strTmp, "%s", atm->chem);

		memcpy(strTmp, str+12, 4);
		strTmp[4] = '\0';
		sscanf(strTmp, "%s", atm->name);

		atm->altLoc = str[16];

		memcpy(strTmp, str+30, 8);
		atm->x = atof(strTmp);
		memcpy(strTmp, str+38, 8);
		atm->y = atof(strTmp);
		memcpy(strTmp, str+46, 8);
		atm->z = atof(strTmp);

		atm->next = NULL;

		if(res->atm == NULL)
			res->atm = atm;
		else	lastatm->next = atm;

		lastatm = atm;
	}

	if(zipflg == TRUE)
		pclose(fp);
	else	fclose(fp);


	return chain - chain0 + 1;
}

/*
refer to PDB format description version 2.2 and Appendix5(formulas and molecular weights
for standard residues) */





/* redundant definition starts here, somewhat obsolete */

static const char	stdaa[] =	"A"	/* ALA */
					"R"	/* ARG */
					"N"	/* ASN */
					"D"	/* ASP */
					"C"	/* CYS */

					"Q"	/* GLN */
					"E"	/* GLU */
					"G"	/* GLY */
					"H"	/* HIS */
					"I"	/* ILE */

					"L"	/* LEU */
					"K"	/* LYS */
					"M"	/* MET */
					"F"	/* PHE */
					"S"	/* SER */

					"T"	/* THR */
					"W"	/* TRP */
					"Y"	/* TYR */
					"V"	/* VAL */
					"P";	/* PRO */
			


/* The number of atoms for the corresponding residues above is as follows. */

static int	resatm[] = {	5,
				11,
				8,
				8,
				6,

				9,
				9,
				4,
				10,
				8,

				8,
				9,
				8,
				11,
				6,

				7,
				14,
				12,
				7,
				7
			};


/* The description of all atoms for the corresponding residues above is as follows. */

static char	*stdatm[] = {	"N","CA","C","O","CB",
				"N","CA","C","O","CB","CG","CD","NE","CZ","NH1","NH2",
				"N","CA","C","O","CB","CG","OD1","ND2",
				"N","CA","C","O","CB","CG","OD1","OD2",
				"N","CA","C","O","CB","SG",

				"N","CA","C","O","CB","CG","CD","OE1","NE2",
				"N","CA","C","O","CB","CG","CD","OE1","OE2",
				"N","CA","C","O",
				"N","CA","C","O","CB","CG","ND1","CD2","CE1","NE2",
				"N","CA","C","O","CB","CG1","CG2","CD1",

				"N","CA","C","O","CB","CG","CD1","CD2",
				"N","CA","C","O","CB","CG","CD","CE","NZ",
				"N","CA","C","O","CB","CG","SD","CE",
				"N","CA","C","O","CB","CG","CD1","CD2","CE1","CE2","CZ",
				"N","CA","C","O","CB","OG",

				"N","CA","C","O","CB","OG1","CG2",
				"N","CA","C","O","CB","CG","CD1","CD2","NE1","CE2","CE3","CZ2","CZ3","CH2",
				"N","CA","C","O","CB","CG","CD1","CD2","CE1","CE2","CZ","OH",
				"N","CA","C","O","CB","CG1","CG2",
				"N","CA","C","O","CB","CG","CD"
			};
