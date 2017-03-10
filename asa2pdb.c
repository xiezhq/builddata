#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "pdb.h"
#include "cnctarea2.h"
#include "asa2pdb.h"
#include "open_file.h"
/*
#include "vdw.h"
*/
#include "memory.h"



void	pdbArray(CHAIN *chain, int nchain, CHAIN *pdb, RES *residue, ATM *atom)
{
	CHAIN	*chaint;
	RES	*res;
	ATM	*atm;

	
	for(chaint = chain + nchain; chain < chaint; chain ++)
	{
		pdb->id = chain->id;
		pdb->nres = chain->nres;
		pdb->res = residue;
		
		for(res = chain->res; res != NULL; res = res->next)
		{
			residue->name = res->name;
			residue->seq = res->seq;
			residue->iCode = res->iCode;
			residue->natm = res->natm;
			residue->atm = atom;

			for(atm = res->atm; atm != NULL; atm = atm->next)
			{
				memcpy(atom->name, atm->name, strlen(atm->name) + 1);
				atom->altLoc = atm->altLoc;
				atom->x = atm->x;
				atom->y = atm->y;
				atom->z = atm->z;

				atom ++;
			}

			residue ++;
		}

		pdb ++;
	}
}


/* get hybrid status of atoms and the corresponding atomic radii

HYBRID	*hybrid: all types of hybrids (the last element of hybrid array is the default value) 

return: number of types of hybrids of atoms
*/
int	hybrid2rad(char *file, HYBRID *hybrid)
{
	FILE	*fp;
	char	str[100];
	int	i;


	fp = open_file(file, "r");

	i = 0;
	for(fgets(str, 100, fp); !feof(fp); fgets(str, 100, fp))
	{
		sscanf(str, "%s %lf", hybrid[i].hybrid, &(hybrid[i].r));
		i ++;
	}

	fclose(fp);

	return i;
}


/* read from file the parameters for atomic VDW radii defined by Gerstein etc. */
int	atm2hybrid(char *file, RESRAD *radi)
{
	FILE	*fp;
	char	str[100];
	RESRAD	*rad;
	ATMRAD	*atm, *lastatm;
	int	i;


	fp = open_file(file, "r");

	rad = radi;

	rad --;
	lastatm = NULL;
	for(fgets(str, 100, fp); memcmp(str, "END", 3) != 0; fgets(str, 100, fp))
	{
		if(str[0] == '#')
			continue;

		/* skip blank line */
		for(i = 0; str[i] != '\n'; i ++)
		{
			if(!isspace(str[i]))
				break;
		}
		if(str[i] == '\n')
			continue;

		/* check and process residue-line */
		if(isalpha(str[0]))
		{
			rad ++;

			sscanf(str, "%s %d", rad->namet, &(rad->natm));
			rad->name = t2sRes(rad->namet);
			rad->atm = NULL;
			continue;
		}

		/* process atom-line */
		atm = (ATMRAD *)memalloc(sizeof(ATMRAD));
		sscanf(str, "%s %s", atm->name, atm->hybrid);
		atm->next = NULL;

		if(rad->atm == NULL)
			rad->atm = atm;
		else	lastatm->next = atm;

		lastatm = atm;
	}

	fclose(fp);

	return rad - radi + 1;
}


void	getRad(HYBRID *hybrid, int nhybrid, RESRAD *rad, int nresrad)
{
	RESRAD	*radt;
	ATMRAD	*atm;
	int	i;


	nhybrid = nhybrid - 1;
	for(radt = rad + nresrad; rad < radt; rad ++)
	{
		for(atm = rad->atm; atm != NULL; atm = atm->next)
		{
			for(i = 0; i < nhybrid; i ++)
			{
				if(strcmp(atm->hybrid, hybrid[i].hybrid) == 0)
				{
					atm->r = hybrid[i].r;

					break;
				}
			}

			if(i == nhybrid)
			{
				/*
				atm->r = RAD_DEFAULT;
				*/
				atm->r = hybrid[i].r;

				printf("#Warning: hybrid status of the atom below wasn't found, using the defualt value for its vdw radius instead:\n"
					"%s(%c) %s\n",
					rad->namet, rad->name, atm->name);
			}
		}
	}
}


/* assign atomic radii to atoms in structure */
void	rad2pdb(CHAIN *chain, int nchain, RESRAD *radi, int nrad)
{
	CHAIN	*chaint;
	RES	*res;
	ATM	*atm;
	RESRAD	*rad, *radt;
	ATMRAD	*atmx;


	/* The radt below is not the specific residue but the collection of some special atoms.
	Refer to parameter file for atomic VDW radii.
	*/
	radt = radi + nrad - 1;

	for(chaint = chain + nchain; chain < chaint; chain ++)
	{
		for(res = chain->res; res != NULL; res = res->next)
		{
			/*
			printf("res %c %d%c\n", res->name, res->seq, res->iCode);
			*/
			for(rad = radi; rad < radt; rad ++)
			{
				/*
				printf("res->name %c rad->name %c\n", res->name, rad->name);
				*/

				if(res->name != rad->name || res->name == 'X')
					continue;

				for(atm = res->atm; atm != NULL; atm = atm->next)
				{
					for(atmx = rad->atm; atmx != NULL; atmx = atmx->next)
					{
						if(strcmp(atm->name, atmx->name) == 0)
						{
							atm->r = atmx->r;
							break;
						}
					}

					if(atmx == NULL)
					{
						/* assign radii to special atoms */
						for(atmx = radt->atm; atmx != NULL; atmx = atmx->next)
						{
							if(strcmp(atm->name, atmx->name) == 0)
							{
								atm->r = atmx->r;
								break;
							}
						}

						if(atmx == NULL)
						{
							/* assign radii based on chemical symbol */
							for(atmx = radt->atm; atmx != NULL; atmx = atmx->next)
							{
								if(*atm->chem == *atmx->name
								&& memcmp(atmx->name + 1, "DUM", 4) == 0)
								{
									printf("#Warning: assign radii based on chemical symbol of atom (%c %d%c %s)\n",
										chain->id, res->seq, res->iCode, atm->name);

									atm->r = atmx->r;
									break;
								}
							}

							/* assign default value */
							if(atmx == NULL)
							{
								printf("#Warning: unknown atom (%c %d%c %s) for assignment of atomic VDW radii\n"
									"\t %.3f was used by default (you can change the value of RAD_DEFAULT in asa2pdb.h)\n",
									chain->id, res->seq, res->iCode, atm->name, RAD_DEFAULT);
							}
						}
					}
				}
				/*
				printf("res->name %c rad->name %c\n", res->name, rad->name);
				*/

				break;
			}

			if(rad == radt)
			{
				printf("#Warning: assign default values of atomic radius for atoms in unknown residue (%c %d%c)\n",
					chain->id, res->seq, res->iCode);

				/* assign radii value to special atoms */
				for(atm = res->atm; atm != NULL; atm = atm->next)
				{
					for(atmx = radt->atm; atmx != NULL; atmx = atmx->next)
					{
						if(strcmp(atm->name, atmx->name) == 0)
						{
							atm->r = atmx->r;
							break;
						}
					}

					if(atmx == NULL)
					{
						/* assign radii based on chemical sybol */
						for(atmx = radt->atm; atmx != NULL; atmx = atmx->next)
						{
							if(*atm->chem == *atmx->name
							&& memcmp(atmx->name + 1, "DUM", 4) == 0)
							{
								atm->r = atmx->r;
								break;
							}
						}

						/* assign default value */
						if(atmx == NULL)
						{
							printf("#Warning: unknown atom (%c %d%c %s) for assignment of atomic VDW radii\n"
								"\t %.3f was used by default (you can change the value of RAD_DEFAULT in asa2pdb.h)\n",
								chain->id, res->seq, res->iCode, atm->name, RAD_DEFAULT);
						}
					}
				}
			}
		}
	}
}


void	rad2pdbVer2(CHAIN *pdb, int nchain, RESRAD *rad, int nresrad)
{
	int	i, j;
	CHAIN	*pdbt;
	RES	*res;
	ATM	*atm;
	RESRAD	*radt;
	ATMRAD	*atmr;


	for(pdbt = pdb + nchain; pdb < pdbt; pdb ++)
	{
		res = pdb->res;
		for(i = 0; i < pdb->nres; i ++)
		{
			for(radt = rad + nresrad; rad < radt; rad ++)
			{
				if(res->name != rad->name)
					continue;

				atm = res->atm;
				for(j = 0; j < res->natm; j ++)
				{
					for(atmr = rad->atm; atmr != NULL; atmr = atmr->next)
					{
						if(strcmp(atm->name, atmr->name) == 0)
						{
							atm->r = atmr->r;
							break;
						}
					}

					if(atmr == NULL)
					{
						printf("#Warning: no radius was found for:\n"
							"chain %c residue %d%c atom %s\n",
							pdb->id, res->seq, res->iCode, atm->name);
					}
				}

				break;
			}

			if(rad >= radt)
			{
				printf("#Warning: no residue in residue-list with atomic VDW radii matched:\n"
					"chain %c residue %d%c\n",
					pdb->id, res->seq, res->iCode); 
			}
		}
	}
}



void	freeRad(RESRAD *rad, int n)
{
	RESRAD	*radt;
	ATMRAD	*atm, *nextatm;


	for(radt = rad + n; rad < radt; rad ++)
	{
		for(atm = rad->atm; atm != NULL; atm = nextatm)
		{
			nextatm = atm->next;

			free((void *)atm);
		}
	}
}


/* count the number of elements in protein unit, such as residues and atoms.

chain0	full details for each atom, residue and chain in computed protein unit
nchain	number of chains
nres	number of residues

return: number of atoms
*/
int	nele(CHAIN *chain, int nchain, int *nres)
{
	int	natm;
	CHAIN	*chaint;
	RES	*res;


	*nres = 0;
	natm = 0;
	for(chaint = chain + nchain; chain < chaint; chain ++)
	{
		(*nres) += chain->nres;

		for(res = chain->res; res != NULL; res = res->next)
		{
			natm += res->natm;
		}
	}

	return natm;
}


/* calculate the solvent accessible area for each atom in the computed part of protein

CHAIN	*chain:	full details for each atom, residue and chain in computed protein unit
int	nchain:	number of chains in computed protein
*/
void	asa2pdb(CHAIN *chain, int nchain, int natm)
{
	CHAIN	*chaint;
	RES	*res;
	ATM	*atm;
	double	probe[natm];
	double	x[natm], y[natm], z[natm];
	double	radx[natm];
	double	accss[natm];
	int	i;


	i = 0;
	for(chaint = chain + nchain; chain < chaint; chain ++)
	{
		for(res = chain->res; res != NULL; res = res->next)
		{
			for(atm = res->atm; atm != NULL; atm = atm->next)
			{
				probe[i] = atm->r + Rwat;
				x[i] = atm->x;
				y[i] = atm->y;
				z[i] = atm->z;
				radx[i] = atm->r;

				i ++;
			}
		}
	}

	if(i != natm)
	{
		printf("Error: atom mismatch in asa calculation (natm (%d) -> i (%d))\n", natm, i);
	}

	/*
	printf("# The radius of probe molecule (water): %.1f\n", Rwat);
	*/
	cnctarea(ERRORPAR, probe, x, y, z, radx, natm, accss);

	i = 0;
	for(chain = chaint - nchain; chain < chaint; chain ++)
	{
		for(res = chain->res; res != NULL; res = res->next)
		{
			res->asa = 0.0;
			for(atm = res->atm; atm != NULL; atm = atm->next)
			{
				atm->accss = accss[i];
				res->asa += accss[i];

				i ++;
			}
		}
	}
}


void	asa2pdbVer2(ATM *atm, int natm, double *accss)
{
	int	i;
	double	probe[natm];
	double	x[natm];
	double	y[natm];
	double	z[natm];
	double	radx[natm];


	for(i = 0; i < natm; i ++)
	{
		probe[i] = atm[i].r + Rwat;
		x[i] = atm[i].x;
		y[i] = atm[i].y;
		z[i] = atm[i].z;
		radx[i] = atm[i].r;
	}

	cnctarea(ERRORPAR, probe, x, y, z, radx, natm, accss);
}


/* All atoms are defined as sidechain atom except four atoms, N, CA, C, O. 
The only exception is GLY since it has no normal sidechain. So the CA atom
is defined as the sidechain of GLY for convenience of further calculation.
*/
void	asa4sidechain(CHAIN *chain, int nchain)
{
	CHAIN	*chaint;
	RES	*res;
	ATM	*atm;


	for(chaint = chain + nchain; chain < chaint; chain ++)
	{
		for(res = chain->res; res != NULL; res = res->next)
		{
			res->asa = 0.0;

			if(res->name == 'G')
			{
				for(atm = res->atm; atm != NULL; atm = atm->next)
				{
					if(memcmp(atm->name, "N", 2) == 0
					|| memcmp(atm->name, "C", 2) == 0
					|| memcmp(atm->name, "O", 2) == 0)
						continue;

					res->asa += atm->accss;
				}
			}
			else
			{
				for(atm = res->atm; atm != NULL; atm = atm->next)
				{
					if(memcmp(atm->name, "N", 2) == 0
					|| memcmp(atm->name, "CA", 3) == 0
					|| memcmp(atm->name, "C", 2) == 0
					|| memcmp(atm->name, "O", 2) == 0)
						continue;

					res->asa += atm->accss;
				}
			}
		}
	}
}


int	rdStdAsa(char *file, STDASA *stdAsa)
{
	FILE	*fp;
	char	str[100];
	STDASA	*stdAsai;
	int	i;
	char	name[10];


	fp = open_file(file, "r");

	stdAsai = stdAsa;

	for(fgets(str, 100, fp); memcmp(str, "END", 3) != 0; fgets(str, 100, fp))
	{
		if(str[0] == '#') continue;

		for(i = 0; str[i] != '\n'; i ++)
		{
			if(!isspace(str[i])) break;
		}
		if(str[i] == '\n') continue;

		sscanf(str, "%s", name);
		stdAsa->name = t2sRes(name);
		stdAsa->asa = atof(str+7);
		stdAsa ++;
	}

	fclose(fp);

	return stdAsa - stdAsai;
}


void	acc(CHAIN *chain, int nchain, STDASA *stdAsa, int nstdAsa)
{
	CHAIN	*chaint;
	RES	*res;
	STDASA	*stdAsai, *stdAsat;


	stdAsai = stdAsa;
	stdAsat = stdAsa + nstdAsa;

	for(chaint = chain + nchain; chain < chaint; chain ++)
	{
		for(res = chain->res; res != NULL; res = res->next)
		{
			for(stdAsa = stdAsai; stdAsa < stdAsat; stdAsa ++)
			{
				if(res->name == stdAsa->name)
				{
					res->acc = res->asa / stdAsa->asa;
					break;
				}
			}

			if(stdAsa == stdAsat)
			{
				printf("#Warning: no standard ASA value for residue %c (%c %4d%c)\n"
					"relative accessibility of 9.0 was assigned to this residue by default\n",
					res->name, chain->id, res->seq, res->iCode); 

				res->acc = 9.0;
			}
		}
	}
}


/* calculate relative solvent accessibility for N-terminal residue of each chain */
void	accN(CHAIN *chain, int nchain, STDASA *stdAsa, int nstdAsa)
{
	CHAIN	*chaint;
	RES	*res;
	STDASA	*stdAsai, *stdAsat;


	stdAsai = stdAsa;
	stdAsat = stdAsa + nstdAsa;

	for(chaint = chain + nchain; chain < chaint; chain ++)
	{
		res = chain->res;

		for(stdAsa = stdAsai; stdAsa < stdAsat; stdAsa ++)
		{
			if(res->name == stdAsa->name)
			{
				res->acc = res->asa / stdAsa->asa;
				break;
			}
		}

		if(stdAsa == stdAsat)
		{
			printf("#Warning: no standard ASA value for residue %c (%c %4d%c)\n"
				"relative accessibility of 9.0 was assigned to this residue by default\n",
				res->name, chain->id, res->seq, res->iCode); 

			res->acc = 9.0;
		}
	}
}


/* calculate the relative solvent accessibility for C-terminal residue of each chain */
void	accC(CHAIN *chain, int nchain, STDASA *stdAsa, int nstdAsa)
{
	CHAIN	*chaint;
	RES	*res;
	STDASA	*stdAsai, *stdAsat;


	stdAsai = stdAsa;
	stdAsat = stdAsa + nstdAsa;

	for(chaint = chain + nchain; chain < chaint; chain ++)
	{
		if((res = chain->res) == NULL) continue;

		while(res->next != NULL) res = res->next;

		for(stdAsa = stdAsai; stdAsa < stdAsat; stdAsa ++)
		{
			if(res->name == stdAsa->name)
			{
				res->acc = res->asa / stdAsa->asa;
				break;
			}
		}

		if(stdAsa == stdAsat)
		{
			printf("#Warning: no standard ASA value for residue %c (%c %4d%c)\n"
				"relative accessibility of 9.0 was assigned to this residue by default\n",
				res->name, chain->id, res->seq, res->iCode); 

			res->acc = 9.0;
		}
	}
}
