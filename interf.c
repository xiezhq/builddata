#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include "pdb.h"
#include "memory.h"
#include "max.h"
#include "interf.h"
#include "asa2pdb.h"
#include "open_file.h"
#include "calzscore.h"
#include "property2res.h"


static void	usage(char *prog)
{
	printf("Usage: %s pdbfile4complex numberOfChainsInProtein1 numberOfChainsInProtein2 "
		"chainsInProtein1 chainsInProtein2 cutoff4surf cutoff4interf atm2hybridFile hybrid2radFile "
		"stdAsaFile outputFile\n"
		"pdbfile4complex: PDB file of protein complex\n"
		"numberOfChainsInProtein1: number of chains in protein1\n"
		"numberOfChainsInProtein2: number of chains in protein2\n"
		"chainsInProtein1: list of chains of protein1 in complex\n"
		"chainsInProtein2: list of chains of protein2 in complex\n"
		"cutoff4surf: residues are defined as surface residues when the residue has a relative accessibility greater than cutoff4surf\n"
		"cutoff4interf: residues are defined as interface residues when the monomer loses at least cutoff4interf A**2 the ASA upon association\n"
		"atm2hybridFile: describes the hybrids of all atoms for all normal residues\n"
		"hybrid2radFile: radii for all types of hybrid atoms\n"
		"stdAsaFile: file containing the value of ASA of normal residues in tripeptide\n"
		"outputFile: output file\n"

		"Note for chains_in_protein1 and chains_in_protein2:\n"
		"The '_' symbol must be presented as the chain identifier if the chain identifier is ' ' in PDB file.\n"
		"For example:\n"
		"interf pdb1i9R.ent 2 3 HL ABC 0.05 1.0 ProtOr.resi-defs.mod-by-xie.dat ProtOr.atom-defs.mod-by-xie.dat glyXgly.ProtOr.dat pdb1i9R.acc\n"
		"interf pdb1avx.ent 1 1 A B 0.05 1.0 ProtOr.resi-defs.mod-by-xie.dat ProtOr.atom-defs.mod-by-xie.dat glyXgly.ProtOr.dat pdb1avx.acc\n",
		prog);
}


static void	getPar(char *argv[], PAR *par)
{
	int	i;
	int	fileNameLen;
	CHAINID	*chainid, *lastchainid;


	par->nchain1 = atoi(argv[2]);
	par->nchain2 = atoi(argv[3]);
	par->surfCutoff = atof(argv[6]);
	par->interfCutoff = atof(argv[7]);
	fileNameLen = strlen(argv[8]);
	memcpy(par->atm2hybridFile, argv[8], fileNameLen + 1);
	/*
	*(par->atm2hybridFile + fileNameLen) = '\0';
	*/
	fileNameLen = strlen(argv[9]);
	memcpy(par->hybrid2radFile, argv[9], fileNameLen + 1);
	/*
	*(par->hybrid2radFile + fileNameLen) = '\0';
	*/

	lastchainid = NULL;
	par->chainId1 = NULL;
	for(i = 0; i < par->nchain1; i ++)
	{
		chainid = (CHAINID *)memalloc(sizeof(CHAINID));

		chainid->id = argv[4][i];
		chainid->next = NULL;

		if(par->chainId1 == NULL)
			par->chainId1 = chainid;
		else	lastchainid->next = chainid;

		lastchainid = chainid;
	}

	par->chainId2 = NULL;
	for(i = 0; i < par->nchain2; i ++)
	{
		chainid = (CHAINID *)memalloc(sizeof(CHAINID));

		chainid->id = argv[5][i];
		chainid->next = NULL;

		if(par->chainId2 == NULL)
			par->chainId2 = chainid;
		else	lastchainid->next = chainid;

		lastchainid = chainid;
	}
}


static void	freePar(PAR *par)
{
	CHAINID	*id, *nextid;


	for(id = par->chainId1; id != NULL; id = nextid)
	{
		nextid = id->next;

		free((void *)id);
	}

	for(id = par->chainId2; id != NULL; id = nextid)
	{
		nextid = id->next;

		free((void *)id);
	}
}


/*
check the identifiers from PDB file (CHAIN * chain) and parameter file (PAR *par), respectively,
to see if match each other

Note:
. it doesn't check the chains originated from HETATM record in PDB file.
*/
static void	checkChain4PdbAndPar(CHAIN *chain, int nchain, PAR *par)
{
	CHAIN	*chaint;
	CHAINID	*id;
	/*
	char	idchain[nchain], idid[par->nchain1 + par->nchain2];
	*/
	char	*idchain, *idid;
	int	i, j;
	int	nid;


	idchain = (char *)memalloc(nchain * sizeof(char));
	idid = (char *)memalloc((par->nchain1 + par->nchain2) * sizeof(char));

	i = 0;
	for(chaint = chain + nchain; chain < chaint; chain ++)
	{
		idchain[i] = chain->id;
		i ++;
	}

	j = 0;
	for(id = par->chainId1; id != NULL; id = id->next)
	{
		idid[j] = id->id;
		j ++;
	}
	for(id = par->chainId2; id != NULL; id = id->next)
	{
		idid[j] = id->id;
		j ++;
	}
	nid = j;

	for(i = 0; i < nchain; i ++)
	{
		for(j = 0; j < nid; j ++)
		{
			if(idchain[i] == idid[j])
			{
				nchain --;
				nid --;
				if(i < nchain)
					idchain[i] = idchain[nchain];
				if(j < nid)
					idid[j] = idid[nid];
				i --;
				break;
			}
		}
	}

	if(nchain > 0)
	{
		printf("#Warning: found chain(s) not belonging to any one of two monomers: ");
		for(i = 0; i < nchain; i ++)
		{
			printf("%c", idchain[i]);
		}
		printf("\n");
	}
	free((void *)idchain);

	if(nid > 0)
	{
		printf("#Error: found chain(s) not existing in original PDB file: ");
		for(j = 0; j < nid; j ++)
		{
			printf("%c", idid[j]);
		}
		printf("\n");

		printf("Program exit!\n");
		exit(1);
	}
	free((void *)idid);
}


/*
#define	TRI	3

static void	stdASAsitu(CHAIN *chain, int nchain)
{
	int	i;
	double	*probe, *x, *y, *z, *radx;
	CHAIN	*chaint;
	RES	*res;
	ATM	*atm;


	for(chaint = chain + n; chain < chaint; chain ++)
	{
		i = 0;
		for(res = chain->res; res != NULL; res = res->next)
		{
			
		}
	}
}
*/


static int	divComplex(CHAIN *chain, int nchain, CHAIN *chain1, CHAINID *idInit)
{
	int	nres1;
	CHAIN	*chaint;
	CHAINID	*id;
	RES	*res, *res1, *lastres1;
	ATM	*atm, *atm1, *lastatm1;


	lastres1 = NULL;
	lastatm1 = NULL;
	nres1 = 0;
	for(chaint = chain + nchain; chain < chaint; chain ++)
	{
		for(id = idInit; id != NULL; id = id->next)
		{
			if(id->id != chain->id)
				continue;

			chain1->id = chain->id;
			chain1->nres = chain->nres;
			chain1->het = chain->het;
			chain1->res = NULL;

			for(res = chain->res; res != NULL; res = res->next)
			{
				nres1 ++;

				res1 = (RES *)memalloc(sizeof(RES));
				res1->name = res->name;
				res1->seq = res->seq;
				res1->iCode = res->iCode;
				res1->natm = res->natm;
				res1->atm = NULL;
				res1->next = NULL;

				if(chain1->res == NULL)
					chain1->res = res1;
				else	lastres1->next = res1;

				for(atm = res->atm; atm != NULL; atm = atm->next)
				{
					atm1 = (ATM *)memalloc(sizeof(ATM));
					memcpy(atm1->name, atm->name, strlen(atm->name)+1);
					atm1->altLoc = atm->altLoc;
					atm1->x = atm->x;
					atm1->y = atm->y;
					atm1->z = atm->z;
					atm1->r = atm->r;
					atm1->next = NULL;

					if(res1->atm == NULL)
						res1->atm = atm1;
					else	lastatm1->next = atm1;

					lastatm1 = atm1;
				}

				lastres1 = res1;
			}

			chain1 ++;
			break;
		}
	}

	return nres1;
}


/* get the value of the lost ASA for each atom and residue upon the association of two protein units

CHAIN *dAccss:
dAccss.res.atm.accss	accss_unbound - accss_bound for one atom
dAccss.res.asa		asa_unbound - asa_bound for one residue

Must guaraintee before call deltaAccssVer1():
1. Each chain in bound form has the same set of residues as that of its corresponding chain
	in unbound form.
&
2. Each residue in bound form has the same set of atoms as that of its corresponding residue
	in unbound form.
*/
static void	deltaAccssVer1(CHAIN *chain, int nchain, CHAIN *chain1, int nchain1, CHAIN *dAccss)
{
	CHAIN	*chaint;
	CHAIN	*chain1i, *chain1t;
	RES	*res, *res1;
	RES	*dres, *lastdres;
	ATM	*atm, *atm1;
	ATM	*datm, *lastdatm;


	lastdres = NULL;
	lastdatm = NULL;
	chain1i = chain1;
	chain1t = chain1 + nchain1;
	for(chaint = chain + nchain; chain < chaint; chain ++)
	{
		for(chain1 = chain1i; chain1 < chain1t; chain1 ++)
		{
			if(chain->id == chain1->id
			&& chain->nres == chain1->nres
			&& chain->het == chain1->het)
			{
				break;
			}
		}

		if(chain1 == chain1t)
			continue;

		/* Must guaraintee:
		two chains to compare have same set of residues
		&
		two residues to compare have same set of atoms
		*/
		dAccss->id = chain->id;
		dAccss->nres = chain->nres;
		dAccss->res = NULL;
		for(res = chain->res, res1 = chain1->res; res != NULL; res = res->next, res1 = res1->next)
		{
			if(res->name != res1->name || res->natm != res1->natm)
			{
				printf("#Error in deltaAccssVer1: the compared chains have different set of residues "
					"(%c %d%c %c : %c %d%c %c)\n",
					chain->id, res->seq, res->iCode, res->name,
					chain1->id, res1->seq, res1->iCode, res1->name);
				exit(1);
			}

			dres = (RES *)memalloc(sizeof(RES));
			dres->name = res->name;
			dres->seq = res->seq;
			dres->iCode = res->iCode;
			dres->natm = res->natm;
			dres->asa = res1->asa - res->asa;
			dres->atm = NULL;
			dres->next = NULL;

			if(dAccss->res == NULL)
				dAccss->res = dres;
			else	lastdres->next = dres;

			for(atm = res->atm, atm1 = res1->atm; atm != NULL; atm = atm->next, atm1 = atm1->next)
			{
				if(strcmp(atm->name, atm1->name) !=0)
				{
					printf("#Error in deltaAccssVer1: the compared residues have different set of atoms "
						"(%c %d%c %c %s : %c %d%c %c %s)\n",
						chain->id, res->seq, res->iCode, res->name, atm->name,
						chain1->id, res1->seq, res1->iCode, res1->name, atm1->name);
					exit(1);
				}

				datm = (ATM *)memalloc(sizeof(ATM));
				memcpy(datm->name, atm->name, strlen(atm->name)+1);
				datm->altLoc = atm->altLoc;

				/*
				datm->x = atm->x;
				datm->y = atm->y;
				datm->z = atm->z;
				datm->r = atm->r;
				*/

				datm->accss = atm1->accss - atm->accss;
				datm->next = NULL;

				if(dres->atm == NULL)
					dres->atm = datm;
				else	lastdatm->next = datm;

				/*
				dres->asa += datm->accss;
				*/

				lastdatm = datm;
			}

			lastdres = dres;
		}
		
		dAccss ++;
	}
}


static void	getArray(CHAIN *chain, int nchain, double *asa, double *acc)
{
	CHAIN	*chaint;
	RES	*res;
	/*
	ATM	*atm;
	*/


	for(chaint = chain + nchain; chain < chaint; chain ++)
	{
		for(res = chain->res; res != NULL; res = res->next)
		{
			/*
			for(atm = res->atm; atm != NULL; atm = atm->next)
			{
				*accss++ = atm->accss;
			}
			*/

			*asa++ = res->asa;
			*acc++ = res->acc;
		}
	}
}


static void	diff(double *data1, double *data2, int n, double *diff)
{
	int	i;


	for(i = 0; i < n; i ++)
	{
		*diff++ = *data2++ - *data1++;
	}
}


static int	defSurf(CHAIN  *chain, int nchain, PAR *par, SURF *surf)
{
	CHAIN	*chaint;
	RES	*res;
	SURF	*surfi;


	surfi = surf;
	for(chaint = chain + nchain; chain < chaint; chain ++)
	{
		for(res = chain->res; res != NULL; res = res->next)
		{
			if(res->acc > par->surfCutoff)
			{
				surf->name = res->name;
				surf->chainid = chain->id;
				surf->seq = res->seq;
				surf->iCode = res->iCode;
				surf->asa = res->asa;
				surf->acc = res->acc;
				surf ++;
			}
		}
	}

	return surf - surfi;
}


static void	surfAsa2array(SURF *surf, int nsurf, double *asa)
{
	SURF	*surft;


	for(surft = surf + nsurf; surf < surft; surf ++)
	{
		*asa++ = surf->asa;
	}
}


static void	surfAcc2array(SURF *surf, int nsurf, double *acc)
{
	SURF	*surft;


	for(surft = surf + nsurf; surf < surft; surf ++)
	{
		*acc++ = surf->acc;
	}
}


static int	defInterf(CHAIN *chain, int nchain, SURF *surf, int nsurf, PAR *par)
{
	CHAIN	*chaint, *chaini;
	RES	*res;
	SURF	*surft;
	int	i;


	i = 0;
	chaini = chain;
	chaint = chain + nchain;
	for(surft = surf + nsurf; surf < surft; surf ++)
	{
		 for(chain = chaini; chain < chaint; chain ++)
		{
			if(chain->id != surf->chainid)
				continue;

			for(res = chain->res; res != NULL; res = res->next)
			{
				if(res->seq != surf->seq || res->iCode != surf->iCode)
					continue;

				if(res->asa > par->interfCutoff)
				{
					surf->intflg = TRUE;
					i ++;
				}
				else	surf->intflg = FALSE;

				break;
			}

			break;
		}
	}

	return i;
}


static void	outputSurf(int nres, int ninterf,
			SURF *surf, int nsurf, double *asa, double *zasa, double *acc, double *zacc,
			char *file)
{
	FILE	*fp;
	int	i;
	SURF	*surfi, *surft;
	char	*chtmp;
	char	log_path[NAME_MAX];
	time_t	time0;
	char	time1[100];


	/* create all the non-existing parent directories specified by the path of file */
	if((chtmp = strrchr(file, '/')) != NULL)
	{
		memcpy(log_path, "mkdir -p ", strlen("mkdir -p ") + 1);
		strcat(log_path, file);
		*strrchr(log_path, '/') = '\0';
		system(log_path);
	}

	fp = open_file(file, "w");

	time(&time0);
	ctime_r(&time0, time1);
	if(chtmp == NULL)
		chtmp = file;
	else	chtmp ++;
	fprintf(fp,
		"# %s was produced by interf at %s\n"
		"# Monomer: nsurf/nres = %d/%d (%.3f), ninterf/nsurf = %d/%d (%.3f)\n"
		"# nsurf: number of surface residues\n"
		"# nres: number of total residues\n"
 		"# ninterf: number of interface residues\n",
		chtmp, time1,
		nsurf, nres, ((double)nsurf)/nres,
		ninterf, nsurf, ((double)ninterf)/nsurf);
	fprintf(fp,"\n"
		"# The records below list all surface residues.\n"
		"# Records format\n"
		"# %-10s %s\n" 
		"# %2d%8c Chain identifier\n"
		"# %2d%8c Residue name\n"
		"# %2d -%3d%3c Residue sequence number\n"
		"# %2d%8c Code for insertion of residues\n"
		"# %2d -%3d%3c Solvent accessible surface area (SASA) of residue\n"
		"# %2d -%3d%3c Relative accessibility of residue\n"
		"# %2d -%3d%3c Z-score for SASA of residue\n"
		"# %2d -%3d%3c Z-score for relative accessibility of residue\n"
		"# %2d -%2d%4c Flag for interface residue (1 if interface residue, else 0)\n"
		"# %2d -%2d%4c Flag for polarity of residue (1 for polar residue, -1 for apolar residue, else 0)\n"
		"# %2d -%2d%4c Flag for hydrophility of residue (1 for hydrophilic residue, -1 for hydrophobic residue, else 0)\n"
		"# %2d -%2d%4c Flag for charged residue (1 for positively charged residue, -1 for negatively charged residue, else 0)\n",
		"COLUMNS", "DEFINITIONS",
		1, ' ',
		2, ' ',
		4, 7, ' ',
		8, ' ',
		10, 17, ' ',
		19, 26, ' ',
		28, 35, ' ',
		37, 44, ' ',
		46, 47, ' ',
		49, 50, ' ',
		52, 53, ' ',
		55, 56, ' ');

	for(surft = surf + nsurf; surf < surft; surf ++)
	{
		/*
		printf("%c %4d%c %8.2f %8.3f %8.3f %8.3f %d\n", 
			surf->chainid, surf->seq, surf->iCode,
			surf->asa, surf->acc, *zasa++, *zacc++, surf->intflg);
		*/

		fprintf(fp, "%c%c %4d%c %8.2f %8.3f %8.3f %8.3f %2d %2d %2d %2d\n", 
			surf->chainid, surf->name, surf->seq, surf->iCode,
			surf->asa, surf->acc, *zasa++, *zacc++, surf->intflg,
			surf->pol, surf->philic, surf->e);

		/* Just for test */
		/*
		fprintf(fp, "%c%c %4d%c %8.2f %8.3f %d\n", 
			surf->chainid, surf->name, surf->seq, surf->iCode,
			surf->asa, surf->acc, surf->intflg);
		*/
		/* Just for test */
	}

	fclose(fp);
}


static void	outputZscore(CHAIN *chain, int nchain, double *zacc, int nres)
{
	CHAIN	*chaint;
	RES	*res;
	int	i;


	i = 0;
	for(chaint = chain + nchain; chain < chaint; chain ++)
	{
		for(res = chain->res; res != NULL; res = res->next)
		{
			printf("%c %4d%c %8.2f %8.3f %8.2f\n", 
				chain->id, res->seq, res->iCode, res->asa, res->acc, *zacc++);
			/*
			printf("i = %d\n", i);
			*/
			i ++;
		}
	}

	if(i != nres)
	{
		printf("Error: number of residues doesn't match: %4d-%4d\n", i, nres);
	}
}


int	main(int argc, char *argv[])
{
	PAR	par;
	int	nchain0;
		/* number of chains in complex */

	int	nres0;
		/* number of residues in complex */

	int	natm0;
		/* number of atoms in complex */

	CHAIN	chain0[CHAIN_MAX];
		/* complex */

	CHAIN	*chain01, *chain02;

	int	nchain01, nchain02;	
		/* number of chains in two different protein units */

	int	nres01, nres02;
		/* number of residues in two different protein units (group chain01 and group chain02) */

	int	natm01, natm02;
		/* number of atoms in two different protein units */

	HYBRID	hybrid[HYBRID_MAX];
	int	nhybrid; /* number of types of atomic hybrid */
	RESRAD	rad[N_RESRAD]; /* radii of atoms in residues or het groups */
	int	nresrad; /* number of standard residues, including normal het groups  */

	STDASA	stdAsa[50];
		/* standard asa for each normal residue in tripeptide */
	int	nstdAsa;
	STDASA	stdAsaN[50], stdAsaC[50];
		/* standard asa for each normal residue in terminal dipeptide */
	int	nstdAsaN, nstdAsaC;

	double	*asa01, *asa02; /* ASA for each residue */
	double	*acc01, *acc02; /* relative accessibility for each residue */
	double	*accss01, *accss02; /* ASA for each atom */
	double	*zacc0, *zacc01, *zacc02; /* zscore */
	double	*zasa01;

	SURF	*surf01; /* surface of monomer1 */
	int	nsurf01; /* number of surface residues in monomer1 */
	int	ninterf01; /* number of interface residues in monomer1 */

	CHAIN	*dAccssChain01, *dAccssChain02;
	double	*dAsa0, *dAcc0, *dAccss0;
	double	*zdAsa0, *zdAccss0;

	char	wat; /* single name for water molecule */
	

	if(argc != 12)
	{
		usage(argv[0]);
		exit(1);
	}


	/* read in pdb file */
	printf("#Reading %s\n", argv[1]);
	/*
	if((nchain0 = pdbATOM(argv[1], chain0)) == 0)
	*/
	if((nchain0 = pdbHETATOM(argv[1], chain0)) == 0)
	{
		printf("#Error: number of chains in %s = 0 !\n", argv[1]);
		exit(1);
	}


	/* test code */
	/*
	coordIdentity(chain0, nchain0);
	exit(0);
	*/
	/* test code */


	/* get the parameters of analysis */
	getPar(argv, &par);


	nresrad = atm2hybrid(par.atm2hybridFile, rad);

	/* remove nucleic acid chain */
	nchain0 = rmnuchain(chain0, nchain0);

	/* remove chains containing water molecules */
	/*
	nchain0 = rmWaterChain(chain0, nchain0);
	*/

	/* remove chains appearing in PDB HETATM records */
	/*
	nchain0 = rmhetchain(chain0, nchain0);
	*/

	/* remove only water molecules */
	wat = 'e';
	rmres(chain0, nchain0, wat);
	wat = 'f';
	rmres(chain0, nchain0, wat);


	/* check the chains, residues and atoms in structure */
	nchain0 = checkATOM(chain0, nchain0, rad, nresrad);

	/* check the number and identifiers of chains */
	checkChain4PdbAndPar(chain0, nchain0, &par);

	nhybrid = hybrid2rad(par.hybrid2radFile, hybrid);
	getRad(hybrid, nhybrid, rad, nresrad);

	/* add radii of atoms to structure */
	rad2pdb(chain0, nchain0, rad, nresrad);
	freeRad(rad, nresrad);


	/* get values of solvent accessible surface area (ASA)

	NOTE: cavity should be treated carefully because the algorithm to calculate ASA 
	doesn't remove the contribution of cavity from total ASA. Furthermore, the ASA
	may not mean what you think if all het-groups are removed from PDB structure.
	So the ASA values calculated here are always greater than or equal to expected values!
	*/

	/* calculate solvent accessible surface area */
	natm0 = nele(chain0, nchain0, &nres0);
	asa2pdb(chain0, nchain0, natm0);

	/* asa of sidechain */
	/*
	asa4sidechain(chain0, nchain0);
	*/


	/* define surface residues */

	nstdAsa = rdStdAsa(argv[10], stdAsa);
	nstdAsaN = rdStdAsa("Xgly.ProtOr.dat", stdAsaN);
	nstdAsaC = rdStdAsa("glyX.ProtOr.dat", stdAsaC);

	/* calculate the relative solvent accessibility of residue */
	acc(chain0, nchain0, stdAsa, nstdAsa);
	accN(chain0, nchain0, stdAsaN, nstdAsaN);
	accC(chain0, nchain0, stdAsaC, nstdAsaC);

	/* In an ideal model, the inaccessible residues should have an ASA value of zero, namely, ASA/stdASA = 0.
	But it is also normal to define residues as surface residues if they satisfy ASA/stdASA < surfCutoff,
	where, stdAsa is the ASA in extended conformation for each residue */

	/* monomer1 */
	nchain01 = par.nchain1;
	chain01 = (CHAIN *)memalloc(nchain01 *sizeof(CHAIN));
	nres01 = divComplex(chain0, nchain0, chain01, par.chainId1);

	natm01 = nele(chain01, nchain01, &nres01);
	asa2pdb(chain01, nchain01, natm01);

	/* asa of sidechain */
	/*
	asa4sidechain(chain01, nchain01);
	*/

	acc(chain01, nchain01, stdAsa, nstdAsa);
	accN(chain01, nchain01, stdAsaN, nstdAsaN);
	accC(chain01, nchain01, stdAsaC, nstdAsaC);
	/*
	asa01 = (double *)memalloc(nres01 * sizeof(double));
	acc01 = (double *)memalloc(nres01 * sizeof(double));
	getArray(chain01, nchain01, asa01, acc01);
	zasa01 = (double *)memalloc(nres01 * sizeof(double));
	calzscore(asa01, nres01, zasa01);
	free((void *)asa01);
	zacc01 = (double *)memalloc(nres01 * sizeof(double));
	calzscore(acc01, nres01, zacc01);
	free((void *)acc01);
	outputZscore(chain01, nchain01, zacc01, nres01);
	free((void *)zacc01);
	*/

	/* analysis of surface and interface start here */
	/* Note: The chains in complex must be included in either monomer1 or monomer2. */

	/* define and analyze surface */
	surf01 = (SURF *)memalloc(nres01 * sizeof(SURF));
	nsurf01 = defSurf(chain01, nchain01, &par, surf01);

	/* remove cavity residues from the residue group early assingned to surface */
	/* Check solvent vector, angle < 110^o */
	/*
	nsurf01 = checkCavity(chain01, nchain01, surf01, nsurf01);
	*/



	/* get more info related to different properties of amino acids to surface residues */
	/*
	rdPar4res(par4resFile, par4res);
	*/

	/* polarity of residue */
	/*
	pol2surf(surf01, nsurf01, par4res);
	*/


	/* hydrophility of residue */
	/*
	hydrophility2surf(surf01, nsurf01, par4res);
	*/



	/* zscore of solvent accessible surface area for each surface residue */
	asa01 = (double *)memalloc(nsurf01 * sizeof(double));
	surfAsa2array(surf01, nsurf01, asa01);
	zasa01 = (double *)memalloc(nsurf01 * sizeof(double));
	calzscore(asa01, nsurf01, zasa01);

	/* zscore of relative solvent accessibility for each surface residue */
	acc01 = (double *)memalloc(nsurf01 * sizeof(double));
	surfAcc2array(surf01, nsurf01, acc01);
	zacc01 = (double *)memalloc(nsurf01 * sizeof(double));
	calzscore(acc01, nsurf01, zacc01);


	/* define interface */
	/* Note: any interface residue must belong to surface. */
	dAccssChain01 = (CHAIN *)memalloc(nchain01 * sizeof(CHAIN));
	deltaAccssVer1(chain0, nchain0, chain01, nchain01, dAccssChain01);
	ninterf01 = defInterf(dAccssChain01, nchain01, surf01, nsurf01, &par);

	printf("#nchain4complex=%d nres4complex=%d natm4complex=%d chains4monomer1=%s chains4monomer2=%s\n"
		"#monomer1: nres=%d natm=%d nsurf=%d ninterf=%d\n",
		nchain0, nres0, natm0, argv[4], argv[5],
		nres01, natm01, nsurf01, ninterf01);

	freeChain(dAccssChain01, nchain01);

	property2surf(surf01, nsurf01);

	outputSurf(nres01, ninterf01, surf01, nsurf01, asa01, zasa01, acc01, zacc01, argv[11]);
	free((void *)surf01);
	free((void *)asa01);
	free((void *)zasa01);
	free((void *)acc01);
	free((void *)zacc01);


	/* analysis of surface and interface end here*/



	/*
	elestatic();

	vdw();

	hb();
	*/


	/* polar and apolar interaction */
	/*
	diffPolarApolar();
	*/


	freePar(&par);

	freeChain(chain0, nchain0);
	freeChain(chain01, nchain01);

	return 0;
}
