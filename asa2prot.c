#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "pdb.h"
#include "asa2pdb.h"
#include "max.h"
#include "memory.h"
#include "open_file.h"


static void	usage(char *prog)
{
	printf("Usage: %s atm2hybridFile hybrid2radFile standardASAfile pdblikeFile outputFile\n", prog);
}


static void	output(CHAIN *chain, int nchain, char *file)
{
	time_t	time0;
	char	time1[100];
	char	*chtmp;
	FILE	*fp;
	CHAIN	*chaini, *chaint;
	RES	*res;
	ATM	*atm;
	char	log_path[NAME_MAX];


	/* create all the non-existing parent directories specified by the path of file */
	if((chtmp = strrchr(file, '/')) != NULL)
	{
		memcpy(log_path, "mkdir -p ", strlen("mkdir -p ") + 1);
		strcat(log_path, file);
		*strrchr(log_path, '/') = '\0';
		system(log_path);
	}

	if((fp = open_file(file, "w")) == NULL)
	{
		printf("#Error: open %s error\n", file);
		exit(1);
	}

	time(&time0);
	ctime_r(&time0, time1);
	if(chtmp == NULL)
		chtmp = file;
	else	chtmp ++;
	fprintf(fp, "# %s was produced by asa2pep at %s\n", chtmp, time1);

	chaini = chain;
	chaint = chaini + nchain;

	/*
	for(; chain < chaint; chain ++)
	{
		for(res = chain->res; res != NULL; res = res->next)
		{
			fprintf(fp, "%c", res->name);
		}
	}
	fprintf(fp, "\n");
	*/

	for(chain = chaini; chain < chaint; chain ++)
	{
		/*
		fprintf(fp, "%-7c %8d\n", chain->id, chain->nres);
		*/
		for(res = chain->res; res != NULL; res = res->next)
		{
			/* only output the value for X in tripeptide (dipeptide),
			GLY-X-GLY (GLY-X) or ALA-X-ALA (ALA-X) */
			/*
			if(res->seq != 2)
				continue;

			fprintf(fp, "%c %8d %8.2f\n", res->name, res->natm, res->asa);
			*/

			/* only output the value for X in dipeptide, X-GLY or X-ALA */
			/*
			if(res->seq != 1)
				continue;

			fprintf(fp, "%c %8d %8.2f\n", res->name, res->natm, res->asa);
			*/

			/* output the values for all residues */
			fprintf(fp, "ASA %4d%c %c %c %8d %8.2f %8.3f\n",
				res->seq, res->iCode, chain->id, res->name, res->natm, res->asa, res->acc);

			/*
			for(atm = res->atm; atm != NULL; atm = atm->next)
			{
				fprintf(fp, "%-7s %8.2f %8.2f\n", atm->name, atm->r, atm->accss);
			}
			*/
		}
	}

	fclose(fp);
}


int	main(int argc, char *argv[])
{
	int	nchain;
	CHAIN	chain[CHAIN_MAX];
	HYBRID	hybrid[HYBRID_MAX];
	int	nhybrid;
	RESRAD	rad[N_RESRAD];
	int	nresrad;
	int	nres, natm;
	STDASA	stdAsa[50];
	int	nstdAsa;

int	i;

	if(argc != 6)
	{
		usage(argv[0]);
		exit(1);
	}

	printf("#Reading %s\n", argv[4]);
	if((nchain = pdbATOM(argv[4], chain)) <= 0)
	{
		printf("Error: number of chains: %d\n", nchain);
		exit(1);
	}


	nresrad = atm2hybrid(argv[1], rad);

	nchain = checkATOM(chain, nchain, rad, nresrad);

	nhybrid = hybrid2rad(argv[2], hybrid);

	getRad(hybrid, nhybrid, rad, nresrad);
	rad2pdb(chain, nchain, rad, nresrad);
	freeRad(rad, nresrad);

	natm = nele(chain, nchain, &nres);
	asa2pdb(chain, nchain, natm);

	/* calculate the relative accessibility for all residues */
	nstdAsa = rdStdAsa(argv[3], stdAsa);
	acc(chain, nchain, stdAsa, nstdAsa);

	/* correct the value of relative accessibility for N-terminal residue */
	nstdAsa = rdStdAsa("Xgly.ProtOr.dat", stdAsa);

	accN(chain, nchain, stdAsa, nstdAsa);
	
	/* correct the value of relative accessibility for C-terminal residue */
	nstdAsa = rdStdAsa("glyX.ProtOr.dat", stdAsa);

	accC(chain, nchain, stdAsa, nstdAsa);

	output(chain, nchain, argv[5]);
	freeChain(chain, nchain);

	return 0;
}
