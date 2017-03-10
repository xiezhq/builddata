#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "pdb.h"
#include "asa2pdb.h"
#include "max.h"


static void	usage(char *prog)
{
	printf("Usage: %s atm2hybridFile hybrid2radFile pdblikeFile\n", prog);
}


static void	output(CHAIN *chain, int nchain)
{
	CHAIN	*chaini, *chaint;
	RES	*res;
	ATM	*atm;
	char	*name;


	chaini = chain;
	chaint = chaini + nchain;

	/*
	for(; chain < chaint; chain ++)
	{
		for(res = chain->res; res != NULL; res = res->next)
		{
			printf("%c", res->name);
		}
	}
	printf("\n");
	*/

	for(chain = chaini; chain < chaint; chain ++)
	{
		/*
		printf("%-7c %8d\n", chain->id, chain->nres);
		*/
		for(res = chain->res; res != NULL; res = res->next)
		{
			/* only output the value for X in tripeptide (dipeptide),
			GLY-X-GLY (GLY-X) or ALA-X-ALA (ALA-X) */
			/*
			if(res->seq != 2)
				continue;

			printf("%3s %3d %8.2f\n", name = s2tRes(res->name), res->natm, res->asa);
			*/

			/* only output the value for X in dipeptide, X-GLY or X-ALA */
			if(res->seq != 1)
				continue;

			printf("%3s %3d %8.2f\n", name = s2tRes(res->name), res->natm, res->asa);


			/*
			printf("ASA %4d%c %c %c %8d %8.2f %8.3f\n",
				res->seq, res->iCode, chain->id, res->name, res->natm, res->asa, res->acc);
			*/

			free((void *)name);

			/*
			for(atm = res->atm; atm != NULL; atm = atm->next)
			{
				printf("%-7s %8.2f %8.2f\n", atm->name, atm->r, atm->accss);
			}
			*/
		}
	}
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


	if(argc != 4)
	{
		usage(argv[0]);
		exit(1);
	}

	if((nchain = pdbATOM(argv[3], chain)) <= 0)
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

	/* asa of sidechains for each residue */
	asa4sidechain(chain, nchain);
	/* asa of sidechains for each residue */

	output(chain, nchain);
	freeChain(chain, nchain);

	return 0;
}
