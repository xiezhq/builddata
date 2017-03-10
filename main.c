#include <stdio.h>
#include <stdlib.h>
#include "pdb.h"
#include "memory.h"
#include "max.h"



static void	usage(char *prog)
{
	printf("Usage: %s pdbFile\n", prog);
}


int	main(int argc, char *argv[])
{
	CHAIN	*chain, *chain0, *chaint;
	RES	*res;
	ATM	*atm;
	int	nchain;
	char	*chp;


	if(argc != 2)
	{
		usage(argv[0]);
		exit(1);
	}


	chain0 = (CHAIN *)memalloc(CHAIN_MAX * sizeof(CHAIN));
	

	if((nchain = pdbATOM(argv[1], chain0)) == 0)
	{
		printf("pdbATOM() processes %s error\n", argv[1]);;
		exit(1);
	}
	else	printf("nchain = %d\n", nchain);


	/* test */
	/*
	chaint = chain0 + nchain;
	for(chain = chain0; chain < chaint; chain ++)
	{
		printf("%c %4d\n", chain->id, chain->nres);

		for(res = chain->res; res != NULL; res = res->next)
		{
			printf("%c %c %.3s %4d%c %d\n",
				chain->id,
				res->name, chp = s2tRes(res->name), res->seq, res->iCode, res->natm);
			free((void *)chp);

			for(atm = res->atm; atm != NULL; atm = atm->next)
			{
				printf("%6d %.3s%c %8.3f %8.3f %8.3f\n",
					0, atm->name, atm->altLoc,
					atm->x, atm->y, atm->z);
			}
		}
	}
	*/
	/* test */


	free((void*) chain0);


	return 0;
}
