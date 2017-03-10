#include <stdio.h>
#include "property4res.h"

/* for definition of struct CHAIN */
#include "pdb.h"

/* for definition of struct SURF */
#include "interf.h"


void	property2chain(CHAIN *chain, int nchain)
{
	extern struct stdresidue	stdres[];

	CHAIN	*chaint;
	RES	*res;
	int	i;


	for(chaint = chain + nchain; chain < chaint && chain != NULL; chain ++)
	{
		for(res = chain->res; res != NULL; res = res->next)
		{
			for(i = 0; stdres[i].single != 'X'; i ++)
			{
				if(res->name == stdres[i].single)
				{
					res->pol = stdres[i].pol;
					res->philic = stdres[i].philic;
					res->e = stdres[i].e;
					break;
				}
			}

			if(stdres[i].single == 'X')
			{
				printf("#Warning: no value for unknown residue (%c %d%c %c), "
					"default values were assigned to it\n",
					chain->id, res->seq, res->iCode, res->name);

				res->pol = stdres[i].pol;
				res->philic = stdres[i].philic;
				res->e = stdres[i].e;
			} 
		}
	}
}

void	property2surf(SURF *surf, int nSurf)
{
	extern struct stdresidue	stdres[];

	SURF	*surft;
	int	i;


	for(surft = surf + nSurf; surf < surft && surf != NULL; surf ++)
	{
		for(i = 0; stdres[i].single != 'X'; i ++)
		{
			if(surf->name == stdres[i].single)
			{
				surf->pol = stdres[i].pol;
				surf->philic = stdres[i].philic;
				surf->e = stdres[i].e;
				break;
			}
		}

		if(stdres[i].single == 'X')
		{
			printf("#Warning: no value for unknown residue (%c %d%c %c), "
				"default values were assigned to it\n",
				surf->chainid, surf->seq, surf->iCode, surf->name);

			surf->pol = stdres[i].pol;
			surf->philic = stdres[i].philic;
			surf->e = stdres[i].e;
		}
	}
}
