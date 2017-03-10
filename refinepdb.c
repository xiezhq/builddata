#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include "max.h"
#include "open_file.c"
#include "memoryxie.c"


#define TRUE	1
#define FALSE	0
#define MARKER	"------"
#define EPSILON	0.000001
#define PDB_MAX	50000 /* max of PDB entries to process, refer to statistics on PDB website
 and the index file of PDB entries (for example, resolu.idx) */



typedef struct chainseq
{
	int	len; /* length of the sequence */
	char	id; /* chain identifier */
	char	*seq; /* sequence of the chain */

	struct chainseq	*next; /* next chain */
} CHAINSEQ;


typedef struct pdbseq
{
	char	id[5]; /* identifier of pdb entry */
	int	nchain; /* number of chains in the pdb entry */

	CHAINSEQ	*chain; /* the first chain */
	struct pdbseq	*next;
} PDBSEQ;



/* 
 read resolu.idx from PDB site and check the resolution of structures

 return the list of structures with the resolution better than (or equal to) res_cutoff */
int	rd_resolu(char *file, float res_cutoff, char **pdbid)
{
	FILE	*fp;
	char	str[100];
	int	npdb;
	float	res;


	fp = open_file(file, "r");

	for(fgets(str, 100, fp); memcmp(str, MARKER, 6) != 0; fgets(str, 100, fp))
	{
	}

	npdb = 0;
	for(fgets(str, 100, fp); !feof(fp); fgets(str, 100, fp))
	{
		res = atof(str+7);

		if((res < res_cutoff || fabs(res - res_cutoff) < EPSILON) && res > 0.0)
		{
			memcpy(pdbid[npdb], str, 4);

			pdbid[npdb][0] = tolower(pdbid[npdb][0]);
			pdbid[npdb][1] = tolower(pdbid[npdb][1]);
			pdbid[npdb][2] = tolower(pdbid[npdb][2]);
			pdbid[npdb][3] = tolower(pdbid[npdb][3]);

			npdb ++;
		}
	}


	fclose(fp);

	return npdb;
}


/*
 read rd_pdb_seqres.txt from PDB site and check the length of every chain in pdb file

 return the list of the wanted pdbid */
int	rd_pdb_seqres(char *file, int nres_cutoff, int *nchain, char **pdbid)
{
	FILE	*fp;
	char	str[100];
	int	i, j;
	int	flg;
	char	lastpdb[5];


	fp = open_file(file, "r");

	i = -1;
	lastpdb[0] = '\0';

	for(fgets(str, 100, fp); !feof(fp); fgets(str, 100, fp))
	{
		if(str[0] != '>')
		{
			continue;
		}

		if(memcmp(lastpdb, str+1, 4) != 0)
		{
			i ++;
			flg = TRUE; 
			memcpy(pdbid[i], str+1, 4);
			nchain[i] = 1;

			/*
			pdb[i].nchain = 1;
			pdb[i].chain = (struct chainseq *)memalloc(sizeof(struct chainseq));
			*/
		}
		else if(flg == TRUE)
		{
			/*
			if(flg == FALSE)
				i ++;
			*/

			nchain[i] += 1;

			/*
			pdb[i].nchain ++;
			pdb[i].chain.next = (struct chainseq *)memalloc(sizeof(struct chainseq));
			*/
		}

		/*
		pdb[i].chain.id = str[6];
		*/

		/* all chains must be the peptide chain with at least nres_cutoff amino acids */
		if(flg == TRUE
		&& (memcmp(str+12, "protein ", 8) != 0 || atoi(strstr(str+12, "length:")+7) < nres_cutoff))
		{
			flg = FALSE;
			i --;
		}

		memcpy(lastpdb, str+1, 4);
	}
	fclose(fp);


	return (i+1);
}


/* return the list of the coincided pdb, pdbid and number of chains are stored in **pdbid1
 and *nchain array, respectively.
 pdbid1[i] and nchain[i] must refer to the same pdbid */
int	cmp_pdblist(char **pdbid1, int *nchain, int n1, char **pdbid2, int n2)
{
	int	i, j;


	for(i = 0; i < n1; i ++)
	{
		for(j = 0; j < n2; j ++)
		{
			if(memcmp(pdbid1[i], pdbid2[j], 4) == 0)
			{
				break;
			}
		}

		if(j == n2)
		{
			memcpy(pdbid1[i], pdbid1[n1-1], 4);
			nchain[i] = nchain[n1-1];
			n1 --;
			i --;
		}
	}

	return n1;
}


/* return the list of multimers with at most nchain_cutoff chains */
int	multimermostlist(char **pdbid, int *nchain, int npdb, int nchain_cutoff)
{
	int	i;


	for(i = 0; i < npdb; i ++)
	{
		if(nchain[i] > nchain_cutoff)
		{
			memcpy(pdbid[i], pdbid[npdb-1], 4);
			nchain[i] = nchain[npdb-1];

			i --;
			npdb --;
		}
	}

	return npdb;
}


/* return the list of multimers with the nchain_cutoff chains */
int	multimerequalist(char **pdbid, int *nchain, int npdb, int nchain_cutoff)
{
	int	i;


	for(i = 0; i < npdb; i ++)
	{
		if(nchain[i] != nchain_cutoff)
		{
			memcpy(pdbid[i], pdbid[npdb-1], 4);
			nchain[i] = nchain[npdb-1];

			i --;
			npdb --;
		}
	}

	return npdb;
}


/* return the list of multimers with at least nchain_cutoff chains */
int	multimerleastlist(char **pdbid, int *nchain, int npdb, int nchain_cutoff)
{
	int	i;


	for(i = 0; i < npdb; i ++)
	{
		if(nchain[i] < nchain_cutoff)
		{
			memcpy(pdbid[i], pdbid[npdb-1], 4);
			nchain[i] = nchain[npdb-1];

			i --;
			npdb --;
		}
	}

	return npdb;
}


/* return the pdb entries containing multiple different chains and the number of such pdb entries */
int	checkseqidentity(char *seqfile, char **pdbid, int *nchain, int npdb, PDBSEQ *pdb)
{
	FILE	*fp;
	char	str[100];
	int	i, j;
	char	lastpdbid[10];
	char	pdbskip[10];
	PDBSEQ	*pdbj;
	struct chainseq	*chain, *chain1, *chain2;
	int	pdbflg;


	fp = open_file(seqfile, "r");

	lastpdbid[0] = '\0';
	pdbskip[0] = '\0';
	j = 0;
	pdbflg = FALSE;
	pdbj = NULL;

	for(fgets(str, 100, fp); !feof(fp); fgets(str, 100, fp))
	{
		if(str[0] == '>')
		{
			if(memcmp(lastpdbid, str+1, 4) != 0 && memcmp(pdbskip, str+1, 4) != 0)
			{
				for(i = 0; i < npdb; i ++)
				{
					if(memcmp(str+1, pdbid[i], 4) == 0)
					{
						memcpy(lastpdbid, pdbid[i], 4);

						pdbj = pdb+j;

						memcpy(pdbj->id, lastpdbid, 4);
						pdbj->nchain = nchain[i];

						pdbj->chain = (struct chainseq *)memalloc(sizeof(struct chainseq));
						chain = pdbj->chain;
						chain->next = NULL;

						chain->id = str[6];
						chain->len = atoi(strstr(str+12, "length:")+7);
						chain->seq = (char *)memalloc(chain->len * sizeof(char) + 1);
						*(chain->seq) = '\0';

						j ++;

						break;
					}
				}

				if(i < npdb)
					pdbflg = TRUE;
				else	
				{
					pdbflg = FALSE;
					memcpy(pdbskip, str+1, 4);
				}
			}
			else if(pdbflg == TRUE)
			{
				pdbj->chain->next = (struct chainseq *)memalloc(sizeof(struct chainseq));

				chain = pdbj->chain->next;
				chain->next = NULL;

				chain->id = str[6];
				chain->len = atoi(strstr(str+12, "length:")+7);
				chain->seq = (char *)memalloc(chain->len * sizeof(char) + 2);
				*(chain->seq) = '\0';
			}

			continue;
		}

		if(pdbflg == FALSE)
			continue;

		snprintf(strchr(chain->seq, '\0'), strlen(str), "%s", str);
		/* refer to the manual of snprintf */

		/*
		printf("%.4s %c %d\n", pdbj->id, chain->id, chain->len);
		*/
	}

	fclose(fp);


	/* check the identity of the sequences in the same pdb entry */
	for(i = 0; i < npdb; i ++)
	{
		for(chain1 = (pdb+i)->chain; chain1 != NULL; chain1 = chain->next)
		{
			for(chain2 = chain1->next; chain2 != NULL; chain2 = chain2->next)
			{
				if(chain1->len != chain2->len)
				{
					break;
				}
				else
				{
					if(memcmp(chain1->seq, chain2->seq, chain1->len) != 0)
					{
						break;
					}
				}
			}

			if(chain2 != NULL)
				break;
		}

		/* remove pdb entry, all chains of which are identical */
		if(chain2 == NULL && (pdb+i)->nchain != 1)
		{
			(pdb+i)->nchain = 0;
		}
	}

	return npdb;
}


void	output_pdbid(char **pdbid, int *nchain, int npdb)
{
	int	i;

	printf("pdbid nchain\n");
	for(i = 0; i < npdb; i ++)
	{
		printf("%.4s %5d\n", pdbid[i], nchain[i]);
	}
}


int	divmer(PDBSEQ *pdb, int npdb, PDBSEQ *initpdbs, PDBSEQ *initpdbm)
{
	int	i, j;
	int	npdbs; /* number of single chain PDB entries */
	PDBSEQ	*lastpdbs, *lastpdbm;


	npdbs = 0;
	initpdbs = NULL;
	initpdbm = NULL;

	for(i = 0; i < npdb; i ++)
	{
		if((pdb+i)->nchain == 1)
		{
			npdbs ++;

			if(initpdbs == NULL)
				initpdbs = pdb + i;
			else	lastpdbs->next = pdb + i;

			lastpdbs = pdb + i;
			lastpdbs->next = NULL;
		}
		else
		{
			if(initpdbm == NULL)
				initpdbm = pdb + i;
			else	lastpdbm->next = pdb + i;

			lastpdbm = pdb + i;
			lastpdbm->next = NULL;
		}
	}

	return npdbs;
}


void	output_pdb(PDBSEQ *pdb, int npdb)
{
	int	i;
	CHAINSEQ	*chain;


	printf("%4s %3s\n", "pdb", "nchain");
	for(i = 0; i < npdb; i ++)
	{
		if((pdb+i)->nchain == 0)
			continue;

		printf("%.4s %3d\n", (pdb+i)->id, (pdb+i)->nchain);
		for(chain = (pdb+i)->chain; chain != NULL; chain = chain->next)
		{
			printf("%c %-5d", chain->id, chain->len);
		}
		printf("\n");
	}
}


void	free_pdb(PDBSEQ *pdb, int npdb)
{
	int		i;
	CHAINSEQ	*chain, *nextchain;


	for(i = 0; i < npdb; i ++)
	{
		for(chain = (pdb+i)->chain; chain != NULL; chain = nextchain)
		{
			nextchain = chain->next;

			free((void *)chain->seq);
			free((void *)chain);
		}
	}

	free((void *)pdb);
}


void	usage(char *prg)
{
	printf("Usage: %s resolu_idx_file resolution_cutoff pdb_seqres_file chainlength_cutoff\n", prg);
	/*
	resolu_idx_file:	resolu.idx file from PDB 
	resolution_cutoff:	the resolution of PDB entry must be better than resolution_cutoff
	pdb_seqres_file:	pdb_seqres.txt from PDB
	chainlength_cutoff:	every chain in a PDB entry must have at least seqlength_cutoff amino acids
	*/
}


int	main(int argc, char *argv[])
{
	char	**pdbid1, **pdbid2;
	int	*nchain;
	int	npdb1, npdb2, npdb;
	char	resolufile[NAME_MAX], pdbseqresfile[NAME_MAX];
	float	resolu_cutoff;
	int	seqlength_cutoff;

	int	nchain_cutoff; /* cutoff of number of chains in pdb structure entry */

	PDBSEQ	*pdb, *initpdbs, *initpdbm;
	int	npdbs; /* number of single-chain PDB entries */
	int	npdbsm; /* number of multimers, each chain of which has a PDB entry */


	if(argc != 5)
	{
		usage(argv[0]);
		exit(1);
	}

	memcpy(resolufile, argv[1], strlen(argv[1])+1);
	resolu_cutoff = atof(argv[2]);
	memcpy(pdbseqresfile, argv[3], strlen(argv[3])+1);
	seqlength_cutoff = atoi(argv[4]);

	pdbid1 = chrmatrix(PDB_MAX, 5);

	npdb1 = rd_resolu(resolufile, resolu_cutoff, pdbid1);

	pdbid2 = chrmatrix(PDB_MAX, 5);
	nchain = (int *)memalloc(PDB_MAX * sizeof(int));

	npdb2 = rd_pdb_seqres(pdbseqresfile, seqlength_cutoff, nchain, pdbid2);

	/* get refined list of pdbid, returned by pdbid2 array */
	npdb = cmp_pdblist(pdbid2, nchain, npdb2, pdbid1, npdb1);


	free_chrmatrix(pdbid1);



	/* dimer and monomer */
	nchain_cutoff = 2;
	npdb = multimermostlist(pdbid2, nchain, npdb, nchain_cutoff);


	/* dimer */
	/*
	nchain_cutoff = 2;
	npdb = multimerequalist(pdbid2, nchain, npdb, nchain_cutoff);
	*/

	/* monomer */
	/*
	nchain_cutoff = 1;
	npdb = multimerequalist(pdbid2, nchain, npdb, nchain_cutoff);
	*/


	/*
	output_pdbid(pdbid2, nchain, npdb);
	*/

	pdb = (PDBSEQ *)memalloc(npdb * sizeof(PDBSEQ));

	npdb = checkseqidentity(pdbseqresfile, pdbid2, nchain, npdb, pdb);

	free((void *)nchain);
	free_chrmatrix(pdbid2);


	/* return monomers and dimers stored in pdb and pdbm array, respectively */
	npdbs = divmer(pdb, npdb, initpdbs, initpdbm)

	npdbsm = multimer_cmp_monomer(pdb, initpdbs, initpdbm);

	output_pdb(pdb, npdb);

	free_pdb(pdb, npdb);



	return 0;
}
