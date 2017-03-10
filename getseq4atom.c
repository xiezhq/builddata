#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "max.h"

#include "open_file.c"



#define TRUE	1
#define FALSE	0


int	getseq4atom(char *pdbfile, char chainid, char *seq)
{
	FILE	*pdbfilep, *seqfilep;
	char	str[100];
	char	*seqp;
	char	lastres[6];
	char	*chrp;
	int	altlocflg;
	char	altloc;


	if(chainid == '_')
		chainid = ' ';

	pdbfilep = open_file(pdbfile, "r");
	
	for(fgets(str, 100, pdbfilep);
		memcmp("END   ", str, 6) != 0 && !feof(pdbfilep);
		fgets(str, 100, pdbfilep))
	{
		if(memcmp("ATOM   ", str, 6) == 0 && str[21] == chainid)
			break;
	}

	if(memcmp("ATOM   ", str, 6) != 0 || str[21] != chainid)
	{
		fclose(pdbfilep);

		seq[0] = '\0';

		return -1;
	}

	altlocflg = FALSE;
	altloc = ' ';
	lastres[0] = '\0';
	seqp = seq;

	for(;memcmp("TER   ", str, 6) != 0
		&&
		memcmp("END   ", str, 6) != 0
		&& !feof(pdbfilep);
		fgets(str, 100, pdbfilep))
	{
		/*
		if(str[21] != chainid || memcmp("ATOM  ", str, 6) != 0)
		{
			continue;
		}
		*/

		if(memcmp("ATOM  ", str, 6) != 0)
		{
			continue;
		}

		if(str[21] != chainid)
		{
			printf("%s", str);
			printf("Error! Found unknown chain #%c# inserted into the chain #%c#\n",
				str[21], chainid);

			exit(1);
		}


		/*
		if(chainid != ' ' && chainid != '_')
		{
			if(lastres[0] != '\0' && memcmp("TER   ", str, 6) == 0)
			{
				break;
			}

			if(str[21] != chainid || memcmp("ATOM  ", str, 6) != 0)
			{
				continue;
			}
		}
		else
		{
			if(memcmp("ATOM  ", str, 6) != 0)
			{
				continue;
			}
		}
		*/



		if(!isspace(str[16]) && altlocflg == FALSE)
		{
			altlocflg = TRUE;
			altloc = str[16];
		}
		if(!isspace(str[16]) && str[16] != altloc)
			continue;

		if(memcmp(lastres, str + 22, 5) != 0)
		{
			if(str[26] != ' ')
			{
				printf("%s", str);
			}

			memcpy(seqp, str + 17, 3);
			seqp += 3;
			memcpy(lastres, str + 22, 5);
		}
	}
	*seqp = '\0';

	fclose(pdbfilep);


	return 0;
}


void	usage(char *prg)
{
	printf("Usage: %s pdb_list_file pdb_seq4atm_file\n", prg);
	/*
	 pdb_list_file		list of pdb identifiers
	 pdb_seq4atm_file	output containing all sequences specified by pdb_list_file in FASTA format,
				the sequences are produced from ATOM record
}


int	main(int argc, char *argv[])
{
	char	pdblistfile[NAME_MAX], pdbseqfile[NAME_MAX];
	char	**pdblist;
	char	*seq;
	char	chainid;
	int	npdb;
	i;


	if(argc != 3)
	{
		usage(argv[0]);
		exit(1);
	}

	memcpy(pdblistfile, argv[1], strlen(argv[1])+1);

	pdblist = (char **)chrmatrix(PDB_MAX, 5);

	npdb = rd_pdblist(pdblistfile, pdblist);



	seq = (char *)memalloc(RES_MAX * sizeof(char));

	for(i = 0; i < npdb; i ++)
	{
		for(j = 0; j < nchain[i]; j ++)
		{
			getseq4atom(pdbfile, chainid, seq)
		}
	}

	free((void *)seq);

	memcpy(pdbseqfile, argv[2], strlen(argv[2])+1);


	return 0;
}
