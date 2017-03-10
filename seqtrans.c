#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "stdres.h"
#include "max.h"


static int 	ThreeToSingle(char *seq)
{
	int	i, ii, j;
	int	resnum;

	
	resnum = strlen(seq) / 3;

	if(resnum >= RES_MAX)
	{
		printf("Too many residues(>= %d) in peptide chain, please increse RES_MAX\n",
			RES_MAX);
		exit(1);	
	}

	for(i = ii = 0; ii < resnum; i += 3, ii ++)
	{
		for(j = 0; j < 23; j ++)
		{
			if(memcmp(seq + i, stdres[j].three, 3) == 0)
			{
				seq[ii] = stdres[j].single;
				break;
			}
		}
		if(j == 23)
		{
			seq[ii] = 'X';
		} /* assign X as the single-letter abbreviation for non-standard group */
	}
	seq[ii] = '\0';

	return 0;
}


static int	SingleToThree(char *seq)
{
}


int	SeqTrans(char *seq, int resflg)
{
	if(resflg == 3)
	{
		if(ThreeToSingle(seq) == 0)
			return 0;
		else	return 1;
	}
	else if(resflg == 1)
	{
		if(SingleToThree(seq) == 0)
			return 0;
		else	return 1;
	}
	else
	{
		printf("Please specify the type of abbreviation of the original sequence!\n");
		return 1;
	}
}


/* All lines of texts should be shorter than 80 characters.
   For details, please refer to FASTA format description */
int	WrtSeqFasta(char *seq, char *seqid, char *description, FILE *out)
{
	int	LengthOfLine;

	if(fprintf(out, ">%s, %s\n", seqid, description) >= 80)
	{
		printf("The first line of text is too long(>= 80 chars),\n");
		printf("please check the following line:\n");
		printf(">%s %s\n", seqid, description);
	}
	else if(strlen(seq) < 1)
	{
		printf("The length of sequence to write out is ZERO(residue)!\n");

		return -1;
	}

	for(LengthOfLine = 0; *seq != '\0';)
	{
		LengthOfLine += fprintf(out, "%c", *seq++);
		if(LengthOfLine % 60 == 0)
			fprintf(out, "\n");
	}
	fprintf(out, "\n");

	return 0;
}
