/* ckopen - open file; check for success */

/*
#include <stdinc.h>
*/
#include <stdio.h>
#include <stdlib.h>


/*
void *ckalloc(int amount);
FILE *ckopen(char *name, char *mode);
The lines above were added by Zhiqun Xie */

#include "cnctarea2.h"
/*This line was added by Zhiqun Xie */


FILE *ckopen(char *name, char *mode)
{
	FILE *fp;

	if ((fp = fopen(name, mode)) == NULL)
		printf("Cannot open %s.\n", name);
	return(fp);
}


/* ckalloc - allocate space; check for success */

void *ckalloc(int amount)
{
	void *p;

	if ((p = (void *) calloc( (unsigned) amount, 1)) == NULL)	{
		printf("Ran out of memory.\n");
                printf("There may be errors as follows:\n");
                printf("1) Not enough memory.\n");
                printf("2) The ARRAY may be overrode.\n");
                printf("3) The wild pointers.\n");
                exit(-1);
	}
	return(p);
}
