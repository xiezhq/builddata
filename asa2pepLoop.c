#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "open_file.h"
#include "max.h"


/*
#define CMD "asa2pep par/fmr74.resi-defs.mod-by-xie.dat par/fmr74.atom-defs.mod-by-xie.dat "
*/
#define CMD "asa2pep par/ProtOr.resi-defs.mod-by-xie.dat par/ProtOr.atom-defs.dat "


void	usage(char *prog)
{
	printf("Usage: %s pdbListFile\n", prog);
}


int	main(int argc, char *argv[])
{
	FILE	*fp;
	char	str[100];
	char	cmd[FILE_NAME_LEN];


	if(argc != 2)
	{
		usage(argv[0]);
		exit(1);
	}

	fp = open_file(argv[1], "r");

	for(fgets(str, 100, fp); !feof(fp); fgets(str, 100, fp))
	{
		memcpy(cmd, CMD, strlen(CMD));
		memcpy(cmd + strlen(CMD), str, strlen(str) + 1);

		/* asa2pep par/ProtOr.resi-defs.mod-by-xie.dat par/ProtOr.atom-defs.dat pdbfile */
		printf("# %s", cmd);
		system(cmd);
	}

	fclose(fp);
	return 0;
}
