#include <stdio.h>
#include <string.h>


static	void	usage(char *prog)
{
	printf("Usage: %s file\n", prog);
}


int	main(int argc, char *argv[])
{
	FILE	*fp;
	char	str[1000];


	if(argc != 2)
	{
		usage(argv[0]);
		exit(1);
	}

	if((fp = fopen(argv[1], "r")) == NULL)
	{
		printf("Error: fail to open %s\n", argv[1]);
		exit(1);
	}

	for(fgets(str, 1000, fp); !feof(fp); fgets(str, 1000, fp))
	{
		if(str[0] == '#') continue;

		printf("%s", str);
	}

	fclose(fp);

	return 0;
}
