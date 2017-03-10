#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "open_file.h"


FILE	*open_file(const char *name, char *mode)
{
	FILE	*fp;
	int	mode_selector;

	if( memcmp("r", mode, 1) == 0 )
		mode_selector = 0;
	else if( memcmp("w", mode, 1) == 0 )
		mode_selector = 1;
	else if( memcmp("a", mode, 1) == 0 )
		mode_selector = 2;
	else if( memcmp("r+", mode, 2) == 0 )
		mode_selector = 3;
	else if( memcmp("w+", mode, 2) ==0 )
		mode_selector = 4;
	else if( memcmp("a+", mode, 2) == 0 )
		mode_selector = 5;
	else
	{
		printf("mode to open %s error!\n", name);
		perror("fopen");
		/*
		return NULL;
		*/
		exit(1);
	}

	switch( mode_selector )
	{
		case 0:
			if( (fp = fopen(name, "r")) == NULL )
			{
				printf("open %s error!\n", name);
				perror("fopen");
				exit(1);
			}break;
		case 1:
			if( (fp = fopen(name, "w")) == NULL )
			{
				printf("open %s error!\n", name);
				perror("fopen");
				exit(1);
			} break;
		case 2:
			if( (fp = fopen(name, "a")) == NULL )
			{
				printf("open %s error!\n", name);
				perror("fopen");
				exit(1);
			} break;
		case 3:
			if( (fp = fopen(name, "r+")) == NULL )
			{
				printf("open %s error!\n", name);
				perror("fopen");
				exit(1);
			} break;
		case 4:
			if( (fp = fopen(name, "w+")) == NULL )
			{
				printf("open %s error!\n", name);
				perror("fopen");
				exit(1);
			} break;
		case 5:
			if( (fp = fopen(name, "a+")) == NULL )
			{
				printf("open %s error!\n", name);
				perror("fopen");
				exit(1);
			} break;
		default:
			{
				printf("open %s error (unknown fopen() mode %s)\n", name, mode);
				perror("fopen");
				exit(1);
			}
	}

	return fp;
}
