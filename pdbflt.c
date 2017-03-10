/* written by Xie ZhiQun in July, 2004 */


/* printf(), fprintf(), perror(), FILE *, fopen(), fclose(), fflush(),
	fgets(), feof() */
#include <stdio.h>

/* malloc(), free(), exit() and system() */
#include <stdlib.h>

/* getopt() */
#include <unistd.h>

/* ftw() */
#include <ftw.h>
#include <sys/types.h>
#include <sys/stat.h>

/* isprint(), isalnum() */
#include <ctype.h>

/* memcmp(), strcat(), strlen() */
#include <string.h>

/* usage(), help(), dir_search(), list_search() and many macros */
#include "pdbflt.h"

/* root directory where all output files should be put */
static char	odir[NAME_MAX];


static void	usage(const char *prog)
{
	fprintf(stderr,
		"Usage: %s [-h] <-d pdb_directory> <-o output_directory> [-l PDB_chain_list_file]\n",
		prog);
	fflush(stderr);
}


static void	help(const char *prog)
{
	printf("Usage: %s [options]\n"
		"Options:\n"
		"-d <DIR>	PDB directory holding all current PDB coordinate files according to the two\n"
		"		letter organization\n"
		"		All PDB coordinate files in <DIR> directory and its subdirectories are  to be\n"
		"		processed if option l not specified\n"
		"-o <ODIR>	Output directory holding all output files\n"
		"-l <LIST>	Process only the specified PDB chains which match the entries listed in the\n"
		"		list file <LIST>\n"
		"		One entry for each line in <LIST>\n"
		"-h		Show this message and exit\n",
		prog);
}


/* NameOfFile(char *ofile, char *suffix):

	'NameOfFile' constructs the name of 'ofile' based on the two-letter organization of
	PDB directory structure. It receives a pointer to a null-terminated character string
	containing the name of a '???file' with suffix '???'. When 'NameOfFile' returns,
	it translates the name of the '???file' into the name of the corresponding 'ofile' with
	suffix 'dom', and places the name in the same character string pointed by the pointer
	'ofile'.

	Arguments:
	ofile		Character pointer to '???file' with suffix '???' before the first call to
			'NameOfFile'. The data stored in the same space are turned into the name
			of the 'ofile' with suffix 'dom' when 'NameOfFile' returns.
	suffix		Character pointer to target suffix to be translated. The corresponding
			root directory is also to be translated such as PDB, ATM, DSP and DOM
			directories etc, the name of which is same as the suffix of 'ofile'.

	Note:
	The parent directory of the ATM directory holding all 'ofile's should be the same directory
	as that of the PDB directory. In fact, PDB, ATM, DSP and DOM directories are all in the same
	directory.
*/
static void	NameOfFile(char *ofile, char *suffix)
{
	char	*p;
	int	i;
	extern char	odir[NAME_MAX];


	p = ofile;

	if((i = strlen(odir) - 1) > NAME_MAX)
	{
		printf("pathname (%s) is too long(%d), please increase NAME_MAX\n", odir, i);
		exit(1);
	}
	if(odir[i] != '/')
	{
		odir[i + 1] = '/';
		odir[i + 2] = '\0';
	}
	if(i + strlen(ofile) > NAME_MAX)
	{
		printf("pathname (%s) is too long, please increase NAME_MAX\n", ofile);
		exit(1);
	}
	if((p = strrchr(ofile, '/')) == NULL)
	{
		strcat(odir, ofile);
	}
	else	strcat(odir, p - 2); /* pdb/divided/xx/ directory, not pdb/all/ directory */

	memcpy(strstr(odir, ".ent"), suffix, strlen(suffix) + 1);
	memcpy(ofile, odir, strlen(odir) + 1);

	/* restore the odir */
	odir[i + 1] = '\0';
}


static int	pdbflt(const char *file, const struct stat *statp, int type)
{
	switch(type)
	{
		case FTW_F:
		{
			switch(statp->st_mode & S_IFMT)
			{
				case S_IFREG:
				{
					char	cmd[NAME_MAX] = {COM};
					char	seqfile[NAME_MAX];

					if(strrchr(file, '.') == NULL)
					{
						printf("%s: unkown file type\n", file);

						return 0;
					}
					else if(strstr(file, ".ent") == NULL)
					{
						printf("%s: expect a PDB coordinate file with suffix "
							".ent\n",
							file);

						return 0;
					}

					memcpy(seqfile, file, strlen(file) + 1);
					NameOfFile(seqfile, ".seq");

					/* full cmd: rd_pdb -i pdbfile seqfile */
					strcat(cmd, file);
					strcat(cmd, " ");
					strcat(cmd, seqfile);

					printf("Launching cmd: %s\n", cmd);

					if(system(cmd) < 0)
					{
						printf("system(%s) error in pdbflt()\n", cmd);
						perror("system()");
						return -1;
					}
				}
					break;
				default:
					printf("%s isn't a regular file\n", file);
					break;
			}
		}
			break;
		case FTW_D:
			printf("%s is a directory\n", file);
			break;
		case FTW_DNR:
			printf("can't read %s\n", file);
			break;
		case FTW_NS:
			printf("stat() error for %s\n", file);
			break;
		default:
			printf("pdbflt() error: unknown type for %s\n", file);
			break;
	}

	return 0;
}


/* dir_search(const char *path):

	'dir_search' recursively descends the directory hierarchy rooted in 'path'.
	For each object in the hierachy, 'dir_search' calls user-defined function
	'pdbflt', passing it a pointer to a null-terminated character string containing
	the name of the object, a pointer to a 'stat' structure containing information
	about the object, and an integer (refer to ftw(3C) in Irix Unix).

	If 'dir_search' fails, it returns -1; otherwise it returns zero.

	argument:

	path	root directory containning all current PDB coordinate files divided
		according to two letter organization
*/
static int	dir_search(const char *path)
{
	int	depth = PDB_DIR_DEPTH;

	if(ftw(path, pdbflt, depth))
	{
		printf("%s: ftw() error in dir_search()\n", path);
		perror("ftw()");
		return -1;
	}

	return 0;
}


/* list_search(const char *path, const char *file):

	'list_search' recursively searches the directory hierarchy rooted in 'path' for
	PDB coordinate files based on the PDB code[chainidentifier] entries listed in 'file'.
	For each PDB coordinate file, 'list_search' calls 'NameOfFile' to get the name of
	the corresponding PDBlike coordinate file, and then 'system' to execute a pdb filting
	program (rd_pdb) to translate a PDB coordinate file into a PDBlike coordinate file.

	If 'list_search' fails, it returns -1; otherwise it returns zero.

	Arguments:

	path	root directory containing all current PDB coordinate files divided
		according to two letter organization
	file	List file holding PDB chains to be processed
*/
static int	list_search(const char *path, const char *file)
{
	FILE	*fp;
	char	str[LINE_LEN_MAX];
	char	pdbfile[NAME_MAX], seqfile[NAME_MAX];
	int	i;


	if((fp = open_file(file, "r")) == NULL)
	{
		printf("open %s error in list_search()\n", file);

		fclose(fp);
		return -1;
	}

	for(fgets(str, LINE_LEN_MAX, fp); !feof(fp) && memcmp("END", str, 3) != 0; fgets(str, LINE_LEN_MAX, fp))
	{
		char	command[NAME_MAX] = {COM};

		for(i = 0; (i < 5) && (isprint(str[i]) != 0); i++)
		{
			/*
			if(str[i] == ' ' || str[i] == '-' || str[i] == '\n')
			*/
			if(isalnum(str[i]) == 0)
				break;
		}
		if(i < 4)	continue;

		str[i] = '\0';
		
		printf("Processing %s entry in %s\n", str, file);

		/* build command: rd_pdb -i -c chainid pdbfile seqfile */
		if(i == 5)
		{
			strcat(command, "-c ");
			strncat(command, str + 4, 1);
			strcat(command, " ");
		}

		memcpy(pdbfile, path, strlen(path) + 1);

		/* use unix style path symbol, '/' */
		if(pdbfile[strlen(pdbfile) - 1] != '/')
		{
			strcat(pdbfile, "/");
		}

		/* concatenate the name of two-letter subdirectory into path name */
		strncat(pdbfile, str + 1, 2);

		/* get the name of pdbfile, path/pdbxxxx?.ent */
		strcat(pdbfile, "/pdb");
		strncat(pdbfile, str, 4);

		/* for normal PDB data */
		/*
		strcat(pdbfile, ".ent");
		*/

		/* for compression format of PDB data */
		strcat(pdbfile, ".ent.Z");

		/* get the name of seqfile */
		memcpy(seqfile, pdbfile, strlen(pdbfile) + 1);
		NameOfFile(seqfile, ".seq");

		/* concatenate filename into command name */
		strcat(strcat(strcat(command, pdbfile), " "), seqfile);

		printf("Launching command: %s\n", command);

		if(system(command) < 0)
		{
			printf("system(%s) in list_search() error\n", command);

			fclose(fp);
			return -1;
		}
	}

	fclose(fp);

	return 0;
}


int	main(int argc, char *argv[])
{
	char	*optstring = "d:o:hl:";
	int	opt;
	int	dflg = 0,
		lflg = 0,
		oflg = 0;
	OPTION	*option;

	extern char	odir[NAME_MAX];


	if(argc > 7 || argc < 2)
	{
		usage(argv[0]);
		return -1;
	}

	if((option = (OPTION *)malloc(sizeof(OPTION))) == NULL)
	{
		printf("malloc() in main() in %s fails\n", argv[0]);
		perror("malloc()");
		exit(1);
	}

	/* parse command line options */
	while((opt = getopt(argc, argv, optstring)) != -1)
	{
		switch(opt)
		{
			case 'd':
				option->dir = optarg;
				dflg ++;
				break;
			case 'o':
				memcpy(odir, optarg, strlen(optarg) + 1);
				oflg ++;
				break;
			case 'l':
				option->list = optarg;
				lflg ++;
				break;
			case 'h':
				help(argv[0]);
				exit(1);
			default:
				usage(argv[0]);
				exit(1);
		}
	}

	if(dflg && oflg)
	{
		if(!lflg)
		{
			if(dir_search(option->dir))
			{
				fprintf(stderr, "dir_search() error in main()\n");
				exit(1);
			}
		}
		else if(lflg)
		{
			if(list_search(option->dir, option->list))
			{
				fprintf(stderr, "list_search() error in main()\n");
				exit(1);
			}
		}
	}
	else
	{
		usage(argv[0]);
		exit(1);
	}

	free(option);

	return 0;
}
