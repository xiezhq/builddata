/* written by Xie ZhiQun in July, 2004 */


/* printf(), fprintf(), perror(), FILE *, fopen(), fclose(), fflush(),
fgets(), feof()
*/
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
#include "searchpdb.h"

#include "pdb.h"
#include "max.h"


/* root directory where all output files should be placed 

for example: if we have /usr/local/database/pdb/data/structures/divided/pdb/??,
we can look /usr/local/database/pdb/data/structures/divided/ as the PDB root directory,
where, ?? is the two-letter directory holding classified PDB files.
If you like, then you can set any directory where you are allowed to write as your own
PDB root directory holding your data files derived from original PDB files. For example,
you can create a /home/xie/data/pdb/ directory as your own PDB root
directory where all two-letter name directories placed, so you will have many
directories as /home/xie/data/pdb/00, /home/xie/data/pdb/01, /home/xie/data/pdb/ab, etc.
*/
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
		"		All PDB coordinate files in <DIR> directory and its subdirectories are to be\n"
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
	int	itmp;


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

	itmp = strlen(suffix);
	memcpy(strstr(odir, ".ent"), suffix, itmp + 1);
	itmp = strlen(odir);
	memcpy(ofile, odir, itmp + 1);

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
					char	cmd[NAME_MAX] = {CMD};
					char	seqfile[NAME_MAX];

					if(strrchr(file, '.') == NULL)
					{
						printf("%s: unkown file type\n", file);

						return 0;
					}
					else if(strstr(file, ".ent") == NULL)
					{
						printf("%s: expect a PDB coordinate file with suffix "
							".ent or .ent.Z\n",
							file);

						return 0;
					}

					memcpy(seqfile, file, strlen(file) + 1);
					NameOfFile(seqfile, SUFFIX);

					/* full command:
					rd_pdb -i pdbfile seqfile
					*/
					/*
					strcat(cmd, file);
					strcat(cmd, " ");
					strcat(cmd, seqfile);
					*/

					/* full command:
					asa2pep par/ProtOr.resi.defs.mod-by-xie.dat par/ProtOr.atoms.defs.dat glyXgly.ProtOr.dat pdbfile outputfile
					*/
					sprintf(cmd + strlen(cmd), "%s %s", file, seqfile);

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
	the corresponding PDBlike coordinate file, and then 'system' to execute a pdb processing
	program to process the PDB coordinate file.

	If 'list_search' fails, it returns -1; otherwise it returns zero.

	Arguments:

	path	root directory containing all current PDB coordinate files divided
		according to two letter organization
	file	list file holding PDB chains to be processed
*/
static int	list_search(const char *path, const char *file)
{
	FILE	*fp;
	char	str[LINE_LEN_MAX];
	char	pdbfile[NAME_MAX], seqfile[NAME_MAX];
	int	i;
	int	intmp;
	char	cmd[NAME_MAX] = {CMD};


	if((fp = fopen(file, "r")) == NULL)
	{
		printf("open %s error in list_search()\n", file);

		fclose(fp);
		return -1;
	}

	for(fgets(str, LINE_LEN_MAX, fp); !feof(fp) && memcmp("END", str, 3) != 0; fgets(str, LINE_LEN_MAX, fp))
	{
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

		memcpy(pdbfile, path, strlen(path) + 1);

		/* use unix style path symbol, '/' */
		intmp = strlen(pdbfile);
		if(pdbfile[intmp - 1] != '/')
		{
			pdbfile[intmp] = '/';
			pdbfile[intmp+1] = '\0';
		}

		/* concatenate the name of two-letter subdirectory into path name */
		str2lower(str, 4);
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
		NameOfFile(seqfile, SUFFIX);

		/* full command:
		asa2pep par/ProtOr.resi.defs.mod-by-xie.dat par/ProtOr.atoms.defs.dat glyXgly.ProtOr.dat pdbfile outputfile
		*/
		sprintf(cmd + strlen(CMD), "%s %s", pdbfile, seqfile);

		printf("Launching command: %s\n", cmd);

		if(system(cmd) < 0)
		{
			printf("system(%s) in list_search() error\n", cmd);

			fclose(fp);
			return -1;
		}
	}

	fclose(fp);

	return 0;
}


/* chains2[0] = '\0' if monomer, else return actual chains in two parts in complex

*str		the 5th collumn in zlablist2005.txt
*chains1	chains in the first part of complex
*chains2	chains in the second part of complex
*/
static void	getChain(char *str, char *chains1, char *chains2)
{
	int	intmp;
	char	*chtmp;


	if((chtmp = strchr(str, ':')) == NULL)
	{
		if(isspace(*str))
		{
			chains1[0] = ' ';
			chains1[1] = '\0';
		}
		else
		{
			if((intmp = strlen(str)) > CHAIN_MAX)
			{
				printf("chains1=%s, too many chains, please increase CHAIN_MAX\n", str);
				exit(1);
			}
			else
			{
				for(; *str != '\n';)
				{
					if(isspace(*str))
					{
						str ++;
						continue;
					}

					*chains1++ = *str++;
				}
			}

			*chains1 = '\0';
		}
		chains2[0] = '\0';
	}
	else
	{
		*chtmp = '\0';
		if((intmp = strlen(str)) > CHAIN_MAX)
		{
			printf("chains1=%s, too many chains, please increase CHAIN_MAX\n", str);
			exit(1);
		}
		else
		{
			for(; *str != '\0';)
			{
				if(isspace(*str))
				{
					str ++;
					continue;
				}

				*chains1++ = *str++;
			}

			*chains1 = '\0';
		}

		chtmp ++;
		if((intmp = strlen(chtmp)) > CHAIN_MAX)
		{
			printf("chains2=%s, too many chains, please increase CHAIN_MAX\n", chtmp);
			exit(1);
		}
		else
		{
			for(; *chtmp != '\n';)
			{
				if(isspace(*chtmp))
				{
					chtmp ++;
					continue;
				}

				*chains2++ = *chtmp++;
			}

			*chains2 = '\0';
		}
	}
}


#undef	CMD
#define CMD	"interf "
#define PAR4CMD	" 0.05 1.0 par/ProtOr.resi-defs.mod-by-xie.dat par/ProtOr.atom-defs.dat glyXgly.ProtOr.dat"

static int	listSearch(const char *path, const char *file)
{
	FILE	*fp;
	char	str[LINE_LEN_MAX];
	char	pdbfile[NAME_MAX], seqfile[NAME_MAX];
	int	i;
	int	intmp;
	char	*chtmp;
	char	cmd[NAME_MAX] = {CMD};
	char	chains1[CHAIN_MAX], chains2[CHAIN_MAX];


	if((fp = fopen(file, "r")) == NULL)
	{
		printf("open %s error in listSearch()\n", file);

		fclose(fp);
		return -1;
	}

	for(fgets(str, LINE_LEN_MAX, fp); !feof(fp) && memcmp("END", str, 3) != 0; fgets(str, LINE_LEN_MAX, fp))
	{
		for(i = 0; i < 5 && isprint(str[i]); i++)
		{
			/*
			if(str[i] == ' ' || str[i] == '-' || str[i] == '\n')
			*/
			if(isalnum(str[i]) == 0)
				break;
		}
		if(i < 4)	continue;

		getChain(str + 5, chains1, chains2);

		/* process only complexes in zlablist2005.txt this time */
		if(chains2[0] == '\0')
			continue;

		str[i] = '\0';
		
		printf("Processing %s entry in %s\n", str, file);

		memcpy(pdbfile, path, strlen(path) + 1);

		/* use unix style path symbol, '/' */
		intmp = strlen(pdbfile);
		if(pdbfile[intmp - 1] != '/')
		{
			pdbfile[intmp] = '/';
			pdbfile[intmp+1] = '\0';
		}

		/* concatenate the name of two-letter subdirectory into path name */
		str2lower(str, 4);
		strncat(pdbfile, str + 1, 2);

		/* get the name of pdbfile, path/pdbxxxx.ent */
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
		NameOfFile(seqfile, SUFFIX);

		/* full command:
		interf pdbfile nchain01 nchain02 chains01 chains02 cutoff4surf cutoff4interf \
		par/ProtOr.resi-defs.mod-by-xie.dat par/ProtOr.atom-defs.dat glyXgly.ProtOr.dat
		*/
		sprintf(cmd + strlen(CMD), "%s %d %d %s %s %s %s",
			pdbfile, strlen(chains1), strlen(chains2), chains1, chains2, PAR4CMD, seqfile);

		printf("Launching command: %s\n", cmd);

		if(system(cmd) < 0)
		{
			printf("system(%s) in listSearch() error\n", cmd);

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
			/*
			if(list_search(option->dir, option->list))
			{
				fprintf(stderr, "list_search() error in main()\n");
				exit(1);
			}
			*/
			if(listSearch(option->dir, option->list))
			{
				fprintf(stderr, "listSearch() error in main()\n");
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
