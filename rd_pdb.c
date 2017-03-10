/*

The program was written by Xie Zhi-Qun in 1997, and some modifications have been made
since that time.

It attempts to retrieve the required informations from the primary structure section(SEQRES),
secondary structure section(HELIX, SHEET, TURN), connectivity annotation section(SSBOND),
crystallographic and coordinate transformation section(CRYST1, ORIGXn, SCALEn, MTRIXn,
TVECT), coordinate section(ATOM) in a PDB file, respectively.
Refer to http://www.rcsb.org/pdb/ for the PDB File Format Contents Guide Version 2.2 (20
December 1996).


Note:
Current version does only include the codes to process the peptide chain(s), and any information
about DNA/RNA will be simply striped off.


Revision history for this file (A revision history record was added here on June 27, 2004):
 
	Revision 4, July 19, 2004:

	4.1	Modified rd_chain() to reserve the amino acid residues with incomplete backbone
		(the original rd_chain() always removes the incomplete residues)
	

	Revision 3, July 18, 2004:

	3.1	Removed strcpy() and strncmp(), used memcpy() or setting value and memcmp()
		instead, respectively, in order to improve performance

	3.2	Revised the codes in order to deal with the coordinate file without SEQRES record

	3.3	Revised the codes so that the program no longer searched SEQRES record for chain
		information in order to be able to retrieve all peptide chains in a PDB coordinate
		file


	Revision 2, July 13, 2004:

	2.1	Revised the program options and arguments, from "rd_pdb 4-letter_pdbid_chainname [-i]"
		to "rd_pdb [-i] [-c <chainname>] pdb_file"

	2.2	Add program options parser getopt()

	2.3	Add usage(), help()


	Revision 1, June 27, 2004:

	1.1	Some minor modifications were made to deal with SIGATM and SIGUIJ record

	1.2	Add codes to deal with the MODEL/ENDMDL pair records in order to deal with the
		structures determined by NMR method

	1.3	Some modifications were made to improve the performance of the program

*/

/* Attention! Especially for HELIX, SHEET and SSBOND records!

The sequence number of residues is a new sequence of digits starting from 1 in PDBlike atomic
coordinate file created by program.

The residue(s) without complete backbone atoms in raw PDB file was(were) removed in ATOM record
and reserved in HELIX, SHEET and SSBOND records in PDBlike file. The new sequence number(s) of
the residue(s) is(are) to set to 0(zero) in HELIX, SHEET and SSBOND records in PDBlike file. */


#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include "max.h"
#include "bool.h"
/*
#include "dos.h"
*/

#ifdef	DOS
/*
#define	INPATH	"..\\pdb\\"
#define	OUTPATH	"..\\atom\\"
*/
#define	INPATH	".\\"
#define	OUTPATH	".\\"
#else
#define INPATH	"../pdb/"
#define OUTPATH	"../atom/"
/*
#define INPATH	"./"
#define	OUTPATH	"./"
*/
#endif


/* deal with the codes printing redundant messages */
#define	print_warning_message
#undef	print_warning_message


#define	SINGLECHAIN
#undef	SINGLECHAIN
/*
export HELIX, SHEET and SSBOND in one chain if defined(SINGLECHAIN),
else export these records in all chains */


#define	BREAKDIST	2.5
/*
maximum allowed peptide bond length. if distance is
 greater, a polypeptide chain interruption is assumed. */


/* max number of characters in the path name followed by a file name */
#define NAME_MAX	255




typedef struct option
{
	int	ignore;
	char	chain;		/* chain identifier */
} OPTION;


extern FILE	*open_file(const char *name, char *model);

extern int	SeqTrans(char *seq, int resflg);

extern int	WrtSeqFasta(char *seq, char *seqid, char *description, FILE *out);


/*
Retrieve the crystallographic data, including CRYST1 record, ORIGXn records and SCALEn records
*/
static void	rd_cryst(void);

/*
Retrieve the specified peptide chain named by "chain" from the PDB source file pointed to
by "in", export the coordinates file pointed to by "out", the coordinates file is in the
directory specified by "outp".
*/
static int	rd_chain(FILE *in, FILE *out, FILE *log,
			char chain, char *outp, int ignore_nonstd_res);


static void	rd_cryst(void)
{
}


void	usage(const char *prog);

void	help(const char *prog);


int	main(int argc, char *argv[])
/*
int	main_rd_pdb(int argc, char *argv[])
*/
{
	FILE	*in,			/* input file, file name suffix = ent/ent.Z */
		*out,			/* output file, filename suffix = seq */
		*log;			/* log_file, filename suffix = log */
	int	s, i, j,
		chainser,		/* chain serial number */
		chainnum;		/* number of chain */

	char	seq[CHAIN_MAX][RES_MAX * 3];
	char	*p = seq[0], *strp;

	/*
	long	input_fppos;
	*/

	char	chain[CHAIN_MAX],
		chainname, resname[4],
		chainid,		/* chain identifier	*/
		*inp, *outp, *logp,
		out_path[NAME_MAX], logfile[NAME_MAX],
		str[PDB_LINE_WIDTH];

	OPTION	*option;
	char	*prog = argv[0];
	char	*pdbfile;
	char	*seqfile;
	int	stringlen;
	int	opt = 0;
	char	*optstring = "hic:";

	extern char	*optarg;
	extern int	optind;


	int	record_flg = 0;

	int	pdbflg = 0;


	if(argc < 3 || argc > 6)
	{
		usage(prog);
		return 1;
	}

	if((option = (OPTION *)malloc(sizeof(OPTION))) == NULL)
	{
		printf("malloc() in main() in %s fails\n\n", prog);
		exit(1);
	}

	/* initialize the structure pointed by option */
	option->ignore = FALSE;
	option->chain = '\0';

	/* parse command line options */
	while((opt = getopt(argc, argv, optstring)) != -1)
	{
		switch(opt)
		{
			case 'i':
				option->ignore = TRUE;
				break;
			case 'c':
				if(strlen(optarg) != 1)
				{
					usage(prog);
					exit(1);
				}
				else
					option->chain = *optarg;
				break;
			case 'h':
				help(prog);
				exit(1);
			default:
				usage(prog);
				exit(1);
		}
	}

	/* shift past option arguments:

		The program requires at least (argc - optind) non-option arguments.
		The argv[optind] is the first non-option argument.
		Input file procedes output file in arguments list.
	*/ 
	argc -= optind;
	argv += optind;

	if(argc < 2)
	{
		printf("%s: missing required PDB coodinate file (input) or/and PDBlike coordinate file "
			"(output) argument\n\n", prog);

		usage(prog);
		exit(1);
	}

	pdbfile = *argv;

	/* check the type of the input file, namely pdbfile */
	if((inp = strrchr(pdbfile, '.')) == NULL)
	{
		printf("%s: unknown file type\n\n", pdbfile);
		exit(1);
	}
	else if(strstr(pdbfile, ".ent") == NULL)
	{
		printf("%s: expect a PDB coordinate file with suffix .ent or .pdb\n\n", pdbfile);
		exit(1);
	}

	/* shift past pdbfile argument */
	argc --;
	argv ++;

	if(argc < 1)
	{
		printf("%s: missing required PDB SEQRES sequence file (output) argument\n\n", prog);
		usage(prog);
		exit(1);
	}

	seqfile = *argv;

	/* shift past seqfile  argument */
	argc --;

	if(argc)
	{
		printf("%s: too many arguments\n\n", prog);
		usage(prog);
		exit(1);
	}

	if(memcmp(inp, ".Z", 2) == 0 || memcmp(inp, ".gz", 3) == 0)
	{
		char	cmd[NAME_MAX];

		pdbflg ++;

		strcat(memcpy(cmd, "zcat ", 6), pdbfile);

		if((in = popen(cmd, "r")) == NULL)
		{
			printf("%s fails to open %s\n", prog, pdbfile);
			exit(1);
		}
	}
	else if((in = open_file(pdbfile, "r")) == NULL)
	{
		printf("%s fails to open %s\n", prog, pdbfile);
		exit(1);
	}

	/* get the name of LOG file which always follows seqfile */
	strcat(memcpy(logfile, seqfile, strlen(seqfile) + 1), ".log");

	/* create all the non-existing parent directories specified by the path of seqfile */
	if(strrchr(logfile, '/') != NULL)
	{
		char	log_path[NAME_MAX];

		memcpy(log_path, "mkdir -p ", strlen("mkdir -p ") + 1);
		strcat(log_path, logfile);
		*strrchr(log_path, '/') = '\0';
		system(log_path);
	}

	if((log = open_file(logfile, "w")) == NULL)
	{
		printf("%s fails to open %s\n\n", prog, logfile);
		exit(1);
	}
	
	fprintf(log, "Reading %s ...\n", pdbfile);

	if(fgets(str, PDB_LINE_WIDTH, in) == NULL)
	{
		printf("Read %s error!\n\n", pdbfile);
		fprintf(log, "Read %s error!\n\n", pdbfile);
		exit(1);
	}


	/* reading primary structure section starts here */

	/* SEQRES record */

	for(fgets(str, PDB_LINE_WIDTH, in);
		memcmp("SEQRES", str, 6) != 0
		&& memcmp("HELIX ", str, 6) != 0 && memcmp("SHEET ", str, 6) != 0
		&& memcmp("SSBOND", str, 6) != 0 && memcmp("CRYST1", str, 6) != 0
		&& memcmp("MODEL ", str, 6) != 0 && memcmp("ATOM  ", str, 6) != 0;
		fgets(str, PDB_LINE_WIDTH, in));

	for(i = record_flg = 0; memcmp("SEQRES", str, 6) == 0; fgets(str, PDB_LINE_WIDTH, in))
	{
		if(strstr(str,"UNK") != NULL)
			fprintf(log, "Warning: unknown residue UNK, unknown residue!\n");
		else if(strstr(str, "ASX") != NULL)
			fprintf(log, "Warning: ambiguous residue ASX, ASP/ASN!\n");
		else if(strstr(str, "GLX") != NULL)
			fprintf(log, "Warning: ambiguous residue GLX, GLU/GLN!\n");

		sscanf(str+8, "%d", &chainser);

		if(isalpha(str[21]) == 0 || isalpha(str[20]) == 0 || isalpha(str[19]) == 0)
		{
			fprintf(log, "Nucleic acid chain [%c] in %s\n", chain[i], pdbfile);
			/*
			continue;
			*/
			exit(1);
		}

		if(1 == chainser)
		{
			if(0 == i)
			{
				seq[0][0] = '\0';

				p = seq[0];

				chain[i] = str[11];
				i ++;
				record_flg ++;

			}
			else
			{
				*p = '\0';

				seq[i][0] = '\0';

				p = seq[i];
				chain[i] = str[11];
				i ++;
			}
		}

		strp = str + 19;
		for(int jj = 0; jj < 13 && isalpha(*strp) != 0; jj ++)
		{
			memcpy(p, strp, 3);
			p += 3;

			strp += 4;
		}

	}
	*p = '\0';

	chainnum = i;

	if(record_flg < 1)
	{
		fprintf(log, "No SEQRES record found in %s\n", pdbfile);
	}
	else if(chainnum < 1)
	{
		fprintf(log,
			"No peptide chain found\n"
			"The program log was exported to %s and program exits.\n",
			logfile);

		exit(1);
	}
	/*
	else if(chainnum == 2)
	*/
	else if(chainnum > CHAIN_MAX)
	{
		printf("Error! Too many chains (%d chains) in protein structure\n", chainnum);
		printf("Please increase CHAIN_MAX (%d).\n", CHAIN_MAX);
		exit(1);
	}
	else
	{
	/* put the codes here which output the sequence presented in SEQRES record */
	/* The sequences must be output in FASTA format.	*/

		int	resflg = 3; /* original sequence is in three-letter abbreviation format */

		for(i = 0; i < chainnum; i ++)
		{
			SeqTrans(seq[i], resflg);
		}

		/* create all the non-existing parent directories specified by seqfile */
		/*
		if(strrchr(seqfile, '/') != NULL)
		{
			memcpy(out_path, "mkdir -p ", strlen("mkdir -p ") + 1);
			strcat(out_path, seqfile);
			*strrchr(out_path, '/') = '\0';
			system(out_path);
		}
		*/

		if((out = fopen(seqfile, "w")) == NULL)
		{
			fprintf(stderr, "%s fails to open %s\n\n", prog, seqfile);
			perror("fopen");
			exit(1);
		}
		

		char 	pdbid[NAME_MAX];
		char	*p;

		if((p = strrchr(pdbfile, '/')) == NULL)
			memcpy(pdbid, pdbfile, strlen(pdbfile) + 1);
		else	memcpy(pdbid, p + 1, strlen(p+1) + 1);
		
		char	chainame[2];
		for(i = 0; i < chainnum; i ++)
		{
			chainame[0] = chain[i];
			chainame[1] = '\0';
			if(WrtSeqFasta(seq[i], pdbid, chainame, out) != 0)
			{
				printf("length of chain #%c#:%d\n", chainname, strlen(seq[i]) / 3);
			}

		}
	}

	fprintf(log, "Found %d peptide chain(s) in SEQRES record in the %s\n", chainnum, pdbfile);

	if(pdbflg > 0)
		pclose(in);
	else	fclose(in);

	fclose(out);
	fclose(log);

	free(option);

	/* reading primary structure section ends here */

	return 0;
}


void	usage(const char *prog)
{
	printf("%s [-h] [-i] [-c <chain_identifier>] PDBfile PDB_SEQRES_sequence_file\n", prog); 
}


void	help(const char *prog)
{
	printf("Usage: %s [options] pdbfile seqfile\n"
		"Retrieve the atomic coordinates contained in the PDB coordinate file, pdbfile,\n"
		"and save the sequence in seqfile\n"
		"Options:\n"
		"-c CHAIN	Specify identifier CHAIN of peptide chain to process\n"
		"		Process all peptide chains in the PDB if missing the option\n"
		"-h		Show this message and exit\n"
		"-i		Ignore warning messages and continue process PDB file even if any amino\n"
		"		acid residue missing atom(s)\n"
		"\n"
		"pdbfile	PDB coordinate file\n"
		"seqfile	sequence file containing PDB SEQRES records in FASTA format\n"
		"\n"
		"Examples:\n"
		"%s -i pdbpath/pdb1ysc.ent seqpath/pdb1ysc.seq\n"
		"Process pdbpath/pdb1ysc.ent and ignore the warning messages though some\n"
		"amino acid resides(LYS 1) missed some atoms, and then output seqpath/pdb1ysc.seq\n"
		"\n"
		"%s -c B pdbpath/pdb9ldb.ent\n"
		"Process B chain in pdbpath/pdb9ldb.ent and then output seqpath/pdb9ldb.seq\n"
		"\n"
		"%s pdbpath/pdb9ldb.ent\n"
		"Process all chains in pdbpath/pdb9ldb.ent and then output seqpath/pdb9ldb.seq\n",
		prog, prog, prog, prog);
}
