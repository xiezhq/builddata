#ifndef PDBFLT_H
#define PDBFLT_H


/* max number of characters in the path name followed by a file name */
#define	NAME_MAX	255


/* max number of characters in one line in a regular file */
#define LINE_LEN_MAX	255


/* for system(command) */
#define COM	"rd_pdb -i "


/* for ftw() */
#define PDB_DIR_DEPTH	2


/* for getopt() */
extern char	*optarg;

typedef struct option
{
	char	*dir,
		*list;
} OPTION;


extern FILE	*open_file(const char *file, char *mode);


#endif /* PDBFLT_H */
