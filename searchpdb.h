#ifndef PDBFLT_H
#define PDBFLT_H


/* max number of characters in the path name followed by a file name */
#define	NAME_MAX	255


/* max number of characters in one line in a regular file */
#define LINE_LEN_MAX	255


/*
#define CMD	"rd_pdb -i "
for system(command)
*/

#define CMD	"asa2pep par/ProtOr.resi-defs.mod-by-xie.dat par/ProtOr.atom-defs.dat glyXgly.ProtOr.dat "
/* full command:
asa2pep par/ProtOr.resi-defs.mod-by-xie.dat par/ProtOr.atom-defs.dat glyXgly.ProtOr.dat pdbfile outputfile
*/

/* full command:
interf pdbfile nchain01 nchain02 chains01 chains02 cutoff4surf cutoff4interf \
par/ProtOr.resi-defs.mod-by-xie.dat par/ProtOr.atom-defs.dat glyXgly.ProtOr.dat
*/
/*
#define CMD	"interf pdbfile nchain01 nchain02 chains01 chains02 cutoff4surf cutoff4interf \
par/ProtOr.resi-defs.mod-by-xie.dat par/ProtOr.atom-defs.dat glyXgly.ProtOr.dat "
#define CMD	"interf "
*/



/* for substitution of suffix of file */
/*
#define SUFFIX	".log"
#define SUFFIX	".asa"
*/
#define SUFFIX	".acc"


/* for ftw() */
#define PDB_DIR_DEPTH	2


/* for getopt() */
extern char	*optarg;

typedef struct option
{
	char	*dir,
		*list;
} OPTION;

#endif /* PDBFLT_H */
