#ifndef PROPERTY4RES_H
#define PROPERTY4RES_H

/*
#define N_STDRES	24
#define N_STDRES	36
#define N_STDRES	42
*/
/* number of standard residues + an unknown residue */


/*
# COLUMNS:
# ResidueName NumberOfAtomsInStandardResidue Polarity Hydrophility Charge
#
# Polarity: 1 for polar residue, -1 for nonpolar residue, 0 for others
# Hydrophility: 1 for hydrophilic residue, -1 for hydrophobic residue, 0 for special residue
# Charge: 1 for positively charged residue, -1 for negatively charged residue, 0 for others
#
ALA   5  -1 -1  0
ARG  11   1  1  1
ASN   8   1  1  0
ASP   8   1  1 -1
CYS   6   1  0  0
GLN   9   1  1  0
GLU   9   1  1 -1
GLY   4  -1  0  0
HIS  10   1  1  1
ILE   8  -1 -1  0
LEU   8  -1 -1  0
LYS   9   1  1  1
MET   8  -1 -1  0
PHE  11  -1 -1  0
PRO   7  -1  0  0
SER   6   1  1  0
THR   7   1  1  0
TRP  14  -1 -1  0
TYR  12   1 -1  0
VAL   7  -1 -1  0
END
*/


static struct stdresidue
{
	char	*three;
	char	single;
	int	natm;
	int	pol;
	int	philic;
	int	e;
} stdres[] =
{
	{"ALA", 'A', 5,  -1, -1,  0},
	{"ARG", 'R', 11,  1,  1,  1},
	{"ASN", 'N', 8,   1,  1,  0},
	{"ASP", 'D', 8,   1,  1, -1},
	{"CYS", 'C', 6,   1,  0,  0},

	{"GLN", 'Q', 9,   1,  1,  0},
	{"GLU", 'E', 9,   1,  1, -1},
	{"GLY", 'G', 4,  -1,  0,  0},
	{"HIS", 'H', 10,  1,  1,  1},
	{"ILE", 'I', 8,  -1, -1,  0},

	{"LEU", 'L', 8,   1, -1,  0},
	{"LYS", 'K', 9,   1,  1,  1},
	{"MET", 'M', 8,  -1, -1,  0},
	{"PHE", 'F', 11, -1, -1,  0},
	{"PRO", 'P', 7,  -1,  0,  0},

	{"SER", 'S', 7,   1,  1,  0},
	{"THR", 'T', 7,   1,  1,  0},
	{"TRP", 'W', 14, -1, -1,  0},
	{"TYR", 'Y', 12,  1, -1,  0},
	{"VAL", 'V', 7,  -1, -1,  0},

	{"ASX", 'B', 8,   1,  1,  0},
	{"GLX", 'Z', 9,   1,  1,  0},
	{"CSS", 'c', 6,   0,  0,  0},

	{"ZN",  'a', 1,   1,  1,  1},
	{"HEM", 'b', 39,  0,  0,  0},
	{"ACE", 'd', 3,   0,  0,  0},
	{"HOH", 'e', 1,   1,  1,  0},
	{"WAT", 'e', 1,   1,  1,  0},
	{"DOD", 'f', 1,   1,  1,  0},

	{"A",   '1', 0,   0,  0,  0},
	{"+A",  '1', 0,   0,  0,  0},
	{"C",   '2', 0,   0,  0,  0},
	{"+C",  '2', 0,   0,  0,  0},
	{"G",   '3', 0,   0,  0,  0},
	{"+G",  '3', 0,   0,  0,  0},
	{"I",   '4', 0,   0,  0,  0},
	{"+I",  '4', 0,   0,  0,  0},
	{"T",   '5', 0,   0,  0,  0},
	{"+T",  '5', 0,   0,  0,  0},
	{"U",   '6', 0,   0,  0,  0},
	{"+U",  '6', 0,   0,  0,  0},
	{"UNK", 'X', 8,   0,  0,  0}
};



static struct stdNucleicAcid
{
	char	*three;
	char	single;
	int	natm;
} stdna[] =
{
	"  A", 'A', 0,
	" +A", 'A', 0,
	"  C", 'C', 0,
	" +C", 'C', 0,
	"  G", 'G', 0,
	" +G", 'G', 0,
	"  I", 'I', 0,
	" +I", 'I', 0,
	"  T", 'T', 0,
	" +T", 'T', 0,
	"  U", 'U', 0,
	" +U", 'U', 0,
	"UNK", 'X', 0
};

#endif
