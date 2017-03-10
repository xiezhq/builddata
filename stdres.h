#ifndef STDRES_H
#define STDRES_H


/*
 * The file was writen by Xie Zhiqun in 1997 and rewriten in 1998,1999,2000.
 * It is its primary purpose to list 20 familiar amino acid residues.
 * If you modify the file, please let me know the modification.
 * Send the modification or question to xiezhq@yahoo.com or xiezhq@dna.sibc.ac.cn
 */


/* refer to PDB format description version 2.2 and Appendix5(formulas and molecular weights
	for standard residues) */

static struct stdresidue
{
	char	*three;
	char	single;
	int	atmnum;
} stdres[] =
{
	"ALA", 'A', 5,
	"ARG", 'R', 11,
	"ASN", 'N', 8,
	"ASP", 'D', 8,
	"CYS", 'C', 6,

	"GLN", 'Q', 9,
	"GLU", 'E', 9,
	"GLY", 'G', 4,
	"HIS", 'H', 10,
	"ILE", 'I', 8,

	"LEU", 'L', 8,
	"LYS", 'K', 9,
	"MET", 'M', 8,
	"PHE", 'F', 11,
	"PRO", 'P', 6,

	"SER", 'S', 7,
	"THR", 'T', 14,
	"TRP", 'W', 12,
	"TYR", 'Y', 7,
	"VAL", 'V', 7,

	"ASX", 'B', 8,
	"GLX", 'Z', 9,
	"UNK", 'X', 8 };

/* --------------------------------------------------------- */

static const char	stdaa[] =	"A"	/* ALA */
					"R"	/* ARG */
					"N"	/* ASN */
					"D"	/* ASP */
					"C"	/* CYS */

					"Q"	/* GLN */
					"E"	/* GLU */
					"G"	/* GLY */
					"H"	/* HIS */
					"I"	/* ILE */

					"L"	/* LEU */
					"K"	/* LYS */
					"M"	/* MET */
					"F"	/* PHE */
					"S"	/* SER */

					"T"	/* THR */
					"W"	/* TRP */
					"Y"	/* TYR */
					"V"	/* VAL */
					"P";	/* PRO */
			


/* The number of atoms for the corresponding residues above is as follows. */

static int	resatm[] = {	5,
				11,
				8,
				8,
				6,

				9,
				9,
				4,
				10,
				8,

				8,
				9,
				8,
				11,
				6,

				7,
				14,
				12,
				7,
				7	};


/* The description of all atoms for the corresponding residues above is as follows. */

static char	*stdatm[] = {	"N","CA","C","O","CB",
				"N","CA","C","O","CB","CG","CD","NE","CZ","NH1","NH2",
				"N","CA","C","O","CB","CG","OD1","ND2",
				"N","CA","C","O","CB","CG","OD1","OD2",
				"N","CA","C","O","CB","SG",

				"N","CA","C","O","CB","CG","CD","OE1","NE2",
				"N","CA","C","O","CB","CG","CD","OE1","OE2",
				"N","CA","C","O",
				"N","CA","C","O","CB","CG","ND1","CD2","CE1","NE2",
				"N","CA","C","O","CB","CG1","CG2","CD1",

				"N","CA","C","O","CB","CG","CD1","CD2",
				"N","CA","C","O","CB","CG","CD","CE","NZ",
				"N","CA","C","O","CB","CG","SD","CE",
				"N","CA","C","O","CB","CG","CD1","CD2","CE1","CE2","CZ",
				"N","CA","C","O","CB","OG",

				"N","CA","C","O","CB","OG1","CG2",
				"N","CA","C","O","CB","CG","CD1","CD2","NE1","CE2","CE3","CZ2","CZ3","CH2",
				"N","CA","C","O","CB","CG","CD1","CD2","CE1","CE2","CZ","OH",
				"N","CA","C","O","CB","CG1","CG2",
				"N","CA","C","O","CB","CG","CD"	};


#endif /* STDRES_H */
