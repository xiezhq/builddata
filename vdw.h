#ifndef VDW_H
#define VDW_H

/* The packing density in proteins: standard radii and volumes. Tsai J, Taylor R, Chothia C and Gerstein M.
JMB, 1999, 290: 253-266. */

/* standard van der Waals radii of the various atomic groups

- All data are adapted from the beforementioned paper.
*/


#define N_HYBRID	20

static struct atomhybrid
{
	char	*hybrid;
		/* name of atomic group defined by status of hybrid of molecular orbit,
		refer to JMB, 1999, 290:253-266.

		Nomenclature:
		The atomic groups found in proteins are given labels of the general form XnHmx, where
		X indicates the chemical nature of the non-hydrogen atoms; n, their valence; Hm, the
		number (m) of hydrogen atoms attached to the non-hydrogen atom; and x, the atom type
		assigned according to their volumes: s, b, or u (small, big, or unique, respectively).
		*/

	double	r;
		/* VDW radius of atomic group defined by Gerstein, JMB, 1999, 290:253-266.
		In Gerstein's definition, atomic group subsume a heavy-atom and its covalently attached 
		hydrogen atoms into one moiety because the positions of hydrogen atoms in protein structures
		are generally not known.
		*/
} atmhybrid[] = /* N_HYBRID */
{
	"C3H0s",  1.61,
	"C3H0b",  1.61,
	"C4H1s",  1.88,
	"C4H1b",  1.88,
	"C3H1s",  1.76,

	"C3H1b",  1.76,
	"C3H2u",  1.76,
	"C4H2s",  1.88,
	"C4H2b",  1.88,
	"C4H3u",  1.88,

	"N3H0u",  1.64,
	"N3H1s",  1.64,
	"N3H1b",  1.64,
	"N4H3u",  1.64,
	"N3H2u",  1.64,

	"O1H0u",  1.42,
	"O2H1u",  1.46,
	"O2H2u",  1.46,
	"S2H0u",  1.77,
	"S2H1u",  1.77,

	"P4",  0.00,
	"O1N3HH", 1.60,
	"Z2ION",  0.74,
	"FE",  1.70,
	"DEFAULT", 1.80
};

#endif
