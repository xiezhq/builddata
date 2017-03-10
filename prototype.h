#ifndef PROTOTYPE_H
#define PROTOTYPE_H


#define TRUE	1
#define FALSE	0
#define true	1
#define false	0



/* parameters for different contacts type */

#define atm2atm_contact		0.0865052
#define hbond_contact		10.91327
#define N2C_dist		2.5
#define N2C_dist_sqr		6.25
#define peptidebond_contact	28.2662
/*
#define peptidebond_contact	16.9852
*/
#define backboneweight		5.0
#define sstweight		1.0
#define DIST2			16.0

/*
atm2atm_contact		atom-atom contacts, distance between two atoms is 3.4 angstroms,
 			atm2atm_contact = 1 / (3.4 * 3.4) = 0.0865052

hbond_contact		a hydrogen bond in sheet corresponds to 15 atom-atom contacts with an average length
			of 2.8 angstroms which is the average length for hydrogen bond(hydrogen Donor-Acceptor),
			15 * atm2atm_contact_2.8A = 15 * (1 / (2.8*2.8)) = 1.91327
			http://www-structure.llnl.gov/xray/tutorial/protein_structure.htm

N2C_dist		maximum allowed peptide bond length. if distance is
			greater, a polypeptide chain interruption is assumed.
			refer to DSSP program

N2C_dist_sqr		square of N2C_dist

peptidebond_contact	contacts for covalent bond N-C with an average length of 1.33 angstroms,
			28.2662 = 50 * (1.0 / (1.33 * 1.33))
			http://www-structure.llnl.gov/xray/tutorial/protein_structure.htm

backboneweight		weight for backbone-backbone atoms contacts,
			5 * atm2atm_contact

sstweight		weight for contacts occuring in secondary structure environment,
			2 * atm2atm_contact

DIST2			distance cutoff(unit is expressed in angstroms**2), their atom-atom contact is ignored
			if the square of distance between two atoms is larger than DIST2.
			4.0 * 4.0 = 16.0 angstrom**2


		Mainchain bond lengths
(http://www-structure.llnl.gov/xray/tutorial/engtabl.html)
--------------------------------------------------------------
Bond           | X-PLOR labelling              | Value | sigma 
--------------------------------------------------------------
C-N            | C-NH1        (except Pro)     | 1.329 | 0.014
               | C-N          (Pro)            | 1.341 | 0.016
               |                               |       |      
C-O            | C-O                           | 1.231 | 0.020
               |                               |       |      
Calpha-C       | CH1E-C       (except Gly)     | 1.525 | 0.021
               | CH2G*-C      (Gly)            | 1.516 | 0.018
               |                               |       |      
Calpha-Cbeta   | CH1E-CH3E    (Ala)            | 1.521 | 0.033
               | CH1E-CH1E    (Ile,Thr,Val)    | 1.540 | 0.027
               | CH1E-CH2E    (the rest)       | 1.530 | 0.020
               |                               |       |      
N-Calpha       | NH1-CH1E     (except Gly,Pro) | 1.458 | 0.019
               | NH1-CH2G*    (Gly)            | 1.451 | 0.016
               | N-CH1E       (Pro)            | 1.466 | 0.015
--------------------------------------------------------------
reference: 
Engh R A & Huber R (1991). Accurate bond and angle parameters for X-ray protein structure refinement.
Acta Cryst., A47, 392-400.


				Types of Interactions
		(http://pps98.man.poznan.pl/ppscore/section7/interact.html)
-------------------------------------------------------------------------------------------------------------------------
interaction	| Nature	| "Bond" Length		|"Bond" Strength| Stabilisation /	| Example  
	   	|	    	|	("optimum")	|		| Free Energy("optimum")|
-------------------------------------------------------------------------------------------------------------------------
	   	|	    	|	angstrom	|  kcal/mol	|	kcal/mol	|
-------------------------------------------------------------------------------------------------------------------------
Atomic Bond	|Covalent	|  1.0-1.6		|	>50	|			|Peptide Bond
-------------------------------------------------------------------------------------------------------------------------
Ion Pair	|		|1.8-4.0(3.5)		|		|0.1-0.5 solvent exposed|Positive Arg NH,
(Salt Bridge)	|Electrostatic	|opposite charges	|1-6		|			|Lys NZ, His NE & ND,
		|		|			|		|			|Amino terminus N.
		|		|3.5-10.0 like charges	|		|3-6 buried destability	|Negative Asp OD,
		|		|			|		|			|Glu OE, Carboxy terminus O  
-------------------------------------------------------------------------------------------------------------------------
hydrogen bond	|hydrogen bond	|2.6-3.5(2.8)		|2-10		|0.5-2.0		|Any hydrogen donor N
                |               |                       |               |solvent <-> buried	|to any acceptor O 
-------------------------------------------------------------------------------------------------------------------------
hydrophobic	|Entropy	|			|2-3		|0.0-2.5 (1.5)		|Side chain Leu, Ile, Val,
		|		|			|		|per methylene group	|Phe, Met, Trp, Ala, Cys,
		|		|			|		|			|Pro, Tyr 
-------------------------------------------------------------------------------------------------------------------------
van der Waals	|Dispersion/	|2.8-4.0		|<1		|			|All atoms
		|Repulsion	|			|		|			|(e.g. aliphatic hydrogen)
-------------------------------------------------------------------------------------------------------------------------
For more contents refer to http://pps98.man.poznan.pl/ppscore/section7/interact.html
*/



/* max length of filename with path */
#define filenamelen	256


#define dssplinelen	140
/* newline_character + width of dssp line, namely, 1 + 139 */


/* define the width of lines in atomic coordinate file created by program rd_pdb */
#ifndef LINE_WIDTH
#define LINE_WIDTH	89 /* line_width + 2 */
#endif



/* critical quantity relative to proteins */
#ifndef PROTEIN_
#define PROTEIN_

#define MIN_GIBBS		-9.9e6
#define gap_dist		10
#define term_segmentsize_min	10
#define segmentsize_min		20
#define dom_max			10
#define segment_max		20
#define DOM_MIN			34
#define DOM_MIN_DOUBLE		68
/*
#define DOM_MIN			40
#define DOM_MIN_DOUBLE		80
*/
#define HELIX_MAX		100
#define STRAND_MAX		100
#define SS_MAX			100
#define HB_MAX			500
#define RES_MAX			1000
#define ATM_MAX			10000

/*
MIN_GIBBS		minimal unfolding gibbs energy of an unfolded protein, unit is expressed in cal/mol

gap_dist		gaps occur at the N- and C-terminus and where sequencial C-alpha atoms
			have a distance larger than gap_dist(angstroms)

term_segmentsize_min	size of N-terminus of which residues are integrated into a whole, namely, an undivided element

segmentsize_min		minimal segments length on the pieces assigned to one or the orther domain,
			but there is an exception to N- and C-terminus since 10 residues close to
			N- and C-terminus are integrated as an element respectively.

dom_max			10 = 1000 / 100. where 100 is the size of domains, the distribution of which
			is the greatest in all domains with the different size.

segment_max		max number of segments in one domain
			The number is 6 in the paper of Holm L & Sander C..
			Refer to Fig4(f) in
			Holm L and Sander C. Parser for protein folding units. Proteins, 1994, 19:256-268.
			However, larger value is safer for running program when definition of domains is not correct.

DOM_MIN			minimum size in residues of the accepted domain
			40 residues is the usual size adopted by many authors.
			34 residues is the size of the minimized protein obtained in wet experiments.
			Refer to Proc Natl Acad Sci USA, 1997, 94:10010-10011

DOM_MIN_DOUBLE		MIN * 2

HELIX_MAX		max number of helix in a peptide chain

STRAND_MAX		max number of strand in a peptide chain

SS_MAX			max number of disulfide bonds in a peptide chain

HB_MAX			max number of hydrogen bonds in a protein

RES_MAX			max number of residues in a pepetide chain

ATM_MAX			max number of atoms in a pepetide chain */
#endif /* define PROTEIN_ */



#define solvent_radius 1.4
/* default is the radius of water molecule */


#ifndef STRUCT_H
#include "struct.h"
#endif


/* to include the definition of structure FILE in stdio.h */
#include <stdio.h>


/* open a file */
extern FILE	*open_file(const char *name, char *model);


/* allocate and free storage */

extern void	*memalloc(unsigned n);

extern float	**fmatrix(unsigned n, unsigned m);

extern double	**dmatrix(unsigned n, unsigned m);

extern void	free_fmatrix(float **matrix);

extern void	free_dmatrix(double **matrix);


/* input data from datafile to matrix */
extern void	get_data(double **matrix, char *datafile, unsigned n);


/* get necessary informations from atomic ordinates file */
extern struct element	*getprotein(const char *file,
				unsigned *nele, unsigned *nres, unsigned *natm,
				unsigned *ss1, unsigned *ss2, unsigned *nss);


extern void	rd_dssp(char *dsspfile, unsigned *nhb, struct hbond *headhb);

extern void	constructmatrix_ca(struct element *ele, unsigned nele, double **matrix);

extern void	contact(struct element *ele, unsigned nele, double **ele_contact);

extern void	contact_ele_reord(struct element *ele, unsigned nele, unsigned *seq,
		double **contact_reord);

extern void	addhb2matrix(struct hbond *headhb, unsigned nhb, struct element *hele, unsigned nele,
		double **matrix);

extern void	eigen(double **matrix, double *eigenvector, double *init, unsigned n);

extern void	schmidt(double *vector, unsigned n);

extern void	reorder_list(double *eigenvector, unsigned *sequence, unsigned n);

extern void	reorder_array(double *eigenvector, unsigned *sequence, unsigned n);

extern void	addinfo2seq(const unsigned *seq, unsigned nele, const struct element *hele,
		struct ele2atm *e2a,
		unsigned nss, unsigned *ss1, unsigned *ss2,

		/* addresses pointed by following pointers to be accessed by partition() and freenergy() */
		double *x, double *y, double *z, double *radx, double *probe,
		unsigned *resid_atm,
		unsigned *ap,
		unsigned *natmres,
		struct unfoldres *unfold,
		struct res2pdbres *resnum,
		unsigned *sslen);

extern double	protein(const struct ele2atm *he2a, unsigned nres, unsigned natm, unsigned nss);

extern btreenode	*partition(const struct element *hele, struct ele2atm *heade2a, struct ele2atm *terme2a,
				struct ele2atm *he2a, struct ele2atm *te2a,
				unsigned nres, unsigned natm, unsigned nss, double **contact_reord);


extern void	cnctarea(double errorpar, double *probe, double *x, double *y,
			double *z, double *radx, int natms, double *accss);

/*
extern double	freenergy(double *probe, double *x, double *y, double *z, double *radx,
			unsigned natm,
			const unsigned *ap, struct unfoldres *unfold,
			unsigned *ss_len, unsigned nss,
			unsigned nres, unsigned *natmres, double T);
*/

/* calculate unfolding Gibbs free energy for each structure unit */
extern double	freenergy(const struct ele2atm *he2a, unsigned natm, unsigned nss, unsigned nres, double T,
			double *asa_dom);


extern double	intracontact(const double *x, const double *y, const double *z, const unsigned *resid_atm, unsigned natm,
			double dist2,
			double *contact);


extern double	intercontact(const double *x1, const double *y1, const double *z1, const unsigned *resid_atm1, unsigned natm1,
			const double *x2, const double *y2, const double *z2, const unsigned *resid_atm2, unsigned natm2,
			double dist2,
			double *contact);


extern double	intracontact_addup(const struct ele2atm *he2a, unsigned nele,
			double **matrix,
			double *contact);


extern double	intercontact_addup(const struct ele2atm *he2a1, unsigned nele1,
			const struct ele2atm *he2a2, unsigned nele2,
			double **matrix,
			double *contact);

extern void	asa_addup(const double *asa_ap, const double *asa_pol, unsigned nres,
			double *asa_sum);

#endif /* PROTOTYPE_H */
