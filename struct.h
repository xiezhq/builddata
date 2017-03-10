#ifndef STRUCT_H
#define STRUCT_H



struct res2pdbres
{
	unsigned	resnum;
	unsigned	pdbresnum;
	char		pdbins;
};


struct resnum2ins
{
	unsigned	resnum;
	char		resnum_ins[7];
};


struct hbond
{
	unsigned	res1, res2;
	char		resnum_ins1[7], resnum_ins2[7];
};


#ifndef ELE2ATM_
#define ELE2ATM_
struct unfoldres
{
	double	ap, pol, asa, /* ASA for amino acids in unfolded conformations */
		bu_ex, ex_u, bb; /* conformational entropies for amino acids */
};

struct ele2atm
{
	double		*probe, *x, *y, *z, *radx;

	double		*asa_ap, *asa_pol;

	unsigned	*resid_atm;
			/* sequence number of residue containning the current atom */

	unsigned	*ap, *natmres,
			/*
			*elenum,
			*/
			nres,
			/* number of residues in an element */

			count_res, count_atm, count_ss,
			/* total of the residues, atoms, disulfide bonds included in from the first
			to the current element respectively, the sequence of elements refered to is
			the sequence reordered according to eigenvectors */

			*sslen,

			row;
			/* subscription for row of element contact matrix */

	struct unfoldres	*unfold;
	struct res2pdbres	*resnum;
};
#endif


#ifndef CHAIN_
#define CHAIN_
struct atom
{
	char		name[5];
			/* atom name */

	/*
	unsigned	resid_atm;
	*/
			/* sequence number of residue containning the current atom */

	double		x, y, z;
			/* orthogonal coordinates */

	struct atom	*next;
};

struct residue
{
	char	name[4];
		/* residue name */

	unsigned	resnum;
			/* sequence number of residue */

	char		resnum_ins[7];
	char		sstdef;

	unsigned	pdbresnum;
	char		pdbins;

	struct atom	*atm;
	struct residue	*prior;
	struct residue	*next;
};

struct element
{
	unsigned	nres;
			/* number of residues in an element */

	struct residue	*res;
};
#endif


#ifndef CONF_
#define CONF_
typedef struct
{
        char    name[5];
        double  bu_ex, ex_u, bb;
} CONF;
#endif

#ifndef BTNODE_
#define BTNODE_
typedef struct btnode
{
	unsigned	ndom;
	struct ele2atm	*start, *end;
	struct btnode	*left, *right;
} btreenode;
#endif /* binary tree node for binary partition */

#ifndef SLNODE_
#define SLNODE_
struct lnode
{
	double	vector;
	int	seqnum;

	struct lnode	*next;
};
#endif /* single linked list */

#ifndef DLNODE_
#define DLNODE_
struct dlnode
{
	double	vector;
	int	seqnum;

	struct dlnode	*prior, *next;
};
#endif /* double linked list */
#endif /* STRUCT_H */
