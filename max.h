#ifndef MAX_H
#define MAX_H


/*
 * The file was writen by Xie Zhiqun in 1997 and modified in 1998,1999,2000.
 * It is its primary purpose to read in and filter the PDB file.
 * If you modify the file, please let me know the modification.
 * Send the modification or question to xiezhq@yahoo.com or xiezhq@dna.sibc.ac.cn
 */

#define NAME_MAX 255	/* max length of file name associated with path */
#define FILE_NAME_LEN	255 /* max length of file name associated with path */

#define	RES_MAX 1000	/* max number of the residues in proteins */
#define RES_MAX2 1000000 /* RES_MAX * RES_MAX */


#define	ATM_MAX 24000	/* max number of the atoms in proteins */
/* if (ATM_MAX > 24708) a stack error will happen when a program(rd_pdb) run on windows platform!
 * but not happen on UNIX platform! why? I don't know!
 */


#define ATM_MAX_OF_RESIDUE 14	/* maxinum of number of atoms in residues */
#define SSNUM 500	/* max number of the disulfide bonds in proteins	*/
#define DOM_MIN 30	/* min number of the residues for a domain in proteins */
#define CHAIN_MAX 20	/* max number of the chain in proteins */
#define SSTNUM 200	/* max of number of the secondary structure units in proteins */
#define HELIX_MAX 500	/* max number of the helixs in a protein */
#define STRAND_MAX 500	/* max of number of the strands in a protein */

#define	PDB_LINE_WIDTH 80 + 1 + 1	/* The real width(bytes) of each line in PDB files is 80 bytes. */


#endif /* MAX_H */
