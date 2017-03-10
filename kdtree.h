#ifndef KDTREE_H
#define KDTREE_H

#include "bool.h"

/* compile option */
#define DEBUG
#define MEDIAN

/*
#undef	DEBUG
*/
/*
#undef	MEDIAN
*/

#define NDIM	3 /* number of dimensions of each data point */


/* Some things to manage simple lists - structure that begin
 * with a pointer to the next element in the list.
 */
struct slList
{
	struct slList	*next;
};

/* an element on a doubly linked list */
struct dlNode
{
	void		*val;
	double		dist;
	struct dlNode	*prev;
	struct dlNode	*next;
};

/* a doubly linked list */
struct dlList
{
	struct dlNode	*head;
	struct dlNode	*tail;
};


#define X	0 /* X dimension */
#define Y	1 /* Y dimension */
#define Z	2 /* Z dimension */
/*
#undef X
#undef Y
#undef Z
*/


#ifdef MEDIAN
#define LEFT	0 /* left child */
#define MID	1 /* median */
#define RIGHT	2 /* right child */

struct kdLeaf
{
	double		coord[NDIM]; /* coordinates of the first point in leaf */
	struct kdLeaf	*next; /* next leaf in list */
	int		hit; /* flag */
};

struct kdBranch
{
	struct kdLeaf	*leaf; /* extra info for leaves on tree */
	int		dim; /* split dimension */
	struct kdBranch	*lo; /* pointer to children with lower coordinates */
	struct kdBranch	*hi; /* pointer to children with higher coordinates */
};
#else
struct kdLeaf
{
	double		coord[NDIM]; /* coordinates of the first point in leaf */
	struct kdLeaf	*next; /* next leaf in list */
	bool		hit; /* flag */
};

struct kdBranch
{
	struct kdLeaf	*leaf; /* extra info for leaves on tree */
	int		dim; /* split dimension */
	double		cutCoord; /* where to split space based on median value */
	struct kdBranch	*lo; /* pointer to children with lower coordinates */
	struct kdBranch	*hi; /* pointer to children with higher coordinates */
};
#endif

struct kdTree
/* the whole tree */
{
	struct kdBranch	*root; /* pointer to root of kd-tree */
};


#define dlEnd(node)	(node->next == NULL)
/* True if node past end */



/*
 * Prototype of functions
 */
extern struct kdTree	*kdTreeMake(struct kdLeaf *leafList, int nodeCount, int ndim);
extern void	traverseKdTree(struct kdTree *tree);

extern void	dlListInit(struct dlList *lists);


extern void	nearSearch_(struct kdTree *tree, struct kdLeaf *target, int ndimension,
			double *hr[2], double distMax,
			int nNear,
			struct dlList *nearList);

/*
 * data		two dimension array, each row is for coordinates of each data point
 * npoint	number of data points in input
 * ndimension	number of dimensions for each data point
 * target	array, coordinate of each dimension for target point
 * hr		hyperrectangle region to search
 * distMax	cutoff of distance to search
 * nNear	value: "nNear" nearest neighbors search if nNear > 0, else near neighbors search.
 * 		nNear is the cutoff for number of neighbors in "nNear" nearest neighbors search.
 * nearList	results, near points and their distance to target
 */
extern void	nearSearch(double **data, int npoint, int ndimension,
			double *target, double *hr[2], double distMax,
			int nNear,
			struct dlList *nearList);

extern void	outPutNearList(struct dlList *list);

extern void	freeNodeList(struct dlNode *node);
extern void	freeKdTree(struct kdTree *tree);

#endif
