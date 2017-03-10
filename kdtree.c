#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include "kdtree.h"
#include "memory.h"
#include "bool.h"


#ifdef DEBUG
static int	computationCount = 0;
#endif

/* number of dimensions of each data point */
static int	ndim = 0;

static int	dimCmp = 0; /* dimension for comparison and sorting */


/* Return the least value */
static double	min(double a, double b)
{
	if(a < b)
		return a;
	else	return b;
}


/* Initialize list to be empty */
void dlListInit(struct dlList *dl)
{
	dl->head = NULL;
	dl->tail = NULL;
}

/* Removes a node from list */
static void	dlRemove(struct dlList *list, struct dlNode *node)
{
	struct dlNode	*before = node->prev;
	struct dlNode	*after = node->next;

	if(before == NULL)
	{
		after->prev = before;
		list->head = after;
	}
	else if(after == NULL)
	{
		before->next = after;
		list->tail = before;
	}
	else
	{
		before->next = after;
		after->prev = before;
	}
	node->prev = NULL;
	node->next = NULL;
}

/*
static void	dlInsertBetween(struct dlNode *before, struct dlNode *after, struct dlNode *newNode)
{
	before->next = newNode;
	newNode->prev = before;
	newNode->next = after;
	after->prev = newNode;
}
*/

/* Add a node to tail of list */
static void	dlAddTail(struct dlList *list, struct dlNode *newNode)
{
	struct dlNode	*tail = list->tail;

	newNode->prev = tail;
	newNode->next = NULL;

	if(tail == NULL)
	{
		list->head = newNode;
	}
	else
	{
		tail->next = newNode;
	}
	list->tail = newNode;
}


/* compare to sort based on dimCmp dimension */
static int	kdLeafCmp(const void *va, const void *vb)
{
	extern int	dimCmp;

	const struct kdLeaf	*a = *((struct kdLeaf **)va);
	const struct kdLeaf	*b = *((struct kdLeaf **)vb);

	if(a->coord[dimCmp] > b->coord[dimCmp])
		return 1;
	else if(a->coord[dimCmp] < b->coord[dimCmp])
		return -1;
	else	return 0;
}

/* compare to sort based on x coordinate */
/*
static int	kdLeafCmpX(const void *va, const void *vb)
{
	const struct kdLeaf	*a = *((struct kdLeaf **)va);
	const struct kdLeaf	*b = *((struct kdLeaf **)vb);

	if(a->coord[X] > b->coord[X])
		return 1;
	else if(a->coord[X] < b->coord[X])
		return -1;
	else	return 0;
}
*/

/* compare to sort based on y coordinate */
/*
static int	kdLeafCmpY(const void *va, const void *vb)
{
	const struct kdLeaf	*a = *((struct kdLeaf **)va);
	const struct kdLeaf	*b = *((struct kdLeaf **)vb);

	if(a->coord[Y] > b->coord[Y])
		return 1;
	else if(a->coord[Y] < b->coord[Y])
		return -1;
	else	return 0;
}
*/

/* compare to sort based on z coordinate */
/*
static int	kdLeafCmpZ(const void *va, const void *vb)
{
	const struct kdLeaf	*a = *((struct kdLeaf **)va);
	const struct kdLeaf	*b = *((struct kdLeaf **)vb);

	if(a->coord[Z] > b->coord[Z])
		return 1;
	else if(a->coord[Z] < b->coord[Z])
		return -1;
	else	return 0;
}
*/


/* Helper structure for sorting dlNodes preserving order */
struct dlSorter
{
	struct dlNode	*node;
};

/* Node comparison pointer, just used by dlSortNodes and helpers. */
static int	(*compareFunc)(const void *elem1, const void *elem2);

/* Compare two dlSorters indirectly, by calling compareFunc. */
static int	dlNodeCmp(const void *elem1, const void *elem2)
{
	struct dlSorter	*a = (struct dlSorter *)elem1;
	struct dlSorter	*b = (struct dlSorter *)elem2;

	return compareFunc(&a->node->val, &b->node->val);
}

/* Count up elements in list */
/*
static int	slCount(void *list)
{
	struct slList	*pt = (struct slList *)list;
	int		len = 0;

	while(pt != NULL)
	{
		len += 1;
		pt = pt->next;
	}

	return len;
}

static int	dlCount(struct dlList *list)
{
	return slCount(list->head) - 1;
}
*/
/* Return length of list */
static int	dlCount(struct dlList *list)
{
	struct dlNode	*node;
	int		len = 0;


	for(node = list->head; node != NULL; node = node->next)
	{
		len += 1;
	}

	return len;
}


/* Sort a singly linked list with Qsort and a temporary array.
 * The arguments to the compare function in real, non-void, life
 * are pointers to pointers of the type that is in the val field of
 * the nodes of the list.
 */
static void	dlSort(struct dlList *list,
		int (*compare)(const void *elem1, const void *elem2))
{
	int	len = dlCount(list);

	/*
	printf("len=%d\n", len);
	*/

	if(len > 1)
	{
	/* Move val's onto an array, sort, and then back into list. */
		struct dlSorter	*sorter, *s;
		struct dlNode	*node;
		int		i;

		sorter = (struct dlSorter *)memalloc(len * sizeof(sorter[0]));
		for(i = 0, node = list->head; i < len; i ++, node = node->next)
		{
			s = &sorter[i];
			s->node = node;
		}
		compareFunc = compare;
		qsort(sorter, len, sizeof(sorter[0]), dlNodeCmp);
		dlListInit(list);
		for(i = 0; i < len; i ++)
		{
			dlAddTail(list, sorter[i].node);
		}

		freeMem(sorter);
	}
}


/* Determine the cutting dimension with maximum variance (most spread).
 * Return 0, 1, 2 for X, Y, Z dimension respctively.
 */
static int	dimension(struct dlList **lists)
{
	struct kdLeaf	*leafHead, *leafTail;
	double		*d;
	double		max;
	int		i, dim;


	d = (double *)memalloc(sizeof(double) * ndim);
	for(i = 0; i < ndim; i ++)
	{
		leafHead = (struct kdLeaf *)(lists[i]->head->val);
		leafTail = (struct kdLeaf *)(lists[i]->tail->val);
		d[i] = fabs(leafHead->coord[i] - leafTail->coord[i]);
	}

	max = 0.0;
	dim = 0;
	for(i = 0; i < ndim; i ++)
	{
		if(d[i] > max)
		{
			max = d[i];
			dim = i;
		}
	}

	freeMem(d);

	return dim;
}

static struct dlList	**dlListMatrix(int n, int m)
{
	int		i;
	struct dlList	**matrix;

	if((matrix = (struct dlList **)malloc(n * sizeof(struct dlList *))) == NULL)
	{
		printf("malloc(%d * %d) in dlListMatrix(%d, %d) caused error!\n",
			n, sizeof(struct dlList*), n, m);
		exit(1);
	}
	if((matrix[0] = (struct dlList *)malloc(n * m * sizeof(struct dlList))) == NULL)
	{
		printf("malloc(%d * %d * %d) in dlListMatrix(%d, %d) error!\n",
			n, m, sizeof(struct dlList), n, m);
		exit(1);
	}
	for(i = 1; i < n; i ++)
	{
		matrix[i] = matrix[i-1] + m;
	}

	return matrix;
}

static struct dlNode	**dlNodeMatrix(int n, int m)
{
	int		i;
	struct dlNode	**matrix;

	if((matrix = (struct dlNode **)malloc(n * sizeof(struct dlList *))) == NULL)
	{
		printf("malloc(%d * %d) in dlNodeMatrix(%d, %d) caused error!\n",
			n, sizeof(struct dlNode *), n, m);
		exit(1);
	}
	if((matrix[0] = (struct dlNode *)malloc(n * m * sizeof(struct dlNode))) == NULL)
	{
		printf("malloc(%d * %d * %d) in dlNodeMatrix(%d, %d) error!\n",
			n, m, sizeof(struct dlNode), n, m);
		exit(1);
	}
	for(i = 1; i < n; i ++)
	{
		matrix[i] = matrix[i-1] + m;
	}

	return matrix;
}


#ifdef MEDIAN

/* clear hit flags of all leafs in list */
static void	clearHits(struct dlList *list)
{
	struct dlNode	*node;
	struct kdLeaf	*leaf;


	for(node = list->head; node != NULL; node = node->next)
	{
		leaf = node->val;
		leaf->hit = RIGHT;
	}
}

/* peel off members of oldList that not hit onto new list */
static int	splitList(struct dlList *oldList, struct dlList *newList)
{
	struct dlNode	*node, *next;
	struct kdLeaf	*leaf;
	int		newCount = 0;


	dlListInit(newList);
	for(node = oldList->head; node != NULL; node = next)
	{
		next = node->next;
		leaf = node->val;
		if(leaf->hit != LEFT)
		{
			dlRemove(oldList, node);
			
			if(leaf->hit == RIGHT)
			{
				dlAddTail(newList, node);
				newCount ++;
			}
		}
	}

	return newCount;
}

static struct kdLeaf	*medianVal(struct dlList *list, int medianIx)
{
	struct dlNode	*node = list->head;
	struct kdLeaf	*leaf = NULL;
	int		i;

	for(i = 0; i < medianIx; i ++)
	{
		leaf = node->val;
		leaf->hit = LEFT;
		node = node->next;
	}
	leaf = node->val;
	leaf->hit = MID;

	return leaf;
}

/* build up kd-tree recursively */
static struct kdBranch	*kdBuild(int nodeCount, struct dlList **lists)
{
	struct kdBranch	*branch;


	branch = (struct kdBranch *)memalloc(sizeof(*branch));

	if(nodeCount < 1)
	{
		branch = NULL;
	}
	else if(nodeCount == 1)
	{
		struct kdLeaf	*leaf;

		leaf = lists[0]->head->val;

		/* connect leaf to tree */
		branch->leaf = leaf;
		branch->dim = 0;
		branch->lo = branch->hi = NULL;
	}
	else
	{
		int		newCount;
		struct dlList	**newLists;
		int		dim;
		int		i;


		/* subdivide lists along median */
		clearHits(lists[0]);
		dim = dimension(lists);
		branch->dim = dim;

		/* How to deal with the points with same position in qsorted array??? */
		/* It seems that it does not matter!!! */
		branch->leaf = medianVal(lists[dim], nodeCount/2);
		/* How to deal with the points with same position in qsorted array??? */
		/* It seems that it does not matter!!! */
		
		newLists = (struct dlList **)dlListMatrix(ndim, 1);

		newCount = splitList(lists[0], newLists[0]);
		for(i = 1; i < ndim; i ++)
		{
			splitList(lists[i], newLists[i]);
		}

		/* recurse on each side */
		branch->lo = kdBuild(nodeCount - newCount - 1, lists);
		branch->hi = kdBuild(newCount, newLists);

		freeMatrix((void *)newLists);
	}

	return branch;
}

struct kdTree	*kdTreeMake(struct kdLeaf *leafList, int nodeCount, int ndimension)
{
	/*
	struct kdLeaf	*leaf;
	*/
	struct kdTree	*tree;
	struct dlList	*dimList, **lists;
	struct dlNode	**dimNodes;
	int		i, j;
	extern int	dimCmp;
	extern int	ndim;


	/* set value for globe variable */
	ndim = ndimension;

	/* build lists sorted in each dimension. This will let us quickly find medians
	 * while constructing the kd-tree.
	 */
	dimList = (struct dlList *)memalloc(sizeof(*dimList) * ndim);
	for(i = 0; i < ndim; i ++)
	{
		dlListInit(&dimList[i]);
	}
	dimNodes = dlNodeMatrix(nodeCount, ndim);
	/*
	for(i = 0, leaf = leafList; leaf != NULL; leaf = leaf->next, i ++)
	{
		for(j = 0; j < ndim; j ++)
		{
			dimNodes[i][j].val = leaf;
			dlAddTail(&dimList[j], &dimNodes[i][j]);
		}
	}
	*/

	for(i = 0; i < nodeCount; i ++)
	{
		for(j = 0; j < ndim; j ++)
		{
			dimNodes[i][j].val = &leafList[i];
			dlAddTail(&dimList[j], &dimNodes[i][j]);
		}
	}

	lists = (struct dlList **)memalloc(ndim * sizeof(*lists));
	for(i = 0; i < ndim; i ++)
	{
		dimCmp = i;
		dlSort(&dimList[i], kdLeafCmp);
		lists[i] = &dimList[i];
	}

	tree = (struct kdTree *)memalloc(sizeof(*tree));
	tree->root = kdBuild(nodeCount, lists);

	freeMatrix((void **)dimNodes);
	freeMem(dimList);
	freeMem(lists);

	return tree;
}

#else

/* peel off members of oldList that not hit onto new list */
static int	splitList(struct dlList *oldList, struct dlList *newList)
{
	struct dlNode	*node, *next;
	struct kdLeaf	*leaf;
	int		newCount = 0;


	dlListInit(newList);
	for(node = oldList->head; node != NULL; node = next)
	{
		next = node->next;
		leaf = node->val;
		if(!leaf->hit)
		{
			dlRemove(oldList, node);
			dlAddTail(newList, node);
			newCount ++;
		}
	}

	return newCount;
}

/* clear hit flags of all leafs in list */
static void	clearHits(struct dlList *list)
{
	struct dlNode	*node;
	struct kdLeaf	*leaf;


	for(node = list->head; node != NULL; node = node->next)
	{
		leaf = node->val;
		leaf->hit = FALSE;
	}
}

/*
static int	dimension(struct dlList **lists, int nodeCount)
{
	struct kdLeaf	*leafHead, *leafTail;
	double		dx, dy, dz;
	double		max;
	int		dim;

	leafHead = (struct kdLeaf *)(lists[0]->head->val);
	leafTail = (struct kdLeaf *)(lists[0]->tail->val);
	dx = fabs(leafHead->coord[X] - leafTail->coord[X]);
	leafHead = (struct kdLeaf *)(lists[1]->head->val);
	leafTail = (struct kdLeaf *)(lists[1]->tail->val);
	dy = fabs(leafHead->coord[Y] - leafTail->coord[Y]);
	leafHead = (struct kdLeaf *)(lists[2]->head->val);
	leafTail = (struct kdLeaf *)(lists[2]->tail->val);
	dz = fabs(leafHead->coord[Z] - leafTail->coord[Z]);

	max = 0.0;
	dim = X;
	if(dx > dy)
	{
		max = dx;
		dim = X;
	}
	else
	{
		max = dy;
		dim = Y;
	}
	if(dz > max)
	{
		dim = Z;
	}

	printf("dim=%d dx %f dy %f dz %f nodeCount=%d\n", dim, dx, dy, dz, nodeCount);

	return dim;
}
*/

/* Return coordinate value of median point (pivot) on given dimension
 * Mark points up to median as hit.
 */
static double	medianVal(struct dlList *list, int medianIx, int dim)
{
	struct dlNode	*node = list->head;
	struct kdLeaf	*leaf = NULL;
	int		i;

	for(i = 0; i < medianIx; i ++)
	{
		leaf = node->val;
		leaf->hit = TRUE;
		node = node->next;
	}

	/* median point belongs to other child */
	leaf = node->val;
	if(dim == X)
		return leaf->coord[X];
	else if(dim == Y)
		return leaf->coord[Y];
	else	return leaf->coord[Z];
}

static struct kdBranch	*kdBuild(int nodeCount, struct dlList **lists)
{
	struct kdBranch	*branch;


	branch = (struct kdBranch *)memalloc(sizeof(*branch));

	if(nodeCount == 1)
	{
		struct kdLeaf	*leaf;

		leaf = lists[0]->head->val;

		/* connect leaf to tree */
		branch->leaf = leaf;
		branch->dim = INT_MAX;
		branch->cutCoord = (double)INT_MAX;
		branch->lo = branch->hi = NULL;
	}
	else
	{
		int		newCount;
		struct dlList	*newLists[NDIM];
		struct dlList	newX, newY, newZ;
		int		dim;


		/* subdivide lists along median */
		newLists[X] = &newX;
		newLists[Y] = &newY;
		newLists[Z] = &newZ;
		clearHits(lists[0]);
		/*
		dim = dimension(lists, nodeCount);
		*/
		dim = dimension(lists);
		branch->dim = dim;
		branch->leaf = NULL;
		branch->cutCoord = medianVal(lists[dim], nodeCount/2, dim); 

		/* How to deal with the points with same position in qsorted array??? */
		/* It seems that it does not matter!!! */
		/*
		branch->leaf = medianVal(lists[dim], nodeCount/2, dim);
		*/
		/* How to deal with the points with same position in qsorted array??? */
		/* It seems that it does not matter!!! */
		

		newCount = splitList(lists[X], newLists[X]);
		splitList(lists[Y], newLists[Y]);
		splitList(lists[Z], newLists[Z]);

		/* recurse on each side */
		branch->lo = kdBuild(nodeCount - newCount, lists);
		branch->hi = kdBuild(newCount, newLists);
	}

	return branch;
}

struct kdTree	*kdTreeMake(struct kdLeaf *leafList, int nodeCount, int ndimension)
{
	/*
	struct kdLeaf	*leaf;
	*/
	struct kdTree	*tree;
	struct dlList	xList, yList, zList, *lists[3];
	struct dlNode	*xNodes, *yNodes, *zNodes;
	int		i;
	extern int	dimCmp;
	extern int	ndim;


	/* set value for globe variable */
	ndim = ndimension;

	/* build lists sorted in each dimension. This will let us quickly find medians
	 * while constructing the kd-tree.
	 */
	dlListInit(&xList);
	dlListInit(&yList);
	dlListInit(&zList);
	xNodes = (struct dlNode *)memalloc(nodeCount * sizeof(*xNodes));
	yNodes = (struct dlNode *)memalloc(nodeCount * sizeof(*yNodes));
	zNodes = (struct dlNode *)memalloc(nodeCount * sizeof(*zNodes));
	/*
	for(i = 0, leaf = leafList; leaf != NULL; leaf = leaf->next, i ++)
	{
		xNodes[i].val = yNodes[i].val = zNodes[i].val = leaf;
		dlAddTail(&xList, &xNodes[i]);
		dlAddTail(&yList, &yNodes[i]);
		dlAddTail(&zList, &zNodes[i]);
	}
	*/
	for(i = 0; i < nodeCount; i ++)
	{
		xNodes[i].val = yNodes[i].val = zNodes[i].val = &leafList[i];
		dlAddTail(&xList, &xNodes[i]);
		dlAddTail(&yList, &yNodes[i]);
		dlAddTail(&zList, &zNodes[i]);
	}

	dimCmp = X;
	dlSort(&xList, kdLeafCmp);
	dimCmp = Y;
	dlSort(&yList, kdLeafCmp);
	dimCmp = Z;
	dlSort(&zList, kdLeafCmp);

	lists[X] = &xList;
	lists[Y] = &yList;
	lists[Z] = &zList;

	tree = (struct kdTree *)memalloc(sizeof(*tree));
	tree->root = kdBuild(nodeCount, lists);

	freeMem(xNodes);
	freeMem(yNodes);
	freeMem(zNodes);

	return tree;
}
#endif


static void	freeBranch(struct kdBranch *b)
{
	if(b->lo != NULL)
		freeBranch(b->lo);
	if(b->hi != NULL)
		freeBranch(b->hi);

	freeMem(b);
}

void	freeKdTree(struct kdTree *tree)
{
	if(tree->root != NULL)
		freeBranch(tree->root);

	freeMem(tree);
}


static void	visitLeaf(struct kdLeaf *leaf)
{
	int	i;

	if(leaf == NULL)	return;

	for(i = 0; i < ndim; i ++)
	{
		printf("%.3f ", leaf->coord[i]);
	}
}

#ifdef MEDIAN
static void	visitBranch(struct kdBranch *branch)
{
	if(branch == NULL)	return;

	printf("(");

	printf("%d ", branch->dim);
	visitLeaf(branch->leaf);
	visitBranch(branch->lo);
	visitBranch(branch->hi);
	
	printf(")");
}

#else
static void	*visitBranch(struct kdBranch *branch)
{
	void	*p1, *p2;


	if(branch == NULL)	return NULL;

	printf("(");

	p1 = visitBranch(branch->lo);
	p2 = visitBranch(branch->hi);	
	
	if(p1 == NULL && p2 == NULL)
	{
		visitLeaf(branch->leaf);
	}
	else
	{
		printf("dim %d cutCoord %.3f", branch->dim, branch->cutCoord);
	}

	printf(")");

	return branch;
}
#endif

void	traverseKdTree(struct kdTree *tree)
{
	if(tree == NULL)
	{
		printf("Empty tree!\n");
	}
	else
	{
		visitBranch(tree->root);
		printf("\n");
	}
}


static void	cutHyperrectangle(double **hr, double cutCoord, int dim, double **hrLo, double **hrHi)
{
	int	i;

	for(i = 0; i < ndim; i ++)
	{
		if(i == dim)
		{
			hrLo[1][i] = hrHi[0][i] = cutCoord;
		}
		else
		{
			hrLo[1][i] = hr[1][i];
			hrHi[0][i] = hr[0][i];
		}
	}

	/* splitting hyperplane does not intersect with hyperrectangle. */
	if(cutCoord < hr[0][dim])
	{
		hrHi[0][dim] = hr[0][dim];
	}
	if(cutCoord > hr[1][dim])
	{
		hrLo[1][dim] = hr[1][dim];
	}
}

/* Return TRUE if input is ilegal hyperrectangle, else FALSE */
static int	ilegalHyperrectangle(double *hr[2], int dim)
{
	if(hr[0] == NULL || hr[1] == NULL)
		return TRUE;

	if(hr[0][dim] > hr[1][dim] || hr[1][dim] < hr[0][dim])
		return TRUE;
	else	return FALSE;
}

/* Return TRUE if input is out of space of hyperrectangle, else FALSE */
static int	outOfRectangle(double *pos, double *hr[2])
{
	int	i;

	for(i = 0; i < ndim; i ++)
	{
		if(pos[i] < hr[0][i] || pos[i] > hr[1][i])
		{
			return TRUE;
		}
	}

	return FALSE;
}


/* Return TRUE if hyperrectangle intersects with hypersphere defined by "target" and "distMax", else FALSE.
 * To check to see if a hyperrectangle "hr" intersects with a hypersphere radius "r" centered at point
 * "t", we find the point "p" in "hr" which is closest to "t". Write "hr_min_i" as the minimum extreme of
 * "hr" in the "ith" dimension and "hr_max_i" as the maximum extrem. "p_i", the "ith" component of this
 * closest point is computed thus:
 *
 *	  {  hr_min_i	if(t_i <= hr_min_i)
 * p_i = <|  t_i	if(hr_min_i < t_i < hr_max_i)
 *	  {  hr_max_i	if(t_i >= hr_max_i) 
 *
 * The objects intersect only if the distance between "p" and "t" is less than or equal to "r".
 * Refer to "Andrew W. Moore. An introductory tutorial on kd-trees. Extract from Andrew Moores's
 * PhD Thesis: Efficient Memory-based Learning for Robot Control, PhD. Thesis. 1991".
*/
static int	hrIntersectSphere(struct kdLeaf *target,
				double *hr[2],
				double distMax)
{
	int	i;
	double	*p, *t;
	double	*d, dist;

	p = (double *)memalloc(sizeof(double) * ndim);
	t = target->coord;
	for(i = 0; i < ndim; i ++)
	{
		if(t[i] <= hr[0][i])
		{
			p[i] = hr[0][i];
		}
		else if(t[i] >= hr[1][i])
		{
			p[i] = hr[1][i];
		}
		else	p[i] = t[i];
	}

	d = (double *)memalloc(sizeof(double) * ndim);
	for(i = 0; i < ndim; i ++)
	{
		d[i] = p[i] - t[i];
	}
	freeMem(p);
	dist = 0.0;
	for(i = 0; i < ndim; i ++)
	{
		dist += d[i] * d[i];
	}
	freeMem(d);

	if(dist > distMax)
		return FALSE;
	else	return TRUE;
}

/* Return maximum value of distance in nearList */
static double	maxDist(struct dlList *nearList)
{
	struct dlNode	*node;
	double		max;

	node = nearList->head;
	if(node == NULL)
	{
		return 0.0;
	}

	max = node->dist;
	for(node = node->next; node != NULL; node = node->next)
	{
		if(node->dist > max)
			max = node->dist;
	}

	return max;
}

/* Remove "node" and all nodes after "node" */
/*
static void	dlRemoveTail(struct dlList *list, struct dlNode *node)
{
	list->tail = node->prev;
	node->prev->next = NULL;
	node->prev = NULL;
}
*/

void	freeNodeList(struct dlNode *node)
{
	struct dlNode	*nextnode;

	for(; node != NULL; node = nextnode)
	{
		nextnode = node->next;
		freeMem((void *)node);
	}
}

static void	dlAddNode(struct dlList *list, struct kdLeaf *leaf, double dist)
{
	struct dlNode	*node;

	node = (struct dlNode *)memalloc(sizeof(*node));
	node->val = leaf;
	node->dist = dist;
	dlAddTail(list, node);
}

static void	dlReplaceNode(struct dlList *nearList, int *nearCount, int nNear,
		struct kdLeaf *leaf, double dist, double distMax)
{
	int		distMaxCount;
	struct dlNode	*node;

	distMaxCount = 0;
	for(node = nearList->head; node != NULL; node = node->next)
	{
		if(!(node->dist < distMax || node->dist > distMax))
		{
			distMaxCount ++;
		}
	}

	if(*nearCount - distMaxCount + 1 < nNear)
	{
		dlAddNode(nearList, leaf, dist);
		(*nearCount) ++;
	}
	else
	/* Remove all (distMaxCount) invalid (distance == distMaxt) nodes */
	{
		/* Remove (distMaxCount-1) nodes */
		int	count = 0;
		for(node = nearList->head; node != NULL; node = node->next)
		{
			if(!(node->dist < distMax || node->dist > distMax))
			{
				count ++;
				if(count < distMaxCount)
				{
					dlRemove(nearList, node);
					freeMem((void *)node);
				}
				else	break;
			}
		}

		/* Replace one invalid node */
		if(node != NULL)
		{
			node->val = leaf;
			node->dist = dist;
			*nearCount -= count - 1;
		}

		/* check number of valid nodes */
		if(*nearCount != nNear)
		{
			printf("#Error! nearCount != nNear (nearCount %d count %d distMaxCount %d nNear %d)\n",
				*nearCount, count, distMaxCount, nNear);
			exit(1);
		}
	}
}

/* Examine leaf in kdtree */
static double	checkLeaf(
			struct kdLeaf *leaf,
			struct kdLeaf *target,
			double *hr[2],
			double distMax,
			int *nearCount,
			int nNear,
			struct dlList *nearList)
{
	double		*pos, *dpos;
	int		i;
	double		dist;

	pos = leaf->coord;

	/* check whether point lies within hyperrectangle. */
	if(outOfRectangle(pos, hr))
		return (double)INT_MAX;

	dpos = (double *)memalloc(sizeof(double) * ndim);
	for(i = 0; i < ndim; i ++)
	{
		dpos[i] = pos[i] - target->coord[i];
	}
	dist = 0.0;
	for(i = 0; i < ndim; i ++)
	{
		dist += dpos[i] * dpos[i];
	}
	freeMem(dpos);

#ifdef DEBUG
	computationCount ++;
	printf("computationCount=%d###\n", computationCount);
#endif

	/* check distance cutoff */
	if(!(dist > distMax))
	{
		if(*nearCount < nNear)
		{
			dlAddNode(nearList, leaf, dist);
			(*nearCount) ++;
		}
		else if(dist < distMax && nNear > 0)
		{
			dlReplaceNode(nearList, nearCount, nNear, leaf, dist, distMax);
		}
		else
		{
			dlAddNode(nearList, leaf, dist);
			(*nearCount) ++;
		}
	}

	return min(dist, distMax);
}


#ifdef MEDIAN
/*
 * Return distance between near point and target as return value and all near points in "nearList".
 *
 * *branch:	kd-tree branch to be searched
 * *target:	speficied target point
 * *hr[2]:	hr[0] and hr[1] for minimum and maximum of coordinates on each dimension of hyperrectangle,
 *		respectively
 * distMax:	cutoff for squre of distance, program searches near points within distance distMax of target
 * *nearList:	nearList points found by program, linked list
 */
static double	nearNeighborSearch(
				struct kdBranch *branch,
				struct kdLeaf *target,
				double *hr[2],
				double distMax,
				int *nearCount,
				int nNear,
				struct dlList *nearList)
{
	extern int	ndim;

	int	dim;
	double	cutCoord;
	double	dist;


	/* ilegal kd-tree or hyperrectangle */
	if(branch == NULL || hr == NULL)
	{
		return (double)INT_MAX;
	}

	/* test the kd-tree */
	printf("\n");
	visitBranch(branch);
	printf("\n");
	/* test the kd-tree */

	dim = branch->dim;
	cutCoord = branch->leaf->coord[dim];

	if(branch->lo == branch->hi && branch->hi == NULL)
	{
		dist = checkLeaf(branch->leaf, target, hr, distMax, nearCount, nNear, nearList);
	}
	else
	{
		double		*hrLo[2];
		double		*hrHi[2];
		double		**hrNearer, **hrFurther;
		struct kdBranch	*branchNearer, *branchFurther;
		int		i;


		hrLo[0] = hr[0];
		hrHi[1] = hr[1];
		hrLo[1] = (double *)memalloc(sizeof(double) * ndim);
		hrHi[0] = (double *)memalloc(sizeof(double) * ndim);

		cutHyperrectangle(hr, cutCoord, dim, hrLo, hrHi);

		/* print parent and sub-hyperrectangles */
		printf("dim=%d cutCoord=%f\n", dim, cutCoord);
		for(i = 0; i < ndim; i ++)
		{
			printf("hr %f %f\n", hr[0][i], hr[1][i]);
		}
		for(i = 0; i < ndim; i ++)
		{
			printf("hrLo %f %f\n", hrLo[0][i], hrLo[1][i]);
		}
		for(i = 0; i < ndim; i ++)
		{
			printf("hrHi %f %f\n", hrHi[0][i], hrHi[1][i]);
		}
		/* print parent and child hyperrectangles */

		if(target->coord[dim] < cutCoord)
		{
			branchNearer = branch->lo;
			hrNearer = hrLo;
			branchFurther = branch->hi;
			hrFurther = hrHi;
		}
		else
		{
			branchNearer = branch->hi;
			hrNearer = hrHi;
			branchFurther = branch->lo;
			hrFurther = hrLo;
		}

		if(!ilegalHyperrectangle(hrNearer, dim))
		{
			dist = nearNeighborSearch(branchNearer, target, hrNearer, distMax,
					nearCount, nNear, nearList);
		}

		if(*nearCount >= nNear && nNear > 0)
		{
			distMax = maxDist(nearList);
		}

		if(!ilegalHyperrectangle(hrFurther, dim) && hrIntersectSphere(target, hrFurther, distMax))
		{
			/* check median point (node) through the splitting hyperplane */
			dist = checkLeaf(branch->leaf, target, hrFurther, distMax, nearCount, nNear, nearList);

			if(*nearCount >= nNear && nNear > 0)
			{
				distMax = maxDist(nearList);
			}

			dist = nearNeighborSearch(branchFurther, target, hrFurther, distMax,
					nearCount, nNear, nearList);
		}

		freeMem(hrLo[1]);
		freeMem(hrHi[0]);
	}

	return min(dist, distMax);
}

#else
/*
 * Return distance between near point and target as return value, all near points and their distances
 * to target are stored in "nearList".
 *
 * *branch:	kd-tree branch to be searched
 * *target:	speficied target point
 * *hr[2]:	hr[0] and hr[1] for minimum and maximum of coordinates on each dimension of hyperrectangle,
 *		respectively
 * distMax:	cutoff for squre of distance, program searches near points within distance distMax of target
 * *nearList:	nearList points found by program, linked list
 */
static double	nearNeighborSearch(
				struct kdBranch *branch,
				struct kdLeaf *target,
				double *hr[2],
				double distMax,
				int *nearCount,
				int nNear,
				struct dlList *nearList)
{
	extern int	ndim;

	int	dim = branch->dim;
	double	cutCoord = branch->cutCoord;
	double	dist;


	/* ilegal kd-tree or hyperrectangle */
	if(branch == NULL || hr == NULL)
	{
		return (double)INT_MAX;
	}

	/* test the kd-tree */
	printf("\n");
	visitBranch(branch);
	printf("\n");
	/* test the kd-tree */

	if(branch->lo == branch->hi && branch->hi == NULL)
	{
		dist = checkLeaf(branch->leaf, target, hr, distMax, nearCount, nNear, nearList);
	}
	else
	{
		double		*hrLo[2];
		double		*hrHi[2];
		double		**hrNearer, **hrFurther;
		struct kdBranch	*branchNearer, *branchFurther;


		hrLo[0] = hr[0];
		hrHi[1] = hr[1];
		hrLo[1] = (double *)memalloc(sizeof(double) * ndim);
		hrHi[0] = (double *)memalloc(sizeof(double) * ndim);

		cutHyperrectangle(hr, cutCoord, dim, hrLo, hrHi);

		/* print parent and sub-hyperrectangles */
		printf("dim=%d cutCoord=%f\n", dim, cutCoord);
		for(int i = 0; i < ndim; i ++)
		{
			printf("hr %f %f\n", hr[0][i], hr[1][i]);
		}
		for(int i = 0; i < ndim; i ++)
		{
			printf("hrLo %f %f\n", hrLo[0][i], hrLo[1][i]);
		}
		for(int i = 0; i < ndim; i ++)
		{
			printf("hrHi %f %f\n", hrHi[0][i], hrHi[1][i]);
		}
		/* print parent and child hyperrectangles */

		if(target->coord[dim] < cutCoord)
		{
			branchNearer = branch->lo;
			hrNearer = hrLo;
			branchFurther = branch->hi;
			hrFurther = hrHi;
		}
		else
		{
			branchNearer = branch->hi;
			hrNearer = hrHi;
			branchFurther = branch->lo;
			hrFurther = hrLo;
		}

		if(!ilegalHyperrectangle(hrNearer, dim))
		{
			dist = nearNeighborSearch(branchNearer, target, hrNearer, distMax,
					nearCount, nNear, nearList);
		}

		if(*nearCount >= nNear && nNear > 0)
		{
			distMax = maxDist(nearList);
		}

		if(!ilegalHyperrectangle(hrFurther, dim) && hrIntersectSphere(target, hrFurther, distMax))
		{
			dist = nearNeighborSearch(branchFurther, target, hrFurther, distMax,
					nearCount, nNear, nearList);
		}

		freeMem(hrLo[1]);
		freeMem(hrHi[0]);
	}

	return min(dist, distMax);
}
#endif


void	outPutNearList(struct dlList *list)
{
	struct dlNode	*node;
	struct kdLeaf	*leaf;

	printf("near:\n");
	for(node = list->head; node != NULL; node = node->next)
	{
		leaf = (struct kdLeaf *)(node->val);
		printf("%f %f %f dist-square %f\n", leaf->coord[0], leaf->coord[1], leaf->coord[2], node->dist);
	}
}

void	nearSearch_(struct kdTree *tree, struct kdLeaf *target, int ndimension, double *hr[2], double distMax,
		int nNear,
		struct dlList *nearList)
{
	extern int	ndim;

	int	nearCount = 0;

	ndim = ndimension;

	nearNeighborSearch(tree->root, target, hr, distMax, &nearCount, nNear, nearList);
}

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
void	nearSearch(double **data, int npoint, int ndimension,
		double *target, double *hr[2], double distMax,
		int nNear,
		struct dlList *nearList)
{
	extern int	ndim;
	int		i, j;
	struct kdLeaf	*leaf, targetLeaf;
	struct kdTree	*kdTree;
	int		nearCount = 0;


	ndim = ndimension;

	if(ndim > NDIM)
	{
		printf("#Error! Too many dimenstions (%d) in data space, "
			"please increase the value of NDIM in kdtree.h\n", NDIM);
		exit(1);
	}

	leaf = (struct kdLeaf *)memalloc(sizeof(*leaf) * npoint);
	for(i = 0; i < npoint; i ++)
	{
		for(j = 0; j < ndim; j ++)
		{
			leaf[i].coord[j] = data[i][j];
		}
	}

	kdTree = kdTreeMake(leaf, npoint, ndim);

#ifdef DEBUG
	traverseKdTree(kdTree);
#endif

	dlListInit(nearList);

	for(i = 0; i < ndim; i ++)
	{
		targetLeaf.coord[i] = target[i];
	}

	nearNeighborSearch(kdTree->root, &targetLeaf, hr, distMax, &nearCount, nNear, nearList);

#ifndef DEBUG
	outPutNearList(nearList);
	freeNodeList(nearList->head);
	freeKdTree(kdTree);
#endif
}


#ifdef DEBUG
int	main()
{
	struct kdLeaf	leaf[4];
	int		i;
	int		npoint;
	struct kdTree	*tree;

	/*
	struct kdLeaf	target;
	*/
	double		*hr[2];
	extern int	ndim;
	double		distMax; /* sqaure of distance cutoff */
	struct dlList	nearList;

	npoint = 4;
	ndim = 3;

/*
	leaf[0].coord[0] = 2; leaf[0].coord[1] = 5; leaf[0].coord[2] = 2.5;
	leaf[1].coord[0] = 2; leaf[1].coord[1] = 5; leaf[1].coord[2] = 2.5;
*/
	leaf[0].coord[0] = 8; leaf[0].coord[1] = 8; leaf[0].coord[2] = 3.8;
	leaf[1].coord[0] = 5; leaf[1].coord[1] = 3; leaf[1].coord[2] = 6.3;
	leaf[2].coord[0] = 8; leaf[2].coord[1] = 8; leaf[2].coord[2] = 3.8;
/*
	leaf[3].coord[0] = 2; leaf[3].coord[1] = 5; leaf[3].coord[2] = 2.5;
*/
	leaf[3].coord[0] = 9; leaf[3].coord[1] = 9; leaf[3].coord[2] = 8.9;

	for(i = 0; i < 3; i ++)
	{
		leaf[i].next = &(leaf[i+1]);

		printf("i=%d %8.3f %8.3f %8.3f leaf=%x next=%x\n",
			i, leaf[i].coord[0], leaf[i].coord[1], leaf[i].coord[2], &leaf[i], leaf[i].next);
	}
	leaf[i].next = NULL;

	printf("i=%d %8.3f %8.3f %8.3f leaf=%x next=%x\n",
		i, leaf[i].coord[0], leaf[i].coord[1], leaf[i].coord[2], &leaf[i], leaf[i].next);

	tree = kdTreeMake(leaf, npoint, ndim);

	traverseKdTree(tree);


	/*
	target.coord[0] = 2;
	target.coord[1] = 5;
	target.coord[2] = 2;
	*/

	hr[0] = (double *)memalloc(sizeof(double) * ndim);
	hr[1] = (double *)memalloc(sizeof(double) * ndim);
	printf("hyperrectangle:\n");
	for(i = 0; i < ndim; i ++)
	{
		hr[0][i] = -99.0;
		hr[1][i] = 99.0;
		printf("%f %f\n", hr[0][i], hr[1][i]);
	}

	distMax = 969.0;
	dlListInit(&nearList);

	{
	int	nNear;
	double	**data;
	double	target[10];
	int	j;

	npoint = 4;
	data = dmatrix(10, 10);
	for(i = 0; i < npoint; i ++)
	{
		for(j = 0; j < ndim; j ++)
		{
			data[i][j] = leaf[i].coord[j];
		}
	}
	target[0] = 20;
	target[1] = 5;
	target[2] = 2;

	nNear = 0;
	/*
	nearSearch_(tree, &target, ndim, hr, distMax, nNear, &nearList);
	*/
	nearSearch(data, npoint, ndim, target, hr, distMax, nNear, &nearList);

	printf(	"target:\n"
		"%f %f %f distMax %f\n",
		target[0], target[1], target[2], distMax);

	free_dmatrix(data);
	}

	freeMem(hr[0]);
	freeMem(hr[1]);

	/*
	printf(	"target:\n"
		"%f %f %f distMax %f\n",
		target.coord[0], target.coord[1], target.coord[2], distMax);
	*/

	outPutNearList(&nearList);

	freeNodeList(nearList.head);
	freeKdTree(tree);

	return 0;
}
#endif
