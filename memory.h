#ifndef MEMORY_H
#define MEMORY_H


extern void	*memalloc(int n);

extern void	freeMem(void *p);

extern void	freeMatrix(void **p);

extern void	free_fmatrix(float **matrix);

extern void	free_dmatrix(double **matrix);

extern void	free_chrmatrix(char **matrix);

extern void	free_imatrix(int **matrix);

extern void	**matrix(int n, int m, int size);

extern int	**imatrix(int n, int m);

extern float	**fmatrix(int n, int m);

extern double	**dmatrix(int n, int m);

extern char	**chrmatrix(int n, int m);

#endif
