#include <stdio.h>
#include <stdlib.h>
#include "memory.h"


void	*memalloc(int n)
{
	void	*p=NULL;

	if((p = malloc(n)) == NULL)
	{
		perror("malloc");
		exit(0);
	}
	return p;
}

void	freeMem(void *p)
{
	if(p != NULL)
	{
		free(p);
	}
}

void	freeMatrix(void **matrix)
{
	freeMem(matrix[0]);
	freeMem(matrix);
}

void	free_fmatrix(float **matrix)
{
	freeMem(matrix[0]);
	freeMem(matrix);
}

void	free_dmatrix(double **matrix)
{
	freeMem(matrix[0]);
	freeMem(matrix);
}


void	free_chrmatrix(char **matrix)
{
	freeMem(matrix[0]);
	freeMem(matrix);
}

void	free_imatrix(int **matrix)
{
	freeMem(matrix[0]);
	freeMem(matrix);
}

int	**imatrix(int n, int m)
{
	int	i;
	int	**matrix;

	if((matrix = (int **)malloc(n * sizeof(int *))) == NULL)
	{
		perror("malloc(n * sizeof(int *))");
		exit(1);
	}
	if((matrix[0] = (int *)malloc(n * m * sizeof(int))) == NULL)
	{
		perror("malloc(n * m * sizeof(int))");
		exit(1);
	}
	for(i = 1; i < n; i++)
		matrix[i] = matrix[i-1] + m;

	return matrix;
}

float	**fmatrix(int n, int m)
{
	int	i;
	float	**matrix;

	if((matrix = (float **)malloc(n * sizeof(float *))) == NULL)
	{
		perror("malloc");
		exit(0);
	}
	if((matrix[0] = (float *)malloc(n * m * sizeof(float))) == NULL)
	{
		perror("malloc");
		exit(0);
	}
	for(i = 1; i < n; i++)
		matrix[i] = matrix[i-1] + m;

	return matrix;
}

double	**dmatrix(int n, int m)
{
	int	i;
	double		**matrix;

	if((matrix = (double **)malloc(n * sizeof(double *))) == NULL)
	{
		perror("malloc");
		exit(1);
	}
	if((matrix[0] = (double *)malloc(n * m * sizeof(double))) == NULL)
	{
		perror("malloc");
		exit(1);
	}
	for(i = 1; i < n; i++)
		matrix[i] = matrix[i-1] + m;

	return matrix;
}


char	**chrmatrix(int n, int m)
{
	int	i;
	char	**matrix;

	if((matrix = (char **)malloc(n * sizeof(char *))) == NULL)
	{
		perror("malloc");
		exit(1);
	}
	if((matrix[0] = (char *)malloc(n * m * sizeof(char))) == NULL)
	{
		perror("malloc");
		exit(1);
	}
	for(i = 1; i < n; i ++)
		matrix[i] = matrix[i-1] + m;

	return matrix;
}
