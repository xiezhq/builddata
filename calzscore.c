#include <math.h>
#include "memory.h"
#include "calzscore.h"


void	calzscore(double *data, int n, double *zscore)
{
	double	sum1;
	double	mean;
	double	diff;
	double	sum2;
	
	double	stddev;
	int	i;


	sum1 = 0.0;
	sum2 = 0.0;

	for(i = 0; i < n; i ++)
	{
		sum1 += data[i];
	}

	mean = sum1 / ((double) n);

	for(i = 0; i < n; i ++)
	{
		diff = data[i] - mean; 
		sum2 += diff * diff;
	}

	stddev = sqrt(sum2 / ((double ) n));

	for(i = 0; i < n; i ++)
	{
		zscore[i] = (data[i] - mean) / stddev;
	}
}
