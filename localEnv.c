#include <stdio.h>
#include <stdlib.h>
#include <math.h>


/* arithmetic mean */
double	arithMean(double *data, int n)
{
	int	i;
	double	mean;

	
	mean = 0.0;
	for(i = 0; i < n; i < ++)
	{
		mean += data[i];
	}

	return (mean / n);
}

/* arithmetic mean for multiple dimensional data */
void	arithMeanMultipleDim(double **data, int nPoint, int nDim, double *mean)
{
	int	i, j;

	for(i = 0; i < nDim; i ++)
	{
		mean[i] = 0.0;
		for(j = 0; j < nPoint; j ++)
		{
			mean[i] += data[i][j];
		}
		mean[i] /= nPoint;
	}
}

/* geometric mean */
double	geomMean(double *data, int n)
{
	int	i;
	double	mean;


	mean = data[0];
	for(i = 1; i < n; i ++)
	{
		mean *= data[i];
	}

	return pow(mean, 1.0/n);
}

/* geometric mean for multiple dimensional data */
void	geomMeanMultipleDim(double **data, int nPoint, int nDim, double *mean)
{
	int	i, j;

	for(i = 0; i < nDim; i ++)
	{
		mean[i] = data[i][j];
		for(j = 1; j < nPoint; j ++)
		{
			mean[i] *= data[i][j];
		}
		mean[i] = pow(mean[i], 1.0/nPoint);
	}
}

/* squre of standard deviation */
double	stdDev(double *data, int n)
{
	int	i;
	double	mean;
	double	s;
	double	tmp;

	mean = arithMean(data, n);
	s = n * mean * mean;
	tmp = 0.0;
	for(i = 0; i < n; i ++)
	{
		s += data[i] * data[i];
		tmp += data[i];
	}
	s -= 2 * mean * tmp;

	return s / (n-1);
}

/* squre of standard deviation */
double	stdDev_(double *data, int n)
{
	int	i;
	double	mean;
	double	d;
	double	s;

	mean = arithMean(data, n);
	s = 0.0;
	for(i = 0; i < n; i ++)
	{
		d = data[i] - mean;
		s += d * d;
	}

	return s / (n-1);
}


/* Pearson correlation coefficient between two vectors (variables) */
double	correlation(double *x, double *y, int n)
{
	int	i;
	double	xmean, ymean, xs, ys;
	double	r;

	xmean = arithMean(x, n);
	ymean = arithMean(y, n);
	r = 0.0;
	for(i = 0; i < n; i ++)
	{
		r += (x[i] - xmean) * (y[i] - ymean);
	}
	xs = stdDev_(x, n);
	ys = stdDev_(y, n);

	return r / ((n-1) * xs * ys);
}
