/*
#include <stdinc.h>
*/


#include "cnctarea2.h"
/* added by Zhiqun Xie */


void sortag(double *a, int n, int *tag)
{
      int tg, i, j, ij, k, m, l, il[16], iu[16];
      double t,tt;

      for(i = 0; i < n; i ++)
        tag[i] = i + 1;
      m=1;
      i=1;
      j=n;

label5:      if (i >= j)	 goto label70;

label10:   k=i;
      ij=(j+i)/2;
      t=a[ij - 1];
      if (a[i-1] <= t)	 goto label20;

      a[ij - 1]= a[i - 1];
      a[i - 1]=t;
      t=a[ij - 1];
      tg=tag[ij - 1];
      tag[ij - 1]=tag[i - 1];
      tag[i - 1]=tg;
label20:   l=j;
      if (a[j - 1] >= t) goto label40;

      a[ij - 1]=a[j - 1];
      a[j - 1]=t;
      t=a[ij - 1];
      tg=tag[ij - 1];
      tag[ij - 1]=tag[j - 1];
      tag[j - 1]=tg;
      if (a[i-1] <= t)	goto label40;

      a[ij - 1]=a[i - 1];
      a[i - 1]=t;
      t=a[ij - 1];
      tg=tag[ij - 1];
      tag[ij - 1]=tag[i - 1];
      tag[i - 1]=tg;
      goto label40;

label30:   a[l - 1]=a[k - 1];
      a[k - 1]=tt;
      tg=tag[l - 1];
      tag[l - 1]=tag[k - 1];
      tag[k - 1]=tg;

label40:   l=l-1;
      if (a[l - 1] > t) goto label40;
      tt=a[l - 1];

label50:   k=k+1;
      if (a[k - 1] < t) goto label50;
      if (k <= l) goto label30;
      if (l-i <= j-k) goto label60;
      il[m - 1]=i;
      iu[m - 1]=l;
      i=k;
      m=m+1;
      goto label80;

label60:   il[m - 1]=k;
      iu[m - 1]=j;
      j=l;
      m=m+1;
      goto label80;

label70:   m=m-1;
      if (m == 0) return;
      i=il[m - 1];
      j=iu[m - 1];

label80:   if (j-i >= 1) goto label10;
      if (i == 1) goto label5;
      i=i-1;

label90:   i=i+1;
      if (i == j) goto label70;
      t=a[i];
      if (a[i-1] <= t) goto label90;
      tg=tag[i];
      k=i;

label100:  a[k]=a[k - 1];
      tag[k]=tag[k - 1];
      k=k-1;
      if (t < a[k - 1]) goto label100;
      a[k]=t;
      tag[k]=tg;
      goto label90;
}
