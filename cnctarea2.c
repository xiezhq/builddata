/***********************************************************************
*
* a paralleliped spanned by the molecule is divided into cubes of 2*rmax,
* where rmax is max (rad atom + water radius):
*
***********************************************************************/






/*-------------------
The comments below are added by Zhiqun Xie in 2005.
---------------------
extern void cnctarea(double errorpar, double *probe, double *x, double *y, double *z,
		      double *radx, int natms, double *accss);

errorpar	z-spacing factor (slice thickness), default errorpar = 0.05 angstrom
probe		probe = atomic radius + radius of solvent molecule, default = Rvdw + Rwater
x, y, z		coordinates for each atom
radx		van der Walls radius of atom
natms		number of atoms in the given protein unit
accss		accessible area for each atom



Note: Lee & Rechards algorithm, refer to:

Lee B & Richards FM. The interpretation of protein structures: estimation of static accessibility.
J Mol Biol, 1971, 55:379-400.

readme.txt of NACCESS.
Hubbard SJ & Thornton JM. NACCESS. Department of Biochemistry and Molecular Biology,
University College London, 1993.

J Biol Chem, 2004, 279(22):23061-23072 for the explanation of errorpar.

-------------------
The comments above are added by Zhiqun Xie in 2005.
-----------------*/







/******************
#include <stdinc.h>
#include <extfunc.h>
*******************/
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

/*
#define maxcubes 3000
*/
#define maxcubes 8000 /* increase from 3000 to 8000 by xie in November, 2005*/

#define maxintersect 80000
#define maxatomscube 800
#define pi 3.141592653
#define pix2 2.0*pi
#define water_radii 1.4


/*
void cnctarea(double errorpar, double *probe, double *x, double *y, double *z,
	      double *radx, int natms, double *accss);
*/

void cnctarea(double errorpar, double *probe, double *x, double *y, double *z,
	      double *radx, int natms, double *accss)
{
	int    idim,jidim,kjidim,i,j,k,l,kji,n,ir,io,kk,jj,mkji,nm,m,in;
	double xr,yr,zr,area,rr,rrx2,rrsq,rsec2r,rsecr,rsecn,b,alpha,
	       beta,tf,arcsum,t,tt,parea,ti,rsec2n;
	double *rad, *radsq;
	int    *cube;
	double zgrid, zres;
	double arci[maxintersect], arcf[maxintersect];
	double dx[maxintersect],dy[maxintersect],d[maxintersect];
	double dsq[maxintersect];
	double xmin, ymin, zmin, xmax, ymax, zmax, rmax;
	int nzp, karc;
	int inov[maxintersect], tag[maxintersect], **iatm,itab[maxcubes];

/* various initializations */

	iatm = (int **) ckalloc(maxatomscube * sizeof(int *));
	for(i = 0; i < maxatomscube; i ++)
		iatm[i] = (int *) ckalloc(maxcubes * sizeof(int));
	rad = (double *) ckalloc(natms * sizeof(double));
	radsq = (double *) ckalloc(natms * sizeof(double));
	cube = (int *) ckalloc(natms * sizeof(int));

/* assign defaults, if not set */

	if (errorpar == 0.0) errorpar = 0.05;

/* number of integration layers along z-axis: */

	nzp = (int) (1./errorpar+1.0);

	xmin= 9999.9;
	ymin= 9999.9;
	zmin= 9999.9;
	xmax=-9999.9;
	ymax=-9999.9;
	zmax=-9999.9;

	rmax=0.0;

	karc=maxintersect;

/* find bounding box for atoms */

	for(i = 0; i < natms; i ++)	{
		rad[i]=radx[i]+water_radii;
		radsq[i]=rad[i]*rad[i];

		if (rad[i] > rmax) rmax=rad[i];

		if (x[i] < -1000.0 || x[i] > 1000.0) continue;
		
		if (xmin > x[i]) xmin=x[i];
		if (ymin > y[i]) ymin=y[i];
		if (zmin > z[i]) zmin=z[i];
		if (xmax < x[i]) xmax=x[i];
		if (ymax < y[i]) ymax=y[i];
		if (zmax < z[i]) zmax=z[i];
	}

	rmax=rmax*2.0;

/* cubicals containing the atoms are setup. the dimension of an edge equals
 the radius of the largest atom sphere the cubes have a single index */

	idim=(int) ((xmax-xmin)/rmax+1.);
	if (idim < 3) idim=3;
	jidim=(int) ((ymax-ymin)/rmax+1.);
	if (jidim < 3) jidim=3;
	jidim=idim*jidim;
	kjidim=(int) ((zmax-zmin)/rmax+1.);
	if (kjidim < 3) kjidim=3;
	kjidim=jidim*kjidim;

	if (kjidim > maxcubes) {
		printf("idim = %d, jidim = %d, kjidim = %d\n", idim, jidim, kjidim);
		printf("kjidim = %d, maxcubes = %d\n", kjidim, maxcubes);
/*****************************************************************************
		alert("area", 4, "too many cubes - increase maxats", 32, 1);
*****************************************************************************/
	}

/* prepare upto ncube cubes each containing upto maxatomscube atoms.
 the cube index is kji. the number of atoms in each cube is in itab;
 the cube index for each atom is in cube; */

	for(l = 0; l < maxcubes; l ++) itab[l]=0;

	for(l = 0; l < natms; l ++)	{
        	i= (int ) ((x[l]-xmin)/rmax+1.);
       	 	j= (int ) ((y[l]-ymin)/rmax);
       	 	k= (int ) ((z[l]-zmin)/rmax);
        	kji=k*jidim+j*idim+i;
        	n=itab[kji - 1]+1;

		if(n > maxatomscube)	{
			fprintf(stdout, "exceeds maxcube %d\n", n);
/*********************************************************************************************
			alert("area", 4, "too many atoms per cube - increase maxatomscube",
				47, 1);
*********************************************************************************************/
		}

		itab[kji - 1]=n;
		iatm[n - 1][kji - 1]=l+1;
		cube[l]=kji;
	}

/* process each atom */

	for(ir = 0; ir < natms; ir ++)	{
		if (x[ir] < -1000.0 || x[ir] > 1000.0) continue;
        	kji=cube[ir];
		io=0;
		area=0.0;
        	xr=x[ir];
        	yr=y[ir];
        	zr=z[ir];
/*
        	rr=rad(ir);
*/
		rr = probe[ir];
        	rrx2=rr*2.;
		rrsq = rr * rr;

/* find the mkji cubes neighboring the kji cube */

		for(kk = 0; kk < 3; kk ++)	{
			k = kk - 1;
			for(jj = 0; jj < 3; jj ++)	{
				j = jj - 1;
				for(i = 0; i < 3; i ++)	{
					mkji = kji + k * jidim + j * idim + i - 1;
					if(mkji < 1)	continue;
					if(mkji > kjidim) break;
					nm = itab[mkji - 1];
					if(nm < 1)	continue;

/* record the atoms in inov that neighbor atom ir */

					for(m = 0; m < nm; m ++)	{
						in = iatm[m][mkji - 1];
						if(in == ir + 1)	continue;
						io ++;
						if(io > maxintersect)	{
/*******************************************************************************************
						  alert("area",4,"too many intersections",
							22, 1);
*******************************************************************************************/
						}
						dx[io - 1]=xr-x[in - 1];
						dy[io - 1]=yr-y[in - 1];
						dsq[io - 1]=dx[io-1]*dx[io-1]+dy[io-1]*dy[io-1];
                				d[io - 1]=sqrt(dsq[io - 1]);
						inov[io - 1]=in;
					}
				}
			}
		}

		if(io == 0)	{
			area=pix2*rrx2;
			goto label1;
		}

/* z resolution determined */

		zres=rrx2/nzp;
		zgrid=z[ir]-rr-zres/2.;

/* section atom spheres perpendicular to the z axis */

		for(i = 0; i < nzp; i ++)	{
        		zgrid=zgrid+zres;

/* find the radius of the circle of intersection of the ir sphere
 on the current z-plane */

			rsec2r=rrsq-(zgrid-zr)*(zgrid-zr);
			rsecr=sqrt(rsec2r);
			for(k = 0; k < karc; k ++)
				arci[k] = 0.0;
			karc = 0;
			for(j = 0; j < io; j ++)	{
				in=inov[j];

/* find radius of circle locus */

				rsec2n=radsq[in-1]-(zgrid-z[in-1])*(zgrid-z[in-1]);
				if (rsec2n <= 0.0)	continue;
				rsecn=sqrt(rsec2n);

/* find intersections of n.circles with ir circles in section */

				if (d[j] >= rsecr+rsecn) continue;

/* do the circles intersect, or is one circle completely inside the other? */

				b=rsecr-rsecn;
				if (d[j] <= fabs(b))	{
					if (b <= 0.0) goto label2;
					continue;
				}

/* if the circles intersect, find the points of intersection */

				karc=karc+1;

				if(karc >= maxintersect)	{
/*************************************************************************************
                                     alert("area",4,"too many intersections", 22, 1);
*************************************************************************************/
				}

/* initial and final arc endpoints are found for the ir circle intersected
 by a neighboring circle contained in the same plane. the initial endpoint
 of the enclosed arc is stored in arci, and the final arc in arcf
 law of cosines */

				alpha=acos((dsq[j]+rsec2r-rsec2n)/(2.*d[j]*rsecr));

/* alpha is the angle between a line containing a point of intersection and
 the reference circle center and the line containing both circle centers */

				beta=atan2(dy[j],dx[j])+pi;

/* beta is the angle between the line containing both circle centers and the
 x-axis */

				ti=  (beta-alpha);
				tf= beta+alpha;
				if (ti < 0.0)ti=ti+pix2;
				if (tf > pix2)tf=tf-pix2;
				arci[karc - 1]=ti;
				if (tf < ti)	{

/* if the arc crosses zero, then it is broken into two segments.
 the first ends at pix2 and the second begins at zero */

					arcf[karc - 1]=pix2;
					karc=karc+1;
				}
				arcf[karc - 1]=tf;
			}

/* find the accssible contact surface area for the sphere ir on
 this section */

			if (karc == 0)	{
				arcsum=pix2;
				goto label3;
			}

/* the arc endpoints are sorted on the value of the initial arc endpoint */

			sortag(arci,karc,tag);

/* calculate the accssible area */

			arcsum=arci[0];
			t=arcf[tag[0] - 1];

			if(karc != 1)	{
				for(k = 1; k < karc; k ++)	{
					if (t < arci[k]) arcsum=arcsum+arci[k]-t;
					tt=arcf[tag[k] - 1];
					if (tt > t) t=tt;
				}
			}
			arcsum=arcsum+pix2-t;

/* the area/radius is equal to the accssible arc length x the section thickness. */

label3:    		parea=arcsum*zres;

/* add the accssible area for this atom in this section to the area for this
 atom for all the section encountered thus far */

			area=area+parea;
label2:			continue;
		}
      
/*	The solvent shell */

label1:		accss[ir]=area*rr;
	}

	free((void *) rad);
	free((void *) radsq);
	free((void *) cube);
	for(i = 0; i < maxatomscube; i ++)
		free((void *) iatm[i]);
	free((void **) iatm);
}
