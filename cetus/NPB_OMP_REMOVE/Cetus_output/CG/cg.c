#include <stdlib.h>
/*
Copyright (C) 1991-2012 Free Software Foundation, Inc.
   This file is part of the GNU C Library.

   The GNU C Library is free software; you can redistribute it andor
   modify it under the terms of the GNU Lesser General Public
   License as published by the Free Software Foundation; either
   version 2.1 of the License, or (at your option) any later version.

   The GNU C Library is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public
   License along with the GNU C Library; if not, see
   <http:www.gnu.org/licenses/>. 
*/
/*
This header is separate from features.h so that the compiler can
   include it implicitly at the start of every compilation.  It must
   not itself include <features.h> or any other header that includes
   <features.h> because the implicit include comes before any feature
   test macros that may be defined in a source file before it first
   explicitly includes a system header.  GCC knows the name of this
   header in order to preinclude it. 
*/
/* We do support the IEC 559 math functionality, real and complex.  */
/*
wchar_t uses ISOIEC 10646 (2nd ed., published 2011-03-15) /
   Unicode 6.0. 
*/
/* We do not support C11 <threads.h>.  */
/* ------------------------------------------------------------------------- */
/*                                                                          */
/*  This benchmark is an OpenMP C version of the NPB CG code. This OpenMP   */
/*  C version is developed by the Center for Manycore Programming at Seoul  */
/*  National University and derived from the OpenMP Fortran versions in     */
/*  "NPB3.3-OMP" developed by NAS.                                          */
/*                                                                          */
/*  Permission to use, copy, distribute and modify this software for any    */
/*  purpose with or without fee is hereby granted. This software is         */
/*  provided "as is" without express or implied warranty.                   */
/*                                                                          */
/*  Information on NPB 3.3, including the technical report, the original    */
/*  specifications, source code, results and information on how to submit   */
/*  new results, is available at:                                           */
/*                                                                          */
/*           http:www.nas.nasa.govSoftware/NPB/                          */
/*                                                                          */
/*  Send comments or suggestions for this OpenMP C version to               */
/*  cmp@aces.snu.ac.kr                                                      */
/*                                                                          */
/*          Center for Manycore Programming                                 */
/*          School of Computer Science and Engineering                      */
/*          Seoul National University                                       */
/*          Seoul 151-744, Korea                                            */
/*                                                                          */
/*          E-mail:  cmp@aces.snu.ac.kr                                     */
/*                                                                          */
/* ------------------------------------------------------------------------- */
/* ------------------------------------------------------------------------- */
/* Authors: Sangmin Seo, Jungwon Kim, Jun Lee, Jeongho Nah, Gangwon Jo,     */
/*          and Jaejin Lee                                                  */
/* ------------------------------------------------------------------------- */
/* --------------------------------------------------------------------- */
/* program cg */
/* --------------------------------------------------------------------- */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "globals.h"
#include "randdp.h"
#include "timers.h"
#include "print_results.h"
/* --------------------------------------------------------------------- */
/* common main_int_mem */
static int colidx[((150000*(15+1))*(15+1))];
static int rowstr[(150000+1)];
static int iv[((((150000*(15+1))*(15+1))+1)+150000)];
static int arow[(150000+1)];
static int acol[(150000*(15+1))];
/* common main_flt_mem */
static double v[((150000*(15+1))*(15+1))];
static double aelt[(150000*(15+1))];
static double a[((150000*(15+1))*(15+1))];
static double x[(150000+2)];
static double z[(150000+2)];
static double p[(150000+2)];
static double q[(150000+2)];
static double r[(150000+2)];
/* commontinof */
static int myid, num_threads, ilow, ihigh;
static int last_n[(1024+1)];
/* common partit_size */
static int naa;
static int nzz;
static int firstrow;
static int lastrow;
static int firstcol;
static int lastcol;
/* commonurando */
static double amult;
static double tran;
/* commontimers */
static logical timeron;
/* --------------------------------------------------------------------- */
/* --------------------------------------------------------------------- */
static void conj_grad(int colidx[], int rowstr[], double x[], double z[], double a[], double p[], double q[], double r[], double * rnorm);
static void makea(int n, int nz, double a[], int colidx[], int rowstr[], int firstrow, int lastrow, int firstcol, int lastcol, int arow[], int acol[][(15+1)], double aelt[][(15+1)], double v[], int iv[]);
static void sparse(double a[], int colidx[], int rowstr[], int n, int nz, int nozer, int arow[], int acol[][(15+1)], double aelt[][(15+1)], int firstrow, int lastrow, int last_n[], double v[], int iv[], int nzloc[], double rcond, double shift);
static void sprnvc(int n, int nz, int nn1, double v[], int iv[]);
static int icnvrt(double x, int ipwr2);
static void vecset(int n, double v[], int iv[], int * nzv, int i, double val);
/* --------------------------------------------------------------------- */
int main(int argc, char * argv[])
{
	int i, j, k, it;
	double zeta;
	double rnorm;
	double norm_temp1, norm_temp2;
	double t, mflops, tmax;
	char Class;
	logical verified;
	double zeta_verify_value, epsilon, err;
	char * t_names[3];
	FILE * fp;
	#pragma cetus private(i) 
	#pragma loop name main#0 
	for (i=0; i<3; i ++ )
	{
		timer_clear(i);
	}
	if ((fp=fopen("timer.flag", "r"))!=((void * )0))
	{
		timeron=true;
		t_names[0]="init";
		t_names[1]="benchmk";
		t_names[2]="conjgd";
		fclose(fp);
	}
	else
	{
		timeron=false;
	}
	timer_start(0);
	firstrow=0;
	lastrow=(150000-1);
	firstcol=0;
	lastcol=(150000-1);
	Class='C';
	zeta_verify_value=28.973605592845;
	printf("\n\n NAS Parallel Benchmarks (NPB3.3-OMP-C) - CG Benchmark\n\n");
	printf(" Size: %11d\n", 150000);
	printf(" Iterations:                  %5d\n", 75);
	printf(" Number of available threads: %5d\n", omp_get_max_threads());
	printf("\n");
	naa=150000;
	nzz=((150000*(15+1))*(15+1));
	/* --------------------------------------------------------------------- */
	/* Inialize random number generator */
	/* --------------------------------------------------------------------- */
	tran=3.14159265E8;
	amult=1.220703125E9;
	zeta=randlc( & tran, amult);
	/* --------------------------------------------------------------------- */
	/*   */
	/* --------------------------------------------------------------------- */
	makea(naa, nzz, a, colidx, rowstr, firstrow, lastrow, firstcol, lastcol, arow, (int (* )[(15+1)])((void * )acol), (double (* )[(15+1)])((void * )aelt), v, iv);
	/* --------------------------------------------------------------------- */
	/* Note: as a result of the above call to makea: */
	/*    values of j used in indexing rowstr go from 0 --> lastrow-firstrow */
	/*    values of colidx which are col indexes go from firstcol --> lastcol */
	/*    So: */
	/*    Shift the col index vals from actual (firstcol --> lastcol )  */
	/*    to local, i.e., (0 --> lastcol-firstcol) */
	/* --------------------------------------------------------------------- */
	#pragma cetus parallel 
	#pragma cetus private(j, k) 
	#pragma omp parallel private(j, k)
	{
		int * reduce = (int * )malloc(((150000*(15+1))*(15+1))*sizeof (int));
		int reduce_span_0;
		for (reduce_span_0=0; reduce_span_0<((150000*(15+1))*(15+1)); reduce_span_0 ++ )
		{
			reduce[reduce_span_0]=0;
		}
		#pragma loop name main#1 
		#pragma cetus for  
		#pragma omp for
		for (j=0; j<((lastrow-firstrow)+1); j ++ )
		{
			#pragma cetus private(k) 
			#pragma loop name main#1#0 
			for (k=rowstr[j]; k<rowstr[j+1]; k ++ )
			{
				reduce[k]=(reduce[k]-firstcol);
			}
		}
		#pragma cetus critical  
		#pragma omp critical
		{
			for (reduce_span_0=0; reduce_span_0<((150000*(15+1))*(15+1)); reduce_span_0 ++ )
			{
				colidx[reduce_span_0]+=reduce[reduce_span_0];
			}
		}
	}
	/* --------------------------------------------------------------------- */
	/* set starting vector to (1, 1, .... 1) */
	/* --------------------------------------------------------------------- */
	#pragma cetus private(i) 
	#pragma loop name main#2 
	#pragma cetus parallel 
	#pragma omp parallel for private(i)
	for (i=0; i<(150000+1); i ++ )
	{
		x[i]=1.0;
	}
	#pragma cetus private(j) 
	#pragma loop name main#3 
	#pragma cetus parallel 
	#pragma omp parallel for if((10000<((7L+(-6L*firstcol))+(6L*lastcol)))) private(j)
	for (j=0; j<((lastcol-firstcol)+1); j ++ )
	{
		q[j]=0.0;
		z[j]=0.0;
		r[j]=0.0;
		p[j]=0.0;
	}
	zeta=0.0;
	/* --------------------------------------------------------------------- */
	/* ----> */
	/* Do one iteration untimed to init all code and data page tables */
	/* ---->                    (then reinit, start timing, to niter its) */
	/* --------------------------------------------------------------------- */
	#pragma cetus private(it, j, norm_temp1, norm_temp2) 
	#pragma loop name main#4 
	for (it=1; it<=1; it ++ )
	{
		/* --------------------------------------------------------------------- */
		/* The call to the conjugate gradient routine: */
		/* --------------------------------------------------------------------- */
		conj_grad(colidx, rowstr, x, z, a, p, q, r,  & rnorm);
		/* --------------------------------------------------------------------- */
		/* zeta = shift + 1(x.z) */
		/* So, first: (x.z) */
		/* Also, find norm of z */
		/* So, first: (z.z) */
		/* --------------------------------------------------------------------- */
		norm_temp1=0.0;
		norm_temp2=0.0;
		#pragma cetus private(j) 
		#pragma loop name main#4#0 
		#pragma cetus reduction(+: norm_temp1, norm_temp2) 
		#pragma cetus parallel 
		#pragma omp parallel for if((10000<((5L+(-4L*firstcol))+(4L*lastcol)))) private(j) reduction(+: norm_temp1, norm_temp2)
		for (j=0; j<((lastcol-firstcol)+1); j ++ )
		{
			norm_temp1=(norm_temp1+(x[j]*z[j]));
			norm_temp2=(norm_temp2+(z[j]*z[j]));
		}
		norm_temp2=(1.0/sqrt(norm_temp2));
		/* --------------------------------------------------------------------- */
		/* Normalize z to obtain x */
		/* --------------------------------------------------------------------- */
		#pragma cetus private(j) 
		#pragma loop name main#4#1 
		#pragma cetus parallel 
		#pragma omp parallel for if((10000<((4L+(-3L*firstcol))+(3L*lastcol)))) private(j)
		for (j=0; j<((lastcol-firstcol)+1); j ++ )
		{
			x[j]=(norm_temp2*z[j]);
		}
	}
	/* end of do one iteration untimed */
	/* --------------------------------------------------------------------- */
	/* set starting vector to (1, 1, .... 1) */
	/* --------------------------------------------------------------------- */
	#pragma cetus private(i) 
	#pragma loop name main#5 
	#pragma cetus parallel 
	#pragma omp parallel for private(i)
	for (i=0; i<(150000+1); i ++ )
	{
		x[i]=1.0;
	}
	zeta=0.0;
	timer_stop(0);
	printf(" Initialization time = %15.3f seconds\n", timer_read(0));
	timer_start(1);
	/* --------------------------------------------------------------------- */
	/* ----> */
	/* Main Iteration for inverse power method */
	/* ----> */
	/* --------------------------------------------------------------------- */
	#pragma cetus private(it, j, norm_temp1, norm_temp2) 
	#pragma cetus lastprivate(zeta) 
	#pragma loop name main#6 
	for (it=1; it<=75; it ++ )
	{
		/* --------------------------------------------------------------------- */
		/* The call to the conjugate gradient routine: */
		/* --------------------------------------------------------------------- */
		if (timeron)
		{
			timer_start(2);
		}
		conj_grad(colidx, rowstr, x, z, a, p, q, r,  & rnorm);
		if (timeron)
		{
			timer_stop(2);
		}
		/* --------------------------------------------------------------------- */
		/* zeta = shift + 1(x.z) */
		/* So, first: (x.z) */
		/* Also, find norm of z */
		/* So, first: (z.z) */
		/* --------------------------------------------------------------------- */
		norm_temp1=0.0;
		norm_temp2=0.0;
		#pragma cetus private(j) 
		#pragma loop name main#6#0 
		#pragma cetus reduction(+: norm_temp1, norm_temp2) 
		#pragma cetus parallel 
		#pragma omp parallel for if((10000<((5L+(-4L*firstcol))+(4L*lastcol)))) private(j) reduction(+: norm_temp1, norm_temp2)
		for (j=0; j<((lastcol-firstcol)+1); j ++ )
		{
			norm_temp1=(norm_temp1+(x[j]*z[j]));
			norm_temp2=(norm_temp2+(z[j]*z[j]));
		}
		norm_temp2=(1.0/sqrt(norm_temp2));
		zeta=(110.0+(1.0/norm_temp1));
		if (it==1)
		{
			printf("\n   iteration           ||r||                 zeta\n");
		}
		printf("    %5d       %20.14E%20.13f\n", it, rnorm, zeta);
		/* --------------------------------------------------------------------- */
		/* Normalize z to obtain x */
		/* --------------------------------------------------------------------- */
		#pragma cetus private(j) 
		#pragma loop name main#6#1 
		#pragma cetus parallel 
		#pragma omp parallel for if((10000<((4L+(-3L*firstcol))+(3L*lastcol)))) private(j)
		for (j=0; j<((lastcol-firstcol)+1); j ++ )
		{
			x[j]=(norm_temp2*z[j]);
		}
	}
	/* end of main iter inv pow meth */
	timer_stop(1);
	/* --------------------------------------------------------------------- */
	/* End of timed section */
	/* --------------------------------------------------------------------- */
	t=timer_read(1);
	printf(" Benchmark completed\n");
	epsilon=1.0E-10;
	if (Class!='U')
	{
		err=(fabs(zeta-zeta_verify_value)/zeta_verify_value);
		if (err<=epsilon)
		{
			verified=true;
			printf(" VERIFICATION SUCCESSFUL\n");
			printf(" Zeta is    %20.13E\n", zeta);
			printf(" Error is   %20.13E\n", err);
		}
		else
		{
			verified=false;
			printf(" VERIFICATION FAILED\n");
			printf(" Zeta                %20.13E\n", zeta);
			printf(" The correct zeta is %20.13E\n", zeta_verify_value);
		}
	}
	else
	{
		verified=false;
		printf(" Problem size unknown\n");
		printf(" NO VERIFICATION PERFORMED\n");
	}
	if (t!=0.0)
	{
		mflops=(((((double)((2*75)*150000))*(((3.0+((double)(15*(15+1))))+(25.0*(5.0+((double)(15*(15+1))))))+3.0))/t)/1000000.0);
	}
	else
	{
		mflops=0.0;
	}
	print_results("CG", Class, 150000, 0, 0, 75, t, mflops, "          floating point", verified, "3.3.1", "06 Nov 2019", "icc", "$(CC)", "-lm", "-I../common", "-g -Wall -O3 -qopenmp -mcmodel=medium", "-O3 -qopenmp -mcmodel=medium", "randdp");
	/* --------------------------------------------------------------------- */
	/* More timers */
	/* --------------------------------------------------------------------- */
	if (timeron)
	{
		tmax=timer_read(1);
		if (tmax==0.0)
		{
			tmax=1.0;
		}
		printf("  SECTION   Time (secs)\n");
		#pragma cetus private(i, t) 
		#pragma loop name main#7 
		for (i=0; i<3; i ++ )
		{
			t=timer_read(i);
			if (i==0)
			{
				printf("  %8s:%9.3f\n", t_names[i], t);
			}
			else
			{
				printf("  %8s:%9.3f  (%6.2f%%)\n", t_names[i], t, (t*100.0)/tmax);
				if (i==2)
				{
					t=(tmax-t);
					printf("    --> %8s:%9.3f  (%6.2f%%)\n", "rest", t, (t*100.0)/tmax);
				}
			}
		}
	}
	return 0;
}

/* --------------------------------------------------------------------- */
/* Floaging point arrays here are named as in NPB1 spec discussion of  */
/* CG algorithm */
/* --------------------------------------------------------------------- */
static void conj_grad(int colidx[], int rowstr[], double x[], double z[], double a[], double p[], double q[], double r[], double * rnorm)
{
	int j, k;
	int cgit, cgitmax = 25;
	double d, sum, rho, rho0, alpha, beta, suml;
	rho=0.0;
	sum=0.0;
	/* --------------------------------------------------------------------- */
	/* Initialize the CG algorithm: */
	/* --------------------------------------------------------------------- */
	#pragma cetus private(j) 
	#pragma loop name conj_grad#0 
	#pragma cetus parallel 
	#pragma omp parallel for if((10000<(7L+(6L*naa)))) private(j)
	for (j=0; j<(naa+1); j ++ )
	{
		q[j]=0.0;
		z[j]=0.0;
		r[j]=x[j];
		p[j]=r[j];
	}
	/* --------------------------------------------------------------------- */
	/* rho = r.r */
	/* Now, obtain the norm of r: First, sum squares of r elements locally... */
	/* --------------------------------------------------------------------- */
	#pragma cetus private(j) 
	#pragma loop name conj_grad#1 
	#pragma cetus reduction(+: rho) 
	#pragma cetus parallel 
	#pragma omp parallel for if((10000<((4L+(-3L*firstcol))+(3L*lastcol)))) private(j) reduction(+: rho)
	for (j=0; j<((lastcol-firstcol)+1); j ++ )
	{
		rho=(rho+(r[j]*r[j]));
	}
	/* --------------------------------------------------------------------- */
	/* ----> */
	/* The conj grad iteration loop */
	/* ----> */
	/* --------------------------------------------------------------------- */
	#pragma cetus private(alpha, beta, cgit, d, j, k, rho0, suml) 
	#pragma loop name conj_grad#2 
	/* #pragma cetus reduction(+: z[j])  */
	for (cgit=1; cgit<=cgitmax; cgit ++ )
	{
		/* --------------------------------------------------------------------- */
		/* Save a temporary of rho and initialize reduction variables */
		/* --------------------------------------------------------------------- */
		rho0=rho;
		d=0.0;
		rho=0.0;
		/* --------------------------------------------------------------------- */
		/* q = A.p */
		/* The partition submatrix-vector multiply: use workspace w */
		/* --------------------------------------------------------------------- */
		/*  */
		/* NOTE: this version of the multiply is actually (slightly: maybe %5)  */
		/*       faster on the sp2 on 16 nodes than is the unrolled-by-2 version  */
		/*       below.   On the Cray t3d, the reverse is true, i.e., the  */
		/*       unrolled-by-two version is some 10% faster.   */
		/*       The unrolled-by-8 version below is significantly faster */
		/*       on the Cray t3d - overall speed of code is 1.5 times faster. */
		#pragma cetus private(j, k, suml) 
		#pragma loop name conj_grad#2#0 
		#pragma cetus parallel 
		#pragma omp parallel for private(j, k, suml)
		for (j=0; j<((lastrow-firstrow)+1); j ++ )
		{
			suml=0.0;
			#pragma cetus private(k) 
			#pragma loop name conj_grad#2#0#0 
			/* #pragma cetus reduction(+: suml)  */
			for (k=rowstr[j]; k<rowstr[j+1]; k ++ )
			{
				suml=(suml+(a[k]*p[colidx[k]]));
			}
			q[j]=suml;
		}
		/*
		
		    for (j = 0; j < lastrow - firstrow + 1; j++) {
			      int i = rowstr[j];
			      int iresidue = (rowstr[j+1] - i) % 2;
			      double sum1 = 0.0;
			      double sum2 = 0.0;
			      if (iresidue == 1)
			        sum1 = sum1 + a[i]p[colidx[i]];
			      for (k = i + iresidue; k <= rowstr[j+1] - 2; k += 2) {
				        sum1 = sum1 + a[k]  *p[colidx[k]];
				        sum2 = sum2 + a[k+1]*p[colidx[k+1]];
			      }
			      q[j] = sum1 + sum2;
		    }
		   
		*/
		/*
		
		    for (j = 0; j < lastrow - firstrow + 1; j++) {
			      int i = rowstr[j]; 
			      int iresidue = (rowstr[j+1] - i) % 8;
			      suml = 0.0;
			      for (k = i; k <= i + iresidue - 1; k++) {
				        suml = suml + a[k]p[colidx[k]];
			      }
			      for (k = i + iresidue; k <= rowstr[j+1] - 8; k += 8) {
				        suml = suml + a[k  ]*p[colidx[k  ]]
				                  + a[k+1]*p[colidx[k+1]]
				                  + a[k+2]*p[colidx[k+2]]
				                  + a[k+3]*p[colidx[k+3]]
				                  + a[k+4]*p[colidx[k+4]]
				                  + a[k+5]*p[colidx[k+5]]
				                  + a[k+6]*p[colidx[k+6]]
				                  + a[k+7]*p[colidx[k+7]];
			      }
			      q[j] = suml;
		    }
		   
		*/
		/* --------------------------------------------------------------------- */
		/* Obtain p.q */
		/* --------------------------------------------------------------------- */
		#pragma cetus private(j) 
		#pragma loop name conj_grad#2#1 
		#pragma cetus reduction(+: d) 
		#pragma cetus parallel 
		#pragma omp parallel for if((10000<((4L+(-3L*firstcol))+(3L*lastcol)))) private(j) reduction(+: d)
		for (j=0; j<((lastcol-firstcol)+1); j ++ )
		{
			d=(d+(p[j]*q[j]));
		}
		/* --------------------------------------------------------------------- */
		/* Obtain alpha = rho (p.q) */
		/* --------------------------------------------------------------------- */
		alpha=(rho0/d);
		/* --------------------------------------------------------------------- */
		/* Obtain z = z + alphap */
		/* and    r = r - alphaq */
		/* --------------------------------------------------------------------- */
		#pragma cetus private(j) 
		#pragma loop name conj_grad#2#2 
		#pragma cetus reduction(+: rho) 
		#pragma cetus parallel 
		#pragma omp parallel for if((10000<((6L+(-5L*firstcol))+(5L*lastcol)))) private(j) reduction(+: rho)
		for (j=0; j<((lastcol-firstcol)+1); j ++ )
		{
			z[j]=(z[j]+(alpha*p[j]));
			r[j]=(r[j]-(alpha*q[j]));
			/* --------------------------------------------------------------------- */
			/* rho = r.r */
			/* Now, obtain the norm of r: First, sum squares of r elements locally.. */
			/* --------------------------------------------------------------------- */
			rho=(rho+(r[j]*r[j]));
		}
		/* --------------------------------------------------------------------- */
		/* Obtain beta: */
		/* --------------------------------------------------------------------- */
		beta=(rho/rho0);
		/* --------------------------------------------------------------------- */
		/* p = r + betap */
		/* --------------------------------------------------------------------- */
		#pragma cetus private(j) 
		#pragma loop name conj_grad#2#3 
		#pragma cetus parallel 
		#pragma omp parallel for if((10000<((4L+(-3L*firstcol))+(3L*lastcol)))) private(j)
		for (j=0; j<((lastcol-firstcol)+1); j ++ )
		{
			p[j]=(r[j]+(beta*p[j]));
		}
	}
	/* end of do cgit=1,cgitmax */
	/* --------------------------------------------------------------------- */
	/* Compute residual norm explicitly:  ||r|| = ||x - A.z|| */
	/* First, form A.z */
	/* The partition submatrix-vector multiply */
	/* --------------------------------------------------------------------- */
	#pragma cetus private(j, k, suml) 
	#pragma loop name conj_grad#3 
	#pragma cetus parallel 
	#pragma omp parallel for private(j, k, suml)
	for (j=0; j<((lastrow-firstrow)+1); j ++ )
	{
		suml=0.0;
		#pragma cetus private(k) 
		#pragma loop name conj_grad#3#0 
		/* #pragma cetus reduction(+: suml)  */
		for (k=rowstr[j]; k<rowstr[j+1]; k ++ )
		{
			suml=(suml+(a[k]*z[colidx[k]]));
		}
		r[j]=suml;
	}
	/* --------------------------------------------------------------------- */
	/* At this point, r contains A.z */
	/* --------------------------------------------------------------------- */
	#pragma cetus private(j, suml) 
	#pragma loop name conj_grad#4 
	#pragma cetus reduction(+: sum) 
	#pragma cetus parallel 
	#pragma omp parallel for if((10000<((5L+(-4L*firstcol))+(4L*lastcol)))) private(j, suml) reduction(+: sum)
	for (j=0; j<((lastcol-firstcol)+1); j ++ )
	{
		suml=(x[j]-r[j]);
		sum=(sum+(suml*suml));
	}
	( * rnorm)=sqrt(sum);
}

/* --------------------------------------------------------------------- */
/* generate the test problem for benchmark 6 */
/* makea generates a sparse matrix with a */
/* prescribed sparsity distribution */
/*  */
/* parameter    type        usage */
/*  */
/* input */
/*  */
/* n            i           number of colsrows of matrix */
/* nz           i           nonzeros as declared array size */
/* rcond        r8         condition number */
/* shift        r8         main diagonal shift */
/*  */
/* output */
/*  */
/* a            r8         array for nonzeros */
/* colidx       i           col indices */
/* rowstr       i           row pointers */
/*  */
/* workspace */
/*  */
/* iv, arow, acol i */
/* aelt           r8 */
/* --------------------------------------------------------------------- */
static void makea(int n, int nz, double a[], int colidx[], int rowstr[], int firstrow, int lastrow, int firstcol, int lastcol, int arow[], int acol[][(15+1)], double aelt[][(15+1)], double v[], int iv[])
{
	int iouter, ivelt, nzv, nn1;
	int ivc[(15+1)];
	double vc[(15+1)];
	/* --------------------------------------------------------------------- */
	/* nonzer is approximately  (int(sqrt(nnzan))); */
	/* --------------------------------------------------------------------- */
	int work;
	/* --------------------------------------------------------------------- */
	/* nn1 is the smallest power of two not less than n */
	/* --------------------------------------------------------------------- */
	nn1=1;
	do
	{
		nn1=(2*nn1);
	}while(nn1<n);
	
	/* --------------------------------------------------------------------- */
	/* Generate nonzero positions and save for the use in sparse. */
	/* --------------------------------------------------------------------- */
	num_threads=omp_get_num_threads();
	myid=omp_get_thread_num();
	if (num_threads>1024)
	{
		if (myid==0)
		{
			printf(" Warning: num_threads%6d exceeded an internal limit%6d\n", num_threads, 1024);
		}
		num_threads=1024;
	}
	work=(((n+num_threads)-1)/num_threads);
	ilow=(work*myid);
	ihigh=(ilow+work);
	if (ihigh>n)
	{
		ihigh=n;
	}
	#pragma cetus private(iouter, ivelt) 
	#pragma loop name makea#0 
	for (iouter=0; iouter<ihigh; iouter ++ )
	{
		nzv=15;
		sprnvc(n, nzv, nn1, vc, ivc);
		if (iouter>=ilow)
		{
			vecset(n, vc, ivc,  & nzv, iouter+1, 0.5);
			arow[iouter]=nzv;
			#pragma cetus private(ivelt) 
			#pragma loop name makea#0#0 
			#pragma cetus parallel 
			#pragma omp parallel for if((10000<(1L+(4L*nzv)))) private(ivelt)
			for (ivelt=0; ivelt<nzv; ivelt ++ )
			{
				acol[iouter][ivelt]=(ivc[ivelt]-1);
				aelt[iouter][ivelt]=vc[ivelt];
			}
		}
	}
	/* --------------------------------------------------------------------- */
	/* ... make the sparse matrix from list of elements with duplicates */
	/*     (v and iv are used as  workspace) */
	/* --------------------------------------------------------------------- */
	sparse(a, colidx, rowstr, n, nz, 15, arow, acol, aelt, firstrow, lastrow, last_n, v,  & iv[0],  & iv[nz], 0.1, 110.0);
}

/* --------------------------------------------------------------------- */
/* rows range from firstrow to lastrow */
/* the rowstr pointers are defined for nrows = lastrow-firstrow+1 values */
/* --------------------------------------------------------------------- */
static void sparse(double a[], int colidx[], int rowstr[], int n, int nz, int nozer, int arow[], int acol[][(15+1)], double aelt[][(15+1)], int firstrow, int lastrow, int last_n[], double v[], int iv[], int nzloc[], double rcond, double shift)
{
	int nrows;
	/* --------------------------------------------------- */
	/* generate a sparse matrix from a list of */
	/* [col, row, element] tri */
	/* --------------------------------------------------- */
	int i, j, j1, j2, nza, k, kk, nzrow, jcol;
	double size, scale, ratio, va;
	logical cont40;
	/* --------------------------------------------------------------------- */
	/* how many rows of result */
	/* --------------------------------------------------------------------- */
	nrows=((lastrow-firstrow)+1);
	j1=(ilow+1);
	j2=(ihigh+1);
	/* --------------------------------------------------------------------- */
	/* ...count the number of triples in each row */
	/* --------------------------------------------------------------------- */
	#pragma cetus private(j) 
	#pragma loop name sparse#0 
	#pragma cetus parallel 
	#pragma omp parallel for if((10000<((1L+(3L*ihigh))+(-3L*ilow)))) private(j)
	for (j=j1; j<j2; j ++ )
	{
		rowstr[j]=0;
	}
	#pragma cetus private(i, j, nza) 
	#pragma loop name sparse#1 
	/* #pragma cetus reduction(+: rowstr[j])  */
	for (i=0; i<n; i ++ )
	{
		#pragma cetus private(j, nza) 
		#pragma loop name sparse#1#0 
		/* #pragma cetus reduction(+: rowstr[j])  */
		for (nza=0; nza<arow[i]; nza ++ )
		{
			j=acol[i][nza];
			if ((j>=ilow)&&(j<ihigh))
			{
				j=(j+1);
				rowstr[j]=(rowstr[j]+arow[i]);
			}
		}
	}
	if (myid==0)
	{
		rowstr[0]=0;
		j1=0;
	}
	#pragma cetus private(j) 
	#pragma loop name sparse#2 
	for (j=(j1+1); j<j2; j ++ )
	{
		rowstr[j]=(rowstr[j]+rowstr[j-1]);
	}
	if (myid<num_threads)
	{
		last_n[myid]=rowstr[j2-1];
	}
	nzrow=0;
	if (myid<num_threads)
	{
		#pragma cetus private(i) 
		#pragma loop name sparse#3 
		#pragma cetus reduction(+: nzrow) 
		#pragma cetus parallel 
		#pragma omp parallel for if((10000<(1L+(3L*myid)))) private(i) reduction(+: nzrow)
		for (i=0; i<myid; i ++ )
		{
			nzrow=(nzrow+last_n[i]);
		}
	}
	if (nzrow>0)
	{
		#pragma cetus private(j) 
		#pragma loop name sparse#4 
		#pragma cetus parallel 
		#pragma omp parallel for if((10000<((4L+(3L*ihigh))+(-3L*j1)))) private(j)
		for (j=j1; j<j2; j ++ )
		{
			rowstr[j]=(rowstr[j]+nzrow);
		}
	}
	nza=(rowstr[nrows]-1);
	/* --------------------------------------------------------------------- */
	/* ... rowstr(j) now is the location of the first nonzero */
	/*     of row j of a */
	/* --------------------------------------------------------------------- */
	if (nza>nz)
	{
		printf("Space for matrix elements exceeded in sparse\n");
		printf("nza, nzmax = %d, %d\n", nza, nz);
		exit(1);
	}
	/* --------------------------------------------------------------------- */
	/* ... preload data pages */
	/* --------------------------------------------------------------------- */
	#pragma cetus private(j, k) 
	#pragma loop name sparse#5 
	for (j=ilow; j<ihigh; j ++ )
	{
		#pragma cetus private(k) 
		#pragma loop name sparse#5#0 
		#pragma cetus parallel 
		#pragma omp parallel for if((10000<((1L+(4L*rowstr[(1L+j)]))+(-4L*rowstr[j])))) private(k)
		for (k=rowstr[j]; k<rowstr[j+1]; k ++ )
		{
			v[k]=0.0;
			iv[k]=( - 1);
		}
		nzloc[j]=0;
	}
	/* --------------------------------------------------------------------- */
	/* ... generate actual values by summing duplicates */
	/* --------------------------------------------------------------------- */
	size=1.0;
	ratio=pow(rcond, 1.0/((double)n));
	#pragma cetus private(i, j, jcol, k, kk, nza, nzrow, scale, va) 
	#pragma loop name sparse#6 
	/* #pragma cetus reduction(+: nzloc[j])  */
	for (i=0; i<n; i ++ )
	{
		#pragma cetus private(j, jcol, k, kk, nza, nzrow, scale, va) 
		#pragma loop name sparse#6#0 
		/* #pragma cetus reduction(+: nzloc[j])  */
		for (nza=0; nza<arow[i]; nza ++ )
		{
			j=acol[i][nza];
			if ((j<ilow)||(j>=ihigh))
			{
				continue;
			}
			scale=(size*aelt[i][nza]);
			#pragma cetus private(jcol, k, kk, nzrow, va) 
			#pragma loop name sparse#6#0#0 
			/* #pragma cetus reduction(+: nzloc[j])  */
			for (nzrow=0; nzrow<arow[i]; nzrow ++ )
			{
				jcol=acol[i][nzrow];
				va=(aelt[i][nzrow]*scale);
				/* -------------------------------------------------------------------- */
				/* ... add the identity rcond to the generated matrix to bound */
				/*     the smallest eigenvalue from below by rcond */
				/* -------------------------------------------------------------------- */
				if ((jcol==j)&&(j==i))
				{
					va=((va+rcond)-shift);
				}
				cont40=false;
				#pragma cetus private(kk) 
				#pragma loop name sparse#6#0#0#0 
				/* #pragma cetus reduction(+: nzloc[j])  */
				for (k=rowstr[j]; k<rowstr[j+1]; k ++ )
				{
					if (iv[k]>jcol)
					{
						/* ---------------------------------------------------------------- */
						/* ... insert colidx here orderly */
						/* ---------------------------------------------------------------- */
						#pragma cetus private(kk) 
						#pragma loop name sparse#6#0#0#0#0 
						for (kk=(rowstr[j+1]-2); kk>=k; kk -- )
						{
							if (iv[kk]>( - 1))
							{
								v[kk+1]=v[kk];
								iv[kk+1]=iv[kk];
							}
						}
						iv[k]=jcol;
						v[k]=0.0;
						cont40=true;
						break;
					}
					else
					{
						if (iv[k]==( - 1))
						{
							iv[k]=jcol;
							cont40=true;
							break;
						}
						else
						{
							if (iv[k]==jcol)
							{
								/* -------------------------------------------------------------- */
								/* ... mark the duplicated entry */
								/* -------------------------------------------------------------- */
								nzloc[j]=(nzloc[j]+1);
								cont40=true;
								break;
							}
						}
					}
				}
				if (cont40==false)
				{
					printf("internal error in sparse: i=%d\n", i);
					exit(1);
				}
				v[k]=(v[k]+va);
			}
		}
		size=(size*ratio);
	}
	/* --------------------------------------------------------------------- */
	/* ... remove empty entries and generate final results */
	/* --------------------------------------------------------------------- */
	#pragma cetus private(j) 
	#pragma loop name sparse#7 
	for (j=(ilow+1); j<ihigh; j ++ )
	{
		nzloc[j]=(nzloc[j]+nzloc[j-1]);
	}
	if (myid<num_threads)
	{
		last_n[myid]=nzloc[ihigh-1];
	}
	nzrow=0;
	if (myid<num_threads)
	{
		#pragma cetus private(i) 
		#pragma loop name sparse#8 
		#pragma cetus reduction(+: nzrow) 
		#pragma cetus parallel 
		#pragma omp parallel for if((10000<(1L+(3L*myid)))) private(i) reduction(+: nzrow)
		for (i=0; i<myid; i ++ )
		{
			nzrow=(nzrow+last_n[i]);
		}
	}
	if (nzrow>0)
	{
		#pragma cetus private(j) 
		#pragma loop name sparse#9 
		#pragma cetus parallel 
		#pragma omp parallel for if((10000<((1L+(3L*ihigh))+(-3L*ilow)))) private(j)
		for (j=ilow; j<ihigh; j ++ )
		{
			nzloc[j]=(nzloc[j]+nzrow);
		}
	}
	#pragma cetus private(j, j1, j2, k, nza) 
	#pragma loop name sparse#10 
	for (j=0; j<nrows; j ++ )
	{
		if (j>0)
		{
			j1=(rowstr[j]-nzloc[j-1]);
		}
		else
		{
			j1=0;
		}
		j2=(rowstr[j+1]-nzloc[j]);
		nza=rowstr[j];
		#pragma cetus private(k) 
		#pragma loop name sparse#10#0 
		#pragma cetus parallel 
		#pragma omp parallel for if((10000<(((1L+(-4L*j1))+(-4L*nzloc[j]))+(4L*rowstr[(1L+j)])))) private(k)
		for (k=j1; k<j2; k ++ )
		{
			a[k]=v[((-1*j1)+k)+rowstr[j]];
			colidx[k]=iv[((-1*j1)+k)+rowstr[j]];
		}
		if ((((-1+(-1*j1))+(-1*nzloc[j]))+rowstr[1+j])>=0)
		{
			nza+=(((-1*j1)+(-1*nzloc[j]))+rowstr[1+j]);
		}
	}
	#pragma cetus private(j) 
	#pragma loop name sparse#11 
	#pragma cetus parallel 
	#pragma omp parallel for if((10000<(1L+(3L*nrows)))) private(j)
	for (j=1; j<(nrows+1); j ++ )
	{
		rowstr[j]=(rowstr[j]-nzloc[j-1]);
	}
	nza=(rowstr[nrows]-1);
}

/* --------------------------------------------------------------------- */
/* generate a sparse n-vector (v, iv) */
/* having nzv nonzeros */
/*  */
/* mark(i) is set to 1 if position i is nonzero. */
/* mark is all zero on entry and is reset to all zero before exit */
/* this corrects a performance bug found by John G. Lewis, caused by */
/* reinitialization of mark on every one of the n calls to sprnvc */
/* --------------------------------------------------------------------- */
static void sprnvc(int n, int nz, int nn1, double v[], int iv[])
{
	int nzv, ii, i;
	double vecelt, vecloc;
	nzv=0;
	while (nzv<nz)
	{
		logical was_gen = false;
		vecelt=randlc( & tran, amult);
		/* --------------------------------------------------------------------- */
		/* generate an integer between 1 and n in a portable manner */
		/* --------------------------------------------------------------------- */
		vecloc=randlc( & tran, amult);
		i=(icnvrt(vecloc, nn1)+1);
		if (i>n)
		{
			continue;
		}
		/* --------------------------------------------------------------------- */
		/* was this integer generated already? */
		/* --------------------------------------------------------------------- */
		#pragma loop name sprnvc#0 
		for (ii=0; ii<nzv; ii ++ )
		{
			if (iv[ii]==i)
			{
				was_gen=true;
				break;
			}
		}
		if (was_gen)
		{
			continue;
		}
		v[nzv]=vecelt;
		iv[nzv]=i;
		nzv=(nzv+1);
	}
}

/* --------------------------------------------------------------------- */
/* scale a double precision number x in (0,1) by a power of 2 and chop it */
/* --------------------------------------------------------------------- */
static int icnvrt(double x, int ipwr2)
{
	return ((int)(ipwr2*x));
}

/* --------------------------------------------------------------------- */
/* set ith element of sparse vector (v, iv) with */
/* nzv nonzeros to val */
/* --------------------------------------------------------------------- */
static void vecset(int n, double v[], int iv[], int * nzv, int i, double val)
{
	int k;
	logical set;
	set=false;
	#pragma cetus private(k) 
	#pragma loop name vecset#0 
	for (k=0; k<( * nzv); k ++ )
	{
		if (iv[k]==i)
		{
			v[k]=val;
			set=true;
		}
	}
	if (set==false)
	{
		v[ * nzv]=val;
		iv[ * nzv]=i;
		( * nzv)=(( * nzv)+1);
	}
}
