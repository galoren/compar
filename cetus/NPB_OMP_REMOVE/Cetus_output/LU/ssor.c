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
/*  This benchmark is an OpenMP C version of the NPB LU code. This OpenMP   */
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
#include <stdio.h>
#include "applu.incl"
#include "timers.h"
/* --------------------------------------------------------------------- */
/* Thread synchronization for pipeline operation */
/* --------------------------------------------------------------------- */
/* commonthreadinfo1 */
int isync[(162+1)];
/* commonthreadinfo2 */
int mthreadnum, iam;
/* --------------------------------------------------------------------- */
/* to perform pseudo-time stepping SSOR iterations */
/* for five nonlinear pde's. */
/* --------------------------------------------------------------------- */
void ssor(int niter)
{
	/* --------------------------------------------------------------------- */
	/* local variables */
	/* --------------------------------------------------------------------- */
	int i, j, k, m, n;
	int istep;
	double tmp, tmp2, tv[162][162][5];
	double delunm[5];
	/* --------------------------------------------------------------------- */
	/* begin pseudo-time stepping iterations */
	/* --------------------------------------------------------------------- */
	tmp=(1.0/(omega*(2.0-omega)));
	/* --------------------------------------------------------------------- */
	/* initialize a,b,c,d to zero (guarantees that page tables have been */
	/* formed, if applicable on given architecture, before timestepping). */
	/* --------------------------------------------------------------------- */
	#pragma cetus private(i, j, m, n) 
	#pragma loop name ssor#0 
	#pragma cetus parallel 
	#pragma omp parallel for if((10000<((((((1L+(3L*jend))+(-3L*jst))+((168L*iend)*jend))+((-168L*iend)*jst))+((-168L*ist)*jend))+((168L*ist)*jst)))) private(i, j, m, n)
	for (j=jst; j<jend; j ++ )
	{
		#pragma cetus private(i, m, n) 
		#pragma loop name ssor#0#0 
		for (i=ist; i<iend; i ++ )
		{
			#pragma cetus private(m, n) 
			#pragma loop name ssor#0#0#0 
			for (n=0; n<5; n ++ )
			{
				#pragma cetus private(m) 
				#pragma loop name ssor#0#0#0#0 
				for (m=0; m<5; m ++ )
				{
					a[j][i][n][m]=0.0;
					b[j][i][n][m]=0.0;
					c[j][i][n][m]=0.0;
					d[j][i][n][m]=0.0;
				}
			}
		}
	}
	#pragma cetus private(i, j, m, n) 
	#pragma loop name ssor#1 
	#pragma cetus parallel 
	#pragma omp parallel for if((10000<((((((((667L+(-336L*iend))+(336L*ist))+(-333L*jend))+(333L*jst))+((168L*iend)*jend))+((-168L*iend)*jst))+((-168L*ist)*jend))+((168L*ist)*jst)))) private(i, j, m, n)
	for (j=(jend-1); j>=jst; j -- )
	{
		#pragma cetus private(i, m, n) 
		#pragma loop name ssor#1#0 
		for (i=(iend-1); i>=ist; i -- )
		{
			#pragma cetus private(m, n) 
			#pragma loop name ssor#1#0#0 
			for (n=0; n<5; n ++ )
			{
				#pragma cetus private(m) 
				#pragma loop name ssor#1#0#0#0 
				for (m=0; m<5; m ++ )
				{
					au[j][i][n][m]=0.0;
					bu[j][i][n][m]=0.0;
					cu[j][i][n][m]=0.0;
					du[j][i][n][m]=0.0;
				}
			}
		}
	}
	#pragma cetus private(i) 
	#pragma loop name ssor#2 
	for (i=1; i<=11; i ++ )
	{
		timer_clear(i);
	}
	/* --------------------------------------------------------------------- */
	/* compute the steady-state residuals */
	/* --------------------------------------------------------------------- */
	rhs();
	/* --------------------------------------------------------------------- */
	/* compute the L2 norms of newton iteration residuals */
	/* --------------------------------------------------------------------- */
	l2norm(162, 162, 162, nx0, ny0, nz0, ist, iend, jst, jend, rsd, rsdnm);
	#pragma cetus private(i) 
	#pragma loop name ssor#3 
	for (i=1; i<=11; i ++ )
	{
		timer_clear(i);
	}
	timer_start(1);
	/* --------------------------------------------------------------------- */
	/* the timestep loop */
	/* --------------------------------------------------------------------- */
	#pragma cetus private(i, j, k, m) 
	#pragma loop name ssor#4 
	/* #pragma cetus reduction(+: u[k][j][i][m])  */
	for (istep=1; istep<=niter; istep ++ )
	{
		if ((((istep%20)==0)||(istep==itmax))||(istep==1))
		{
			if (niter>1)
			{
				printf(" Time step %4d\n", istep);
			}
		}
		/* --------------------------------------------------------------------- */
		/* perform SSOR iteration */
		/* --------------------------------------------------------------------- */
		if (timeron)
		{
			timer_start(5);
		}
		tmp2=dt;
		#pragma cetus private(i, j, k, m) 
		#pragma loop name ssor#4#0 
		#pragma cetus parallel 
		#pragma omp parallel for if((10000<(((((((((((((-5L+(-6L*jend))+(6L*jst))+(3L*nz))+((-36L*iend)*jend))+((36L*iend)*jst))+((36L*ist)*jend))+((-36L*ist)*jst))+((3L*jend)*nz))+((-3L*jst)*nz))+(((18L*iend)*jend)*nz))+(((-18L*iend)*jst)*nz))+(((-18L*ist)*jend)*nz))+(((18L*ist)*jst)*nz)))) private(i, j, k, m)
		for (k=1; k<(nz-1); k ++ )
		{
			#pragma cetus private(i, j, m) 
			#pragma loop name ssor#4#0#0 
			for (j=jst; j<jend; j ++ )
			{
				#pragma cetus private(i, m) 
				#pragma loop name ssor#4#0#0#0 
				for (i=ist; i<iend; i ++ )
				{
					#pragma cetus private(m) 
					#pragma loop name ssor#4#0#0#0#0 
					for (m=0; m<5; m ++ )
					{
						rsd[k][j][i][m]=(tmp2*rsd[k][j][i][m]);
					}
				}
			}
		}
		if (timeron)
		{
			timer_stop(5);
		}
		mthreadnum=0;
		mthreadnum=(omp_get_num_threads()-1);
		if (mthreadnum>(jend-jst))
		{
			mthreadnum=(jend-jst);
		}
		iam=0;
		iam=omp_get_thread_num();
		if (iam<=mthreadnum)
		{
			isync[iam]=0;
		}
		#pragma cetus private(k) 
		#pragma loop name ssor#4#1 
		for (k=1; k<(nz-1); k ++ )
		{
			/* --------------------------------------------------------------------- */
			/* form the lower triangular part of the jacobian matrix */
			/* --------------------------------------------------------------------- */
			if (timeron)
			{
				timer_start(6);
			}
			jacld(k);
			if (timeron)
			{
				timer_stop(6);
			}
			/* --------------------------------------------------------------------- */
			/* perform the lower triangular solution */
			/* --------------------------------------------------------------------- */
			if (timeron)
			{
				timer_start(7);
			}
			blts(162, 162, 162, nx, ny, nz, k, omega, rsd, a, b, c, d, ist, iend, jst, jend, nx0, ny0);
			if (timeron)
			{
				timer_stop(7);
			}
		}
		#pragma cetus private(k) 
		#pragma loop name ssor#4#2 
		for (k=(nz-2); k>0; k -- )
		{
			/* --------------------------------------------------------------------- */
			/* form the strictly upper triangular part of the jacobian matrix */
			/* --------------------------------------------------------------------- */
			if (timeron)
			{
				timer_start(8);
			}
			jacu(k);
			if (timeron)
			{
				timer_stop(8);
			}
			/* --------------------------------------------------------------------- */
			/* perform the upper triangular solution */
			/* --------------------------------------------------------------------- */
			if (timeron)
			{
				timer_start(9);
			}
			buts(162, 162, 162, nx, ny, nz, k, omega, rsd, tv, du, au, bu, cu, ist, iend, jst, jend, nx0, ny0);
			if (timeron)
			{
				timer_stop(9);
			}
		}
		/* --------------------------------------------------------------------- */
		/* update the variables */
		/* --------------------------------------------------------------------- */
		if (timeron)
		{
			timer_start(10);
		}
		tmp2=tmp;
		#pragma cetus private(i, j, k, m) 
		#pragma loop name ssor#4#3 
		#pragma cetus parallel 
		#pragma omp parallel for if((10000<(((((((((((((-5L+(-6L*jend))+(6L*jst))+(3L*nz))+((-36L*iend)*jend))+((36L*iend)*jst))+((36L*ist)*jend))+((-36L*ist)*jst))+((3L*jend)*nz))+((-3L*jst)*nz))+(((18L*iend)*jend)*nz))+(((-18L*iend)*jst)*nz))+(((-18L*ist)*jend)*nz))+(((18L*ist)*jst)*nz)))) private(i, j, k, m)
		for (k=1; k<(nz-1); k ++ )
		{
			#pragma cetus private(i, j, m) 
			#pragma loop name ssor#4#3#0 
			for (j=jst; j<jend; j ++ )
			{
				#pragma cetus private(i, m) 
				#pragma loop name ssor#4#3#0#0 
				for (i=ist; i<iend; i ++ )
				{
					#pragma cetus private(m) 
					#pragma loop name ssor#4#3#0#0#0 
					for (m=0; m<5; m ++ )
					{
						u[k][j][i][m]=(u[k][j][i][m]+(tmp2*rsd[k][j][i][m]));
					}
				}
			}
		}
		if (timeron)
		{
			timer_stop(10);
		}
		/* --------------------------------------------------------------------- */
		/* compute the max-norms of newton iteration corrections */
		/* --------------------------------------------------------------------- */
		if ((istep%inorm)==0)
		{
			if (timeron)
			{
				timer_start(11);
			}
			l2norm(162, 162, 162, nx0, ny0, nz0, ist, iend, jst, jend, rsd, delunm);
			if (timeron)
			{
				timer_stop(11);
			}
			/*
			
			      if ( ipr == 1 ) {
				        printf(" \n RMS-norm of SSOR-iteration correction "
				               "for first pde  = %12.5E\n"
				               " RMS-norm of SSOR-iteration correction "
				               "for second pde = %12.5E\n"
				               " RMS-norm of SSOR-iteration correction "
				               "for third pde  = %12.5E\n"
				               " RMS-norm of SSOR-iteration correction "
				               "for fourth pde = %12.5E\n",
				               " RMS-norm of SSOR-iteration correction "
				               "for fifth pde  = %12.5E\n", 
				               delunm[0], delunm[1], delunm[2], delunm[3], delunm[4]); 
			      } else if ( ipr == 2 ) {
				        printf("(%5d,%15.6f)\n", istep, delunm[4]);
			      }
			     
			*/
		}
		/* --------------------------------------------------------------------- */
		/* compute the steady-state residuals */
		/* --------------------------------------------------------------------- */
		rhs();
		/* --------------------------------------------------------------------- */
		/* compute the max-norms of newton iteration residuals */
		/* --------------------------------------------------------------------- */
		if (((istep%inorm)==0)||(istep==itmax))
		{
			if (timeron)
			{
				timer_start(11);
			}
			l2norm(162, 162, 162, nx0, ny0, nz0, ist, iend, jst, jend, rsd, rsdnm);
			if (timeron)
			{
				timer_stop(11);
			}
			/*
			
			      if ( ipr == 1 ) {
				        printf(" \n RMS-norm of steady-state residual for "
				               "first pde  = %12.5E\n"
				               " RMS-norm of steady-state residual for "
				               "second pde = %12.5E\n"
				               " RMS-norm of steady-state residual for "
				               "third pde  = %12.5E\n"
				               " RMS-norm of steady-state residual for "
				               "fourth pde = %12.5E\n"
				               " RMS-norm of steady-state residual for "
				               "fifth pde  = %12.5E\n", 
				               rsdnm[0], rsdnm[1], rsdnm[2], rsdnm[3], rsdnm[4]);
			      }
			     
			*/
		}
		/* --------------------------------------------------------------------- */
		/* check the newton-iteration residuals against the tolerance levels */
		/* --------------------------------------------------------------------- */
		if (((((rsdnm[0]<tolrsd[0])&&(rsdnm[1]<tolrsd[1]))&&(rsdnm[2]<tolrsd[2]))&&(rsdnm[3]<tolrsd[3]))&&(rsdnm[4]<tolrsd[4]))
		{
			/* if (ipr == 1 ) { */
				printf(" \n convergence was achieved after %4d pseudo-time steps\n", istep);
			/* } */
			break;
		}
	}
	timer_stop(1);
	maxtime=timer_read(1);
}
