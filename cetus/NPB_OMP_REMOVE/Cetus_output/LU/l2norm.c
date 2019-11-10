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
#include <math.h>
#include "applu.incl"
/* --------------------------------------------------------------------- */
/* to compute the l2-norm of vector v. */
/* --------------------------------------------------------------------- */
/* --------------------------------------------------------------------- */
/* To improve cache performance, second two dimensions padded by 1  */
/* for even number sizes only.  Only needed in v. */
/* --------------------------------------------------------------------- */
void l2norm(int ldx, int ldy, int ldz, int nx0, int ny0, int nz0, int ist, int iend, int jst, int jend, double v[][(((ldy/2)*2)+1)][(((ldx/2)*2)+1)][5], double sum[5])
{
	/* --------------------------------------------------------------------- */
	/* local variables */
	/* --------------------------------------------------------------------- */
	double sum_local[5];
	int i, j, k, m;
	#pragma cetus private(m) 
	#pragma loop name l2norm#0 
	#pragma cetus parallel 
	/*
	Disabled due to low profitability: #pragma omp parallel for private(m)
	*/
	for (m=0; m<5; m ++ )
	{
		sum[m]=0.0;
	}
	#pragma cetus private(m) 
	#pragma loop name l2norm#1 
	#pragma cetus parallel 
	/*
	Disabled due to low profitability: #pragma omp parallel for private(m)
	*/
	for (m=0; m<5; m ++ )
	{
		sum_local[m]=0.0;
	}
	#pragma cetus parallel 
	#pragma cetus private(i, j, k, m) 
	#pragma omp parallel if((10000<(((((((((((((129L+(-6L*jend))+(6L*jst))+(3L*nz0))+((-36L*iend)*jend))+((36L*iend)*jst))+((36L*ist)*jend))+((-36L*ist)*jst))+((3L*jend)*nz0))+((-3L*jst)*nz0))+(((18L*iend)*jend)*nz0))+(((-18L*iend)*jst)*nz0))+(((-18L*ist)*jend)*nz0))+(((18L*ist)*jst)*nz0)))) private(i, j, k, m)
	{
		double * reduce = (double * )malloc(5*sizeof (double));
		int reduce_span_0;
		for (reduce_span_0=0; reduce_span_0<5; reduce_span_0 ++ )
		{
			reduce[reduce_span_0]=0;
		}
		#pragma loop name l2norm#2 
		#pragma cetus for  
		#pragma omp for
		for (k=1; k<(nz0-1); k ++ )
		{
			#pragma cetus private(i, j, m) 
			#pragma loop name l2norm#2#0 
			/* #pragma cetus reduction(+: sum_local[m])  */
			for (j=jst; j<jend; j ++ )
			{
				#pragma cetus private(i, m) 
				#pragma loop name l2norm#2#0#0 
				/* #pragma cetus reduction(+: sum_local[m])  */
				for (i=ist; i<iend; i ++ )
				{
					#pragma cetus private(m) 
					#pragma loop name l2norm#2#0#0#0 
					for (m=0; m<5; m ++ )
					{
						reduce[m]=(reduce[m]+(v[k][j][i][m]*v[k][j][i][m]));
					}
				}
			}
		}
		#pragma cetus critical  
		#pragma omp critical
		{
			for (reduce_span_0=0; reduce_span_0<5; reduce_span_0 ++ )
			{
				sum_local[reduce_span_0]+=reduce[reduce_span_0];
			}
		}
	}
	#pragma cetus private(m) 
	#pragma loop name l2norm#3 
	#pragma cetus parallel 
	/*
	Disabled due to low profitability: #pragma omp parallel for private(m)
	*/
	for (m=0; m<5; m ++ )
	{
		sum[m]+=sum_local[m];
	}
	/* end parallel */
	#pragma cetus private(m) 
	#pragma loop name l2norm#4 
	#pragma cetus parallel 
	/*
	Disabled due to low profitability: #pragma omp parallel for private(m)
	*/
	for (m=0; m<5; m ++ )
	{
		sum[m]=sqrt(sum[m]/(((nx0-2)*(ny0-2))*(nz0-2)));
	}
}
