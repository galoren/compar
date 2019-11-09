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
/*  This benchmark is an OpenMP C version of the NPB BT code. This OpenMP   */
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
#include "header.h"
/* --------------------------------------------------------------------- */
/* this function computes the norm of the difference between the */
/* computed solution and the exact solution */
/* --------------------------------------------------------------------- */
void error_norm(double rms[5])
{
	int i, j, k, m, d;
	double xi, eta, zeta, u_exact[5], add;
	double rms_local[5];
	#pragma cetus private(m) 
	#pragma loop name error_norm#0 
	#pragma cetus parallel 
	/*
	Disabled due to low profitability: #pragma omp parallel for private(m)
	*/
	for (m=0; m<5; m ++ )
	{
		rms[m]=0.0;
	}
	#pragma cetus private(m) 
	#pragma loop name error_norm#1 
	#pragma cetus parallel 
	/*
	Disabled due to low profitability: #pragma omp parallel for private(m)
	*/
	for (m=0; m<5; m ++ )
	{
		rms_local[m]=0.0;
	}
	#pragma cetus private(add, eta, i, j, k, m, xi, zeta) 
	#pragma loop name error_norm#2 
	/* #pragma cetus reduction(+: rms_local[m])  */
	for (k=0; k<=(grid_points[2]-1); k ++ )
	{
		zeta=(((double)k)*dnzm1);
		#pragma cetus private(add, eta, i, j, m, xi) 
		#pragma loop name error_norm#2#0 
		/* #pragma cetus reduction(+: rms_local[m])  */
		for (j=0; j<=(grid_points[1]-1); j ++ )
		{
			eta=(((double)j)*dnym1);
			#pragma cetus private(add, i, m, xi) 
			#pragma loop name error_norm#2#0#0 
			/* #pragma cetus reduction(+: rms_local[m])  */
			for (i=0; i<=(grid_points[0]-1); i ++ )
			{
				xi=(((double)i)*dnxm1);
				exact_solution(xi, eta, zeta, u_exact);
				#pragma cetus private(add, m) 
				#pragma loop name error_norm#2#0#0#0 
				#pragma cetus parallel 
				/*
				Disabled due to low profitability: #pragma omp parallel for private(add, m)
				*/
				for (m=0; m<5; m ++ )
				{
					add=(u[k][j][i][m]-u_exact[m]);
					rms_local[m]=(rms_local[m]+(add*add));
				}
			}
		}
	}
	#pragma cetus private(m) 
	#pragma loop name error_norm#3 
	#pragma cetus parallel 
	/*
	Disabled due to low profitability: #pragma omp parallel for private(m)
	*/
	for (m=0; m<5; m ++ )
	{
		rms[m]+=rms_local[m];
	}
	#pragma cetus private(d, m) 
	#pragma loop name error_norm#4 
	#pragma cetus parallel 
	/*
	Disabled due to low profitability: #pragma omp parallel for private(d, m)
	*/
	for (m=0; m<5; m ++ )
	{
		#pragma cetus private(d) 
		#pragma loop name error_norm#4#0 
		for (d=0; d<3; d ++ )
		{
			rms[m]=(rms[m]/((double)(grid_points[d]-2)));
		}
		rms[m]=sqrt(rms[m]);
	}
}

void rhs_norm(double rms[5])
{
	int i, j, k, d, m;
	double add;
	double rms_local[5];
	#pragma cetus private(m) 
	#pragma loop name rhs_norm#0 
	#pragma cetus parallel 
	/*
	Disabled due to low profitability: #pragma omp parallel for private(m)
	*/
	for (m=0; m<5; m ++ )
	{
		rms[m]=0.0;
	}
	#pragma cetus private(m) 
	#pragma loop name rhs_norm#1 
	#pragma cetus parallel 
	/*
	Disabled due to low profitability: #pragma omp parallel for private(m)
	*/
	for (m=0; m<5; m ++ )
	{
		rms_local[m]=0.0;
	}
	#pragma cetus parallel 
	#pragma cetus private(add, i, j, k, m) 
	#pragma omp parallel if((10000<(((((((-43L+(92L*grid_points[0L]))+(86L*grid_points[1L]))+(89L*grid_points[2L]))+((-46L*grid_points[0L])*grid_points[1L]))+((-46L*grid_points[0L])*grid_points[2L]))+((-43L*grid_points[1L])*grid_points[2L]))+(((23L*grid_points[0L])*grid_points[1L])*grid_points[2L])))) private(add, i, j, k, m)
	{
		double * reduce = (double * )malloc(5*sizeof (double));
		int reduce_span_0;
		for (reduce_span_0=0; reduce_span_0<5; reduce_span_0 ++ )
		{
			reduce[reduce_span_0]=0;
		}
		#pragma loop name rhs_norm#2 
		#pragma cetus for  
		#pragma omp for
		for (k=1; k<=(grid_points[2]-2); k ++ )
		{
			#pragma cetus private(add, i, j, m) 
			#pragma loop name rhs_norm#2#0 
			/* #pragma cetus reduction(+: rms_local[m])  */
			for (j=1; j<=(grid_points[1]-2); j ++ )
			{
				#pragma cetus private(add, i, m) 
				#pragma loop name rhs_norm#2#0#0 
				/* #pragma cetus reduction(+: rms_local[m])  */
				for (i=1; i<=(grid_points[0]-2); i ++ )
				{
					#pragma cetus private(add, m) 
					#pragma loop name rhs_norm#2#0#0#0 
					for (m=0; m<5; m ++ )
					{
						add=rhs[k][j][i][m];
						reduce[m]=(reduce[m]+(add*add));
					}
				}
			}
		}
		#pragma cetus critical  
		#pragma omp critical
		{
			for (reduce_span_0=0; reduce_span_0<5; reduce_span_0 ++ )
			{
				rms_local[reduce_span_0]+=reduce[reduce_span_0];
			}
		}
	}
	#pragma cetus private(m) 
	#pragma loop name rhs_norm#3 
	#pragma cetus parallel 
	/*
	Disabled due to low profitability: #pragma omp parallel for private(m)
	*/
	for (m=0; m<5; m ++ )
	{
		rms[m]+=rms_local[m];
	}
	#pragma cetus private(d, m) 
	#pragma loop name rhs_norm#4 
	#pragma cetus parallel 
	/*
	Disabled due to low profitability: #pragma omp parallel for private(d, m)
	*/
	for (m=0; m<5; m ++ )
	{
		#pragma cetus private(d) 
		#pragma loop name rhs_norm#4#0 
		for (d=0; d<3; d ++ )
		{
			rms[m]=(rms[m]/((double)(grid_points[d]-2)));
		}
		rms[m]=sqrt(rms[m]);
	}
}
