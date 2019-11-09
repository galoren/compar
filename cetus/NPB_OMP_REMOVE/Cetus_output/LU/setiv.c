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
#include "applu.incl"
/* --------------------------------------------------------------------- */
/*  */
/* set the initial values of independent variables based on tri-linear */
/* interpolation of boundary values in the computational space. */
/*  */
/* --------------------------------------------------------------------- */
void setiv()
{
	/* --------------------------------------------------------------------- */
	/* local variables */
	/* --------------------------------------------------------------------- */
	int i, j, k, m;
	double xi, eta, zeta;
	double pxi, peta, pzeta;
	double ue_1jk[5], ue_nx0jk[5], ue_i1k[5];
	double ue_iny0k[5], ue_ij1[5], ue_ijnz[5];
	#pragma cetus private(eta, i, j, k, m, peta, pxi, pzeta, xi, zeta) 
	#pragma loop name setiv#0 
	for (k=1; k<(nz-1); k ++ )
	{
		zeta=(((double)k)/(nz-1));
		#pragma cetus private(eta, i, j, m, peta, pxi, pzeta, xi) 
		#pragma loop name setiv#0#0 
		for (j=1; j<(ny-1); j ++ )
		{
			eta=(((double)j)/(ny0-1));
			#pragma cetus private(i, m, peta, pxi, pzeta, xi) 
			#pragma loop name setiv#0#0#0 
			for (i=1; i<(nx-1); i ++ )
			{
				xi=(((double)i)/(nx0-1));
				exact(0, j, k, ue_1jk);
				exact(nx0-1, j, k, ue_nx0jk);
				exact(i, 0, k, ue_i1k);
				exact(i, ny0-1, k, ue_iny0k);
				exact(i, j, 0, ue_ij1);
				exact(i, j, nz-1, ue_ijnz);
				#pragma cetus private(m, peta, pxi, pzeta) 
				#pragma loop name setiv#0#0#0#0 
				#pragma cetus parallel 
				/*
				Disabled due to low profitability: #pragma omp parallel for private(m, peta, pxi, pzeta)
				*/
				for (m=0; m<5; m ++ )
				{
					pxi=(((1.0-xi)*ue_1jk[m])+(xi*ue_nx0jk[m]));
					peta=(((1.0-eta)*ue_i1k[m])+(eta*ue_iny0k[m]));
					pzeta=(((1.0-zeta)*ue_ij1[m])+(zeta*ue_ijnz[m]));
					u[k][j][i][m]=((((((pxi+peta)+pzeta)-(pxi*peta))-(peta*pzeta))-(pzeta*pxi))+((pxi*peta)*pzeta));
				}
			}
		}
	}
}
