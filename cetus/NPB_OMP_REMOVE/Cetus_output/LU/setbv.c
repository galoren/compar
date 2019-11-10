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
/* set the boundary values of dependent variables */
/* --------------------------------------------------------------------- */
void setbv()
{
	/* --------------------------------------------------------------------- */
	/* local variables */
	/* --------------------------------------------------------------------- */
	int i, j, k, m;
	double temp1[5], temp2[5];
	/* --------------------------------------------------------------------- */
	/* set the dependent variable values along the top and bottom faces */
	/* --------------------------------------------------------------------- */
	#pragma cetus private(i, j, m) 
	#pragma loop name setbv#0 
	for (j=0; j<ny; j ++ )
	{
		#pragma cetus private(i, m) 
		#pragma loop name setbv#0#0 
		for (i=0; i<nx; i ++ )
		{
			exact(i, j, 0, temp1);
			exact(i, j, nz-1, temp2);
			#pragma cetus private(m) 
			#pragma loop name setbv#0#0#0 
			#pragma cetus parallel 
			/*
			Disabled due to low profitability: #pragma omp parallel for private(m)
			*/
			for (m=0; m<5; m ++ )
			{
				u[0][j][i][m]=temp1[m];
				u[nz-1][j][i][m]=temp2[m];
			}
		}
	}
	/* --------------------------------------------------------------------- */
	/* set the dependent variable values along north and south faces */
	/* --------------------------------------------------------------------- */
	#pragma cetus private(i, k, m) 
	#pragma loop name setbv#1 
	for (k=0; k<nz; k ++ )
	{
		#pragma cetus private(i, m) 
		#pragma loop name setbv#1#0 
		for (i=0; i<nx; i ++ )
		{
			exact(i, 0, k, temp1);
			exact(i, ny-1, k, temp2);
			#pragma cetus private(m) 
			#pragma loop name setbv#1#0#0 
			#pragma cetus parallel 
			/*
			Disabled due to low profitability: #pragma omp parallel for private(m)
			*/
			for (m=0; m<5; m ++ )
			{
				u[k][0][i][m]=temp1[m];
				u[k][ny-1][i][m]=temp2[m];
			}
		}
	}
	/* --------------------------------------------------------------------- */
	/* set the dependent variable values along east and west faces */
	/* --------------------------------------------------------------------- */
	#pragma cetus private(j, k, m) 
	#pragma loop name setbv#2 
	for (k=0; k<nz; k ++ )
	{
		#pragma cetus private(j, m) 
		#pragma loop name setbv#2#0 
		for (j=0; j<ny; j ++ )
		{
			exact(0, j, k, temp1);
			exact(nx-1, j, k, temp2);
			#pragma cetus private(m) 
			#pragma loop name setbv#2#0#0 
			#pragma cetus parallel 
			/*
			Disabled due to low profitability: #pragma omp parallel for private(m)
			*/
			for (m=0; m<5; m ++ )
			{
				u[k][j][0][m]=temp1[m];
				u[k][j][nx-1][m]=temp2[m];
			}
		}
	}
}
