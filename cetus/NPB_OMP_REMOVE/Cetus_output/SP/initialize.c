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
/*  This benchmark is an OpenMP C version of the NPB SP code. This OpenMP   */
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
#include "header.h"
/* --------------------------------------------------------------------- */
/* This subroutine initializes the field variable u using  */
/* tri-linear transfinite interpolation of the boundary values      */
/* --------------------------------------------------------------------- */
void initialize()
{
	int i, j, k, m, ix, iy, iz;
	double xi, eta, zeta, Pface[2][3][5], Pxi, Peta, Pzeta, temp[5];
	/* --------------------------------------------------------------------- */
	/*  Later (in compute_rhs) we compute 1u for every element. A few of  */
	/*  the corner elements are not used, but it convenient (and faster)  */
	/*  to compute the whole thing with a simple loop. Make sure those  */
	/*  values are nonzero by initializing the whole thing here.  */
	/* --------------------------------------------------------------------- */
	#pragma cetus private(i, j, k) 
	#pragma loop name initialize#0 
	#pragma cetus parallel 
	#pragma omp parallel for if((10000<(((1L+(3L*grid_points[2L]))+((3L*grid_points[1L])*grid_points[2L]))+(((7L*grid_points[0L])*grid_points[1L])*grid_points[2L])))) private(i, j, k)
	for (k=0; k<=(grid_points[2]-1); k ++ )
	{
		#pragma cetus private(i, j) 
		#pragma loop name initialize#0#0 
		for (j=0; j<=(grid_points[1]-1); j ++ )
		{
			#pragma cetus private(i) 
			#pragma loop name initialize#0#0#0 
			for (i=0; i<=(grid_points[0]-1); i ++ )
			{
				u[k][j][i][0]=1.0;
				u[k][j][i][1]=0.0;
				u[k][j][i][2]=0.0;
				u[k][j][i][3]=0.0;
				u[k][j][i][4]=1.0;
			}
		}
	}
	/* --------------------------------------------------------------------- */
	/* first store the "interpolated" values everywhere on the grid     */
	/* --------------------------------------------------------------------- */
	#pragma cetus private(Peta, Pxi, Pzeta, eta, i, ix, iy, iz, j, k, m, xi, zeta) 
	#pragma loop name initialize#1 
	for (k=0; k<=(grid_points[2]-1); k ++ )
	{
		zeta=(((double)k)*dnzm1);
		#pragma cetus private(Peta, Pxi, Pzeta, eta, i, ix, iy, iz, j, m, xi) 
		#pragma loop name initialize#1#0 
		for (j=0; j<=(grid_points[1]-1); j ++ )
		{
			eta=(((double)j)*dnym1);
			#pragma cetus private(Peta, Pxi, Pzeta, i, ix, iy, iz, m, xi) 
			#pragma loop name initialize#1#0#0 
			for (i=0; i<=(grid_points[0]-1); i ++ )
			{
				xi=(((double)i)*dnxm1);
				#pragma cetus private(Pxi, ix) 
				#pragma loop name initialize#1#0#0#0 
				for (ix=0; ix<2; ix ++ )
				{
					Pxi=((double)ix);
					exact_solution(Pxi, eta, zeta,  & Pface[ix][0][0]);
				}
				#pragma cetus private(Peta, iy) 
				#pragma loop name initialize#1#0#0#1 
				for (iy=0; iy<2; iy ++ )
				{
					Peta=((double)iy);
					exact_solution(xi, Peta, zeta,  & Pface[iy][1][0]);
				}
				#pragma cetus private(Pzeta, iz) 
				#pragma loop name initialize#1#0#0#2 
				for (iz=0; iz<2; iz ++ )
				{
					Pzeta=((double)iz);
					exact_solution(xi, eta, Pzeta,  & Pface[iz][2][0]);
				}
				#pragma cetus private(Peta, Pxi, Pzeta, m) 
				#pragma loop name initialize#1#0#0#3 
				#pragma cetus parallel 
				/*
				Disabled due to low profitability: #pragma omp parallel for private(Peta, Pxi, Pzeta, m)
				*/
				for (m=0; m<5; m ++ )
				{
					Pxi=((xi*Pface[1][0][m])+((1.0-xi)*Pface[0][0][m]));
					Peta=((eta*Pface[1][1][m])+((1.0-eta)*Pface[0][1][m]));
					Pzeta=((zeta*Pface[1][2][m])+((1.0-zeta)*Pface[0][2][m]));
					u[k][j][i][m]=((((((Pxi+Peta)+Pzeta)-(Pxi*Peta))-(Pxi*Pzeta))-(Peta*Pzeta))+((Pxi*Peta)*Pzeta));
				}
			}
		}
	}
	/* --------------------------------------------------------------------- */
	/* now store the exact values on the boundaries         */
	/* --------------------------------------------------------------------- */
	/* --------------------------------------------------------------------- */
	/* west face                                                   */
	/* --------------------------------------------------------------------- */
	xi=0.0;
	i=0;
	#pragma cetus private(eta, j, k, m, zeta) 
	#pragma loop name initialize#2 
	for (k=0; k<=(grid_points[2]-1); k ++ )
	{
		zeta=(((double)k)*dnzm1);
		#pragma cetus private(eta, j, m) 
		#pragma loop name initialize#2#0 
		for (j=0; j<=(grid_points[1]-1); j ++ )
		{
			eta=(((double)j)*dnym1);
			exact_solution(xi, eta, zeta, temp);
			#pragma cetus private(m) 
			#pragma loop name initialize#2#0#0 
			#pragma cetus parallel 
			/*
			Disabled due to low profitability: #pragma omp parallel for private(m)
			*/
			for (m=0; m<5; m ++ )
			{
				u[k][j][i][m]=temp[m];
			}
		}
	}
	/* --------------------------------------------------------------------- */
	/* east face                                                       */
	/* --------------------------------------------------------------------- */
	xi=1.0;
	i=(grid_points[0]-1);
	#pragma cetus private(eta, j, k, m, zeta) 
	#pragma loop name initialize#3 
	for (k=0; k<=(grid_points[2]-1); k ++ )
	{
		zeta=(((double)k)*dnzm1);
		#pragma cetus private(eta, j, m) 
		#pragma loop name initialize#3#0 
		for (j=0; j<=(grid_points[1]-1); j ++ )
		{
			eta=(((double)j)*dnym1);
			exact_solution(xi, eta, zeta, temp);
			#pragma cetus private(m) 
			#pragma loop name initialize#3#0#0 
			#pragma cetus parallel 
			/*
			Disabled due to low profitability: #pragma omp parallel for private(m)
			*/
			for (m=0; m<5; m ++ )
			{
				u[k][j][i][m]=temp[m];
			}
		}
	}
	/* --------------------------------------------------------------------- */
	/* south face                                                  */
	/* --------------------------------------------------------------------- */
	eta=0.0;
	j=0;
	#pragma cetus private(i, k, m, xi, zeta) 
	#pragma loop name initialize#4 
	for (k=0; k<=(grid_points[2]-1); k ++ )
	{
		zeta=(((double)k)*dnzm1);
		#pragma cetus private(i, m, xi) 
		#pragma loop name initialize#4#0 
		for (i=0; i<=(grid_points[0]-1); i ++ )
		{
			xi=(((double)i)*dnxm1);
			exact_solution(xi, eta, zeta, temp);
			#pragma cetus private(m) 
			#pragma loop name initialize#4#0#0 
			#pragma cetus parallel 
			/*
			Disabled due to low profitability: #pragma omp parallel for private(m)
			*/
			for (m=0; m<5; m ++ )
			{
				u[k][j][i][m]=temp[m];
			}
		}
	}
	/* --------------------------------------------------------------------- */
	/* north face                                     */
	/* --------------------------------------------------------------------- */
	eta=1.0;
	j=(grid_points[1]-1);
	#pragma cetus private(i, k, m, xi, zeta) 
	#pragma loop name initialize#5 
	for (k=0; k<=(grid_points[2]-1); k ++ )
	{
		zeta=(((double)k)*dnzm1);
		#pragma cetus private(i, m, xi) 
		#pragma loop name initialize#5#0 
		for (i=0; i<=(grid_points[0]-1); i ++ )
		{
			xi=(((double)i)*dnxm1);
			exact_solution(xi, eta, zeta, temp);
			#pragma cetus private(m) 
			#pragma loop name initialize#5#0#0 
			#pragma cetus parallel 
			/*
			Disabled due to low profitability: #pragma omp parallel for private(m)
			*/
			for (m=0; m<5; m ++ )
			{
				u[k][j][i][m]=temp[m];
			}
		}
	}
	/* --------------------------------------------------------------------- */
	/* bottom face                                        */
	/* --------------------------------------------------------------------- */
	zeta=0.0;
	k=0;
	#pragma cetus private(eta, i, j, m, xi) 
	#pragma loop name initialize#6 
	for (j=0; j<=(grid_points[1]-1); j ++ )
	{
		eta=(((double)j)*dnym1);
		#pragma cetus private(i, m, xi) 
		#pragma loop name initialize#6#0 
		for (i=0; i<=(grid_points[0]-1); i ++ )
		{
			xi=(((double)i)*dnxm1);
			exact_solution(xi, eta, zeta, temp);
			#pragma cetus private(m) 
			#pragma loop name initialize#6#0#0 
			#pragma cetus parallel 
			/*
			Disabled due to low profitability: #pragma omp parallel for private(m)
			*/
			for (m=0; m<5; m ++ )
			{
				u[k][j][i][m]=temp[m];
			}
		}
	}
	/* --------------------------------------------------------------------- */
	/* top face      */
	/* --------------------------------------------------------------------- */
	zeta=1.0;
	k=(grid_points[2]-1);
	#pragma cetus private(eta, i, j, m, xi) 
	#pragma loop name initialize#7 
	for (j=0; j<=(grid_points[1]-1); j ++ )
	{
		eta=(((double)j)*dnym1);
		#pragma cetus private(i, m, xi) 
		#pragma loop name initialize#7#0 
		for (i=0; i<=(grid_points[0]-1); i ++ )
		{
			xi=(((double)i)*dnxm1);
			exact_solution(xi, eta, zeta, temp);
			#pragma cetus private(m) 
			#pragma loop name initialize#7#0#0 
			#pragma cetus parallel 
			/*
			Disabled due to low profitability: #pragma omp parallel for private(m)
			*/
			for (m=0; m<5; m ++ )
			{
				u[k][j][i][m]=temp[m];
			}
		}
	}
}

void lhsinit(int ni, int nj)
{
	int j, m;
	/* --------------------------------------------------------------------- */
	/* zap the whole left hand side for starters */
	/* set all diagonal values to 1. This is overkill, but convenient */
	/* --------------------------------------------------------------------- */
	#pragma cetus private(j, m) 
	#pragma loop name lhsinit#0 
	#pragma cetus parallel 
	#pragma omp parallel for if((10000<(1L+(49L*nj)))) private(j, m)
	for (j=1; j<=nj; j ++ )
	{
		#pragma cetus private(m) 
		#pragma loop name lhsinit#0#0 
		for (m=0; m<5; m ++ )
		{
			lhs[j][0][m]=0.0;
			lhsp[j][0][m]=0.0;
			lhsm[j][0][m]=0.0;
			lhs[j][ni][m]=0.0;
			lhsp[j][ni][m]=0.0;
			lhsm[j][ni][m]=0.0;
		}
		lhs[j][0][2]=1.0;
		lhsp[j][0][2]=1.0;
		lhsm[j][0][2]=1.0;
		lhs[j][ni][2]=1.0;
		lhsp[j][ni][2]=1.0;
		lhsm[j][ni][2]=1.0;
	}
}

void lhsinitj(int nj, int ni)
{
	int i, m;
	/* --------------------------------------------------------------------- */
	/* zap the whole left hand side for starters */
	/* set all diagonal values to 1. This is overkill, but convenient */
	/* --------------------------------------------------------------------- */
	#pragma cetus private(i, m) 
	#pragma loop name lhsinitj#0 
	#pragma cetus parallel 
	#pragma omp parallel for if((10000<(1L+(49L*ni)))) private(i, m)
	for (i=1; i<=ni; i ++ )
	{
		#pragma cetus private(m) 
		#pragma loop name lhsinitj#0#0 
		for (m=0; m<5; m ++ )
		{
			lhs[0][i][m]=0.0;
			lhsp[0][i][m]=0.0;
			lhsm[0][i][m]=0.0;
			lhs[nj][i][m]=0.0;
			lhsp[nj][i][m]=0.0;
			lhsm[nj][i][m]=0.0;
		}
		lhs[0][i][2]=1.0;
		lhsp[0][i][2]=1.0;
		lhsm[0][i][2]=1.0;
		lhs[nj][i][2]=1.0;
		lhsp[nj][i][2]=1.0;
		lhsm[nj][i][2]=1.0;
	}
}
