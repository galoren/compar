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
	/* Later (in compute_rhs) we compute 1u for every element. A few of  */
	/* the corner elements are not used, but it convenient (and faster)  */
	/* to compute the whole thing with a simple loop. Make sure those  */
	/* values are nonzero by initializing the whole thing here.  */
	/* --------------------------------------------------------------------- */
	#pragma cetus private(Peta, Pxi, Pzeta, eta, i, ix, iy, iz, j, k, m, xi, zeta) 
	#pragma loop name initialize#0 
	for (k=0; k<=(grid_points[2]-1); k ++ )
	{
		#pragma cetus private(i, j, m) 
		#pragma loop name initialize#0#0 
		#pragma cetus parallel 
		#pragma omp parallel for if((10000<((1L+(3L*grid_points[1L]))+((18L*grid_points[0L])*grid_points[1L])))) private(i, j, m)
		for (j=0; j<=(grid_points[1]-1); j ++ )
		{
			#pragma cetus private(i, m) 
			#pragma loop name initialize#0#0#0 
			for (i=0; i<=(grid_points[0]-1); i ++ )
			{
				#pragma cetus private(m) 
				#pragma loop name initialize#0#0#0#0 
				for (m=0; m<5; m ++ )
				{
					u[k][j][i][m]=1.0;
				}
			}
		}
		/* --------------------------------------------------------------------- */
		/* first store the "interpolated" values everywhere on the grid     */
		/* --------------------------------------------------------------------- */
		#pragma cetus private(Peta, Pxi, Pzeta, eta, i, ix, iy, iz, j, k, m, xi, zeta) 
		#pragma loop name initialize#0#1 
		for (k=0; k<=(grid_points[2]-1); k ++ )
		{
			zeta=(((double)k)*dnzm1);
			#pragma cetus private(Peta, Pxi, Pzeta, eta, i, ix, iy, iz, j, m, xi) 
			#pragma loop name initialize#0#1#0 
			for (j=0; j<=(grid_points[1]-1); j ++ )
			{
				eta=(((double)j)*dnym1);
				#pragma cetus private(Peta, Pxi, Pzeta, i, ix, iy, iz, m, xi) 
				#pragma loop name initialize#0#1#0#0 
				for (i=0; i<=(grid_points[0]-1); i ++ )
				{
					xi=(((double)i)*dnxm1);
					#pragma cetus private(ix) 
					#pragma loop name initialize#0#1#0#0#0 
					for (ix=0; ix<2; ix ++ )
					{
						exact_solution((double)ix, eta, zeta,  & Pface[ix][0][0]);
					}
					#pragma cetus private(iy) 
					#pragma loop name initialize#0#1#0#0#1 
					for (iy=0; iy<2; iy ++ )
					{
						exact_solution(xi, (double)iy, zeta,  & Pface[iy][1][0]);
					}
					#pragma cetus private(iz) 
					#pragma loop name initialize#0#1#0#0#2 
					for (iz=0; iz<2; iz ++ )
					{
						exact_solution(xi, eta, (double)iz,  & Pface[iz][2][0]);
					}
					#pragma cetus private(Peta, Pxi, Pzeta, m) 
					#pragma loop name initialize#0#1#0#0#3 
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
		i=0;
		xi=0.0;
		#pragma cetus private(eta, j, k, m, zeta) 
		#pragma loop name initialize#0#2 
		for (k=0; k<=(grid_points[2]-1); k ++ )
		{
			zeta=(((double)k)*dnzm1);
			#pragma cetus private(eta, j, m) 
			#pragma loop name initialize#0#2#0 
			for (j=0; j<=(grid_points[1]-1); j ++ )
			{
				eta=(((double)j)*dnym1);
				exact_solution(xi, eta, zeta, temp);
				#pragma cetus private(m) 
				#pragma loop name initialize#0#2#0#0 
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
		i=(grid_points[0]-1);
		xi=1.0;
		#pragma cetus private(eta, j, k, m, zeta) 
		#pragma loop name initialize#0#3 
		for (k=0; k<=(grid_points[2]-1); k ++ )
		{
			zeta=(((double)k)*dnzm1);
			#pragma cetus private(eta, j, m) 
			#pragma loop name initialize#0#3#0 
			for (j=0; j<=(grid_points[1]-1); j ++ )
			{
				eta=(((double)j)*dnym1);
				exact_solution(xi, eta, zeta, temp);
				#pragma cetus private(m) 
				#pragma loop name initialize#0#3#0#0 
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
		j=0;
		eta=0.0;
		#pragma cetus private(i, k, m, xi, zeta) 
		#pragma loop name initialize#0#4 
		for (k=0; k<=(grid_points[2]-1); k ++ )
		{
			zeta=(((double)k)*dnzm1);
			#pragma cetus private(i, m, xi) 
			#pragma loop name initialize#0#4#0 
			for (i=0; i<=(grid_points[0]-1); i ++ )
			{
				xi=(((double)i)*dnxm1);
				exact_solution(xi, eta, zeta, temp);
				#pragma cetus private(m) 
				#pragma loop name initialize#0#4#0#0 
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
		j=(grid_points[1]-1);
		eta=1.0;
		#pragma cetus private(i, k, m, xi, zeta) 
		#pragma loop name initialize#0#5 
		for (k=0; k<=(grid_points[2]-1); k ++ )
		{
			zeta=(((double)k)*dnzm1);
			#pragma cetus private(i, m, xi) 
			#pragma loop name initialize#0#5#0 
			for (i=0; i<=(grid_points[0]-1); i ++ )
			{
				xi=(((double)i)*dnxm1);
				exact_solution(xi, eta, zeta, temp);
				#pragma cetus private(m) 
				#pragma loop name initialize#0#5#0#0 
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
		k=0;
		zeta=0.0;
		#pragma cetus private(eta, i, j, m, xi) 
		#pragma loop name initialize#0#6 
		for (j=0; j<=(grid_points[1]-1); j ++ )
		{
			eta=(((double)j)*dnym1);
			#pragma cetus private(i, m, xi) 
			#pragma loop name initialize#0#6#0 
			for (i=0; i<=(grid_points[0]-1); i ++ )
			{
				xi=(((double)i)*dnxm1);
				exact_solution(xi, eta, zeta, temp);
				#pragma cetus private(m) 
				#pragma loop name initialize#0#6#0#0 
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
		k=(grid_points[2]-1);
		zeta=1.0;
		#pragma cetus private(eta, i, j, m, xi) 
		#pragma loop name initialize#0#7 
		for (j=0; j<=(grid_points[1]-1); j ++ )
		{
			eta=(((double)j)*dnym1);
			#pragma cetus private(i, m, xi) 
			#pragma loop name initialize#0#7#0 
			for (i=0; i<=(grid_points[0]-1); i ++ )
			{
				xi=(((double)i)*dnxm1);
				exact_solution(xi, eta, zeta, temp);
				#pragma cetus private(m) 
				#pragma loop name initialize#0#7#0#0 
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
	/* end parallel */
}

void lhsinit(double lhs[][3][5][5], int ni)
{
	int i, m, n;
	/* --------------------------------------------------------------------- */
	/* zero the whole left hand side for starters */
	/* set all diagonal values to 1. This is overkill, but convenient */
	/* --------------------------------------------------------------------- */
	i=0;
	#pragma cetus private(m, n) 
	#pragma loop name lhsinit#0 
	for (n=0; n<5; n ++ )
	{
		#pragma cetus private(m) 
		#pragma loop name lhsinit#0#0 
		#pragma cetus parallel 
		/*
		Disabled due to low profitability: #pragma omp parallel for private(m)
		*/
		for (m=0; m<5; m ++ )
		{
			lhs[i][0][n][m]=0.0;
			lhs[i][1][n][m]=0.0;
			lhs[i][2][n][m]=0.0;
		}
		lhs[i][1][n][n]=1.0;
	}
	i=ni;
	#pragma cetus private(m, n) 
	#pragma loop name lhsinit#1 
	for (n=0; n<5; n ++ )
	{
		#pragma cetus private(m) 
		#pragma loop name lhsinit#1#0 
		#pragma cetus parallel 
		/*
		Disabled due to low profitability: #pragma omp parallel for private(m)
		*/
		for (m=0; m<5; m ++ )
		{
			lhs[i][0][n][m]=0.0;
			lhs[i][1][n][m]=0.0;
			lhs[i][2][n][m]=0.0;
		}
		lhs[i][1][n][n]=1.0;
	}
}
