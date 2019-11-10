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
/* this function performs the solution of the approximate factorization */
/* step in the x-direction for all five matrix components */
/* simultaneously. The Thomas algorithm is employed to solve the */
/* systems for the x-lines. Boundary conditions are non-periodic */
/* --------------------------------------------------------------------- */
void x_solve()
{
	int i, j, k, i1, i2, m;
	double ru1, fac1, fac2;
	if (timeron)
	{
		timer_start(6);
	}
	#pragma cetus private(fac1, fac2, i, i1, i2, j, k, m, ru1) 
	#pragma loop name x_solve#0 
	for (k=1; k<=nz2; k ++ )
	{
		lhsinit(nx2+1, ny2);
		/* --------------------------------------------------------------------- */
		/* Computes the left hand side for the three x-factors   */
		/* --------------------------------------------------------------------- */
		/* --------------------------------------------------------------------- */
		/* first fill the lhs for the u-eigenvalue                    */
		/* --------------------------------------------------------------------- */
		#pragma cetus private(i, j, ru1) 
		#pragma loop name x_solve#0#0 
		for (j=1; j<=ny2; j ++ )
		{
			#pragma cetus private(i, ru1) 
			#pragma loop name x_solve#0#0#0 
			#pragma cetus parallel 
			#pragma omp parallel for if((10000<(1L+(5L*grid_points[0L])))) private(i, ru1)
			for (i=0; i<=(grid_points[0]-1); i ++ )
			{
				ru1=(c3c4*rho_i[k][j][i]);
				cv[i]=us[k][j][i];
				rhon[i]=(((((dx2+(con43*ru1))>(dx5+(c1c5*ru1))) ? (dx2+(con43*ru1)) : (dx5+(c1c5*ru1)))>(((dxmax+ru1)>dx1) ? (dxmax+ru1) : dx1)) ? (((dx2+(con43*ru1))>(dx5+(c1c5*ru1))) ? (dx2+(con43*ru1)) : (dx5+(c1c5*ru1))) : (((dxmax+ru1)>dx1) ? (dxmax+ru1) : dx1));
			}
			#pragma cetus private(i) 
			#pragma loop name x_solve#0#0#1 
			#pragma cetus parallel 
			#pragma omp parallel for if((10000<(1L+(7L*nx2)))) private(i)
			for (i=1; i<=nx2; i ++ )
			{
				lhs[j][i][0]=0.0;
				lhs[j][i][1]=((( - dttx2)*cv[i-1])-(dttx1*rhon[i-1]));
				lhs[j][i][2]=(1.0+(c2dttx1*rhon[i]));
				lhs[j][i][3]=((dttx2*cv[i+1])-(dttx1*rhon[i+1]));
				lhs[j][i][4]=0.0;
			}
		}
		/* --------------------------------------------------------------------- */
		/* add fourth order dissipation                              */
		/* --------------------------------------------------------------------- */
		#pragma cetus private(i, j) 
		#pragma loop name x_solve#0#1 
		#pragma cetus parallel 
		#pragma omp parallel for if((10000<(1L+(10L*ny2)))) private(i, j)
		for (j=1; j<=ny2; j ++ )
		{
			i=1;
			lhs[j][i][2]=(lhs[j][i][2]+comz5);
			lhs[j][i][3]=(lhs[j][i][3]-comz4);
			lhs[j][i][4]=(lhs[j][i][4]+comz1);
			lhs[j][i+1][1]=(lhs[j][i+1][1]-comz4);
			lhs[j][i+1][2]=(lhs[j][i+1][2]+comz6);
			lhs[j][i+1][3]=(lhs[j][i+1][3]-comz4);
			lhs[j][i+1][4]=(lhs[j][i+1][4]+comz1);
		}
		#pragma cetus private(i, j) 
		#pragma loop name x_solve#0#2 
		#pragma cetus parallel 
		#pragma omp parallel for if((10000<((1L+(-39L*ny2))+((7L*ny2)*grid_points[0L])))) private(i, j)
		for (j=1; j<=ny2; j ++ )
		{
			#pragma cetus private(i) 
			#pragma loop name x_solve#0#2#0 
			for (i=3; i<=(grid_points[0]-4); i ++ )
			{
				lhs[j][i][0]=(lhs[j][i][0]+comz1);
				lhs[j][i][1]=(lhs[j][i][1]-comz4);
				lhs[j][i][2]=(lhs[j][i][2]+comz6);
				lhs[j][i][3]=(lhs[j][i][3]-comz4);
				lhs[j][i][4]=(lhs[j][i][4]+comz1);
			}
		}
		#pragma cetus private(i, j) 
		#pragma loop name x_solve#0#3 
		#pragma cetus parallel 
		#pragma omp parallel for if((10000<(1L+(10L*ny2)))) private(i, j)
		for (j=1; j<=ny2; j ++ )
		{
			i=(grid_points[0]-3);
			lhs[j][i][0]=(lhs[j][i][0]+comz1);
			lhs[j][i][1]=(lhs[j][i][1]-comz4);
			lhs[j][i][2]=(lhs[j][i][2]+comz6);
			lhs[j][i][3]=(lhs[j][i][3]-comz4);
			lhs[j][i+1][0]=(lhs[j][i+1][0]+comz1);
			lhs[j][i+1][1]=(lhs[j][i+1][1]-comz4);
			lhs[j][i+1][2]=(lhs[j][i+1][2]+comz5);
		}
		/* --------------------------------------------------------------------- */
		/* subsequently, fill the other factors (u+c), (u-c) by adding to  */
		/* the first   */
		/* --------------------------------------------------------------------- */
		#pragma cetus private(i, j) 
		#pragma loop name x_solve#0#4 
		#pragma cetus parallel 
		#pragma omp parallel for if((10000<((1L+(3L*ny2))+((12L*nx2)*ny2)))) private(i, j)
		for (j=1; j<=ny2; j ++ )
		{
			#pragma cetus private(i) 
			#pragma loop name x_solve#0#4#0 
			for (i=1; i<=nx2; i ++ )
			{
				lhsp[j][i][0]=lhs[j][i][0];
				lhsp[j][i][1]=(lhs[j][i][1]-(dttx2*speed[k][j][i-1]));
				lhsp[j][i][2]=lhs[j][i][2];
				lhsp[j][i][3]=(lhs[j][i][3]+(dttx2*speed[k][j][i+1]));
				lhsp[j][i][4]=lhs[j][i][4];
				lhsm[j][i][0]=lhs[j][i][0];
				lhsm[j][i][1]=(lhs[j][i][1]+(dttx2*speed[k][j][i-1]));
				lhsm[j][i][2]=lhs[j][i][2];
				lhsm[j][i][3]=(lhs[j][i][3]-(dttx2*speed[k][j][i+1]));
				lhsm[j][i][4]=lhs[j][i][4];
			}
		}
		/* --------------------------------------------------------------------- */
		/* FORWARD ELIMINATION   */
		/* --------------------------------------------------------------------- */
		/* --------------------------------------------------------------------- */
		/* perform the Thomas algorithm; first, FORWARD ELIMINATION      */
		/* --------------------------------------------------------------------- */
		#pragma cetus private(fac1, i, i1, i2, j, m) 
		#pragma loop name x_solve#0#5 
		#pragma cetus parallel 
		#pragma omp parallel for if((10000<((1L+(-79L*ny2))+((41L*ny2)*grid_points[0L])))) private(fac1, i, i1, i2, j, m)
		for (j=1; j<=ny2; j ++ )
		{
			#pragma cetus private(fac1, i, i1, i2, m) 
			#pragma loop name x_solve#0#5#0 
			for (i=0; i<=(grid_points[0]-3); i ++ )
			{
				i1=(i+1);
				i2=(i+2);
				fac1=(1.0/lhs[j][i][2]);
				lhs[j][i][3]=(fac1*lhs[j][i][3]);
				lhs[j][i][4]=(fac1*lhs[j][i][4]);
				#pragma cetus private(m) 
				#pragma loop name x_solve#0#5#0#0 
				for (m=0; m<3; m ++ )
				{
					rhs[k][j][i][m]=(fac1*rhs[k][j][i][m]);
				}
				lhs[j][i1][2]=(lhs[j][i1][2]-(lhs[j][i1][1]*lhs[j][i][3]));
				lhs[j][i1][3]=(lhs[j][i1][3]-(lhs[j][i1][1]*lhs[j][i][4]));
				#pragma cetus private(m) 
				#pragma loop name x_solve#0#5#0#1 
				for (m=0; m<3; m ++ )
				{
					rhs[k][j][i1][m]=(rhs[k][j][i1][m]-(lhs[j][i1][1]*rhs[k][j][i][m]));
				}
				lhs[j][i2][1]=(lhs[j][i2][1]-(lhs[j][i2][0]*lhs[j][i][3]));
				lhs[j][i2][2]=(lhs[j][i2][2]-(lhs[j][i2][0]*lhs[j][i][4]));
				#pragma cetus private(m) 
				#pragma loop name x_solve#0#5#0#2 
				for (m=0; m<3; m ++ )
				{
					rhs[k][j][i2][m]=(rhs[k][j][i2][m]-(lhs[j][i2][0]*rhs[k][j][i][m]));
				}
			}
		}
		/* --------------------------------------------------------------------- */
		/* The last two rows in this grid block are a bit different,  */
		/* since they for (not have two more rows available for the */
		/* elimination of off-diagonal entries */
		/* --------------------------------------------------------------------- */
		#pragma cetus private(fac1, fac2, i, i1, j, m) 
		#pragma loop name x_solve#0#6 
		#pragma cetus parallel 
		#pragma omp parallel for if((10000<(1L+(40L*ny2)))) private(fac1, fac2, i, i1, j, m)
		for (j=1; j<=ny2; j ++ )
		{
			i=(grid_points[0]-2);
			i1=(grid_points[0]-1);
			fac1=(1.0/lhs[j][i][2]);
			lhs[j][i][3]=(fac1*lhs[j][i][3]);
			lhs[j][i][4]=(fac1*lhs[j][i][4]);
			#pragma cetus private(m) 
			#pragma loop name x_solve#0#6#0 
			for (m=0; m<3; m ++ )
			{
				rhs[k][j][i][m]=(fac1*rhs[k][j][i][m]);
			}
			lhs[j][i1][2]=(lhs[j][i1][2]-(lhs[j][i1][1]*lhs[j][i][3]));
			lhs[j][i1][3]=(lhs[j][i1][3]-(lhs[j][i1][1]*lhs[j][i][4]));
			#pragma cetus private(m) 
			#pragma loop name x_solve#0#6#1 
			for (m=0; m<3; m ++ )
			{
				rhs[k][j][i1][m]=(rhs[k][j][i1][m]-(lhs[j][i1][1]*rhs[k][j][i][m]));
			}
			/* --------------------------------------------------------------------- */
			/* scale the last row immediately  */
			/* --------------------------------------------------------------------- */
			fac2=(1.0/lhs[j][i1][2]);
			#pragma cetus private(m) 
			#pragma loop name x_solve#0#6#2 
			for (m=0; m<3; m ++ )
			{
				rhs[k][j][i1][m]=(fac2*rhs[k][j][i1][m]);
			}
		}
		/* --------------------------------------------------------------------- */
		/* for (the u+c and the u-c factors                  */
		/* --------------------------------------------------------------------- */
		#pragma cetus private(fac1, i, i1, i2, j, m) 
		#pragma loop name x_solve#0#7 
		#pragma cetus parallel 
		#pragma omp parallel for if((10000<((1L+(-49L*ny2))+((26L*ny2)*grid_points[0L])))) private(fac1, i, i1, i2, j, m)
		for (j=1; j<=ny2; j ++ )
		{
			#pragma cetus private(fac1, i, i1, i2, m) 
			#pragma loop name x_solve#0#7#0 
			for (i=0; i<=(grid_points[0]-3); i ++ )
			{
				i1=(i+1);
				i2=(i+2);
				m=3;
				fac1=(1.0/lhsp[j][i][2]);
				lhsp[j][i][3]=(fac1*lhsp[j][i][3]);
				lhsp[j][i][4]=(fac1*lhsp[j][i][4]);
				rhs[k][j][i][m]=(fac1*rhs[k][j][i][m]);
				lhsp[j][i1][2]=(lhsp[j][i1][2]-(lhsp[j][i1][1]*lhsp[j][i][3]));
				lhsp[j][i1][3]=(lhsp[j][i1][3]-(lhsp[j][i1][1]*lhsp[j][i][4]));
				rhs[k][j][i1][m]=(rhs[k][j][i1][m]-(lhsp[j][i1][1]*rhs[k][j][i][m]));
				lhsp[j][i2][1]=(lhsp[j][i2][1]-(lhsp[j][i2][0]*lhsp[j][i][3]));
				lhsp[j][i2][2]=(lhsp[j][i2][2]-(lhsp[j][i2][0]*lhsp[j][i][4]));
				rhs[k][j][i2][m]=(rhs[k][j][i2][m]-(lhsp[j][i2][0]*rhs[k][j][i][m]));
				m=4;
				fac1=(1.0/lhsm[j][i][2]);
				lhsm[j][i][3]=(fac1*lhsm[j][i][3]);
				lhsm[j][i][4]=(fac1*lhsm[j][i][4]);
				rhs[k][j][i][m]=(fac1*rhs[k][j][i][m]);
				lhsm[j][i1][2]=(lhsm[j][i1][2]-(lhsm[j][i1][1]*lhsm[j][i][3]));
				lhsm[j][i1][3]=(lhsm[j][i1][3]-(lhsm[j][i1][1]*lhsm[j][i][4]));
				rhs[k][j][i1][m]=(rhs[k][j][i1][m]-(lhsm[j][i1][1]*rhs[k][j][i][m]));
				lhsm[j][i2][1]=(lhsm[j][i2][1]-(lhsm[j][i2][0]*lhsm[j][i][3]));
				lhsm[j][i2][2]=(lhsm[j][i2][2]-(lhsm[j][i2][0]*lhsm[j][i][4]));
				rhs[k][j][i2][m]=(rhs[k][j][i2][m]-(lhsm[j][i2][0]*rhs[k][j][i][m]));
			}
		}
		/* --------------------------------------------------------------------- */
		/* And again the last two rows separately */
		/* --------------------------------------------------------------------- */
		#pragma cetus private(fac1, i, i1, j, m) 
		#pragma loop name x_solve#0#8 
		#pragma cetus parallel 
		#pragma omp parallel for if((10000<(1L+(22L*ny2)))) private(fac1, i, i1, j, m)
		for (j=1; j<=ny2; j ++ )
		{
			i=(grid_points[0]-2);
			i1=(grid_points[0]-1);
			m=3;
			fac1=(1.0/lhsp[j][i][2]);
			lhsp[j][i][3]=(fac1*lhsp[j][i][3]);
			lhsp[j][i][4]=(fac1*lhsp[j][i][4]);
			rhs[k][j][i][m]=(fac1*rhs[k][j][i][m]);
			lhsp[j][i1][2]=(lhsp[j][i1][2]-(lhsp[j][i1][1]*lhsp[j][i][3]));
			lhsp[j][i1][3]=(lhsp[j][i1][3]-(lhsp[j][i1][1]*lhsp[j][i][4]));
			rhs[k][j][i1][m]=(rhs[k][j][i1][m]-(lhsp[j][i1][1]*rhs[k][j][i][m]));
			m=4;
			fac1=(1.0/lhsm[j][i][2]);
			lhsm[j][i][3]=(fac1*lhsm[j][i][3]);
			lhsm[j][i][4]=(fac1*lhsm[j][i][4]);
			rhs[k][j][i][m]=(fac1*rhs[k][j][i][m]);
			lhsm[j][i1][2]=(lhsm[j][i1][2]-(lhsm[j][i1][1]*lhsm[j][i][3]));
			lhsm[j][i1][3]=(lhsm[j][i1][3]-(lhsm[j][i1][1]*lhsm[j][i][4]));
			rhs[k][j][i1][m]=(rhs[k][j][i1][m]-(lhsm[j][i1][1]*rhs[k][j][i][m]));
			/* --------------------------------------------------------------------- */
			/* Scale the last row immediately */
			/* --------------------------------------------------------------------- */
			rhs[k][j][i1][3]=(rhs[k][j][i1][3]/lhsp[j][i1][2]);
			rhs[k][j][i1][4]=(rhs[k][j][i1][4]/lhsm[j][i1][2]);
		}
		/* --------------------------------------------------------------------- */
		/* BACKSUBSTITUTION  */
		/* --------------------------------------------------------------------- */
		#pragma cetus private(i, i1, j, m) 
		#pragma loop name x_solve#0#9 
		#pragma cetus parallel 
		#pragma omp parallel for if((10000<(1L+(16L*ny2)))) private(i, i1, j, m)
		for (j=1; j<=ny2; j ++ )
		{
			i=(grid_points[0]-2);
			i1=(grid_points[0]-1);
			#pragma cetus private(m) 
			#pragma loop name x_solve#0#9#0 
			for (m=0; m<3; m ++ )
			{
				rhs[k][j][i][m]=(rhs[k][j][i][m]-(lhs[j][i][3]*rhs[k][j][i1][m]));
			}
			rhs[k][j][i][3]=(rhs[k][j][i][3]-(lhsp[j][i][3]*rhs[k][j][i1][3]));
			rhs[k][j][i][4]=(rhs[k][j][i][4]-(lhsm[j][i][3]*rhs[k][j][i1][4]));
		}
		/* --------------------------------------------------------------------- */
		/* The first three factors */
		/* --------------------------------------------------------------------- */
		#pragma cetus private(i, i1, i2, j, m) 
		#pragma loop name x_solve#0#10 
		#pragma cetus parallel 
		#pragma omp parallel for if((10000<((1L+(-61L*ny2))+((16L*ny2)*grid_points[0L])))) private(i, i1, i2, j, m)
		for (j=1; j<=ny2; j ++ )
		{
			#pragma cetus private(i, i1, i2, m) 
			#pragma loop name x_solve#0#10#0 
			for (i=(grid_points[0]-3); i>=0; i -- )
			{
				i1=(i+1);
				i2=(i+2);
				#pragma cetus private(m) 
				#pragma loop name x_solve#0#10#0#0 
				for (m=0; m<3; m ++ )
				{
					rhs[k][j][i][m]=((rhs[k][j][i][m]-(lhs[j][i][3]*rhs[k][j][i1][m]))-(lhs[j][i][4]*rhs[k][j][i2][m]));
				}
				/* ------------------------------------------------------------------- */
				/* And the remaining two */
				/* ------------------------------------------------------------------- */
				rhs[k][j][i][3]=((rhs[k][j][i][3]-(lhsp[j][i][3]*rhs[k][j][i1][3]))-(lhsp[j][i][4]*rhs[k][j][i2][3]));
				rhs[k][j][i][4]=((rhs[k][j][i][4]-(lhsm[j][i][3]*rhs[k][j][i1][4]))-(lhsm[j][i][4]*rhs[k][j][i2][4]));
			}
		}
	}
	if (timeron)
	{
		timer_stop(6);
	}
	/* --------------------------------------------------------------------- */
	/* Do the block-diagonal inversion           */
	/* --------------------------------------------------------------------- */
	ninvr();
}
