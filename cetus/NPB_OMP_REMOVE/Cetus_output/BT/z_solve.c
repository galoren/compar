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
#include "work_lhs.h"
#include "timers.h"
/* --------------------------------------------------------------------- */
/* Performs line solves in Z direction by first factoring */
/* the block-tridiagonal matrix into an upper triangular matrix,  */
/* and then performing back substitution to solve for the unknow */
/* vectors of each line.   */
/*  */
/* Make sure we treat elements zero to cell_size in the direction */
/* of the sweep. */
/* --------------------------------------------------------------------- */
void z_solve()
{
	int i, j, k, m, n, ksize;
	/* --------------------------------------------------------------------- */
	/* --------------------------------------------------------------------- */
	if (timeron)
	{
		timer_start(8);
	}
	/* --------------------------------------------------------------------- */
	/* --------------------------------------------------------------------- */
	/* --------------------------------------------------------------------- */
	/* This function computes the left hand side for the three z-factors    */
	/* --------------------------------------------------------------------- */
	ksize=(grid_points[2]-1);
	/* --------------------------------------------------------------------- */
	/* Compute the indices for storing the block-diagonal matrix; */
	/* determine c (labeled f) and s jacobians */
	/* --------------------------------------------------------------------- */
	#pragma cetus private(i, j, k, m, n) 
	#pragma loop name z_solve#0 
	for (j=1; j<=(grid_points[1]-2); j ++ )
	{
		#pragma cetus private(i, k, m, n) 
		#pragma loop name z_solve#0#0 
		for (i=1; i<=(grid_points[0]-2); i ++ )
		{
			#pragma cetus private(k) 
			#pragma cetus lastprivate(tmp1, tmp2, tmp3) 
			#pragma loop name z_solve#0#0#0 
			#pragma cetus parallel 
			#pragma omp parallel for if((10000<(56L+(55L*ksize)))) private(k) lastprivate(tmp1, tmp2, tmp3)
			for (k=0; k<=ksize; k ++ )
			{
				tmp1=(1.0/u[k][j][i][0]);
				tmp2=(tmp1*tmp1);
				tmp3=(tmp1*tmp2);
				fjac[k][0][0]=0.0;
				fjac[k][1][0]=0.0;
				fjac[k][2][0]=0.0;
				fjac[k][3][0]=1.0;
				fjac[k][4][0]=0.0;
				fjac[k][0][1]=(( - (u[k][j][i][1]*u[k][j][i][3]))*tmp2);
				fjac[k][1][1]=(u[k][j][i][3]*tmp1);
				fjac[k][2][1]=0.0;
				fjac[k][3][1]=(u[k][j][i][1]*tmp1);
				fjac[k][4][1]=0.0;
				fjac[k][0][2]=(( - (u[k][j][i][2]*u[k][j][i][3]))*tmp2);
				fjac[k][1][2]=0.0;
				fjac[k][2][2]=(u[k][j][i][3]*tmp1);
				fjac[k][3][2]=(u[k][j][i][2]*tmp1);
				fjac[k][4][2]=0.0;
				fjac[k][0][3]=(( - ((u[k][j][i][3]*u[k][j][i][3])*tmp2))+(c2*qs[k][j][i]));
				fjac[k][1][3]=((( - c2)*u[k][j][i][1])*tmp1);
				fjac[k][2][3]=((( - c2)*u[k][j][i][2])*tmp1);
				fjac[k][3][3]=(((2.0-c2)*u[k][j][i][3])*tmp1);
				fjac[k][4][3]=c2;
				fjac[k][0][4]=(((((c2*2.0)*square[k][j][i])-(c1*u[k][j][i][4]))*u[k][j][i][3])*tmp2);
				fjac[k][1][4]=((( - c2)*(u[k][j][i][1]*u[k][j][i][3]))*tmp2);
				fjac[k][2][4]=((( - c2)*(u[k][j][i][2]*u[k][j][i][3]))*tmp2);
				fjac[k][3][4]=((c1*(u[k][j][i][4]*tmp1))-(c2*(qs[k][j][i]+((u[k][j][i][3]*u[k][j][i][3])*tmp2))));
				fjac[k][4][4]=((c1*u[k][j][i][3])*tmp1);
				njac[k][0][0]=0.0;
				njac[k][1][0]=0.0;
				njac[k][2][0]=0.0;
				njac[k][3][0]=0.0;
				njac[k][4][0]=0.0;
				njac[k][0][1]=((( - c3c4)*tmp2)*u[k][j][i][1]);
				njac[k][1][1]=(c3c4*tmp1);
				njac[k][2][1]=0.0;
				njac[k][3][1]=0.0;
				njac[k][4][1]=0.0;
				njac[k][0][2]=((( - c3c4)*tmp2)*u[k][j][i][2]);
				njac[k][1][2]=0.0;
				njac[k][2][2]=(c3c4*tmp1);
				njac[k][3][2]=0.0;
				njac[k][4][2]=0.0;
				njac[k][0][3]=(((( - con43)*c3c4)*tmp2)*u[k][j][i][3]);
				njac[k][1][3]=0.0;
				njac[k][2][3]=0.0;
				njac[k][3][3]=(((con43*c3)*c4)*tmp1);
				njac[k][4][3]=0.0;
				njac[k][0][4]=(((((( - (c3c4-c1345))*tmp3)*(u[k][j][i][1]*u[k][j][i][1]))-(((c3c4-c1345)*tmp3)*(u[k][j][i][2]*u[k][j][i][2])))-((((con43*c3c4)-c1345)*tmp3)*(u[k][j][i][3]*u[k][j][i][3])))-((c1345*tmp2)*u[k][j][i][4]));
				njac[k][1][4]=(((c3c4-c1345)*tmp2)*u[k][j][i][1]);
				njac[k][2][4]=(((c3c4-c1345)*tmp2)*u[k][j][i][2]);
				njac[k][3][4]=((((con43*c3c4)-c1345)*tmp2)*u[k][j][i][3]);
				njac[k][4][4]=(c1345*tmp1);
			}
			/* --------------------------------------------------------------------- */
			/* now jacobians set, so form left hand side in z direction */
			/* --------------------------------------------------------------------- */
			lhsinit(lhs, ksize);
			#pragma cetus private(k) 
			#pragma cetus lastprivate(tmp1, tmp2) 
			#pragma loop name z_solve#0#0#1 
			#pragma cetus parallel 
			#pragma omp parallel for if((10000<(-78L+(79L*ksize)))) private(k) lastprivate(tmp1, tmp2)
			for (k=1; k<=(ksize-1); k ++ )
			{
				tmp1=(dt*tz1);
				tmp2=(dt*tz2);
				lhs[k][0][0][0]=(((( - tmp2)*fjac[k-1][0][0])-(tmp1*njac[k-1][0][0]))-(tmp1*dz1));
				lhs[k][0][1][0]=((( - tmp2)*fjac[k-1][1][0])-(tmp1*njac[k-1][1][0]));
				lhs[k][0][2][0]=((( - tmp2)*fjac[k-1][2][0])-(tmp1*njac[k-1][2][0]));
				lhs[k][0][3][0]=((( - tmp2)*fjac[k-1][3][0])-(tmp1*njac[k-1][3][0]));
				lhs[k][0][4][0]=((( - tmp2)*fjac[k-1][4][0])-(tmp1*njac[k-1][4][0]));
				lhs[k][0][0][1]=((( - tmp2)*fjac[k-1][0][1])-(tmp1*njac[k-1][0][1]));
				lhs[k][0][1][1]=(((( - tmp2)*fjac[k-1][1][1])-(tmp1*njac[k-1][1][1]))-(tmp1*dz2));
				lhs[k][0][2][1]=((( - tmp2)*fjac[k-1][2][1])-(tmp1*njac[k-1][2][1]));
				lhs[k][0][3][1]=((( - tmp2)*fjac[k-1][3][1])-(tmp1*njac[k-1][3][1]));
				lhs[k][0][4][1]=((( - tmp2)*fjac[k-1][4][1])-(tmp1*njac[k-1][4][1]));
				lhs[k][0][0][2]=((( - tmp2)*fjac[k-1][0][2])-(tmp1*njac[k-1][0][2]));
				lhs[k][0][1][2]=((( - tmp2)*fjac[k-1][1][2])-(tmp1*njac[k-1][1][2]));
				lhs[k][0][2][2]=(((( - tmp2)*fjac[k-1][2][2])-(tmp1*njac[k-1][2][2]))-(tmp1*dz3));
				lhs[k][0][3][2]=((( - tmp2)*fjac[k-1][3][2])-(tmp1*njac[k-1][3][2]));
				lhs[k][0][4][2]=((( - tmp2)*fjac[k-1][4][2])-(tmp1*njac[k-1][4][2]));
				lhs[k][0][0][3]=((( - tmp2)*fjac[k-1][0][3])-(tmp1*njac[k-1][0][3]));
				lhs[k][0][1][3]=((( - tmp2)*fjac[k-1][1][3])-(tmp1*njac[k-1][1][3]));
				lhs[k][0][2][3]=((( - tmp2)*fjac[k-1][2][3])-(tmp1*njac[k-1][2][3]));
				lhs[k][0][3][3]=(((( - tmp2)*fjac[k-1][3][3])-(tmp1*njac[k-1][3][3]))-(tmp1*dz4));
				lhs[k][0][4][3]=((( - tmp2)*fjac[k-1][4][3])-(tmp1*njac[k-1][4][3]));
				lhs[k][0][0][4]=((( - tmp2)*fjac[k-1][0][4])-(tmp1*njac[k-1][0][4]));
				lhs[k][0][1][4]=((( - tmp2)*fjac[k-1][1][4])-(tmp1*njac[k-1][1][4]));
				lhs[k][0][2][4]=((( - tmp2)*fjac[k-1][2][4])-(tmp1*njac[k-1][2][4]));
				lhs[k][0][3][4]=((( - tmp2)*fjac[k-1][3][4])-(tmp1*njac[k-1][3][4]));
				lhs[k][0][4][4]=(((( - tmp2)*fjac[k-1][4][4])-(tmp1*njac[k-1][4][4]))-(tmp1*dz5));
				lhs[k][1][0][0]=((1.0+((tmp1*2.0)*njac[k][0][0]))+((tmp1*2.0)*dz1));
				lhs[k][1][1][0]=((tmp1*2.0)*njac[k][1][0]);
				lhs[k][1][2][0]=((tmp1*2.0)*njac[k][2][0]);
				lhs[k][1][3][0]=((tmp1*2.0)*njac[k][3][0]);
				lhs[k][1][4][0]=((tmp1*2.0)*njac[k][4][0]);
				lhs[k][1][0][1]=((tmp1*2.0)*njac[k][0][1]);
				lhs[k][1][1][1]=((1.0+((tmp1*2.0)*njac[k][1][1]))+((tmp1*2.0)*dz2));
				lhs[k][1][2][1]=((tmp1*2.0)*njac[k][2][1]);
				lhs[k][1][3][1]=((tmp1*2.0)*njac[k][3][1]);
				lhs[k][1][4][1]=((tmp1*2.0)*njac[k][4][1]);
				lhs[k][1][0][2]=((tmp1*2.0)*njac[k][0][2]);
				lhs[k][1][1][2]=((tmp1*2.0)*njac[k][1][2]);
				lhs[k][1][2][2]=((1.0+((tmp1*2.0)*njac[k][2][2]))+((tmp1*2.0)*dz3));
				lhs[k][1][3][2]=((tmp1*2.0)*njac[k][3][2]);
				lhs[k][1][4][2]=((tmp1*2.0)*njac[k][4][2]);
				lhs[k][1][0][3]=((tmp1*2.0)*njac[k][0][3]);
				lhs[k][1][1][3]=((tmp1*2.0)*njac[k][1][3]);
				lhs[k][1][2][3]=((tmp1*2.0)*njac[k][2][3]);
				lhs[k][1][3][3]=((1.0+((tmp1*2.0)*njac[k][3][3]))+((tmp1*2.0)*dz4));
				lhs[k][1][4][3]=((tmp1*2.0)*njac[k][4][3]);
				lhs[k][1][0][4]=((tmp1*2.0)*njac[k][0][4]);
				lhs[k][1][1][4]=((tmp1*2.0)*njac[k][1][4]);
				lhs[k][1][2][4]=((tmp1*2.0)*njac[k][2][4]);
				lhs[k][1][3][4]=((tmp1*2.0)*njac[k][3][4]);
				lhs[k][1][4][4]=((1.0+((tmp1*2.0)*njac[k][4][4]))+((tmp1*2.0)*dz5));
				lhs[k][2][0][0]=(((tmp2*fjac[k+1][0][0])-(tmp1*njac[k+1][0][0]))-(tmp1*dz1));
				lhs[k][2][1][0]=((tmp2*fjac[k+1][1][0])-(tmp1*njac[k+1][1][0]));
				lhs[k][2][2][0]=((tmp2*fjac[k+1][2][0])-(tmp1*njac[k+1][2][0]));
				lhs[k][2][3][0]=((tmp2*fjac[k+1][3][0])-(tmp1*njac[k+1][3][0]));
				lhs[k][2][4][0]=((tmp2*fjac[k+1][4][0])-(tmp1*njac[k+1][4][0]));
				lhs[k][2][0][1]=((tmp2*fjac[k+1][0][1])-(tmp1*njac[k+1][0][1]));
				lhs[k][2][1][1]=(((tmp2*fjac[k+1][1][1])-(tmp1*njac[k+1][1][1]))-(tmp1*dz2));
				lhs[k][2][2][1]=((tmp2*fjac[k+1][2][1])-(tmp1*njac[k+1][2][1]));
				lhs[k][2][3][1]=((tmp2*fjac[k+1][3][1])-(tmp1*njac[k+1][3][1]));
				lhs[k][2][4][1]=((tmp2*fjac[k+1][4][1])-(tmp1*njac[k+1][4][1]));
				lhs[k][2][0][2]=((tmp2*fjac[k+1][0][2])-(tmp1*njac[k+1][0][2]));
				lhs[k][2][1][2]=((tmp2*fjac[k+1][1][2])-(tmp1*njac[k+1][1][2]));
				lhs[k][2][2][2]=(((tmp2*fjac[k+1][2][2])-(tmp1*njac[k+1][2][2]))-(tmp1*dz3));
				lhs[k][2][3][2]=((tmp2*fjac[k+1][3][2])-(tmp1*njac[k+1][3][2]));
				lhs[k][2][4][2]=((tmp2*fjac[k+1][4][2])-(tmp1*njac[k+1][4][2]));
				lhs[k][2][0][3]=((tmp2*fjac[k+1][0][3])-(tmp1*njac[k+1][0][3]));
				lhs[k][2][1][3]=((tmp2*fjac[k+1][1][3])-(tmp1*njac[k+1][1][3]));
				lhs[k][2][2][3]=((tmp2*fjac[k+1][2][3])-(tmp1*njac[k+1][2][3]));
				lhs[k][2][3][3]=(((tmp2*fjac[k+1][3][3])-(tmp1*njac[k+1][3][3]))-(tmp1*dz4));
				lhs[k][2][4][3]=((tmp2*fjac[k+1][4][3])-(tmp1*njac[k+1][4][3]));
				lhs[k][2][0][4]=((tmp2*fjac[k+1][0][4])-(tmp1*njac[k+1][0][4]));
				lhs[k][2][1][4]=((tmp2*fjac[k+1][1][4])-(tmp1*njac[k+1][1][4]));
				lhs[k][2][2][4]=((tmp2*fjac[k+1][2][4])-(tmp1*njac[k+1][2][4]));
				lhs[k][2][3][4]=((tmp2*fjac[k+1][3][4])-(tmp1*njac[k+1][3][4]));
				lhs[k][2][4][4]=(((tmp2*fjac[k+1][4][4])-(tmp1*njac[k+1][4][4]))-(tmp1*dz5));
			}
			/* --------------------------------------------------------------------- */
			/* --------------------------------------------------------------------- */
			/* --------------------------------------------------------------------- */
			/* performs guaussian elimination on this cell. */
			/*  */
			/* assumes that unpacking routines for non-first cells  */
			/* preload C' and rhs' from previous cell. */
			/*  */
			/* assumed send happens outside this routine, but that */
			/* c'(KMAX) and rhs'(KMAX) will be sent to next cell. */
			/* --------------------------------------------------------------------- */
			/* --------------------------------------------------------------------- */
			/* outer most do loops - sweeping in i direction */
			/* --------------------------------------------------------------------- */
			/* --------------------------------------------------------------------- */
			/* multiply c[0][j][i] by b_inverse and copy back to c */
			/* multiply rhs(0) by b_inverse(0) and copy to rhs */
			/* --------------------------------------------------------------------- */
			binvcrhs(lhs[0][1], lhs[0][2], rhs[0][j][i]);
			/* --------------------------------------------------------------------- */
			/* begin inner most do loop */
			/* do all the elements of the cell unless last  */
			/* --------------------------------------------------------------------- */
			#pragma cetus private(k) 
			#pragma loop name z_solve#0#0#2 
			for (k=1; k<=(ksize-1); k ++ )
			{
				/* ------------------------------------------------------------------- */
				/* subtract Alhs_vector(k-1) from lhs_vector(k) */
				/*  */
				/* rhs(k) = rhs(k) - Arhs(k-1) */
				/* ------------------------------------------------------------------- */
				matvec_sub(lhs[k][0], rhs[k-1][j][i], rhs[k][j][i]);
				/* ------------------------------------------------------------------- */
				/* B(k) = B(k) - C(k-1)A(k) */
				/* matmul_sub(AA,i,j,k,c,CC,i,j,k-1,c,BB,i,j,k) */
				/* ------------------------------------------------------------------- */
				matmul_sub(lhs[k][0], lhs[k-1][2], lhs[k][1]);
				/* ------------------------------------------------------------------- */
				/* multiply c[k][j][i] by b_inverse and copy back to c */
				/* multiply rhs[0][j][i] by b_inverse[0][j][i] and copy to rhs */
				/* ------------------------------------------------------------------- */
				binvcrhs(lhs[k][1], lhs[k][2], rhs[k][j][i]);
			}
			/* --------------------------------------------------------------------- */
			/* Now finish up special cases for last cell */
			/* --------------------------------------------------------------------- */
			/* --------------------------------------------------------------------- */
			/* rhs(ksize) = rhs(ksize) - Arhs(ksize-1) */
			/* --------------------------------------------------------------------- */
			matvec_sub(lhs[ksize][0], rhs[ksize-1][j][i], rhs[ksize][j][i]);
			/* --------------------------------------------------------------------- */
			/* B(ksize) = B(ksize) - C(ksize-1)A(ksize) */
			/* matmul_sub(AA,i,j,ksize,c, */
			/* $              CC,i,j,ksize-1,c,BB,i,j,ksize) */
			/* --------------------------------------------------------------------- */
			matmul_sub(lhs[ksize][0], lhs[ksize-1][2], lhs[ksize][1]);
			/* --------------------------------------------------------------------- */
			/* multiply rhs(ksize) by b_inverse(ksize) and copy to rhs */
			/* --------------------------------------------------------------------- */
			binvrhs(lhs[ksize][1], rhs[ksize][j][i]);
			/* --------------------------------------------------------------------- */
			/* --------------------------------------------------------------------- */
			/* --------------------------------------------------------------------- */
			/* back solve: if last cell, then generate U(ksize)=rhs(ksize) */
			/* else assume U(ksize) is loaded in un pack backsub_info */
			/* so just use it */
			/* after u(kstart) will be sent to next cell */
			/* --------------------------------------------------------------------- */
			#pragma cetus private(k, m, n) 
			#pragma loop name z_solve#0#0#3 
			for (k=(ksize-1); k>=0; k -- )
			{
				#pragma cetus private(m, n) 
				#pragma loop name z_solve#0#0#3#0 
				#pragma cetus parallel 
				/*
				Disabled due to low profitability: #pragma omp parallel for private(m, n)
				*/
				for (m=0; m<5; m ++ )
				{
					#pragma cetus private(n) 
					#pragma loop name z_solve#0#0#3#0#0 
					for (n=0; n<5; n ++ )
					{
						rhs[k][j][i][m]=(rhs[k][j][i][m]-(lhs[k][2][n][m]*rhs[k+1][j][i][n]));
					}
				}
			}
		}
	}
	if (timeron)
	{
		timer_stop(8);
	}
}
