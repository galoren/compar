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
/*  */
/* Performs line solves in X direction by first factoring */
/* the block-tridiagonal matrix into an upper triangular matrix,  */
/* and then performing back substitution to solve for the unknow */
/* vectors of each line.   */
/*  */
/* Make sure we treat elements zero to cell_size in the direction */
/* of the sweep. */
/*  */
/* --------------------------------------------------------------------- */
void x_solve()
{
	int i, j, k, m, n, isize;
	/* --------------------------------------------------------------------- */
	/* --------------------------------------------------------------------- */
	if (timeron)
	{
		timer_start(6);
	}
	/* --------------------------------------------------------------------- */
	/* --------------------------------------------------------------------- */
	/* --------------------------------------------------------------------- */
	/* This function computes the left hand side in the xi-direction */
	/* --------------------------------------------------------------------- */
	isize=(grid_points[0]-1);
	/* --------------------------------------------------------------------- */
	/* determine a (labeled f) and n jacobians */
	/* --------------------------------------------------------------------- */
	#pragma cetus private(i, j, k, m, n) 
	#pragma loop name x_solve#0 
	for (k=1; k<=(grid_points[2]-2); k ++ )
	{
		#pragma cetus private(i, j, m, n) 
		#pragma loop name x_solve#0#0 
		for (j=1; j<=(grid_points[1]-2); j ++ )
		{
			#pragma cetus private(i) 
			#pragma cetus lastprivate(tmp1, tmp2, tmp3) 
			#pragma loop name x_solve#0#0#0 
			#pragma cetus parallel 
			#pragma omp parallel for if((10000<(56L+(55L*isize)))) private(i) lastprivate(tmp1, tmp2, tmp3)
			for (i=0; i<=isize; i ++ )
			{
				tmp1=rho_i[k][j][i];
				tmp2=(tmp1*tmp1);
				tmp3=(tmp1*tmp2);
				/* ------------------------------------------------------------------- */
				/*  */
				/* ------------------------------------------------------------------- */
				fjac[i][0][0]=0.0;
				fjac[i][1][0]=1.0;
				fjac[i][2][0]=0.0;
				fjac[i][3][0]=0.0;
				fjac[i][4][0]=0.0;
				fjac[i][0][1]=(( - ((u[k][j][i][1]*tmp2)*u[k][j][i][1]))+(c2*qs[k][j][i]));
				fjac[i][1][1]=((2.0-c2)*(u[k][j][i][1]/u[k][j][i][0]));
				fjac[i][2][1]=(( - c2)*(u[k][j][i][2]*tmp1));
				fjac[i][3][1]=(( - c2)*(u[k][j][i][3]*tmp1));
				fjac[i][4][1]=c2;
				fjac[i][0][2]=(( - (u[k][j][i][1]*u[k][j][i][2]))*tmp2);
				fjac[i][1][2]=(u[k][j][i][2]*tmp1);
				fjac[i][2][2]=(u[k][j][i][1]*tmp1);
				fjac[i][3][2]=0.0;
				fjac[i][4][2]=0.0;
				fjac[i][0][3]=(( - (u[k][j][i][1]*u[k][j][i][3]))*tmp2);
				fjac[i][1][3]=(u[k][j][i][3]*tmp1);
				fjac[i][2][3]=0.0;
				fjac[i][3][3]=(u[k][j][i][1]*tmp1);
				fjac[i][4][3]=0.0;
				fjac[i][0][4]=((((c2*2.0)*square[k][j][i])-(c1*u[k][j][i][4]))*(u[k][j][i][1]*tmp2));
				fjac[i][1][4]=(((c1*u[k][j][i][4])*tmp1)-(c2*(((u[k][j][i][1]*u[k][j][i][1])*tmp2)+qs[k][j][i])));
				fjac[i][2][4]=((( - c2)*(u[k][j][i][2]*u[k][j][i][1]))*tmp2);
				fjac[i][3][4]=((( - c2)*(u[k][j][i][3]*u[k][j][i][1]))*tmp2);
				fjac[i][4][4]=(c1*(u[k][j][i][1]*tmp1));
				njac[i][0][0]=0.0;
				njac[i][1][0]=0.0;
				njac[i][2][0]=0.0;
				njac[i][3][0]=0.0;
				njac[i][4][0]=0.0;
				njac[i][0][1]=(((( - con43)*c3c4)*tmp2)*u[k][j][i][1]);
				njac[i][1][1]=((con43*c3c4)*tmp1);
				njac[i][2][1]=0.0;
				njac[i][3][1]=0.0;
				njac[i][4][1]=0.0;
				njac[i][0][2]=((( - c3c4)*tmp2)*u[k][j][i][2]);
				njac[i][1][2]=0.0;
				njac[i][2][2]=(c3c4*tmp1);
				njac[i][3][2]=0.0;
				njac[i][4][2]=0.0;
				njac[i][0][3]=((( - c3c4)*tmp2)*u[k][j][i][3]);
				njac[i][1][3]=0.0;
				njac[i][2][3]=0.0;
				njac[i][3][3]=(c3c4*tmp1);
				njac[i][4][3]=0.0;
				njac[i][0][4]=(((((( - ((con43*c3c4)-c1345))*tmp3)*(u[k][j][i][1]*u[k][j][i][1]))-(((c3c4-c1345)*tmp3)*(u[k][j][i][2]*u[k][j][i][2])))-(((c3c4-c1345)*tmp3)*(u[k][j][i][3]*u[k][j][i][3])))-((c1345*tmp2)*u[k][j][i][4]));
				njac[i][1][4]=((((con43*c3c4)-c1345)*tmp2)*u[k][j][i][1]);
				njac[i][2][4]=(((c3c4-c1345)*tmp2)*u[k][j][i][2]);
				njac[i][3][4]=(((c3c4-c1345)*tmp2)*u[k][j][i][3]);
				njac[i][4][4]=(c1345*tmp1);
			}
			/* --------------------------------------------------------------------- */
			/* now jacobians set, so form left hand side in x direction */
			/* --------------------------------------------------------------------- */
			lhsinit(lhs, isize);
			#pragma cetus private(i) 
			#pragma cetus lastprivate(tmp1, tmp2) 
			#pragma loop name x_solve#0#0#1 
			#pragma cetus parallel 
			#pragma omp parallel for if((10000<(-78L+(79L*isize)))) private(i) lastprivate(tmp1, tmp2)
			for (i=1; i<=(isize-1); i ++ )
			{
				tmp1=(dt*tx1);
				tmp2=(dt*tx2);
				lhs[i][0][0][0]=(((( - tmp2)*fjac[i-1][0][0])-(tmp1*njac[i-1][0][0]))-(tmp1*dx1));
				lhs[i][0][1][0]=((( - tmp2)*fjac[i-1][1][0])-(tmp1*njac[i-1][1][0]));
				lhs[i][0][2][0]=((( - tmp2)*fjac[i-1][2][0])-(tmp1*njac[i-1][2][0]));
				lhs[i][0][3][0]=((( - tmp2)*fjac[i-1][3][0])-(tmp1*njac[i-1][3][0]));
				lhs[i][0][4][0]=((( - tmp2)*fjac[i-1][4][0])-(tmp1*njac[i-1][4][0]));
				lhs[i][0][0][1]=((( - tmp2)*fjac[i-1][0][1])-(tmp1*njac[i-1][0][1]));
				lhs[i][0][1][1]=(((( - tmp2)*fjac[i-1][1][1])-(tmp1*njac[i-1][1][1]))-(tmp1*dx2));
				lhs[i][0][2][1]=((( - tmp2)*fjac[i-1][2][1])-(tmp1*njac[i-1][2][1]));
				lhs[i][0][3][1]=((( - tmp2)*fjac[i-1][3][1])-(tmp1*njac[i-1][3][1]));
				lhs[i][0][4][1]=((( - tmp2)*fjac[i-1][4][1])-(tmp1*njac[i-1][4][1]));
				lhs[i][0][0][2]=((( - tmp2)*fjac[i-1][0][2])-(tmp1*njac[i-1][0][2]));
				lhs[i][0][1][2]=((( - tmp2)*fjac[i-1][1][2])-(tmp1*njac[i-1][1][2]));
				lhs[i][0][2][2]=(((( - tmp2)*fjac[i-1][2][2])-(tmp1*njac[i-1][2][2]))-(tmp1*dx3));
				lhs[i][0][3][2]=((( - tmp2)*fjac[i-1][3][2])-(tmp1*njac[i-1][3][2]));
				lhs[i][0][4][2]=((( - tmp2)*fjac[i-1][4][2])-(tmp1*njac[i-1][4][2]));
				lhs[i][0][0][3]=((( - tmp2)*fjac[i-1][0][3])-(tmp1*njac[i-1][0][3]));
				lhs[i][0][1][3]=((( - tmp2)*fjac[i-1][1][3])-(tmp1*njac[i-1][1][3]));
				lhs[i][0][2][3]=((( - tmp2)*fjac[i-1][2][3])-(tmp1*njac[i-1][2][3]));
				lhs[i][0][3][3]=(((( - tmp2)*fjac[i-1][3][3])-(tmp1*njac[i-1][3][3]))-(tmp1*dx4));
				lhs[i][0][4][3]=((( - tmp2)*fjac[i-1][4][3])-(tmp1*njac[i-1][4][3]));
				lhs[i][0][0][4]=((( - tmp2)*fjac[i-1][0][4])-(tmp1*njac[i-1][0][4]));
				lhs[i][0][1][4]=((( - tmp2)*fjac[i-1][1][4])-(tmp1*njac[i-1][1][4]));
				lhs[i][0][2][4]=((( - tmp2)*fjac[i-1][2][4])-(tmp1*njac[i-1][2][4]));
				lhs[i][0][3][4]=((( - tmp2)*fjac[i-1][3][4])-(tmp1*njac[i-1][3][4]));
				lhs[i][0][4][4]=(((( - tmp2)*fjac[i-1][4][4])-(tmp1*njac[i-1][4][4]))-(tmp1*dx5));
				lhs[i][1][0][0]=((1.0+((tmp1*2.0)*njac[i][0][0]))+((tmp1*2.0)*dx1));
				lhs[i][1][1][0]=((tmp1*2.0)*njac[i][1][0]);
				lhs[i][1][2][0]=((tmp1*2.0)*njac[i][2][0]);
				lhs[i][1][3][0]=((tmp1*2.0)*njac[i][3][0]);
				lhs[i][1][4][0]=((tmp1*2.0)*njac[i][4][0]);
				lhs[i][1][0][1]=((tmp1*2.0)*njac[i][0][1]);
				lhs[i][1][1][1]=((1.0+((tmp1*2.0)*njac[i][1][1]))+((tmp1*2.0)*dx2));
				lhs[i][1][2][1]=((tmp1*2.0)*njac[i][2][1]);
				lhs[i][1][3][1]=((tmp1*2.0)*njac[i][3][1]);
				lhs[i][1][4][1]=((tmp1*2.0)*njac[i][4][1]);
				lhs[i][1][0][2]=((tmp1*2.0)*njac[i][0][2]);
				lhs[i][1][1][2]=((tmp1*2.0)*njac[i][1][2]);
				lhs[i][1][2][2]=((1.0+((tmp1*2.0)*njac[i][2][2]))+((tmp1*2.0)*dx3));
				lhs[i][1][3][2]=((tmp1*2.0)*njac[i][3][2]);
				lhs[i][1][4][2]=((tmp1*2.0)*njac[i][4][2]);
				lhs[i][1][0][3]=((tmp1*2.0)*njac[i][0][3]);
				lhs[i][1][1][3]=((tmp1*2.0)*njac[i][1][3]);
				lhs[i][1][2][3]=((tmp1*2.0)*njac[i][2][3]);
				lhs[i][1][3][3]=((1.0+((tmp1*2.0)*njac[i][3][3]))+((tmp1*2.0)*dx4));
				lhs[i][1][4][3]=((tmp1*2.0)*njac[i][4][3]);
				lhs[i][1][0][4]=((tmp1*2.0)*njac[i][0][4]);
				lhs[i][1][1][4]=((tmp1*2.0)*njac[i][1][4]);
				lhs[i][1][2][4]=((tmp1*2.0)*njac[i][2][4]);
				lhs[i][1][3][4]=((tmp1*2.0)*njac[i][3][4]);
				lhs[i][1][4][4]=((1.0+((tmp1*2.0)*njac[i][4][4]))+((tmp1*2.0)*dx5));
				lhs[i][2][0][0]=(((tmp2*fjac[i+1][0][0])-(tmp1*njac[i+1][0][0]))-(tmp1*dx1));
				lhs[i][2][1][0]=((tmp2*fjac[i+1][1][0])-(tmp1*njac[i+1][1][0]));
				lhs[i][2][2][0]=((tmp2*fjac[i+1][2][0])-(tmp1*njac[i+1][2][0]));
				lhs[i][2][3][0]=((tmp2*fjac[i+1][3][0])-(tmp1*njac[i+1][3][0]));
				lhs[i][2][4][0]=((tmp2*fjac[i+1][4][0])-(tmp1*njac[i+1][4][0]));
				lhs[i][2][0][1]=((tmp2*fjac[i+1][0][1])-(tmp1*njac[i+1][0][1]));
				lhs[i][2][1][1]=(((tmp2*fjac[i+1][1][1])-(tmp1*njac[i+1][1][1]))-(tmp1*dx2));
				lhs[i][2][2][1]=((tmp2*fjac[i+1][2][1])-(tmp1*njac[i+1][2][1]));
				lhs[i][2][3][1]=((tmp2*fjac[i+1][3][1])-(tmp1*njac[i+1][3][1]));
				lhs[i][2][4][1]=((tmp2*fjac[i+1][4][1])-(tmp1*njac[i+1][4][1]));
				lhs[i][2][0][2]=((tmp2*fjac[i+1][0][2])-(tmp1*njac[i+1][0][2]));
				lhs[i][2][1][2]=((tmp2*fjac[i+1][1][2])-(tmp1*njac[i+1][1][2]));
				lhs[i][2][2][2]=(((tmp2*fjac[i+1][2][2])-(tmp1*njac[i+1][2][2]))-(tmp1*dx3));
				lhs[i][2][3][2]=((tmp2*fjac[i+1][3][2])-(tmp1*njac[i+1][3][2]));
				lhs[i][2][4][2]=((tmp2*fjac[i+1][4][2])-(tmp1*njac[i+1][4][2]));
				lhs[i][2][0][3]=((tmp2*fjac[i+1][0][3])-(tmp1*njac[i+1][0][3]));
				lhs[i][2][1][3]=((tmp2*fjac[i+1][1][3])-(tmp1*njac[i+1][1][3]));
				lhs[i][2][2][3]=((tmp2*fjac[i+1][2][3])-(tmp1*njac[i+1][2][3]));
				lhs[i][2][3][3]=(((tmp2*fjac[i+1][3][3])-(tmp1*njac[i+1][3][3]))-(tmp1*dx4));
				lhs[i][2][4][3]=((tmp2*fjac[i+1][4][3])-(tmp1*njac[i+1][4][3]));
				lhs[i][2][0][4]=((tmp2*fjac[i+1][0][4])-(tmp1*njac[i+1][0][4]));
				lhs[i][2][1][4]=((tmp2*fjac[i+1][1][4])-(tmp1*njac[i+1][1][4]));
				lhs[i][2][2][4]=((tmp2*fjac[i+1][2][4])-(tmp1*njac[i+1][2][4]));
				lhs[i][2][3][4]=((tmp2*fjac[i+1][3][4])-(tmp1*njac[i+1][3][4]));
				lhs[i][2][4][4]=(((tmp2*fjac[i+1][4][4])-(tmp1*njac[i+1][4][4]))-(tmp1*dx5));
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
			/* c'(IMAX) and rhs'(IMAX) will be sent to next cell */
			/* --------------------------------------------------------------------- */
			/* --------------------------------------------------------------------- */
			/* outer most do loops - sweeping in i direction */
			/* --------------------------------------------------------------------- */
			/* --------------------------------------------------------------------- */
			/* multiply c[k][j][0] by b_inverse and copy back to c */
			/* multiply rhs(0) by b_inverse(0) and copy to rhs */
			/* --------------------------------------------------------------------- */
			binvcrhs(lhs[0][1], lhs[0][2], rhs[k][j][0]);
			/* --------------------------------------------------------------------- */
			/* begin inner most do loop */
			/* do all the elements of the cell unless last  */
			/* --------------------------------------------------------------------- */
			#pragma cetus private(i) 
			#pragma loop name x_solve#0#0#2 
			for (i=1; i<=(isize-1); i ++ )
			{
				/* ------------------------------------------------------------------- */
				/* rhs(i) = rhs(i) - Arhs(i-1) */
				/* ------------------------------------------------------------------- */
				matvec_sub(lhs[i][0], rhs[k][j][i-1], rhs[k][j][i]);
				/* ------------------------------------------------------------------- */
				/* B(i) = B(i) - C(i-1)A(i) */
				/* ------------------------------------------------------------------- */
				matmul_sub(lhs[i][0], lhs[i-1][2], lhs[i][1]);
				/* ------------------------------------------------------------------- */
				/* multiply c[k][j][i] by b_inverse and copy back to c */
				/* multiply rhs[k][j][0] by b_inverse[k][j][0] and copy to rhs */
				/* ------------------------------------------------------------------- */
				binvcrhs(lhs[i][1], lhs[i][2], rhs[k][j][i]);
			}
			/* --------------------------------------------------------------------- */
			/* rhs(isize) = rhs(isize) - Arhs(isize-1) */
			/* --------------------------------------------------------------------- */
			matvec_sub(lhs[isize][0], rhs[k][j][isize-1], rhs[k][j][isize]);
			/* --------------------------------------------------------------------- */
			/* B(isize) = B(isize) - C(isize-1)A(isize) */
			/* --------------------------------------------------------------------- */
			matmul_sub(lhs[isize][0], lhs[isize-1][2], lhs[isize][1]);
			/* --------------------------------------------------------------------- */
			/* multiply rhs() by b_inverse() and copy to rhs */
			/* --------------------------------------------------------------------- */
			binvrhs(lhs[isize][1], rhs[k][j][isize]);
			/* --------------------------------------------------------------------- */
			/* back solve: if last cell, then generate U(isize)=rhs(isize) */
			/* else assume U(isize) is loaded in un pack backsub_info */
			/* so just use it */
			/* after u(istart) will be sent to next cell */
			/* --------------------------------------------------------------------- */
			#pragma cetus private(i, m, n) 
			#pragma loop name x_solve#0#0#3 
			for (i=(isize-1); i>=0; i -- )
			{
				#pragma cetus private(m, n) 
				#pragma loop name x_solve#0#0#3#0 
				#pragma cetus parallel 
				/*
				Disabled due to low profitability: #pragma omp parallel for private(m, n)
				*/
				for (m=0; m<5; m ++ )
				{
					#pragma cetus private(n) 
					#pragma loop name x_solve#0#0#3#0#0 
					for (n=0; n<5; n ++ )
					{
						rhs[k][j][i][m]=(rhs[k][j][i][m]-(lhs[i][2][n][m]*rhs[k][j][i+1][n]));
					}
				}
			}
		}
	}
	if (timeron)
	{
		timer_stop(6);
	}
}
