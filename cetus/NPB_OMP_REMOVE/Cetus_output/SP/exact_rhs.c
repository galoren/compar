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
/* compute the right hand side based on exact solution */
/* --------------------------------------------------------------------- */
void exact_rhs()
{
	double dtemp[5], xi, eta, zeta, dtpp;
	int m, i, j, k, ip1, im1, jp1, jm1, km1, kp1;
	/* --------------------------------------------------------------------- */
	/* initialize                                   */
	/* --------------------------------------------------------------------- */
	#pragma cetus private(i, j, k, m) 
	#pragma loop name exact_rhs#0 
	#pragma cetus parallel 
	#pragma omp parallel for if((10000<(((1L+(3L*grid_points[2L]))+((3L*grid_points[1L])*grid_points[2L]))+(((18L*grid_points[0L])*grid_points[1L])*grid_points[2L])))) private(i, j, k, m)
	for (k=0; k<=(grid_points[2]-1); k ++ )
	{
		#pragma cetus private(i, j, m) 
		#pragma loop name exact_rhs#0#0 
		for (j=0; j<=(grid_points[1]-1); j ++ )
		{
			#pragma cetus private(i, m) 
			#pragma loop name exact_rhs#0#0#0 
			for (i=0; i<=(grid_points[0]-1); i ++ )
			{
				#pragma cetus private(m) 
				#pragma loop name exact_rhs#0#0#0#0 
				for (m=0; m<5; m ++ )
				{
					forcing[k][j][i][m]=0.0;
				}
			}
		}
	}
	/* --------------------------------------------------------------------- */
	/* xi-direction flux differences                       */
	/* --------------------------------------------------------------------- */
	#pragma cetus private(dtpp, eta, i, im1, ip1, j, k, m, xi, zeta) 
	#pragma loop name exact_rhs#1 
	for (k=1; k<=(grid_points[2]-2); k ++ )
	{
		zeta=(((double)k)*dnzm1);
		#pragma cetus private(dtpp, eta, i, im1, ip1, j, m, xi) 
		#pragma loop name exact_rhs#1#0 
		for (j=1; j<=(grid_points[1]-2); j ++ )
		{
			eta=(((double)j)*dnym1);
			#pragma cetus private(dtpp, i, m, xi) 
			#pragma loop name exact_rhs#1#0#0 
			for (i=0; i<=(grid_points[0]-1); i ++ )
			{
				xi=(((double)i)*dnxm1);
				exact_solution(xi, eta, zeta, dtemp);
				#pragma cetus private(m) 
				#pragma loop name exact_rhs#1#0#0#0 
				#pragma cetus parallel 
				/*
				Disabled due to low profitability: #pragma omp parallel for private(m)
				*/
				for (m=0; m<5; m ++ )
				{
					ue[i][m]=dtemp[m];
				}
				dtpp=(1.0/dtemp[0]);
				#pragma cetus private(m) 
				#pragma loop name exact_rhs#1#0#0#1 
				#pragma cetus parallel 
				/*
				Disabled due to low profitability: #pragma omp parallel for private(m)
				*/
				for (m=1; m<5; m ++ )
				{
					buf[i][m]=(dtpp*dtemp[m]);
				}
				cuf[i]=(buf[i][1]*buf[i][1]);
				buf[i][0]=((cuf[i]+(buf[i][2]*buf[i][2]))+(buf[i][3]*buf[i][3]));
				q[i]=(0.5*(((buf[i][1]*ue[i][1])+(buf[i][2]*ue[i][2]))+(buf[i][3]*ue[i][3])));
			}
			#pragma cetus private(i, im1, ip1) 
			#pragma loop name exact_rhs#1#0#1 
			#pragma cetus parallel 
			#pragma omp parallel for if((10000<(-17L+(9L*grid_points[0L])))) private(i, im1, ip1)
			for (i=1; i<=(grid_points[0]-2); i ++ )
			{
				im1=(i-1);
				ip1=(i+1);
				forcing[k][j][i][0]=((forcing[k][j][i][0]-(tx2*(ue[ip1][1]-ue[im1][1])))+(dx1tx1*((ue[ip1][0]-(2.0*ue[i][0]))+ue[im1][0])));
				forcing[k][j][i][1]=(((forcing[k][j][i][1]-(tx2*(((ue[ip1][1]*buf[ip1][1])+(c2*(ue[ip1][4]-q[ip1])))-((ue[im1][1]*buf[im1][1])+(c2*(ue[im1][4]-q[im1]))))))+(xxcon1*((buf[ip1][1]-(2.0*buf[i][1]))+buf[im1][1])))+(dx2tx1*((ue[ip1][1]-(2.0*ue[i][1]))+ue[im1][1])));
				forcing[k][j][i][2]=(((forcing[k][j][i][2]-(tx2*((ue[ip1][2]*buf[ip1][1])-(ue[im1][2]*buf[im1][1]))))+(xxcon2*((buf[ip1][2]-(2.0*buf[i][2]))+buf[im1][2])))+(dx3tx1*((ue[ip1][2]-(2.0*ue[i][2]))+ue[im1][2])));
				forcing[k][j][i][3]=(((forcing[k][j][i][3]-(tx2*((ue[ip1][3]*buf[ip1][1])-(ue[im1][3]*buf[im1][1]))))+(xxcon2*((buf[ip1][3]-(2.0*buf[i][3]))+buf[im1][3])))+(dx4tx1*((ue[ip1][3]-(2.0*ue[i][3]))+ue[im1][3])));
				forcing[k][j][i][4]=(((((forcing[k][j][i][4]-(tx2*((buf[ip1][1]*((c1*ue[ip1][4])-(c2*q[ip1])))-(buf[im1][1]*((c1*ue[im1][4])-(c2*q[im1]))))))+((0.5*xxcon3)*((buf[ip1][0]-(2.0*buf[i][0]))+buf[im1][0])))+(xxcon4*((cuf[ip1]-(2.0*cuf[i]))+cuf[im1])))+(xxcon5*((buf[ip1][4]-(2.0*buf[i][4]))+buf[im1][4])))+(dx5tx1*((ue[ip1][4]-(2.0*ue[i][4]))+ue[im1][4])));
			}
			/* --------------------------------------------------------------------- */
			/* Fourth-order dissipation                          */
			/* --------------------------------------------------------------------- */
			#pragma cetus private(i, m) 
			#pragma loop name exact_rhs#1#0#2 
			#pragma cetus parallel 
			/*
			Disabled due to low profitability: #pragma omp parallel for private(i, m)
			*/
			for (m=0; m<5; m ++ )
			{
				i=1;
				forcing[k][j][i][m]=(forcing[k][j][i][m]-(dssp*(((5.0*ue[i][m])-(4.0*ue[i+1][m]))+ue[i+2][m])));
				i=2;
				forcing[k][j][i][m]=(forcing[k][j][i][m]-(dssp*((((( - 4.0)*ue[i-1][m])+(6.0*ue[i][m]))-(4.0*ue[i+1][m]))+ue[i+2][m])));
			}
			#pragma cetus private(i, m) 
			#pragma loop name exact_rhs#1#0#3 
			#pragma cetus parallel 
			#pragma omp parallel for if((10000<(-107L+(18L*grid_points[0L])))) private(i, m)
			for (i=3; i<=(grid_points[0]-4); i ++ )
			{
				#pragma cetus private(m) 
				#pragma loop name exact_rhs#1#0#3#0 
				for (m=0; m<5; m ++ )
				{
					forcing[k][j][i][m]=(forcing[k][j][i][m]-(dssp*((((ue[i-2][m]-(4.0*ue[i-1][m]))+(6.0*ue[i][m]))-(4.0*ue[i+1][m]))+ue[i+2][m])));
				}
			}
			#pragma cetus private(i, m) 
			#pragma loop name exact_rhs#1#0#4 
			#pragma cetus parallel 
			/*
			Disabled due to low profitability: #pragma omp parallel for private(i, m)
			*/
			for (m=0; m<5; m ++ )
			{
				i=(grid_points[0]-3);
				forcing[k][j][i][m]=(forcing[k][j][i][m]-(dssp*(((ue[i-2][m]-(4.0*ue[i-1][m]))+(6.0*ue[i][m]))-(4.0*ue[i+1][m]))));
				i=(grid_points[0]-2);
				forcing[k][j][i][m]=(forcing[k][j][i][m]-(dssp*((ue[i-2][m]-(4.0*ue[i-1][m]))+(5.0*ue[i][m]))));
			}
		}
	}
	/* --------------------------------------------------------------------- */
	/* eta-direction flux differences              */
	/* --------------------------------------------------------------------- */
	#pragma cetus private(dtpp, eta, i, j, jm1, jp1, k, m, xi, zeta) 
	#pragma loop name exact_rhs#2 
	for (k=1; k<=(grid_points[2]-2); k ++ )
	{
		zeta=(((double)k)*dnzm1);
		#pragma cetus private(dtpp, eta, i, j, jm1, jp1, m, xi) 
		#pragma loop name exact_rhs#2#0 
		for (i=1; i<=(grid_points[0]-2); i ++ )
		{
			xi=(((double)i)*dnxm1);
			#pragma cetus private(dtpp, eta, j, m) 
			#pragma loop name exact_rhs#2#0#0 
			for (j=0; j<=(grid_points[1]-1); j ++ )
			{
				eta=(((double)j)*dnym1);
				exact_solution(xi, eta, zeta, dtemp);
				#pragma cetus private(m) 
				#pragma loop name exact_rhs#2#0#0#0 
				#pragma cetus parallel 
				/*
				Disabled due to low profitability: #pragma omp parallel for private(m)
				*/
				for (m=0; m<5; m ++ )
				{
					ue[j][m]=dtemp[m];
				}
				dtpp=(1.0/dtemp[0]);
				#pragma cetus private(m) 
				#pragma loop name exact_rhs#2#0#0#1 
				#pragma cetus parallel 
				/*
				Disabled due to low profitability: #pragma omp parallel for private(m)
				*/
				for (m=1; m<5; m ++ )
				{
					buf[j][m]=(dtpp*dtemp[m]);
				}
				cuf[j]=(buf[j][2]*buf[j][2]);
				buf[j][0]=((cuf[j]+(buf[j][1]*buf[j][1]))+(buf[j][3]*buf[j][3]));
				q[j]=(0.5*(((buf[j][1]*ue[j][1])+(buf[j][2]*ue[j][2]))+(buf[j][3]*ue[j][3])));
			}
			#pragma cetus private(j, jm1, jp1) 
			#pragma loop name exact_rhs#2#0#1 
			#pragma cetus parallel 
			#pragma omp parallel for if((10000<(-17L+(9L*grid_points[1L])))) private(j, jm1, jp1)
			for (j=1; j<=(grid_points[1]-2); j ++ )
			{
				jm1=(j-1);
				jp1=(j+1);
				forcing[k][j][i][0]=((forcing[k][j][i][0]-(ty2*(ue[jp1][2]-ue[jm1][2])))+(dy1ty1*((ue[jp1][0]-(2.0*ue[j][0]))+ue[jm1][0])));
				forcing[k][j][i][1]=(((forcing[k][j][i][1]-(ty2*((ue[jp1][1]*buf[jp1][2])-(ue[jm1][1]*buf[jm1][2]))))+(yycon2*((buf[jp1][1]-(2.0*buf[j][1]))+buf[jm1][1])))+(dy2ty1*((ue[jp1][1]-(2.0*ue[j][1]))+ue[jm1][1])));
				forcing[k][j][i][2]=(((forcing[k][j][i][2]-(ty2*(((ue[jp1][2]*buf[jp1][2])+(c2*(ue[jp1][4]-q[jp1])))-((ue[jm1][2]*buf[jm1][2])+(c2*(ue[jm1][4]-q[jm1]))))))+(yycon1*((buf[jp1][2]-(2.0*buf[j][2]))+buf[jm1][2])))+(dy3ty1*((ue[jp1][2]-(2.0*ue[j][2]))+ue[jm1][2])));
				forcing[k][j][i][3]=(((forcing[k][j][i][3]-(ty2*((ue[jp1][3]*buf[jp1][2])-(ue[jm1][3]*buf[jm1][2]))))+(yycon2*((buf[jp1][3]-(2.0*buf[j][3]))+buf[jm1][3])))+(dy4ty1*((ue[jp1][3]-(2.0*ue[j][3]))+ue[jm1][3])));
				forcing[k][j][i][4]=(((((forcing[k][j][i][4]-(ty2*((buf[jp1][2]*((c1*ue[jp1][4])-(c2*q[jp1])))-(buf[jm1][2]*((c1*ue[jm1][4])-(c2*q[jm1]))))))+((0.5*yycon3)*((buf[jp1][0]-(2.0*buf[j][0]))+buf[jm1][0])))+(yycon4*((cuf[jp1]-(2.0*cuf[j]))+cuf[jm1])))+(yycon5*((buf[jp1][4]-(2.0*buf[j][4]))+buf[jm1][4])))+(dy5ty1*((ue[jp1][4]-(2.0*ue[j][4]))+ue[jm1][4])));
			}
			/* --------------------------------------------------------------------- */
			/* Fourth-order dissipation                       */
			/* --------------------------------------------------------------------- */
			#pragma cetus private(j, m) 
			#pragma loop name exact_rhs#2#0#2 
			#pragma cetus parallel 
			/*
			Disabled due to low profitability: #pragma omp parallel for private(j, m)
			*/
			for (m=0; m<5; m ++ )
			{
				j=1;
				forcing[k][j][i][m]=(forcing[k][j][i][m]-(dssp*(((5.0*ue[j][m])-(4.0*ue[j+1][m]))+ue[j+2][m])));
				j=2;
				forcing[k][j][i][m]=(forcing[k][j][i][m]-(dssp*((((( - 4.0)*ue[j-1][m])+(6.0*ue[j][m]))-(4.0*ue[j+1][m]))+ue[j+2][m])));
			}
			#pragma cetus private(j, m) 
			#pragma loop name exact_rhs#2#0#3 
			#pragma cetus parallel 
			#pragma omp parallel for if((10000<(-107L+(18L*grid_points[1L])))) private(j, m)
			for (j=3; j<=(grid_points[1]-4); j ++ )
			{
				#pragma cetus private(m) 
				#pragma loop name exact_rhs#2#0#3#0 
				for (m=0; m<5; m ++ )
				{
					forcing[k][j][i][m]=(forcing[k][j][i][m]-(dssp*((((ue[j-2][m]-(4.0*ue[j-1][m]))+(6.0*ue[j][m]))-(4.0*ue[j+1][m]))+ue[j+2][m])));
				}
			}
			#pragma cetus private(j, m) 
			#pragma loop name exact_rhs#2#0#4 
			#pragma cetus parallel 
			/*
			Disabled due to low profitability: #pragma omp parallel for private(j, m)
			*/
			for (m=0; m<5; m ++ )
			{
				j=(grid_points[1]-3);
				forcing[k][j][i][m]=(forcing[k][j][i][m]-(dssp*(((ue[j-2][m]-(4.0*ue[j-1][m]))+(6.0*ue[j][m]))-(4.0*ue[j+1][m]))));
				j=(grid_points[1]-2);
				forcing[k][j][i][m]=(forcing[k][j][i][m]-(dssp*((ue[j-2][m]-(4.0*ue[j-1][m]))+(5.0*ue[j][m]))));
			}
		}
	}
	/* --------------------------------------------------------------------- */
	/* zeta-direction flux differences                       */
	/* --------------------------------------------------------------------- */
	#pragma cetus private(dtpp, eta, i, j, k, km1, kp1, m, xi, zeta) 
	#pragma loop name exact_rhs#3 
	for (j=1; j<=(grid_points[1]-2); j ++ )
	{
		eta=(((double)j)*dnym1);
		#pragma cetus private(dtpp, i, k, km1, kp1, m, xi, zeta) 
		#pragma loop name exact_rhs#3#0 
		for (i=1; i<=(grid_points[0]-2); i ++ )
		{
			xi=(((double)i)*dnxm1);
			#pragma cetus private(dtpp, k, m, zeta) 
			#pragma loop name exact_rhs#3#0#0 
			for (k=0; k<=(grid_points[2]-1); k ++ )
			{
				zeta=(((double)k)*dnzm1);
				exact_solution(xi, eta, zeta, dtemp);
				#pragma cetus private(m) 
				#pragma loop name exact_rhs#3#0#0#0 
				#pragma cetus parallel 
				/*
				Disabled due to low profitability: #pragma omp parallel for private(m)
				*/
				for (m=0; m<5; m ++ )
				{
					ue[k][m]=dtemp[m];
				}
				dtpp=(1.0/dtemp[0]);
				#pragma cetus private(m) 
				#pragma loop name exact_rhs#3#0#0#1 
				#pragma cetus parallel 
				/*
				Disabled due to low profitability: #pragma omp parallel for private(m)
				*/
				for (m=1; m<5; m ++ )
				{
					buf[k][m]=(dtpp*dtemp[m]);
				}
				cuf[k]=(buf[k][3]*buf[k][3]);
				buf[k][0]=((cuf[k]+(buf[k][1]*buf[k][1]))+(buf[k][2]*buf[k][2]));
				q[k]=(0.5*(((buf[k][1]*ue[k][1])+(buf[k][2]*ue[k][2]))+(buf[k][3]*ue[k][3])));
			}
			#pragma cetus private(k, km1, kp1) 
			#pragma loop name exact_rhs#3#0#1 
			#pragma cetus parallel 
			#pragma omp parallel for if((10000<(-17L+(9L*grid_points[2L])))) private(k, km1, kp1)
			for (k=1; k<=(grid_points[2]-2); k ++ )
			{
				km1=(k-1);
				kp1=(k+1);
				forcing[k][j][i][0]=((forcing[k][j][i][0]-(tz2*(ue[kp1][3]-ue[km1][3])))+(dz1tz1*((ue[kp1][0]-(2.0*ue[k][0]))+ue[km1][0])));
				forcing[k][j][i][1]=(((forcing[k][j][i][1]-(tz2*((ue[kp1][1]*buf[kp1][3])-(ue[km1][1]*buf[km1][3]))))+(zzcon2*((buf[kp1][1]-(2.0*buf[k][1]))+buf[km1][1])))+(dz2tz1*((ue[kp1][1]-(2.0*ue[k][1]))+ue[km1][1])));
				forcing[k][j][i][2]=(((forcing[k][j][i][2]-(tz2*((ue[kp1][2]*buf[kp1][3])-(ue[km1][2]*buf[km1][3]))))+(zzcon2*((buf[kp1][2]-(2.0*buf[k][2]))+buf[km1][2])))+(dz3tz1*((ue[kp1][2]-(2.0*ue[k][2]))+ue[km1][2])));
				forcing[k][j][i][3]=(((forcing[k][j][i][3]-(tz2*(((ue[kp1][3]*buf[kp1][3])+(c2*(ue[kp1][4]-q[kp1])))-((ue[km1][3]*buf[km1][3])+(c2*(ue[km1][4]-q[km1]))))))+(zzcon1*((buf[kp1][3]-(2.0*buf[k][3]))+buf[km1][3])))+(dz4tz1*((ue[kp1][3]-(2.0*ue[k][3]))+ue[km1][3])));
				forcing[k][j][i][4]=(((((forcing[k][j][i][4]-(tz2*((buf[kp1][3]*((c1*ue[kp1][4])-(c2*q[kp1])))-(buf[km1][3]*((c1*ue[km1][4])-(c2*q[km1]))))))+((0.5*zzcon3)*((buf[kp1][0]-(2.0*buf[k][0]))+buf[km1][0])))+(zzcon4*((cuf[kp1]-(2.0*cuf[k]))+cuf[km1])))+(zzcon5*((buf[kp1][4]-(2.0*buf[k][4]))+buf[km1][4])))+(dz5tz1*((ue[kp1][4]-(2.0*ue[k][4]))+ue[km1][4])));
			}
			/* --------------------------------------------------------------------- */
			/* Fourth-order dissipation */
			/* --------------------------------------------------------------------- */
			#pragma cetus private(k, m) 
			#pragma loop name exact_rhs#3#0#2 
			#pragma cetus parallel 
			/*
			Disabled due to low profitability: #pragma omp parallel for private(k, m)
			*/
			for (m=0; m<5; m ++ )
			{
				k=1;
				forcing[k][j][i][m]=(forcing[k][j][i][m]-(dssp*(((5.0*ue[k][m])-(4.0*ue[k+1][m]))+ue[k+2][m])));
				k=2;
				forcing[k][j][i][m]=(forcing[k][j][i][m]-(dssp*((((( - 4.0)*ue[k-1][m])+(6.0*ue[k][m]))-(4.0*ue[k+1][m]))+ue[k+2][m])));
			}
			#pragma cetus private(k, m) 
			#pragma loop name exact_rhs#3#0#3 
			#pragma cetus parallel 
			#pragma omp parallel for if((10000<(-107L+(18L*grid_points[2L])))) private(k, m)
			for (k=3; k<=(grid_points[2]-4); k ++ )
			{
				#pragma cetus private(m) 
				#pragma loop name exact_rhs#3#0#3#0 
				for (m=0; m<5; m ++ )
				{
					forcing[k][j][i][m]=(forcing[k][j][i][m]-(dssp*((((ue[k-2][m]-(4.0*ue[k-1][m]))+(6.0*ue[k][m]))-(4.0*ue[k+1][m]))+ue[k+2][m])));
				}
			}
			#pragma cetus private(k, m) 
			#pragma loop name exact_rhs#3#0#4 
			#pragma cetus parallel 
			/*
			Disabled due to low profitability: #pragma omp parallel for private(k, m)
			*/
			for (m=0; m<5; m ++ )
			{
				k=(grid_points[2]-3);
				forcing[k][j][i][m]=(forcing[k][j][i][m]-(dssp*(((ue[k-2][m]-(4.0*ue[k-1][m]))+(6.0*ue[k][m]))-(4.0*ue[k+1][m]))));
				k=(grid_points[2]-2);
				forcing[k][j][i][m]=(forcing[k][j][i][m]-(dssp*((ue[k-2][m]-(4.0*ue[k-1][m]))+(5.0*ue[k][m]))));
			}
		}
	}
	/* --------------------------------------------------------------------- */
	/* now change the sign of the forcing function,  */
	/* --------------------------------------------------------------------- */
	#pragma cetus private(i, j, k, m) 
	#pragma loop name exact_rhs#4 
	#pragma cetus parallel 
	#pragma omp parallel for if((10000<(((((((-137L+(72L*grid_points[0L]))+(66L*grid_points[1L]))+(69L*grid_points[2L]))+((-36L*grid_points[0L])*grid_points[1L]))+((-36L*grid_points[0L])*grid_points[2L]))+((-33L*grid_points[1L])*grid_points[2L]))+(((18L*grid_points[0L])*grid_points[1L])*grid_points[2L])))) private(i, j, k, m)
	for (k=1; k<=(grid_points[2]-2); k ++ )
	{
		#pragma cetus private(i, j, m) 
		#pragma loop name exact_rhs#4#0 
		for (j=1; j<=(grid_points[1]-2); j ++ )
		{
			#pragma cetus private(i, m) 
			#pragma loop name exact_rhs#4#0#0 
			for (i=1; i<=(grid_points[0]-2); i ++ )
			{
				#pragma cetus private(m) 
				#pragma loop name exact_rhs#4#0#0#0 
				for (m=0; m<5; m ++ )
				{
					forcing[k][j][i][m]=(( - 1.0)*forcing[k][j][i][m]);
				}
			}
		}
	}
}
