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
#include <stdio.h>
#include "applu.incl"
void pintgr()
{
	/* --------------------------------------------------------------------- */
	/* local variables */
	/* --------------------------------------------------------------------- */
	int i, j, k;
	int ibeg, ifin, ifin1;
	int jbeg, jfin, jfin1;
	double phi1[(162+2)][(162+2)];
	double phi2[(162+2)][(162+2)];
	double frc1, frc2, frc3;
	/* --------------------------------------------------------------------- */
	/* set up the sub-domains for integeration in each processor */
	/* --------------------------------------------------------------------- */
	ibeg=ii1;
	ifin=ii2;
	jbeg=ji1;
	jfin=ji2;
	ifin1=(ifin-1);
	jfin1=(jfin-1);
	#pragma cetus private(i, j, k) 
	#pragma loop name pintgr#0 
	#pragma cetus parallel 
	#pragma omp parallel for if((10000<((((((1L+(-3L*ji1))+(3L*ji2))+((6L*ii1)*ji1))+((-6L*ii1)*ji2))+((-6L*ii2)*ji1))+((6L*ii2)*ji2)))) private(i, j, k)
	for (j=jbeg; j<jfin; j ++ )
	{
		#pragma cetus private(i, k) 
		#pragma loop name pintgr#0#0 
		for (i=ibeg; i<ifin; i ++ )
		{
			k=ki1;
			phi1[j][i]=(0.4*(u[k][j][i][4]-((0.5*(((u[k][j][i][1]*u[k][j][i][1])+(u[k][j][i][2]*u[k][j][i][2]))+(u[k][j][i][3]*u[k][j][i][3])))/u[k][j][i][0])));
			k=(ki2-1);
			phi2[j][i]=(0.4*(u[k][j][i][4]-((0.5*(((u[k][j][i][1]*u[k][j][i][1])+(u[k][j][i][2]*u[k][j][i][2]))+(u[k][j][i][3]*u[k][j][i][3])))/u[k][j][i][0])));
		}
	}
	frc1=0.0;
	#pragma cetus private(i, j) 
	#pragma loop name pintgr#1 
	#pragma cetus reduction(+: frc1) 
	#pragma cetus parallel 
	#pragma omp parallel for if((10000<((((((1L+(3L*ii1))+(-3L*ii2))+((3L*ii1)*ji1))+((-3L*ii1)*ji2))+((-3L*ii2)*ji1))+((3L*ii2)*ji2)))) private(i, j) reduction(+: frc1)
	for (j=jbeg; j<jfin1; j ++ )
	{
		#pragma cetus private(i) 
		#pragma loop name pintgr#1#0 
		/* #pragma cetus reduction(+: frc1)  */
		for (i=ibeg; i<ifin1; i ++ )
		{
			frc1=(frc1+(((((((phi1[j][i]+phi1[j][i+1])+phi1[j+1][i])+phi1[j+1][i+1])+phi2[j][i])+phi2[j][i+1])+phi2[j+1][i])+phi2[j+1][i+1]));
		}
	}
	frc1=((dxi*deta)*frc1);
	#pragma cetus private(i, k) 
	#pragma loop name pintgr#2 
	#pragma cetus parallel 
	#pragma omp parallel for if((10000<((((((1L+(-3L*ki1))+(3L*ki2))+((3L*ii1)*ki1))+((-3L*ii1)*ki2))+((-3L*ii2)*ki1))+((3L*ii2)*ki2)))) private(i, k)
	for (k=ki1; k<ki2; k ++ )
	{
		#pragma cetus private(i) 
		#pragma loop name pintgr#2#0 
		for (i=ibeg; i<ifin; i ++ )
		{
			phi1[k][i]=(0.4*(u[k][jbeg][i][4]-((0.5*(((u[k][jbeg][i][1]*u[k][jbeg][i][1])+(u[k][jbeg][i][2]*u[k][jbeg][i][2]))+(u[k][jbeg][i][3]*u[k][jbeg][i][3])))/u[k][jbeg][i][0])));
		}
	}
	#pragma cetus private(i, k) 
	#pragma loop name pintgr#3 
	#pragma cetus parallel 
	#pragma omp parallel for if((10000<((((((1L+(-3L*ki1))+(3L*ki2))+((3L*ii1)*ki1))+((-3L*ii1)*ki2))+((-3L*ii2)*ki1))+((3L*ii2)*ki2)))) private(i, k)
	for (k=ki1; k<ki2; k ++ )
	{
		#pragma cetus private(i) 
		#pragma loop name pintgr#3#0 
		for (i=ibeg; i<ifin; i ++ )
		{
			phi2[k][i]=(0.4*(u[k][jfin-1][i][4]-((0.5*(((u[k][jfin-1][i][1]*u[k][jfin-1][i][1])+(u[k][jfin-1][i][2]*u[k][jfin-1][i][2]))+(u[k][jfin-1][i][3]*u[k][jfin-1][i][3])))/u[k][jfin-1][i][0])));
		}
	}
	frc2=0.0;
	#pragma cetus private(i, k) 
	#pragma loop name pintgr#4 
	#pragma cetus reduction(+: frc2) 
	#pragma cetus parallel 
	#pragma omp parallel for if((10000<((((((1L+(3L*ii1))+(-3L*ii2))+((3L*ii1)*ki1))+((-3L*ii1)*ki2))+((-3L*ii2)*ki1))+((3L*ii2)*ki2)))) private(i, k) reduction(+: frc2)
	for (k=ki1; k<(ki2-1); k ++ )
	{
		#pragma cetus private(i) 
		#pragma loop name pintgr#4#0 
		/* #pragma cetus reduction(+: frc2)  */
		for (i=ibeg; i<ifin1; i ++ )
		{
			frc2=(frc2+(((((((phi1[k][i]+phi1[k][i+1])+phi1[k+1][i])+phi1[k+1][i+1])+phi2[k][i])+phi2[k][i+1])+phi2[k+1][i])+phi2[k+1][i+1]));
		}
	}
	frc2=((dxi*dzeta)*frc2);
	#pragma cetus private(j, k) 
	#pragma loop name pintgr#5 
	#pragma cetus parallel 
	#pragma omp parallel for if((10000<((((((1L+(-3L*ki1))+(3L*ki2))+((3L*ji1)*ki1))+((-3L*ji1)*ki2))+((-3L*ji2)*ki1))+((3L*ji2)*ki2)))) private(j, k)
	for (k=ki1; k<ki2; k ++ )
	{
		#pragma cetus private(j) 
		#pragma loop name pintgr#5#0 
		for (j=jbeg; j<jfin; j ++ )
		{
			phi1[k][j]=(0.4*(u[k][j][ibeg][4]-((0.5*(((u[k][j][ibeg][1]*u[k][j][ibeg][1])+(u[k][j][ibeg][2]*u[k][j][ibeg][2]))+(u[k][j][ibeg][3]*u[k][j][ibeg][3])))/u[k][j][ibeg][0])));
		}
	}
	#pragma cetus private(j, k) 
	#pragma loop name pintgr#6 
	#pragma cetus parallel 
	#pragma omp parallel for if((10000<((((((1L+(-3L*ki1))+(3L*ki2))+((3L*ji1)*ki1))+((-3L*ji1)*ki2))+((-3L*ji2)*ki1))+((3L*ji2)*ki2)))) private(j, k)
	for (k=ki1; k<ki2; k ++ )
	{
		#pragma cetus private(j) 
		#pragma loop name pintgr#6#0 
		for (j=jbeg; j<jfin; j ++ )
		{
			phi2[k][j]=(0.4*(u[k][j][ifin-1][4]-((0.5*(((u[k][j][ifin-1][1]*u[k][j][ifin-1][1])+(u[k][j][ifin-1][2]*u[k][j][ifin-1][2]))+(u[k][j][ifin-1][3]*u[k][j][ifin-1][3])))/u[k][j][ifin-1][0])));
		}
	}
	frc3=0.0;
	#pragma cetus private(j, k) 
	#pragma loop name pintgr#7 
	#pragma cetus reduction(+: frc3) 
	#pragma cetus parallel 
	#pragma omp parallel for if((10000<((((((1L+(3L*ji1))+(-3L*ji2))+((3L*ji1)*ki1))+((-3L*ji1)*ki2))+((-3L*ji2)*ki1))+((3L*ji2)*ki2)))) private(j, k) reduction(+: frc3)
	for (k=ki1; k<(ki2-1); k ++ )
	{
		#pragma cetus private(j) 
		#pragma loop name pintgr#7#0 
		/* #pragma cetus reduction(+: frc3)  */
		for (j=jbeg; j<jfin1; j ++ )
		{
			frc3=(frc3+(((((((phi1[k][j]+phi1[k][j+1])+phi1[k+1][j])+phi1[k+1][j+1])+phi2[k][j])+phi2[k][j+1])+phi2[k+1][j])+phi2[k+1][j+1]));
		}
	}
	frc3=((deta*dzeta)*frc3);
	frc=(0.25*((frc1+frc2)+frc3));
	/* printf("\n\n     surface integral = %12.5E\n\n\n", frc); */
}
