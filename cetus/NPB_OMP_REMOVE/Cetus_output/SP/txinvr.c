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
/* block-diagonal matrix-vector multiplication                   */
/* --------------------------------------------------------------------- */
void txinvr()
{
	int i, j, k;
	double t1, t2, t3, ac, ru1, uu, vv, ww, r1, r2, r3, r4, r5, ac2inv;
	if (timeron)
	{
		timer_start(11);
	}
	#pragma cetus private(ac, ac2inv, i, j, k, r1, r2, r3, r4, r5, ru1, t1, t2, t3, uu, vv, ww) 
	#pragma loop name txinvr#0 
	#pragma cetus parallel 
	#pragma omp parallel for if((10000<(((1L+(3L*nz2))+((3L*ny2)*nz2))+(((21L*nx2)*ny2)*nz2)))) private(ac, ac2inv, i, j, k, r1, r2, r3, r4, r5, ru1, t1, t2, t3, uu, vv, ww)
	for (k=1; k<=nz2; k ++ )
	{
		#pragma cetus private(ac, ac2inv, i, j, r1, r2, r3, r4, r5, ru1, t1, t2, t3, uu, vv, ww) 
		#pragma loop name txinvr#0#0 
		for (j=1; j<=ny2; j ++ )
		{
			#pragma cetus private(ac, ac2inv, i, r1, r2, r3, r4, r5, ru1, t1, t2, t3, uu, vv, ww) 
			#pragma loop name txinvr#0#0#0 
			for (i=1; i<=nx2; i ++ )
			{
				ru1=rho_i[k][j][i];
				uu=us[k][j][i];
				vv=vs[k][j][i];
				ww=ws[k][j][i];
				ac=speed[k][j][i];
				ac2inv=(ac*ac);
				r1=rhs[k][j][i][0];
				r2=rhs[k][j][i][1];
				r3=rhs[k][j][i][2];
				r4=rhs[k][j][i][3];
				r5=rhs[k][j][i][4];
				t1=((c2/ac2inv)*(((((qs[k][j][i]*r1)-(uu*r2))-(vv*r3))-(ww*r4))+r5));
				t2=((bt*ru1)*((uu*r1)-r2));
				t3=(((bt*ru1)*ac)*t1);
				rhs[k][j][i][0]=(r1-t1);
				rhs[k][j][i][1]=(( - ru1)*((ww*r1)-r4));
				rhs[k][j][i][2]=(ru1*((vv*r1)-r3));
				rhs[k][j][i][3]=(( - t2)+t3);
				rhs[k][j][i][4]=(t2+t3);
			}
		}
	}
	if (timeron)
	{
		timer_stop(11);
	}
}
