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
#include <stdlib.h>
#include "applu.incl"
void domain()
{
	/* --------------------------------------------------------------------- */
	/* local variables */
	/* --------------------------------------------------------------------- */
	nx=nx0;
	ny=ny0;
	nz=nz0;
	/* --------------------------------------------------------------------- */
	/* check the sub-domain size */
	/* --------------------------------------------------------------------- */
	if (((nx<4)||(ny<4))||(nz<4))
	{
		printf("     SUBDOMAIN SIZE IS TOO SMALL - \n""     ADJUST PROBLEM SIZE OR NUMBER OF PROCESSORS\n""     SO THAT NX, NY AND NZ ARE GREATER THAN OR EQUAL\n""     TO 4 THEY ARE CURRENTLY%3d%3d%3d\n", nx, ny, nz);
		exit(1);
	}
	if (((nx>162)||(ny>162))||(nz>162))
	{
		printf("     SUBDOMAIN SIZE IS TOO LARGE - \n""     ADJUST PROBLEM SIZE OR NUMBER OF PROCESSORS\n""     SO THAT NX, NY AND NZ ARE LESS THAN OR EQUAL TO \n""     ISIZ1, ISIZ2 AND ISIZ3 RESPECTIVELY.  THEY ARE\n""     CURRENTLYi%4d%4d%4d\n", nx, ny, nz);
		exit(1);
	}
	/* --------------------------------------------------------------------- */
	/* set up the start and end in i and j extents for all processors */
	/* --------------------------------------------------------------------- */
	ist=1;
	iend=(nx-1);
	jst=1;
	jend=(ny-1);
	ii1=1;
	ii2=(nx0-1);
	ji1=1;
	ji2=(ny0-2);
	ki1=2;
	ki2=(nz0-1);
}
