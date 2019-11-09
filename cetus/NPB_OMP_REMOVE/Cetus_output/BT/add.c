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
#include "timers.h"
/* --------------------------------------------------------------------- */
/* addition of update to the vector u */
/* --------------------------------------------------------------------- */
void add()
{
	int i, j, k, m;
	if (timeron)
	{
		timer_start(11);
	}
	#pragma cetus private(i, j, k, m) 
	#pragma loop name add#0 
	#pragma cetus parallel 
	#pragma omp parallel for if((10000<(((((((-137L+(72L*grid_points[0L]))+(66L*grid_points[1L]))+(69L*grid_points[2L]))+((-36L*grid_points[0L])*grid_points[1L]))+((-36L*grid_points[0L])*grid_points[2L]))+((-33L*grid_points[1L])*grid_points[2L]))+(((18L*grid_points[0L])*grid_points[1L])*grid_points[2L])))) private(i, j, k, m)
	for (k=1; k<=(grid_points[2]-2); k ++ )
	{
		#pragma cetus private(i, j, m) 
		#pragma loop name add#0#0 
		for (j=1; j<=(grid_points[1]-2); j ++ )
		{
			#pragma cetus private(i, m) 
			#pragma loop name add#0#0#0 
			for (i=1; i<=(grid_points[0]-2); i ++ )
			{
				#pragma cetus private(m) 
				#pragma loop name add#0#0#0#0 
				for (m=0; m<5; m ++ )
				{
					u[k][j][i][m]=(u[k][j][i][m]+rhs[k][j][i][m]);
				}
			}
		}
	}
	if (timeron)
	{
		timer_stop(11);
	}
}
