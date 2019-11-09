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
#include "applu.incl"
void setcoeff()
{
	/* --------------------------------------------------------------------- */
	/* local variables */
	/* --------------------------------------------------------------------- */
	/* --------------------------------------------------------------------- */
	/* set up coefficients */
	/* --------------------------------------------------------------------- */
	dxi=(1.0/(nx0-1));
	deta=(1.0/(ny0-1));
	dzeta=(1.0/(nz0-1));
	tx1=(1.0/(dxi*dxi));
	tx2=(1.0/(2.0*dxi));
	tx3=(1.0/dxi);
	ty1=(1.0/(deta*deta));
	ty2=(1.0/(2.0*deta));
	ty3=(1.0/deta);
	tz1=(1.0/(dzeta*dzeta));
	tz2=(1.0/(2.0*dzeta));
	tz3=(1.0/dzeta);
	/* --------------------------------------------------------------------- */
	/* diffusion coefficients */
	/* --------------------------------------------------------------------- */
	dx1=0.75;
	dx2=dx1;
	dx3=dx1;
	dx4=dx1;
	dx5=dx1;
	dy1=0.75;
	dy2=dy1;
	dy3=dy1;
	dy4=dy1;
	dy5=dy1;
	dz1=1.0;
	dz2=dz1;
	dz3=dz1;
	dz4=dz1;
	dz5=dz1;
	/* --------------------------------------------------------------------- */
	/* fourth difference dissipation */
	/* --------------------------------------------------------------------- */
	dssp=(((((dx1>dy1) ? dx1 : dy1)>dz1) ? ((dx1>dy1) ? dx1 : dy1) : dz1)/4.0);
	/* --------------------------------------------------------------------- */
	/* coefficients of the exact solution to the first pde */
	/* --------------------------------------------------------------------- */
	ce[0][0]=2.0;
	ce[0][1]=0.0;
	ce[0][2]=0.0;
	ce[0][3]=4.0;
	ce[0][4]=5.0;
	ce[0][5]=3.0;
	ce[0][6]=0.5;
	ce[0][7]=0.02;
	ce[0][8]=0.01;
	ce[0][9]=0.03;
	ce[0][10]=0.5;
	ce[0][11]=0.4;
	ce[0][12]=0.3;
	/* --------------------------------------------------------------------- */
	/* coefficients of the exact solution to the second pde */
	/* --------------------------------------------------------------------- */
	ce[1][0]=1.0;
	ce[1][1]=0.0;
	ce[1][2]=0.0;
	ce[1][3]=0.0;
	ce[1][4]=1.0;
	ce[1][5]=2.0;
	ce[1][6]=3.0;
	ce[1][7]=0.01;
	ce[1][8]=0.03;
	ce[1][9]=0.02;
	ce[1][10]=0.4;
	ce[1][11]=0.3;
	ce[1][12]=0.5;
	/* --------------------------------------------------------------------- */
	/* coefficients of the exact solution to the third pde */
	/* --------------------------------------------------------------------- */
	ce[2][0]=2.0;
	ce[2][1]=2.0;
	ce[2][2]=0.0;
	ce[2][3]=0.0;
	ce[2][4]=0.0;
	ce[2][5]=2.0;
	ce[2][6]=3.0;
	ce[2][7]=0.04;
	ce[2][8]=0.03;
	ce[2][9]=0.05;
	ce[2][10]=0.3;
	ce[2][11]=0.5;
	ce[2][12]=0.4;
	/* --------------------------------------------------------------------- */
	/* coefficients of the exact solution to the fourth pde */
	/* --------------------------------------------------------------------- */
	ce[3][0]=2.0;
	ce[3][1]=2.0;
	ce[3][2]=0.0;
	ce[3][3]=0.0;
	ce[3][4]=0.0;
	ce[3][5]=2.0;
	ce[3][6]=3.0;
	ce[3][7]=0.03;
	ce[3][8]=0.05;
	ce[3][9]=0.04;
	ce[3][10]=0.2;
	ce[3][11]=0.1;
	ce[3][12]=0.3;
	/* --------------------------------------------------------------------- */
	/* coefficients of the exact solution to the fifth pde */
	/* --------------------------------------------------------------------- */
	ce[4][0]=5.0;
	ce[4][1]=4.0;
	ce[4][2]=3.0;
	ce[4][3]=2.0;
	ce[4][4]=0.1;
	ce[4][5]=0.4;
	ce[4][6]=0.3;
	ce[4][7]=0.05;
	ce[4][8]=0.04;
	ce[4][9]=0.03;
	ce[4][10]=0.1;
	ce[4][11]=0.3;
	ce[4][12]=0.2;
}
