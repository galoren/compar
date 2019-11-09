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
/* --------------------------------------------------------------------- */
/* program SP */
/* --------------------------------------------------------------------- */
#include <stdio.h>
#include <stdlib.h>
#include "header.h"
#include "print_results.h"
/* commonglobal */
int grid_points[3], nx2, ny2, nz2;
logical timeron;
/* commonconstants */
double tx1, tx2, tx3, ty1, ty2, ty3, tz1, tz2, tz3, dx1, dx2, dx3, dx4, dx5, dy1, dy2, dy3, dy4, dy5, dz1, dz2, dz3, dz4, dz5, dssp, dt, ce[5][13], dxmax, dymax, dzmax, xxcon1, xxcon2, xxcon3, xxcon4, xxcon5, dx1tx1, dx2tx1, dx3tx1, dx4tx1, dx5tx1, yycon1, yycon2, yycon3, yycon4, yycon5, dy1ty1, dy2ty1, dy3ty1, dy4ty1, dy5ty1, zzcon1, zzcon2, zzcon3, zzcon4, zzcon5, dz1tz1, dz2tz1, dz3tz1, dz4tz1, dz5tz1, dnxm1, dnym1, dnzm1, c1c2, c1c5, c3c4, c1345, conz1, c1, c2, c3, c4, c5, c4dssp, c5dssp, dtdssp, dttx1, bt, dttx2, dtty1, dtty2, dttz1, dttz2, c2dttx1, c2dtty1, c2dttz1, comz1, comz4, comz5, comz6, c3c4tx3, c3c4ty3, c3c4tz3, c2iv, con43, con16;
/* commonfields */
double u[162][(((162/2)*2)+1)][(((162/2)*2)+1)][5];
double us[162][(((162/2)*2)+1)][(((162/2)*2)+1)];
double vs[162][(((162/2)*2)+1)][(((162/2)*2)+1)];
double ws[162][(((162/2)*2)+1)][(((162/2)*2)+1)];
double qs[162][(((162/2)*2)+1)][(((162/2)*2)+1)];
double rho_i[162][(((162/2)*2)+1)][(((162/2)*2)+1)];
double speed[162][(((162/2)*2)+1)][(((162/2)*2)+1)];
double square[162][(((162/2)*2)+1)][(((162/2)*2)+1)];
double rhs[162][(((162/2)*2)+1)][(((162/2)*2)+1)][5];
double forcing[162][(((162/2)*2)+1)][(((162/2)*2)+1)][5];
/* commonwork_1d */
double cv[162];
double rhon[162];
double rhos[162];
double rhoq[162];
double cuf[162];
double q[162];
double ue[162][5];
double buf[162][5];
/* commonwork_lhs */
double lhs[(((162/2)*2)+1)][(((162/2)*2)+1)][5];
double lhsp[(((162/2)*2)+1)][(((162/2)*2)+1)][5];
double lhsm[(((162/2)*2)+1)][(((162/2)*2)+1)][5];
int main(int argc, char * argv[])
{
	int i, niter, step, n3;
	double mflops, t, tmax, trecs[(15+1)];
	logical verified;
	char Class;
	char * t_names[(15+1)];
	/* --------------------------------------------------------------------- */
	/* Read input file (if it exists), else take */
	/* defaults from parameters */
	/* --------------------------------------------------------------------- */
	FILE * fp;
	if ((fp=fopen("timer.flag", "r"))!=((void * )0))
	{
		timeron=true;
		t_names[1]="total";
		t_names[2]="rhsx";
		t_names[3]="rhsy";
		t_names[4]="rhsz";
		t_names[5]="rhs";
		t_names[6]="xsolve";
		t_names[7]="ysolve";
		t_names[8]="zsolve";
		t_names[9]="redist1";
		t_names[10]="redist2";
		t_names[14]="tzetar";
		t_names[13]="ninvr";
		t_names[12]="pinvr";
		t_names[11]="txinvr";
		t_names[15]="add";
		fclose(fp);
	}
	else
	{
		timeron=false;
	}
	printf("\n\n NAS Parallel Benchmarks (NPB3.3-OMP-C) - SP Benchmark\n\n");
	if ((fp=fopen("inputsp.data", "r"))!=((void * )0))
	{
		int result;
		printf(" Reading from input file inputsp.data\n");
		result=fscanf(fp, "%d",  & niter);
		while (fgetc(fp)!='\n')
		{
			;
		}
		result=fscanf(fp, "%lf",  & dt);
		while (fgetc(fp)!='\n')
		{
			;
		}
		result=fscanf(fp, "%d%d%d",  & grid_points[0],  & grid_points[1],  & grid_points[2]);
		fclose(fp);
	}
	else
	{
		printf(" No input file inputsp.data. Using compiled defaults\n");
		niter=400;
		dt=6.7E-4;
		grid_points[0]=162;
		grid_points[1]=162;
		grid_points[2]=162;
	}
	printf(" Size: %4dx%4dx%4d\n", grid_points[0], grid_points[1], grid_points[2]);
	printf(" Iterations: %4d    dt:  %11.7f\n", niter, dt);
	printf(" Number of available threads: %5d\n", omp_get_max_threads());
	printf("\n");
	if (((grid_points[0]>162)||(grid_points[1]>162))||(grid_points[2]>162))
	{
		printf(" %d, %d, %d\n", grid_points[0], grid_points[1], grid_points[2]);
		printf(" Problem size too big for compiled array sizes\n");
		return 0;
	}
	nx2=(grid_points[0]-2);
	ny2=(grid_points[1]-2);
	nz2=(grid_points[2]-2);
	set_constants();
	#pragma cetus private(i) 
	#pragma loop name main#0 
	for (i=1; i<=15; i ++ )
	{
		timer_clear(i);
	}
	exact_rhs();
	initialize();
	/* --------------------------------------------------------------------- */
	/* do one time step to touch all code, and reinitialize */
	/* --------------------------------------------------------------------- */
	adi();
	initialize();
	#pragma cetus private(i) 
	#pragma loop name main#1 
	for (i=1; i<=15; i ++ )
	{
		timer_clear(i);
	}
	timer_start(1);
	#pragma cetus private(step) 
	#pragma loop name main#2 
	for (step=1; step<=niter; step ++ )
	{
		if (((step%20)==0)||(step==1))
		{
			printf(" Time step %4d\n", step);
		}
		adi();
	}
	timer_stop(1);
	tmax=timer_read(1);
	verify(niter,  & Class,  & verified);
	if (tmax!=0.0)
	{
		n3=((grid_points[0]*grid_points[1])*grid_points[2]);
		t=(((grid_points[0]+grid_points[1])+grid_points[2])/3.0);
		mflops=((((((881.174*((double)n3))-(4683.91*(t*t)))+(11484.5*t))-19272.4)*((double)niter))/(tmax*1000000.0));
	}
	else
	{
		mflops=0.0;
	}
	print_results("SP", Class, grid_points[0], grid_points[1], grid_points[2], niter, tmax, mflops, "          floating point", verified, "3.3.1", "08 Nov 2019", "icc", "$(CC)", "-lm", "-I../common", "-g -Wall -O3 -qopenmp -mcmodel=medium", "-O3 -qopenmp -mcmodel=medium", "(none)");
	/* --------------------------------------------------------------------- */
	/* More timers */
	/* --------------------------------------------------------------------- */
	if (timeron)
	{
		#pragma cetus private(i) 
		#pragma loop name main#3 
		for (i=1; i<=15; i ++ )
		{
			trecs[i]=timer_read(i);
		}
		if (tmax==0.0)
		{
			tmax=1.0;
		}
		printf("  SECTION   Time (secs)\n");
		#pragma cetus private(i, t) 
		#pragma loop name main#4 
		for (i=1; i<=15; i ++ )
		{
			printf("  %-8s:%9.3f  (%6.2f%%)\n", t_names[i], trecs[i], (trecs[i]*100.0)/tmax);
			if (i==5)
			{
				t=((trecs[2]+trecs[3])+trecs[4]);
				printf("    --> %8s:%9.3f  (%6.2f%%)\n", "sub-rhs", t, (t*100.0)/tmax);
				t=(trecs[5]-t);
				printf("    --> %8s:%9.3f  (%6.2f%%)\n", "rest-rhs", t, (t*100.0)/tmax);
			}
			else
			{
				if (i==8)
				{
					t=((trecs[8]-trecs[9])-trecs[10]);
					printf("    --> %8s:%9.3f  (%6.2f%%)\n", "sub-zsol", t, (t*100.0)/tmax);
				}
				else
				{
					if (i==10)
					{
						t=(trecs[9]+trecs[10]);
						printf("    --> %8s:%9.3f  (%6.2f%%)\n", "redist", t, (t*100.0)/tmax);
					}
				}
			}
		}
	}
	return 0;
}
