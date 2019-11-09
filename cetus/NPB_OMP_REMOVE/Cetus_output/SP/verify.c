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
#include <stdio.h>
#include <math.h>
#include "header.h"
/* --------------------------------------------------------------------- */
/* verification routine                          */
/* --------------------------------------------------------------------- */
void verify(int no_time_steps, char * Class, logical * verified)
{
	double xcrref[5], xceref[5], xcrdif[5], xcedif[5];
	double epsilon, xce[5], xcr[5], dtref = 0.0;
	int m;
	/* --------------------------------------------------------------------- */
	/* tolerance level */
	/* --------------------------------------------------------------------- */
	epsilon=1.0E-8;
	/* --------------------------------------------------------------------- */
	/* compute the error norm and the residual norm, and exit if not printing */
	/* --------------------------------------------------------------------- */
	error_norm(xce);
	compute_rhs();
	rhs_norm(xcr);
	#pragma cetus private(m) 
	#pragma loop name verify#0 
	#pragma cetus parallel 
	/*
	Disabled due to low profitability: #pragma omp parallel for private(m)
	*/
	for (m=0; m<5; m ++ )
	{
		xcr[m]=(xcr[m]/dt);
	}
	( * Class)='U';
	( * verified)=true;
	#pragma cetus private(m) 
	#pragma loop name verify#1 
	#pragma cetus parallel 
	/*
	Disabled due to low profitability: #pragma omp parallel for private(m)
	*/
	for (m=0; m<5; m ++ )
	{
		xcrref[m]=1.0;
		xceref[m]=1.0;
	}
	/* --------------------------------------------------------------------- */
	/* reference data for 12X12X12 grids after 100 time steps,  */
	/* with DT = 1.50e-02 */
	/* --------------------------------------------------------------------- */
	if ((((grid_points[0]==12)&&(grid_points[1]==12))&&(grid_points[2]==12))&&(no_time_steps==100))
	{
		( * Class)='S';
		dtref=0.015;
		/* --------------------------------------------------------------------- */
		/* Reference values of RMS-norms of residual. */
		/* --------------------------------------------------------------------- */
		xcrref[0]=0.02747031545133948;
		xcrref[1]=0.010360746705285417;
		xcrref[2]=0.016235745065095532;
		xcrref[3]=0.015840557224455615;
		xcrref[4]=0.03484904060936246;
		/* --------------------------------------------------------------------- */
		/* Reference values of RMS-norms of solution error. */
		/* --------------------------------------------------------------------- */
		xceref[0]=2.7289258557377225E-5;
		xceref[1]=1.0364446640837285E-5;
		xceref[2]=1.615479828716647E-5;
		xceref[3]=1.57507049944801E-5;
		xceref[4]=3.417766618339053E-5;
		/* --------------------------------------------------------------------- */
		/* reference data for 36X36X36 grids after 400 time steps,  */
		/* with DT = 1.5e-03 */
		/* --------------------------------------------------------------------- */
	}
	else
	{
		if ((((grid_points[0]==36)&&(grid_points[1]==36))&&(grid_points[2]==36))&&(no_time_steps==400))
		{
			( * Class)='W';
			dtref=0.0015;
			/* --------------------------------------------------------------------- */
			/* Reference values of RMS-norms of residual. */
			/* --------------------------------------------------------------------- */
			xcrref[0]=0.001893253733584;
			xcrref[1]=1.717075447775E-4;
			xcrref[2]=2.778153350936E-4;
			xcrref[3]=2.887475409984E-4;
			xcrref[4]=0.003143611161242;
			/* --------------------------------------------------------------------- */
			/* Reference values of RMS-norms of solution error. */
			/* --------------------------------------------------------------------- */
			xceref[0]=7.542088599534E-5;
			xceref[1]=6.512852253086E-6;
			xceref[2]=1.049092285688E-5;
			xceref[3]=1.128838671535E-5;
			xceref[4]=1.212845639773E-4;
			/* --------------------------------------------------------------------- */
			/* reference data for 64X64X64 grids after 400 time steps,  */
			/* with DT = 1.5e-03 */
			/* --------------------------------------------------------------------- */
		}
		else
		{
			if ((((grid_points[0]==64)&&(grid_points[1]==64))&&(grid_points[2]==64))&&(no_time_steps==400))
			{
				( * Class)='A';
				dtref=0.0015;
				/* --------------------------------------------------------------------- */
				/* Reference values of RMS-norms of residual. */
				/* --------------------------------------------------------------------- */
				xcrref[0]=2.4799822399300195;
				xcrref[1]=1.1276337964368832;
				xcrref[2]=1.5028977888770492;
				xcrref[3]=1.421781621169518;
				xcrref[4]=2.129211303513828;
				/* --------------------------------------------------------------------- */
				/* Reference values of RMS-norms of solution error. */
				/* --------------------------------------------------------------------- */
				xceref[0]=1.090014029782055E-4;
				xceref[1]=3.734395176928209E-5;
				xceref[2]=5.009278540654163E-5;
				xceref[3]=4.767109393952825E-5;
				xceref[4]=1.3621613399213E-4;
				/* --------------------------------------------------------------------- */
				/* reference data for 102X102X102 grids after 400 time steps, */
				/* with DT = 1.0e-03 */
				/* --------------------------------------------------------------------- */
			}
			else
			{
				if ((((grid_points[0]==102)&&(grid_points[1]==102))&&(grid_points[2]==102))&&(no_time_steps==400))
				{
					( * Class)='B';
					dtref=0.001;
					/* --------------------------------------------------------------------- */
					/* Reference values of RMS-norms of residual. */
					/* --------------------------------------------------------------------- */
					xcrref[0]=69.03293579998;
					xcrref[1]=30.95134488084;
					xcrref[2]=41.03336647017;
					xcrref[3]=38.64769009604;
					xcrref[4]=56.43482272596;
					/* --------------------------------------------------------------------- */
					/* Reference values of RMS-norms of solution error. */
					/* --------------------------------------------------------------------- */
					xceref[0]=0.009810006190188;
					xceref[1]=0.00102282790567;
					xceref[2]=0.001720597911692;
					xceref[3]=0.001694479428231;
					xceref[4]=0.01847456263981;
					/* --------------------------------------------------------------------- */
					/* reference data for 162X162X162 grids after 400 time steps, */
					/* with DT = 0.67e-03 */
					/* --------------------------------------------------------------------- */
				}
				else
				{
					if ((((grid_points[0]==162)&&(grid_points[1]==162))&&(grid_points[2]==162))&&(no_time_steps==400))
					{
						( * Class)='C';
						dtref=6.7E-4;
						/* --------------------------------------------------------------------- */
						/* Reference values of RMS-norms of residual. */
						/* --------------------------------------------------------------------- */
						xcrref[0]=588.1691581829;
						xcrref[1]=245.4417603569;
						xcrref[2]=329.3829191851;
						xcrref[3]=308.1924971891;
						xcrref[4]=459.7223799176;
						/* --------------------------------------------------------------------- */
						/* Reference values of RMS-norms of solution error. */
						/* --------------------------------------------------------------------- */
						xceref[0]=0.2598120500183;
						xceref[1]=0.02590888922315;
						xceref[2]=0.0513288641632;
						xceref[3]=0.04806073419454;
						xceref[4]=0.5483377491301;
						/* --------------------------------------------------------------------- */
						/* reference data for 408X408X408 grids after 500 time steps, */
						/* with DT = 0.3e-03 */
						/* --------------------------------------------------------------------- */
					}
					else
					{
						if ((((grid_points[0]==408)&&(grid_points[1]==408))&&(grid_points[2]==408))&&(no_time_steps==500))
						{
							( * Class)='D';
							dtref=3.0E-4;
							/* --------------------------------------------------------------------- */
							/* Reference values of RMS-norms of residual. */
							/* --------------------------------------------------------------------- */
							xcrref[0]=10446.96216887;
							xcrref[1]=3204.427762578;
							xcrref[2]=4648.680733032;
							xcrref[3]=4238.923283697;
							xcrref[4]=7588.412036136;
							/* --------------------------------------------------------------------- */
							/* Reference values of RMS-norms of solution error. */
							/* --------------------------------------------------------------------- */
							xceref[0]=5.089471423669;
							xceref[1]=0.5323514855894;
							xceref[2]=1.187051008971;
							xceref[3]=1.083734951938;
							xceref[4]=11.64108338568;
							/* --------------------------------------------------------------------- */
							/* reference data for 1020X1020X1020 grids after 500 time steps, */
							/* with DT = 0.1e-03 */
							/* --------------------------------------------------------------------- */
						}
						else
						{
							if ((((grid_points[0]==1020)&&(grid_points[1]==1020))&&(grid_points[2]==1020))&&(no_time_steps==500))
							{
								( * Class)='E';
								dtref=1.0E-4;
								/* --------------------------------------------------------------------- */
								/* Reference values of RMS-norms of residual. */
								/* --------------------------------------------------------------------- */
								xcrref[0]=62553.87422609;
								xcrref[1]=14953.17020012;
								xcrref[2]=23475.95750586;
								xcrref[3]=20910.99783534;
								xcrref[4]=47704.12841218;
								/* --------------------------------------------------------------------- */
								/* Reference values of RMS-norms of solution error. */
								/* --------------------------------------------------------------------- */
								xceref[0]=67.42735164909;
								xceref[1]=5.390656036938;
								xceref[2]=16.80647196477;
								xceref[3]=15.36963126457;
								xceref[4]=157.5330146156;
							}
							else
							{
								( * verified)=false;
							}
						}
					}
				}
			}
		}
	}
	/* --------------------------------------------------------------------- */
	/* verification test for residuals if gridsize is one of  */
	/* the defined grid sizes above (class .ne. 'U') */
	/* --------------------------------------------------------------------- */
	/* --------------------------------------------------------------------- */
	/* Compute the difference of solution values and the known reference values. */
	/* --------------------------------------------------------------------- */
	#pragma cetus private(m) 
	#pragma loop name verify#2 
	#pragma cetus parallel 
	/*
	Disabled due to low profitability: #pragma omp parallel for private(m)
	*/
	for (m=0; m<5; m ++ )
	{
		xcrdif[m]=fabs((xcr[m]-xcrref[m])/xcrref[m]);
		xcedif[m]=fabs((xce[m]-xceref[m])/xceref[m]);
	}
	/* --------------------------------------------------------------------- */
	/* Output the comparison of computed results to known cases. */
	/* --------------------------------------------------------------------- */
	if (( * Class)!='U')
	{
		printf(" Verification being performed for class %c\n",  * Class);
		printf(" accuracy setting for epsilon = %20.13E\n", epsilon);
		( * verified)=(fabs(dt-dtref)<=epsilon);
		if ( ! ( * verified))
		{
			( * Class)='U';
			printf(" DT does not match the reference value of %15.8E\n", dtref);
		}
	}
	else
	{
		printf(" Unknown class\n");
	}
	if (( * Class)!='U')
	{
		printf(" Comparison of RMS-norms of residual\n");
	}
	else
	{
		printf(" RMS-norms of residual\n");
	}
	#pragma cetus private(m) 
	#pragma loop name verify#3 
	for (m=0; m<5; m ++ )
	{
		if (( * Class)=='U')
		{
			printf("          %2d%20.13E\n", m+1, xcr[m]);
		}
		else
		{
			if (xcrdif[m]<=epsilon)
			{
				printf("          %2d%20.13E%20.13E%20.13E\n", m+1, xcr[m], xcrref[m], xcrdif[m]);
			}
			else
			{
				( * verified)=false;
				printf(" FAILURE: %2d%20.13E%20.13E%20.13E\n", m+1, xcr[m], xcrref[m], xcrdif[m]);
			}
		}
	}
	if (( * Class)!='U')
	{
		printf(" Comparison of RMS-norms of solution error\n");
	}
	else
	{
		printf(" RMS-norms of solution error\n");
	}
	#pragma cetus private(m) 
	#pragma loop name verify#4 
	for (m=0; m<5; m ++ )
	{
		if (( * Class)=='U')
		{
			printf("          %2d%20.13E\n", m+1, xce[m]);
		}
		else
		{
			if (xcedif[m]<=epsilon)
			{
				printf("          %2d%20.13E%20.13E%20.13E\n", m+1, xce[m], xceref[m], xcedif[m]);
			}
			else
			{
				( * verified)=false;
				printf(" FAILURE: %2d%20.13E%20.13E%20.13E\n", m+1, xce[m], xceref[m], xcedif[m]);
			}
		}
	}
	if (( * Class)=='U')
	{
		printf(" No reference values provided\n");
		printf(" No verification performed\n");
	}
	else
	{
		if ( * verified)
		{
			printf(" Verification Successful\n");
		}
		else
		{
			printf(" Verification failed\n");
		}
	}
}
