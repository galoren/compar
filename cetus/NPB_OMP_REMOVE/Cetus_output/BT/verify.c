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
	/* reference data for 12X12X12 grids after 60 time steps, with DT = 1.0e-02 */
	/* --------------------------------------------------------------------- */
	if ((((grid_points[0]==12)&&(grid_points[1]==12))&&(grid_points[2]==12))&&(no_time_steps==60))
	{
		( * Class)='S';
		dtref=0.01;
		/* --------------------------------------------------------------------- */
		/* Reference values of RMS-norms of residual. */
		/* --------------------------------------------------------------------- */
		xcrref[0]=0.17034283709541312;
		xcrref[1]=0.012975252070034096;
		xcrref[2]=0.032527926989486054;
		xcrref[3]=0.0264364212751668;
		xcrref[4]=0.1921178413174443;
		/* --------------------------------------------------------------------- */
		/* Reference values of RMS-norms of solution error. */
		/* --------------------------------------------------------------------- */
		xceref[0]=4.997691334581158E-4;
		xceref[1]=4.519566678296193E-5;
		xceref[2]=7.397376517292135E-5;
		xceref[3]=7.382123863243973E-5;
		xceref[4]=8.926963098749145E-4;
		/* --------------------------------------------------------------------- */
		/* reference data for 24X24X24 grids after 200 time steps,  */
		/* with DT = 0.8e-3 */
		/* --------------------------------------------------------------------- */
	}
	else
	{
		if ((((grid_points[0]==24)&&(grid_points[1]==24))&&(grid_points[2]==24))&&(no_time_steps==200))
		{
			( * Class)='W';
			dtref=8.0E-4;
			/* --------------------------------------------------------------------- */
			/* Reference values of RMS-norms of residual. */
			/* --------------------------------------------------------------------- */
			xcrref[0]=112.5590409344;
			xcrref[1]=11.80007595731;
			xcrref[2]=27.10329767846;
			xcrref[3]=24.69174937669;
			xcrref[4]=263.8427874317;
			/* --------------------------------------------------------------------- */
			/* Reference values of RMS-norms of solution error. */
			/* --------------------------------------------------------------------- */
			xceref[0]=4.419655736008;
			xceref[1]=0.4638531260002;
			xceref[2]=1.011551749967;
			xceref[3]=0.9235878729944;
			xceref[4]=10.18045837718;
			/* --------------------------------------------------------------------- */
			/* reference data for 64X64X64 grids after 200 time steps,  */
			/* with DT = 0.8e-3 */
			/* --------------------------------------------------------------------- */
		}
		else
		{
			if ((((grid_points[0]==64)&&(grid_points[1]==64))&&(grid_points[2]==64))&&(no_time_steps==200))
			{
				( * Class)='A';
				dtref=8.0E-4;
				/* --------------------------------------------------------------------- */
				/* Reference values of RMS-norms of residual. */
				/* --------------------------------------------------------------------- */
				xcrref[0]=108.06346714637264;
				xcrref[1]=11.319730901220813;
				xcrref[2]=25.974354511582465;
				xcrref[3]=23.66562254467891;
				xcrref[4]=252.78963211748345;
				/* --------------------------------------------------------------------- */
				/* Reference values of RMS-norms of solution error. */
				/* --------------------------------------------------------------------- */
				xceref[0]=4.2348416040525025;
				xceref[1]=0.443902824969957;
				xceref[2]=0.9669248013634565;
				xceref[3]=0.8830206303976548;
				xceref[4]=9.737990177082928;
				/* --------------------------------------------------------------------- */
				/* reference data for 102X102X102 grids after 200 time steps, */
				/* with DT = 3.0e-04 */
				/* --------------------------------------------------------------------- */
			}
			else
			{
				if ((((grid_points[0]==102)&&(grid_points[1]==102))&&(grid_points[2]==102))&&(no_time_steps==200))
				{
					( * Class)='B';
					dtref=3.0E-4;
					/* --------------------------------------------------------------------- */
					/* Reference values of RMS-norms of residual. */
					/* --------------------------------------------------------------------- */
					xcrref[0]=1423.3597229287254;
					xcrref[1]=99.33052259015024;
					xcrref[2]=356.46025644535285;
					xcrref[3]=324.8544795908409;
					xcrref[4]=3270.7541254659363;
					/* --------------------------------------------------------------------- */
					/* Reference values of RMS-norms of solution error. */
					/* --------------------------------------------------------------------- */
					xceref[0]=52.96984714093686;
					xceref[1]=4.463289611567067;
					xceref[2]=13.122573342210174;
					xceref[3]=12.006925323559145;
					xceref[4]=124.59576151035986;
					/* --------------------------------------------------------------------- */
					/* reference data for 162X162X162 grids after 200 time steps, */
					/* with DT = 1.0e-04 */
					/* --------------------------------------------------------------------- */
				}
				else
				{
					if ((((grid_points[0]==162)&&(grid_points[1]==162))&&(grid_points[2]==162))&&(no_time_steps==200))
					{
						( * Class)='C';
						dtref=1.0E-4;
						/* --------------------------------------------------------------------- */
						/* Reference values of RMS-norms of residual. */
						/* --------------------------------------------------------------------- */
						xcrref[0]=6239.8116551764615;
						xcrref[1]=507.93239190423964;
						xcrref[2]=1542.3530093013596;
						xcrref[3]=1330.238792929119;
						xcrref[4]=11604.087428436455;
						/* --------------------------------------------------------------------- */
						/* Reference values of RMS-norms of solution error. */
						/* --------------------------------------------------------------------- */
						xceref[0]=164.62008369091265;
						xceref[1]=11.497107903824313;
						xceref[2]=41.20744620746151;
						xceref[3]=37.08765105969417;
						xceref[4]=362.11053051841265;
						/* --------------------------------------------------------------------- */
						/* reference data for 408x408x408 grids after 250 time steps, */
						/* with DT = 0.2e-04 */
						/* --------------------------------------------------------------------- */
					}
					else
					{
						if ((((grid_points[0]==408)&&(grid_points[1]==408))&&(grid_points[2]==408))&&(no_time_steps==250))
						{
							( * Class)='D';
							dtref=2.0E-5;
							/* --------------------------------------------------------------------- */
							/* Reference values of RMS-norms of residual. */
							/* --------------------------------------------------------------------- */
							xcrref[0]=25331.88551738;
							xcrref[1]=2346.39371698;
							xcrref[2]=6294.554366904;
							xcrref[3]=5352.56537603;
							xcrref[4]=39058.64038618;
							/* --------------------------------------------------------------------- */
							/* Reference values of RMS-norms of solution error. */
							/* --------------------------------------------------------------------- */
							xceref[0]=310.0009377557;
							xceref[1]=24.24086324913;
							xceref[2]=77.82212022645;
							xceref[3]=68.35623860116;
							xceref[4]=606.5737200368;
							/* --------------------------------------------------------------------- */
							/* reference data for 1020x1020x1020 grids after 250 time steps, */
							/* with DT = 0.4e-05 */
							/* --------------------------------------------------------------------- */
						}
						else
						{
							if ((((grid_points[0]==1020)&&(grid_points[1]==1020))&&(grid_points[2]==1020))&&(no_time_steps==250))
							{
								( * Class)='E';
								dtref=4.0E-6;
								/* --------------------------------------------------------------------- */
								/* Reference values of RMS-norms of residual. */
								/* --------------------------------------------------------------------- */
								xcrref[0]=97953.72484517;
								xcrref[1]=9739.814511521;
								xcrref[2]=24676.06342965;
								xcrref[3]=20924.1957286;
								xcrref[4]=139213.8856939;
								/* --------------------------------------------------------------------- */
								/* Reference values of RMS-norms of solution error. */
								/* --------------------------------------------------------------------- */
								xceref[0]=432.7562208414;
								xceref[1]=36.99051964887;
								xceref[2]=108.9845040954;
								xceref[3]=94.62517622043;
								xceref[4]=776.5512765309;
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
	/* the defined grid sizes above (Class != 'U') */
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
