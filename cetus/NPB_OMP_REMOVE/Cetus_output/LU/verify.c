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
#include <math.h>
#include "applu.incl"
/* --------------------------------------------------------------------- */
/* verification routine                          */
/* --------------------------------------------------------------------- */
void verify(double xcr[5], double xce[5], double xci, char * Class, logical * verified)
{
	double xcrref[5], xceref[5], xciref;
	double xcrdif[5], xcedif[5], xcidif;
	double epsilon, dtref = 0.0;
	int m;
	/* --------------------------------------------------------------------- */
	/* tolerance level */
	/* --------------------------------------------------------------------- */
	epsilon=1.0E-8;
	( * Class)='U';
	( * verified)=true;
	#pragma cetus private(m) 
	#pragma loop name verify#0 
	#pragma cetus parallel 
	/*
	Disabled due to low profitability: #pragma omp parallel for private(m)
	*/
	for (m=0; m<5; m ++ )
	{
		xcrref[m]=1.0;
		xceref[m]=1.0;
	}
	xciref=1.0;
	if ((((nx0==12)&&(ny0==12))&&(nz0==12))&&(itmax==50))
	{
		( * Class)='S';
		dtref=0.5;
		/* --------------------------------------------------------------------- */
		/* Reference values of RMS-norms of residual, for the (12X12X12) grid, */
		/* after 50 time steps, with  DT = 5.0e-01 */
		/* --------------------------------------------------------------------- */
		xcrref[0]=0.016196343210976703;
		xcrref[1]=0.002197674516482132;
		xcrref[2]=0.0015179927653399185;
		xcrref[3]=0.0015029584435994323;
		xcrref[4]=0.03426407315589646;
		/* --------------------------------------------------------------------- */
		/* Reference values of RMS-norms of solution error,  */
		/* for the (12X12X12) grid, */
		/* after 50 time steps, with  DT = 5.0e-01 */
		/* --------------------------------------------------------------------- */
		xceref[0]=6.422331995796092E-4;
		xceref[1]=8.414434204734792E-5;
		xceref[2]=5.8588269616485187E-5;
		xceref[3]=5.847422259515735E-5;
		xceref[4]=0.0013103347914111294;
		/* --------------------------------------------------------------------- */
		/* Reference value of surface integral, for the (12X12X12) grid, */
		/* after 50 time steps, with DT = 5.0e-01 */
		/* --------------------------------------------------------------------- */
		xciref=7.841892886593708;
	}
	else
	{
		if ((((nx0==33)&&(ny0==33))&&(nz0==33))&&(itmax==300))
		{
			( * Class)='W';
			/* SPEC95fp size */
			dtref=0.0015;
			/* --------------------------------------------------------------------- */
			/* Reference values of RMS-norms of residual, for the (33x33x33) grid, */
			/* after 300 time steps, with  DT = 1.5e-3 */
			/* --------------------------------------------------------------------- */
			xcrref[0]=12.36511638192;
			xcrref[1]=1.317228477799;
			xcrref[2]=2.550120713095;
			xcrref[3]=2.326187750252;
			xcrref[4]=28.26799444189;
			/* --------------------------------------------------------------------- */
			/* Reference values of RMS-norms of solution error,  */
			/* for the (33X33X33) grid, */
			/* --------------------------------------------------------------------- */
			xceref[0]=0.4867877144216;
			xceref[1]=0.05064652880982;
			xceref[2]=0.0928181810196;
			xceref[3]=0.08570126542733;
			xceref[4]=1.084277417792;
			/* --------------------------------------------------------------------- */
			/* Reference value of surface integral, for the (33X33X33) grid, */
			/* after 300 time steps, with  DT = 1.5e-3 */
			/* --------------------------------------------------------------------- */
			xciref=11.61399311023;
		}
		else
		{
			if ((((nx0==64)&&(ny0==64))&&(nz0==64))&&(itmax==250))
			{
				( * Class)='A';
				dtref=2.0;
				/* --------------------------------------------------------------------- */
				/* Reference values of RMS-norms of residual, for the (64X64X64) grid, */
				/* after 250 time steps, with  DT = 2.0e+00 */
				/* --------------------------------------------------------------------- */
				xcrref[0]=779.0210760668937;
				xcrref[1]=63.40276525969287;
				xcrref[2]=194.9924972729248;
				xcrref[3]=178.45301160418538;
				xcrref[4]=1838.4760349464248;
				/* --------------------------------------------------------------------- */
				/* Reference values of RMS-norms of solution error,  */
				/* for the (64X64X64) grid, */
				/* after 250 time steps, with  DT = 2.0e+00 */
				/* --------------------------------------------------------------------- */
				xceref[0]=29.964085685471943;
				xceref[1]=2.819457636500335;
				xceref[2]=7.347341269877474;
				xceref[3]=6.713922568777705;
				xceref[4]=70.71531568839258;
				/* --------------------------------------------------------------------- */
				/* Reference value of surface integral, for the (64X64X64) grid, */
				/* after 250 time steps, with DT = 2.0e+00 */
				/* --------------------------------------------------------------------- */
				xciref=26.030925604886278;
			}
			else
			{
				if ((((nx0==102)&&(ny0==102))&&(nz0==102))&&(itmax==250))
				{
					( * Class)='B';
					dtref=2.0;
					/* --------------------------------------------------------------------- */
					/* Reference values of RMS-norms of residual, for the (102X102X102) grid, */
					/* after 250 time steps, with  DT = 2.0e+00 */
					/* --------------------------------------------------------------------- */
					xcrref[0]=3553.2672969982737;
					xcrref[1]=262.1475079531069;
					xcrref[2]=883.3372185095219;
					xcrref[3]=778.1277473942527;
					xcrref[4]=7308.796959254531;
					/* --------------------------------------------------------------------- */
					/* Reference values of RMS-norms of solution error, for the (102X102X102)  */
					/* grid, after 250 time steps, with  DT = 2.0e+00 */
					/* --------------------------------------------------------------------- */
					xceref[0]=114.01176380212709;
					xceref[1]=8.109896365542157;
					xceref[2]=28.480597317698308;
					xceref[3]=25.90539456783294;
					xceref[4]=260.54907504857414;
					/* --------------------------------------------------------------------- */
					/* Reference value of surface integral, for the (102X102X102) grid, */
					/* after 250 time steps, with DT = 2.0e+00 */
					/* --------------------------------------------------------------------- */
					xciref=47.88716270330823;
				}
				else
				{
					if ((((nx0==162)&&(ny0==162))&&(nz0==162))&&(itmax==250))
					{
						( * Class)='C';
						dtref=2.0;
						/* --------------------------------------------------------------------- */
						/* Reference values of RMS-norms of residual, for the (162X162X162) grid, */
						/* after 250 time steps, with  DT = 2.0e+00 */
						/* --------------------------------------------------------------------- */
						xcrref[0]=10376.698032353785;
						xcrref[1]=892.2124588010086;
						xcrref[2]=2562.3881458266087;
						xcrref[3]=2191.9434385783143;
						xcrref[4]=17807.80572610612;
						/* --------------------------------------------------------------------- */
						/* Reference values of RMS-norms of solution error, for the (162X162X162)  */
						/* grid, after 250 time steps, with  DT = 2.0e+00 */
						/* --------------------------------------------------------------------- */
						xceref[0]=215.98639971694928;
						xceref[1]=15.57895592398636;
						xceref[2]=54.13188630772078;
						xceref[3]=48.22626431540454;
						xceref[4]=455.90291004325036;
						/* --------------------------------------------------------------------- */
						/* Reference value of surface integral, for the (162X162X162) grid, */
						/* after 250 time steps, with DT = 2.0e+00 */
						/* --------------------------------------------------------------------- */
						xciref=66.64045535721813;
						/* --------------------------------------------------------------------- */
						/* Reference value of surface integral, for the (162X162X162) grid, */
						/* after 250 time steps, with DT = 2.0e+00 */
						/* --------------------------------------------------------------------- */
						xciref=66.64045535721813;
					}
					else
					{
						if ((((nx0==408)&&(ny0==408))&&(nz0==408))&&(itmax==300))
						{
							( * Class)='D';
							dtref=1.0;
							/* --------------------------------------------------------------------- */
							/* Reference values of RMS-norms of residual, for the (408X408X408) grid, */
							/* after 300 time steps, with  DT = 1.0e+00 */
							/* --------------------------------------------------------------------- */
							xcrref[0]=48684.17937025;
							xcrref[1]=4696.371050071;
							xcrref[2]=12181.14549776;
							xcrref[3]=10338.01493461;
							xcrref[4]=71423.98413817;
							/* --------------------------------------------------------------------- */
							/* Reference values of RMS-norms of solution error, for the (408X408X408)  */
							/* grid, after 300 time steps, with  DT = 1.0e+00 */
							/* --------------------------------------------------------------------- */
							xceref[0]=375.2393004482;
							xceref[1]=30.84128893659;
							xceref[2]=94.34276905469;
							xceref[3]=82.30686681928;
							xceref[4]=700.262063621;
							/* --------------------------------------------------------------------- */
							/* Reference value of surface integral, for the (408X408X408) grid, */
							/* after 300 time steps, with DT = 1.0e+00 */
							/* --------------------------------------------------------------------- */
							xciref=83.34101392503;
						}
						else
						{
							if ((((nx0==1020)&&(ny0==1020))&&(nz0==1020))&&(itmax==300))
							{
								( * Class)='E';
								dtref=0.5;
								/* --------------------------------------------------------------------- */
								/* Reference values of RMS-norms of residual,  */
								/* for the (1020X1020X1020) grid, */
								/* after 300 time steps, with  DT = 0.5e+00 */
								/* --------------------------------------------------------------------- */
								xcrref[0]=209964.1687874;
								xcrref[1]=21304.03143165;
								xcrref[2]=53192.28789371;
								xcrref[3]=45097.61639833;
								xcrref[4]=293236.000659;
								/* --------------------------------------------------------------------- */
								/* Reference values of RMS-norms of solution error,  */
								/* for the (1020X1020X1020)  */
								/* grid, after 300 time steps, with  DT = 0.5e+00 */
								/* --------------------------------------------------------------------- */
								xceref[0]=480.0572578333;
								xceref[1]=42.21993400184;
								xceref[2]=121.0851906824;
								xceref[3]=104.788898677;
								xceref[4]=836.3028257389;
								/* --------------------------------------------------------------------- */
								/* Reference value of surface integral, for the (1020X1020X1020) grid, */
								/* after 300 time steps, with DT = 0.5e+00 */
								/* --------------------------------------------------------------------- */
								xciref=95.12163272273;
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
	#pragma loop name verify#1 
	#pragma cetus parallel 
	/*
	Disabled due to low profitability: #pragma omp parallel for private(m)
	*/
	for (m=0; m<5; m ++ )
	{
		xcrdif[m]=fabs((xcr[m]-xcrref[m])/xcrref[m]);
		xcedif[m]=fabs((xce[m]-xceref[m])/xceref[m]);
	}
	xcidif=fabs((xci-xciref)/xciref);
	/* --------------------------------------------------------------------- */
	/* Output the comparison of computed results to known cases. */
	/* --------------------------------------------------------------------- */
	if (( * Class)!='U')
	{
		printf("\n Verification being performed for class %c\n",  * Class);
		printf(" Accuracy setting for epsilon = %20.13E\n", epsilon);
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
	#pragma loop name verify#2 
	for (m=0; m<5; m ++ )
	{
		if (( * Class)=='U')
		{
			printf("          %2d  %20.13E\n", m+1, xcr[m]);
		}
		else
		{
			if (xcrdif[m]<=epsilon)
			{
				printf("          %2d  %20.13E%20.13E%20.13E\n", m+1, xcr[m], xcrref[m], xcrdif[m]);
			}
			else
			{
				( * verified)=false;
				printf(" FAILURE: %2d  %20.13E%20.13E%20.13E\n", m+1, xcr[m], xcrref[m], xcrdif[m]);
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
	#pragma loop name verify#3 
	for (m=0; m<5; m ++ )
	{
		if (( * Class)=='U')
		{
			printf("          %2d  %20.13E\n", m+1, xce[m]);
		}
		else
		{
			if (xcedif[m]<=epsilon)
			{
				printf("          %2d  %20.13E%20.13E%20.13E\n", m+1, xce[m], xceref[m], xcedif[m]);
			}
			else
			{
				( * verified)=false;
				printf(" FAILURE: %2d  %20.13E%20.13E%20.13E\n", m+1, xce[m], xceref[m], xcedif[m]);
			}
		}
	}
	if (( * Class)!='U')
	{
		printf(" Comparison of surface integral\n");
	}
	else
	{
		printf(" Surface integral\n");
	}
	if (( * Class)=='U')
	{
		printf("              %20.13E\n", xci);
	}
	else
	{
		if (xcidif<=epsilon)
		{
			printf("              %20.13E%20.13E%20.13E\n", xci, xciref, xcidif);
		}
		else
		{
			( * verified)=false;
			printf(" FAILURE:     %20.13E%20.13E%20.13E\n", xci, xciref, xcidif);
		}
	}
	if (( * Class)=='U')
	{
		printf(" No reference values provided\n");
		printf("No verification performed\n");
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
