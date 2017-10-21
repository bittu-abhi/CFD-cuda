#include "lid.h"
#include <stdio.h>

__device__ void TDMA(int points, double *D, int *coef,int y_level,int flagD, int flagXY, point *pt, int flagPSIW)
{
	double mod_c[500];
	double mod_d[500];
	for (int i = 0; i < points; ++i)
	{
		if (i==0)
		{
			mod_c[i] = coef[2]/coef[1];
			mod_d[i]=D[i]/coef[1];
		}
		else if( i!=0 || i!=points-1 || i!=1 || i!=points-2)
		{
			mod_c[i] = coef[8]/(coef[7]-coef[6]*mod_c[i-1]);
			mod_d[i] = (D[i]-coef[6]*mod_d[i-1])/(coef[7]-coef[6]*mod_c[i-1]);
		}
		else if(i == 1 || i == points-2)
		{
			mod_c[i] = coef[5]/(coef[4]-coef[3]*mod_c[i-1]);
			mod_d[i] = (D[i]-coef[3]*mod_d[i-1])/(coef[4]-coef[3]*mod_c[i-1]);	
		}	
	}
	mod_d[points-1]= (D[points-1]-coef[9]*mod_d[points-2])/(coef[10]-coef[9]*mod_c[points-2]);

	for (int i = points-1; i >= 0; i=i-1)
	{
		if (i==points-1)
		{
			if(flagD==0 && flagPSIW == 0)
				pt[(flagXY==0)?(i+points*y_level):(i*points+y_level)].delW[flagXY]=mod_d[i];
			if (flagD==1 && flagPSIW == 0)
				pt[(flagXY==0)?(i+points*y_level):(i*points+y_level)].delW2[flagXY]=mod_d[i];
			if (flagPSIW == 1)
				pt[(flagXY==0)?(i+points*y_level):(i*points+y_level)].delpsi[flagXY]=mod_d[i];
			
		}
		else
		{
			if(flagD==0 && flagPSIW==0)
				pt[(flagXY==0)?(i+points*y_level):(i*points+y_level)].delW[flagXY]=mod_d[i]-mod_c[i]*pt[(flagXY==0)?(i+1+points*y_level):((i+1)*points+y_level)].delW[flagXY];
			if(flagD!=0 && flagPSIW==0)
				pt[(flagXY==0)?(i+points*y_level):(i*points+y_level)].delW2[flagXY]=mod_d[i]-mod_c[i]*pt[(flagXY==0)?(i+1+points*y_level):((i+1)*points+y_level)].delW2[flagXY];
			if (flagPSIW==1)
				pt[(flagXY==0)?(i+points*y_level):(i*points+y_level)].delpsi[flagXY]=mod_d[i]-mod_c[i]*pt[(flagXY==0)?(i+1+points*y_level):((i+1)*points+y_level)].delpsi[flagXY];
		}
	}
}

__global__ void compact1D(point *pt, int points, double delta, int flagXY, int flagPSIW)
{
	int coef[12] = {0,2,4,1,4,1,1,3,1,4,2,0};
	int y_level = blockIdx.x;
	double D[500];
	for (int i = 0; i < points; ++i)
	{
		if(flagPSIW==0)
		{
			if (i==0)
				D[i] = (-5*pt[(flagXY==0)?(i+points*y_level):(i*points+y_level)].w+4*pt[(flagXY==0)?(i+1+points*y_level):((i+1)*points+y_level)].w+pt[(flagXY==0)?(i+2+points*y_level):((i+2)*points+y_level)].w)/delta;
			else if(i==1)
				D[i] = 3/delta*(pt[(flagXY==0)?(i+1+points*y_level):((i+1)*points+y_level)].w-pt[(flagXY==0)?(i-1+points*y_level):((i-1)*points+y_level)].w);
			else if (i== points-1)
				D[i] = -1*(-5*pt[(flagXY==0)?(i+points*y_level):(i*points+y_level)].w+4*pt[(flagXY==0)?(i-1+points*y_level):((i-1)*points+y_level)].w+pt[(flagXY==0)?(i-2+points*y_level):((i-2)*points+y_level)].w)/delta;
			else if (i==points-2)
				D[i] = -3/delta*(pt[(flagXY==0)?(i-1+points*y_level):((i-1)*points+y_level)].w-pt[(flagXY==0)?(i+1+points*y_level):((i+1)*points+y_level)].w);
			else
				D[i] = (-pt[(flagXY==0)?(i-2+points*y_level):((i-2)*points+y_level)].w-28*pt[(flagXY==0)?(i-1+points*y_level):((i-1)*points+y_level)].w+28*pt[(flagXY==0)?(i+1+points*y_level):((i+1)*points+y_level)].w+pt[(flagXY==0)?(i+2+points*y_level):((i+2)*points+y_level)].w)/(12*delta);
		}
		else
		{
			if (i==0)
				D[i] = (-5*pt[(flagXY==0)?(i+points*y_level):(i*points+y_level)].psi+4*pt[(flagXY==0)?(i+1+points*y_level):((i+1)*points+y_level)].psi+pt[(flagXY==0)?(i+2+points*y_level):((i+2)*points+y_level)].psi)/delta;
			else if(i==1)
				D[i] = 3/delta*(pt[(flagXY==0)?(i+1+points*y_level):((i+1)*points+y_level)].psi-pt[(flagXY==0)?(i-1+points*y_level):((i-1)*points+y_level)].psi);
			else if (i== points-1)
				D[i] = -1*(-5*pt[(flagXY==0)?(i+points*y_level):(i*points+y_level)].psi+4*pt[(flagXY==0)?(i-1+points*y_level):((i-1)*points+y_level)].psi+pt[(flagXY==0)?(i-2+points*y_level):((i-2)*points+y_level)].psi)/delta;
			else if (i==points-2)
				D[i] = -3/delta*(pt[(flagXY==0)?(i-1+points*y_level):((i-1)*points+y_level)].psi-pt[(flagXY==0)?(i+1+points*y_level):((i+1)*points+y_level)].psi);
			else
				D[i] = (-pt[(flagXY==0)?(i-2+points*y_level):((i-2)*points+y_level)].psi-28*pt[(flagXY==0)?(i-1+points*y_level):((i-1)*points+y_level)].psi+28*pt[(flagXY==0)?(i+1+points*y_level):((i+1)*points+y_level)].psi+pt[(flagXY==0)?(i+2+points*y_level):((i+2)*points+y_level)].psi)/(12*delta);
		}
		/*if(y_level==points-1 && i==points-1)
		{
			printf("%lf delpsix=%lf delpsiy=%lf delwx=%lf delwy=%lf delw2x=%5.14lf delw2y=%5.14lf \n",D[i],pt[i+y_level*points].delpsi[0],pt[i+y_level*points].delpsi[1],pt[i+y_level*points].delW[0],pt[i+y_level*points].delW[1],pt[i+y_level*points].delW2[0],pt[i+y_level*points].delW2[1]);
		}*/		
		
	}
	TDMA(points,D,coef,y_level,0,flagXY,pt,flagPSIW);
}

__global__ void compact2D(point *pt, int points, double delta, int flagXY)
{
	int coef[12] = {0,1,11,1,10,1,2/11,1,2/11,11,1,0};
	int y_level = blockIdx.x;
	double D[500];
	for (int i = 0; i < points; ++i)
	{
		if (i==0)
			D[i] = (13*pt[(flagXY==0)?(i+points*y_level):(i*points+y_level)].w-27*pt[(flagXY==0)?(i+1+points*y_level):((i+1)*points+y_level)].w+15*pt[(flagXY==0)?(i+2+points*y_level):((i+2)*points+y_level)].w-pt[(flagXY==0)?(i+3+points*y_level):((i+3)*points+y_level)].w)/(delta*delta);
		else if(i==1 || i==points-2)
			D[i] = 12/(delta*delta)*(pt[(flagXY==0)?(i+1+points*y_level):((i+1)*points+y_level)].w-2*pt[(flagXY==0)?(i+points*y_level):(i*points+y_level)].w+pt[(flagXY==0)?(i-1+points*y_level):((i-1)*points+y_level)].w);
		else if (i== points-1)
			D[i] = 1 /(delta*delta)*(13*pt[(flagXY==0)?(i+points*y_level):(i*points+y_level)].w-27*pt[(flagXY==0)?(i-1+points*y_level):((i-1)*points+y_level)].w+15*pt[(flagXY==0)?(i-2+points*y_level):((i-2)*points+y_level)].w-pt[(flagXY==0)?(i-3+points*y_level):((i-3)*points+y_level)].w);
		else
			D[i] = 3/(11*4*delta*delta)*(pt[(flagXY==0)?(i-2+points*y_level):((i-2)*points+y_level)].w-2*pt[(flagXY==0)?(i+points*y_level):(i*points+y_level)].w+pt[(flagXY==0)?(i+2+points*y_level):((i+2)*points+y_level)].w)+12/(11*delta*delta)*(pt[(flagXY==0)?(i-1+points*y_level):((i-1)*points+y_level)].w-2*pt[(flagXY==0)?(i+points*y_level):(i*points+y_level)].w+pt[(flagXY==0)?(i+1+points*y_level):((i+1)*points+y_level)].w);
	}
	TDMA(points,D,coef,y_level,1,flagXY,pt,0);
}