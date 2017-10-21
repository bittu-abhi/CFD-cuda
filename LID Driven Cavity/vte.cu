#include "lid.h"
#include <stdio.h>
#ifndef CUR
#define CUR curr[y+x*(points)]
#endif

using namespace std;

__global__ void timeW(point *curr, double deltat, double re, int points)
{
	int x = blockIdx.x;
	int y = threadIdx.x;
	/*
	//This part of the code requires some back of the copy formulations hence I am currently using the forward time scheme

	//LDDRK scheme is used
	double k[4];
	//coefficiens of LDDRK scheme
	double c[]={1,0.5,0.162997,0.0407574};
	double beta[]={0,0,0,c[1]};
	beta[2] = c[2]/beta[3];
	beta[1] = c[3]/(beta[3]*beta[2]);

	k[0] = deltat*curr[y+x*sizeof(point)].w;
	k[1] = deltat*

	*/
	CUR.w=CUR.w-deltat*((CUR.delW[0]*CUR.delpsi[1]-CUR.delW[1]*CUR.delpsi[0])-(1/re)*(CUR.delW2[0]+CUR.delW2[1]));
	if(x==points-1 && y == points-1)
		printf("w=%5.14lf psi=%5.14lf delW2X=%5.14lf delW2Y=%5.14lf\n",CUR.w,CUR.psi,CUR.delW2[0],CUR.delW2[1]);
}