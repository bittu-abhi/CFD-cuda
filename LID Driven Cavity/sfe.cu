#include "lid.h"
#include <stdio.h>


__global__ void streamfunc(double deltax, double deltay, point *pt,int points, double convergence, int *address)
{
	int i = threadIdx.x;
	int x = blockIdx.x;
	double temp = pt[i+x*points].psi;
	if(x!=0 && x!=points-1)
	{
		if(i!=0 && i!=points-1)
		{
			pt[i+x*points].psi=1/(2/(deltax*deltax)+2/(deltay*deltay))*(pt[i+x*points].w+(pt[i+1+x*points].psi+pt[i-1+x*points].psi)/(deltax*deltax)+(pt[i+(x+1)*points].psi+pt[i+(x-1)*points].psi)/(deltay*deltay));
		}
	}
	else if (x==0 || x==points-1 || i==0 || i==points-1)
	{
			pt[i+x*points].psi=0;
	}
	pt[i+x*points].error = pt[i+x*points].psi-temp;
	if((pt[i+x*points].psi-temp)<convergence)
		*address = 1;
	//if(i==points-1 && x==points-1)
	//		printf("w=%5.14lf act=%5.14lf error=%5.14lf %d %d\n",pt[i+x*points].w,pt[i+x*points].psi,pt[i+x*points].error, x, i );

}

__global__ void updateboundary(point *pt, int points, double deltay, double deltax)
{
	int x = threadIdx.x;
	int y = blockIdx.x;	
	if(x==0)	
		pt[x+y*points].w=2*(pt[x+y*points].psi-pt[x+1+y*points].psi)/(deltax*deltax);
	else if(x==points-1)
		pt[x+y*points].w=2*(pt[x+y*points].psi-pt[x-1+y*points].psi)/(deltax*deltax);
	if(y==0)
		pt[x+y*points].w=2*(pt[x+y*points].psi-pt[x+(y+1)*points].psi)/(deltay*deltay);
	else if(y==points-1)
		pt[x+y*points].w= -2/deltay+2*(pt[x+y*points].psi-pt[x+(y-1)*points].psi)/(deltay*deltay);
}	
