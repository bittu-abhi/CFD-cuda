#include "sbli.h"

point::point()
{}

void sbli(point *pt, int points, double timesteps, double delta)
{
	int x = blockIdx.x;
	int y = blockIdx.y;
	int state = threadIdx.x;

	//CPU variables
	point *ptr = new point[points*points];
	int flag=0;

	//GPU variables
	point *d_pt;
	int *d_flag;

	cudaMalloc((void **)&d_pt,points*points*sizeof(double));
	cudaMalloc((void **)&d_flag,sizeof(int));

	cudaStream_t stream1;
	cudaStreamCreate(&stream1);

	cudaMemcpy(d_pt,&ptr[0],points*points*sizeof(double),cudaMemcpyHostToDevice);
	cudaMemcpy(d_flag,&flag[0],sizeof(int),cudaMemcpyHostToDevice);

	WENO<<<points,points,4>>>(pt,points,0);
	WENO<<<points,points,4>>>(pt,points,1);
	pt[x+y*points]
}
