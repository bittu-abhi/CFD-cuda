#include "lid.h"
#include <iostream>

using namespace std;

__global__ void initial(point *pt , int points)
{
	int x = blockIdx.x;
	int y = threadIdx.x;
	pt[y+x*sizeof(point)].psi = 0;
	pt[y+x*sizeof(point)].w = 0;
	for (int i = 0; i < 2; ++i)
	{
		pt[y+x*sizeof(point)].delW[i]=0;
		pt[y+x*sizeof(point)].delW2[i]=0;
		pt[y+x*sizeof(point)].delpsi[i]=0;
	}
}

void LID(double deltax, double deltay, double deltat,int points, double convergence, int iterations, double reynold, double timesteps)
{
	//CPU variables
	point *ptr = new point[points*points];
	int flag=0;

	//Gpu variables
	point *d_ptr;
	int *d_flag;
	cudaStream_t stream1, stream2,stream3,stream4,stream5,stream6;
	cudaStreamCreate(&stream1);
	cudaStreamCreate(&stream2);
	cudaStreamCreate(&stream3);
	cudaStreamCreate(&stream4);
	cudaStreamCreate(&stream5);
	cudaStreamCreate(&stream6);

	//Memory allocation on GPU
	cudaMalloc((void **)&d_ptr,points*points*sizeof(point));
	cudaMalloc((void **)&d_flag,sizeof(int));
	//Copying memory from CPU to GPU
	cudaMemcpy(d_ptr,&ptr[0],points*points*sizeof(point),cudaMemcpyHostToDevice);
	cudaMemcpy(d_flag,&flag,sizeof(int),cudaMemcpyHostToDevice);
	cout<<"Memory copy on the GPU : done"<<endl;
	cout<<endl;

	initial<<<points,points>>>(d_ptr,points);
	cudaDeviceSynchronize();

	cout<<"initialisation done..."<<endl;
	cout<<endl;

	int iter = 0;

	for (double t = 0; t < timesteps*deltat; t=t+deltat)
	{
		if((int)(t/deltat)%10000==0)
		{
			cudaMemcpyAsync(&ptr[0],d_ptr,points*points*sizeof(point),cudaMemcpyDeviceToHost,stream4);
			visu(ptr,t,deltat,points,deltax,deltay);
		}
		updateboundary<<<points,points,0,stream1>>>(d_ptr,points,deltay,deltax);
		cout<<"vorticity boundary updated"<<endl;
		//Streamfunction evalation
		//Flag set for convergence (1 means error is greater than required threshold)

		//Calculates the total number of iteratios
		iter = 0;
		flag = 0;
		cudaMemcpy(d_flag,&flag,sizeof(int),cudaMemcpyHostToDevice);

		while(flag==0 && iter<iterations)
		{
			streamfunc<<<points,points>>>(deltax,deltay,d_ptr,points,convergence,d_flag);
			cudaDeviceSynchronize();
			iter=iter+1;
			cudaMemcpyAsync(&flag,d_flag,sizeof(int),cudaMemcpyDeviceToHost);
		}
		cout<<endl;
		/*cudaDeviceSynchronize();
		cout<<"hehe"<<endl;
		secondW<<<points,points,1,stream1>>>(d_ptr,points,deltax,0);
		secondW<<<points,points,1,stream2>>>(d_ptr,points,deltay,1);*/
		//1st derivative of w
		compact1D<<<points,1,0,stream1>>>(d_ptr,points,deltax,0,0);
		compact1D<<<points,1,0,stream2>>>(d_ptr,points,deltay,1,0);
		//2nd derivative of w
		compact2D<<<points,1,0,stream3>>>(d_ptr,points,deltax,0);
		compact2D<<<points,1,0,stream4>>>(d_ptr,points,deltay,1);
		//1st derivative of psi
		compact1D<<<points,1,0,stream5>>>(d_ptr,points,deltax,0,1);
		compact1D<<<points,1,0,stream6>>>(d_ptr,points,deltay,1,1);
		cudaDeviceSynchronize();
		timeW<<<points,points>>>(d_ptr,deltat,reynold,points);
		cudaDeviceSynchronize();
		cout<<"t="<<t<<endl;
	}


	cout<<endl;
	cudaMemcpyAsync(&ptr[0],d_ptr,points*points*sizeof(point),cudaMemcpyDeviceToHost);

	cout<<"Copying final values on the CPU from GPU : done"<<endl;
	cout<<endl;
	visu(ptr,timesteps*deltat,deltat,points,deltax,deltay);
	cudaStreamDestroy(stream1);
	cudaStreamDestroy(stream2);
	cudaStreamDestroy(stream3);	
	cudaStreamDestroy(stream4);	
	cudaStreamDestroy(stream5);	
	cudaStreamDestroy(stream6);	

	cudaFree(d_ptr);
	delete[] ptr;
}