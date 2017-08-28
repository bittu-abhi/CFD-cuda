#include <iostream>
#include <fstream>
#include <cuda.h>
#include <cuda_runtime.h>
#include <omp.h>
#include "ausmPlus.h"

using namespace std;

cell::cell(float *initial)
{
	//omp_set_nested(1);
	//omp_set_num_threads(4);
	//#pragma omp parallel for
	stateVar[0]=initial[0];
	stateVar[1]=0.0;
	stateVar[2]=0.0;
	stateVar[3]=(gammma-1)*(initial[3]-0.5*(powf(stateVar[1],2)+powf(stateVar[2],2))/stateVar[0]);
}

cell::cell(){}

__global__ void evaluate(cell *domain,float deltat)
{
	int x=blockIdx.x;
	int y=threadIdx.x;
	float vol=0.5*abs((domain[x].nodes[0][0]-domain[x].nodes[1][0])*(domain[x].nodes[0][1]+domain[x].nodes[1][1])+\
		(domain[x].nodes[1][0]-domain[x].nodes[2][0])*(domain[x].nodes[1][1]+domain[x].nodes[2][1])+\
		(domain[x].nodes[2][0]-domain[x].nodes[0][0])*(domain[x].nodes[2][1]+domain[x].nodes[0][1]));
	if(domain[x].flag==0 || domain[x].flag==4 || domain[x].flag==2)
	{	
		
		domain[x].stateVar[y]-=deltat*(domain[x].convflux[0][y]+domain[x].convflux[1][y]+domain[x].convflux[2][y]\
		-(domain[x].diffflux[0][y]+domain[x].diffflux[1][y]+domain[x].diffflux[2][y]))/vol;

		if(y==1)
			domain[x].stateVar[y]-=deltat*(domain[x].presflux[0][0]+domain[x].presflux[1][0]+domain[x].presflux[2][0])/vol;
		if(y==2)
		domain[x].stateVar[y]-=deltat*(domain[x].presflux[0][1]+domain[x].presflux[1][1]+domain[x].presflux[2][1])/vol;
	}
}

__global__ void Boundary(cell *domain,float *initial)
{
	int x=blockIdx.x;
	int y=threadIdx.x;
	int note=-8;
	int next=-8;

	//Outlet element evaluation
	if(domain[x].flag==1)
	{
		domain[x].stateVar[y]=initial[y];
	}		
}

void ausmplus(float *initial,float timesteps, float deltat)
{
	//GPU variables
	cell *d_domain;
	float *d_node,*d_boundary,*d_initial;
	float *d_R,*d_k,*d_gammma,*d_mu;

	//Store the values of nodes for faster access
	float *nodes=new float[26233*3];

	//Values of boundary elements
	float *boundary=new float[2200*2];

	//Alocate memory to domain
	cell *domain=new cell[50266];

	cout<<"Allocation on the host PC : done"<<endl;
	cout<<endl;

	//Open mesh files
	fstream myfile1,myfile2,myfile3;
	myfile1.open("Elements.txt",ios::in);
	myfile2.open("Nodes.txt",ios::in);
	myfile3.open("boundary.txt",ios::in);

	//Fill the array nodes for access in GPU
	if(myfile2.is_open())
	{
		for(int i=0;i<26233*3;i+=3)
		{
			myfile2>>nodes[i]>>nodes[i+1]>>nodes[i+2];
		}
	}
	else
	{
		cout<<"Could not open Nodes.txt"<<endl;
	}

	//Fill the array boundary for access in GPU
	if(myfile3.is_open())
	{
		for(int i=0;i<2200*2;i+=2)
		{
			myfile3>>boundary[i]>>boundary[i+1];	
		}
	}
	else
	{
		cout<<"Could not open boundary.txt"<<endl;
	}
	//Feed the file just once for access in GPU
	if(myfile1.is_open())
	{
		for(int i=0;i<50266;i++)
		{	
			domain[i]=cell(initial);
			myfile1>>domain[i].nodes[0][2]>>domain[i].nodes[1][2]>>domain[i].nodes[2][2];
		}
	}
	else
	{
		cout<<"Could not open Elements.txt"<<endl;	
	}
	myfile1.close();
	myfile2.close();
	myfile3.close();

	cout<<"Initialisation : done"<<endl;
	cout<<endl;

	cudaStream_t stream1, stream2,stream3;
	cudaStreamCreate(&stream1);
	cudaStreamCreate(&stream2);
	cudaStreamCreate(&stream3);

	cudaMalloc((void **)&d_domain,50266*sizeof(cell));
	cudaMalloc((void **)&d_node,26233*3*sizeof(float));
	cudaMalloc((void **)&d_boundary,2200*2*sizeof(float));
	cudaMalloc((void **)&d_initial,4*sizeof(float));
	cudaMalloc((void **)&d_R,sizeof(float));
	cudaMalloc((void **)&d_k,sizeof(float));
	cudaMalloc((void **)&d_gammma,sizeof(float));
	cudaMalloc((void **)&d_mu,sizeof(float));

	cout<<"Allocation on the GPU : done"<<endl;
	cout<<endl;

	cudaMemcpyAsync(d_domain,&domain[0],50266*sizeof(cell),cudaMemcpyHostToDevice,stream1);
	cudaMemcpyAsync(d_node,&nodes[0],26233*3*sizeof(float),cudaMemcpyHostToDevice,stream2);
	cudaMemcpyAsync(d_boundary,&boundary[0],2200*2*sizeof(float),cudaMemcpyHostToDevice,stream3);
	cudaMemcpyAsync(d_initial,&initial[0],4*sizeof(float),cudaMemcpyHostToDevice,stream3);
	cudaMemcpyAsync(d_R,&R,sizeof(float),cudaMemcpyHostToDevice,stream3);
	cudaMemcpyAsync(d_k,&k,sizeof(float),cudaMemcpyHostToDevice,stream3);
	cudaMemcpyAsync(d_gammma,&gammma,sizeof(float),cudaMemcpyHostToDevice,stream3);
	cudaMemcpyAsync(d_mu,&mu,sizeof(float),cudaMemcpyHostToDevice,stream3);

	cout<<"Memory copy on the GPU : done"<<endl;
	cout<<endl;

	set_nodes<<<50266,3>>>(d_node,d_domain,d_boundary,d_initial,d_gammma);
	set_neighbour<<<50266,3>>>(d_domain);
	calculate_norm<<<50266,3>>>(d_domain);
	read_values<<<50266,3>>>(d_domain);

	cudaDeviceSynchronize();

	cout<<"Initialisation on the GPU : done"<<endl;
	cout<<endl;

	//Euler first order method

	for (float t = 0; t < timesteps*deltat; t+=deltat)
	{
		pressureFlux<<<50266,3,0,stream1>>>(d_domain,d_R,d_gammma);
		convectiveflux<<<50266,3,0,stream2>>>(d_domain,d_R,d_gammma);
		diffusiveFlux<<<50266,3,0,stream3>>>(d_domain,d_R,d_gammma,d_mu,300,d_k);
		cudaDeviceSynchronize();
		if((int)(t/deltat)%1000000==0)
		{
			cudaMemcpyAsync(&domain[0],d_domain,50266*sizeof(cell),cudaMemcpyDeviceToHost,stream3);
			visual(domain,t);
		}
		evaluate<<<50266,4>>>(d_domain,deltat);
		//Boundary<<<50266,3,0,stream1>>>(d_domain,d_initial);
		cudaDeviceSynchronize();
		cout<<"time = "<<t<<endl;
		cout<<endl;
		read_values<<<50266,3,0,stream1>>>(d_domain);
		cudaDeviceSynchronize();
	}

	cout<<endl;
	cudaMemcpy(&domain[0],d_domain,50266*sizeof(cell),cudaMemcpyDeviceToHost);

	cout<<"Copying final values on the CPU from GPU : done"<<endl;
	cout<<endl;

	visual(domain,deltat*timesteps);
	cudaStreamDestroy(stream1);
	cudaStreamDestroy(stream2);
	cudaStreamDestroy(stream3);	

	cudaFree(d_initial);
	cudaFree(d_boundary);
	cudaFree(d_node);
	cudaFree(d_domain);

	delete[] nodes;
	delete[] boundary;
	delete[] domain;
}
