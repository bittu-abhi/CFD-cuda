#include <iostream>
#include <fstream>
#include <cuda.h>
#include <cuda_runtime.h>
#include <omp.h>
#include "ausmPlus.h"
#include <stdio.h>
#include <algorithm>

using namespace std;

cell::cell(double *initial)
{
	//omp_set_nested(1);
	//omp_set_num_threads(4);
	//#pragma omp parallel for
	stateVar[0]=initial[0];
	stateVar[1]=0.0000;
	stateVar[2]=0.0000;
	stateVar[3]=101325.0000/(gammma-1.0000);	
}

cell::cell()
{}

__global__ void evaluate(cell *domain,double deltat)
{
	int x=blockIdx.x;
	int y=threadIdx.x;
	double vol=abs(0.500f*(domain[x].nodes[0][0]-domain[x].nodes[2][0])*(domain[x].nodes[1][1]-domain[x].nodes[3][1])+\
	0.500f*(domain[x].nodes[3][0]-domain[x].nodes[1][0])*(domain[x].nodes[0][1]-domain[x].nodes[2][1]));
	if(domain[x].flag==0 || domain[x].flag==4)
	{	
		
		domain[x].stateVar[y]-=deltat*(domain[x].convflux[0][y]+domain[x].convflux[1][y]+domain[x].convflux[2][y]+domain[x].convflux[3][y]\
		-(domain[x].diffflux[0][y]+domain[x].diffflux[1][y]+domain[x].diffflux[2][y]+domain[x].diffflux[3][y]))/vol;

		if(y==1)
			domain[x].stateVar[y]-=deltat*(domain[x].presflux[0][0]+domain[x].presflux[1][0]+domain[x].presflux[2][0]+domain[x].presflux[3][0])/vol;
		if(y==2)
		domain[x].stateVar[y]-=deltat*(domain[x].presflux[0][1]+domain[x].presflux[1][1]+domain[x].presflux[2][1]+domain[x].presflux[3][1])/vol;
	}
}

__global__ void Boundary(cell *domain,double *initial)
{
	int x=blockIdx.x;
	int y=threadIdx.x;
	int note=-8;
	int neigh1=-9,next=-8;

	//Outlet element evaluation
	if(domain[x].flag==1)
	{
		domain[x].stateVar[y]=initial[y];
	}
	if(domain[x].flag==2)
	{
		if(domain[x].face[y]<1 || domain[x].face[y]>42830)
		{
			note=y;
		}
			next=(int)domain[x].face[(y+2)%4]-1;
		if(y==note)
		{
			for (int i = 0; i < 4; ++i)
			{
				if(domain[x].nodes[(y+2)%4][0]==domain[next].nodes[i][0] && domain[x].nodes[(y+2)%4][1]==domain[next].nodes[i][1] )
				{
					if(domain[x].nodes[(y+3)%4][0]==domain[next].nodes[(i+1)%4][0] && domain[x].nodes[(y+3)%4][1]==domain[next].nodes[(i+1)%4][1] )
					{
						neigh1=i;
					}
					else if(domain[x].nodes[(y+3)%4][0]==domain[next].nodes[(i+4-1)%4][0] && domain[x].nodes[(y+3)%4][1]==domain[next].nodes[(i+4-1)%4][1])
					{
						neigh1=(i-1+4)%4;
					}
				}
			}
			domain[x].stateVar[0]=initial[0];
			domain[x].stateVar[1]=domain[next].stateVar[1];//-domain[next].temp_var[(neigh1+2)%4][1];
			domain[x].stateVar[2]=domain[next].stateVar[2];//-domain[next].temp_var[(neigh1+2)%4][2];
			domain[x].stateVar[3]=domain[next].stateVar[3];//-domain[next].temp_var[(neigh1+2)%4][3];
		}
	}
}

void ausmplus(double *initial,double timesteps, double deltat)
{
	//GPU variables
	cell *d_domain;
	double *d_node,*d_boundary,*d_initial;
	double *d_R,*d_k,*d_gammma,*d_mu;

	//Store the values of nodes for faster access
	double *nodes=new double[43351*3];

	//Values of boundary elements
	double *boundary=new double[1040*2];

	//Alocate memory to domain
	cell *domain=new cell[42830];

	cout<<"Allocation on the host PC : done"<<endl;
	cout<<endl;

	//Open mesh files
	fstream myfile1,myfile2,myfile3;
	myfile1.open("Elements.txt",ios::in);
	myfile2.open("Nodes.txt",ios::in);
	myfile3.open("boundary.txt",ios::in);

	double initial1[]={initial[0],initial[1],initial[2],initial[3]};

	//Fill the array nodes for access in GPU
	if(myfile2.is_open())
	{
		for(int i=0;i<43351*3;i+=3)
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
		for(int i=0;i<1040*2;i+=2)
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
		for(int i=0;i<42830;i++)
		{	
			domain[i]=cell(initial);
			myfile1>>domain[i].nodes[0][2]>>domain[i].nodes[1][2]>>domain[i].nodes[2][2]>>domain[i].nodes[3][2];
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

	cudaStream_t stream1, stream2,stream3,stream4;
	cudaStreamCreate(&stream1);
	cudaStreamCreate(&stream2);
	cudaStreamCreate(&stream3);
	cudaStreamCreate(&stream4);

	cudaMalloc((void **)&d_domain,42830*sizeof(cell));
	cudaMalloc((void **)&d_node,43351*3*sizeof(double));
	cudaMalloc((void **)&d_boundary,1040*2*sizeof(double));
	cudaMalloc((void **)&d_initial,4*sizeof(double));
	cudaMalloc((void **)&d_R,sizeof(double));
	cudaMalloc((void **)&d_k,sizeof(double));
	cudaMalloc((void **)&d_gammma,sizeof(double));
	cudaMalloc((void **)&d_mu,sizeof(double));

	cout<<"Allocation on the GPU : done"<<endl;
	cout<<endl;

	cudaMemcpyAsync(d_domain,&domain[0],42830*sizeof(cell),cudaMemcpyHostToDevice,stream1);
	cudaMemcpyAsync(d_node,&nodes[0],43351*3*sizeof(double),cudaMemcpyHostToDevice,stream2);
	cudaMemcpyAsync(d_boundary,&boundary[0],1040*2*sizeof(double),cudaMemcpyHostToDevice,stream3);
	cudaMemcpyAsync(d_initial,&initial1[0],4*sizeof(double),cudaMemcpyHostToDevice,stream3);
	cudaMemcpyAsync(d_R,&R,sizeof(double),cudaMemcpyHostToDevice,stream3);
	cudaMemcpyAsync(d_k,&k,sizeof(double),cudaMemcpyHostToDevice,stream3);
	cudaMemcpyAsync(d_gammma,&gammma,sizeof(double),cudaMemcpyHostToDevice,stream3);
	cudaMemcpyAsync(d_mu,&mu,sizeof(double),cudaMemcpyHostToDevice,stream3);

	cout<<"Memory copy on the GPU : done"<<endl;
	cout<<endl;

	set_nodes<<<42830,4>>>(d_node,d_domain,d_boundary,d_initial);
	set_neighbour<<<42830,4>>>(d_domain);
	calculate_norm<<<42830,4>>>(d_domain);
	read_values<<<42830,4>>>(d_domain);

	cudaDeviceSynchronize();

	cout<<"Initialisation on the GPU : done"<<endl;
	cout<<endl;

	//Euler first order method
	for (double t = 0.0; t < timesteps*deltat; t+=deltat)
	{
		pressureFlux<<<42830,4,0,stream1>>>(d_domain,d_R,d_gammma);
		convectiveflux<<<42830,4,0,stream2>>>(d_domain,d_R,d_gammma);
		diffusiveFlux<<<42830,4,0,stream3>>>(d_domain,d_R,d_gammma,d_mu,300,d_k);
		cudaDeviceSynchronize();
		evaluate<<<42830,4>>>(d_domain,deltat);
		Boundary<<<42830,4>>>(d_domain,d_initial);
		cout<<endl;
		read_values<<<42830,4>>>(d_domain);
		cudaDeviceSynchronize();
		cout<<"time = "<<t<<endl;
		if((int)(t/deltat)%1000==0)
		{
			cudaMemcpyAsync(&domain[0],d_domain,42830*sizeof(cell),cudaMemcpyDeviceToHost,stream4);
			visual(domain,t);
		}	
	}

	cout<<endl;
	cudaMemcpy(&domain[0],d_domain,42830*sizeof(cell),cudaMemcpyDeviceToHost);

	cout<<"Copying final values on the CPU from GPU : done"<<endl;
	cout<<endl;

	visual(domain,deltat*timesteps);
	cudaStreamDestroy(stream1);
	cudaStreamDestroy(stream2);
	cudaStreamDestroy(stream3);	
	cudaStreamDestroy(stream4);	

	cudaFree(d_initial);
	cudaFree(d_boundary);
	cudaFree(d_node);
	cudaFree(d_domain);

	delete[] nodes;
	delete[] boundary;
	delete[] domain;
}