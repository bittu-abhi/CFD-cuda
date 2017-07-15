#include <iostream>
#include <fstream>
#include <cuda.h>
#include <cuda_runtime.h>
#include <omp.h>
#include "ausmPlus.h"

using namespace std;

cell::cell(double *state)
{
	omp_set_nested(1);
	omp_set_num_threads(4);
	#pragma omp parallel for
	for(int i=0;i<4;i++)
	{
		stateVar[i]=state[i];
	}
}

cell::cell(){}

__global__ void evaluate(cell *domain,double deltat)
{
	int x=blockIdx.x;
	int y=threadIdx.x;
	double vol=0.5*(domain[x].nodes[0][0]-domain[x].nodes[2][0])*(domain[x].nodes[1][1]-domain[x].nodes[3][1])+\
	0.5*(domain[x].nodes[3][0]-domain[x].nodes[1][0])*(domain[x].nodes[0][1]-domain[x].nodes[2][1]);
	if(domain[x].flag==0 || domain[x].flag==4)
	{
		domain[x].stateVar[y]=domain[x].stateVar[y]-(domain[x].convflux[0][y]+domain[x].convflux[1][y]+domain[x].convflux[2][y]+domain[x].convflux[3][y]\
			-(domain[x].diffflux[0][y]+domain[x].diffflux[1][y]+domain[x].diffflux[2][y]+domain[x].diffflux[3][y]))/vol*deltat;
		if(y==1)
			domain[x].stateVar[y]+=(domain[x].presflux[0][0]+domain[x].presflux[1][0]+domain[x].presflux[2][0]+domain[x].presflux[3][0])/vol*deltat;
		if(y==2)
			domain[x].stateVar[y]+=(domain[x].presflux[0][1]+domain[x].presflux[1][1]+domain[x].presflux[2][1]+domain[x].presflux[3][1])/vol*deltat;
	}
}

__global__ void Boundary(cell *domain,double *initial)
{
	int x=blockIdx.x;
	int y=threadIdx.x;
	//Inlet element evaluation
	if(domain[x].flag==1)
	{
		domain[x].stateVar[y]=initial[y];
	}
	//Outlet element evaluation
	else if(domain[x].flag==2)
	{
		if(threadIdx.x>0)
		{
			domain[x].stateVar[y]=2*domain[x].stateVar[y]-domain[(int)domain[x].face[3][0]].stateVar[y];
		}
		else
			domain[x].stateVar[y]=initial[0];
	}
	//Farfield element evaluation
	else if(domain[x].flag==3)
	{
		/* nothing to change :-) */
	}
}

void ausmplus(double *initial,double timesteps, double deltat)
{
	//OpenMP flags
	omp_set_nested(1);
	omp_set_num_threads(2*omp_get_num_procs());

	//GPU variables
	cell *d_domain;
	double *d_node,*d_boundary,*d_initial;
	double *d_R,*d_k,*d_gammma,*d_mu;

	//Store the values of nodes for faster access
	double *nodes=new double[25452*3];

	//Values of boundary elements
	double *boundary=new double[599*2];

	//Alocate memory to domain
	cell *domain=new cell[25000];

	cout<<"Allocation on the host PC : done"<<endl;
	cout<<endl;

	//Open mesh files
	fstream myfile1,myfile2,myfile3;
	myfile1.open("Elements.txt",ios::in);
	myfile2.open("Nodes.txt",ios::in);
	myfile3.open("boundary.txt",ios::in);

	//Fill the array nodes for access in GPU
	for(int i=0;i<25452*3;i+=3)
	{
		myfile2>>nodes[i]>>nodes[i+1]>>nodes[i+2];
	}

	//Fill the array boundary for access in GPU
	for(int i=0;i<599*2;i+=2)
	{
		myfile3>>boundary[i]>>boundary[i+1];
	}

	//Feed the file just once for access in GPU
	for(int i=0;i<25000;i++)
	{
		domain[i]=cell(initial);
		myfile1>>domain[i].nodes[0][2]>>domain[i].nodes[1][2]>>domain[i].nodes[2][2]>>domain[i].nodes[3][2];
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
	cudaStreamCreateWithFlags(&stream4,cudaStreamNonBlocking);
	cudaMalloc((void **)&d_domain,25000*sizeof(cell));
	cudaMalloc((void **)&d_node,25452*3*sizeof(double));
	cudaMalloc((void **)&d_boundary,599*2*sizeof(double));
	cudaMalloc((void **)&d_initial,4*sizeof(double));
	cudaMalloc((void **)&d_R,sizeof(double));
	cudaMalloc((void **)&d_k,sizeof(double));
	cudaMalloc((void **)&d_gammma,sizeof(double));
	cudaMalloc((void **)&d_mu,sizeof(double));

	cout<<"Allocation on the GPU : done"<<endl;
	cout<<endl;

	cudaMemcpyAsync(d_domain,&domain[0],25000*sizeof(cell),cudaMemcpyHostToDevice,stream1);
	cudaMemcpyAsync(d_node,&nodes[0],25452*3*sizeof(double),cudaMemcpyHostToDevice,stream2);
	cudaMemcpyAsync(d_boundary,&boundary[0],599*2*sizeof(double),cudaMemcpyHostToDevice,stream3);
	cudaMemcpyAsync(d_initial,&initial[0],4*sizeof(double),cudaMemcpyHostToDevice,stream3);
	cudaMemcpyAsync(d_R,&R,sizeof(double),cudaMemcpyHostToDevice,stream3);
	cudaMemcpyAsync(d_k,&k,sizeof(double),cudaMemcpyHostToDevice,stream3);
	cudaMemcpyAsync(d_gammma,&gammma,sizeof(double),cudaMemcpyHostToDevice,stream3);
	cudaMemcpyAsync(d_mu,&mu,sizeof(double),cudaMemcpyHostToDevice,stream3);

	cout<<"Memory copy on the GPU : done"<<endl;
	cout<<endl;

	set_nodes<<<25000,4>>>(d_node,d_domain,d_boundary);
	set_neighbour<<<25000,4>>>(d_domain);
	
	cudaDeviceSynchronize();

	cout<<"Initialisation on the GPU : done"<<endl;
	cout<<endl;

	//Euler first order method
	for (double t = 0; t < timesteps*deltat; t+=deltat)
	{
		pressureFlux<<<25000,4,0,stream1>>>(d_domain,d_R,d_gammma);
		convectiveflux<<<25000,4,0,stream2>>>(d_domain,d_R,d_gammma);
		diffusiveFlux<<<25000,4,0,stream3>>>(d_domain,d_R,d_gammma,d_mu,300,d_k);
		cudaDeviceSynchronize();
		evaluate<<<25000,4>>>(d_domain,deltat);
		Boundary<<<25000,4,0,stream4>>>(d_domain,d_initial);
		cudaDeviceSynchronize();
		cout<<"time = "<<t<<endl;
	}

	cout<<endl;
	cudaMemcpy(&domain[0],d_domain,25000*sizeof(cell),cudaMemcpyDeviceToHost);
	cudaDeviceSynchronize();
	cout<<"Copying final values on the CPU from GPU : done"<<endl;
	cout<<endl;

	visual(domain);
	cudaStreamDestroy(stream1);
	cudaStreamDestroy(stream2);
	cudaStreamDestroy(stream3);	
	cudaStreamDestroy(stream4);	
	cudaFree(d_node);
	cudaFree(d_domain);

	delete[] nodes;
	delete[] boundary;

/*
	//For verification of the correct neighbours and nodes
	for (int i = 0; i < 25000; ++i)
	{
		//cout<<domain[i].diffflux[0][1]<<","<<domain[i].diffflux[1][1]<<","<<domain[i].diffflux[2][1]<<","<<domain[i].diffflux[3][1]<<endl;
		cout<<"("<<domain[i].nodes[0][0]<<","<<domain[i].nodes[0][1]<<")"<<","<<"("<<domain[i].nodes[1][0]<<","<<domain[i].nodes[1][1]<<")";
		cout<<"("<<domain[i].nodes[2][0]<<","<<domain[i].nodes[2][1]<<")"<<","<<"("<<domain[i].nodes[3][0]<<","<<domain[i].nodes[3][1]<<")"<<endl;
		//cout<<domain[i].face[0][0]<<","<<domain[i].face[1][0]<<","<<domain[i].face[2][0]<<","<<domain[i].face[3][0]<<"		"<<domain[i].flag<<endl;
		cout<<endl;
	} 
*/
	delete[] domain;
}
