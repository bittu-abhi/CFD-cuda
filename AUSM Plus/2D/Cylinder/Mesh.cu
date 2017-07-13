#include "ausmPlus.h"

__global__ void set_nodes(double *node, cell *domain, double *boundary)
{
	int x=blockIdx.x;
	int y=threadIdx.x;
	int flag1=0;
	int temp=3*((int)(domain[x].nodes[y][2])-1);
	domain[x].nodes[y][1]=node[temp+1];
	domain[x].nodes[y][0]=node[temp];
	if(domain[x].nodes[y][2]>126 && domain[x].nodes[y][2]<176)
	{
		domain[x].flag=1;	
	}
	else if(domain[x].nodes[y][2]>=25402 && domain[x].nodes[y][2]<=25452)
	{
		domain[x].flag=2;
	}
	for(int i=0;i<599*2;i++)
	{
		if(domain[x].nodes[0][2]==boundary[i] || domain[x].nodes[1][2]==boundary[i] || domain[x].nodes[2][2]==boundary[i] || domain[x].nodes[3][2]==boundary[i])
		{
			flag1=1;
			break;
		}
	}
	if(flag1==0)
		domain[x].flag=0;

	for (int i = 0; i < 4; ++i)
	{
		if(abs(domain[x].nodes[i][0]-(20+sqrt(4-pow(domain[x].nodes[i][1]-15,2))))<0.002 || abs(domain[x].nodes[i][0]-(20-sqrt(4-pow(domain[x].nodes[i][1]-15,2))))<0.002)
			domain[x].flag=4;
	}
	if(flag1==1 && domain[x].flag!=1 && domain[x].flag!=2 && domain[x].flag!=4 && domain[x].flag!=0)
			domain[x].flag=3;
}

__global__ void set_neighbour(cell *domain)
{
	int x=blockIdx.x;
	int y=threadIdx.x;
	int i,flag1=0,flag2=0;
	for (i = 0; i < 25000; i+=1)
	{
		for (int j = 0; j < 4; j+=1)
		{
			if(domain[i].nodes[j][2]==domain[x].nodes[y][2] && i!=x)
				flag1=1;
			if(domain[i].nodes[j][2]==domain[x].nodes[(y+1)%4][2] && i!=x)
				flag2=1;
		}
		if(flag1==1 && flag2==1)
		{	
			domain[x].face[y][0]=i;
			domain[x].face[y][1]=domain[x].nodes[y][2];
			domain[x].face[y][2]=domain[x].nodes[(y+1)%4][2];
			break;
		}
		flag1=0;
		flag2=0;
	}
}