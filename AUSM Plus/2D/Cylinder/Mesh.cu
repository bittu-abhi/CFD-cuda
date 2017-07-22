#include "ausmPlus.h"
#include <cuda_runtime_api.h>
#include <math.h>

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
	for(int i=0;i<600*2;i++)
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
		if(abs(domain[x].nodes[i][0]-(20+sqrt(4-pow(domain[x].nodes[i][1]-15,2))))<0.005 || abs(domain[x].nodes[i][0]-(20-sqrt(4-pow(domain[x].nodes[i][1]-15,2))))<0.005)
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
			if(domain[i].nodes[j][0]==domain[x].nodes[y][0] && i!=x && domain[i].nodes[j][1]==domain[x].nodes[y][1])
				flag1=1;
			if( domain[i].nodes[j][0]==domain[x].nodes[(y+1)%4][0] && i!=x && domain[i].nodes[j][1]==domain[x].nodes[(y+1)%4][1])
				flag2=1;
		}
		if(flag1==1 && flag2==1)
		{	
			domain[x].face[y]=i;
			break;
		}
		flag1=0;
		flag2=0;
	}
}

__global__ void calculate_norm(cell *domain)
{
	int x=blockIdx.x;
	int y=threadIdx.x;

	//Now to determine if the normal is pointing outward, and if not, then change accordingly
	double cen_cord[2];
	cen_cord[0]=0.25*(domain[x].nodes[0][0]+domain[x].nodes[1][0]+domain[x].nodes[2][0]+domain[x].nodes[3][0]);
	cen_cord[1]=0.25*(domain[x].nodes[0][1]+domain[x].nodes[1][1]+domain[x].nodes[2][1]+domain[x].nodes[3][1]);

	//construct the face
	double m,c;
	//if(domain[x].nodes[(y+1)%4][0]-domain[x].nodes[y][0]!=0)
		m=(domain[x].nodes[(y+1)%4][1]-domain[x].nodes[y][1])/(domain[x].nodes[(y+1)%4][0]-domain[x].nodes[y][0]);

	c=domain[x].nodes[y][1]-m*domain[x].nodes[y][0];

	//A perpendicular line passing through the centre of the element
	if(m!=0 && !isinf(m))
	{
		double req_m=-1/m;
		double req_c=cen_cord[1]-req_m*cen_cord[0];

		//Intersection of this line with the face would give a point on the face. Now using this point as (x1,y2), we would
		//always get a vector pointing outward from the face,regardless of the way the nodes are number(clockwise or anticlockwise)
		double req_x=(c-req_c)/(req_m-m);
		double req_y=m*req_x+c;
		
		double dino=sqrt(pow(req_x-cen_cord[0],2)+pow(req_y-cen_cord[1],2));

		domain[x].norms[y][0]=(req_x-cen_cord[0])/dino;
		domain[x].norms[y][1]=(req_y-cen_cord[1])/dino;
	}
	else if(m==0)
	{
		domain[x].norms[y][0]=0;
		if(domain[x].nodes[(y+1)%4][1]<cen_cord[1])
			domain[x].norms[y][1]=-1;
		else
			domain[x].norms[y][1]=1;
	}
	else
	{
		domain[x].norms[y][1]=0;
		if(domain[x].nodes[(y+1)%4][0]<cen_cord[0])
			domain[x].norms[y][0]=-1;
		else
			domain[x].norms[y][0]=1;
	}
}
