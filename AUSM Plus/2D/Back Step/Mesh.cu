#include "ausmPlus.h"
#include <cuda_runtime_api.h>
#include <math.h>
#include <stdio.h>

__global__ void set_nodes(double *node, cell *domain, double *boundary, double *initial)
{
	int x=blockIdx.x;
	int y=threadIdx.x;
	int flag1=0;
	int temp=3*((int)(domain[x].nodes[y][2])-1);
	domain[x].nodes[y][1]=node[temp+1];
	domain[x].nodes[y][0]=node[temp];
	
	for(int i=0;i<1040*2;i++)
	{
		if(domain[x].nodes[y][2]==boundary[i])
		{
			flag1=1;
			break;
		}
	}

	if(flag1==0)
		domain[x].flag=0;

	if(domain[x].nodes[y][2]>521 && domain[x].nodes[y][2]<1041)
	{
		domain[x].flag=2;
	}
	
	if((domain[x].nodes[y][2]>66 && domain[x].nodes[y][2]<522) || (domain[x].nodes[y][2]<1041 && domain[x].nodes[y][2]>641) || domain[x].nodes[y][2]==522)
			domain[x].flag=4;

	if(domain[x].nodes[y][2]>=1 && domain[x].nodes[y][2]<67)
	{
		domain[x].flag=1;
		domain[x].stateVar[0]=initial[0];
		domain[x].stateVar[1]=initial[1];
		domain[x].stateVar[2]=initial[2];
		domain[x].stateVar[3]=initial[3];	
	}
}

__global__ void set_neighbour(cell *domain)
{
	int x=blockIdx.x;
	int y=threadIdx.x;
	int i,flag1=0,flag2=0;
	for (i = 0; i < 42830; i++)
	{
		for (int j = 0; j < 4; j+=1)
		{
			if(abs(domain[i].nodes[j][0]-domain[x].nodes[y][0])<0.0001 && i!=x && abs(domain[i].nodes[j][1]-domain[x].nodes[y][1])<0.0001)
				flag1=1;
			if(abs(domain[i].nodes[j][0]-domain[x].nodes[(y+1)%4][0])<0.0001 && i!=x && abs(domain[i].nodes[j][1]-domain[x].nodes[(y+1)%4][1])<0.0001)
				flag2=1;
		}
		if(flag1==1 && flag2==1)
		{	
			domain[x].face[y]=i+1;
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
	cen_cord[0]=0.250000*(domain[x].nodes[0][0]+domain[x].nodes[1][0]+domain[x].nodes[2][0]+domain[x].nodes[3][0]);
	cen_cord[1]=0.250000*(domain[x].nodes[0][1]+domain[x].nodes[1][1]+domain[x].nodes[2][1]+domain[x].nodes[3][1]);

	//construct the face
	double m,c;
	m=(domain[x].nodes[(y+1)%4][1]-domain[x].nodes[y][1])/(domain[x].nodes[(y+1)%4][0]-domain[x].nodes[y][0]);

	c=domain[x].nodes[y][1]-m*domain[x].nodes[y][0];

	//A perpendicular line passing through the centre of the element
	if(m!=0.0000 && !isinf(m))
	{
		double req_m=-1/m;
		double req_c=cen_cord[1]-req_m*cen_cord[0];

		//Intersection of this line with the face would give a point on the face. Now using this point as (x1,y2), we would
		//always get a vector pointing outward from the face,regardless of the way the nodes are number(clockwise or anticlockwise)
		double req_x=(c-req_c)/(req_m-m);
		double req_y=m*req_x+c;
		
		domain[x].norms[y][0]=(req_x-cen_cord[0]);
		domain[x].norms[y][1]=(req_y-cen_cord[1]);

		double dino=sqrt(pow((req_x-cen_cord[0]),2)+pow((req_y-cen_cord[1]),2));
		
		domain[x].norms[y][0]/=dino;
		domain[x].norms[y][1]/=dino;		

	}
	else if(m==0.0000)
	{
		domain[x].norms[y][0]=0;
		if(domain[x].nodes[y][1]<cen_cord[1])
			domain[x].norms[y][1]=-1.000;
		else 
			domain[x].norms[y][1]=1.000;
	}
	else
	{
		domain[x].norms[y][1]=0;
		if(domain[x].nodes[y][0]<cen_cord[0])
			domain[x].norms[y][0]=-1.000;
		else
			domain[x].norms[y][0]=1.0000;
	}
}

__global__ void read_values(cell *domain)
{
	int x=blockIdx.x;
	int y=threadIdx.x;
	int faces=(int)domain[x].face[y]-1;
	int note=-10;

	if(faces<0 || faces >42829)
	{
		note=y;
	}
	if(y!=note)
	{
		for (int i = 0; i < 4; ++i)
		{
			domain[x].temp_var[y][i]=domain[faces].stateVar[i];
		}
	}
	else
	{
		if(domain[x].flag==4)
		{
			domain[x].temp_var[note][0]=1.225;
			domain[x].temp_var[note][1]=-1.0000*domain[x].stateVar[1];
			domain[x].temp_var[note][2]=-1.0000*domain[x].stateVar[2];
			domain[x].temp_var[note][3]=domain[x].stateVar[3];
			//printf("%lf %lf %lf %lf %d %d %d\n",domain[x].stateVar[1],domain[x].temp_var[y][1],domain[x].stateVar[2],domain[x].temp_var[y][2],note, x,y);
		}
	}
}