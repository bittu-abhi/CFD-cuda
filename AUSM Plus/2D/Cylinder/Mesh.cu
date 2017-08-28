#include "ausmPlus.h"
#include <cuda_runtime_api.h>
#include <math.h>
#include <stdio.h>

__global__ void set_nodes(float *node, cell *domain, float *boundary,float *initial,float *gammma)
{
	int x=blockIdx.x;
	int y=threadIdx.x;
	int temp=3*((int)(domain[x].nodes[y][2])-1);
	domain[x].nodes[y][1]=node[temp+1];
	domain[x].nodes[y][0]=node[temp];
	
	//Inflow
	if(domain[x].nodes[y][2]>500 && domain[x].nodes[y][2]<601 || domain[x].nodes[y][2]==1)
	{
		if((domain[x].nodes[(y+1)%3][2]>500 && domain[x].nodes[(y+1)%3][2]<601 || domain[x].nodes[(y+1)%3][2]==1)||(domain[x].nodes[(y-1+3)%3][2]>500 && domain[x].nodes[(y-1+3)%3][2]<601 || domain[x].nodes[(y+1)%3][2]==1))
		{
			domain[x].flag=1;
			domain[x].stateVar[0]=initial[0];
			domain[x].stateVar[1]=initial[1];
			domain[x].stateVar[2]=initial[2];
			domain[x].stateVar[3]=(101325.000-0.5000*(powf(initial[1],2.0000)+powf(initial[2],2.0000))/initial[0])/(*gammma-1.0000)+\
		0.5000*(powf(initial[1],2.0000)+powf(initial[2],2.0000))/initial[0];	
		}
	}
	
	//Outflow
	if(domain[x].nodes[y][2]>200 && domain[x].nodes[y][2]<302)
	{
		if((domain[x].nodes[(y+1)%3][2]>200 && domain[x].nodes[(y+1)%3][2]<302)|| (domain[x].nodes[(y-1+3)%3][2]>200 && domain[x].nodes[(y-1+3)%3][2]<302))
		{
			domain[x].flag=2;
		}
	}

	//Assumed Farfield
	if((domain[x].nodes[y][2]>301 && domain[x].nodes[y][2]<501) || (domain[x].nodes[y][2]>1 && domain[x].nodes[y][2]<201))
	{
		if(((domain[x].nodes[(y+1)%3][2]>301 && domain[x].nodes[(y+1)%3][2]<501) || (domain[x].nodes[(y+1)%3][2]>1 && domain[x].nodes[(y+1)%3][2]<201))||((domain[x].nodes[(y-1+3)%3][2]>301 && domain[x].nodes[(y-1+3)%3][2]<501) || (domain[x].nodes[(y-1+3)%3][2]>1 && domain[x].nodes[(y-1+3)%3][2]<201)))
		{
			domain[x].flag=2;
			/*domain[x].stateVar[0]=initial[0];
			domain[x].stateVar[1]=initial[1];
			domain[x].stateVar[2]=initial[2];
			domain[x].stateVar[3]=(101325.000-0.5000*(powf(initial[1],2.0000)+powf(initial[2],2.0000))/initial[0])/(*gammma-1.0000)+\
			0.5000*(powf(initial[1],2.0000)+powf(initial[2],2.0000))/initial[0];*/
		}	
	}

	//Cylinder Wall
	if(domain[x].nodes[y][2]>600 && domain[x].nodes[y][2]<2201)
	{
		if((domain[x].nodes[(y+1)%3][2]>600 && domain[x].nodes[(y+1)%3][2]<2201)||(domain[x].nodes[(y-1+3)%3][2]>600 && domain[x].nodes[(y-1+3)%3][2]<2201))
		{
			domain[x].flag=4;
		}
	}
	
	//Fluid Region
	if(domain[x].flag!=1 && domain[x].flag!=2 && domain[x].flag!=3 && domain[x].flag!=4)
		domain[x].flag=0;

	if ((domain[x].nodes[y][2]==301 && domain[x].nodes[(y+1)%3][2]==302) || (domain[x].nodes[y][2]==301 && domain[x].nodes[(y+3-1)%3][2]==302))
		domain[x].flag=2;

	if ((domain[x].nodes[y][2]==200 && domain[x].nodes[(y+1)%3][2]==201) || 
		(domain[x].nodes[y][2]==200 && domain[x].nodes[(y+3-1)%3][2]==201))
		domain[x].flag=2;
}

__global__ void set_neighbour(cell *domain)
{
	int x=blockIdx.x;
	int y=threadIdx.x;
	int i,flag1=0,flag2=0;
	for (i = 0; i < 50266; i+=1)
	{
		for (int j = 0; j <3; j+=1)
		{
			if(domain[i].nodes[j][0]==domain[x].nodes[y][0] && i!=x && domain[i].nodes[j][1]==domain[x].nodes[y][1])
				flag1=1;
			if( domain[i].nodes[j][0]==domain[x].nodes[(y+1)%3][0] && i!=x && domain[i].nodes[j][1]==domain[x].nodes[(y+1)%3][1])
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
	float cen_cord[2];
	cen_cord[0]=(1.0/3.0)*(domain[x].nodes[0][0]+domain[x].nodes[1][0]+domain[x].nodes[2][0]);
	cen_cord[1]=(1.0/3.0)*(domain[x].nodes[0][1]+domain[x].nodes[1][1]+domain[x].nodes[2][1]);

	//construct the face
	float m,c;
	m=(domain[x].nodes[(y+1)%3][1]-domain[x].nodes[y][1])/(domain[x].nodes[(y+1)%3][0]-domain[x].nodes[y][0]);

	c=domain[x].nodes[y][1]-m*domain[x].nodes[y][0];

	//A perpendicular line passing through the centre of the element
	if(m!=0.0000 && !isinf(m))
	{
		float req_m=-1/m;
		float req_c=cen_cord[1]-req_m*cen_cord[0];

		//Intersection of this line with the face would give a point on the face. Now using this point as (x1,y2), we would
		//always get a vector pointing outward from the face,regardless of the way the nodes are number(clockwise or anticlockwise)
		float req_x=(c-req_c)/(req_m-m);
		float req_y=m*req_x+c;
		
		domain[x].norms[y][0]=(req_x-cen_cord[0]);
		domain[x].norms[y][1]=(req_y-cen_cord[1]);

		float dino=sqrt(powf((req_x-cen_cord[0]),2)+powf((req_y-cen_cord[1]),2));
		
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

	if(faces<0 || faces >50266)
	{
		note=y;
	}
	if(y!=note)
	{
		for (int i = 0; i < 4; ++i)
		{
			domain[x].temp_var[y][i]=domain[faces].stateVar[i];
		}
		if(note!=y && domain[x].flag==2)
		{
			domain[x].temp_var[y][0]=1.225;
			domain[x].temp_var[y][1]=domain[(int)domain[x].face[y]-1].stateVar[1];
			domain[x].temp_var[y][2]=domain[(int)domain[x].face[y]-1].stateVar[2];
			domain[x].temp_var[y][3]=domain[(int)domain[x].face[y]-1].stateVar[3];
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
		}
		if(domain[x].flag==2)
		{
			domain[x].temp_var[y][0]=domain[x].stateVar[0];
			domain[x].temp_var[y][1]=domain[x].stateVar[1];
			domain[x].temp_var[y][2]=domain[x].stateVar[2];
			domain[x].temp_var[y][3]=domain[x].stateVar[3];
		}
	}
}
