#include "ausmPlus.h"
#include <cmath>
#include <algorithm>

__global__ void pressureFlux(cell *domain, double *R, double *gammma)
{
	int y=threadIdx.x;
	int x=blockIdx.x;
	if(domain[x].flag==0 || domain[x].flag==4)
	{
		//Calculating the critical speed of sound for all the four sides/faces and the element itself
		double a_s[2];
		//Element
		a_s[0]=sqrt(2*(gammma[0]-1)/(gammma[0]+1)*(domain[x].stateVar[3]+(gammma[0]-1)*(domain[x].stateVar[3]\
			-0.5*(pow(domain[x].stateVar[1],2)+pow(domain[x].stateVar[2],2))/domain[x].stateVar[0])));
		//Side/face
		if(domain[x].flag!=4)
			a_s[1]=sqrt(2*(gammma[0]-1)/(gammma[0]+1)*(domain[(int)domain[x].face[y][0]].stateVar[3]+(gammma[0]-1)*(domain[(int)domain[x].face[y][0]].stateVar[3]\
			-0.5*(pow(domain[(int)domain[x].face[y][0]].stateVar[1],2)+pow(domain[(int)domain[x].face[y][0]].stateVar[2],2))/domain[(int)domain[x].face[y][0]].stateVar[0])));
		else
		{
			a_s[1]=a_s[0];
		}

		//speed for the boundary calculation
		a_s[0]=pow(a_s[0],2)/max(a_s[0],abs(sqrt(pow(domain[x].stateVar[1]/domain[x].stateVar[0],2)+pow(domain[x].stateVar[2]/domain[x].stateVar[0],2))));
		if(domain[x].flag!=4)
			a_s[1]=pow(a_s[0],2)/max(a_s[0],abs(sqrt(pow(domain[(int)domain[x].face[y][0]].stateVar[1]/domain[(int)domain[x].face[y][0]].stateVar[0],2)\
			+pow(domain[(int)domain[x].face[y][0]].stateVar[2]/domain[(int)domain[x].face[y][0]].stateVar[0],2))));
		else
		{
			a_s[1]=a_s[0];
		}

		int i1,i2;
		for (int i = 0; i < 4; ++i)
		{
			if(domain[x].nodes[i][2]==domain[x].face[y][0])
				i1=i;
			if(domain[x].nodes[(i+1)%4][2]==domain[x].face[(y+1)%4][0])
				i2=i+1;
		}
		//Speed of sound at facial interface
		double a_mid=min(a_s[0],a_s[1]);
		//Mach number of incoming and outgoing waves
		double machplus=(domain[x].stateVar[1]/domain[x].stateVar[0]*(domain[x].nodes[i1][0]-domain[x].nodes[i2][0])+\
		domain[x].stateVar[2]/domain[x].stateVar[0]*(domain[x].nodes[i2][1]-domain[x].nodes[i1][1]))/a_mid/sqrt(pow(domain[x].nodes[i1][0]-domain[x].nodes[i2][0],2)\
			+pow(domain[x].nodes[i2][1]-domain[x].nodes[i1][1],2));
		double machminus;
		if(domain[x].flag!=4)
		{
			machminus=-(domain[(int)domain[x].face[y][0]].stateVar[1]/domain[(int)domain[x].face[y][0]].stateVar[0]*(domain[x].nodes[i1][0]-domain[x].nodes[i2][0])+\
			domain[(int)domain[x].face[y][0]].stateVar[2]/domain[(int)domain[x].face[y][0]].stateVar[0]*(domain[x].nodes[i2][1]-domain[x].nodes[i1][1]))/a_mid/\
			sqrt(pow(domain[x].nodes[i1][0]-domain[x].nodes[i2][0],2)+pow(domain[x].nodes[i2][1]-domain[x].nodes[i1][1],2));
		}
			
		else
		{
			machminus=machplus;
		}

		//Pressure Fluxes
		double pressplus=domain[x].stateVar[3]+(gammma[0]-1)*(domain[x].stateVar[3]-0.5*(pow(domain[x].stateVar[1],2)+pow(domain[x].stateVar[2],2))/domain[x].stateVar[0]);
		double presminus;
		if(domain[x].flag!=4)
			presminus=domain[(int)domain[x].face[y][0]].stateVar[3]+(gammma[0]-1)*(domain[(int)domain[x].face[y][0]].stateVar[3]\
			-0.5*(pow(domain[(int)domain[x].face[y][0]].stateVar[1],2)+pow(domain[(int)domain[x].face[y][0]].stateVar[2],2))/domain[(int)domain[x].face[y][0]].stateVar[0]);
		else
		{
			presminus=domain[x].stateVar[3]+(gammma[0]-1)*(domain[x].stateVar[3]-0.5*(pow(domain[x].stateVar[1],2)+pow(domain[x].stateVar[2],2))/domain[x].stateVar[0]);
		}

		if(abs(machplus)>=1)
			pressplus*=0.5*(1+machplus/abs(machplus));
		else
			pressplus*=0.25*pow((machplus+1),2)*(2-machplus)+3/16*machplus*pow((pow(machplus,2)-1),2);
		if(abs(machminus)>=1)
			presminus*=0.5*(1-machminus/abs(machminus));
		else
			presminus*=0.25*pow((machminus-1),2)*(2+machminus)-3/16*machminus*pow((pow(machminus,2)-1),2);
		
		domain[x].presflux[y][0]=(pressplus+presminus)*(domain[x].nodes[i1][0]-domain[x].nodes[i2][0]);
		domain[x].presflux[y][1]=(pressplus+presminus)*(domain[x].nodes[i2][1]-domain[x].nodes[i1][1]);
	}
}