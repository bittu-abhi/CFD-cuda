#include "ausmPlus.h"
#include <cmath>
#include <algorithm>

__global__ void pressureFlux(cell *domain, double *R, double *gammma)
{
	int y=threadIdx.x;
	int x=blockIdx.x;
	int faces=(int)domain[x].face[y][0];
	int note;
	if(domain[x].flag==0 || domain[x].flag==4)
	{
		//Calculating the critical speed of sound for all the four sides/faces and the element itself
		double a_s[2];

		if(domain[x].face[y][0]<1 || domain[x].face[y][0]>26000)
			note=y;

		//Element
		a_s[0]=sqrt(2*(gammma[0]-1)/(gammma[0]+1)*(domain[x].stateVar[3]+(gammma[0]-1)*(domain[x].stateVar[3]\
			-0.5*(pow(domain[x].stateVar[1],2)+pow(domain[x].stateVar[2],2))/domain[x].stateVar[0])));
		//Side/face
		if(domain[x].flag!=4)
			a_s[1]=sqrt(2*(gammma[0]-1)/(gammma[0]+1)*(domain[faces].stateVar[3]+(gammma[0]-1)*(domain[faces].stateVar[3]\
			-0.5*(pow(domain[faces].stateVar[1],2)+pow(domain[faces].stateVar[2],2))/domain[faces].stateVar[0])));
		else
		{
			a_s[1]=a_s[0];
		}

		//speed for the boundary calculation
		a_s[0]=pow(a_s[0],2)/max(a_s[0],abs(sqrt(pow(domain[x].stateVar[1]/domain[x].stateVar[0],2)+pow(domain[x].stateVar[2]/domain[x].stateVar[0],2))));
		if(domain[x].flag!=4)
			a_s[1]=pow(a_s[0],2)/max(a_s[0],abs(sqrt(pow(domain[faces].stateVar[1]/domain[faces].stateVar[0],2)+pow(domain[faces].stateVar[2]/domain[faces].stateVar[0],2))));
		else
		{
			a_s[1]=a_s[0];
		}

		int i1,i2;
		if(domain[x].flag!=4)
		{
			for (int i = 0; i < 4; ++i)
			{
				for (int j = 0; j < 4; ++j)
				{
					if(domain[x].nodes[i][0]==domain[faces].nodes[j][0] && domain[x].nodes[i][1]==domain[faces].nodes[j][1])
					{
						if(domain[x].nodes[(i+1)%4][0]==domain[faces].nodes[(j+1)%4][0] && domain[x].nodes[(i+1)%4][1]==domain[faces].nodes[(j+1)%4][1])
						{
							i1=i;
							i2=(i+1)%4;
						}
					}
				}
			}
		}
		else
		{
			for (int i = 0; i < 4; ++i)
			{
				if(domain[x].nodes[i][2]==domain[x].face[note][1] && domain[x].nodes[(i+1)%4][2]==domain[x].face[note][2])
				{
					i1=i;
					i2=(i+1)%4;
				}
			}
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
			machminus=-(domain[faces].stateVar[1]/domain[faces].stateVar[0]*(domain[x].nodes[i1][0]-domain[x].nodes[i2][0])+\
			domain[faces].stateVar[2]/domain[faces].stateVar[0]*(domain[x].nodes[i2][1]-domain[x].nodes[i1][1]))/a_mid/\
			sqrt(pow(domain[x].nodes[i1][0]-domain[x].nodes[i2][0],2)+pow(domain[x].nodes[i2][1]-domain[x].nodes[i1][1],2));
		}
			
		else
		{
			machminus=machplus;
		}
		
		//Pressure Fluxes
		double pressplus=(gammma[0]-1)*(domain[x].stateVar[3]-0.5*(pow(domain[x].stateVar[1],2)+pow(domain[x].stateVar[2],2))/domain[x].stateVar[0]);
		double presminus;
		if(domain[x].flag!=4)
			presminus=(gammma[0]-1)*(domain[faces].stateVar[3]-0.5*(pow(domain[faces].stateVar[1],2)+pow(domain[faces].stateVar[2],2))/domain[faces].stateVar[0]);
		else
		{
			presminus=(gammma[0]-1)*(domain[x].stateVar[3]-0.5*(pow(domain[x].stateVar[1],2)+pow(domain[x].stateVar[2],2))/domain[x].stateVar[0]);
		}
		double plus,minus;
		if(abs(machplus)>=1)
			plus=0.5*(1+machplus/abs(machplus));
		else
			plus=0.25*pow((machplus+1),2)*(2-machplus)+3/16*machplus*pow((pow(machplus,2)-1),2);
		if(abs(machminus)>=1)
			minus=0.5*(1-machminus/abs(machminus));
		else
			minus=0.25*pow((machminus-1),2)*(2+machminus)-3/16*machminus*pow((pow(machminus,2)-1),2);

		domain[x].presflux[y][0]=(pressplus*plus+presminus*minus)*(domain[x].nodes[i1][0]-domain[x].nodes[i2][0]);
		domain[x].presflux[y][1]=(pressplus*plus+presminus*minus)*(domain[x].nodes[i2][1]-domain[x].nodes[i1][1]);
	}
}
