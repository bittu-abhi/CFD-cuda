#include "ausmPlus.h"
#include <math.h>

__global__ void convectiveflux(cell *domain, double *R, double *gammma)
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
			a_s[1]=pow(a_s[0],2)/max(a_s[0],abs(sqrt(pow(domain[faces].stateVar[1]/domain[faces].stateVar[0],2)\
			+pow(domain[faces].stateVar[2]/domain[faces].stateVar[0],2))));
		else
			a_s[1]=a_s[0];
		
		//Speed of sound at facial interface
		double a_mid=min(a_s[0],a_s[1]);

		//Pressure Fluxes
		double pressplus=domain[x].stateVar[3]+(*gammma-1)*(domain[x].stateVar[3]-0.5*(pow(domain[x].stateVar[1],2)+pow(domain[x].stateVar[2],2))/domain[x].stateVar[0]);
		double presminus;
		if(domain[x].flag!=4)
			presminus=(*gammma-1)*(domain[faces].stateVar[3]-0.5*(pow(domain[faces].stateVar[1],2)+pow(domain[faces].stateVar[2],2))/domain[faces].stateVar[0]);
		else
			presminus=(gammma[0]-1)*(domain[x].stateVar[3]-0.5*(pow(domain[x].stateVar[1],2)+pow(domain[x].stateVar[2],2))/domain[x].stateVar[0]);

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

		//Machnumber of contravarient velocity
		double machplus=(domain[x].stateVar[1]/domain[x].stateVar[0]*(domain[x].nodes[i1][0]-domain[x].nodes[i2][0])+\
		domain[x].stateVar[2]/domain[x].stateVar[0]*(domain[x].nodes[i2][1]-domain[x].nodes[i1][1]))/a_mid/sqrt(pow(domain[x].nodes[i1][0]-domain[x].nodes[i2][0],2)\
			+pow(domain[x].nodes[i2][1]-domain[x].nodes[i1][1],2));
		double machminus;
		if(domain[x].flag!=4)
		{
			machminus=(domain[faces].stateVar[1]/domain[faces].stateVar[0]*(domain[x].nodes[i1][0]-domain[x].nodes[i2][0])+\
			domain[faces].stateVar[2]/domain[faces].stateVar[0]*(domain[x].nodes[i2][1]-domain[x].nodes[i1][1]))/a_mid/\
			sqrt(pow(domain[x].nodes[i1][0]-domain[x].nodes[i2][0],2)+pow(domain[x].nodes[i2][1]-domain[x].nodes[i1][1],2));
		}
			
		else
		{
			machminus=machplus;
		}
		double split_mach_plus,split_mach_minus;

		if(abs(machplus)>=1)
			split_mach_plus=0.5*(machplus+abs(machplus));
		else
			split_mach_plus=0.5*pow(machplus+1.0,2.0)+1/8*pow(pow(machplus,2.0)-1.0,2.0);
		if(abs(machminus)>=1)
			split_mach_minus=0.5*(machminus-abs(machminus));
		else
			split_mach_minus=-0.5*pow(machminus-1.0,2.0)-1/8*pow(pow(machminus,2.0)-1.0,2.0);

		if(domain[x].flag!=4)
		{
			for (int i = 0; i < 4; ++i)
			{
				domain[x].convflux[y][i]=a_mid*(0.5*(split_mach_plus+abs(split_mach_plus))*domain[x].stateVar[i]+0.5*(split_mach_minus-abs(split_mach_minus))\
					*domain[faces].stateVar[i]);
			}
			domain[x].convflux[y][3]+=a_mid*(0.5*(split_mach_plus+abs(split_mach_plus))*pressplus+0.5*(split_mach_minus-abs(split_mach_minus))*presminus);
		}
		else
		{
			for (int i = 0; i < 4; ++i)
			{
				domain[x].convflux[y][i]=a_mid*(0.5*(split_mach_plus+abs(split_mach_plus))*domain[x].stateVar[i]+0.5*(split_mach_minus-abs(split_mach_minus))\
					*domain[faces].stateVar[i]);
				if(i==1 || i==2)
					domain[x].convflux[y][i]=0;
				
			}
			domain[x].convflux[y][3]+=a_mid*(0.5*(split_mach_plus+abs(split_mach_plus))*pressplus+0.5*(split_mach_minus-abs(split_mach_minus))*presminus);
		}
	}
}