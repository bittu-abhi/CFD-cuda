#include "ausmPlus.h"
#include <math.h>
#include <stdio.h>

__global__ void convectiveflux(cell *domain, double *R, double *gammma)
{
	int y=threadIdx.x;
	int x=blockIdx.x;
	int ourFlag=(int)domain[x].flag;
	if(ourFlag==0 || ourFlag==4)
	{
		//Calculating the critical speed of sound for all the four sides/faces and the element itself
		double a_s[2];

		//Element
		a_s[0]=sqrt(2*(gammma[0]-1)/(gammma[0]+1)*(domain[x].stateVar[3]+(gammma[0]-1)*(domain[x].stateVar[3]\
			-0.5*(pow(domain[x].stateVar[1],2)+pow(domain[x].stateVar[2],2))/domain[x].stateVar[0]))/domain[x].stateVar[0]);
		//Side/face
		a_s[1]=sqrt(2*(gammma[0]-1)/(gammma[0]+1)*(domain[x].temp_var[y][3]+(gammma[0]-1)*(domain[x].temp_var[y][3]\
			-0.5*(pow(domain[x].temp_var[y][1],2)+pow(domain[x].temp_var[y][2],2))/domain[x].temp_var[y][0]))/domain[x].temp_var[y][0]);

		//speed for the boundary calculation
		a_s[0]=pow(a_s[0],2)/max(a_s[0],abs(sqrt(pow(domain[x].stateVar[1]/domain[x].stateVar[0],2)+pow(domain[x].stateVar[2]/domain[x].stateVar[0],2))));
		a_s[1]=pow(a_s[1],2)/max(a_s[1],abs(sqrt(pow(domain[x].temp_var[y][1]/domain[x].temp_var[y][0],2)+pow(domain[x].temp_var[y][2]/domain[x].temp_var[y][0],2))));

		//Speed of sound at facial interface
		double a_mid=min(a_s[0],a_s[1]);

		//Pressure
		double press=domain[x].stateVar[3]+(*gammma-1)*(domain[x].stateVar[3]-0.5*(pow(domain[x].stateVar[1],2)+pow(domain[x].stateVar[2],2))/domain[x].stateVar[0]);

		//Machnumber of contravarient velocity(V=u*nx+v*ny)
		double machplus=(domain[x].stateVar[1]/domain[x].stateVar[0]*domain[x].norms[y][0]+domain[x].stateVar[2]/domain[x].stateVar[0]*domain[x].norms[y][1])/a_mid;
		double machminus=(domain[x].temp_var[y][1]/domain[x].temp_var[y][0]*domain[x].norms[y][0]+domain[x].temp_var[y][2]/domain[x].temp_var[y][0]*domain[x].norms[y][1])/a_mid;

		double split_mach_plus,split_mach_minus;

		if(abs(machplus)>=1)
			split_mach_plus=0.5*(machplus+abs(machplus));
		else
			split_mach_plus=0.5*pow(machplus+1.0,2.0)+1/8*pow(pow(machplus,2.0)-1.0,2.0);
		if(abs(machminus)>=1)
			split_mach_minus=0.5*(machminus-abs(machminus));
		else
			split_mach_minus=-0.5*pow(machminus-1.0,2.0)-1/8*pow(pow(machminus,2.0)-1.0,2.0);

		for (int i = 0; i < 4; ++i)
		{
			domain[x].convflux[y][i]=a_mid*(0.5*(split_mach_plus+abs(split_mach_plus))*domain[x].stateVar[i]+0.5*(split_mach_minus-abs(split_mach_minus))\
				*domain[x].temp_var[y][i])*sqrt(pow(domain[x].nodes[y][0]-domain[x].nodes[(y+1)%4][0],2)+pow(domain[x].nodes[y][1]-domain[x].nodes[(y+1)%4][1],2));
		}
		domain[x].convflux[y][3]+=a_mid*(0.5*(split_mach_plus+abs(split_mach_plus))*press+0.5*(split_mach_minus-abs(split_mach_minus))*press)\
		*sqrt(pow(domain[x].nodes[y][0]-domain[x].nodes[(y+1)%4][0],2)+pow(domain[x].nodes[y][1]-domain[x].nodes[(y+1)%4][1],2));
	}
}
