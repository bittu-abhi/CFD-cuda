#include "ausmPlus.h"
#include <math.h>

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
		double pressplus=(*gammma-1)*(domain[x].stateVar[3]-0.5*(pow(domain[x].stateVar[1],2)+pow(domain[x].stateVar[2],2))/domain[x].stateVar[0]);
		double pressminus=(*gammma-1)*(domain[x].temp_var[y][3]-0.5*(pow(domain[x].temp_var[y][1],2)+pow(domain[x].temp_var[y][2],2))/domain[x].temp_var[y][0]);

		//Machnumber of contravarient velocity(V=u*nx+v*ny)
		double mach_one=(domain[x].stateVar[1]/domain[x].stateVar[0]*domain[x].norms[y][0]+domain[x].stateVar[2]/domain[x].stateVar[0]*domain[x].norms[y][1])/a_mid;
		double mach_two=(domain[x].temp_var[y][1]/domain[x].temp_var[y][0]*domain[x].norms[y][0]+domain[x].temp_var[y][2]/domain[x].temp_var[y][0]*domain[x].norms[y][1])/a_mid;

		double split_mach_one,split_mach_two;
		if(abs(mach_one)>=1)
			split_mach_one=0.5*(mach_one+abs(mach_one));
		else
			split_mach_one=(0.25*pow(mach_one+1,2))+((1.0/8.0)*(pow(pow(mach_one,2.0)-1.0,2.0)));

		if(abs(mach_two)>=1)
			split_mach_two=0.5*(mach_two-abs(mach_two));
		else
			split_mach_two=-0.25*pow(mach_two-1.0,2.0)-1.0/8.0*pow(pow(mach_two,2.0)-1.0,2.0);
	
			double split_mach=split_mach_one+split_mach_two;
		
		for (int i = 0; i < 4; ++i)
		{
			domain[x].convflux[y][i]=a_mid*(0.5*(split_mach+abs(split_mach))*domain[x].stateVar[i]+0.5*(split_mach-abs(split_mach))\
				*domain[x].temp_var[y][i])*sqrt(pow(domain[x].nodes[y][0]-domain[x].nodes[(y+1)%4][0],2)+pow(domain[x].nodes[y][1]-domain[x].nodes[(y+1)%4][1],2));
		}
		domain[x].convflux[y][3]+=a_mid*(0.5*(split_mach+abs(split_mach))*pressplus+0.5*(split_mach-abs(split_mach))*pressminus)\
		*sqrt(pow(domain[x].nodes[y][0]-domain[x].nodes[(y+1)%4][0],2)+pow(domain[x].nodes[y][1]-domain[x].nodes[(y+1)%4][1],2));
	}
}
