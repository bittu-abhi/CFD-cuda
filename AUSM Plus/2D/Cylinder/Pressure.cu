#include "ausmPlus.h"
#include <math.h>
#include <algorithm>
#include <stdio.h>

__global__ void pressureFlux(cell *domain, double *R, double *gammma)
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
		//Mach number of incoming and outgoing waves
		double machplus=(domain[x].stateVar[1]/domain[x].stateVar[0]*domain[x].norms[y][0]+domain[x].stateVar[2]/domain[x].stateVar[0]*domain[x].norms[y][1])/a_mid;
		double machminus=(domain[x].temp_var[y][1]/domain[x].temp_var[y][0]*domain[x].norms[y][0]+domain[x].temp_var[y][2]/domain[x].temp_var[y][0]*domain[x].norms[y][1])/a_mid;
		
		//Pressure Fluxes
		double pressplus=(gammma[0]-1)*(domain[x].stateVar[3]-0.5*(pow(domain[x].stateVar[1],2)+pow(domain[x].stateVar[2],2))/domain[x].stateVar[0]);
		double presminus=(gammma[0]-1)*(domain[x].temp_var[y][3]-0.5*(pow(domain[x].temp_var[y][1],2)+pow(domain[x].temp_var[y][2],2))/domain[x].temp_var[y][0]);

		double plus,minus;
		if(abs(machplus)>=1)
			plus=0.5*(1+machplus/abs(machplus));
		else
			plus=0.25*pow((machplus+1),2)*(2-machplus)+3/16*machplus*pow((pow(machplus,2)-1),2);
		if(abs(machminus)>=1)
			minus=0.5*(1-machminus/abs(machminus));
		else
			minus=0.25*pow((machminus-1),2)*(2+machminus)-3/16*machminus*pow((pow(machminus,2)-1),2);

		domain[x].presflux[y][0]=(pressplus*plus+presminus*minus)*domain[x].norms[y][0]*sqrt(pow(domain[x].nodes[(y+1)%4][0]-domain[x].nodes[y][0],2)\
		+pow(domain[x].nodes[(y+1)%4][1]-domain[x].nodes[y][1],2));
		domain[x].presflux[y][1]=(pressplus*plus+presminus*minus)*domain[x].norms[y][1]*sqrt(pow(domain[x].nodes[(y+1)%4][0]-domain[x].nodes[y][0],2)\
		+pow(domain[x].nodes[(y+1)%4][1]-domain[x].nodes[y][1],2));

		//if(domain[x].flag==4)
		//	printf("%lf %lf %d %d\n",domain[x].presflux[y][0],domain[x].presflux[y][1],x,y );
	}
}
