#include "ausmPlus.h"
#include <math.h>
#include <algorithm>

__global__ void pressureFlux(cell *domain, double *R, double *gammma)
{
	int y=threadIdx.x;
	int x=blockIdx.x;
	int faces=(int)domain[x].face[y];
	int note=-1;
	int ourFlag=(int)domain[x].flag;
	if(ourFlag==0 || ourFlag==4)
	{

		//Calculating the critical speed of sound for all the four sides/faces and the element itself
		double a_s[2];
		if(domain[x].face[y]<1)
		{
			note=y;
		}
		
		//Element
		a_s[0]=sqrt(2*(gammma[0]-1)/(gammma[0]+1)*(domain[x].stateVar[3]+(gammma[0]-1)*(domain[x].stateVar[3]\
			-0.5*(pow(domain[x].stateVar[1],2)+pow(domain[x].stateVar[2],2))/domain[x].stateVar[0])));
		//Side/face
		if(ourFlag!=4 || (ourFlag==4 && y!=note))
			a_s[1]=sqrt(2*(gammma[0]-1)/(gammma[0]+1)*(domain[faces].stateVar[3]+(gammma[0]-1)*(domain[faces].stateVar[3]\
			-0.5*(pow(domain[faces].stateVar[1],2)+pow(domain[faces].stateVar[2],2))/domain[faces].stateVar[0])));
		else
		{
			a_s[1]=a_s[0];
		}
		
		//speed for the boundary calculation
		a_s[0]=pow(a_s[0],2)/max(a_s[0],abs(sqrt(pow(domain[x].stateVar[1]/domain[x].stateVar[0],2)+pow(domain[x].stateVar[2]/domain[x].stateVar[0],2))));
		if(ourFlag!=4 || (ourFlag==4 && y!=note))
			a_s[1]=pow(a_s[1],2)/max(a_s[1],abs(sqrt(pow(domain[faces].stateVar[1]/domain[faces].stateVar[0],2)+pow(domain[faces].stateVar[2]/domain[faces].stateVar[0],2))));
		else
		{
			a_s[1]=a_s[0];
		}

		//Speed of sound at facial interface
		double a_mid=min(a_s[0],a_s[1]);
		//Mach number of incoming and outgoing waves
		double machplus=(domain[x].stateVar[1]/domain[x].stateVar[0]*domain[x].norms[y][0]+domain[x].stateVar[2]/domain[x].stateVar[0]*domain[x].norms[y][1])/a_mid;
		double machminus;
		if(ourFlag!=4 || (ourFlag==4 && y!=note))
		{
			machminus=(domain[faces].stateVar[1]/domain[faces].stateVar[0]*domain[x].norms[y][0]+domain[faces].stateVar[2]/domain[faces].stateVar[0]*domain[x].norms[y][1])/a_mid;
		}
		else
		{
			machminus=-machplus;
		}
		
		//Pressure Fluxes
		double pressplus=(gammma[0]-1)*(domain[x].stateVar[3]-0.5*(pow(domain[x].stateVar[1],2)+pow(domain[x].stateVar[2],2))/domain[x].stateVar[0]);
		double presminus;
		if(ourFlag!=4 || (ourFlag==4 && y!=note))
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

		domain[x].presflux[y][0]=(pressplus*plus+presminus*minus)*domain[x].norms[y][0]*sqrt(pow(domain[x].nodes[(y+1)%4][0]-domain[x].nodes[y][0],2)\
		+pow(domain[x].nodes[(y+1)%4][1]-domain[x].nodes[y][1],2));
		domain[x].presflux[y][1]=(pressplus*plus+presminus*minus)*domain[x].norms[y][1]*sqrt(pow(domain[x].nodes[(y+1)%4][0]-domain[x].nodes[y][0],2)\
		+pow(domain[x].nodes[(y+1)%4][1]-domain[x].nodes[y][1],2));
	}
}
