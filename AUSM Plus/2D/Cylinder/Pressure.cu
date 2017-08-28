#include "ausmPlus.h"
#include <math.h>
#include <algorithm>
#include <stdio.h>

__global__ void pressureFlux(cell *domain, float *R, float *gammma)
{
	int y=threadIdx.x;
	int x=blockIdx.x;
	int ourFlag=(int)domain[x].flag;
	if(ourFlag==0 || ourFlag==4 || ourFlag==2)
	{


		//Calculating the critical speed of sound for all the four sides/faces and the element itself
		float a_s[2];

		//Element
		a_s[0]=sqrt(2*(gammma[0]-1)/(gammma[0]+1)*(domain[x].stateVar[3]+(gammma[0]-1)*(domain[x].stateVar[3]\
			-0.5*(powf(domain[x].stateVar[1],2)+powf(domain[x].stateVar[2],2))/domain[x].stateVar[0]))/domain[x].stateVar[0]);
		//Side/face
		a_s[1]=sqrt(2*(gammma[0]-1)/(gammma[0]+1)*(domain[x].temp_var[y][3]+(gammma[0]-1)*(domain[x].temp_var[y][3]\
			-0.5*(powf(domain[x].temp_var[y][1],2)+powf(domain[x].temp_var[y][2],2))/domain[x].temp_var[y][0]))/domain[x].temp_var[y][0]);
		
		//speed for the boundary calculation
		a_s[0]=powf(a_s[0],2)/max(a_s[0],abs(sqrt(powf(domain[x].stateVar[1]/domain[x].stateVar[0],2)+powf(domain[x].stateVar[2]/domain[x].stateVar[0],2))));
		a_s[1]=powf(a_s[1],2)/max(a_s[1],abs(sqrt(powf(domain[x].temp_var[y][1]/domain[x].temp_var[y][0],2)+powf(domain[x].temp_var[y][2]/domain[x].temp_var[y][0],2))));

		//Speed of sound at facial interface
		float a_mid=min(a_s[0],a_s[1]);

		//Mach number of incoming and outgoing waves
		float mach_one=(domain[x].stateVar[1]/domain[x].stateVar[0]*domain[x].norms[y][0]+domain[x].stateVar[2]/domain[x].stateVar[0]*domain[x].norms[y][1])/a_mid;
		float mach_two=(domain[x].temp_var[y][1]/domain[x].temp_var[y][0]*domain[x].norms[y][0]+domain[x].temp_var[y][2]/domain[x].temp_var[y][0]*domain[x].norms[y][1])/a_mid;
		
		//Pressure Fluxes
		float press_one=(*gammma-1)*(domain[x].stateVar[3]-0.5*(powf(domain[x].stateVar[1],2)+powf(domain[x].stateVar[2],2))/domain[x].stateVar[0]);
		float press_two=(*gammma-1)*(domain[x].temp_var[y][3]-0.5*(powf(domain[x].temp_var[y][1],2)+powf(domain[x].temp_var[y][2],2))/domain[x].temp_var[y][0]);

		float one,two;
		if(abs(mach_one)>=1)
			one=0.5*(1+mach_one/abs(mach_one));
		else
			one=0.25*powf((mach_one+1),2)*(2-mach_one)+3/16*mach_one*powf((powf(mach_one,2)-1),2);
	
		if(abs(mach_two)>=1)
			two=0.5*(1-mach_two/abs(mach_two));
		else
			two=0.25*powf((mach_two-1),2)*(2+mach_two)-3/16*mach_two*powf((powf(mach_two,2)-1),2);

		domain[x].presflux[y][0]=(press_one*one+press_two*two)*domain[x].norms[y][0]*sqrt(powf(domain[x].nodes[y][0]-domain[x].nodes[(y+1)%3][0],2)\
			+powf(domain[x].nodes[y][1]-domain[x].nodes[(y+1)%3][1],2));
		domain[x].presflux[y][1]=(press_one*one+press_two*two)*domain[x].norms[y][1]*sqrt(powf(domain[x].nodes[y][0]-domain[x].nodes[(y+1)%3][0],2)\
			+powf(domain[x].nodes[y][1]-domain[x].nodes[(y+1)%3][1],2));

		/*if(abs((1.0/3.0)*(domain[x].nodes[0][0]+domain[x].nodes[1][0]+domain[x].nodes[2][0])+0.524158619)<0.001 && abs((1.0/3.0)*(domain[x].nodes[0][1]+domain[x].nodes[1][1]+domain[x].nodes[2][1])-0.8526501336)<0.001)
		{
			printf("pressure  %5.14lf %5.14lf	 %d %d %d\n",domain[x].presflux[y][0],domain[x].presflux[y][1],domain[x].flag,x+1,y);
			//printf("mach(%5.14lf %5.14lf) pressure(%5.14lf %5.14lf) %d %d\n",one,two,press_one,press_two,x+1,y );

		}*/
	}
}
