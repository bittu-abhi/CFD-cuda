#include "ausmPlus.h"
#include <stdio.h>
	
__global__ void diffusiveFlux(cell *domain,float *R, float *gammma, float *mu,float wall_temp,float *k)
{
	int x=blockIdx.x;
	int y=threadIdx.x;
	int note=-10;
	int faces=(int)domain[x].face[y]-1;
	int ourFlag=(int)domain[x].flag;
	float delu_delx=0.0,delv_delx=0.0,delu_dely=0.0,delv_dely=0.0;
	if(ourFlag==0 || ourFlag==4 || ourFlag==2)
	{
		float x_cord[]={0,0},y_cord[]={0,0};
		
		if(faces<0 || faces>50266)
		{
			note=y;
		}

		int i1,i2;
		if(ourFlag==4 && y==note)
		{
			i1=note;
			i2=(note+1)%3;
			x_cord[1]=(1.0/2.0)*(domain[x].nodes[i1][0]+domain[x].nodes[i2][0]);
			y_cord[1]=(1.0/2.0)*(domain[x].nodes[i1][1]+domain[x].nodes[i2][1]);
		}

		for (int i = 0; i < 3; ++i)
		{
			if(ourFlag!=4 || (ourFlag==4 && y!=note))
			{
				//x_cordinate of the elements
				x_cord[0]+=(1.0/3.0)*(domain[x].nodes[i][0]);
				x_cord[1]+=(1.0/3.0)*(domain[faces].nodes[i][0]);
				//Y coordinate of the elements
				y_cord[0]+=(1.0/3.0)*(domain[x].nodes[i][1]);
				y_cord[1]+=(1.0/3.0)*(domain[faces].nodes[i][1]);
			}
			else
			{
				//x_cordinate of the elements
				x_cord[0]+=(1.0/3.0)*(domain[x].nodes[i][0]);
				//Y coordinate of the elements
				y_cord[0]+=(1.0/3.0)*(domain[x].nodes[i][1]);
			}
		}

		if(ourFlag==2 && y==note)
		{
			i1=note;
			i2=(note+1)%3;
			x_cord[1]=(1.0/2.0)*(domain[x].nodes[i1][0]+domain[x].nodes[i2][0]);
			y_cord[1]=(1.0/2.0)*(domain[x].nodes[i1][1]+domain[x].nodes[i2][1]);
		}
		else if(ourFlag==2 && y!=note)
		{
			for (int i = 0; i < 3; ++i)
			{
				//x_cordinate of the elements
				x_cord[0]+=(1.0/3.0)*(domain[x].nodes[i][0]);
				x_cord[1]+=(1.0/3.0)*(domain[faces].nodes[i][0]);
				//Y coordinate of the elements
				y_cord[0]+=(1.0/3.0)*(domain[x].nodes[i][1]);
				y_cord[1]+=(1.0/3.0)*(domain[faces].nodes[i][1]);
			}
		}

		if(abs(x_cord[1]-x_cord[0])<=0.001)
		{
			delu_delx=0.0;
			delv_delx=0.0;
		}
		else
		{
			delu_delx=(domain[x].temp_var[y][1]/domain[x].temp_var[y][0]-domain[x].stateVar[1]/domain[x].stateVar[0])/(x_cord[1]-x_cord[0]);
			delv_delx=(domain[x].temp_var[y][2]/domain[x].temp_var[y][0]-domain[x].stateVar[2]/domain[x].stateVar[0])/(x_cord[1]-x_cord[0]);
		}
		if(abs(y_cord[1]-y_cord[0])<=0.001)
		{
			delu_dely=0.0;
			delv_dely=0.0;
		}
		else
		{
			delu_dely=(domain[x].temp_var[y][1]/domain[x].temp_var[y][0]-domain[x].stateVar[1]/domain[x].stateVar[0])/(y_cord[1]-y_cord[0]);
			delv_dely=(domain[x].temp_var[y][2]/domain[x].temp_var[y][0]-domain[x].stateVar[2]/domain[x].stateVar[0])/(y_cord[1]-y_cord[0]);
		}

		float tau_xx=2*mu[0]*(delu_delx-1/3*(delu_delx+delv_dely));
		float tau_yy=2*mu[0]*(delv_dely-1/3*(delu_delx+delv_dely));
		float tau_xy=mu[0]*(delu_dely+delv_delx);

		float temp[2];
		temp[0]=(gammma[0]-1)/R[0]*(domain[x].stateVar[3]-0.5*(powf(domain[x].stateVar[1],2)+powf(domain[x].stateVar[2],2))/domain[x].stateVar[0])/domain[x].stateVar[0];
		if(ourFlag!=4 || (ourFlag==4 && y!=note))
			temp[1]=(gammma[0]-1)/R[0]*(domain[x].temp_var[y][3]-0.5*(powf(domain[x].temp_var[y][1],2)\
				+powf(domain[x].temp_var[y][2],2))/domain[x].temp_var[y][0])/domain[x].temp_var[y][0];
		else
		{
			temp[1]=wall_temp;
		}

		float delT_delx,delT_dely;	
		if(abs(x_cord[1]-x_cord[0])<=0.001)
			delT_delx=0;
		else
			delT_delx=(temp[1]-temp[0])/(x_cord[1]-x_cord[0]);
		if(abs(y_cord[1]-y_cord[0])<=0.001)
			delT_dely=0;
		else
			delT_dely=(temp[1]-temp[0])/(y_cord[1]-y_cord[0]);

		float thetaX=domain[x].stateVar[1]/domain[x].stateVar[0]*tau_xx+domain[x].stateVar[2]/domain[x].stateVar[0]*tau_xy+k[0]*delT_delx;
		float thetaY=domain[x].stateVar[1]/domain[x].stateVar[0]*tau_xy+domain[x].stateVar[2]/domain[x].stateVar[0]*tau_yy+k[0]*delT_dely;

		domain[x].diffflux[y][0]=0;
		domain[x].diffflux[y][1]=(tau_xx*domain[x].norms[y][0]+tau_xy*domain[x].norms[y][1])\
		*sqrt(powf(domain[x].nodes[y][0]-domain[x].nodes[(y+1)%3][0],2)+powf(domain[x].nodes[y][1]-domain[x].nodes[(y+1)%3][1],2));
		domain[x].diffflux[y][2]=(tau_xy*domain[x].norms[y][0]+tau_yy*domain[x].norms[y][1])\
		*sqrt(powf(domain[x].nodes[y][0]-domain[x].nodes[(y+1)%3][0],2)+powf(domain[x].nodes[y][1]-domain[x].nodes[(y+1)%3][1],2));
		domain[x].diffflux[y][3]=(thetaX*domain[x].norms[y][0]+thetaY*domain[x].norms[y][1])\
		*sqrt(powf(domain[x].nodes[y][0]-domain[x].nodes[(y+1)%3][0],2)+powf(domain[x].nodes[y][1]-domain[x].nodes[(y+1)%3][1],2));

		/*if(abs((1.0/3.0)*(domain[x].nodes[0][0]+domain[x].nodes[1][0]+domain[x].nodes[2][0])+0.524158619)<0.001 && abs((1.0/3.0)*(domain[x].nodes[0][1]+domain[x].nodes[1][1]+domain[x].nodes[2][1])-0.8526501336)<0.001)
		{
			printf("diffusive %5.14lf %5.14lf %5.14lf %5.14lf	 %d %d %d\n",domain[x].diffflux[y][0],domain[x].diffflux[y][1],domain[x].diffflux[y][2],domain[x].diffflux[y][3],domain[x].flag,x+1,y);
			printf("%5.14lf %5.14lf %5.14lf %5.14lf\n",x_cord[0],x_cord[1],y_cord[0],y_cord[1]);
		}*/
	}
}
