#include "ausmPlus.h"
#include <stdio.h>
	
__global__ void diffusiveFlux(cell *domain,double *R, double *gammma, double *mu,double wall_temp,double *k)
{
	int x=blockIdx.x;
	int y=threadIdx.x;
	int note=-10;
	int faces=(int)domain[x].face[y]-1;
	int ourFlag=(int)domain[x].flag;
	double delu_delx=0.0,delv_delx=0.0,delu_dely=0.0,delv_dely=0.0;
	if(ourFlag==0 || ourFlag==4)
	{
		double x_cord[]={0,0},y_cord[]={0,0};
		
		if(faces<0 || faces>42860)
		{
			note=y;
		}

		int i1,i2;
		if(ourFlag==4 && y==note)
		{
			i1=note;
			i2=(note+1)%4;
			x_cord[1]=0.5*(domain[x].nodes[i1][0]+domain[x].nodes[i2][0]);
			y_cord[1]=0.5*(domain[x].nodes[i1][1]+domain[x].nodes[i2][1]);
		}

		for (int i = 0; i < 4; ++i)
		{
			if(ourFlag!=4 || (ourFlag==4 && y!=note))
			{
				//x_cordinate of the elements
				x_cord[0]+=0.25*(domain[x].nodes[i][0]);
				x_cord[1]+=0.25*(domain[faces].nodes[i][0]);
				//Y coordinate of the elements
				y_cord[0]+=0.25*(domain[x].nodes[i][1]);
				y_cord[1]+=0.25*(domain[faces].nodes[i][1]);
			}
			else
			{
				//x_cordinate of the elements
				x_cord[0]+=0.25*(domain[x].nodes[i][0]);
				//Y coordinate of the elements
				y_cord[0]+=0.25*(domain[x].nodes[i][1]);
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

		double tau_xx=2*mu[0]*(delu_delx-1/2*(delu_delx+delv_dely));
		double tau_yy=2*mu[0]*(delv_dely-1/2*(delu_delx+delv_dely));
		double tau_xy=mu[0]*(delu_dely+delv_delx);

		double temp[2];
		temp[0]=(gammma[0]-1)/R[0]*(domain[x].stateVar[3]-0.5*(pow(domain[x].stateVar[1],2)+pow(domain[x].stateVar[2],2))/domain[x].stateVar[0])/domain[x].stateVar[0];
		if(ourFlag!=4 || (ourFlag==4 && y!=note))
			temp[1]=(gammma[0]-1)/R[0]*(domain[x].temp_var[y][3]-0.5*(pow(domain[x].temp_var[y][1],2)\
				+pow(domain[x].temp_var[y][2],2))/domain[x].temp_var[y][0])/domain[x].temp_var[y][0];
		else
		{
			temp[1]=wall_temp;
		}

		double delT_delx,delT_dely;	
		if(abs(x_cord[1]-x_cord[0])<=0.001)
			delT_delx=0;
		else
			delT_delx=(temp[1]-temp[0])/(x_cord[1]-x_cord[0]);
		if(abs(y_cord[1]-y_cord[0])<=0.001)
			delT_dely=0;
		else
			delT_dely=(temp[1]-temp[0])/(y_cord[1]-y_cord[0]);

		double thetaX=domain[x].stateVar[1]/domain[x].stateVar[0]*tau_xx+domain[x].stateVar[2]/domain[x].stateVar[0]*tau_xy+k[0]*delT_delx;
		double thetaY=domain[x].stateVar[1]/domain[x].stateVar[0]*tau_xy+domain[x].stateVar[2]/domain[x].stateVar[0]*tau_yy+k[0]*delT_dely;

		domain[x].diffflux[y][0]=0;
		domain[x].diffflux[y][1]=(tau_xx*domain[x].norms[y][0]+tau_xy*domain[x].norms[y][1])\
		*sqrt(pow(domain[x].nodes[y][0]-domain[x].nodes[(y+1)%4][0],2)+pow(domain[x].nodes[y][1]-domain[x].nodes[(y+1)%4][1],2));
		domain[x].diffflux[y][2]=(tau_xy*domain[x].norms[y][0]+tau_yy*domain[x].norms[y][1])\
		*sqrt(pow(domain[x].nodes[y][0]-domain[x].nodes[(y+1)%4][0],2)+pow(domain[x].nodes[y][1]-domain[x].nodes[(y+1)%4][1],2));
		domain[x].diffflux[y][3]=(thetaX*domain[x].norms[y][0]+thetaY*domain[x].norms[y][1])\
		*sqrt(pow(domain[x].nodes[y][0]-domain[x].nodes[(y+1)%4][0],2)+pow(domain[x].nodes[y][1]-domain[x].nodes[(y+1)%4][1],2));

		/*if(abs(0.25*(domain[x].nodes[0][0]+domain[x].nodes[1][0]+domain[x].nodes[2][0]+domain[x].nodes[3][0])-0.375)<0.000001 && abs(0.25*(domain[x].nodes[0][1]+domain[x].nodes[1][1]+domain[x].nodes[2][1]+domain[x].nodes[3][1])-29.875)<0.00001)
			printf("upper %5.14lf %5.14lf %5.14lf %5.14lf %5.14lf %d %d %d\n",domain[x].diffflux[y][2],delu_dely,(domain[x].temp_var[y][1]/domain[x].temp_var[y][0]-domain[x].stateVar[1]/domain[x].stateVar[0])/(y_cord[1]-y_cord[0]),y_cord[1]-y_cord[0],(domain[x].temp_var[y][1]/domain[x].temp_var[y][0]-domain[x].stateVar[1]/domain[x].stateVar[0]),x,y,note );
	
		if(abs(0.25*(domain[x].nodes[0][0]+domain[x].nodes[1][0]+domain[x].nodes[2][0]+domain[x].nodes[3][0])-0.3757366911)<0.000001 && abs(0.25*(domain[x].nodes[0][1]+domain[x].nodes[1][1]+domain[x].nodes[2][1]+domain[x].nodes[3][1])-13.8545675278)<0.00001)
			printf("lower %5.14lf %5.14lf %5.14lf %5.14lf %5.14lf %d %d %d\n",domain[x].diffflux[y][2],delu_dely,(domain[x].temp_var[y][1]/domain[x].temp_var[y][0]-domain[x].stateVar[1]/domain[x].stateVar[0])/(y_cord[1]-y_cord[0]), y_cord[1]-y_cord[0],(domain[x].temp_var[y][1]/domain[x].temp_var[y][0]-domain[x].stateVar[1]/domain[x].stateVar[0]),x,y,note );*/
	}
}	