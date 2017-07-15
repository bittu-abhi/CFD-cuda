#include "ausmPlus.h"
	
__global__ void diffusiveFlux(cell *domain,double *R, double *gammma, double *mu,double wall_temp,double *k)
{
	int x=blockIdx.x;
	int y=threadIdx.x;
	int note,faces=(int)domain[x].face[y][0];
	if(domain[x].flag==0 || domain[x].flag==4)
	{
		double x_cord[]={0,0},y_cord[]={0,0};
		if(domain[x].face[y][0]<1 || domain[x].face[y][0]>26000)
			note=y;
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
		for (int i = 0; i < 4; ++i)
		{
			if(domain[x].flag!=4)
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
		if(domain[x].flag==4)
		{
			x_cord[1]=0.5*(domain[x].nodes[i1][0]+domain[x].nodes[i2][0]);
			y_cord[1]=0.5*(domain[x].nodes[i1][1]+domain[x].nodes[i2][1]);
		}
		double delu_delx,delv_delx,delu_dely,delv_dely;
		if(y!=note)
		{
			if(abs(x_cord[1]-x_cord[0])<=0.001)
			{
				delu_delx=0;
				delv_delx=0;
			}
			else
			{
				delu_delx=(domain[faces].stateVar[1]/domain[faces].stateVar[0]-domain[x].stateVar[1]/domain[x].stateVar[0])/(x_cord[1]-x_cord[0]);
				delv_delx=(domain[faces].stateVar[2]/domain[faces].stateVar[0]-domain[x].stateVar[2]/domain[x].stateVar[0])/(x_cord[1]-x_cord[0]);
			}
			if(abs(y_cord[1]-y_cord[0])<=0.001)
			{
				delu_dely=0;
				delv_dely=0;
			}
			else
			{
				delu_dely=(domain[faces].stateVar[1]/domain[faces].stateVar[0]-domain[x].stateVar[1]/domain[x].stateVar[0])/(y_cord[1]-y_cord[0]);
				delv_dely=(domain[faces].stateVar[2]/domain[faces].stateVar[0]-domain[x].stateVar[2]/domain[x].stateVar[0])/(y_cord[1]-y_cord[0]);
			}
		}
		else
		{
			if(abs(x_cord[1]-x_cord[0])<=0.001)
			{
				delu_delx=0;
				delv_delx=0;
			}
			else
			{
				delu_delx=(0-domain[x].stateVar[1]/domain[x].stateVar[0])/(x_cord[1]-x_cord[0]);
				delv_delx=(0-domain[x].stateVar[2]/domain[x].stateVar[0])/(x_cord[1]-x_cord[0]);
			}
			if(abs(y_cord[1]-y_cord[0])<=0.001)
			{
				delu_dely=0;
				delv_dely=0;
			}
			else
			{
				delu_dely=(0-domain[x].stateVar[1]/domain[x].stateVar[0])/(y_cord[1]-y_cord[0]);
				delv_dely=(0-domain[x].stateVar[2]/domain[x].stateVar[0])/(y_cord[1]-y_cord[0]);
			}
		}

		double tau_xx=2*mu[0]*(delu_delx-1/3*(delu_delx+delv_dely));
		double tau_yy=2*mu[0]*(delv_dely-1/3*(delu_delx+delv_dely));
		double tau_xy=mu[0]*(delu_dely+delv_delx);

		double temp[2];
		temp[0]=(gammma[0]-1)/R[0]*(domain[x].stateVar[3]-0.5*(pow(domain[x].stateVar[1],2)+pow(domain[x].stateVar[2],2))/domain[x].stateVar[0])/domain[x].stateVar[0];
		if(domain[x].flag!=4)
			temp[1]=(gammma[0]-1)/R[0]*(domain[faces].stateVar[3]-0.5*(pow(domain[faces].stateVar[1],2)\
				+pow(domain[faces].stateVar[2],2))/domain[faces].stateVar[0])/domain[faces].stateVar[0];
		else
		{
			temp[1]=wall_temp;
		}

		double delT_delx,delT_dely;
		if((x_cord[1]-x_cord[0])==0)
			delT_delx=0;
		else
			delT_delx=(temp[1]-temp[0])/(x_cord[1]-x_cord[0]);
		if((y_cord[1]-y_cord[0])==0)
			delT_dely=0;
		else
			delT_dely=(temp[1]-temp[0])/(y_cord[1]-y_cord[0]);

		double thetaX=domain[x].stateVar[1]/domain[x].stateVar[0]*tau_xx+domain[x].stateVar[2]/domain[x].stateVar[0]*tau_xy+k[0]*delT_delx;
		double thetaY=domain[x].stateVar[1]/domain[x].stateVar[0]*tau_xy+domain[x].stateVar[2]/domain[x].stateVar[0]*tau_yy+k[0]*delT_dely;

		if(domain[x].flag!=4)
		{
			domain[x].diffflux[y][0]=0;
			domain[x].diffflux[y][1]=(tau_xx*(domain[x].nodes[i1][0]-domain[x].nodes[i2][0])+tau_xy*(domain[x].nodes[i2][1]-domain[x].nodes[i1][1]));
			domain[x].diffflux[y][2]=(tau_xy*(domain[x].nodes[i1][0]-domain[x].nodes[i2][0])+tau_yy*(domain[x].nodes[i2][1]-domain[x].nodes[i1][1]));
			domain[x].diffflux[y][3]=(thetaX*(domain[x].nodes[i1][0]-domain[x].nodes[i2][0])+thetaY*(domain[x].nodes[i2][1]-domain[x].nodes[i1][1]));
		}
		else
		{
			domain[x].diffflux[y][0]=0;
			domain[x].diffflux[y][1]=(tau_xx*(domain[x].nodes[i1][0]-domain[x].nodes[i2][0])+tau_xy*(domain[x].nodes[i2][1]-domain[x].nodes[i1][1]));
			domain[x].diffflux[y][2]=(tau_xy*(domain[x].nodes[i1][0]-domain[x].nodes[i2][0])+tau_yy*(domain[x].nodes[i2][1]-domain[x].nodes[i1][1]));
			domain[x].diffflux[y][3]=(thetaX*(domain[x].nodes[i1][0]-domain[x].nodes[i2][0])+thetaY*(domain[x].nodes[i2][1]-domain[x].nodes[i1][1]));
		}
	}
}
