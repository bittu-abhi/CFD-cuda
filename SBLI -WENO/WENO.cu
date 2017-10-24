#include <iostream>
#include <cuda.h>
#include <cuda_runtime.h>
#include <math.h>
#include "sbli.h"

void WENO(point *pt,int points,int flagXY)
{
	double cons[3] = {1/10,3/5,3/10};
	double beta[4][3]={{0},{0},{0}};
	int x = blockIdx.x;
	int y = blockIdx.y;
	int state = threadIdx.x;

	beta[state][0] = 13/12*pow((pt[(flagXY==0):(x-2+y*points)?(x+(y-2)*points)].stateVar[state]-2*pt[(flagXY==0):(x-1+y*points)?(x+(y-1)*points)].stateVar[state]+pt[(flagXY==0):(x+y*points)?(x+(y)*points)].stateVar[state]),2)+1/4*pow((pt[(flagXY==0):(x-2+y*points)?(x+(y-2)*points)].stateVar[state]-4*pt[(flagXY==0):(x-1+y*points)?(x+(y-1)*points)].stateVar[state]+3*pt[(flagXY==0):(x+y*points)?(x+(y)*points)].stateVar[state]),2);

	beta[state][1] = 13/12*pow((pt[(flagXY==0):(x-1+y*points)?(x+(y-1)*points)].stateVar[state]-2*pt[(flagXY==0):(x+y*points)?(x+(y)*points)].stateVar[state]+pt[(flagXY==0):(x+1+y*points)?(x+(y+1)*points)].stateVar[state]),2)+1/4*pow((pt[(flagXY==0):(x-1+y*points)?(x+(y-1)*points)].stateVar[state]-pt[(flagXY==0):(x+1+y*points)?(x+(y+1)*points)].stateVar[state]),2);

	beta[state][2] = 13/12*pow((pt[(flagXY==0):(x+y*points)?(x+(y)*points)].stateVar[state]-2*pt[(flagXY==0):(x+1+y*points)?(x+(y+1)*points)].stateVar[state]+pt[(flagXY==0):(x+2+y*points)?(x+(y+2)*points)].stateVar[state]),2)+1/4*pow((3*pt[(flagXY==0):(x+y*points)?(x+(y)*points)].stateVar[state]-4*pt[(flagXY==0):(x+1+y*points)?(x+(y+1)*points)].stateVar[state]+pt[(flagXY==0):(x+2+y*points)?(x+(y+2)*points)].stateVar[state]),2);

	double wtilda[4][3]={{0},{0},{0}};

	wtilda[state][0] = (1/10)/pow((beta[state][0]+0.000001),2);
	wtilda[state][1] = (3/5)/pow((beta[state][1]+0.000001),2);
	wtilda[state][2] = (3/10)/pow((beta[state][2]+0.000001),2);

	double deno = wtilda[state][0]+wtilda[state][1]+wtilda[state][2];

	double weight[4][3] = {{wtilda[state][0]/deno},{wtilda[state][1]/deno},{wtilda[state][2]/deno}};

	double recon_poly[4][3] =  {{0},{0},{0}};

	recon_poly[state][0] = 1/3*(pt[(flagXY==0):(x-2+y*points)?(x+(y-2)*points)].stateVar[state])-7/6*(pt[(flagXY==0):(x-1+y*points)?(x+(y-1)*points)].stateVar[state])+11/6*(pt[(flagXY==0):(x+y*points)?(x+(y)*points)].stateVar[state]);

	recon_poly[state][1] = -1/6*(pt[(flagXY==0):(x-1+y*points)?(x+(y-1)*points)].stateVar[state])+5/6*(pt[(flagXY==0):(x+y*points)?(x+(y)*points)].stateVar[state])+1/3*(pt[(flagXY==0):(x+1+y*points)?(x+(y+1)*points)].stateVar[state]);

	recon_poly[state][2] = 1/3*(pt[(flagXY==0):(x+y*points)?(x+(y)*points)].stateVar[state])+5/6*(pt[(flagXY==0):(x+1+y*points)?(x+(y+1)*points)].stateVar[state])-1/6*(pt[(flagXY==0):(x+2+y*points)?(x+(y+2)*points)].stateVar[state]);

	pt[x+y*points].interface[state][flagXY] = weight[state][0]*recon_poly[state][0] + weight[state][1]*recon_poly[state][1] + weight[state][2]*recon_poly[state][2]; //It can be returned. Same for i-1/2
}

