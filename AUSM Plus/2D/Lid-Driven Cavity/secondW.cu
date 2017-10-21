#include "lid.h"
#include <stdio.h>

__global__ void secondW(point *pt, int points, double delta, int flagXY )
{
	int y_level=blockIdx.x;
	int x = threadIdx.x;

	if(x!=0 && y_level!=0 && x!=points-1 && y_level!=points-1)
	{
		pt[(x+y_level*points)].delW[flagXY] = (pt[(flagXY==0)?(x+1+y_level*points):((y_level+1)*points+x)].w-pt[(flagXY==0)?(x-1+y_level*points):((y_level-1)*points+x)].w)/(2*delta);

		pt[(x+y_level*points)].delpsi[flagXY] = (pt[(flagXY==0)?(x+1+y_level*points):((y_level+1)*points+x)].psi-pt[(flagXY==0)?(x-1+y_level*points):((y_level-1)*points+x)].psi)/(2*delta);


		pt[(x+y_level*points)].delW2[flagXY]=(pt[(flagXY==0)?(x-1+points*y_level):((y_level-1)*points+x)].w-2*pt[(flagXY==0)?(x+points*y_level):(y_level*points+x)].w+pt[(flagXY==0)?(x+1+points*y_level):((y_level+1)*points+x)].w)/(delta*delta);
	}
	if(x==points-1)
	{

		if(y_level!=(points-1) && y_level!=0)
		{	
			pt[(x+y_level*points)].delW[flagXY] = (flagXY==0)?((pt[(x+points*y_level)].w-pt[((x-1)+points*y_level)].w)/(delta)): (pt[((y_level+1)*points+x)].w-pt[((y_level-1)*points+x)].w)/(2*delta);

			pt[(x+y_level*points)].delpsi[flagXY] = (flagXY==0)?((pt[(x+points*y_level)].psi-pt[((x-1)+points*y_level)].psi)/(delta)): (pt[((y_level+1)*points+x)].psi-pt[((y_level
				-1)*points+x)].psi)/(2*delta);

			pt[(x+y_level*points)].delW2[flagXY] = (flagXY==0)?(2*(pt[(x+points*y_level)].w-5*pt[((x-1)+points*y_level)].w+4*pt[((x-2)+points*y_level)].w-pt[((x-3)+points*y_level)].w)/(delta*delta)): (pt[((y_level-1)*points+x)].w-2*pt[(y_level*points+x)].w+pt[((y_level+1)*points+x)].w)/(delta*delta);
		}
	}
	if(x==0)
	{
		if(y_level!=(points-1) && y_level!=0)
		{
			pt[x+y_level*points].delW[flagXY] = (flagXY==0)?((pt[(x+1+points*y_level)].w-pt[((x)+points*y_level)].w)/(delta)): (pt[((y_level+1)*points+x)].w-pt[((y_level-1)*points+x)].w)/(2*delta);

			pt[(x+y_level*points)].delpsi[flagXY] = (flagXY==0)?((pt[(x+1+points*y_level)].psi-pt[((x)+points*y_level)].psi)/(delta)): (pt[((y_level+1)*points+x)].psi-pt[((y_level-1)*points+x)].psi)/(2*delta);

			pt[(x+y_level*points)].delW2[flagXY] = (flagXY==0)?(2*(pt[(x+points*y_level)].w-5*pt[((x+1)+points*y_level)].w+4*pt[((x+2)+points*y_level)].w-pt[((x+3)+points*y_level)].w)/(delta*delta)): (pt[((y_level-1)*points+x)].w-2*pt[(y_level*points+x)].w+pt[((y_level+1)*points+x)].w)/(delta*delta);
		}
	}
	if(y_level==0)
	{
		if(x!=(points-1) && x!=0)
		{
			pt[(x+y_level*points)].delW[flagXY] = (flagXY==0)?((pt[(x+1+points*y_level)].w-pt[((x-1)+points*y_level)].w)/(2*delta)): (pt[((y_level+1)*points+x)].w-pt[((y_level)*points+x)].w)/(delta);

			pt[(x+y_level*points)].delpsi[flagXY] = (flagXY==0)?((pt[(x+1+points*y_level)].psi-pt[((x-1)+points*y_level)].psi)/(2*delta)): (pt[((y_level+1)*points+x)].psi-pt[((y_level)*points+x)].psi)/(delta);

			pt[(x+y_level*points)].delW2[flagXY] = (flagXY==0)?((pt[(x+1+points*y_level)].w+pt[((x-1)+points*y_level)].w-2*pt[((x)+points*y_level)].w)/(delta*delta)): (2*pt[((y_level)*points+x)].w-5*pt[((y_level+1)*points+x)].w+4*pt[((y_level+2)*points+x)].w-pt[((y_level+3)*points+x)].w)/(delta*delta);
		}	
	}
	if(y_level==(points-1))
	{
		if(x!=(points-1) && x!=0)
		{
			pt[(x+y_level*points)].delW[flagXY] = (flagXY==0)?((pt[(x+1+points*y_level)].w-pt[((x-1)+points*y_level)].w)/(2*delta)): (pt[((y_level)*points+x)].w-pt[((y_level-1)*points+x)].w)/(delta);

			pt[(x+y_level*points)].delpsi[flagXY] = (flagXY==0)?((pt[(x+1+points*y_level)].psi-pt[((x-1)+points*y_level)].psi)/(2*delta)): (pt[((y_level)*points+x)].psi-pt[((y_level-1)*points+x)].psi)/(delta);

			pt[(x+y_level*points)].delW2[flagXY] = (flagXY==0)?((pt[(x+1+points*y_level)].w+pt[((x-1)+points*y_level)].w-2*pt[((x)+points*y_level)].w)/(delta*delta)): (2*pt[((y_level)*points+x)].w-5*pt[((y_level-1)*points+x)].w+4*pt[((y_level-2)*points+x)].w-pt[((y_level-3)*points+x)].w)/(delta*delta);
		}
	}

	//Four Cornor points
	if(y_level==(points-1) && x==points-1)
	{
		pt[(x+y_level*points)].delW[flagXY] = (flagXY==0)?((pt[(x+points*y_level)].w-pt[((x-1)+points*y_level)].w)/(delta)): (pt[((y_level)*points+x)].w-pt[((y_level-1)*points+x)].w)/(delta);

		pt[(x+y_level*points)].delpsi[flagXY] = (flagXY==0)?((pt[(x+points*y_level)].psi-pt[((x-1)+points*y_level)].psi)/(delta)): (pt[((y_level)*points+x)].psi-pt[((y_level-1)*points+x)].psi)/(delta);

		pt[(x+y_level*points)].delW2[flagXY] = (flagXY==0)?(2*(pt[(x+points*y_level)].w-5*pt[((x-1)+points*y_level)].w+4*pt[((x-2)+points*y_level)].w-pt[((x-3)+points*y_level)].w)/(delta*delta)): (2*pt[((y_level)*points+x)].w-5*pt[((y_level-1)*points+x)].w+4*pt[((y_level-2)*points+x)].w-pt[((y_level-3)*points+x)].w)/(delta*delta);
	}
	if(y_level==(points-1) && x==0)
	{
		pt[(x+y_level*points)].delW[flagXY] = (flagXY==0)?((pt[((x+1)+points*y_level)].w-pt[((x)+points*y_level)].w)/(delta)): (pt[((y_level)*points+x)].w-pt[((y_level-1)*points+x)].w)/(delta);

		pt[(x+y_level*points)].delpsi[flagXY] = (flagXY==0)?((pt[((x+1)+points*y_level)].psi-pt[((x+0)+points*y_level)].psi)/(delta)): (pt[((y_level)*points+x)].psi-pt[((y_level-1)*points+x)].psi)/(delta);

		pt[(x+y_level*points)].delW2[flagXY] = (flagXY==0)?(2*(pt[(x+points*y_level)].w-5*pt[((x+1)+points*y_level)].w+4*pt[((x+2)+points*y_level)].w-pt[((x+3)+points*y_level)].w)/(delta*delta)): (2*pt[((y_level)*points+x)].w-5*pt[((y_level-1)*points+x)].w+4*pt[((y_level-2)*points+x)].w-pt[((y_level-3)*points+x)].w)/(delta*delta);
	}
	if(y_level==0 && x==(points-1))
	{
		pt[(x+y_level*points)].delW[flagXY] = (flagXY==0)?((pt[(x+points*y_level)].w-pt[((x-1)+points*y_level)].w)/(delta)): (pt[((y_level+1)*points+x)].w-pt[((y_level)*points+x)].w)/(delta);

		pt[(x+y_level*points)].delpsi[flagXY] = (flagXY==0)?((pt[(x+points*y_level)].psi-pt[((x-1)+points*y_level)].psi)/(delta)): (pt[((y_level+1)*points+x)].psi-pt[((y_level)*points+x)].psi)/(delta);

		pt[(x+y_level*points)].delW2[flagXY] = (flagXY==0)?(2*(pt[(x+points*y_level)].w-5*pt[((x-1)+points*y_level)].w+4*pt[((x-2)+points*y_level)].w-pt[((x-3	)+points*y_level)].w)/(delta*delta)): (2*pt[((y_level)*points+x)].w-5*pt[((y_level+1)*points+x)].w+4*pt[((y_level+2)*points+x)].w-pt[((y_level+3)*points+x)].w)/(delta*delta);
	}
	if(y_level==0 && x==0)
	{
		pt[(x+y_level*points)].delW[flagXY] = (flagXY==0)?((pt[((x+1)+points*y_level)].w-pt[((x)+points*y_level)].w)/(delta)): (pt[((y_level+1)*points+x)].w-pt[((y_level)*points+x)].w)/(delta);

		pt[(x+y_level*points)].delpsi[flagXY] = (flagXY==0)?((pt[((x+1)+points*y_level)].psi-pt[((x)+points*y_level)].psi)/(delta)): (pt[((y_level+1)*points+x)].psi-pt[((y_level)*points+x)].psi)/(delta);

		pt[(x+y_level*points)].delW2[flagXY] = (flagXY==0)?(2*(pt[(x+points*y_level)].w-5*pt[((x+1)+points*y_level)].w+4*pt[((x+2)+points*y_level)].w-pt[((x+3)+points*y_level)].w)/(delta*delta)): (2*pt[((y_level)*points+x)].w-5*pt[((y_level+1)*points+x)].w+4*pt[((y_level+2)*points+x)].w-pt[((y_level+3)*points+x)].w)/(delta*delta);
	}

	//if(y_level==(points-1) && x==(points-1))
	//		printf("product=%5.14lf, delWX=%5.14lf, delWY=%5.14lf, flagXY=%d  %d %d\n",(flagXY==0)?((pt[(x+points*y_level)].w-pt[((x-1)+points*y_level)].w)/(delta)): (pt[((y_level)*points+x)].w-pt[((y_level-1)*points+x)].w)/(delta),pt[(x+points*y_level)].delW[0],pt[(x+points*y_level)].delW[1],flagXY,x,y_level);
}