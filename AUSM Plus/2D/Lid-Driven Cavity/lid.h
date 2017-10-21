#ifndef lid_H
#define lid_H

class point
{
public:
	//Voricity values
	double w;

	//stream functions values
	double psi;

	//first and second derivative values after getting calculated from compact scheme
	//First column is for derivative wrt x, second for derivative wrt y
	double delW[2];
	double delW2[2];
	double delpsi[2];

	double error;
	point();

};

//Streamfunction equation solver using jacobi iteration
__global__ void streamfunc(double deltax, double deltay, point *pt,int points, double convergence, int *address);

//Complete flowchart of the lid driven cavity is implemented 
void LID(double deltax, double deltay, double deltat,int points, double convergence, int iterations, double reynold, double timesteps);

//Forward time stepping of omega
__global__ void timeW(point *curr, double deltat, double re, int points);

//calculations of coefficients using 6th order comapct scheme for 1st derivative
__global__ void compact1D(point *pt, int points, double delta, int flagXY, int flagPSIW);

//calculations of coefficients using 6th order comapct scheme for 2nd derivative
__global__ void compact2D(point *pt, int points, double delta, int flagXY);

//Visualisation file writer, to be opened in paraview
void visu(point *ptr,double t,double deltat,int points,double deltax, double deltay);

__global__ void updateboundary(point *pt, int points, double deltay, double deltax);
//Uses a second order discretisation for liznearisation of the PDEs
__global__ void secondW(point *pt, int points, double delta,int flagXY );

#endif