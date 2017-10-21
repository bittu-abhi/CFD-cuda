#include <iostream>
#include <cuda.h>
#include <cuda_runtime.h>
#include <omp.h>
#include <algorithm>
#include <math.h>
#include "lid.h"

using namespace std;

point::point()
{}

int main()
{
	//number of points
	int points = 500;

	//dimension in x and y direction
	double x_dim = 1;
	double y_dim = 1;

	//delta values in different directions and time
	double deltax = x_dim/points;
	double deltay = y_dim/points;
	double deltat = 0.0001;

	//Total time steps
	double timesteps = 1000000;
	//Convergence cirterion
	double convergence = 0.000001;

	//total number of iterations to be performed
	int iter = 5000;

	//Velocity of the plate
	double u = 1;

	//density of the fluid
	double density = 1.225;
	
	//viscosity of the fluid
	double viscosity = 1.225*pow(10,-3);
	
	//Reynolds number calculated according to the dimensions of the cavity
	double reynold = max(x_dim,y_dim)*density*u/viscosity;

	cout<<"Reynolds number="<<reynold<<endl;

	LID(deltax, deltay, deltat, points, convergence, iter, reynold, timesteps);

	return 0;
}