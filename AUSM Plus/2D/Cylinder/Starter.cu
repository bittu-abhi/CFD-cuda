#include <iostream>
#include <fstream>
#include <cmath>
#include "ausmPlus.h"

using namespace std;

double gammma;
double mu;
double k;
double R;	

int main()
{
	gammma=1.4;
	mu=pow(1.798,-5);
	k=0.0251;
	R=286.9;
	
	double initial[4];
	//Rho
	initial[0]=1.225;
	//Rho*U
	initial[1]=2;
	//Rho *V
	initial[2]=0;
	//Rho*E, E is the internal energy including the kinetic energy(i.e. total intenal energy)
	initial[3]=101325.;
	//Time steps and delta_t
	double timesteps=1;
	double deltat=0.0001;
	ausmplus(initial,timesteps,deltat);

	return 0;
}