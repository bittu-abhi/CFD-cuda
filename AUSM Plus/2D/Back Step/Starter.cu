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
	mu=pow(1.789,-5);
	k=0.6065;
	R=186.9;
	
	double initial[4];
	//Rho
	initial[0]=1.225;
	//Rho*U
	initial[1]=200.000*initial[0];
	//Rho *V
	initial[2]=0.0000;
	//Rho*E, E is the internal energy including the kinetic energy(i.e. total intenal energy)
	initial[3]=(101325.000-0.5000*(pow(initial[1],2.0000)+pow(initial[2],2.0000))/initial[0])/(gammma-1.0000)+\
	0.5000*(pow(initial[1],2.0000)+pow(initial[2],2.0000))/initial[0];
	//Time steps and delta_t
	double timesteps=1000000;
	double deltat=0.00001;
	ausmplus(initial,timesteps,deltat);

	return 0;
}