#ifndef ausmPlus_H
#define ausmPlus_H

class cell
{
public:
	//stateVar is the array containing the state variable rho, rhoU, rhoV, rhoE
	//E is the total energy = internal energy plus kinetic energy
	double stateVar[4];

	//Co-ordinates of the four nodes, First column is the X, then is Y, third is the node number
	double nodes[4][3];

	//Type of cell, 0-fluid, 1-Inlet,2-Outflow, 3-Farfield, 4-Cylinder wall
	int flag;

	//Contains the respective faces connected to and the nodes common tho them
	//First row is the id of the next element,, next two are the nodes 
	double face[4][3];

	//ALl the fluxes for the paricular element
	//The convective flux
	double convflux[4][4];
	//The diffusive flux
	double diffflux[4][4];
	//The pressure flux
	double presflux[4][2];

	//constructor(s)
	cell(double *state);
	cell();
};

__global__ void set_nodes(double *node, cell *domain, double *boundary);

__global__ void set_neighbour(cell *domain);

__global__ void pressureFlux(cell *domain, double *R, double *gammma);

__global__ void convectiveflux(cell *domain, double *R, double *gammma);

__global__ void diffusiveFlux(cell *domain,double *R, double *gammma, double *mu,double wall_temp,double *k);

void ausmplus(double *initial,double timesteps, double deltat);

void visual(cell *domain);

extern double gammma;
extern double mu;
extern double k;
extern double R;

#endif