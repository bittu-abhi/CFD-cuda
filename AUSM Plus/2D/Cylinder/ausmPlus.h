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
	double face[4];

	//Contains the outward normal vector to the corresponding face of the element(nx,ny),
	double norms[4][2];

	//ALl the fluxes for the paricular element
	//The convective flux
	double convflux[4][4];
	//The diffusive flux
	double diffflux[4][4];
	//The pressure flux
	double presflux[4][2];

	//To store the values of the next face element so that there is no race conditions.Columns are the 
	//state variables. Rows are faces.
	double temp_var[4][4];

	//constructor(s)
	cell(double *state);
	cell();
};

//Set the nodes of each element-co-ordinates, and the node number
__global__ void set_nodes(double *node, cell *domain, double *boundary);

//Set the element number of the faes surrounding the current element
__global__ void set_neighbour(cell *domain);

//Calcualte the pressure flux using AUSM+ scheme
__global__ void pressureFlux(cell *domain, double *R, double *gammma);

//Calcualte the convective flux using AUSM+ scheme
__global__ void convectiveflux(cell *domain, double *R, double *gammma);

//Calcualte the pressure flux using forward divided difference scheme
__global__ void diffusiveFlux(cell *domain,double *R, double *gammma, double *mu,double wall_temp,double *k);

//Calcualte the facial norms of the four faces
__global__ void calculate_norm(cell *domain);

//Gather the state variables of the different cell values surrounding the current element to avoid race conditions(if any)
__global__ void read_values(cell *domain);

//AUSM+ scheme
void ausmplus(double *initial,double timesteps, double deltat);

//Write the values of final state variables, fluxes in different .csv files
void visual(cell *domain);

extern double gammma;
extern double mu;
extern double k;
extern double R;

#endif
