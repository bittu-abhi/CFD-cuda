# CFD using CUDA
Cuda based CFD codes on different bodies (currently its a cylinder of radius 2 center at (20,15). 
The scheme for convective fluxes is AUSM+, while for diffusive fluxes is simply divided difference.

At present I am having oscillation in the solution. Fluid is assumed to be air. Domain is a rectangle of size 
100x30.
