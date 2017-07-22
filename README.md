# CFD using CUDA
Cuda based CFD codes on different bodies. The goal is to learn parallel processing and computational fluid dynamics using GPU, Since GPU uses a massive parallel architecture. I have not completely optimised any code. 

(Currently its a cylinder of radius 2 center at (20,15).The scheme for convective fluxes is AUSM+, while for diffusive fluxes is simply divided difference.At present I am having oscillation in the solution. Fluid is assumed to be air. Domain is a rectangle of size 100x30.For this I have 25000 mesh elements on which 25000 grid block operte using 4 threads each giving a total of 100000 threads.)
