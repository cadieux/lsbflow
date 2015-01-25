# lsbflow
Pseudo-spectral Navier-Stokes Solver for Laminar Separation Bubble Flows

LSBFLOW is a Fourier-Chebyshev pseudo-spectral collocated Navier-Stokes solver written in Fortran, and parallelized with MPI for use on a Linux cluster. Its purpose is to simulation boundary layer flows over flat plates with pressure gradients imposed through suction and blowing at the wall or ceiling of the domain. Large eddy simulation capability with the following subgrid-scale models is included: Smagorinsky, Dynamic Smagorinsky, sigma, interscale energy transfer, and truncated Navier-Stokes.

For more information on the methods and models used, please look at the numerous helpful comments in the code. I will add a how to guide if anyone is interested - just email me. To compile, just type "make" in a terminal inside the "src" folder. To run, from the "run" folder type "mpiexec -np 4 lsb.x > output.txt" and replace the number "4" by however many mpi processes you'd like to split the domain into.

Pre-requisites:
* gfortran (only minor changes req'd for intel or pg compilers)
* MPI (tested on mpich2 and openmpi
* FFTW (tested on version 3.3.3)
* LAPACK & BLAS (any modern version should work)
