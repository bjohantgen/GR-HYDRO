# GR-HYDRO
A general relativistic hydrocode written as the final project for Relativity during Fall 2026. The project currently runs Newtonian 1-D radial hydro. The structure is as follows:

### main
This file contains the code that is run during the course of the simulation. The config parameters are put into this file and the output comes from this file.

### physics
This file holds all of the physics functions. This includes wave speeds, energies, transformations of the primitive and conserved variables and flux and source terms.

### simulation
This file contains the algorithm used to solve Euler's equations, the functions that initializes the simulation grid, the functions that apply boundary conditions, the functions that apply the timestep, the function to write the output files, and the function that calculates the performance metrics.

### test_conds
This file holds the initial conditions of famous test problems, that can be applied as a diagnostic after changing the code.