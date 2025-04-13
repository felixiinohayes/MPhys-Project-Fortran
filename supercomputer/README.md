Project overview:
Band-structure computations for BiTeI in the topological insulator phase with degrees of broken periodicity. 

2D.f90:
 - Lattice periodicity broken in 2 dimensions.

cube.f90:
- Lattice periodicity broken in 3 dimensions (real space).

Dependencies: 
Compilers: 
- mpif90, gfortran
  
Libraries:
- openmpi,lapack, openblas, arpack (and parpack)

Configuration:
1. Ensure arpack (and therefore parpack) is installed
2. Edit Makefile ARPACK_LIB variable to arpack lib path.
3. I believe these are normally system-dependent but if there are any optmisation flags you think we should use, feel free to add - although numerical precision is of a decent priority for the arnoldi iteration so do NOT put floating point optimisations (e.g -ffast-math)

Execution:
1. make
2. make run

Expected behaviour:
- make : Makefile compiles both 2D.f90 and cube.f90 into executables
- make run: executes mpirun on both programs.
- The cube program may finish without runtime errors but have a parpack error which you should see in console (in this case something is obviously wrong but this should NOT happen).
