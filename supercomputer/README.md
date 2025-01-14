Project overview:
Band-structure computations for BiTeI in the topological insulator phase with degrees of broken periodicity. 

1D.f90: 
 - Periodicity broken in 1 dimension.
  
2D.f90:
 - Periodicity broken in 2 dimensions.

cube.f90:
- Periodicity broken in 3 dimensions.

Dependencies: 
Compilers: 
- mpif90, gfortran
  
Libraries:
- openmpi,lapack, openblas, arpack (and parpack)

Expected behaviour:
- 

INSTRUCTIONS:

1. Ensure arpack is installed
2. Edit Makefile ARPACK_LIB variable to arpack installation path.
3. Edit Makefile NUMPROCS variable.
4. make
5. make run
