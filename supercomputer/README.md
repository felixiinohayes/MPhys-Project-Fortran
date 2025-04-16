Project overview:
Band-structure computations for BiTeI in the topological insulator phase with degrees of broken periodicity.

DEPENDENCIES:
Compilers:
- mpif90, gfortran

Libraries:
- openmpi, lapack, openblas, arpack (and parpack)

INSTRUCTIONS:
1. Edit ARPACK_LIB variable in the Makefile to arpack lib path.
2. Optional: I believe these are normally system-dependent but if there are any optimisation flags you think we should use, feel free to add - do NOT put floating point optimisations (e.g -ffast-math)
3. Run "make run" and it should compile and run everything.
4. The results will be saved in the "data" directory. Please send us back this directory as well as the "logs" directory.
