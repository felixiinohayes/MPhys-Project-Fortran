# Project Overview

Band-structure computations for BiTeI in the topological insulator phase, incorporating broken periodicity.

## Dependencies

**Compilers:**
- `mpif90`
- `gfortran`

**Libraries:**
- `openmpi`
- `lapack`
- `openblas`
- `arpack` (including `parpack`)

## Instructions

- **Compilation**
1. Edit the `ARPACK_LIB` variable in the `Makefile` to point to the correct ARPACK library path on your system.
2. *(Optional)* You may add system-specific optimisation flags if desired. **Do not** include floating-point optimisations such as `-ffast-math`, as they may compromise numerical accuracy.
3. Run `make all` to produce the binaries in the `bin` directory.

- **Execution**
1. Use the provided SGE job scripts in the `jobs` directory. These scripts compile and execute the binary on the cluster.
2. To submit jobs:

qsub jobs/run_2D.sge
qsub jobs/run_cube_B0.00.sge
qsub jobs/run_cube_B0.05.sge

4. Output logs are saved in the `logs` directory, results are stored in the `data` directory. Please send us back both of these.
