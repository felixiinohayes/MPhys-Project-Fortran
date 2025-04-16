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

1. Edit the `ARPACK_LIB` variable in the `Makefile` to point to the correct ARPACK library path on your system.
2. *(Optional)* You may add system-specific optimisation flags if desired. **Do not** include floating-point optimisations such as `-ffast-math`, as they may compromise numerical accuracy.
3. Run `make run` to compile and execute the code.
4. Results will be saved in the `data` directory. Please send us back both the `data` and `logs` directories.
