ARPACK_LIB = /Users/felixiinohayes/Programming/arpack-ng/build/lib
NUMPROCS = 4

FFLAGS = -fallow-argument-mismatch -ffree-line-length-0 -O3 -cpp
FC = mpif90
LIB = -llapack -lblas -ld_classic -larpack -lparpack -L${ARPACK_LIB}

main1 = cube
main2 = 2D

all: ${main1} ${main2}

${main1}: ${main1}.f90
	${FC} ${FFLAGS} -o ${main1} ${main1}.f90 ${LIB}

${main2}: ${main2}.f90
	${FC} ${FFLAGS} -o ${main2} ${main2}.f90 ${LIB}

run: ${main1} ${main2}
	@echo "Running ${main2} once"
	mpirun -n ${NUMPROCS} ./${main2}
	@echo "Running ${main1} with B=0.0"
	@${MAKE} run-cube B=0.0
	@echo "Running ${main1} with B=0.2"
	@${MAKE} run-cube B=0.2

run-cube:
	${FC} ${FFLAGS} -DB=${B} -o ${main1} ${main1}.f90 ${LIB}
	mpirun -n ${NUMPROCS} ./${main1}

clean:
	rm -rf *.mod *.o *.c ${main1} ${main2}
