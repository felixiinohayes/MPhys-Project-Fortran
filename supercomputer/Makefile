ARPACK_LIB = /opt/homebrew/Cellar/arpack/3.9.1_1/lib

FFLAGS = -fallow-argument-mismatch -ffree-line-length-0 -O3
FC = mpif90
LIB = -llapack -lblas -ld_classic -larpack -lparpack -L${ARPACK_LIB}

main1 = cube
main2 = 2D
numprocs= 1

all: ${main1} ${main2}

${main1}: ${main1}.f90
	${FC} ${FFLAGS} -o ${main1} ${main1}.f90 ${LIB}

${main2}: ${main2}.f90
	${FC} ${FFLAGS} -o ${main2} ${main2}.f90 ${LIB}

run: ${main1} ${main2}
	mpirun -n ${numprocs} ./${main1}
	mpirun -n ${numprocs} ./${main2}

clean:
	rm -rf *.mod *.o *.c ${main1} ${main2}
