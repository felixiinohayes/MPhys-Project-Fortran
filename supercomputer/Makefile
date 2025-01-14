FFLAGS = -fallow-argument-mismatch -ffree-line-length-0 -O3
FC = mpif90
LIB = -L/Users/felixiinohayes/Programming/arpack-ng/build/lib -llapack -lblas -lparpack

main1 = cube
main2 = 2D

all: ${main1} ${main2}

${main1}: ${main1}.f90
	${FC} ${FFLAGS} -o ${main1} ${main1}.f90 ${LIB}

${main2}: ${main2}.f90
	${FC} ${FFLAGS} -o ${main2} ${main2}.f90 ${LIB}

run: ${main1} ${main2}
	./${main1}
	./${main2}

clean:
	rm -rf *.mod *.o *.c ${main1} ${main2}
