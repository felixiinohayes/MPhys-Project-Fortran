ARPACK_LIB = /opt/homebrew/Cellar/arpack/3.9.1_1/lib

FFLAGS = -fallow-argument-mismatch -ffree-line-length-0 -O3 -cpp
FC = mpif90
LIB = -llapack -lblas -larpack -lparpack -L${ARPACK_LIB}

SRC = src
BIN = bin
main1 = cube_shift
main2 = 2D

all: ${BIN}/${main1} ${BIN}/${main2}

${BIN}/${main1}: ${SRC}/${main1}.f90
	${FC} ${FFLAGS} -o $@ $< ${LIB}

${BIN}/${main2}: ${SRC}/${main2}.f90
	${FC} ${FFLAGS} -o $@ $< ${LIB}

clean:
	rm -rf *.mod *.o ${BIN}/* logs/*