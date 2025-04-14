ARPACK_LIB = /opt/homebrew/Cellar/arpack/3.9.1_1/lib
NUMPROCS = 4
#2D: nb = 60, NUMPROCS = 80 

FFLAGS = -fallow-argument-mismatch -ffree-line-length-0 -O3 -cpp
FC = mpif90
LIB = -llapack -lblas -larpack -lparpack -L${ARPACK_LIB}

main1 = cube_mem
main2 = 2D

all: ${main1} ${main2}

${main1}: ${main1}.f90 
	${FC} ${FFLAGS} -o ${main1} ${main1}.f90 ${LIB}

${main2}: ${main2}.f90 
	${FC} ${FFLAGS} -o ${main2} ${main2}.f90 ${LIB}

run: ${main1} ${main2}
#	@echo "Running ${main2}"
#	@timestamp=$$(date +"%Y%m%d_%H%M%S"); \
#	mpirun -n ${NUMPROCS} ./${main2}
#--------------------------------
#---Run B=0.05
	@echo "Running ${main1} with B=0.05 E=0.2"
	@${MAKE} run-cube B=0.05 E=0.2
	@echo "Running ${main1} with B=0.05 E=0.0"
	@${MAKE} run-cube B=0.05 E=0.0
	@echo "Running ${main1} with B=0.05 E=-0.2"
	@${MAKE} run-cube B=0.05 E=-0.2
#--------------------------------
#---Run B=0.0
	@echo "Running ${main1} with B=0 E=0.2"
	@${MAKE} run-cube B=0 E=0.2
	@echo "Running ${main1} with B=0 E=0.0"
	@${MAKE} run-cube B=0 E=0.0
	@echo "Running ${main1} with B=0 E=-0.2"
	@${MAKE} run-cube B=0 E=-0.2
#--------------------------------

run-cube:
	${FC} ${FFLAGS} -DBVAL=${B} -DEVAL=${E} -o ${main1} ${main1}.f90 ${LIB}
	@timestamp=$$(date +"%Y%m%d_%H%M%S"); \
	mpirun -n ${NUMPROCS} ./${main1}

clean:
	rm -rf *.mod *.csv *.o *.c ${main1} ${main2}
	rm -rf logs/*.log
