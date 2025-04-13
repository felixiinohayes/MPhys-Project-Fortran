ARPACK_LIB = /opt/homebrew/Cellar/arpack/3.9.1_1/lib
NUMPROCS = 6 #135
#2D: nb = 75, NUMPROCS = 90 -> 3.5GB per process

FFLAGS = -fallow-argument-mismatch -ffree-line-length-0 -O3 -cpp
FC = mpif90
LIB = -llapack -lblas -larpack -lparpack -L${ARPACK_LIB}

main1 = cube_mem
main2 = 2D

all: ${main1} ${main2}

mem_usage.o: mem/mem_usage.c
	gcc -c -o mem_usage.o mem/mem_usage.c

${main1}: ${main1}.f90 mem_usage.o
	${FC} ${FFLAGS} -o ${main1} ${main1}.f90 mem_usage.o ${LIB}

${main2}: ${main2}.f90 mem_usage.o
	${FC} ${FFLAGS} -o ${main2} ${main2}.f90 mem_usage.o ${LIB}

run: ${main1} ${main2}
	@echo "Running ${main2}"
	@timestamp=$$(date +"%Y%m%d_%H%M%S"); \
	mpirun -n ${NUMPROCS} ./${main2} > logs/2D_$${timestamp}.log 2>&1
	@echo "Running ${main1} with B=0.05"
	@${MAKE} run-cube B=0.05
	@echo "Running ${main1} with B=0"
	@${MAKE} run-cube B=0.0

run-cube:
	${FC} ${FFLAGS} -DBVAL=${B} -o ${main1} ${main1}.f90 mem_usage.o ${LIB}
	@timestamp=$$(date +"%Y%m%d_%H%M%S"); \
	mpirun -n ${NUMPROCS} ./${main1} > logs/cube_${B}_$${timestamp}.log 2>&1

clean:
	rm -rf *.mod *.csv *.o *.c ${main1} ${main2}
	rm -rf logs/*.log
