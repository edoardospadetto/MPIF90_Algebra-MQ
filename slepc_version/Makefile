default: ex1
export SLEPC_DIR=/home/edoardo/lib/slepc-3.14.2
export PETSC_ARCH=arch-linux-c-debug
export PETSC_DIR=/home/edoardo/lib/petsc-3.14.4
include ${SLEPC_DIR}/lib/slepc/conf/slepc_common

ex1: slepc_distributed_matrix_diag.o 
	-${FLINKER} -o ex1 slepc_distributed_matrix_diag.o ${SLEPC_EPS_LIB} 
	${RM} slepc_distributed_matrix_diag.o
	mpirun -np 1 --use-hwthread-cpus  ./ex1 1
	mpirun -np 2 --use-hwthread-cpus  ./ex1 2
	mpirun -np 4 --use-hwthread-cpus  ./ex1 4
	mpirun -np 6 --use-hwthread-cpus  ./ex1 6
	mpirun -np 8 --use-hwthread-cpus  ./ex1 8
	
	
	gnuplot -c plot.gnuplot
