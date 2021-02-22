default: ex1
export SLEPC_DIR=/home/costafi/slepc-3.14.2
export PETSC_ARCH=arch-linux-c-debug
export PETSC_DIR=/home/costafi/petsc
include ${SLEPC_DIR}/lib/slepc/conf/slepc_common

ex1: slepc_distributed_matrix_diag.o 
	-${FLINKER} -o ex1 slepc_distributed_matrix_diag.o ${SLEPC_EPS_LIB} 
	${RM} slepc_distributed_matrix_diag.o
	mpirun -np 4 --use-hwthread-cpus  ./ex1
	#gnuplot -c plot.gnuplot
