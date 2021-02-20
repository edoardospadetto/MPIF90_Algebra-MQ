mpif90 -g distributed_matrix_diag.F90 -L /home/edoardo/lib/scalapack-2.1.0 -lscalapack  -llapack -lblas  -o ex1
mpirun -np 4  ./ex1