mpif90 -g distributed_matrix_diag.f90 -L/usr/local/lib -lscalapack -o ex1

mpirun --mca btl_vader_backing_directory /tmp -np 4 -oversubscribe ./ex1