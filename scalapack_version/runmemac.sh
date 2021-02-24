mpif90 -g distributed_matrix_diag.f90 -L/usr/local/lib -lscalapack -o ex1
mpirun --mca btl_vader_backing_directory /tmp -np 4 -oversubscribe ./ex1
#mpirun -np 2 ./ex1
gnuplot gnuplot/plot_gs.gnu
gnuplot gnuplot/plot_time.gnu 
#open times.pdf
#open gs_H.pdf
#open gs_I.pdf
