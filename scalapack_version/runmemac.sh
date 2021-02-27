mpif90 -g distributed_matrix_diag.f90 -L/usr/local/lib -lscalapack -o ./bin/ex1
#mpirun --mca btl_vader_backing_directory /tmp -np 4 -oversubscribe ./bin/ex1 2 1 
mpirun -np 2 ./bin/ex1 2 1 
#gnuplot gnuplot/plot_gs.gnu
#gnuplot gnuplot/plot_time.gnu 
#open times.pdf
#open gs_H.pdf
#open gs_I.pdf
