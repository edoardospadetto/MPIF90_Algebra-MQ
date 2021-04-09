echo "Compiling.."
mpif90 -g distributed_matrix_diag.F90 -L /home/edoardo/lib/scalapack-2.1.0 -lscalapack -lblas -llapack  -o ./bin/ex1
rm *mod
echo "Run .."
#mpirun -np 1  ./bin/ex1 1 1 
#mpirun -np 2  ./bin/ex1 1 2 
mpirun -np 4  ./bin/ex1 2 2 
mpirun  -np 6  ./bin/ex1 3 2 
mpirun  -np 8   ./bin/ex1 4 2 



echo "Plot.."
gnuplot "./plot_gs.gnu"
gnuplot "./plot_time.gnu"


