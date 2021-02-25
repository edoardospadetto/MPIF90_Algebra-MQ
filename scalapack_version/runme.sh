echo "Compiling.."
mpif90 -g distributed_matrix_diag.F90 -L /home/edoardo/lib/scalapack-2.1.0 -lscalapack  -llapack -lblas  -o ./bin/ex1
rm *mod
echo "Run .."
mpirun -np 4  ./bin/ex1



echo "Plot.."
gnuplot "./plot_gs.gnu"
gnuplot "./plot_time.gnu"


