echo "Compiling.."
mpif90 -g distributed_matrix_diag.F90 -L /home/costafi/scalapack-2.1.0 -lscalapack  -llapack -lblas  -o ex1
rm *mod
echo "Run .."
mpirun -np 4 --use-hwthread-cpus  ./ex1