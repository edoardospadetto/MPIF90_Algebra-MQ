mpif90 exampl1.f -L /home/edoardo/Downloads/scalapack-2.1.0 -lscalapack  -llapack -lblas  -L /home/edoardo/Downloads/BLACS/LIB -lblacs -lblacsmpi -lblacsF77init  -o ex1
mpirun -np 6 --oversubscribe ./ex1
