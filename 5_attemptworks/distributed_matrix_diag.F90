#include "modules/debugqic.F90"
#include "modules/matrixqic.F90"
#include "schandler.F90"



program test_scalapack
    use mpi
    use matrixqic
    use schandler
    implicit none
    integer :: dim
    integer :: ii,jj
    integer :: nprocs_mpi, ierr, rank 
    character*(MPI_MAX_PROCESSOR_NAME) :: hostname
    integer :: namesize
    double precision :: eigvals(4) 


    double complex, dimension(:,:), allocatable :: m,z,another
    
	
	
    
    	!######
    	
    
    
    !print *, 'Insert the dimensions of the matrix: '
    dim = 16
    allocate(m(dim,dim))
    allocate(z(dim,dim))
    allocate(another(dim,dim))
   
  
    call ddzm(m)
     
    
     
    
 
    

	
	 
   

 

end program test_scalapack
