#include "modules/debugqic.F90"
#include "modules/matrixqic.F90"
#include "schandler.F90"



program test_scalapack
    use matrixqic
    use schandler
    implicit none
    integer :: dim
    integer :: ii,jj


    double complex, dimension(:,:), allocatable :: m,z

   
    
    
    external blacs_pinfo, blacs_setup
    !print *, 'Insert the dimensions of the matrix: '
    dim = 4

   
    
    allocate(m(dim,dim))
    allocate(z(dim,dim))
    m = rghcm(dim) 
    call pzm(m)
 
    

	
	 
   call ddzm(dim,m,z)


end program test_scalapack
