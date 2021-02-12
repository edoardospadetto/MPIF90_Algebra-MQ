#include "modules/debug_module.F90"
#include "modules/matrix_interface.F90"
#include "modules/scalapack_interface.F90"
#include "modules/hamiltonians.F90"



program test_scalapack
    use mpi
    use matrix_interface
    use scalapack_interface
    use hamiltonians
    
    implicit none
    !) BLACS
    integer :: iam, nprocs, nprow, npcol, myrow, mycol, context
    !)  MATRICES
    integer :: N, sizeg
    PARAMETER (N= 4) 
    PARAMETER (sizeg= 2**N) 
    double complex, dimension(sizeg, sizeg) :: M,H,L
    double precision , dimension(sizeg) :: eigvaltest, w
    
    integer :: maxn
    PARAMETER  (MAXN = 100)
    integer :: lda
    PARAMETER (lda = maxn)
    
    double complex, dimension(lda,lda) :: A,Z
    integer, dimension(9) :: desca,  descz
    integer :: info , nb
    
    
    real*8 :: lambda
    
    
    nb = 4
    lambda = 1

    
  
    
	!1) INITIALIZE BLACS
	NPROW= 2
	NPCOL=2
	
	CALL BLACS_PINFO( IAM, NPROCS )
	IF( ( NPROCS.LT.1 ) ) THEN
		CALL BLACS_SETUP( IAM, NPROW*NPCOL )
	END IF
	
	!###############
	
	CALL BLACS_GET( -1, 0, CONTEXT )
	CALL BLACS_GRIDINIT( CONTEXT, 'R', NPROW, NPCOL )
	CALL BLACS_GRIDINFO( CONTEXT, NPROW, NPCOL, MYROW, MYCOL )
	
    	!##############
    	
    	
    	IF(MYROW .eq. -1) THEN
		WRITE(*,*) " ERROR, blacs context not valid "
		CALL BLACS_EXIT(0) 
		STOP
	END IF 
	
	!BUILD GLOBAL MATRIX
	!m = rghcm(sizeg)
	
	
	A= dcmplx(0.d0,0.d0)
	
	
	!Prepare A,B,and C
	call DESCINIT( DESCA, sizeg, sizeg, nb, nb, 0, 0, context, lda, info)
	call transverse_field_ising_model_hamiltonian(context, lambda,N, A, descA)
	
	
    
	
	!IF (IAM .eq. 0 ) then
	!	print*, "EIGENVALUES" 
		!L = matmul(M,H)
		!call pzm (m+h)
	        !call eigz(L, size(L, dim=1), eigvaltest)
		!print*, eigvaltest
	!END IF
    
    !print *, 'Insert the dimensions of the matrix: '
   
   
    call DESCINIT( DESCZ, sizeg, sizeg, nb, nb, 0, 0, context, lda, info)
    
    call ddzm(A,descA, Z, descz, W)
    
    IF (IAM .eq. 0 ) then
		print*, "EIGENVALUES"      
		print*, w
    END IF
    
     
    
     
    
 
    
	CALL BLACS_GRIDEXIT( CONTEXT )
	CALL BLACS_EXIT(0)



	
	 
   

 

end program test_scalapack
