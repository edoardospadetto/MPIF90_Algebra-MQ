#include "modules/debugqic.F90"
#include "modules/matrixqic.F90"
#include "schandler.F90"



program test_scalapack
    use mpi
    use matrixqic
    use schandler
    implicit none
    !) BLACS
    integer :: iam, nprocs, nprow, npcol, myrow, mycol, context
    integer :: dim
    !)  MATRICES
    integer :: sizeg
    PARAMETER (sizeg = 20 )
    double complex, dimension(sizeg, sizeg) :: M
    double precision , dimension(sizeg) :: eigvaltest, w
    
    integer :: maxn
    PARAMETER  (MAXN = 100)
    integer :: lda
    PARAMETER (lda = maxn)
    
    double complex, dimension(lda,lda) :: A,Z
    integer, dimension(9) :: desca, descz
    integer :: info 
    !)
    
  
    
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
	m = rghcm(sizeg)
	call  build_matrix(M,iam,nprow, npcol, myrow, mycol, context, A, desca)
	
	IF (IAM .eq. 0 ) then
		print*, "EIGENVALUES" 
	        call eigz(M, size(M, dim=1), eigvaltest)
		print*, eigvaltest
	END IF
    
    !print *, 'Insert the dimensions of the matrix: '
   
   
    call DESCINIT( DESCZ, sizeg, sizeg, 4, 4, 0, 0, context, lda, info)
    call ddzm(A,desca, Z, descz, W)
    
    IF (IAM .eq. 0 ) then
		print*, "EIGENVALUES"      
		print*, w
    END IF
    
     
    
     
    
 
    
	CALL BLACS_GRIDEXIT( CONTEXT )
	CALL BLACS_EXIT(0)



	
	 
   

 

end program test_scalapack
