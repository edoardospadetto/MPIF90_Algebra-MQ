
module schandler
	implicit none
	contains

	subroutine ddzm(dim,m,z)
	    implicit none
	    integer ,intent(in) :: dim
	    double complex, dimension(:,:), intent(in) :: m
	    double complex, dimension(:,:), intent(inout) :: z
	    double complex, dimension(:), allocatable ::  rwork
	    integer :: ii,jj 
	    integer ia, ja, iz,jz,dlen
	    integer :: lwork, lrwork
	    integer :: info, lda
	    
	    double precision, dimension(:), allocatable :: work
	    double precision, dimension(:), allocatable :: eigvl
	    integer, dimension(:), allocatable :: desca, descz
	    INTEGER :: CONTEXT, I, IAM, MYCOL, MYROW, N, NB, NPCOL, NPROCS, NPROW
	    
	  
	    
	    
	    external blacs_pinfo, blacs_setup

	   

	    
	    ia = 1
	    ja = 1
	    iz = 1
	    jz = 1
	    dlen = 50
	    lda = dim
	    lwork = 1000
	    lrwork = 2*dim + 2*dim -2

	    !    Set up the problem
	
	    NB = 1
	    NPROW = 2
	    NPCOL = 2
	    MYROW = -1


	    allocate(eigvl(dim))
	    allocate(desca(dlen))
	    allocate(descz(dlen))
	    allocate(work(lwork))
	    allocate(rwork(lrwork))


	
	  	 
		 
	!     Initialize the BLACS

	      CALL BLACS_PINFO( IAM, NPROCS )
	      IF( ( NPROCS.LT.1 ) ) THEN
		 CALL BLACS_SETUP( IAM, NPROW*NPCOL )
	      END IF

	      print*, IAM

	     ! Initialize a single BLACS context
	  
	      CALL BLACS_GET( -1, 0, CONTEXT )
	      CALL BLACS_GRIDINIT( CONTEXT, 'R', NPROW, NPCOL )
	      CALL BLACS_GRIDINFO( CONTEXT, NPROW, NPCOL, MYROW, MYCOL )
	      
	      IF(MYROW .eq. -1) THEN
			WRITE(*,*) " ERROR, blacs context not valid "
			CALL BLACS_EXIT(0) 
			STOP
	      ELSE
			!INIT Descriptors 
			CALL DESCINIT( DESCA, dim, dim, NB, NB, 0, 0, CONTEXT, LDA, INFO )
			CALL DESCINIT( DESCZ, dim, dim, NB, NB, 0, 0, CONTEXT, LDA, INFO )

			!DEBUG PRINT 
			!CALL PSLAPRNT( dim, dim, M, 1, 1, DESCA, 0, 0, 'A', 6, WORK )
			!DIAGONALIZE 
			CALL pzheev( 'V', 'U', dim, m, 1,1, DESCA, eigvl, Z, 1, 1, DESCZ, WORK, LWORK, RWORK, LRWORK, INFO )

			print*, eigvl
	      END IF
	      	CALL BLACS_EXIT(0)
       
	    deallocate(desca)
	    deallocate(descz)
	    deallocate(work)
	    deallocate(rwork)
	    
	
	end subroutine

end module schandler 
