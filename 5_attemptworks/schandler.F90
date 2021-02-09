
module schandler
use matrixqic 
implicit none
contains

subroutine ddzm(M)
	use mpi
	use matrixqic
	implicit none
	

	double complex, dimension(:,:), intent(inout) :: m
	double complex, dimension(:), allocatable ::  rwork
	double complex, dimension(:,:) , allocatable :: A, z
	integer :: ii,jj 
	integer ia, ja, iz,jz,dlen
	integer LNCOLSA, LNROWSA
	integer :: lwork, lrwork
	integer :: info
	integer :: dim

	double precision, dimension(:), allocatable :: work
	double precision, dimension(:), allocatable :: eigvl, eigvltest
	integer, dimension(:), allocatable :: desca, descz
	INTEGER :: CONTEXT, I, IAM, MYCOL, MYROW, N, NB, NPCOL, NPROCS, NPROW,ierr


	external blacs_pinfo, blacs_setup
	integer numroc



	!CHECK if input matrix is square
	if (size(M,dim =1) .eq. size(M,dim =1)) then 
		dim =  size(M,dim =1)
	
	
	else if (IAM .eq. 0) then
		write(*,*) "Error" 
		STOP
	end if
	
	!################
	
	
	!###############
	ia = 1
	ja = 1
	iz = 1
	jz = 1
	dlen = 50
	
	lwork = 1000
	lrwork = 2*dim + 2*dim -2

	!Set up the problem

	NB = 2
	NPROW = 2
	NPCOL = 2
	MYROW = -1

	
	
	! Calculate the number of rows, LNROWSA, and the number of columns, LNCOLSA, for the local matrices A:
	! Calculating the block sizes:
	!print*, NPROW, NPCOL
	!CALL BLOCKSET( NB,NB, N, NPROW, NPCOL )

	
	

	
	

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
	print*,"number of processes",  NPROCS
	print*, IAM
	!###############



	! Initialize a single BLACS context

	CALL BLACS_GET( -1, 0, CONTEXT )
	CALL BLACS_GRIDINIT( CONTEXT, 'R', NPROW, NPCOL )
	CALL BLACS_GRIDINFO( CONTEXT, NPROW, NPCOL, MYROW, MYCOL )
	
	!########################
	!TEST Generate the matrix and share it
	PRINT*, "---------------------------------------"
	if(iam .eq. 0) then 
        m = rghcm(dim)
	call MPI_BCAST(m,dim*dim, MPI_DOUBLE_COMPLEX, 0 ,MPI_COMM_WORLD,ierr )
		 
	end if
	
	call MPI_BCAST(m,dim*dim, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD,ierr )
	!print*, m(:,dim)
	
	
	!######################
	PRINT*,"---------------------------------------"


	IF(MYROW .eq. -1) THEN
		WRITE(*,*) " ERROR, blacs context not valid "
		CALL BLACS_EXIT(0) 
		STOP
	ELSE
		!Set local matrices
		LNROWSA = NUMROC(dim,NB,MYROW,0,NPROW)
		LNCOLSA = NUMROC(dim,NB,MYCOL,0,NPCOL)
		

		ALLOCATE(A(LNROWSA,LNCOLSA),Z(LNROWSA,LNCOLSA))
		A = 0.0d0
		
		CALL DESCINIT( DESCA, dim, dim, NB, NB, 0, 0, CONTEXT, LNROWSA, INFO )
		CALL DESCINIT( DESCZ, dim, dim, NB, NB, 0, 0, CONTEXT, LNROWSA, INFO )
		
		DO ii=1,dim
		DO jj=1,dim
		
			CALL PZELSET( A, ii, jj, DESCA, M(ii,jj) )	
			
		ENDDO
		ENDDO
		
		
		!CALL PZLAPRNT( dim, dim, A, 1, 1, DESCA, 0, 0, 'A', 6, WORK )
		!DIAGONALIZE 
		
		CALL pzheev( 'V', 'U', dim, A, 1,1, DESCA, eigvl, Z, 1, 1, DESCZ, WORK, LWORK, RWORK, LRWORK, INFO )
		
	END IF
	
	IF (IAM .eq. 0 ) then
		print*, "EIGENVALUES" 
		print*, " "
	        allocate(eigvltest(dim))
	        call eigz(M, size(M, dim=1), eigvltest)
		print*, eigvl
		print*, " "
		print*, eigvltest
	END IF
	
	
	CALL BLACS_GRIDEXIT( CONTEXT )
	CALL BLACS_EXIT(0)


	!deallocate(desca)
	!deallocate(descz)
	!deallocate(work)
	!deallocate(rwork)


end subroutine

end module schandler 
