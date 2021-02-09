
module schandler
use matrixqic 
implicit none
contains

!input M , --> A, desca
subroutine build_matrix(M,iam,nprows, npcols, myrow, mycol, context, A, desca)
	use mpi
	implicit none
	
	!input
	double complex , dimension(:,:), intent(INOUT) :: M, A
	integer, intent(INOUT) ::  iam, nprows, npcols, myrow, mycol, context
	integer, dimension(:), intent(INOUT) ::  desca
	!vars 
	integer :: sizeG, lda
	integer :: ierr
	integer :: nb
	integer :: ii, jj 
	PARAMETER (nb = 4)
	
	
	
	!Check M
	if( size(M, dim=1) .eq. size(M, dim = 2)) then 
		sizeG =   size(M, dim=1)
	else 
		print*, "not a square matrix "
		STOP
	end if
	
	lda = size(A,dim = 1)	
	
	!USE M from the first process
	if(iam .eq. 0) then 
		call MPI_BCAST(m,sizeg*sizeg, MPI_DOUBLE_COMPLEX, 0 ,MPI_COMM_WORLD,ierr)
	end if
		call MPI_BCAST(m,sizeg*sizeg, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD,ierr )

		
	CALL DESCINIT( DESCA, sizeg, sizeg, nb, nb, 0, 0, context, lda, ierr )
	
	if (ierr .ne. 0 ) then 
		print*, "error, building matrix "
		STOP
	end if
	
	!Distribute M 
	DO ii=1,sizeg
		DO jj=1,sizeg
			CALL PZELSET( A, ii, jj, DESCA, M(ii,jj) )	
		ENDDO
	ENDDO
		
	!CALL PZLAPRNT( dim, dim, A, 1, 1, DESCA, 0, 0, 'A', 6, WORK )



end subroutine 


!input A, desca
subroutine ddzm(A, desca, Z, descz, W)
	
	use matrixqic
	implicit none
	
	!input 
	double complex , dimension(:,:), intent(INOUT) :: A, Z
	integer , dimension(:), intent(INOUT):: desca, descz
	double precision, dimension(:), intent(INOUT)::W
	
	
	double complex, dimension(:), allocatable ::  rwork
	integer :: lwork, lrwork
	integer :: info
	double precision, dimension(:), allocatable :: work


	!external blacs_pinfo, blacs_setup

	
	lwork = 1000
	lrwork = 4*size(w) -2 ! dopo proviamo 

	allocate(work(lwork))
	allocate(rwork(lrwork))

	
	CALL pzheev( 'V', 'U', size(w) , A, 1,1, DESCA, W, Z, 1, 1, DESCZ, WORK, LWORK, RWORK, LRWORK, INFO )
		
	
	deallocate(work)
	deallocate(rwork)


end subroutine

!subroutine zdistrmm

!end subroutine

end module schandler 






