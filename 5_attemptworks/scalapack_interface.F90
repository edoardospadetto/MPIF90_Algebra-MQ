!##################################################
!Module to handle matrices in a distribute fashion
!#################################################
module scalapack_interface
use matrix_interface
implicit none
contains

	
	!compute kroeneker product of a matrix with and Identity matrices of the different dimensions
	! I idb x A x I ida
	!exploits properties of identitiy matrices.
	!It does not work.
	function getbigmatdiff(idb, inputmat,ida) result(bigmat)
	integer ::  N, ii,jj, kk, idim, ida, idb
	double complex , dimension(:,:) :: inputmat
	double complex, dimension(size(inputmat,dim = 1)*idb*ida,size(inputmat,dim = 1)*ida*idb) :: bigmat
	idim = size(inputmat,dim = 1)
	bigmat = 0.0
	do ii = 0, idb-1
	    do jj = 1, idim*ida
		do kk = 1,idim*ida
		   if (jj-(jj/ida) .eq. kk-(kk/ida) ) then
		    bigmat(jj+(idim*ida)*ii, kk+(idim*ida)*ii) = inputmat(jj/ida, kk/ida ) 
		   end if 
		  
		end do
	    end do
	end do

	end function

!input M , --> A, desca
subroutine build_matrix(M,iam, A, desca)
	use mpi
	implicit none
	
	!input
	double complex , dimension(:,:), intent(INOUT) :: M, A
	integer, intent(INOUT) ::  iam
	integer, dimension(:), intent(INOUT) ::  desca
	!vars 
	integer :: sizeG, lda
	integer :: ierr
	integer :: ii, jj 
	
	
	
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
		print*, "------------------"
		call MPI_BCAST(m,sizeg*sizeg, MPI_DOUBLE_COMPLEX, 0 ,MPI_COMM_WORLD,ierr)
	else 
		call MPI_BCAST(m,sizeg*sizeg, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD,ierr )
	end if
		
	!Print to show matrices
	!do ii =0 , nprows*npcols
    	!if(iam .eq. ii) then
    	!print*, "******"
    	!call pzm(M)
    	!end if 
    	!call mpi_barrier(MPI_COMM_WORLD, ierr) 
    	!end do

		
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
		
	
end subroutine 

!simple interface to print distributed matrix
subroutine printmat(A, descA)
	integer , dimension(:):: desca 
	double complex, dimensioN(:,:):: A
	double complex, dimension(30) :: work
	CALL PZLAPRNT( descA(3),descA(4), A, 1, 1, DESCA, 0, 0, 'A', 6, WORK )
end subroutine 


!diagonalize matrix
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

        
        lwork = 1000
        lrwork = 4*size(w) -2 ! dopo proviamo 

        allocate(work(lwork))

        allocate(rwork(lrwork))

        CALL pzheev( 'V', 'U', size(w) , A, 1,1, DESCA, W, Z, 1, 1, DESCZ, WORK, LWORK, RWORK, LRWORK, INFO )

        deallocate(work)
        deallocate(rwork)


end subroutine

!matmul
subroutine dmatmul(A,descA,B,descB,C,descC)
	use mpi 
	implicit none
        double complex, dimension(:,:), intent(INOUT) :: A, B, C
        integer, dimension(:), intent(INOUT) :: descA, descB, descC 
 
        
        if (descA(4) .ne. descB(3)) then 
                write(*,*) "invalid matrix dimensions"
                STOP
        end if
        if (descC(4) .ne. descB(4))  then 
                 write(*,*) "invalid matrix dimensions"
                 STOP
        end if
        if (descC(3) .ne. descA(3)) then 
                 write(*,*) "invalid matrix dimensions"
                 STOP
        end if
       
 	!do ii =0 , 3
    	!if(iam .eq. ii) then
    	!print*, iam , ia, ja, ib ,jb
    	!end if 
    	!call mpi_barrier(MPI_COMM_WORLD, info) 
    	!end do
    	
 	C = dcmplx(0.d0,0.d0)
        call pzgemm('N', 'N',descA(3), descB(4), descA(4), dcmplx(1.0,0.0), & 
        A, 1, 1, descA, B,  1, 1, descB, dcmplx(1.0,0.0), C, 1,1, descC)
       
     
        
        end subroutine


end module scalapack_interface




