!##################################################
!Module to handle matrices in a distribute fashion
!#################################################
module scalapack_interface
use matrix_interface
use debug_module
implicit none
contains


!------------------------------------------------------------------------------------------
!Compute the sum for distribute matrices, safely.
subroutine dsum(A, descA, B, descB, C,descC)
	integer, dimension(:):: descA, descB , descC
	double complex, dimension(:,:):: A,B,C
	double complex :: alpha, beta
	
	integer ii , jj
	
	call breakifn("Wrong dimensions" , all( descA(3:4) .eq. descB(3:4)) , .true.)
	call breakifn("Wrong dimensions" , all( descA(3:4) .eq. descC(3:4)) , .true.)
	
	if(all((descA .eq. descB) .and. (descB .eq. descC) )) then 
		C = A+B
		
	!else
	!	do ii = 1 , descA(3)
	!		do jj = 1, descA(4)
	!		call pzelget('A','I',A,ii, jj ,descA, alpha)
	!		call pzelget('A','I',B,ii, jj ,descB, beta)
	!		 
	!		
	!		call pzelset(C, ii ,jj , descC, alpha + beta)
	!		end do 
	!	end do 
	end if
	

end subroutine


!------------------------------------------------------------------------------------------
!input M , --> A, desca
subroutine build_matrix(M,iam, A, desca)
	use mpi
	use matrix_interface
	implicit none
	
	!input
	double complex , dimension(:,:), intent(INOUT) :: M, A
	integer, intent(INOUT) ::  iam
	integer, dimension(:), intent(INOUT) ::  desca
	!vars 
	integer :: sizeG 
	integer :: ierr
	integer :: ii, jj 
	
	
	
	!Check M
	if( size(M, dim=1) .eq. size(M, dim = 2)) then 
		sizeG =   size(M, dim=1)
	else 
		print*, "not a square matrix "
		STOP
	end if
	

	
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

!------------------------------------------------------------------------------------------
!simple interface to print distributed matrix
subroutine printmat(A, descA)
	integer , dimension(:):: desca 
	double complex, dimensioN(:,:):: A
	double complex, dimension(30) :: work
	CALL PZLAPRNT( descA(3),descA(4), A, 1, 1, DESCA, 0, 0, 'A', 6, WORK )
end subroutine 

!------------------------------------------------------------------------------------------
!diagonalize matrix
subroutine ddzm(A, desca, Z, descz, W)
        
        use matrix_interface
        implicit none
        
        !input 
        double complex , dimension(:,:), intent(INOUT) :: A, Z
        integer , dimension(:), intent(INOUT):: desca, descz
        double precision, dimension(:), intent(INOUT)::W
        double complex, dimension(:), allocatable ::  rwork
        integer :: lwork, lrwork
        integer :: info
        double precision, dimension(:), allocatable :: work

        
        lwork = 100000
        lrwork = 4*size(w) -2 ! dopo proviamo 

        allocate(work(lwork))

        allocate(rwork(lrwork))

        CALL pzheev( 'V', 'U', size(w) , A, 1,1, DESCA, W, Z, 1, 1, DESCZ, WORK, LWORK, RWORK, LRWORK, INFO )

        deallocate(work)
        deallocate(rwork)


end subroutine
!------------------------------------------------------------------------------------------

!diagonalize matrix
subroutine ddzm2(A, desca, Z, descz, W)
        
        use matrix_interface
        implicit none
        
        !input 
        double complex , dimension(:,:), intent(INOUT) :: A, Z
        integer , dimension(:), intent(INOUT):: desca, descz
        double precision, dimension(:), intent(INOUT)::W
        double complex, dimension(:), allocatable ::  rwork
        integer :: lwork, lrwork,liwork
        integer :: info
        double precision, dimension(:), allocatable :: work
	integer, dimension(:), allocatable :: iwork
	
        
        lwork = 100000
        lrwork = 30000!4*size(w) -2 ! dopo proviamo 
        liwork = 1000
        allocate(work(lwork))
	allocate(iwork(liwork))
        allocate(rwork(lrwork))

        CALL pzheevr( 'V','A' ,'U', size(w) , A, 1,1, DESCA,&
        -100.d0,100.d0,0,100,size(w),size(w) ,W, Z, 1, 1, DESCZ, WORK, LWORK, &
        	 RWORK, LRWORK,IWORK,LIWORK, INFO )

        deallocate(work)
        deallocate(rwork)


end subroutine
!-----------


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
!-----------------------------------------------------------
!Must be tested
	subroutine distributedvsfull(A,descA, B)
		implicit none
		double complex , dimension(:,:) ::A
		double complex, dimension (:,:), allocatable :: Bprime
		double complex , dimension(:,:) ::B
		integer , dimension(:) :: descA
		integer, dimension(9) :: descB
		double complex :: alpha,tot
		integer :: ii , jj 
		
		descB = descA
		allocate(Bprime(size(A, dim = 1), size(A, dim= 2)))
		do ii = 1 , descA(3)
			do jj = 1, descA(4)
			
			call pzelset(Bprime, ii ,jj , descB, B(ii,jj))
			
			end do 
		end do 
		print*, sum(A-Bprime)
		
	end subroutine


!------------------------------------------------------
 !compute kroeneker product of a matrix with and Identity matrices of the same dimensions
! prod{1_i-1} I_{sizeA} x A x prod{i+1_N} I_{sizeA}
!exploits properties of identitiy matrices.
!The result is stored in a scalapack distributed matrix
!input mat is not distributed
!MUST be tested

subroutine getbigmat(inputmat, index, N ,bigmat,descBigMat)
integer :: index, N, ii,jj, kk1, kk2,idim
integer :: helpidx(2), blocksize, dimleft, dimright
integer , dimension(:) :: descBigMat
double complex , dimension(:,:) :: inputmat
double complex, dimension( size(inputmat,dim = 1)**N, size(inputmat,dim = 1)**N) :: bigmat

call breakifn("Invalid rows number" ,(size(inputmat,dim = 1)**N .eq. size(bigmat,dim =1)), .true.)
call breakifn("Invalid columns number" ,( size(inputmat,dim = 1)**N .eq. size(bigmat,dim =2)), .true.)

idim = size(inputmat,dim = 1)
blocksize = (idim**(N-index))


dimright= idim**(index-1)

bigmat = dcmplx(0.d0,0.d0)
do ii = 1, idim
    do jj = 1, idim
        do kk1 = 0, dimright-1
            helpidx(1)= (ii)+idim*kk1
            helpidx(2)= (jj)+idim*kk1
            do kk2 = 1, blocksize
                !bigmat((helpidx(1)-1)*blocksize+kk2 ,(helpidx(2)-1)*blocksize+kk2)  &
                !   = inputmat(ii,jj)
                CALL PZELSET( bigmat, (helpidx(1)-1)*blocksize+kk2, (helpidx(2)-1)*blocksize+kk2, descBigMat, inputmat(ii,jj) )
            end do
        end do
    end do
end do
end subroutine

end module scalapack_interface




