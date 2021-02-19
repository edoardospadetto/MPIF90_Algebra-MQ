!##################################################
!Module to handle matrices in a distribute fashion
!#################################################
module scalapack_interface

	use matrix_interface
	use debug_module
	
	implicit none
	
contains

	!Computes the sum for distribute matrices, safely.
	subroutine dsum(A, descA, B, descB, C, descC)
		
		integer, dimension(:) :: descA, descB, descC
		double complex, dimension(:,:) :: A, B, C
		double complex :: alpha, beta
		integer ii, jj
		
		call breakifn("Wrong dimensions", all( descA(3:4) .eq. descB(3:4)), .true.)
		call breakifn("Wrong dimensions", all( descA(3:4) .eq. descC(3:4)), .true.)
		
		if (all((descA .eq. descB) .and. (descB .eq. descC))) then 
			C = A + B	
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

	! Distributed matrix builder
	subroutine build_matrix(M, iam, A, desca)

		use mpi
		use matrix_interface

		implicit none
		
		! input
		double complex , dimension(:,:), intent(INOUT) :: M, A
		integer, intent(INOUT) ::  iam
		integer, dimension(:), intent(INOUT) ::  desca
		! variables
		integer :: sizeG 
		integer :: ierr
		integer :: ii, jj 
		
		! Check if M is square
		if( size(M, dim=1) .eq. size(M, dim=2)) then 
			sizeG = size(M, dim=1)
		else 
			print*, "M is not a square matrix!"
			STOP
		end if
		
		! Use M from the first process
		if(iam .eq. 0) then 
			print*, "IAM = 0!"
			call MPI_BCAST(m,sizeg*sizeg, MPI_DOUBLE_COMPLEX, 0 ,MPI_COMM_WORLD,ierr)
		else 
			call MPI_BCAST(m,sizeg*sizeg, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD,ierr )
		end if
			
		! Print matrices
		!do ii = 0, nprows*npcols
			!if (iam .eq. ii) then
				!print*, "******"
				!call pzm(M)
			!end if 
			!call mpi_barrier(MPI_COMM_WORLD, ierr) 
		!end do

		if (ierr .ne. 0) then 
			print*, "Error building matrix!"
			STOP
		end if
		
		! Distribute M 
		DO ii = 1, sizeg
			DO jj = 1, sizeg
				CALL PZELSET(A, ii, jj, DESCA, M(ii,jj))	
			END DO
		END DO

	end subroutine 

	! Simple interface to print a distributed matrix
	subroutine printmat(A, descA)

		integer, dimension(:) :: desca 
		double complex, dimensioN(:,:) :: A
		double complex, dimension(30) :: work

		CALL PZLAPRNT(descA(3), descA(4), A, 1, 1, DESCA, 0, 0, 'A', 6, WORK)

	end subroutine 

	! Distributed matrix diagonalization 
	subroutine ddzm(A, desca, Z, descz, W)
			
		use matrix_interface

		implicit none
		
		double complex, dimension(:,:), intent(INOUT) :: A, Z
		integer, dimension(:), intent(INOUT) :: desca, descz
		double precision, dimension(:), intent(INOUT) :: W
		double complex, dimension(:), allocatable :: rwork
		integer :: lwork, lrwork
		integer :: info
		double precision, dimension(:), allocatable :: work

		lwork = 100000
		lrwork = 4*size(w) - 2 ! dopo proviamo 

		allocate(work(lwork))
		allocate(rwork(lrwork))

		CALL pzheev('V', 'U', size(w), A, 1, 1, DESCA, W, Z, 1, 1, DESCZ, WORK, LWORK, RWORK, LRWORK, INFO)
	
		deallocate(work)
		deallocate(rwork)

	end subroutine

	! Distributed matrix multiplication
	subroutine dmatmul(A, descA, B, descB, C, descC)
		
		use mpi 

		implicit none

		double complex, dimension(:,:), intent(INOUT) :: A, B, C
		integer, dimension(:), intent(INOUT) :: descA, descB, descC 

		if (descA(4) .ne. descB(3)) then 
			write(*,*) "Invalid matrix dimensions!"
			STOP
		end if
		if (descC(4) .ne. descB(4))  then 
			write(*,*) "Invalid matrix dimensions!"
			STOP
		end if
		if (descC(3) .ne. descA(3)) then 
			write(*,*) "Invalid matrix dimensions!"
			STOP
		end if
			
		C = dcmplx(0.d0,0.d0)
		
		call pzgemm('N', 'N', descA(3), descB(4), descA(4), dcmplx(1.0,0.0), & 
			        & A, 1, 1, descA, B, 1, 1, descB, dcmplx(1.0,0.0), C, 1, 1, descC)
		
	end subroutine

	! -----------------------------------------------------------
	! SUBROUTINES/FUNCTIONS TO BE TESTED OR NOT WORKING
	! -----------------------------------------------------------

	!diagonalize matrix
	!it doesn't work
	subroutine ddzm2(A, desca, Z, descz, W)
			
		use matrix_interface

		implicit none

		double complex, dimension(:,:), intent(INOUT) :: A, Z
		integer, dimension(:), intent(INOUT) :: desca, descz
		double precision, dimension(:), intent(INOUT) :: W
		double complex, dimension(:), allocatable :: rwork
		integer :: lwork, lrwork, liwork
		integer :: info
		double precision, dimension(:), allocatable :: work
		integer, dimension(:), allocatable :: iwork

		lwork = 100000
		lrwork = 30000 !4*size(w) -2 ! dopo proviamo 
		liwork = 1000

		allocate(work(lwork))
		allocate(iwork(liwork))
		allocate(rwork(lrwork))

		CALL pzheevr('V', 'A', 'U', size(w), A, 1, 1, DESCA, -100.d0, 100.d0, 0, 100, size(w), &
		             & size(w), W, Z, 1, 1, DESCZ, WORK, LWORK, RWORK, LRWORK, IWORK, LIWORK, INFO)

		deallocate(work)
		deallocate(rwork)

	end subroutine

	!Must be tested
	subroutine distributedvsfull(A, descA, B)

		implicit none

		double complex , dimension(:,:) :: A
		double complex, dimension (:,:), allocatable :: Bprime
		double complex , dimension(:,:) :: B
		integer , dimension(:) :: descA
		integer, dimension(9) :: descB
		double complex :: alpha, tot
		integer :: ii, jj 
		
		descB = descA

		allocate(Bprime(size(A, dim = 1), size(A, dim= 2)))

		do ii = 1 , descA(3)
			do jj = 1, descA(4)
				call pzelset(Bprime, ii ,jj , descB, B(ii,jj))
			end do 
		end do

		print*, sum(A-Bprime)
		
	end subroutine
		
	! Get column
	function getcol(A, descA, col) result(res)
		
		double complex, dimension(:,:) :: A
		integer:: dimr, dimc, descA(9), col, ii
		double complex :: alpha
		!parameter(dimr = descA(3), dimc = descA(4))
		double complex, dimension(descA(3)) :: res
		
		do ii = 1, descA(3)
			call pzelget('A', ' ', res(ii), A, ii, col, descA)
		end do
		
		return

	end function

end module scalapack_interface