include "modules/debug_module.F90"
include "modules/matrix_interface.F90"
include "modules/scalapack_interface.F90"
include "modules/hamiltonians.F90"
include "modules/utils.F90"


program test_scalapack
    
    use mpi
    use matrix_interface
    use scalapack_interface
    use hamiltonians
    use utils
    
    implicit none

    integer :: iam, nprocs, nprow, npcol, myrow, mycol, context
    integer :: N, sizeg
    !PARAMETER (N = 8) 
    !PARAMETER (sizeg = 2**N) 
    double complex, dimension(:,:), allocatable :: M, H, L
    double precision, dimension(:), allocatable :: eigvaltest, w
    integer :: maxn
    PARAMETER (MAXN = 500)
    integer :: lda, ii
    PARAMETER (lda = maxn)
    double complex, dimension(lda,lda) :: A,Z
    integer, dimension(9) :: desca, descz
    integer :: info, nb
    real*8 :: couplings(3)
    real*8 :: lambda
    !complex*16 :: testmatrix(2,2)
    !double complex, dimension(N) :: reselx
    !double complex, dimension(2**N) :: statex
    
    couplings = (/1.d0,1.d0,1.d0/)
    nb = 4
    
    ! ---- INITIALIZE BLACS -------------------------------------------------------
	NPROW = 2
	NPCOL = 2
	
	CALL BLACS_PINFO(IAM, NPROCS)
	IF((NPROCS.LT.1)) THEN
		CALL BLACS_SETUP(IAM, NPROW*NPCOL)
	END IF

	CALL BLACS_GET(-1, 0, CONTEXT)
	CALL BLACS_GRIDINIT(CONTEXT, 'R', NPROW, NPCOL)
	CALL BLACS_GRIDINFO(CONTEXT, NPROW, NPCOL, MYROW, MYCOL)
		
    IF(MYROW .eq. -1) THEN
		WRITE(*,*) "ERROR, blacs context not valid."
		CALL BLACS_EXIT(0) 
		STOP
	END IF 
    ! -----------------------------------------------------------------------------

    ! ---- COMPUTATIONS -----------------------------------------------------------
    do N = 2, 8

        if (iam .eq. 0 ) then
            open(unit=73, file='eig_'//trim(str_i(N))//'.txt', action="write")
        end if

        sizeg = 2**N 
        allocate(M(sizeg, sizeg), H(sizeg, sizeg), L(sizeg, sizeg), eigvaltest(sizeg), w(sizeg))
        
        do ii = 1, 20

            print *, "---- N:", N, "--------------- lambda:", lambda, "----"
            
            !lambda = (ii-20)*(1.d0/19.d0)
            lambda = 0.15*(ii-1)

            A = dcmplx(0.d0,0.d0)
        
            call DESCINIT(DESCA, sizeg, sizeg, nb, nb, 0, 0, context, lda, info)
            !call heisenbergmodel_hamiltonian(context, lambda, couplings,N, A, descA)
            call transverse_field_ising_model_hamiltonian(context, lambda, N, A, descA)
        
            call DESCINIT(DESCZ, sizeg, sizeg, nb, nb, 0, 0, context, lda, info)
            
            call ddzm(A, descA, Z, descz, W)
            
            IF (IAM .eq. 0 ) then
                !print*, ii 
                !write(*, ('B64')) maxval()  
                !print*, statex
                write(73,*) lambda, w(1)/(N-1) !(ii-20)*dble(1.d0/19.d0), w(1)/((N-1))
                !reselx = getavgvalue(statex, 'y', N)
                !print*, reselx
            END IF 

        end do 

        deallocate(M, H, L, eigvaltest, w) 

    end do
    ! -----------------------------------------------------------------------------

    ! ---- EXIT BLACS -------------------------------------------------------------
    CALL BLACS_GRIDEXIT(CONTEXT)
	CALL BLACS_EXIT(0)
    ! -----------------------------------------------------------------------------
	
	close(73) 

end program test_scalapack