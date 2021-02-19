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
	
    ! blacs variables
    integer :: iam, nprocs, nprow, npcol, myrow, mycol, context
    ! system variables
    integer :: N, sizeg
    PARAMETER (N = 7) 
    PARAMETER (sizeg = 2**N) 
    double complex, dimension(sizeg, sizeg) :: M, H, L
    double precision , dimension(sizeg) :: eigvaltest, w
    real*8 :: lambda
    real*8 :: couplings(3)
    !complex*16 :: testmatrix(2,2)
    double complex, dimension(N) :: reselx
    double complex, dimension(2**N) :: statex
    ! scalapack variables
    integer :: maxn
    PARAMETER  (MAXN = 500)
    integer :: lda, ii
    PARAMETER (lda = maxn)
    double complex, dimension(lda,lda) :: A,Z
    integer, dimension(9) :: desca,  descz
    integer :: info , nb
    
    couplings = (/1.d0,1.d0,1.d0/)
    nb = 4
    
    ! ---- initialize BLACS -------------------------------------------------------
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
		WRITE(*,*) "ERROR, blacs context not valid!"
		CALL BLACS_EXIT(0) 
		STOP
    END IF 
    ! -----------------------------------------------------------------------------

    ! ---- Do the computations ----------------------------------------------------
	if (iam .eq. 0) then
		open(unit=22, file="eig_2.txt", action="write" , status="old")
    end if
    
    do ii = 1, 20

        lambda = ii*0.15 !(ii-20)*(1.d0/19.d0)

        A= dcmplx(0.d0, 0.d0)
        
        call DESCINIT(DESCA, sizeg, sizeg, nb, nb, 0, 0, context, lda, info)
        !call heisenbergmodel_hamiltonian(context, lambda,couplings,N, A, descA)
        call transverse_field_ising_model_hamiltonian(context, lambda,N, A, descA)
    
        call DESCINIT(DESCZ, sizeg, sizeg, nb, nb, 0, 0, context, lda, info)
        
        call ddzm(A,descA, Z, descz, W)
        
        !statex = getcol(Z,descZ,1)
        !call printmat(Z,descZ)
        IF (IAM .eq. 0) then
            !print*, ii 
            !write(*, ('B64')) maxval()  
            !print*, statex
            write (22, *) lambda, w(1)/(N-1) !(ii-20)*dble(1.d0/19.d0), w(1)/(N-1)
            !reselx = getavgvalue(statex, 'y', N)
            !print*, reselx 
        END IF

    end do 

    close(22) 

    ! ---- exit BLACS -------------------------------------------------------------
    CALL BLACS_GRIDEXIT(CONTEXT)
    CALL BLACS_EXIT(0)
    ! -----------------------------------------------------------------------------

end program test_scalapack
