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
    integer :: N, sizeg, lda, ii, info, nb, iter 
    !PARAMETER (MAXN = 1024)
    !PARAMETER (lda = maxn)
    double complex, dimension(:,:), allocatable :: M, H, L
    double precision, dimension(:), allocatable :: eigvaltest, w
    double complex, dimension(:,:), allocatable :: A, Z
    integer, dimension(9) :: desca, descz
    real*8 :: couplings(3), lambda, start, finish, start_all, finish_all 
    character(1) :: which_model
    !double complex, dimension(:), allocatable:: firsteigenstate 
  
    which_model = 'H'
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
    if (iam .eq. 0) then 
        open(unit=22, file='times.txt', action="write")
    end if 
    
    if(iam == 0) then
        start_all = MPI_Wtime()
    end if 

    do N = 2, 10

        if (iam .eq. 0) then 
            open(unit=73, file='./results/eig_'//trim(which_model)//'_'//trim(str_i(N))//'.txt', action="write")
        end if 

        sizeg = 2**N 
        allocate(M(sizeg,sizeg), H(sizeg,sizeg), L(sizeg,sizeg), eigvaltest(sizeg), w(sizeg))!, firsteigenstate(sizeg))
        
        do ii = 1, 21

            lambda = 0.15*(ii-1)

            if (iam .eq. 0) then 
                print *, "---- N:", N, "--------------- lambda:", lambda, "----"
            end if 

            lda = 2**N
            allocate(A(lda,lda))
            allocate(Z(lda,lda))

            A = dcmplx(0.d0,0.d0)

            call DESCINIT(DESCA, sizeg, sizeg, nb, nb, 0, 0, context, lda, info)

            if (which_model .eq. 'H') then 
                call heisenbergmodel_hamiltonian(context, lambda, couplings,N, A, descA, lda)
            else if (which_model .eq. 'I') then 
                call transverse_field_ising_model_hamiltonian(context, lambda, N, A, descA)
            end if
        
            call DESCINIT(DESCZ, sizeg, sizeg, nb, nb, 0, 0, context, lda, info)
            
            if (ii == 10 .and. iam == 0) then
                start = MPI_Wtime()
            end if
            
            !call ddzm2(A, descA, Z, descz, W)
            call ddzm2(A, descA, Z, descz, W)

            if (ii == 10 .and. iam == 0) then
                finish = MPI_Wtime()
            end if 
            

            if (iam .eq. 0) then
                write(73,*) lambda, w(1)/(N-1) 
                !write(73,*) lambda*(real(N)/real(N-1)), w(1)/(N-1) 
            end if 
            
            !####temp 
            !firsteigenstate =  getcol(Z,descZ,1)
            !    if (iam .eq. 0) then
            !	print*, firsteigenstate
            !end if 
            !#######

            deallocate(A, Z)

        end do 

        if (iam .eq. 0) then 
            write(22,*) N, finish - start 
        end if

        deallocate(M, H, L, eigvaltest, w) 

    end do

    if (iam == 0) then
        finish_all = MPI_Wtime()
        print *, 'ELAPSED TIME:', finish_all-start_all
    end if 

    if (iam .eq. 0) then 
        close(22)
        close(73)  
    end if 
    ! -----------------------------------------------------------------------------

    ! ---- EXIT BLACS -------------------------------------------------------------
    CALL BLACS_GRIDEXIT(CONTEXT)
    CALL BLACS_EXIT(0)
    ! -----------------------------------------------------------------------------

end program test_scalapack