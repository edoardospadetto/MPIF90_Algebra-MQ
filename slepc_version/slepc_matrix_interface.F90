#include <slepc/finclude/slepceps.h>

module matrix_interface
    use slepceps
    contains 

    !compute eigenvalues and eigenvector with zheev,
    subroutine eigz(A,eps,nu,xr,xi,rank,kr,error,time)

        Mat            A
        EPS            eps
        EPSType        tname
        PetscReal      tol, error, t1, t2, t3, time
        PetscScalar    kr, ki
        Vec            xr, xi
        PetscInt       nu,i,nev
        PetscInt       maxit, its, nconv
        PetscMPIInt    rank
        PetscErrorCode ierr


    !   ** Set operators
        call EPSSetOperators(eps,A,PETSC_NULL_MAT,ierr)
        call EPSSetProblemType(eps,EPS_HEP,ierr)
   
        call EPSSetWhichEigenpairs(eps,EPS_SMALLEST_REAL,ierr)
   
        call PetscTime(t1,ierr)
        call EPSSetUp(eps,ierr)
        call PetscTime(t2,ierr)
        call EPSSolve(eps,ierr)
        !if (rank .eq. 0) then
        !        write(*,*) 'Error?', ierr
        !end if
        call PetscTime(t3,ierr)
        time = t3 - t2

        call EPSGetDimensions(eps,nev,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,ierr)
        call EPSGetIterationNumber(eps,its,ierr)
        call EPSGetType(eps,tname,ierr)
        call EPSGetTolerances(eps,tol,maxit,ierr)

        !if (rank .eq. 0) then
        !write(*,*) ' Number of iterations of the method:',its
        !write(*,*) ' Solution method:', tname
        !write(*,*) ' Number of requested eigenvalues:',nev
        !write(*,*) ' Stopping condition: tol=',tol,', maxit=',maxit
        !endif

    !   ** Get number of converged eigenpairs
        call EPSGetConverged(eps,nconv,ierr)

    !   ** Display eigenvalues and relative errors
        if (nconv.gt.0) then

            ! Get eigenvalue
            i = 0
            call EPSGetEigenvalue(eps,i,kr,ki,ierr)

            !Compute the relative error associated to each eigenpair
            call EPSComputeError(eps,i,EPS_ERROR_RELATIVE,error,ierr)

        else
            print *, 'Error'
        endif

    end subroutine eigz

end module 