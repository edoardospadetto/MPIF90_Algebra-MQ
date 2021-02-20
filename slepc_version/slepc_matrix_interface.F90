#include <slepc/finclude/slepceps.h>

module matrix_interface
    use slepceps
    contains 

    !compute eigenvalues and eigenvector with zheev,
    subroutine eigz(A,eps,nu,xr,xi,rank,kr,error)
        Mat            A
        EPS            eps
        EPSType        tname
        PetscReal      tol, error
        PetscScalar    kr, ki
        Vec            xr, xi
        PetscInt       nu,i,j
        PetscInt       nev, maxit, its, nconv
        PetscMPIInt    rank
        PetscErrorCode ierr
        PetscInt       ncv,mpd

        nev = 1
        ncv = 2*nev
        mpd = ncv

    !     ** Set operators. In this case, it is a standard eigenvalue problem
        call EPSSetOperators(eps,A,PETSC_NULL_MAT,ierr)
        call EPSSetProblemType(eps,EPS_HEP,ierr)

        call EPSSolve(eps,ierr)
        call EPSGetIterationNumber(eps,its,ierr)
        call EPSGetType(eps,tname,ierr)
        call EPSGetTolerances(eps,tol,maxit,ierr)

        !if (rank .eq. 0) then
        !write(*,*) ' Number of iterations of the method:',its
        !write(*,*) ' Solution method:', tname
        !write(*,*) ' Number of requested eigenvalues:',nev
        !write(*,*) ' Stopping condition: tol=',tol,', maxit=',maxit
        !endif

    !     ** Get number of converged eigenpairs
        call EPSGetConverged(eps,nconv,ierr)

    !     ** Display eigenvalues and relative errors
        if (nconv.gt.0) then
            i = 0
            j = 1
            call EPSGetEigenpair(eps,i,kr,ki,xr,xi,ierr)

            if(PetscRealPart(kr) .gt. 0) then
                call EPSGetEigenpair(eps,j,kr,ki,xr,xi,ierr)
            end if

            !Compute the relative error associated to each eigenpair
            call EPSComputeError(eps,i,EPS_ERROR_RELATIVE,error,ierr)

        else
            print *, 'Error'
        endif

    end subroutine eigz

end module 