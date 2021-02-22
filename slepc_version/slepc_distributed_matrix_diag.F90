#include "slepc_hamiltonians.F90"
#include "slepc_matrix_interface.F90"
#include <slepc/finclude/slepceps.h>

program main

      use slepceps
      use hamiltonians
      use matrix_interface
      implicit none


!  Variables:
!     A     operator matrix
!     eps   eigenproblem solver context

      Mat            A
      EPS            eps
      Vec            xr, xi
      PetscReal      error,time
      PetscScalar    kr
      PetscInt       nu,n
      PetscMPIInt    rank
      PetscErrorCode ierr
      PetscScalar    couplings(3)
      PetscScalar    lambda

      integer :: ii,jj,l
      real :: from,to

      call SlepcInitialize(PETSC_NULL_CHARACTER,ierr)
      call MPI_Comm_rank(PETSC_COMM_WORLD,rank,ierr)

      nu = 7
      n = 2**nu
      couplings = (/1.d0,1.d0,1.d0/)

      l = 20
      from = 0
      to = 3.d0

      if(rank .eq. 0 ) then
            open(unit = 22, file="data.txt", action="write" , status="old")
      end if


      do ii = 2,10
            do jj = 0,l

                  nu = ii
                  n = 2**nu

                  ! Create matrix and set values
                  call MatCreate(PETSC_COMM_WORLD,A,ierr)
                  call MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,n,n,ierr)
                  call MatSetUp(A,ierr)
            
                  lambda = from + (to-from)/l*jj 
                  !lambda = 0.d0
                  call heisenbergmodel_hamiltonian(A,nu,lambda,couplings,rank)
                  call MatCreateVecs(A,xr,xi,ierr)

                  !if( jj .eq. int(l/2) ) then
                  !      call MatView(A,PETSC_VIEWER_STDOUT_WORLD,ierr)
                  !end if

                  ! Diagonalize
                  kr = 0.d0
                  call EPSCreate(PETSC_COMM_WORLD,eps,ierr)
                  call eigz(A,eps,nu,xr,xi,rank,kr,error,time)

                  if (rank .eq. 0) then
                        print*, "N :", nu," ---- lambda :",PetscRealPart(lambda)
                        write(22,*) nu,PetscRealPart(lambda),PetscRealPart(kr)/(nu-1),error,time
                  endif

                  call MatDestroy(A,ierr)
                  call VecDestroy(xr,ierr)
                  call VecDestroy(xi,ierr)
                  call EPSDestroy(eps,ierr)

            end do
            !write(22,*)
      end do
      
      call SlepcFinalize(ierr)

end