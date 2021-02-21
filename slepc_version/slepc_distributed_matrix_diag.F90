#include "slepc_hamiltonians.F90"
#include "slepc_matrix_interface.F90"
#include <slepc/finclude/slepceps.h>

program main

      use slepceps
      use hamiltonians
      use matrix_interface
      implicit none

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Declarations
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!  Variables:
!     A     operator matrix
!     eps   eigenproblem solver context

      Mat            A
      EPS            eps
      Vec            xr, xi
      PetscReal      error
      PetscScalar    kr
      PetscInt       nu,n
      PetscMPIInt    rank
      PetscErrorCode ierr
      PetscScalar    couplings(3)
      PetscScalar    lambda
      PetscBool      flg

      integer :: ii,l
      real :: from,to

      call SlepcInitialize(PETSC_NULL_CHARACTER,ierr)
      call MPI_Comm_rank(PETSC_COMM_WORLD,rank,ierr)

      nu = 7
      n = 2**nu
      couplings = (/1.d0,0.d0,0.d0/)


      l = 40
      from = -2.d0
      to = 2.d0

      call PetscOptionsGetInt(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,  &
    &                        '-n',n,flg,ierr)

      if (rank .eq. 0) then
        write(*,100) nu
      endif
100  format (/'1-D Laplacian Eigenproblem, n =',I3,' units')


      if(rank .eq. 0 ) then
            open(unit = 22, file="data.txt", action="write" , status="old")
      end if

      ! Create matrix and set values

      call MatCreate(PETSC_COMM_WORLD,A,ierr)
      call MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,n,n,ierr)
      call MatSetFromOptions(A,ierr)
      call MatSetUp(A,ierr)

      do ii = 0,l
            
            lambda = from + (to-from)/l*ii     
            call heisenbergmodel_hamiltonian(A,nu,lambda,couplings,rank)
            call MatCreateVecs(A,xr,xi,ierr)

            ! Diagonalize

            call EPSCreate(PETSC_COMM_WORLD,eps,ierr)
            call eigz(A,eps,nu,xr,xi,rank,kr,error)

            if (rank .eq. 0) then
                  write(22,*) PetscRealPart(lambda),PetscRealPart(kr)/(nu-1),error
            endif
        
      end do

!     ** Free work space

      call EPSDestroy(eps,ierr)
      call MatDestroy(A,ierr)
      call VecDestroy(xr,ierr)
      call VecDestroy(xi,ierr)

      call SlepcFinalize(ierr)

end