#include "slepc_hamiltonians.F90"
#include "slepc_matrix_interface.F90"
#include <slepc/finclude/slepceps.h>

program main

      use slepceps
      use hamiltonians
      use matrix_interface
      implicit none


      ! Variables:

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


      integer :: ii,jj,l,nprocs, stat
      real :: from,to
      real*8 :: start_build, finish_build, start, finish 
      character*8 :: args(2)
      
      
       call get_command_argument(1,args(1))
       read(args(1),*,iostat=stat)  nprocs
   
    
	    if (stat .ne. 0 ) then 
	    	print *, "error"
	    	STOP 
	    end if
	    
  
    

      ! Initialize communicator
      call SlepcInitialize(PETSC_NULL_CHARACTER,ierr)
      call MPI_Comm_rank(PETSC_COMM_WORLD,rank,ierr)

      ! Number of particles (ca be modified in the following do loop)
      nu = 7
      ! Dimension of the problem
      n = 2**nu
      ! Couplings of the Heisenberg model
      couplings = (/2d0,2d0,1d0/)
      ! Field strenght
      lambda = 0.d0

      l = 50
      from = 0
      to = 3.d0

      if(rank .eq. 0 ) then
            open(unit = 22, file="results/data.txt", action="write" )
      end if
      
       if(rank .eq. 0 ) then
            open(unit=23, file='./results/times'//trim(args(1))//'.txt', action="write")
      end if



      do ii = 2,12
            do jj = 0,l

                  nu = ii
                  n = 2**nu

                  ! Create matrix and set values
                  call MatCreate(PETSC_COMM_WORLD,A,ierr)
                  
                  call MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,n,n,ierr)
                  !call MatSetType(A,MATSBAIJ,ierr)
                  call MatSetUp(A,ierr)
            
                  lambda = from + (to-from)/l*jj 
                  !lambda = 0d0
                  ! Store values in the distributed matrix
                  
                   if (jj == l .and. rank == 0) then
		       start_build = MPI_Wtime()
		   end if
		    
                  
                  call heisenbergmodel_hamiltonian(A,nu,lambda,couplings,rank)
                  call MatCreateVecs(A,xr,xi,ierr)

		 if (jj == l .and. rank == 0) then
                	finish_build = MPI_Wtime()
            	end if
            
                  ! Diagonalize
                  kr = 0.d0
                  
                  if (jj == l .and. rank == 0) then
                	start = MPI_Wtime()
                  end if
                  
                  call EPSCreate(PETSC_COMM_WORLD,eps,ierr)
                  call eigz(A,eps,nu,xr,xi,rank,kr,error,time)
                  
                  if (jj == l .and. rank == 0) then
                	finish = MPI_Wtime()
            	  end if

                  if (rank .eq. 0) then
                        print*, "N :", nu," ---- lambda :",PetscRealPart(lambda)
                        write(22,*) nu,PetscRealPart(lambda),PetscRealPart(kr)/(nu-1),error,time
                  endif
                  
                 

                  ! Deallocate
                  call MatDestroy(A,ierr)
                  call VecDestroy(xr,ierr)
                  call VecDestroy(xi,ierr)
                  call EPSDestroy(eps,ierr)

            end do
            
              if (rank .eq. 0) then 
           		 write(23,*) nu, finish_build - start_build, finish - start 
            		 !N , Build time, diagonalize time
       		  end if
            !write(22,*)
      end do
      
      call SlepcFinalize(ierr)

end
