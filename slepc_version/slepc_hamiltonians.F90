
#include "slepc/finclude/slepcbv.h"
#include "slepc/finclude/slepcds.h"
#include "slepc/finclude/slepcrg.h"
#include "slepc/finclude/slepcfn.h"
#include "slepc/finclude/slepceps.h"
#include "slepc/finclude/slepcpep.h"

module hamiltonians
    use slepceps
    implicit none 
    
    contains 
    
    subroutine hamfield(A,n,lambda)
          
        Mat            :: A
        PetscErrorCode :: ierr
        PetscInt       :: n
        PetscScalar    :: lambda

        integer:: ii,jj, temp
        
        do ii = 0, 2**n-1
            temp = 0
            do jj = 0, n-1 
                temp = temp - 2*mod(ii/(2**jj), 2) +1
            end do 
            if(temp .ne.0) then
                call MatSetValue(A,ii,ii,dcmplx(dble(temp*lambda),0.0),ADD_VALUES,ierr)
            end if
        end do
            
    end subroutine

    
    subroutine haminteractionx(A,n,J)

        Mat            :: A
        PetscErrorCode :: ierr
        PetscInt       :: n
        PetscScalar    :: J
        integer :: ii , jj , kk, res

        do ii = 0, 2**n-1
            do jj = 0, 2**n-1
                res = 0 
                do kk = 1, n-1  
                    if ((2**(kk-1)+2**kk) .eq. xor(jj,ii)) then
                        res = res + 1
                    end if    
                end do 

               ! if( xor(jj,ii) .eq. (2**(n-1)+1) ) then
               !     res = res + 1
               ! end if

                if(res .ne.0) then
                    call MatSetValue(A,ii,jj,-J*dcmplx(real(res),0.0),ADD_VALUES,ierr) 
                end if
            end do 
        end do 
    
    end subroutine

    
    subroutine haminteractiony(A,n,J)

        Mat            :: A
        PetscErrorCode :: ierr
        PetscInt       :: n
        PetscScalar    :: J
        integer :: ii , jj , kk, res, testa, testb

        do ii = 0, 2**n-1
            do jj = 0, 2**n-1
                res = 0 
                do kk = 1, n-1    
                    if ((2**(kk-1)+2**kk) .eq. xor(jj,ii)) then         
                        testa = mod(ii/2**kk,2)
                        testb = mod(ii/2**(kk-1),2) 
                        res   = mod(not(XOR(TESTA,TESTB)),2)*2+1 +res     
                    end if 
                end do 

               ! if( xor(jj,ii) .eq. (2**(n-1)+1) ) then
               !     testa = mod(ii/2**0,2)
               !     testb = mod(ii/2**(n-1),2) 
               !     res   = mod(not(XOR(TESTA,TESTB)),2)*2+1 +res     
               ! end if

                if(res .ne.0) then
                    call MatSetValue(A,ii,jj,-J*dcmplx(real(res),0.0),ADD_VALUES,ierr) 
                end if
            end do 
        end do 
    
    end subroutine
    
    
    subroutine haminteractionz(A,n,J)

        Mat            :: A
        PetscErrorCode :: ierr
        PetscInt       :: n
        PetscScalar    :: J
        integer :: ii, kk, res, testa, testb

        do ii = 0, 2**N-1
            res = 0 
            do kk = 1, N-1     
                testa = mod(ii/2**kk,2)
                testb = mod(ii/2**(kk-1),2) 
                res =   -2*abs(testa-testb)+1 +res              
            end do 

            !testa = mod(ii/2**0,2)
            !testb = mod(ii/2**(n-1),2) 
            !res =   -2*abs(testa-testb)+1 +res 

            if(res .ne.0) then
                call MatSetValue(A,ii,ii,-J*dcmplx(dble(res),0.0),ADD_VALUES,ierr) 
            end if
        end do 
    
    end subroutine
       
    
    subroutine heisenbergmodel_hamiltonian(A,nu,lambda,couplings,rank)
        implicit none
        Mat            :: A
        PetscErrorCode :: ierr
        PetscInt       :: nu
        PetscScalar    :: couplings(3)
        PetscScalar    :: lambda  
        PetscMPIInt    :: rank 
           
        ! Initialize to zero
        call MatZeroEntries(A, ierr)

        
        if(rank .eq. 0) then

            ! X interaction
            call haminteractionx(A,nu,couplings(1))
    
            ! Y interaction
            call haminteractiony(A,nu,couplings(2))
            
            ! Z interaction
            call haminteractionz(A,nu,couplings(3))

            ! Field
            call hamfield(A,nu,lambda)

        end if

        ! Assembly and View
        call MatSetOption(A,MAT_IGNORE_ZERO_ENTRIES,PETSC_TRUE,ierr)
        !call MatSetOption(A,MAT_IGNORE_LOWER_TRIANGULAR,PETSC_TRUE,ierr)
        call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr)
        call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr)
        !call MatView(A,PETSC_VIEWER_STDOUT_WORLD,ierr)

    end subroutine
    
    
end module
