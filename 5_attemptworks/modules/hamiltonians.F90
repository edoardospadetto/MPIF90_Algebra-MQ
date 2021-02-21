module hamiltonians
    
    !Module containing all subroutine to get a full hamiltonian 
    !the resulting matrix must be a scalapack distributed one

    use debug_module
    use matrix_interface
    use scalapack_interface

    implicit none 

    integer :: nb_for_hamiltonians = 4 
    integer :: lda_for_hamiltonians = 1000

contains 

    ! Field term 
    subroutine hamfield(N, lambda, context, pt, descpt)
        
        integer:: N, ii, jj, temp
        complex*16, dimension(:,:), intent(INOUT) :: pt
        real*8 :: lambda
        integer :: context, info
        integer, dimension(9) :: descpt 
        
        call breakifn("Invalid columns/Rows number!", ((2**N .eq. descpt(4)) .and. (2**N .eq. descpt(3))), .true.)
        
        pt = 0.0

        do ii = 0, 2**N-1
            temp = 0
            do jj = 0, N-1 
                temp = temp - 2*mod(ii/(2**jj), 2) + 1
            end do 
            call pzelset(pt, ii+1, ii+1, descpt, dcmplx(dble(temp*lambda),0.0))
        end do
     
    end subroutine

    ! Interaction term along x
    subroutine haminteractionx(A, descA, N)

        integer, dimension(9) :: descA 
        double complex, dimension(:,:) :: A
        integer :: ii, jj, kk, ll, res, N

        A = dcmplx(0.d0,0.d0)
        do ii = 0, 2**N-1
            do jj = 0, 2**N-1
                res = 0 
                do kk = 1, N-1  
                    if ((2**(kk-1)+2**kk) .eq. xor(jj,ii)) then 
                        res = res + 1
                    end if 
                end do 
                call pzelset(A, ii+1, jj+1, descA, dcmplx(real(res),0.0))
            end do 
        end do 

    end subroutine

    ! Interaction term along y
    subroutine haminteractiony(A, descA, N)

        integer, dimension(9) :: descA 
        double complex, dimension(:,:) :: A
        integer :: ii, jj, kk, ll, res, N, testa, testb

        A = dcmplx(0.d0,0.d0)
        do ii = 0, 2**N-1
            do jj = 0, 2**N-1
                res = 0 
                do kk = 1, N-1  
                    if ((2**(kk-1)+2**kk) .eq. xor(jj,ii)) then 
                        testa = mod(ii/2**kk, 2)
                        testb = mod(ii/2**(kk-1), 2) 
                        res = MOD(not(XOR(TESTA,TESTB)), 2)*2 + 1 + res
                    end if 
                end do 
                call pzelset(A, ii+1, jj+1, descA, dcmplx(real(res),0.0))
            end do 
        end do 

    end subroutine

    ! Interaction term along z
    subroutine haminteractionz(A, descA, N)

        integer, dimension(9) :: descA 
        double complex, dimension(:,:) :: A
        integer :: ii, jj, kk, ll, res, N, testa, testb

        A = dcmplx(0.d0,0.d0)
        do ii = 0, 2**N-1
            res = 0 
            do kk = 1, N-1  
                testa = mod(ii/2**kk,2)
                testb = mod(ii/2**(kk-1), 2) 
                res = -2*abs(testa-testb) + 1 + res
            end do 
            call pzelset(A, ii+1, ii+1, descA, dcmplx(dble(res),0.0))
        end do 

    end subroutine

    ! Transverse field Ising model hamiltonian
    subroutine transverse_field_ising_model_hamiltonian(context, lambda, N, hamiltonian, descHamiltonian)

        implicit none

        integer :: N, context, info
        real*8 :: lambda
        integer, dimension(9) :: descA, deschamiltonian
        complex*16, dimension(lda_for_hamiltonians,lda_for_hamiltonians):: A, hamiltonian

        hamiltonian = dcmplx(0.d0,0.d0)
        
        CALL DESCINIT(DESCA, 2**N, 2**N, nb_for_hamiltonians, nb_for_hamiltonians, 0, 0, context, lda_for_hamiltonians, info)

        A = dcmplx(0.d0,0.d0)

        call haminteractionx(A, descA, N)
        call dsum(A, descA, hamiltonian, descHamiltonian, hamiltonian, descHamiltonian) 

        A = dcmplx(0.d0,0.d0)

        call hamfield(N, lambda, context , A, descA)
        call dsum(A, descA, hamiltonian, descHamiltonian, hamiltonian, descHamiltonian) 

        ! call printmat(hamiltonian, descHamiltonian)
    
    end subroutine

    ! Heisenberg model hamiltonian
    subroutine heisenbergmodel_hamiltonian(context, lambda, couplings, N, hamiltonian, descHamiltonian,lda)
        
        implicit none
        
        integer :: N, context, info,lda
        real*8 :: lambda, couplings(3)
        integer, dimension(9) :: descA, deschamiltonian
        complex*16, dimension(lda,lda) :: A, hamiltonian

        CALL DESCINIT(descA, 2**N, 2**N, nb_for_hamiltonians, nb_for_hamiltonians, 0, 0, context, lda, info)

        A = dcmplx(0.d0,0.d0)

        call haminteractionx(A, descA, N)
        call dsum(-1*couplings(1)*A, descA, hamiltonian, descHamiltonian, hamiltonian, descHamiltonian)

        A = dcmplx(0.d0,0.d0)
        call haminteractiony(A, descA, N)
        call dsum(-1*couplings(2)*A, descA, hamiltonian, descHamiltonian, hamiltonian, descHamiltonian) 

        A = dcmplx(0.d0,0.d0)
        call haminteractionz(A, descA, N)
        call dsum(-1*couplings(3)*A, descA, hamiltonian, descHamiltonian, hamiltonian, descHamiltonian) 

        A = dcmplx(0.d0,0.d0)
        call hamfield(N, lambda, context, A, descA)
        call dsum(A, descA, hamiltonian, descHamiltonian, hamiltonian, descHamiltonian) 

        !call printmat(hamiltonian, descHamiltonian)

    end subroutine

end module
