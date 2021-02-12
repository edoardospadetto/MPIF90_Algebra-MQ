module hamiltonians
!Module containing all subroutine to get a full hamiltonian 
!the resulting matrix must be a scalapack distributed one
use debug_module
use matrix_interface
use scalapack_interface



implicit none 

integer :: nb_for_hamiltonians = 4 
integer :: lda_for_hamiltonians = 200


contains 
!TRANSVERSE FIELD ISING MODEL
!------------------------------------------------------------

!computes field term of the hamiltonian
subroutine hamfield(N,lambda,context,pt,descpt)
       
       integer:: N ,ii
       complex*16, dimension(:, :), intent(INOUT) :: pt
       complex*16,dimension(2,2):: pauliz
       real*8 ::  lambda
       
       integer :: context,info
       integer , dimension(9) :: descA, descpt
       complex*16, dimension(lda_for_hamiltonians,lda_for_hamiltonians) :: A
       
 
       
       call breakifn("Invalid columns/Rows number", ((2**N .eq. descpt(4)) .and. (2**N .eq. descpt(3))) , .true.)
       
       pauliz = 0.0
       pauliz(1,1)%re=1.0
       pauliz(1,2)%re=0.0
       pauliz(2,1)%re=0.0
       pauliz(2,2)%re=-1.0
       pt = 0.0
       CALL DESCINIT( DESCA, 2**N, 2**N, nb_for_hamiltonians, nb_for_hamiltonians, 0, 0, context, lda_for_hamiltonians, info)
       
       do ii = 1,N
       	    call getbigmat(pauliz,ii,N,A,descA)
            call dsum(A,descA,pt,descpt,pt,descpt) 
       end do
       pt = lambda*pt
       
 end subroutine

subroutine haminteraction(A,descA,N)
integer , dimension(9) :: descA 
double complex, dimension(:,:) :: A
integer :: ii , jj , kk, ll ,res, N
A = dcmplx(0.d0,0.d0)
do ii = 0, 2**N-1
	do jj = 0, 2**N-1
		res = 0 
	 	do kk = 1, N-1
		    if ((2**(kk-1)+2**kk) .eq. xor(jj,ii)) then 
		        res = res + 1
		    end if 
                end do 
                call pzelset(A,ii+1,jj+1,descA, dcmplx(real(res),0.0) )
        end do 
end do 

end subroutine

subroutine hamintbroken(N, context , pt, descpt)
      implicit none
      integer:: N,ii
      complex*16, dimension(2**N, 2**N) :: pt
      complex*16,dimension(2,2):: paulix
    
      integer , dimension(9) :: descA, descB ,descC ,descpt
      complex*16, dimension(lda_for_hamiltonians,lda_for_hamiltonians):: A, B,C
      integer :: context, info 
 
      
      paulix = 0.0
      paulix(1,1)%re=0.0
      paulix(1,2)%re=1.0
      paulix(2,1)%re=1.0
      paulix(2,2)%re=0.0
      pt = 0.0
 
      CALL DESCINIT( DESCA, 2**N, 2**N, nb_for_hamiltonians, nb_for_hamiltonians, 0, 0, context, lda_for_hamiltonians, info)
      CALL DESCINIT( DESCB, 2**N, 2**N, nb_for_hamiltonians, nb_for_hamiltonians, 0, 0, context, lda_for_hamiltonians, info)  
      CALL DESCINIT( DESCC, 2**N, 2**N, nb_for_hamiltonians, nb_for_hamiltonians, 0, 0, context, lda_for_hamiltonians, info) 
      
      do ii = 1,N-1
      	    A = dcmplx(0.d0,0.d0)
      	    B = dcmplx(0.d0,0.d0)
      	    C = dcmplx(0.d0,0.d0)
      	    call getbigmat(paulix,ii+1,N,A,descA)
      	    call getbigmat(paulix,ii,N,B,descB)
      	    call dmatmul(A,descA,B,descB,C,descC)
            call dsum(C,descC,pt,descpt,pt,descpt) 
     end do
     	!call printmat(pt, descpt)
    
 end subroutine








subroutine transverse_field_ising_model_hamiltonian(context, lambda,N, hamiltonian, descHamiltonian)
	implicit none
	integer:: N, context, info
	real*8:: lambda

	integer , dimension(9) :: descA, deschamiltonian
	complex*16, dimension(lda_for_hamiltonians,lda_for_hamiltonians):: A, hamiltonian

       !TRANSVERSE FIELD ISING MODEL 
       CALL DESCINIT( DESCA, 2**N, 2**N, nb_for_hamiltonians, nb_for_hamiltonians, 0, 0, context, lda_for_hamiltonians, info)
       call haminteraction( A, descA, N)
       call dsum(A,descA,hamiltonian,descHamiltonian,hamiltonian,descHamiltonian) 
       !A = dcmplx(0.d0,0.d0)
       !call hamintbroken(N, context , A, descA)
       
       !print*, sum(hamiltonian-A)
       
       !call dsum(A,descA,hamiltonian,descHamiltonian,hamiltonian,descHamiltonian) 
       A = dcmplx(0.d0,0.d0)
       call hamfield(N,lambda, context , A, descA)
       call dsum(A,descA,hamiltonian,descHamiltonian,hamiltonian,descHamiltonian) 
      
       
end subroutine





!------------------------------------------------------------


end module
