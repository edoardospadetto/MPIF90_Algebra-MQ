module hamiltonians
!Module containing all subroutine to get a full hamiltonian 
!the resulting matrix must be a scalapack distributed one
use debug_module
use matrix_interface
use scalapack_interface



implicit none 

contains 
!TRANSVERSE FIELD ISING MODEL
!------------------------------------------------------------

!computes field term of the hamiltonian
subroutine hamfield(N,lambda,context,pt,descpt)
       
       integer:: N ,ii
       complex*16, dimension(:, :), intent(INOUT) :: pt
       complex*16,dimension(2,2):: pauliz
       real*8 ::  lambda
       
       integer :: lda ,context,info
       PARAMETER (LDA = 100)
       integer , dimension(9) :: descA, descpt
       complex*16, dimension(lda,lda) :: A
       
       
       !nb
       integer :: nb = 4
       
       call breakifn("Invalid columns/Rows number", ((2**N .eq. descpt(4)) .and. (2**N .eq. descpt(3))) , .true.)
       
       pauliz = 0.0
       pauliz(1,1)%re=1.0
       pauliz(1,2)%re=0.0
       pauliz(2,1)%re=0.0
       pauliz(2,2)%re=-1.0
       pt = 0.0
       CALL DESCINIT( DESCA, 2**N, 2**N, nb, nb, 0, 0, context, lda, info)
       
       do ii = 1,N
       	    call getbigmat(pauliz,ii,N,A,descA)
            call dsum(A,descA,pt,descpt,pt,descpt) 
       end do
       pt = lambda*pt
 end subroutine

!computes interaction term of the hamiltonian
subroutine hamint(N, context , pt, descpt)
      implicit none
      integer:: N,ii
      complex*16, dimension(2**N, 2**N) :: pt
      complex*16,dimension(2,2):: paulix
      integer :: lda 
      PARAMETER (LDA = 100)
      integer , dimension(9) :: descA, descB ,descC ,descpt
      complex*16, dimension(lda,lda):: A, B,C
      integer :: context, info 
      integer :: nb 

      nb = 4
      
      paulix = 0.0
      paulix(1,1)%re=0.0
      paulix(1,2)%re=1.0
      paulix(2,1)%re=1.0
      paulix(2,2)%re=0.0
      pt = 0.0
 
      CALL DESCINIT( DESCA, 2**N, 2**N, nb, nb, 0, 0, context, lda, info)
      CALL DESCINIT( DESCB, 2**N, 2**N, nb, nb, 0, 0, context, lda, info)  
      CALL DESCINIT( DESCC, 2**N, 2**N, nb, nb, 0, 0, context, lda, info) 
      
      do ii = 1,N-1
      	    call getbigmat(paulix,ii+1,N,A,descA)
      	    call getbigmat(paulix,ii,N, B,descB)
      	    call dmatmul(A,descA,B,descB,C,descC)
            call dsum(C,descC,pt,descpt,pt,descpt) 
     end do
    
 end subroutine




!compute full hamiltonian, summing up interaction and
!field term
subroutine transverse_field_ising_model_hamiltonian(context, lambda,N, hamiltonian, descHamiltonian)
	implicit none
	integer:: N, context, info
	real*8:: lambda
	integer :: lda 
	PARAMETER (LDA = 100)
	integer , dimension(9) :: descA, deschamiltonian
	complex*16, dimension(lda,lda):: A, hamiltonian
        integer :: nb 
        nb = 4
       !TRANSVERSE FIELD ISING MODEL 
       CALL DESCINIT( DESCA, 2**N, 2**N, nb, nb, 0, 0, context, lda, info)
  
       call hamint(N, context , A, descA)
       !print*, "done interaction"
       call dsum(A,descA,hamiltonian,descHamiltonian,hamiltonian,descHamiltonian) 
       !print*, "done update"
       !call hamfield(N,lambda, context , A, descA)
       !print*, "done field"
       !call dsum(A,descA,hamiltonian,descHamiltonian,hamiltonian,descHamiltonian) 
       !print*, "done update"
       
end subroutine



!------------------------------------------------------------


end module
