!computes interaction term of the hamiltonian
!do not use
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
      	
      	    call dmatmul(A,descA,B,descB,C,descC)
            call dsum(C,descC,pt,descpt,pt,descpt) 
     end do
     	!call printmat(pt, descpt)
    
 end subroutine
 
 
 !compute kroeneker product of a matrix with and Identity matrices of the same dimensions
! prod{1_i-1} I_{sizeA} x A x prod{i+1_N} I_{sizeA}
!exploits properties of identitiy matrices.
!The result is stored in a scalapack distributed matrix
!input mat is not distributed
!It does not work !

subroutine getbigmat(inputmat, index, N ,bigmat,descBigMat)
integer :: index, N, ii,jj, kk1, kk2,idim
integer :: helpidx(2), blocksize
integer , dimension(:) :: descBigMat
double complex , dimension(:,:) :: inputmat
double complex, dimension( size(inputmat,dim = 1)**N, size(inputmat,dim = 1)**N) :: bigmat

call breakifn("Invalid rows number" ,(size(inputmat,dim = 1)**N .eq. size(bigmat,dim =1)), .true.)
call breakifn("Invalid columns number" ,( size(inputmat,dim = 1)**N .eq. size(bigmat,dim =2)), .true.)

idim = size(inputmat,dim = 1)
blocksize = (idim**(N-index))

bigmat = 0.0
do ii = 1, idim
    do jj = 1, idim
        do kk1 = 0,idim**(index-1)-1
            helpidx(1)= (ii)+idim*kk1
            helpidx(2)= (jj)+idim*kk1
            do kk2 = 1, blocksize
                !bigmat((helpidx(1)-1)*blocksize+kk2 ,(helpidx(2)-1)*blocksize+kk2)  &
                !   = inputmat(ii,jj)
                CALL PZELSET( bigmat, (helpidx(1)-1)*blocksize+kk2, (helpidx(2)-1)*blocksize+kk2, descBigMat, inputmat(ii,jj) )
            end do
        end do
    end do
end do
end subroutine
