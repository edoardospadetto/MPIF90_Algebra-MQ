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


subroutine hamintsigmax(N, context , pt, descpt)
      implicit none
      integer:: N,ii
      complex*16, dimension(2**N, 2**N) :: pt
      complex*16,dimension(2,2):: paulix
    
      integer , dimension(9) :: descA, descB ,descC ,descpt
      complex*16, dimension(lda_for_hamiltonians,lda_for_hamiltonians):: A, B,C
      integer :: context, info 
 
      
      paulix = dcmplx(0.d0,0d0)
     
      paulix(1,2)=dcmplx(1.d0,0.d0)
      paulix(2,1)=dcmplx(1.d0,1.d0)
      
      pt = dcmplx(0.d0,0d0)
 
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
     !A = dcmplx(0.d0,0.d0)
     !B = dcmplx(0.d0,0.d0)
     !C = dcmplx(0.d0,0.d0)
     !call getbigmat(paulix,N,N,A,descA)
     !call getbigmat(paulix,1,N,B,descB)
     !call dmatmul(A,descA,B,descB,C,descC)
     !call dsum(C,descC,pt,descpt,pt,descpt) 
     	!call printmat(pt, descpt)
    
 end subroutine


subroutine hamintsigmaz(N, context , pt, descpt)
      implicit none
      integer:: N,ii
      complex*16, dimension(2**N, 2**N) :: pt
      complex*16,dimension(2,2):: pauliz
    
      integer , dimension(9) :: descA, descB ,descC ,descpt
      complex*16, dimension(lda_for_hamiltonians,lda_for_hamiltonians):: A, B,C
      integer :: context, info 
 
        pauliz = dcmplx(0.d0,0.d0)
	pauliz(1,1)%re=1.0
	pauliz(1,2)%re=0.0
	pauliz(2,1)%re=0.0
	pauliz(2,2)%re=-1.0
        pt = 0.0
 
      CALL DESCINIT( DESCA, 2**N, 2**N, nb_for_hamiltonians, nb_for_hamiltonians, 0, 0, context, lda_for_hamiltonians, info)
      CALL DESCINIT( DESCB, 2**N, 2**N, nb_for_hamiltonians, nb_for_hamiltonians, 0, 0, context, lda_for_hamiltonians, info)  
      CALL DESCINIT( DESCC, 2**N, 2**N, nb_for_hamiltonians, nb_for_hamiltonians, 0, 0, context, lda_for_hamiltonians, info) 
      
      do ii = 1,N-1
      	    A = dcmplx(0.d0,0.d0)
      	    B = dcmplx(0.d0,0.d0)
      	    C = dcmplx(0.d0,0.d0)
      	    call getbigmat(pauliz,ii+1,N,A,descA)
      	    call getbigmat(pauliz,ii,N,B,descB)
      	    call dmatmul(A,descA,B,descB,C,descC)
            call dsum(C,descC,pt,descpt,pt,descpt) 
     end do
     !A = dcmplx(0.d0,0.d0)
     !B = dcmplx(0.d0,0.d0)
     !C = dcmplx(0.d0,0.d0)
     !call getbigmat(pauliz,N,N,A,descA)
     !call getbigmat(pauliz,1,N,B,descB)
     !call dmatmul(A,descA,B,descB,C,descC)
     !call dsum(C,descC,pt,descpt,pt,descpt) 
     	!call printmat(pt, descpt)
    
 end subroutine


subroutine hamintsigmay(N, context , pt, descpt)
      implicit none
      integer:: N,ii
      complex*16, dimension(2**N, 2**N) :: pt
      complex*16,dimension(2,2):: pauliy
    
      integer , dimension(9) :: descA, descB ,descC ,descpt
      complex*16, dimension(lda_for_hamiltonians,lda_for_hamiltonians):: A, B,C
      integer :: context, info 
 
      
      pauliy = dcmplx(0.d0,0.d0)
      pauliy(1,2)%im=-1.0
      pauliy(2,1)%im=1.0
      pt = 0.0
 
      CALL DESCINIT( DESCA, 2**N, 2**N, nb_for_hamiltonians, nb_for_hamiltonians, 0, 0, context, lda_for_hamiltonians, info)
      CALL DESCINIT( DESCB, 2**N, 2**N, nb_for_hamiltonians, nb_for_hamiltonians, 0, 0, context, lda_for_hamiltonians, info)  
      CALL DESCINIT( DESCC, 2**N, 2**N, nb_for_hamiltonians, nb_for_hamiltonians, 0, 0, context, lda_for_hamiltonians, info) 
      
      do ii = 1,N-1
      	    A = dcmplx(0.d0,0.d0)
      	    B = dcmplx(0.d0,0.d0)
      	    C = dcmplx(0.d0,0.d0)
      	    call getbigmat(pauliy,ii+1,N,A,descA)
      	    call getbigmat(pauliy,ii,N,B,descB)
      	    call dmatmul(A,descA,B,descB,C,descC)
            call dsum(C,descC,pt,descpt,pt,descpt) 
     end do
     !A = dcmplx(0.d0,0.d0)
     !B = dcmplx(0.d0,0.d0)
     !C = dcmplx(0.d0,0.d0)
     !call getbigmat(pauliy,N,N,A,descA)
     !call getbigmat(pauliy,1,N,B,descB)
     !call dmatmul(A,descA,B,descB,C,descC)
     !call dsum(C,descC,pt,descpt,pt,descpt) 
     	!call printmat(pt, descpt)
    
 end subroutine


subroutine haminteractionz(A,descA,N)
integer , dimension(9) :: descA 
double complex, dimension(:,:) :: A
integer :: ii , jj , kk, ll ,res, N, testa, testb
A = dcmplx(0.d0,0.d0)
do ii = 0, 2**N-1
	do jj = 0, 2**N-1
		res = 0 
	 	do kk = 1, N-1  
	 		
		    if (0 .eq. xor(jj,ii)) then ! if (1 .eq. xor(jj,ii)) then 
		        if (ii .ne.jj ) then 
		        STOP 
		        end if 
		    	testa = mod(ii/2**kk,2)
		        testb= mod(ii/2**(kk-1),2) 
		        res =   -2*abs(testa-testb)+1 +res
		        
		    end if 
                end do 
                call pzelset(A,ii+1,jj+1,descA, dcmplx(dble(res),0.0) )
        end do 
end do 

end subroutine


!------------------------------------------------------
 !compute kroeneker product of a matrix with and Identity matrices of the same dimensions
! prod{1_i-1} I_{sizeA} x A x prod{i+1_N} I_{sizeA}
!exploits properties of identitiy matrices.
!The result is stored in a scalapack distributed matrix
!input mat is not distributed
!MUST be tested

subroutine getbigmat(inputmat, index, N ,bigmat,descBigMat)
integer :: index, N, ii,jj, kk1, kk2,idim
integer :: helpidx(2), blocksize, dimleft, dimright
integer , dimension(:) :: descBigMat
double complex , dimension(:,:) :: inputmat
double complex, dimension( size(inputmat,dim = 1)**N, size(inputmat,dim = 1)**N) :: bigmat

call breakifn("Invalid rows number" ,(size(inputmat,dim = 1)**N .eq. size(bigmat,dim =1)), .true.)
call breakifn("Invalid columns number" ,( size(inputmat,dim = 1)**N .eq. size(bigmat,dim =2)), .true.)

idim = size(inputmat,dim = 1)
blocksize = (idim**(N-index))


dimright= idim**(index-1)

bigmat = dcmplx(0.d0,0.d0)
do ii = 1, idim
    do jj = 1, idim
        do kk1 = 0, dimright-1
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

end module scalapack_interface


