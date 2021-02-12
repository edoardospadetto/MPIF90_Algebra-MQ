!###################################################
!generic module to handle matrices
!###################################################

module matrix_interface
	use debug_module
contains 


	!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	!generate random complex double matrix hermitian
	FUNCTION rghcm(tsize) RESULT(gen)
	!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	IMPLICIT NONE
	INTEGER, INTENT(IN) :: tsize
	INTEGER ::  ii, jj
	DOUBLE precision :: a, b
	DOUBLE COMPLEX, DIMENSION(tsize,tsize) :: gen

	do ii= 1, tsize
	  do jj = 1, ii
	      call RANDOM_NUMBER(a)
	      call RANDOM_NUMBER(b)
	      gen(ii,jj) = cmplx(a,b,kind(0D0))
	      if (ii .NE. jj) then
		  gen(jj,ii) = cmplx(a,-b,kind(0D0))
	      end if
	  end do
	end do
	END FUNCTION rghcm
	
	!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	!print a complexmatrix*8 on terminal
	SUBROUTINE pzm(mat)
	!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	   IMPLICIT NONE
	   INTEGER :: ii,jj
	   DOUBLE COMPLEX, DIMENSION(:,:), INTENT(IN) :: mat
	   INTEGER, DIMENSION(SIZE(SHAPE(mat))) :: sizem
	   sizem=SHAPE(mat)
	   print*, "ok"
	   PRINT*, ""
	   DO ii=1,sizem(1)
		write(*, "(*(' 'sf6.2xspf6.2x'i ':x))")  mat(ii,:)
	     PRINT*, ""
	   END DO


	END SUBROUTINE pzm
	
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !compute eigenvalues and eigenvector with zheev,
    subroutine eigz(A,dimA,ev)
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    	   use debug_module
           external zheev
           double complex, dimension(:,:), intent(INOUT) :: A
           double precision, dimension(:), intent(INOUT):: ev
           integer ::dimA,ii
           integer :: LWMAX, LWORK,INFO
           parameter(LWMAX = 1000000)
           complex*16 :: WORK(LWMAX)
           double precision :: RWORK(3*dimA-2)

           call breakifn("size eigenvalue not equal dimA  .. eigc:matrixqic ", size(ev, dim=1) .eq. dimA)
           call breakifn("size eigenvalue not equal dimA  .. eigc:matrixqic ", size(A, dim=1) .eq. dimA)
           call breakifn("size eigenvalue not equal dimA  .. eigc:matrixqic ", size(A, dim=2) .eq. dimA)

           LWORK = -1


            call zheev ( 'V','U',dimA,A,dimA,ev,WORK,LWORK,RWORK,INFO)

            !PRINT*, INT(WORK(1))
            LWORK = MIN(LWMAX, INT(WORK(1)))
            call zheev ( 'V','U',dimA,A,dimA,ev,WORK,LWORK,RWORK,INFO)
            call breakifn("bhoo", info .eq. 0, .TRUE.)



            !print*, "eig stuff computed"



    end subroutine eigz



end module 
