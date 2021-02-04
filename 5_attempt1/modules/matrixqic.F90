

module matrixqic
    use debugqic

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
  !generate random complex double matrix
  FUNCTION rgcm(sizea,sizeb) RESULT(gen)
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: sizea,sizeb
    INTEGER ::  ii, jj
    DOUBLE precision :: a, b
    DOUBLE COMPLEX, DIMENSION(sizea,sizeb) :: gen

    do ii= 1, sizea
        do jj = 1, sizeb
            call RANDOM_NUMBER(a)
            call RANDOM_NUMBER(b)
            gen(ii,jj) = cmplx(a,b,kind(0D0))
        end do
    end do
END FUNCTION rgcm

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!generate random double diagonal matrix
FUNCTION rgddm(sizem) RESULT(gen)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: sizem
  INTEGER ::  ii, jj
  DOUBLE precision :: a
  DOUBLE precision, DIMENSION(sizem,sizem) :: gen

  do ii= 1, sizem
      do jj = 1, sizem
          gen(ii,jj) =0.0
      end do
      call RANDOM_NUMBER(a)
      gen(ii,ii) = a
  end do
END FUNCTION rgddm


    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !generates random real matrix
    FUNCTION rgrm(sizerow,sizecol) RESULT(gen)
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: sizerow,  sizecol
      REAL, DIMENSION(sizerow,sizecol) :: gen
      call RANDOM_NUMBER(gen)
  END FUNCTION rgrm

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !print a matrix on terminal
    SUBROUTINE prm(mat)
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      IMPLICIT NONE
      !inout vars
      REAL, DIMENSION(:,:), INTENT(IN) :: mat
      !tmp vars
      INTEGER*2 :: ii,jj
      INTEGER, DIMENSION(SIZE(SHAPE(mat))) :: sizem

      sizem=SHAPE(mat)

      PRINT*, ""
      DO ii=1,sizem(1)
        DO jj=1,sizem(2)
          WRITE(*, fmt="(f0.2, tr2)", advance="no") mat(ii,jj)
        END DO
        PRINT*, ""
      END DO

  END SUBROUTINE prm

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !read dimension from file.
    FUNCTION read_dim(path) RESULT(dim)
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    character(len = *) :: path
    integer :: dim
    logical :: file_exist = .FALSE.
    INQUIRE(FILE=path, EXIST=file_exist)

    call breakifn(path // "does not exists",file_exist, .TRUE.)

    open(unit = 3, file = path , status = "old" )

    read(3, fmt= "(i8)") dim

    close(3)

    END FUNCTION

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!print a complexmatrix*8 on terminal
SUBROUTINE pzm(mat)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   IMPLICIT NONE
   INTEGER :: ii,jj
   DOUBLE COMPLEX, DIMENSION(:,:), INTENT(IN) :: mat
   INTEGER, DIMENSION(SIZE(SHAPE(mat))) :: sizem
   sizem=SHAPE(mat)

   PRINT*, ""
   DO ii=1,sizem(1)
        write(*, "(*(' 'sf6.2xspf6.2x'i ':x))")  mat(ii,:)
     PRINT*, ""
   END DO


END SUBROUTINE pzm

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!print a complexmatrix on terminal
SUBROUTINE pcm(mat)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   IMPLICIT NONE
   INTEGER :: ii,jj
   COMPLEX, DIMENSION(:,:), INTENT(IN) :: mat
   INTEGER, DIMENSION(SIZE(SHAPE(mat))) :: sizem
   sizem=SHAPE(mat)

   PRINT*, ""
   DO ii=1,sizem(1)
        write(*, "(*(' 'sf6.2xspf6.2x'i ':x))")  mat(ii,:)
     PRINT*, ""
   END DO


END SUBROUTINE pcm

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!print a real matrix on terminal
SUBROUTINE pdm(mat)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   IMPLICIT NONE
   INTEGER :: ii,jj
   DOUBLE PRECISION, DIMENSION(:,:), INTENT(IN) :: mat
   INTEGER, DIMENSION(SIZE(SHAPE(mat))) :: sizem
   sizem=SHAPE(mat)

   PRINT*, ""
   DO ii=1,sizem(1)
        write(*, *)  mat(ii,:)
     PRINT*, ""
   END DO


END SUBROUTINE pdm




    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !compute eigenvalues and eigenvector with zheev, 
    subroutine eigz(A,dimA,ev)
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

            PRINT*, INT(WORK(1))
            LWORK = MIN(LWMAX, INT(WORK(1)))
            call zheev ( 'V','U',dimA,A,dimA,ev,WORK,LWORK,RWORK,INFO)
            call breakifn("bhoo", info .eq. 0, .TRUE.)



            print*, "eig stuff computed"



    end subroutine eigz





end module
