


MODULE matmulqic
    use debugqic
    !variable to enable debug
    LOGICAL :: debug_var = .FALSE.

  CONTAINS

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !enable debug for the function declared in this module
    SUBROUTINE set_debug(mybool)
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        LOGICAL, INTENT(IN) :: mybool
        debug_var = mybool
    END SUBROUTINE set_debug






    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !perform matrix multiplication
    SUBROUTINE matmul1(aa,bb,result)
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      IMPLICIT NONE
      !inout vars
      REAL, DIMENSION(:,:), INTENT(IN) :: aa,bb
      REAL, DIMENSION(:,:), ALLOCATABLE, INTENT(OUT) :: result
      !temp vars
      INTEGER*2 :: ii,jj,kk
      INTEGER, DIMENSION(SIZE(SHAPE(aa))) :: sizea
      INTEGER, DIMENSION(SIZE(SHAPE(bb))) :: sizeb

      sizea=SHAPE(aa)
      sizeb=SHAPE(bb)
      !allocate memory to the result of multiplication
      ALLOCATE (result(sizea(1),sizeb(2)))

      !check equality between rows and columns using the debug subroutine
      call breakifn("Invalid Matrix Shapes", (sizea(2) .EQ. sizeb(1)), debug_var)

      !compute the multiplication using the loop
        DO ii=1,sizea(1)
          DO jj=1,sizeb(2)
            result(ii,jj)=0
            DO kk=1,sizea(2)
              result(ii,jj) = result(ii,jj)+ aa(ii,kk)*bb(kk,jj)
            END DO
          END DO
        END DO
    END SUBROUTINE matmul1



    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !perform matrix multiplication with different order
    SUBROUTINE matmul2(aa,bb,result)
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      !everything is equal to matmul 1 except the order of the nested for loops
      IMPLICIT NONE
      INTEGER*2 :: ii,jj,kk
      REAL, DIMENSION(:,:), INTENT(IN) :: aa,bb
      REAL, DIMENSION(:,:), ALLOCATABLE, INTENT(OUT) :: result

      INTEGER, DIMENSION(SIZE(SHAPE(aa))) :: sizea
      INTEGER, DIMENSION(SIZE(SHAPE(bb))) :: sizeb

      sizea=SHAPE(aa)
      sizeb=SHAPE(bb)

      ALLOCATE (result(sizea(1),sizeb(2)))

      call breakifn("Invalid Matrix Shapes", (sizea(2) .EQ. sizeb(1)), debug_var)

        DO jj=1,sizeb(2)
          DO ii=1,sizea(1)
            result(ii,jj)=0
            DO kk=1,sizea(2)
              result(ii,jj) = result(ii,jj)+ aa(ii,kk)*bb(kk,jj)
            END DO
          END DO
        END DO


    END SUBROUTINE matmul2

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !perform matrix multiplication with different order
    SUBROUTINE matmul3(aa,bb,result)
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

     ! I introduced another variable, aa1 that's the transpose of aa. I use this in the loop to perform 2 loops over columns and one over rows.

      IMPLICIT NONE
      INTEGER*2 :: ii,jj,kk
      REAL, DIMENSION(:,:), INTENT(IN) :: aa,bb
      REAL, DIMENSION(:,:), ALLOCATABLE, INTENT(OUT) :: result
      REAL, DIMENSION(:,:), ALLOCATABLE :: aa1

      INTEGER, DIMENSION(SIZE(SHAPE(aa))) :: sizea
      INTEGER, DIMENSION(SIZE(SHAPE(bb))) :: sizeb

      sizea=SHAPE(aa)
      sizeb=SHAPE(bb)

      ALLOCATE (result(sizea(1),sizeb(2)))
      call breakifn("Invalid Matrix Shapes", (sizea(2) .EQ. sizeb(1)), debug_var)

        ALLOCATE(aa1(sizea(2), sizea(1)))
        aa1=TRANSPOSE(aa)
        DO jj=1,sizeb(2)
          DO ii=1,sizea(1)
            result(ii,jj)=0
            DO kk=1,sizea(2)
              result(ii,jj) = result(ii,jj)+ aa1(kk,ii)*bb(kk,jj)
            END DO
          END DO
        END DO


    END SUBROUTINE matmul3

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !perform matrix multiplication with different order
    SUBROUTINE matmul4(aa,bb,result)
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !as matmul3 but 2 loops on rows and one on columns transposing bb

      IMPLICIT NONE
      INTEGER*2 :: ii,jj,kk
      REAL, DIMENSION(:,:), INTENT(IN) :: aa,bb
      REAL, DIMENSION(:,:), ALLOCATABLE, INTENT(OUT) :: result
      REAL, DIMENSION(:,:), ALLOCATABLE :: bb1

      INTEGER, DIMENSION(SIZE(SHAPE(aa))) :: sizea
      INTEGER, DIMENSION(SIZE(SHAPE(bb))) :: sizeb

      sizea=SHAPE(aa)
      sizeb=SHAPE(bb)

      ALLOCATE (result(sizea(1),sizeb(2)))

      call breakifn("Invalid Matrix Shapes", (sizea(2) .EQ. sizeb(1)), debug_var)


        ALLOCATE(bb1(sizeb(2), sizeb(1)))
        bb1=TRANSPOSE(bb)
        DO ii=1,sizea(1)
          DO jj=1,sizeb(2)
            result(ii,jj)=0
            DO kk=1,sizea(2)
              result(ii,jj) = result(ii,jj)+ aa(ii,kk)*bb1(jj,kk)
            END DO
          END DO
        END DO

    END SUBROUTINE matmul4

END MODULE matmulqic
