module printer
    use debugqic
contains

   !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !add line to file outfile as INT   REAL
    SUBROUTINE wfio(outfile, myint, myreal)
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        real , intent(IN) :: myreal
        integer , intent(IN):: myint
        character(len = *), intent(IN):: outfile
        logical::  file_exist = .FALSE.
        character(len = 3) :: stat = "new"

        inquire(FILE=outfile, EXIST=file_exist)

        if (file_exist) then
            stat = "old"
        end if

        open(unit = 2, file=outfile, status = stat, access = 'append')

        write(2, "(I5,E22.3)") myint,myreal

        close(2)


    END SUBROUTINE

    !###############################################
    !compute and output eigenvectors real part and probabilyty density function
    subroutine ozgrid(x1,x2,vec,time, filename)
    !###############################################
        use matrixqic
        implicit none
        double complex, dimension(:) :: vec
        character*3 :: filename
        character*20 :: numb
        integer :: which, L,ii
        double precision :: x1, x2,eps, helpr,helpim
        integer::time

        L = size(vec, 1)
        eps=  (x2-x1)/real(L)

        write(numb,"('_',I8.8)") time
        numb = trim(adjustl(numb))

        vec = vec/sqrt(norm2(vec%re)**2+norm2(vec%im)**2)
        !output data
        do ii = 1,L
            helpr = (vec(ii)%re)/sqrt(eps)
            helpim = (vec(ii)%im)/sqrt(eps)
            call wdddo(trim("./results/")//trim(filename)//trim(numb)//".txt" , x1 + (x2-x1)/real(L)*ii, helpr, helpim)
        end do
    end subroutine


    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !add line to file outfile as real   REAL
    SUBROUTINE wffo(outfile, myr, myreal)
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        real , intent(IN) :: myreal,myr

        character(len = *), intent(IN):: outfile
        logical::  file_exist = .FALSE.
        character(len = 3) :: stat = "old"

        inquire(FILE=outfile, EXIST=file_exist)

        if (.not. file_exist) then
            stat = "new"
        else if (file_exist) then
            stat = "old"
        end if



        open(unit = 2, file=outfile, status = stat, access = 'append')

        write(2, "(E22.3,E22.3)") myr,myreal

        close(2)
    END SUBROUTINE

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !add line to file outfile as real   REAL
    SUBROUTINE wddo(outfile, myr, myreal)
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        double precision , intent(IN) :: myreal,myr

        character(len = *), intent(IN):: outfile
        logical::  file_exist = .FALSE.
        character(len = 3) :: stat = "old"

        inquire(FILE=outfile, EXIST=file_exist)

        if (.not. file_exist) then
            stat = "new"
        else if (file_exist) then
            stat = "old"
        end if



        open(unit = 2, file=outfile, status = stat, access = 'append')

        write(2, "(E22.3,E22.3)") myr,myreal

        close(2)
    END SUBROUTINE

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !add line to file outfile as real   REAL
    SUBROUTINE wdddo(outfile, myr, myreal,myrealll)
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        double precision , intent(IN) :: myreal,myr,myrealll

        character(len = *), intent(IN):: outfile
        logical::  file_exist = .FALSE.
        character(len = 3) :: stat = "old"

        inquire(FILE=outfile, EXIST=file_exist)

        if (.not. file_exist) then
            stat = "new"
        else if (file_exist) then
            stat = "old"
        end if



        open(unit = 2, file=outfile, status = stat, access = 'append')

        write(2, "(E22.3,E22.3,E22.3)") myr,myreal,myrealll

        close(2)
    END SUBROUTINE

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !add line to file outfile as real   REAL
    SUBROUTINE widdo(outfile,myi, myr, myreal)
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        real*8 , intent(IN) :: myreal,myr
        integer, intent(IN):: myi
        character(len = *), intent(IN):: outfile
        logical::  file_exist = .FALSE.
        character(len = 3) :: stat = "old"

        inquire(FILE=outfile, EXIST=file_exist)

        if (.not. file_exist) then
            stat = "new"
        else if (file_exist) then
            stat = "old"
        end if


        open(unit = 2, file=outfile, status = stat, access = 'append')

        write(2, "(I5,E22.3,E22.3)") myi,myr,myreal

        close(2)


    END SUBROUTINE


    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !add line to file outfile as REAL
    SUBROUTINE wdo(outfile, myreal)
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        double precision , intent(IN) :: myreal
        character(len = *), intent(IN):: outfile
        logical::  file_exist = .FALSE.
        character(len = 3) :: stat = "new"

        inquire(FILE=outfile, EXIST=file_exist)

        if (file_exist) then
            stat = "old"
        else
            stat="new"
        end if

        open(unit = 2, file=outfile, status = stat, access = 'append')

        write(2, "(E22.3)") myreal

        close(2)


    END SUBROUTINE
end module
