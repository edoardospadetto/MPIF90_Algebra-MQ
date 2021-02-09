
    module debugqic
    contains
        subroutine breakifn(d_string, cndtion, debug)
            implicit none
            character(len = *), intent(IN) :: d_string

            logical, intent(IN), optional :: debug
            logical, intent(IN) :: cndtion
            logical :: de
            ! if debug is enabled continue otherwise noop
            ! set debug to its right value
            de=.false.
            if( present(debug)) then
                de=debug
            end if

            if (de) then
                !if error occurs stop the program and print
                if (.not. cndtion) then
                    write (*,*) d_string
                    write(*,*) "Execution terminated"
                    STOP

                end if
            end if
        end subroutine breakifn

    end module debugqic
