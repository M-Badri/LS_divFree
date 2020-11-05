module mod_write


contains



    subroutine write_out_scalar_field (x1d, y1d, x, y, u, fname)

        implicit none
        real (kind=8), intent(in) :: x(:), y(:), u(:), x1d(:), y1d(:)
        character(len=*), intent (in) :: fname
        integer :: i, j, n, nunit
        character(len=100) :: fname_c

        call execute_command_line ("mkdir -p outputs" )
        fname_c = trim (adjustl ("./outputs/" // fname))
        open(newunit = nunit, file = fname_c)
        write (nunit, *) "ZONE I= ",  size(x1d), " ,J=", size(y1d) , " ,f=point"

        n = 0
        do j = 1, size(y1d)
            do i = 1, size(x1d)
                n = n + 1
                write(nunit,*) x1d(i), x1d(j), u(n)
            end do
        end do
        close (nunit)

    end subroutine write_out_scalar_field



end module mod_write
