module mod_tests


contains



    subroutine get_test_1_data (x1d, y1d, xd, yd, ud, vd)

        implicit none
        real(kind = 8), allocatable, dimension(:), intent (in out) :: xd, yd, ud, vd
        real(kind = 8), allocatable, dimension(:), intent (in out) :: x1d, y1d
        integer, parameter :: nx = 8, ny = 8
        real (kind = 8), parameter  :: pi = acos(-1.d0)
        real (kind = 8), parameter  :: lx = 2.d0 * pi, ly = 2.d0 * pi
        integer :: i, j, n
        integer :: nd

        nd = nx * ny
        allocate (xd(nd), yd(nd), ud(nd), vd(nd))
        allocate (x1d(nx), y1d(ny))

        do i = 1, nx
            x1d(i) = (i - 0.5d0) * (lx / real(nx, kind = 8) )
        end do
        do j = 1, ny
            y1d(j) = (j - 0.5d0) * (ly / real(ny, kind = 8) )
        end do

        n = 0
        do j = 1, ny
            do i = 1, nx
                n = n + 1
                xd(n) = (i - 0.5d0) * (lx / real(nx, kind = 8) )
                yd(n) = (j - 0.5d0) * (ly / real(ny, kind = 8) )
                ud(n) = xd(n)**2 + yd(n)**2
                vd(n) =  dsin (xd(n)) * dcos (yd(n))
            end do
        end do

    end subroutine get_test_1_data



    subroutine get_test_1_interplation_coordinates (xd, yd, xi, yi, ui_e, vi_e)

        implicit none
        real (kind = 8), intent(in)  :: xd(:), yd(:)
        real (kind = 8), allocatable, intent(in out) :: xi(:), yi(:), ui_e(:), vi_e(:)
        real (kind = 8), allocatable :: xt(:), yt(:)
        integer :: ni, i, j, m

        allocate (xt(size(xd)-1))
        allocate (yt(size(yd)-1))
        xt(:) = 0.5d0 * &
                (xd(lbound(xd,1):ubound(xd,1)-1) + xd(lbound(xd,1)+1:ubound(xd,1)))
        yt(:) = 0.5d0 * &
                (yd(lbound(yd,1):ubound(yd,1)-1) + yd(lbound(yd,1)+1:ubound(yd,1)))

        ni = size(xt) * size(yt)
        allocate (xi(ni))
        allocate (yi(ni))
        allocate (vi_e(ni))
        allocate (ui_e(ni))

        m = 0
        do j = 1, size(yt)
            do i = 1, size(xt)
                m = m + 1
                xi(m) = xt(i)
                yi(m) = yt(j)
                ui_e(m) = xi(m)**2 + yi(m)**2 !dcos (xi(i)) * dsin (yi(j))
            end do
        end do

    end subroutine get_test_1_interplation_coordinates



    subroutine get_test_2_data (x1d, y1d, xd, yd, ud, vd)

        implicit none
        real(kind = 8), allocatable, dimension(:), intent (in out) :: xd, yd, ud, vd
        real(kind = 8), allocatable, dimension(:), intent (in out) :: x1d, y1d
        integer, parameter :: nx = 128,  ny = 128
        integer :: nd_xlo = 57, nd_xhi = 60
        integer :: nd_ylo = 14, nd_yhi = 17
        real (kind = 8), parameter  :: pi = acos(-1.d0)
        real (kind = 8), parameter  :: lx = 2.d0 * pi, ly = 2.d0 * pi
        integer :: i, j, n, nd

        nd = (nd_xhi - nd_xlo + 1) * (nd_yhi - nd_ylo + 1)
        if(allocated(xd))  deallocate(xd)
        if(allocated(yd))  deallocate(yd)
        if(allocated(ud))  deallocate(ud)
        if(allocated(vd))  deallocate(vd)
        if(allocated(x1d)) deallocate(x1d)
        if(allocated(y1d)) deallocate(y1d)
        allocate (xd(nd), yd(nd), ud(nd), vd(nd))
        allocate (x1d(nd_xhi - nd_xlo + 1), y1d(nd_yhi - nd_ylo + 1))

        n = 0
        do i = nd_xlo, nd_xhi
            n = n + 1
            x1d(n) = (i - 0.5d0) * (lx / real(nx, kind = 8) )
        end do
        n = 0
        do j = nd_ylo, nd_yhi
            n = n + 1
            y1d(n) = (j - 0.5d0) * (ly / real(ny, kind = 8) )
        end do

        n = 0
        do j = nd_ylo , nd_yhi
            do i = nd_xlo , nd_xhi
                n = n + 1
                xd(n) = (i - 0.5d0) * (lx / real(nx, kind = 8) )

                yd(n) = (j - 0.5d0) * (ly / real(ny, kind = 8) )
                ud(n) = -dcos (xd(n)) * dsin (yd(n))
                vd(n) =  dsin (xd(n)) * dcos (yd(n))
            end do
        end do

    end subroutine get_test_2_data



    subroutine get_test_2_interplation_coordinates (xd, yd, xi, yi, ui_e, vi_e)

        implicit none
        real (kind = 8), intent(in)  :: xd(:), yd(:)
        real (kind = 8), allocatable, intent(in out) :: xi(:), yi(:), ui_e(:), vi_e(:)
        real (kind = 8), allocatable :: xt(:), yt(:)
        integer :: ni, i, j, m

        allocate (xt(size(xd)-1))
        allocate (yt(size(yd)-1))
        xt(:) = 0.5d0 * &
                (xd(lbound(xd,1):ubound(xd,1)-1) + xd(lbound(xd,1)+1:ubound(xd,1)))
        yt(:) = 0.5d0 * &
                (yd(lbound(yd,1):ubound(yd,1)-1) + yd(lbound(yd,1)+1:ubound(yd,1)))

        ni = size(xt) * size(yt)
        if(allocated(xi))   deallocate (xi)
        if(allocated(yi))   deallocate (yi)
        if(allocated(vi_e)) deallocate (vi_e)
        if(allocated(ui_e)) deallocate (ui_e)
        allocate (xi(ni))
        allocate (yi(ni))
        allocate (vi_e(ni))
        allocate (ui_e(ni))

        m = 0
        do j = 1, size(yt)
            do i = 1, size(xt)
                m = m + 1
                xi(m) = xt(i)
                yi(m) = yt(j)
                ui_e(m) = -dcos (xi(m)) * dsin (yi(m))
                vi_e(m) =  dsin (xi(m)) * dcos (yi(m))
            end do
        end do

    end subroutine get_test_2_interplation_coordinates



end module mod_tests
