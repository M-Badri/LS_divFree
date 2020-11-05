module mod_tests


contains



    subroutine get_test_1 (xd, yd, ud, vd)

        implicit none
        real(kind = 8), allocatable, dimension(:), intent (in out) :: xd, yd, ud, vd
        integer, parameter :: nx = 8, ny = 8
        real (kind = 8), parameter  :: pi = acos(-1.d0)
        real (kind = 8), parameter  :: lx = 2.d0 * pi, ly = 2.d0 * pi
        real (kind = 8) :: x, y
        integer :: i, j, n
        integer :: nd

        nd = nx * ny
        allocate (xd(nd), yd(nd), ud(nd), vd(nd))

        n = 0
        do j = 1, ny
            do i = 1, nx
                n = n + 1
                xd(n) = (i - 0.5d0) * (lx / real(nx-1, kind = 8) )
                yd(n) = (j - 0.5d0) * (ly / real(ny-1, kind = 8) )
                ud(n) = -dcos (xd(n)) * dsin (yd(n))
                vd(n) =  dsin (xd(n)) * dcos (yd(n))
            end do
        end do

    end subroutine get_test_1



end module mod_tests
