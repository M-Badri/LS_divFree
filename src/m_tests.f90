module m_tests
    use iso_fortran_env, only: R8 => real64, I4 => int32


    type t_test_data

        integer(kind = I4) :: n1, n2
        integer(kind = I4) :: n_lvls
        real(kind = R8) :: dx_c, dy_c
        real(kind = R8), allocatable, dimension(:) :: ll_corner_x, ll_corner_y
        real(kind = R8), allocatable, dimension(:,:) :: x1d, y1d
        real(kind = R8), allocatable, dimension(:,:) :: x, y
        real(kind = R8), allocatable, dimension(:,:) :: u, v

    end type t_test_data


contains



    function test_data_factory (test_data_num) result(td)

        implicit none
        class(t_test_data), allocatable :: td
        integer(kind = I4), intent(in) :: test_data_num


        selectcase (test_data_num)
            case (1)
                td = test_data_1 ()
        end select

    end function test_data_factory



    function test_data_1 () result (td1)

        implicit none
        class(t_test_data), allocatable :: td1
        real(kind = R8), parameter :: pi = dacos(-1.0d0)
        real(kind = R8) :: lx = 2.0d0*pi, ly = 2.0d0*pi
        integer(kind = I4) :: nx = 16, ny = 16
        integer(kind = I4) :: i, j, n, m
        real(kind = R8) :: dx, dy

        allocate(td1)
        td1%n_lvls = 1
        td1%n1 = 4
        td1%n2 = 4

        td1%dx_c = lx / real(nx, kind = R8)
        td1%dy_c = ly / real(ny, kind = R8)

        allocate (td1%ll_corner_x(td1%n_lvls))
        allocate (td1%ll_corner_y(td1%n_lvls))
        td1%ll_corner_x(:) = ((nx / 2) * td1%dx_c)
        td1%ll_corner_y(:) = ((ny / 2) * td1%dy_c)

        allocate (td1%x1d(td1%n_lvls, td1%n1))
        allocate (td1%y1d(td1%n_lvls, td1%n2))
        do n = 1, td1%n_lvls
            dx = td1%dx_c / 2**(n-1)
            dy = td1%dy_c / 2**(n-1)
            do j = 1, td1%n2
                td1%y1d(n, j) = (real(j, kind=R8) - 0.5d0) * dy
            end do
            do i = 1, td1%n1
                td1%x1d(n, i) = (real(i, kind=R8) - 0.5d0) * dx
            end do
        end do

        allocate (td1%u(td1%n_lvls, td1%n1 * td1%n2))
        allocate (td1%v(td1%n_lvls, td1%n1 * td1%n2))
        allocate (td1%x(td1%n_lvls, td1%n1 * td1%n2))
        allocate (td1%y(td1%n_lvls, td1%n1 * td1%n2))
        do m = 1, td1%n_lvls
            n = 0
            do j = 1, td1%n2
                do i = 1, td1%n1
                    n = n + 1
                    td1%x(m, n) =  td1%x1d(m, i)
                    td1%y(m, n) =  td1%y1d(m, j)
                    td1%u(m, n) =  dcos (td1%x(m, n)) * dsin (td1%y(m, n))
                    td1%v(m, n) = -dsin (td1%x(m, n)) * dcos (td1%y(m, n))
                end do
            end do
        end do

    end function test_data_1



end module m_tests
