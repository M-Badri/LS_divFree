module m_interpolations

    use iso_fortran_env, only: R8 => real64, I4 => int32
    use m_least_square

    private
    public :: t_scalar_mls, t_vector_mls



    type a_of_p
        real(kind = R8), dimension(:,:), pointer :: a
    end type a_of_p



    type :: t_interpolations

        private
        integer(kind = I4) :: nlvls
        integer(kind = I4) :: nc, nchlds
        type(a_of_p),    allocatable :: arr_of_a_inv(:,:)
        type(a_of_p),    allocatable :: arr_of_b_mat(:,:)
        real(kind = R8), allocatable :: xld(:), yld(:)
        real(kind = R8), allocatable :: uld(:), vld(:)
        real(kind = R8)              :: dx_c, dy_c

    contains

        procedure, private :: fill_stencil_matrises
        procedure, private :: dx_at_level
        procedure, private :: dy_at_level

    end type t_interpolations



    type, extends(t_interpolations) :: t_scalar_mls

    contains

        procedure, nopass :: constructor => scalar_constructor
        procedure :: get_interpolated => get_s_interpolated
        procedure :: get_momentum_matrices

    end type t_scalar_mls



    type, extends(t_interpolations) :: t_vector_mls

    contains

        procedure, nopass :: constructor => vector_constructor
        procedure :: get_interpolated => get_v_interpolated

    end type t_vector_mls



contains



    function scalar_constructor (order, nlvls, nx, ny, dx_c, dy_c) result(intrp)

        implicit none
        class(t_scalar_mls), allocatable :: intrp
        integer(kind = I4), intent(in) :: order, nlvls, nx, ny
        real(kind = R8), intent(in) :: dx_c, dy_c

        integer (kind = I4) :: m, n, nd

        allocate (t_scalar_mls :: intrp)

        if (order == 2) then
            intrp%nc = 6
        else
            stop "Wrong order for interpolation!"
        end if

        intrp%nchlds = 4
        intrp%nlvls = nlvls
        intrp%dx_c = dx_c
        intrp%dy_c = dy_c

        nd = 9
        allocate (intrp%xld(nd))
        allocate (intrp%yld(nd))
        allocate (intrp%uld(nd))

        allocate (intrp%arr_of_a_inv(intrp%nchlds, nlvls))
        allocate (intrp%arr_of_b_mat(intrp%nchlds, nlvls))
        do n = 1, nlvls
            do m = 1, intrp%nchlds
                allocate(intrp%arr_of_a_inv(m, n)%a(intrp%nc, intrp%nc))
                allocate(intrp%arr_of_b_mat(m, n)%a(intrp%nc, nd))
            end do
        end do

        call intrp%get_momentum_matrices(nd)

    end function scalar_constructor



    subroutine get_momentum_matrices(this, nd)

        implicit none
        class(t_scalar_mls) :: this
        integer(kind = I4), intent(in) :: nd
        integer(kind = I4) :: m, n
        real(kind = R8), allocatable :: p(:,:), w(:,:), a(:,:), a_inv(:,:)
        real(kind = R8), allocatable :: xtmp(:), ytmp(:), b(:,:)
        real(kind = R8) :: dx, dy
        real(kind = R8) :: xi(4), yi(4)

        do n = 1, this%nlvls
            call this%fill_stencil_matrises (n)
            dx = this%dx_at_level (n) / 2.0d0
            dy = this%dx_at_level (n) / 2.0d0
!            xi(1:2:1) = -dx
!            xi(3:4:1) =  dx
!            yi(1:3:2) =  dy
!            yi(2:4:2) = -dy
            xi(1:3:2) = -dx
            xi(2:4:2) =  dx
            yi(1:2:1) = -dy
            yi(3:4:1) =  dy

            do m = 1, this%nchlds
                xtmp = this%xld - xi(m)
                ytmp = this%yld - yi(m)
                p = mls_p_matrix_scalar (nd, this%nc, xtmp, ytmp)
                w = mls_w_matrix_scalar (nd, xtmp, ytmp, 0.0d0, 0.0d0)
                b = matmul (transpose(p), w)
                a = matmul (b, p)
                this%arr_of_a_inv(m, n)%a(:,:) = inverse_matrix (a)
                this%arr_of_b_mat(m, n)%a(:,:) = b
            end do
        end do

    end subroutine get_momentum_matrices



    subroutine fill_stencil_matrises(this, nlvl)

        implicit none
        class(t_interpolations), intent(in out) :: this
        integer(kind = I4), intent(in) :: nlvl
        real(kind = R8) :: dx, dy

        dx = this%dx_at_level (nlvl)
        dy = this%dy_at_level (nlvl)

        this%xld(1:7:3) = -dx
        this%xld(2:8:3) =  0.0d0
        this%xld(3:9:3) =  dx

        this%yld(1:3:1) = -dy
        this%yld(4:6:1) =  0.0d0
        this%yld(7:9:1) =  dy

    end subroutine fill_stencil_matrises



    function dx_at_level (this, lvl) result(dx)

        implicit none
        class(t_interpolations), intent(in) :: this
        integer(kind = I4), intent(in) :: lvl
        real(kind = R8) :: dx

        dx = this%dx_c / (2.d0**(lvl-1))
    end function dx_at_level



    function dy_at_level (this, lvl) result(dy)

        implicit none
        class(t_interpolations), intent(in) :: this
        integer(kind = I4), intent(in) :: lvl
        real(kind = R8) :: dy

        dy = this%dy_c / (2.d0**(lvl-1))
    end function dy_at_level



    function vector_constructor (order, nlvls, nx, ny, dx_c, dy_c) result(intrp)

        implicit none
        class(t_interpolations), allocatable :: intrp
        integer(kind = I4), intent(in) :: order, nlvls, nx, ny
        real(kind = R8), intent(in) :: dx_c, dy_c

        integer (kind = I4) :: nc, n

        allocate (t_scalar_mls :: intrp)

        if (order == 2) then
            nc = 6
        else
            stop "Wrong order for interpolation!"
        end if

!        allocate (intrp%arr_of_a_inv(nlvls))
!        do n = 1, nlvls
!            allocate(intrp%arr_of_a_inv(4, n)%a(nx,ny))
!        end do

    end function vector_constructor



    subroutine get_s_interpolated (this, lvl, n1, n2, xd, yd, u, u_f)

        implicit none
        class(t_scalar_mls), intent(in) :: this
        integer(kind = I4),  intent(in) :: lvl, n1, n2
        real(kind = R8), intent(out)    :: u_f(:,:)
        real(kind = R8), intent(in)     :: u(:)
        real(kind = R8), intent(in)     :: xd(:), yd(:)!, llcor_x, llcor_y
        real(kind = R8), allocatable    :: xt(:), yt(:), ut(:)
        integer (kind = I4) :: n, i, j, m
        real(kind = R8) :: dxc, dyc
        real(kind = R8), allocatable :: xc(:), yc(:), temp(:), c(:), uc(:)

!        print*, this%dx_at_level(lvl), lvl, llcor_y
!        print*, xd, size(u_f, 2)
!        print*, " "
!        print*, yd(1:n1*n2:n2)
!        allocate(xt, source = xd(1:n1:1))
!        allocate(yt, source = yd(1:n1*n2:n2))
        allocate(ut(9))
        allocate(xt(9))
        allocate(yt(9))
        allocate(xc(this%nchlds))
        allocate(yc(this%nchlds))
        allocate(uc(this%nchlds))
!        xt(:) = xd(1:n1:1)
!        yt(:) = yd(1:n2*n1:n2)
!        print*, xt
!        print*, yt
!        print*,""
        n = 0
        do i = 2, n1 - 1
            do j = 2, n2 - 1
                n = (n1*(i-1)) + j
                ut(1:3) = u(n-n1-1:n-n1+1)
                ut(4:6) = u(n-1:n+1)
                ut(7:9) = u(n+n1-1:n+n1+1)
                xt(1:3) = xd(n-n1-1:n-n1+1)
                xt(4:6) = xd(n-1:n+1)
                xt(7:9) = xd(n+n1-1:n+n1+1)
                yt(1:3) = yd(n-n1-1:n-n1+1)
                yt(4:6) = yd(n-1:n+1)
                yt(7:9) = yd(n+n1-1:n+n1+1)
                dxc = 0.5d0 * this%dx_at_level(lvl)
                dyc = 0.5d0 * this%dy_at_level(lvl)
                xc(1:3:2) =  - dxc
                xc(2:4:2) =  + dxc
                yc(1:2:1) =  - dyc
                yc(3:4:1) =  + dyc

                do m = 1, this%nchlds
!                    xt(:) = xt(:) - xc(n)
!                    yt(:) = yt(:) - yc(n)
                    c = matmul (this%arr_of_a_inv(m, lvl)%a, &
                        matmul (this%arr_of_b_mat(m, lvl)%a, ut))
!                    temp = [1.d0, xc(m), yc(m), xc(m)**2, xc(m)*yc(m), yc(m)**2]
                    temp = [1.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0]
                    uc(m) = dot_product (c, temp)
                    print*, m, uc(m)- dcos(xt(5)+xc(m)) * dsin(yt(5)+yc(m))
                    print*," "

                end do

            end do
        end do


    end subroutine get_s_interpolated



    subroutine get_v_interpolated (this, u, v, vel_f)

        implicit none
        class(t_vector_mls) :: this
        real (8) :: u(:), v(:), vel_f(:)

    end subroutine get_v_interpolated



end module m_interpolations
