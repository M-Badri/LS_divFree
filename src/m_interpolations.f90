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
        procedure :: get_scalar_momentum_matrices

    end type t_scalar_mls



    type, extends(t_interpolations) :: t_vector_mls

    contains

        procedure, nopass :: constructor => vector_constructor
        procedure :: get_interpolated => get_v_interpolated
        procedure :: get_vector_momentum_matrices

    end type t_vector_mls



contains



    function scalar_constructor (order, nlvls, nx, ny, dx_c, dy_c) result(intrp)

        implicit none
        class(t_scalar_mls), allocatable :: intrp
        integer(kind = I4), intent(in) :: order, nlvls, nx, ny
        real(kind = R8), intent(in) :: dx_c, dy_c

        integer (kind = I4) :: m, n, nd

        allocate (t_scalar_mls :: intrp)

        selectcase (order)
            case (2)
                intrp%nc = 6

            case (3)
                intrp%nc = 8

            case default
                stop "Wrong order for interpolation!"
        end select

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

        call intrp%get_scalar_momentum_matrices(nd)

    end function scalar_constructor



    subroutine get_scalar_momentum_matrices(this, nd)

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

    end subroutine get_scalar_momentum_matrices



    subroutine get_s_interpolated (this, lvl, n1, n2, xd, yd, u, u_f)

        implicit none
        class(t_scalar_mls), intent(in) :: this
        integer(kind = I4),  intent(in) :: lvl, n1, n2
        real(kind = R8), intent(out)    :: u_f(:,:)
        real(kind = R8), intent(in)     :: u(:)
        real(kind = R8), intent(in)     :: xd(:), yd(:)

        real(kind = R8), allocatable    :: ut(:), xt(:), yt(:)
        real(kind = R8), allocatable    :: xc(:), yc(:), c(:), uc(:)
        real(kind = R8), allocatable    :: temp(:), temp2(:)
        real(kind = R8)    :: dxc, dyc
        integer(kind = I4) :: n, i, j, m

        allocate(ut(9))
        allocate(xt(9))
        allocate(yt(9))
        allocate(xc(this%nchlds))
        allocate(yc(this%nchlds))
        allocate(uc(this%nchlds))

        selectcase (this%nc)

            case(6)
                temp  = [1.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0]
                temp2 = [0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0]

            case(8)
                temp  = [1.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0]
                temp2 = [0.0d0, 0.0d0, 0.0d0, 2.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0]

        end select

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

                    c = matmul (this%arr_of_a_inv(m, lvl)%a, &
                        matmul (this%arr_of_b_mat(m, lvl)%a, ut))
                    uc(m) = dot_product (c, temp)

!                    print*, m, uc(m) - dcos(xt(5)+xc(m)) * dsin(yt(5)+yc(m))
                    if(this%nc == 6) then
                        print*, m, dot_product (c, temp) &
                                 , + dcos(xt(5)+xc(m)) * dsin(yt(5)+yc(m))
                    end if
                end do

            end do
        end do

    end subroutine get_s_interpolated



    function vector_constructor (order, nlvls, nx, ny, dx_c, dy_c) result(intrp)

        implicit none
        class(t_vector_mls), allocatable :: intrp
        integer(kind = I4), intent(in) :: order, nlvls, nx, ny
        real(kind = R8), intent(in) :: dx_c, dy_c

        integer (kind = I4) :: m, n, nd

        allocate (t_vector_mls :: intrp)

        selectcase (order)
            case (2)
                intrp%nc = 12

            case (3)
                intrp%nc = 16

            case default
                stop "Wrong order for interpolation!"
        end select

        intrp%nchlds = 4
        intrp%nlvls = nlvls
        intrp%dx_c = dx_c
        intrp%dy_c = dy_c

        nd = 9
        allocate (intrp%xld(nd))
        allocate (intrp%yld(nd))
        allocate (intrp%uld(nd))
        allocate (intrp%vld(nd))

        allocate (intrp%arr_of_a_inv(intrp%nchlds, nlvls))
        allocate (intrp%arr_of_b_mat(intrp%nchlds, nlvls))
        do n = 1, nlvls
            do m = 1, intrp%nchlds
                allocate(intrp%arr_of_a_inv(m, n)%a(intrp%nc, intrp%nc))
                allocate(intrp%arr_of_b_mat(m, n)%a(intrp%nc, 3*nd))
            end do
        end do

        call intrp%get_vector_momentum_matrices(3*nd)

    end function vector_constructor



    subroutine get_vector_momentum_matrices(this, nd)

        implicit none
        class(t_vector_mls) :: this
        integer(kind = I4), intent(in) :: nd
        integer(kind = I4) :: m, n
        real(kind = R8), allocatable :: p(:,:), w(:,:), a(:,:), a_inv(:,:)
        real(kind = R8), allocatable :: xtmp(:), ytmp(:), b(:,:)
        real(kind = R8) :: dx, dy
        real(kind = R8) :: xi(4), yi(4)

        do n = 1, this%nlvls
            call this%fill_stencil_matrises (n)
            dx = this%dx_at_level (n) / 2.0d0
            dy = this%dy_at_level (n) / 2.0d0
            xi(1:3:2) = -dx
            xi(2:4:2) =  dx
            yi(1:2:1) = -dy
            yi(3:4:1) =  dy

            do m = 1, this%nchlds
                xtmp = this%xld - xi(m)
                ytmp = this%yld - yi(m)
                p = mls_p_matrix_vector (nd/3, this%nc, xtmp, ytmp)
                w = mls_w_matrix_vector (nd/3, xtmp, ytmp, 0.0d0, 0.0d0)
                b = matmul (transpose(p), w)
                a = matmul (b, p)
                this%arr_of_a_inv(m, n)%a(:,:) = inverse_matrix (a)
!                print*, size(b,1), size(b,2),nd
                this%arr_of_b_mat(m, n)%a(:,:) = b
            end do
        end do

    end subroutine get_vector_momentum_matrices



    subroutine get_v_interpolated (this, lvl, n1, n2, xd, yd, u, v, u_f)

        implicit none
        class(t_vector_mls), intent(in) :: this
        integer(kind = I4),  intent(in) :: lvl, n1, n2
        real(kind = R8), intent(out)    :: u_f(:,:)
        real(kind = R8), intent(in)     :: u(:), v(:)
        real(kind = R8), intent(in)     :: xd(:), yd(:)

        real(kind = R8), allocatable    :: ut(:), xt(:), yt(:)
        real(kind = R8), allocatable    :: xc(:), yc(:), c(:), uc(:)
        real(kind = R8), allocatable    :: temp(:)!, temp2(:)
        real(kind = R8)    :: dxc, dyc
        integer(kind = I4) :: n, i, j, m

        allocate(ut(27))
        allocate(xt(9))
        allocate(yt(9))
        allocate(xc(this%nchlds))
        allocate(yc(this%nchlds))
        allocate(uc(this%nchlds))

        selectcase (this%nc)

            case(12)
                temp  = [1.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, &
                         1.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0]
!                temp2 = [0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0]

            case(8)
                temp  = [1.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0]
!                temp2 = [0.0d0, 0.0d0, 0.0d0, 2.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0]

        end select

        n = 0
        do i = 2, n1 - 1
            do j = 2, n2 - 1
                n = (n1*(i-1)) + j
                ut(1:7:3)   = u(n-n1-1:n-n1+1)
                ut(2:8:3)   = v(n-n1-1:n-n1+1)
                ut(3:9:3)   = 0.0d0
                ut(10:16:3) = u(n-1:n+1)
                ut(11:17:3) = v(n-1:n+1)
                ut(12:18:3) = 0.0d0
                ut(19:25:3) = u(n+n1-1:n+n1+1)
                ut(20:26:3) = v(n+n1-1:n+n1+1)
                ut(21:27:3) = 0.0d0
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


!                    print*, size(this%arr_of_b_mat(m, lvl)%a, 1)
!                    print*, size(this%arr_of_b_mat(m, lvl)%a, 2)
                    c = matmul (this%arr_of_a_inv(m, lvl)%a, &
                        matmul (this%arr_of_b_mat(m, lvl)%a, ut))
                    uc(m) = dot_product (c(6:12), temp(6:12))

                    print*, m, uc(m) , dcos(xt(5)+xc(m)) * dsin(yt(5)+yc(m)), &
                                       -dsin(xt(5)+xc(m)) * dcos(yt(5)+yc(m))
!                    if(this%nc == 8) then
!                        print*, m, dot_product (c, temp2) &
!                                 + dcos(xt(5)+xc(m)) * dsin(yt(5)+yc(m))
!                    end if
!                    print*," "

                end do

            end do
        end do

    end subroutine get_v_interpolated



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



end module m_interpolations
