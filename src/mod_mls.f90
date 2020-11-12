module mod_least_square

    use mod_math

contains



    function ls_matrix_scalar (xd, yd) result (b)

        implicit none
        integer, parameter :: nc = 6
        real(kind = 8), intent (in) :: xd(:), yd(:)

        integer :: i, nd
        real(kind = 8), allocatable ::  b(:,:)

        nd = size (xd)
        allocate (b(nd, nc))

        do i = 1, nd
            b(i, :) = [1.d0, xd(i), yd(i), xd(i)**2,  xd(i)*yd(i), yd(i)**2]
        end do

    end function ls_matrix_scalar



    subroutine get_interpolation_scalar (xi, yi, coef, ui)

        implicit none
        real(kind = 8), intent(in) :: xi(:), yi(:), coef(:)
        real(kind = 8), allocatable, intent(in out) :: ui(:)
        real(kind = 8), allocatable :: a(:)
        integer :: n

        allocate (ui(size(xi)))
        allocate (a(size(coef)))

        do n = 1, size(ui)
            a(:) = [1.d0, xi(n), yi(n), xi(n)**2, xi(n)*yi(n), yi(n)**2]
            ui(n) = 0.0d0
            ui(n) = dot_product(coef, a)
        end do

    end subroutine get_interpolation_scalar



    subroutine get_mls_scalar (nd, xd, yd, ud, xi, yi, ui)

        implicit none
        integer, parameter :: nc = 6
        integer, intent(in) :: nd
        real(kind = 8), intent(in) :: xd(nd), yd(nd), ud(nd), xi(:), yi(:)
        real(kind = 8), intent(in out), allocatable  :: ui(:)

        integer :: ni, n
        real(kind = 8) :: p(nd,nc), w(nd,nd), b(nc,nd), a(nc,nc), a_inv(nc,nc)
        real(kind = 8) :: c(nc), temp(nc)

        ni = size(xi)
        allocate(ui(ni))
        do n = 1, ni
            p = mls_p_matrix_scalar (nd, nc, xd, yd)
            w = mls_w_matrix_scalar (nd, xd, yd, xi(n), yi(n))
            b = matmul (transpose(p), w)
            a = matmul (b, p)
            a_inv = inverse_matrix (a)
            c = matmul (a_inv, matmul(b, ud))
            temp = [1.d0, xi(n), yi(n), xi(n)**2, xi(n)*yi(n), yi(n)**2]
            ui(n) = dot_product (c, temp)
        end do

    end subroutine get_mls_scalar



    function mls_p_matrix_scalar (nd, nc, xd, yd) result (b)

        implicit none
        integer, intent(in) :: nd, nc
        real(kind = 8), intent(in) :: xd(nd), yd(nd)

        integer :: i
        real(kind = 8) :: b(nd,nc)

        do i = 1, nd
            b(i, :) = [1.d0, xd(i), yd(i), xd(i)**2, xd(i)*yd(i), yd(i)**2]
        end do

    end function mls_p_matrix_scalar



    function mls_w_matrix_scalar (nd, xd, yd, xi, yi) result (w)

        implicit none
        integer, intent(in) :: nd
        real(kind = 8), intent(in)  :: xd(nd), yd(nd), xi, yi

        integer :: n
        real(kind = 8) :: w(nd,nd), r(nd)
        real(kind = 8) :: maxdis

        do n = 1, nd
            r(n) = norm2 ([xi-xd(n), yi-yd(n)])
        end do

        maxdis = maxval (r(:))
        r(:) = r(:) / (maxdis * 1.05d0)

        w(:,:) = 0.0d0
        do n = 1, nd
            w(n,n) = 1.d0 - (6.d0*r(n)**2) + (8.d0*r(n)**3) - (3.d0*r(n)**4)
        end do

    end function mls_w_matrix_scalar



    subroutine get_mls_vector (nd, xd, yd, ud, vd, xi, yi, ui, vi)

        implicit none
        integer, parameter :: nc = 12
        integer, intent(in) :: nd
        real(kind = 8), intent(in) :: xd(nd), yd(nd), ud(nd), vd(nd), xi(:), yi(:)
        real(kind = 8), intent(in out), allocatable  :: ui(:), vi(:)

        integer :: ni, n
        real(kind = 8) :: p(3*nd,nc), w(3*nd,3*nd), b(nc,3*nd), a(nc,nc), a_inv(nc,nc)
        real(kind = 8) :: c(nc), temp(nc)
        real(kind = 8), allocatable :: vec_d(:)

        allocate (vec_d(3*nd))
        do n = 1, nd
            vec_d(3*n-2) = ud(n)
            vec_d(3*n-1) = vd(n)
            vec_d(3*n)   = 0.0d0
        end do

        ni = size(xi)
        allocate(ui(ni))
        allocate(vi(ni))

        do n = 1, ni
            p = mls_p_matrix_vector (nd, nc, xd, yd)
            w = mls_w_matrix_vector (nd, xd, yd, xi(n), yi(n))
            b = matmul (transpose(p), w)
            a = matmul (b, p)
            a_inv = inverse_matrix (a)
            c = matmul (a_inv, matmul(b, vec_d))

            !> Calculates U and V
            temp = [1.0d0, xi(n), yi(n), xi(n)**2, yi(n)**2, xi(n)*yi(n), &
                    1.0d0, xi(n), yi(n), xi(n)**2, yi(n)**2, xi(n)*yi(n)  ]
            !> U
            ui(n) = dot_product (c(1:nc/2), temp(1:nc/2))

            !> V
            vi(n) = dot_product (c(nc/2+1:nc), temp(nc/2+1:nc))
        end do

    end subroutine get_mls_vector



    function mls_p_matrix_vector (nd, nc, xd, yd) result (p)

        implicit none
        integer, intent(in) :: nd, nc
        real(kind = 8), intent(in) :: xd(nd), yd(nd)

        integer :: i
        real(kind = 8) :: p(3*nd,nc)

        do i = 1, nd
            p(3*i-2, :) = [1.0d0, xd(i), yd(i), xd(i)**2, yd(i)**2, xd(i)*yd(i), &
                           0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0]
            p(3*i-1, :) = [0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, &
                           1.0d0, xd(i), yd(i), xd(i)**2, yd(i)**2, xd(i)*yd(i)]
            p(3*i,   :) = [0.0d0, 1.0d0, 0.0d0, 2.0d0*xd(i), 0.0d0, yd(i), &
                           0.0d0, 0.0d0, 1.0d0, 0.0d0, 2.0d0*yd(i), xd(i) ]
        end do

    end function mls_p_matrix_vector



    function mls_w_matrix_vector (nd, xd, yd, xi, yi) result (w)

        implicit none
        integer, intent(in) :: nd
        real(kind = 8), intent(in)  :: xd(nd), yd(nd), xi, yi

        integer :: n
        real(kind = 8) :: w(3*nd,3*nd), r(nd)
        real(kind = 8) :: maxdis, wtemp

        do n = 1, nd
            r(n) = norm2 ([xi-xd(n), yi-yd(n)])
        end do

        maxdis = maxval (r(:))
        r(:) = r(:) / (maxdis * 1.05d0)

        w(:,:) = 0.0d0
        do n = 1, nd
            wtemp = 1.d0 - (6.d0*r(n)**2) + (8.d0*r(n)**3) - (3.d0*r(n)**4)
            w(3*n-2, 3*n-2) = wtemp
            w(3*n-1, 3*n-1) = wtemp
            w(3*n,   3*n)   = wtemp
        end do

    end function mls_w_matrix_vector



end module mod_least_square
