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
        real(kind = 8), intent(in out) :: ui(:)
        real(kind = 8), intent(in) :: xd(nd), yd(nd), ud(nd), xi(:), yi(:)

        integer :: ni, n
        real(kind = 8) :: p(nd,nc), w(nd,nd), b(nc,nd), a(nc,nc), a_inv(nc,nc)
        real(kind = 8) :: c(nc), temp(nc)

        ni = size(xi)
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



end module mod_least_square
