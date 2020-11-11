module mod_least_square

    use mod_math

contains



    function ls_matrix_scalar (xd, yd) result (b)

        implicit none
        integer, parameter :: nc = 6
        real(kind = 8), intent (in) :: xd(:)    , yd(:)
        real(kind = 8), allocatable ::  b(:,:)
        integer :: i, nd

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
        real(kind = 8) :: p(nd,nc), w(nd,nd), b(nc,nd), a(nc,nc), a_inv(nc,nc)
        integer :: ni, n
        real(8) :: ui(1), c(nc), tt(6)

        ni = size(xi)
        do n = 1, ni
            p = mls_p_matrix_scalar (nd, nc, xd, yd)
            w = mls_w_matrix_scalar (nd, xd, yd, xi(n), yi(n))
            a = matmul (transpose(p), matmul(w, p))
            a_inv = inverse_matrix (a)
            b = matmul (transpose(p), w)
!            phi = matmul (transpose(p), matmul(a_inv, b))
            c = matmul (a_inv, matmul(b, ud))
!            ui = matmul (transpose(p), c)
        end do

        ui = 0.d0
        tt = [1.d0, xi, yi, xi**2, xi*yi, yi**2]
        do n = 1, nc
            ui = ui + c(n)* tt(n)
        end do

    end subroutine get_mls_scalar



    function mls_p_matrix_scalar (nd, nc, xd, yd) result (b)

        implicit none
        integer, intent(in) :: nd, nc
        real(kind = 8), intent(in) :: xd(nd), yd(nd)
        real(kind = 8) :: b(nd,nc)
        integer :: i

        do i = 1, nd
            b(i, :) = [1.d0, xd(i), yd(i), xd(i)**2, xd(i)*yd(i), yd(i)**2]
        end do

    end function mls_p_matrix_scalar



    function mls_w_matrix_scalar (nd, xd, yd, xi, yi) result (w)

        implicit none
        integer, intent(in) :: nd
        real(kind = 8), intent(in)  :: xd(nd), yd(nd), xi, yi
        real(kind = 8) :: w(nd,nd), r(nd)
        real(kind = 8) :: maxdis
        integer :: n

        maxdis = 0.0d0
        do n = 1, nd
            r(n) = norm2 ([xi-xd(n), yi-yd(n)])
            maxdis = max (r(n), maxdis)
        end do

        maxdis = maxval (r(:))
        r(:) = r(:) / (maxdis * 1.05d0)

        w(:,:) = 0.0d0
        do n = 1, nd
            w(n,n) = 1.d0 - (6.d0*r(n)**2) + (8.d0*r(n)**3) - (3.d0*r(n)**4)
        end do

    end function mls_w_matrix_scalar



end module mod_least_square
