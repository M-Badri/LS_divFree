module mod_least_square


contains



    function get_scalar_ls_matrix (xd, yd, ud) result (b)

        implicit none
        integer, parameter :: nc = 6
        real(kind = 8), intent (in) :: xd(:), yd(:), ud(:)
        real(kind = 8), allocatable :: r(:), c(:), b(:,:)
        integer :: i, j, nd, n

        nd = size (ud)
        allocate (b(nd, nc))

        do n = 1, nd
            b(n, :) = [1.d0, xd(n), yd(n), xd(n)**2,  xd(n)*yd(n), yd(n)**2]
        end do

!        a = matmul (transpose (at), a)

    end function get_scalar_ls_matrix



end module mod_least_square
