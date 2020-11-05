module mod_least_square


contains

    function get_scalar_ls_matrix (xd, yd, ud) result (a)

        implicit none
        real(kind = 8), intent (in) :: xd(:), yd(:), ud(:)
        real(kind = 8), allocatable :: a (:,:)
        integer :: i, j, nd

        nd = size (ud)
        allocate (a(nd, nd))
    end function get_scalar_ls_matrix

end module mod_least_square
