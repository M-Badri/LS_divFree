program simplified_2

    use mod_tests
    use mod_least_square
    use mod_write

    implicit none
    real(kind = 8), allocatable, dimension (:)   :: xd, yd, ud, vd
    real(kind = 8), allocatable, dimension (:)   :: x1d, y1d
    real(kind = 8), allocatable, dimension (:,:) :: a

    call get_test_1 (x1d, y1d, xd, yd, ud, vd)
    a = get_scalar_ls_matrix (xd, yd, ud)

    call write_out_scalar_field (x1d, y1d, xd, yd, ud, "intial_field.dat")


end program simplified_2


