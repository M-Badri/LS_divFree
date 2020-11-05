program simplified_2

    use mod_tests
    use mod_math
    use mod_least_square
    use mod_write

    implicit none
    real(kind = 8), allocatable, dimension (:)   :: xd, yd, ud, vd, btu, coef
    real(kind = 8), allocatable, dimension (:)   :: x1d, y1d
    real(kind = 8), allocatable, dimension (:,:) :: b, btb, btb_i

    call get_test_1 (x1d, y1d, xd, yd, ud, vd)
    call write_out_scalar_field (x1d, y1d, xd, yd, ud, "intial_field.dat")

    b = get_scalar_ls_matrix (xd, yd, ud)
    btb = matmul (transpose (b), b)
    btu = matmul (transpose (b), ud)
    btb_i = inverse_matrix (size(btb, 1), size(btb, 2), btb)
    coef =  matmul (btb_i, btu)

end program simplified_2


