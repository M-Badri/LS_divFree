program simplified_2

    use mod_tests
    use mod_math
    use mod_least_square
    use mod_write
    USE LA_PRECISION, ONLY: WP => DP
    USE F95_LAPACK, ONLY: LA_GESVD, LA_GESV


    implicit none
    real(kind = 8), allocatable, dimension (:)   :: xd, yd, ud, vd, btu, coef
    real(kind = 8), allocatable, dimension (:)   :: x1d, y1d
    real(kind = 8), allocatable, dimension (:,:) :: b, btb, btb_inv
    real(kind = 8), allocatable, dimension (:)   :: xi, yi, ui_e, vi_e
    real(kind = 8), allocatable, dimension (:)   :: ui_scalar
    real(kind = 8), allocatable :: ui(:)
    integer ::  nd

!    call get_test_1_data (x1d, y1d, xd, yd, ud, vd)
!    call write_out_scalar_field (x1d, y1d, xd, yd, ud, "intial_field.dat")

    call get_test_2_data (x1d, y1d, xd, yd, ud, vd)

    b = ls_matrix_scalar (xd, yd)
    btb = matmul (transpose (b), b)
    btu = matmul (transpose (b), ud)

    btb_inv = inverse_matrix (btb)

!    b_i = inverse_matrix (b)
!    btb_inv_bt = matmul (btb_inv, transpose(b))

     coef =  matmul (btb_inv, btu)

    call get_test_2_interplation_coordinates (x1d, y1d, xi, yi, ui_e, vi_e)
    call get_interpolation_scalar (xi, yi, coef, ui_scalar)

    nd = size(ud)
    allocate(ui(1))
    call get_mls_scalar (nd, xd, yd, ud, [xi(1)], [yi(1)], ui)
    print*, ui - ui_e(1)

end program simplified_2


