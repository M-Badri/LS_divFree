program simplified_2

    use mod_tests
    use m_math
    use m_least_square
    use mod_write
    use m_tests
    use m_interpolations

    implicit none
    real(kind = 8), allocatable, dimension (:)   :: xd, yd, ud, vd, btu, coef
    real(kind = 8), allocatable, dimension (:)   :: x1d, y1d
    real(kind = 8), allocatable, dimension (:,:) :: b, btb, btb_inv
    real(kind = 8), allocatable, dimension (:)   :: xi, yi, ui_e, vi_e
    real(kind = 8), allocatable, dimension (:)   :: ui_scalar
    real(kind = 8), allocatable :: ui(:), vi(:)
    integer ::  nd, i

!    call get_test_1_data (x1d, y1d, xd, yd, ud, vd)
!    call write_out_scalar_field (x1d, y1d, xd, yd, ud, "intial_field.dat")

!    call get_test_2_data (x1d, y1d, xd, yd, ud, vd)
!
!    b = ls_matrix_scalar (xd, yd)
!    btb = matmul (transpose (b), b)
!    btu = matmul (transpose (b), ud)
!
!    btb_inv = inverse_matrix (btb)
!
!!    b_i = inverse_matrix (b)
!!    btb_inv_bt = matmul (btb_inv, transpose(b))
!
!     coef =  matmul (btb_inv, btu)
!
!
!    !> Tests the scalar MLS for u
!    call get_test_2_interplation_coordinates (x1d, y1d, xi, yi, ui_e, vi_e)
!    call get_mls_scalar (size(ud), xd, yd, ud, xi, yi, ui)
!    print*, ui(:) - ui_e(:)
!    print*, ""
!
!    !> Tests the scalar MLS for v
!    call get_test_2_interplation_coordinates (x1d, y1d, xi, yi, ui_e, vi_e)
!    call get_mls_scalar (size(ud), xd, yd, vd, xi, yi, vi)
!    print*,  vi(:) - vi_e(:)
!
!    print*, "****************************************************************************"
!    print*, "****************************************************************************"
!    print*, " "
!
!    !> Tests the vector MLS for u and v
!    if (allocated(vi)) deallocate(vi)
!    if (allocated(ui)) deallocate(ui)
!    call get_test_2_interplation_coordinates (x1d, y1d, xi, yi, ui_e, vi_e)
!    call get_mls_vector (size(vd), xd, yd, ud, vd, xi, yi, ui, vi)
!    print*,  ui(:) - ui_e(:)
!    print*, " "
!    print*,  vi(:) - vi_e(:)
!    print*, "****************************************************************************"
!    print*, "****************************************************************************"

    block
        real(8), allocatable :: cdomain(:,:), u_fine(:,:)
        integer, allocatable :: ngrids(:), ngrids_box(:)
        real(8), parameter :: pi = acos(-1.0d0)
        type(t_test_data) :: td
        class(t_scalar_mls), allocatable :: s_intrp
        class(t_vector_mls), allocatable :: v_intrp
!        type(tp_one_lvl_test) :: tst
!        type(t_test_data) :: tt
        cdomain = reshape( [ 0.0d0, 0.0d0, 2.0d0*pi, 2.0d0*pi], [ 2, 2 ])
        ngrids = [128, 128]
        ngrids_box =[16, 16]
!        call tst%make_test (ngrids, cdomain, ngrids_box)
!        call tst%get_data_for_test (1, box_data)
!        print*, box_data%x
!        print*, box_data%y

        td = test_data_factory (test_data_num = 1)

!        s_intrp = interpolation_factory ("scalar", 2, td%n_lvls, &
!                                        td%n1, td%n2, td%dx_c, td%dy_c)

!        s_intrp = s_intrp%constructor (2, td%n_lvls, &
!                                        td%n1, td%n2, td%dx_c, td%dy_c)

        v_intrp = v_intrp%constructor (2, td%n_lvls, &
                                        td%n1, td%n2, td%dx_c, td%dy_c)

        s_intrp = s_intrp%constructor (2, td%n_lvls, &
                                        td%n1, td%n2, td%dx_c, td%dy_c)

        allocate (u_fine(4, td%n1*td%n2))
        do i = 1, td%n_lvls
            call v_intrp%get_interpolated ( i, td%n1, td%n2, &
                                            td%x(i,:), td%y(i,:), &
                                            td%u(i,:), td%v(i,:),  u_fine )
print*, "############"
            call s_intrp%get_interpolated ( i, td%n1, td%n2, &
                                            td%x(i,:), td%y(i,:), &
                                            td%u(i,:),  u_fine )
        end do

    end block


end program simplified_2
