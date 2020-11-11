module mod_math


contains



    function inverse_matrix (a, tol) result (a_inv)

        implicit none
        real(kind = 8), optional   :: tol
        real(kind = 8), intent(in) :: a(:,:)
        real(kind = 8) :: a_inv (size(a,1),size(a,2))
        real(kind = 8) :: a_copy(size(a,1),size(a,2))
        real(kind = 8) :: s(min(size(a,1), size(a,2)))
        real(kind = 8) :: ss(size(s),size(s))
        real(kind = 8) :: u(size(a,1),size(a,1))
        real(kind = 8) :: vt(size(a,2),size(a,2))
        real(kind = 8) :: my_tol
        real(kind = 8), allocatable :: work(:)
        integer(kind = 4) :: m, n, i
        integer(kind = 4) :: info
        integer(kind = 4) :: lda
        integer(kind = 4) :: ldu
        integer(kind = 4) :: ldv
        integer(kind = 4) :: ldvt
        integer(kind = 4) :: lwork
        character(len=1)  :: jobu, jobvt

        if (present(tol)) then
            my_tol = tol
        else
            my_tol = 10.d0 ** (-12.0d0)
        end if

        m = size(a, 1)
        n = size(a, 2)
        jobu  = "a"
        jobvt = "a"
        ldu = m
        ldv = n
        lda = max(m, 1)
        ldvt = n

        lwork = 5 * max( 3 * min(m,n) + max(m,n), 5 * min(m,n))
        allocate (work(lwork))

        a_copy(:,:) = a(:,:)
        call dgesvd ( jobu, jobvt, m, n, a_copy, lda, s, u, ldu, vt, ldvt,  work, &
                      lwork, info )

        ss = 0.d0
        do i = 1, size(s)
            if (s(i) >= my_tol) then
                ss(i,i) = 1.d0 / s(i)
            else
                ss(i,i) = 0.0d0
            end if
        end do

        a_inv = matmul(transpose(vt), matmul(ss, transpose(u)))

    end function inverse_matrix



end module mod_math
