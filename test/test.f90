program test

    use :: multilinear_interpolation_mod

    implicit none
    type(table_type) :: xt(3)
    real(rp), allocatable :: f(:)
    real(rp) :: x(3)
    real(rp) :: value, value_exact
    integer :: i, j, k, icount
    real(rp) :: coeffs(3) = [0.3_rp, 7._rp, 13._rp]

    allocate(xt(1)%params(5))
    xt(1)%params = [100._rp, 500._rp, 1000._rp, 1500._rp, 2000._rp]
    allocate(xt(2)%params(4))
    xt(2)%params = [15._rp, 25._rp, 35._rp, 45._rp]
    allocate(xt(3)%params(3))
    xt(3)%params = [7._rp, 9._rp, 11._rp]

    allocate(f(size(xt(1)%params)*size(xt(2)%params)*size(xt(3)%params)))

    !!!Write example table
    open(1, file = "example_table_3param.dat ")
    write(1,'(4a10)') "xt(1)", "xt(2)", "xt(3)", "f"
    icount = 0
    do i = 1, size(xt(1)%params)
        do j = 1, size(xt(2)%params)
            do k = 1, size(xt(3)%params)
                icount = icount + 1
                f(icount) = coeffs(1) * xt(1)%params(i) + &
                            coeffs(2) * xt(2)%params(j) + &
                            coeffs(3) * xt(3)%params(k)

                write(1,'(4f10.4)') xt(1)%params(i), xt(2)%params(j), xt(3)%params(k), f(icount)
            end do
        end do
    end do

    close(1)

    x = [575.39_rp, 39.72_rp, 7.27_rp]
    write(*,'(1a25,3f10.2)') "Input Params: ", x
    write(*,*) "----------------------------------------------------"
    value = multilinear_interpolation(xt, f, x)
    write(*,'(1a25,1f10.2)') "Interpolated value: ", value
    value_exact = coeffs(1) * x(1) + coeffs(2) * x(2) + coeffs(3) * x(3)
    write(*,'(1a25,1f10.2)') "Exact value: ", value_exact
    write(*,'(1a25,1l10)') "Interpolation correct: ", abs(value_exact - value) < small
    write(*,*) "----------------------------------------------------"

    do i = 1, size(xt)
        deallocate(xt(i)%params)
    end do

end program
