program test

    use :: multilinear_interpolation_mod

    implicit none
    type(table_type) :: table
    real(rp) :: x(3)
    real(rp) :: value, value_exact
    integer :: i, j, k, icount
    real(rp) :: coeffs(3) = [0.3_rp, 7._rp, 13._rp]
    real(rp) :: rand1, rand2, rand3
    character(13) :: x_(3)
    logical :: passed

    allocate(table%x(3))
    allocate(table%x(1)%params(5))
    table%x(1)%params = [100._rp, 500._rp, 1000._rp, 1500._rp, 2000._rp]
    allocate(table%x(2)%params(4))
    table%x(2)%params = [15._rp, 25._rp, 35._rp, 45._rp]
    allocate(table%x(3)%params(3))
    table%x(3)%params = [7._rp, 9._rp, 11._rp]

    allocate(table%f(size(table%x(1)%params)*size(table%x(2)%params)*size(table%x(3)%params)))

    !!!Write example table
    open(1, file = "example_table_3param.dat ")
    write(1,'(1a5,4a12)') "No", "xt(1)", "xt(2)", "xt(3)", "f"
    icount = 0
    do i = 1, size(table%x(1)%params)
        do j = 1, size(table%x(2)%params)
            do k = 1, size(table%x(3)%params)
                icount = icount + 1
                table%f(icount) = coeffs(1) * table%x(1)%params(i) + &
                            coeffs(2) * table%x(2)%params(j) + &
                            coeffs(3) * table%x(3)%params(k)

                write(1,'(1i5,4f12.4)') icount, table%x(1)%params(i), table%x(2)%params(j), table%x(3)%params(k), table%f(icount)
            end do
        end do
    end do

    close(1)

    open(152, file = "ForInterp.log")
    passed = .true.
    do i = 1, 100
        call random_number(rand1)
        call random_number(rand2)
        call random_number(rand3)
        rand1 = rand1 * (1.5_rp * table%x(1)%params(size(table%x(1)%params)) - 0._rp) + 0._rp
        rand2 = rand2 * (1.5_rp * table%x(2)%params(size(table%x(2)%params)) - 0._rp) + 0._rp
        rand3 = rand3 * (1.5_rp * table%x(3)%params(size(table%x(3)%params)) - 0._rp) + 0._rp
        x = [rand1, rand2, rand3]
        write(152,'(1a25,3f15.2)') "Input Params: ", x
        do j = 1, size(x)
            if(x(j) < table%x(j)%params(1)) then
                x_(j) = "out of lower"
            else if(x(j) > table%x(j)%params(size(table%x(j)%params))) then
                x_(j) = "out of upper"
            else
                x_(j) = "in"
            end if
        end do
        write(152,'(1a25,3a15)') "", (trim(x_(j)), j = 1, size(x))
        write(152,*) "---------------------------------------------------------------------"
        value = table%interpolation(x)
        write(152,'(1a25,1f15.2)') "Interpolated value: ", value
        value_exact = coeffs(1) * x(1) + coeffs(2) * x(2) + coeffs(3) * x(3)
        write(152,'(1a25,1f15.2)') "Exact value: ", value_exact
        write(152,'(1a25,1l15)') "Interpolation correct: ", abs(value_exact - value) < small
        write(152,*) "_____________________________________________________________________"
        passed = passed .and. abs(value_exact - value) < small
    end do
    write(152,*) ""
    write(152,*) "All tests passed =", passed
    write(*,*) "All tests passed =", passed
    close(152)


end program
