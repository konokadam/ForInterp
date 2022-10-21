program test

    use :: multilinear_interpolation_mod

    implicit none
    type(table_type) :: xt(3)
    real(rp), allocatable :: f(:)
    real(rp) :: x(3)
    real(rp) :: value

    allocate(xt(1)%params(5))
    xt(1)%params = [13._rp, 17._rp, 23._rp, 29._rp, 31._rp]
    allocate(xt(2)%params(4))
    xt(2)%params = [2._rp, 6._rp, 11._rp, 12._rp]
    allocate(xt(3)%params(3))
    xt(3)%params = [7._rp, 9._rp, 11._rp]

    allocate(f(size(xt(1)%params)*size(xt(2)%params)*size(xt(3)%params)))
    f(1) = 0.4754_rp
    f(2) = 0.6112_rp
    f(3) = 0.7471_rp
    f(4) = 1.4263_rp
    f(5) = 1.8338_rp
    f(6) = 2.2413_rp
    f(7) = 2.6149_rp
    f(8) = 3.362_rp
    f(9) = 4.1091_rp
    f(10) = 2.8526_rp
    f(11) = 3.6677_rp
    f(12) = 4.4827_rp
    f(13) = 0.6217_rp
    f(14) = 0.7993_rp
    f(15) = 0.977_rp
    f(16) = 1.8652_rp
    f(17) = 2.3981_rp
    f(18) = 2.931_rp
    f(19) = 3.4195_rp
    f(20) = 4.3965_rp
    f(21) = 5.3735_rp
    f(22) = 3.7304_rp
    f(23) = 4.7962_rp
    f(24) = 5.862_rp
    f(25) = 0.8411_rp
    f(26) = 1.0815_rp
    f(27) = 1.3218_rp
    f(28) = 2.5235_rp
    f(29) = 3.2445_rp
    f(30) = 3.9655_rp
    f(31) = 4.6264_rp
    f(32) = 5.9482_rp
    f(33) = 7.2701_rp
    f(34) = 5.047_rp
    f(35) = 6.489_rp
    f(36) = 7.931_rp
    f(37) = 1.0606_rp
    f(38) = 1.3636_rp
    f(39) = 1.6666_rp
    f(40) = 3.1818_rp
    f(41) = 4.0909_rp
    f(42) = 5_rp
    f(43) = 5.8333_rp
    f(44) = 7.5_rp
    f(45) = 9.1666_rp
    f(46) = 6.3636_rp
    f(47) = 8.1818_rp
    f(48) = 10_rp
    f(49) = 1.1337_rp
    f(50) = 1.4576_rp
    f(51) = 1.7816_rp
    f(52) = 3.4012_rp
    f(53) = 4.373_rp
    f(54) = 5.3448_rp
    f(55) = 6.2356_rp
    f(56) = 8.0172_rp
    f(57) = 9.7988_rp
    f(58) = 6.8025_rp
    f(59) = 8.746_rp
    f(60) = 10.6896_rp

    x = [34.39_rp, 2.75_rp, 1.37_rp]
    write(*,*) "Input Params: ", x
    write(*,*) "----------------------------------------------------"
    value = multilinear_interpolation(xt, f, x)
    write(*,*) "Interpolated value: ", value
    write(*,*) "----------------------------------------------------"

end program
