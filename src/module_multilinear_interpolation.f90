module precision

    use, intrinsic :: iso_fortran_env, only: real32, real64, real128    !!!single, double, quadrable
    implicit none
    private
    integer, parameter, public :: rp = real64                    !.. working precision for reals
    real(rp), public :: small = 1.e-6_rp

end module


!!!Author H. Emrah Konokman
!! General, multi-linear interpolation
!  dintrv is borrowed from https://github.com/jacobwilliams/finterp
module multilinear_interpolation_mod

    use :: precision

    implicit none
    type :: param_type
        real(rp), allocatable :: params(:)
    end type

    type :: table_type
        type(param_type), allocatable :: x(:)   !!!Table parameters
        real(rp), allocatable :: f(:)           !!!Table line values
    contains
        procedure :: interpolation => multilinear_interpolation
        final :: delete_table
    end type

contains

!*****************************************************************************************
!>
!  Returns the indices in `xt` that bound `x`, to use for interpolation.
!  If outside the range, then the indices are returned that can
!  be used for extrapolation.
!  Precisely,
!
!```fortran
!         if            x < xt(1)   then ileft=1,   iright=2,    mflag=-1
!         if   xt(i) <= x < xt(i+1) then ileft=i,   iright=i+1,  mflag=0
!         if   xt(n) <= x           then ileft=n-1, iright=n,    mflag=1
!```
!
!### History
!
!  * interv written by carl de boor [5]
!  * dintrv author: amos, d. e., (snla) : date written 800901
!  * revision date 820801
!  * Jacob Williams, 2/24/2015 : updated to free-form Fortran.
!  * Jacob Williams, 2/17/2016 : additional refactoring (eliminated GOTOs).
!  * Jacob Williams, 2/22/2016 : modified bspline-fortran `dintrv` routine for
!    linear interpolation/extrapolation use.
!  * Jacob Williams, 10/9/2019 : added optional `inearest` output.

    pure subroutine dintrv(xt,x,ilo,ileft,iright,mflag,inearest)

        implicit none

        real(rp),dimension(:),intent(in) :: xt       !! a knot or break point vector
        real(rp),intent(in)              :: x        !! argument
        integer,intent(inout)            :: ilo      !! an initialization parameter which must be set
        !! to 1 the first time the array `xt` is
        !! processed by dintrv. `ilo` contains information for
        !! efficient processing after the initial call and `ilo`
        !! must not be changed by the user.  each dimension
        !! requires a distinct `ilo` parameter.
        integer,intent(out)              :: ileft    !! left index
        integer,intent(out)              :: iright   !! right index
        integer,intent(out)              :: mflag    !! signals when `x` lies out of bounds
        integer,intent(out),optional     :: inearest !! nearest index

        integer :: ihi, istep, imid, n

        n = size(xt)

        if (n==1) then
            ! this is only allowed for nearest interpolation
            if (present(inearest)) then
                inearest = 1
                return
            end if
        end if

        ihi = ilo + 1
        if ( ihi>=n ) then
            if ( x>=xt(n) ) then
                mflag = 1
                ileft = n-1
                iright= n
                if (present(inearest)) inearest = n
                return
            end if
            if ( n<=1 ) then
                mflag = -1
                ileft = 1
                iright= 2
                if (present(inearest)) inearest = 1
                return
            end if
            ilo = n - 1
            ihi = n
        endif

        if ( x>=xt(ihi) ) then

            ! now x >= xt(ilo). find upper bound
            istep = 1
            do
                ilo = ihi
                ihi = ilo + istep
                if ( ihi>=n ) then
                    if ( x>=xt(n) ) then
                        mflag = 1
                        ileft = n-1
                        iright= n
                        if (present(inearest)) inearest = n
                        return
                    end if
                    ihi = n
                elseif ( x>=xt(ihi) ) then
                    istep = istep*2
                    cycle
                endif
                exit
            end do

        else

            if ( x>=xt(ilo) ) then
                mflag = 0
                ileft = ilo
                iright= ilo+1
                if (present(inearest)) then
                    if ( abs(x-xt(ileft)) <= abs(x-xt(iright)) ) then
                        inearest = ileft
                    else
                        inearest = iright
                    end if
                end if
                return
            end if
            ! now x <= xt(ihi). find lower bound
            istep = 1
            do
                ihi = ilo
                ilo = ihi - istep
                if ( ilo<=1 ) then
                    ilo = 1
                    if ( x<xt(1) ) then
                        mflag = -1
                        ileft = 1
                        iright= 2
                        if (present(inearest)) inearest = 1
                        return
                    end if
                elseif ( x<xt(ilo) ) then
                    istep = istep*2
                    cycle
                endif
                exit
            end do

        endif

        ! now xt(ilo) <= x < xt(ihi). narrow the interval
        do
            imid = (ilo+ihi)/2
            if ( imid==ilo ) then
                mflag = 0
                ileft = ilo
                iright= ilo+1
                if (present(inearest)) then
                    if ( abs(x-xt(ileft)) <= abs(x-xt(iright)) ) then
                        inearest = ileft
                    else
                        inearest = iright
                    end if
                end if
                return
            end if
            ! note. it is assumed that imid = ilo in case ihi = ilo+1
            if ( x<xt(imid) ) then
                ihi = imid
            else
                ilo = imid
            endif
        end do

    end subroutine dintrv

    subroutine delete_table(self)

        type(table_type) :: self

        integer :: i

        if(allocated(self%x)) then
            do i = 1, size(self%x)
                if(allocated(self%x(i)%params)) deallocate(self%x(i)%params)
            end do
            deallocate(self%x)
        end if
        if(allocated(self%f)) deallocate(self%f)

    end subroutine

    !!!Interpolation function
    !! If input parameter is out of lower or upper limits then extrapolation is performed.
    pure function multilinear_interpolation(self, x) result(value)

        class(table_type), intent(in) :: self     !!!Table parameters
        real(rp), intent(in) :: x(:)              !!!Input parameters to interpolate
        real(rp) :: value                         !!!Interpolated value

        integer :: ix(size(self%x),2)
        integer :: mflag
        integer :: i, j, num_param, num_f
        integer :: indx, iopp
        integer :: ilo
        real(rp) :: vol_opp, vol_opp_sum
        integer :: indexf, fr
        real(rp) :: sign

        num_param = size(x)
        num_f = size(self%f)
        sign = 1._rp

        !!!Contribution of the bound value is proportional to the opposite volume
        !!!     *------------------------------*
        !!!     |                   |          |
        !!!     |     opposite      |          |
        !!!     |      volume       |Interp.val|
        !!!     |-------------------x----------|
        !!!     |                   |          |
        !!!     |                   |          |
        !!!     *------------------------------*bound value

        !!!Find lower and upper bounds of the table
        do i = 1, num_param
            ilo = 1
            call dintrv(self%x(i)%params, x(i), ilo, ix(i,1), ix(i,2), mflag)
        end do

        !!!Determine bound and opposite bound indices and calculate opposite volumes
        value = 0._rp
        vol_opp_sum = 0._rp
        do i = 1, 2**num_param

            vol_opp = 1._rp
            indexf = 0
            fr = num_f
            do j = 1, num_param
                !!!indx: Index of ix (from dintrv) to get the table bound value (lower/lef or upper/right bound):
                indx = abs(mod(int((i-1)/2**(num_param-j))+1,2) - 2)
                !!!iopp: Index of ix (from dintrv) to get the table parameters to calculate opposite volume
                iopp = 3 - indx
                !!!Opposite volume of the bound
                if(iopp == 1 .and. x(j) > self%x(j)%params(size(self%x(j)%params))) sign = -1._rp
                if(iopp == 2 .and. x(j) < self%x(j)%params(1)) sign = -1._rp
                vol_opp = vol_opp * sign * abs(self%x(j)%params(ix(j,iopp)) - x(j))
                sign = 1._rp
                !!!indexf: Index of the bound value (f)
                fr = fr/size(self%x(j)%params)
                indexf = indexf + (ix(j,indx) - 1) * fr
            end do

            indexf = indexf + 1

            !!!Contribution of the bound value is proportional to opposite volume
            value = value + vol_opp * self%f(indexf)
            vol_opp_sum = vol_opp_sum + vol_opp
        end do

        value = value/vol_opp_sum

    end function

end module
