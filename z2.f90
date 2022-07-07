program main
        use iso_fortran_env
        implicit none
        integer, parameter :: dp = REAL64

        integer, parameter :: N = 50
        integer, dimension(N,N,N,N,4) :: latt

        real(dp) :: beta
        real(dp) :: dbeta
        real(dp) :: actionS
        integer :: i1
        integer :: NT

        ! coldstart
        latt = 1

        NT = 50
        dbeta = 1.0/NT

        ! cold -> hot
        do i1 = 1, NT
                beta = 1. - i1*dbeta
                call sweep(N, beta, latt, actionS)
                print*, beta, actionS
        end do

        ! hot -> cold
        do i1 = 1, NT
                beta = i1*dbeta
                call sweep(N, beta, latt, actionS)
                print*, beta, actionS
        end do

        contains

        subroutine moveup(xvec, d, N)
                implicit none
                integer, intent(in) :: d
                integer, intent(in) :: N
                integer, dimension(4), intent(inout) :: xvec
                xvec(d) = xvec(d) + 1
                if(xvec(d) > N) then
                        xvec(d) = xvec(d) - N
                else
                end if
        end subroutine moveup

        subroutine movedown(xvec, d, N)
                implicit none
                integer, intent(in) :: d
                integer, intent(in) :: N
                integer, dimension(4), intent(inout) :: xvec
                xvec(d) = xvec(d) - 1
                if(xvec(d) < 1) then
                        xvec(d) = xvec(d) + N
                else
                end if
        end subroutine movedown

        subroutine sweep(N, beta, lat, actionS)
                
                use iso_fortran_env
                implicit none
                integer, parameter :: dp = REAL64

                integer, intent(in) :: N
                real(dp), intent(in) :: beta
                integer, dimension(N,N,N,N,4), intent(inout) :: lat
                real(dp), intent(out) :: actionS

                integer :: i1, i2, i3, i4, d, dperp
                real(dp) :: r_num
                real(dp), dimension(N,N,N,N,4):: rand_arr

                integer, dimension(4) :: xvec
                real(dp) :: staple
                real(dp) :: staplesum
                real(dp) :: bplus
                real(dp) :: bminus

                call random_seed()
                call random_number(rand_arr)

                actionS = 0.
                do i1 = 1, N
                do i2 = 1, N
                do i3 = 1, N
                do i4 = 1, N
                do d = 1, 4

                xvec = [ i1, i2, i3, i4 ]

                staplesum = 0.
                do dperp = 1, 4
                if (dperp /= d) then
                call movedown(xvec, dperp, N)
                staple = lat(xvec(1), xvec(2), xvec(3), xvec(4), dperp)
                staple = staple*lat(xvec(1), xvec(2), xvec(3), xvec(4), d)
                call moveup(xvec, d, N)
                staple = staple*lat(xvec(1), xvec(2), xvec(3), xvec(4), dperp)
                call moveup(xvec, dperp, N)
                staplesum = staplesum + staple

                staple = lat(xvec(1), xvec(2), xvec(3), xvec(4), dperp)
                call moveup(xvec, dperp, N)
                call movedown(xvec, d, N)
                staple = staple*lat(xvec(1), xvec(2), xvec(3), xvec(4), d)
                call movedown(xvec, dperp, N)
                staple = staple*lat(xvec(1), xvec(2), xvec(3), xvec(4), dperp)
                staplesum = staplesum + staple
                else
                end if
                end do

                bplus = exp(beta * staplesum)
                bminus = exp(-beta * staplesum)
                bplus = bplus / (bplus + bminus)

                r_num = rand_arr(i1, i2, i3, i4, d)
                !call random_number(r_num)

                if (r_num <= bplus) then
                        lat(i1, i2, i3, i4, d) = 1
                        actionS = actionS + staplesum/N**4/24.
                else
                        lat(i1, i2, i3, i4, d) = -1
                        actionS = actionS - staplesum/N**4/24.
                end if

                end do
                end do
                end do
                end do
                end do

                actionS = 1. - actionS

        end subroutine sweep

end program main
