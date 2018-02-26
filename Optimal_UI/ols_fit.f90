SUBROUTINE ols_fit(xin,yin,outvec)
!Simple 1-d linear fit: y= ax+b
! Adaptation of LLSQ, which solves a linear least squares problem matching a line to data.
! LLSQ is by John Burkardt

implicit none
real*8, dimension(:), intent(in) :: xin,yin
real*8, dimension(:), intent(out) :: outvec
real*8:: a,b,bot,top,xbar,ybar
integer :: i,n

n = size(xin)
!
!  Special case.
!
if ( n == 1 ) then
a = 0.0D+00
b = yin(1)
return
end if
!
!  Average X and Y.
!
xbar = sum ( xin(1:n) ) / real ( n, kind = 8 )
ybar = sum ( yin(1:n) ) / real ( n, kind = 8 )
!
!  Compute Beta.
!
top = dot_product ( xin(1:n) - xbar, yin(1:n) - ybar )
bot = dot_product ( xin(1:n) - xbar, xin(1:n) - xbar )

a = top / bot

b = ybar - a * xbar

do i=1,n
outvec(i) = a * xin(i) + b
enddo

END SUBROUTINE ols_fit
