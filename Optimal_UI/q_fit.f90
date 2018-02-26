SUBROUTINE q_fit(xin,yin,outvec)
!Simple 1-d linear fit: y= ax+b
! Adaptation of LLSQ, which solves a linear least squares problem matching a line to data.
! LLSQ is by John Burkardt

implicit none

real*8, dimension(:), intent(in) :: xin,yin
real*8, dimension(:), intent(out) :: outvec
real*8:: A_hat, X_hat, x1,xn,y1,yn
real*8:: a,b,c,bot,top,xbar,ybar
real*8, dimension(:), allocatable :: y_tilde,x_tilde
integer :: i,n

n = size(xin)
allocate(x_tilde(n),y_tilde(n))

x1 = xin(1)
xn = xin(n)
y1 = yin(1)
yn = yin(n)

A_hat = (y1-yn)/(x1-xn)
X_hat = (x1**2.0d0 - xn**2.0d0)/(x1-xn)

do i=1,n
y_tilde(i) = yin(i)-A_hat*xin(i)
x_tilde(i) = xin(i)**2.0d0 - X_hat*xin(i)
enddo

!
!  Average X_tilde and Y_tilde.
!
xbar = sum ( x_tilde(1:n) ) / real ( n, kind = 8 )
ybar = sum ( y_tilde(1:n) ) / real ( n, kind = 8 )
!
!  Compute Beta.
!
top = dot_product ( x_tilde(1:n) - xbar, y_tilde(1:n) - ybar )
bot = dot_product ( x_tilde(1:n) - xbar, x_tilde(1:n) - xbar )

a = top / bot

c = ybar - a * xbar

b = A_hat - a*X_hat

do i=1,n
outvec(i) = a*xin(i)**2.0d0 + b*xin(i) + c
if (isnan(outvec(i))) then
print*,'outvec',i,' is not a number'
pause
endif
enddo

deallocate(x_tilde,y_tilde)


END SUBROUTINE q_fit
