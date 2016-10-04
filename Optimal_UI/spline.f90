!   spline.f90
!   Optimal_UI
!   From http://ww2.odu.edu/~agodunov/computing/programs/book2/Ch01/spline.f90
!   Copyright 2016 Larry Warren. All rights reserved.
!
SUBROUTINE spline (xvec, yvec, b, c, d, n)
!======================================================================
!  Calculate the coefficients b(i), c(i), and d(i), i=1,2,...,n
!  for cubic spline interpolation
!  s(x) = y(i) + b(i)*(x-x(i)) + c(i)*(x-x(i))**2 + d(i)*(x-x(i))**3
!  for  x(i) <= x <= x(i+1)
!  Alex G: January 2010
!----------------------------------------------------------------------
!  input..
!  xvec = the arrays of data abscissas (in strictly increasing order)
!  yvec = the arrays of data ordinates
!  n = size of the arrays xveci() and yi() (n>=2)
!  output..
!  b, c, d  = arrays of spline coefficients
!  comments ...
!  spline.f90 program is based on fortran version of program spline.f
!  the accompanying function fspline can be used for interpolation
!======================================================================
implicit none
integer:: n
double precision, dimension(n):: xvec, yvec, b, c, d
integer:: i, j, gap
double precision:: h

gap = n-1
! check input
if ( n < 2 ) return
if ( n < 3 ) then
b(1) = (yvec(2)-yvec(1))/(xvec(2)-xvec(1))   ! linear interpolation
c(1) = 0.
d(1) = 0.
b(2) = b(1)
c(2) = 0.
d(2) = 0.
return
end if
!
! step 1: preparation
!
d(1) = xvec(2) - xvec(1)
c(2) = (yvec(2) - yvec(1))/d(1)
do i = 2, gap
d(i) = xvec(i+1) - xvec(i)
b(i) = 2.0d0*(d(i-1) + d(i))
c(i+1) = (yvec(i+1) - yvec(i))/d(i)
c(i) = c(i+1) - c(i)
end do
!
! step 2: end conditions
!
b(1) = -d(1)
b(n) = -d(n-1)
c(1) = 0.0d0
c(n) = 0.0d0
if(n /= 3) then
c(1) = c(3)/(xvec(4)-xvec(2)) - c(2)/(xvec(3)-xvec(1))
c(n) = c(n-1)/(xvec(n)-xvec(n-2)) - c(n-2)/(xvec(n-1)-xvec(n-3))
c(1) = c(1)*d(1)**2/(xvec(4)-xvec(1))
c(n) = -c(n)*d(n-1)**2/(xvec(n)-xvec(n-3))
end if
!
! step 3: forward elimination
!
do i = 2, n
h = d(i-1)/b(i-1)
b(i) = b(i) - h*d(i-1)
c(i) = c(i) - h*c(i-1)
end do
!
! step 4: back substitution
!
c(n) = c(n)/b(n)
do j = 1, gap
i = n-j
c(i) = (c(i) - d(i)*c(i+1))/b(i)
end do
!
! step 5: compute spline coefficients
!
b(n) = (yvec(n) - yvec(gap))/d(gap) + d(gap)*(c(gap) + 2.0d0*c(n))
do i = 1, gap
b(i) = (yvec(i+1) - yvec(i))/d(i) - d(i)*(c(i+1) + 2.0d0*c(i))
d(i) = (c(i+1) - c(i))/d(i)
c(i) = 3.0d0*c(i)
end do
c(n) = 3.0*c(n)
d(n) = d(n-1)
END SUBROUTINE spline

double precision function ispline(u, xvec, yvec, b, c, d, n)
!======================================================================
! function ispline evaluates the cubic spline interpolation at point z
! ispline = y(i)+b(i)*(u-x(i))+c(i)*(u-x(i))**2+d(i)*(u-x(i))**3
! where  x(i) <= u <= x(i+1)
!----------------------------------------------------------------------
! input..
! u       = the abscissa at which the spline is to be evaluated
! xvec, yvec    = the arrays of given data points
! b, c, d = arrays of spline coefficients computed by spline
! n       = the number of data points
! output:
! ispline = interpolated value at point u
!=======================================================================
implicit none
integer:: n
double precision::  u, xvec(n), yvec(n), b(n), c(n), d(n)
integer:: i, j, k
double precision:: dx

! if u is ouside the xvec() interval take a boundary value (left or right)
if(u <= xvec(1)) then
ispline = yvec(1)
return
end if
if(u >= xvec(n)) then
ispline = yvec(n)
return
end if

!*
!  binary search for for i, such that xvec(i) <= u <= xvec(i+1)
!*
i = 1
j = n+1
do while (j > i+1)
k = (i+j)/2
if(u < xvec(k)) then
j=k
else
i=k
end if
end do
!*
!  evaluate spline interpolation
!*
dx = u - xvec(i)
ispline = yvec(i) + dx*(b(i) + dx*(c(i) + dx*d(i)))
end function ispline