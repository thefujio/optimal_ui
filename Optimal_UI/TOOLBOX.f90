MODULE TOOLBOX
  USE PARAM
  implicit none
  
  !private
  public:: restrict,unrestrict
  double precision :: yval,ypval,linmat
CONTAINS

  FUNCTION repmat(A, M, N)
  !Replicate A M times in dim1, N times in dim2
  !Mimics Matlab's repmat for 2-D arrays
  integer, INTENT(IN) :: M, N
  real*8, INTENT(IN) :: A(:,:)
  real*8 :: repmat(SIZE(A,1)*M, SIZE(A,2)*N)
  integer :: n1A, n2A
  integer :: i, j, n11, n12, n21, n22

  n1A = SIZE(A,1); n2A = SIZE(A,2)

  do i = 1,M,1
    n11 = 1+( i-1 )*n1A
    n12 = n11+n1A-1

    do j = 1,N,1
      n21 = 1+( j-1 )*n2A
      n22 = n21+n2A-1
      repmat( n11:n12, n21:n22) = A
    enddo ! j loop
  enddo ! i loop
  END FUNCTION repmat

  SUBROUTINE linspace(linmat, a, b, n)
  ! returns a linearly spaced vector with n points in [a, b] with the option
  implicit none
  integer, intent(in) :: n
  Real*8, intent(in) :: a
  Real*8, intent(in) :: b
  Real*8, intent(out) :: linmat(n)
  Real*8 :: dx
  integer :: i

  dx = (b-a) / (real(n,8)-1.00d0)

  linmat = [(i*dx+a, i = 0, n-1)]

END SUBROUTINE linspace

SUBROUTINE grid(x,xmin,xmax,s)
  ! Purpose: Generate grid x on [xmin,xmax] using spacing parameter s set as follows:
  ! s=1		linear spacing
  ! s>1		left skewed grid spacing with power s
  ! 0<s<1		right skewed grid spacing with power s
  ! s<0		geometric spacing with distances changing by a factor -s^(1/(n-1)), (>1 grow, <1 shrink)
  ! s=-1		logarithmic spacing with distances changing by a factor (xmax-xmin+1)^(1/(n-1))
  ! s=0		logarithmic spacing with distances changing by a factor (xmax/xmin)^(1/(n-1)), only if xmax,xmin>0
  implicit none
	real*8, DIMENSION(:), INTENT(OUT) :: x
	real*8, INTENT(IN) :: xmin,xmax,s
	real*8 :: c ! growth rate of grid subintervals for logarithmic spacing
	integer :: n,i
	n=size(x)
	forall(i=1:n) x(i)=(real(i,8)-1.00d0)/real(n-1,8)
	if (s>0.00d0) then
		x=x**s*(xmax-xmin)+xmin
		if (s==1.00d0) then
!			PRINT '(a,i8,a,f6.3,a,f6.3,a)', 'Using ',n,' equally spaced grid points over domain [',xmin,',',xmax,']'
		else
!			PRINT '(a,i8,a,f6.3,a,f6.3,a,f6.3,a)', 'Using ',n,' skewed spaced grid points with power ',s,' over domain [',xmin,',',xmax,']'
		endif
	else
		if (s==-1.000d0) then
			c=xmax-xmin+1.00d0
!		ELSEIF (s==0.0_WP) THEN
!			IF (xmin>0.0_WP) THEN
!				c=xmax/xmin
!			ELSE
!				STOP 'grid: can not use logarithmic spacing for nonpositive values'
!			END IF
		else
			c=-s
		endif
!		PRINT '(a,i8,a,f6.3,a,f6.3,a,f6.3,a)', 'Using ',n,' logarithmically spaced grid points with growth rate ',c,' over domain [',xmin,',',xmax,']'
		x=((xmax-xmin)/(c-1.00d0))*(c**x)-((xmax-c*xmin)/(c-1.00d0))
	endif
END SUBROUTINE grid

! 'restrict' function provides a transformation of y from the real line
! to the interval [x1,x2] 
FUNCTION restrict(y,x1,x2, spread)
	REAL*8 :: restrict, spread
    REAL*8, INTENT(IN) :: x1, x2, y
	!print *, 'restrict'
    restrict=(x2-x1)*atan(y/spread)/pie+(x2+x1)*0.500d0
	!print *, 'ok restrict'
END FUNCTION restrict

! 'unrestrict' function provides a transformation of y from the interval [x1,x2]
! to the real line 
FUNCTION unrestrict(y,x1,x2, spread)
    REAL*8 :: unrestrict, spread
    REAL*8, INTENT(IN) :: x1, x2, y
	!print *, 'unrestrict'
    unrestrict=spread*tan(pie*(2.00d0*y-x2-x1)/(2.00d0*(x2-x1)))
	!print *, 'ok unrestrict'
END FUNCTION unrestrict

END MODULE TOOLBOX