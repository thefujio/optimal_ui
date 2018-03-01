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

SUBROUTINE nearest_interp ( nd, xd, yd, ni, xi, yi )

!*****************************************************************************80
!
!! NEAREST_INTERP evaluates the nearest neighbor interpolant.
!  Discussion:
!    The nearest neighbor interpolant L(ND,XD,YD)(X) is the piecewise
!    constant function which interpolates the data (XD(I),YD(I)) for I = 1
!    to ND.
!  Licensing:
!    This code is distributed under the GNU LGPL license.
!  Modified:
!    04 September 2012
!  Author:
!    John Burkardt
!
!  Parameters:
!    Input, integer ( kind = 4 ) ND, the number of data points.
!    ND must be at least 1.
!    Input, real ( kind = 8 ) XD(ND), the data points.
!    Input, real ( kind = 8 ) YD(ND), the data values.
!    Input, integer ( kind = 4 ) NI, the number of interpolation points.
!    Input, real ( kind = 8 ) XI(NI), the interpolation points.
!    Output, real ( kind = 8 ) YI(NI), the interpolated values.
!
  implicit none

  integer ( kind = 4 ) nd
  integer ( kind = 4 ) ni

  real ( kind = 8 ) d
  real ( kind = 8 ) d2
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real ( kind = 8 ) xd(nd)
  real ( kind = 8 ) xi(ni)
  real ( kind = 8 ) yd(nd)
  real ( kind = 8 ) yi(ni)

  do i = 1, ni
    k = 1
    d = abs ( xi(i) - xd(k) )
    do j = 2, nd
      d2 = abs ( xi(i) - xd(j) )
      if ( d2 < d ) then
        k = j
        d = d2
      end if
    end do
    yi(i) = yd(k)
  end do
  return
END SUBROUTINE nearest_interp

subroutine pwl_value_1d ( nd, xd, yd, ni, xi, yi )

!*****************************************************************************80
!
!! PWL_VALUE_1D evaluates the piecewise linear interpolant.
!
!  Discussion:
!
!    The piecewise linear interpolant L(ND,XD,YD)(X) is the piecewise
!    linear function which interpolates the data (XD(I),YD(I)) for I = 1
!    to ND.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 September 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ND, the number of data points.
!    ND must be at least 1.
!
!    Input, real ( kind = 8 ) XD(ND), the data points.
!
!    Input, real ( kind = 8 ) YD(ND), the data values.
!
!    Input, integer ( kind = 4 ) NI, the number of interpolation points.
!
!    Input, real ( kind = 8 ) XI(NI), the interpolation points.
!
!    Output, real ( kind = 8 ) YI(NI), the interpolated values.
!
  implicit none

  integer ( kind = 4 ) nd
  integer ( kind = 4 ) ni

  integer ( kind = 4 ) i
  integer ( kind = 4 ) k
  real ( kind = 8 ) t
  real ( kind = 8 ) xd(nd)
  real ( kind = 8 ) yd(nd)
  real ( kind = 8 ) xi(ni)
  real ( kind = 8 ) yi(ni)

  yi(1:ni) = 0.0D+00

  if ( nd == 1 ) then
    yi(1:ni) = yd(1)
    return
  end if

  do i = 1, ni

    if ( xi(i) <= xd(1) ) then

      t = ( xi(i) - xd(1) ) / ( xd(2) - xd(1) )
      yi(i) = ( 1.0D+00 - t ) * yd(1) + t * yd(2)

    else if ( xd(nd) <= xi(i) ) then

      t = ( xi(i) - xd(nd-1) ) / ( xd(nd) - xd(nd-1) )
      yi(i) = ( 1.0D+00 - t ) * yd(nd-1) + t * yd(nd)

    else

      do k = 2, nd

        if ( xd(k-1) <= xi(i) .and. xi(i) <= xd(k) ) then

          t = ( xi(i) - xd(k-1) ) / ( xd(k) - xd(k-1) )
          yi(i) = ( 1.0D+00 - t ) * yd(k-1) + t * yd(k)
          exit

        end if

      end do

    end if

  end do
  
  return
end subroutine pwl_value_1d

SUBROUTINE ergodic(p,s)
! Purpose: Compute ergodic distribution s of Markov transition matrix p
IMPLICIT NONE
REAL*8, DIMENSION(:,:), INTENT(IN) :: p
REAL*8, DIMENSION(:), INTENT(OUT) :: s
REAL*8, DIMENSION(size(s),size(s)) :: ip,VL,VR
REAL*8, DIMENSION(size(s)) :: WR,WI,DW
INTEGER :: m,uw,w1(1),LWKOPT,INFO
INTEGER, PARAMETER :: LWORK = 50000
REAL*8 :: ds,DUMMY(1,1),WORK(LWORK)


m=size(s)
IF (size(p,dim=1)/=m .or. size(p,dim=2)/=m) THEN
  PRINT '(a,i3,a,i3)', 'sd: p must be a square matrix of size ',m,' x ',m
  STOP 'program terminated by sd'
END IF
ip=p
CALL DGEEV('V','V',m,ip,m,WR,WI,VL,m,VR,m,WORK,LWORK,INFO)
LWKOPT = WORK(1)
DW=abs(sqrt(WR*WR+WI*WI)-1)
w1=minloc(DW)
uw=count(DW<1000*epsilon(DW))
IF (uw<1) PRINT '(a)', 'Warning: No unitary eigenvalue is found. Stationary distribution of Markov chain does not exist.'
IF (uw>1) PRINT '(a)', 'Warning: More than one unitary eigenvalue is found. Stationary distribution of Markov chain is not unique.'
IF (uw<1 .or. uw>1) PRINT *, 'Using eigenvalue ', WR(w1(1)),'+i',WI(w1(1))
s=vl(:,w1(1))/sum(vl(:,w1(1)))
IF (any(s<0)) THEN
  PRINT '(a)', 'The stationary distribution of Markov chain has negative values. Rebalancing...'
  ds=sum(s,mask=s<0)/count(s>=0)
  WHERE(s<0)
  s=0
  ELSEWHERE
  s=s+ds
  ENDWHERE
END IF

END SUBROUTINE ergodic

END MODULE TOOLBOX