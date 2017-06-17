SUBROUTINE stadist(m,pimat,pss,initpss)
  !*************************************************************************
  !  stadist.f90
  !  
  !  PURPOSE:
  !   Given a square matrix pimat of dimension (m,m), computes the stationary
  !   distribution and stores it in the vector pss(m).
  !
  !
  !  Program written by: J. Fisher and M. Gervais
  !   Date:    April 2005
  !   Updated: May 2005
  !*************************************************************************
  implicit none
  
  !dummy arguments declarations
  integer, intent(in)                   :: m
  double precision, intent(in)          :: pimat(m,m)
  double precision, intent(out)         :: pss(m)
  double precision, intent(in), dimension(m), optional:: initpss
  
  !Local variables declarations
  integer                         :: i,j,iter
  integer, parameter              :: niter = 1000
  double precision, dimension (m) :: pn,pnp1
  double precision                :: psstot,temp
  double precision, parameter     :: tol = 0.1d-12
  double precision, parameter     :: zero = 0.0d0, one = 1.0d0
  
  pn   = zero
  pnp1 = zero
  temp = zero
  pss  = zero
  
  !Making sure that the rows of pimat sum to 1
  do i=1,m
!    temp = SUM(pimat(i,:))
    if (dabs(SUM(pimat(i,:))-one) .GT. tol) then
      write(*,'(''row '',i3,'' does not sum to 1'')')  i
      write(*,'(''SUM(mat('',i3,'',:) = '',f14.10)') i,SUM(pimat(i,:))
      pause
    end if
  end do
  
  !SOLVING FOR THE STATIONARY DISTRIBUTION OF PIMAT
  if (PRESENT(initpss)) then
    pnp1 = initpss
  else
    pnp1 = one/dble(m)
  end if
  temp = one
  iter = 1
  do while (iter.le.niter .and. dabs(temp)>tol)
    pn = pnp1
    temp = zero
    do j=1,m
      pnp1(j) = zero
      do i=1,m
        pnp1(j) = pnp1(j) + pn(i)*pimat(i,j)
      end do
      if (dabs(pnp1(j)-pn(j)).gt.temp) temp = pnp1(j)-pn(j)
    end do
    iter = iter+1
  end do
!  write (*,fmt='(''supnorm(pnp1-pn) = '',f10.8)') temp
  
  !Making sure pn has converged
  if (iter.ge.niter) then
    write (*,'(''Stationary distribution did not converge after '',&
               &i4,'' iterations:'',/,5x,''supnorm = '',f14.12)') iter,temp
  else
    write (*,'(''Stationary distribution converged after '',&
               &i3,'' iterations: supnorm = '',f14.12)') iter,temp
  end if
  
  pss = pnp1
  psstot = SUM(pss)
  
  if (dabs(psstot-one).gt.tol) then
    print *, 'THE STATIONARY DISTRIBUTION DOES NOT ADD UP TO ONE'
    print *, 'psstot = ',psstot
    STOP
  end if
  
  !LAST CARD OF SUBROUTINE stadist
  
  end SUBROUTINE stadist
