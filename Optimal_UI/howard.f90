SUBROUTINE howard(J1)
  !****************************************************************************
  !  howard.f90
  !  
  !  PURPOSE:
  !   Iterate on value functions given policy rules
  !
  !
  !  Program written by: M. Gervais
  !   Date:    Feb 2016
  !   Updated: Feb 2016
  !****************************************************************************
  USE PARAM
  USE UTILITY
  implicit none
  
  !Dummy arguments declarations
  real(8), dimension(nx,ny,nz), intent(inout):: J1
!  real(8), dimension(nx,ny), intent(in)   :: RU
!  real(8), dimension(ny)   , intent(in)   :: U
!  real(8), dimension(nx,ns), intent(in)   :: w
!  real(8), dimension(nx,ns,ns), intent(in):: dprime
!  integer, dimension(nx,ns,ns), intent(in):: Vp
  
  !Local variables declarations
  integer:: iter
  integer:: is,ix,iy,iz
  integer:: isp,ixp,iyp,izp
  real(8), dimension(nx,ny,nz):: J0
  real(8):: EJ,dp
  real(8):: norm0,norm1
  
  
  norm1 = zero
  do iter=1,nupdate
    norm0 = norm1
    J0 = J1
    J1 = zero
    do ix=1,nx
      do is=1,ns
        iy = iyfun(is)
        iz = izfun(is)
        if (dabs(J0(ix,iy,iz)-hell)<high_tol) CYCLE
        EJ = zero
        do isp=1,ns
          iyp = iyfun(isp)   !this gives the aggregate state in is
          izp = izfun(isp)   !this gives the idiosyncratic state in is
          ixp = iVprime(ix,is,isp) !this is the index of V' in state s'
!          Vp = x(ixp) !this is V' in state s
          dp = dprime(ix,is,isp) !This is d'
!          EV = EV + betta*ps(is,isp)*(dp*U(iyp) + (1-dp)*(Vp+lambda*R(ixp,iyp)))
          EJ = EJ + betta*ps(is,isp)*((one-dp)*(one-lambda*Ptilde(ixp,iyp))*J0(ixp,iyp,izp))
        end do
        J1(ix,iy,iz) = y(iy)+z(iz)-w(ix,iy,iz) + EJ
      end do
    end do
    norm1 = MAXVAL(dabs(J1-J0))
    print*, 'howard iteration: ', iter
    if (norm1<high_tol) then
      EXIT
    else if (iter>1 .and. norm1>norm0) then
      write(*,"('Norm went up: No Howard this iteration')")
      EXIT
    end if
  end do
  
  RETURN
END SUBROUTINE howard
