SUBROUTINE howardfull(J1,U1)
  !****************************************************************************
  !  howardfull.f90
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
  real(8), dimension(ny)      , intent(inout):: U1
  
  !Local variables declarations
  integer:: iter,is,ix,iy,iz
  integer:: isp,ixp,iyp,izp
  real(8), dimension(nx,ny,nz):: J0
  real(8), dimension(ny)      :: U0
  real(8):: dp,EJ
  real(8):: norm0,norm1
  
  
  norm1 = zero
  do iter=1,nupdate
    U0 = U1
    do iy=1,ny
      U1(iy) = Ufunc(b) + betta*DOT_PRODUCT(py(iy,:),U0+RU)
    end do
    
    norm0 = norm1    
    J0 = J1
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
    if (norm1<high_tol) then
      EXIT
    else if (iter>1 .and. norm1>norm0) then
      write(*,"('Norm went up: No Howard this iteration')")
      EXIT
    end if
  end do
  
  RETURN
END SUBROUTINE howardfull
