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
  real(8), dimension(nx,ns), intent(inout):: J1
  real(8), dimension(ny)   , intent(inout):: U1
!  real(8), dimension(ny)   , intent(in)   :: RU
!  real(8), dimension(nx,ns), intent(in)   :: w
!  real(8), dimension(nx,ns,ns), intent(in):: dprime
  
  !Local variables declarations
  integer:: iter,is,ix,iy,iz
  integer:: isp,ixp,iyp
  real(8), dimension(nx,ns):: J0
  real(8), dimension(ny)   :: U0
  real(8):: dp
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
        EJ = zero
        do isp=1,ns
          iyp = iyfun(isp)   !this gives the aggregate state in is
!          izp = izfun(isp)   !this gives the idiosyncratic state in is
          ixp = Vprime(ix,is,isp) !this is the index of V' in state s'
!          Vp = x(ixp) !this is V' in state s
          dp = dprime(ix,is,isp) !This is d'
!          EV = EV + betta*ps(is,isp)*(dp*U(iyp) + (1-dp)*(Vp+lambda*R(ixp,iyp)))
          EJ = EJ + betta*ps(is,isp)*((one-dp)*(one-lambda*Ptilde(ixp,iyp))*J0(ixp,isp))
        end do
        J1(ic) = y(iy)+z(iz)-w(ix,is) + EJ
      end do
    end do
    norm1 = MAXVAL(dabs(J1-J0))
    if (norm1<high_tol) then
      J = J1
      EXIT
    else if (iter>1 .and. norm1>norm0) then
      write(*,"('Norm went up: No Howard this iteration')")
      EXIT
    end if
  end do
  
  J = J1
  
  RETURN
END SUBROUTINE howardfull
