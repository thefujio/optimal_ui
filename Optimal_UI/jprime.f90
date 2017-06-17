SUBROUTINE JPRIME(U1,J1)
  !****************************************************************************
  !  JPRIME.f90
  !  
  !  PURPOSE:
  !   Finds the value functions J1 and U1 and policy rules given a tax code 
  !   (tau,b) and current value functions J0 and U0
  !     
  !  
  !  Program written by: M. Gervais
  !   Date:    Feb 2016
  !   Updated: Feb 2016
  !****************************************************************************
  USE IOOP
  USE PARAM
  USE UTILITY
  implicit none
  
  !Dummy arguments declarations
!  double precision, intent(in):: tau,b
  double precision, dimension(ny,ne)      , intent(inout):: U1
  double precision, dimension(nx,ny,nz), intent(inout):: J1

  INTERFACE
    SUBROUTINE howard(J1)
      USE PARAM
      USE UTILITY
      implicit none
      !Dummy arguments declarations
      real(8), dimension(nx,ny,nz), intent(inout):: J1
    END SUBROUTINE howard

    SUBROUTINE howard_full(J1,U1)
      USE PARAM
      USE UTILITY
      implicit none
      !Dummy arguments declarations
      real(8), dimension(nx,ns), intent(inout):: J1
      real(8), dimension(ny,ne)   , intent(inout):: U1
    END SUBROUTINE howard_full
  END INTERFACE
  
  !Local variables declarations
  integer:: i,ii,ie,iu,ix,iy,iz,is,ic
  integer:: ixp,iyp,izp,isp,iup
  real(8), dimension(nu):: ExpU
  real(8), dimension(nx):: ret
  real(8), dimension(ny,ne):: U0
  real(8), dimension(nx,ny,nz):: J0
  integer, dimension(nx,ny,nz):: iJ1
  real(8), dimension(nx,ny)   :: Jtilde
  real(8), dimension(nc):: Jtemp,wage
  real(8):: EV,EJ,V,Vp,dp,c
  
  
  !Setting current value functions and resetting U1 and J1
  U0 = U1
  J0 = J1
  U1 = zero
  J1 = zero

  !Below I assume that there is initial idiosyncratic uncertainty. If there
  !I need to take a stand on what the value of getting into
  !a new match is for firms. For example, in the paper all new matches
  !have z=0/iz=2 when nz=3. It may be just as well to flip a (fair) coin once the match
  !is created, i.e. integrate over possible states for z (assign pzss*J0()).
 	
  !Given J, find Tightness and JFP
  if (nz==1) then
    Jtilde = J0(:,:,1)
  else if (nz==2) then !starts at the exp. value
    Jtilde = zero
    do iz=1,nz
      Jtilde(:,:) = Jtilde(:,:) + pzss(iz)*J0(:,:,iz)
    end do
  else if (nz==3) then
    Jtilde = J0(:,:,2) !no initial uncertainty, starts at midpoint
    !Jtilde(:,:) = Jtilde(:,:) + pzss(iz)*J0(:,:,iz) !Flip a coin
  end if
  
  theta = zero
  P = zero
  do iy=1,ny
    do ix=1,nx
      if (Jtilde(ix,iy)>kappa) then
        theta(ix,iy) = (((Jtilde(ix,iy)/kappa)**(gamma))-1)**(1/gamma)
        P(ix,iy) = theta(ix,iy)*(one+theta(ix,iy)**gamma)**(-one/gamma)
      else
        theta(ix,iy) = zero
        P(ix,iy) = zero
      end if
    end do
  end do
 	
  !Given U and P, find optimal market for unemployed
  RU = zero
  MU = 0
  PUtilde = zero
  do ie=1,ne
    do iy=1,ny
      do ix=1,nx
        ret(ix) = P(ix,iy)*(x(ix)-U0(iy,ie))
      end do
      RU(iy,ie) = MAXVAL(ret)
      MU(iy,ie) = MAXLOC(ret,DIM=1)
      PUtilde(iy,ie) = P(MU(iy,ie),iy)
    end do
  end do
  !Update the value of unemployment
  do i=1,nu
    ExpU(i) = zero
    do ii=1,nu
      ExpU(i) = ExpU(i) + pus(i,ii)*(U0(iuyfun(ii),iuefun(ii))+RU(iuyfun(ii),iuefun(ii)))
    end do
    U1(iuyfun(i),iuefun(i)) = Ufunc(bvec(iuefun(i))) + ExpU(i)
  end do
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! could check convergence of U here !!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  !Given V (some x) and P, find optimal OJS market
  R = zero
  Ptilde = zero
  do iy=1,ny
    do ix=1,nx
      M(ix,iy) = ix
      ret = zero
      V = x(ix)
      do ixp=ix,nx !can start at ix as ret is negative below ix
        ret(ixp) = P(ixp,iy)*(x(ixp)-V)
      end do
      R(ix,iy) = MAXVAL(ret)
        if (R(ix,iy)>0.0d0) then
        M(ix,iy) = MAXLOC(ret,DIM=1)
        endif
      Ptilde(ix,iy) = P(M(ix,iy),iy)
    end do
  end do
	
  !Given R and U, construct dprime(x,y)
  dprimevec = delta
  do ix=1,nx
    do iy=1,ny
      if (U1(iy,1)>x(ix)+lambda*R(ix,iy)) then
        dprimevec(ix,iy,:) = 1.0d0
      end if
    end do
  end do
	
  !Find optimal contract and update J
  !Here cont is a matrix with a contract in each row, 
  !  i.e. each row has the index of vprime in state s, s=1,ns
  iJ1 = 0
  iVprime = 0
  w = zero
  do is=1,ns
    iy = iyfun(is)
    iz = izfun(is)
    do ix=1,nx
      do ic=1,nc !nc is the number of possible contracts
        EV = zero
        EJ = zero
        do isp=1,ns
          iyp = iyfun(isp)   !this gives the aggregate state in is'
          izp = izfun(isp)   !this gives the idiosyncratic state in is'
          ixp = cont(ic,isp) !this is the index of V' in state is'
          Vp = x(ixp)        !this is V' in state s'
          dp = dprimevec(ixp,iyp,1) !This is dprime
          EV = EV + betta*ps(is,isp)*(dp*U1(iyp,1) + (1-dp)*(Vp+lambda*R(ixp,iyp)))
          EJ = EJ + betta*ps(is,isp)*((one-dp)*(one-lambda*Ptilde(ixp,iyp))*J0(ixp,iyp,izp))
        end do
        c = Um1(x(ix)-EV)
        if (dabs(c-hell)<high_tol) then
          wage(ic) = zero
          Jtemp(ic) = hell
        else
          wage(ic) = c/(one-tau)
          Jtemp(ic) = y(iy)+z(iz)-wage(ic) + EJ
        end if
      end do
      J1(ix,iy,iz)  = MAXVAL(Jtemp)
      iJ1(ix,iy,iz) = MAXLOC(Jtemp,DIM=1)
      w(ix,iy,iz)   = wage(iJ1(ix,iy,iz))
      wind(ix,iy,iz) = 1
      iVprime(ix,is,:) = cont(iJ1(ix,iy,iz),:)
      do isp=1,ns
        iyp = iyfun(isp)
        ixp = cont(iJ1(ix,iy,iz),isp)
        dprime(ix,is,isp) = dprimevec(ixp,iyp,1)
      end do
    end do
  end do
 	
  if (do_howard) then 
    call howard(J1)
  end if

RETURN
END SUBROUTINE JPRIME
