SUBROUTINE SDI(J1,U1)
  !****************************************************************************
  !  SDI.f90
  !  
  !  PURPOSE:
  !   Finds the statoinary distribution given a tax code (tau,b)
  !     
  !  
  !  Program written by: M. Gervais
  !   Date:    Feb 2016
  !   Updated: Feb 2016
  !****************************************************************************
!  USE INTERFACES
  USE IOOP
  USE PARAM
  implicit none
  
  INTERFACE
    SUBROUTINE stadist(m,pimat,pss,initpss)
      implicit none
      !dummy arguments declarations
      integer, intent(in)                   :: m
      double precision, intent(in)          :: pimat(m,m)
      double precision, intent(out)         :: pss(m)
      double precision, intent(in), dimension(m), optional:: initpss
    END SUBROUTINE stadist
  END INTERFACE
  
  !Dummy arguments declarations
  double precision, dimension(ny,ne)      , intent(inout):: U1
  double precision, dimension(nx,ny,nz), intent(inout):: J1
  
  !Local variables declarations
  integer:: is,ix,iy,iz,iu
  integer:: isp,ixp,iyp,izp,iup,ixpojs
  real(8), dimension(ns*nx+nu):: muinit
  real(8):: Pojs,dp,PU
  
  pimat = zero
  !Building transition matrix
    !pimat is (nx+1)*ny + nx*nz by (nx+1)*ny + nx*nz, or
    !pimat is nx*ns+ny by nx*ns+ny
    !is=1,ix=1
    !...
    !is=1,ix=nx
    !is=2,ix=1
    !... and so on...
    !is=ns,ix=nx
    !iy=1,ix=U1 (unemployed)
    !...
    !iy=ny,ix=Uny (unemployed)
 	
  !dealing with employed individuals
  do is=1,ns
    iy = iyfun(is)
    do ix=1,nx
      do isp=1,ns
        iyp = iyfun(isp)
        ixp = iVprime(ix,is,isp)
        dp  = dprime(ix,is,isp)
        ixpojs = M(ixp,iyp) !this is the market in which (ixp,iyp) searches
        Pojs   = lambda*Ptilde(ixpojs,iyp) !this is the prob of success of OJS
        !transits to contract Vprime
        pimat((is-1)*nx+ix,(isp-1)*nx+ixp) = ps(is,isp)*(1.0d0-dp)*(1.0d0-Pojs)
        !transits to OJS market
        pimat((is-1)*nx+ix,(isp-1)*nx+ixpojs) = &
            pimat((is-1)*nx+ix,(isp-1)*nx+ixpojs) + ps(is,isp)*(1.0d0-dp)*Pojs
        !transits to unemployment (and eligible for UI)
        pimat((is-1)*nx+ix,ns*nx+iufun(iyp,1)) = pimat((is-1)*nx+ix,ns*nx+iufun(iyp,1)) + ps(is,isp)*dp
      end do
    end do
  end do

  !dealing with unemployed individuals
  do iu=1,nu
    do iup=1,nu
      ixp = MU(iuyfun(iup),iuefun(iup)) !this is the market in which U(iyp) searches
      PU = PUtilde(iuyfun(iup),iuefun(iup)) !this is the prob of success of search
      !transits back to unemployment
      pimat(ns*nx+iu,ns*nx+iup) = pus(iu,iup)*(one-PU)
      do izp=1,nz
        !transits to market ixp in state isp=(iyp,izp)
        isp = isfun(iuyfun(iup),izp)
        pimat(ns*nx+iu,(isp-1)*nx+ixp) = pus(iu,iup)*Pztilde(izp)*PU
      end do
    end do
  end do
  
  
  !setting initial distribution to all unemployed
  muinit = zero
  muinit(ns*nx+1:) = one/dble(nu)
  
  call stadist(ns*nx+nu,pimat,muss,muinit)

  !Aggregate Statistics
  ee=zero
  eu=zero
  ue=zero
  do is=1,ns
    do ix=1,nx
      do isp=1,ns
        eu = eu + ps(is,isp)*dprimevec(ix,iyfun(isp))*muss((is-1)*nx+ix)
        ee = ee + &
        ps(is,isp)*(1.0d0-dprimevec(ix,iyfun(isp)))*lambda*Ptilde(ix,iyfun(isp))*muss((is-1)*nx+ix)
      end do
    end do
  end do

  do iu=1,nu
    do iup=1,nu
      ue = ue + pus(iu,iup)*PUtilde(iuyfun(iup),iuefun(iup))*muss(ns*nx+iu)
    end do
  end do

  unemp=SUM(muss(ns*nx+1:))
  em=SUM(muss(1:ns*nx))
  UEflow=ue/unemp
  EEflow=ee/em
  EUflow=eu/em

  tot_wage = zero
  do is=1,ns
    iy = iyfun(is)
    iz = izfun(is)
    do ix=1,nx
      tot_wage = tot_wage + muss((is-1)*nx+ix)*w(ix,iy,iz)
    end do
  end do
  avg_wage = tot_wage/em
  tranfers = b*unemp

  welfare = zero
  do iu=1,nu
    welfare = welfare + U1(iuyfun(iu),iuefun(iu))*muss(ns*nx+iu)
  enddo

  do is=1,ns
    do ix=1,nx
    welfare = welfare + x(ix)*muss((is-1)*nx+ix)
    enddo
  enddo

  print*,'job-finding probability: ',UEflow
  print*,'job-to-job flow: ',EEflow
  print*,'job destruction rate: ',EUflow
  print*,'unemployment rate: ', unemp
  print*,'average wage: ', avg_wage
  print*,'welfare: ', welfare
  RETURN
END SUBROUTINE SDI
