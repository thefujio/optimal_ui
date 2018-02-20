SUBROUTINE SDIFINE(J1,U1)
  !****************************************************************************
  !  SDI.f90
  !  
  !  PURPOSE:
  !   Finds the statoinary distribution given a tax code (tau,b)
  !   Calculate aggregate moments
  !  
  !  Program written by: M. Gervais
  !   Date:    Feb 2016
  !   Updated: Feb 2016
  !****************************************************************************
!  USE INTERFACES
  USE IOOP
  USE PARAM
  USE UTILITY
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
  double precision, dimension(ny,ne)      , intent(in):: U1
  double precision, dimension(nx,ny,nz), intent(in):: J1
  
  !Local variables declarations
  integer:: is,ix,iy,iz,iu,ie
  integer:: isp,ixp,iyp,izp,iup,ixpojs
  real(8), dimension(ne):: tempe
  real(8), dimension(ns*ndist+nu):: muinit
  real(8):: Pojs,dp,PU, ubenmeasure, unobenmeasure
  
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
    iz = izfun(is)
    do ix=1,ndist
      do isp=1,ns
        iyp = iyfun(isp)
        ixp = iVprimefine(ix,is,isp)
        dp  = dprimefine(ix,is,isp)
        ixpojs = Mfine(ixp,iyp) !this is the market in which (ixp,iyp) searches
        Pojs   = lambda*Ptildefine(ixpojs,iyp) !this is the prob of success of OJS
        !print*, 'pojs', Pojs
        !print*, 'dp', dp
        !print*, 'ixpojs', ixpojs
        !pause
        !transits to contract Vprime
        pimat((is-1)*ndist+ix,(isp-1)*ndist+ixp) = pimat((is-1)*ndist+ix,(isp-1)*ndist+ixp) + &
            ps(is,isp)*(1.0d0-dp)*(1.0d0-Pojs)
        !transits to OJS market
        pimat((is-1)*ndist+ix,(isp-1)*ndist+ixpojs) = pimat((is-1)*ndist+ix,(isp-1)*ndist+ixpojs) + &
            ps(is,isp)*(1.0d0-dp)*Pojs
        !transits to unemployment (and eligible for UI)
        pimat((is-1)*ndist+ix,ns*ndist+iufun(iyp,windfine(ix,iy,iz))) = pimat((is-1)*ndist+ix,ns*ndist + &
            iufun(iyp,windfine(ix,iy,iz))) + ps(is,isp)*dp
      end do
    end do
  end do

  !dealing with unemployed individuals
  do iu=1,nu
    do iup=1,nu
      ixp = MUfine(iuyfun(iup),iuefun(iup)) !this is the market in which U(iyp) searches
      PU = PUtildefine(iuyfun(iup),iuefun(iup)) !this is the prob of success of search
      !transits back to unemployment
      pimat(ns*ndist+iu,ns*ndist+iup) = pimat(ns*ndist+iu,ns*ndist+iup) + pus(iu,iup)*(one-PU)
      do izp=1,nz
        !transits to market ixp in state isp=(iyp,izp)
        isp = isfun(iuyfun(iup),izp)
        pimat(ns*ndist+iu,(isp-1)*ndist+ixp) = pimat(ns*ndist+iu,(isp-1)*ndist+ixp) + pus(iu,iup)*Pztilde(izp)*PU
      end do
    end do
  end do

  !print*, pimat(201,:)
  !print*, sum(pimat(201,:))
  !pause

  !setting initial distribution to all unemployed
  muinit = zero
  muinit(ns*ndist+1:) = one/dble(nu)

  call stadist(ns*ndist+nu,pimat,muss,muinit)

  !Aggregate Statistics: flows
  ee=zero
  eu=zero
  ue=zero
  submktval=zero
  submktmeasure = zero
  submktwgt = zero
  do is=1,ns
    do ix=1,ndist
      do isp=1,ns
        eu = eu + ps(is,isp)*dprimevecfine(ix,iyfun(isp),windfine(ix,iyfun(is),izfun(is)))*muss((is-1)*ndist+ix)
        ee = ee + &
        ps(is,isp)*(1.0d0-dprimevecfine(ix,iyfun(isp),windfine(ix,iyfun(is),izfun(is))))*lambda*Ptildefine(ix,iyfun(isp))*muss((is-1)*ndist+ix)
        submktwgt = submktwgt + xfine(Mfine(ix,iyfun(isp)))*ps(is,isp)* &
        (1.0d0-dprimevecfine(ix,iyfun(isp),windfine(ix,iyfun(is),izfun(is))))*lambda*Ptildefine(ix,iyfun(isp))*muss((is-1)*ndist+ix)
        submktmeasure = submktmeasure + ps(is,isp)* &
        (1.0d0-dprimevecfine(ix,iyfun(isp),windfine(ix,iyfun(is),izfun(is))))*lambda*Ptildefine(ix,iyfun(isp))*muss((is-1)*ndist+ix)
      end do
    end do
  end do

  do iu=1,nu
    do iup=1,nu
      ue = ue + pus(iu,iup)*PUtildefine(iuyfun(iup),iuefun(iup))*muss(ns*ndist+iu)
      submktwgt = submktwgt + xfine(MUfine(iuyfun(iup),iuefun(iup)))*pus(iu,iup)*PUtildefine(iuyfun(iup),iuefun(iup))*muss(ns*ndist+iu)
      submktmeasure = submktmeasure + pus(iu,iup)*PUtildefine(iuyfun(iup),iuefun(iup))*muss(ns*ndist+iu)
    end do
  end do
  submktval = submktwgt/submktmeasure
  unemp=SUM(muss(ns*ndist+1:))
  em=SUM(muss(1:ns*ndist))
  UEflow=ue/unemp
  EEflow=ee/em
  EUflow=eu/em

  !Aggregate Statistics: stocks
  tot_output = zero
  tot_wage = zero
  tot_util = zero
  tot_vf = zero
  do is=1,ns
    iy = iyfun(is)
    iz = izfun(is)
    do ix=1,ndist
      tot_output = tot_output + muss((is-1)*ndist+ix)*(y(iyfun(is))+z(izfun(is)))
      tot_wage = tot_wage + muss((is-1)*ndist+ix)*wfine(ix,iy,iz)
      tot_util = tot_util + muss((is-1)*ndist+ix)*Ufunc(wfine(ix,iy,iz))
      tot_vf = tot_vf + muss((is-1)*ndist+ix)*xfine(ix)
    end do
  end do
  avg_wage = tot_wage/em

  if (ne>2) then
    do ie=1,ne
      tempe(ie) = DABS(e(ie)-avg_wage)
    end do
    avg_wage_ind = MINLOC(tempe, DIM=1)
  endif

  transfers = zero
  welfare = zero
  uval = zero
  uwgt = zero
  umeasure = zero
  do iu=1,nu
  transfers = transfers + muss(ns*ndist+iu)*bvec(iuefun(iu))
  tot_util = tot_util + muss(ns*ndist+iu)*Ufunc(bvec(iuefun(iu)))
  tot_vf = tot_vf + muss(ns*ndist+iu)*Ufine(iuyfun(iu),iuefun(iu))
  welfare = welfare + muss(ns*ndist+iu)*Ufine(iuyfun(iu),iuefun(iu))
  uwgt = uwgt + muss(ns*ndist+iu)*Ufine(iuyfun(iu),iuefun(iu))
  umeasure = umeasure + muss(ns*ndist+iu)
  end do

  avg_benefit = transfers/unemp

  uval = uwgt/umeasure
  unobenval = zero
  ubenval = zero
  do iu=1,nu
    ie = iuyfun(iu)
    ie = iuefun(iu)
    if (ie < ne) then
      ubenval = ubenval + muss(ns*ndist+iu)*Ufine(iy,ie)
      ubenmeasure = ubenmeasure + muss(ns*ndist+iu)
    else
      unobenval = unobenval + muss(ns*ndist+iu)*Ufine(iy,ie)
      unobenmeasure = unobenmeasure + muss(ns*ndist+iu)
    end if
  end do

  if (nu>1) then
    ubenval = ubenval/ubenmeasure
    unobenval = unobenval/unobenmeasure
    else
    ubenval = Ufine(1,1)
    unobenval = Ufine(1,1)
  endif

  do is=1,ns
    do ix=1,ndist
    welfare = welfare + xfine(ix)*muss((is-1)*ndist+ix)
    end do
  end do

  print*,'job-finding probability: ',UEflow
  print*,'job-to-job flow: ',EEflow
  print*,'job destruction rate: ',EUflow
  print*,'unemployment rate: ', unemp
  print*,'output: ', tot_output
  print*,'average wage: ', avg_wage
  print*,'welfare: ', welfare
  print*,'vf     : ', tot_vf
  print*,'inst. util.: ', tot_util


  !Grid Output
  !uval, submktval calculated in loops above
  ceval = welfare
  trval = avg_benefit
  rrval = avg_benefit/avg_wage !avg b/ avg wage
  taxval = tau
  jfpval = UEflow
  grosswageval = avg_wage
  netwageval = avg_wage*(1.0d0-tau)
  urateval = unemp
  uuval = 1.0d0-UEflow !1.0d0 - UEflow The current charts are the percent of population who are in uu,ee, not the flow prob!
  eeval = EEflow


  RETURN
END SUBROUTINE SDIFINE
