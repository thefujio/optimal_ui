  Subroutine agg_moments
  !****************************************************************************
  !  GBD.f90
  !
  !  PURPOSE:
  !   Returns the government budget deficit in a long run kind of fashion
  !
  !
  !  Program written by: M. Gervais
  !   Date:    Feb 2016
  !   Updated: Feb 2016
  !****************************************************************************
  USE IOOP
  USE PARAM
  implicit none

  !Dummy arguments declarations

  !Local variables declarations
  integer:: is,ix,iy,iu,ie

  do iu=1,nu
  uval = uval + U1(iuyfun(iu),iuefun(iu))*muss(ns*nx+iu)
  end do

  wagebill = zero
  do is=1,ns
  iy = iyfun(is)
  iz = izfun(is)
  do ix=1,nx
  wagebill = wagebill + muss((is-1)*nx+ix)*w(ix,iy,iz)
  end do
  end do
  !Adds up obligations (b+hp for eligible, and hp for ineligible)
  uibill = zero
  urate = zero
  do iu=1,nu
  urate = urate + muss(ns*nx+iu)
  ie = iuefun(iu)
  uibill = uibill + uivec(ie)*muss(ns*nx+iu)
  end do

  GBD = tax*wagebill - uibill

  !Grid Output
    ceval = welfare
    taxval = tau
    jfpval = UEflow
    uval = 0.0d0
    submktval = 0.0d0
    grosswageval = avgwage
    netwageval = avgwage*(1.0d0-tax)
    urateval = urate
    uuval = 1.0d0 - UEflow
    eeval = EEflow

  end subroutine agg_moments
