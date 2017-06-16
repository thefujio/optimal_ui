REAL(8) FUNCTION GBD(tax,uivec)
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
  real(8):: tax
  real(8), dimension(*)   :: uivec
  
  !Local variables declarations
  integer:: is,ix,iy,iz,iu,ie
  real(8):: wagebill,uibill,urate
  
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
  
  RETURN
END FUNCTION GBD
