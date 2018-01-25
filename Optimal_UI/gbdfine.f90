REAL(8) FUNCTION GBDFINE(tax,uivec)
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
    do ix=1,ndist
      wagebill = wagebill + muss((is-1)*ndist+ix)*wfine(ix,iy,iz)
    end do
  end do
  !Adds up obligations (b+hp for eligible, and hp for ineligible)
  uibill = zero
  urate = zero
  do iu=1,nu
    urate = urate + muss(ns*ndist+iu)
    ie = iuefun(iu)
    uibill = uibill + uivec(ie)*muss(ns*ndist+iu)
  end do
  
  GBDFINE = tax*wagebill - uibill
  
  RETURN
END FUNCTION GBDFINE
