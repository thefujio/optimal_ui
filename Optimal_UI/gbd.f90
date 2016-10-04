REAL(8) FUNCTION GBD(tax,ui)
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
  real(8):: tax,ui
  
  !Local variables declarations
  integer:: is,ix,iy,iz
  real(8):: wagebill,urate
  
  wagebill = zero
  do is=1,ns
    iy = iyfun(is)
    iz = izfun(is)
    do ix=1,nx
      wagebill = wagebill + muss((is-1)*nx+ix)*w(ix,iy,iz)
    end do
  end do
  
  urate = zero
  do iy=1,ny
    urate = urate + muss(ns*nx+iy)
  end do
  
  GBD = tax*wagebill - ui*urate
  
  RETURN
END FUNCTION GBD
