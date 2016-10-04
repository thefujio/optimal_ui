SUBROUTINE centergrid(xl,xu)
!
!   refine.f90
!   Optimal_UI
!
!   Created by Larry Warren on 9/28/16.
!   Copyright 2016 Larry Warren. All rights reserved.
!
USE PARAM
USE UTILITY
implicit none
!Dummy Variable Declarations
double precision, intent(inout)::xl,xu

!Local variables declarations
!integer:: is,ix,iy,iz
!integer:: isp,ixp,iyp,izp
!real(8), dimension(nx,ny,nz):: J0
!real(8):: norm0,norm1
integer:: ND1, NT0
real(8):: gridstep

ND1=count(dprime(:,1,1) == 1.0d0)
print*, ND1
NT0=count(theta(:,ny) == 0.0d0)
print*,NT0
gridstep = (xu-xl)/real(nx-1,8)
print*, 'gridstep= ',gridstep
if (ND1 >= 5) then
xl = xl + 0.15d0*gridstep*real(ND1,8)
else if (ND1<=2) then
xl = xl - 3.0d0*gridstep
else
xl = xl
end if

if (NT0 >= 5) then
xu = xu - 0.15d0*gridstep*real(NT0,8)
else if (NT0<=2) then
xu = xu + 3.0d0*gridstep
else
xu = xu
end if
print*, xl, xu
RETURN
END SUBROUTINE centergrid
