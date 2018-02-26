SUBROUTINE centergrid(xgrid,xl,xu)
!
!   centergrid.f90
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
double precision, dimension(nx), intent(out)::xgrid
!Local variables declarations
 integer:: i
!integer:: isp,ixp,iyp,izp
!real(8), dimension(nx,ny,nz):: J0
!real(8):: norm0,norm1
real(8):: gridstep
integer:: ND1, NT0
gridstep = (xu-xl)/real(nx-1,8)
ND1=count(dprimevec(:,1,1) == 1.0d0)
if (want_print) then
print*, 'Submkts w/delta=1:',ND1
endif
NT0=count(theta(:,ny) == 0.0d0)
if (want_print) then
print*, 'Submkts w/theta=0:',NT0
endif
if (ND1 >= 20) then
xl = xl + 0.05d0*gridstep*real(ND1,8)
else if (ND1<=15) then
xl = xl - 2.0d0*(real(nx,8)/30.0d0)*gridstep
end if

if (NT0 >= 25) then
xu = xu - 0.05d0*gridstep*real(NT0,8)
else if (NT0<=20) then
xu = xu + 0.50d0*(real(nx,8)/30.0d0)*gridstep
end if

xgrid(1:nx) = (/ ( &
((xu-xl)/(real(nx-1,8)))*(real(i-1,8)) + xl, i=1,nx &
) /)
if (want_print) then
print*, 'gridstep= ',gridstep
print*, 'xmin = : ', xl, ' xmax = : ', xu
endif
RETURN
END SUBROUTINE centergrid
