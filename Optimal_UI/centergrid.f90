SUBROUTINE centergrid(xgrid,xl,xu,U1)
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
double precision, dimension(ny), intent(in)::U1
!Local variables declarations
 integer:: i
!integer:: isp,ixp,iyp,izp
!real(8), dimension(nx,ny,nz):: J0
!real(8):: norm0,norm1
real(8):: gridstep
integer:: ND1, NT0
gridstep = (xu-xl)/real(nx-1,8)
ND1=count(dprimevec(:,1) == 1.0d0)
print*, 'Submkts w/delta=1:',ND1
NT0=count(theta(:,ny) == 0.0d0)
print*, 'Submkts w/theta=0:',NT0
print*, dprimevec(:,1)

if (ND1 >= 5) then
xl = xl + 0.1d0*gridstep*real(ND1,8)
else if (ND1<=2) then
xl = xl - 3.0d0*gridstep
end if
!
!if (NT0 >= 5) then
!xu = 0.015d0*gridstep*real(NT0,8)
!else if (NT0<=2) then
!xu = xu + 1.0d0*gridstep
!end if

xgrid(1:nx) = (/ ( &
((xu-xl)/(real(nx-1,8)))*(real(i-1,8)) + xl, i=1,nx &
) /)

print*, 'gridstep= ',gridstep
print*, 'xmin = : ', xl, ' xmax = : ', xu

RETURN
END SUBROUTINE centergrid
