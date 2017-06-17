SUBROUTINE getwage(Jinit,Uinit,indx,inds,indc,inde,wageval,Jtemp)
!
!   getwage.f90
!   Optimal_UI
!
!   For case with wage-dependent UI benefit, get wage associated with a contract V'(s')
!   to satisfy PK constraint for a given V
!   Created by Larry Warren on 10/6/16.
!   Copyright 2016 Larry Warren. All rights reserved.
!
USE PARAM
USE UTILITY
implicit none

!Dummy arguments declarations
real(8), dimension(nx,ny,nz), intent(in):: Jinit
real(8), dimension(ny,ne), intent(in):: Uinit
integer, intent(in):: indx
integer, intent(in):: inds
integer, intent(in):: indc
integer, intent(out):: inde
real(8), intent(out):: wageval
real(8), intent(out):: Jtemp

!Local variables declarations
integer:: j,indy,indz,indyp,indzp,indxp,indsp,witers
real(8)::wage0,wage1
real(8), dimension(ne)      :: tempe
real(8), dimension(nx,ny)   :: Jtilde
real(8):: EV0,EV1,EJ0,EJ1,Vp,c,dp,werror

indy = iyfun(inds)
indz = izfun(inds)
inde = ne
EV0 = zero
EJ0 = zero
wage0 = zero
!initially set to d at ne (ineligible or wage=0)
dp = dprimevec(indx,indy,ne)
witers=1
werror=one
do while(werror>1.0d-4 .and. witers<30)
  do indsp=1,ns
    indyp = iyfun(indsp)   !this gives the aggregate state in is'
    indzp = izfun(indsp)   !this gives the idiosyncratic state in is'
    indxp = cont(indc,indsp) !this is the index of V' in state is'
    Vp = x(indxp)        !this is V' in state s'
    dp = dprimevec(indxp,indyp,inde) !This is dprime

    !Should I interpolate U1 based on w?
    EV1 = EV1 + betta*ps(inds,indsp)*(dp*Uinit(indyp,inde) + (one-dp)*(Vp+lambda*R(indxp,indyp)))
    EJ1 = EJ1 + betta*ps(inds,indsp)*((one-dp)*(one-lambda*Ptilde(indxp,indyp))*Jinit(indxp,indyp,indzp))
  end do
  c = Um1(x(indx)-EV1)
  if (dabs(c-hell)<high_tol) then
    wage1 = 0.0d0
    Jtemp = hell
    inde = ne
  else
    wage1 = c/(one-tau)
    do j=1,ne
      tempe(j) = DABS(e(j)-wage1)
    end do
    inde = MINLOC(tempe, DIM=1)
    Jtemp = y(indy)+z(indz)-wage1 + EJ1
  end if
  werror = DABS(wage1-wage0)
  EV0=EV1
  EJ0=EJ1
  witers=witers+1
  wage0=wage1
  wageval = wage1
end do
RETURN
END SUBROUTINE getwage