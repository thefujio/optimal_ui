!
!   write_output.f90
!   Optimal_UI
!
!   Created by Larry Warren on 10/3/16.
!   Copyright 2016 Larry Warren. All rights reserved.
!
SUBROUTINE write_output(U1,J1)

!  USE INTERFACES
USE IOOP
USE PARAM
implicit none

!Dummy arguments declarations
double precision, dimension(ny)      , intent(in):: U1
double precision, dimension(nx,ny,nz), intent(in):: J1



! CURRENT VALUE FUNCTIONS AND POLICY RULES
!real(8), dimension(nx,ny)   :: theta,P
!  real(8), dimension(nx,ny,nz):: J
!real(8), dimension(nx,ny,nz):: w
!real(8), dimension(nx,ns,ns):: dprime
!integer, dimension(nx,ns,ns):: iVprime
!real(8), dimension(nx,ny)   :: R,Ptilde
!integer, dimension(nx,ny)   :: M
!  real(8), dimension(ny)      :: U
!real(8), dimension(ny)      :: RU,PUtilde
!integer, dimension(ny)      :: MU

!Local variables declarations
integer:: i,ix,iy,iz,is,ic
integer:: ixp,iyp,izp,isp,ixpojs
real(8), dimension(nx,ny)   :: Jtilde
real(8):: Pojs,dp,PU

if (nz==1) then
Jtilde = J1(:,:,1)
else if (nz==2) then
Jtilde = zero
do iz=1,nz
Jtilde(:,:) = Jtilde(:,:) + pzss(iz)*J1(:,:,iz)
end do
else if (nz==3) then
Jtilde = J1(:,:,2)
end if


!do is=1,ns
!iy = iyfun(is)
!do ix=1,nx
!do isp=1,ns
!iyp = iyfun(isp)
!ixp = iVprime(ix,is,isp)
!dp  = dprime(ix,is,isp)
!ixpojs = M(ixp,iyp) !this is the market in which (ixp,iyp) searches
!Pojs   = lambda*Ptilde(ixpojs,iyp) !this is the prob of success of OJS
!end do
!end do
!end do

!dealing with unemployed individuals
do iy=1,ny
do iyp=1,ny
ixp = MU(iyp) !this is the market in which U(iyp) searches
PU = PUtilde(iyp) !this is the prob of success of search
do izp=1,nz
!transits to market ixp in state isp=(iyp,izp)
isp = isfun(iyp,izp)
end do
end do
end do

Call wri1file(ny,U1,root_dir//out_dir//"U.txt")
Call wri1file(ny,PUtilde,root_dir//out_dir//"PUtilde.txt")
Call wri1file(ny,RU,root_dir//out_dir//"RU.txt")
Call int1file(ny,MU,root_dir//out_dir//"MU.txt")

Call wri2file(nx,ny,P,root_dir//out_dir//"P.txt")
Call wri2file(nx,ny,theta,root_dir//out_dir//"theta.txt")
Call wri2file(nx,ny,R,root_dir//out_dir//"R.txt")
Call wri2file(nx,ny,Ptilde,root_dir//out_dir//"Ptilde.txt")
Call wri2file(nx,ny,Jtilde,root_dir//out_dir//"Jtilde.txt")
Call int2file(nx,ny,M,root_dir//out_dir//"M.txt")

!Call wri3file(nx,ny,nz,J1,root_dir//out_dir//"J.txt")
!Call wri3file(nx,ny,nz,w,root_dir//out_dir//"w.txt")

open(105,file=root_dir//out_dir//'empfuns.csv',status='replace')
  write(105,40)
  40 format('V,','J(y1;z1),','J(y1;z2),','J(y2;z1),','J(y2;z2),','w(y1;z1),','w(y1;z2),','w(y2;z1),','w(y2;z2),','P(y1),','P(y2),','th(y1),','th(y2),', &
    'R(y1),','R(y2),','Ptilde(y1),','Ptilde(y2),','Jtilde(y1),','Jtilde(y2)')
    do i=1,nx
      write(105,'(<1+ns*2+ny*5>(f15.4,","))') x(i),J1(i,1,1),J1(i,1,2),J1(i,2,1),J1(i,2,2),w(i,1,1),w(i,1,2),w(i,2,1),w(i,2,2),P(i,1),P(i,2),theta(i,1), &
      theta(i,2),R(i,1),R(i,2),Ptilde(i,1),Ptilde(i,2),Jtilde(i,1),Jtilde(i,2)
    enddo
close(105)

END SUBROUTINE write_output