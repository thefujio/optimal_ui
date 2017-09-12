!
!   write_output.f90
!   Optimal_UI
!
!   Created by Larry Warren on 10/3/16.
!   Copyright 2016 Larry Warren. All rights reserved.
!
SUBROUTINE write_output(J1,U1)

!  USE INTERFACES
USE IOOP
USE PARAM
implicit none

!Dummy arguments declarations
double precision, dimension(ny,ne)      , intent(in):: U1
double precision, dimension(nx,ny,nz), intent(in):: J1



! CURRENT VALUE FUNCTIONS AND POLICY RULES (for reference)
!real(8), dimension(nx,ny)   :: theta,P
!real(8), dimension(nx,ny,nz):: J
!real(8), dimension(nx,ny,nz):: w
!real(8), dimension(nx,ns,ns):: dprime
!integer, dimension(nx,ns,ns):: iVprime
!real(8), dimension(nx,ny)   :: R,Ptilde
!integer, dimension(nx,ny)   :: M
!real(8), dimension(ny)      :: U
!real(8), dimension(ny)      :: RU,PUtilde
!integer, dimension(ny)      :: MU

!Local variables declarations
integer:: i,ix,iy,iz,iu,is
real(8), dimension(nx,ny)   :: Jtilde
real(8), dimension(nx+nu):: muvec

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


Call wri1file(ny,U1,root_dir//out_dir//"U.txt")
Call wri1file(ny,PUtilde,root_dir//out_dir//"PUtilde.txt")
Call wri1file(ny,RU,root_dir//out_dir//"RU.txt")
Call int1file(ny,MU,root_dir//out_dir//"MU.txt")

Call wri2file(nx,ny,P,root_dir//out_dir//"P.txt")
Call wri2file(nx,ny,theta,root_dir//out_dir//"theta.txt")
Call wri2file(nx,ny,R,root_dir//out_dir//"R.txt")
Call wri2file(nx,ny,Ptilde,root_dir//out_dir//"Ptilde.txt")
Call wri2file(nx,ny,Jtilde,root_dir//out_dir//"Jtilde.txt")
Call wri2file(nx,ny,dprimevec,root_dir//out_dir//"dprime.txt")
Call int2file(nx,ny,M,root_dir//out_dir//"M.txt")

Call wri3file(nx,ny,nz,J1,root_dir//out_dir//"J.txt")
!Call wri3file(nx,ny,nz,w,root_dir//out_dir//"w.txt")
if (nz==2 .and. ny==2) then
  open(105,file=root_dir//out_dir//'empfuns.csv',status='replace')
    write(105,40)
    40 format('V,','J(y1;z1),','J(y1;z2),','J(y2;z1),','J(y2;z2),','w(y1;z1),','w(y1;z2),','w(y2;z1),','w(y2;z2),','dprime(y1),','dprime(y2),', &
      'P(y1),','P(y2),','th(y1),','th(y2),','R(y1),','R(y2),','Ptilde(y1),','Ptilde(y2),','Jtilde(y1),','Jtilde(y2),','M(y1),','M(y2)')
      do i=1,nx
        write(105,'(<1+ns*2+ny*6>(f15.4,","),2(I6,","))') x(i),J1(i,1,1),J1(i,1,2),J1(i,2,1),J1(i,2,2),w(i,1,1),w(i,1,2),w(i,2,1),w(i,2,2),dprimevec(i,1,1),dprimevec(i,2,1),&
        P(i,1),P(i,2),theta(i,1),theta(i,2),R(i,1),R(i,2),Ptilde(i,1),Ptilde(i,2),Jtilde(i,1),Jtilde(i,2),M(i,1),M(i,2)
      enddo
  close(105)
elseif (nz==1 .and. ny==2) then
  open(105,file=root_dir//out_dir//'empfuns.csv',status='replace')
    write(105,50)
    50 format('V,','J(y1),','J(y2),','w(y1),','w(y2),','dprime(y1),','dprime(y2),', &
      'P(y1),','P(y2),','th(y1),','th(y1),','R(y1),','R(y1),','Ptilde(y),','Ptilde(y2),','Jtilde(y1),','Jtilde(y2),','M(y1),','M(y2),')
      do i=1,nx
        write(105,'(<17>(f15.4,","),2(I6,","))') x(i),J1(i,1,1),J1(i,2,1),w(i,1,1),w(i,2,1),dprimevec(i,1,1),dprimevec(i,2,1),&
        P(i,1),P(i,2),theta(i,1),theta(i,2),R(i,1),R(i,2),Ptilde(i,1),Ptilde(i,2),Jtilde(i,1),Jtilde(i,2),M(i,1),M(i,2)
      enddo
  close(105)

elseif (nz==2 .and. ny==1) then
  open(105,file=root_dir//out_dir//'empfuns.csv',status='replace')
    write(105,60)
    60 format('V,','J(y1;z1),','J(y1;z2),','w(y1;z1),','w(y1;z2),','dprime(y1),', &
      'P(y1),','th(y1),','R(y1),','Ptilde(y1),','Jtilde(y1),','M(y1),')
      do i=1,nx
        write(105,'(<11>(f15.4,","),1(I6,","))') x(i),J1(i,1,1),J1(i,1,2),w(i,1,1),w(i,1,2),dprimevec(i,1,1),&
        P(i,1),theta(i,1),R(i,1),Ptilde(i,1),Jtilde(i,1),M(i,1)
      enddo
  close(105)

elseif (nz==3 .and. ny==1) then
  open(105,file=root_dir//out_dir//'empfuns.csv',status='replace')
    write(105,70)
    70 format('V,','J(y1;z1),','J(y1;z2),','J(y1;z3),','w(y1;z1),','w(y1;z2),','w(y1;z3),','dprime(y1),', &
      'P(y1),','th(y1),','R(y1),','Ptilde(y1),','Jtilde(y1),','M(y1),')
    do i=1,nx
      write(105,'(<13>(f15.4,","),1(I6,","))') x(i),J1(i,1,1),J1(i,1,2),J1(i,1,3),w(i,1,1),w(i,1,2),w(i,1,3),dprimevec(i,1,1),&
      P(i,1),theta(i,1),R(i,1),Ptilde(i,1),Jtilde(i,1),M(i,1)
    enddo
  close(105)
elseif (nz==3 .and. ny==2) then
  open(105,file=root_dir//out_dir//'empfuns.csv',status='replace')
    write(105,75)
    75 format('V,','J(y1;z1),','J(y1;z2),','J(y1;z3),','J(y2;z1),','J(y2;z2),','J(y2;z3),','w(y1;z1),','w(y1;z2),','w(y1;z3),','w(y2;z1),','w(y2;z2),','w(y2;z3),','dprime(y1),','dprime(y2),', &
      'P(y1),','P(y2),','th(y1),','th(y2),','R(y1),','R(y2),','Ptilde(y1),','Ptilde(y2),','Jtilde(y1),','Jtilde(y2),','M(y1),','M(y2),')
    do i=1,nx
      write(105,'(<25>(f15.4,","),2(I6,","))') x(i),J1(i,1,1),J1(i,1,2),J1(i,1,3),J1(i,2,1),J1(i,2,2),J1(i,2,3),&
      w(i,1,1),w(i,1,2),w(i,1,3),w(i,2,1),w(i,2,2),w(i,2,3),dprimevec(i,1,1),dprimevec(i,2,1), &
      P(i,1),P(i,2),theta(i,1),theta(i,2),R(i,1),R(i,2),Ptilde(i,1),Ptilde(i,2),Jtilde(i,1),Jtilde(i,2),M(i,1),M(i,2)
    enddo
  close(105)
else
  open(105,file=root_dir//out_dir//'empfuns.csv',status='replace')
    write(105,80)
    80 format('V,','J,','w,','dprime,', &
    'P,','th,','R,','Ptilde(y1),','Jtilde,','M,')
    do i=1,nx
      write(105,'(<9>(f15.4,","),1(I6,","))') x(i),J1(i,1,1),w(i,1,1),dprimevec(i,1,1),&
      P(i,1),theta(i,1),R(i,1),Ptilde(i,1),Jtilde(i,1),M(i,1)
    enddo
  close(105)
endif

if (nz==3 .and. ny==1) then
  open(107,file=root_dir//out_dir//'vprime.csv',status='replace')
    write(107,'((A6,","),<ns*ns>(I6,","))') 'is',1,1,1,2,2,2,3,3,3
    write(107,'((A6,","),<ns*ns>(I6,","))') 'isp',1,2,3,1,2,3,1,2,3
    write(107,'((A6,","),<ns*ns>(I6,","))') 'y(isp)', iyfun(1),iyfun(2),iyfun(3),iyfun(1),iyfun(2),iyfun(3),iyfun(1),iyfun(2),iyfun(3)
    write(107,'((A6,","),<ns*ns>(I6,","))') 'z(isp)', izfun(1), izfun(2),izfun(3),izfun(1),izfun(2),izfun(3),izfun(1),izfun(2),izfun(3)
    do i=1,nx
      write(107,'((f15.4,","),<ns*ns>(I6,","))') x(i),iVprime(i,1,1),iVprime(i,1,2),iVprime(i,1,3),iVprime(i,2,1),iVprime(i,2,2),iVprime(i,2,3),iVprime(i,3,1),iVprime(i,3,2),iVprime(i,3,3)
    enddo
  close(107)
else if (nz==2 .and. ny==1) then
  open(107,file=root_dir//out_dir//'vprime.csv',status='replace')
    write(107,'((A6,","),<2*ns>(I6,","))') 'is',1,1,2,2
    write(107,'((A6,","),<2*ns>(I6,","))') 'isp',1,2,1,2
    write(107,'((A6,","),<2*ns>(I6,","))') 'y(isp)', iyfun(1), iyfun(2),iyfun(1), iyfun(2)
    write(107,'((A6,","),<2*ns>(I6,","))') 'z(isp)', izfun(1), izfun(2), izfun(1), izfun(2)
    do i=1,nx
      write(107,'((f15.4,","),<2*ns>(I6,","))') x(i),iVprime(i,1,1),iVprime(i,1,2),iVprime(i,2,1),iVprime(i,2,2)
    enddo
  close(107)
else if (nz==1 .and. ny==2) then
  open(107,file=root_dir//out_dir//'vprime.csv',status='replace')
    write(107,'((A6,","),<2*ns>(I6,","))') 'is',1,1,2,2
    write(107,'((A6,","),<2*ns>(I6,","))') 'isp',1,2,1,2
    write(107,'((A6,","),<2*ns>(I6,","))') 'y(isp)', iyfun(1), iyfun(2),iyfun(1), iyfun(2)
    write(107,'((A6,","),<2*ns>(I6,","))') 'z(isp)', izfun(1), izfun(2), izfun(1), izfun(2)
    do i=1,nx
      write(107,'((f15.4,","),<2*ns>(I6,","))') x(i),iVprime(i,1,1),iVprime(i,1,2),iVprime(i,2,1),iVprime(i,2,2)
    enddo
  close(107)

else if(nz==1 .and. ny==1) then
  open(107,file=root_dir//out_dir//'vprime.csv',status='replace')
    write(107,'((A6,","),<2*ns>(I6,","))') 'is',1
    write(107,'((A6,","),<2*ns>(I6,","))') 'isp',1
    write(107,'((A6,","),<2*ns>(I6,","))') 'y(isp)', iyfun(1)
    write(107,'((A6,","),<2*ns>(I6,","))') 'z(isp)', izfun(1)
    do i=1,nx
      write(107,'((f15.4,","),<2*ns>(I6,","))') x(i),iVprime(i,1,1)
    enddo
  close(107)
endif

if (ne>2) then
  open(109,file=root_dir//out_dir//'ufuns.csv',status='replace')
    write(109,90)
    90 format('Past Wage,','U(y1),','U(y2),','wtilde(y1;z1),','wtilde(y1;z2),','wtilde(y2;z1),','wtilde(y2;z2),', &
'PUtilde(y1),','PUtilde(y2),','RU(y1),','RU(y2),','MU(y1),','MU(y2)')
    do i=1,ne
      write(109,'(<1+ny*nz+ny*3>(f15.4,","),2(I6,","))') e(i),U1(1,i),U1(2,i),w(MU(1,i),1,1),w(MU(1,i),1,2),w(MU(2,i),2,1),w(MU(2,i),2,2), &
          PUtilde(1,i),PUtilde(2,i),RU(1,i),RU(2,i),MU(1,i),MU(2,i)
    enddo
  close(109)
endif

!Stationary Distribution
Call wri1file(ns*nx+nu,muss,root_dir//out_dir//"muss.txt")
!Collapsed to Rows of length nx+1: U is 1st element of muvec
muvec = zero
do iu=1,nu
  muvec(iu) = muss(ns*nx+iu)
enddo

do is=1,ns
  do ix=1,nx
  muvec(nu+ix) = muvec(nu+ix) + muss((is-1)*nx+ix)
  enddo
enddo

open(106,file=root_dir//out_dir//'dist.csv',status='replace')
  write(106,30)
  30 format('V,','dist')
  do iu=1,nu
    write(106,'(2(f15.4,","))') U1(iuyfun(iu),iuefun(iu)),muvec(iu)
  enddo
  do ix=1,nx
    write(106,'(2(f15.4,","))') x(ix),muvec(ix+nu)
  enddo
close(106)

!other parameters, moments
params(1) = betta
params(2) = eta
params(3) = sigma
params(4) = gamma
params(5) = lambda
params(6) = delta
params(7) = kappa
params(8) = nx
params(9) = rr
params(10)= tau
params(11)= hp
params(12)= em
params(13)= unemp
params(14)= UEflow
params(15)= EEflow
params(16)= EUflow
params(17)= tot_wage
params(18)= avg_wage
params(19)= avg_benefit
params(20)= transfers
params(21)= welfare
call wri2file(nparams,1,params,root_dir//out_dir//"params.txt")

END SUBROUTINE write_output