!
!   write_output_fine.f90
!   Optimal_UI
!
!   Created by Larry Warren on 10/3/16.
!   Copyright 2016 Larry Warren. All rights reserved.
!
SUBROUTINE write_output_fine(J1,U1)

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

! Value and Policy Functions on fine grid
! real(8), dimension(ndist,ny,ne) :: dprimevecfine
! real(8), dimension(ndist,ny,nz) :: wfine
! integer, dimension(ndist,ny,nz) :: windfine
! real(8), dimension(ndist,ns,ns) :: dprimefine
! integer, dimension(ndist,ns,ns) :: iVprimefine
! real(8), dimension(ndist,ny)    :: Rfine,Ptildefine
! integer, dimension(ndist,ny)    :: Mfine
! !real(8), dimension(ny,ne)       :: RUfine,PUtildefine
! integer, dimension(ny,ne)       :: MUfine


!Local variables declarations
integer:: i,ix,iy,iz,iu,is
real(8), dimension(ndist,ny)   :: Jtildefine
real(8), dimension(ndist+nu):: muvec

if (nz==1) then
  Jtildefine = Jfine(:,:,1)
else if (nz==2) then
  Jtildefine = zero
  do iz=1,nz
    Jtildefine(:,:) = Jtildefine(:,:) + pzss(iz)*Jfine(:,:,iz)
  end do
else if (nz==3) then
  Jtildefine = Jfine(:,:,2)
end if


Call wri1file(ny,U1,root_dir//out_dir//"U.txt")
Call wri1file(ny,PUtilde,root_dir//out_dir//"PUtilde.txt")
Call wri1file(ny,RU,root_dir//out_dir//"RU.txt")
Call int1file(ny,MU,root_dir//out_dir//"MU.txt")
Call wri2file(nx,ny,P,root_dir//out_dir//"P.txt")
Call wri2file(nx,ny,theta,root_dir//out_dir//"theta.txt")
Call wri2file(nx,ny,R,root_dir//out_dir//"R.txt")
Call wri2file(nx,ny,Ptilde,root_dir//out_dir//"Ptilde.txt")
!Call wri2file(nx,ny,Jtilde,root_dir//out_dir//"Jtilde.txt")
Call wri2file(nx,ny,dprimevec,root_dir//out_dir//"dprime.txt")
Call int2file(nx,ny,M,root_dir//out_dir//"M.txt")
Call wri3file(nx,ny,nz,J1,root_dir//out_dir//"J.txt")


Call wri1file(ny,PUtildefine,root_dir//out_dir//"PUtildefine.txt")
Call wri1file(ny,RUfine,root_dir//out_dir//"RUfine.txt")
Call int1file(ny,MUfine,root_dir//out_dir//"MUfine.txt")
Call wri2file(ny,ne,Ufine,root_dir//out_dir//"Ufine.txt")
Call wri2file(ndist,ny,Pfine,root_dir//out_dir//"Pfine.txt")
Call wri2file(ndist,ny,thetafine,root_dir//out_dir//"thetafine.txt")
Call wri2file(ndist,ny,Rfine,root_dir//out_dir//"Rfine.txt")
Call wri2file(ndist,ny,Ptildefine,root_dir//out_dir//"Ptildefine.txt")
Call wri2file(ndist,ny,Jtildefine,root_dir//out_dir//"Jtildefine.txt")
Call wri2file(ndist,ny,dprimevecfine,root_dir//out_dir//"dprimefine.txt")
Call int2file(ndist,ny,Mfine,root_dir//out_dir//"Mfine.txt")

Call wri1file(ndist,Xfine,root_dir//out_dir//"Xfine.txt")

Call wri3file(ndist,ny,nz,Jfine,root_dir//out_dir//"Jfine.txt")
!Call wri3file(nx,ny,nz,w,root_dir//out_dir//"w.txt")
if (nz==2 .and. ny==2) then
  open(105,file=root_dir//out_dir//'empfuns_fine.csv',status='replace')
    write(105,40)
    40 format('V,','J(y1;z1),','J(y1;z2),','J(y2;z1),','J(y2;z2),','w(y1;z1),','w(y1;z2),','w(y2;z1),','w(y2;z2),','dprime(y1),','dprime(y2),', &
      'P(y1),','P(y2),','th(y1),','th(y2),','R(y1),','R(y2),','Ptilde(y1),','Ptilde(y2),','M(y1),','M(y2)')
      do i=1,ndist
        write(105,'(<1+ns*2+ny*5>(f15.4,","),2(I6,","))') xfine(i),Jfine(i,1,1),Jfine(i,1,2),Jfine(i,2,1),Jfine(i,2,2),wfine(i,1,1),wfine(i,1,2),&
        wfine(i,2,1),wfine(i,2,2),dprimevecfine(i,1,1),dprimevecfine(i,2,1),P(i,1),P(i,2),theta(i,1),theta(i,2),Rfine(i,1),Rfine(i,2),&
        Ptildefine(i,1),Ptildefine(i,2),Mfine(i,1),Mfine(i,2)
      enddo
  close(105)
elseif (nz==1 .and. ny==2) then
  open(105,file=root_dir//out_dir//'empfuns_fine.csv',status='replace')
    write(105,50)
    50 format('V,','J(y1),','J(y2),','w(y1),','w(y2),','dprime(y1),','dprime(y2),','P(y1),','P(y2),','th(y1),','th(y1),',&
      'R(y1),','R(y1),','Ptilde(y),','Ptilde(y2),','M(y1),','M(y2),')
      do i=1,ndist
        write(105,'(<15>(f15.4,","),2(I6,","))') xfine(i),Jfine(i,1,1),Jfine(i,2,1),wfine(i,1,1),wfine(i,2,1),dprimevecfine(i,1,1),dprimevecfine(i,2,1),&
        Pfine(i,1),Pfine(i,2),thetafine(i,1),thetafine(i,2),Rfine(i,1),Rfine(i,2),Ptildefine(i,1),Ptildefine(i,2),&
        Mfine(i,1),Mfine(i,2)
      enddo
  close(105)

elseif (nz==2 .and. ny==1) then
  open(105,file=root_dir//out_dir//'empfuns_fine.csv',status='replace')
    write(105,60)
    60 format('V,','J(y1;z1),','J(y1;z2),','w(y1;z1),','w(y1;z2),','dprime(y1),', &
      'P(y1),','th(y1),','R(y1),','Ptilde(y1),','M(y1),')
      do i=1,ndist
        write(105,'(<10>(f15.4,","),1(I6,","))') xfine(i),Jfine(i,1,1),Jfine(i,1,2),wfine(i,1,1),wfine(i,1,2),dprimevecfine(i,1,1),&
        Pfine(i,1),thetafine(i,1),Rfine(i,1),Ptildefine(i,1),Mfine(i,1)
      enddo
  close(105)

elseif (nz==3 .and. ny==1) then
  open(105,file=root_dir//out_dir//'empfuns_fine.csv',status='replace')
    write(105,70)
    70 format('V,','J(y1;z1),','J(y1;z2),','J(y1;z3),','w(y1;z1),','w(y1;z2),','w(y1;z3),','dprime(y1),', &
      'P(y1),','th(y1),','R(y1),','Ptilde(y1),','M(y1),')
    do i=1,ndist
      write(105,'(<12>(f15.4,","),1(I6,","))') xfine(i),Jfine(i,1,1),Jfine(i,1,2),Jfine(i,1,3),wfine(i,1,1),wfine(i,1,2),wfine(i,1,3),dprimevecfine(i,1,1),&
      Pfine(i,1),thetafine(i,1),Rfine(i,1),Ptildefine(i,1),Mfine(i,1)
    enddo
  close(105)
elseif (nz==3 .and. ny==2) then
  open(105,file=root_dir//out_dir//'empfuns_fine.csv',status='replace')
    write(105,75)
    75 format('V,','J(y1;z1),','J(y1;z2),','J(y1;z3),','J(y2;z1),','J(y2;z2),','J(y2;z3),','w(y1;z1),','w(y1;z2),','w(y1;z3),',&
      'w(y2;z1),','w(y2;z2),','w(y2;z3),','dprime(y1),','dprime(y2),', 'P(y1),','P(y2),','th(y1),','th(y2),','R(y1),','R(y2),',&
      'Ptilde(y1),','Ptilde(y2),','M(y1),','M(y2),')
    do i=1,ndist
      write(105,'(<23>(f15.4,","),2(I6,","))') xfine(i),Jfine(i,1,1),Jfine(i,1,2),Jfine(i,1,3),Jfine(i,2,1),Jfine(i,2,2),Jfine(i,2,3),&
      wfine(i,1,1),wfine(i,1,2),wfine(i,1,3),wfine(i,2,1),wfine(i,2,2),wfine(i,2,3),dprimevecfine(i,1,1),dprimevecfine(i,2,1), &
      Pfine(i,1),Pfine(i,2),thetafine(i,1),thetafine(i,2),Rfine(i,1),Rfine(i,2),Ptildefine(i,1),Ptildefine(i,2),Mfine(i,1),Mfine(i,2)
    enddo
  close(105)
else
  open(105,file=root_dir//out_dir//'empfuns_fine.csv',status='replace')
    write(105,80)
    80 format('V,','J,','w,','dprime,', &
    'P,','th,','R,','Ptilde(y1),','M,')
    do i=1,ndist
      write(105,'(<8>(f15.4,","),1(I6,","))') xfine(i),Jfine(i,1,1),wfine(i,1,1),dprimevecfine(i,1,1),&
      Pfine(i,1),thetafine(i,1),Rfine(i,1),Ptildefine(i,1),Mfine(i,1)
    enddo
  close(105)
endif

if (nz==3 .and. ny==1) then
  open(107,file=root_dir//out_dir//'vprime_fine.csv',status='replace')
    write(107,'((A6,","),<ns*ns>(I6,","))') 'is',1,1,1,2,2,2,3,3,3
    write(107,'((A6,","),<ns*ns>(I6,","))') 'isp',1,2,3,1,2,3,1,2,3
    write(107,'((A6,","),<ns*ns>(I6,","))') 'y(isp)', iyfun(1),iyfun(2),iyfun(3),iyfun(1),iyfun(2),&
    iyfun(3),iyfun(1),iyfun(2),iyfun(3)
    write(107,'((A6,","),<ns*ns>(I6,","))') 'z(isp)', izfun(1), izfun(2),izfun(3),izfun(1),izfun(2),&
    izfun(3),izfun(1),izfun(2),izfun(3)
    do i=1,ndist
      write(107,'((f15.4,","),<ns*ns>(I6,","))') xfine(i),iVprimefine(i,1,1),iVprimefine(i,1,2),iVprimefine(i,1,3),iVprimefine(i,2,1),iVprimefine(i,2,2),&
      iVprimefine(i,2,3),iVprimefine(i,3,1),iVprimefine(i,3,2),iVprimefine(i,3,3)
    enddo
  close(107)
else if (nz==2 .and. ny==1) then
  open(107,file=root_dir//out_dir//'vprime_fine.csv',status='replace')
    write(107,'((A6,","),<2*ns>(I6,","))') 'is',1,1,2,2
    write(107,'((A6,","),<2*ns>(I6,","))') 'isp',1,2,1,2
    write(107,'((A6,","),<2*ns>(I6,","))') 'y(isp)', iyfun(1), iyfun(2),iyfun(1), iyfun(2)
    write(107,'((A6,","),<2*ns>(I6,","))') 'z(isp)', izfun(1), izfun(2), izfun(1), izfun(2)
    do i=1,ndist
      write(107,'((f15.4,","),<2*ns>(I6,","))') xfine(i),iVprimefine(i,1,1),iVprimefine(i,1,2),iVprimefine(i,2,1),iVprimefine(i,2,2)
    enddo
  close(107)
else if (nz==1 .and. ny==2) then
  open(107,file=root_dir//out_dir//'vprime_fine.csv',status='replace')
    write(107,'((A6,","),<2*ns>(I6,","))') 'is',1,1,2,2
    write(107,'((A6,","),<2*ns>(I6,","))') 'isp',1,2,1,2
    write(107,'((A6,","),<2*ns>(I6,","))') 'y(isp)', iyfun(1), iyfun(2),iyfun(1), iyfun(2)
    write(107,'((A6,","),<2*ns>(I6,","))') 'z(isp)', izfun(1), izfun(2), izfun(1), izfun(2)
    do i=1,ndist
      write(107,'((f15.4,","),<2*ns>(I6,","))') xfine(i),iVprimefine(i,1,1),iVprimefine(i,1,2),iVprimefine(i,2,1),iVprimefine(i,2,2)
    enddo
  close(107)

else if(nz==1 .and. ny==1) then
  open(107,file=root_dir//out_dir//'vprime_fine.csv',status='replace')
    write(107,'((A6,","),<2*ns>(I6,","))') 'is',1
    write(107,'((A6,","),<2*ns>(I6,","))') 'isp',1
    write(107,'((A6,","),<2*ns>(I6,","))') 'y(isp)', iyfun(1)
    write(107,'((A6,","),<2*ns>(I6,","))') 'z(isp)', izfun(1)
    do i=1,ndist
      write(107,'((f15.4,","),<2*ns>(I6,","))') xfine(i),iVprimefine(i,1,1)
    enddo
  close(107)
endif

if (ne>2) then
  open(109,file=root_dir//out_dir//'ufuns_fine.csv',status='replace')
    write(109,90)
    90 format('Past Wage,','U(y1),','U(y2),','wtilde(y1;z1),','wtilde(y1;z2),','wtilde(y2;z1),','wtilde(y2;z2),', &
'PUtilde(y1),','PUtilde(y2),','RU(y1),','RU(y2),','MU(y1),','MU(y2)')
    do i=1,ne
      write(109,'(<1+ny*nz+ny*3>(f15.4,","),2(I6,","))') e(i),U1(1,i),U1(2,i),wfine(MUfine(1,i),1,1),wfine(MUfine(1,i),1,2),wfine(MUfine(2,i),2,1),&
      wfine(MUfine(2,i),2,2),PUtildefine(1,i),PUtildefine(2,i),RUfine(1,i),RUfine(2,i),MUfine(1,i),MUfine(2,i)
    enddo
  close(109)
endif

!Stationary Distribution
Call wri1file(ns*nx+nu,muss,root_dir//out_dir//"muss.txt")
!Collapsed to Rows of length nx+1: U is 1st element of muvec
muvec = zero
do iu=1,nu
  muvec(iu) = muss(ns*ndist+iu)
enddo

do is=1,ns
  do ix=1,ndist
  muvec(nu+ix) = muvec(nu+ix) + muss((is-1)*nx+ix)
  enddo
enddo

open(106,file=root_dir//out_dir//'dist.csv',status='replace')
  write(106,30)
  30 format('V,','dist')
  do iu=1,nu
    write(106,'(2(f15.4,","))') U1(iuyfun(iu),iuefun(iu)),muvec(iu)
  enddo
  do ix=1,ndist
    write(106,'(2(f15.4,","))') xfine(ix),muvec(ix+nu)
  enddo
close(106)

!Print details of the parameters and environment to detail:
write (detail,20) nx,ny,nz,ne,y,pyss,z,pzss,hp,bval,psi,tau,betta,kappa,delta,lambda,sigma,gamma,em,unemp,&
UEflow,EEflow,EUflow,tot_output,avg_wage,avg_benefit,transfers,rrval,welfare
20 format ('   PARAMETERS'/&
'----------------'/&
'nx       = ',I6/&
'ny       = ',I6/&
'nz       = ',I6/&
'ne       = ',I6/&
'y        = ',(<ny>(f12.8,","))/&
'y sdist. = ',(<ny>(f12.8,","))/&
'z        = ',(<nz>(f12.8,","))/&
'z sdist. = ',(<nz>(f12.8,","))/&
'hp       = ',f12.8/&
'b        = ',f12.8/&
'psi      = ',f12.8/&
'tau      = ',f12.8/&
'betta    = ',f12.8/&
'kappa    = ',f12.8/&
'delta    = ',f12.8/&
'lambda   = ',f12.8/&
'sigma    = ',f12.8/&
'gamma    = ',f12.8/&
'   MOMENTS'/&
'----------------'/&

'em           = ',f12.8/&
'unemp        = ',f12.8/&
'UEflow       = ',f12.8/&
'EEflow       = ',f12.8/&
'EUflow       = ',f12.8/&
'tot_output   = ',f12.8/&
'avg_wage     = ',f12.8/&
'avg_benefit  = ',f12.8/&
'transfers    = ',f12.8/&
'rr           = ',f12.8/&
'welfare      = ',f12.8///)


END SUBROUTINE write_output_fine