SUBROUTINE refine(J1,U1)
!
!   refine.f90
!   Optimal_UI
!
!   Created by Larry Warren on 10/31/17.
!   Copyright 2017 Larry Warren. All rights reserved.
!
USE PARAM
USE UTILITY
USE TOOLBOX
implicit none
!Dummy Variable Declarations
double precision, dimension(ny,ne), intent(inout):: U1
double precision, dimension(nx,ny,nz), intent(inout):: J1

INTERFACE
  SUBROUTINE movavg(invec,outvec,pd)
    real*8, dimension(:), intent(in) :: invec
    integer, intent(in) ::pd
    real*8, dimension(:), intent(out) :: outvec
  END SUBROUTINE movavg

  SUBROUTINE spline(xvec,yvec,b,c,d,n)
    integer:: n
    double precision, dimension(n):: xvec, yvec
    double precision, dimension(n):: b, c, d
  END SUBROUTINE spline

  double precision function ispline(u, xvec, yvec, b, c, d, n)
    implicit none
    integer:: n
    double precision:: u
    double precision, dimension(n):: xvec, yvec, b, c, d
  end function ispline

END INTERFACE

real(8), dimension(nx,ny,nz) :: wma
real(8), dimension(nx,ny,nz) :: windreal,windma
real(8), dimension(nx,ns,ns) :: dprimema
real(8), dimension(nx,ns,ns) :: Vprime, Vprimema
real(8), dimension(nx,ny)    :: Ptildema
real(8), dimension(nx,ny)    :: Mma
!real(8), dimension(ny,ne)    :: RUfine,PUtildefine
real(8), dimension(nx)       :: wvec, wvectemp, Vprimevec, dvec, Mvec, Ptildevec, Rvec, thetavec, Pvec, Jtvec
real(8), dimension(nx)       :: splineb,splinec,splined
real(8), dimension(ndist)    :: tempM, tempiV, tempiw, retfine
real(8), dimension(nu)       :: ExpUfine
real(8), dimension(ndist,ny) :: Mfinereal
integer, dimension(ny,ne)    :: NDind
integer, dimension(ny)       :: NTind

integer :: i,ii,j,jj,ix,ixx,iy,iz,ie,is,isp, smp
!Create fine grid on X
call linspace(xfine,xmin,xmax,ndist)

!smoothing parameter for ma
smp = 19

!To interpolate:
!M(ix,iy)
do iy=1,ny
  splineb = 0.0d0
  splinec = 0.0d0
  splined = 0.0d0
  do ix=1,nx
    Mvec(ix) = x(M(ix,iy))
  end do
  call spline(x,Mvec,splineb,splinec,splined,nx)
  do ix=1,ndist
  Mfinereal(ix,iy) = ispline(xfine(ix), x, Mvec, splineb, splinec, splined, nx)
  enddo
enddo

do iy=1,ny
  do ix=1,ndist
    tempM = ABS(Mfinereal(ix,iy)-xfine)
    Mfine(ix,iy) = MINLOC(tempM, DIM=1)
  end do
end do

!Ptilde(ix,iy)
do iy=1,ny
  splineb = 0.0d0
  splinec = 0.0d0
  splined = 0.0d0
  Ptildevec = Ptilde(:,iy)
  call spline(x,Ptildevec,splineb,splinec,splined,nx)
  do ix=1,ndist
    Ptildefine(ix,iy) = ispline(xfine(ix), x, Ptildevec, splineb, splinec, splined, nx)
  enddo
enddo

!R(ix,iy)
do iy=1,ny
  splineb = 0.0d0
  splinec = 0.0d0
  splined = 0.0d0
  Rvec = R(:,iy)
  call spline(x,Rvec,splineb,splinec,splined,nx)
  do ix=1,ndist
    Rfine(ix,iy) = ispline(xfine(ix), x, Rvec, splineb, splinec, splined, nx)
  enddo
enddo

!theta(ix,iy)
do iy=1,ny
  splineb = 0.0d0
  splinec = 0.0d0
  splined = 0.0d0
  thetavec = theta(:,iy)
  call spline(x,thetavec,splineb,splinec,splined,nx)
  do ix=1,ndist
    thetafine(ix,iy) = ispline(xfine(ix), x, thetavec, splineb, splinec, splined, nx)
  enddo
enddo

!P(ix,iy)
do iy=1,ny
  splineb = 0.0d0
  splinec = 0.0d0
  splined = 0.0d0
  Pvec = P(:,iy)
  call spline(x,Pvec,splineb,splinec,splined,nx)
  do ix=1,ndist
    Pfine(ix,iy) = ispline(xfine(ix), x, Pvec, splineb, splinec, splined, nx)
  enddo
enddo


do ie=1,ne
  do iy=1,ny
    do ix=1,ndist
      retfine(ix) = Pfine(ix,iy)*(xfine(ix)-U1(iy,ie))
    end do
    RUfine(iy,ie) = MAXVAL(retfine)
    MUfine(iy,ie) = MAXLOC(retfine,DIM=1)
    PUtildefine(iy,ie) = Pfine(MUfine(iy,ie),iy)
  end do
end do

!Recalculate U1 on xfine grid
  do i=1,nu
    ExpUfine(i) = zero
    do ii=1,nu
      ExpUfine(i) = ExpUfine(i) + betta*pus(i,ii)*(U1(iuyfun(ii),iuefun(ii))+RUfine(iuyfun(ii),iuefun(ii)))
    end do
    Ufine(iuyfun(i),iuefun(i)) = Ufunc(bvec(iuefun(i))) + ExpUfine(i)
  end do

!dprimevec(ix,iy,ie) is either delta or 1 - so don't spline, rather set to nearest gridpt value
do iy=1,ny
  do ie=1,ne
    dvec = dprimevec(:,iy,ie)
    call nearest_interp (nx,x,dvec,ndist,xfine,dprimevecfine(:,iy,ie))
  enddo
enddo

!w and V'

!only smooth where d=delta and theta postitive?
do ie = 1,ne
  do iy = 1,ny
    NDind(iy,ie)=count(dprimevec(:,iy,ie) == 1.0d0)
  end do
end do

do iy = 1,ny
  NTind(iy)=count(theta(:,ny) <= 0.0d0)
end do

do iz=1,nz
  do iy=1,ny
    wma(:,iy,iz) = w(:,iy,iz)
    call movavg(w(NDind(iy,1):nx-NTind(iy),iy,iz),wma(NDind(iy,1):nx-NTind(iy),iy,iz),smp)
  end do
end do

do jj=1,ns
  do ii=1,ns
    do ix=1,nx
      Vprime(ix,ii,jj) = x(iVprime(ix,ii,jj))
    end do
    Vprimema(:,ii,jj) = Vprime(:,ii,jj)
    iy=iyfun(ii)
    call movavg(Vprime(NDind(iy,1):nx-NTind(iy),ii,jj),Vprimema(NDind(iy,1):nx-NTind(iy),ii,jj),smp)
  end do
end do

!Spline interpolate on fine grid:
!spline (xvec, yvec, b, c, d, n)
do iz=1,nz
  do iy=1,ny
  splineb = 0.0d0
  splinec = 0.0d0
  splined = 0.0d0
  wvec = wma(:,iy,iz)
  call spline(x,wvec,splineb,splinec,splined,nx)
    do ix=1,ndist
    wfine(ix,iy,iz) = ispline(xfine(ix), x, wvec, splineb, splinec, splined, nx)
    enddo
  enddo
enddo


do is=1,ns
  do isp=1,ns
  splineb = 0.0d0
  splinec = 0.0d0
  splined = 0.0d0
  Vprimevec = Vprimema(:,is,isp)
  call spline(x,Vprimevec,splineb,splinec,splined,nx)
    do ix=1,ndist
    Vprimefine(ix,is,isp) = ispline(xfine(ix), x, Vprimevec, splineb, splinec, splined, nx)
    enddo
  enddo
enddo

do is=1,ns
  do isp=1,ns
    do ix=1,ndist
      tempiV = ABS(Vprimefine(ix,is,isp)-xfine)
      iVprimefine(ix,is,isp) = minloc(tempiV, DIM=1)
    end do
  end do
end do

!Need to interpolate wind, which is sometimes continuous and sometimes binary. Best option is just nearest-neighbor.
do iz=1,nz
  do iy=1,ny
    wvectemp = real(wind(:,iy,iz),8)
    call nearest_interp (nx,x,wvectemp,ndist,xfine,tempiw)
    windfine(:,iy,iz) = NINT(tempiw)
  end do
end do

do is=1,ns
  do isp=1,ns
    do ix=1,ndist
      iy = iyfun(is)
      iz = izfun(is)
      dprimefine(ix,is,isp) = dprimevecfine(iVprimefine(ix,is,isp),iyfun(isp),windfine(ix,iy,iz))
    end do
  end do
end do

!J(ix,iy,iz)
Jfine = zero
do iz=1,nz
  do iy=1,ny
  splineb = 0.0d0
  splinec = 0.0d0
  splined = 0.0d0
  Jtvec = J1(:,iy,iz)
  call spline(x,Jtvec,splineb,splinec,splined,nx)
    do ix=1,ndist
    Jfine(ix,iy,iz) = ispline(xfine(ix), x, Jtvec, splineb, splinec, splined, nx)
    enddo
  enddo
enddo

END SUBROUTINE refine
