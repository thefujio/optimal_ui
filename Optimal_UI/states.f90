  SUBROUTINE STATES(J1,U1)
  !****************************************************************************
  !  states.f90 - Entry point of console application.
  !
  !  PURPOSE:
  !   Sets up grid for PVU, bvec, states and transition matrices for y, z, etc.
  !   Originially resided in MAIN.f90
  !
  !  MODULES:
  !   PARAM      - Contains all the parameters of the model
  !   UTILITY    - Contains all the functions related to the utility function
  !   IOOP       - I/O routines
  !   INTERFACES - Contains the interface statement of some functions/subroutines
  !
  !
  !
  !  Program written by: M. Gervais and edited by Larry Warren
  !   Date:    Jan 2016
  !   Updated: June 2017
  !****************************************************************************
  USE PARAM
  USE UTILITY
  USE IOOP
  implicit none

  INTERFACE
    SUBROUTINE stadist(m,pimat,pss,initpss)
      implicit none
      !dummy arguments declarations
      integer, intent(in)                   :: m
      double precision, intent(in)          :: pimat(m,m)
      double precision, intent(out)         :: pss(m)
      double precision, intent(in), dimension(m), optional:: initpss
    END SUBROUTINE stadist
  END INTERFACE

  !Dummy arguments declarations
  double precision, dimension(ny,ne)      , intent(inout):: U1
  double precision, dimension(nx,ny,nz), intent(inout):: J1

  integer:: iter,i,jj                      !Generic indexes
  integer:: is,ix,iy,iz,ie
  integer:: isp,ixp,iyp,izp,iep

  nc = nx**ns
  allocate(cont(nc,ns))
  call contract(cont,nx,ns,nc)

  !Transition matrix and Aggregate TFP
  if(ny==1) then
    y = one
    py = one
  else if (ny==2) then
    py(1,1) = (recession_length-one)/recession_length
    py(1,2) = one - py(1,1)
    py(2,2) = (boom_length-one)/boom_length
    py(2,1) = one - py(2,2)
  call stadist(ny,py,pyss)
    y(1) = 0.955d0
    y(2) = (one-pyss(1)*y(1))/pyss(2) !This makes expected output 1 in the LR
    !y = one
  else
    write (*,'(3x,''This code cannot handle ny>2: Quitting'')')
    STOP
  end if

  !Transition matrix for Idiosyncratic Productivity
  if(nz==1) then
    z = zero
    pz = one
    pzss = one
    Pztilde = one
  else if (nz==2) then
    pz(1,1) = 0.25d0
    pz(1,2) = one - pz(1,1)
    pz(2,2) = 0.25d0
    pz(2,1) = one - pz(2,2)
    call stadist(nz,pz,pzss)
    z(1) =-0.2d0
    z(2) = 0.2d0
    Pztilde = pzss
  else if (nz==3) then
    pz = zero
    pz(1,1) = 0.75d0
    pz(1,3) = one - pz(1,1)
    pz(2,1) = 0.5d0
    pz(2,3) = one - pz(2,1)
    pz(3,1) = 0.25d0
    pz(3,3) = one - pz(3,1)
    call stadist(nz,pz,pzss)
    z(1) =-0.2d0
    z(2) = zero
    z(3) = 0.2d0
    Pztilde = zero
    Pztilde(2) = one
  else if (nz>3) then
    write (*,'(3x,''This code cannot handle nz>3: Quitting'')')
    STOP
  end if
  !Transition matrix and shocks over all states
  i=1
  do iy=1,ny
    do iz=1,nz
      iyfun(i) = iy
      izfun(i) = iz
      jj=1
      do iyp=1,ny
        do izp=1,nz
          ps(i,jj) = py(iy,iyp)*pz(iz,izp)
          jj=jj+1
        end do
      end do
      isfun(iy,iz) = i
      i=i+1
    end do
  end do

  !Transition matrix for UI eligibility
  ! if ne>1, gridpt ne is ineligible.
  if(ne==1) then
    pe = one
  else if (ne==2) then
    pe(1,1) = one - psi
    pe(1,2) = psi
    pe(2,2) = one
    pe(2,1) = zero
  else if (ne>2) then
    do ie=1,ne-1
      pe(ie,ie) = one-psi
      pe(ie,ne) = psi
    end do
    pe(ne,ne) = one
  end if

  !Transition process for unemployed over agg. state and eligibility
  i=1
  do iy=1,ny
    do ie=1,ne
      iuyfun(i) = iy
      iuefun(i) = ie
      jj=1
      do iyp=1,ny
        do iep=1,ne
          pus(i,jj) = py(iy,iyp)*pe(ie,iep)
          jj=jj+1
        end do
      end do
      iufun(iy,ie) = i
      i=i+1
    end do
  end do

  !Grid of unemployment benefits (e) (past wage):
  !replacement rate is roughly 46% in U.S.
  if (ne==1) then
    bvec = hp+bval
  elseif (ne==2) then
    bvec(1) = hp+bval
    bvec(2) = hp
    print*, 'bvec:', bvec
  elseif (ne>2) then
    emin = 0.9d0*MINVAL(y)+MINVAL(z)
    emax = MAXVAL(y)+MAXVAL(z)
    e(1:ne-1) = (/ ( &
    ((emax-emin)/(real(ne-2,8)))*(real(i-1,8)) + emin, i=1,ne-1) /)
    e(ne) = zero
    call wri2file(ne,1,e,root_dir//out_dir//"egrid.txt")

    !Set unemployment benefit: e=1:ne-1 is eligible for UI ne = ineligible (wage=0)
    do ie=1,ne-1
      bvec(ie) = hp+bval*e(ie)
    end do
    bvec(ne) = hp
    print*, 'bvec:', bvec
  endif
  !Grid on PVU (x)
  xmin = Ufunc(hp+bval)/(one-betta)
  xmax = Ufunc(MAXVAL(y)+MAXVAL(z))/(one-betta)
  print*,'maximum output is: ',MAXVAL(y)+MAXVAL(z)

  x(1:nx) = (/ ( &
  ((xmax-xmin)/(real(nx-1,8)))*(real(i-1,8)) + xmin, i=1,nx &
  ) /)
  call wri2file(nx,1,x,root_dir//out_dir//"xgrid.txt")

  !Setting initial conditions for value functions

  U1 = -18.929393939237876d0
  do iz=1,nz
    do iy=1,ny
      do ix=1,nx
        J1(ix,iy,iz)=(one/(one-betta))*&
        (one-(((one-sigma)*(one-betta)*x(ix)+one)**(one/(one-sigma))))
      end do
    end do
  end do
  !If code has been run already, use existing solution to speed up:
    !call readfile(nx,1,x,root_dir//out_dir//"x.txt")
    !call read3file(nx,ny,nz,J1,root_dir//out_dir//"J.txt")
  !
  END SUBROUTINE STATES