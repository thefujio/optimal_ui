PROGRAM MAIN
  !****************************************************************************
  !  MAIN.f90 - Entry point of console application.
  !  
  !  PURPOSE:
  !   Finds a BR allocation
  !   
  !  
  !  MODULES:
  !   PARAM      - Contains all the parameters of the model
  !   UTILITY    - Contains all the functions related to the utility function
  !   IOOP       - I/O routines
  !   INTERFACES - Contains the interface statement of some functions/subroutines
  !   Many others...
  !   
  !   
  !  Program written by: M. Gervais
  !   Date:    Jan 2016
  !   Updated: Jan 2016
  !****************************************************************************
  USE PARAM
  USE UTILITY
  USE IOOP
!  USE omp_lib
!  USE isnan_int
  implicit none
  
  INTERFACE
    SUBROUTINE VFI(J1,U1)
      USE IOOP
      USE PARAM
      USE UTILITY
      implicit none
      !Dummy arguments declarations
      double precision, dimension(ny,ne)      , intent(inout):: U1
      double precision, dimension(nx,ny,nz), intent(inout):: J1
      END SUBROUTINE VFI

    SUBROUTINE write_output(U1,J1)
    !  USE INTERFACES
      USE IOOP
      USE PARAM
      implicit none
      !Dummy arguments declarations
      double precision, dimension(ny,ne)      , intent(in):: U1
      double precision, dimension(nx,ny,nz), intent(in):: J1
    END SUBROUTINE write_output


    SUBROUTINE SDI(J1,U1)
      USE IOOP
      USE PARAM
      implicit none
      double precision, dimension(ny,ne), intent(in)   :: U1
      double precision, dimension(nx,ny,nz), intent(in):: J1

    END SUBROUTINE SDI
    
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
    END FUNCTION GBD

    SUBROUTINE stadist(m,pimat,pss,initpss)
    implicit none
    !dummy arguments declarations
    integer, intent(in)                   :: m
    double precision, intent(in)          :: pimat(m,m)
    double precision, intent(out)         :: pss(m)
    double precision, intent(in), dimension(m), optional:: initpss
    END SUBROUTINE stadist
  END INTERFACE


  !Variables declarations
  real(8):: time_begin,time_end            !Used to keep track of time
!  real(8):: start_iter, end_iter           !Used to keep track of time
  
  real(8), dimension(ny,ne)      :: U
  real(8), dimension(nx,ny,nz):: J
  
  integer:: iter,i,jj                      !Generic indexes
  integer:: is,ix,iy,iz,ie
  integer:: isp,ixp,iyp,izp,iep
  real(8):: fl,fu
  real(8):: bd
!  real(8):: norm
  
  ! TAX CODE
  real(8):: taul,tauu
  
  
  !Timing of the program
  call CPU_Time(time_begin)
  ! CHANGING IMSL HANDLING OF ERRORS
!  CALL erset(4, 1, 0) !this prevents BCONF from stopping when the max # of iter is reached
!  CALL erset(0, 0, 0) !this prevents BCONF from stopping OR printing when the max # of iter is reached
  
  
  !Open the output file where details are printed
  open(unit=detail,  file=root_dir//out_dir//"detail.txt",  status='replace')
  
  
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
  if(ne==1) then
  pe = one
  else if (ne==2) then
  pe(1,1) = one - psi
  pe(1,2) = psi
  pe(2,2) = one
  pe(2,1) = zero
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
  !Grid on PVU (x)
  xmin = Ufunc(hp+bmin)/(one-betta)
  xmax = Ufunc(MAXVAL(y)+MAXVAL(z))/(one-betta)
  print*,'maximum output is: ',MAXVAL(y)+MAXVAL(z)

  x(1:nx) = (/ ( &
                ((xmax-xmin)/(real(nx-1,8)))*(real(i-1,8)) + xmin, i=1,nx &
                ) /)
  call wri2file(nx,1,x,root_dir//out_dir//"xgrid.txt")
  
  !Setting initial conditions for value functions
  U = xmin
  do iz=1,nz
    do iy=1,ny
      do ix=1,nx
        J(ix,iy,iz)=(one/(one-betta))*&
          (one-(((one-sigma)*(one-betta)*x(ix)+one)**(one/(one-sigma))))
      end do
    end do
  end do
  
!  call readfile(nx,1,x,root_dir//out_dir//"x.txt")
!  call readfile(nx,1,J,root_dir//out_dir//"jfunc.txt")
   U = -18.929393939237876d0
  !Set unemployment benefit: e=1 is eligible for UI
   b = bmin
  
  !Bisection on tax rate
  taul = 0.0001d0
!  taul = 0.037890682220459d0
  tauu = 0.05d0
  tau = taul
  call vfi(J,U)
  call sdi(J,U)
  fl = gbd(taul,b)
  write (*,'(5x,''Budget Deficit = '',f10.6)') fl

  tau = tauu
  call vfi(J,U)
  call sdi(J,U)
  fu = gbd(tauu,b)
  write (*,'(5x,''Budget Deficit = '',f10.6)') fu
  if (fl*fu>zero) then
    write (*,'(3x,''Stop: Bisection not bracketed'')')
    STOP
  end if

  !Bisection loop
  do iter=1,niter
    tau = half*(taul+tauu)
    call vfi(J,U)
    call sdi(J,U)
    bd = gbd(tau,b)
    if (bd>0) then
      tauu=tau
    else
      taul=tau
    end if
    if (dabs(taul-tauu) < tol .or. dabs(bd) < tol) EXIT
  end do
  if (iter.ge.niter) then
    write (*,'(3x,''Bisection did not converge after '',i6,'' iterations '')') iter
    write (*,'(5x,''Budget Deficit = '',f6.5)') bd
  else
    write (*,'(3x,''Bisection converged after '',i6,'' iterations '')') iter
    write (*,'(5x,''tau = '',f6.5)') tau
    write (*,'(5x,''Budget Deficit = '',f10.8)') bd
  end if

  params(1) = betta
  params(2) = eta
  params(3) = sigma
  params(4) = gamma
  params(5) = lambda
  params(6) = delta
  params(7) = kappa
  params(8) = nx
  params(9) = b
  params(10)= hp
  call wri2file(nparams,1,params,root_dir//out_dir//"params.txt")
  call write_output(U,J)
  call CPU_Time(time_end)
  !Timing of the program
  write(*,"('TOTAL EXECUTION TIME (seconds): ',f8.2)") &
  (time_end-time_begin)/8.0d0
  
  STOP
END PROGRAM MAIN
