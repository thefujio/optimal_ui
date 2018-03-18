  SUBROUTINE calfun(n,xvec,fval)
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
  USE TOOLBOX
  USE POWELL

  !  USE omp_lib
  !  USE isnan_int
  implicit none

  INTERFACE
    SUBROUTINE STATES(J1,U1)
      USE PARAM
      USE UTILITY
      USE IOOP
      implicit none
      !Dummy arguments declarations
      double precision, dimension(ny,ne)      , intent(inout):: U1
      double precision, dimension(nx,ny,nz), intent(inout):: J1
    END SUBROUTINE STATES

    SUBROUTINE VFI(J1,U1)
      USE IOOP
      USE PARAM
      USE UTILITY
      implicit none
      !Dummy arguments declarations
      double precision, dimension(ny,ne)      , intent(inout):: U1
      double precision, dimension(nx,ny,nz), intent(inout):: J1
    END SUBROUTINE VFI

    SUBROUTINE REFINE(J1,U1)
      USE PARAM
      USE UTILITY
      USE TOOLBOX
      implicit none
      !Dummy Variable Declarations
      double precision, dimension(ny,ne), intent(inout):: U1
      double precision, dimension(nx,ny,nz), intent(inout):: J1
    END SUBROUTINE REFINE

    SUBROUTINE write_output_fine(J1,U1)
      USE IOOP
      USE PARAM
      implicit none
      !Dummy arguments declarations
      double precision, dimension(ny,ne)      , intent(in):: U1
      double precision, dimension(nx,ny,nz), intent(in):: J1   
    END SUBROUTINE write_output_fine

    SUBROUTINE SDIFINE(J1,U1)
      USE IOOP
      USE PARAM
      implicit none
      double precision, dimension(ny,ne), intent(in)   :: U1
      double precision, dimension(nx,ny,nz), intent(in):: J1
    END SUBROUTINE SDIFINE

    REAL(8) FUNCTION GBDFINE(tax,uivec)
      USE IOOP
      USE PARAM
      implicit none
      !Dummy arguments declarations
      real(8):: tax
      real(8), dimension(*)   :: uivec
    END FUNCTION GBDFINE

    SUBROUTINE stadist(m,pimat,pss,initpss)
      implicit none
      !dummy arguments declarations
      integer, intent(in)                   :: m
      double precision, intent(in)          :: pimat(m,m)
      double precision, intent(out)         :: pss(m)
      double precision, intent(in), dimension(m), optional:: initpss
    END SUBROUTINE stadist
  END INTERFACE

  integer, intent(IN):: n
  real*8, intent(IN), dimension(n):: xvec
  real*8, intent(OUT):: fval
  real(8), dimension(ny,ne)     :: U
  real(8), dimension(nx,ny,nz)  :: J
  real(8), dimension(biter)     :: bdvec,tauvec
  integer:: iter,tauind                      !Generic indexes
  real(8):: fl,fu
  real(8):: bd,jfploss,urloss,j_jloss

  ! TAX CODE - upper and lower bound tax rates for bisection
  real(8):: taul,tauu,tau0,tau1,f0,f1

  if (transform == 1) then
    kappa = restrict(xvec(1),lb(1),ub(1),spread)
    delta = restrict(xvec(2),lb(2),ub(2),spread)
    lambda = restrict(xvec(3),lb(3),ub(3),spread)
  else
    kappa = xvec(1)
    delta = xvec(2)
    lambda = xvec(3)
  endif

  !Set up grids and transition matrices
  call states(J,U)
  !Bisection on tax rate
  taul = 0.015d0
  tauu = 0.055d0
  !Evaluate at endpoints taul,tauu
  100 tau = taul
  call vfi(J,U)
  call refine(J,U)
  call sdifine(J,U)
  fl = gbdfine(taul,bvec)
  write (*,'(5x,''Tax lower = '',f10.6)') tau
  write (*,'(5x,''Budget Deficit = '',f10.6)') fl
  tau0 = tau
  f0 = fl

  tau = tauu
  call vfi(J,U)
  call refine(J,U)
  call sdifine(J,U)
  fu = gbdfine(tauu,bvec)
  write (*,'(5x,''Tax upper = '',f10.6)') tau
  write (*,'(5x,''Budget Deficit = '',f10.6)') fu
  tau1 = tau
  f1=fu
  if (fl*fu>zero) then
    write (*,'(3x,''Stop: Root not bracketed'')')
    !PAUSE
  end if

  tauvec=1.0d0
  bdvec=1.0d0
  !Secant loop
  do iter=1,biter
    !tau = half*(taul+tauu)
    tau = tau1 - f1*((tau1-tau0)/(f1-f0))
    call vfi(J,U)
    call refine(J,U)
    call sdifine(J,U)
    bd = gbdfine(tau,bvec)
!    if (bd>0.0d0) then
!      tauu=tau
!    else
!      taul=tau
!    end if

    tau0 = tau1
    f0 = f1
    tau1 = tau
    f1 = bd

    write (*,'(5x,''tau = '',f10.8)') tau
    write (*,'(5x,''Budget Deficit = '',f10.8)') bd
    write (*,'(5x,''iteration = '',i4)') iter

    tauvec(iter)=tau
    bdvec(iter)=bd

    if (dabs(tau1-tau0) < tol .or. dabs(bd) < bis_tol) EXIT

  end do !end of iter loop

  if (iter.ge.biter) then
    write (*,'(3x,''Bisection did not converge after '',i6,'' iterations, choosing closest value '')') iter
      tauind = minloc(abs(bdvec),DIM=1)
      write (*,'(5x,''Recalculating at range of closest tau = '',f10.8)') tauvec(tauind)
      write (*,'(5x,''With Budget Deficit = '',f10.8)') bdvec(tauind)
      tauu = tauvec(tauind) + 0.005d0
      taul = tauvec(tauind) - 0.005d0
      GO TO 100
      !tau = tauvec(tauind)
      !call vfi(J,U)
      !call refine(J,U)
      !call sdifine(J,U)
      !bd = gbdfine(tau,bvec)
      !write (*,'(5x,''tau = '',f10.8)') tau
      !write (*,'(5x,''Budget Deficit = '',f10.8)') bd
  else
    write (*,'(3x,''Bisection converged after '',i6,'' iterations '')') iter
    write (*,'(5x,''tau = '',f10.8)') tau
    write (*,'(5x,''Budget Deficit = '',f10.8)') bd
  end if

  bdval = bd
  call write_output(J,U)
  call write_output_fine(J,U)
  deallocate(cont)

  !Calibration Output
  jfploss =  min(1.0d8, abs(100.0d0*(UEflow - jfptarget)/jfptarget))
  urloss = min(1.0d8, abs(100.0d0*(unemp - urtarget)/urtarget))
  j_jloss = min(1.0d8, abs(100.0d0*(EEflow - j_jtarget)/j_jtarget))

  funcerror = (jfploss**2.0d0 + urloss**2.0d0 + j_jloss**2.0d0)
  if (isnan(real(funcerror,8))) then
    funcerror = 1.00d8
  endif
  fval = funcerror
  print*, 'kappa = :', kappa, 'delta = ', delta, 'lambda = ', lambda
  print*, 'distance from calibration targets: ', fval

  write (calibout,14) kappa, delta, lambda, jfploss, urloss, j_jloss, funcerror
  14 format ('   Loss Function'/&
  '----------------'/&
  "kappa            ",f18.8/&
  "delta            ",f18.8/&
  "lambda           ",f18.8/&
  "jfploss          ",f18.8/&
  "urloss           ",f18.8/&
  "j2jloss          ",f18.8/&
  "funcerror        ",f18.8///)

  END SUBROUTINE calfun
