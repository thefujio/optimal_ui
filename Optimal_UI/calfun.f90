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

    SUBROUTINE write_output(J1,U1)
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

    REAL(8) FUNCTION GBD(tax,uivec)
      USE IOOP
      USE PARAM
      implicit none
      !Dummy arguments declarations
      real(8):: tax
      real(8), dimension(*)   :: uivec
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

  integer, intent(IN):: n
  real*8, intent(IN), dimension(n):: xvec
  real*8, intent(OUT):: fval
  real(8), dimension(ny,ne)     :: U
  real(8), dimension(nx,ny,nz)  :: J

  integer:: iter                      !Generic indexes
  real(8):: fl,fu
  real(8):: bd,jfploss,seploss,j_jloss

  ! TAX CODE - upper and lower bound tax rates for bisection
  real(8):: taul,tauu

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
  taul = 0.0000001d0
  tauu = 0.065d0
  !Evaluate at endpoints taul,tauu
  tau = taul
  call vfi(J,U)
  call sdi(J,U)
  fl = gbd(taul,bvec)
  write (*,'(5x,''Budget Deficit = '',f10.6)') fl

  tau = tauu
  call vfi(J,U)
  call sdi(J,U)
  fu = gbd(tauu,bvec)
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
    bd = gbd(tau,bvec)
    if (bd>0) then
      tauu=tau
    else
      taul=tau
    end if
    if (dabs(taul-tauu) < tol) EXIT !removed .or. dabs(bd) < tol
  end do
  if (iter.ge.niter) then
    write (*,'(3x,''Bisection did not converge after '',i6,'' iterations '')') iter
    write (*,'(5x,''Budget Deficit = '',f6.5)') bd
  else
    write (*,'(3x,''Bisection converged after '',i6,'' iterations '')') iter
    write (*,'(5x,''tau = '',f6.5)') tau
    write (*,'(5x,''Budget Deficit = '',f10.8)') bd
  end if

  call write_output(J,U)
  deallocate(cont)
  !Calibration Output
  jfploss =  min(1.0d8, abs(100.0d0*(UEflow - jfptarget)/jfptarget))
  seploss = min(1.0d8, abs(100.0d0*(EUflow - septarget)/septarget))
  j_jloss = min(1.0d8, abs(100.0d0*(EEflow - j_jtarget)/j_jtarget))

  funcerror = (jfploss + seploss + j_jloss)
  if (isnan(real(funcerror,8))) then
    funcerror = 1.00d8
  endif
  fval = funcerror
  print*, 'kappa = :', kappa, 'delta = ', delta, 'lambda = ', lambda
  print*, 'distance from calibration targets: ', fval

  write (calibout,14) kappa, delta, lambda, jfploss, seploss, j_jloss, funcerror
  14 format ('   Loss Function'/&
  '----------------'/&
  "kappa            ",f18.8/&
  "delta            ",f18.8/&
  "lambda           ",f18.8/&
  "jfploss          ",f18.8/&
  "seploss          ",f18.8/&
  "j2jloss          ",f18.8/&
  "funcerror        ",f18.8///)

  END SUBROUTINE calfun
