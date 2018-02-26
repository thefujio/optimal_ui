  SUBROUTINE tau_eval(n,xvec,tax,fval)
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
  real*8, intent(IN):: tax
  real*8, intent(OUT):: fval
  real(8), dimension(ny,ne)     :: U
  real(8), dimension(nx,ny,nz)  :: J

  integer:: iter                      !Generic indexes
  real(8):: f
  real(8):: bd,jfploss,urloss,j_jloss

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
  !Set tax rate
  tau = tax
  call vfi(J,U)
  call refine(J,U)
  call sdifine(J,U)
  f = gbdfine(tau,bvec)
  write (*,'(5x,''Tax Rate = '',f10.6)') tau
  write (*,'(5x,''Budget Deficit = '',f10.6)') f
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

  END SUBROUTINE tau_eval
