SUBROUTINE VFI(J1,U1)
  !****************************************************************************
  !  VFI.f90
  !  
  !  PURPOSE:
  !   Finds the value functions and policy rules given a tax code (tau,b)
  !     
  !  
  !  Program written by: M. Gervais
  !   Date:    Feb 2016
  !   Updated: Feb 2016
  !****************************************************************************
!  USE INTERFACES
  USE IOOP
  USE PARAM
  USE UTILITY
  implicit none
  
  !Dummy arguments declarations
  double precision, dimension(ny,ne)      , intent(inout):: U1
  double precision, dimension(nx,ny,nz), intent(inout):: J1

  INTERFACE
    SUBROUTINE JPRIME(U1,J1)
      USE IOOP
      USE PARAM
      USE UTILITY
      implicit none
      !Dummy arguments declarations
      double precision, dimension(ny,ne)      , intent(inout):: U1
      double precision, dimension(nx,ny,nz), intent(inout):: J1
    END SUBROUTINE JPRIME

    SUBROUTINE centergrid(xgrid,xl,xu)
      USE PARAM
      USE UTILITY
      implicit none
      !Dummy arguments declarations
      double precision, intent(inout):: xl
      double precision, intent(inout):: xu
      double precision, dimension(nx), intent(out)::xgrid
    END SUBROUTINE centergrid

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

  !Local variables declarations
  integer:: iter,iz,iy,ix
  real(8):: norm, normx
  double precision, dimension(ny,ne)   :: U0
  double precision, dimension(nx,ny,nz):: J0
  double precision, dimension(nx)      :: xnew,Jvec,splineb,splinec,splined


  do iter=1,niter
    U0 = U1
    J0 = J1
    if (ne>2) then
      call jprime_wage(U1,J1)
    else
      call jprime(U1,J1)
    endif
    norm = MAXVAL(dabs(J1-J0))
    if (want_print) then
    write(*,'(5x,''norm at iteration '',i3,'' = '',f15.8)') iter,norm
    write(*,'(5x,''xmin: '',f11.6,'' xmax: '',f11.6)') xmin, xmax
    endif
    if (iter>2 .AND. want_howard) do_howard   = .true.

    !adjust space of grid of submarkets
    !print*, dprime(:,1,1)
    !print*, theta
    if (want_print) then
      print*, 'U: ', U1
    endif

    if (norm > 1.0d-2) then
    call centergrid(xnew,xmin,xmax)
    !interpolate J on the new submarket grid
    do iz=1,nz
      do iy=1,ny
      splineb = 0.0d0
      splinec = 0.0d0
      splined = 0.0d0
      Jvec = J1(:,iy,iz)
      call spline(x,Jvec,splineb,splinec,splined,nx)
        do ix=1,nx
        J1(ix,iy,iz) = ispline(xnew(ix), x, Jvec, splineb, splinec, splined, nx)
        enddo
      !print*, 'Jdiff', iy, iz,': ', J1(:,iy,iz)-Jvec
      enddo
    enddo
    normx = MAXVAL(dabs(xnew-x))
    if (want_print) then
      print*, 'norm of Xnew-X:', normx
    endif
    x = xnew
    endif


    if (norm<tol) EXIT
  end do
  write (*,'(3x,''Value function converged after '',i3,'' iterations: norm = '',f11.8)') iter,norm
  
  do_howard = .false.
 	
  RETURN
END SUBROUTINE VFI
