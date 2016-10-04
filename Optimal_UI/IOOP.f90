MODULE IOOP
!*************************************************************************
!  ioop.f90
!  
!   
!  PURPOSE:
!   Reading and writing m by n matrices
!
!
!  Program written by: J. Fisher and M. Gervais
!   Date:    May 2005
!   Updated: July 2006
!*************************************************************************
  implicit none
  private
  public:: wri1file, wri2file, int1file, int2file, int2file1, readfile, &
           readintfile, readintfile1, &
           wriscalar, readscalar, &
           wri3file,read3file, &
           wri4file,read4file, &
           int3file,readint3file, &
           wrilogical, readlogical

CONTAINS

    !**********************************************************************
    ! SUBROUTINE wri1file
    !**********************************************************************
    SUBROUTINE wri1file(n,x,filename)
    implicit none

    !dummy arguments declarations
    character(*)                                  :: filename
    integer, intent(in)                           :: n
    double precision, intent(in), dimension (n) :: x

    !Local variables declarations
    character (len=*), parameter  :: newline = '()'
    integer                       :: j


    open (unit=99,file=filename)
    do j=1,n
    write (99,'(f25.15,$)') x(j)
    end do

    !    open (unit=99,file=filename)
    !    do i=1,m
    !      write (99,'(2X,150f20.10,2X)'), (x(i,j), j=1, n)
    !    end do

    close(99)

    END SUBROUTINE wri1file

  !**********************************************************************
  ! SUBROUTINE wri2file
  !**********************************************************************
  SUBROUTINE wri2file(m,n,x,filename)
    implicit none
    
    !dummy arguments declarations
    character(*)                                  :: filename
    integer, intent(in)                           :: m,n
    double precision, intent(in), dimension (m,n) :: x
    
    !Local variables declarations
    character (len=*), parameter  :: newline = '()'
    integer                       :: i,j
    
    
    open (unit=99,file=filename)  
    do i=1,m
      do j=1,n
        write (99,'(f25.15,$)') x(i,j)
      end do
      write (99,newline)
    end do
    
!    open (unit=99,file=filename)  
!    do i=1,m
!      write (99,'(2X,150f20.10,2X)'), (x(i,j), j=1, n)
!    end do
    
    close(99)
    
  END SUBROUTINE wri2file
  
  
  !**********************************************************************
  ! SUBROUTINE wri3file
  !**********************************************************************
  SUBROUTINE wri3file(m,n,o,x,filename)
    implicit none
    
    !dummy arguments declarations
    character(*)                                   :: filename
    integer, intent(in)                            :: m,n,o
    double precision, intent(in), dimension (m,n,o):: x
    
    !Local variables declarations
    integer                       :: i,j,k
    
    open (unit=99,file=filename)
    do i=1,m
      do j=1,n
        do k=1,o
          write (99,'(f25.15)') x(i,j,k)
        end do
      end do
    end do
    
    close(99)
  END SUBROUTINE wri3file
  
  
  !**********************************************************************
  ! SUBROUTINE wri4file
  !**********************************************************************
  SUBROUTINE wri4file(m,n,o,p,x,filename)
    implicit none
    
    !dummy arguments declarations
    character(*)                                      :: filename
    integer, intent(in)                               :: m,n,o,p
    double precision, intent(in), dimension (m,n,o,p) :: x
    
    !Local variables declarations
    character (len=*), parameter  :: newline = '()'
    integer                       :: i,j,k,l
    
    open (unit=99,file=filename)  
    do i=1,m
      do j=1,n
        do k=1,o
          do l=1,p
            write (99,'(f25.15)') x(i,j,k,l)
          end do
        end do
      end do
    end do
    
    close(99)
  END SUBROUTINE wri4file
  
  
  !**********************************************************************
  ! SUBROUTINE wriscalar
  !**********************************************************************
  SUBROUTINE wriscalar(x,filename)
    implicit none
    
    !dummy arguments declarations
    character(*)                :: filename
    double precision, intent(in):: x
    
    !Local variables declarations
    
    open (unit=99,file=filename)  
    write (99,'(f25.15)') x
    close(99)
  END SUBROUTINE wriscalar
  
  
  !***********************************************************************
  ! SUBROUTINE readfile
  !***********************************************************************
  SUBROUTINE readfile(m,n,x,filename)
    implicit none
    
    !dummy arguments declarations
    character(*)                                   :: filename
    integer, intent(in)                            :: m,n
    double precision, intent(out), dimension (m,n) :: x
    
    !Local variables declarations
    character (len=*), parameter  :: newline = '()'
    integer                       :: i,j
    
    
    open (unit=99,file=filename)
    
    do i=1,m
      read (99,'(5000f25.15)') (x(i,j), j=1, n)
    end do

    close(99)
    
  END SUBROUTINE readfile
  
  
  !***********************************************************************
  ! SUBROUTINE read3file
  !***********************************************************************
  SUBROUTINE read3file(m,n,o,x,filename)
    implicit none
    
    !dummy arguments declarations
    character(*)                                    :: filename
    integer, intent(in)                             :: m,n,o
    double precision, intent(out), dimension (m,n,o):: x
    
    !Local variables declarations
    character (len=*), parameter  :: newline = '()'
    integer                       :: i,j,k
    
    open (unit=99,file=filename)
    
    do i=1,m
      do j=1,n
        do k=1,o
          read (99,'(f25.15)') x(i,j,k)
        end do
      end do
    end do
    
    close(99)
    
  END SUBROUTINE read3file
  
  
  !***********************************************************************
  ! SUBROUTINE read4file
  !***********************************************************************
  SUBROUTINE read4file(m,n,o,p,x,filename)
    implicit none
    
    !dummy arguments declarations
    character(*)                                       :: filename
    integer, intent(in)                                :: m,n,o,p
    double precision, intent(out), dimension (m,n,o,p) :: x
    
    !Local variables declarations
    character (len=*), parameter  :: newline = '()'
    integer                       :: i,j,k,l
    
    
    open (unit=99,file=filename)
    
    do i=1,m
      do j=1,n
        do k=1,o
          do l=1,p
            read (99,'(f25.15)') x(i,j,k,l)
          end do
        end do
      end do
    end do
    
    close(99)
    
  END SUBROUTINE read4file
  
  
  !***********************************************************************
  ! SUBROUTINE readscalar
  !***********************************************************************
  SUBROUTINE readscalar(x,filename)
    implicit none
    
    !dummy arguments declarations
    character(*)                 :: filename
    double precision, intent(out):: x
    
    !Local variables declarations
    
    open (unit=99,file=filename)
    read (99,'(f25.15)') x
    close(99)
  END SUBROUTINE readscalar
  
  
  !**********************************************************************
  ! SUBROUTINE wrilogical
  !**********************************************************************
  SUBROUTINE wrilogical(x,filename)
    implicit none
    
    !dummy arguments declarations
    character(*)       :: filename
    logical, intent(in):: x
    
    !Local variables declarations
    
    open (unit=99,file=filename)  
    write (99,'(L1)') x
    close(99)
  END SUBROUTINE wrilogical
  
  
  !***********************************************************************
  ! SUBROUTINE readlogical
  !***********************************************************************
  SUBROUTINE readlogical(x,filename)
    implicit none
    
    !dummy arguments declarations
    character(*)        :: filename
    logical, intent(out):: x
    
    !Local variables declarations
    
    open (unit=99,file=filename)
    read (99,'(L1)') x
    close(99)
  END SUBROUTINE readlogical

  !**********************************************************************
  ! SUBROUTINE int1file
  !**********************************************************************
  SUBROUTINE int1file(n,x,filename)
  implicit none

  !dummy arguments declarations
  character(*)                          :: filename
  integer, intent(in)                   :: n
  integer, intent(in), dimension (n)  :: x

  !Local variables declarations
  character (len=*), parameter  :: newline = '()'
  integer                       :: j


  open (unit=99,file=filename)
  do j=1,n
  write (99,'(i10, 5x, $)') x(j)
  end do

  close(99)

  !LAST CARD OF SUBROUTINE int1file

  END SUBROUTINE int1file


  !**********************************************************************
  ! SUBROUTINE int2file
  !**********************************************************************
  SUBROUTINE int2file(m,n,x,filename)
    implicit none
    
    !dummy arguments declarations
    character(*)                          :: filename
    integer, intent(in)                   :: m,n
    integer, intent(in), dimension (m,n)  :: x
    
    !Local variables declarations
    character (len=*), parameter  :: newline = '()'
    integer                       :: i,j
    
    
    open (unit=99,file=filename)  
    do i=1,m
      do j=1,n
        write (99,'(i10, 5x, $)') x(i,j)
      end do
      write (99,newline)
    end do
    
    close(99)
    
    !LAST CARD OF SUBROUTINE int2file
    
  END SUBROUTINE int2file
  
  
  !**********************************************************************
  ! SUBROUTINE int3file
  !**********************************************************************
  SUBROUTINE int3file(m,n,o,x,filename)
    implicit none
    
    !dummy arguments declarations
    character(*)                          :: filename
    integer, intent(in)                   :: m,n,o
    integer, intent(in), dimension (m,n,o):: x
    
    !Local variables declarations
    integer                       :: i,j,k
    
    open (unit=99,file=filename)  
    do i=1,m
      do j=1,n
        do k=1,o
          write (99,'(i10)') x(i,j,k)
        end do
      end do
    end do
    
    close(99)
  END SUBROUTINE int3file
  
  
  !**********************************************************************
  ! SUBROUTINE readint3file
  !**********************************************************************
  SUBROUTINE readint3file(m,n,o,x,filename)
    implicit none
    
    !dummy arguments declarations
    character(*)                           :: filename
    integer, intent(in)                    :: m,n,o
    integer, intent(out), dimension (m,n,o):: x
    
    !Local variables declarations
    integer                       :: i,j,k
    
    open (unit=99,file=filename)  
    do i=1,m
      do j=1,n
        do k=1,o
          read (99,'(i10)') x(i,j,k)
        end do
      end do
    end do
    
    close(99)
  END SUBROUTINE readint3file
  
  
  !**********************************************************************
  ! SUBROUTINE int2file1
  !**********************************************************************
  SUBROUTINE int2file1(x,filename)
    implicit none
    
    !dummy arguments declarations
    character(*)       :: filename
    integer, intent(in):: x
    
    
    open (unit=99,file=filename)  
    write (99,'(i10, 5x, $)') x
    
    close(99)
    
    !LAST CARD OF SUBROUTINE int2file1
    
  END SUBROUTINE int2file1
  
  
  !**********************************************************************
  ! SUBROUTINE readintfile
  !**********************************************************************
  SUBROUTINE readintfile(m,n,x,filename)
    implicit none
    
    !dummy arguments declarations
    character(*)                          :: filename
    integer, intent(in)                   :: m,n
    integer, intent(out), dimension (m,n) :: x
    
    !Local variables declarations
    character (len=*), parameter  :: newline = '()'
    integer                       :: i,j
    
    
    open (unit=99,file=filename)  
    do i=1,m
!      do j=1,n
!        read (99,'(i10, 5x, $)') x(i,j)
      read (99,'(1500i10)') (x(i,j), j=1, n)
!      end do
!      write (99,newline)
    end do
    
    close(99)
    
    !LAST CARD OF SUBROUTINE int2file
    
  END SUBROUTINE readintfile
  
  
  !**********************************************************************
  ! SUBROUTINE readintfile1
  !**********************************************************************
  SUBROUTINE readintfile1(x,filename)
    implicit none
    
    !dummy arguments declarations
    character(*)        :: filename
    integer, intent(out):: x
    
    
    open (unit=99,file=filename)  
    read (99,'(i10, 5x, $)') x
    
    close(99)
    
    !LAST CARD OF SUBROUTINE readintfile1
  END SUBROUTINE readintfile1
END MODULE IOOP