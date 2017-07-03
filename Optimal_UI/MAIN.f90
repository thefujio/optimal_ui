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
  USE TOOLBOX
!  USE POWELL
!  USE omp_lib
!  USE isnan_int
  implicit none
  !Variables declarations
  real(8):: time_begin,time_end            !Used to keep track of time
!  real(8):: start_iter, end_iter           !Used to keep track of time
  real(8), dimension(dims):: paramvec
!  real(8):: rhobeg, rhoend
  integer:: i !,iprint, maxfun
!  real(8), dimension((2*dims+5)*(2*dims+dims)+3*dims*(dims+5)/2):: wspace
!Timing of the program
  call CPU_Time(time_begin)
  ! CHANGING IMSL HANDLING OF ERRORS
!  CALL erset(4, 1, 0) !this prevents BCONF from stopping when the max # of iter is reached
!  CALL erset(0, 0, 0) !this prevents BCONF from stopping OR printing when the max # of iter is reached

  !Open the output file where details are printed
  open(unit=detail,  file=root_dir//out_dir//"detail.txt",  status='replace')
  open(unit=calibout,  file=root_dir//out_dir//"calibout.txt",  status='replace')
  if (transform == 1) then
    do i=1,dims
      paramvec(i) = unrestrict(guess(i),lb(i),ub(i),spread)
    enddo
  else
    do i=1,dims
      paramvec(i)       =   guess(i)
    enddo
    print*,paramvec
  endif


  call linspace(bgrid,0.4d0,0.89d0,gridpoints)

  print *, "Run bisection method to find tau for each rr in grid"
  !For non-calibration testing:
  do i=1,gridpoints
    bval = bgrid(i)
    Call calfun(dims,paramvec,funcerror)
    rrgrid(i) = rrval
    cegrid(i) = ceval
    taxgrid(i) = taxval
    jfpgrid(i) = jfpval
    uvalgrid(i) = uval
    submktgrid(i) = submktval
    grosswagegrid(i) = grosswageval
    netwagegrid(i) = netwageval
    urategrid(i) = urateval
    uugrid(i) = uuval
    eegrid(i) = eeval
    print*,'gridpoint ', i ,' completed'
  end do

  open(unit=gridout,  file=root_dir//out_dir//"gridout.txt",  status='replace')
  write(gridout,15)
  15 format('b,','RR,','CE,','Tax,','JFP(U),','U VF,','Avg. Open Submkt,','Gross W,','Net W,','Urate,','UU,','EE,')
    do i=1,gridpoints
      write(gridout,'(<12>(f15.4,","))') bgrid(i),rrgrid(i),cegrid(i),taxgrid(i),jfpgrid(i),uvalgrid(i),submktgrid(i),grosswagegrid(i),netwagegrid(i),urategrid(i),uugrid(i),eegrid(i)
    enddo
  close(gridout)


  call CPU_Time(time_end)
  !Timing of the program
  write(*,"('TOTAL EXECUTION TIME (seconds): ',f8.2)") &
  (time_end-time_begin)/8.0d0
  STOP

END PROGRAM MAIN

