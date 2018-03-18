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
  INTERFACE
    SUBROUTINE ols_fit(xin,yin,outvec)
      implicit none
      real*8, dimension(:), intent(in) :: xin,yin
      real*8, dimension(:), intent(out) :: outvec
    END SUBROUTINE ols_fit

    SUBROUTINE tau_eval(n,xvec,tax,fval)
      implicit none
      integer, intent(IN):: n
      real*8, intent(IN), dimension(n):: xvec
      real*8, intent(IN):: tax
      real*8, intent(OUT):: fval
    END SUBROUTINE
  END INTERFACE

  !Variables declarations
  real(8):: time_begin,time_end            !Used to keep track of time
!  real(8):: start_iter, end_iter           !Used to keep track of time
  real(8), dimension(dims):: paramvec
!  real(8):: rhobeg, rhoend
  integer:: i,j !,iprint, maxfun
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

  !0.233333333d0
  call linspace(bgrid,0.210d0,0.270d0,gridpoints)
  !call linspace(hpgrid,0.43d0,0.78d0,gridpoints)

  !call linspace(psigrid,1.0d0,0.0d0,gridpoints)
  !call linspace(durgrid,1.0d0,36.0d0,gridpoints-1)
  !durgrid = (/ 6.0d0, (6.0d0*0.95d0), (6.0d0*1.05d0) /)
  !do i=1,gridpoints
  !psigrid(i) = 1.0d0/durgrid(i)
  !enddo
  call linspace(yygrid,0.00d0,0.10d0,bgridpoints)
  !bgrid = 0.25d0
  !print*, psigrid
  !pause
  print *, "Run secant method to find tau for each rr in grid"
  open(unit=gridout,  file=root_dir//out_dir//"gridout.txt",  status='replace')
  write(gridout,15)
  15 format('Psi,','RR,','CE,','Tax,','JFP(U),','U VF,','Avg. Open Submkt,','Gross W,','Net W,',&
  'Urate,','UU,','EE,','b,','transfers,','Utility,','VF','U VF ben,','U VF noben,')
  !For non-calibration testing:
  do j=1,bgridpoints
    do i=1,gridpoints
      !hp = hpgrid(i)
      yyval = yygrid(j)
      bval = bgrid(i)
      psi = psigrid(i)
      !
      Call calfun(dims,paramvec,funcerror)
      !
      rrgrid(i) = rrval
      cegrid(i) = ceval
      taxgrid(i) = taxval
      bdgrid(i) = bdval
      jfpgrid(i) = jfpval
      uvalgrid(i) = uval
      submktgrid(i) = submktval
      grosswagegrid(i) = grosswageval
      netwagegrid(i) = netwageval
      urategrid(i) = urateval
      uugrid(i) = uuval
      eegrid(i) = eeval
      trgrid(i) = trval
      utilgrid(i)=tot_util
      vfgrid(i)=tot_vf
      uvalbengrid(i) = ubenval
      uvalnobengrid(i) = unobenval
      print*,'gridpoint ', i ,' completed'
    end do

  do i=1,gridpoints
    write(gridout,'(<20>(f15.4,","))') yyval, psigrid(i),rrgrid(i),cegrid(i),taxgrid(i),jfpgrid(i),uvalgrid(i), &
      submktgrid(i), grosswagegrid(i),netwagegrid(i),urategrid(i),uugrid(i),eegrid(i),bgrid(i),trgrid(i), &
      utilgrid(i),vfgrid(i),uvalbengrid(i),uvalnobengrid(i),bdgrid(i)
  enddo
  enddo !end j grid
  close(gridout)

!  !Smooth the tax and run at each point in grid:
  call ols_fit(bgrid,taxgrid,taxgrid_smooth)
  print*,'taxgrid_smooth',taxgrid_smooth
  print *, "Run for final smoothed tau for each rr in grid"
  open(unit=gridout_smooth,  file=root_dir//out_dir//"gridout_smooth.txt",  status='replace')
  write(gridout_smooth,19)
  19 format('Psi,','RR,','CE,','Tax,','JFP(U),','U VF,','Avg. Open Submkt,','Gross W,','Net W,',&
  'Urate,','UU,','EE,','b,','transfers,','Utility,','VF','U VF ben,','U VF noben,')
  !For non-calibration testing:
  do j=1,bgridpoints
    do i=1,gridpoints
      !hp = hpgrid(i)
      yyval = yygrid(j)
      bval = bgrid(i)
      psi = psigrid(i)
      tau = taxgrid_smooth(i)
      !
      Call tau_eval(dims,paramvec,tau,funcerror)
      !
      rrgrid(i) = rrval
      cegrid(i) = ceval
      taxgrid(i) = taxval
      bdgrid(i) = bdval
      jfpgrid(i) = jfpval
      uvalgrid(i) = uval
      submktgrid(i) = submktval
      grosswagegrid(i) = grosswageval
      netwagegrid(i) = netwageval
      urategrid(i) = urateval
      uugrid(i) = uuval
      eegrid(i) = eeval
      trgrid(i) = trval
      utilgrid(i)=tot_util
      vfgrid(i)=tot_vf
      uvalbengrid(i) = ubenval
      uvalnobengrid(i) = unobenval
      print*,'gridpoint ', i ,' completed'
    end do

  do i=1,gridpoints
    write(gridout_smooth,'(<20>(f15.4,","))') yyval, psigrid(i),rrgrid(i),cegrid(i),taxgrid(i),jfpgrid(i),uvalgrid(i), &
      submktgrid(i), grosswagegrid(i),netwagegrid(i),urategrid(i),uugrid(i),eegrid(i),bgrid(i),trgrid(i), &
      utilgrid(i),vfgrid(i),uvalbengrid(i),uvalnobengrid(i),bdgrid(i)
  enddo
  enddo !end j grid
  close(gridout_smooth)


  call CPU_Time(time_end)
  !Timing of the program
  write(*,"('TOTAL EXECUTION TIME (seconds): ',f8.2)") &
  (time_end-time_begin)/8.0d0
  STOP

END PROGRAM MAIN

