MODULE PARAM
  !****************************************************************************
  !  PARAM.f90
  !  
  !  PURPOSE:
  !   Sets all parameters and global variables for model solution
  !  
  !  
  !  Program written by: M. Gervais
  !   Date:    Jan 2016
  !   Updated: Feb 2016
  !****************************************************************************

  implicit none
  
  ! CONSTANTS
  real(8), parameter:: eps   = 1.0d-10
  real(8), parameter:: zero  = 0.00d0
  real(8), parameter:: half  = 0.50d0
  real(8), parameter:: one   = 1.00d0
  real(8), parameter:: two   = 2.00d0
  real(8), parameter:: three = 3.00d0
  real(8), parameter:: hell  =-1.00d7
  
  
  ! DIRECTORIES AND FILES
  character(LEN=*), parameter:: root_dir = "/Users/fujio/Fortran/Optimal_UI/Optimal_UI/"
  character(LEN=*), parameter:: out_dir  = "output/"
  integer,          parameter:: detail = 7
  integer,          parameter:: calibout = 8
  integer,          parameter:: gridout = 9


  integer, parameter:: nparams = 21        !For storage purposes
  real(8), dimension(nparams):: params     !To store values of parameters
  ! LOGICAL VARIABLES TO CONTROL INITIAL VALUES AND REFINEMENTS
  logical:: use_old_v   = .false.
  logical:: do_howard   = .false.
  logical:: want_howard = .false.
  logical:: do_refine   = .false.
  logical:: want_refine = .true.
  logical:: simul_only  = .false.
  logical:: want_print  = .false.

  !Grid over b
  integer, parameter        :: gridpoints = 1
  integer, parameter        :: bgridpoints = 1
  real(8)                   :: bval
  real(8)                   :: rrval, ceval, taxval, jfpval, uval, uwgt, umeasure, submktval, submktwgt, submktmeasure, grosswageval, netwageval, urateval, uuval, eeval
  real(8), dimension(gridpoints) :: psigrid, durgrid, bgrid, rrgrid, cegrid, taxgrid, jfpgrid, uvalgrid, submktgrid, grosswagegrid, netwagegrid, urategrid, uugrid, eegrid
  !Flag for whether optimization routine used is going to be unconstrained (1) or bound-constrained (0)
  integer,          parameter:: transform = 0
  !Calibration Targets:
  real(8), parameter        :: jfptarget = 0.343d0
  real(8), parameter        :: urtarget = 0.0537d0
  real(8), parameter        :: j_jtarget = 0.0223d0
  real(8)                   :: funcerror

  !Initial Guess on Calibrated Parameters
  integer, parameter        :: dims = 3
  !real(8), dimension(dims), parameter   :: guess     = (/ 0.13100000d0, 0.02280000d0, 1.00000000d0 /) !kappa, delta, lambda
  real(8), dimension(dims), parameter   :: guess     = (/ 0.08683515d0, 0.01946000d0, 0.41697092d0 /)
  real(8), dimension(dims), parameter   :: lb        = (/ 0.02000000d0, 0.01000000d0, 0.15000000d0 /)
  real(8), dimension(dims), parameter   :: ub        = (/ 0.20000000d0, 0.03000000d0, 1.00000000d0 /)
  real(8)                               :: kappa, delta, lambda
  !Spread for restr/unrestr function - must be larger than 1!
  real(8), parameter:: spread = 2.000d0
  real(8), parameter:: pie = 3.14159265d0

  ! GRID ON PV
  integer, parameter        :: nx = 30
  real(8)                   :: xmin,xmax
  real(8), dimension(nx)    :: x
    
  ! STOCHASTIC AGGREGATE PRODUCTIVITY
  integer, parameter        :: ny = 1
  !NBER avg lenght (month) of recession (1945-2009)
  real(8), parameter        :: recession_length = 11.1d0
  !NBER avg lenght (month) of booms (1945-2009)
  real(8), parameter        :: boom_length = 58.4d0
  real(8), dimension(ny)    :: y,pyss
  real(8), dimension(ny,ny) :: py
  
  ! STOCHASTIC IDIOSYNCRATIC PRODUCTIVITY
  integer, parameter        :: nz = 3
  real(8), dimension(nz)    :: z,pzss,Pztilde
  real(8), dimension(nz,nz) :: pz
  
  ! INDEX CONVERSION
  integer, parameter        :: ns = ny*nz
  integer, dimension(ns)    :: iyfun,izfun
  integer, dimension(ny,nz) :: isfun
  real(8), dimension(ns)    :: psss
  real(8), dimension(ns,ns) :: ps

  !Stochastic Expiration of UI benefits
  ! GRID ON UI benefit (past wage+ineligible)
  !integer, parameter        :: ne = 31
  integer, parameter        :: ne=1
  real(8)                   :: emin,emax
  real(8), dimension(ne)    :: e

  real(8)                   :: psi
  real(8), dimension(ne,ne) :: pe
  !Index Conversion for unemployment transmat
  integer, parameter        :: nu = ny*ne
  integer, dimension(nu)    :: iuyfun,iuefun
  integer, dimension(ny,ne) :: iufun
  real(8), dimension(ne)    :: bvec
  real(8), dimension(nu,nu) :: pus

  !About 46% of avg wage
  real(8), parameter:: bmin = 0.250d0
  !Home production plus UI = about 2/3 of avg. wage
  real(8), parameter:: hp = 0.430d0
  ! CONTRACT SPACE
  integer:: nc
  integer, allocatable:: cont(:,:)
  
  ! CURRENT VALUE FUNCTIONS AND POLICY RULES
  real(8), dimension(nx,ny)   :: theta,P
  real(8), dimension(nx,ny,ne):: dprimevec
!  real(8), dimension(nx,ny,nz):: J
  real(8), dimension(nx,ny,nz):: w
  integer, dimension(nx,ny,nz):: wind
  real(8), dimension(nx,ns,ns):: dprime
  integer, dimension(nx,ns,ns):: iVprime
  real(8), dimension(nx,ny)   :: R,Ptilde
  integer, dimension(nx,ny)   :: M
!  real(8), dimension(ny)      :: U
  real(8), dimension(ny,ne)   :: RU,PUtilde
  integer, dimension(ny,ne)   :: MU
  
  !STATIONARY DISTRIBUTION
  real(8), dimension(ns*nx+nu,ns*nx+nu):: pimat
  real(8), dimension(ns*nx+nu):: muss
  
  ! MODEL PARAMETERS
  ! INDIVIDUALS / PREFERENCES
  real(8), parameter:: betta = 0.995942407351067d0 !one/(1.05d0**(1/12))
  real(8), parameter:: eta   = zero
  real(8), parameter:: sigma = 2.000d0
  
  ! MATCHING TECHNOLOGY, SEARCH, SEPARATION, COST
  real(8), parameter:: gamma  = 0.5000d0
  !real(8), parameter:: lambda = 1.000d0
  !real(8), parameter:: delta  = 0.026d0
  !real(8), parameter:: kappa  = 0.029d0
  
  !INITIAL CONDITIONS AND POLICY PARAMETERS
  real(8):: tau,rr

  !Aggregate Statistics
  real(8):: ee,eu,ue,unemp,em,EUflow,EEflow,UEflow
  real(8):: tot_wage,avg_wage,avg_benefit,transfers,welfare
  integer:: avg_wage_ind

  !PROGRAMMING PARAMETERS
  ! FLOW CONTROLS
  integer, parameter:: niter = 500
  integer, parameter:: nupdate = 30
  
  ! TOLERANCE LEVEL
  real(8):: tol = 1.0d-8
  real(8), parameter:: high_tol = 1.0d-6
  real(8), parameter:: low_tol  = 1.0d-10
  real(8), parameter:: errrel   = 1.0d-12
  
END MODULE PARAM
