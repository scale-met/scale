!-------------------------------------------------------------------------------
!> Program SCALE-LES (a launcher of main routine)
!!
!! @par Description
!!          SCALE: Scalable Computing by Advanced Library and Environment
!!          Numerical model for LES-scale weather
!!
!! @author Team SCALE
!!
!<
!-------------------------------------------------------------------------------
program scaleles_launcher
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use dc_log, only: &
     LogInit
  use gtool_file, only: &
     FileCloseAll
  use scale_precision
  use scale_stdio
  use scale_prof
  use scale_process, only: &
     PRC_setup,         &
     PRC_MPIstart,      &
     PRC_MPIsplit,      &
     PRC_MPIfinish,     &
     PRC_MPIstop,       &
     MASTER_LOG,        &
     MASTER_COMM_WORLD, &
     max_depth
  use scale_fileio, only: &
     FILEIO_setup
  use scale_comm, only: &
     COMM_setup
  use mod_les_driver
  !
  !-----------------------------------------------------------------------------
  implicit none
  !-----------------------------------------------------------------------------
  !
  !++ included parameters
  !
  !-----------------------------------------------------------------------------
  !
  !++ parameters & variables
  !
  !=============================================================================

  integer :: myrank_local                 ! rank number in local communicator
  integer :: nmax_local                   ! # of total process in local communicator
  integer :: nmax_parent                  ! # of total process in parent world
  integer :: nmax_child                   ! # of total process in child world
  integer :: inter_parent                 ! inter communicator with parent
  integer :: inter_child                  ! inter communicator with child
  integer :: NUM_DOMAIN                   ! number of domains
  integer :: PRC_DOMAINS(max_depth)       ! # of total process in each domain
  integer :: MY_COMM_WORLD                ! assigned local communicator

  logical :: LOG_SPLIT = .false.      ! flag of log-output for mpi splitting

  character(len=H_LONG) :: CONF_FILES(max_depth)  ! names of configulation files
  character(len=H_LONG) :: fname_launch           ! config file for launcher
  character(len=H_LONG) :: fname_local            ! config file for local domain

  integer :: ierr
  integer :: LNC_FID_CONF

  namelist / PARAM_LAUNCHER / &
     NUM_DOMAIN,     &
     PRC_DOMAINS,    &
     CONF_FILES,     &
     LOG_SPLIT
  !-----------------------------------------------------------

  NUM_DOMAIN       = 1
  PRC_DOMAINS      = 0
  CONF_FILES(:)    = ""

  ! start MPI
  call PRC_MPIstart

  !--- read from argument
  if ( COMMAND_ARGUMENT_COUNT() < 1 ) then
     if (MASTER_LOG) write(*,*) 'WARNING: No config file specified!'
     if (MASTER_LOG) write(*,*) '         Default values are used.'
  else
     call get_command_argument(1,fname_launch)
  endif
  CONF_FILES(1) = fname_launch

  !--- open config file till end
  LNC_FID_CONF = IO_get_available_fid()
  open( LNC_FID_CONF,                &
        file   = trim(fname_launch), &
        form   = 'formatted',        &
        status = 'old',              &
        iostat = ierr                )
  if ( ierr /= 0 ) then
     if (MASTER_LOG) write(*,*)
     if (MASTER_LOG) write(*,*) 'WARNING: Failed to open config file! :', trim(fname_launch)
     if (MASTER_LOG) write(*,*) '         Default values are used.'
     if (MASTER_LOG) write(*,*)
  end if

  !--- read PARAM
  rewind(LNC_FID_CONF)
  read(LNC_FID_CONF,nml=PARAM_LAUNCHER,iostat=ierr)

  if( ierr > 0 ) then !--- fatal error
     if (MASTER_LOG) write(*,*) ' xxx Not appropriate names in namelist PARAM_LAUNCHER . Check!'
     call PRC_MPIstop
  endif

  close( LNC_FID_CONF )

  if ( MASTER_LOG ) write(*,*) '*** Start Launch System for SCALE-LES'
  if ( MASTER_LOG ) write ( *, '(1X,A,I2)' ) "TOTAL DOMAIN NUMBER = ", NUM_DOMAIN


  ! split MPI communicator
  call PRC_MPIsplit(     &
      MASTER_COMM_WORLD, & ! [in ]
      NUM_DOMAIN,        & ! [in ]
      PRC_DOMAINS,       & ! [in ]
      CONF_FILES,        & ! [in ]
      LOG_SPLIT,         & ! [in ]
      MY_COMM_WORLD,     & ! [out]
      inter_parent,      & ! [out]
      inter_child,       & ! [out]
      fname_local        ) ! [out]


  ! start main routine
  call scaleles ( MY_COMM_WORLD,   &
                  inter_parent,    &
                  inter_child,     &
                  fname_local      )


  ! stop MPI
  call PRC_MPIfinish

  stop
  !=============================================================================
  !
end program scaleles_launcher
