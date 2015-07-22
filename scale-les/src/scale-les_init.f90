!-------------------------------------------------------------------------------
!> Program SCALE-LES init (a launcher of main routine)
!!
!! @par Description
!!          This program makes initial data for ideal test cases
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2012-04-08 (H.Yashiro)  [mod] merge all init programs
!! @li      2012-06-13 (Y.Sato)     [mod] add HBINW option
!!
!<
!-------------------------------------------------------------------------------
program scaleles_init_launcher
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
     PRC_BULKsetup,     &
     MASTER_LOG,        &
     MASTER_COMM_WORLD, &
     MASTER_nmax,       &
     max_depth
  use scale_fileio, only: &
     FILEIO_setup
  use scale_comm, only: &
     COMM_setup
  use mod_init_driver
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

  integer :: NUM_BULKJOB                  ! number of bulk jobs
  integer :: NUM_DOMAIN                   ! number of domains
  integer :: PRC_BULKJOB(max_depth)       ! # of total process in each bulk job
  integer :: PRC_DOMAINS(max_depth)       ! # of total process in each domain
  integer :: BULK_COMM_WORLD              ! split communicator for each bulk job
  integer :: MY_COMM_WORLD                ! assigned local communicator
  integer :: inter_parent                 ! inter communicator with parent
  integer :: inter_child                  ! inter communicator with child
  integer :: inter_parent_null            ! inter communicator with parent
  integer :: inter_child_null             ! inter communicator with child

  logical :: ABORT_ALL_JOBS = .false.     ! flag of abort all jobs or not
  logical :: LOG_SPLIT = .false.          ! flag of log-output for mpi splitting

  character(len=H_LONG) :: CONF_FILES(max_depth)  ! names of configulation files
  character(len=H_LONG) :: CONF_BULKS(max_depth)  ! names of configulation files (dummy)
  character(len=H_LONG) :: fname_launch           ! config file for launcher
  character(len=H_LONG) :: fname_local            ! config file for local domain
  character(len=H_LONG) :: fname_bulks            ! config file for dummy use
  character(len=H_LONG) :: fname_scaleles         ! config file for scaleles

  integer :: check, nprocs
  integer :: ierr
  integer :: LNC_FID_CONF

  namelist / PARAM_LAUNCHER / &
     NUM_BULKJOB,     &
     NUM_DOMAIN,      &
     PRC_DOMAINS,     &
     CONF_FILES,      &
     ABORT_ALL_JOBS,  &
     LOG_SPLIT
  !-----------------------------------------------------------

  NUM_BULKJOB      = 1
  NUM_DOMAIN       = 1
  PRC_BULKJOB(:)   = 0
  PRC_DOMAINS(:)   = 0
  CONF_FILES(:)    = ""
  CONF_BULKS(:)    = ""

  ! start MPI
  call PRC_MPIstart
  if ( MASTER_LOG ) write(*,*) '*** Start Launch System for SCALE-LES'

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

  !--- read PARAM_LAUNCHER
  rewind(LNC_FID_CONF)
  read(LNC_FID_CONF,nml=PARAM_LAUNCHER,iostat=ierr)
  if( ierr > 0 ) then !--- fatal error
     if (MASTER_LOG) write(*,*) ' xxx Not appropriate names in namelist PARAM_LAUNCHER. Check!'
     call PRC_MPIstop
  endif
  close( LNC_FID_CONF )
  fname_local = CONF_FILES(1)

  !--- communicator split for bulk jobs
  check = mod(MASTER_nmax,NUM_BULKJOB)
  if( check /= 0 ) then !--- fatal error
     if (MASTER_LOG) write(*,*) ' xxx Total Num of Processes must be divisible by NUM_BULKJOB. Check!'
     if (MASTER_LOG) write(*,*) '     Total Num of Processes: ', MASTER_nmax
     if (MASTER_LOG) write(*,*) '                NUM_BULKJOB: ', NUM_BULKJOB
     call PRC_MPIstop
  endif

  nprocs = MASTER_nmax / NUM_BULKJOB
  PRC_BULKJOB(1:NUM_BULKJOB) = nprocs
  if ( MASTER_LOG ) write ( *, '(1X,A,I4)' ) "TOTAL BULK JOB NUMBER   = ", NUM_BULKJOB
  if ( MASTER_LOG ) write ( *, '(1X,A,I5)' ) "PROCESS NUM of EACH JOB = ", nprocs
  if ( MASTER_LOG ) write ( *, '(1X,A,I4)' ) "TOTAL DOMAIN NUMBER     = ", NUM_DOMAIN
  if ( MASTER_LOG ) write ( *, * ) "Flag of ABORT ALL JOBS  = ", ABORT_ALL_JOBS

  call PRC_MPIsplit(     &
      MASTER_COMM_WORLD, & ! [in ]
      NUM_BULKJOB,       & ! [in ]
      PRC_BULKJOB,       & ! [in ]
      CONF_BULKS,        & ! [in ] dummy
      LOG_SPLIT,         & ! [in ]
      .true.,            & ! [in ] flag bulk_split
      BULK_COMM_WORLD,   & ! [out]
      inter_parent_null, & ! [out] dummy
      inter_child_null,  & ! [out] dummy
      fname_bulks        ) ! [out]

  call PRC_BULKsetup(    &
      BULK_COMM_WORLD,   & ! [in]
      ABORT_ALL_JOBS     ) ! [in]


  !--- communicator split for nesting domains
  call PRC_MPIsplit(     &
      BULK_COMM_WORLD,   & ! [in ]
      NUM_DOMAIN,        & ! [in ]
      PRC_DOMAINS,       & ! [in ]
      CONF_FILES,        & ! [in ]
      LOG_SPLIT,         & ! [in ]
      .false.,           & ! [in ] flag bulk_split
      MY_COMM_WORLD,     & ! [out]
      inter_parent,      & ! [out]
      inter_child,       & ! [out]
      fname_local        ) ! [out]


  ! start main routine
  if ( NUM_BULKJOB > 1 ) then
     fname_scaleles = trim(fname_bulks)//"/"//trim(fname_local)
  else
     fname_scaleles = fname_local
  endif
  call scaleles_init ( MY_COMM_WORLD,     &
                       inter_parent_null, &
                       inter_child_null,  &
                       fname_scaleles     )


  ! stop MPI
  call PRC_MPIfinish

  stop
  !=============================================================================
end program scaleles_init_launcher
