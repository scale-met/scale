!-------------------------------------------------------------------------------
!> Program SCALE-RM (a launcher of main routine)
!!
!! @par Description
!!          SCALE: Scalable Computing by Advanced Library and Environment
!!          Limited area model for regional weather, regional climate, and large-Eddy Simulation (LES)
!!
!! @author Team SCALE
!!
!<
!-------------------------------------------------------------------------------
program scalerm_pp
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_io
  use scale_prof

  use scale_prc, only: &
     PRC_DOMAIN_nlim,     &
     PRC_COMM_NULL,       &
     PRC_abort,           &
     PRC_MPIstart,        &
     PRC_MPIfinish,       &
     PRC_MPIsplit_bulk,   &
     PRC_MPIsplit_nest,   &
     PRC_UNIVERSAL_setup, &
     PRC_GLOBAL_setup,    &
     PRC_GLOBAL_ROOT,     &
     PRC_ERRHANDLER_setup
  use scale_fpm, only: &
     FPM_Init
  use mod_rm_prep, only: &
     rm_prep
  use mod_rm_driver, only: &
     rm_driver
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

  logical               :: EXECUTE_PREPROCESS           = .true.  ! execute preprocess tools?
  logical               :: EXECUTE_MODEL                = .false. ! execute main model?
  integer               :: NUM_BULKJOB                  = 1       ! number of bulk jobs
  integer               :: NUM_BULKJOB_ONCE             = 1       ! number of bulk jobs for one iteration
  integer               :: NUM_ITERATION_BULK           = 1       ! number of iteration for bulk job
  integer               :: NUM_DOMAIN                   = 1       ! number of domains
  integer               :: NUM_FAIL_TOLERANCE           = 1       ! tolerance number of failure processes
  integer               :: FREQ_FAIL_CHECK              = 0       ! FPM polling frequency per DT (0: no polling)
  integer               :: PRC_DOMAINS(PRC_DOMAIN_nlim) = 0       ! number of total process in each domain
  character(len=H_LONG) :: CONF_FILES (PRC_DOMAIN_nlim) = ""      ! name of configulation files
  logical               :: ABORT_ALL_JOBS               = .false. ! abort all jobs or not?
  logical               :: LOG_SPLIT                    = .false. ! log-output for mpi splitting?
  logical               :: COLOR_REORDER                = .true.  ! coloring reorder for mpi splitting?
  logical               :: FAILURE_PRC_MANAGE           = .false. ! use failure process management?

  namelist / PARAM_LAUNCHER / &
!      EXECUTE_PREPROCESS, &
!      EXECUTE_MODEL,      &
     NUM_BULKJOB,        &
     NUM_ITERATION_BULK, &
     NUM_DOMAIN,         &
     NUM_FAIL_TOLERANCE, &
     FREQ_FAIL_CHECK,    &
     PRC_DOMAINS,        &
     CONF_FILES,         &
     ABORT_ALL_JOBS,     &
     LOG_SPLIT,          &
     COLOR_REORDER,      &
     FAILURE_PRC_MANAGE

  integer               :: universal_comm                   ! universal communicator
  integer               :: universal_nprocs                 ! number of procs in universal communicator
  integer               :: universal_myrank                 ! my rank         in universal communicator
  logical               :: universal_master                 ! master process  in universal communicator?
  character(len=H_LONG) :: universal_cnf_fname              ! config file for launcher

  integer               :: global_comm                      ! communicator for each member
  integer               :: global_nprocs                    ! number of procs in global communicator
  integer               :: PRC_BULKJOB(PRC_DOMAIN_nlim) = 0 ! number of procs in each bulk job = global_nprocs
  integer               :: ID_BULKJOB                       ! bulk job ID

  logical               :: use_fpm = .false.                ! switch for fpm module

  integer               :: local_comm                       ! assigned local communicator
  integer               :: ID_DOMAIN                        ! domain ID
  integer               :: intercomm_parent                 ! inter communicator with parent
  integer               :: intercomm_child                  ! inter communicator with child
  character(len=H_LONG) :: local_cnf_fname                  ! config file for local domain

  integer :: itr
  integer :: fid, ierr
  !-----------------------------------------------------------

  ! start MPI
  call PRC_MPIstart( universal_comm ) ! [OUT]

  call PRC_UNIVERSAL_setup( universal_comm,   & ! [IN]
                            universal_nprocs, & ! [OUT]
                            universal_myrank, & ! [OUT]
                            universal_master  ) ! [OUT]

  if( universal_master ) write(*,*) '*** Start Launch System for SCALE-RM'

  !--- read launcher config

  universal_cnf_fname = IO_ARG_getfname( universal_master )

  fid = IO_CNF_open( universal_cnf_fname, & ! [IN]
                     universal_master     ) ! [IN]

  ! set default
  CONF_FILES(1) = universal_cnf_fname

  ! read namelist
  rewind(fid)
  read(fid,nml=PARAM_LAUNCHER,iostat=ierr)
  if ( ierr < 0 ) then !--- missing
     ! keep default setting (no members, no nesting)
  elseif( ierr > 0 ) then !--- fatal error
     if( universal_master ) write(*,*) 'xxx Not appropriate names in namelist PARAM_LAUNCHER. Check!'
     call PRC_abort
  endif

  close(fid)
  FREQ_FAIL_CHECK = 0 ! force 0, coz no time integrations

  if (      EXECUTE_PREPROCESS &
       .OR. EXECUTE_MODEL      ) then
     if( universal_master ) write(*,*) "*** Execute preprocess? : ", EXECUTE_PREPROCESS
     if( universal_master ) write(*,*) "*** Execute model?      : ", EXECUTE_MODEL
  else
     if( universal_master ) write(*,*) 'xxx No execution. please check PARAM_LAUNCHER. STOP'
     call PRC_abort
  endif

  !--- split for bulk jobs

  if ( NUM_BULKJOB == 1 ) NUM_ITERATION_BULK = 1
  NUM_BULKJOB_ONCE = ceiling( real(NUM_BULKJOB) / NUM_ITERATION_BULK )
  if ( mod(universal_nprocs,NUM_BULKJOB_ONCE) /= 0 ) then !--- fatal error
     if( universal_master ) write(*,*) 'xxx Total Num of Processes must be divisible by NUM_BULKJOB. Check!'
     if( universal_master ) write(*,*) 'xxx Total Num of Processes           = ', universal_nprocs
     if( universal_master ) write(*,*) 'xxx NUM_BULKJOB                      = ', NUM_BULKJOB
     if( universal_master ) write(*,*) 'xxx NUM_ITERATION_BULK               = ', NUM_ITERATION_BULK
     if( universal_master ) write(*,*) 'xxx NUM_BULKJOB / NUM_ITERATION_BULK = ', NUM_BULKJOB_ONCE
     call PRC_abort
  endif

  global_nprocs = universal_nprocs / NUM_BULKJOB_ONCE
  PRC_BULKJOB(1:NUM_BULKJOB) = global_nprocs
  if ( NUM_BULKJOB > 1 ) then
     if( universal_master ) write(*,'(1x,A,I5)') "*** TOTAL # of BULK JOBS             = ", NUM_BULKJOB
     if( universal_master ) write(*,'(1x,A,I5)') "*** # of BULK JOB for each iteration = ", NUM_BULKJOB_ONCE
     if( universal_master ) write(*,'(1x,A,I5)') "*** # of PROCESS of each JOB         = ", global_nprocs

     if ( FAILURE_PRC_MANAGE ) then
        if( universal_master ) write(*,'(1x,A)') "*** Available: Failure Process Management"
        use_fpm = .true.                    !--- available only in bulk job
        if ( NUM_FAIL_TOLERANCE <= 0 ) then !--- fatal error
           if( universal_master ) write(*,*) 'xxx Num of Failure Processes must be positive number. Check!'
           if( universal_master ) write(*,*) 'xxx NUM_FAIL_TOLERANCE = ', NUM_FAIL_TOLERANCE
           call PRC_abort
        endif

        if ( NUM_FAIL_TOLERANCE > NUM_BULKJOB ) then !--- fatal error
           write(*,*) 'xxx NUM_FAIL_TOLERANCE is bigger than # of NUM_BLUKJOB'
           write(*,*) '    set to be: NUM_FAIL_TOLERANCE <= NUM_BLUKJOB'
           call PRC_abort
        endif
     endif
  endif

  ! communicator split for bulk/ensemble
  call PRC_MPIsplit_bulk( universal_comm,   & ! [IN]
                          NUM_BULKJOB_ONCE, & ! [IN]
                          PRC_BULKJOB(:),   & ! [IN]
                          LOG_SPLIT,        & ! [IN]
                          global_comm,      & ! [OUT]
                          ID_BULKJOB        ) ! [OUT]

  call PRC_GLOBAL_setup( ABORT_ALL_JOBS, & ! [IN]
                         global_comm     ) ! [IN]

  !--- split for nesting

  if ( NUM_DOMAIN > 1 ) then
     if( universal_master ) write(*,'(1x,A,I5)') "*** TOTAL DOMAIN NUMBER     = ", NUM_DOMAIN
     if( universal_master ) write(*,'(1x,A,L5)') "*** Flag of ABORT ALL JOBS  = ", ABORT_ALL_JOBS
  endif

  ! communicator split for nesting domains
  call PRC_MPIsplit_nest( global_comm,      & ! [IN]
                          NUM_DOMAIN,       & ! [IN]
                          PRC_DOMAINS(:),   & ! [IN]
                          LOG_SPLIT,        & ! [IN]
                          COLOR_REORDER,    & ! [IN]
                          local_comm,       & ! [OUT]
                          ID_DOMAIN,        & ! [OUT]
                          intercomm_parent, & ! [OUT]
                          intercomm_child   ) ! [OUT]

  !--- initialize FPM module & error handler
  call FPM_Init( NUM_FAIL_TOLERANCE, & ! [IN]
                 FREQ_FAIL_CHECK,    & ! [IN]
                 universal_comm,     & ! [IN]
                 global_comm,        & ! [IN]
                 local_comm,         & ! [IN]
                 NUM_BULKJOB,        & ! [IN]
                 PRC_GLOBAL_ROOT,    & ! [IN]
                 use_fpm             ) ! [IN]

  call PRC_ERRHANDLER_setup( use_fpm, universal_master )

  !--- start main routine

  do itr = 1, NUM_ITERATION_BULK

     if ( ID_BULKJOB > NUM_BULKJOB ) exit

     call IO_set_universalrank( universal_myrank, & ! [IN]
                                ID_BULKJOB,       & ! [IN]
                                ID_DOMAIN         ) ! [IN]

     if ( NUM_BULKJOB > 1 ) then
        write(local_cnf_fname,'(I4.4,2A)') ID_BULKJOB, "/", trim(CONF_FILES(ID_DOMAIN))
     else
        local_cnf_fname = trim(CONF_FILES(ID_DOMAIN))
     endif

     if ( EXECUTE_PREPROCESS ) then
        call rm_prep( local_comm,     & ! [IN]
                      PRC_COMM_NULL,  & ! [IN]
                      PRC_COMM_NULL,  & ! [IN]
                      local_cnf_fname ) ! [IN]
     endif

     if ( EXECUTE_MODEL ) then
        call rm_driver( local_comm,       & ! [IN]
                        intercomm_parent, & ! [IN]
                        intercomm_child,  & ! [IN]
                        local_cnf_fname   ) ! [IN]
     endif

     ID_BULKJOB = ID_BULKJOB + NUM_BULKJOB_ONCE

  end do

  ! stop MPI
  call PRC_MPIfinish

  if( universal_master ) write(*,*) '*** End   Launch System for SCALE-RM'

  stop
end program scalerm_pp
