!-----------------------------------------------------------------------------------------
!> post-process for scale high-performance (under the nickname of POPSCA)
!!
!! @par Description
!!      convert from netcdf to grads format (with combine and slice)
!!
!! @author Team SCALE
!!
!! @par History
!! @li  2015-02-03 (R.Yoshida)  original
!<
!-----------------------------------------------------------------------------------------
program netcdf2grads_h_launcher
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use gtool_file, only: &
     FileCloseAll
  use scale_precision
  use scale_stdio
  use scale_prof

  use scale_process, only: &
     PRC_DOMAIN_nlim,     &
     PRC_MPIstart,        &
     PRC_MPIstop,         &
     PRC_MPIfinish,       &
     PRC_MPIsplit,        &
     PRC_UNIVERSAL_setup, &
     PRC_GLOBAL_setup
  use mod_net2g_comm
  use mod_netcdf2grads_h
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

  integer               :: NUM_BULKJOB                  = 1       ! number of bulk jobs
  integer               :: NUM_DOMAIN                   = 1       ! number of domains
  integer               :: PRC_DOMAINS(PRC_DOMAIN_nlim) = 0       ! number of total process in each domain
  character(len=H_LONG) :: CONF_FILES (PRC_DOMAIN_nlim) = ""      ! name of configulation files
  logical               :: ABORT_ALL_JOBS               = .false. ! abort all jobs or not?
  logical               :: LOG_SPLIT                    = .false. ! log-output for mpi splitting?
  logical               :: COLOR_REORDER                = .true.  ! coloring reorder for mpi splitting?

  namelist / PARAM_LAUNCHER / &
     NUM_BULKJOB,     &
     NUM_DOMAIN,      &
     PRC_DOMAINS,     &
     CONF_FILES,      &
     ABORT_ALL_JOBS,  &
     LOG_SPLIT,       &
     COLOR_REORDER

  integer               :: universal_comm                         ! universal communicator
  integer               :: universal_nprocs                       ! number of procs in universal communicator
  logical               :: universal_master                       ! master process  in universal communicator?
  character(len=H_LONG) :: universal_cnf_fname                    ! config file for launcher

  integer               :: global_comm                            ! communicator for each member
  integer               :: global_nprocs                          ! number of procs in global communicator
  integer               :: PRC_BULKJOB(PRC_DOMAIN_nlim) = 0       ! = global_nprocs
  character(len=H_LONG) :: dummy1     (PRC_DOMAIN_nlim) = ""
  integer               :: intercomm_parent_null                  ! NULL inter communicator with parent
  integer               :: intercomm_child_null                   ! NULL inter communicator with child
  character(len=H_LONG) :: bulk_prefix                            ! dirname of each member

  integer               :: local_comm       ! assigned local communicator
  integer               :: intercomm_parent ! inter communicator with parent
  integer               :: intercomm_child  ! inter communicator with child
  character(len=H_LONG) :: local_cnf_fname  ! config file for local domain

  integer :: fid, ierr
  !-----------------------------------------------------------

  ! start MPI
  call PRC_MPIstart( universal_comm ) ! [OUT]

  call PRC_UNIVERSAL_setup( universal_comm,   & ! [IN]
                            universal_nprocs, & ! [OUT]
                            universal_master  ) ! [OUT]

  if( universal_master ) write(*,*) '*** Start Launch System for SCALE (POPSCA)'

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
     call PRC_MPIstop
  endif

  close(fid)

  !--- split for bulk jobs

  if ( mod(universal_nprocs,NUM_BULKJOB) /= 0 ) then !--- fatal error
     if( universal_master ) write(*,*) 'xxx Total Num of Processes must be divisible by NUM_BULKJOB. Check!'
     if( universal_master ) write(*,*) 'xxx Total Num of Processes = ', universal_nprocs
     if( universal_master ) write(*,*) 'xxx            NUM_BULKJOB = ', NUM_BULKJOB
     call PRC_MPIstop
  endif

  global_nprocs = universal_nprocs / NUM_BULKJOB
  PRC_BULKJOB(1:NUM_BULKJOB) = global_nprocs
  if ( NUM_BULKJOB > 1 ) then
     if( universal_master ) write(*,'(1x,A,I5)') "*** TOTAL BULK JOB NUMBER   = ", NUM_BULKJOB
     if( universal_master ) write(*,'(1x,A,I5)') "*** PROCESS NUM of EACH JOB = ", global_nprocs
  endif

  ! communicator split for bulk/ensemble
  call PRC_MPIsplit( universal_comm,        & ! [IN]
                     NUM_BULKJOB,           & ! [IN]
                     PRC_BULKJOB(:),        & ! [IN]
                     dummy1 (:),            & ! [IN]  dummy
                     LOG_SPLIT,             & ! [IN]
                     .true.,                & ! [IN]  flag bulk_split
                     .false.,               & ! [IN]  no reordering
                     global_comm,           & ! [OUT]
                     intercomm_parent_null, & ! [OUT] null
                     intercomm_child_null,  & ! [OUT] null
                     bulk_prefix            ) ! [OUT] dir name instead of file name

  call PRC_GLOBAL_setup( ABORT_ALL_JOBS, & ! [IN]
                         global_comm     ) ! [IN]

  !--- split for nesting

  if ( NUM_DOMAIN > 1 ) then
     if( universal_master ) write(*,'(1x,A,I5)') "*** TOTAL DOMAIN NUMBER     = ", NUM_DOMAIN
     if( universal_master ) write(*,'(1x,A,L5)') "*** Flag of ABORT ALL JOBS  = ", ABORT_ALL_JOBS
  endif

  ! communicator split for nesting domains
  call PRC_MPIsplit( global_comm,      & ! [IN]
                     NUM_DOMAIN,       & ! [IN]
                     PRC_DOMAINS(:),   & ! [IN]
                     CONF_FILES (:),   & ! [IN]
                     LOG_SPLIT,        & ! [IN]
                     .false.,          & ! [IN] flag bulk_split
                     COLOR_REORDER,    & ! [IN]
                     local_comm,       & ! [OUT]
                     intercomm_parent, & ! [OUT]
                     intercomm_child,  & ! [OUT]
                     local_cnf_fname   ) ! [OUT]

  !--- start main routine

  if ( NUM_BULKJOB > 1 ) then
     local_cnf_fname = trim(bulk_prefix)//"/"//trim(local_cnf_fname)
  endif

  call popsca( local_comm,     & ! [IN]
               local_cnf_fname ) ! [IN]

  ! stop MPI
  call comm_finalize

  stop
end program netcdf2grads_h_launcher
