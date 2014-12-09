program driver
  !
  !-----------------------------------------------------------
  ! 2014/06/13: Original (Ryuji Yoshida)
  !
  !-----------------------------------------------------------

  use mpi !external
  use mod_vars
  use mod_prof
  use mod_comm

  implicit none

  integer :: i, j, k, l, t
  integer :: value
  integer :: ierr

  namelist / PARAM_INTERCOMM / &
     NX,             &
     NY,             &
     NZ,             &
     NT,             &
     DOMAIN_NUM,     &
     IAM_PARENT,     &
     IAM_CHILD,      &
     LAUNCH_PRC
  !-----------------------------------------------------------

  call get_command_argument(1,fnamelist)

  ! default settings
  NX         = 8
  NY         = 8
  NZ         = 40
  NT         = 5
  DOMAIN_NUM = 0
  LAUNCH_PRC = 0
  IAM_PARENT = .false.
  IAM_CHILD  = .false.

  value      = 1

  ! read namelist
  open  ( FID_NML, file=trim(fnamelist), status='old', delim='apostrophe' )
  read  ( FID_NML, nml=PARAM_INTERCOMM )
  close ( FID_NML )

  ! initialize MPI
  call MPI_INIT( ierr )
  call MPI_COMM_RANK( MPI_COMM_WORLD, myproc, ierr )
  call MPI_COMM_SIZE( MPI_COMM_WORLD, nprocs, ierr )
  call prof_initialize
  call prof_start ( 'init' )

  ! open log file
  LFID = DOMAIN_NUM*1000 + myproc
  write(domnum,'(I2.2)') DOMAIN_NUM
  write(prcnum,'(I3.3)') myproc
  flogfile = "LOG.d"//domnum//".pe"//prcnum
  open  ( LFID, file=trim(flogfile), status='replace' )

  write ( LFID, * ) "---------------------------------------------------"
  write ( LFID, * ) "     START MINI TEST OF INTER COMMUNICATION"
  write ( LFID, * ) "---------------------------------------------------"
  write ( LFID, * ) ""
  write ( LFID, * ) "NAMELIST: DOMAIN_NUM = ", DOMAIN_NUM
  write ( LFID, * ) "NAMELIST: NX = ", NX
  write ( LFID, * ) "NAMELIST: NY = ", NY
  write ( LFID, * ) "NAMELIST: NZ = ", NZ
  write ( LFID, * ) "NAMELIST: NT = ", NT
  write ( LFID, * ) "NAMELIST: IAM_PARENT = ", IAM_PARENT
  write ( LFID, * ) "NAMELIST: IAM_CHILD = ", IAM_CHILD
  if ( IAM_PARENT ) then
     write ( LFID, * ) "NAMELIST: LAUNCH_PRC = ", LAUNCH_PRC
  endif
  write ( LFID, * ) ""
  write ( LFID, * ) "---------------------------------------------------"
  write ( LFID, * ) "MPI INITIALIZE: TOTAL PROCESSES = ", nprocs
  write ( LFID, * ) "MPI INITIALIZE: MY PROCESS NUMBER = ", myproc

  ! memory allocate
  allocate ( var  ( NX, NY, NZ ) )
  allocate ( dummy( 1,  1,  1  ) )

  ! launch child process
  write ( LFID, * ) ""
  write ( LFID, * ) "---------------------------------------------------"
  call prof_start ( 'launch' )
  if ( IAM_PARENT ) then
     call comm_spawn ( DOMAIN_NUM, LAUNCH_PRC )
  endif
  if ( IAM_CHILD ) then
     call comm_get_parent
  endif
  call prof_end ( 'launch' )

  ! initialize variables
  if ( IAM_PARENT ) then
     do k = 1, nz
     do j = 1, ny
     do i = 1, nx
        var(i,j,k) = float( myproc+1 )*1000 + float( DOMAIN_NUM )*100 + float(value)
     enddo
     enddo
     enddo
  else
     var(:,:,:) = 0.0D0
  endif

  ! setup inter-domain communication
  write ( LFID, * ) ""
  write ( LFID, * ) "---------------------------------------------------"
  if ( IAM_PARENT ) then
     call prof_start ( 'make yp' )
     call comm_makeyp ( nprocs, LAUNCH_PRC, parent )
     call prof_end ( 'make yp' )
  endif
  if ( IAM_CHILD ) then
     call prof_start ( 'make yp' )
     call comm_makeyp ( PARENT_SIZE, nprocs, child )
     call prof_end ( 'make yp' )
  endif


  sleep_loop = DOMAIN_NUM

  if ( IAM_PARENT ) then
     call prof_start ( 'barrier' )
     call comm_barrier ( parent )
     call prof_end ( 'barrier' )
  endif
  if ( IAM_CHILD ) then
     call prof_start ( 'barrier' )
     call comm_barrier ( child )
     call prof_end ( 'barrier' )
  endif

  call prof_end ( 'init' )
  call prof_start ( 'main' )


  !----------------------------------------------------------->>>
  do t = 1, NT
  write ( LFID, * ) ""
  write ( LFID, '(1X,A,I2,A,I2)' ) "time step: ", t, "/", NT
  write ( LFID, * ) "---------------------------------------------------"
  write ( LFID, '(1X,A,F9.3,A,F9.3)' ) &
        "before comm; var check: max = ", maxval(var), "  min = ", minval(var)

  ! issue send or recv
  if ( IAM_PARENT ) then
     call prof_start ( 'send' )
     call comm_send ( var, parent )
     call prof_end ( 'send' )
  endif
  if ( IAM_CHILD ) then
     call prof_start ( 'recv' )
     call comm_recv ( child )
     call prof_end ( 'recv' )
  endif

  ! dummy main routine
  do l = 1, sleep_loop
     write ( LFID, '(3X,A,I2,A,I2)' ) "DO DUMMY JOB: ", l, "/", sleep_loop
     call sleep(1)
  enddo

  ! update a value
  if ( IAM_PARENT ) then
     value = value + 1
     do k = 1, nz
     do j = 1, ny
     do i = 1, nx
        var(i,j,k) = float( myproc+1 )*1000 + float( DOMAIN_NUM )*100 + float(value)
     enddo
     enddo
     enddo
  endif

  ! issue wait
  if ( IAM_PARENT ) then
     call prof_start ( 'wait' )
     call comm_wait ( dummy, parent )
     call prof_end ( 'wait' )
  endif
  if ( IAM_CHILD ) then
     call prof_start ( 'wait' )
     call comm_wait ( var, child )
     call prof_end ( 'wait' )
  endif

  ! check variables
  write ( LFID, '(1X,A,F9.3,A,F9.3)' ) &
        "after comm; var check: max = ", maxval(var), "  min = ", minval(var)

  enddo
  !-----------------------------------------------------------<<<

  call prof_end ( 'main' )
  call prof_start ( 'finaliz' )

  ! finalize the process
  write ( LFID, * ) ""
  write ( LFID, * ) "---------------------------------------------------"
  write ( LFID, * ) "finalizing the process"

  if ( IAM_PARENT ) then
     call prof_start ( 'barrier' )
     call comm_barrier ( parent )
     call prof_end ( 'barrier' )
     call prof_start ( 'free' )
     call comm_free ( parent )
     call prof_end ( 'free' )
  endif
  if ( IAM_CHILD ) then
     call prof_start ( 'barrier' )
     call comm_barrier ( child )
     call prof_end ( 'barrier' )
     call prof_start ( 'free' )
     call comm_free ( child )
     call prof_end ( 'free' )
  endif
  call prof_end ( 'finaliz' )

  ! report profile
  write ( LFID, * ) ""
  write ( LFID, * ) "---------------------------------------------------"
  call prof_finalize


  call MPI_FINALIZE(ierr)
  write ( LFID, * ) ""
  write ( LFID, * ) "---------------------------------------------------"
  write ( LFID, * ) "finished"
  close ( LFID )
  stop

  !-----------------------------------------------------------
end program driver
