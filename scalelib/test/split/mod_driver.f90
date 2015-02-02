module mod_driver
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

  !--- Public
  !-----------------------------------------------------------
  public :: main

  !--- Private
  !-----------------------------------------------------------
  !
  !-----------------------------------------------------------
  contains
  !-----------------------------------------------------------


  !-----------------------------------------------------------
  subroutine main ( &
     myproc_main, &
     nprocs_main, &
     parent_prc,  &
     child_prc,   &
     fnamelist    )
    integer, intent(in) :: myproc_main
    integer, intent(in) :: nprocs_main
    integer, intent(in) :: parent_prc
    integer, intent(in) :: child_prc
    character(len=128), intent(in) :: fnamelist

  integer :: i, j, k, l, t
  integer :: value
  integer :: ierr

  namelist / PARAM_INTERCOMM / &
     NX,             &
     NY,             &
     NZ,             &
     NT,             &
     DOMAIN_ID,      &
     IAM_PARENT,     &
     IAM_CHILD
  !-----------------------------------------------------------

  ! default settings
  NX         = 8
  NY         = 8
  NZ         = 40
  NT         = 5
  DOMAIN_ID  = 0
  IAM_PARENT = .false.
  IAM_CHILD  = .false.

  ! initialize local communication
  myproc = myproc_main
  nprocs = nprocs_main

  value      = 1

  ! read namelist
  open  ( FID_NML, file=trim(fnamelist), status='old', delim='apostrophe' )
  read  ( FID_NML, nml=PARAM_INTERCOMM )
  close ( FID_NML )

  call prof_initialize
  call prof_start ( 'init' )

  ! open log file
  LFID = DOMAIN_ID*1000 + myproc
  write(domnum,'(I2.2)') DOMAIN_ID
  write(prcnum,'(I3.3)') myproc
  flogfile = "LOG.d"//domnum//".pe"//prcnum
  open  ( LFID, file=trim(flogfile), status='replace' )

  write ( LFID, * ) "---------------------------------------------------"
  write ( LFID, * ) "     START MINI TEST OF INTER COMMUNICATION"
  write ( LFID, * ) "---------------------------------------------------"
  write ( LFID, * ) ""
  write ( LFID, * ) "NAMELIST: DOMAIN_ID = ", DOMAIN_ID
  write ( LFID, * ) "NAMELIST: NX = ", NX
  write ( LFID, * ) "NAMELIST: NY = ", NY
  write ( LFID, * ) "NAMELIST: NZ = ", NZ
  write ( LFID, * ) "NAMELIST: NT = ", NT
  write ( LFID, * ) "NAMELIST: IAM_PARENT = ", IAM_PARENT
  write ( LFID, * ) "NAMELIST: IAM_CHILD = ", IAM_CHILD
  write ( LFID, * ) ""
  write ( LFID, * ) "---------------------------------------------------"
  write ( LFID, * ) "INITIALIZE: TOTAL PROCESSES = ", nprocs
  write ( LFID, * ) "INITIALIZE: MY PROCESS NUMBER = ", myproc
  write ( LFID, * ) "INITIALIZE: MY PROCESS NUMBER in GLOBAL = ", myproc_glb
  if ( IAM_PARENT ) then
     write ( LFID, * ) "INITIALIZE: CHILD PRC = ", child_prc
     write ( LFID, * ) "INITIALIZE: INTERCOMM_CHILD = ", INTERCOMM_CHILD
  endif
  if ( IAM_CHILD ) then
     write ( LFID, * ) "INITIALIZE: PARENT PRC = ", parent_prc
     write ( LFID, * ) "INITIALIZE: INTERCOMM_PARENT = ", INTERCOMM_PARENT
  endif

  ! memory allocate
  allocate ( var  ( NX, NY, NZ ) )
  allocate ( dummy( 1,  1,  1  ) )

  ! initialize variables
  if ( IAM_PARENT ) then
     do k = 1, nz
     do j = 1, ny
     do i = 1, nx
        var(i,j,k) = float( myproc+1 )*1000 + float( DOMAIN_ID )*100 + float(value)
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
     call comm_makeyp ( nprocs, child_prc, parent )
     call prof_end ( 'make yp' )
  endif
  if ( IAM_CHILD ) then
     call prof_start ( 'make yp' )
     call comm_makeyp ( parent_prc, nprocs, child )
     call prof_end ( 'make yp' )
  endif


  sleep_loop = DOMAIN_ID

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
!     call sleep(1)
  enddo

  ! update a value
  if ( IAM_PARENT ) then
     value = value + 1
     do k = 1, nz
     do j = 1, ny
     do i = 1, nx
        var(i,j,k) = float( myproc+1 )*1000 + float( DOMAIN_ID )*100 + float(value)
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


  !call MPI_FINALIZE(ierr)
  write ( LFID, * ) ""
  write ( LFID, * ) "---------------------------------------------------"
  write ( LFID, * ) "main routine finished"
  close ( LFID )
  return

  end subroutine main
  !-----------------------------------------------------------

  !-----------------------------------------------------------
end module mod_driver
