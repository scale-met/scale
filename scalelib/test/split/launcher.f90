program launcher
  !
  !-----------------------------------------------------------
  ! 2014/12/09: Original (Ryuji Yoshida)
  !
  !-----------------------------------------------------------

  use mpi !external
  use mod_vars
  use mod_prof
  use mod_comm
  use mod_driver

  implicit none

  integer :: i
  integer :: total_prc, nprc, color, key
  integer :: myproc_main, nprocs_main
  integer :: parent_prc, child_prc
  integer :: itag, ierr
  integer :: newerr, stat
  logical :: do_create_p(max_domain)
  logical :: do_create_c(max_domain)
  character(len=128) :: fname_main = "run.conf"

  namelist / PARAM_LAUNCHER / &
     NUM_DOMAIN,     &
     PRC_DOMAINS,    &
     CONF_FILES
  !-----------------------------------------------------------

  call get_command_argument(1,fnamelist)

  ! default settings
  NUM_DOMAIN  = 1
  PRC_DOMAINS = 0
  CONF_FILES  = "run.d01.conf"

  ! read namelist
  open  ( FID_NML, file=trim(fnamelist), status='old', delim='apostrophe' )
  read  ( FID_NML, nml=PARAM_LAUNCHER )
  close ( FID_NML )

  ! initialize MPI
  call MPI_INIT( ierr )
  call MPI_COMM_RANK( MPI_COMM_WORLD, myproc_glb, ierr )
  call MPI_COMM_SIZE( MPI_COMM_WORLD, nprocs_glb, ierr )
  GLB_COMM_WORLD = MPI_COMM_WORLD

  if ( myproc_glb == master ) then
     masternode = .true.
  else
     masternode = .false.
  endif

  ! open log file
  if ( masternode ) then
     MFID = 999
     flogfile = "LOG.launcher.masternode"
     open  ( MFID, file=trim(flogfile), status='replace' )

     write ( MFID, * ) "---------------------------------------------------"
     write ( MFID, * ) "     LAUNCH MINI TEST OF INTER COMMUNICATION"
     write ( MFID, * ) "---------------------------------------------------"
     write ( MFID, * ) ""
     write ( MFID, * ) "NAMELIST: NUM_DOMAIN = ", NUM_DOMAIN
     do i = 1, NUM_DOMAIN
        write ( MFID, '(1X,A,I1,A,I3)' ) "NAMELIST: PRC_DOMAINS(",i,") = ", PRC_DOMAINS(i)
        write ( MFID, '(1X,A,I1,A,A)' ) "NAMELIST: CONF_FILES(",i,")  = ", trim(CONF_FILES(i))
     enddo
  endif

  ! communicator management
  if ( masternode ) write ( MFID, * ) ""
  if ( masternode ) write ( MFID, * ) "---------------------------------------------------"

  if ( NUM_DOMAIN == 1 ) then ! single domain run
     LOC_COMM_WORLD = MPI_COMM_WORLD
     parent_prc = 0
     child_prc  = 0
     myproc_main = myproc_glb
     nprocs_main = nprocs_glb
     fname_main  = CONF_FILES(1)

  elseif ( NUM_DOMAIN > 1 ) then ! multi domain run

     total_prc = 0
     do i = 1, NUM_DOMAIN
        total_prc = total_prc + PRC_DOMAINS(i)
     enddo
     if ( total_prc /= nprocs_glb ) then
        if ( masternode ) write ( MFID, * ) ""
        if ( masternode ) write ( MFID, * ) "ERROR: PROCESS NUMBER is INCONSISTENT"
        if ( masternode ) write ( MFID, * ) " REQ. TOTAL PRC = ", total_prc, "  NPROCS = ", nprocs_glb
        call comm_abort
     endif

     allocate ( COLOR_LIST(0:nprocs_glb-1) )
     allocate ( KEY_LIST  (0:nprocs_glb-1) )

     ! make a process table
     color = 1
     key   = 0
     nprc  = PRC_DOMAINS(color) - 1
     PRC_ROOT(:) = -999
     do i = 0, nprocs_glb-1
        if ( key == 0 ) PRC_ROOT(color) = i
        COLOR_LIST(i) = color
        KEY_LIST(i)   = key
        if ( masternode ) write ( MFID, '(1X,4(A,I3))' ) "PE:", i, "   COLOR:", COLOR_LIST(i), &
                                "   KEY:", KEY_LIST(i), "   PRC_ROOT:", PRC_ROOT(color)
        key = key + 1
        if ( key >= nprc ) then
           color = color + 1
           key   = 0
           nprc  = PRC_DOMAINS(color)
        endif
     enddo

     ! split comm_world
     my_color = COLOR_LIST(myproc_glb) ! equal to domain number
     my_key = KEY_LIST(myproc_glb)     ! equal to process number in the local communicator
     call MPI_COMM_SPLIT(MPI_COMM_WORLD, my_color, &
                         my_key, LOC_COMM_WORLD, ierr)

     parent_prc = 0
     child_prc  = 0
     myproc_main = my_key
     nprocs_main = PRC_DOMAINS(my_color)
     fname_main  = CONF_FILES(my_color)

     ! sync communicator
     if ( masternode ) write ( MFID, * ) ""
     if ( masternode ) write ( MFID, * ) "---------------------------------------------------"


     do_create_p(:) = .false.
     do_create_c(:) = .false.

     ! set parent-child relationship
     do i = 1, NUM_DOMAIN-1
        if ( masternode ) write ( MFID, '(1X,A,I2)' ) "relationship: ", i
        if ( masternode ) write ( MFID, '(1X,A,I2,A,I2)' ) "--- parent color = ", i, "  child color = ", i+1
        if ( my_color == i ) then
           flag_parent  = .true.
           do_create_p(i) = .true.
           child_prc   = PRC_DOMAINS(my_color+1)
        elseif ( my_color == i+1 ) then
           flag_child   = .true.
           do_create_c(i) = .true.
           parent_prc  = PRC_DOMAINS(my_color-1)
        endif
     enddo

     ! create inter-commnunicator
     do i = 1, NUM_DOMAIN-1
        itag = i*100
        if ( do_create_p(i) ) then ! as a parent
              call MPI_INTERCOMM_CREATE(LOC_COMM_WORLD, root,                 &
                                        MPI_COMM_WORLD, PRC_ROOT(my_color+1), &
                                        itag, INTERCOMM_CHILD, ierr)
        elseif ( do_create_c(i) ) then ! as a child
              call MPI_INTERCOMM_CREATE(LOC_COMM_WORLD, root,                 &
                                        MPI_COMM_WORLD, PRC_ROOT(my_color-1), &
                                        itag, INTERCOMM_PARENT, ierr)
        endif

        call MPI_BARRIER(MPI_COMM_WORLD, ierr)
     enddo

  else
     if ( masternode ) write ( MFID, * ) ""
     if ( masternode ) write ( MFID, * ) "ERROR: REQUESTED DOMAIN NUMBER IS NOT ACCEPTABLE"
     call comm_abort
  endif


  ! start main routine
  if ( masternode ) write ( MFID, * ) ""
  if ( masternode ) write ( MFID, * ) "---------------------------------------------------"
  if ( masternode ) write ( MFID, * ) "---> go into main routine"

  call main ( myproc_main,  &
              nprocs_main,  &
              parent_prc,   &
              child_prc,    &
              fname_main    )

  if ( masternode ) write ( MFID, * ) ""
  if ( masternode ) write ( MFID, * ) "<--- return from main routine"


  ! finalizing
  if ( masternode ) write ( MFID, * ) "---------------------------------------------------"
  if ( masternode ) write ( MFID, * ) "free the splitted comm world"

  call MPI_COMM_FREE(LOC_COMM_WORLD, ierr)
  call MPI_FINALIZE(ierr)

  if ( masternode ) write ( MFID, * ) ""
  if ( masternode ) write ( MFID, * ) "finished"

  close ( MFID )
  stop

  !-----------------------------------------------------------
end program launcher
