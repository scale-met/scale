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
     PRC_setup,     &
     PRC_MPIstart,  &
     PRC_MPIsplit,  &
     PRC_MPIfinish, &
     PRC_MPIstop,   &
     GLOBAL_nmax,   &
     GLOBAL_LOG,    &
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

  integer :: myrank_local             ! rank number in local communicator
  integer :: nmax_local               ! # of total process in local communicator
  integer :: nmax_parent              ! # of total process in parent world
  integer :: nmax_child               ! # of total process in child world
  integer :: icomm_parent             ! inter communicator with parent
  integer :: icomm_child              ! inter communicator with child
  integer :: NUM_DOMAIN               ! number of domains
  integer :: PRC_DOMAINS(max_depth)   ! # of total process in each domain
  integer :: PRC_ORDER(max_depth)     ! reordered number of process
  integer :: COLOR_DOMAINS(max_depth) ! # of color in each domain

  logical :: flag_parent              ! flag of "I am parent domain"
  logical :: flag_child               ! flag of "I am child domain"
  logical :: LOG_SPLIT = .false.      ! flag of log-output for mpi splitting

  character(len=H_LONG) :: CONF_FILES(max_depth)  ! names of configulation files
  character(len=H_LONG) :: fname_launch           ! config file for launcher
  character(len=H_LONG) :: fname_local            ! config file for local domain

  integer :: i, j, is, ie
  integer :: previous, shift
  integer :: ierr
  integer :: LNC_FID_CONF

  namelist / PARAM_LAUNCHER / &
     NUM_DOMAIN,     &
     PRC_DOMAINS,    &
     COLOR_DOMAINS,  &
     CONF_FILES,     &
     LOG_SPLIT
  !-----------------------------------------------------------

  NUM_DOMAIN       = 1
  PRC_DOMAINS      = 0
  CONF_FILES(:)    = ""
  COLOR_DOMAINS(:) = -1

  ! start MPI
  call PRC_MPIstart

  !--- read from argument
  if ( COMMAND_ARGUMENT_COUNT() < 1 ) then
     if (GLOBAL_LOG) write(*,*) 'WARNING: No config file specified!'
     if (GLOBAL_LOG) write(*,*) '         Default values are used.'
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
     if (GLOBAL_LOG) write(*,*)
     if (GLOBAL_LOG) write(*,*) 'WARNING: Failed to open config file! :', trim(fname_launch)
     if (GLOBAL_LOG) write(*,*) '         Default values are used.'
     if (GLOBAL_LOG) write(*,*)
  end if

  !--- read PARAM
  rewind(LNC_FID_CONF)
  read(LNC_FID_CONF,nml=PARAM_LAUNCHER,iostat=ierr)

  if( ierr > 0 ) then !--- fatal error
     if (GLOBAL_LOG) write(*,*) ' xxx Not appropriate names in namelist PARAM_LAUNCHER . Check!'
     call PRC_MPIstop
  endif

  close( LNC_FID_CONF )

  !--- make color order
  is = 1
  ie = NUM_DOMAIN
  previous = -999
  shift    = 1
  PRC_ORDER(:) = PRC_DOMAINS(:)
  call sort_ascd( PRC_ORDER(is:ie), is, ie )
print *, PRC_ORDER(is:ie)
print *, "-----------------"
print *, PRC_DOMAINS(is:ie)
print *, "-----------------"
  do i = is, ie
     do j = is, ie
print *, PRC_DOMAINS(i), PRC_ORDER(j)
        if ( PRC_DOMAINS(i) .eq. PRC_ORDER(j) ) then
print *, i, j
           COLOR_DOMAINS(i) = j - 1
           if ( previous .eq. COLOR_DOMAINS(i) ) then
              COLOR_DOMAINS(i) = COLOR_DOMAINS(i) + shift
              shift = shift + 1
           else
              previous = COLOR_DOMAINS(i)
              shift    = 1
           endif
           exit
        endif
     enddo
  enddo

do i = is, ie
if ( GLOBAL_LOG ) write (*,*) "[DEBUG] ", i, PRC_DOMAINS(i), COLOR_DOMAINS(i)
enddo
stop

  if ( NUM_DOMAIN == 1 ) then
     PRC_DOMAINS(1)   = GLOBAL_nmax
     COLOR_DOMAINS(1) = 0
  endif

  if ( GLOBAL_LOG ) write(*,*) '*** Start Launch System for SCALE-LES'
  if ( GLOBAL_LOG ) write ( *, * ) "NUM_DOMAIN = ", NUM_DOMAIN
  do i = 1, NUM_DOMAIN
     if ( GLOBAL_LOG ) write ( *, '(1X,A,I1,A,I5)' ) "PRC_DOMAINS(",i,") = ", PRC_DOMAINS(i)
     if ( GLOBAL_LOG ) write ( *, '(1X,A,I1,A,I3)' ) "COLOR_DOMAINS(",i,") = ", COLOR_DOMAINS(i)
     if ( GLOBAL_LOG ) write ( *, '(1X,A,I1,A,A)' ) "CONF_FILES(",i,")  = ", trim(CONF_FILES(i))
  enddo

  ! split MPI communicator
  call PRC_MPIsplit(    &
      NUM_DOMAIN,       & ! [in ]
      PRC_DOMAINS,      & ! [in ]
      COLOR_DOMAINS,    & ! [in ]
      CONF_FILES,       & ! [in ]
      LOG_SPLIT,        & ! [in ]
      nmax_parent,      & ! [out]
      nmax_child,       & ! [out]
      myrank_local,     & ! [out]
      nmax_local,       & ! [out]
      flag_parent,      & ! [out]
      flag_child,       & ! [out]
      icomm_parent,     & ! [out]
      icomm_child,      & ! [out]
      fname_local       ) ! [out]

  ! start main routine
  call scaleles ( myrank_local,    &
                  nmax_local,      &
                  flag_parent,     &
                  flag_child,      &
                  nmax_parent,     &
                  nmax_child,      &
                  icomm_parent,    &
                  icomm_child,     &
                  fname_local      )

  ! stop MPI
  call PRC_MPIfinish

  stop
  !=============================================================================
  !
  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> quicksort (ascending order)
  recursive subroutine sort_ascd(a, top, bottom)
    implicit none
    integer, intent(inout) :: a(:)
    integer, intent(in)    :: top, bottom
    integer :: i, j, cnt, trg
    !---------------------------------------------------------------------------
    cnt = a( (top+bottom) / 2 )
    i = top; j = bottom
    do
       do while ( a(i) > cnt ) !ascending evaluation
          i = i + 1
       enddo
       do while ( cnt > a(j) ) !ascending evaluation
          j = j - 1
       enddo
       if ( i >= j ) exit
       trg = a(i);  a(i) = a(j);  a(j) = trg
       i = i + 1
       j = j - 1
    enddo
    if ( top < i-1    ) call sort_ascd( a, top, i-1    )
    if ( j+1 < bottom ) call sort_ascd( a, j+1, bottom )
    return
  end subroutine sort_ascd
  !-----------------------------------------------------------------------------
  !
end program scaleles_launcher
