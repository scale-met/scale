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
  integer :: F_PRC_DOMAINS(max_depth)   ! # of total process in each domain (final)
  integer :: PRC_ORDER(max_depth)     ! reordered number of process
  integer :: COLOR_DOMAINS(max_depth) ! # of color in each domain
  integer :: F_COLOR_DOMAINS(max_depth) ! # of color in each domain (final)
  integer :: PARENT_COLOR(max_depth)  ! parent color number
  integer :: F_PARENT_COLOR(max_depth)  ! parent color number (final)
  integer :: CHILD_COLOR(max_depth)   ! child color number
  integer :: F_CHILD_COLOR(max_depth)   ! child color number (final)
  integer :: PARENT_PRC(max_depth)  ! parent color number
  integer :: CHILD_PRC(max_depth)   ! child color number
  integer :: F_PARENT_PRC(max_depth)  ! parent color number (final)
  integer :: F_CHILD_PRC(max_depth)   ! child color number (final)
  integer :: REL_PARENT_COLOR(max_depth)  ! parent color number (relationship)
  integer :: REL_CHILD_COLOR(max_depth)   ! child color number (relationship)
  integer :: REL_PARENT_PRC(max_depth)  ! parent color number (relationship)
  integer :: REL_CHILD_PRC(max_depth)   ! child color number (relationship)
  integer :: REF_COL2DOM(0:max_depth) ! refering from # of color to # of domain
  integer :: F_REF_COL2DOM(0:max_depth) ! refering from # of color to # of domain
  integer :: REF_DOMAIN(max_depth) ! refering from # of color to # of domain
  integer :: REF_DOM2ORDER(max_depth) ! refering from # of domain to # of order

  logical :: HAVE_PARENT(max_depth) = .false.      ! flag of have parent
  logical :: HAVE_CHILD(max_depth) = .false.      ! flag of have child

  logical :: flag_parent              ! flag of "I am parent domain"
  logical :: flag_child               ! flag of "I am child domain"
  logical :: LOG_SPLIT = .false.      ! flag of log-output for mpi splitting

  character(len=H_LONG) :: CONF_FILES(max_depth)  ! names of configulation files
  character(len=H_LONG) :: F_CONF_FILES(max_depth)  ! names of configulation files (final)
  character(len=H_LONG) :: fname_launch           ! config file for launcher
  character(len=H_LONG) :: fname_local            ! config file for local domain

  integer :: dnum
  integer :: id_parent  ! parent domain number
  integer :: id_child   ! child domain number

  integer :: i, j, is, ie
  integer :: touch(max_depth)
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
  F_CONF_FILES(:)  = ""
  COLOR_DOMAINS(:) = -1
  F_COLOR_DOMAINS(:) = -1
  PARENT_COLOR(:) = -1
  F_PARENT_COLOR(:) = -1
  CHILD_COLOR(:) = -1
  F_CHILD_COLOR(:) = -1
  PARENT_PRC(:) = -1
  CHILD_PRC(:) = -1
  F_PARENT_PRC(:) = -1
  F_CHILD_PRC(:) = -1

     REL_PARENT_COLOR(:) = -1
     REL_CHILD_COLOR(:)  = -1
     REL_PARENT_PRC(:)   = -1
     REL_CHILD_PRC(:)    = -1

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
  !    domain num is counted from 1
  !    color num  is counted from 0
  if ( NUM_DOMAIN > 1 ) then ! multiple domain case
  is = 1
  ie = NUM_DOMAIN
  touch(:) = -1
  PRC_ORDER(:) = PRC_DOMAINS(:)
  call sort_ascd( PRC_ORDER(is:ie), is, ie )
  do i = is, ie
  do j = is, ie
     if ( PRC_DOMAINS(i) .eq. PRC_ORDER(j) .and. touch(j) < 0 ) then
        COLOR_DOMAINS(i) = j - 1            ! domain_num --> color_num
        REF_COL2DOM(COLOR_DOMAINS(i)) = i   ! color_num  --> domain_num
        touch(j) = 1
        exit
     endif
  enddo
  enddo
  do i = is, ie
     id_parent = i - 1
     id_child  = i + 1
     if ( 1 <= id_parent .and. id_parent <= NUM_DOMAIN ) then
        PARENT_COLOR(i) = COLOR_DOMAINS(id_parent)
        PARENT_PRC(i)   = PRC_DOMAINS(id_parent)
     else
        PARENT_COLOR(i) = -1
        PARENT_PRC(i)   = -1
     endif
     if ( 1 <= id_child  .and. id_child  <= NUM_DOMAIN ) then
        CHILD_COLOR(i) = COLOR_DOMAINS(id_child)
        CHILD_PRC(i)   = PRC_DOMAINS(id_child)
     else
        CHILD_COLOR(i) = -1
        CHILD_PRC(i)   = -1
     endif
     if ( GLOBAL_LOG .and. LOG_SPLIT ) write( *, '(1X,A,I2,1X,A,I2,2(2X,A,I2,1X,A,I5,A))' ) &
                                       "DOMAIN: ", i, "MY_COL: ", COLOR_DOMAINS(i), &
                                       "PARENT: COL= ", PARENT_COLOR(i), "PRC(", PARENT_PRC(i), ")", &
                                       "CHILD: COL= ", CHILD_COLOR(i),  "PRC(", CHILD_PRC(i), ")"
  enddo

  !--- reorder following color order
  do i = is, ie
     dnum = REF_COL2DOM(i-1)
     REF_DOMAIN(i)      = dnum
     REF_DOM2ORDER(dnum) = i
     F_PRC_DOMAINS(i)   = PRC_DOMAINS(dnum)
     F_COLOR_DOMAINS(i) = COLOR_DOMAINS(dnum)
     F_CONF_FILES(i)    = CONF_FILES(dnum)
     F_PARENT_COLOR(i)  = PARENT_COLOR(dnum)
     F_CHILD_COLOR(i)   = CHILD_COLOR(dnum)
     F_PARENT_PRC(i)    = PARENT_PRC(dnum)
     F_CHILD_PRC(i)     = CHILD_PRC(dnum)
!     F_REF_COL2DOM(F_COLOR_DOMAINS(i)) = REF_DOMAIN(i)
!     if ( REF_DOMAIN(i) /= NUM_DOMAIN ) then
!        HAVE_PARENT(i) = .true.
!     endif
!     if ( REF_DOMAIN(i) /= 1) then
!        HAVE_CHILD(i) =  .true.
!     endif
!if ( GLOBAL_LOG ) print *, "REF_DOMAIN", i, REF_DOMAIN(i)
  enddo

  !--- set relationship
  HAVE_PARENT(:) = .false.
  HAVE_CHILD(:)  = .false.
  do i = is, ie-1
     id_parent = REF_DOM2ORDER(i)
     id_child  = REF_DOM2ORDER(i+1)
     REL_PARENT_COLOR(i) = F_PARENT_COLOR(id_child)
     REL_CHILD_COLOR(i)  = F_CHILD_COLOR(id_parent)
     REL_PARENT_PRC(i)   = F_PARENT_PRC(id_child)
     REL_CHILD_PRC(i)    = F_CHILD_PRC(id_parent)
     HAVE_PARENT(id_child) = .true.
     HAVE_CHILD(id_parent) = .true.
!if ( GLOBAL_LOG ) print *, i, id_parent, id_child
!if ( GLOBAL_LOG ) print *, i, REL_PARENT_COLOR(i), REL_CHILD_COLOR(i), REL_PARENT_PRC(i), REL_CHILD_PRC(i)
  enddo

!if ( GLOBAL_LOG ) print *, HAVE_PARENT(:)
!if ( GLOBAL_LOG ) print *, HAVE_CHILD(:)
!do i = is, ie
!if ( GLOBAL_LOG ) write (*,*) "[DEBUG] ", i, PRC_DOMAINS(i), COLOR_DOMAINS(i)
!enddo
!stop

  else ! single domain case
     F_PRC_DOMAINS(1)   = GLOBAL_nmax
     F_COLOR_DOMAINS(1) = 0
     F_CONF_FILES(1)    = CONF_FILES(1)
     HAVE_PARENT(:) = .false.
     HAVE_CHILD(:)  = .false.
  endif

  if ( GLOBAL_LOG ) write(*,*) '*** Start Launch System for SCALE-LES'
  if ( GLOBAL_LOG ) write ( *, '(1X,A,I2)' ) "TOTAL DOMAIN NUMBER = ", NUM_DOMAIN
  do i = 1, NUM_DOMAIN
     if ( GLOBAL_LOG ) write ( *, * ) ""
     if ( GLOBAL_LOG ) write ( *, '(1X,A,I2,A,I5)' ) "ORDER (",i,") -> DOMAIN: ", REF_DOMAIN(i)
     if ( GLOBAL_LOG ) write ( *, '(1X,A,I1,A,I5)' ) "NUM PRC_DOMAINS(",i,")  = ", F_PRC_DOMAINS(i)
     if ( GLOBAL_LOG ) write ( *, '(1X,A,I1,A,I3)' ) "MY COLOR_DOMAINS(",i,") = ", F_COLOR_DOMAINS(i)
!     if ( GLOBAL_LOG ) write ( *, '(1X,A,I1,A,I3)' ) "PARENT_COLOR(",i,")  = ", F_PARENT_COLOR(i)
!     if ( GLOBAL_LOG ) write ( *, '(1X,A,I1,A,I3)' ) "CHILD_COLOR(",i,")   = ", F_CHILD_COLOR(i)
!     if ( GLOBAL_LOG ) write ( *, '(1X,A,I1,A,I3)' ) "REF_COL2DOM(",COLOR_DOMAINS(i),")  = ", REF_COL2DOM(COLOR_DOMAINS(i))
     if ( GLOBAL_LOG ) write ( *, '(1X,A,I1,A,A)'  ) "CONF_FILES(",i,")    = ", trim(F_CONF_FILES(i))
  enddo
  if ( GLOBAL_LOG ) write ( *, * ) ""

do i = 1, NUM_DOMAIN-1
 if ( GLOBAL_LOG ) write ( *, '(1X,A,I2)' ) "relationship: ", i
 if ( GLOBAL_LOG ) write ( *, '(1X,A,I2,A,I2)' ) &
                   "--- parent color = ", REL_PARENT_COLOR(i), "  child color = ", REL_CHILD_COLOR(i)
 if ( GLOBAL_LOG ) write ( *, '(1X,A,I5,A,I5)' ) &
                   "--- parent prc = ", REL_PARENT_PRC(i), "  child prc = ", REL_CHILD_PRC(i)
enddo


  ! split MPI communicator
  call PRC_MPIsplit(    &
      NUM_DOMAIN,       & ! [in ]
      F_PRC_DOMAINS,      & ! [in ]
      F_COLOR_DOMAINS,    & ! [in ]
      REL_PARENT_COLOR,    & ! [in ]
      REL_CHILD_COLOR,    & ! [in ]
      REL_PARENT_PRC,    & ! [in ]
      REL_CHILD_PRC,    & ! [in ]
      HAVE_PARENT,    & ! [in ]
      HAVE_CHILD,    & ! [in ]
!!!!      REF_COL2DOM,      & ! [in ]
      F_CONF_FILES,       & ! [in ]
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

print *, "call-scale-les", myrank_local, flag_parent, flag_child, nmax_parent, nmax_child, icomm_parent, icomm_child, fname_local

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
