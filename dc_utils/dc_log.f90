!-------------------------------------------------------------------------------
!> module DC_Log
!!
!! @par Description
!!          Handring Log
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2012-07-12 (S.Nishizawa) [new]
!!
!<
!-------------------------------------------------------------------------------
module dc_log
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: LogInit
  public :: LogFinalize
  public :: Log

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters
  !
  integer, parameter :: LOG_NONE  = 0
  integer, parameter :: LOG_ERROR = 1
  integer, parameter :: LOG_WARN  = 2
  integer, parameter :: LOG_INFO  = 3
  integer, parameter :: LOG_DEBUG = 4

  integer, parameter :: STDERR = 0
  integer, parameter :: STDOUT = 6

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  integer, public, parameter :: LOG_LMSG = 1024
  integer, public            :: LOG_fid  = STDOUT

  !-----------------------------------------------------------------------------
  !
  !++ Private variables
  !
  integer, private :: LOG_ilevel = LOG_INFO

  logical, private :: LOG_master = .true.
  logical, private :: LOG_opened = .false.

contains

  subroutine LogInit( &
       fid_conf, & ! (in)
       fid_log,  & ! (in) optional
       master    & ! (in) optional
       )
    implicit none
    integer, intent(in) :: fid_conf
    integer, intent(in), optional :: fid_log
    logical, intent(in), optional :: master

    character(len=5)   :: LOG_LEVEL = 'I'
    character(len=100) :: LOG_FILE = "LOG_"

    integer :: ierr

    namelist / PARAM_DC_LOG / &
         LOG_LEVEL, &
         LOG_FILE
  !-----------------------------------------------------------------------------

    if ( present(master) ) LOG_master = master

    call date_and_time(LOG_FILE(4:11), LOG_FILE(12:21))

    !--- read PARAM
    rewind(fid_conf)
    read(fid_conf,nml=PARAM_DC_LOG,iostat=ierr)

    if ( ierr > 0 ) then
       call Log('E', 'xxx Not appropriate names in namelist PARAM_DC_LOG. Check!')
    end if

    if ( present(fid_log) ) then
       LOG_FID = fid_log
    else
       open( LOG_FID,                & ! (out)
            file   = trim(LOG_FILE), & ! (in)
            form   = 'formatted',    & ! (in)
            iostat = ierr            ) ! (in)
       if ( ierr /= 0 ) then
          call Log('E', 'xxx File open error! :' // trim(LOG_FILE))
       end if
       LOG_opened = .true.
    end if

    select case (trim(LOG_LEVEL))
    case ('E', 'e', 'ERROR', 'error')
       LOG_ilevel = LOG_ERROR
    case ('W', 'w', 'WARN', 'warn')
       LOG_ilevel = LOG_WARN
    case ('I', 'i', 'INFO', 'info')
       LOG_ilevel = LOG_INFO
    case ('D', 'd', 'DEBUG', 'debug')
       LOG_ilevel = LOG_DEBUG
    case default
       call Log('E', 'xxx LOG_LEVEL is invalid. Check!')
    end select
       
    return
  end subroutine LogInit

  subroutine LogFinalize

    if ( LOG_opened ) close(LOG_fid)
    LOG_opened = .false.

  end subroutine LogFinalize

  subroutine Log( &
       type,   & ! (in)
       message & ! (in)
       )
    implicit none
    character(len=*), intent(in) :: type
    character(len=*), intent(in) :: message

    select case (trim(type))
    case ('E', 'e')
       if ( LOG_ilevel >= LOG_ERROR ) call LogPut(message)
       write(STDERR,*) message
       call abort
    case ('W', 'w')
       if ( LOG_ilevel >= LOG_WARN ) call LogPut(message)
    case ('I', 'i')
       if ( LOG_ilevel >= LOG_INFO ) call LogPut(message)
    case ('D', 'd')
       if ( LOG_ilevel >= LOG_DEBUG ) call LogPut(message)
    case default
       write(STDERR,*) 'BUG: wrong log level'
       call abort
    end select

    return
  end subroutine Log

!  private
  subroutine LogPut( &
       message & ! (in)
       )
    character(len=*) :: message

    if ( LOG_master ) write(LOG_fid, *) trim(message)

    return
  end subroutine LogPut

end module dc_log
