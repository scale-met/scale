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
  public :: Log_nml

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
  integer, public, parameter :: LOG_LMSG    = 4096
  integer, public            :: LOG_fid
  integer, public            :: LOG_fid_nml

  !-----------------------------------------------------------------------------
  !
  !++ Private variables
  !
  integer, private :: LOG_ilevel = LOG_INFO

  logical, private :: LOG_master
#if defined(PGI) || defined(SX)
  logical, public :: LOG_master_nml
#else
  logical, private :: LOG_master_nml
#endif
  logical, private :: LOG_opened     = .false.

contains

  subroutine LogInit( &
       fid_conf,  & ! (in)
       fid_log,   & ! (in) optional
       master,    & ! (in) optional
       fid_nml,   & ! (in) optional
       master_nml & ! (in) optional
       )
    implicit none
    integer, intent(in) :: fid_conf
    integer, intent(in), optional :: fid_log
    logical, intent(in), optional :: master
    integer, intent(in), optional :: fid_nml
    logical, intent(in), optional :: master_nml

    character(len=5)   :: LOG_LEVEL = 'I'
    character(len=100) :: LOG_FILE  = "LOG_"

    integer :: ierr

    namelist / PARAM_DC_LOG / &
         LOG_LEVEL, &
         LOG_FILE

    character(len=LOG_LMSG) :: message
  !-----------------------------------------------------------------------------

    LOG_fid = STDOUT
    LOG_fid_nml = STDOUT
    LOG_master = .true.

    if ( present(master)     ) LOG_master     = master
    LOG_master_nml = LOG_master
    if ( present(master_nml) ) LOG_master_nml = master_nml

    call date_and_time(LOG_FILE(4:11), LOG_FILE(12:21))

    !--- read PARAM
    rewind(fid_conf)
    read(fid_conf,nml=PARAM_DC_LOG,iostat=ierr)

    if ( ierr > 0 ) then
       call Log('E', 'xxx Not appropriate names in namelist PARAM_DC_LOG. Check!')
    end if

    if ( present(fid_log) ) then
       LOG_FID  = fid_log
       LOG_FILE = ''
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

    LOG_fid_nml = LOG_FID
    if ( present(fid_nml) ) then
       LOG_fid_nml = fid_nml
    end if

#if defined(PGI) || defined(SX)
    if ( LOG_master_nml ) write(LOG_fid_nml,nml=PARAM_DC_LOG)
#else
    write(message,nml=PARAM_DC_LOG)
    call Log_nml('I',message)
#endif

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
       write(STDERR,*) trim(message)
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

  subroutine Log_nml( &
       type,   & ! (in)
       message & ! (in)
       )
    implicit none
    character(len=*), intent(in) :: type
    character(len=*), intent(in) :: message

    select case (trim(type))
    case ('E', 'e')
       if ( LOG_ilevel >= LOG_ERROR ) call LogPut_nml(message)
       write(STDERR,*) trim(message)
       call abort
    case ('W', 'w')
       if ( LOG_ilevel >= LOG_WARN ) call LogPut_nml(message)
    case ('I', 'i')
       if ( LOG_ilevel >= LOG_INFO ) call LogPut_nml(message)
    case ('D', 'd')
       if ( LOG_ilevel >= LOG_DEBUG ) call LogPut_nml(message)
    case default
       write(STDERR,*) 'BUG: wrong log level'
       call abort
    end select

    return
  end subroutine Log_nml

!  private
  subroutine LogPut_nml( &
       message & ! (in)
       )
    character(len=*) :: message

    if ( LOG_master_nml ) write(LOG_fid_nml, *) trim(message)

    return
  end subroutine LogPut_nml

end module dc_log
