!-------------------------------------------------------------------------------
!> module History
!!
!! @par Description
!!          History output
!!
!! @author H.Tomita and SCALE developpers
!!
!! @par History
!! @li      2011-11-11 (H.Yashiro) [new] Imported from SCALE-LES ver.2
!!
!<
!-------------------------------------------------------------------------------
module mod_history
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedures
  !
  public :: HISTRY_setup
  public :: HISTRY_put
  public :: HISTRY_write
  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private procedures
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !

  character(len=IO_FILECHR), private, save :: HISTORY_BASENAME      = 'history'

  integer, private, parameter :: HIST_limit = 1000 !> number limit for history item

  integer,                   private, allocatable, save :: HIST_nmax
  character(len=FIO_HSHORT), private, allocatable, save :: HIST_item(HIST_limit)
  character(len=FIO_HMID),   private, allocatable, save :: HIST_desc(HIST_limit)
  character(len=FIO_HSHORT), private, allocatable, save :: HIST_unit(HIST_limit)
  real(8) ,                  private, allocatable, save :: HIST_kmax(HIST_limit)

  !-----------------------------------------------------------------------------  
contains
  !-----------------------------------------------------------------------------
  subroutine HISTRY_setup
    use mod_stdio, only: &
       IO_FID_CONF, &
       IO_FID_LOG,  &
       IO_L
    use mod_process, only: &
       PRC_MPIstop
    use mod_time, only: &
       TIME_ymdhms2sec
    implicit none

    real(8)                  :: HISTORY_DEFAULT_TINTERVAL = 1.D0
    real(8)                  :: HISTORY_DEFAULT_TAVERAGE  = 1.D0
    character(len=IO_SYSCHR) :: HISTORY_DEFAULT_TUNIT     = "MIN"
    character(len=IO_SYSCHR) :: HISTORY_DATATYPE          = "REAL8"

    NAMELIST / PARAM_HISTORY / &
       HISTORY_BASENAME,          &
       HISTORY_DEFAULT_TINTERVAL, &
       HISTORY_DEFAULT_TAVERAGE,  &
       HISTORY_DEFAULT_TUNIT,     &
       HISTORY_DATATYPE

    character(len=FIO_HSHORT) :: ITEM  !> name of history item
    real(8),                  :: TINT  !> time interval to output
    real(8),                  :: TAVG  !> time average  to output
    character(len=FIO_HSHORT) :: TUNIT !> time unit

    NAMELIST / HISTORYITEM / &
       ITEM, &
       TINT, &
       TAVG, &
       TUNIT

    integer :: n
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[HISTORY]/Categ[IO]'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_HISTORY,iostat=ierr)

    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_HISTORY. Check!'
       call PRC_MPIstop
    endif
    if( IO_L ) write(IO_FID_LOG,nml=PARAM_HISTORY)

    ! listup history request
    rewind(IO_FID_CONF)
    do n = 1, HIST_limit
       read(IO_FID_CONF,nml=HISTORYITEM,iostat=ierr)
       if( ierr /= 0 ) exit
    enddo
    HIST_nmax = n - 1

    if ( HIST_nmax > 0 ) then
       if( IO_L ) write(IO_FID_LOG,*) '*** NUmber of registered history item:', HIST_nmax
    else
       if( IO_L ) write(IO_FID_LOG,*) '*** No history file specified.'
    endif

    rewind(IO_FID_CONF)
    do n = 1, HIST_nmax
       ! set default
       ITEM  = 'unknown'
       TINT  = HISTORY_DEFAULT_TINTERVAL
       TAVG  = HISTORY_DEFAULT_TAVERAGE
       TUNIT = HISTORY_DEFAULT_TUNIT

       read(IO_FID_CONF,nml=HISTORYITEM,iostat=ierr)
       if( ierr /= 0 ) exit

       HIST_item(n) = ITEM
       call TIME_ymdhms2sec( HIST_tintsec(n), TINT, TUNIT )
       call TIME_ymdhms2sec( HIST_tavgsec(n), TAVG, TUNIT )
    enddo

    return
  end subroutine HISTRY_setup

  !-----------------------------------------------------------------------------
  subroutine HISTRY_reg( &
      itemid, &
      item,   &
      desc,   &
      ktype   )
    implicit none

    !---------------------------------------------------------------------------


    return
  end subroutine HISTRY_reg

  !-----------------------------------------------------------------------------
  subroutine HISTRY_put( &
      var,   &
      item,  &
      desc,  &
      ktype, )
    implicit none

    !---------------------------------------------------------------------------


    return
  end subroutine HISTRY_put

end module mod_history
!-------------------------------------------------------------------------------
