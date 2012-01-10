!-------------------------------------------------------------------------------
!> module TIME
!!
!! @par Description
!!          Calendar date and time module
!!
!! @author H.Tomita and SCALE developpers
!!
!! @par History
!! @li      2011-11-11 (H.Yashiro) [new]
!!
!<
!-------------------------------------------------------------------------------
module mod_time
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use mod_stdio, only: &
     IO_SYSCHR
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: TIME_setup
  public :: TIME_checkstate
  public :: TIME_advance
  public :: TIME_date2sec
  public :: TIME_sec2date
  public :: TIME_ymdhms2sec
  public :: TIME_rapstart
  public :: TIME_rapend
  public :: TIME_rapreport
  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  real(8), public, save :: TIME_DTSEC

  real(8), public, save :: TIME_DTSEC_ATMOS_DYN
  integer, public, save :: TIME_NSTEP_ATMOS_DYN
  real(8), public, save :: TIME_DTSEC_ATMOS_PHY_TB
  real(8), public, save :: TIME_DTSEC_ATMOS_PHY_MP
  real(8), public, save :: TIME_DTSEC_ATMOS_PHY_RD
  real(8), public, save :: TIME_DTSEC_ATMOS_RESTART
  real(8), public, save :: TIME_DTSEC_OCEAN

  real(8), public, save :: TIME_NOWSEC
  integer, public, save :: TIME_NOWDATE(6)
  real(8), public, save :: TIME_NOWMS
  integer, public, save :: TIME_NOWSTEP

  logical, public, save :: TIME_DOATMOS_step
  logical, public, save :: TIME_DOATMOS_DYN
  logical, public, save :: TIME_DOATMOS_PHY_TB
  logical, public, save :: TIME_DOATMOS_PHY_MP
  logical, public, save :: TIME_DOATMOS_PHY_RD
  logical, public, save :: TIME_DOATMOS_restart
  logical, public, save :: TIME_DOOCEAN_step
  logical, public, save :: TIME_DOend

  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  public :: TIME_rapid
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  real(8), private,      save :: TIME_STARTSEC
  integer, private,      save :: TIME_STARTDATE(6) = (/ 0000, 01, 01, 00, 00, 00 /)
  real(8), private,      save :: TIME_STARTMS      = 0.D0

  real(8), private,      save :: TIME_ENDSEC
  integer, private,      save :: TIME_ENDDATE(6)
  real(8), private,      save :: TIME_ENDMS

  integer, private,      save :: TIME_NSTEP

  real(8), private,      save :: TIME_RES_ATMOS_DYN     = 0.D0
  real(8), private,      save :: TIME_RES_ATMOS_PHY_TB  = 0.D0
  real(8), private,      save :: TIME_RES_ATMOS_PHY_MP  = 0.D0
  real(8), private,      save :: TIME_RES_ATMOS_PHY_RD  = 0.D0
  real(8), private,      save :: TIME_RES_ATMOS_RESTART = 0.D0
  real(8), private,      save :: TIME_RES_OCEAN         = 0.D0

  real(8), private, parameter :: TIME_DOY  = 365.D0
  real(8), private, parameter :: TIME_DOM(12) = (/ 31,28,31,30,31,30,31,31,30,31,30,31 /)
  real(8), private, parameter :: TIME_HOUR =  24.D0
  real(8), private, parameter :: TIME_MIN  =  60.D0
  real(8), private, parameter :: TIME_SEC  =  60.D0

  integer,                  private, parameter :: TIME_rapnlimit = 100
  integer,                  private,      save :: TIME_rapnmax   = 0
  character(len=IO_SYSCHR), private,      save :: TIME_rapname(TIME_rapnlimit)
  real(8),                  private,      save :: TIME_raptstr(TIME_rapnlimit)
  real(8),                  private,      save :: TIME_rapttot(TIME_rapnlimit)
  integer,                  private,      save :: TIME_rapnstr(TIME_rapnlimit)
  integer,                  private,      save :: TIME_rapnend(TIME_rapnlimit)

  real(8), private, parameter :: eps = 1.D-10 !> epsilon for timesec

  !-----------------------------------------------------------------------------
contains

  !-----------------------------------------------------------------------------
  !> Setup time
  !-----------------------------------------------------------------------------
  subroutine TIME_setup
    use mod_stdio, only: &
       IO_FID_CONF, &
       IO_FID_LOG,  &
       IO_L
    use mod_process, only: &
       PRC_MPIstop
    implicit none

    real(8)                  :: TIME_DURATION              = 1.D0
    character(len=IO_SYSCHR) :: TIME_DURATION_UNIT         = "MIN"
    real(8)                  :: TIME_DT                    = 300.D0
    character(len=IO_SYSCHR) :: TIME_DT_UNIT               = "MSEC"
    real(8)                  :: TIME_DT_ATMOS_DYN          = 0.D0
    character(len=IO_SYSCHR) :: TIME_DT_ATMOS_DYN_UNIT     = "SEC"
    real(8)                  :: TIME_DT_ATMOS_PHY_TB       = 1.D0
    character(len=IO_SYSCHR) :: TIME_DT_ATMOS_PHY_TB_UNIT  = "SEC"
    real(8)                  :: TIME_DT_ATMOS_PHY_MP       = 5.D0
    character(len=IO_SYSCHR) :: TIME_DT_ATMOS_PHY_MP_UNIT  = "SEC"
    real(8)                  :: TIME_DT_ATMOS_PHY_RD       = 1.D0
    character(len=IO_SYSCHR) :: TIME_DT_ATMOS_PHY_RD_UNIT  = "MIN"
    real(8)                  :: TIME_DT_ATMOS_RESTART      = 1.D0
    character(len=IO_SYSCHR) :: TIME_DT_ATMOS_RESTART_UNIT = "MIN"
    real(8)                  :: TIME_DT_OCEAN              = 1.D0
    character(len=IO_SYSCHR) :: TIME_DT_OCEAN_UNIT         = "HOUR"

    NAMELIST / PARAM_TIME / &
       TIME_STARTDATE,             &
       TIME_STARTMS,               &
       TIME_DURATION,              &
       TIME_DURATION_UNIT,         &
       TIME_DT,                    &
       TIME_DT_UNIT,               &
       TIME_DT_ATMOS_DYN,          &
       TIME_DT_ATMOS_DYN_UNIT,     &
       TIME_NSTEP_ATMOS_DYN,       &
       TIME_DT_ATMOS_PHY_TB,       &
       TIME_DT_ATMOS_PHY_TB_UNIT,  &
       TIME_DT_ATMOS_PHY_MP,       &
       TIME_DT_ATMOS_PHY_MP_UNIT,  &
       TIME_DT_ATMOS_PHY_RD,       &
       TIME_DT_ATMOS_PHY_RD_UNIT,  &
       TIME_DT_ATMOS_RESTART,      &
       TIME_DT_ATMOS_RESTART_UNIT, &
       TIME_DT_OCEAN,              &
       TIME_DT_OCEAN_UNIT

    real(8) :: TIME_DURATIONSEC

    integer :: ierr
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[TIME]/Categ[COMMON]'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_TIME,iostat=ierr)

    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_TIME. Check!'
       call PRC_MPIstop
    endif
    if( IO_L ) write(IO_FID_LOG,nml=PARAM_TIME)

    !--- calculate time
    call TIME_date2sec( TIME_STARTSEC, TIME_STARTDATE(:), TIME_STARTMS )

    call TIME_ymdhms2sec( TIME_DURATIONSEC, TIME_DURATION, TIME_DURATION_UNIT )

    TIME_ENDSEC = TIME_STARTSEC + TIME_DURATIONSEC

    call TIME_sec2date( TIME_ENDDATE(:), TIME_ENDMS, TIME_ENDSEC )

    call TIME_ymdhms2sec( TIME_DTSEC, TIME_DT, TIME_DT_UNIT )

    TIME_NSTEP = int( TIME_DURATIONSEC / TIME_DTSEC )

    TIME_NOWSEC     = TIME_STARTSEC
    TIME_NOWDATE(:) = TIME_STARTDATE(:)
    TIME_NOWMS      = TIME_STARTMS

    TIME_NOWSTEP    = 0

    if( IO_L ) write(IO_FID_LOG,'(1x,A,I4.4,A,I2.2,A,I2.2,A,I2.2,A,I2.2,A,I2.2,A,F6.3)') '*** START Date     : ', &
               TIME_STARTDATE(1),'/',TIME_STARTDATE(2),'/',TIME_STARTDATE(3),' ', &
               TIME_STARTDATE(4),':',TIME_STARTDATE(5),':',TIME_STARTDATE(6),' +', &
               TIME_STARTMS
    if( IO_L ) write(IO_FID_LOG,'(1x,A,I4.4,A,I2.2,A,I2.2,A,I2.2,A,I2.2,A,I2.2,A,F6.3)') '*** END   Date     : ', &
               TIME_ENDDATE(1),'/',TIME_ENDDATE(2),'/',TIME_ENDDATE(3),' ', &
               TIME_ENDDATE(4),':',TIME_ENDDATE(5),':',TIME_ENDDATE(6),' +', &
               TIME_ENDMS
    if( IO_L ) write(IO_FID_LOG,'(1x,A,F12.3)') '*** delta t (sec.) :', TIME_DTSEC
    if( IO_L ) write(IO_FID_LOG,*) '*** No. of steps   :', TIME_NSTEP

    if ( TIME_DTSEC <= 0.D0 ) then
       write(*,*) ' xxx Delta t <= 0.D0 is not accepted. Check!'
       call PRC_MPIstop
    endif

    !--- calculate intervals for atmosphere
    call TIME_ymdhms2sec( TIME_DTSEC_ATMOS_DYN,     TIME_DT_ATMOS_DYN,     TIME_DT_ATMOS_DYN_UNIT     )
    call TIME_ymdhms2sec( TIME_DTSEC_ATMOS_PHY_TB,  TIME_DT_ATMOS_PHY_TB,  TIME_DT_ATMOS_PHY_TB_UNIT  )
    call TIME_ymdhms2sec( TIME_DTSEC_ATMOS_PHY_MP,  TIME_DT_ATMOS_PHY_MP,  TIME_DT_ATMOS_PHY_MP_UNIT  )
    call TIME_ymdhms2sec( TIME_DTSEC_ATMOS_PHY_RD,  TIME_DT_ATMOS_PHY_RD,  TIME_DT_ATMOS_PHY_RD_UNIT  )
    call TIME_ymdhms2sec( TIME_DTSEC_ATMOS_RESTART, TIME_DT_ATMOS_RESTART, TIME_DT_ATMOS_RESTART_UNIT )
    call TIME_ymdhms2sec( TIME_DTSEC_OCEAN,         TIME_DT_OCEAN,         TIME_DT_OCEAN_UNIT         )

    TIME_DTSEC_ATMOS_DYN     = max( TIME_DTSEC_ATMOS_DYN,     TIME_DTSEC )
    TIME_DTSEC_ATMOS_PHY_TB  = max( TIME_DTSEC_ATMOS_PHY_TB,  TIME_DTSEC_ATMOS_DYN*TIME_NSTEP_ATMOS_DYN )
    TIME_DTSEC_ATMOS_PHY_MP  = max( TIME_DTSEC_ATMOS_PHY_MP,  TIME_DTSEC_ATMOS_DYN*TIME_NSTEP_ATMOS_DYN )
    TIME_DTSEC_ATMOS_PHY_RD  = max( TIME_DTSEC_ATMOS_PHY_RD,  TIME_DTSEC_ATMOS_DYN*TIME_NSTEP_ATMOS_DYN )
    TIME_DTSEC_ATMOS_RESTART = max( TIME_DTSEC_ATMOS_RESTART, TIME_DTSEC_ATMOS_DYN*TIME_NSTEP_ATMOS_DYN )
    TIME_DTSEC_OCEAN         = max( TIME_DTSEC_OCEAN,         TIME_DTSEC_ATMOS_DYN*TIME_NSTEP_ATMOS_DYN )

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** Time interval for atmospheric processes (sec.)'
    if( IO_L ) write(IO_FID_LOG,'(1x,A,F10.3)') '*** Dynamics (time)                  :', TIME_DTSEC_ATMOS_DYN
    if( IO_L ) write(IO_FID_LOG,'(1x,A,I4)')    '***          (step)                  :', TIME_NSTEP_ATMOS_DYN
    if( IO_L ) write(IO_FID_LOG,'(1x,A,F10.3)') '*** Physics,Turbulence               :', TIME_DTSEC_ATMOS_PHY_TB
    if( IO_L ) write(IO_FID_LOG,'(1x,A,F10.3)') '*** Physics,Cloud Microphysics       :', TIME_DTSEC_ATMOS_PHY_MP
    if( IO_L ) write(IO_FID_LOG,'(1x,A,F10.3)') '*** Physics,Radiation                :', TIME_DTSEC_ATMOS_PHY_RD
    if( IO_L ) write(IO_FID_LOG,*) '*** Time interval for ocean process (sec.)'
    if( IO_L ) write(IO_FID_LOG,'(1x,A,F10.3)') '*** SST update                       :', TIME_DTSEC_OCEAN
    if( IO_L ) write(IO_FID_LOG,*) '*** Time interval for Restart (sec.)'
    if( IO_L ) write(IO_FID_LOG,'(1x,A,F10.3)') '*** Prognostic Variables, Atmosphere :', TIME_DTSEC_ATMOS_RESTART

    return
  end subroutine TIME_setup

  !-----------------------------------------------------------------------------
  !> Check state of this time
  !-----------------------------------------------------------------------------
  subroutine TIME_checkstate
    use mod_stdio, only: &
       IO_FID_LOG, &
       IO_L
    implicit none

    TIME_DOATMOS_step   = .false.
    TIME_DOATMOS_DYN    = .false.
    TIME_DOATMOS_PHY_TB = .false.
    TIME_DOATMOS_PHY_MP = .false.
    TIME_DOATMOS_PHY_RD = .false.
    TIME_DOOCEAN_step   = .false.

    TIME_RES_ATMOS_DYN     = TIME_RES_ATMOS_DYN     + TIME_DTSEC
    TIME_RES_ATMOS_PHY_TB  = TIME_RES_ATMOS_PHY_TB  + TIME_DTSEC
    TIME_RES_ATMOS_PHY_MP  = TIME_RES_ATMOS_PHY_MP  + TIME_DTSEC
    TIME_RES_ATMOS_PHY_RD  = TIME_RES_ATMOS_PHY_RD  + TIME_DTSEC
    TIME_RES_OCEAN         = TIME_RES_OCEAN         + TIME_DTSEC
 
    if ( TIME_RES_ATMOS_DYN - TIME_DTSEC_ATMOS_DYN > -eps ) then
       TIME_DOATMOS_step  = .true.
       TIME_DOATMOS_DYN   = .true.
       TIME_RES_ATMOS_DYN = TIME_RES_ATMOS_DYN - TIME_DTSEC_ATMOS_DYN
    endif
    if ( TIME_RES_ATMOS_PHY_TB - TIME_DTSEC_ATMOS_PHY_TB > -eps ) then
       TIME_DOATMOS_step     = .true.
       TIME_DOATMOS_PHY_TB   = .true.
       TIME_RES_ATMOS_PHY_TB = TIME_RES_ATMOS_PHY_TB - TIME_DTSEC_ATMOS_PHY_TB
    endif
    if ( TIME_RES_ATMOS_PHY_MP - TIME_DTSEC_ATMOS_PHY_MP > -eps ) then
       TIME_DOATMOS_step     = .true.
       TIME_DOATMOS_PHY_MP   = .true.
       TIME_RES_ATMOS_PHY_MP = TIME_RES_ATMOS_PHY_MP - TIME_DTSEC_ATMOS_PHY_MP
    endif
    if ( TIME_RES_ATMOS_PHY_RD - TIME_DTSEC_ATMOS_PHY_RD > -eps ) then
       TIME_DOATMOS_step     = .true.
       TIME_DOATMOS_PHY_RD   = .true.
       TIME_RES_ATMOS_PHY_RD = TIME_RES_ATMOS_PHY_RD - TIME_DTSEC_ATMOS_PHY_RD
    endif

    if ( TIME_RES_OCEAN - TIME_DTSEC_OCEAN > -eps ) then
       TIME_DOOCEAN_step = .true.
       TIME_RES_OCEAN    = TIME_RES_OCEAN - TIME_DTSEC_OCEAN
    endif

    call TIME_sec2date( TIME_NOWDATE(:), TIME_NOWMS, TIME_NOWSEC )

    if( IO_L ) write(IO_FID_LOG,'(1x,A,I4.4,A,I2.2,A,I2.2,A,I2.2,A,I2.2,A,I2.2,A,F6.3,A,I6)') '*** TIME: ', &
               TIME_NOWDATE(1),'/',TIME_NOWDATE(2),'/',TIME_NOWDATE(3),' ', &
               TIME_NOWDATE(4),':',TIME_NOWDATE(5),':',TIME_NOWDATE(6),' +', &
               TIME_NOWMS,' STEP:',TIME_NOWSTEP

  end subroutine TIME_checkstate

  !-----------------------------------------------------------------------------
  !> Advance the time & eval. timing to restart&stop
  !-----------------------------------------------------------------------------
  subroutine TIME_advance
    use mod_stdio, only: &
       IO_FID_LOG, &
       IO_L
    implicit none

    TIME_DOend = .false.

    TIME_NOWSEC  = TIME_NOWSEC + TIME_DTSEC
    TIME_NOWSTEP = TIME_NOWSTEP + 1

    if( TIME_NOWSEC - TIME_ENDSEC > -eps ) TIME_DOend = .true.

    TIME_DOATMOS_restart = .false.

    TIME_RES_ATMOS_RESTART = TIME_RES_ATMOS_RESTART + TIME_DTSEC

    if ( TIME_RES_ATMOS_RESTART - TIME_DTSEC_ATMOS_RESTART > -eps ) then
       TIME_DOATMOS_restart   = .true.
       TIME_RES_ATMOS_RESTART = TIME_RES_ATMOS_RESTART - TIME_DTSEC_ATMOS_RESTART
    elseif( TIME_DOend ) then
       TIME_DOATMOS_restart   = .true.
    endif

  end subroutine TIME_advance

  !-----------------------------------------------------------------------------
  !> convert date to second
  !
  !@todo fit to gregorian calendar
  !-----------------------------------------------------------------------------
  subroutine TIME_date2sec( &
     second,   &
     datetime, &
     microsec  )
    implicit none

    real(8), intent(out) :: second
    integer, intent( in) :: datetime(6)
    real(8), intent( in) :: microsec

    integer :: m
    !---------------------------------------------------------------------------

    second = 0.D0
    second = second + datetime(1) * TIME_SEC * TIME_MIN * TIME_HOUR * TIME_DOY
    if ( datetime(2) > 1 ) then
       do m = 1, datetime(2)-1
          second = second + TIME_DOM(m) * TIME_SEC * TIME_MIN * TIME_HOUR
       enddo
    endif
    second = second + ( datetime(3)-1 ) * TIME_SEC * TIME_MIN * TIME_HOUR
    second = second + datetime(4) * TIME_SEC * TIME_MIN
    second = second + datetime(5) * TIME_SEC
    second = second + datetime(6)
    second = second + microsec

    return
  end subroutine TIME_date2sec

  !-----------------------------------------------------------------------------
  !> convert date to second
  !
  !@todo fit to gregorian calendar
  !-----------------------------------------------------------------------------
  subroutine TIME_sec2date( &
     datetime, &
     microsec, &
     second    )
    implicit none

    integer, intent(out) :: datetime(6)
    real(8), intent(out) :: microsec
    real(8), intent( in) :: second

    real(8) :: temp
    real(8) :: nsec, nmin, nhour, nday

    integer :: m
    !---------------------------------------------------------------------------

    temp  = dble( int( second,kind=8 ) )
    microsec = second - temp
    nsec  = temp

    temp  = mod( nsec, TIME_SEC )
    datetime(6) = temp
    nmin  = ( nsec-temp ) / TIME_SEC

    temp  = mod( nmin, TIME_MIN )
    datetime(5) = temp
    nhour = ( nmin-temp ) / TIME_MIN

    temp  = mod( nhour, TIME_HOUR )
    datetime(4) = temp
    nday  = ( nhour-temp ) / TIME_HOUR

    temp  = mod( nday, TIME_DOY )
    datetime(1) = ( nday-temp ) / TIME_DOY
    nday  = temp

    do m = 1, 12
       if ( nday <= TIME_DOM(m) ) then
          datetime(2) = m
          datetime(3) = nday + 1
          exit
       endif
       nday = nday - TIME_DOM(m)
    enddo

    return
  end subroutine TIME_sec2date

  !-----------------------------------------------------------------------------
  !> convert several units to second
  !
  !@todo fit to gregorian calendar
  !-----------------------------------------------------------------------------
  subroutine TIME_ymdhms2sec( &
     second, &
     value,  &
     unit    )
    use mod_process, only: &
       PRC_MPIstop
    implicit none

    real(8),          intent(out) :: second
    real(8),          intent( in) :: value
    character(len=*), intent( in) :: unit
    !---------------------------------------------------------------------------

    select case(trim(unit))
    case('MSEC')
       second = value * 1.D-3
    case('SEC')
       second = value
    case('MIN')
       second = value * TIME_SEC
    case('HOUR')
       second = value * TIME_SEC * TIME_MIN
    case('DAY')
       second = value * TIME_SEC * TIME_MIN * TIME_HOUR
    case default
       write(*,*) ' xxx Unsupported UNIT:', trim(unit), value
       call PRC_MPIstop
    endselect

    return
  end subroutine TIME_ymdhms2sec

  !-----------------------------------------------------------------------------
  function TIME_rapid( rapname ) result(id)
    implicit none

    character(len=*), intent(in) :: rapname

    integer :: id
    !---------------------------------------------------------------------------

    do id = 1, TIME_rapnlimit
       if( trim(rapname) == trim(TIME_rapname(id)) ) return
    enddo

    TIME_rapnmax     = TIME_rapnmax + 1
    id               = TIME_rapnmax
    TIME_rapname(id) = trim(rapname)

  end function TIME_rapid

  !-----------------------------------------------------------------------------
  subroutine TIME_rapstart( rapname )
    use mod_process, only: &
       PRC_MPItime
    implicit none

    character(len=*), intent(in) :: rapname

    integer :: id
    !---------------------------------------------------------------------------

    id = TIME_rapid( rapname )

    TIME_raptstr(id) = PRC_MPItime()
    TIME_rapnstr(id) = TIME_rapnstr(id) + 1

    return
  end subroutine TIME_rapstart

  !-----------------------------------------------------------------------------
  subroutine TIME_rapend( rapname )
    use mod_process, only: &
       PRC_MPItime
    implicit none

    character(len=*), intent(in) :: rapname

    integer :: id
    !---------------------------------------------------------------------------

    id = TIME_rapid( rapname )

    TIME_rapttot(id) = TIME_rapttot(id) + ( PRC_MPItime()-TIME_raptstr(id) )
    TIME_rapnend(id) = TIME_rapnend(id) + 1

    return
  end subroutine TIME_rapend

  !-----------------------------------------------------------------------------
  subroutine TIME_rapreport
    use mod_stdio, only: &
       IO_FID_LOG,  &
       IO_L
    implicit none

    integer :: id
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** Computational Time Report'

    do id = 1, TIME_rapnmax

       if ( TIME_rapnstr(id) /= TIME_rapnend(id) ) then
           write(*,*) '*** Computational Time Report'
       endif

       if( IO_L ) write(IO_FID_LOG,'(1x,A,I3.3,A,A,A,F10.3,A,I5)') &
                                '*** ID=',id,':',TIME_rapname(id), &
                                ' T=',TIME_rapttot(id),' N=',TIME_rapnstr(id)
    enddo

    return
  end subroutine TIME_rapreport

end module mod_time
