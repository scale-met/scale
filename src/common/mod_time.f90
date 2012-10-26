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
  use mod_stdio, only : &
     IO_FID_LOG, &
     IO_L,       &
     IO_SYSCHR
  use dc_types, only : &
      DP
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
  !++ included parameters
  !
  include 'inc_precision.h'
  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  real(DP), public, save :: TIME_DTSEC                !< time interval of model       [sec]

  real(DP), public, save :: TIME_DTSEC_ATMOS_DYN      !< time interval of dynamics     [sec]
  integer, public, save :: TIME_NSTEP_ATMOS_DYN = 20 !< small step of dynamics
  real(DP), public, save :: TIME_DTSEC_ATMOS_PHY_TB   !< time interval of turbulence   [sec]
  real(DP), public, save :: TIME_DTSEC_ATMOS_PHY_MP   !< time interval of microphysics [sec]
  real(DP), public, save :: TIME_DTSEC_ATMOS_PHY_RD   !< time interval of radiation    [sec]
  real(DP), public, save :: TIME_DTSEC_ATMOS_RESTART  !< time interval of restart      [sec]
  real(DP), public, save :: TIME_DTSEC_OCEAN          !< time interval of ocean        [sec]
  real(DP), public, save :: TIME_DTSEC_OCEAN_RESTART  !< time interval of restart      [sec]

  real(DP), public, save :: TIME_NOWSEC               !< current time [sec]
  integer, public, save :: TIME_NOWSTEP              !< current step [number]

  integer, public, save :: TIME_NOWDATE(6)           !< current time [YYYY MM DD HH MM SS]
  real(DP), public, save :: TIME_NOWMS                !< subsecond part of current time [millisec]
  real(DP), public, save :: TIME_NOWSECL              !< current time [sec]

  logical, public, save :: TIME_DOATMOS_step         !< execute atmospheric component in this step?
  logical, public, save :: TIME_DOATMOS_DYN          !< execute dynamics?
  logical, public, save :: TIME_DOATMOS_PHY_TB       !< execute physics(turbulence)?
  logical, public, save :: TIME_DOATMOS_PHY_MP       !< execute physics(microphysics)?
  logical, public, save :: TIME_DOATMOS_PHY_RD       !< execute physics(radiation)?
  logical, public, save :: TIME_DOATMOS_restart      !< execute restart output?
  logical, public, save :: TIME_DOOCEAN_step         !< execute ocean component in this step?
  logical, public, save :: TIME_DOOCEAN_restart      !< execute restart output?
  logical, public, save :: TIME_DOend                !< finish program in this step?

  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  private :: TIME_rapid

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  integer, private,      save :: TIME_STARTDATE(6) = (/ 2000, 01, 01, 00, 00, 00 /)
  real(DP), private,      save :: TIME_STARTMS      = 0.0_DP                           !< [millisec]
  real(DP), private,      save :: TIME_STARTSECL

  integer, private,      save :: TIME_ENDDATE(6)
  real(DP), private,      save :: TIME_ENDMS
  real(DP), private,      save :: TIME_ENDSECL

  integer, private,      save :: TIME_NSTEP

  real(DP), private,      save :: TIME_RES_ATMOS_DYN     = 0.0_DP
  real(DP), private,      save :: TIME_RES_ATMOS_PHY_TB  = 0.E0_DP
  real(DP), private,      save :: TIME_RES_ATMOS_PHY_MP  = 0.E0_DP
  real(DP), private,      save :: TIME_RES_ATMOS_PHY_RD  = 0.E0_DP
  real(DP), private,      save :: TIME_RES_ATMOS_RESTART = 0.E0_DP
  real(DP), private,      save :: TIME_RES_OCEAN         = 0.E0_DP
  real(DP), private,      save :: TIME_RES_OCEAN_RESTART = 0.E0_DP

  real(DP), private, parameter :: TIME_DOY     = 365.E0_DP
  real(DP), private, parameter :: TIME_DOM(12) = (/ 31,28,31,30,31,30,31,31,30,31,30,31 /)
  real(DP), private, parameter :: TIME_HOUR    =  24.E0_DP
  real(DP), private, parameter :: TIME_MIN     =  60.E0_DP
  real(DP), private, parameter :: TIME_SEC     =  60.E0_DP

  integer,                  private, parameter :: TIME_rapnlimit = 100
  integer,                  private,      save :: TIME_rapnmax   = 0
  character(len=IO_SYSCHR), private,      save :: TIME_rapname(TIME_rapnlimit)
  real(DP),                  private,      save :: TIME_raptstr(TIME_rapnlimit)
  real(DP),                  private,      save :: TIME_rapttot(TIME_rapnlimit)
  integer,                  private,      save :: TIME_rapnstr(TIME_rapnlimit)
  integer,                  private,      save :: TIME_rapnend(TIME_rapnlimit)

  real(DP), private, parameter :: eps = 1.E-10_DP !> epsilon for timesec

  !-----------------------------------------------------------------------------
contains

  !-----------------------------------------------------------------------------
  !> Setup time
  !-----------------------------------------------------------------------------
  subroutine TIME_setup
    use mod_stdio, only: &
       IO_FID_CONF
    use mod_process, only: &
       PRC_MPIstop
    implicit none

    real(DP)                  :: TIME_DURATION              = 60.0E0_DP
    character(len=IO_SYSCHR) :: TIME_DURATION_UNIT         = "MIN"
    real(DP)                  :: TIME_DT                    =  0.6E0_DP
    character(len=IO_SYSCHR) :: TIME_DT_UNIT               = "SEC"
    real(DP)                  :: TIME_DT_ATMOS_DYN          =  0.03E0_DP
    character(len=IO_SYSCHR) :: TIME_DT_ATMOS_DYN_UNIT     = "SEC"
    real(DP)                  :: TIME_DT_ATMOS_PHY_TB       =  0.6E0_DP
    character(len=IO_SYSCHR) :: TIME_DT_ATMOS_PHY_TB_UNIT  = "SEC"
    real(DP)                  :: TIME_DT_ATMOS_PHY_MP       =  0.6E0_DP
    character(len=IO_SYSCHR) :: TIME_DT_ATMOS_PHY_MP_UNIT  = "SEC"
    real(DP)                  :: TIME_DT_ATMOS_PHY_RD       =  0.6E0_DP
    character(len=IO_SYSCHR) :: TIME_DT_ATMOS_PHY_RD_UNIT  = "SEC"
    real(DP)                  :: TIME_DT_ATMOS_RESTART      = 60.0E0_DP
    character(len=IO_SYSCHR) :: TIME_DT_ATMOS_RESTART_UNIT = "SEC"
    real(DP)                  :: TIME_DT_OCEAN              = 60.0E0_DP
    character(len=IO_SYSCHR) :: TIME_DT_OCEAN_UNIT         = "MIN"
    real(DP)                  :: TIME_DT_OCEAN_RESTART      = 60.0E0_DP
    character(len=IO_SYSCHR) :: TIME_DT_OCEAN_RESTART_UNIT = "SEC"

    NAMELIST / PARAM_TIME / &
       TIME_STARTDATE,             &
       TIME_STARTMS,               &
       TIME_DURATION,              &
       TIME_DURATION_UNIT,         &
       TIME_DT,                    &
       TIME_DT_UNIT,               &
       TIME_DT_ATMOS_DYN,          &
       TIME_DT_ATMOS_DYN_UNIT,     &
       TIME_DT_ATMOS_PHY_TB,       &
       TIME_DT_ATMOS_PHY_TB_UNIT,  &
       TIME_DT_ATMOS_PHY_MP,       &
       TIME_DT_ATMOS_PHY_MP_UNIT,  &
       TIME_DT_ATMOS_PHY_RD,       &
       TIME_DT_ATMOS_PHY_RD_UNIT,  &
       TIME_DT_ATMOS_RESTART,      &
       TIME_DT_ATMOS_RESTART_UNIT, &
       TIME_DT_OCEAN,              &
       TIME_DT_OCEAN_UNIT,         &
       TIME_DT_OCEAN_RESTART,      &
       TIME_DT_OCEAN_RESTART_UNIT

    real(DP) :: TIME_DURATIONSEC
    real(DP) :: temp

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
    call TIME_date2sec( TIME_STARTSECL, TIME_STARTDATE(:) )
    TIME_STARTMS = TIME_STARTMS * 1.E-3_DP

    call TIME_ymdhms2sec( TIME_DURATIONSEC, TIME_DURATION, TIME_DURATION_UNIT )

    temp  = dble( int( TIME_STARTMS+TIME_DURATIONSEC,kind=8 ) )
    TIME_ENDMS   = TIME_STARTMS   + TIME_DURATIONSEC - temp
    TIME_ENDSECL = TIME_STARTSECL + temp

    call TIME_sec2date( TIME_ENDDATE(:), TIME_ENDSECL )

    call TIME_ymdhms2sec( TIME_DTSEC, TIME_DT, TIME_DT_UNIT )

    TIME_NSTEP_ATMOS_DYN = int ( TIME_DTSEC / TIME_DTSEC_ATMOS_DYN )

    TIME_NSTEP = int( TIME_DURATIONSEC / TIME_DTSEC )

    TIME_NOWSECL    = TIME_STARTSECL
    TIME_NOWDATE(:) = TIME_STARTDATE(:)
    TIME_NOWMS      = TIME_STARTMS

    TIME_NOWSEC     = TIME_STARTSECL + TIME_STARTMS
    TIME_NOWSTEP    = 1

    if( IO_L ) write(IO_FID_LOG,'(1x,A,I4.4,A,I2.2,A,I2.2,A,I2.2,A,I2.2,A,I2.2,A,F6.3)') '*** START Date     : ', &
               TIME_STARTDATE(1),'/',TIME_STARTDATE(2),'/',TIME_STARTDATE(3),' ', &
               TIME_STARTDATE(4),':',TIME_STARTDATE(5),':',TIME_STARTDATE(6),' +', &
               TIME_STARTMS
    if( IO_L ) write(IO_FID_LOG,'(1x,A,I4.4,A,I2.2,A,I2.2,A,I2.2,A,I2.2,A,I2.2,A,F6.3)') '*** END   Date     : ', &
               TIME_ENDDATE(1),'/',TIME_ENDDATE(2),'/',TIME_ENDDATE(3),' ', &
               TIME_ENDDATE(4),':',TIME_ENDDATE(5),':',TIME_ENDDATE(6),' +', &
               TIME_ENDMS
    if( IO_L ) write(IO_FID_LOG,'(1x,A,F10.3)') '*** delta t (sec.) :', TIME_DTSEC
    if( IO_L ) write(IO_FID_LOG,*) '*** No. of steps   :', TIME_NSTEP

    if ( TIME_DTSEC <= 0.E0_DP ) then
       write(*,*) ' xxx Delta t <= 0.E0_DP is not accepted. Check!'
       call PRC_MPIstop
    endif

    !--- calculate intervals for atmosphere
    call TIME_ymdhms2sec( TIME_DTSEC_ATMOS_DYN,     TIME_DT_ATMOS_DYN,     TIME_DT_ATMOS_DYN_UNIT     )
    call TIME_ymdhms2sec( TIME_DTSEC_ATMOS_PHY_TB,  TIME_DT_ATMOS_PHY_TB,  TIME_DT_ATMOS_PHY_TB_UNIT  )
    call TIME_ymdhms2sec( TIME_DTSEC_ATMOS_PHY_MP,  TIME_DT_ATMOS_PHY_MP,  TIME_DT_ATMOS_PHY_MP_UNIT  )
    call TIME_ymdhms2sec( TIME_DTSEC_ATMOS_PHY_RD,  TIME_DT_ATMOS_PHY_RD,  TIME_DT_ATMOS_PHY_RD_UNIT  )
    call TIME_ymdhms2sec( TIME_DTSEC_ATMOS_RESTART, TIME_DT_ATMOS_RESTART, TIME_DT_ATMOS_RESTART_UNIT )
    call TIME_ymdhms2sec( TIME_DTSEC_OCEAN,         TIME_DT_OCEAN,         TIME_DT_OCEAN_UNIT         )
    call TIME_ymdhms2sec( TIME_DTSEC_OCEAN_RESTART, TIME_DT_OCEAN_RESTART, TIME_DT_OCEAN_RESTART_UNIT )

    TIME_DTSEC_ATMOS_DYN     = max( TIME_DTSEC_ATMOS_DYN,     TIME_DTSEC          /TIME_NSTEP_ATMOS_DYN )
    TIME_DTSEC_ATMOS_PHY_TB  = max( TIME_DTSEC_ATMOS_PHY_TB,  TIME_DTSEC_ATMOS_DYN*TIME_NSTEP_ATMOS_DYN )
    TIME_DTSEC_ATMOS_PHY_MP  = max( TIME_DTSEC_ATMOS_PHY_MP,  TIME_DTSEC_ATMOS_DYN*TIME_NSTEP_ATMOS_DYN )
    TIME_DTSEC_ATMOS_PHY_RD  = max( TIME_DTSEC_ATMOS_PHY_RD,  TIME_DTSEC_ATMOS_DYN*TIME_NSTEP_ATMOS_DYN )
    TIME_DTSEC_ATMOS_RESTART = max( TIME_DTSEC_ATMOS_RESTART, TIME_DTSEC_ATMOS_DYN*TIME_NSTEP_ATMOS_DYN )
    TIME_DTSEC_OCEAN         = max( TIME_DTSEC_OCEAN,         TIME_DTSEC_ATMOS_DYN*TIME_NSTEP_ATMOS_DYN )
    TIME_DTSEC_OCEAN_RESTART = max( TIME_DTSEC_OCEAN_RESTART, TIME_DTSEC_ATMOS_DYN*TIME_NSTEP_ATMOS_DYN )

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** Time interval for atmospheric processes (sec.)'
    if( IO_L ) write(IO_FID_LOG,'(1x,A,F10.3)') '*** Dynamics (time)                  :', TIME_DTSEC_ATMOS_DYN
    if( IO_L ) write(IO_FID_LOG,'(1x,A,I4)')    '***          (step)                  :', TIME_NSTEP_ATMOS_DYN
    if( IO_L ) write(IO_FID_LOG,'(1x,A,F10.3)') '*** Physics, Turbulence              :', TIME_DTSEC_ATMOS_PHY_TB
    if( IO_L ) write(IO_FID_LOG,'(1x,A,F10.3)') '*** Physics, Cloud Microphysics      :', TIME_DTSEC_ATMOS_PHY_MP
    if( IO_L ) write(IO_FID_LOG,'(1x,A,F10.3)') '*** Physics, Radiation               :', TIME_DTSEC_ATMOS_PHY_RD
    if( IO_L ) write(IO_FID_LOG,*) '*** Time interval for ocean process (sec.)'
    if( IO_L ) write(IO_FID_LOG,'(1x,A,F10.3)') '*** SST update                       :', TIME_DTSEC_OCEAN
    if( IO_L ) write(IO_FID_LOG,*) '*** Time interval for Restart (sec.)'
    if( IO_L ) write(IO_FID_LOG,'(1x,A,F10.3)') '*** Prognostic Variables, Atmosphere :', TIME_DTSEC_ATMOS_RESTART
    if( IO_L ) write(IO_FID_LOG,'(1x,A,F10.3)') '*** Ocean Variables                  :', TIME_DTSEC_OCEAN_RESTART

    return
  end subroutine TIME_setup

  !-----------------------------------------------------------------------------
  !> Check state of this time
  !-----------------------------------------------------------------------------
  subroutine TIME_checkstate
    implicit none
    !---------------------------------------------------------------------------

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

    call TIME_sec2date( TIME_NOWDATE(:), TIME_NOWSEC )

    if( IO_L ) write(IO_FID_LOG,'(1x,A,I4.4,A,I2.2,A,I2.2,A,I2.2,A,I2.2,A,I2.2,A,F6.3,A,I6)') '*** TIME: ', &
               TIME_NOWDATE(1),'/',TIME_NOWDATE(2),'/',TIME_NOWDATE(3),' ', &
               TIME_NOWDATE(4),':',TIME_NOWDATE(5),':',TIME_NOWDATE(6),' +', &
               TIME_NOWMS,' STEP:',TIME_NOWSTEP

  end subroutine TIME_checkstate

  !-----------------------------------------------------------------------------
  !> Advance the time & eval. timing to restart&stop
  !-----------------------------------------------------------------------------
  subroutine TIME_advance
    use mod_process, only: &
       PRC_master, &
       PRC_myrank
    implicit none

    logical :: exists

    real(DP) :: temp
    !---------------------------------------------------------------------------

    TIME_DOend = .false.

    temp  = dble( int( TIME_NOWMS+TIME_DTSEC,kind=8 ) )
    TIME_NOWMS   = TIME_NOWMS + TIME_DTSEC - temp
    TIME_NOWMS   = nint( TIME_NOWMS * 1.E5_DP ) * 1.D-5
    TIME_NOWSECL = TIME_NOWSECL + temp

    TIME_NOWSEC  = TIME_NOWSECL + TIME_NOWMS
    TIME_NOWSTEP = TIME_NOWSTEP + 1

    if (       TIME_NOWSECL - TIME_ENDSECL > -eps &
         .AND. TIME_NOWMS   - TIME_ENDMS   > -eps ) then
      TIME_DOend = .true.
    endif

    ! QUIT file control
    if ( PRC_myrank == PRC_master ) then ! master node
       inquire(file='QUIT', exist=exists)

       if( exists ) TIME_DOend = .true.
    endif

    TIME_DOATMOS_restart = .false.
    TIME_DOOCEAN_restart = .false.

    TIME_RES_ATMOS_RESTART = TIME_RES_ATMOS_RESTART + TIME_DTSEC
    TIME_RES_OCEAN_RESTART = TIME_RES_OCEAN_RESTART + TIME_DTSEC

    if ( TIME_RES_ATMOS_RESTART - TIME_DTSEC_ATMOS_RESTART > -eps ) then
       TIME_DOATMOS_restart   = .true.
       TIME_RES_ATMOS_RESTART = TIME_RES_ATMOS_RESTART - TIME_DTSEC_ATMOS_RESTART
    elseif( TIME_DOend ) then
       TIME_DOATMOS_restart   = .true.
    endif

    if ( TIME_RES_OCEAN_RESTART - TIME_DTSEC_OCEAN_RESTART > -eps ) then
       TIME_DOOCEAN_restart   = .true.
       TIME_RES_OCEAN_RESTART = TIME_RES_OCEAN_RESTART - TIME_DTSEC_OCEAN_RESTART
    elseif( TIME_DOend ) then
       TIME_DOOCEAN_restart   = .true.
    endif

  end subroutine TIME_advance

  !-----------------------------------------------------------------------------
  !> convert date to second
  !@todo fit to gregorian calendar
  !-----------------------------------------------------------------------------
  subroutine TIME_date2sec( &
     second,  &
     datetime )
    implicit none

    real(DP), intent(out) :: second
    integer, intent( in) :: datetime(6)

    integer :: m
    !---------------------------------------------------------------------------

    second = 0.E0_DP
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

    return
  end subroutine TIME_date2sec

  !-----------------------------------------------------------------------------
  !> convert date to second
  !
  !@todo fit to gregorian calendar
  !-----------------------------------------------------------------------------
  subroutine TIME_sec2date( &
     datetime, &
     second    )
    implicit none

    integer, intent(out) :: datetime(6)
    real(DP), intent( in) :: second

    real(DP) :: temp
    real(DP) :: nmin, nhour, nday

    integer :: m
    !---------------------------------------------------------------------------

    temp  = mod( second, TIME_SEC )
    datetime(6) = int( temp )
    nmin  = ( second-temp ) / TIME_SEC

    temp  = mod( nmin, TIME_MIN )
    datetime(5) = int( temp )
    nhour = ( nmin-temp ) / TIME_MIN

    temp  = mod( nhour, TIME_HOUR )
    datetime(4) = int( temp )
    nday  = ( nhour-temp ) / TIME_HOUR

    temp  = mod( nday, TIME_DOY )
    datetime(1) = int( ( nday-temp ) / TIME_DOY )
    nday  = temp

    do m = 1, 12
       if ( nday <= TIME_DOM(m) ) then
          datetime(2) = m
          datetime(3) = int( nday ) + 1
          exit
       endif
       nday = nday - TIME_DOM(m)
    enddo

    return
  end subroutine TIME_sec2date

  !-----------------------------------------------------------------------------
  !> convert several units to second
  !@todo fit to gregorian calendar
  !-----------------------------------------------------------------------------
  subroutine TIME_ymdhms2sec( &
     second, &
     value,  &
     unit    )
    use mod_process, only: &
       PRC_MPIstop
    implicit none

    real(DP),          intent(out) :: second
    real(DP),          intent( in) :: value
    character(len=*), intent( in) :: unit
    !---------------------------------------------------------------------------

    select case(trim(unit))
    case('MSEC')
       second = value * 1.E-3_DP
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
    character (len=IO_SYSCHR) :: trapname
    !---------------------------------------------------------------------------

    trapname = trim(rapname)

    do id = 1, TIME_rapnmax
       if( trapname == TIME_rapname(id) ) return
    enddo

    TIME_rapnmax     = TIME_rapnmax + 1
    id               = TIME_rapnmax
    TIME_rapname(id) = trapname

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
    use mod_stdio, only : &
       IO_LOG_SUPPRESS, &
       IO_LOG_ALLNODE
    use mod_process, only: &
       PRC_master, &
       PRC_myrank, &
       PRC_MPItimestat
    implicit none

    real(DP) :: avgvar(TIME_rapnlimit)
    real(DP) :: maxvar(TIME_rapnlimit)
    real(DP) :: minvar(TIME_rapnlimit)
    integer :: maxidx(TIME_rapnlimit)
    integer :: minidx(TIME_rapnlimit)

    integer :: id
    !---------------------------------------------------------------------------

    do id = 1, TIME_rapnmax
       if ( TIME_rapnstr(id) /= TIME_rapnend(id) ) then
           write(*,*) '*** Mismatch Report',id,TIME_rapname(id),TIME_rapnstr(id),TIME_rapnend(id)
       endif
    enddo

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** Computational Time Report'

    if ( IO_LOG_ALLNODE ) then ! report for each node

       do id = 1, TIME_rapnmax
          if( IO_L ) write(IO_FID_LOG,'(1x,A,I3.3,A,A,A,F10.3,A,I7)') &
                     '*** ID=',id,' : ',TIME_rapname(id),' T=',TIME_rapttot(id),' N=',TIME_rapnstr(id)
       enddo

    else

       call PRC_MPItimestat( avgvar(1:TIME_rapnmax), &
                             maxvar(1:TIME_rapnmax), &
                             minvar(1:TIME_rapnmax), &
                             maxidx(1:TIME_rapnmax), &
                             minidx(1:TIME_rapnmax), &
                             TIME_rapttot(1:TIME_rapnmax) )

       do id = 1, TIME_rapnmax
          if( IO_L ) write(IO_FID_LOG,'(1x,A,I3.3,A,A,A,F10.3,A,F10.3,A,I5,A,A,F10.3,A,I5,A,A,I7)') &
                     '*** ID=',id,' : ',TIME_rapname(id), &
                     ' T(avg)=',avgvar(id), &
                     ', T(max)=',maxvar(id),'[',maxidx(id),']', &
                     ', T(min)=',minvar(id),'[',minidx(id),']', &
                     ' N=',TIME_rapnstr(id)
       enddo

       if ( IO_LOG_SUPPRESS ) then ! report to STDOUT
          if ( PRC_myrank == PRC_master ) then ! master node
             write(*,*) '*** Computational Time Report'
             do id = 1, TIME_rapnmax
                write(*,'(1x,A,I3.3,A,A,A,F10.3,A,F10.3,A,I5,A,A,F10.3,A,I5,A,A,I7)') &
                     '*** ID=',id,' : ',TIME_rapname(id), &
                     ' T(avg)=',avgvar(id), &
                     ', T(max)=',maxvar(id),'[',maxidx(id),']', &
                     ', T(min)=',minvar(id),'[',minidx(id),']', &
                     ' N=',TIME_rapnstr(id)
             enddo
          endif
       endif
    endif

    return
  end subroutine TIME_rapreport

end module mod_time

