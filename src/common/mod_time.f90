!-------------------------------------------------------------------------------
!> module TIME
!!
!! @par Description
!!          Time management and timer module
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2011-11-11 (H.Yashiro)  [new]
!! @li      2013-01-29 (H.Yashiro)  [mod] exclude calendar tools
!!
!<
module mod_time
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use dc_types, only: &
     DP
  use mod_stdio, only: &
     IO_FID_LOG, &
     IO_L,       &
     IO_SYSCHR
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ included parameters
  !
  include 'inc_precision.h'

  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: TIME_setup
  public :: TIME_checkstate
  public :: TIME_advance

  public :: TIME_rapstart
  public :: TIME_rapend
  public :: TIME_rapreport

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  real(DP), public, save :: TIME_DTSEC                !< time interval of model              [sec]

  real(DP), public, save :: TIME_DTSEC_ATMOS_DYN      !< time interval of dynamics           [sec]
  integer,  public, save :: TIME_NSTEP_ATMOS_DYN      !< small step of dynamics
  real(DP), public, save :: TIME_DTSEC_ATMOS_PHY_SF   !< time interval of surface flux       [sec]
  real(DP), public, save :: TIME_DTSEC_ATMOS_PHY_TB   !< time interval of turbulence         [sec]
  real(DP), public, save :: TIME_DTSEC_ATMOS_PHY_MP   !< time interval of microphysics       [sec]
  real(DP), public, save :: TIME_DTSEC_ATMOS_PHY_RD   !< time interval of radiation          [sec]
  real(DP), public, save :: TIME_DTSEC_ATMOS_RESTART  !< time interval of atmosphere restart [sec]
  real(DP), public, save :: TIME_DTSEC_OCEAN          !< time interval of ocean step         [sec]
  real(DP), public, save :: TIME_DTSEC_OCEAN_RESTART  !< time interval of ocean restart      [sec]
  real(DP), public, save :: TIME_DTSEC_LAND           !< time interval of land step          [sec]
  real(DP), public, save :: TIME_DTSEC_LAND_RESTART   !< time interval of land restart       [sec]
  real(DP), public, save :: TIME_DTSEC_CPL            !< time interval of coupler calc.      [sec]
  real(DP), public, save :: TIME_DTSEC_CPL_RESTART    !< time interval of coupler restart    [sec]

  integer,  public, save :: TIME_NOWDATE(6)           !< current time [YYYY MM DD HH MM SS]
  real(DP), public, save :: TIME_NOWMS                !< subsecond part of current time [millisec]
  integer,  public, save :: TIME_NOWDAY               !< absolute day of current time [day]
  real(DP), public, save :: TIME_NOWSEC               !< subday part  of current time [sec]
  real(DP), public, save :: TIME_NOWDAYSEC            !< second of current time [sec]
  integer,  public, save :: TIME_NOWSTEP              !< current step [number]

  logical,  public, save :: TIME_DOATMOS_step         !< execute atmospheric component in this step?
  logical,  public, save :: TIME_DOATMOS_DYN          !< execute dynamics?
  logical,  public, save :: TIME_DOATMOS_PHY_SF       !< execute physics(surface flux)?
  logical,  public, save :: TIME_DOATMOS_PHY_TB       !< execute physics(turbulence)?
  logical,  public, save :: TIME_DOATMOS_PHY_MP       !< execute physics(microphysics)?
  logical,  public, save :: TIME_DOATMOS_PHY_RD       !< execute physics(radiation)?
  logical,  public, save :: TIME_DOATMOS_restart      !< execute atmosphere restart output?
  logical,  public, save :: TIME_DOOCEAN_step         !< execute ocean component in this step?
  logical,  public, save :: TIME_DOOCEAN_restart      !< execute ocean restart output?
  logical,  public, save :: TIME_DOLAND_step          !< execute land component in this step?
  logical,  public, save :: TIME_DOLAND_restart       !< execute land restart output?
  logical,  public, save :: TIME_DOCPL_calc           !< execute coupler component in this step?
  logical,  public, save :: TIME_DOCPL_restart        !< execute coupler restart output?
  logical,  public, save :: TIME_DOend                !< finish program in this step?

  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  private :: TIME_rapid

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  integer,  private, save :: TIME_STARTDATE(6) = (/ 0, 1, 1, 0, 0, 0 /)
  real(DP), private, save :: TIME_STARTMS      = 0.0_DP !< [millisec]
  integer,  private, save :: TIME_STARTDAY
  real(DP), private, save :: TIME_STARTSEC

  integer,  private, save :: TIME_ENDDATE(6)
  real(DP), private, save :: TIME_ENDMS
  integer,  private, save :: TIME_ENDDAY
  real(DP), private, save :: TIME_ENDSEC

  integer,  private, save :: TIME_NSTEP

  real(DP), private, save :: TIME_RES_ATMOS_DYN     = 0.0_DP
  real(DP), private, save :: TIME_RES_ATMOS_PHY_SF  = 0.0_DP
  real(DP), private, save :: TIME_RES_ATMOS_PHY_TB  = 0.0_DP
  real(DP), private, save :: TIME_RES_ATMOS_PHY_MP  = 0.0_DP
  real(DP), private, save :: TIME_RES_ATMOS_PHY_RD  = 0.0_DP
  real(DP), private, save :: TIME_RES_ATMOS_RESTART = 0.0_DP
  real(DP), private, save :: TIME_RES_OCEAN         = 0.0_DP
  real(DP), private, save :: TIME_RES_OCEAN_RESTART = 0.0_DP
  real(DP), private, save :: TIME_RES_LAND          = 0.0_DP
  real(DP), private, save :: TIME_RES_LAND_RESTART  = 0.0_DP
  real(DP), private, save :: TIME_RES_CPL           = 0.0_DP
  real(DP), private, save :: TIME_RES_CPL_RESTART   = 0.0_DP

  integer,                  private, parameter :: TIME_rapnlimit = 100
  integer,                  private,      save :: TIME_rapnmax   = 0
  character(len=IO_SYSCHR), private,      save :: TIME_rapname(TIME_rapnlimit)
  real(DP),                 private,      save :: TIME_raptstr(TIME_rapnlimit)
  real(DP),                 private,      save :: TIME_rapttot(TIME_rapnlimit)
  integer,                  private,      save :: TIME_rapnstr(TIME_rapnlimit)
  integer,                  private,      save :: TIME_rapnend(TIME_rapnlimit)

  real(DP), private, parameter :: eps = 1.E-10_DP !> epsilon for timesec

  !-----------------------------------------------------------------------------
contains

  !-----------------------------------------------------------------------------
  !> Setup
  subroutine TIME_setup( &
       flg_init ) ! (in)
    use mod_stdio, only: &
       IO_FID_CONF
    use mod_process, only: &
       PRC_MPIstop
    use mod_const, only: &
       UNDEF => CONST_UNDEF
    use mod_calendar, only: &
       CALENDAR_date2daysec,    &
       CALENDAR_daysec2date,    &
       CALENDAR_adjust_daysec,  &
       CALENDAR_combine_daysec, &
       CALENDAR_unit2sec
    implicit none

    logical, intent(in), optional :: flg_init
    logical :: flgi = .false.

    real(DP)                 :: TIME_DURATION
    character(len=IO_SYSCHR) :: TIME_DURATION_UNIT         = "SEC"
    real(DP)                 :: TIME_DT
    character(len=IO_SYSCHR) :: TIME_DT_UNIT               = "SEC"

    real(DP)                 :: TIME_DT_ATMOS_DYN
    character(len=IO_SYSCHR) :: TIME_DT_ATMOS_DYN_UNIT     = "SEC"
    real(DP)                 :: TIME_DT_ATMOS_PHY_SF
    character(len=IO_SYSCHR) :: TIME_DT_ATMOS_PHY_SF_UNIT  = ""
    real(DP)                 :: TIME_DT_ATMOS_PHY_TB
    character(len=IO_SYSCHR) :: TIME_DT_ATMOS_PHY_TB_UNIT  = ""
    real(DP)                 :: TIME_DT_ATMOS_PHY_MP
    character(len=IO_SYSCHR) :: TIME_DT_ATMOS_PHY_MP_UNIT  = ""
    real(DP)                 :: TIME_DT_ATMOS_PHY_RD
    character(len=IO_SYSCHR) :: TIME_DT_ATMOS_PHY_RD_UNIT  = ""
    real(DP)                 :: TIME_DT_ATMOS_RESTART
    character(len=IO_SYSCHR) :: TIME_DT_ATMOS_RESTART_UNIT = ""

    real(DP)                 :: TIME_DT_OCEAN
    character(len=IO_SYSCHR) :: TIME_DT_OCEAN_UNIT         = ""
    real(DP)                 :: TIME_DT_OCEAN_RESTART
    character(len=IO_SYSCHR) :: TIME_DT_OCEAN_RESTART_UNIT = ""

    real(DP)                 :: TIME_DT_LAND
    character(len=IO_SYSCHR) :: TIME_DT_LAND_UNIT          = ""
    real(DP)                 :: TIME_DT_LAND_RESTART
    character(len=IO_SYSCHR) :: TIME_DT_LAND_RESTART_UNIT  = ""

    real(DP)                 :: TIME_DT_CPL
    character(len=IO_SYSCHR) :: TIME_DT_CPL_UNIT           = ""
    real(DP)                 :: TIME_DT_CPL_RESTART
    character(len=IO_SYSCHR) :: TIME_DT_CPL_RESTART_UNIT   = ""

    NAMELIST / PARAM_TIME / &
       TIME_STARTDATE,             &
       TIME_STARTMS,               &
       TIME_DURATION,              &
       TIME_DURATION_UNIT,         &
       TIME_DT,                    &
       TIME_DT_UNIT,               &
       TIME_DT_ATMOS_DYN,          &
       TIME_DT_ATMOS_DYN_UNIT,     &
       TIME_DT_ATMOS_PHY_SF,       &
       TIME_DT_ATMOS_PHY_SF_UNIT,  &
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
       TIME_DT_OCEAN_RESTART_UNIT, &
       TIME_DT_LAND,               &
       TIME_DT_LAND_UNIT,          &
       TIME_DT_LAND_RESTART,       &
       TIME_DT_LAND_RESTART_UNIT,  &
       TIME_DT_CPL,                &
       TIME_DT_CPL_UNIT,           &
       TIME_DT_CPL_RESTART,        &
       TIME_DT_CPL_RESTART_UNIT

    real(DP) :: TIME_DURATIONSEC

    integer :: ierr
    !---------------------------------------------------------------------------

    if ( present(flg_init) ) flgi = flg_init

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[TIME]/Categ[COMMON]'

    TIME_DURATION         = UNDEF
    TIME_DT               = UNDEF
    TIME_DT_ATMOS_DYN     = UNDEF
    TIME_DT_ATMOS_PHY_SF  = UNDEF
    TIME_DT_ATMOS_PHY_TB  = UNDEF
    TIME_DT_ATMOS_PHY_MP  = UNDEF
    TIME_DT_ATMOS_PHY_RD  = UNDEF
    TIME_DT_ATMOS_RESTART = UNDEF
    TIME_DT_OCEAN         = UNDEF
    TIME_DT_OCEAN_RESTART = UNDEF
    TIME_DT_LAND          = UNDEF
    TIME_DT_LAND_RESTART  = UNDEF
    TIME_DT_CPL           = UNDEF
    TIME_DT_CPL_RESTART   = UNDEF

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


    if ( .not. flgi ) then

       ! check time setting
       if ( TIME_DT == UNDEF ) then
          if(IO_L) write(IO_FID_LOG,*) 'xxx Not found TIME_DT.'
          call PRC_MPIstop
       end if

       if ( TIME_DURATION == UNDEF ) then
          if(IO_L) write(IO_FID_LOG,*) 'xxx Not found TIME_DURATION.'
          call PRC_MPIstop
       end if

       ! DYN
       if ( TIME_DT_ATMOS_DYN == UNDEF ) then
          if(IO_L) write(IO_FID_LOG,*) '*** Not found TIME_DT_ATMOS_DYN. TIME_DT is used.'
          TIME_DT_ATMOS_DYN = TIME_DT
       end if
       if ( TIME_DT_ATMOS_DYN_UNIT == '' ) then
          if(IO_L) write(IO_FID_LOG,*) '*** Not found TIME_DT_ATMOS_DYN_UNIT. TIME_DT_UNIT is used.'
          TIME_DT_ATMOS_DYN_UNIT = TIME_DT_UNIT
       end if
       ! PHY_SF
       if ( TIME_DT_ATMOS_PHY_SF == UNDEF ) then
          if(IO_L) write(IO_FID_LOG,*) '*** Not found TIME_DT_ATMOS_PHY_SF. TIME_DT is used.'
          TIME_DT_ATMOS_PHY_SF = TIME_DT
       end if
       if ( TIME_DT_ATMOS_PHY_SF_UNIT == '' ) then
          if(IO_L) write(IO_FID_LOG,*) '*** Not found TIME_DT_ATMOS_PHY_SF_UNIT. TIME_DT_UNIT is used.'
          TIME_DT_ATMOS_PHY_SF_UNIT = TIME_DT_UNIT
       end if
       ! PHY_TB
       if ( TIME_DT_ATMOS_PHY_TB == UNDEF ) then
          if(IO_L) write(IO_FID_LOG,*) '*** Not found TIME_DT_ATMOS_PHY_TB. TIME_DT is used.'
          TIME_DT_ATMOS_PHY_TB = TIME_DT
       end if
       if ( TIME_DT_ATMOS_PHY_TB_UNIT == '' ) then
          if(IO_L) write(IO_FID_LOG,*) '*** Not found TIME_DT_ATMOS_PHY_TB_UNIT. TIME_DT_UNIT is used.'
          TIME_DT_ATMOS_PHY_TB_UNIT = TIME_DT_UNIT
       end if
       ! PHY_MP
       if ( TIME_DT_ATMOS_PHY_MP == UNDEF ) then
          if(IO_L) write(IO_FID_LOG,*) '*** Not found TIME_DT_ATMOS_PHY_MP. TIME_DT is used.'
          TIME_DT_ATMOS_PHY_MP = TIME_DT
       end if
       if ( TIME_DT_ATMOS_PHY_MP_UNIT == '' ) then
          if(IO_L) write(IO_FID_LOG,*) '*** Not found TIME_DT_ATMOS_PHY_MP_UNIT. TIME_DT_UNIT is used.'
          TIME_DT_ATMOS_PHY_MP_UNIT = TIME_DT_UNIT
       end if
       ! PHY_RD
       if ( TIME_DT_ATMOS_PHY_RD == UNDEF ) then
          if(IO_L) write(IO_FID_LOG,*) '*** Not found TIME_DT_ATMOS_PHY_RD. TIME_DT is used.'
          TIME_DT_ATMOS_PHY_RD = TIME_DT
       end if
       if ( TIME_DT_ATMOS_PHY_RD_UNIT == '' ) then
          if(IO_L) write(IO_FID_LOG,*) '*** Not found TIME_DT_ATMOS_PHY_RD_UNIT. TIME_DT_UNIT is used.'
          TIME_DT_ATMOS_PHY_RD_UNIT = TIME_DT_UNIT
       end if
       ! ATMOS RESTART
       if ( TIME_DT_ATMOS_RESTART == UNDEF ) then
          if(IO_L) write(IO_FID_LOG,*) '*** Not found TIME_DT_ATMOS_RESTART. TIME_DURATION is used.'
          TIME_DT_ATMOS_RESTART = TIME_DURATION
       end if
       if ( TIME_DT_ATMOS_RESTART_UNIT == '' ) then
          if(IO_L) write(IO_FID_LOG,*) '*** Not found TIME_DT_ATMOS_RESTART_UNIT. TIME_DURATION_UNIT is used.'
          TIME_DT_ATMOS_RESTART_UNIT = TIME_DURATION_UNIT
       end if
       ! OCEAN
       if ( TIME_DT_OCEAN == UNDEF ) then
          if(IO_L) write(IO_FID_LOG,*) '*** Not found TIME_DT_OCEAN. TIME_DT is used.'
          TIME_DT_OCEAN = TIME_DT
       end if
       if ( TIME_DT_OCEAN_UNIT == '' ) then
          if(IO_L) write(IO_FID_LOG,*) '*** Not found TIME_DT_OCEAN_UNIT. TIME_DT_UNIT is used.'
          TIME_DT_OCEAN_UNIT = TIME_DT_UNIT
       end if
       ! OCEAN RESTART
       if ( TIME_DT_OCEAN_RESTART == UNDEF ) then
          if(IO_L) write(IO_FID_LOG,*) '*** Not found TIME_DT_OCEAN_RESTART. TIME_DURATION is used.'
          TIME_DT_OCEAN_RESTART = TIME_DURATION
       end if
       if ( TIME_DT_OCEAN_RESTART_UNIT == '' ) then
          if(IO_L) write(IO_FID_LOG,*) '*** Not found TIME_DT_OCEAN_RESTART_UNIT. TIME_DURATION_UNIT is used.'
          TIME_DT_OCEAN_RESTART_UNIT = TIME_DURATION_UNIT
       end if
       ! LAND
       if ( TIME_DT_LAND == UNDEF ) then
          if(IO_L) write(IO_FID_LOG,*) '*** Not found TIME_DT_LAND. TIME_DT is used.'
          TIME_DT_LAND = TIME_DT
       end if
       if ( TIME_DT_LAND_UNIT == '' ) then
          if(IO_L) write(IO_FID_LOG,*) '*** Not found TIME_DT_LAND_UNIT. TIME_DT_UNIT is used.'
          TIME_DT_LAND_UNIT = TIME_DT_UNIT
       end if
       ! LAND RESTART
       if ( TIME_DT_LAND_RESTART == UNDEF ) then
          if(IO_L) write(IO_FID_LOG,*) '*** Not found TIME_DT_LAND_RESTART. TIME_DURATION is used.'
          TIME_DT_LAND_RESTART = TIME_DURATION
       end if
       if ( TIME_DT_LAND_RESTART_UNIT == '' ) then
          if(IO_L) write(IO_FID_LOG,*) '*** Not found TIME_DT_LAND_RESTART_UNIT. TIME_DURATION_UNIT is used.'
          TIME_DT_LAND_RESTART_UNIT = TIME_DURATION_UNIT
       end if
       ! COUPLER
       if ( TIME_DT_CPL == UNDEF ) then
          if(IO_L) write(IO_FID_LOG,*) '*** Not found TIME_DT_CPL. TIME_DT is used.'
          TIME_DT_CPL = TIME_DT
       end if
       if ( TIME_DT_CPL_UNIT == '' ) then
          if(IO_L) write(IO_FID_LOG,*) '*** Not found TIME_DT_CPL_UNIT. TIME_DT_UNIT is used.'
          TIME_DT_CPL_UNIT = TIME_DT_UNIT
       end if
       ! CPL RESTART
       if ( TIME_DT_CPL_RESTART == UNDEF ) then
          if(IO_L) write(IO_FID_LOG,*) '*** Not found TIME_DT_CPL_RESTART. TIME_DURATION is used.'
          TIME_DT_CPL_RESTART = TIME_DURATION
       end if
       if ( TIME_DT_CPL_RESTART_UNIT == '' ) then
          if(IO_L) write(IO_FID_LOG,*) '*** Not found TIME_DT_CPL_RESTART_UNIT. TIME_DURATION_UNIT is used.'
          TIME_DT_CPL_RESTART_UNIT = TIME_DURATION_UNIT
       end if

    end if

    !--- calculate time
    TIME_STARTMS = TIME_STARTMS * 1.E-3_DP
    call CALENDAR_date2daysec( TIME_STARTDAY,     & ! [OUT]
                               TIME_STARTSEC,     & ! [OUT]
                               TIME_STARTDATE(:), & ! [IN]
                               TIME_STARTMS       ) ! [IN]

    TIME_NOWDATE(:) = TIME_STARTDATE(:)
    TIME_NOWMS      = TIME_STARTMS
    TIME_NOWDAY     = TIME_STARTDAY
    TIME_NOWSEC     = TIME_STARTSEC
    TIME_NOWDAYSEC  = CALENDAR_combine_daysec( TIME_NOWDAY, TIME_NOWSEC )

    TIME_ENDDAY = TIME_STARTDAY

    if ( flgi ) then
       TIME_ENDSEC = TIME_STARTSEC
    else
       call CALENDAR_unit2sec( TIME_DURATIONSEC, TIME_DURATION, TIME_DURATION_UNIT )
       TIME_ENDSEC = TIME_STARTSEC + TIME_DURATIONSEC
    end if

    call CALENDAR_adjust_daysec( TIME_ENDDAY, TIME_ENDSEC ) ! [INOUT]

    call CALENDAR_daysec2date( TIME_ENDDATE(:), & ! [OUT]
                               TIME_ENDMS,      & ! [OUT]
                               TIME_ENDDAY,     & ! [IN]
                               TIME_ENDSEC      ) ! [IN]

    if ( .not. flgi ) then

       call CALENDAR_unit2sec( TIME_DTSEC, TIME_DT, TIME_DT_UNIT )

       TIME_NSTEP   = int( TIME_DURATIONSEC / TIME_DTSEC )
       TIME_NOWSTEP = 1

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

    !if( IO_L ) write(IO_FID_LOG,*) '*** now:', TIME_NOWDAYSEC, TIME_NOWDAY, TIME_NOWSEC

       !--- calculate intervals for atmosphere
       call CALENDAR_unit2sec( TIME_DTSEC_ATMOS_DYN,     TIME_DT_ATMOS_DYN,     TIME_DT_ATMOS_DYN_UNIT     )
       call CALENDAR_unit2sec( TIME_DTSEC_ATMOS_PHY_SF,  TIME_DT_ATMOS_PHY_SF,  TIME_DT_ATMOS_PHY_SF_UNIT  )
       call CALENDAR_unit2sec( TIME_DTSEC_ATMOS_PHY_TB,  TIME_DT_ATMOS_PHY_TB,  TIME_DT_ATMOS_PHY_TB_UNIT  )
       call CALENDAR_unit2sec( TIME_DTSEC_ATMOS_PHY_MP,  TIME_DT_ATMOS_PHY_MP,  TIME_DT_ATMOS_PHY_MP_UNIT  )
       call CALENDAR_unit2sec( TIME_DTSEC_ATMOS_PHY_RD,  TIME_DT_ATMOS_PHY_RD,  TIME_DT_ATMOS_PHY_RD_UNIT  )
       call CALENDAR_unit2sec( TIME_DTSEC_ATMOS_RESTART, TIME_DT_ATMOS_RESTART, TIME_DT_ATMOS_RESTART_UNIT )
       call CALENDAR_unit2sec( TIME_DTSEC_OCEAN,         TIME_DT_OCEAN,         TIME_DT_OCEAN_UNIT         )
       call CALENDAR_unit2sec( TIME_DTSEC_OCEAN_RESTART, TIME_DT_OCEAN_RESTART, TIME_DT_OCEAN_RESTART_UNIT )
       call CALENDAR_unit2sec( TIME_DTSEC_LAND,          TIME_DT_LAND,          TIME_DT_LAND_UNIT          )
       call CALENDAR_unit2sec( TIME_DTSEC_LAND_RESTART,  TIME_DT_LAND_RESTART,  TIME_DT_LAND_RESTART_UNIT  )
       call CALENDAR_unit2sec( TIME_DTSEC_CPL,           TIME_DT_CPL,           TIME_DT_CPL_UNIT           )
       call CALENDAR_unit2sec( TIME_DTSEC_CPL_RESTART,   TIME_DT_CPL_RESTART,   TIME_DT_CPL_RESTART_UNIT   )

       TIME_NSTEP_ATMOS_DYN = max( int( TIME_DTSEC / TIME_DTSEC_ATMOS_DYN ), 1 )

       TIME_DTSEC_ATMOS_DYN     = max( TIME_DTSEC_ATMOS_DYN,     TIME_DTSEC          /TIME_NSTEP_ATMOS_DYN )
       TIME_DTSEC_ATMOS_PHY_SF  = max( TIME_DTSEC_ATMOS_PHY_SF,  TIME_DTSEC_ATMOS_DYN*TIME_NSTEP_ATMOS_DYN )
       TIME_DTSEC_ATMOS_PHY_TB  = max( TIME_DTSEC_ATMOS_PHY_TB,  TIME_DTSEC_ATMOS_DYN*TIME_NSTEP_ATMOS_DYN )
       TIME_DTSEC_ATMOS_PHY_MP  = max( TIME_DTSEC_ATMOS_PHY_MP,  TIME_DTSEC_ATMOS_DYN*TIME_NSTEP_ATMOS_DYN )
       TIME_DTSEC_ATMOS_PHY_RD  = max( TIME_DTSEC_ATMOS_PHY_RD,  TIME_DTSEC_ATMOS_DYN*TIME_NSTEP_ATMOS_DYN )
       TIME_DTSEC_ATMOS_RESTART = max( TIME_DTSEC_ATMOS_RESTART, TIME_DTSEC_ATMOS_DYN*TIME_NSTEP_ATMOS_DYN )
       TIME_DTSEC_OCEAN         = max( TIME_DTSEC_OCEAN,         TIME_DTSEC_ATMOS_DYN*TIME_NSTEP_ATMOS_DYN )
       TIME_DTSEC_OCEAN_RESTART = max( TIME_DTSEC_OCEAN_RESTART, TIME_DTSEC_ATMOS_DYN*TIME_NSTEP_ATMOS_DYN )
       TIME_DTSEC_LAND          = max( TIME_DTSEC_LAND,          TIME_DTSEC_ATMOS_DYN*TIME_NSTEP_ATMOS_DYN )
       TIME_DTSEC_LAND_RESTART  = max( TIME_DTSEC_LAND_RESTART,  TIME_DTSEC_ATMOS_DYN*TIME_NSTEP_ATMOS_DYN )
       TIME_DTSEC_CPL           = max( TIME_DTSEC_CPL,           TIME_DTSEC_ATMOS_DYN*TIME_NSTEP_ATMOS_DYN )
       TIME_DTSEC_CPL_RESTART   = max( TIME_DTSEC_CPL_RESTART,   TIME_DTSEC_ATMOS_DYN*TIME_NSTEP_ATMOS_DYN )

       if( IO_L ) write(IO_FID_LOG,*)
       if( IO_L ) write(IO_FID_LOG,*) '*** Time interval for atmospheric processes (sec.)'
       if( IO_L ) write(IO_FID_LOG,'(1x,A,F10.3)') '*** Dynamics (time)             : ', TIME_DTSEC_ATMOS_DYN
       if( IO_L ) write(IO_FID_LOG,'(1x,A,I4)')    '***          (step)             : ', TIME_NSTEP_ATMOS_DYN
       if( IO_L ) write(IO_FID_LOG,'(1x,A,F10.3)') '*** Physics, Surface Flux       : ', TIME_DTSEC_ATMOS_PHY_SF
       if( IO_L ) write(IO_FID_LOG,'(1x,A,F10.3)') '*** Physics, Turbulence         : ', TIME_DTSEC_ATMOS_PHY_TB
       if( IO_L ) write(IO_FID_LOG,'(1x,A,F10.3)') '*** Physics, Cloud Microphysics : ', TIME_DTSEC_ATMOS_PHY_MP
       if( IO_L ) write(IO_FID_LOG,'(1x,A,F10.3)') '*** Physics, Radiation          : ', TIME_DTSEC_ATMOS_PHY_RD
       if( IO_L ) write(IO_FID_LOG,*) '*** Time interval for ocean process (sec.)'
       if( IO_L ) write(IO_FID_LOG,'(1x,A,F10.3)') '*** Ocean update                : ', TIME_DTSEC_OCEAN
       if( IO_L ) write(IO_FID_LOG,*) '*** Time interval for land process (sec.)'
       if( IO_L ) write(IO_FID_LOG,'(1x,A,F10.3)') '*** Land  update                : ', TIME_DTSEC_LAND
       if( IO_L ) write(IO_FID_LOG,*) '*** Time interval for coupler process (sec.)'
       if( IO_L ) write(IO_FID_LOG,'(1x,A,F10.3)') '*** Coupler interval            : ', TIME_DTSEC_CPL
       if( IO_L ) write(IO_FID_LOG,*) '*** Time interval for Restart (sec.)'
       if( IO_L ) write(IO_FID_LOG,'(1x,A,F10.3)') '*** Atmospheric Variables       : ', TIME_DTSEC_ATMOS_RESTART
       if( IO_L ) write(IO_FID_LOG,'(1x,A,F10.3)') '*** Ocean Variables             : ', TIME_DTSEC_OCEAN_RESTART
       if( IO_L ) write(IO_FID_LOG,'(1x,A,F10.3)') '*** Land Variables              : ', TIME_DTSEC_LAND_RESTART
       if( IO_L ) write(IO_FID_LOG,'(1x,A,F10.3)') '*** Coupler Variables           : ', TIME_DTSEC_CPL_RESTART

    end if

    return
  end subroutine TIME_setup

  !-----------------------------------------------------------------------------
  !> Evaluate component execution
  subroutine TIME_checkstate
    implicit none
    !---------------------------------------------------------------------------

    TIME_DOATMOS_step      = .false.
    TIME_DOATMOS_DYN       = .false.
    TIME_DOATMOS_PHY_SF    = .false.
    TIME_DOATMOS_PHY_TB    = .false.
    TIME_DOATMOS_PHY_MP    = .false.
    TIME_DOATMOS_PHY_RD    = .false.
    TIME_DOOCEAN_step      = .false.
    TIME_DOLAND_step       = .false.
    TIME_DOCPL_calc        = .false.

    TIME_RES_ATMOS_DYN    = TIME_RES_ATMOS_DYN    + TIME_DTSEC
    TIME_RES_ATMOS_PHY_SF = TIME_RES_ATMOS_PHY_SF + TIME_DTSEC
    TIME_RES_ATMOS_PHY_TB = TIME_RES_ATMOS_PHY_TB + TIME_DTSEC
    TIME_RES_ATMOS_PHY_MP = TIME_RES_ATMOS_PHY_MP + TIME_DTSEC
    TIME_RES_ATMOS_PHY_RD = TIME_RES_ATMOS_PHY_RD + TIME_DTSEC
    TIME_RES_OCEAN        = TIME_RES_OCEAN        + TIME_DTSEC
    TIME_RES_LAND         = TIME_RES_LAND         + TIME_DTSEC
    TIME_RES_CPL          = TIME_RES_CPL          + TIME_DTSEC

    if ( TIME_RES_ATMOS_DYN - TIME_DTSEC_ATMOS_DYN > -eps ) then
       TIME_DOATMOS_step  = .true.
       TIME_DOATMOS_DYN   = .true.
       TIME_RES_ATMOS_DYN = TIME_RES_ATMOS_DYN - TIME_DTSEC_ATMOS_DYN
    endif
    if ( TIME_RES_ATMOS_PHY_SF - TIME_DTSEC_ATMOS_PHY_SF > -eps ) then
       TIME_DOATMOS_step     = .true.
       TIME_DOATMOS_PHY_SF   = .true.
       TIME_RES_ATMOS_PHY_SF = TIME_RES_ATMOS_PHY_SF - TIME_DTSEC_ATMOS_PHY_SF
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
    if ( TIME_RES_LAND  - TIME_DTSEC_LAND  > -eps ) then
       TIME_DOLAND_step  = .true.
       TIME_RES_LAND     = TIME_RES_LAND  - TIME_DTSEC_LAND
    endif

    if ( TIME_RES_CPL   - TIME_DTSEC_CPL   > -eps ) then
       TIME_DOCPL_calc   = .true.
       TIME_RES_CPL      = TIME_RES_CPL   - TIME_DTSEC_CPL
    endif

    if( IO_L ) write(IO_FID_LOG,'(1x,A,I4.4,A,I2.2,A,I2.2,A,I2.2,A,I2.2,A,I2.2,A,F6.3,A,I6,A,I6)') &
               '*** TIME: ', TIME_NOWDATE(1),'/',TIME_NOWDATE(2),'/',TIME_NOWDATE(3),' ', &
                             TIME_NOWDATE(4),':',TIME_NOWDATE(5),':',TIME_NOWDATE(6),' +', &
                             TIME_NOWMS,' STEP:',TIME_NOWSTEP, '/', TIME_NSTEP

    return
  end subroutine TIME_checkstate

  !-----------------------------------------------------------------------------
  !> Advance the time & evaluate restart & stop
  subroutine TIME_advance
    use mod_process, only: &
       PRC_master, &
       PRC_myrank
    use mod_calendar, only: &
       CALENDAR_daysec2date,   &
       CALENDAR_adjust_daysec, &
       CALENDAR_combine_daysec
    implicit none

    logical :: exists
    !---------------------------------------------------------------------------

    TIME_DOend = .false.

    TIME_NOWSEC = TIME_NOWSEC + TIME_DTSEC

    call CALENDAR_adjust_daysec( TIME_NOWDAY, TIME_NOWSEC ) ! [INOUT]

    TIME_NOWDAYSEC = CALENDAR_combine_daysec( TIME_NOWDAY, TIME_NOWSEC )

    call CALENDAR_daysec2date( TIME_NOWDATE(:), & ! [OUT]
                               TIME_NOWMS,      & ! [OUT]
                               TIME_NOWDAY,     & ! [IN]
                               TIME_NOWSEC      ) ! [IN]

    TIME_NOWSTEP = TIME_NOWSTEP + 1

    if (       TIME_NOWDAY - TIME_ENDDAY >= 0   &
         .AND. TIME_NOWSEC - TIME_ENDSEC > -eps ) then
       TIME_DOend = .true.
    endif

    ! QUIT file control
    if ( PRC_myrank == PRC_master ) then ! master node
       inquire(file='QUIT', exist=exists)

       if( exists ) TIME_DOend = .true.
    endif

    TIME_DOATMOS_restart = .false.
    TIME_DOOCEAN_restart = .false.
    TIME_DOLAND_restart  = .false.
    TIME_DOCPL_restart   = .false.

    TIME_RES_ATMOS_RESTART = TIME_RES_ATMOS_RESTART + TIME_DTSEC
    TIME_RES_OCEAN_RESTART = TIME_RES_OCEAN_RESTART + TIME_DTSEC
    TIME_RES_LAND_RESTART  = TIME_RES_LAND_RESTART  + TIME_DTSEC
    TIME_RES_CPL_RESTART   = TIME_RES_CPL_RESTART   + TIME_DTSEC

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

    if ( TIME_RES_LAND_RESTART  - TIME_DTSEC_LAND_RESTART  > -eps ) then
       TIME_DOLAND_restart    = .true.
       TIME_RES_LAND_RESTART  = TIME_RES_LAND_RESTART  - TIME_DTSEC_LAND_RESTART
    elseif( TIME_DOend ) then
       TIME_DOLAND_restart    = .true.
    endif

    if ( TIME_RES_CPL_RESTART  - TIME_DTSEC_CPL_RESTART  > -eps ) then
       TIME_DOCPL_restart    = .true.
       TIME_RES_CPL_RESTART  = TIME_RES_CPL_RESTART  - TIME_DTSEC_CPL_RESTART
    elseif( TIME_DOend ) then
       TIME_DOCPL_restart    = .true.
    endif

    return
  end subroutine TIME_advance

  !-----------------------------------------------------------------------------
  !> Get item ID or register item
  function TIME_rapid( rapname ) result(id)
    implicit none

    character(len=*), intent(in) :: rapname !< name of item

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

    TIME_rapnstr(id) = 0
    TIME_rapnend(id) = 0
    TIME_rapttot(id) = 0.0_DP

  end function TIME_rapid

  !-----------------------------------------------------------------------------
  !> Start raptime
  subroutine TIME_rapstart( rapname )
    use mod_process, only: &
       PRC_MPItime
    implicit none

    character(len=*), intent(in) :: rapname !< name of item

    integer :: id
    !---------------------------------------------------------------------------

    id = TIME_rapid( rapname )

    TIME_raptstr(id) = PRC_MPItime()
    TIME_rapnstr(id) = TIME_rapnstr(id) + 1

    !if( IO_L ) write(IO_FID_LOG,*) rapname, TIME_rapnstr(id)

#ifdef _FAPP_
call START_COLLECTION( rapname )
#endif

    return
  end subroutine TIME_rapstart

  !-----------------------------------------------------------------------------
  !> Save raptime
  subroutine TIME_rapend( rapname )
    use mod_process, only: &
       PRC_MPItime
    implicit none

    character(len=*), intent(in) :: rapname !< name of item

    integer :: id
    !---------------------------------------------------------------------------

    id = TIME_rapid( rapname )

    TIME_rapttot(id) = TIME_rapttot(id) + ( PRC_MPItime()-TIME_raptstr(id) )
    TIME_rapnend(id) = TIME_rapnend(id) + 1

#ifdef _FAPP_
call STOP_COLLECTION( rapname )
#endif

    return
  end subroutine TIME_rapend

  !-----------------------------------------------------------------------------
  !> Report raptime
  subroutine TIME_rapreport
    use mod_stdio, only: &
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
    integer  :: maxidx(TIME_rapnlimit)
    integer  :: minidx(TIME_rapnlimit)

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

       call PRC_MPItimestat( avgvar      (1:TIME_rapnmax), &
                             maxvar      (1:TIME_rapnmax), &
                             minvar      (1:TIME_rapnmax), &
                             maxidx      (1:TIME_rapnmax), &
                             minidx      (1:TIME_rapnmax), &
                             TIME_rapttot(1:TIME_rapnmax)  )

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
!-------------------------------------------------------------------------------
