!-------------------------------------------------------------------------------
!> module ADMIN TIME
!!
!! @par Description
!!          Time management for SCALE-RM
!!
!! @author Team SCALE
!!
!<
module mod_admin_time
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_stdio
  use scale_prof
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: ADMIN_TIME_setup
  public :: ADMIN_TIME_checkstate
  public :: ADMIN_TIME_advance

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  real(DP), public :: TIME_DTSEC_ATMOS_RESTART  !< time interval of atmosphere restart [sec]
  real(DP), public :: TIME_DTSEC_OCEAN_RESTART  !< time interval of ocean restart      [sec]
  real(DP), public :: TIME_DTSEC_LAND_RESTART   !< time interval of land restart       [sec]
  real(DP), public :: TIME_DTSEC_URBAN_RESTART  !< time interval of urban restart      [sec]
  real(DP), public :: TIME_DTSEC_RESUME         !< time interval for resume            [sec]

  integer,  public :: TIME_DSTEP_ATMOS_RESTART  !< interval of atmosphere restart [step]
  integer,  public :: TIME_DSTEP_OCEAN_RESTART  !< interval of ocean restart      [step]
  integer,  public :: TIME_DSTEP_LAND_RESTART   !< interval of land restart       [step]
  integer,  public :: TIME_DSTEP_URBAN_RESTART  !< interval of urban restart      [step]
  integer,  public :: TIME_DSTEP_RESUME         !< interval for resume            [step]

  logical,  public :: TIME_DOATMOS_step         !< execute atmosphere component in this step?
  logical,  public :: TIME_DOATMOS_DYN          !< execute dynamics in this step?
  logical,  public :: TIME_DOATMOS_PHY_CP       !< execute physics  in this step? (cumulus     )
  logical,  public :: TIME_DOATMOS_PHY_MP       !< execute physics  in this step? (microphysics)
  logical,  public :: TIME_DOATMOS_PHY_RD       !< execute physics  in this step? (radiation   )
  logical,  public :: TIME_DOATMOS_PHY_SF       !< execute physics  in this step? (surface flux)
  logical,  public :: TIME_DOATMOS_PHY_TB       !< execute physics  in this step? (turbulence  )
  logical,  public :: TIME_DOATMOS_PHY_CH       !< execute physics  in this step? (chemistry   )
  logical,  public :: TIME_DOATMOS_PHY_AE       !< execute physics  in this step? (aerosol     )
  logical,  public :: TIME_DOATMOS_restart      !< execute atmosphere restart output in this step?
  logical,  public :: TIME_DOOCEAN_step         !< execute ocean      component      in this step?
  logical,  public :: TIME_DOOCEAN_restart      !< execute ocean      restart output in this step?
  logical,  public :: TIME_DOLAND_step          !< execute land       component      in this step?
  logical,  public :: TIME_DOLAND_restart       !< execute land       restart output in this step?
  logical,  public :: TIME_DOURBAN_step         !< execute urban      component      in this step?
  logical,  public :: TIME_DOURBAN_restart      !< execute urban      restart output in this step?
  logical,  public :: TIME_DOresume             !< resume in this step?
  logical,  public :: TIME_DOend                !< finish program in this step?

  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  integer,  private :: TIME_STARTDATE(6) = (/ -999, 1, 1, 0, 0, 0 /)
  real(DP), private :: TIME_STARTMS      = 0.0_DP !< [millisec]
  integer,  private :: TIME_STARTDAY
  real(DP), private :: TIME_STARTSEC

  integer,  private :: TIME_ENDDATE(6)
  real(DP), private :: TIME_ENDMS
  integer,  private :: TIME_ENDDAY
  real(DP), private :: TIME_ENDSEC

  integer,  private :: TIME_RES_ATMOS_DYN     = 0
  integer,  private :: TIME_RES_ATMOS_PHY_CP  = 0
  integer,  private :: TIME_RES_ATMOS_PHY_MP  = 0
  integer,  private :: TIME_RES_ATMOS_PHY_RD  = 0
  integer,  private :: TIME_RES_ATMOS_PHY_SF  = 0
  integer,  private :: TIME_RES_ATMOS_PHY_TB  = 0
  integer,  private :: TIME_RES_ATMOS_PHY_CH  = 0
  integer,  private :: TIME_RES_ATMOS_PHY_AE  = 0
  integer,  private :: TIME_RES_ATMOS_RESTART = 0
  integer,  private :: TIME_RES_OCEAN         = 0
  integer,  private :: TIME_RES_OCEAN_RESTART = 0
  integer,  private :: TIME_RES_LAND          = 0
  integer,  private :: TIME_RES_LAND_RESTART  = 0
  integer,  private :: TIME_RES_URBAN         = 0
  integer,  private :: TIME_RES_URBAN_RESTART = 0
  integer,  private :: TIME_RES_RESUME

  real(DP), private :: TIME_WALLCLOCK_START             ! Start time of wall clock             [sec]
  real(DP), private :: TIME_WALLCLOCK_LIMIT   = -1.0_DP ! Elapse time limit of wall clock time [sec]
  real(DP), private :: TIME_WALLCLOCK_SAFE    =  0.9_DP ! Safety coefficient for elapse time limit

  real(DP), private, parameter :: eps = 1.E-6_DP !> epsilon for timesec

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine ADMIN_TIME_setup( &
       setup_TimeIntegration )
    use gtool_file, only: &
       FileGetDatainfo
    use scale_process, only: &
       PRC_myrank,  &
       PRC_MPIstop, &
       PRC_MPItime
    use scale_const, only: &
       UNDEF8 => CONST_UNDEF8
    use scale_calendar, only: &
       CALENDAR_date2daysec,    &
       CALENDAR_daysec2date,    &
       CALENDAR_adjust_daysec,  &
       CALENDAR_combine_daysec, &
       CALENDAR_unit2sec,       &
       CALENDAR_CFunits2sec,    &
       CALENDAR_date2char
    use scale_time, only: &
       TIME_DTSEC,                 &
       TIME_NOWDATE,               &
       TIME_NOWMS,                 &
       TIME_NOWDAY,                &
       TIME_NOWSEC,                &
       TIME_NOWDAYSEC,             &
       TIME_NOWSTEP,               &
       TIME_NSTEP,                 &
       TIME_DTSEC_ATMOS_DYN,       &
       NSTEP_DYN => TIME_NSTEP_ATMOS_DYN, &
       TIME_DTSEC_ATMOS_PHY_CP,    &
       TIME_DTSEC_ATMOS_PHY_MP,    &
       TIME_DTSEC_ATMOS_PHY_RD,    &
       TIME_DTSEC_ATMOS_PHY_SF,    &
       TIME_DTSEC_ATMOS_PHY_TB,    &
       TIME_DTSEC_ATMOS_PHY_CH,    &
       TIME_DTSEC_ATMOS_PHY_AE,    &
       TIME_DTSEC_OCEAN,           &
       TIME_DTSEC_LAND,            &
       TIME_DTSEC_URBAN,           &
       TIME_DTSEC_WALLCLOCK_CHECK, &
       TIME_DSTEP_ATMOS_DYN,       &
       TIME_DSTEP_ATMOS_PHY_CP,    &
       TIME_DSTEP_ATMOS_PHY_MP,    &
       TIME_DSTEP_ATMOS_PHY_RD,    &
       TIME_DSTEP_ATMOS_PHY_SF,    &
       TIME_DSTEP_ATMOS_PHY_TB,    &
       TIME_DSTEP_ATMOS_PHY_CH,    &
       TIME_DSTEP_ATMOS_PHY_AE,    &
       TIME_DSTEP_OCEAN,           &
       TIME_DSTEP_LAND,            &
       TIME_DSTEP_URBAN,           &
       TIME_DSTEP_WALLCLOCK_CHECK, &
       TIME_OFFSET_YEAR,           &
       TIME_STARTDAYSEC
    use mod_atmos_vars, only: &
       RESTART_IN_BASENAME => ATMOS_RESTART_IN_BASENAME
    implicit none

    logical, intent(in) :: setup_TimeIntegration

    real(DP)               :: TIME_DURATION                = UNDEF8
    character(len=H_SHORT) :: TIME_DURATION_UNIT           = "SEC"
    real(DP)               :: TIME_DT                      = UNDEF8
    character(len=H_SHORT) :: TIME_DT_UNIT                 = "SEC"

    real(DP)               :: TIME_DT_ATMOS_DYN            = UNDEF8
    character(len=H_SHORT) :: TIME_DT_ATMOS_DYN_UNIT       = "SEC"
    integer                :: TIME_NSTEP_ATMOS_DYN         = -1
    real(DP)               :: TIME_DT_ATMOS_PHY_CP         = UNDEF8
    character(len=H_SHORT) :: TIME_DT_ATMOS_PHY_CP_UNIT    = ""
    real(DP)               :: TIME_DT_ATMOS_PHY_MP         = UNDEF8
    character(len=H_SHORT) :: TIME_DT_ATMOS_PHY_MP_UNIT    = ""
    real(DP)               :: TIME_DT_ATMOS_PHY_RD         = UNDEF8
    character(len=H_SHORT) :: TIME_DT_ATMOS_PHY_RD_UNIT    = ""
    real(DP)               :: TIME_DT_ATMOS_PHY_SF         = UNDEF8
    character(len=H_SHORT) :: TIME_DT_ATMOS_PHY_SF_UNIT    = ""
    real(DP)               :: TIME_DT_ATMOS_PHY_TB         = UNDEF8
    character(len=H_SHORT) :: TIME_DT_ATMOS_PHY_TB_UNIT    = ""
    real(DP)               :: TIME_DT_ATMOS_PHY_CH         = UNDEF8
    character(len=H_SHORT) :: TIME_DT_ATMOS_PHY_CH_UNIT    = ""
    real(DP)               :: TIME_DT_ATMOS_PHY_AE         = UNDEF8
    character(len=H_SHORT) :: TIME_DT_ATMOS_PHY_AE_UNIT    = ""
    real(DP)               :: TIME_DT_ATMOS_RESTART        = UNDEF8
    character(len=H_SHORT) :: TIME_DT_ATMOS_RESTART_UNIT   = ""

    real(DP)               :: TIME_DT_OCEAN                = UNDEF8
    character(len=H_SHORT) :: TIME_DT_OCEAN_UNIT           = ""
    real(DP)               :: TIME_DT_OCEAN_RESTART        = UNDEF8
    character(len=H_SHORT) :: TIME_DT_OCEAN_RESTART_UNIT   = ""

    real(DP)               :: TIME_DT_LAND                 = UNDEF8
    character(len=H_SHORT) :: TIME_DT_LAND_UNIT            = ""
    real(DP)               :: TIME_DT_LAND_RESTART         = UNDEF8
    character(len=H_SHORT) :: TIME_DT_LAND_RESTART_UNIT    = ""

    real(DP)               :: TIME_DT_URBAN                = UNDEF8
    character(len=H_SHORT) :: TIME_DT_URBAN_UNIT           = ""
    real(DP)               :: TIME_DT_URBAN_RESTART        = UNDEF8
    character(len=H_SHORT) :: TIME_DT_URBAN_RESTART_UNIT   = ""

    real(DP)               :: TIME_DT_RESUME               = UNDEF8
    character(len=H_SHORT) :: TIME_DT_RESUME_UNIT          = ""

    real(DP)               :: TIME_DT_WALLCLOCK_CHECK      = UNDEF8
    character(len=H_SHORT) :: TIME_DT_WALLCLOCK_CHECK_UNIT = ""

    NAMELIST / PARAM_TIME / &
       TIME_STARTDATE,               &
       TIME_STARTMS,                 &
       TIME_DURATION,                &
       TIME_DURATION_UNIT,           &
       TIME_DT,                      &
       TIME_DT_UNIT,                 &
       TIME_DT_ATMOS_DYN,            &
       TIME_DT_ATMOS_DYN_UNIT,       &
       TIME_NSTEP_ATMOS_DYN,         &
       TIME_DT_ATMOS_PHY_CP,         &
       TIME_DT_ATMOS_PHY_CP_UNIT,    &
       TIME_DT_ATMOS_PHY_MP,         &
       TIME_DT_ATMOS_PHY_MP_UNIT,    &
       TIME_DT_ATMOS_PHY_RD,         &
       TIME_DT_ATMOS_PHY_RD_UNIT,    &
       TIME_DT_ATMOS_PHY_SF,         &
       TIME_DT_ATMOS_PHY_SF_UNIT,    &
       TIME_DT_ATMOS_PHY_TB,         &
       TIME_DT_ATMOS_PHY_TB_UNIT,    &
       TIME_DT_ATMOS_PHY_CH,         &
       TIME_DT_ATMOS_PHY_CH_UNIT,    &
       TIME_DT_ATMOS_PHY_AE,         &
       TIME_DT_ATMOS_PHY_AE_UNIT,    &
       TIME_DT_ATMOS_RESTART,        &
       TIME_DT_ATMOS_RESTART_UNIT,   &
       TIME_DT_OCEAN,                &
       TIME_DT_OCEAN_UNIT,           &
       TIME_DT_OCEAN_RESTART,        &
       TIME_DT_OCEAN_RESTART_UNIT,   &
       TIME_DT_LAND,                 &
       TIME_DT_LAND_UNIT,            &
       TIME_DT_LAND_RESTART,         &
       TIME_DT_LAND_RESTART_UNIT,    &
       TIME_DT_URBAN,                &
       TIME_DT_URBAN_UNIT,           &
       TIME_DT_URBAN_RESTART,        &
       TIME_DT_URBAN_RESTART_UNIT,   &
       TIME_DT_WALLCLOCK_CHECK,      &
       TIME_DT_WALLCLOCK_CHECK_UNIT, &
       TIME_DT_RESUME,               &
       TIME_DT_RESUME_UNIT,          &
       TIME_WALLCLOCK_LIMIT,         &
       TIME_WALLCLOCK_SAFE

    integer              :: dateday
    real(DP)             :: datesec
    real(DP)             :: cftime
    character(len=H_MID) :: cfunits

    real(DP)          :: TIME_DURATIONSEC
    character(len=27) :: startchardate
    character(len=27) :: endchardate

    integer :: ierr
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '++++++ Module[TIME] / Categ[COMMON] / Origin[SCALElib]'

    TIME_NSTEP_ATMOS_DYN = -1

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_TIME,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_TIME. Check!'
       call PRC_MPIstop
    endif
    if( IO_LNML ) write(IO_FID_LOG,nml=PARAM_TIME)

    ! check time setting
    if ( setup_TimeIntegration ) then
       if ( TIME_DT == UNDEF8 ) then
          write(*,*) 'xxx Not found TIME_DT.'
          call PRC_MPIstop
       endif
       if ( TIME_DURATION == UNDEF8 ) then
          write(*,*) 'xxx Not found TIME_DURATION.'
          call PRC_MPIstop
       endif

       ! DYN
       if ( TIME_DT_ATMOS_DYN == UNDEF8 ) then
          if ( TIME_NSTEP_ATMOS_DYN < 0 ) then
             if( IO_L ) write(IO_FID_LOG,*) '*** Not found TIME_DT_ATMOS_DYN. TIME_DT is used.'
             TIME_DT_ATMOS_DYN = TIME_DT
          endif
       endif
       if ( TIME_DT_ATMOS_DYN_UNIT == '' ) then
          if( IO_L ) write(IO_FID_LOG,*) '*** Not found TIME_DT_ATMOS_DYN_UNIT. TIME_DT_UNIT is used.'
          TIME_DT_ATMOS_DYN_UNIT = TIME_DT_UNIT
       endif
       ! PHY_CP
       if ( TIME_DT_ATMOS_PHY_CP == UNDEF8 ) then
          if( IO_L ) write(IO_FID_LOG,*) '*** Not found TIME_DT_ATMOS_PHY_CP. TIME_DT is used.'
          TIME_DT_ATMOS_PHY_CP = TIME_DT
       endif
       if ( TIME_DT_ATMOS_PHY_CP_UNIT == '' ) then
          if( IO_L ) write(IO_FID_LOG,*) '*** Not found TIME_DT_ATMOS_PHY_CP_UNIT. TIME_DT_UNIT is used.'
          TIME_DT_ATMOS_PHY_CP_UNIT = TIME_DT_UNIT
       endif
       ! PHY_MP
       if ( TIME_DT_ATMOS_PHY_MP == UNDEF8 ) then
          if( IO_L ) write(IO_FID_LOG,*) '*** Not found TIME_DT_ATMOS_PHY_MP. TIME_DT is used.'
          TIME_DT_ATMOS_PHY_MP = TIME_DT
       endif
       if ( TIME_DT_ATMOS_PHY_MP_UNIT == '' ) then
          if( IO_L ) write(IO_FID_LOG,*) '*** Not found TIME_DT_ATMOS_PHY_MP_UNIT. TIME_DT_UNIT is used.'
          TIME_DT_ATMOS_PHY_MP_UNIT = TIME_DT_UNIT
       endif
       ! PHY_RD
       if ( TIME_DT_ATMOS_PHY_RD == UNDEF8 ) then
          if( IO_L ) write(IO_FID_LOG,*) '*** Not found TIME_DT_ATMOS_PHY_RD. TIME_DT is used.'
          TIME_DT_ATMOS_PHY_RD = TIME_DT
       endif
       if ( TIME_DT_ATMOS_PHY_RD_UNIT == '' ) then
          if( IO_L ) write(IO_FID_LOG,*) '*** Not found TIME_DT_ATMOS_PHY_RD_UNIT. TIME_DT_UNIT is used.'
          TIME_DT_ATMOS_PHY_RD_UNIT = TIME_DT_UNIT
       endif
       ! PHY_SF
       if ( TIME_DT_ATMOS_PHY_SF == UNDEF8 ) then
          if( IO_L ) write(IO_FID_LOG,*) '*** Not found TIME_DT_ATMOS_PHY_SF. TIME_DT is used.'
          TIME_DT_ATMOS_PHY_SF = TIME_DT
       endif
       if ( TIME_DT_ATMOS_PHY_SF_UNIT == '' ) then
          if( IO_L ) write(IO_FID_LOG,*) '*** Not found TIME_DT_ATMOS_PHY_SF_UNIT. TIME_DT_UNIT is used.'
          TIME_DT_ATMOS_PHY_SF_UNIT = TIME_DT_UNIT
       endif
       ! PHY_TB
       if ( TIME_DT_ATMOS_PHY_TB == UNDEF8 ) then
          if( IO_L ) write(IO_FID_LOG,*) '*** Not found TIME_DT_ATMOS_PHY_TB. TIME_DT is used.'
          TIME_DT_ATMOS_PHY_TB = TIME_DT
       endif
       if ( TIME_DT_ATMOS_PHY_TB_UNIT == '' ) then
          if( IO_L ) write(IO_FID_LOG,*) '*** Not found TIME_DT_ATMOS_PHY_TB_UNIT. TIME_DT_UNIT is used.'
          TIME_DT_ATMOS_PHY_TB_UNIT = TIME_DT_UNIT
       endif
       ! PHY_CH
       if ( TIME_DT_ATMOS_PHY_CH == UNDEF8 ) then
          if( IO_L ) write(IO_FID_LOG,*) '*** Not found TIME_DT_ATMOS_PHY_CH. TIME_DT is used.'
          TIME_DT_ATMOS_PHY_CH = TIME_DT
       endif
       if ( TIME_DT_ATMOS_PHY_CH_UNIT == '' ) then
          if( IO_L ) write(IO_FID_LOG,*) '*** Not found TIME_DT_ATMOS_PHY_CH_UNIT. TIME_DT_UNIT is used.'
          TIME_DT_ATMOS_PHY_CH_UNIT = TIME_DT_UNIT
       endif
       ! PHY_AE
       if ( TIME_DT_ATMOS_PHY_AE == UNDEF8 ) then
          if( IO_L ) write(IO_FID_LOG,*) '*** Not found TIME_DT_ATMOS_PHY_AE. TIME_DT is used.'
          TIME_DT_ATMOS_PHY_AE = TIME_DT
       endif
       if ( TIME_DT_ATMOS_PHY_AE_UNIT == '' ) then
          if( IO_L ) write(IO_FID_LOG,*) '*** Not found TIME_DT_ATMOS_PHY_AE_UNIT. TIME_DT_UNIT is used.'
          TIME_DT_ATMOS_PHY_AE_UNIT = TIME_DT_UNIT
       endif
       ! ATMOS RESTART
       if ( TIME_DT_ATMOS_RESTART == UNDEF8 ) then
          if( IO_L ) write(IO_FID_LOG,*) '*** Not found TIME_DT_ATMOS_RESTART. TIME_DURATION is used.'
          TIME_DT_ATMOS_RESTART = TIME_DURATION
       endif
       if ( TIME_DT_ATMOS_RESTART_UNIT == '' ) then
          if( IO_L ) write(IO_FID_LOG,*) '*** Not found TIME_DT_ATMOS_RESTART_UNIT. TIME_DURATION_UNIT is used.'
          TIME_DT_ATMOS_RESTART_UNIT = TIME_DURATION_UNIT
       endif
       ! OCEAN
       if ( TIME_DT_OCEAN == UNDEF8 ) then
          if( IO_L ) write(IO_FID_LOG,*) '*** Not found TIME_DT_OCEAN. TIME_DT is used.'
          TIME_DT_OCEAN = TIME_DT
       endif
       if ( TIME_DT_OCEAN_UNIT == '' ) then
          if( IO_L ) write(IO_FID_LOG,*) '*** Not found TIME_DT_OCEAN_UNIT. TIME_DT_UNIT is used.'
          TIME_DT_OCEAN_UNIT = TIME_DT_UNIT
       endif
       ! OCEAN RESTART
       if ( TIME_DT_OCEAN_RESTART == UNDEF8 ) then
          if( IO_L ) write(IO_FID_LOG,*) '*** Not found TIME_DT_OCEAN_RESTART. TIME_DURATION is used.'
          TIME_DT_OCEAN_RESTART = TIME_DURATION
       endif
       if ( TIME_DT_OCEAN_RESTART_UNIT == '' ) then
          if( IO_L ) write(IO_FID_LOG,*) '*** Not found TIME_DT_OCEAN_RESTART_UNIT. TIME_DURATION_UNIT is used.'
          TIME_DT_OCEAN_RESTART_UNIT = TIME_DURATION_UNIT
       endif
       ! LAND
       if ( TIME_DT_LAND == UNDEF8 ) then
          if( IO_L ) write(IO_FID_LOG,*) '*** Not found TIME_DT_LAND. TIME_DT is used.'
          TIME_DT_LAND = TIME_DT
       endif
       if ( TIME_DT_LAND_UNIT == '' ) then
          if( IO_L ) write(IO_FID_LOG,*) '*** Not found TIME_DT_LAND_UNIT. TIME_DT_UNIT is used.'
          TIME_DT_LAND_UNIT = TIME_DT_UNIT
       endif
       ! LAND RESTART
       if ( TIME_DT_LAND_RESTART == UNDEF8 ) then
          if( IO_L ) write(IO_FID_LOG,*) '*** Not found TIME_DT_LAND_RESTART. TIME_DURATION is used.'
          TIME_DT_LAND_RESTART = TIME_DURATION
       endif
       if ( TIME_DT_LAND_RESTART_UNIT == '' ) then
          if( IO_L ) write(IO_FID_LOG,*) '*** Not found TIME_DT_LAND_RESTART_UNIT. TIME_DURATION_UNIT is used.'
          TIME_DT_LAND_RESTART_UNIT = TIME_DURATION_UNIT
       endif
       ! URBAN
       if ( TIME_DT_URBAN == UNDEF8 ) then
          if( IO_L ) write(IO_FID_LOG,*) '*** Not found TIME_DT_URBAN. TIME_DT is used.'
          TIME_DT_URBAN = TIME_DT
       endif
       if ( TIME_DT_URBAN_UNIT == '' ) then
          if( IO_L ) write(IO_FID_LOG,*) '*** Not found TIME_DT_URBAN_UNIT. TIME_DT_UNIT is used.'
          TIME_DT_URBAN_UNIT = TIME_DT_UNIT
       endif
       ! URBAN RESTART
       if ( TIME_DT_URBAN_RESTART == UNDEF8 ) then
          if( IO_L ) write(IO_FID_LOG,*) '*** Not found TIME_DT_URBAN_RESTART. TIME_DURATION is used.'
          TIME_DT_URBAN_RESTART = TIME_DURATION
       endif
       if ( TIME_DT_URBAN_RESTART_UNIT == '' ) then
          if( IO_L ) write(IO_FID_LOG,*) '*** Not found TIME_DT_URBAN_RESTART_UNIT. TIME_DURATION_UNIT is used.'
          TIME_DT_URBAN_RESTART_UNIT = TIME_DURATION_UNIT
       endif
       ! Resume
       if ( TIME_DT_RESUME == UNDEF8 ) then
          TIME_DT_RESUME = TIME_DURATION
       endif
       if ( TIME_DT_RESUME_UNIT == '' ) then
          if( IO_L ) write(IO_FID_LOG,*) '*** Not found TIME_DT_RESUME_UNIT. TIME_DURATION_UNIT is used.'
          TIME_DT_RESUME_UNIT = TIME_DURATION_UNIT
       endif
    endif

    !--- calculate time
    if ( TIME_STARTDATE(1) == -999 ) then
       if ( RESTART_IN_BASENAME /= '' ) then ! read start time from the restart data
          call FileGetDatainfo( RESTART_IN_BASENAME, & ! [IN]
                                'DENS',              & ! [IN]
                                PRC_myrank,          & ! [IN]
                                0,                   & ! [IN] step
                                time_start = cftime, & ! [OUT]
                                time_units = cfunits ) ! [OUT]

          dateday = 0
          datesec = CALENDAR_CFunits2sec( cftime, cfunits, 0 )

          call CALENDAR_adjust_daysec( dateday, datesec )

          call CALENDAR_daysec2date( TIME_STARTDATE, & ! [OUT]
                                     TIME_STARTMS,   & ! [OUT]
                                     dateday,        & ! [IN]
                                     datesec,        & ! [IN]
                                     0               ) ! [IN]
       else
          TIME_STARTDATE = (/ 0, 1, 1, 0, 0, 0 /)
          TIME_STARTMS = 0.0_DP
       endif
    else
       TIME_STARTMS = TIME_STARTMS * 1.E-3_DP
    endif

    TIME_OFFSET_YEAR = TIME_STARTDATE(1)

    call CALENDAR_date2daysec( TIME_STARTDAY,     & ! [OUT]
                               TIME_STARTSEC,     & ! [OUT]
                               TIME_STARTDATE(:), & ! [IN]
                               TIME_STARTMS,      & ! [IN]
                               TIME_OFFSET_YEAR   ) ! [IN]

    call CALENDAR_date2char( startchardate,     & ! [OUT]
                             TIME_STARTDATE(:), & ! [IN]
                             TIME_STARTMS       ) ! [IN]

    TIME_STARTDAYSEC  = CALENDAR_combine_daysec( TIME_STARTDAY, TIME_STARTSEC )

    TIME_NOWDATE(:)   = TIME_STARTDATE(:)
    TIME_NOWMS        = TIME_STARTMS
    TIME_NOWDAY       = TIME_STARTDAY
    TIME_NOWSEC       = TIME_STARTSEC
    TIME_NOWDAYSEC    = CALENDAR_combine_daysec( TIME_NOWDAY, TIME_NOWSEC )

    TIME_ENDDAY       = TIME_STARTDAY

    if ( setup_TimeIntegration ) then
       call CALENDAR_unit2sec( TIME_DURATIONSEC, TIME_DURATION, TIME_DURATION_UNIT )
       TIME_ENDSEC = TIME_STARTSEC + TIME_DURATIONSEC
    else
       TIME_ENDSEC = TIME_STARTSEC
    endif

    call CALENDAR_adjust_daysec( TIME_ENDDAY, TIME_ENDSEC ) ! [INOUT]

    call CALENDAR_daysec2date( TIME_ENDDATE(:), & ! [OUT]
                               TIME_ENDMS,      & ! [OUT]
                               TIME_ENDDAY,     & ! [IN]
                               TIME_ENDSEC,     & ! [IN]
                               TIME_OFFSET_YEAR ) ! [IN]

    call CALENDAR_date2char( endchardate,     & ! [OUT]
                             TIME_ENDDATE(:), & ! [IN]
                             TIME_ENDMS       ) ! [IN]

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** Date/time setting ***'
    if( IO_L ) write(IO_FID_LOG,'(1x,A,A)') '*** START Date     : ', startchardate
    if( IO_L ) write(IO_FID_LOG,'(1x,A,A)') '*** END   Date     : ', endchardate

    if ( setup_TimeIntegration ) then

       call CALENDAR_unit2sec( TIME_DTSEC, TIME_DT, TIME_DT_UNIT )

       TIME_NSTEP   = int( TIME_DURATIONSEC / TIME_DTSEC )
       TIME_NOWSTEP = 1

       if( IO_L ) write(IO_FID_LOG,'(1x,A,F10.3)') '*** delta t (sec.) :', TIME_DTSEC
       if( IO_L ) write(IO_FID_LOG,'(1x,A,I10)'  ) '*** No. of steps   :', TIME_NSTEP

       !--- calculate intervals for atmosphere
       if ( TIME_DT_ATMOS_DYN /= UNDEF8 ) then
          call CALENDAR_unit2sec( TIME_DTSEC_ATMOS_DYN,     TIME_DT_ATMOS_DYN,     TIME_DT_ATMOS_DYN_UNIT     )
       else
          TIME_DTSEC_ATMOS_DYN = TIME_DTSEC / TIME_NSTEP_ATMOS_DYN
       endif
       call CALENDAR_unit2sec( TIME_DTSEC_ATMOS_PHY_CP,  TIME_DT_ATMOS_PHY_CP,  TIME_DT_ATMOS_PHY_CP_UNIT  )
       call CALENDAR_unit2sec( TIME_DTSEC_ATMOS_PHY_MP,  TIME_DT_ATMOS_PHY_MP,  TIME_DT_ATMOS_PHY_MP_UNIT  )
       call CALENDAR_unit2sec( TIME_DTSEC_ATMOS_PHY_RD,  TIME_DT_ATMOS_PHY_RD,  TIME_DT_ATMOS_PHY_RD_UNIT  )
       call CALENDAR_unit2sec( TIME_DTSEC_ATMOS_PHY_SF,  TIME_DT_ATMOS_PHY_SF,  TIME_DT_ATMOS_PHY_SF_UNIT  )
       call CALENDAR_unit2sec( TIME_DTSEC_ATMOS_PHY_TB,  TIME_DT_ATMOS_PHY_TB,  TIME_DT_ATMOS_PHY_TB_UNIT  )
       call CALENDAR_unit2sec( TIME_DTSEC_ATMOS_PHY_CH,  TIME_DT_ATMOS_PHY_CH,  TIME_DT_ATMOS_PHY_CH_UNIT  )
       call CALENDAR_unit2sec( TIME_DTSEC_ATMOS_PHY_AE,  TIME_DT_ATMOS_PHY_AE,  TIME_DT_ATMOS_PHY_AE_UNIT  )
       call CALENDAR_unit2sec( TIME_DTSEC_ATMOS_RESTART, TIME_DT_ATMOS_RESTART, TIME_DT_ATMOS_RESTART_UNIT )
       call CALENDAR_unit2sec( TIME_DTSEC_OCEAN,         TIME_DT_OCEAN,         TIME_DT_OCEAN_UNIT         )
       call CALENDAR_unit2sec( TIME_DTSEC_OCEAN_RESTART, TIME_DT_OCEAN_RESTART, TIME_DT_OCEAN_RESTART_UNIT )
       call CALENDAR_unit2sec( TIME_DTSEC_LAND,          TIME_DT_LAND,          TIME_DT_LAND_UNIT          )
       call CALENDAR_unit2sec( TIME_DTSEC_LAND_RESTART,  TIME_DT_LAND_RESTART,  TIME_DT_LAND_RESTART_UNIT  )
       call CALENDAR_unit2sec( TIME_DTSEC_URBAN,         TIME_DT_URBAN,         TIME_DT_URBAN_UNIT         )
       call CALENDAR_unit2sec( TIME_DTSEC_URBAN_RESTART, TIME_DT_URBAN_RESTART, TIME_DT_URBAN_RESTART_UNIT )
       call CALENDAR_unit2sec( TIME_DTSEC_RESUME,        TIME_DT_RESUME,        TIME_DT_RESUME_UNIT        )

       TIME_NSTEP_ATMOS_DYN = max( nint( TIME_DTSEC / TIME_DTSEC_ATMOS_DYN ), 1 )
       NSTEP_DYN = TIME_NSTEP_ATMOS_DYN

       TIME_DTSEC_ATMOS_DYN     = max( TIME_DTSEC_ATMOS_DYN,     TIME_DTSEC          /TIME_NSTEP_ATMOS_DYN )
       TIME_DTSEC_ATMOS_PHY_CP  = max( TIME_DTSEC_ATMOS_PHY_CP,  TIME_DTSEC_ATMOS_DYN*TIME_NSTEP_ATMOS_DYN )
       TIME_DTSEC_ATMOS_PHY_MP  = max( TIME_DTSEC_ATMOS_PHY_MP,  TIME_DTSEC_ATMOS_DYN*TIME_NSTEP_ATMOS_DYN )
       TIME_DTSEC_ATMOS_PHY_RD  = max( TIME_DTSEC_ATMOS_PHY_RD,  TIME_DTSEC_ATMOS_DYN*TIME_NSTEP_ATMOS_DYN )
       TIME_DTSEC_ATMOS_PHY_SF  = max( TIME_DTSEC_ATMOS_PHY_SF,  TIME_DTSEC_ATMOS_DYN*TIME_NSTEP_ATMOS_DYN )
       TIME_DTSEC_ATMOS_PHY_TB  = max( TIME_DTSEC_ATMOS_PHY_TB,  TIME_DTSEC_ATMOS_DYN*TIME_NSTEP_ATMOS_DYN )
       TIME_DTSEC_ATMOS_PHY_CH  = max( TIME_DTSEC_ATMOS_PHY_CH,  TIME_DTSEC_ATMOS_DYN*TIME_NSTEP_ATMOS_DYN )
       TIME_DTSEC_ATMOS_PHY_AE  = max( TIME_DTSEC_ATMOS_PHY_AE,  TIME_DTSEC_ATMOS_DYN*TIME_NSTEP_ATMOS_DYN )
       TIME_DTSEC_ATMOS_RESTART = max( TIME_DTSEC_ATMOS_RESTART, TIME_DTSEC_ATMOS_DYN*TIME_NSTEP_ATMOS_DYN )
       TIME_DTSEC_OCEAN         = max( TIME_DTSEC_OCEAN,         TIME_DTSEC_ATMOS_DYN*TIME_NSTEP_ATMOS_DYN )
       TIME_DTSEC_OCEAN_RESTART = max( TIME_DTSEC_OCEAN_RESTART, TIME_DTSEC_ATMOS_DYN*TIME_NSTEP_ATMOS_DYN )
       TIME_DTSEC_LAND          = max( TIME_DTSEC_LAND,          TIME_DTSEC_ATMOS_DYN*TIME_NSTEP_ATMOS_DYN )
       TIME_DTSEC_LAND_RESTART  = max( TIME_DTSEC_LAND_RESTART,  TIME_DTSEC_ATMOS_DYN*TIME_NSTEP_ATMOS_DYN )
       TIME_DTSEC_URBAN         = max( TIME_DTSEC_URBAN,         TIME_DTSEC_ATMOS_DYN*TIME_NSTEP_ATMOS_DYN )
       TIME_DTSEC_URBAN_RESTART = max( TIME_DTSEC_URBAN_RESTART, TIME_DTSEC_ATMOS_DYN*TIME_NSTEP_ATMOS_DYN )
       TIME_DTSEC_RESUME        = max( TIME_DTSEC_RESUME,        TIME_DTSEC_ATMOS_DYN*TIME_NSTEP_ATMOS_DYN )

       TIME_DSTEP_ATMOS_DYN     = nint( TIME_DTSEC_ATMOS_DYN     / TIME_DTSEC * TIME_NSTEP_ATMOS_DYN )
       TIME_DSTEP_ATMOS_PHY_CP  = nint( TIME_DTSEC_ATMOS_PHY_CP  / TIME_DTSEC )
       TIME_DSTEP_ATMOS_PHY_MP  = nint( TIME_DTSEC_ATMOS_PHY_MP  / TIME_DTSEC )
       TIME_DSTEP_ATMOS_PHY_RD  = nint( TIME_DTSEC_ATMOS_PHY_RD  / TIME_DTSEC )
       TIME_DSTEP_ATMOS_PHY_SF  = nint( TIME_DTSEC_ATMOS_PHY_SF  / TIME_DTSEC )
       TIME_DSTEP_ATMOS_PHY_TB  = nint( TIME_DTSEC_ATMOS_PHY_TB  / TIME_DTSEC )
       TIME_DSTEP_ATMOS_PHY_CH  = nint( TIME_DTSEC_ATMOS_PHY_CH  / TIME_DTSEC )
       TIME_DSTEP_ATMOS_PHY_AE  = nint( TIME_DTSEC_ATMOS_PHY_AE  / TIME_DTSEC )
       TIME_DSTEP_OCEAN         = nint( TIME_DTSEC_OCEAN         / TIME_DTSEC )
       TIME_DSTEP_LAND          = nint( TIME_DTSEC_LAND          / TIME_DTSEC )
       TIME_DSTEP_URBAN         = nint( TIME_DTSEC_URBAN         / TIME_DTSEC )
       TIME_DSTEP_ATMOS_RESTART = nint( TIME_DTSEC_ATMOS_RESTART / TIME_DTSEC )
       TIME_DSTEP_OCEAN_RESTART = nint( TIME_DTSEC_OCEAN_RESTART / TIME_DTSEC )
       TIME_DSTEP_LAND_RESTART  = nint( TIME_DTSEC_LAND_RESTART  / TIME_DTSEC )
       TIME_DSTEP_URBAN_RESTART = nint( TIME_DTSEC_URBAN_RESTART / TIME_DTSEC )
       TIME_DSTEP_RESUME        = nint( TIME_DTSEC_RESUME        / TIME_DTSEC )

       TIME_RES_RESUME = TIME_DSTEP_RESUME - 1

       if ( abs( real(TIME_NSTEP_ATMOS_DYN,kind=DP)*TIME_DTSEC_ATMOS_DYN &
               - real(TIME_DSTEP_ATMOS_DYN,kind=DP)*TIME_DTSEC            ) > eps ) then
          write(*,*) 'xxx delta t(ATMOS_DYN) must be a multiple of delta t ', &
                     TIME_DTSEC_ATMOS_DYN, real(TIME_DSTEP_ATMOS_DYN,kind=DP)*TIME_DTSEC
          call PRC_MPIstop
       endif
       if ( abs(TIME_DTSEC_ATMOS_PHY_CP-real(TIME_DSTEP_ATMOS_PHY_CP,kind=DP)*TIME_DTSEC) > eps ) then
          write(*,*) 'xxx delta t(ATMOS_PHY_CP) must be a multiple of delta t ', &
                     TIME_DTSEC_ATMOS_PHY_CP, real(TIME_DSTEP_ATMOS_PHY_CP,kind=DP)*TIME_DTSEC
          call PRC_MPIstop
       endif
       if ( abs(TIME_DTSEC_ATMOS_PHY_MP-real(TIME_DSTEP_ATMOS_PHY_MP,kind=DP)*TIME_DTSEC) > eps ) then
          write(*,*) 'xxx delta t(ATMOS_PHY_MP) must be a multiple of delta t ', &
                     TIME_DTSEC_ATMOS_PHY_MP, real(TIME_DSTEP_ATMOS_PHY_MP,kind=DP)*TIME_DTSEC
          call PRC_MPIstop
       endif
       if ( abs(TIME_DTSEC_ATMOS_PHY_RD-real(TIME_DSTEP_ATMOS_PHY_RD,kind=DP)*TIME_DTSEC) > eps ) then
          write(*,*) 'xxx delta t(ATMOS_PHY_RD) must be a multiple of delta t ', &
                     TIME_DTSEC_ATMOS_PHY_RD, real(TIME_DSTEP_ATMOS_PHY_RD,kind=DP)*TIME_DTSEC
          call PRC_MPIstop
       endif
       if ( abs(TIME_DTSEC_ATMOS_PHY_SF-real(TIME_DSTEP_ATMOS_PHY_SF,kind=DP)*TIME_DTSEC) > eps ) then
          write(*,*) 'xxx delta t(ATMOS_PHY_SF) must be a multiple of delta t ', &
                     TIME_DTSEC_ATMOS_PHY_SF, real(TIME_DSTEP_ATMOS_PHY_SF,kind=DP)*TIME_DTSEC
          call PRC_MPIstop
       endif
       if ( abs(TIME_DTSEC_ATMOS_PHY_TB-real(TIME_DSTEP_ATMOS_PHY_TB,kind=DP)*TIME_DTSEC) > eps ) then
          write(*,*) 'xxx delta t(ATMOS_PHY_TB) must be a multiple of delta t ', &
                     TIME_DTSEC_ATMOS_PHY_TB, real(TIME_DSTEP_ATMOS_PHY_TB,kind=DP)*TIME_DTSEC
          call PRC_MPIstop
       endif
       if ( abs(TIME_DTSEC_ATMOS_PHY_CH-real(TIME_DSTEP_ATMOS_PHY_CH,kind=DP)*TIME_DTSEC) > eps ) then
          write(*,*) 'xxx delta t(ATMOS_PHY_CH) must be a multiple of delta t ', &
                     TIME_DTSEC_ATMOS_PHY_CH, real(TIME_DSTEP_ATMOS_PHY_CH,kind=DP)*TIME_DTSEC
          call PRC_MPIstop
       endif
       if ( abs(TIME_DTSEC_ATMOS_PHY_AE-real(TIME_DSTEP_ATMOS_PHY_AE,kind=DP)*TIME_DTSEC) > eps ) then
          write(*,*) 'xxx delta t(ATMOS_PHY_AE) must be a multiple of delta t ', &
                     TIME_DTSEC_ATMOS_PHY_AE, real(TIME_DSTEP_ATMOS_PHY_AE,kind=DP)*TIME_DTSEC
          call PRC_MPIstop
       endif
       if ( abs(TIME_DTSEC_OCEAN-real(TIME_DSTEP_OCEAN,kind=DP)*TIME_DTSEC) > eps ) then
          write(*,*) 'xxx delta t(OCEAN) must be a multiple of delta t ', &
                     TIME_DTSEC_OCEAN, real(TIME_DSTEP_OCEAN,kind=DP)*TIME_DTSEC
          call PRC_MPIstop
       endif
       if ( abs(TIME_DTSEC_LAND-real(TIME_DSTEP_LAND,kind=DP)*TIME_DTSEC) > eps ) then
          write(*,*) 'xxx delta t(LAND) must be a multiple of delta t ', &
                     TIME_DTSEC_LAND, real(TIME_DSTEP_LAND,kind=DP)*TIME_DTSEC
          call PRC_MPIstop
       endif
       if ( abs(TIME_DTSEC_URBAN-real(TIME_DSTEP_URBAN,kind=DP)*TIME_DTSEC) > eps ) then
          write(*,*) 'xxx delta t(URBAN) must be a multiple of delta t ', &
                     TIME_DTSEC_URBAN, real(TIME_DSTEP_URBAN,kind=DP)*TIME_DTSEC
          call PRC_MPIstop
       endif
       if ( abs(TIME_DTSEC_ATMOS_RESTART-real(TIME_DSTEP_ATMOS_RESTART,kind=DP)*TIME_DTSEC) > eps ) then
          write(*,*) 'xxx delta t(ATMOS_RESTART) must be a multiple of delta t ', &
                     TIME_DTSEC_ATMOS_RESTART, real(TIME_DSTEP_ATMOS_RESTART,kind=DP)*TIME_DTSEC
          call PRC_MPIstop
       endif
       if ( abs(TIME_DTSEC_OCEAN_RESTART-real(TIME_DSTEP_OCEAN_RESTART,kind=DP)*TIME_DTSEC) > eps ) then
          write(*,*) 'xxx delta t(OCEAN_RESTART) must be a multiple of delta t ', &
                     TIME_DTSEC_OCEAN_RESTART, real(TIME_DSTEP_OCEAN_RESTART,kind=DP)*TIME_DTSEC
          call PRC_MPIstop
       endif
       if ( abs(TIME_DTSEC_LAND_RESTART-real(TIME_DSTEP_LAND_RESTART,kind=DP)*TIME_DTSEC) > eps ) then
          write(*,*) 'xxx delta t(LAND_RESTART) must be a multiple of delta t ', &
                     TIME_DTSEC_LAND_RESTART, real(TIME_DSTEP_LAND_RESTART,kind=DP)*TIME_DTSEC
          call PRC_MPIstop
       endif
       if ( abs(TIME_DTSEC_URBAN_RESTART-real(TIME_DSTEP_URBAN_RESTART,kind=DP)*TIME_DTSEC) > eps ) then
          write(*,*) 'xxx delta t(URBAN_RESTART) must be a multiple of delta t ', &
                     TIME_DTSEC_URBAN_RESTART, real(TIME_DSTEP_URBAN_RESTART,kind=DP)*TIME_DTSEC
          call PRC_MPIstop
       endif
       if ( abs(TIME_DTSEC_URBAN_RESTART-real(TIME_DSTEP_URBAN_RESTART,kind=DP)*TIME_DTSEC) > eps ) then
          write(*,*) 'xxx delta t(URBAN_RESTART) must be a multiple of delta t ', &
                     TIME_DTSEC_URBAN_RESTART, real(TIME_DSTEP_URBAN_RESTART,kind=DP)*TIME_DTSEC
          call PRC_MPIstop
       endif
       if ( abs(TIME_DTSEC_RESUME-real(TIME_DSTEP_RESUME,kind=DP)*TIME_DTSEC) > eps ) then
          write(*,*) 'xxx delta t(RESUME) must be a multiple of delta t ', &
                     TIME_DTSEC_RESUME, real(TIME_DSTEP_RESUME,kind=DP)*TIME_DTSEC
          call PRC_MPIstop
       endif

       if( IO_L ) write(IO_FID_LOG,*)
       if( IO_L ) write(IO_FID_LOG,*)                     '*** Time interval for each processes (sec.)'
       if( IO_L ) write(IO_FID_LOG,*)                     '*** Atmosphere'
       if( IO_L ) write(IO_FID_LOG,'(1x,A,F10.3)')        '*** Dynamics (time)             : ', TIME_DTSEC_ATMOS_DYN
       if( IO_L ) write(IO_FID_LOG,'(1x,A,I10,A,I8,A)')   '***          (step)             : ', TIME_NSTEP_ATMOS_DYN, &
                                                          ' (step interval=', TIME_DSTEP_ATMOS_DYN,    ')'
       if( IO_L ) write(IO_FID_LOG,'(1x,A,F10.3,A,I8,A)') '*** Physics, Cumulus            : ', TIME_DTSEC_ATMOS_PHY_CP, &
                                                          ' (step interval=', TIME_DSTEP_ATMOS_PHY_CP, ')'
       if( IO_L ) write(IO_FID_LOG,'(1x,A,F10.3,A,I8,A)') '*** Physics, Cloud Microphysics : ', TIME_DTSEC_ATMOS_PHY_MP, &
                                                          ' (step interval=', TIME_DSTEP_ATMOS_PHY_MP, ')'
       if( IO_L ) write(IO_FID_LOG,'(1x,A,F10.3,A,I8,A)') '*** Physics, Radiation          : ', TIME_DTSEC_ATMOS_PHY_RD, &
                                                          ' (step interval=', TIME_DSTEP_ATMOS_PHY_RD, ')'
       if( IO_L ) write(IO_FID_LOG,'(1x,A,F10.3,A,I8,A)') '*** Physics, Surface Flux       : ', TIME_DTSEC_ATMOS_PHY_SF, &
                                                          ' (step interval=', TIME_DSTEP_ATMOS_PHY_SF, ')'
       if( IO_L ) write(IO_FID_LOG,'(1x,A,F10.3,A,I8,A)') '*** Physics, Turbulence         : ', TIME_DTSEC_ATMOS_PHY_TB, &
                                                          ' (step interval=', TIME_DSTEP_ATMOS_PHY_TB, ')'
       if( IO_L ) write(IO_FID_LOG,'(1x,A,F10.3,A,I8,A)') '*** Physics, Chemistry          : ', TIME_DTSEC_ATMOS_PHY_CH, &
                                                          ' (step interval=', TIME_DSTEP_ATMOS_PHY_CH, ')'
       if( IO_L ) write(IO_FID_LOG,'(1x,A,F10.3,A,I8,A)') '*** Physics, Aerosol            : ', TIME_DTSEC_ATMOS_PHY_AE, &
                                                          ' (step interval=', TIME_DSTEP_ATMOS_PHY_AE, ')'
       if( IO_L ) write(IO_FID_LOG,'(1x,A,F10.3,A,I8,A)') '*** Ocean                       : ', TIME_DTSEC_OCEAN, &
                                                          ' (step interval=', TIME_DSTEP_OCEAN,        ')'
       if( IO_L ) write(IO_FID_LOG,'(1x,A,F10.3,A,I8,A)') '*** Land                        : ', TIME_DTSEC_LAND, &
                                                          ' (step interval=', TIME_DSTEP_LAND,         ')'
       if( IO_L ) write(IO_FID_LOG,'(1x,A,F10.3,A,I8,A)') '*** Urban                       : ', TIME_DTSEC_URBAN, &
                                                          ' (step interval=', TIME_DSTEP_URBAN,        ')'
       if( IO_L ) write(IO_FID_LOG,*)
       if( IO_L ) write(IO_FID_LOG,*)                     '*** Time interval for restart (sec.)'
       if( IO_L ) write(IO_FID_LOG,'(1x,A,F10.3,A,I8,A)') '*** Atmospheric Variables       : ', TIME_DTSEC_ATMOS_RESTART, &
                                                          ' (step interval=', TIME_DSTEP_ATMOS_RESTART, ')'
       if( IO_L ) write(IO_FID_LOG,'(1x,A,F10.3,A,I8,A)') '*** Ocean Variables             : ', TIME_DTSEC_OCEAN_RESTART, &
                                                          ' (step interval=', TIME_DSTEP_OCEAN_RESTART, ')'
       if( IO_L ) write(IO_FID_LOG,'(1x,A,F10.3,A,I8,A)') '*** Land Variables              : ', TIME_DTSEC_LAND_RESTART,  &
                                                          ' (step interval=', TIME_DSTEP_LAND_RESTART,  ')'
       if( IO_L ) write(IO_FID_LOG,'(1x,A,F10.3,A,I8,A)') '*** Urban Variables             : ', TIME_DTSEC_URBAN_RESTART, &
                                                          ' (step interval=', TIME_DSTEP_URBAN_RESTART, ')'
       if( IO_L ) write(IO_FID_LOG,'(1x,A,F10.3,A,I8,A)') '*** Resume                      : ', TIME_DTSEC_RESUME, &
                                                          ' (step interval=', TIME_DSTEP_RESUME,        ')'
    else
       TIME_DTSEC = 1.0_RP
    endif

    ! WALLCLOCK TERMINATOR SETUP
    if ( TIME_WALLCLOCK_LIMIT > 0.0_DP ) then
       if( IO_L ) write(IO_FID_LOG,*)
       if( IO_L ) write(IO_FID_LOG,*) '*** Wall clock time limit of execution is specified.'

       if ( TIME_DT_WALLCLOCK_CHECK == UNDEF8 ) then
          if( IO_L ) write(IO_FID_LOG,*) '*** Not found TIME_DT_WALLCLOCK_CHECK. largest time step interval is used.'
          TIME_DTSEC_WALLCLOCK_CHECK = max( TIME_DTSEC_ATMOS_DYN*TIME_NSTEP_ATMOS_DYN, &
                                            TIME_DTSEC_ATMOS_PHY_CP,                   &
                                            TIME_DTSEC_ATMOS_PHY_MP,                   &
                                            TIME_DTSEC_ATMOS_PHY_RD,                   &
                                            TIME_DTSEC_ATMOS_PHY_SF,                   &
                                            TIME_DTSEC_ATMOS_PHY_TB,                   &
                                            TIME_DTSEC_ATMOS_PHY_CH,                   &
                                            TIME_DTSEC_ATMOS_PHY_AE,                   &
                                            TIME_DTSEC_OCEAN,                          &
                                            TIME_DTSEC_LAND,                           &
                                            TIME_DTSEC_URBAN                           )
       else
          if ( TIME_DT_WALLCLOCK_CHECK_UNIT == '' ) then
             if( IO_L ) write(IO_FID_LOG,*) '*** Not found TIME_DT_WALLCLOCK_CHECK_UNIT. TIME_DURATION_UNIT is used.'
             TIME_DT_WALLCLOCK_CHECK_UNIT = TIME_DURATION_UNIT
          endif
          call CALENDAR_unit2sec( TIME_DTSEC_WALLCLOCK_CHECK, TIME_DT_WALLCLOCK_CHECK, TIME_DT_WALLCLOCK_CHECK_UNIT )
          TIME_DTSEC_WALLCLOCK_CHECK = max( TIME_DTSEC_WALLCLOCK_CHECK, TIME_DTSEC_ATMOS_DYN*TIME_NSTEP_ATMOS_DYN )
       endif

       TIME_DSTEP_WALLCLOCK_CHECK = int( TIME_DTSEC_WALLCLOCK_CHECK / TIME_DTSEC )

       TIME_WALLCLOCK_SAFE  = max( min( TIME_WALLCLOCK_SAFE, 1.0_DP ), 0.0_DP )
       TIME_WALLCLOCK_START = PRC_MPItime()

       if( IO_L ) write(IO_FID_LOG,'(1x,A,F10.1,A)')      '*** This job stops after ', &
                                                          TIME_WALLCLOCK_LIMIT * TIME_WALLCLOCK_SAFE, ' seconds.'
       if( IO_L ) write(IO_FID_LOG,'(1x,A,F10.3,A,I8,A)') '*** Step interval for check     : ', TIME_DTSEC_WALLCLOCK_CHECK, &
                                                          ' (step interval=', TIME_DSTEP_WALLCLOCK_CHECK, ')'
    endif

    return
  end subroutine ADMIN_TIME_setup

  !-----------------------------------------------------------------------------
  !> Evaluate component execution
  subroutine ADMIN_TIME_checkstate
    use scale_process, only: &
       PRC_MPItime
    use scale_calendar, only: &
       CALENDAR_date2char
    use scale_time, only: &
       TIME_NOWDATE,            &
       TIME_NOWMS,              &
       TIME_NOWSTEP,            &
       TIME_NSTEP,              &
       TIME_DSTEP_ATMOS_DYN,    &
       TIME_DSTEP_ATMOS_PHY_CP, &
       TIME_DSTEP_ATMOS_PHY_MP, &
       TIME_DSTEP_ATMOS_PHY_RD, &
       TIME_DSTEP_ATMOS_PHY_SF, &
       TIME_DSTEP_ATMOS_PHY_TB, &
       TIME_DSTEP_ATMOS_PHY_CH, &
       TIME_DSTEP_ATMOS_PHY_AE, &
       TIME_DSTEP_OCEAN,        &
       TIME_DSTEP_LAND,         &
       TIME_DSTEP_URBAN
    implicit none

    real(DP)          :: WALLCLOCK_elapse
    character(len=27) :: nowchardate
    !---------------------------------------------------------------------------

    TIME_DOATMOS_step     = .false.
    TIME_DOATMOS_DYN      = .false.
    TIME_DOATMOS_PHY_CP   = .false.
    TIME_DOATMOS_PHY_MP   = .false.
    TIME_DOATMOS_PHY_RD   = .false.
    TIME_DOATMOS_PHY_SF   = .false.
    TIME_DOATMOS_PHY_TB   = .false.
    TIME_DOATMOS_PHY_CH   = .false.
    TIME_DOATMOS_PHY_AE   = .false.
    TIME_DOOCEAN_step     = .false.
    TIME_DOLAND_step      = .false.
    TIME_DOURBAN_step     = .false.
    TIME_DOresume         = .false.

    TIME_RES_ATMOS_DYN    = TIME_RES_ATMOS_DYN    + 1
    TIME_RES_ATMOS_PHY_CP = TIME_RES_ATMOS_PHY_CP + 1
    TIME_RES_ATMOS_PHY_MP = TIME_RES_ATMOS_PHY_MP + 1
    TIME_RES_ATMOS_PHY_RD = TIME_RES_ATMOS_PHY_RD + 1
    TIME_RES_ATMOS_PHY_SF = TIME_RES_ATMOS_PHY_SF + 1
    TIME_RES_ATMOS_PHY_TB = TIME_RES_ATMOS_PHY_TB + 1
    TIME_RES_ATMOS_PHY_CH = TIME_RES_ATMOS_PHY_CH + 1
    TIME_RES_ATMOS_PHY_AE = TIME_RES_ATMOS_PHY_AE + 1
    TIME_RES_OCEAN        = TIME_RES_OCEAN        + 1
    TIME_RES_LAND         = TIME_RES_LAND         + 1
    TIME_RES_URBAN        = TIME_RES_URBAN        + 1
    TIME_RES_RESUME       = TIME_RES_RESUME       + 1

    if ( TIME_RES_ATMOS_DYN    == TIME_DSTEP_ATMOS_DYN ) then
       TIME_DOATMOS_step     = .true.
       TIME_DOATMOS_DYN      = .true.
       TIME_RES_ATMOS_DYN    = 0
    endif
    if ( TIME_RES_ATMOS_PHY_CP == TIME_DSTEP_ATMOS_PHY_CP ) then
       TIME_DOATMOS_step     = .true.
       TIME_DOATMOS_PHY_CP   = .true.
       TIME_RES_ATMOS_PHY_CP = 0
    endif
    if ( TIME_RES_ATMOS_PHY_MP == TIME_DSTEP_ATMOS_PHY_MP ) then
       TIME_DOATMOS_step     = .true.
       TIME_DOATMOS_PHY_MP   = .true.
       TIME_RES_ATMOS_PHY_MP = 0
    endif
    if ( TIME_RES_ATMOS_PHY_RD == TIME_DSTEP_ATMOS_PHY_RD ) then
       TIME_DOATMOS_step     = .true.
       TIME_DOATMOS_PHY_RD   = .true.
       TIME_RES_ATMOS_PHY_RD = 0
    endif
    if ( TIME_RES_ATMOS_PHY_SF == TIME_DSTEP_ATMOS_PHY_SF ) then
       TIME_DOATMOS_step     = .true.
       TIME_DOATMOS_PHY_SF   = .true.
       TIME_RES_ATMOS_PHY_SF = 0
    endif
    if ( TIME_RES_ATMOS_PHY_TB == TIME_DSTEP_ATMOS_PHY_TB ) then
       TIME_DOATMOS_step     = .true.
       TIME_DOATMOS_PHY_TB   = .true.
       TIME_RES_ATMOS_PHY_TB = 0
    endif
    if ( TIME_RES_ATMOS_PHY_CH == TIME_DSTEP_ATMOS_PHY_CH ) then
       TIME_DOATMOS_step     = .true.
       TIME_DOATMOS_PHY_CH   = .true.
       TIME_RES_ATMOS_PHY_CH = 0
    endif
    if ( TIME_RES_ATMOS_PHY_AE == TIME_DSTEP_ATMOS_PHY_AE ) then
       TIME_DOATMOS_step     = .true.
       TIME_DOATMOS_PHY_AE   = .true.
       TIME_RES_ATMOS_PHY_AE = 0
    endif

    if ( TIME_RES_OCEAN  == TIME_DSTEP_OCEAN  ) then
       TIME_DOOCEAN_step     = .true.
       TIME_RES_OCEAN        = 0
    endif
    if ( TIME_RES_LAND   == TIME_DSTEP_LAND   ) then
       TIME_DOLAND_step      = .true.
       TIME_RES_LAND         = 0
    endif
    if ( TIME_RES_URBAN  == TIME_DSTEP_URBAN  ) then
       TIME_DOURBAN_step     = .true.
       TIME_RES_URBAN        = 0
    endif
    if ( TIME_RES_RESUME == TIME_DSTEP_RESUME ) then
       TIME_DOresume         = .true.
       TIME_RES_RESUME       = 0
    endif

    call CALENDAR_date2char( nowchardate,     & ! [OUT]
                             TIME_NOWDATE(:), & ! [IN]
                             TIME_NOWMS       ) ! [IN]

    if ( TIME_WALLCLOCK_LIMIT > 0.0_DP ) then
       WALLCLOCK_elapse = PRC_MPItime() - TIME_WALLCLOCK_START
       if( IO_L ) write(IO_FID_LOG,'(1x,3A,I7,A,I7,A,F10.1)') '*** TIME: ', nowchardate,' STEP:',TIME_NOWSTEP, '/', TIME_NSTEP, &
                                                              ' WCLOCK:', WALLCLOCK_elapse
    else
       if( IO_L ) write(IO_FID_LOG,'(1x,3A,I7,A,I7)') '*** TIME: ', nowchardate,' STEP:',TIME_NOWSTEP, '/', TIME_NSTEP
    endif

    return
  end subroutine ADMIN_TIME_checkstate

  !-----------------------------------------------------------------------------
  !> Advance the time & evaluate restart & stop
  subroutine ADMIN_TIME_advance
    use scale_process, only: &
       PRC_IsMaster, &
       PRC_MPItime
    use scale_calendar, only: &
       CALENDAR_daysec2date,   &
       CALENDAR_adjust_daysec, &
       CALENDAR_combine_daysec
    use scale_time, only: &
       TIME_DTSEC,                &
       TIME_NOWDATE,              &
       TIME_NOWMS,                &
       TIME_NOWDAY,               &
       TIME_NOWSEC,               &
       TIME_NOWDAYSEC,            &
       TIME_NOWSTEP,              &
       TIME_NSTEP,                &
       TIME_OFFSET_YEAR,          &
       TIME_DSTEP_WALLCLOCK_CHECK
    implicit none

    real(DP) :: WALLCLOCK_elapse
    logical  :: exists
    !---------------------------------------------------------------------------

    TIME_DOend = .false.

    TIME_NOWSTEP = TIME_NOWSTEP + 1
    TIME_NOWDAY  = TIME_STARTDAY
    TIME_NOWSEC  = TIME_STARTSEC + real(TIME_NOWSTEP-1,kind=DP) * TIME_DTSEC

    ! reallocate day & sub-day
    call CALENDAR_adjust_daysec( TIME_NOWDAY, TIME_NOWSEC ) ! [INOUT]

    call CALENDAR_daysec2date( TIME_NOWDATE(:), & ! [OUT]
                               TIME_NOWMS,      & ! [OUT]
                               TIME_NOWDAY,     & ! [IN]
                               TIME_NOWSEC,     & ! [IN]
                               TIME_OFFSET_YEAR ) ! [IN]

    TIME_NOWDAYSEC = CALENDAR_combine_daysec( TIME_NOWDAY, TIME_NOWSEC )

    if ( TIME_NOWSTEP > TIME_NSTEP ) then
       TIME_DOend = .true.
    endif

    if (       TIME_WALLCLOCK_LIMIT > 0.0_DP                     &
         .AND. mod(TIME_NOWSTEP-1,TIME_DSTEP_WALLCLOCK_CHECK) == 0 ) then
       WALLCLOCK_elapse = PRC_MPItime() - TIME_WALLCLOCK_START

       if ( WALLCLOCK_elapse > TIME_WALLCLOCK_LIMIT * TIME_WALLCLOCK_SAFE ) then
          if( IO_L ) write(IO_FID_LOG,*) '*** Elapse time limit is detected. Termination operation starts.'
          TIME_DOend = .true.
       endif
    endif

    ! QUIT file control
    if ( PRC_IsMaster ) then ! master node
       inquire(file='QUIT', exist=exists)

       if( exists ) then
          if( IO_L ) write(IO_FID_LOG,*) '*** QUIT file is found. Termination operation starts.'
          TIME_DOend = .true.
       endif
    endif

    TIME_DOATMOS_restart = .false.
    TIME_DOOCEAN_restart = .false.
    TIME_DOLAND_restart  = .false.
    TIME_DOURBAN_restart = .false.

    TIME_RES_ATMOS_RESTART = TIME_RES_ATMOS_RESTART + 1
    TIME_RES_OCEAN_RESTART = TIME_RES_OCEAN_RESTART + 1
    TIME_RES_LAND_RESTART  = TIME_RES_LAND_RESTART  + 1
    TIME_RES_URBAN_RESTART = TIME_RES_URBAN_RESTART + 1

    if ( TIME_RES_ATMOS_RESTART == TIME_DSTEP_ATMOS_RESTART ) then
       TIME_DOATMOS_restart   = .true.
       TIME_RES_ATMOS_RESTART = 0
    elseif( TIME_DOend ) then
       TIME_DOATMOS_restart   = .true.
    endif

    if ( TIME_RES_OCEAN_RESTART == TIME_DSTEP_OCEAN_RESTART ) then
       TIME_DOOCEAN_restart   = .true.
       TIME_RES_OCEAN_RESTART = 0
    elseif( TIME_DOend ) then
       TIME_DOOCEAN_restart   = .true.
    endif

    if ( TIME_RES_LAND_RESTART  == TIME_DSTEP_LAND_RESTART  ) then
       TIME_DOLAND_restart    = .true.
       TIME_RES_LAND_RESTART  = 0
    elseif( TIME_DOend ) then
       TIME_DOLAND_restart    = .true.
    endif

    if ( TIME_RES_URBAN_RESTART == TIME_DSTEP_URBAN_RESTART ) then
       TIME_DOURBAN_restart   = .true.
       TIME_RES_URBAN_RESTART = 0
    elseif( TIME_DOend ) then
       TIME_DOURBAN_restart   = .true.
    endif

    return
  end subroutine ADMIN_TIME_advance

end module mod_admin_time
