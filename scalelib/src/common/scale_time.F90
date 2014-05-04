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
module scale_time
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
  public :: TIME_setup
  public :: TIME_checkstate
  public :: TIME_advance

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  real(DP), public :: TIME_DTSEC                !< time interval of model              [sec]

  real(DP), public :: TIME_DTSEC_ATMOS_DYN      !< time interval of dynamics              [sec]
  integer,  public :: TIME_NSTEP_ATMOS_DYN      !< small step of dynamics
  real(DP), public :: TIME_DTSEC_ATMOS_PHY_CP   !< time interval of physics(cumulus     ) [sec]
  real(DP), public :: TIME_DTSEC_ATMOS_PHY_MP   !< time interval of physics(microphysics) [sec]
  real(DP), public :: TIME_DTSEC_ATMOS_PHY_RD   !< time interval of physics(radiation   ) [sec]
  real(DP), public :: TIME_DTSEC_ATMOS_PHY_SF   !< time interval of physics(surface flux) [sec]
  real(DP), public :: TIME_DTSEC_ATMOS_PHY_TB   !< time interval of physics(turbulence  ) [sec]
  real(DP), public :: TIME_DTSEC_ATMOS_PHY_CH   !< time interval of physics(chemistry   ) [sec]
  real(DP), public :: TIME_DTSEC_ATMOS_PHY_AE   !< time interval of physics(aerosol     ) [sec]
  real(DP), public :: TIME_DTSEC_ATMOS_RESTART  !< time interval of atmosphere restart    [sec]
  real(DP), public :: TIME_DTSEC_OCEAN          !< time interval of ocean step            [sec]
  real(DP), public :: TIME_DTSEC_OCEAN_RESTART  !< time interval of ocean restart         [sec]
  real(DP), public :: TIME_DTSEC_LAND           !< time interval of land step             [sec]
  real(DP), public :: TIME_DTSEC_LAND_RESTART   !< time interval of land restart          [sec]
  real(DP), public :: TIME_DTSEC_CPL            !< time interval of coupler calc.         [sec]
  real(DP), public :: TIME_DTSEC_CPL_RESTART    !< time interval of coupler restart       [sec]

  integer,  public :: TIME_NOWDATE(6)           !< current time [YYYY MM DD HH MM SS]
  real(DP), public :: TIME_NOWMS                !< subsecond part of current time [millisec]
  integer,  public :: TIME_NOWDAY               !< absolute day of current time [day]
  real(DP), public :: TIME_NOWSEC               !< subday part  of current time [sec]
  real(DP), public :: TIME_NOWDAYSEC            !< second of current time [sec]
  integer,  public :: TIME_NOWSTEP              !< current step [number]

  logical,  public :: TIME_DOATMOS_step         !< execute atmospheric component in this step?
  logical,  public :: TIME_DOATMOS_DYN          !< execute dynamics?
  logical,  public :: TIME_DOATMOS_PHY_CP       !< execute physics(cumulus     )?
  logical,  public :: TIME_DOATMOS_PHY_MP       !< execute physics(microphysics)?
  logical,  public :: TIME_DOATMOS_PHY_RD       !< execute physics(radiation   )?
  logical,  public :: TIME_DOATMOS_PHY_SF       !< execute physics(surface flux)?
  logical,  public :: TIME_DOATMOS_PHY_TB       !< execute physics(turbulence  )?
  logical,  public :: TIME_DOATMOS_PHY_CH       !< execute physics(chemistry   )?
  logical,  public :: TIME_DOATMOS_PHY_AE       !< execute physics(aerosol     )?
  logical,  public :: TIME_DOATMOS_restart      !< execute atmosphere restart output?
  logical,  public :: TIME_DOOCEAN_step         !< execute ocean component in this step?
  logical,  public :: TIME_DOOCEAN_restart      !< execute ocean restart output?
  logical,  public :: TIME_DOLAND_step          !< execute land component in this step?
  logical,  public :: TIME_DOLAND_restart       !< execute land restart output?
  logical,  public :: TIME_DOCPL_calc           !< execute coupler component in this step?
  logical,  public :: TIME_DOCPL_restart        !< execute coupler restart output?
  logical,  public :: TIME_DOend                !< finish program in this step?

  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  integer,  private :: TIME_STARTDATE(6) = (/ 0, 1, 1, 0, 0, 0 /)
  real(DP), private :: TIME_STARTMS      = 0.0_DP !< [millisec]
  integer,  private :: TIME_STARTDAY
  real(DP), private :: TIME_STARTSEC

  integer,  private :: TIME_ENDDATE(6)
  real(DP), private :: TIME_ENDMS
  integer,  private :: TIME_ENDDAY
  real(DP), private :: TIME_ENDSEC

  integer,  private :: TIME_NSTEP

  real(DP), private :: TIME_RES_ATMOS_DYN     = 0.0_DP
  real(DP), private :: TIME_RES_ATMOS_PHY_CP  = 0.0_DP
  real(DP), private :: TIME_RES_ATMOS_PHY_MP  = 0.0_DP
  real(DP), private :: TIME_RES_ATMOS_PHY_RD  = 0.0_DP
  real(DP), private :: TIME_RES_ATMOS_PHY_SF  = 0.0_DP
  real(DP), private :: TIME_RES_ATMOS_PHY_TB  = 0.0_DP
  real(DP), private :: TIME_RES_ATMOS_PHY_CH  = 0.0_DP
  real(DP), private :: TIME_RES_ATMOS_PHY_AE  = 0.0_DP
  real(DP), private :: TIME_RES_ATMOS_RESTART = 0.0_DP
  real(DP), private :: TIME_RES_OCEAN         = 0.0_DP
  real(DP), private :: TIME_RES_OCEAN_RESTART = 0.0_DP
  real(DP), private :: TIME_RES_LAND          = 0.0_DP
  real(DP), private :: TIME_RES_LAND_RESTART  = 0.0_DP
  real(DP), private :: TIME_RES_CPL           = 0.0_DP
  real(DP), private :: TIME_RES_CPL_RESTART   = 0.0_DP

  real(DP), private, parameter :: eps = 1.E-10_DP !> epsilon for timesec

  !-----------------------------------------------------------------------------
contains

  !-----------------------------------------------------------------------------
  !> Setup
  subroutine TIME_setup( &
       setup_TimeIntegration )
    use scale_process, only: &
       PRC_MPIstop
    use scale_const, only: &
       UNDEF => CONST_UNDEF
    use scale_calendar, only: &
       CALENDAR_date2daysec,    &
       CALENDAR_daysec2date,    &
       CALENDAR_adjust_daysec,  &
       CALENDAR_combine_daysec, &
       CALENDAR_unit2sec
    implicit none

    logical, intent(in) :: setup_TimeIntegration

    real(DP)               :: TIME_DURATION
    character(len=H_SHORT) :: TIME_DURATION_UNIT         = "SEC"
    real(DP)               :: TIME_DT
    character(len=H_SHORT) :: TIME_DT_UNIT               = "SEC"

    real(DP)               :: TIME_DT_ATMOS_DYN
    character(len=H_SHORT) :: TIME_DT_ATMOS_DYN_UNIT     = "SEC"
    real(DP)               :: TIME_DT_ATMOS_PHY_CP
    character(len=H_SHORT) :: TIME_DT_ATMOS_PHY_CP_UNIT  = ""
    real(DP)               :: TIME_DT_ATMOS_PHY_MP
    character(len=H_SHORT) :: TIME_DT_ATMOS_PHY_MP_UNIT  = ""
    real(DP)               :: TIME_DT_ATMOS_PHY_RD
    character(len=H_SHORT) :: TIME_DT_ATMOS_PHY_RD_UNIT  = ""
    real(DP)               :: TIME_DT_ATMOS_PHY_SF
    character(len=H_SHORT) :: TIME_DT_ATMOS_PHY_SF_UNIT  = ""
    real(DP)               :: TIME_DT_ATMOS_PHY_TB
    character(len=H_SHORT) :: TIME_DT_ATMOS_PHY_TB_UNIT  = ""
    real(DP)               :: TIME_DT_ATMOS_PHY_CH
    character(len=H_SHORT) :: TIME_DT_ATMOS_PHY_CH_UNIT  = ""
    real(DP)               :: TIME_DT_ATMOS_PHY_AE
    character(len=H_SHORT) :: TIME_DT_ATMOS_PHY_AE_UNIT  = ""
    real(DP)               :: TIME_DT_ATMOS_RESTART
    character(len=H_SHORT) :: TIME_DT_ATMOS_RESTART_UNIT = ""

    real(DP)               :: TIME_DT_OCEAN
    character(len=H_SHORT) :: TIME_DT_OCEAN_UNIT         = ""
    real(DP)               :: TIME_DT_OCEAN_RESTART
    character(len=H_SHORT) :: TIME_DT_OCEAN_RESTART_UNIT = ""

    real(DP)               :: TIME_DT_LAND
    character(len=H_SHORT) :: TIME_DT_LAND_UNIT          = ""
    real(DP)               :: TIME_DT_LAND_RESTART
    character(len=H_SHORT) :: TIME_DT_LAND_RESTART_UNIT  = ""

    real(DP)               :: TIME_DT_CPL
    character(len=H_SHORT) :: TIME_DT_CPL_UNIT           = ""
    real(DP)               :: TIME_DT_CPL_RESTART
    character(len=H_SHORT) :: TIME_DT_CPL_RESTART_UNIT   = ""

    NAMELIST / PARAM_TIME / &
       TIME_STARTDATE,             &
       TIME_STARTMS,               &
       TIME_DURATION,              &
       TIME_DURATION_UNIT,         &
       TIME_DT,                    &
       TIME_DT_UNIT,               &
       TIME_DT_ATMOS_DYN,          &
       TIME_DT_ATMOS_DYN_UNIT,     &
       TIME_DT_ATMOS_PHY_CP,       &
       TIME_DT_ATMOS_PHY_CP_UNIT,  &
       TIME_DT_ATMOS_PHY_MP,       &
       TIME_DT_ATMOS_PHY_MP_UNIT,  &
       TIME_DT_ATMOS_PHY_RD,       &
       TIME_DT_ATMOS_PHY_RD_UNIT,  &
       TIME_DT_ATMOS_PHY_SF,       &
       TIME_DT_ATMOS_PHY_SF_UNIT,  &
       TIME_DT_ATMOS_PHY_TB,       &
       TIME_DT_ATMOS_PHY_TB_UNIT,  &
       TIME_DT_ATMOS_PHY_CH,       &
       TIME_DT_ATMOS_PHY_CH_UNIT,  &
       TIME_DT_ATMOS_PHY_AE,       &
       TIME_DT_ATMOS_PHY_AE_UNIT,  &
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

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[TIME]/Categ[COMMON]'

    TIME_DURATION         = UNDEF
    TIME_DT               = UNDEF
    TIME_DT_ATMOS_DYN     = UNDEF
    TIME_DT_ATMOS_PHY_CP  = UNDEF
    TIME_DT_ATMOS_PHY_MP  = UNDEF
    TIME_DT_ATMOS_PHY_RD  = UNDEF
    TIME_DT_ATMOS_PHY_SF  = UNDEF
    TIME_DT_ATMOS_PHY_TB  = UNDEF
    TIME_DT_ATMOS_PHY_CH  = UNDEF
    TIME_DT_ATMOS_PHY_AE  = UNDEF
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

    ! check time setting
    if ( setup_TimeIntegration ) then
       if ( TIME_DT == UNDEF ) then
          if( IO_L ) write(IO_FID_LOG,*) 'xxx Not found TIME_DT.'
          call PRC_MPIstop
       endif
       if ( TIME_DURATION == UNDEF ) then
          if( IO_L ) write(IO_FID_LOG,*) 'xxx Not found TIME_DURATION.'
          call PRC_MPIstop
       endif

       ! DYN
       if ( TIME_DT_ATMOS_DYN == UNDEF ) then
          if( IO_L ) write(IO_FID_LOG,*) '*** Not found TIME_DT_ATMOS_DYN. TIME_DT is used.'
          TIME_DT_ATMOS_DYN = TIME_DT
       endif
       if ( TIME_DT_ATMOS_DYN_UNIT == '' ) then
          if( IO_L ) write(IO_FID_LOG,*) '*** Not found TIME_DT_ATMOS_DYN_UNIT. TIME_DT_UNIT is used.'
          TIME_DT_ATMOS_DYN_UNIT = TIME_DT_UNIT
       endif
       ! PHY_CP
       if ( TIME_DT_ATMOS_PHY_CP == UNDEF ) then
          if( IO_L ) write(IO_FID_LOG,*) '*** Not found TIME_DT_ATMOS_PHY_CP. TIME_DT is used.'
          TIME_DT_ATMOS_PHY_CP = TIME_DT
       endif
       if ( TIME_DT_ATMOS_PHY_CP_UNIT == '' ) then
          if( IO_L ) write(IO_FID_LOG,*) '*** Not found TIME_DT_ATMOS_PHY_CP_UNIT. TIME_DT_UNIT is used.'
          TIME_DT_ATMOS_PHY_CP_UNIT = TIME_DT_UNIT
       endif
       ! PHY_MP
       if ( TIME_DT_ATMOS_PHY_MP == UNDEF ) then
          if( IO_L ) write(IO_FID_LOG,*) '*** Not found TIME_DT_ATMOS_PHY_MP. TIME_DT is used.'
          TIME_DT_ATMOS_PHY_MP = TIME_DT
       endif
       if ( TIME_DT_ATMOS_PHY_MP_UNIT == '' ) then
          if( IO_L ) write(IO_FID_LOG,*) '*** Not found TIME_DT_ATMOS_PHY_MP_UNIT. TIME_DT_UNIT is used.'
          TIME_DT_ATMOS_PHY_MP_UNIT = TIME_DT_UNIT
       endif
       ! PHY_RD
       if ( TIME_DT_ATMOS_PHY_RD == UNDEF ) then
          if( IO_L ) write(IO_FID_LOG,*) '*** Not found TIME_DT_ATMOS_PHY_RD. TIME_DT is used.'
          TIME_DT_ATMOS_PHY_RD = TIME_DT
       endif
       if ( TIME_DT_ATMOS_PHY_RD_UNIT == '' ) then
          if( IO_L ) write(IO_FID_LOG,*) '*** Not found TIME_DT_ATMOS_PHY_RD_UNIT. TIME_DT_UNIT is used.'
          TIME_DT_ATMOS_PHY_RD_UNIT = TIME_DT_UNIT
       endif
       ! PHY_SF
       if ( TIME_DT_ATMOS_PHY_SF == UNDEF ) then
          if( IO_L ) write(IO_FID_LOG,*) '*** Not found TIME_DT_ATMOS_PHY_SF. TIME_DT is used.'
          TIME_DT_ATMOS_PHY_SF = TIME_DT
       endif
       if ( TIME_DT_ATMOS_PHY_SF_UNIT == '' ) then
          if( IO_L ) write(IO_FID_LOG,*) '*** Not found TIME_DT_ATMOS_PHY_SF_UNIT. TIME_DT_UNIT is used.'
          TIME_DT_ATMOS_PHY_SF_UNIT = TIME_DT_UNIT
       endif
       ! PHY_TB
       if ( TIME_DT_ATMOS_PHY_TB == UNDEF ) then
          if( IO_L ) write(IO_FID_LOG,*) '*** Not found TIME_DT_ATMOS_PHY_TB. TIME_DT is used.'
          TIME_DT_ATMOS_PHY_TB = TIME_DT
       endif
       if ( TIME_DT_ATMOS_PHY_TB_UNIT == '' ) then
          if( IO_L ) write(IO_FID_LOG,*) '*** Not found TIME_DT_ATMOS_PHY_TB_UNIT. TIME_DT_UNIT is used.'
          TIME_DT_ATMOS_PHY_TB_UNIT = TIME_DT_UNIT
       endif
       ! PHY_CH
       if ( TIME_DT_ATMOS_PHY_CH == UNDEF ) then
          if( IO_L ) write(IO_FID_LOG,*) '*** Not found TIME_DT_ATMOS_PHY_CH. TIME_DT is used.'
          TIME_DT_ATMOS_PHY_CH = TIME_DT
       endif
       if ( TIME_DT_ATMOS_PHY_CH_UNIT == '' ) then
          if( IO_L ) write(IO_FID_LOG,*) '*** Not found TIME_DT_ATMOS_PHY_CH_UNIT. TIME_DT_UNIT is used.'
          TIME_DT_ATMOS_PHY_CH_UNIT = TIME_DT_UNIT
       endif
       ! PHY_AE
       if ( TIME_DT_ATMOS_PHY_AE == UNDEF ) then
          if( IO_L ) write(IO_FID_LOG,*) '*** Not found TIME_DT_ATMOS_PHY_AE. TIME_DT is used.'
          TIME_DT_ATMOS_PHY_AE = TIME_DT
       endif
       if ( TIME_DT_ATMOS_PHY_AE_UNIT == '' ) then
          if( IO_L ) write(IO_FID_LOG,*) '*** Not found TIME_DT_ATMOS_PHY_AE_UNIT. TIME_DT_UNIT is used.'
          TIME_DT_ATMOS_PHY_AE_UNIT = TIME_DT_UNIT
       endif
       ! ATMOS RESTART
       if ( TIME_DT_ATMOS_RESTART == UNDEF ) then
          if( IO_L ) write(IO_FID_LOG,*) '*** Not found TIME_DT_ATMOS_RESTART. TIME_DURATION is used.'
          TIME_DT_ATMOS_RESTART = TIME_DURATION
       endif
       if ( TIME_DT_ATMOS_RESTART_UNIT == '' ) then
          if( IO_L ) write(IO_FID_LOG,*) '*** Not found TIME_DT_ATMOS_RESTART_UNIT. TIME_DURATION_UNIT is used.'
          TIME_DT_ATMOS_RESTART_UNIT = TIME_DURATION_UNIT
       endif
       ! OCEAN
       if ( TIME_DT_OCEAN == UNDEF ) then
          if( IO_L ) write(IO_FID_LOG,*) '*** Not found TIME_DT_OCEAN. TIME_DT is used.'
          TIME_DT_OCEAN = TIME_DT
       endif
       if ( TIME_DT_OCEAN_UNIT == '' ) then
          if( IO_L ) write(IO_FID_LOG,*) '*** Not found TIME_DT_OCEAN_UNIT. TIME_DT_UNIT is used.'
          TIME_DT_OCEAN_UNIT = TIME_DT_UNIT
       endif
       ! OCEAN RESTART
       if ( TIME_DT_OCEAN_RESTART == UNDEF ) then
          if( IO_L ) write(IO_FID_LOG,*) '*** Not found TIME_DT_OCEAN_RESTART. TIME_DURATION is used.'
          TIME_DT_OCEAN_RESTART = TIME_DURATION
       endif
       if ( TIME_DT_OCEAN_RESTART_UNIT == '' ) then
          if( IO_L ) write(IO_FID_LOG,*) '*** Not found TIME_DT_OCEAN_RESTART_UNIT. TIME_DURATION_UNIT is used.'
          TIME_DT_OCEAN_RESTART_UNIT = TIME_DURATION_UNIT
       endif
       ! LAND
       if ( TIME_DT_LAND == UNDEF ) then
          if( IO_L ) write(IO_FID_LOG,*) '*** Not found TIME_DT_LAND. TIME_DT is used.'
          TIME_DT_LAND = TIME_DT
       endif
       if ( TIME_DT_LAND_UNIT == '' ) then
          if( IO_L ) write(IO_FID_LOG,*) '*** Not found TIME_DT_LAND_UNIT. TIME_DT_UNIT is used.'
          TIME_DT_LAND_UNIT = TIME_DT_UNIT
       endif
       ! LAND RESTART
       if ( TIME_DT_LAND_RESTART == UNDEF ) then
          if( IO_L ) write(IO_FID_LOG,*) '*** Not found TIME_DT_LAND_RESTART. TIME_DURATION is used.'
          TIME_DT_LAND_RESTART = TIME_DURATION
       endif
       if ( TIME_DT_LAND_RESTART_UNIT == '' ) then
          if( IO_L ) write(IO_FID_LOG,*) '*** Not found TIME_DT_LAND_RESTART_UNIT. TIME_DURATION_UNIT is used.'
          TIME_DT_LAND_RESTART_UNIT = TIME_DURATION_UNIT
       endif
       ! COUPLER
       if ( TIME_DT_CPL == UNDEF ) then
          if( IO_L ) write(IO_FID_LOG,*) '*** Not found TIME_DT_CPL. TIME_DT is used.'
          TIME_DT_CPL = TIME_DT
       endif
       if ( TIME_DT_CPL_UNIT == '' ) then
          if( IO_L ) write(IO_FID_LOG,*) '*** Not found TIME_DT_CPL_UNIT. TIME_DT_UNIT is used.'
          TIME_DT_CPL_UNIT = TIME_DT_UNIT
       endif
       ! CPL RESTART
       if ( TIME_DT_CPL_RESTART == UNDEF ) then
          if( IO_L ) write(IO_FID_LOG,*) '*** Not found TIME_DT_CPL_RESTART. TIME_DURATION is used.'
          TIME_DT_CPL_RESTART = TIME_DURATION
       endif
       if ( TIME_DT_CPL_RESTART_UNIT == '' ) then
          if( IO_L ) write(IO_FID_LOG,*) '*** Not found TIME_DT_CPL_RESTART_UNIT. TIME_DURATION_UNIT is used.'
          TIME_DT_CPL_RESTART_UNIT = TIME_DURATION_UNIT
       endif
    endif

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
                               TIME_ENDSEC      ) ! [IN]

    if ( setup_TimeIntegration ) then

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

       !--- calculate intervals for atmosphere
       call CALENDAR_unit2sec( TIME_DTSEC_ATMOS_DYN,     TIME_DT_ATMOS_DYN,     TIME_DT_ATMOS_DYN_UNIT     )
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
       call CALENDAR_unit2sec( TIME_DTSEC_CPL,           TIME_DT_CPL,           TIME_DT_CPL_UNIT           )
       call CALENDAR_unit2sec( TIME_DTSEC_CPL_RESTART,   TIME_DT_CPL_RESTART,   TIME_DT_CPL_RESTART_UNIT   )

       TIME_NSTEP_ATMOS_DYN = max( int( TIME_DTSEC / TIME_DTSEC_ATMOS_DYN ), 1 )

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
       TIME_DTSEC_CPL           = max( TIME_DTSEC_CPL,           TIME_DTSEC_ATMOS_DYN*TIME_NSTEP_ATMOS_DYN )
       TIME_DTSEC_CPL_RESTART   = max( TIME_DTSEC_CPL_RESTART,   TIME_DTSEC_ATMOS_DYN*TIME_NSTEP_ATMOS_DYN )

       if( IO_L ) write(IO_FID_LOG,*)
       if( IO_L ) write(IO_FID_LOG,*) '*** Time interval for atmospheric processes (sec.)'
       if( IO_L ) write(IO_FID_LOG,'(1x,A,F10.3)') '*** Dynamics (time)             : ', TIME_DTSEC_ATMOS_DYN
       if( IO_L ) write(IO_FID_LOG,'(1x,A,I4)')    '***          (step)             : ', TIME_NSTEP_ATMOS_DYN
       if( IO_L ) write(IO_FID_LOG,'(1x,A,F10.3)') '*** Physics, Cumulus            : ', TIME_DTSEC_ATMOS_PHY_CP
       if( IO_L ) write(IO_FID_LOG,'(1x,A,F10.3)') '*** Physics, Cloud Microphysics : ', TIME_DTSEC_ATMOS_PHY_MP
       if( IO_L ) write(IO_FID_LOG,'(1x,A,F10.3)') '*** Physics, Radiation          : ', TIME_DTSEC_ATMOS_PHY_RD
       if( IO_L ) write(IO_FID_LOG,'(1x,A,F10.3)') '*** Physics, Surface Flux       : ', TIME_DTSEC_ATMOS_PHY_SF
       if( IO_L ) write(IO_FID_LOG,'(1x,A,F10.3)') '*** Physics, Turbulence         : ', TIME_DTSEC_ATMOS_PHY_TB
       if( IO_L ) write(IO_FID_LOG,'(1x,A,F10.3)') '*** Physics, Chemistry          : ', TIME_DTSEC_ATMOS_PHY_CH
       if( IO_L ) write(IO_FID_LOG,'(1x,A,F10.3)') '*** Physics, Aerosol            : ', TIME_DTSEC_ATMOS_PHY_AE
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

    endif

    ! only for register
    call PROF_rapstart('Debug')
    call PROF_rapend  ('Debug')

    return
  end subroutine TIME_setup

  !-----------------------------------------------------------------------------
  !> Evaluate component execution
  subroutine TIME_checkstate
    implicit none
    !---------------------------------------------------------------------------

    TIME_DOATMOS_step      = .false.
    TIME_DOATMOS_DYN       = .false.
    TIME_DOATMOS_PHY_CP    = .false.
    TIME_DOATMOS_PHY_MP    = .false.
    TIME_DOATMOS_PHY_RD    = .false.
    TIME_DOATMOS_PHY_SF    = .false.
    TIME_DOATMOS_PHY_TB    = .false.
    TIME_DOATMOS_PHY_CH    = .false.
    TIME_DOATMOS_PHY_AE    = .false.
    TIME_DOOCEAN_step      = .false.
    TIME_DOLAND_step       = .false.
    TIME_DOCPL_calc        = .false.

    TIME_RES_ATMOS_DYN    = TIME_RES_ATMOS_DYN    + TIME_DTSEC
    TIME_RES_ATMOS_PHY_CP = TIME_RES_ATMOS_PHY_CP + TIME_DTSEC
    TIME_RES_ATMOS_PHY_MP = TIME_RES_ATMOS_PHY_MP + TIME_DTSEC
    TIME_RES_ATMOS_PHY_RD = TIME_RES_ATMOS_PHY_RD + TIME_DTSEC
    TIME_RES_ATMOS_PHY_SF = TIME_RES_ATMOS_PHY_SF + TIME_DTSEC
    TIME_RES_ATMOS_PHY_TB = TIME_RES_ATMOS_PHY_TB + TIME_DTSEC
    TIME_RES_ATMOS_PHY_CH = TIME_RES_ATMOS_PHY_CH + TIME_DTSEC
    TIME_RES_ATMOS_PHY_AE = TIME_RES_ATMOS_PHY_AE + TIME_DTSEC
    TIME_RES_OCEAN        = TIME_RES_OCEAN        + TIME_DTSEC
    TIME_RES_LAND         = TIME_RES_LAND         + TIME_DTSEC
    TIME_RES_CPL          = TIME_RES_CPL          + TIME_DTSEC

    if ( TIME_RES_ATMOS_DYN - TIME_DTSEC_ATMOS_DYN > -eps ) then
       TIME_DOATMOS_step  = .true.
       TIME_DOATMOS_DYN   = .true.
       TIME_RES_ATMOS_DYN = TIME_RES_ATMOS_DYN - TIME_DTSEC_ATMOS_DYN
    endif
    if ( TIME_RES_ATMOS_PHY_CP - TIME_DTSEC_ATMOS_PHY_CP > -eps ) then
       TIME_DOATMOS_step     = .true.
       TIME_DOATMOS_PHY_CP   = .true.
       TIME_RES_ATMOS_PHY_CP = TIME_RES_ATMOS_PHY_CP - TIME_DTSEC_ATMOS_PHY_CP
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
    if ( TIME_RES_ATMOS_PHY_CH - TIME_DTSEC_ATMOS_PHY_CH > -eps ) then
       TIME_DOATMOS_step     = .true.
       TIME_DOATMOS_PHY_CH   = .true.
       TIME_RES_ATMOS_PHY_CH = TIME_RES_ATMOS_PHY_CH - TIME_DTSEC_ATMOS_PHY_CH
    endif
    if ( TIME_RES_ATMOS_PHY_AE - TIME_DTSEC_ATMOS_PHY_AE > -eps ) then
       TIME_DOATMOS_step     = .true.
       TIME_DOATMOS_PHY_AE   = .true.
       TIME_RES_ATMOS_PHY_AE = TIME_RES_ATMOS_PHY_AE - TIME_DTSEC_ATMOS_PHY_AE
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
    use scale_process, only: &
       PRC_master, &
       PRC_myrank
    use scale_calendar, only: &
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
  !> generate time label
  subroutine TIME_gettimelabel( &
       timelabel )
    implicit none

    character(len=15), intent(out) :: timelabel

    integer :: n
    !---------------------------------------------------------------------------

    timelabel = ''
    write(timelabel(1:15), '(F15.3)') TIME_NOWDAYSEC
    do n = 1, 15
       if ( timelabel(n:n) == ' ' ) timelabel(n:n) = '0'
    enddo

    return
  end subroutine TIME_gettimelabel

end module scale_time
!-------------------------------------------------------------------------------
