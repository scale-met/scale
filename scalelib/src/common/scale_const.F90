!-------------------------------------------------------------------------------
!> module CONSTANT
!!
!! @par Description
!!          Physical constants module
!!
!! @author Team SCALE
!!
!<
#include "scalelib.h"
module scale_const
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_io
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: CONST_setup

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  real(RP), public            :: CONST_PI      = 3.14159265358979_RP !< pi
  real(RP), public            :: CONST_D2R                           !< degree to radian
  real(RP), public            :: CONST_EPS     = 1.E-16_RP           !< small number
  real(RP), public            :: CONST_EPS1    = 0.99999999999999_RP !< small number
  real(RP), public            :: CONST_HUGE    = 1.E+30_RP           !< huge  number
  !$acc declare create(CONST_PI, CONST_EPS, CONST_EPS1, CONST_HUGE)

  integer,  public, parameter :: CONST_UNDEF2  = -32768              !< undefined value (INT2)
  real(SP), public, parameter :: CONST_UNDEF4  = -9.9999E30          !< undefined value (REAL4)
  real(DP), public, parameter :: CONST_UNDEF8  = -9.9999D30          !< undefined value (REAL8)
  real(RP), public            :: CONST_UNDEF

  ! adopted constants
  real(RP), public            :: CONST_RADIUS  = 6.37122E+6_RP       !< radius of the planet [m]
  real(RP), public            :: CONST_OHM     = 7.2920E-5_RP        !< angular velocity of the planet [1/s]
  real(RP), public            :: CONST_GRAV    = 9.80665_RP          !< standard acceleration of gravity [m/s2]

  ! physical constants
  real(RP), public, parameter :: CONST_STB     = 5.670373E-8_RP      !< Stefan-Boltzman constant [W/m2/K4]
  real(RP), public, parameter :: CONST_KARMAN  = 0.4_RP              !< von Karman constant
  real(RP), public, parameter :: CONST_R       = 8.31436_RP          !< universal gas constant [J/mol/K]

  ! dry air constants
  real(RP), public            :: CONST_Mdry    =  28.966_RP          !< mass weight (dry air)                     [g/mol]
  real(RP), public            :: CONST_Rdry    =  287.04_RP          !< specific gas constant (dry air)           [J/kg/K]
  real(RP), public            :: CONST_CPdry   = 1004.64_RP          !< specific heat (dry air,constant pressure) [J/kg/K]
  real(RP), public            :: CONST_CVdry                         !< specific heat (dry air,constant volume)   [J/kg/K]
  real(RP), public            :: CONST_LAPS    = 6.5E-3_RP           !< lapse rate of ISA                         [K/m]
  real(RP), public            :: CONST_LAPSdry                       !< dry adiabatic lapse rate                  [K/m]

  ! water constants
  real(RP), public            :: CONST_Mvap    =  18.016_RP          !< mass weight (water vapor)                      [g/mol]
  real(RP), public, parameter :: CONST_Rvap    =  461.50_RP          !< specific gas constant (water vapor)            [J/kg/K]
  real(RP), public, parameter :: CONST_CPvap   = 1846.00_RP          !< specific heat (water vapor, constant pressure) [J/kg/K]
  real(RP), public            :: CONST_CVvap                         !< specific heat (water vapor, constant volume)   [J/kg/K]
  real(RP), public, parameter :: CONST_CL      = 4218.0_RP           !< specific heat (liquid water)                   [J/kg/K]
  real(RP), public, parameter :: CONST_CI      = 2106.0_RP           !< specific heat (ice)                            [J/kg/K]

  real(RP), public            :: CONST_EPSvap                        !< Rdry / Rvap
  real(RP), public            :: CONST_EPSTvap                       !< 1 / epsilon - 1

  real(RP), public, parameter :: CONST_EMELT   = 3.4E5_RP
  real(RP), public, parameter :: CONST_TMELT   = 273.15_RP

  real(RP), public, parameter :: CONST_LHV0    = 2.501E+6_RP         !< latent heat of vaporizaion at 0C [J/kg]
  real(RP), public            :: CONST_LHV00                         !< latent heat of vaporizaion at 0K [J/kg]
  real(RP), public, parameter :: CONST_LHS0    = 2.834E+6_RP         !< latent heat of sublimation at 0C [J/kg]
  real(RP), public            :: CONST_LHS00                         !< latent heat of sublimation at 0K [J/kg]
  real(RP), public            :: CONST_LHF0                          !< latent heat of fusion      at 0C [J/kg]
  real(RP), public            :: CONST_LHF00                         !< latent heat of fusion      at 0K [J/kg]
  real(RP), public, parameter :: CONST_PSAT0   =  610.78_RP          !< saturate pressure of water vapor at 0C [Pa]
  real(RP), public, parameter :: CONST_DWATR   = 1000.0_RP           !< density of water [kg/m3]
  real(RP), public, parameter :: CONST_DICE    =  916.8_RP           !< density of ice   [kg/m3]

  real(RP), public            :: CONST_SOUND                         !< speed of sound (dry air at 0C) [m/s]

  real(RP), public            :: CONST_Pstd    = 101325.0_RP         !< standard pressure [Pa]
  real(RP), public            :: CONST_PRE00   = 100000.0_RP         !< pressure reference [Pa]
  real(RP), public            :: CONST_Tstd    = 288.15_RP           !< standard temperature (15C) [K]
  real(RP), public, parameter :: CONST_TEM00   = 273.15_RP           !< temperature reference (0C) [K]
  real(RP), public, parameter :: CONST_PPM     = 1.E-6_RP            !< parts par million

  integer,  public            :: CONST_I_LW    = 1                   !< long-wave radiation index
  integer,  public            :: CONST_I_SW    = 2                   !< short-wave radiation index

  character(len=H_SHORT), public :: CONST_THERMODYN_TYPE = 'EXACT'   !< internal energy type

  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  logical, private :: initialized = .false.
  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine CONST_setup
    use scale_prc, only: &
       PRC_abort
    implicit none

    real(RP) :: CONST_SmallPlanetFactor = 1.0_RP !< factor for small planet

    namelist / PARAM_CONST / &
       CONST_RADIUS,           &
       CONST_OHM,              &
       CONST_GRAV,             &
       CONST_Rdry,             &
       CONST_CPdry,            &
       CONST_LAPS,             &
       CONST_Pstd,             &
       CONST_PRE00,            &
       CONST_Tstd,             &
       CONST_THERMODYN_TYPE,   &
       CONST_SmallPlanetFactor

    integer :: ierr
    !---------------------------------------------------------------------------

    if ( initialized ) return
    initialized = .true.

    LOG_NEWLINE
    LOG_INFO("CONST_setup",*) 'Setup'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_CONST,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("CONST_setup",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("CONST_setup",*) 'Not appropriate names in namelist PARAM_CONST. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_CONST)

    if    ( RP == SP ) then
       CONST_UNDEF = real(CONST_UNDEF4,kind=RP)
    elseif( RP == DP ) then
       CONST_UNDEF = real(CONST_UNDEF8,kind=RP)
    else
       LOG_ERROR("CONST_setup",*) 'unsupported precision: ', RP
       call PRC_abort
    endif

    CONST_PI      = 4.0_RP * atan( 1.0_RP )
    CONST_D2R     = CONST_PI / 180.0_RP
    CONST_EPS     =          epsilon(0.0_RP)
    CONST_EPS1    = 1.0_RP - epsilon(0.0_RP)
    CONST_HUGE    =             huge(0.0_RP)

    CONST_RADIUS  = CONST_RADIUS / CONST_SmallPlanetFactor
    CONST_OHM     = CONST_OHM    * CONST_SmallPlanetFactor

    CONST_CVdry   = CONST_CPdry - CONST_Rdry
    CONST_LAPSdry = CONST_GRAV / CONST_CPdry

    CONST_CVvap   = CONST_CPvap - CONST_Rvap
    CONST_EPSvap  = CONST_Rdry / CONST_Rvap
    CONST_EPSTvap = 1.0_RP / CONST_EPSvap - 1.0_RP

    CONST_LHF0    = CONST_LHS0 - CONST_LHV0

    CONST_LHV00   = CONST_LHV0 - ( CONST_CPvap - CONST_CL ) * CONST_TEM00
    CONST_LHS00   = CONST_LHS0 - ( CONST_CPvap - CONST_CI ) * CONST_TEM00
    CONST_LHF00   = CONST_LHF0 - ( CONST_CL    - CONST_CI ) * CONST_TEM00

    CONST_SOUND = sqrt( CONST_CPdry * CONST_Rdry / ( CONST_CPdry - CONST_Rdry ) * CONST_TEM00 )

    LOG_NEWLINE
    LOG_INFO("CONST_setup",*) 'List of constants '
    LOG_INFO_CONT(*) 'PI                                                : PI      = ', CONST_PI
    LOG_INFO_CONT(*) 'Small number                                      : EPS     = ', CONST_EPS
    LOG_INFO_CONT(*) 'Small number (1-EPS)                              : EPS1    = ', CONST_EPS1
    LOG_INFO_CONT(*) 'Huge  number                                      : HUGE    = ', CONST_HUGE
    LOG_INFO_CONT(*) 'undefined number (INT2)                           : UNDEF2  = ', CONST_UNDEF2
    LOG_INFO_CONT(*) 'undefined number (REAL,general use)               : UNDEF   = ', CONST_UNDEF
    LOG_INFO_CONT(*) 'undefined number (REAL4)                          : UNDEF4  = ', CONST_UNDEF4
    LOG_INFO_CONT(*) 'undefined number (REAL8)                          : UNDEF8  = ', CONST_UNDEF8

    LOG_INFO_CONT(*) 'radius of the planet                          [m] : RADIUS  = ', CONST_RADIUS
    LOG_INFO_CONT(*) 'angular velocity of the planet              [1/s] : OHM     = ', CONST_OHM
    LOG_INFO_CONT(*) 'standard acceleration of gravity           [m/s2] : GRAV    = ', CONST_GRAV

    LOG_INFO_CONT(*) 'Stefan-Boltzman constant                [W/m2/K4] : STB     = ', CONST_STB
    LOG_INFO_CONT(*) 'von Karman constant                               : KARMAN  = ', CONST_KARMAN
    LOG_INFO_CONT(*) 'universal gas constant                  [J/mol/K] : R       = ', CONST_R

    LOG_INFO_CONT(*) 'mass weight (dry air)                     [g/mol] : Mdry    = ', CONST_Mdry
    LOG_INFO_CONT(*) 'specific gas constant (dry air)          [J/kg/K] : Rdry    = ', CONST_Rdry
    LOG_INFO_CONT(*) 'specific heat (dry air, const. pressure) [J/kg/K] : CPdry   = ', CONST_CPdry
    LOG_INFO_CONT(*) 'specific heat (dry air, const. volume)   [J/kg/K] : Cvdry   = ', CONST_CVdry
    LOG_INFO_CONT(*) 'lapse rate of ISA                           [K/m] : LAPS    = ', CONST_LAPS
    LOG_INFO_CONT(*) 'dry adiabatic lapse rate                    [K/m] : LAPSdry = ', CONST_LAPSdry

    LOG_INFO_CONT(*) 'mass weight (water vapor)                 [g/mol] : Mvap    = ', CONST_Mvap
    LOG_INFO_CONT(*) 'specific gas constant (water vapor)      [J/kg/K] : Rvap    = ', CONST_Rvap
    LOG_INFO_CONT(*) 'specific heat (vapor, const. pressure)   [J/kg/K] : CPvap   = ', CONST_CPvap
    LOG_INFO_CONT(*) 'specific heat (vapor, const. volume)     [J/kg/K] : CVvap   = ', CONST_CVvap
    LOG_INFO_CONT(*) 'specific heat (liquid water)             [J/kg/K] : CL      = ', CONST_CL
    LOG_INFO_CONT(*) 'specific heat (ice)                      [J/kg/K] : CI      = ', CONST_CI
    LOG_INFO_CONT(*) 'Rdry / Rvap                                       : EPSvap  = ', CONST_EPSvap
    LOG_INFO_CONT(*) '1 / EPSvap - 1                                    : EPSTvap = ', CONST_EPSTvap

    LOG_INFO_CONT(*) 'latent heat of vaporizaion at 0C           [J/kg] : LHV0    = ', CONST_LHV0
    LOG_INFO_CONT(*) 'latent heat of sublimation at 0C           [J/kg] : LHS0    = ', CONST_LHS0
    LOG_INFO_CONT(*) 'latent heat of fusion      at 0C           [J/kg] : LHF0    = ', CONST_LHF0
    LOG_INFO_CONT(*) 'latent heat of vaporizaion at 0K           [J/kg] : LHV00   = ', CONST_LHV00
    LOG_INFO_CONT(*) 'latent heat of sublimation at 0K           [J/kg] : LHS00   = ', CONST_LHS00
    LOG_INFO_CONT(*) 'latent heat of fusion      at 0K           [J/kg] : LHF00   = ', CONST_LHF00
    LOG_INFO_CONT(*) 'Thermodynamics calculation type : ', trim(CONST_THERMODYN_TYPE)
    LOG_INFO_CONT(*) 'saturate pressure of water vapor at 0C       [Pa] : PSAT0   = ', CONST_PSAT0
    LOG_INFO_CONT(*) 'density of water                          [kg/m3] : DWATR   = ', CONST_DWATR
    LOG_INFO_CONT(*) 'density of ice                            [kg/m3] : DICE    = ', CONST_DICE

    LOG_INFO_CONT(*) 'speed of sound (dry air at 0C)              [m/s] : SOUND   = ', CONST_SOUND
    LOG_INFO_CONT(*) 'standard pressure                            [Pa] : Pstd    = ', CONST_Pstd
    LOG_INFO_CONT(*) 'pressure reference                           [Pa] : PRE00   = ', CONST_PRE00
    LOG_INFO_CONT(*) 'standard temperature (15C)                    [K] : Tstd    = ', CONST_Tstd
    LOG_INFO_CONT(*) 'temperature reference (0C)                    [K] : TEM00   = ', CONST_TEM00

    return
  end subroutine CONST_setup

end module scale_const
