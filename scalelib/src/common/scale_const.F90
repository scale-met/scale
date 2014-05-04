!-------------------------------------------------------------------------------
!> module CONSTANT
!!
!! @par Description
!!          Physical constants module
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2011-11-11 (H.Yashiro)  [new]
!! @li      2013-08-31 (T.Yamaura)  [add] Stefan-Boltzman constant
!!
!<
module scale_const
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_stdio
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

  integer,  public, parameter :: CONST_UNDEF2  = -32768              !< undefined value (INT2)
  real(SP), public, parameter :: CONST_UNDEF4  = -9.9999E30          !< undefined value (REAL4)
  real(DP), public, parameter :: CONST_UNDEF8  = -9.9999D30          !< undefined value (REAL8)
  real(RP), public            :: CONST_UNDEF

  ! physical constants
  real(RP), public            :: CONST_RADIUS  = 6.37122E+6_RP       !< radius of the planet [m]
  real(RP), public            :: CONST_OHM     = 7.2920E-5_RP        !< Angular velocity of the planet [1/s]
  real(RP), public            :: CONST_GRAV    = 9.80665_RP          !< gravitational constant [m/s2]
  real(RP), public, parameter :: CONST_KARMAN  = 0.4_RP              !< karman constant
  real(RP), public, parameter :: CONST_STB     = 5.67E-8_RP          !< Stefan-Boltzman constant

  real(RP), public            :: CONST_R       = 8.3144621_RP        !< universal gas constant [J/mol/K]

  real(RP), public            :: CONST_Mdry    =   28.97_RP          !< mass weight (dry)                      [g/mol]
  real(RP), public            :: CONST_Rdry    =  287.04_RP          !< specific gas constant   (dry)          [J/kg/K]
  real(RP), public            :: CONST_CPdry   = 1004.64_RP          !< specific heat  (dry,constant pressure) [J/kg/K]
  real(RP), public            :: CONST_CVdry                         !< specific heat  (dry,constant volume)   [J/kg/K]
  real(RP), public            :: CONST_LASPdry                       !< g / Cp, Dry adiabatic lapse rate

  ! [todo] out of date
  real(RP), public            :: CONST_RovCP                         !< R / Cp = kappa (dry)
  real(RP), public            :: CONST_CPovR                         !< 1 / kappa      (dry)
  real(RP), public            :: CONST_RovCV                         !< R / Cv         (dry)
  real(RP), public            :: CONST_CPovCV                        !< Cp / Cv        (dry)
  real(RP), public            :: CONST_CVovCP                        !< Cv / Cp        (dry)

  ! water constants
  real(RP), public            :: CONST_Mvap    =  18.02_RP           !< mass weight  (water vapor) [g/mol]
  real(RP), public, parameter :: CONST_Rvap    = 461.46_RP           !< gas constant (water vapor) [J/kg/K]
  real(RP), public            :: CONST_EPSvap                        !< Rdry / Rvap
  real(RP), public            :: CONST_EPSTvap                       !< 1 / epsilon - 1

  real(RP), public, parameter :: CONST_CPvap   = 1845.60_RP          !< specific heat (water vapor, consant pressure) [J/kg/K]
  real(RP), public            :: CONST_CVvap                         !< specific heat (water vapor, consant volume)   [J/kg/K]
  real(RP), public, parameter :: CONST_CL      = 4218.0_RP           !< specific heat (liquid water)
  real(RP), public, parameter :: CONST_CI      = 2006.0_RP           !< specific heat (ice)

  real(RP), public, parameter :: CONST_EMELT   = 3.4E5_RP
  real(RP), public, parameter :: CONST_TMELT   = 273.15_RP

  real(RP), public, parameter :: CONST_LH0     = 2.5008E6_RP         !< latent heat of vaporizaion at 0 degree
  real(RP), public            :: CONST_LH00
  real(RP), public, parameter :: CONST_LHS0    = 2.8342E6_RP
  real(RP), public            :: CONST_LHS00
  real(RP), public            :: CONST_LHF0
  real(RP), public            :: CONST_LHF00
  real(RP), public, parameter :: CONST_PSAT0   =  610.7_RP           !< saturate pressure of water vapor at 0C
  real(RP), public, parameter :: CONST_DWATR   = 1000.0_RP           !< density of water [kg/m3]
  real(RP), public, parameter :: CONST_DICE    =  916.8_RP           !< density of ice   [kg/m3]

  real(RP), public            :: CONST_Pstd    = 101325.0_RP         !< standard pressure [Pa]
  real(RP), public            :: CONST_PRE00   = 100000.0_RP         !< pressure reference [Pa]
  real(RP), public            :: CONST_Tstd    = 288.15_RP           !< standard temperature (15 degree C) [K]
  real(RP), public, parameter :: CONST_TEM00   = 273.15_RP           !< temperature reference (0 degree C) [K]
  real(RP), public, parameter :: CONST_PPM     = 1.E-6_RP            !< parts par million

  ! tentative
  integer,  public, save      :: CONST_NRAD    = 2                   !< # of radiation categories
  integer,  public, save      :: CONST_I_SW    = 1                   !< short-wave radiation index
  integer,  public, save      :: CONST_I_LW    = 2                   !< long-wave radiation index

  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine CONST_setup
    use scale_process, only: &
       PRC_MPIstop
    implicit none

    namelist / PARAM_CONST / &
       CONST_RADIUS, &
       CONST_OHM, &
       CONST_GRAV, &
       CONST_Rdry, &
       CONST_CPdry, &
       CONST_Pstd, &
       CONST_PRE00, &
       CONST_Tstd

    integer :: ierr
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[CONST]/Categ[COMMON]'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_CONST,iostat=ierr)

    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_CONST. Check!'
       call PRC_MPIstop
    endif
    if( IO_L ) write(IO_FID_LOG,nml=PARAM_CONST)

    if    ( RP == SP ) then
       CONST_UNDEF = real(CONST_UNDEF4,kind=RP)
    elseif( RP == DP ) then
       CONST_UNDEF = real(CONST_UNDEF8,kind=RP)
    else
       write(*,*) 'xxx unsupported precision: ', RP
       call PRC_MPIstop
    endif

    CONST_PI      = 4.E0_RP * atan( 1.0_RP )
    CONST_D2R     = CONST_PI / 180.0_RP
    CONST_EPS     =          epsilon(0.0_RP)
    CONST_EPS1    = 1.0_RP - epsilon(0.0_RP)

    CONST_CVdry   = CONST_CPdry - CONST_Rdry
    CONST_LASPdry = CONST_GRAV / CONST_CPdry

    CONST_CVvap   = CONST_CPvap - CONST_Rvap
    CONST_EPSvap  = CONST_Rdry / CONST_Rvap
    CONST_EPSTvap = 1.E0_RP / CONST_EPSvap - 1.E0_RP

    ! [todo] out of date
    CONST_RovCP   = CONST_Rdry  / CONST_CPdry
    CONST_RovCV   = CONST_Rdry  / CONST_CVdry
    CONST_CPovR   = CONST_CPdry / CONST_Rdry
    CONST_CPovCV  = CONST_CPdry / CONST_CVdry
    CONST_CVovCP  = CONST_CVdry / CONST_CPdry

    CONST_LH00    = CONST_LH0  - ( CONST_CPvap - CONST_CL ) * CONST_TEM00
    CONST_LHS00   = CONST_LHS0 - ( CONST_CPvap - CONST_CI ) * CONST_TEM00
    CONST_LHF0    = CONST_LHS0 - CONST_LH0
    CONST_LHF00   = CONST_LHF0 - ( CONST_CL  - CONST_CI ) * CONST_TEM00

    if( IO_L ) write(IO_FID_LOG,*) '*** PI                                            : PI      = ', CONST_PI
    if( IO_L ) write(IO_FID_LOG,*) '*** Small num                                     : EPS     = ', CONST_EPS
    if( IO_L ) write(IO_FID_LOG,*) '*** Small num (1-EPS)                             : EPS1    = ', CONST_EPS1
    if( IO_L ) write(IO_FID_LOG,*) '*** undefined num(INT2)                           : UNDEF2  = ', CONST_UNDEF2
    if( IO_L ) write(IO_FID_LOG,*) '*** undefined num(REAL,general use)               : UNDEF   = ', CONST_UNDEF
    if( IO_L ) write(IO_FID_LOG,*) '*** undefined num(REAL4)                          : UNDEF4  = ', CONST_UNDEF4
    if( IO_L ) write(IO_FID_LOG,*) '*** undefined num(REAL8)                          : UNDEF8  = ', CONST_UNDEF8

    if( IO_L ) write(IO_FID_LOG,*) '*** radius of the planet [m]                      : RADIUS  = ', CONST_RADIUS
    if( IO_L ) write(IO_FID_LOG,*) '*** angular velocity of the planet [1/s]          : OHM     = ', CONST_OHM
    if( IO_L ) write(IO_FID_LOG,*) '*** gravitational constant [m/s2]                 : GRAV    = ', CONST_GRAV
    if( IO_L ) write(IO_FID_LOG,*) '*** karman constant                               : KARMAN  = ', CONST_KARMAN

    if( IO_L ) write(IO_FID_LOG,*) '*** universal gas constant [J/mol/K]              : R       = ', CONST_R
    if( IO_L ) write(IO_FID_LOG,*) '*** gas constant   (dry)          [J/kg/K]        : Rdry    = ', CONST_Rdry
    if( IO_L ) write(IO_FID_LOG,*) '*** specific heat  (dry,pressure) [J/kg/K]        : CPdry   = ', CONST_CPdry
    if( IO_L ) write(IO_FID_LOG,*) '*** specific heat  (dry,volume)   [J/kg/K]        : Cvdry   = ', CONST_CVdry
    if( IO_L ) write(IO_FID_LOG,*) '*** dry adiabatic lapse rate                      : LASPdry = ', CONST_LASPdry
    if( IO_L ) write(IO_FID_LOG,*) '*** gas constant (water vapor) [J/kg/K]           : Rvap    = ', CONST_Rvap
    if( IO_L ) write(IO_FID_LOG,*) '*** specific heat (water vapor,pressure) [J/kg/K] : CPvap   = ', CONST_CPvap
    if( IO_L ) write(IO_FID_LOG,*) '*** specific heat (water vapor,volume)   [J/kg/K] : CVvap   = ', CONST_CVvap
    if( IO_L ) write(IO_FID_LOG,*) '*** specific heat (liquid water)         [J/kg/K] : CL      = ', CONST_CL
    if( IO_L ) write(IO_FID_LOG,*) '*** specific heat (ice)                  [J/kg/K] : CI      = ', CONST_CI
    if( IO_L ) write(IO_FID_LOG,*) '*** Rdry / Rvap                                   : EPSvap  = ', CONST_EPSvap
    if( IO_L ) write(IO_FID_LOG,*) '*** 1 / EPSvap - 1                                : EPSTvap = ', CONST_EPSTvap

    if( IO_L ) write(IO_FID_LOG,*) '*** standard pressure  [Pa]                       : Pstd    = ', CONST_Pstd
    if( IO_L ) write(IO_FID_LOG,*) '*** pressure reference [Pa]                       : PRE00   = ', CONST_PRE00
    if( IO_L ) write(IO_FID_LOG,*) '*** standard temperature (15C) [K]                : Tstd    = ', CONST_Tstd
    if( IO_L ) write(IO_FID_LOG,*) '*** temperature reference (0C) [K]                : TEM00   = ', CONST_TEM00

    return
  end subroutine CONST_setup

end module scale_const
!-------------------------------------------------------------------------------
