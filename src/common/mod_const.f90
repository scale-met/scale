!-------------------------------------------------------------------------------
!> module CONSTANT
!!
!! @par Description
!!          Physical constants module
!!
!! @author H.Tomita and SCALE developpers
!!
!! @par History
!! @li      2011-11-11 (H.Yashiro) [new]
!!
!<
!-------------------------------------------------------------------------------
module mod_const
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use mod_stdio, only : &
     IO_FID_LOG, &
     IO_L
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
  !++ included parameters
  !
  include 'inc_precision.h'
  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  real(RP), public, save      :: CONST_PI    = 3.14159265358979E0_RP !< pi
  real(RP), public, save      :: CONST_EPS   = 1.E-34_RP             !< small number

  integer, public, parameter :: CONST_UNDEF2 = -32768            !< undefined value (INT2)
  real(4), public, parameter :: CONST_UNDEF4 = -9.9999E30        !< undefined value (REAL4)
  real(8), public, parameter :: CONST_UNDEF8 = -9.9999D30        !< undefined value (REAL8)

  ! physical constants
  real(RP), public, parameter :: CONST_ERADIUS = 6.37122E+6_RP !< earth radius           [m]
  real(RP), public, parameter :: CONST_EOHM    = 7.292E-5_RP   !< Angular velocity of the Earth [1/s]
  real(RP), public, save      :: CONST_GRAV    = 9.80616E0_RP  !< gravitational constant [m/s2]
  real(RP), public, parameter :: CONST_KARMAN  = 0.4E0_RP      !< karman constant

  real(RP), public, parameter :: CONST_Rdry    = 287.04E0_RP   !< gas constant   (dry)
  real(RP), public, parameter :: CONST_CPdry   = 1003.5E0_RP   !< specific heat  (dry,constant pressure)
  real(RP), public, save      :: CONST_CVdry                !< specific heat  (dry,constant volume)
  real(RP), public, save      :: CONST_RovCP                !< R / Cp = kappa (dry)
  real(RP), public, save      :: CONST_CPovR                !< 1 / kappa      (dry)
  real(RP), public, save      :: CONST_RovCV                !< R / Cv         (dry)
  real(RP), public, save      :: CONST_CPovCV               !< Cp / Cv        (dry)
  real(RP), public, save      :: CONST_CVovCP               !< Cv / Cp        (dry)
  real(RP), public, save      :: CONST_LASPdry              !< g / Cp, Dry adiabatic lapse rate

  real(RP), public, parameter :: CONST_Pstd  = 101325.E0_RP    !< standard pressure [Pa]
  real(RP), public, parameter :: CONST_PRE00 = 100000.E0_RP    !< pressure reference [Pa]
  real(RP), public, parameter :: CONST_Tstd  = 288.15E0_RP     !< standard temperature (15 degree C) [K]
  real(RP), public, parameter :: CONST_TEM00 = 273.15E0_RP     !< temperature reference (0 degree C) [K]

  ! water constants
  real(RP), public, parameter :: CONST_Rvap  = 461.5E0_RP      !< gas constant (water vapor)
  real(RP), public, save      :: CONST_EPSvap               !< Rdry / Rvap
  real(RP), public, save      :: CONST_EPSTvap              !< 1 / epsilon - 1

  real(RP), public, parameter :: CONST_CPvap = 1850.E0_RP      !< specific heat (water vapor, consant pressure)
  real(RP), public, save      :: CONST_CVvap                !< specific heat (water vapor, consant volume)
  real(RP), public, parameter :: CONST_CL    = 4218.E0_RP      !< specific heat (liquid water) 
  real(RP), public, parameter :: CONST_CI    = 2006.E0_RP      !< specific heat (ice)

  real(RP), public, parameter :: CONST_EMELT = 3.4E5_RP
  real(RP), public, parameter :: CONST_TMELT = 273.15E0_RP

  real(RP), public, parameter :: CONST_LH0   = 2.5008E6_RP     !< latent heat of vaporizaion at 0 degree
  real(RP), public, save      :: CONST_LH00
  real(RP), public, parameter :: CONST_LHS0  = 2.8342E6_RP
  real(RP), public, save      :: CONST_LHS00
  real(RP), public, save      :: CONST_LHF0
  real(RP), public, save      :: CONST_LHF00
  real(RP), public, parameter :: CONST_PSAT0 = 610.7E0_RP      !< saturate pressure of water vapor at 0C
  real(RP), public, parameter :: CONST_DWATR = 1000.E0_RP      !< density of water
  real(RP), public, parameter :: CONST_DICE  = 916.8E0_RP      !< density of ice

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
  !> Setup constants
  !-----------------------------------------------------------------------------
  subroutine CONST_setup
    use mod_stdio, only: &
       IO_FID_CONF
    use mod_process, only: &
       PRC_MPIstop
    implicit none

    namelist / PARAM_CONST / &
       CONST_GRAV

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

    CONST_PI    = 4.E0_RP * atan( 1.D0 )
    CONST_EPS   = epsilon(0.E0_RP)

    CONST_CVdry   = CONST_CPdry - CONST_Rdry
    CONST_RovCP   = CONST_Rdry  / CONST_CPdry
    CONST_RovCV   = CONST_Rdry  / CONST_CVdry
    CONST_CPovR   = CONST_CPdry / CONST_Rdry
    CONST_CPovCV  = CONST_CPdry / CONST_CVdry
    CONST_CVovCP  = CONST_CVdry / CONST_CPdry
    CONST_LASPdry = CONST_GRAV / CONST_CPdry

    CONST_EPSvap  = CONST_Rdry / CONST_Rvap
    CONST_EPSTvap = 1.E0_RP / CONST_EPSvap - 1.E0_RP
    CONST_CVvap   = CONST_CPvap - CONST_Rvap

    CONST_LH00    = CONST_LH0  - ( CONST_CPvap - CONST_CL ) * CONST_TEM00
    CONST_LHS00   = CONST_LHS0 - ( CONST_CPvap - CONST_CI ) * CONST_TEM00
    CONST_LHF0    = CONST_LHS0 - CONST_LH0
    CONST_LHF00   = CONST_LHF0 - ( CONST_CL  - CONST_CI ) * CONST_TEM00

    if( IO_L ) write(IO_FID_LOG,*) '*** PI                                   : PI      = ', CONST_PI
    if( IO_L ) write(IO_FID_LOG,*) '*** Small num(REAL8)                     : EPS     = ', CONST_EPS
    if( IO_L ) write(IO_FID_LOG,*) '*** undefined num(INT2)                  : UNDEF2  = ', CONST_UNDEF2
    if( IO_L ) write(IO_FID_LOG,*) '*** undefined num(REAL4)                 : UNDEF4  = ', CONST_UNDEF4
    if( IO_L ) write(IO_FID_LOG,*) '*** undefined num(REAL8)                 : UNDEF8  = ', CONST_UNDEF8
    if( IO_L ) write(IO_FID_LOG,*) '*** radius of the Earth [m]              : ERADIUS = ', CONST_ERADIUS
    if( IO_L ) write(IO_FID_LOG,*) '*** angular velocity of the Earth [1/s]  : EOHM    = ', CONST_EOHM
    if( IO_L ) write(IO_FID_LOG,*) '*** gravitational constant [m/s2]        : GRAV    = ', CONST_GRAV
    if( IO_L ) write(IO_FID_LOG,*) '*** karman constant                      : KARMAN  = ', CONST_KARMAN
    if( IO_L ) write(IO_FID_LOG,*) '*** gas constant   (dry)                 : Rdry    = ', CONST_Rdry
    if( IO_L ) write(IO_FID_LOG,*) '*** specific heat  (dry,pressure)        : CPdry   = ', CONST_CPdry
    if( IO_L ) write(IO_FID_LOG,*) '*** specific heat  (dry,volume)          : Cvdry   = ', CONST_CVdry
    if( IO_L ) write(IO_FID_LOG,*) '*** R / Cp = kappa (dry)                 : RovCP   = ', CONST_RovCP
    if( IO_L ) write(IO_FID_LOG,*) '*** 1 / kappa      (dry)                 : CPovR   = ', CONST_CPovR
    if( IO_L ) write(IO_FID_LOG,*) '*** R / Cv         (dry)                 : RovCV   = ', CONST_RovCV
    if( IO_L ) write(IO_FID_LOG,*) '*** Cp / Cv        (dry)                 : CPovCV  = ', CONST_CPovCV
    if( IO_L ) write(IO_FID_LOG,*) '*** Cv / Cp        (dry)                 : CVovCP  = ', CONST_CVovCP
    if( IO_L ) write(IO_FID_LOG,*) '*** dry adiabatic lapse rate             : LASPdry = ', CONST_LASPdry
    if( IO_L ) write(IO_FID_LOG,*) '*** standard pressure  [Pa]              : Pstd    = ', CONST_Pstd
    if( IO_L ) write(IO_FID_LOG,*) '*** pressure reference [Pa]              : PRE00   = ', CONST_PRE00
    if( IO_L ) write(IO_FID_LOG,*) '*** standard temperature (15C) [K]       : Tstd    = ', CONST_Tstd
    if( IO_L ) write(IO_FID_LOG,*) '*** temperature reference (0C) [K]       : TEM00   = ', CONST_TEM00
    if( IO_L ) write(IO_FID_LOG,*) '*** gas constant (water vapor)           : Rvap    = ', CONST_Rvap
    if( IO_L ) write(IO_FID_LOG,*) '*** Rdry / Rvap = epsilon                : EPSvap  = ', CONST_EPSvap
    if( IO_L ) write(IO_FID_LOG,*) '*** 1 / epsilon - 1                      : EPSTvap = ', CONST_EPSTvap
    if( IO_L ) write(IO_FID_LOG,*) '*** specific heat (water vapor,pressure) : CPvap   = ', CONST_CPvap
    if( IO_L ) write(IO_FID_LOG,*) '*** specific heat (water vapor,volume)   : CVvap   = ', CONST_CVvap
    if( IO_L ) write(IO_FID_LOG,*) '*** specific heat (liquid water)         : CL      = ', CONST_CL
    if( IO_L ) write(IO_FID_LOG,*) '*** specific heat (ice)                  : CI      = ', CONST_CI

    return
  end subroutine CONST_setup

end module mod_const
