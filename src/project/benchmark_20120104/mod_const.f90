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
  real(8), public, save      :: CONST_PI    = 3.14159265358979D0 !< pi
  real(8), public, save      :: CONST_EPS   = 1.D-34             !< small number

  integer, public, save      :: CONST_UNDEF2 = -32768            !< undefined value
  real(4), public, save      :: CONST_UNDEF4 = -999.E20          !< undefined value
  real(8), public, save      :: CONST_UNDEF8 = -999.D20          !< undefined value

  ! physical constants
  real(8), public, save      :: CONST_ERADIUS = 6.37122D+6 !< earth radius
  real(8), public, save      :: CONST_GRAV    = 9.80616D0  !< gravitational constant [m/s**2]
  real(8), public, save      :: CONST_KARMAN  = 0.4D0      !< karman constant

  real(8), public, save      :: CONST_Rdry    = 287.04D0   !< gas constant (dry)
  real(8), public, save      :: CONST_CPdry   = 1005.7D0   !< specific heat (dry,constant pressure)
  real(8), public, save      :: CONST_CVdry                !< specific heat (dry,constant volume)
  real(8), public, save      :: CONST_RovCP                !< kappa
  real(8), public, save      :: CONST_CPovR                !< 1 / kappa
  real(8), public, save      :: CONST_CPovCV               !< Cp / Cv
  real(8), public, save      :: CONST_CVovCP               !< Cv / Cp

  real(8), public, parameter :: CONST_Pstd  = 101325.D0    !< standard pressure [Pa]
  real(8), public, parameter :: CONST_Tstd  = 288.15D0     !< standard temperature (15 degree C) [K]
  real(8), public, parameter :: CONST_PRE00 = 100000.D0    ! pressure reference [Pa]
  real(8), public, parameter :: CONST_TEM00 = 273.15D0     ! temperature reference (0 degree C) [K]

  ! water constants
  real(8), public, parameter :: CONST_Rvap  = 461.5D0     ! gas constant (water vapor)
  real(8), public, save      :: CONST_EPSvap              ! RD / Rvap
  real(8), public, save      :: CONST_EPSTvap             ! 1 / epsilon - 1

  real(8), public, parameter :: CONST_CPvap = 1850.D0     ! specific heat (water vapor, consant pressure)
  real(8), public, save      :: CONST_CVvap               ! specific heat (water vapor, consant volume)
  real(8), public, parameter :: CONST_CL    = 4218.D0     ! specific heat (liquid water) 
  real(8), public, parameter :: CONST_CI    = 2006.D0     ! specific heat (ice)

  real(8), public, parameter :: CONST_EMELT = 3.4D5
  real(8), public, parameter :: CONST_TMELT = 273.15D0

  real(8), public, parameter :: CONST_LH0   = 2.5008D6    ! latent heat of vaporizaion at 0 degree
  real(8), public, save      :: CONST_LH00
  real(8), public, parameter :: CONST_LHS0  = 2.8342D6
  real(8), public, save      :: CONST_LHS00
  real(8), public, save      :: CONST_LHF0
  real(8), public, save      :: CONST_LHF00
  real(8), public, parameter :: CONST_PSAT0 = 610.7D0     ! saturate pressure of water vapor at 0C
  real(8), public, parameter :: CONST_DWATR = 1000.D0     ! density of water
  real(8), public, parameter :: CONST_DICE  = 916.8D0     ! density of ice

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
       IO_FID_CONF, &
       IO_FID_LOG,  &
       IO_L
    use mod_process, only: &
       PRC_MPIstop
    implicit none

    NAMELIST / PARAM_CONST / &
       CONST_UNDEF8, &
       CONST_GRAV,   &
       CONST_Rdry,   &
       CONST_CPdry

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

    CONST_PI    = 4.D0 * atan( 1.D0 )
    CONST_EPS   = epsilon(0.D0)

    CONST_CVdry   = CONST_CPdry - CONST_Rdry

    CONST_RovCP   = CONST_Rdry  / CONST_CPdry
    CONST_CPovR   = CONST_CPdry / CONST_Rdry
    CONST_CPovCV  = CONST_CPdry / CONST_CVdry
    CONST_CVovCP  = CONST_CVdry / CONST_CPdry

    CONST_EPSvap  = CONST_Rdry / CONST_Rvap
    CONST_EPSTvap = 1.D0 / CONST_EPSvap - 1.D0

    CONST_CVvap   = CONST_CPvap - CONST_Rvap

    CONST_LH00    = CONST_LH0  - ( CONST_CPvap - CONST_CL ) * CONST_TEM00
    CONST_LHS00   = CONST_LHS0 - ( CONST_CPvap - CONST_CI ) * CONST_TEM00
    CONST_LHF0    = CONST_LHS0 - CONST_LH0
    CONST_LHF00   = CONST_LHF0 - ( CONST_CL  - CONST_CI ) * CONST_TEM00

    if( IO_L ) write(IO_FID_LOG,*) '*** PI          :', CONST_PI
    if( IO_L ) write(IO_FID_LOG,*) '*** EPS(REAL8)  :', CONST_EPS
    if( IO_L ) write(IO_FID_LOG,*) '*** UNDEF(REAL8):', CONST_UNDEF8
    if( IO_L ) write(IO_FID_LOG,*) '*** UNDEF(REAL4):', CONST_UNDEF4
    if( IO_L ) write(IO_FID_LOG,*) '*** UNDEF(INT2) :', CONST_UNDEF2
    if( IO_L ) write(IO_FID_LOG,*) '*** GRAV        :', CONST_GRAV
    if( IO_L ) write(IO_FID_LOG,*) '*** R(dry)      :', CONST_Rdry
    if( IO_L ) write(IO_FID_LOG,*) '*** CP(dry)     :', CONST_CPdry
    if( IO_L ) write(IO_FID_LOG,*) '*** CV(dry)     :', CONST_CVdry
    if( IO_L ) write(IO_FID_LOG,*) '*** Std. Pres.  :', CONST_Pstd

    return
  end subroutine CONST_setup

end module mod_const
