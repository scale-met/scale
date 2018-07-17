!-------------------------------------------------------------------------------
!> module ATMOS admin
!!
!! @par Description
!!          Atmosphere submodel administrator
!!
!! @author Team SCALE
!<
!-------------------------------------------------------------------------------
#include "scalelib.h"
module mod_atmos_admin
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_io
  use scale_prof
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: ATMOS_ADMIN_setup
  public :: ATMOS_ADMIN_getscheme

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  logical,                public :: ATMOS_do          = .true. ! main switch for the model

  character(len=H_SHORT), public :: ATMOS_DYN_TYPE    = 'NONE'
  character(len=H_SHORT), public :: ATMOS_PHY_MP_TYPE = 'NONE'
  character(len=H_SHORT), public :: ATMOS_PHY_AE_TYPE = 'NONE'
  character(len=H_SHORT), public :: ATMOS_PHY_CH_TYPE = 'NONE'
  character(len=H_SHORT), public :: ATMOS_PHY_RD_TYPE = 'NONE'
  character(len=H_SHORT), public :: ATMOS_PHY_SF_TYPE = 'NONE'
  character(len=H_SHORT), public :: ATMOS_PHY_TB_TYPE = 'NONE'
  character(len=H_SHORT), public :: ATMOS_PHY_BL_TYPE = 'NONE'
  character(len=H_SHORT), public :: ATMOS_PHY_CP_TYPE = 'NONE'

  character(len=H_SHORT), public :: ATMOS_PHY_PRECIP_TYPE = 'Upwind-Euler'

  logical,                public :: ATMOS_USE_AVERAGE = .false.

  logical,                public :: ATMOS_sw_dyn
  logical,                public :: ATMOS_sw_phy_mp
  logical,                public :: ATMOS_sw_phy_ae
  logical,                public :: ATMOS_sw_phy_ch
  logical,                public :: ATMOS_sw_phy_rd
  logical,                public :: ATMOS_sw_phy_sf
  logical,                public :: ATMOS_sw_phy_tb
  logical,                public :: ATMOS_sw_phy_bl
  logical,                public :: ATMOS_sw_phy_cp

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
  subroutine ATMOS_ADMIN_setup
    use scale_prc, only: &
       PRC_abort
    implicit none

    namelist / PARAM_ATMOS / &
       ATMOS_do,          &
       ATMOS_DYN_TYPE,    &
       ATMOS_PHY_MP_TYPE, &
       ATMOS_PHY_AE_TYPE, &
       ATMOS_PHY_CH_TYPE, &
       ATMOS_PHY_RD_TYPE, &
       ATMOS_PHY_SF_TYPE, &
       ATMOS_PHY_TB_TYPE, &
       ATMOS_PHY_BL_TYPE, &
       ATMOS_PHY_CP_TYPE, &
       ATMOS_PHY_PRECIP_TYPE, &
       ATMOS_USE_AVERAGE

    integer :: ierr
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("ATMOS_ADMIN_setup",*) 'Setup'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_ATMOS,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("ATMOS_ADMIN_setup",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("ATMOS_ADMIN_setup",*) 'Not appropriate names in namelist PARAM_ATMOS. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_ATMOS)

    !-----< module component check >-----

    LOG_NEWLINE
    LOG_INFO("ATMOS_ADMIN_setup",*) 'Atmosphere model components '

    if ( ATMOS_do ) then
       LOG_INFO_CONT(*) 'Atmosphere model       : ON'
    else
       LOG_INFO_CONT(*) 'Atmosphere model       : OFF'
    endif

    LOG_INFO_CONT(*) 'Dynamics...'

    if    ( ATMOS_DYN_TYPE == 'OFF' ) then
       LOG_INFO_CONT(*) '+ Dynamical core       : OFF'
       LOG_INFO_CONT(*) '+ Advection            : OFF'
       ATMOS_sw_dyn = .false.
    elseif( ATMOS_DYN_TYPE == 'NONE' ) then
       ! The advection is disbled
       ! The tendencies calculated by physical processed are added
       LOG_INFO_CONT(*) '+ Dynamical core       : ON, ', trim(ATMOS_DYN_TYPE)
       LOG_INFO_CONT(*) '+ Advection            : OFF'
       ATMOS_sw_dyn = .true.
    else ! default
       LOG_INFO_CONT(*) '+ Dynamical core       : ON, ', trim(ATMOS_DYN_TYPE)
       LOG_INFO_CONT(*) '+ Advection            : ON'
       ATMOS_sw_dyn = .true.
    endif

    LOG_INFO_CONT(*) 'Physics...'

    if ( ATMOS_PHY_MP_TYPE /= 'OFF' .AND. ATMOS_PHY_MP_TYPE /= 'NONE' ) then
       LOG_INFO_CONT(*) '+ Cloud Microphysics   : ON, ', trim(ATMOS_PHY_MP_TYPE)
       ATMOS_sw_phy_mp = .true.
    else
       LOG_INFO_CONT(*) '+ Cloud Microphysics   : OFF'
       ATMOS_sw_phy_mp = .false.
    endif

    if ( ATMOS_PHY_AE_TYPE /= 'OFF' .AND. ATMOS_PHY_AE_TYPE /= 'NONE' ) then
       LOG_INFO_CONT(*) '+ Aerosol Microphysics : ON, ', trim(ATMOS_PHY_AE_TYPE)
       ATMOS_sw_phy_ae = .true.
    else
       LOG_INFO_CONT(*) '+ Aerosol Microphysics : OFF'
       ATMOS_sw_phy_ae = .false.
    endif

    if ( ATMOS_PHY_CH_TYPE /= 'OFF' .AND. ATMOS_PHY_CH_TYPE /= 'NONE' ) then
       LOG_INFO_CONT(*) '+ Chemistry            : ON, ', trim(ATMOS_PHY_CH_TYPE)
       ATMOS_sw_phy_ch = .true.
    else
       LOG_INFO_CONT(*) '+ Chemistry            : OFF'
       ATMOS_sw_phy_ch = .false.
    endif

    if ( ATMOS_PHY_RD_TYPE /= 'OFF' .AND. ATMOS_PHY_RD_TYPE /= 'NONE' ) then
       LOG_INFO_CONT(*) '+ Radiative transfer   : ON, ', trim(ATMOS_PHY_RD_TYPE)
       ATMOS_sw_phy_rd = .true.
    else
       LOG_INFO_CONT(*) '+ Radiative transfer   : OFF'
       ATMOS_sw_phy_rd = .false.
    endif

    if ( ATMOS_PHY_SF_TYPE /= 'OFF' .AND. ATMOS_PHY_SF_TYPE /= 'NONE' ) then
       LOG_INFO_CONT(*) '+ Surface Flux         : ON, ', trim(ATMOS_PHY_SF_TYPE)
       ATMOS_sw_phy_sf = .true.
    else
       LOG_INFO_CONT(*) '+ Surface Flux         : OFF'
       ATMOS_sw_phy_sf = .false.
    endif

    if ( ATMOS_PHY_TB_TYPE /= 'OFF' .AND. ATMOS_PHY_TB_TYPE /= 'NONE' ) then
       LOG_INFO_CONT(*) '+ Sub-grid Turbulence  : ON, ', trim(ATMOS_PHY_TB_TYPE)
       ATMOS_sw_phy_tb = .true.
    else
       LOG_INFO_CONT(*) '+ Sub-grid Turbulence  : OFF'
       ATMOS_sw_phy_tb = .false.
    endif

    if ( ATMOS_PHY_BL_TYPE /= 'OFF' .AND. ATMOS_PHY_BL_TYPE /= 'NONE' ) then
       LOG_INFO_CONT(*) '+ PBL Turbulence       : ON, ', trim(ATMOS_PHY_BL_TYPE)
       ATMOS_sw_phy_bl = .true.
    else
       LOG_INFO_CONT(*) '+ PBL Turbulence       : OFF'
       ATMOS_sw_phy_bl = .false.
    endif

    if ( ATMOS_PHY_CP_TYPE /= 'OFF' .AND. ATMOS_PHY_CP_TYPE /= 'NONE' ) then
       LOG_INFO_CONT(*) '+ Convection Param.    : ON, ', trim(ATMOS_PHY_CP_TYPE)
       ATMOS_sw_phy_cp = .true.
    else
       LOG_INFO_CONT(*) '+ Convection Param.    : OFF'
       ATMOS_sw_phy_cp = .false.
    endif

    if ( ATMOS_USE_AVERAGE ) then
       LOG_INFO_CONT(*) '+ Use time-averaging value for physics? : YES'
    else
       LOG_INFO_CONT(*) '+ Use time-averaging value for physics? : NO'
    endif

    return
  end subroutine ATMOS_ADMIN_setup

  !-----------------------------------------------------------------------------
  !> Get name of scheme for each component
  subroutine ATMOS_ADMIN_getscheme( &
       component_name, &
       scheme_name     )
    use scale_prc, only: &
       PRC_abort
    implicit none

    character(len=*),       intent(in)  :: component_name
    character(len=H_SHORT), intent(out) :: scheme_name
    !---------------------------------------------------------------------------

    select case(component_name)
    case("DYN")
       scheme_name = ATMOS_DYN_TYPE
    case("PHY_MP")
       scheme_name = ATMOS_PHY_MP_TYPE
    case("PHY_AE")
       scheme_name = ATMOS_PHY_AE_TYPE
    case("PHY_CH")
       scheme_name = ATMOS_PHY_CH_TYPE
    case("PHY_RD")
       scheme_name = ATMOS_PHY_RD_TYPE
    case("PHY_SF")
       scheme_name = ATMOS_PHY_SF_TYPE
    case("PHY_TB")
       scheme_name = ATMOS_PHY_TB_TYPE
    case("PHY_BL")
       scheme_name = ATMOS_PHY_BL_TYPE
    case("PHY_CP")
       scheme_name = ATMOS_PHY_CP_TYPE
    case default
       LOG_ERROR("ATMOS_ADMIN_getscheme",*) 'Unsupported component_name. Check!', trim(component_name)
       call PRC_abort
    end select

    return
  end subroutine ATMOS_ADMIN_getscheme

end module mod_atmos_admin
