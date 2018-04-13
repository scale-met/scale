!-------------------------------------------------------------------------------
!> module LAND / Physics Matsiro model
!!
!! @par Description
!!          matsiro-type land physics module
!!
!! @author Team SCALE
!<
!-------------------------------------------------------------------------------
#include "scalelib.h"
module scale_land_phy_matsiro
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_io
  use scale_prof
  use scale_debug
  use scale_land_grid_cartesC_index
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: LAND_PHY_MATSIRO_setup
  public :: LAND_PHY_MATSIRO

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
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
  subroutine LAND_PHY_MATSIRO_setup( LAND_TYPE )
    use scale_prc, only: &
       PRC_abort
    implicit none

    character(len=*), intent(in) :: LAND_TYPE

    logical :: dummy

    namelist / PARAM_LAND_MATSIRO / &
       dummy

    integer :: ierr
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("LAND_PHY_MATSIRO_setup",*) 'Setup'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_LAND_MATSIRO,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("LAND_PHY_MATSIRO_setup",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("LAND_PHY_MATSIRO_setup",*) 'Not appropriate names in namelist PARAM_LAND_MATSIRO. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_LAND_MATSIRO)

    LOG_ERROR("LAND_PHY_MATSIRO", *) "This scheme is under development"
    call PRC_abort

    return
  end subroutine LAND_PHY_MATSIRO_setup

  !-----------------------------------------------------------------------------
  !> Physical processes for land submodel
  subroutine LAND_PHY_MATSIRO( &
       LAND_TEMP_t,       &
       LAND_WATER_t,      &
       LAND_TEMP,         &
       LAND_WATER,        &
       LAND_WaterLimit,   &
       LAND_ThermalCond,  &
       LAND_HeatCapacity, &
       LAND_WaterDiff,    &
       LAND_SFLX_GH,      &
       LAND_SFLX_prec,    &
       LAND_SFLX_evap,    &
       CDZ,               &
       dt                 )
    use scale_atmos_grid_cartesC_index
    implicit none

    ! arguments
    real(RP), intent(out) :: LAND_TEMP_t      (LKMAX,LIA,LJA)
    real(RP), intent(out) :: LAND_WATER_t     (LKMAX,LIA,LJA)

    real(RP), intent(in)  :: LAND_TEMP        (LKMAX,LIA,LJA)
    real(RP), intent(in)  :: LAND_WATER       (LKMAX,LIA,LJA)
    real(RP), intent(in)  :: LAND_WaterLimit  (LIA,LJA)
    real(RP), intent(in)  :: LAND_ThermalCond (LIA,LJA)
    real(RP), intent(in)  :: LAND_HeatCapacity(LIA,LJA)
    real(RP), intent(in)  :: LAND_WaterDiff   (LIA,LJA)
    real(RP), intent(in)  :: LAND_SFLX_GH     (LIA,LJA)
    real(RP), intent(in)  :: LAND_SFLX_prec   (LIA,LJA)
    real(RP), intent(in)  :: LAND_SFLX_evap   (LIA,LJA)
    real(RP), intent(in)  :: CDZ              (LKMAX)
    real(DP), intent(in)  :: dt
    !---------------------------------------------------------------------------

    LOG_PROGRESS(*) 'land / physics / Matsiro'

    LAND_TEMP_t (:,:,:) = 0.0_RP
    LAND_WATER_t(:,:,:) = 0.0_RP

    return
  end subroutine LAND_PHY_MATSIRO

end module scale_land_phy_matsiro
