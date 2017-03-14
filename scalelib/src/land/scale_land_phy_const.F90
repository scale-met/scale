!-------------------------------------------------------------------------------
!> module LAND / Physics Constant model
!!
!! @par Description
!!          constant land physics module
!!
!! @author Team SCALE
!<
!-------------------------------------------------------------------------------
module scale_land_phy_const
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_stdio
  use scale_prof
  use scale_debug
  use scale_grid_index
  use scale_land_grid_index
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: LAND_PHY_CONST_setup
  public :: LAND_PHY_CONST

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
  subroutine LAND_PHY_CONST_setup( LAND_TYPE )
    implicit none

    character(len=*), intent(in) :: LAND_TYPE
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '++++++ Module[CONST] / Categ[LAND PHY] / Origin[SCALElib]'

    return
  end subroutine LAND_PHY_CONST_setup

  !-----------------------------------------------------------------------------
  !> Physical processes for land submodel
  subroutine LAND_PHY_CONST( &
       TEMP_t,       &
       WATER_t,      &
       TEMP,         &
       WATER,        &
       WaterLimit,   &
       ThermalCond,  &
       HeatCapacity, &
       WaterDiff,    &
       SFLX_GH,      &
       SFLX_prec,    &
       SFLX_evap,    &
       CDZ,          &
       dt            )
    implicit none

    ! arguments
    real(RP), intent(out) :: TEMP_t      (LKMAX,IA,JA)
    real(RP), intent(out) :: WATER_t     (LKMAX,IA,JA)

    real(RP), intent(in)  :: TEMP        (LKMAX,IA,JA)
    real(RP), intent(in)  :: WATER       (LKMAX,IA,JA)
    real(RP), intent(in)  :: WaterLimit  (IA,JA)
    real(RP), intent(in)  :: ThermalCond (IA,JA)
    real(RP), intent(in)  :: HeatCapacity(IA,JA)
    real(RP), intent(in)  :: WaterDiff   (IA,JA)
    real(RP), intent(in)  :: SFLX_GH     (IA,JA)
    real(RP), intent(in)  :: SFLX_prec   (IA,JA)
    real(RP), intent(in)  :: SFLX_evap   (IA,JA)
    real(RP), intent(in)  :: CDZ         (LKMAX)
    real(DP), intent(in)  :: dt
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*) '*** Land physics step: Const'

    TEMP_t (:,:,:) = 0.0_RP
    WATER_t(:,:,:) = 0.0_RP

    return
  end subroutine LAND_PHY_CONST

end module scale_land_phy_const
