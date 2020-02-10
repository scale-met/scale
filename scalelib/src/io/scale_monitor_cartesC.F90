!-------------------------------------------------------------------------------
!> module MONITOR CartesC
!!
!! @par Description
!!          Monitor output module for the cartesianC grid
!!
!! @author Team SCALE
!!
!<
!-------------------------------------------------------------------------------
#include "scalelib.h"
module scale_monitor_cartesc
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use scale_precision
  use scale_io
  use scale_prof
  use scale_atmos_grid_cartesC_index
  use scale_ocean_grid_cartesC_index
  use scale_land_grid_cartesC_index
  use scale_urban_grid_cartesC_index
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedures
  !
  public :: MONITOR_CARTESC_setup

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private procedures
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine MONITOR_CARTESC_setup( &
       dt, &
       ATMOS_do, OCEAN_do, LAND_do, URBAN_do )
    use scale_prc_cartesC, only: &
       PRC_TwoD
    use scale_monitor, only: &
       MONITOR_setup, &
       MONITOR_set_dim
    use scale_atmos_grid_cartesC_real, only: &
       ATMOS_GRID_CARTESC_REAL_AREA,         &
       ATMOS_GRID_CARTESC_REAL_TOTAREA,      &
       ATMOS_GRID_CARTESC_REAL_AREAUY,       &
       ATMOS_GRID_CARTESC_REAL_TOTAREAUY,    &
       ATMOS_GRID_CARTESC_REAL_AREAXV,       &
       ATMOS_GRID_CARTESC_REAL_TOTAREAXV,    &
       ATMOS_GRID_CARTESC_REAL_VOL,          &
       ATMOS_GRID_CARTESC_REAL_TOTVOL,       &
       ATMOS_GRID_CARTESC_REAL_VOLWXY,       &
       ATMOS_GRID_CARTESC_REAL_TOTVOLWXY,    &
       ATMOS_GRID_CARTESC_REAL_VOLZUY,       &
       ATMOS_GRID_CARTESC_REAL_TOTVOLZUY,    &
       ATMOS_GRID_CARTESC_REAL_VOLZXV,       &
       ATMOS_GRID_CARTESC_REAL_TOTVOLZXV,    &
       ATMOS_GRID_CARTESC_REAL_AREAZUY_X,    &
       ATMOS_GRID_CARTESC_REAL_TOTAREAZUY_X, &
       ATMOS_GRID_CARTESC_REAL_AREAZXV_Y,    &
       ATMOS_GRID_CARTESC_REAL_TOTAREAZXV_Y
    use scale_ocean_grid_cartesC_real, only: &
       OCEAN_GRID_CARTESC_REAL_AREA,    &
       OCEAN_GRID_CARTESC_REAL_TOTAREA, &
       OCEAN_GRID_CARTESC_REAL_VOL,     &
       OCEAN_GRID_CARTESC_REAL_TOTVOL
    use scale_land_grid_cartesC_real, only: &
       LAND_GRID_CARTESC_REAL_AREA,    &
       LAND_GRID_CARTESC_REAL_TOTAREA, &
       LAND_GRID_CARTESC_REAL_VOL,     &
       LAND_GRID_CARTESC_REAL_TOTVOL
    use scale_urban_grid_cartesC_real, only: &
       URBAN_GRID_CARTESC_REAL_AREA,    &
       URBAN_GRID_CARTESC_REAL_TOTAREA, &
       URBAN_GRID_CARTESC_REAL_VOL,     &
       URBAN_GRID_CARTESC_REAL_TOTVOL
    implicit none

    real(DP), intent(in) :: dt
    logical,  intent(in) :: ATMOS_do
    logical,  intent(in) :: OCEAN_do
    logical,  intent(in) :: LAND_do
    logical,  intent(in) :: URBAN_do
    !---------------------------------------------------------------------------


    call MONITOR_setup( dt )

    ! atmos
    if ( ATMOS_do ) then
       call MONITOR_set_dim( KA, KS, KE, IA, IS, IE, JA, JS, JE, &
                             "ZXY", 3, &
                             ATMOS_GRID_CARTESC_REAL_AREA(:,:), &
                             ATMOS_GRID_CARTESC_REAL_TOTAREA,    &
                             ATMOS_GRID_CARTESC_REAL_VOL(:,:,:), &
                             ATMOS_GRID_CARTESC_REAL_TOTVOL      )
       call MONITOR_set_dim( KA, KS, KE, IA, IS, IE, JA, JS, JE, &
                             "WXY", 3, &
                             ATMOS_GRID_CARTESC_REAL_AREA(:,:),     &
                             ATMOS_GRID_CARTESC_REAL_TOTAREA,       &
                             ATMOS_GRID_CARTESC_REAL_VOLWXY(:,:,:), &
                             ATMOS_GRID_CARTESC_REAL_TOTVOLWXY      )
       call MONITOR_set_dim( KA, KS, KE, IA, IS, IE, JA, JS, JE, &
                             "ZUY", 3, &
                             ATMOS_GRID_CARTESC_REAL_AREAUY(:,:),   &
                             ATMOS_GRID_CARTESC_REAL_TOTAREAUY,     &
                             ATMOS_GRID_CARTESC_REAL_VOLZUY(:,:,:), &
                             ATMOS_GRID_CARTESC_REAL_TOTVOLZUY      )
       call MONITOR_set_dim( KA, KS, KE, IA, IS, IE, JA, JS, JE, &
                             "ZXV", 3, &
                             ATMOS_GRID_CARTESC_REAL_AREAXV(:,:),   &
                             ATMOS_GRID_CARTESC_REAL_TOTAREAXV,     &
                             ATMOS_GRID_CARTESC_REAL_VOLZXV(:,:,:), &
                             ATMOS_GRID_CARTESC_REAL_TOTVOLZXV      )
       call MONITOR_set_dim( KA, KS, KE, IA, IS, IE, JA, JS, JE, &
                             "XY", 2, &
                             ATMOS_GRID_CARTESC_REAL_AREA(:,:), &
                             ATMOS_GRID_CARTESC_REAL_TOTAREA    )
       if ( PRC_TwoD ) then
          call MONITOR_set_dim( KA, KS, KE, KA, KS, KE, JA, JS, JE, &
                                "ZY-W", 2, &
                                ATMOS_GRID_CARTESC_REAL_AREAZUY_X(:,IS,:), &
                                ATMOS_GRID_CARTESC_REAL_TOTAREAZUY_X(IS)   )
       else
          call MONITOR_set_dim( KA, KS, KE, KA, KS, KE, JA, JS, JE, &
                                "ZY-W", 2, &
                                ATMOS_GRID_CARTESC_REAL_AREAZUY_X(:,IS-1,:), &
                                ATMOS_GRID_CARTESC_REAL_TOTAREAZUY_X(IS-1)   )
       end if
       call MONITOR_set_dim( KA, KS, KE, KA, KS, KE, JA, JS, JE, &
                             "ZY-E", 2, &
                             ATMOS_GRID_CARTESC_REAL_AREAZUY_X(:,IE,:), &
                             ATMOS_GRID_CARTESC_REAL_TOTAREAZUY_X(IE)    )
       call MONITOR_set_dim( KA, KS, KE, KA, KS, KE, IA, IS, IE, &
                             "ZX-S", 2, &
                             ATMOS_GRID_CARTESC_REAL_AREAZXV_Y(:,:,JS-1), &
                             ATMOS_GRID_CARTESC_REAL_TOTAREAZXV_Y(JS-1)   )
       call MONITOR_set_dim( KA, KS, KE, KA, KS, KE, IA, IS, IE, &
                             "ZX-N", 2, &
                             ATMOS_GRID_CARTESC_REAL_AREAZXV_Y(:,:,JE), &
                             ATMOS_GRID_CARTESC_REAL_TOTAREAZXV_Y(JE)   )
    end if

    ! ocean
    if ( OCEAN_do ) then
       call MONITOR_set_dim( OKA, OKS, OKE, OIA, OIS, OIE, OJA, OJS, OJE, &
                             "OXY", 3, &
                             OCEAN_GRID_CARTESC_REAL_AREA(:,:),  &
                             OCEAN_GRID_CARTESC_REAL_TOTAREA,    &
                             OCEAN_GRID_CARTESC_REAL_VOL(:,:,:), &
                             OCEAN_GRID_CARTESC_REAL_TOTVOL      )
    end if

    ! land
    if ( LAND_do ) then
       call MONITOR_set_dim( LKA, LKS, LKE, LIA, LIS, LIE, LJA, LJS, LJE, &
                             "LXY", 3, &
                             LAND_GRID_CARTESC_REAL_AREA(:,:),  &
                             LAND_GRID_CARTESC_REAL_TOTAREA,    &
                             LAND_GRID_CARTESC_REAL_VOL(:,:,:), &
                             LAND_GRID_CARTESC_REAL_TOTVOL      )
    end if

    ! urban
    if ( URBAN_do ) then
       call MONITOR_set_dim( UKA, UKS, UKE, UIA, UIS, UIE, UJA, UJS, UJE, &
                             "UXY", 3, &
                             URBAN_GRID_CARTESC_REAL_AREA(:,:),  &
                             URBAN_GRID_CARTESC_REAL_TOTAREA,    &
                             URBAN_GRID_CARTESC_REAL_VOL(:,:,:), &
                             URBAN_GRID_CARTESC_REAL_TOTVOL      )
    end if

    return
  end subroutine MONITOR_CARTESC_setup

end module scale_monitor_cartesc
