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
module scale_monitor_cartesc
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use scale_precision
  use scale_stdio
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
  subroutine MONITOR_CARTESC_setup( dt )
    use scale_monitor, only: &
       MONITOR_setup, &
       MONITOR_set_dim
    use scale_atmos_grid_cartesC_real, only: &
       ATMOS_GRID_CARTESC_REAL_AREA,      &
       ATMOS_GRID_CARTESC_REAL_TOTAREA,   &
       ATMOS_GRID_CARTESC_REAL_AREAUY,    &
       ATMOS_GRID_CARTESC_REAL_TOTAREAUY, &
       ATMOS_GRID_CARTESC_REAL_AREAXV,    &
       ATMOS_GRID_CARTESC_REAL_TOTAREAXV, &
       ATMOS_GRID_CARTESC_REAL_VOL,       &
       ATMOS_GRID_CARTESC_REAL_TOTVOL,    &
       ATMOS_GRID_CARTESC_REAL_VOLWXY,    &
       ATMOS_GRID_CARTESC_REAL_TOTVOLWXY, &
       ATMOS_GRID_CARTESC_REAL_VOLZUY,    &
       ATMOS_GRID_CARTESC_REAL_TOTVOLZUY, &
       ATMOS_GRID_CARTESC_REAL_VOLZXV,    &
       ATMOS_GRID_CARTESC_REAL_TOTVOLZXV
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

    integer :: ierr
    integer :: n
    !---------------------------------------------------------------------------


    call MONITOR_setup( dt )

    ! atmos
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

    ! ocean
    call MONITOR_set_dim( OKA, OKS, OKE, OIA, OIS, OIE, OJA, OJS, OJE, &
                          "OXY", 3, &
                          OCEAN_GRID_CARTESC_REAL_AREA(:,:),  &
                          OCEAN_GRID_CARTESC_REAL_TOTAREA,    &
                          OCEAN_GRID_CARTESC_REAL_VOL(:,:,:), &
                          OCEAN_GRID_CARTESC_REAL_TOTVOL      )

    ! land
    call MONITOR_set_dim( LKA, LKS, LKE, LIA, LIS, LIE, LJA, LJS, LJE, &
                          "LXY", 3, &
                          LAND_GRID_CARTESC_REAL_AREA(:,:),  &
                          LAND_GRID_CARTESC_REAL_TOTAREA,    &
                          LAND_GRID_CARTESC_REAL_VOL(:,:,:), &
                          LAND_GRID_CARTESC_REAL_TOTVOL      )

    ! urban
    call MONITOR_set_dim( UKA, UKS, UKE, UIA, UIS, UIE, UJA, UJS, UJE, &
                          "UXY", 3, &
                          URBAN_GRID_CARTESC_REAL_AREA(:,:),  &
                          URBAN_GRID_CARTESC_REAL_TOTAREA,    &
                          URBAN_GRID_CARTESC_REAL_VOL(:,:,:), &
                          URBAN_GRID_CARTESC_REAL_TOTVOL      )

    return
  end subroutine MONITOR_CARTESC_setup

end module scale_monitor_cartesc
