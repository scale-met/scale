!-------------------------------------------------------------------------------
!> module ocean / grid / cartesianC / real
!!
!! @par Description
!!          Grid module for cartesian coordinate for ocean
!!
!! @author Team SCALE
!!
!<
!-------------------------------------------------------------------------------
#include "scalelib.h"
module scale_ocean_grid_cartesC_real
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_io
  use scale_prof
  use scale_atmos_grid_cartesC_index
  use scale_ocean_grid_cartesC_index
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: OCEAN_GRID_CARTESC_REAL_setup

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  real(RP), public, allocatable :: OCEAN_GRID_CARTESC_REAL_AREA(:,:)   !< area   of grid cell
  real(RP), public              :: OCEAN_GRID_CARTESC_REAL_TOTAREA     !< total area
  real(RP), public, allocatable :: OCEAN_GRID_CARTESC_REAL_VOL (:,:,:) !< volume of grid cell
  real(RP), public              :: OCEAN_GRID_CARTESC_REAL_TOTVOL      !< total volume

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
  !> Setup area and volume
  subroutine OCEAN_GRID_CARTESC_REAL_setup
    use scale_atmos_grid_cartesC_real, only: &
       ATMOS_GRID_CARTESC_REAL_AREA, &
       ATMOS_GRID_CARTESC_REAL_TOTAREA
    use scale_ocean_grid_cartesC, only: &
       OCEAN_GRID_CARTESC_CDZ
    use scale_file_cartesC, only: &
       FILE_CARTESC_set_coordinates_ocean
    integer :: k, i, j

    ! at this moment, horizontal grid is identical to that of the atmosphere
    allocate( OCEAN_GRID_CARTESC_REAL_AREA(    OIA,OJA) )
    allocate( OCEAN_GRID_CARTESC_REAL_VOL (OKA,OIA,OJA) )
    OCEAN_GRID_CARTESC_REAL_AREA(:,:) = ATMOS_GRID_CARTESC_REAL_AREA(:,:)
    OCEAN_GRID_CARTESC_REAL_TOTAREA   = ATMOS_GRID_CARTESC_REAL_TOTAREA

    do j = 1,   OJA
    do i = 1,   OIA
    do k = OKS, OKE
       OCEAN_GRID_CARTESC_REAL_VOL(k,i,j) = OCEAN_GRID_CARTESC_REAL_AREA(i,j) * OCEAN_GRID_CARTESC_CDZ(k)
    enddo
    enddo
    enddo

    OCEAN_GRID_CARTESC_REAL_TOTVOL = 0.0_RP
    do j = OJS, OJE
    do i = OIS, OIE
    do k = OKS, OKE
       OCEAN_GRID_CARTESC_REAL_TOTVOL = OCEAN_GRID_CARTESC_REAL_TOTVOL + OCEAN_GRID_CARTESC_REAL_VOL(k,i,j)
    end do
    end do
    end do

    call FILE_CARTESC_set_coordinates_ocean( OCEAN_GRID_CARTESC_REAL_VOL(:,:,:) ) ! [IN]

    return
  end subroutine OCEAN_GRID_CARTESC_REAL_setup

end module scale_ocean_grid_cartesC_real
