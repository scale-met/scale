!-------------------------------------------------------------------------------
!> module land / grid / cartesianC / real
!!
!! @par Description
!!          Grid module for cartesian coordinate for land
!!
!! @author Team SCALE
!<
!-------------------------------------------------------------------------------
#include "scalelib.h"
module scale_land_grid_cartesC_real
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_io
  use scale_prof
  use scale_atmos_grid_cartesC_index
  use scale_land_grid_cartesC_index
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: LAND_GRID_CARTESC_REAL_setup
  public :: LAND_GRID_CARTESC_REAL_finalize
  public :: LAND_GRID_CARTESC_REAL_set_areavol

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  real(RP), public, allocatable :: LAND_GRID_CARTESC_REAL_AREA(:,:)   !< area   of grid cell
  real(RP), public              :: LAND_GRID_CARTESC_REAL_TOTAREA     !< total area
  real(RP), public, allocatable :: LAND_GRID_CARTESC_REAL_VOL (:,:,:) !< volume of grid cell
  real(RP), public              :: LAND_GRID_CARTESC_REAL_TOTVOL      !< total volume

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
  !> Setup real grid
  subroutine LAND_GRID_CARTESC_REAL_setup

    ! at this moment, horizontal grid is identical to that of the atmosphere
    allocate( LAND_GRID_CARTESC_REAL_AREA(    LIA,LJA) )
    allocate( LAND_GRID_CARTESC_REAL_VOL (LKA,LIA,LJA) )

    return
  end subroutine LAND_GRID_CARTESC_REAL_setup

  !-----------------------------------------------------------------------------
  !> Finalize
  subroutine LAND_GRID_CARTESC_REAL_finalize
    implicit none

    deallocate( LAND_GRID_CARTESC_REAL_AREA  )
    deallocate( LAND_GRID_CARTESC_REAL_VOL   )

    return
  end subroutine LAND_GRID_CARTESC_REAL_finalize
  !-----------------------------------------------------------------------------

  subroutine LAND_GRID_CARTESC_REAL_set_areavol
    use scale_atmos_grid_cartesC_real, only: &
       ATMOS_GRID_CARTESC_REAL_AREA
    use scale_land_grid_cartesC, only: &
       LAND_GRID_CARTESC_CDZ
    use scale_file_cartesC, only: &
       FILE_CARTESC_set_coordinates_land
    use scale_landuse, only: &
       LANDUSE_fact_land

    integer :: k, i, j

    do j = 1, LJA
    do i = 1, LIA
       LAND_GRID_CARTESC_REAL_AREA(i,j) = ATMOS_GRID_CARTESC_REAL_AREA(i,j) * LANDUSE_fact_land(i,j)
    end do
    end do

    LAND_GRID_CARTESC_REAL_TOTAREA = 0.0_RP
    do j = LJS, LJE
    do i = LIS, LIE
       LAND_GRID_CARTESC_REAL_TOTAREA = LAND_GRID_CARTESC_REAL_TOTAREA + LAND_GRID_CARTESC_REAL_AREA(i,j)
    end do
    end do

    do j = 1,   LJA
    do i = 1,   LIA
    do k = LKS, LKE
       LAND_GRID_CARTESC_REAL_VOL(k,i,j) = LAND_GRID_CARTESC_REAL_AREA(i,j) * LAND_GRID_CARTESC_CDZ(k)
    enddo
    enddo
    enddo

    LAND_GRID_CARTESC_REAL_TOTVOL = 0.0_RP
    do j = LJS, LJE
    do i = LIS, LIE
    do k = LKS, LKE
       LAND_GRID_CARTESC_REAL_TOTVOL = LAND_GRID_CARTESC_REAL_TOTVOL + LAND_GRID_CARTESC_REAL_VOL(k,i,j)
    end do
    end do
    end do

    call FILE_CARTESC_set_coordinates_land( LAND_GRID_CARTESC_REAL_VOL(:,:,:) ) ! [IN]

    return
  end subroutine LAND_GRID_CARTESC_REAL_set_areavol

end module scale_land_grid_cartesC_real
