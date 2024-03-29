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
  public :: OCEAN_GRID_CARTESC_REAL_finalize
  public :: OCEAN_GRID_CARTESC_REAL_set_areavol

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  real(RP), public, allocatable :: OCEAN_GRID_CARTESC_REAL_AREA(:,:)   !< area   of grid cell
  real(RP), public              :: OCEAN_GRID_CARTESC_REAL_TOTAREA     !< total area
  real(RP), public, allocatable :: OCEAN_GRID_CARTESC_REAL_VOL (:,:,:) !< volume of grid cell
  real(RP), public              :: OCEAN_GRID_CARTESC_REAL_TOTVOL      !< total volume
  !$acc declare create(OCEAN_GRID_CARTESC_REAL_TOTAREA,OCEAN_GRID_CARTESC_REAL_TOTVOL)

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
    implicit none

    ! at this moment, horizontal grid is identical to that of the atmosphere
    allocate( OCEAN_GRID_CARTESC_REAL_AREA(    OIA,OJA) )
    allocate( OCEAN_GRID_CARTESC_REAL_VOL (OKA,OIA,OJA) )
    !$acc enter data create(OCEAN_GRID_CARTESC_REAL_AREA,OCEAN_GRID_CARTESC_REAL_VOL)

    return
  end subroutine OCEAN_GRID_CARTESC_REAL_setup

  !-----------------------------------------------------------------------------
  !> Finalize
  subroutine OCEAN_GRID_CARTESC_REAL_finalize
    implicit none

    !$acc exit data delete(OCEAN_GRID_CARTESC_REAL_AREA,OCEAN_GRID_CARTESC_REAL_VOL)
    deallocate( OCEAN_GRID_CARTESC_REAL_AREA )
    deallocate( OCEAN_GRID_CARTESC_REAL_VOL )

    return
  end subroutine OCEAN_GRID_CARTESC_REAL_finalize

  !-----------------------------------------------------------------------------
  subroutine OCEAN_GRID_CARTESC_REAL_set_areavol
    use scale_atmos_grid_cartesC_real, only: &
       ATMOS_GRID_CARTESC_REAL_AREA
    use scale_ocean_grid_cartesC, only: &
       OCEAN_GRID_CARTESC_CDZ
    use scale_file_cartesC, only: &
       FILE_CARTESC_set_coordinates_ocean
    use scale_landuse, only: &
       LANDUSE_fact_ocean
    implicit none

    integer :: k, i, j

    do j = 1,   OJA
    do i = 1,   OIA
       OCEAN_GRID_CARTESC_REAL_AREA(i,j) = ATMOS_GRID_CARTESC_REAL_AREA(i,j) * LANDUSE_fact_ocean(i,j)
    end do
    end do

    OCEAN_GRID_CARTESC_REAL_TOTAREA = 0.0_RP
    do j = OJS, OJE
    do i = OIS, OIE
       OCEAN_GRID_CARTESC_REAL_TOTAREA = OCEAN_GRID_CARTESC_REAL_TOTAREA + OCEAN_GRID_CARTESC_REAL_AREA(i,j)
    end do
    end do

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

    !$acc update device(OCEAN_GRID_CARTESC_REAL_AREA,OCEAN_GRID_CARTESC_REAL_VOL)
    !$acc update device(OCEAN_GRID_CARTESC_REAL_TOTAREA,OCEAN_GRID_CARTESC_REAL_TOTVOL)
    return
  end subroutine OCEAN_GRID_CARTESC_REAL_set_areavol

end module scale_ocean_grid_cartesC_real
