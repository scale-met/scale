!-------------------------------------------------------------------------------
!> module urban / grid / cartesianC / real
!!
!! @par Description
!!          Grid module for cartesian coordinate for urban
!!
!! @author Team SCALE
!!
!<
!-------------------------------------------------------------------------------
#include "scalelib.h"
module scale_urban_grid_cartesC_real
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_io
  use scale_prof
  use scale_atmos_grid_cartesC_index
  use scale_urban_grid_cartesC_index
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: URBAN_GRID_CARTESC_REAL_setup
  public :: URBAN_GRID_CARTESC_REAL_finalize
  public :: URBAN_GRID_CARTESC_REAL_set_areavol

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  real(RP), public, allocatable :: URBAN_GRID_CARTESC_REAL_AREA(:,:)   !< area   of grid cell
  real(RP), public              :: URBAN_GRID_CARTESC_REAL_TOTAREA     !< total area
  real(RP), public, allocatable :: URBAN_GRID_CARTESC_REAL_VOL (:,:,:) !< volume of grid cell
  real(RP), public              :: URBAN_GRID_CARTESC_REAL_TOTVOL      !< total volume

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
  subroutine URBAN_GRID_CARTESC_REAL_setup

    ! at this moment, horizontal grid is identical to that of the atmosphere
    allocate( URBAN_GRID_CARTESC_REAL_AREA(    UIA,UJA) )
    allocate( URBAN_GRID_CARTESC_REAL_VOL (UKA,UIA,UJA) )
    !$acc enter data create(URBAN_GRID_CARTESC_REAL_AREA,URBAN_GRID_CARTESC_REAL_VOL)

    return
  end subroutine URBAN_GRID_CARTESC_REAL_setup

  !-----------------------------------------------------------------------------
  !> Finalize
  subroutine URBAN_GRID_CARTESC_REAL_finalize
    implicit none

    !$acc exit data delete(URBAN_GRID_CARTESC_REAL_AREA,URBAN_GRID_CARTESC_REAL_VOL)
    deallocate( URBAN_GRID_CARTESC_REAL_AREA )
    deallocate( URBAN_GRID_CARTESC_REAL_VOL  )

    return
  end subroutine URBAN_GRID_CARTESC_REAL_finalize

  !-----------------------------------------------------------------------------
  subroutine URBAN_GRID_CARTESC_REAL_set_areavol
    use scale_atmos_grid_cartesC_real, only: &
       ATMOS_GRID_CARTESC_REAL_AREA
    use scale_urban_grid_cartesC, only: &
       URBAN_GRID_CARTESC_CDZ
    use scale_file_cartesC, only: &
       FILE_CARTESC_set_coordinates_urban
    use scale_landuse, only: &
       LANDUSE_fact_urban

    integer :: k, i, j

    do j = 1, UJA
    do i = 1, UIA
       URBAN_GRID_CARTESC_REAL_AREA(i,j) = ATMOS_GRID_CARTESC_REAL_AREA(i,j) * LANDUSE_fact_urban(i,j)
    end do
    end do

    URBAN_GRID_CARTESC_REAL_TOTAREA = 0.0_RP
    do j = UJS, UJE
    do i = UIS, UIE
       URBAN_GRID_CARTESC_REAL_TOTAREA = URBAN_GRID_CARTESC_REAL_TOTAREA + URBAN_GRID_CARTESC_REAL_AREA(i,j)
    end do
    end do

    do j = 1,   UJA
    do i = 1,   UIA
    do k = UKS, UKE
       URBAN_GRID_CARTESC_REAL_VOL(k,i,j) = URBAN_GRID_CARTESC_REAL_AREA(i,j) * URBAN_GRID_CARTESC_CDZ(k)
    enddo
    enddo
    enddo

    URBAN_GRID_CARTESC_REAL_TOTVOL = 0.0_RP
    do j = UJS, UJE
    do i = UIS, UIE
    do k = UKS, UKE
       URBAN_GRID_CARTESC_REAL_TOTVOL = URBAN_GRID_CARTESC_REAL_TOTVOL + URBAN_GRID_CARTESC_REAL_VOL(k,i,j)
    end do
    end do
    end do

    !$acc update device(URBAN_GRID_CARTESC_REAL_AREA,URBAN_GRID_CARTESC_REAL_VOL)

    call FILE_CARTESC_set_coordinates_urban( URBAN_GRID_CARTESC_REAL_VOL(:,:,:) ) ! [IN]

    return
  end subroutine URBAN_GRID_CARTESC_REAL_set_areavol

end module scale_urban_grid_cartesC_real
