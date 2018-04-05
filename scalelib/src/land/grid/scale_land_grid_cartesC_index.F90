!-------------------------------------------------------------------------------
!> module land grid index for the CartesianC grid
!!
!! @par Description
!!          Grid Index module for land
!!
!! @author Team SCALE
!<
!-------------------------------------------------------------------------------
#include "scalelib.h"
module scale_land_grid_cartesC_index
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_stdio
  use scale_prof
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: LAND_GRID_CARTESC_INDEX_setup

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  integer, public :: LKMAX = 1 ! # of computational cells: z for land
  integer, public :: LIMAX = 1 ! # of computational cells: x for land
  integer, public :: LJMAX = 1 ! # of computational cells: y for land

  integer, public :: LKA       ! # of total grids: z for land, local
  integer, public :: LIA       ! # of total grids: x for land, local
  integer, public :: LJA       ! # of total grids: y for land, local

  integer, public :: LKS       ! start point of inner domain: z for land, local
  integer, public :: LKE       ! end   point of inner domain: z for land, local
  integer, public :: LIS       ! start point of inner domain: x for land, local
  integer, public :: LIE       ! end   point of inner domain: x for land, local
  integer, public :: LJS       ! start point of inner domain: y for land, local
  integer, public :: LJE       ! end   point of inner domain: y for land, local

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
  subroutine LAND_GRID_CARTESC_INDEX_setup
    use scale_prc, only: &
       PRC_abort
    use scale_atmos_grid_cartesC_index, only: &
         IMAX, &
         IA, IS, IE, &
         JMAX, &
         JA, JS, JE
    implicit none

    namelist / PARAM_LAND_GRID_CARTESC_INDEX / &
       LKMAX

    integer :: ierr
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_PROGRESS(*) 'Module[CartesC INDEX] / Categ[LAND GRID] / Origin[SCALElib]'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_LAND_GRID_CARTESC_INDEX,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("LAND_GRID_CARTESC_INDEX_setup",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("LAND_GRID_CARTESC_INDEX_setup",*) 'Not appropriate names in namelist PARAM_LAND_GRID_CARTESC_INDEX. Check!'
       call PRC_abort
    endif
    if( IO_NML ) write(IO_FID_NML,nml=PARAM_LAND_GRID_CARTESC_INDEX)

    LKS  = 1
    LKE  = LKMAX
    LKA  = LKMAX

    LOG_NEWLINE
    LOG_INFO("LAND_GRID_CARTESC_INDEX_setup",*) 'Land grid index information '
    LOG_INFO_CONT('(1x,A,I6,A,I6,A,I6)') 'z-axis levels :', LKMAX

    ! at this moment horizontal grid is same as that in atmosphere
    LIMAX = IMAX
    LIA = IA
    LIS = IS
    LIE = IE

    LJMAX = JMAX
    LJA = JA
    LJS = JS
    LJE = JE

    return
  end subroutine LAND_GRID_CARTESC_INDEX_setup

end module scale_land_grid_cartesC_index
