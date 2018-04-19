!-------------------------------------------------------------------------------
!> module ocean grid cartesc index
!!
!! @par Description
!!          Grid Index module for ocean
!!
!! @author Team SCALE
!<
!-------------------------------------------------------------------------------
#include "scalelib.h"
module scale_ocean_grid_cartesC_index
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
  public :: OCEAN_GRID_CARTESC_INDEX_setup

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  integer, public :: OKMAX = -1 ! # of computational cells: z for ocean
  integer, public :: OIMAX = -1 ! # of computational cells: x for ocean
  integer, public :: OJMAX = -1 ! # of computational cells: y for ocean

  integer, public :: OKA       ! # of total grids: z for ocean, local
  integer, public :: OKS       ! start point of inner domain: z for ocean, local
  integer, public :: OKE       ! end   point of inner domain: z for ocean, local

  integer, public :: OIA       ! # of total grids: x for ocean, local
  integer, public :: OIS       ! start point of inner domain: x for ocean, local
  integer, public :: OIE       ! end   point of inner domain: x for ocean, local

  integer, public :: OJA       ! # of total grids: y for ocean, local
  integer, public :: OJS       ! start point of inner domain: y for ocean, local
  integer, public :: OJE       ! end   point of inner domain: y for ocean, local

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
  subroutine OCEAN_GRID_CARTESC_INDEX_setup
    use scale_prc, only: &
         PRC_abort
    use scale_atmos_grid_cartesC_index, only: &
         IMAX, &
         IA, IS, IE, &
         JMAX, &
         JA, JS, JE
    implicit none

    namelist / PARAM_OCEAN_GRID_CARTESC_INDEX / &
       OKMAX

    integer :: ierr
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("OCEAN_GRID_CARTESC_INDEX_setup",*) 'Setup'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_OCEAN_GRID_CARTESC_INDEX,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("OCEAN_GRID_CARTESC_INDEX_setup",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("OCEAN_GRID_CARTESC_INDEX_setup",*) 'Not appropriate names in namelist PARAM_OCEAN_GRID_CARTESC_INDEX. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_OCEAN_GRID_CARTESC_INDEX)

    OKA  = OKMAX
    OKS  = 1
    OKE  = OKMAX

    LOG_NEWLINE
    LOG_INFO("OCEAN_GRID_CARTESC_INDEX_setup",*) 'Ocean grid index information '
    LOG_INFO_CONT('(1x,A,I6,A,I6,A,I6)') 'z-axis levels :', OKMAX


    ! at this moment horizontal grid is same as that in atmosphere
    OIMAX = IMAX
    OIA = IA
    OIS = IS
    OIE = IE

    OJMAX = JMAX
    OJA = JA
    OJS = JS
    OJE = JE

    return
  end subroutine OCEAN_GRID_CARTESC_INDEX_setup

end module scale_ocean_grid_cartesC_index
