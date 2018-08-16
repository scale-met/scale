!-------------------------------------------------------------------------------
!> module urban / grid / icosahedralA / index
!!
!! @par Description
!!          Grid Index module for urban
!!
!! @author Team SCALE
!<
!-------------------------------------------------------------------------------
#include "scalelib.h"
module scale_urban_grid_cartesC_index
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
  public :: URBAN_GRID_CARTESC_INDEX_setup

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  integer, public :: UKMAX = -1 ! # of computational cells: z for urban
  integer, public :: UIMAX = -1 ! # of computational cells: x for urban
  integer, public :: UJMAX = -1 ! # of computational cells: y for urban

  integer, public :: UKA       ! # of total grids: z for urban, local
  integer, public :: UKS       ! start point of inner domain: z for urban, local
  integer, public :: UKE       ! end   point of inner domain: z for urban, local
  integer, public :: UIA       ! # of total grids: x for urban, local
  integer, public :: UIS       ! start point of inner domain: x for urban, local
  integer, public :: UIE       ! end   point of inner domain: x for urban, local
  integer, public :: UJA       ! # of total grids: Y for urban, local
  integer, public :: UJS       ! start point of inner domain: y for urban, local
  integer, public :: UJE       ! end   point of inner domain: y for urban, local

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
  subroutine URBAN_GRID_CARTESC_INDEX_setup
    use scale_prc, only: &
       PRC_abort
    use scale_atmos_grid_cartesC_index, only: &
       IMAX,       &
       IA, IS, IE, &
       JMAX,       &
       JA, JS, JE
    implicit none

    namelist / PARAM_URBAN_GRID_CARTESC_INDEX / &
       UKMAX

    integer :: ierr
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("URBAN_GRID_CARTESC_INDEX_setup",*) 'Setup'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_URBAN_GRID_CARTESC_INDEX,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("URBAN_GRID_CARTESC_INDEX_setup",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("URBAN_GRID_CARTESC_INDEX_setup",*) 'Not appropriate names in namelist PARAM_URBAN_GRID_CARTESC_INDEX. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_URBAN_GRID_CARTESC_INDEX)

    if ( UKMAX < 1 ) then
       LOG_ERROR("URBAN_GRID_CARTESC_INDEX_setup",*) 'UKMAX must be >= 1 ', UKMAX
       call PRC_abort
    end if

    UKA  = UKMAX
    UKS  = 1
    UKE  = UKMAX

    LOG_NEWLINE
    LOG_INFO("URBAN_GRID_CARTESC_INDEX_setup",*) 'Urban grid index information '
    LOG_INFO_CONT('(1x,A,I6,A,I6,A,I6)') 'z-axis levels :', UKMAX

    ! at this moment horizontal grid is same as that in atmosphere
    UIMAX = IMAX
    UIA = IA
    UIS = IS
    UIE = IE

    UJMAX = JMAX
    UJA = JA
    UJS = JS
    UJE = JE

    return
  end subroutine URBAN_GRID_CARTESC_INDEX_setup

end module scale_urban_grid_cartesC_index
