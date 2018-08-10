!-------------------------------------------------------------------------------
!> module ATMOSPHERE driver
!!
!! @par Description
!!          Atmosphere module
!!
!! @author Team SCALE
!!
!<
!-------------------------------------------------------------------------------
module mod_atmos_surface
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_io
  use scale_atmos_grid_icoA_index

  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: ATMOS_SURFACE_GET
  public :: ATMOS_SURFACE_SET

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
  !> Get surface boundary condition
  subroutine ATMOS_SURFACE_GET

    return
  end subroutine ATMOS_SURFACE_GET

  !-----------------------------------------------------------------------------
  !> Set surface boundary condition
  subroutine ATMOS_SURFACE_SET

    return
  end subroutine ATMOS_SURFACE_SET

end module mod_atmos_surface
