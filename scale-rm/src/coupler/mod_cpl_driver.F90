!-------------------------------------------------------------------------------
!> module CPL driver
!!
!! @par Description
!!          Coupler driver
!!
!! @author Team SCALE
!<
!-------------------------------------------------------------------------------
#include "scalelib.h"
module mod_cpl_driver
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
  public :: CPL_driver_setup
  public :: CPL_driver

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
  !> Setup
  subroutine CPL_driver_setup
    implicit none
    !---------------------------------------------------------------------------

    return
  end subroutine CPL_driver_setup

  !-----------------------------------------------------------------------------
  !> CPL calcuration
  subroutine CPL_driver
    implicit none
    !---------------------------------------------------------------------------

    return
  end subroutine CPL_driver

end module mod_cpl_driver
