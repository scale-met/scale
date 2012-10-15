!-------------------------------------------------------------------------------
!> module Atmosphere / Physics Radiation
!!
!! @par Description
!!          Atmospheric radiation transfer process
!!          mstrnX
!!
!! @author H.Tomita and SCALE developpers
!!
!! @par History
!! @li      2012-03-26 (H.Yashiro) [new]
!!
!<
!-------------------------------------------------------------------------------
module mod_atmos_phy_rd
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use mod_stdio, only: &
     IO_FID_LOG,  &
     IO_L
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: ATMOS_PHY_RD_setup
  public :: ATMOS_PHY_RD

  !-----------------------------------------------------------------------------
  !
  !++ included parameters
  !
  include "inc_precision.h"
  include 'inc_index.h'
  include 'inc_tracer.h'

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
  subroutine ATMOS_PHY_RD_setup
    implicit none
    !---------------------------------------------------------------------------

    return
  end subroutine ATMOS_PHY_RD_setup

  !-----------------------------------------------------------------------------
  subroutine ATMOS_PHY_RD
    implicit none
    !---------------------------------------------------------------------------

    return
  end subroutine ATMOS_PHY_RD

end module mod_atmos_phy_rd
