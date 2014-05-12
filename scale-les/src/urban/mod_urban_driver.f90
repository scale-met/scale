!-------------------------------------------------------------------------------
!> module URBAN driver
!!
!! @par Description
!!          Urban module driver
!!
!! @author Team SCALE
!!
!<
!-------------------------------------------------------------------------------
module mod_urban_driver
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_stdio
  use scale_prof
  use scale_grid_index
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: URBAN_driver_setup
  public :: URBAN_driver

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
  subroutine URBAN_driver_setup
    use mod_urban_vars, only: &
       sw_phy => URBAN_sw_phy
    use mod_urban_phy_ucm, only: &
       URBAN_PHY_driver_setup
    implicit none
    !---------------------------------------------------------------------------

    if ( sw_phy ) call URBAN_PHY_driver_setup

    return
  end subroutine URBAN_driver_setup

  !-----------------------------------------------------------------------------
  !> Urban step
  subroutine URBAN_driver
    use mod_urban_vars, only: &
       sw_phy => URBAN_sw_phy, &
       URBAN_vars_history
    use mod_urban_phy_ucm, only: &
       URBAN_PHY_driver_first, &
       URBAN_PHY_driver_final
    implicit none
    !---------------------------------------------------------------------------

    !########## Physics First ##########
    if ( sw_phy ) then
      call PROF_rapstart('URB Physics')
      call URBAN_PHY_driver_first
      call PROF_rapend  ('URB Physics')
    endif

    !########## Physics Final ##########
    if ( sw_phy ) then
      call PROF_rapstart('URB Physics')
      call URBAN_PHY_driver_final
      call PROF_rapend  ('URB Physics')
    endif

    !########## History & Monitor ##########
    call PROF_rapstart('URB History')
    call URBAN_vars_history
    call PROF_rapend  ('URB History')

    return
  end subroutine URBAN_driver

end module mod_urban_driver
