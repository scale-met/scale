!-------------------------------------------------------------------------------
!> module LAND driver
!!
!! @par Description
!!          Land module driver
!!
!! @author Team SCALE
!! @li      2013-08-31 (T.Yamaura)  [new]
!!
!<
!-------------------------------------------------------------------------------
module mod_land
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
  public :: LAND_setup
  public :: LAND_step

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
  subroutine LAND_setup
    use mod_land_vars, only: &
       sw_phy => LAND_sw_phy
    use mod_land_phy_bucket, only: &
       LAND_PHY_driver_setup
    implicit none
    !---------------------------------------------------------------------------

    if ( sw_phy ) call LAND_PHY_driver_setup

    return
  end subroutine LAND_setup

  !-----------------------------------------------------------------------------
  !> Land step
  subroutine LAND_step
    use mod_land_vars, only: &
       sw_phy => LAND_sw_phy, &
       LAND_vars_history
    use mod_land_phy_bucket, only: &
       LAND_PHY_driver_first, &
       LAND_PHY_driver_final
    implicit none
    !---------------------------------------------------------------------------

    !########## Physics First ##########
    if ( sw_phy ) then
      call PROF_rapstart('LND Physics')
      call LAND_PHY_driver_first
      call PROF_rapend  ('LND Physics')
    endif

    !########## Physics Final ##########
    if ( sw_phy ) then
      call PROF_rapstart('LND Physics')
      call LAND_PHY_driver_final
      call PROF_rapend  ('LND Physics')
    endif

    !########## History & Monitor ##########
    call PROF_rapstart('LND History')
    call LAND_vars_history
    call PROF_rapend  ('LND History')

    return
  end subroutine LAND_step

end module mod_land
