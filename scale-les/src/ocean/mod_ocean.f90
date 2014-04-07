!-------------------------------------------------------------------------------
!> module OCEAN driver
!!
!! @par Description
!!          Ocean module driver
!!
!! @author Team SCALE
!!
!<
!-------------------------------------------------------------------------------
module mod_ocean
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
  public :: OCEAN_setup
  public :: OCEAN_step

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
  subroutine OCEAN_setup
    use mod_ocean_vars, only: &
       sw_phy => OCEAN_sw_phy
    use mod_ocean_phy_slab, only: &
       OCEAN_PHY_driver_setup
    implicit none
    !---------------------------------------------------------------------------

    if ( sw_phy ) call OCEAN_PHY_driver_setup

    return
  end subroutine OCEAN_setup

  !-----------------------------------------------------------------------------
  !> Ocean step
  subroutine OCEAN_step
    use mod_ocean_vars, only: &
       sw_phy => OCEAN_sw_phy, &
       OCEAN_vars_history
    use mod_ocean_phy_slab, only: &
       OCEAN_PHY_driver_first, &
       OCEAN_PHY_driver_final
    implicit none
    !---------------------------------------------------------------------------

    !########## Physics First ##########
    if ( sw_phy ) then
      call PROF_rapstart('OCN Physics')
      call OCEAN_PHY_driver_first
      call PROF_rapend  ('OCN Physics')
    endif

    !########## Physics Final ##########
    if ( sw_phy ) then
      call PROF_rapstart('OCN Physics')
      call OCEAN_PHY_driver_final
      call PROF_rapend  ('OCN Physics')
    endif

    !########## History & Monitor ##########
    call PROF_rapstart('OCN History')
    call OCEAN_vars_history
    call PROF_rapend  ('OCN History')

    return
  end subroutine OCEAN_step

end module mod_ocean
