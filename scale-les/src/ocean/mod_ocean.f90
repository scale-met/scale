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
  use mod_stdio, only: &
     IO_FID_LOG,  &
     IO_L
  use mod_time, only: &
     TIME_rapstart, &
     TIME_rapend
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
       sw_phy => OCEAN_sw_phy,  &
       OCEAN_vars_setup,        &
       OCEAN_vars_restart_read
    use mod_ocean_phy, only: &
       OCEAN_PHY_setup
    implicit none
    !---------------------------------------------------------------------------

    call OCEAN_vars_setup

    call OCEAN_vars_restart_read

    if ( sw_phy ) call OCEAN_PHY_setup

    return
  end subroutine OCEAN_setup

  !-----------------------------------------------------------------------------
  !> Ocean step
  subroutine OCEAN_step
    use mod_ocean_vars, only: &
       sw_phy => OCEAN_sw_phy, &
       OCEAN_vars_history
    use mod_ocean_phy, only: &
       OCEAN_PHY
    implicit none
    !---------------------------------------------------------------------------

    !########## Physics ##########
    call TIME_rapstart('OCN Physics')
    if ( sw_phy ) then
       call OCEAN_PHY
    endif
    call TIME_rapend  ('OCN Physics')

    !########## History & Monitor ##########
    call TIME_rapstart('OCN History')
       call OCEAN_vars_history
    call TIME_rapend  ('OCN History')

    return
  end subroutine OCEAN_step

end module mod_ocean
