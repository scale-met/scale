!-------------------------------------------------------------------------------
!> module LAND driver
!!
!! @par Description
!!          Land module driver
!!
!! @author Team SCALE
!!
!<
!-------------------------------------------------------------------------------
module mod_land
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
       sw_phy => LAND_sw_phy,  &
       LAND_vars_setup,        &
       LAND_vars_restart_read
    use mod_land_phy, only: &
       LAND_PHY_setup
    implicit none
    !---------------------------------------------------------------------------

    call LAND_vars_setup

    call LAND_vars_restart_read

    if ( sw_phy ) call LAND_PHY_setup

    return
  end subroutine LAND_setup

  !-----------------------------------------------------------------------------
  !> Land step
  subroutine LAND_step
    use mod_land_vars, only: &
       sw_phy => LAND_sw_phy, &
       LAND_vars_history
    use mod_land_phy, only: &
       LAND_PHY
    implicit none
    !---------------------------------------------------------------------------

    !########## Physics ##########
    call TIME_rapstart('LND Physics')
    if ( sw_phy ) then
       call LAND_PHY
    endif
    call TIME_rapend  ('LND Physics')

    !########## History & Monitor ##########
    call TIME_rapstart('LND History')
       call LAND_vars_history
    call TIME_rapend  ('LND History')

    return
  end subroutine LAND_step

end module mod_land
