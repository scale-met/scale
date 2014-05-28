!-------------------------------------------------------------------------------
!> module LAND driver
!!
!! @par Description
!!          Land model driver
!!
!! @author Team SCALE
!! @li      2013-08-31 (T.Yamaura)  [new]
!<
!-------------------------------------------------------------------------------
module mod_land_driver
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
  public :: LAND_driver_setup
  public :: LAND_driver

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
  subroutine LAND_driver_setup
    use mod_land_admin, only: &
       LAND_TYPE, &
       LAND_sw
    use mod_land_phy_bucket, only: &
       LAND_PHY_driver_setup
    implicit none
    !---------------------------------------------------------------------------

    if( LAND_sw ) call LAND_PHY_driver_setup( LAND_TYPE )

    return
  end subroutine LAND_driver_setup

  !-----------------------------------------------------------------------------
  !> Land step
  subroutine LAND_driver
    use mod_land_admin, only: &
       LAND_sw
    use mod_land_vars, only: &
       LAND_vars_history
    use mod_land_phy_bucket, only: &
       LAND_PHY_driver_first, &
       LAND_PHY_driver_final
    implicit none
    !---------------------------------------------------------------------------

    !########## Physics First ##########
    if ( LAND_sw ) then
      call PROF_rapstart('LND Physics')
      call LAND_PHY_driver_first
      call PROF_rapend  ('LND Physics')
    endif

    !########## Physics Final ##########
    if ( LAND_sw ) then
      call PROF_rapstart('LND Physics')
      call LAND_PHY_driver_final
      call PROF_rapend  ('LND Physics')
    endif

    !########## History & Monitor ##########
    call PROF_rapstart('LND History')
    call LAND_vars_history
    call PROF_rapend  ('LND History')

    return
  end subroutine LAND_driver

end module mod_land_driver
