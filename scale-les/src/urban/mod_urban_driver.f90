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
    use mod_urban_admin, only: &
       URBAN_TYPE, &
       URBAN_sw
    use mod_urban_phy_ucm, only: &
       URBAN_PHY_driver_setup
    use mod_urban_phy_ucm, only: &
       URBAN_SURFACE_SET
    implicit none
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '++++++ Module[DRIVER] / Categ[URBAN] / Origin[SCALE-LES]'

    if( URBAN_sw ) call URBAN_PHY_driver_setup( URBAN_TYPE )

    if( URBAN_sw ) call URBAN_SURFACE_SET

    return
  end subroutine URBAN_driver_setup

  !-----------------------------------------------------------------------------
  !> Urban step
  subroutine URBAN_driver
    use mod_urban_admin, only: &
       URBAN_sw
    use mod_urban_vars, only: &
       URBAN_vars_history
    use mod_urban_phy_ucm, only: &
       URBAN_PHY_driver, &
       URBAN_SURFACE_SET
    implicit none
    !---------------------------------------------------------------------------

    !########## Physics First ##########
    if ( URBAN_sw ) then
      call PROF_rapstart('URB Physics')
      call URBAN_PHY_driver
      call PROF_rapend  ('URB Physics')
    endif

    !########## Physics Final ##########
    if ( URBAN_sw ) then
      call PROF_rapstart('URB Physics')
      call URBAN_SURFACE_SET
      call PROF_rapend  ('URB Physics')
    endif

    !########## History & Monitor ##########
    call PROF_rapstart('URB History')
    call URBAN_vars_history
    call PROF_rapend  ('URB History')

    return
  end subroutine URBAN_driver

end module mod_urban_driver
