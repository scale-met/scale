!-------------------------------------------------------------------------------
!> module OCEAN driver
!!
!! @par Description
!!          Ocean model driver
!!
!! @author Team SCALE
!<
!-------------------------------------------------------------------------------
module mod_ocean_driver
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
  public :: OCEAN_driver_setup
  public :: OCEAN_driver

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
  subroutine OCEAN_driver_setup
    use mod_ocean_admin, only: &
       OCEAN_TYPE, &
       OCEAN_sw
    use mod_ocean_phy_slab, only: &
       OCEAN_PHY_driver_setup
    use mod_ocean_phy_slab, only: &
       OCEAN_PHY_driver_final
    implicit none
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '++++++ Module[DRIVER] / Categ[OCEAN] / Origin[SCALE-LES]'

    if( OCEAN_sw ) call OCEAN_PHY_driver_setup( OCEAN_TYPE )

    if( OCEAN_sw ) call OCEAN_PHY_driver_final

    return
  end subroutine OCEAN_driver_setup

  !-----------------------------------------------------------------------------
  !> Ocean step
  subroutine OCEAN_driver
    use mod_ocean_admin, only: &
       OCEAN_sw
    use mod_ocean_vars, only: &
       OCEAN_vars_history
    use mod_ocean_phy_slab, only: &
       OCEAN_PHY_driver_first, &
       OCEAN_PHY_driver_final
    implicit none
    !---------------------------------------------------------------------------

    !########## Physics First ##########
    if ( OCEAN_sw ) then
      call PROF_rapstart('OCN Physics')
      call OCEAN_PHY_driver_first
      call PROF_rapend  ('OCN Physics')
    endif

    !########## Physics Final ##########
    if ( OCEAN_sw ) then
      call PROF_rapstart('OCN Physics')
      call OCEAN_PHY_driver_final
      call PROF_rapend  ('OCN Physics')
    endif

    !########## History & Monitor ##########
    call PROF_rapstart('OCN History')
    call OCEAN_vars_history
    call PROF_rapend  ('OCN History')

    return
  end subroutine OCEAN_driver

end module mod_ocean_driver
