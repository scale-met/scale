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
  public :: URBAN_SURFACE_SET

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
    implicit none
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '++++++ Module[DRIVER] / Categ[URBAN] / Origin[SCALE-LES]'

    if( URBAN_sw ) call URBAN_PHY_driver_setup( URBAN_TYPE )

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
       URBAN_PHY_driver
    implicit none
    !---------------------------------------------------------------------------

    !########## Physics ##########
    if ( URBAN_sw ) then
      call PROF_rapstart('URB Physics')
      call URBAN_PHY_driver
      call PROF_rapend  ('URB Physics')
    endif

    !########## Put surface boundary to other model ##########
    call URBAN_SURFACE_SET( setup=.false. )

    !########## History & Monitor ##########
    call PROF_rapstart('URB History')
    call URBAN_vars_history
    call PROF_rapend  ('URB History')

    return
  end subroutine URBAN_driver

  !-----------------------------------------------------------------------------
  !> Put surface boundary to other model
  subroutine URBAN_SURFACE_SET( setup )
    use mod_urban_vars, only: &
       URBAN_vars_fillhalo
    use mod_cpl_admin, only: &
       CPL_sw_AtmUrb
    use mod_cpl_vars, only: &
       CPL_putUrb_setup, &
       CPL_putUrb
    implicit none

    logical, intent(in) :: setup
    !---------------------------------------------------------------------------

    if ( CPL_sw_AtmUrb ) then
       call URBAN_vars_fillhalo

       if ( setup ) then
!       call CPL_putUrb_setup( URBAN_TEMP(:,:) ) ! [IN]
       endif

!       call CPL_putUrb( URBAN_TEMP(:,:) ) ! [IN]

    endif

    return
  end subroutine URBAN_SURFACE_SET

end module mod_urban_driver
