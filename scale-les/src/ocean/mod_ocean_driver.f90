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
  use scale_debug
  use scale_grid_index
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: OCEAN_driver_setup
  public :: OCEAN_driver
  public :: OCEAN_SURFACE_SET

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
    use mod_ocean_phy_slab, only: &
       OCEAN_PHY_driver_setup
!    use mod_ocean_frc_nudge, only: &
!       OCEAN_FRC_driver_setup
    implicit none
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '++++++ Module[DRIVER] / Categ[OCEAN] / Origin[SCALE-LES]'

    call OCEAN_PHY_driver_setup

!    if( OCEAN_FRC_sw ) call OCEAN_FRC_driver_setup

    return
  end subroutine OCEAN_driver_setup

  !-----------------------------------------------------------------------------
  !> Ocean step
  subroutine OCEAN_driver
    use scale_time, only: &
       dt => TIME_DTSEC_OCEAN
    use mod_ocean_admin, only: &
       OCEAN_sw
    use mod_ocean_vars, only: &
       OCEAN_TEMP,        &
       OCEAN_TEMP_t,      &
       OCEAN_vars_total,  &
       OCEAN_vars_history
    use mod_ocean_phy_slab, only: &
       OCEAN_PHY_driver
!    use mod_ocean_forcing, only: &
!       OCEAN_forcing
    implicit none

    integer :: i, j
    !---------------------------------------------------------------------------

    !########## Physics ##########
    if ( OCEAN_sw ) then
       call PROF_rapstart('OCN Physics')
       call OCEAN_PHY_driver( update_flag=.true., history_flag=.true. )
       call PROF_rapend  ('OCN Physics')
    endif

    !########## Forcing ##########
!    if ( OCEAN_FORCE_sw ) then
!       call PROF_rapstart('OCN Forcing')
!       call OCEAN_forcing
!       call PROF_rapend  ('OCN Forcing')
!    endif

    !########## Update ##########
    do j = JS, JE
    do i = IS, IE
       OCEAN_TEMP(i,j) = OCEAN_TEMP(i,j) + OCEAN_TEMP_t(i,j) * dt
    enddo
    enddo

    call OCEAN_vars_total

    !########## Put Surface Boundary to coupler ##########
    call OCEAN_SURFACE_SET( setup=.false. )

    !########## History & Monitor ##########
    call PROF_rapstart('OCN History')
    call OCEAN_vars_history
    call PROF_rapend  ('OCN History')

    do j = JS, JE
    do i = IS, IE
       OCEAN_TEMP_t(i,j) = 0.0_RP
    enddo
    enddo

    return
  end subroutine OCEAN_driver

  !-----------------------------------------------------------------------------
  !> Put surface boundary to other model
  subroutine OCEAN_SURFACE_SET( setup )
    use mod_ocean_vars, only: &
       OCEAN_vars_fillhalo, &
       OCEAN_TEMP,          &
       OCEAN_SFC_TEMP,      &
       OCEAN_SFC_albedo,    &
       OCEAN_SFC_Z0
    use mod_cpl_admin, only: &
       CPL_sw_AtmOcn
    use mod_cpl_vars, only: &
       CPL_putOcn_setup, &
       CPL_putOcn
    implicit none

    logical, intent(in) :: setup
    !---------------------------------------------------------------------------

    if ( CPL_sw_AtmOcn ) then
       call OCEAN_vars_fillhalo

       if ( setup ) then
          call CPL_putOcn_setup( OCEAN_SFC_TEMP  (:,:),   & ! [IN]
                                 OCEAN_SFC_albedo(:,:,:), & ! [IN]
                                 OCEAN_SFC_Z0    (:,:)    ) ! [IN]
       endif

       call CPL_putOcn( OCEAN_TEMP(:,:) ) ! [IN]
    endif

    return
  end subroutine OCEAN_SURFACE_SET

end module mod_ocean_driver
