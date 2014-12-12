!-------------------------------------------------------------------------------
!> module URBAN driver
!!
!! @par Description
!!          Urban module driver
!!
!! @author Team SCALE
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
  use scale_debug
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
    use mod_urban_phy_ucm, only: &
       URBAN_PHY_driver_setup
    use mod_urban_vars, only: &
       URBAN_vars_history
    use mod_urban_admin, only: &
       URBAN_sw
    implicit none
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '++++++ Module[DRIVER] / Categ[URBAN] / Origin[SCALE-LES]'

    !########## Get Surface Boundary from coupler ##########
    call URBAN_SURFACE_GET( setup=.true. )

    call URBAN_PHY_driver_setup

    !########## History & Monitor ##########
    if ( URBAN_sw ) then
       call PROF_rapstart('URB History', 1)
       call URBAN_vars_history
       call PROF_rapend  ('URB History', 1)
    endif

    return
  end subroutine URBAN_driver_setup

  !-----------------------------------------------------------------------------
  !> Urban step
  subroutine URBAN_driver
    use scale_time, only: &
       dt => TIME_DTSEC_URBAN
    use mod_urban_admin, only: &
       URBAN_sw
    use mod_urban_vars, only: &
       URBAN_TEMP,        &
       URBAN_TEMP_t,      &
       URBAN_vars_total,  &
       URBAN_vars_history
    use mod_urban_phy_ucm, only: &
       URBAN_PHY_driver
    implicit none

    integer :: i, j
    !---------------------------------------------------------------------------

    !########## Get Surface Boundary from coupler ##########
    call URBAN_SURFACE_GET( setup=.false. )

    !########## Physics ##########
    if ( URBAN_sw ) then
       call PROF_rapstart('URB Physics', 1)
       call URBAN_PHY_driver( update_flag = .true. )
       call PROF_rapend  ('URB Physics', 1)
    endif

    !########## Update ##########
!    do j = JS, JE
!    do i = IS, IE
!       URBAN_TEMP(i,j) = URBAN_TEMP(i,j) + URBAN_TEMP_t(i,j) * dt
!    enddo
!    enddo

    call URBAN_vars_total

    !########## Put Surface Boundary to coupler ##########
    call URBAN_SURFACE_SET( setup=.false. )

    !########## History & Monitor ##########
    call PROF_rapstart('URB History', 1)
    call URBAN_vars_history
    call PROF_rapend  ('URB History', 1)

    do j = JS, JE
    do i = IS, IE
       URBAN_TEMP_t(i,j) = 0.0_RP
    enddo
    enddo

    return
  end subroutine URBAN_driver

  !-----------------------------------------------------------------------------
  !> Get surface boundary
  subroutine URBAN_SURFACE_GET( setup )
    use mod_cpl_admin, only: &
       CPL_sw_AtmUrb
    use mod_cpl_vars, only: &
       CPL_getUrb
    implicit none

    logical, intent(in) :: setup
    !---------------------------------------------------------------------------

    return
  end subroutine URBAN_SURFACE_GET

  !-----------------------------------------------------------------------------
  !> Put surface boundary to other model
  subroutine URBAN_SURFACE_SET( setup )
    use mod_cpl_admin, only: &
       CPL_sw_AtmUrb
    use mod_cpl_vars, only: &
       CPL_putUrb_setup, &
       CPL_putUrb_restart, &
       CPL_putUrb
    implicit none

    logical, intent(in) :: setup
    !---------------------------------------------------------------------------

    return
  end subroutine URBAN_SURFACE_SET

end module mod_urban_driver
