!-------------------------------------------------------------------------------
!> module CPL driver
!!
!! @par Description
!!          Coupler driver
!!
!! @author Team SCALE
!! @li      2013-08-31 (T.Yamaura)  [new]
!<
!-------------------------------------------------------------------------------
module mod_cpl_driver
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_stdio
  use scale_prof
  use scale_debug
  use scale_grid_index
  use scale_tracer
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: CPL_driver_setup
  public :: CPL_driver

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
  subroutine CPL_driver_setup
    use mod_cpl_vars, only: &
       CPL_vars_merge
    use mod_cpl_atmos_ocean_driver, only: &
       CPL_AtmOcn_driver_setup
    use mod_cpl_atmos_land_driver, only: &
       CPL_AtmLnd_driver_setup
    use mod_cpl_atmos_urban_driver, only: &
       CPL_AtmUrb_driver_setup
    implicit none
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '++++++ Module[DRIVER] / Categ[CPL] / Origin[SCALE-LES]'

    call CPL_AtmOcn_driver_setup
    call CPL_AtmLnd_driver_setup
    call CPL_AtmUrb_driver_setup

    call CPL_vars_merge

    return
  end subroutine CPL_driver_setup

  !-----------------------------------------------------------------------------
  !> CPL calcuration
  subroutine CPL_driver
    use mod_cpl_admin, only: &
       CPL_sw_AtmOcn,              &
       CPL_sw_AtmLnd,              &
       CPL_sw_AtmUrb,              &
       CPL_OCN_SFC_TEMP_UPDATE, &
       CPL_LND_SFC_TEMP_UPDATE, &
       CPL_URB_SFC_TEMP_UPDATE
    use mod_cpl_vars, only: &
       CPL_vars_merge
    use mod_cpl_atmos_ocean_driver, only: &
       CPL_AtmOcn_driver
    use mod_cpl_atmos_land_driver, only: &
       CPL_AtmLnd_driver
    use mod_cpl_atmos_urban_driver, only: &
       CPL_AtmUrb_driver
    implicit none
    !---------------------------------------------------------------------------

    !########## Coupler Atoms-Ocean ##########
    if( CPL_sw_AtmOcn ) then
      call PROF_rapstart('CPL Atmos-Ocean', 1)
      call CPL_AtmOcn_driver( CPL_OCN_SFC_TEMP_UPDATE )
      call PROF_rapend  ('CPL Atmos-Ocean', 1)
    endif

    !########## Coupler Atoms-Land ##########
    if( CPL_sw_AtmLnd ) then
      call PROF_rapstart('CPL Atmos-Land', 1)
      call CPL_AtmLnd_driver( CPL_LND_SFC_TEMP_UPDATE )
      call PROF_rapend  ('CPL Atmos-Land', 1)
    endif

    !########## Coupler Atoms-Urban ##########
    if( CPL_sw_AtmUrb ) then
      call PROF_rapstart('CPL Atmos-Urban', 1)
      call CPL_AtmUrb_driver( CPL_URB_SFC_TEMP_UPDATE )
      call PROF_rapend  ('CPL Atmos-Urban', 1)
    endif

    !########## Merge Land/Ocean/Urban with fraction ##########
    call CPL_vars_merge

    return
  end subroutine CPL_driver

end module mod_cpl_driver
