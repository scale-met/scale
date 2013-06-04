!-------------------------------------------------------------------------------
!> module ATMOSPHERE driver
!!
!! @par Description
!!          Atmosphere module
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2011-11-11 (H.Yashiro) [new]
!! @li      2012-02-15 (H.Yashiro) [add] Microphysics
!! @li      2012-03-23 (H.Yashiro) [mod] Cleaning
!!
!<
!-------------------------------------------------------------------------------
module mod_atmos
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
  public :: ATMOS_setup
  public :: ATMOS_step

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
  !> Setup atmosphere
  !-----------------------------------------------------------------------------
  subroutine ATMOS_setup
    use mod_atmos_thermodyn, only: &
       ATMOS_THERMODYN_setup
    use mod_atmos_saturation, only: &
       ATMOS_SATURATION_setup
    use mod_atmos_vars, only: &
       sw_dyn    => ATMOS_sw_dyn,    &
       sw_phy_sf => ATMOS_sw_phy_sf, &
       sw_phy_tb => ATMOS_sw_phy_tb, &
       sw_phy_mp => ATMOS_sw_phy_mp, &
       sw_phy_rd => ATMOS_sw_phy_rd, &
       ATMOS_vars_setup,  &
       ATMOS_vars_restart_read
    use mod_atmos_vars_sf, only: &
       ATMOS_vars_sf_setup
    use mod_atmos_refstate, only: &
       ATMOS_REFSTATE_setup
    use mod_atmos_boundary, only: &
       ATMOS_BOUNDARY_setup
    use mod_atmos_dyn, only: &
       ATMOS_DYN_setup
    use mod_atmos_phy_sf, only: &
       ATMOS_PHY_SF_setup
    use mod_atmos_phy_tb, only: &
       ATMOS_PHY_TB_setup
    use mod_atmos_phy_mp, only: &
       ATMOS_PHY_MP_setup
    use mod_atmos_phy_rd, only: &
       ATMOS_PHY_RD_setup
    implicit none
    !---------------------------------------------------------------------------

    ! setup common tools
    call ATMOS_THERMODYN_setup

    call ATMOS_SATURATION_setup

    ! setup variable container
    call ATMOS_vars_setup

    call ATMOS_vars_sf_setup

    ! read restart
    call ATMOS_vars_restart_read

    call ATMOS_REFSTATE_setup

    call ATMOS_BOUNDARY_setup

    ! setup each components
    if ( sw_dyn    ) call ATMOS_DYN_setup

    if ( sw_phy_sf ) call ATMOS_PHY_SF_setup

    if ( sw_phy_tb ) call ATMOS_PHY_TB_setup

    if ( sw_phy_mp ) call ATMOS_PHY_MP_setup

    if ( sw_phy_rd ) call ATMOS_PHY_RD_setup

    return
  end subroutine ATMOS_setup

  !-----------------------------------------------------------------------------
  !> advance atmospheric state
  !-----------------------------------------------------------------------------
  subroutine ATMOS_step
    use mod_time, only: &
       do_dyn    => TIME_DOATMOS_DYN,    &
       do_phy_sf => TIME_DOATMOS_PHY_SF, &
       do_phy_tb => TIME_DOATMOS_PHY_TB, &
       do_phy_mp => TIME_DOATMOS_PHY_MP, &
       do_phy_rd => TIME_DOATMOS_PHY_RD
    use mod_atmos_vars, only: &
       sw_dyn    => ATMOS_sw_dyn,    &
       sw_phy_sf => ATMOS_sw_phy_sf, &
       sw_phy_tb => ATMOS_sw_phy_tb, &
       sw_phy_mp => ATMOS_sw_phy_mp, &
       sw_phy_rd => ATMOS_sw_phy_rd, &
       ATMOS_vars_history
    use mod_atmos_dyn, only: &
       ATMOS_DYN
    use mod_atmos_phy_sf, only: &
       ATMOS_PHY_SF
    use mod_atmos_phy_tb, only: &
       ATMOS_PHY_TB
    use mod_atmos_phy_mp, only: &
       ATMOS_PHY_MP
    use mod_atmos_phy_rd, only: &
       ATMOS_PHY_RD
    implicit none
    !---------------------------------------------------------------------------

    !########## Dynamics ##########
    call TIME_rapstart('Dynamics')
    if ( sw_dyn .AND. do_dyn ) then
       call ATMOS_DYN
    endif
    call TIME_rapend  ('Dynamics')

    !########## Surface Flux ##########

    call TIME_rapstart('SurfaceFlux')
    if ( sw_phy_sf .AND. do_phy_sf ) then
       call ATMOS_PHY_SF
    endif
    call TIME_rapend  ('SurfaceFlux')

    !########## Turbulence ##########

    call TIME_rapstart('Turbulence')
    if ( sw_phy_tb .AND. do_phy_tb ) then
       call ATMOS_PHY_TB
    endif
    call TIME_rapend  ('Turbulence')

    !########## Microphysics ##########
    call TIME_rapstart('Microphysics')
    if ( sw_phy_mp .AND. do_phy_mp ) then
       call ATMOS_PHY_MP
    endif
    call TIME_rapend  ('Microphysics')

    !########## Radiation ##########
    call TIME_rapstart('Radiation')
    if ( sw_phy_rd .AND. do_phy_rd ) then
       call ATMOS_PHY_RD
    endif
    call TIME_rapend  ('Radiation')

    !########## History&Monitor ##########
    call TIME_rapstart('History')
       call ATMOS_vars_history
    call TIME_rapend  ('History')

    return
  end subroutine ATMOS_step

end module mod_atmos
