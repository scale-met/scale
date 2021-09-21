!-------------------------------------------------------------------------------
!> Module Atmospheric Physics driver
!!
!! @par Description
!!          This module for physical processes
!!
!! @author Team SCALE
!<
!-------------------------------------------------------------------------------
#include "scalelib.h"
module mod_atmos_phy_driver
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use scale_precision
  use scale_io
  use scale_prof
  use scale_atmos_grid_icoA_index
  use scale_tracer

  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: atmos_phy_driver_tracer_setup
  public :: atmos_phy_driver_setup
  public :: atmos_phy_driver
  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private procedures
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Tracer setup
  subroutine atmos_phy_driver_tracer_setup
    use mod_atmos_phy_mp_driver, only: &
       ATMOS_PHY_MP_driver_tracer_setup
    use mod_atmos_phy_bl_driver, only: &
       ATMOS_PHY_BL_driver_tracer_setup
    implicit none
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("atmos_phy_driver_tracer_setup",*) 'Setup'

    call ATMOS_PHY_MP_driver_tracer_setup
    call ATMOS_PHY_BL_driver_tracer_setup

    return
  end subroutine atmos_phy_driver_tracer_setup

  !-----------------------------------------------------------------------------
  !> Setup
  subroutine atmos_phy_driver_setup
    use mod_atmos_phy_mp_driver, only: &
       ATMOS_PHY_MP_driver_setup
    use mod_atmos_phy_rd_driver, only: &
       ATMOS_PHY_RD_driver_setup
    use mod_atmos_phy_sf_driver, only: &
       ATMOS_PHY_SF_driver_setup
    use mod_atmos_phy_bl_driver, only: &
       ATMOS_PHY_BL_driver_setup
    use scale_atmos_solarins, only: &
       ATMOS_SOLARINS_setup
    use scale_time, only: &
       TIME_NOWDATE
    implicit none
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("atmos_phy_driver_setup",*) 'Setup'

    LOG_NEWLINE
    LOG_INFO("atmos_phy_driver_setup",*) 'Setup each atmospheric components ...'

    !---< solar insolation setup >---
    call atmos_solarins_setup( 0.0_RP, 0.0_RP, TIME_NOWDATE(1) )

    ! setup each components
    call ATMOS_PHY_MP_driver_setup
    call ATMOS_PHY_RD_driver_setup
    call ATMOS_PHY_SF_driver_setup
    call ATMOS_PHY_BL_driver_setup

    LOG_NEWLINE
    LOG_INFO("atmos_phy_driver_setup",*) 'Finish setup of each atmospheric components.'

    return
  end subroutine atmos_phy_driver_setup

  !-----------------------------------------------------------------------------
  subroutine atmos_phy_driver
    use mod_atmos_vars, only: &
       ATMOS_vars_calc_diagnostics_fromIcoGrid, &
       ATMOS_vars_calc_diagnostics, &
       ATMOS_vars_calc_prognostics
    use mod_atmos_phy_mp_driver, only: &
       ATMOS_PHY_MP_driver_step, &
       ATMOS_PHY_MP_driver_adjustment
    use mod_atmos_phy_bl_driver, only: &
       ATMOS_PHY_BL_driver_step
    use mod_atmos_phy_sf_driver, only: &
       ATMOS_PHY_SF_driver_step
    use mod_atmos_phy_rd_driver, only: &
       ATMOS_PHY_RD_driver_step
    use mod_atmos_admin, only: &
       ATMOS_sw_phy_mp,   &
       ATMOS_sw_phy_rd,   &
       ATMOS_sw_phy_sf,   &
       ATMOS_sw_phy_bl
    use mod_prgvar, only: &
       prgvar_get_in, &
       prgvar_set_in
    implicit none

    real(RP) :: rhog  (ADM_gall_in,ADM_kall,ADM_lall)
    real(RP) :: rhogvx(ADM_gall_in,ADM_kall,ADM_lall)
    real(RP) :: rhogvy(ADM_gall_in,ADM_kall,ADM_lall)
    real(RP) :: rhogvz(ADM_gall_in,ADM_kall,ADM_lall)
    real(RP) :: rhogw (ADM_gall_in,ADM_kall,ADM_lall)
    real(RP) :: rhoge (ADM_gall_in,ADM_kall,ADM_lall)
    real(RP) :: rhogq (ADM_gall_in,ADM_kall,ADM_lall,QA)

    !---------------------------------------------------------------------------

    call PROF_rapstart  ('__Physics',1)

    call prgvar_get_in( &
       rhog(:,:,:),                                               & ! [OUT]
       rhogvx(:,:,:), rhogvy(:,:,:), rhogvz(:,:,:), rhogw(:,:,:), & ! [OUT]
       rhoge(:,:,:), rhogq(:,:,:,:)                               ) ! [OUT]

    call ATMOS_vars_calc_diagnostics_fromIcoGrid( &
       rhog(:,:,:),                                               & ! [IN]
       rhogvx(:,:,:), rhogvy(:,:,:), rhogvz(:,:,:), rhogw(:,:,:), & ! [IN]
       rhoge(:,:,:), rhogq(:,:,:,:)                               ) ! [IN]


    call atmos_vars_calc_diagnostics


    if( atmos_sw_phy_mp ) then
       call PROF_rapstart  ('___Microphysics',1)
       call atmos_phy_mp_driver_step
       call atmos_phy_mp_driver_adjustment
       call atmos_vars_calc_diagnostics
       call PROF_rapend    ('___Microphysics',1)
    end if
    if( atmos_sw_phy_rd ) then
       call PROF_rapstart  ('___Radiation',1)
       call atmos_phy_rd_driver_step
       call atmos_vars_calc_diagnostics
       call PROF_rapend    ('___Radiation',1)
    end if
    if( atmos_sw_phy_sf ) then
       call PROF_rapstart  ('___Surfacefluxes',1)
       call atmos_phy_sf_driver_step
       call atmos_vars_calc_diagnostics
       call PROF_rapend    ('___Surfacefluxes',1)
    end if
    if( atmos_sw_phy_bl ) then
       call PROF_rapstart  ('___Boundarylayer',1)
       call atmos_phy_bl_driver_step
       call atmos_vars_calc_diagnostics
       call PROF_rapend    ('___Boundarylayer',1)
    end if


    call atmos_vars_calc_prognostics( &
       rhog(:,:,:),                                               & ! [OUT]
       rhogvx(:,:,:), rhogvy(:,:,:), rhogvz(:,:,:), rhogw(:,:,:), & ! [OUT]
       rhoge(:,:,:), rhogq(:,:,:,:)                               ) ! [OUT]


    !--- set the prognostic variables
    call prgvar_set_in( &
       rhog(:,:,:),                                               & ! [IN]
       rhogvx(:,:,:), rhogvy(:,:,:), rhogvz(:,:,:), rhogw(:,:,:), & ! [IN]
       rhoge(:,:,:), rhogq(:,:,:,:)                               ) ! [IN]

    call PROF_rapend  ('__Physics',1)

    return
  end subroutine atmos_phy_driver

end module mod_atmos_phy_driver
