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
!! @li      2013-08-31 (T.Yamaura) [add] Coupler
!!
!<
!-------------------------------------------------------------------------------
module mod_atmos
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
    use scale_time, only: &
       TIME_NOWDATE
    use mod_atmos_vars, only: &
       ATMOS_DYN_TYPE, &
       ATMOS_PHY_SF_TYPE, &
       ATMOS_PHY_TB_TYPE, &
       ATMOS_PHY_MP_TYPE, &
       ATMOS_PHY_RD_TYPE, &
       ATMOS_PHY_AE_TYPE, &
       ATMOS_PHY_CH_TYPE, &
       DENS, &
       MOMZ, &
       MOMX, &
       MOMY, &
       RHOT, &
       QTRC, &
       DENS_tp, &
       MOMZ_tp, &
       MOMX_tp, &
       MOMY_tp, &
       RHOT_tp, &
       QTRC_tp, &
       sw_dyn    => ATMOS_sw_dyn,    &
       sw_phy_sf => ATMOS_sw_phy_sf, &
       sw_phy_tb => ATMOS_sw_phy_tb, &
       sw_phy_mp => ATMOS_sw_phy_mp, &
       sw_phy_rd => ATMOS_sw_phy_rd, &
       sw_phy_ae => ATMOS_sw_phy_ae
    use scale_atmos_solarins, only: &
       ATMOS_SOLARINS_setup
    use scale_atmos_refstate, only: &
       ATMOS_REFSTATE_setup
    use scale_atmos_boundary, only: &
       ATMOS_BOUNDARY_setup
    use mod_atmos_dyn_driver, only: &
       ATMOS_DYN_setup => ATMOS_DYN_driver_setup
    use mod_atmos_phy_sf_driver, only: &
       ATMOS_PHY_SF_setup => ATMOS_PHY_SF_driver_setup
    use mod_atmos_phy_tb_driver, only: &
       ATMOS_PHY_TB_setup => ATMOS_PHY_TB_driver_setup
    use mod_atmos_phy_mp_driver, only: &
       ATMOS_PHY_MP_setup => ATMOS_PHY_MP_driver_setup
    use mod_atmos_phy_rd_driver, only: &
       ATMOS_PHY_RD_setup => ATMOS_PHY_RD_driver_setup
    use mod_atmos_phy_ae_driver, only: &
       ATMOS_PHY_AE_setup => ATMOS_PHY_AE_driver_setup
    implicit none
    !---------------------------------------------------------------------------

    !--- setup solar insolation
    call ATMOS_SOLARINS_setup( TIME_NOWDATE(1) )

    call ATMOS_REFSTATE_setup( DENS, RHOT, QTRC ) ! (in)

    call ATMOS_BOUNDARY_setup( &
         DENS, MOMZ, MOMX, MOMY, RHOT, QTRC ) ! (in)

    ! setup each components
    if ( sw_dyn    ) call ATMOS_DYN_setup   ( ATMOS_DYN_TYPE    )

    if ( sw_phy_sf ) call ATMOS_PHY_SF_setup( ATMOS_PHY_SF_TYPE )

    if ( sw_phy_tb ) call ATMOS_PHY_TB_setup( ATMOS_PHY_TB_TYPE )

    if ( sw_phy_mp ) call ATMOS_PHY_MP_setup( ATMOS_PHY_MP_TYPE )

    if ( sw_phy_rd .or. sw_phy_ae ) call ATMOS_PHY_AE_setup( ATMOS_PHY_AE_TYPE )

    if ( sw_phy_rd ) call ATMOS_PHY_RD_setup( ATMOS_PHY_RD_TYPE )

    !########## initialize tendencies ##########
    DENS_tp(:,:,:)   = 0.0_RP
    MOMZ_tp(:,:,:)   = 0.0_RP
    MOMX_tp(:,:,:)   = 0.0_RP
    MOMY_tp(:,:,:)   = 0.0_RP
    RHOT_tp(:,:,:)   = 0.0_RP
    QTRC_tp(:,:,:,:) = 0.0_RP

    return
  end subroutine ATMOS_setup

  !-----------------------------------------------------------------------------
  !> advance atmospheric state
  !-----------------------------------------------------------------------------
  subroutine ATMOS_step
    use scale_time, only: &
       do_dyn    => TIME_DOATMOS_DYN,    &
       do_phy_sf => TIME_DOATMOS_PHY_SF, &
       do_phy_tb => TIME_DOATMOS_PHY_TB, &
       do_phy_mp => TIME_DOATMOS_PHY_MP, &
       do_phy_rd => TIME_DOATMOS_PHY_RD, &
       do_phy_ae => TIME_DOATMOS_PHY_AE
    use mod_atmos_vars, only: &
       DENS, &
       RHOT, &
       QTRC, &
       DENS_tp, &
       MOMZ_tp, &
       MOMX_tp, &
       MOMY_tp, &
       RHOT_tp, &
       QTRC_tp, &
       sw_dyn    => ATMOS_sw_dyn,    &
       sw_phy_sf => ATMOS_sw_phy_sf, &
       sw_phy_tb => ATMOS_sw_phy_tb, &
       sw_phy_mp => ATMOS_sw_phy_mp, &
       sw_phy_rd => ATMOS_sw_phy_rd, &
       sw_phy_ae => ATMOS_sw_phy_ae, &
       ATMOS_vars_history
    use mod_atmos_dyn_driver, only: &
       ATMOS_DYN => ATMOS_DYN_driver
    use mod_atmos_phy_sf_driver, only: &
       ATMOS_PHY_SF => ATMOS_PHY_SF_driver, &
       ATMOS_PHY_SF_driver_first,  &
       ATMOS_PHY_SF_driver_final
    use mod_atmos_phy_tb_driver, only: &
       ATMOS_PHY_TB => ATMOS_PHY_TB_driver
    use mod_atmos_phy_mp_driver, only: &
       ATMOS_PHY_MP => ATMOS_PHY_MP_driver
    use mod_atmos_phy_rd_driver, only: &
       ATMOS_PHY_RD => ATMOS_PHY_RD_driver
    use mod_atmos_phy_ae_driver, only: &
       ATMOS_PHY_AE => ATMOS_PHY_AE_driver
    use scale_atmos_refstate, only: &
       ATMOS_REFSTATE_update, &
       ATMOS_REFSTATE_UPDATE_FLAG
    use mod_cpl_vars, only: &
       sw_AtmLnd => CPL_sw_AtmLnd, &
       sw_AtmOcn => CPL_sw_AtmOcn
    implicit none
    !---------------------------------------------------------------------------

    !########## Reference State ###########
    if ( ATMOS_REFSTATE_UPDATE_FLAG ) then
       call ATMOS_REFSTATE_update( DENS, RHOT, QTRC ) ! (in)
    endif

    !########## Surface First ##########
    if ( sw_phy_sf ) then
       call PROF_rapstart('ATM SurfaceFlux')
       if ( sw_AtmLnd .or. sw_AtmOcn ) then
          call ATMOS_PHY_SF_driver_first
       else
          call ATMOS_PHY_SF( do_phy_sf, .true. )
       endif
       call PROF_rapend  ('ATM SurfaceFlux')
    endif

    !########## Turbulence ##########
    if ( sw_phy_tb ) then
       call PROF_rapstart('ATM Turbulence')
       call ATMOS_PHY_TB( do_phy_tb, .true. )
       call PROF_rapend  ('ATM Turbulence')
    endif

    !########## Microphysics ##########
    if ( sw_phy_mp ) then
       call PROF_rapstart('ATM Microphysics')
       call ATMOS_PHY_MP( do_phy_mp, .true. )
       call PROF_rapend  ('ATM Microphysics')
    endif

    !########## Aerosol ##########
    if ( sw_phy_ae ) then
       call PROF_rapstart('ATM Aerosol')
       call ATMOS_PHY_AE( do_phy_ae, .true. )
       call PROF_rapend  ('ATM Aerosol')
    endif

    !########## Radiation ##########
    if ( sw_phy_rd ) then
       call PROF_rapstart('ATM Radiation')
       call ATMOS_PHY_RD( do_phy_rd, .true. )
       call PROF_rapend  ('ATM Radiation')
    endif

    !########## Dynamics ##########
    if ( sw_dyn ) then
       call PROF_rapstart('ATM Dynamics')
       call ATMOS_DYN( do_dyn )
       call PROF_rapend  ('ATM Dynamics')
    endif

    !########## Surface Final ##########
    if ( sw_phy_sf ) then
       call PROF_rapstart('ATM SurfaceFlux')
       if ( sw_AtmLnd .or. sw_AtmOcn ) then
          call ATMOS_PHY_SF_driver_final
       endif
       call PROF_rapend  ('ATM SurfaceFlux')
    endif

    !########## History & Monitor ##########
    call PROF_rapstart('ATM History Vars')
    call ATMOS_vars_history
    call PROF_rapend  ('ATM History Vars')

    !########## reset tendencies ##########
    DENS_tp(:,:,:)   = 0.0_RP
    MOMZ_tp(:,:,:)   = 0.0_RP
    MOMX_tp(:,:,:)   = 0.0_RP
    MOMY_tp(:,:,:)   = 0.0_RP
    RHOT_tp(:,:,:)   = 0.0_RP
    QTRC_tp(:,:,:,:) = 0.0_RP

    return
  end subroutine ATMOS_step

end module mod_atmos
