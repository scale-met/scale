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
module mod_atmos_driver
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
  public :: ATMOS_driver_setup
  public :: ATMOS_driver
  public :: ATMOS_SURFACE_SET

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
  subroutine ATMOS_driver_setup
    use scale_time, only: &
       TIME_NOWDATE
    use mod_atmos_vars, only: &
       ATMOS_vars_diagnostics, &
       DENS,              &
       MOMZ,              &
       MOMX,              &
       MOMY,              &
       RHOT,              &
       QTRC,              &
       DENS_tp,           &
       MOMZ_tp,           &
       MOMX_tp,           &
       MOMY_tp,           &
       RHOT_tp,           &
       QTRC_tp
    use scale_atmos_solarins, only: &
       ATMOS_SOLARINS_setup
    use scale_atmos_refstate, only: &
       ATMOS_REFSTATE_setup
    use scale_atmos_boundary, only: &
       ATMOS_BOUNDARY_setup
    use mod_atmos_dyn_driver, only: &
       ATMOS_DYN_driver_setup
    use mod_atmos_phy_mp_driver, only: &
       ATMOS_PHY_MP_driver_setup
    use mod_atmos_phy_ae_driver, only: &
       ATMOS_PHY_AE_driver_setup
    use mod_atmos_phy_ch_driver, only: &
       ATMOS_PHY_CH_driver_setup
    use mod_atmos_phy_rd_driver, only: &
       ATMOS_PHY_RD_driver_setup
    use mod_atmos_phy_sf_driver, only: &
       ATMOS_PHY_SF_driver_setup
    use mod_atmos_phy_tb_driver, only: &
       ATMOS_PHY_TB_driver_setup
    use mod_atmos_phy_cp_driver, only: &
       ATMOS_PHY_CP_driver_setup
    implicit none
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[DRIVER]/Categ[ATMOS]'
    if( IO_L ) write(IO_FID_LOG,*) '*** setup each atmospheric components'

    !--- setup solar insolation
    call ATMOS_SOLARINS_setup( TIME_NOWDATE(1) )

    call ATMOS_REFSTATE_setup( DENS, RHOT, QTRC ) ! (in)

    call ATMOS_BOUNDARY_setup( DENS, MOMZ, MOMX, MOMY, RHOT, QTRC ) ! (in)

    ! setup each components
    call ATMOS_DYN_driver_setup
    call ATMOS_PHY_MP_driver_setup
    call ATMOS_PHY_AE_driver_setup
    call ATMOS_PHY_CH_driver_setup
    call ATMOS_PHY_RD_driver_setup
    call ATMOS_PHY_SF_driver_setup
    call ATMOS_PHY_TB_driver_setup
    call ATMOS_PHY_CP_driver_setup

    !########## Calculate diagnostic variables ##########
    call ATMOS_vars_diagnostics

    !########## Set surface boundary & put to other model ##########
    call ATMOS_SURFACE_SET

    !########## initialize tendencies ##########
    DENS_tp(:,:,:)   = 0.0_RP
    MOMZ_tp(:,:,:)   = 0.0_RP
    MOMX_tp(:,:,:)   = 0.0_RP
    MOMY_tp(:,:,:)   = 0.0_RP
    RHOT_tp(:,:,:)   = 0.0_RP
    QTRC_tp(:,:,:,:) = 0.0_RP

    if( IO_L ) write(IO_FID_LOG,*) '*** setup ATMOS finished.'

    return
  end subroutine ATMOS_driver_setup

  !-----------------------------------------------------------------------------
  !> advance atmospheric state
  subroutine ATMOS_driver
    use scale_time, only: &
       do_dyn    => TIME_DOATMOS_DYN,    &
       do_phy_mp => TIME_DOATMOS_PHY_MP, &
       do_phy_ae => TIME_DOATMOS_PHY_AE, &
       do_phy_ch => TIME_DOATMOS_PHY_CH, &
       do_phy_rd => TIME_DOATMOS_PHY_RD, &
       do_phy_sf => TIME_DOATMOS_PHY_SF, &
       do_phy_tb => TIME_DOATMOS_PHY_TB, &
       do_phy_cp => TIME_DOATMOS_PHY_CP
    use scale_atmos_refstate, only: &
       ATMOS_REFSTATE_UPDATE_FLAG, &
       ATMOS_REFSTATE_update
    use scale_atmos_boundary, only: &
       ATMOS_BOUNDARY_UPDATE_FLAG, &
       ATMOS_BOUNDARY_update
    use mod_atmos_admin, only: &
       ATMOS_sw_dyn,    &
       ATMOS_sw_phy_mp, &
       ATMOS_sw_phy_ae, &
       ATMOS_sw_phy_ch, &
       ATMOS_sw_phy_rd, &
       ATMOS_sw_phy_sf, &
       ATMOS_sw_phy_tb, &
       ATMOS_sw_phy_cp
    use mod_atmos_vars, only: &
       ATMOS_vars_diagnostics, &
       ATMOS_vars_history,     &
       DENS,                   &
       MOMZ,                   &
       MOMX,                   &
       MOMY,                   &
       RHOT,                   &
       QTRC,                   &
       DENS_tp,                &
       MOMZ_tp,                &
       MOMX_tp,                &
       MOMY_tp,                &
       RHOT_tp,                &
       QTRC_tp
    use mod_atmos_dyn_driver, only: &
       ATMOS_DYN_driver
    use mod_atmos_phy_mp_driver, only: &
       ATMOS_PHY_MP_driver
    use mod_atmos_phy_ae_driver, only: &
       ATMOS_PHY_AE_driver
    use mod_atmos_phy_ch_driver, only: &
       ATMOS_PHY_CH_driver
    use mod_atmos_phy_rd_driver, only: &
       ATMOS_PHY_RD_driver
    use mod_atmos_phy_sf_driver, only: &
       ATMOS_PHY_SF_driver
    use mod_atmos_phy_tb_driver, only: &
       ATMOS_PHY_TB_driver
    use mod_atmos_phy_cp_driver, only: &
       ATMOS_PHY_CP_driver
    implicit none
    !---------------------------------------------------------------------------

    !########## Reference State ###########
    if ( ATMOS_REFSTATE_UPDATE_FLAG ) then
       call ATMOS_REFSTATE_update( DENS, RHOT, QTRC ) ! (in)
    endif

    !########## Boundary Condition ###########
    if ( ATMOS_BOUNDARY_UPDATE_FLAG ) then
       call ATMOS_BOUNDARY_update()
    endif

    !########## Microphysics ##########
    if ( ATMOS_sw_phy_mp ) then
       call PROF_rapstart('ATM Microphysics')
       call ATMOS_PHY_MP_driver( update_flag=do_phy_mp, history_flag=.true. )
       call PROF_rapend  ('ATM Microphysics')
    endif

    !########## Aerosol ##########
    if ( ATMOS_sw_phy_ae ) then
       call PROF_rapstart('ATM Aerosol')
       call ATMOS_PHY_AE_driver( update_flag=do_phy_ae, history_flag=.true. )
       call PROF_rapend  ('ATM Aerosol')
    endif

    !########## Chemistry ##########
    if ( ATMOS_sw_phy_ch ) then
       call PROF_rapstart('ATM Chemistry')
       call ATMOS_PHY_CH_driver( update_flag=do_phy_ch, history_flag=.true. )
       call PROF_rapend  ('ATM Chemistry')
    endif

    !########## Radiation ##########
    if ( ATMOS_sw_phy_rd ) then
       call PROF_rapstart('ATM Radiation')
       call ATMOS_PHY_RD_driver( update_flag=do_phy_rd, history_flag=.true. )
       call PROF_rapend  ('ATM Radiation')
    endif

    !########## Surface Flux ##########
    if ( ATMOS_sw_phy_sf ) then
       call PROF_rapstart('ATM SurfaceFlux')
       call ATMOS_PHY_SF_driver( update_flag=do_phy_sf, history_flag=.true. )
       call PROF_rapend  ('ATM SurfaceFlux')
    endif

    !########## Turbulence ##########
    if ( ATMOS_sw_phy_tb ) then
       call PROF_rapstart('ATM Turbulence')
       call ATMOS_PHY_TB_driver( update_flag=do_phy_tb, history_flag=.true. )
       call PROF_rapend  ('ATM Turbulence')
    endif

    !########## Cumulus ##########
    if ( ATMOS_sw_phy_cp ) then
       call PROF_rapstart('ATM Cumulus')
       call ATMOS_PHY_CP_driver( update_flag=do_phy_cp, history_flag=.true. )
       call PROF_rapend  ('ATM Cumulus')
    endif

    !########## Dynamics ##########
    if ( ATMOS_sw_dyn ) then
       call PROF_rapstart('ATM Dynamics')
       call ATMOS_DYN_driver( do_dyn )
       call PROF_rapend  ('ATM Dynamics')
    endif

    !########## Calculate diagnostic variables ##########
    call ATMOS_vars_diagnostics

    !########## Set surface boundary & put to other model ##########
    call ATMOS_SURFACE_SET

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
  end subroutine ATMOS_driver

  !-----------------------------------------------------------------------------
  !> Set surface boundary condition (ATM2CPL)
  subroutine ATMOS_SURFACE_SET
    use scale_const, only: &
       RovCP => CONST_RovCP
    use scale_grid_real, only: &
       REAL_CZ, &
       REAL_FZ, &
       REAL_Z1
    use scale_atmos_bottom, only: &
       BOTTOM_estimate => ATMOS_BOTTOM_estimate
    use mod_atmos_vars, only: &
       DENS, &
       QTRC, &
       TEMP, &
       PRES, &
       W,    &
       U,    &
       V
    use mod_atmos_phy_mp_vars, only: &
       SFLX_rain => ATMOS_PHY_MP_SFLX_rain, &
       SFLX_snow => ATMOS_PHY_MP_SFLX_snow
    use mod_atmos_phy_rd_vars, only: &
       SFLX_LW_dn => ATMOS_PHY_RD_SFLX_LW_dn, &
       SFLX_SW_dn => ATMOS_PHY_RD_SFLX_SW_dn
    use mod_atmos_phy_sf_vars, only: &
       SFC_DENS => ATMOS_PHY_SF_SFC_DENS, &
       SFC_PRES => ATMOS_PHY_SF_SFC_PRES
    use mod_cpl_admin, only: &
       CPL_sw
    use mod_cpl_vars, only: &
       CPL_putATM
    implicit none

    ! works
    real(RP) :: SFC_TEMP(IA,JA)
    !---------------------------------------------------------------------------

    ! update surface density, surface pressure
    call BOTTOM_estimate( DENS    (:,:,:), & ! [IN]
                          PRES    (:,:,:), & ! [IN]
                          REAL_CZ (:,:,:), & ! [IN]
                          REAL_FZ (:,:,:), & ! [IN]
                          REAL_Z1 (:,:),   & ! [IN]
                          SFC_DENS(:,:),   & ! [OUT]
                          SFC_PRES(:,:)    ) ! [OUT]

    if ( CPL_sw ) then
       call CPL_putAtm( TEMP      (KS,:,:),   & ! [IN]
                        PRES      (KS,:,:),   & ! [IN]
                        W         (KS,:,:),   & ! [IN]
                        U         (KS,:,:),   & ! [IN]
                        V         (KS,:,:),   & ! [IN]
                        DENS      (KS,:,:),   & ! [IN]
                        QTRC      (KS,:,:,:), & ! [IN]
                        SFC_PRES  (:,:),      & ! [IN]
                        SFLX_LW_dn(:,:),      & ! [IN]
                        SFLX_SW_dn(:,:),      & ! [IN]
                        SFLX_rain (:,:),      & ! [IN]
                        SFLX_snow (:,:)       ) ! [IN]
    endif

    return
  end subroutine ATMOS_SURFACE_SET

end module mod_atmos_driver
