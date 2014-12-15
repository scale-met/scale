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

  use scale_const, only: &
     I_SW  => CONST_I_SW, &
     I_LW  => CONST_I_LW
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: ATMOS_driver_setup
  public :: ATMOS_driver
  public :: ATMOS_driver_finalize
  public :: ATMOS_SURFACE_GET
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
       TIME_NOWDATE,    &
       TIME_OFFSET_YEAR
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
       RHOQ_tp
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
    if( IO_L ) write(IO_FID_LOG,*) '++++++ Module[DRIVER] / Categ[ATMOS] / Origin[SCALE-LES]'

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** Setup each atmospheric components...'

    !--- setup solar insolation
    call ATMOS_SOLARINS_setup( TIME_NOWDATE(1)+TIME_OFFSET_YEAR )

    call ATMOS_REFSTATE_setup( DENS, RHOT, QTRC ) ! (in)

    call ATMOS_BOUNDARY_setup( DENS, MOMZ, MOMX, MOMY, RHOT, QTRC ) ! (in)

    !########## initialize tendencies ##########
    DENS_tp(:,:,:)   = 0.0_RP
    MOMZ_tp(:,:,:)   = 0.0_RP
    MOMX_tp(:,:,:)   = 0.0_RP
    MOMY_tp(:,:,:)   = 0.0_RP
    RHOT_tp(:,:,:)   = 0.0_RP
    RHOQ_tp(:,:,:,:) = 0.0_RP

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

    !########## History & Monitor ##########
    call PROF_rapstart('ATM History Vars', 1)
    call ATMOS_vars_history
    call PROF_rapend  ('ATM History Vars', 1)

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** Finish setup of each atmospheric components.'

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
       ATMOS_vars_monitor,     &
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
       RHOQ_tp
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

    !########## Lateral/Top Boundary Condition ###########
    if ( ATMOS_BOUNDARY_UPDATE_FLAG ) then
       call ATMOS_BOUNDARY_update( DENS, MOMZ, MOMX, MOMY, RHOT, QTRC ) ! (in)
    endif

    !########## Get Surface Boundary Condition from coupler ##########
    call ATMOS_SURFACE_GET( setup=.false. )

    !########## Dynamics ##########
    if ( ATMOS_sw_dyn ) then
       call PROF_rapstart('ATM Dynamics', 1)
       call ATMOS_DYN_driver( do_dyn )
       call PROF_rapend  ('ATM Dynamics', 1)
    endif

    !########## reset tendencies ##########
    DENS_tp(:,:,:)   = 0.0_RP
    MOMZ_tp(:,:,:)   = 0.0_RP
    MOMX_tp(:,:,:)   = 0.0_RP
    MOMY_tp(:,:,:)   = 0.0_RP
    RHOT_tp(:,:,:)   = 0.0_RP
    RHOQ_tp(:,:,:,:) = 0.0_RP

    !########## Microphysics ##########
    if ( ATMOS_sw_phy_mp ) then
       call PROF_rapstart('ATM Microphysics', 1)
       call ATMOS_PHY_MP_driver( update_flag = do_phy_mp )
       call PROF_rapend  ('ATM Microphysics', 1)
    endif

    !########## Aerosol ##########
    if ( ATMOS_sw_phy_ae ) then
       call PROF_rapstart('ATM Aerosol', 1)
       call ATMOS_PHY_AE_driver( update_flag = do_phy_ae )
       call PROF_rapend  ('ATM Aerosol', 1)
    endif

    !########## Chemistry ##########
    if ( ATMOS_sw_phy_ch ) then
       call PROF_rapstart('ATM Chemistry', 1)
       call ATMOS_PHY_CH_driver( update_flag = do_phy_ch )
       call PROF_rapend  ('ATM Chemistry', 1)
    endif

    !########## Radiation ##########
    if ( ATMOS_sw_phy_rd ) then
       call PROF_rapstart('ATM Radiation', 1)
       call ATMOS_PHY_RD_driver( update_flag = do_phy_rd )
       call PROF_rapend  ('ATM Radiation', 1)
    endif

    !########## Surface Flux ##########
    if ( ATMOS_sw_phy_sf ) then
       call PROF_rapstart('ATM SurfaceFlux', 1)
       call ATMOS_PHY_SF_driver( update_flag = do_phy_sf )
       call PROF_rapend  ('ATM SurfaceFlux', 1)
    endif

    !########## Turbulence ##########
    if ( ATMOS_sw_phy_tb ) then
       call PROF_rapstart('ATM Turbulence', 1)
       call ATMOS_PHY_TB_driver( update_flag = do_phy_tb )
       call PROF_rapend  ('ATM Turbulence', 1)
    endif

    !########## Cumulus ##########
    if ( ATMOS_sw_phy_cp ) then
       call PROF_rapstart('ATM Cumulus', 1)
       call ATMOS_PHY_CP_driver( update_flag = do_phy_cp )
       call PROF_rapend  ('ATM Cumulus', 1)
    endif

    !########## Calculate diagnostic variables ##########
    call ATMOS_vars_diagnostics

    !########## Set Surface Boundary & Put to coupler ##########
    call ATMOS_SURFACE_SET( setup=.false. )

    !########## History & Monitor ##########
    call PROF_rapstart('ATM History Vars', 1)
    call ATMOS_vars_history
    call ATMOS_vars_monitor
    call PROF_rapend  ('ATM History Vars', 1)

    return
  end subroutine ATMOS_driver

  !-----------------------------------------------------------------------------
  !> Finalize
  subroutine ATMOS_driver_finalize
    use scale_atmos_boundary, only: &
       ATMOS_BOUNDARY_UPDATE_FLAG, &
       ATMOS_BOUNDARY_update
    use scale_grid_nest, only: &
       NEST_COMM_disconnect
    use mod_atmos_vars, only: &
       DENS, &
       MOMZ, &
       MOMX, &
       MOMY, &
       RHOT, &
       QTRC

    !########## Lateral/Top Boundary Condition ###########
    if ( ATMOS_BOUNDARY_UPDATE_FLAG ) then
       ! If this run is parent of online nesting, boundary data must be sent
       call ATMOS_BOUNDARY_update( DENS, MOMZ, MOMX, MOMY, RHOT, QTRC, .true. ) ! (in)

       ! Finialize Inter-Communicators
       call NEST_COMM_disconnect
    endif

  end subroutine ATMOS_driver_finalize


  !-----------------------------------------------------------------------------
  !> Get surface boundary condition (CPL2ATM)
  subroutine ATMOS_SURFACE_GET( setup )
    use scale_time, only: &
       TIME_DTSEC
    use scale_history, only: &
       HIST_in
    use scale_topography, only: &
       TOPO_Zsfc
    use scale_grid_real, only: &
       REAL_CZ, &
       REAL_Z1
    use scale_atmos_bottom, only: &
       BOTTOM_estimate => ATMOS_BOTTOM_estimate
    use mod_atmos_vars, only: &
       DENS, &
       PRES
    use mod_atmos_phy_sf_vars, only: &
       SFC_DENS   => ATMOS_PHY_SF_SFC_DENS,   &
       SFC_PRES   => ATMOS_PHY_SF_SFC_PRES,   &
       SFC_TEMP   => ATMOS_PHY_SF_SFC_TEMP,   &
       SFC_albedo => ATMOS_PHY_SF_SFC_albedo, &
       SFC_Z0M    => ATMOS_PHY_SF_SFC_Z0M,    &
       SFC_Z0H    => ATMOS_PHY_SF_SFC_Z0H,    &
       SFC_Z0E    => ATMOS_PHY_SF_SFC_Z0E,    &
       SFLX_MW    => ATMOS_PHY_SF_SFLX_MW,    &
       SFLX_MU    => ATMOS_PHY_SF_SFLX_MU,    &
       SFLX_MV    => ATMOS_PHY_SF_SFLX_MV,    &
       SFLX_SH    => ATMOS_PHY_SF_SFLX_SH,    &
       SFLX_LH    => ATMOS_PHY_SF_SFLX_LH,    &
       SFLX_QTRC  => ATMOS_PHY_SF_SFLX_QTRC
    use mod_cpl_admin, only: &
       CPL_sw
    use mod_cpl_vars, only: &
       CPL_getATM,   &
       CPL_getATM_SF
    implicit none

    logical, intent(in) :: setup

    real(RP) :: Uabs10  (IA,JA) ! 10m absolute wind [m/s]
    real(RP) :: U10     (IA,JA) ! 10m x-wind [m/s]
    real(RP) :: V10     (IA,JA) ! 10m y-wind [m/s]
    real(RP) :: T2      (IA,JA) !  2m Temp   [K]
    real(RP) :: Q2      (IA,JA) !  2m Vapor  [kg/kg]
    real(RP) :: FLX_heat(IA,JA) ! ground heat flux [J/m2/s]
    !---------------------------------------------------------------------------

    if ( setup ) then
       SFLX_MW  (:,:)   = 0.0_RP
       SFLX_MU  (:,:)   = 0.0_RP
       SFLX_MV  (:,:)   = 0.0_RP
       SFLX_SH  (:,:)   = 0.0_RP
       SFLX_LH  (:,:)   = 0.0_RP
       SFLX_QTRC(:,:,:) = 0.0_RP

       if ( CPL_sw ) then
          call CPL_getATM_SF( SFLX_MW  (:,:),  & ! [OUT]
                              SFLX_MU  (:,:),  & ! [OUT]
                              SFLX_MV  (:,:),  & ! [OUT]
                              SFLX_SH  (:,:),  & ! [OUT]
                              SFLX_LH  (:,:),  & ! [OUT]
                              SFLX_QTRC(:,:,:) ) ! [OUT]
       endif
    endif

    if ( CPL_sw ) then
       Uabs10   (:,:)   = 0.0_RP
       U10      (:,:)   = 0.0_RP
       V10      (:,:)   = 0.0_RP
       T2       (:,:)   = 0.0_RP
       Q2       (:,:)   = 0.0_RP
       FLX_heat (:,:)   = 0.0_RP

       call CPL_getATM( SFC_TEMP  (:,:),   & ! [OUT]
                        SFC_albedo(:,:,:), & ! [OUT]
                        SFC_Z0M   (:,:),   & ! [OUT]
                        SFC_Z0H   (:,:),   & ! [OUT]
                        SFC_Z0E   (:,:),   & ! [OUT]
                        Uabs10    (:,:),   & ! [OUT]
                        U10       (:,:),   & ! [OUT]
                        V10       (:,:),   & ! [OUT]
                        T2        (:,:),   & ! [OUT]
                        Q2        (:,:),   & ! [OUT]
                        FLX_heat  (:,:)    ) ! [OUT]

       call HIST_in( SFC_Z0M(:,:), 'SFC_Z0M', 'roughness length (momentum)', 'm'     )
       call HIST_in( SFC_Z0H(:,:), 'SFC_Z0H', 'roughness length (heat)',     'm'     )
       call HIST_in( SFC_Z0E(:,:), 'SFC_Z0E', 'roughness length (vapor)',    'm'     )
       call HIST_in( Uabs10 (:,:), 'Uabs10',  '10m absolute wind',           'm/s'   )
       call HIST_in( U10    (:,:), 'U10',     '10m x-wind',                  'm/s'   )
       call HIST_in( V10    (:,:), 'V10',     '10m y-wind',                  'm/s'   )
       call HIST_in( T2     (:,:), 'T2 ',     '2m temperature',              'K'     )
       call HIST_in( Q2     (:,:), 'Q2 ',     '2m water vapor',              'kg/kg' )

       call HIST_in( FLX_heat(:,:), 'GHFLX', 'ground heat flux (merged)', 'W/m2' )
    endif

    ! calculate surface density, surface pressure
    call BOTTOM_estimate( DENS     (:,:,:), & ! [IN]
                          PRES     (:,:,:), & ! [IN]
                          REAL_CZ  (:,:,:), & ! [IN]
                          TOPO_Zsfc(:,:),   & ! [IN]
                          REAL_Z1  (:,:),   & ! [IN]
                          SFC_DENS (:,:),   & ! [OUT]
                          SFC_PRES (:,:)    ) ! [OUT]

    call HIST_in( SFC_DENS  (:,:),      'SFC_DENS',   'surface atmospheric density',       'kg/m3' )
    call HIST_in( SFC_PRES  (:,:),      'SFC_PRES',   'surface atmospheric pressure',      'Pa'    )

    ! if coupler is disabled, SFC_TEMP, SFC_albedo is set in ATMOS_PHY_SF_vars
    call HIST_in( SFC_TEMP  (:,:),      'SFC_TEMP',   'surface skin temperature (merged)', 'K'     )
    call HIST_in( SFC_albedo(:,:,I_LW), 'SFC_ALB_LW', 'surface albedo (longwave, merged)', '0-1'   )
    call HIST_in( SFC_albedo(:,:,I_SW), 'SFC_ALB_SW', 'surface albedo (shortwave,merged)', '0-1'   )

    return
  end subroutine ATMOS_SURFACE_GET

  !-----------------------------------------------------------------------------
  !> Set surface boundary condition (ATM->CPL)
  subroutine ATMOS_SURFACE_SET( setup )
    use scale_topography, only: &
       TOPO_Zsfc
    use scale_grid_real, only: &
       REAL_CZ, &
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
    use mod_admin_restart, only: &
       RESTART_RUN
    use mod_atmos_phy_mp_vars, only: &
       SFLX_rain => ATMOS_PHY_MP_SFLX_rain, &
       SFLX_snow => ATMOS_PHY_MP_SFLX_snow
    use mod_atmos_phy_rd_vars, only: &
       SFLX_LW_dn => ATMOS_PHY_RD_SFLX_LW_dn, &
       SFLX_SW_dn => ATMOS_PHY_RD_SFLX_SW_dn
    use mod_atmos_phy_sf_vars, only: &
       SFC_DENS   => ATMOS_PHY_SF_SFC_DENS,   &
       SFC_PRES   => ATMOS_PHY_SF_SFC_PRES,   &
       SFC_TEMP   => ATMOS_PHY_SF_SFC_TEMP,   &
       SFC_albedo => ATMOS_PHY_SF_SFC_albedo, &
       SFC_Z0M    => ATMOS_PHY_SF_SFC_Z0M,    &
       SFC_Z0H    => ATMOS_PHY_SF_SFC_Z0H,    &
       SFC_Z0E    => ATMOS_PHY_SF_SFC_Z0E,    &
       SFLX_MW    => ATMOS_PHY_SF_SFLX_MW,    &
       SFLX_MU    => ATMOS_PHY_SF_SFLX_MU,    &
       SFLX_MV    => ATMOS_PHY_SF_SFLX_MV,    &
       SFLX_SH    => ATMOS_PHY_SF_SFLX_SH,    &
       SFLX_LH    => ATMOS_PHY_SF_SFLX_LH,    &
       SFLX_QTRC  => ATMOS_PHY_SF_SFLX_QTRC
    use mod_cpl_admin, only: &
       CPL_sw
    use mod_cpl_vars, only: &
       CPL_putAtm, &
       CPL_putAtm_restart
    implicit none

    logical, intent(in) :: setup

    ! works
    real(RP) :: ATM_PBL (IA,JA)
    !---------------------------------------------------------------------------

    ! update surface density, surface pressure
    call BOTTOM_estimate( DENS     (:,:,:), & ! [IN]
                          PRES     (:,:,:), & ! [IN]
                          REAL_CZ  (:,:,:), & ! [IN]
                          TOPO_Zsfc(:,:),   & ! [IN]
                          REAL_Z1  (:,:),   & ! [IN]
                          SFC_DENS (:,:),   & ! [OUT]
                          SFC_PRES (:,:)    ) ! [OUT]

    ! planetary boundary layer
    ATM_PBL(:,:) = 1.0E+2_RP ! tentative

    if ( CPL_sw ) then
       if ( setup .and. RESTART_RUN ) then
          call CPL_putAtm_restart( SFC_TEMP  (:,:),   & ! [IN]
                                   SFC_albedo(:,:,:), & ! [IN]
                                   SFC_Z0M   (:,:),   & ! [IN]
                                   SFC_Z0H   (:,:),   & ! [IN]
                                   SFC_Z0E   (:,:),   & ! [IN]
                                   SFLX_MW   (:,:),   & ! [IN]
                                   SFLX_MU   (:,:),   & ! [IN]
                                   SFLX_MV   (:,:),   & ! [IN]
                                   SFLX_SH   (:,:),   & ! [IN]
                                   SFLX_LH   (:,:),   & ! [IN]
                                   SFLX_QTRC (:,:,:)  ) ! [IN]
       else
          call CPL_putAtm( TEMP      (KS,:,:),   & ! [IN]
                           PRES      (KS,:,:),   & ! [IN]
                           W         (KS,:,:),   & ! [IN]
                           U         (KS,:,:),   & ! [IN]
                           V         (KS,:,:),   & ! [IN]
                           DENS      (KS,:,:),   & ! [IN]
                           QTRC      (KS,:,:,:), & ! [IN]
                           ATM_PBL   (:,:),      & ! [IN]
                           SFC_PRES  (:,:),      & ! [IN]
                           SFLX_LW_dn(:,:),      & ! [IN]
                           SFLX_SW_dn(:,:),      & ! [IN]
                           SFLX_rain (:,:),      & ! [IN]
                           SFLX_snow (:,:)       ) ! [IN]
       endif
    endif

    return
  end subroutine ATMOS_SURFACE_SET

end module mod_atmos_driver
