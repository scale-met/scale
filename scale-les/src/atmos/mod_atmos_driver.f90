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
  public :: ATMOS_driver_setup1
  public :: ATMOS_driver_setup2
  public :: ATMOS_driver
  public :: ATMOS_driver_firstsend
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
  subroutine ATMOS_driver_setup1
    use scale_time, only: &
       TIME_NOWDATE,    &
       TIME_OFFSET_YEAR
    use mod_atmos_vars, only: &
       ATMOS_vars_diagnostics, &
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
    use mod_atmos_phy_rd_driver, only: &
       ATMOS_PHY_RD_driver_setup
    use mod_atmos_phy_ch_driver, only: &
       ATMOS_PHY_CH_driver_setup
    use mod_atmos_phy_ae_driver, only: &
       ATMOS_PHY_AE_driver_setup
    implicit none
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '++++++ Module[DRIVER] / Categ[ATMOS] / Origin[SCALE-LES]'

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** Setup each atmospheric components 1 ...'

    !--- setup solar insolation
    call ATMOS_SOLARINS_setup( TIME_NOWDATE(1)+TIME_OFFSET_YEAR )

    call PROF_rapstart('ATM_Refstate', 2)
    call ATMOS_REFSTATE_setup( DENS, RHOT, QTRC )
    call PROF_rapend  ('ATM_Refstate', 2)

    call PROF_rapstart('ATM_Boundary', 2)
    call ATMOS_BOUNDARY_setup( DENS, MOMZ, MOMX, MOMY, RHOT, QTRC )
    call PROF_rapend  ('ATM_Boundary', 2)

    !########## Get Surface Boundary Condition ##########
    call PROF_rapstart('ATM_SfcExch', 2)
    call ATMOS_SURFACE_GET
    call PROF_rapend  ('ATM_SfcExch', 2)

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

    !########## Calculate diagnostic variables ##########
    call ATMOS_vars_diagnostics

    !########## Set Surface Boundary Condition ##########
    call PROF_rapstart('ATM_SfcExch', 2)
    call ATMOS_SURFACE_SET( countup = .false. )
    call PROF_rapend  ('ATM_SfcExch', 2)

    return
  end subroutine ATMOS_driver_setup1

  !-----------------------------------------------------------------------------
  !> Setup
  subroutine ATMOS_driver_setup2
    use mod_atmos_vars, only: &
       ATMOS_vars_history
    use mod_atmos_phy_cp_driver, only: &
       ATMOS_PHY_CP_driver_setup
    use mod_atmos_phy_sf_driver, only: &
       ATMOS_PHY_SF_driver_setup
    use mod_atmos_phy_tb_driver, only: &
       ATMOS_PHY_TB_driver_setup
    implicit none
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '++++++ Module[DRIVER] / Categ[ATMOS] / Origin[SCALE-LES]'

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** Setup each atmospheric components 2 ...'

    !########## Get Surface Boundary Condition ##########
    call PROF_rapstart('ATM_SfcExch', 2)
    call ATMOS_SURFACE_GET
    call PROF_rapend  ('ATM_SfcExch', 2)

    ! setup each components
    call ATMOS_PHY_SF_driver_setup
    call ATMOS_PHY_TB_driver_setup
    call ATMOS_PHY_CP_driver_setup

    !########## Set Surface Boundary Condition ##########
    call PROF_rapstart('ATM_SfcExch', 2)
    call ATMOS_SURFACE_SET( countup = .true. )
    call PROF_rapend  ('ATM_SfcExch', 2)

    !########## History & Monitor ##########
    call PROF_rapstart('ATM_History', 1)
    call ATMOS_vars_history
    call PROF_rapend  ('ATM_History', 1)

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** Finish setup of each atmospheric components.'

    return
  end subroutine ATMOS_driver_setup2

  !-----------------------------------------------------------------------------
  !> advance atmospheric state
  subroutine ATMOS_driver
    use mod_admin_time, only: &
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
       call PROF_rapstart('ATM_Refstate', 2)
       call ATMOS_REFSTATE_update( DENS, RHOT, QTRC ) ! [IN]
       call PROF_rapend  ('ATM_Refstate', 2)
    endif

    !########## Get Surface Boundary Condition ##########
    call PROF_rapstart('ATM_SfcExch', 2)
    call ATMOS_SURFACE_GET
    call PROF_rapend  ('ATM_SfcExch', 2)

    !########## Dynamics ##########
    if ( ATMOS_sw_dyn ) then
       call PROF_rapstart('ATM_Dynamics', 1)
       call ATMOS_DYN_driver( do_dyn )
       call PROF_rapend  ('ATM_Dynamics', 1)
    endif

    !########## Lateral/Top Boundary Condition ###########
    if ( ATMOS_BOUNDARY_UPDATE_FLAG ) then
       call PROF_rapstart('ATM_Boundary', 2)
       call ATMOS_BOUNDARY_update( DENS, MOMZ, MOMX, MOMY, RHOT, QTRC ) ! [INOUT]
       call PROF_rapend  ('ATM_Boundary', 2)
    endif

    !########## reset tendencies ##########
    DENS_tp(:,:,:)   = 0.0_RP
    MOMZ_tp(:,:,:)   = 0.0_RP
    MOMX_tp(:,:,:)   = 0.0_RP
    MOMY_tp(:,:,:)   = 0.0_RP
    RHOT_tp(:,:,:)   = 0.0_RP
    RHOQ_tp(:,:,:,:) = 0.0_RP

    !########## Calculate diagnostic variables ##########
    call PROF_rapstart('ATM_History', 1)
    call ATMOS_vars_diagnostics
    call PROF_rapend  ('ATM_History', 1)

    !########## Microphysics ##########
    if ( ATMOS_sw_phy_mp ) then
       call PROF_rapstart('ATM_Microphysics', 1)
       call ATMOS_PHY_MP_driver( update_flag = do_phy_mp )
       call PROF_rapend  ('ATM_Microphysics', 1)
    endif

    !########## Aerosol ##########
    if ( ATMOS_sw_phy_ae ) then
       call PROF_rapstart('ATM_Aerosol', 1)
       call ATMOS_PHY_AE_driver( update_flag = do_phy_ae )
       call PROF_rapend  ('ATM_Aerosol', 1)
    endif

    !########## Chemistry ##########
    if ( ATMOS_sw_phy_ch ) then
       call PROF_rapstart('ATM_Chemistry', 1)
       call ATMOS_PHY_CH_driver( update_flag = do_phy_ch )
       call PROF_rapend  ('ATM_Chemistry', 1)
    endif

    !########## Radiation ##########
    if ( ATMOS_sw_phy_rd ) then
       call PROF_rapstart('ATM_Radiation', 1)
       call ATMOS_PHY_RD_driver( update_flag = do_phy_rd )
       call PROF_rapend  ('ATM_Radiation', 1)
    endif

    !########## Surface Flux ##########
    if ( ATMOS_sw_phy_sf ) then
       call PROF_rapstart('ATM_SurfaceFlux', 1)
       call ATMOS_PHY_SF_driver( update_flag = do_phy_sf )
       call PROF_rapend  ('ATM_SurfaceFlux', 1)
    endif

    !########## Turbulence ##########
    if ( ATMOS_sw_phy_tb ) then
       call PROF_rapstart('ATM_Turbulence', 1)
       call ATMOS_PHY_TB_driver( update_flag = do_phy_tb )
       call PROF_rapend  ('ATM_Turbulence', 1)
    endif

    !########## Cumulus ##########
    if ( ATMOS_sw_phy_cp ) then
       call PROF_rapstart('ATM_Cumulus', 1)
       call ATMOS_PHY_CP_driver( update_flag = do_phy_cp )
       call PROF_rapend  ('ATM_Cumulus', 1)
    endif

    !########## Set Surface Boundary Condition ##########
    call PROF_rapstart('ATM_SfcExch', 2)
    call ATMOS_SURFACE_SET( countup = .true. )
    call PROF_rapend  ('ATM_SfcExch', 2)

    !########## History & Monitor ##########
    call PROF_rapstart('ATM_History', 1)
    call ATMOS_vars_history
    call ATMOS_vars_monitor
    call PROF_rapend  ('ATM_History', 1)

    return
  end subroutine ATMOS_driver

  !-----------------------------------------------------------------------------
  !> First send
  subroutine ATMOS_driver_firstsend
    use scale_atmos_boundary, only: &
       ATMOS_BOUNDARY_firstsend
    use mod_atmos_vars, only: &
       DENS, &
       MOMZ, &
       MOMX, &
       MOMY, &
       RHOT, &
       QTRC
    implicit none
    !---------------------------------------------------------------------------

    ! If this run is parent of online nesting, boundary data must be sent
    call PROF_rapstart('ATM_Boundary', 2)
    call ATMOS_BOUNDARY_firstsend( DENS, MOMZ, MOMX, MOMY, RHOT, QTRC ) ! [IN]
    call PROF_rapend  ('ATM_Boundary', 2)

    return
  end subroutine ATMOS_driver_firstsend

  !-----------------------------------------------------------------------------
  !> Finalize
  subroutine ATMOS_driver_finalize
    use scale_atmos_boundary, only: &
       ATMOS_BOUNDARY_UPDATE_FLAG, &
       ATMOS_BOUNDARY_finalize
    use scale_grid_nest, only: &
       NEST_COMM_disconnect
    use mod_atmos_vars, only: &
       DENS, &
       MOMZ, &
       MOMX, &
       MOMY, &
       RHOT, &
       QTRC
    implicit none
    !---------------------------------------------------------------------------

    !########## Lateral/Top Boundary Condition ###########
    if ( ATMOS_BOUNDARY_UPDATE_FLAG ) then
       ! If this run is parent of online nesting, boundary data must be sent
       call ATMOS_BOUNDARY_finalize

       ! Finialize Inter-Communicators
       call NEST_COMM_disconnect
    endif

    return
  end subroutine ATMOS_driver_finalize

  !-----------------------------------------------------------------------------
  !> Get surface boundary condition
  subroutine ATMOS_SURFACE_GET
    use mod_atmos_phy_sf_vars, only: &
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
       SFLX_GH    => ATMOS_PHY_SF_SFLX_GH,    &
       SFLX_QTRC  => ATMOS_PHY_SF_SFLX_QTRC,  &
       U10        => ATMOS_PHY_SF_U10,        &
       V10        => ATMOS_PHY_SF_V10,        &
       T2         => ATMOS_PHY_SF_T2,         &
       Q2         => ATMOS_PHY_SF_Q2
    use mod_cpl_admin, only: &
       CPL_sw
    use mod_cpl_vars, only: &
       CPL_getSFC_ATM
    implicit none
    !---------------------------------------------------------------------------

    if ( CPL_sw ) then
       call CPL_getSFC_ATM( SFC_TEMP  (:,:),   & ! [OUT]
                            SFC_albedo(:,:,:), & ! [OUT]
                            SFC_Z0M   (:,:),   & ! [OUT]
                            SFC_Z0H   (:,:),   & ! [OUT]
                            SFC_Z0E   (:,:),   & ! [OUT]
                            SFLX_MW   (:,:),   & ! [OUT]
                            SFLX_MU   (:,:),   & ! [OUT]
                            SFLX_MV   (:,:),   & ! [OUT]
                            SFLX_SH   (:,:),   & ! [OUT]
                            SFLX_LH   (:,:),   & ! [OUT]
                            SFLX_GH   (:,:),   & ! [OUT]
                            SFLX_QTRC (:,:,:), & ! [OUT]
                            U10       (:,:),   & ! [OUT]
                            V10       (:,:),   & ! [OUT]
                            T2        (:,:),   & ! [OUT]
                            Q2        (:,:)    ) ! [OUT]
    endif

    return
  end subroutine ATMOS_SURFACE_GET

  !-----------------------------------------------------------------------------
  !> Set surface boundary condition
  subroutine ATMOS_SURFACE_SET( countup )
    use scale_grid_real, only: &
       REAL_CZ, &
       REAL_Z1
    use scale_topography, only: &
       TOPO_Zsfc
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
       SFLX_SW_dn => ATMOS_PHY_RD_SFLX_SW_dn, &
       cosSZA     => ATMOS_PHY_RD_cosSZA
    use mod_cpl_admin, only: &
       CPL_sw
    use mod_cpl_vars, only: &
       CPL_putATM
    implicit none

    ! arguments
    logical, intent(in) :: countup

    ! works
    real(RP) :: SFC_DENS(IA,JA)
    real(RP) :: SFC_PRES(IA,JA)
    real(RP) :: ATM_PBL (IA,JA)
    !---------------------------------------------------------------------------

    if ( CPL_sw ) then
       ! planetary boundary layer
       ATM_PBL(:,:) = 100.0_RP ! tentative

       call BOTTOM_estimate( DENS     (:,:,:), & ! [IN]
                             PRES     (:,:,:), & ! [IN]
                             REAL_CZ  (:,:,:), & ! [IN]
                             TOPO_Zsfc(:,:),   & ! [IN]
                             REAL_Z1  (:,:),   & ! [IN]
                             SFC_DENS (:,:),   & ! [OUT]
                             SFC_PRES (:,:)    ) ! [OUT]

       call CPL_putATM( TEMP      (KS,:,:),   & ! [IN]
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
                        cosSZA    (:,:),      & ! [IN]
                        SFLX_rain (:,:),      & ! [IN]
                        SFLX_snow (:,:),      & ! [IN]
                        countup               ) ! [IN]
    endif

    return
  end subroutine ATMOS_SURFACE_SET

end module mod_atmos_driver
