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
#include "inc_openmp.h"
module mod_atmos_driver
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_stdio
  use scale_prof
  use scale_atmos_grid_cartesC_index
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: ATMOS_driver_config
  public :: ATMOS_driver_setup
  public :: ATMOS_driver_resume1
  public :: ATMOS_driver_resume2
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
  !> Config
  subroutine ATMOS_driver_config
    use mod_atmos_phy_mp_driver, only: &
       ATMOS_PHY_MP_driver_tracer_setup
    use mod_atmos_phy_ae_driver, only: &
       ATMOS_PHY_AE_driver_tracer_setup
    use mod_atmos_phy_ch_driver, only: &
       ATMOS_PHY_CH_driver_tracer_setup
    use mod_atmos_phy_tb_driver, only: &
       ATMOS_PHY_TB_driver_config
    use mod_atmos_phy_bl_driver, only: &
       ATMOS_PHY_BL_driver_tracer_setup
    implicit none
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '++++++ Module[CONFIG] / Categ[ATMOS] / Origin[SCALE-RM]'

    call ATMOS_PHY_MP_driver_tracer_setup
    call ATMOS_PHY_AE_driver_tracer_setup
    call ATMOS_PHY_CH_driver_tracer_setup
    call ATMOS_PHY_TB_driver_config
    call ATMOS_PHY_BL_driver_tracer_setup

    return
  end subroutine ATMOS_driver_config

  !-----------------------------------------------------------------------------
  !> Setup
  subroutine ATMOS_driver_setup
    use scale_time, only: &
       TIME_NOWDATE
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
    use mod_atmos_phy_mp_vars, only: &
       QA_MP, &
       QS_MP, &
       QE_MP
    use mod_atmos_phy_ch_driver, only: &
       ATMOS_PHY_CH_driver_setup
    use mod_atmos_phy_ae_driver, only: &
       ATMOS_PHY_AE_driver_setup
    use mod_atmos_phy_rd_driver, only: &
       ATMOS_PHY_RD_driver_setup
    use mod_atmos_phy_sf_driver, only: &
       ATMOS_PHY_SF_driver_setup
    use mod_atmos_phy_tb_driver, only: &
       ATMOS_PHY_TB_driver_setup
    use mod_atmos_phy_bl_driver, only: &
       ATMOS_PHY_BL_driver_setup
    use mod_atmos_phy_cp_driver, only: &
       ATMOS_PHY_CP_driver_setup
    use scale_atmos_grid_cartesC, only: &
       CZ => ATMOS_GRID_CARTESC_CZ, &
       FZ => ATMOS_GRID_CARTESC_FZ
    use scale_atmos_grid_cartesC_real, only: &
       BASE_LON => ATMOS_GRID_CARTESC_REAL_BASEPOINT_LON, &
       BASE_LAT => ATMOS_GRID_CARTESC_REAL_BASEPOINT_LAT, &
       REAL_CZ  => ATMOS_GRID_CARTESC_REAL_CZ, &
       REAL_FZ  => ATMOS_GRID_CARTESC_REAL_FZ, &
       REAL_PHI => ATMOS_GRID_CARTESC_REAL_PHI
    implicit none
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '++++++ Module[DRIVER] / Categ[ATMOS] / Origin[SCALE-RM]'

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** Setup each atmospheric components ...'

    !--- setup solar insolation
    call ATMOS_SOLARINS_setup( BASE_LON, BASE_LAT, TIME_NOWDATE(1) )

    call PROF_rapstart('ATM_Refstate', 2)
    call ATMOS_REFSTATE_setup( KA, KS, KE, IA, ISB, IEB, JA, JSB, JEB, &
                               CZ(:), FZ(:), REAL_CZ(:,:,:), REAL_FZ(:,:,:), REAL_PHI(:,:,:) )

    call PROF_rapend  ('ATM_Refstate', 2)

    call PROF_rapstart('ATM_Boundary', 2)
    call ATMOS_BOUNDARY_setup( QA_MP, QS_MP, QE_MP )
    call PROF_rapend  ('ATM_Boundary', 2)

    ! setup each components
    call ATMOS_DYN_driver_setup
    call ATMOS_PHY_MP_driver_setup
    call ATMOS_PHY_AE_driver_setup
    call ATMOS_PHY_CH_driver_setup
    call ATMOS_PHY_RD_driver_setup
    call ATMOS_PHY_SF_driver_setup
    call ATMOS_PHY_TB_driver_setup
    call ATMOS_PHY_BL_driver_setup
    call ATMOS_PHY_CP_driver_setup

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** Finish setup of each atmospheric components.'

    return
  end subroutine ATMOS_driver_setup

  !-----------------------------------------------------------------------------
  !> Resume
  subroutine ATMOS_driver_resume1
    use mod_atmos_vars, only: &
       ATMOS_vars_calc_diagnostics,     &
       ATMOS_vars_history_setpres, &
       DENS,                       &
       MOMZ,                       &
       MOMX,                       &
       MOMY,                       &
       RHOT,                       &
       QTRC,                       &
       QV,                         &
       POTT,                       &
       TEMP,                       &
       PRES,                       &
       DENS_tp,                    &
       MOMZ_tp,                    &
       RHOU_tp,                    &
       RHOV_tp,                    &
       RHOT_tp,                    &
       RHOH_p,                     &
       RHOQ_tp,                    &
       MOMX_tp,                    &
       MOMY_tp
    use mod_atmos_phy_mp_vars, only: &
       QS_MP, &
       QE_MP
    use scale_atmos_refstate, only: &
       ATMOS_REFSTATE_resume
    use scale_atmos_boundary, only: &
       ATMOS_BOUNDARY_resume
    use mod_atmos_phy_mp_driver, only: &
       ATMOS_PHY_MP_driver_resume
    use mod_atmos_phy_ch_driver, only: &
       ATMOS_PHY_CH_driver_resume
    use mod_atmos_phy_ae_driver, only: &
       ATMOS_PHY_AE_driver_resume
    use mod_atmos_phy_rd_driver, only: &
       ATMOS_PHY_RD_driver_resume
    use scale_atmos_grid_cartesC, only: &
       CZ   => ATMOS_GRID_CARTESC_CZ,  &
       FZ   => ATMOS_GRID_CARTESC_FZ,  &
       FDZ  => ATMOS_GRID_CARTESC_FDZ, &
       RCDZ => ATMOS_GRID_CARTESC_RCDZ
    use scale_atmos_grid_cartesC_real, only: &
       REAL_CZ  => ATMOS_GRID_CARTESC_REAL_CZ, &
       REAL_FZ  => ATMOS_GRID_CARTESC_REAL_FZ, &
       REAL_PHI => ATMOS_GRID_CARTESC_REAL_PHI
    use scale_time, only: &
       TIME_NOWSEC
    implicit none
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '++++++ Module[DRIVER] / Categ[ATMOS] / Origin[SCALE-RM]'

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** Resume each atmospheric components 1 ...'

    call PROF_rapstart('ATM_Boundary', 2)
    call ATMOS_BOUNDARY_resume( DENS, MOMZ, MOMX, MOMY, RHOT, QTRC(:,:,:,QS_MP:QE_MP) )
    call PROF_rapend  ('ATM_Boundary', 2)

    !########## Get Surface Boundary Condition ##########
    call PROF_rapstart('ATM_SfcExch', 2)
    call ATMOS_SURFACE_GET
    call PROF_rapend  ('ATM_SfcExch', 2)

    !########## initialize tendencies ##########
!OCL XFILL
    DENS_tp(:,:,:)   = 0.0_RP
!OCL XFILL
    MOMZ_tp(:,:,:)   = 0.0_RP
!OCL XFILL
    RHOU_tp(:,:,:)   = 0.0_RP
!OCL XFILL
    RHOV_tp(:,:,:)   = 0.0_RP
!OCL XFILL
    RHOT_tp(:,:,:)   = 0.0_RP
!OCL XFILL
    RHOH_p (:,:,:)   = 0.0_RP
!OCL XFILL
    RHOQ_tp(:,:,:,:) = 0.0_RP
!OCL XFILL
    MOMX_tp(:,:,:)   = 0.0_RP
!OCL XFILL
    MOMY_tp(:,:,:)   = 0.0_RP

    ! setup each components
    call ATMOS_PHY_MP_driver_resume
    call ATMOS_PHY_AE_driver_resume
    call ATMOS_PHY_CH_driver_resume
    call ATMOS_PHY_RD_driver_resume

    !########## Calculate diagnostic variables ##########
    call ATMOS_vars_calc_diagnostics
    call ATMOS_vars_history_setpres

    !########## Set reference state ##########
    call PROF_rapstart('ATM_Refstate', 2)
    call ATMOS_REFSTATE_resume( KA, KS, KE, IA, ISB, IEB, JA, JSB, JEB, &
                                DENS(:,:,:), POTT(:,:,:), TEMP(:,:,:), PRES(:,:,:), QV(:,:,:), &
                                CZ(:), FZ(:), FDZ(:), RCDZ(:),                                 &
                                REAL_CZ(:,:,:), REAL_FZ(:,:,:), REAL_PHI(:,:,:),               &
                                TIME_NOWSEC                                                    )

    call PROF_rapend  ('ATM_Refstate', 2)

    !########## Set Surface Boundary Condition ##########
    call PROF_rapstart('ATM_SfcExch', 2)
    call ATMOS_SURFACE_SET( countup = .false. )
    call PROF_rapend  ('ATM_SfcExch', 2)

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** Finish resume of each atmospheric components 1.'

    return
  end subroutine ATMOS_driver_resume1

  !-----------------------------------------------------------------------------
  !> Setup
  subroutine ATMOS_driver_resume2
    use mod_atmos_vars, only: &
       ATMOS_vars_history, &
       ATMOS_vars_monitor
    use mod_atmos_phy_sf_driver, only: &
       ATMOS_PHY_SF_driver_resume
    use mod_atmos_phy_tb_driver, only: &
       ATMOS_PHY_TB_driver_resume
    use mod_atmos_phy_cp_driver, only: &
       ATMOS_PHY_CP_driver_resume
    use mod_atmos_phy_bl_driver, only: &
       ATMOS_PHY_BL_driver_resume
    implicit none
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '++++++ Module[DRIVER] / Categ[ATMOS] / Origin[SCALE-RM]'

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** Resume each atmospheric components 2 ...'

    !########## Get Surface Boundary Condition ##########
    call PROF_rapstart('ATM_SfcExch', 2)
    call ATMOS_SURFACE_GET
    call PROF_rapend  ('ATM_SfcExch', 2)

    ! setup each components
    call ATMOS_PHY_SF_driver_resume
    call ATMOS_PHY_TB_driver_resume
    call ATMOS_PHY_BL_driver_resume
    call ATMOS_PHY_CP_driver_resume

    !########## Set Surface Boundary Condition ##########
    call PROF_rapstart('ATM_SfcExch', 2)
    call ATMOS_SURFACE_SET( countup = .true. )
    call PROF_rapend  ('ATM_SfcExch', 2)

    !########## History & Monitor ##########
    call PROF_rapstart('ATM_History', 1)
    call ATMOS_vars_history
    call ATMOS_vars_monitor
    call PROF_rapend  ('ATM_History', 1)

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** Finish resume of each atmospheric components 2.'

    return
  end subroutine ATMOS_driver_resume2

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
       do_phy_bl => TIME_DOATMOS_PHY_BL, &
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
       ATMOS_sw_phy_bl, &
       ATMOS_sw_phy_cp
    use mod_atmos_vars, only: &
       ATMOS_vars_history,         &
       ATMOS_vars_calc_diagnostics,&
       ATMOS_vars_history_setpres, &
       ATMOS_vars_monitor,         &
       DENS,                       &
       MOMZ,                       &
       MOMX,                       &
       MOMY,                       &
       RHOT,                       &
       QTRC,                       &
       QV,                         &
       PRES,                       &
       POTT,                       &
       TEMP,                       &
       DENS_tp,                    &
       MOMZ_tp,                    &
       RHOU_tp,                    &
       RHOV_tp,                    &
       RHOT_tp,                    &
       RHOH_p,                     &
       RHOQ_tp,                    &
       MOMX_tp,                    &
       MOMY_tp
    use mod_atmos_dyn_driver, only: &
       ATMOS_DYN_driver
    use mod_atmos_phy_mp_driver, only: &
       ATMOS_PHY_MP_driver_calc_tendency, &
       ATMOS_PHY_MP_driver_adjustment
    use mod_atmos_phy_ae_driver, only: &
       ATMOS_PHY_AE_driver_calc_tendency, &
       ATMOS_PHY_AE_driver_adjustment
    use mod_atmos_phy_ch_driver, only: &
       ATMOS_PHY_CH_driver_calc_tendency
    use mod_atmos_phy_rd_driver, only: &
       ATMOS_PHY_RD_driver
    use mod_atmos_phy_sf_driver, only: &
       ATMOS_PHY_SF_driver
    use mod_atmos_phy_tb_driver, only: &
       ATMOS_PHY_TB_driver
    use mod_atmos_phy_bl_driver, only: &
       ATMOS_PHY_BL_driver_calc_tendency
    use mod_atmos_phy_cp_driver, only: &
       ATMOS_PHY_CP_driver
    use scale_atmos_grid_cartesC, only: &
       CZ   => ATMOS_GRID_CARTESC_CZ,  &
       FZ   => ATMOS_GRID_CARTESC_FZ,  &
       FDZ  => ATMOS_GRID_CARTESC_FDZ, &
       RCDZ => ATMOS_GRID_CARTESC_RCDZ
    use scale_atmos_grid_cartesC_real, only: &
       REAL_CZ  => ATMOS_GRID_CARTESC_REAL_CZ, &
       REAL_FZ  => ATMOS_GRID_CARTESC_REAL_FZ, &
       REAL_PHI => ATMOS_GRID_CARTESC_REAL_PHI
    use scale_time, only: &
       TIME_NOWSEC
    implicit none
    !---------------------------------------------------------------------------

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

    !########## Calculate diagnostic variables ##########
    call ATMOS_vars_calc_diagnostics


    !########## Adjustment ##########
    ! Microphysics
    if ( ATMOS_sw_phy_mp .and. do_phy_mp ) then
       call PROF_rapstart('ATM_Microphysics', 1)
       call ATMOS_PHY_MP_driver_adjustment
       call PROF_rapend  ('ATM_Microphysics', 1)
       call ATMOS_vars_calc_diagnostics
    endif
    ! Aerosol
    if ( ATMOS_sw_phy_ae .and. do_phy_ae ) then
       call PROF_rapstart('ATM_Aerosol', 1)
       call ATMOS_PHY_AE_driver_adjustment
       call PROF_rapend  ('ATM_Aerosol', 1)
       call ATMOS_vars_calc_diagnostics
    endif

    !########## Reference State ###########
    if ( ATMOS_REFSTATE_UPDATE_FLAG ) then
       call PROF_rapstart('ATM_Refstate', 2)
       call ATMOS_REFSTATE_update( KA, KS, KE, IA, ISB, IEB, JA, JSB, JEB, &
                                   DENS(:,:,:), POTT(:,:,:), TEMP(:,:,:), PRES(:,:,:), QV(:,:,:), & ! [IN]
                                   CZ(:), FZ(:), FDZ(:), RCDZ(:),                                 & ! [IN]
                                   REAL_CZ(:,:,:), REAL_FZ(:,:,:), REAL_PHI(:,:,:),               & ! [IN]
                                   TIME_NOWSEC                                                    ) ! [IN]
       call PROF_rapend  ('ATM_Refstate', 2)
    endif


    !########## Set hydrostatic pressure coordinate ##########
    call PROF_rapstart('ATM_History', 1)
    call ATMOS_vars_history_setpres
    call PROF_rapend  ('ATM_History', 1)


    !########## calculate tendency ##########
    ! reset tendencies
!OCL XFILL
    DENS_tp(:,:,:)   = 0.0_RP
!OCL XFILL
    MOMZ_tp(:,:,:)   = 0.0_RP
!OCL XFILL
    RHOU_tp(:,:,:)   = 0.0_RP
!OCL XFILL
    RHOV_tp(:,:,:)   = 0.0_RP
!OCL XFILL
    RHOT_tp(:,:,:)   = 0.0_RP
!OCL XFILL
    RHOH_p (:,:,:)   = 0.0_RP
!OCL XFILL
    RHOQ_tp(:,:,:,:) = 0.0_RP
!OCL XFILL
    MOMX_tp(:,:,:)   = 0.0_RP
!OCL XFILL
    MOMY_tp(:,:,:)   = 0.0_RP

    ! Microphysics
    if ( ATMOS_sw_phy_mp ) then
       call PROF_rapstart('ATM_Microphysics', 1)
       call ATMOS_PHY_MP_driver_calc_tendency( update_flag = do_phy_mp )
       call PROF_rapend  ('ATM_Microphysics', 1)
    endif
    ! Aerosol
    if ( ATMOS_sw_phy_ae ) then
       call PROF_rapstart('ATM_Aerosol', 1)
       call ATMOS_PHY_AE_driver_calc_tendency( update_flag = do_phy_ae )
       call PROF_rapend  ('ATM_Aerosol', 1)
    endif
    ! Chemistry
    if ( ATMOS_sw_phy_ch ) then
       call PROF_rapstart('ATM_Chemistry', 1)
       call ATMOS_PHY_CH_driver_calc_tendency( update_flag = do_phy_ch )
       call PROF_rapend  ('ATM_Chemistry', 1)
    endif
    ! Radiation
    if ( ATMOS_sw_phy_rd ) then
       call PROF_rapstart('ATM_Radiation', 1)
       call ATMOS_PHY_RD_driver( update_flag = do_phy_rd )
       call PROF_rapend  ('ATM_Radiation', 1)
    endif
    ! Surface Flux
    if ( ATMOS_sw_phy_sf ) then
       call PROF_rapstart('ATM_SurfaceFlux', 1)
       call ATMOS_PHY_SF_driver( update_flag = do_phy_sf )
       call PROF_rapend  ('ATM_SurfaceFlux', 1)
    endif
    ! Turbulence
    if ( ATMOS_sw_phy_tb ) then
       call PROF_rapstart('ATM_Turbulence', 1)
       call ATMOS_PHY_TB_driver( update_flag = do_phy_tb )
       call PROF_rapend  ('ATM_Turbulence', 1)
    endif
    ! Pnaletary Boundary layer
    if ( ATMOS_sw_phy_bl ) then
       call PROF_rapstart('ATM_PBL', 1)
       call ATMOS_PHY_BL_driver_calc_tendency( update_flag = do_phy_bl )
       call PROF_rapend  ('ATM_PBL', 1)
    endif
    ! Cumulus
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
  !> Finalize
  subroutine ATMOS_driver_finalize
    use scale_atmos_boundary, only: &
       ATMOS_BOUNDARY_UPDATE_FLAG, &
       ATMOS_BOUNDARY_finalize
    use scale_comm_cartesC_nest, only: &
       NEST_COMM_disconnect => COMM_CARTESC_NEST_disconnect
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
    use scale_atmos_grid_cartesC_real, only: &
       REAL_CZ => ATMOS_GRID_CARTESC_REAL_CZ, &
       REAL_Z1 => ATMOS_GRID_CARTESC_REAL_Z1
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
       SFLX_rain_MP => ATMOS_PHY_MP_SFLX_rain, &
       SFLX_snow_MP => ATMOS_PHY_MP_SFLX_snow
    use mod_atmos_phy_cp_vars, only: &
       SFLX_rain_CP => ATMOS_PHY_CP_SFLX_rain
    use mod_atmos_phy_rd_vars, only: &
       SFLX_rad_dn => ATMOS_PHY_RD_SFLX_downall, &
       cosSZA      => ATMOS_PHY_RD_cosSZA
    use mod_atmos_phy_bl_vars, only: &
       ATM_PBL => ATMOS_PHY_BL_Zi
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

    real(RP) :: SFLX_rain(IA,JA)
    real(RP) :: SFLX_snow(IA,JA)

    integer  :: i,j
    !---------------------------------------------------------------------------

    if ( CPL_sw ) then
       ! sum of rainfall from mp and cp
       !$omp parallel do private(i,j) OMP_SCHEDULE_
       do j = 1, JA
       do i = 1, IA
          SFLX_rain(i,j) = SFLX_rain_MP(i,j) + SFLX_rain_CP(i,j)
          SFLX_snow(i,j) = SFLX_snow_MP(i,j)
       enddo
       enddo

       ! planetary boundary layer
       call BOTTOM_estimate( KA, KS, KE, IA, ISB, IEB, JA, JSB, JEB, &
                             DENS(:,:,:), PRES(:,:,:),                     & ! [IN]
                             REAL_CZ(:,:,:), TOPO_Zsfc(:,:), REAL_Z1(:,:), & ! [IN]
                             SFC_DENS(:,:), SFC_PRES(:,:)                  ) ! [OUT]

       call CPL_putATM( TEMP       (KS,:,:),   & ! [IN]
                        PRES       (KS,:,:),   & ! [IN]
                        W          (KS,:,:),   & ! [IN]
                        U          (KS,:,:),   & ! [IN]
                        V          (KS,:,:),   & ! [IN]
                        DENS       (KS,:,:),   & ! [IN]
                        QTRC       (KS,:,:,:), & ! [IN]
                        ATM_PBL    (:,:),      & ! [IN]
                        SFC_DENS   (:,:),      & ! [IN]
                        SFC_PRES   (:,:),      & ! [IN]
                        SFLX_rad_dn(:,:,:,:),  & ! [IN]
                        cosSZA     (:,:),      & ! [IN]
                        SFLX_rain  (:,:),      & ! [IN]
                        SFLX_snow  (:,:),      & ! [IN]
                        countup                ) ! [IN]
    endif

    return
  end subroutine ATMOS_SURFACE_SET

  subroutine ATMOS_driver_boundary_update

  end subroutine ATMOS_driver_boundary_update
end module mod_atmos_driver
