!-------------------------------------------------------------------------------
!> module ATMOSPHERE driver
!!
!! @par Description
!!          Atmosphere module
!!
!! @author Team SCALE
!!
!<
!-------------------------------------------------------------------------------
#include "scalelib.h"
module mod_atmos_driver
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_io
  use scale_prof
  use scale_atmos_grid_cartesC_index
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: ATMOS_driver_tracer_setup
  public :: ATMOS_driver_setup
  public :: ATMOS_driver_finalize
  public :: ATMOS_driver_calc_tendency
  public :: ATMOS_driver_calc_tendency_from_sflux
  public :: ATMOS_driver_update
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
  !> Tracer setup
  subroutine ATMOS_driver_tracer_setup
    use mod_atmos_phy_mp_driver, only: &
       ATMOS_PHY_MP_driver_tracer_setup, &
       ATMOS_PHY_MP_USER_qhyd2qtrc,      &
       ATMOS_PHY_MP_driver_qhyd2qtrc_onlyqv
    use mod_atmos_phy_mp_vars, only: &
       QA_MP, &
       QS_MP, &
       QE_MP
    use mod_atmos_phy_ae_driver, only: &
       ATMOS_PHY_AE_driver_tracer_setup
    use mod_atmos_phy_ch_driver, only: &
       ATMOS_PHY_CH_driver_tracer_setup
    use mod_atmos_phy_tb_driver, only: &
       ATMOS_PHY_TB_driver_tracer_setup
    use mod_atmos_phy_bl_driver, only: &
       ATMOS_PHY_BL_driver_tracer_setup
    use scale_atmos_hydrometeor, only: &
       ATMOS_HYDROMETEOR_regist, &
       ATMOS_HYDROMETEOR_dry
    use mod_atmos_admin, only: &
       ATMOS_USE_QV
    use mod_atmos_phy_lt_driver, only: &
       ATMOS_PHY_LT_driver_tracer_setup
    implicit none

    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("ATMOS_driver_tracer_setup",*) 'Setup'


    call ATMOS_PHY_MP_driver_tracer_setup
    call ATMOS_PHY_CH_driver_tracer_setup
    call ATMOS_PHY_AE_driver_tracer_setup
    call ATMOS_PHY_TB_driver_tracer_setup
    call ATMOS_PHY_BL_driver_tracer_setup
    call ATMOS_PHY_LT_driver_tracer_setup

    if ( ATMOS_HYDROMETEOR_dry .and. ATMOS_USE_QV ) then
       LOG_INFO("ATMOS_driver_tracer_setup",*) 'Regist QV'
       call ATMOS_HYDROMETEOR_regist( 0, 0,                                    & ! (in)
            (/'QV'/),                                                          & ! (in)
            (/'Ratio of Water Vapor mass to total mass (Specific humidity)'/), & ! (in)
            (/'kg/kg'/),                                                       & ! (in)
            QS_MP                                                              ) ! (out)
       QA_MP = 1
       QE_MP = QS_MP
       ATMOS_PHY_MP_USER_qhyd2qtrc => ATMOS_PHY_MP_driver_qhyd2qtrc_onlyqv

    end if

    return
  end subroutine ATMOS_driver_tracer_setup

  !-----------------------------------------------------------------------------
  !> Setup
  subroutine ATMOS_driver_setup( init )
    use scale_time, only: &
       TIME_NOWDATE
    use scale_atmos_solarins, only: &
       ATMOS_SOLARINS_setup
    use scale_atmos_refstate, only: &
       ATMOS_REFSTATE_setup
    use mod_atmos_bnd_driver, only: &
       ATMOS_BOUNDARY_driver_setup
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
    use mod_atmos_phy_lt_driver, only: &
       ATMOS_PHY_LT_driver_setup
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
    logical, intent(in), optional :: init
    logical :: not_init
    !---------------------------------------------------------------------------

    not_init = .true.
    if ( present(init) ) not_init = .not. init

    LOG_NEWLINE
    LOG_INFO("ATMOS_driver_setup",*) 'Setup'

    LOG_NEWLINE
    LOG_INFO("ATMOS_driver_setup",*) 'Setup each atmospheric components ...'

    !--- setup solar insolation
    call ATMOS_SOLARINS_setup( BASE_LON, BASE_LAT, TIME_NOWDATE(1) )

    call PROF_rapstart('ATM_Refstate', 2)
    call ATMOS_REFSTATE_setup( KA, KS, KE, IA, ISB, IEB, JA, JSB, JEB, &
                               CZ(:), FZ(:), REAL_CZ(:,:,:), REAL_FZ(:,:,:), REAL_PHI(:,:,:) )
    call PROF_rapend  ('ATM_Refstate', 2)

    call PROF_rapstart('ATM_Boundary', 2)
    if ( not_init ) call ATMOS_BOUNDARY_driver_setup
    call PROF_rapend  ('ATM_Boundary', 2)

    ! setup each components
    if ( not_init ) call ATMOS_DYN_driver_setup
    if ( not_init ) call ATMOS_PHY_LT_driver_setup
                    call ATMOS_PHY_MP_driver_setup
    if ( not_init ) call ATMOS_PHY_AE_driver_setup
    if ( not_init ) call ATMOS_PHY_CH_driver_setup
    if ( not_init ) call ATMOS_PHY_RD_driver_setup
    if ( not_init ) call ATMOS_PHY_SF_driver_setup
    if ( not_init ) call ATMOS_PHY_TB_driver_setup
                    call ATMOS_PHY_BL_driver_setup
    if ( not_init ) call ATMOS_PHY_CP_driver_setup

    LOG_NEWLINE
    LOG_INFO("ATMOS_driver_setup",*) 'Finish setup of each atmospheric components.'

    return
  end subroutine ATMOS_driver_setup

  !-----------------------------------------------------------------------------
  !> Calculation tendency
  subroutine ATMOS_driver_calc_tendency( force )
    use mod_atmos_vars, only: &
       DENS_tp, &
       MOMZ_tp, &
       RHOU_tp, &
       RHOV_tp, &
       RHOT_tp, &
       RHOH_p,  &
       RHOQ_tp, &
       MOMX_tp, &
       MOMY_tp
    use mod_atmos_phy_mp_driver, only: &
       ATMOS_PHY_MP_driver_calc_tendency
    use mod_atmos_phy_ch_driver, only: &
       ATMOS_PHY_CH_driver_calc_tendency
    use mod_atmos_phy_ae_driver, only: &
       ATMOS_PHY_AE_driver_calc_tendency
    use mod_atmos_phy_rd_driver, only: &
       ATMOS_PHY_RD_driver_calc_tendency
    use mod_atmos_phy_sf_driver, only: &
       ATMOS_PHY_SF_driver_calc_tendency
    use mod_atmos_phy_tb_driver, only: &
       ATMOS_PHY_TB_driver_calc_tendency
    use mod_atmos_phy_cp_driver, only: &
       ATMOS_PHY_CP_driver_calc_tendency
    use mod_atmos_phy_bl_driver, only: &
       ATMOS_PHY_BL_driver_calc_tendency
    use mod_admin_time, only: &
       do_phy_mp => TIME_DOATMOS_PHY_MP, &
       do_phy_ae => TIME_DOATMOS_PHY_AE, &
       do_phy_ch => TIME_DOATMOS_PHY_CH, &
       do_phy_rd => TIME_DOATMOS_PHY_RD, &
       do_phy_sf => TIME_DOATMOS_PHY_SF, &
       do_phy_tb => TIME_DOATMOS_PHY_TB, &
       do_phy_bl => TIME_DOATMOS_PHY_BL, &
       do_phy_cp => TIME_DOATMOS_PHY_CP, &
       do_phy_lt => TIME_DOATMOS_PHY_LT
    use mod_atmos_admin, only: &
       ATMOS_sw_phy_mp, &
       ATMOS_sw_phy_ae, &
       ATMOS_sw_phy_ch, &
       ATMOS_sw_phy_rd, &
       ATMOS_sw_phy_sf, &
       ATMOS_sw_phy_tb, &
       ATMOS_sw_phy_bl, &
       ATMOS_sw_phy_cp
    use mod_cpl_admin, only: &
       CPL_sw
    implicit none
    logical, intent(in) :: force
    !---------------------------------------------------------------------------

    !########## Get Surface Boundary from coupler ##########
    call ATMOS_SURFACE_GET

    !########## calculate tendency ##########
    ! reset tendencies
    !$omp parallel workshare
    !$acc kernels
!OCL XFILL
    DENS_tp(:,:,:)   = 0.0_RP
    !$acc end kernels
    !$acc kernels
!OCL XFILL
    MOMZ_tp(:,:,:)   = 0.0_RP
    !$acc end kernels
    !$acc kernels
!OCL XFILL
    RHOU_tp(:,:,:)   = 0.0_RP
    !$acc end kernels
    !$acc kernels
!OCL XFILL
    RHOV_tp(:,:,:)   = 0.0_RP
    !$acc end kernels
    !$acc kernels
!OCL XFILL
    RHOT_tp(:,:,:)   = 0.0_RP
    !$acc end kernels
    !$acc kernels
!OCL XFILL
    RHOH_p (:,:,:)   = 0.0_RP
    !$acc end kernels
    !$acc kernels
!OCL XFILL
    RHOQ_tp(:,:,:,:) = 0.0_RP
    !$acc end kernels
    !$acc kernels
!OCL XFILL
    MOMX_tp(:,:,:)   = 0.0_RP
    !$acc end kernels
    !$acc kernels
!OCL XFILL
    MOMY_tp(:,:,:)   = 0.0_RP
    !$acc end kernels
    !$omp end parallel workshare

    ! Microphysics
    if ( ATMOS_sw_phy_mp ) then
       call PROF_rapstart('ATM_Microphysics', 1)
       call ATMOS_PHY_MP_driver_calc_tendency( update_flag = do_phy_mp .or. force )
       call PROF_rapend  ('ATM_Microphysics', 1)
    endif
    ! Aerosol
    if ( ATMOS_sw_phy_ae ) then
       call PROF_rapstart('ATM_Aerosol', 1)
       call ATMOS_PHY_AE_driver_calc_tendency( update_flag = do_phy_ae .or. force )
       call PROF_rapend  ('ATM_Aerosol', 1)
    endif
    ! Chemistry
    if ( ATMOS_sw_phy_ch ) then
       call PROF_rapstart('ATM_Chemistry', 1)
       call ATMOS_PHY_CH_driver_calc_tendency( update_flag = do_phy_ch .or. force )
       call PROF_rapend  ('ATM_Chemistry', 1)
    endif
    ! Radiation
    if ( ATMOS_sw_phy_rd ) then
       call PROF_rapstart('ATM_Radiation', 1)
       call ATMOS_PHY_RD_driver_calc_tendency( update_flag = do_phy_rd .or. force )
       call PROF_rapend  ('ATM_Radiation', 1)
    endif
    ! Turbulence
    if ( ATMOS_sw_phy_tb ) then
       call PROF_rapstart('ATM_Turbulence', 1)
       call ATMOS_PHY_TB_driver_calc_tendency( update_flag = do_phy_tb .or. force )
       call PROF_rapend  ('ATM_Turbulence', 1)
    endif
    ! Cumulus
    if ( ATMOS_sw_phy_cp ) then
       call PROF_rapstart('ATM_Cumulus', 1)
       call ATMOS_PHY_CP_driver_calc_tendency( update_flag = do_phy_cp .or. force )
       call PROF_rapend  ('ATM_Cumulus', 1)
    endif
    if ( .not. CPL_sw ) then
       ! Surface Flux
       if ( ATMOS_sw_phy_sf ) then
          call PROF_rapstart('ATM_SurfaceFlux', 1)
          call ATMOS_PHY_SF_driver_calc_tendency( update_flag = do_phy_sf .or. force )
          call PROF_rapend  ('ATM_SurfaceFlux', 1)
       endif
       ! Planetary Boundary layer
       if ( ATMOS_sw_phy_bl ) then
          call PROF_rapstart('ATM_PBL', 1)
          call ATMOS_PHY_BL_driver_calc_tendency( update_flag = do_phy_bl .or. force )
          call PROF_rapend  ('ATM_PBL', 1)
       endif
    end if

    !########## Set Surface Boundary Condition ##########
    call ATMOS_SURFACE_SET( countup = .true. )

    return
  end subroutine ATMOS_driver_calc_tendency

  !-----------------------------------------------------------------------------
  !> Calculation tendency from surface flux with coupler
  subroutine ATMOS_driver_calc_tendency_from_sflux( force )
    use mod_atmos_phy_sf_driver, only: &
       ATMOS_PHY_SF_driver_calc_tendency
    use mod_atmos_phy_bl_driver, only: &
       ATMOS_PHY_BL_driver_calc_tendency
    use mod_cpl_admin, only: &
       CPL_sw
    use mod_atmos_admin, only: &
       ATMOS_sw_phy_sf, &
       ATMOS_sw_phy_bl
    use mod_admin_time, only: &
       do_phy_sf => TIME_DOATMOS_PHY_SF, &
       do_phy_bl => TIME_DOATMOS_PHY_BL
    implicit none
    logical, intent(in) :: force
    !---------------------------------------------------------------------------

    if ( CPL_sw ) then

       !########## Get Surface Boundary Condition ##########
       call ATMOS_SURFACE_GET

       ! Surface Flux
       if ( ATMOS_sw_phy_sf ) then
          call PROF_rapstart('ATM_SurfaceFlux', 1)
          call ATMOS_PHY_SF_driver_calc_tendency( update_flag = do_phy_sf .or. force )
          call PROF_rapend  ('ATM_SurfaceFlux', 1)
       endif

       ! Planetary Boundary layer
       if ( ATMOS_sw_phy_bl ) then
          call PROF_rapstart('ATM_PBL', 1)
          call ATMOS_PHY_BL_driver_calc_tendency( update_flag = do_phy_bl .or. force )
          call PROF_rapend  ('ATM_PBL', 1)
       endif

    end if

    return
  end subroutine ATMOS_driver_calc_tendency_from_sflux

  !-----------------------------------------------------------------------------
  !> advance atmospheric state
  subroutine ATMOS_driver_update
    use mod_atmos_admin, only: &
       ATMOS_sw_dyn,    &
       ATMOS_sw_phy_mp, &
       ATMOS_sw_phy_ae, &
       ATMOS_sw_phy_lt
    use mod_admin_time, only: &
       do_dyn    => TIME_DOATMOS_DYN,    &
       do_phy_mp => TIME_DOATMOS_PHY_MP, &
       do_phy_ae => TIME_DOATMOS_PHY_AE
    use scale_atmos_refstate, only: &
       ATMOS_REFSTATE_update
    use mod_atmos_vars, only: &
       ATMOS_vars_calc_diagnostics,&
       ATMOS_vars_fillhalo,        &
       DENS,                       &
       TEMP,                       &
       PRES,                       &
       POTT,                       &
       QV
    use mod_atmos_bnd_driver, only: &
       ATMOS_BOUNDARY_driver_update, &
       ATMOS_BOUNDARY_driver_send,   &
       ATMOS_BOUNDARY_UPDATE_FLAG
    use mod_atmos_dyn_driver, only: &
       ATMOS_DYN_driver
    use mod_atmos_phy_mp_driver, only: &
       ATMOS_PHY_MP_driver_adjustment
    use mod_atmos_phy_ae_driver, only: &
       ATMOS_PHY_AE_driver_adjustment
    use mod_atmos_phy_lt_driver, only: &
       ATMOS_PHY_LT_driver_adjustment
    use scale_atmos_grid_cartesC, only: &
       CZ   => ATMOS_GRID_CARTESC_CZ,  &
       FZ   => ATMOS_GRID_CARTESC_FZ,  &
       FDZ  => ATMOS_GRID_CARTESC_FDZ, &
       RCDZ => ATMOS_GRID_CARTESC_RCDZ
    use scale_atmos_grid_cartesC_real, only: &
       REAL_CZ  => ATMOS_GRID_CARTESC_REAL_CZ, &
       REAL_FZ  => ATMOS_GRID_CARTESC_REAL_FZ, &
       REAL_PHI => ATMOS_GRID_CARTESC_REAL_PHI, &
       AREA     => ATMOS_GRID_CARTESC_REAL_AREA
    use scale_time, only: &
       TIME_NOWDAYSEC, &
       TIME_DTSEC
    implicit none

    !---------------------------------------------------------------------------

    !########## Dynamics ##########
    if ( ATMOS_sw_dyn ) then
       if ( ATMOS_BOUNDARY_UPDATE_FLAG ) then
          call PROF_rapstart('ATM_Boundary', 2)
          call ATMOS_BOUNDARY_driver_update( TIME_NOWDAYSEC - TIME_DTSEC * 0.5_DP )
          call PROF_rapend  ('ATM_Boundary', 2)
          call ATMOS_vars_fillhalo
       endif
       call PROF_rapstart('ATM_Dynamics', 1)
       call ATMOS_DYN_driver( do_dyn )
       call PROF_rapend  ('ATM_Dynamics', 1)
    endif


    !########## Calculate diagnostic variables ##########
    call ATMOS_vars_calc_diagnostics


    !########## Adjustment ##########
    ! Microphysics
    if ( ATMOS_sw_phy_mp ) then
       call PROF_rapstart('ATM_Microphysics', 1)
       call ATMOS_PHY_MP_driver_adjustment
       call PROF_rapend  ('ATM_Microphysics', 1)
       call ATMOS_vars_fillhalo
       call ATMOS_vars_calc_diagnostics
    endif
    ! Aerosol
    if ( ATMOS_sw_phy_ae ) then
       call PROF_rapstart('ATM_Aerosol', 1)
       call ATMOS_PHY_AE_driver_adjustment
       call PROF_rapend  ('ATM_Aerosol', 1)
       call ATMOS_vars_fillhalo
       call ATMOS_vars_calc_diagnostics
    endif
    ! Lightning
    if ( ATMOS_sw_phy_lt ) then
       call PROF_rapstart('ATM_Lightning', 1)
       call ATMOS_PHY_LT_driver_adjustment
       call PROF_rapend  ('ATM_Lightning', 1)
       call ATMOS_vars_fillhalo
       ! calc_diagnostics is not necessary
    end if


    !########## Lateral/Top Boundary Condition ###########
    if ( ATMOS_BOUNDARY_UPDATE_FLAG ) then
       call PROF_rapstart('ATM_Boundary', 2)
       call ATMOS_BOUNDARY_driver_update( TIME_NOWDAYSEC )
       call PROF_rapend  ('ATM_Boundary', 2)
       call ATMOS_vars_fillhalo
       call ATMOS_vars_calc_diagnostics
    endif


    !########## Send Lateral/Top Boundary Condition (Online nesting) ###########
    if ( ATMOS_BOUNDARY_UPDATE_FLAG ) then
       call PROF_rapstart('ATM_Boundary', 2)
       call ATMOS_BOUNDARY_driver_send
       call PROF_rapend  ('ATM_Boundary', 2)
    endif


    !########## Reference State ###########
    call PROF_rapstart('ATM_Refstate', 2)
    call ATMOS_REFSTATE_update( KA, KS, KE, IA, IS, IE, ISB, IEB, JA, JS, JE, JSB, JEB, &
                                DENS(:,:,:), POTT(:,:,:), TEMP(:,:,:), PRES(:,:,:), QV(:,:,:), & ! [IN]
                                CZ(:), FZ(:), FDZ(:), RCDZ(:),                                 & ! [IN]
                                REAL_CZ(:,:,:), REAL_FZ(:,:,:), REAL_PHI(:,:,:), AREA(:,:),    & ! [IN]
                                TIME_NOWDAYSEC                                                 ) ! [IN]
    call PROF_rapend  ('ATM_Refstate', 2)

    return
  end subroutine ATMOS_driver_update

  !-----------------------------------------------------------------------------
  !> Finalize
  subroutine ATMOS_driver_finalize( init )
    use scale_atmos_refstate, only: &
       ATMOS_REFSTATE_finalize
    use mod_atmos_bnd_driver, only: &
       ATMOS_BOUNDARY_driver_finalize
    use mod_atmos_dyn_driver, only: &
       ATMOS_DYN_driver_finalize
    use mod_atmos_phy_lt_driver, only: &
       ATMOS_PHY_LT_driver_finalize
    use mod_atmos_phy_mp_driver, only: &
       ATMOS_PHY_MP_driver_finalize
    use mod_atmos_phy_ae_driver, only: &
       ATMOS_PHY_AE_driver_finalize
    use mod_atmos_phy_ch_driver, only: &
       ATMOS_PHY_CH_driver_finalize
    use mod_atmos_phy_rd_driver, only: &
       ATMOS_PHY_RD_driver_finalize
    use mod_atmos_phy_tb_driver, only: &
       ATMOS_PHY_TB_driver_finalize
    use mod_atmos_phy_bl_driver, only: &
       ATMOS_PHY_BL_driver_finalize
    use mod_atmos_phy_cp_driver, only: &
       ATMOS_PHY_CP_driver_finalize
    implicit none
    logical, intent(in), optional :: init
    logical :: not_init
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("ATMOS_driver_finalize",*) 'Finalize'

    not_init = .true.
    if ( present(init) ) not_init = .not. init

    if ( not_init ) call ATMOS_DYN_driver_finalize
    if ( not_init ) call ATMOS_PHY_LT_driver_finalize
                    call ATMOS_PHY_MP_driver_finalize
    if ( not_init ) call ATMOS_PHY_AE_driver_finalize
    if ( not_init ) call ATMOS_PHY_CH_driver_finalize
    if ( not_init ) call ATMOS_PHY_RD_driver_finalize
    !if ( not_init ) call ATMOS_PHY_SF_driver_finalize
    if ( not_init ) call ATMOS_PHY_TB_driver_finalize
                    call ATMOS_PHY_BL_driver_finalize
    if ( not_init )call ATMOS_PHY_CP_driver_finalize

    if ( not_init ) call ATMOS_BOUNDARY_driver_finalize

    call ATMOS_REFSTATE_finalize

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
       SFLX_SHEX  => ATMOS_PHY_SF_SFLX_SHEX,  &
       SFLX_LHEX  => ATMOS_PHY_SF_SFLX_LHEX,  &
       SFLX_QVEX  => ATMOS_PHY_SF_SFLX_QVEX,  &
       SFLX_GH    => ATMOS_PHY_SF_SFLX_GH,    &
       SFLX_QTRC  => ATMOS_PHY_SF_SFLX_QTRC,  &
       SFLX_ENGI  => ATMOS_PHY_SF_SFLX_ENGI,  &
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

    call PROF_rapstart('ATM_SfcExch', 2)

    if ( CPL_sw ) then
       call CPL_getSFC_ATM( SFC_TEMP  (:,:),     & ! [OUT]
                            SFC_albedo(:,:,:,:), & ! [OUT]
                            SFC_Z0M   (:,:),     & ! [OUT]
                            SFC_Z0H   (:,:),     & ! [OUT]
                            SFC_Z0E   (:,:),     & ! [OUT]
                            SFLX_MW   (:,:),     & ! [OUT]
                            SFLX_MU   (:,:),     & ! [OUT]
                            SFLX_MV   (:,:),     & ! [OUT]
                            SFLX_SH   (:,:),     & ! [OUT]
                            SFLX_LH   (:,:),     & ! [OUT]
                            SFLX_SHEX (:,:),     & ! [OUT]
                            SFLX_LHEX (:,:),     & ! [OUT]
                            SFLX_QVEX (:,:),     & ! [OUT]
                            SFLX_GH   (:,:),     & ! [OUT]
                            SFLX_QTRC (:,:,:),   & ! [OUT]
                            SFLX_ENGI (:,:),     & ! [OUT]
                            U10       (:,:),     & ! [OUT]
                            V10       (:,:),     & ! [OUT]
                            T2        (:,:),     & ! [OUT]
                            Q2        (:,:)      ) ! [OUT]
    endif

    call PROF_rapend  ('ATM_SfcExch', 2)

    return
  end subroutine ATMOS_SURFACE_GET

  !-----------------------------------------------------------------------------
  !> Set surface boundary condition
  subroutine ATMOS_SURFACE_SET( countup )
    use scale_atmos_grid_cartesC_real, only: &
       REAL_FZ => ATMOS_GRID_CARTESC_REAL_FZ
    use scale_atmos_bottom, only: &
       BOTTOM_estimate => ATMOS_BOTTOM_estimate
    use mod_atmos_vars, only: &
       DENS, &
       QV,   &
       TEMP, &
       PRES, &
       W,    &
       U,    &
       V,    &
       PREC, &
       PREC_ENGI
    use mod_atmos_phy_sf_vars, only: &
       SFC_TEMP  => ATMOS_PHY_SF_SFC_TEMP
    use mod_atmos_phy_mp_vars, only: &
       SFLX_rain_MP => ATMOS_PHY_MP_SFLX_rain, &
       SFLX_snow_MP => ATMOS_PHY_MP_SFLX_snow, &
       SFLX_ENGI_MP => ATMOS_PHY_MP_SFLX_ENGI
    use mod_atmos_phy_cp_vars, only: &
       SFLX_rain_CP => ATMOS_PHY_CP_SFLX_rain, &
       SFLX_snow_CP => ATMOS_PHY_CP_SFLX_snow, &
       SFLX_ENGI_CP => ATMOS_PHY_CP_SFLX_ENGI
    use mod_atmos_phy_rd_vars, only: &
       SFLX_rad_dn => ATMOS_PHY_RD_SFLX_down, &
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

    real(RP) :: TEMP1(IA,JA)
    real(RP) :: PRES1(IA,JA)
    real(RP) :: W1   (IA,JA)
    real(RP) :: U1   (IA,JA)
    real(RP) :: V1   (IA,JA)
    real(RP) :: DENS1(IA,JA)
    real(RP) :: QV1  (IA,JA)

    integer  :: i,j
    !---------------------------------------------------------------------------

    call PROF_rapstart('ATM_SfcExch', 2)

    ! sum of rainfall from mp and cp
    !$omp parallel do private(i,j) OMP_SCHEDULE_
    !$acc kernels
    do j = JS, JE
    do i = IS, IE
       PREC     (i,j) = SFLX_rain_MP(i,j) + SFLX_rain_CP(i,j) + SFLX_snow_MP(i,j) + SFLX_snow_CP(i,j)
       PREC_ENGI(i,j) = SFLX_ENGI_MP(i,j) + SFLX_ENGI_CP(i,j)
    enddo
    enddo
    !$acc end kernels

    if ( CPL_sw ) then

       !$acc data create(SFC_DENS,SFC_PRES,TEMP1,PRES1,W1,U1,V1,DENS1,QV1)

       ! planetary boundary layer
       call BOTTOM_estimate( KA, KS, KE, IA, IS, IE, JA, JS, JE, &
                             DENS(:,:,:), PRES(:,:,:), QV(:,:,:), & ! [IN]
                             SFC_TEMP(:,:),                       & ! [IN]
                             REAL_FZ(:,:,:),                      & ! [IN]
                             SFC_DENS(:,:), SFC_PRES(:,:)         ) ! [OUT]

       !$omp parallel do
       !$acc kernels
       do j = JS, JE
       do i = IS, IE
          TEMP1(i,j) = TEMP(KS,i,j)
          PRES1(i,j) = PRES(KS,i,j)
          W1   (i,j) = W   (KS,i,j)
          U1   (i,j) = U   (KS,i,j)
          V1   (i,j) = V   (KS,i,j)
          DENS1(i,j) = DENS(KS,i,j)
          QV1  (i,j) = QV  (KS,i,j)
       end do
       end do
       !$acc end kernels

       call CPL_putATM( TEMP1      (:,:),     & ! [IN]
                        PRES1      (:,:),     & ! [IN]
                        W1         (:,:),     & ! [IN]
                        U1         (:,:),     & ! [IN]
                        V1         (:,:),     & ! [IN]
                        DENS1      (:,:),     & ! [IN]
                        QV1        (:,:),     & ! [IN]
                        ATM_PBL    (:,:),     & ! [IN]
                        SFC_DENS   (:,:),     & ! [IN]
                        SFC_PRES   (:,:),     & ! [IN]
                        SFLX_rad_dn(:,:,:,:), & ! [IN]
                        cosSZA     (:,:),     & ! [IN]
                        PREC       (:,:),     & ! [IN]
                        PREC_ENGI  (:,:),     & ! [IN]
                        countup               ) ! [IN]

       !$acc end data

    endif

    call PROF_rapend  ('ATM_SfcExch', 2)

    return
  end subroutine ATMOS_SURFACE_SET

  subroutine ATMOS_driver_boundary_update

  end subroutine ATMOS_driver_boundary_update
end module mod_atmos_driver
