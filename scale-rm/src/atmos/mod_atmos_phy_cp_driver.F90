!-------------------------------------------------------------------------------
!> module ATMOSPHERE / Physics Cumulus
!!
!! @par Description
!!          Cumulus parameterization driver
!!
!! @author Team SCALE
!!
!<
!-------------------------------------------------------------------------------
#include "scalelib.h"
module mod_atmos_phy_cp_driver
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_io
  use scale_prof
  use scale_atmos_grid_cartesC_index
  use scale_tracer
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: ATMOS_PHY_CP_driver_setup
  public :: ATMOS_PHY_CP_driver_calc_tendency

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
  subroutine ATMOS_PHY_CP_driver_setup
    use scale_prc, only: &
       PRC_abort
    use scale_atmos_phy_cp_common, only: &
       ATMOS_PHY_CP_common_setup
    use scale_atmos_phy_cp_kf, only: &
       ATMOS_PHY_CP_kf_setup
    use mod_atmos_admin, only: &
       ATMOS_PHY_CP_TYPE, &
       ATMOS_sw_phy_cp
    use scale_time , only :&
       TIME_DTSEC,             &
       TIME_DTSEC_ATMOS_PHY_CP
    use scale_atmos_grid_cartesC_real, only: &
       ATMOS_GRID_CARTESC_REAL_CZ, &
       ATMOS_GRID_CARTESC_REAL_AREA
    use scale_atmos_hydrometeor, only: &
       ATMOS_HYDROMETEOR_ice_phase
    implicit none

    logical :: warmrain
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("ATMOS_PHY_CP_driver_setup",*) 'Setup'

    if ( ATMOS_sw_phy_cp ) then

       ! setup library component
       call ATMOS_PHY_CP_common_setup
       select case ( ATMOS_PHY_CP_TYPE )
       case ( 'KF' )
          warmrain = ( .not. ATMOS_HYDROMETEOR_ice_phase )
          call ATMOS_PHY_CP_kf_setup( KA, KS, KE, IA, IS, IE, JA, JS, JE, &
                                      ATMOS_GRID_CARTESC_REAL_CZ, ATMOS_GRID_CARTESC_REAL_AREA, &
                                      warmrain                                                  )
       case default
          LOG_ERROR("ATMOS_PHY_CP_driver_setup",*) 'ATMOS_PHY_CP_TYPE (', trim(ATMOS_PHY_CP_TYPE), ') is invalid. Check!'
          call PRC_abort
       end select

    else
       LOG_INFO("ATMOS_PHY_CP_driver_setup",*) 'this component is never called.'
    endif

    return
  end subroutine ATMOS_PHY_CP_driver_setup

  !-----------------------------------------------------------------------------
  !> Driver
  subroutine ATMOS_PHY_CP_driver_calc_tendency( update_flag )
    use scale_statistics, only: &
       STATISTICS_checktotal, &
       STATISTICS_total
    use scale_atmos_grid_cartesC_real, only: &
       ATMOS_GRID_CARTESC_REAL_VOL,       &
       ATMOS_GRID_CARTESC_REAL_TOTVOL,    &
       ATMOS_GRID_CARTESC_REAL_VOLWXY,    &
       ATMOS_GRID_CARTESC_REAL_TOTVOLWXY, &
       ATMOS_GRID_CARTESC_REAL_VOLZUY,    &
       ATMOS_GRID_CARTESC_REAL_TOTVOLZUY, &
       ATMOS_GRID_CARTESC_REAL_VOLZXV,    &
       ATMOS_GRID_CARTESC_REAL_TOTVOLZXV
    use scale_file_history, only: &
       FILE_HISTORY_in
    use scale_time , only :&
       TIME_DTSEC, &
       TIME_DTSEC_ATMOS_PHY_CP
    use scale_atmos_grid_cartesC_real, only: &
       FZ => ATMOS_GRID_CARTESC_REAL_FZ
    use scale_atmos_hydrometeor, only: &
       HYD_NAME, &
       CV_WATER
    use scale_atmos_phy_cp_kf, only: &
       ATMOS_PHY_CP_kf_tendency
    use scale_atmos_phy_cp_common, only: &
       ATMOS_PHY_CP_common_wmean
    use mod_atmos_phy_mp_vars, only: &
       QS_MP, &
       QE_MP
    use mod_atmos_phy_mp_driver, only: &
       ATMOS_PHY_MP_driver_qhyd2qtrc
    use mod_atmos_admin, only: &
       ATMOS_PHY_CP_TYPE
    use mod_atmos_vars, only: &
       DENS   => DENS_av, &
       MOMZ   => MOMZ_av, &
       MOMX   => MOMX_av, &
       MOMY   => MOMY_av, &
       RHOT   => RHOT_av, &
       QTRC   => QTRC_av, &
       DENS_t => DENS_tp, &
       MOMZ_t => MOMZ_tp, &
       MOMX_t => MOMX_tp, &
       MOMY_t => MOMY_tp, &
       RHOT_t => RHOT_tp, &
       RHOQ_t => RHOQ_tp, &
       U,                 &
       V,                 &
       W,                 &
       TEMP,              &
       PRES,              &
       QDRY,              &
       QV,                &
       Rtot,              &
       CPtot
    use scale_atmos_hydrometeor, only: &
       N_HYD
    use mod_atmos_phy_cp_vars, only: &
       DENS_t_CP      => ATMOS_PHY_CP_DENS_t,         &
       RHOT_t_CP      => ATMOS_PHY_CP_RHOT_t,         &
       RHOQV_t_CP     => ATMOS_PHY_CP_RHOQV_t,        &
       RHOHYD_t_CP    => ATMOS_PHY_CP_RHOHYD_t,       &
       MFLX_cloudbase => ATMOS_PHY_CP_MFLX_cloudbase, &
       SFLX_rain      => ATMOS_PHY_CP_SFLX_rain,      &  ! convective rain [kg/m2/s]
       SFLX_ENGI      => ATMOS_PHY_CP_SFLX_ENGI,      &  ! internal energy flux [J/m2/s]
       cloudtop       => ATMOS_PHY_CP_cloudtop,       &  ! cloud top height [m]
       cloudbase      => ATMOS_PHY_CP_cloudbase,      &  ! cloud base height [m]
       cldfrac_dp     => ATMOS_PHY_CP_cldfrac_dp,     &  ! cloud fraction (deep convection) (0-1)
       cldfrac_sh     => ATMOS_PHY_CP_cldfrac_sh,     &  ! cloud fraction (shallow convection) (0-1)
       w0mean         => ATMOS_PHY_CP_w0mean,         &  ! running mean vertical wind velocity [m/s]
       kf_nca         => ATMOS_PHY_CP_kf_nca             ! advection/cumulus convection timescale/dt for KF [step]
    implicit none

    logical, intent(in) :: update_flag

    real(RP) :: RHOQ_t_CP(KA,IA,JA,QS_MP:QE_MP)

    integer  :: k, i, j, iq
    !---------------------------------------------------------------------------

    ! temporal running mean of vertical velocity
    call ATMOS_PHY_CP_common_wmean( KA, KS, KE, IA, IS, IE, JA, JS, JE, &
                                    W(:,:,:),                            & ! [IN]
                                    TIME_DTSEC, TIME_DTSEC_ATMOS_PHY_CP, & ! [IN]
                                    w0mean(:,:,:)                        ) ! [INOUT]
    call FILE_HISTORY_in( w0mean(:,:,:), 'w0mean', 'running mean vertical wind velocity', 'kg/m2/s', fill_halo=.true. )

    if ( update_flag ) then ! update
       select case ( ATMOS_PHY_CP_TYPE )
       case ( 'KF' )
          call ATMOS_PHY_CP_kf_tendency( KA, KS, KE, IA, IS, IE, JA, JS, JE, &
                                         DENS(:,:,:),                              & ! [IN]
                                         U(:,:,:), V(:,:,:),                       & ! [IN]
                                         RHOT(:,:,:), TEMP(:,:,:), PRES(:,:,:),    & ! [IN]
                                         QDRY(:,:,:), QV(:,:,:),                   & ! [IN]
                                         Rtot(:,:,:), CPtot(:,:,:),                & ! [IN]
                                         w0mean(:,:,:),                            & ! [IN]
                                         FZ,                                       & ! [IN]
                                         TIME_DTSEC_ATMOS_PHY_CP,                  & ! [IN]
                                         DENS_t_CP(:,:,:),                         & ! [INOUT]
                                         RHOT_t_CP(:,:,:),                         & ! [INOUT]
                                         RHOQV_t_CP(:,:,:), RHOHYD_t_CP(:,:,:,:),  & ! [INOUT]
                                         kf_nca(:,:),                              & ! [INOUT]
                                         SFLX_rain(:,:),                           & ! [OUT]
                                         cloudtop(:,:), cloudbase(:,:),            & ! [OUT]
                                         cldfrac_dp(:,:,:), cldfrac_sh(:,:,:)      ) ! [OUT]
          !$omp parallel do
          do j = JS, JE
          do i = IS, IE
             SFLX_ENGI(i,j) = SFLX_rain(i,j) * TEMP(KS,i,j) * CV_WATER
          end do
          end do

       end select

!OCL XFILL
       do j  = JSB, JEB
       do i  = ISB, IEB
          MFLX_cloudbase(i,j) = 0.0_RP
       enddo
       enddo

       ! diagnose tendency of number concentration

       call FILE_HISTORY_in( MFLX_cloudbase(:,:),   'CBMFX',     'cloud base mass flux',             'kg/m2/s', fill_halo=.true. )
       call FILE_HISTORY_in( SFLX_rain     (:,:),   'RAIN_CP',   'surface rain rate by CP',          'kg/m2/s', fill_halo=.true. )
       call FILE_HISTORY_in( SFLX_rain     (:,:),   'PREC_CP',   'surface precipitation rate by CP', 'kg/m2/s', fill_halo=.true. )
       call FILE_HISTORY_in( cloudtop      (:,:),   'CUMHGT',    'CP cloud top height',              'm',       fill_halo=.true. )
       call FILE_HISTORY_in( cloudbase     (:,:),   'CUBASE',    'CP cloud base height',             'm',       fill_halo=.true. )
       call FILE_HISTORY_in( cldfrac_dp    (:,:,:), 'CUMFRC_DP', 'CP cloud fraction (deep)',         '1',       fill_halo=.true. )
       call FILE_HISTORY_in( cldfrac_sh    (:,:,:), 'CUMFRC_SH', 'CP cloud fraction (shallow)',      '1',       fill_halo=.true. )
       call FILE_HISTORY_in( kf_nca        (:,:),   'kf_nca',    'advection or cumulus convection timescale for KF', 's', fill_halo=.true. )

       call FILE_HISTORY_in( DENS_t_CP(:,:,:), 'DENS_t_CP', 'tendency DENS in CP', 'kg/m3/s'  , fill_halo=.true. )
       call FILE_HISTORY_in( RHOT_t_CP(:,:,:), 'RHOT_t_CP', 'tendency RHOT in CP', 'K*kg/m3/s', fill_halo=.true. )

       call FILE_HISTORY_in( RHOQV_t_CP(:,:,:), 'QV_t_CP',  'tendency rho*QV in CP', 'kg/m3/s', fill_halo=.true. )
       do iq = 1, N_HYD
          call FILE_HISTORY_in( RHOHYD_t_CP(:,:,:,iq), trim(HYD_NAME(iq))//'_t_CP', &
                                'tendency rho*'//trim(HYD_NAME(iq))//' in CP', 'kg/m3/s', fill_halo=.true. )
       enddo

    endif ! update

    !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
    do j = JSB, JEB
    do i = ISB, IEB
    do k = KS, KE
       DENS_t(k,i,j) = DENS_t(k,i,j) + DENS_t_CP(k,i,j)
       RHOT_t(k,i,j) = RHOT_t(k,i,j) + RHOT_t_CP(k,i,j)
    enddo
    enddo
    enddo

    call ATMOS_PHY_MP_driver_qhyd2qtrc( KA, KS, KE, IA, IS, IE, JA, JS, JE, &
                                        RHOQV_t_CP(:,:,:), RHOHYD_t_CP(:,:,:,:), & ! [IN]
                                        RHOQ_t_CP(:,:,:,QS_MP:QE_MP)             ) ! [OUT]

    do iq = QS_MP, QE_MP
    !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(3)
    do j  = JS, JE
    do i  = IS, IE
    do k  = KS, KE
       RHOQ_t(k,i,j,iq) = RHOQ_t(k,i,j,iq) + RHOQ_t_CP(k,i,j,iq)
    enddo
    enddo
    enddo
    enddo

    if ( STATISTICS_checktotal ) then
       call STATISTICS_total( KA, KS, KE, IA, IS, IE, JA, JS, JE, &
                              DENS_t_CP(:,:,:), 'DENS_t_CP',       &
                              ATMOS_GRID_CARTESC_REAL_VOL(:,:,:),  &
                              ATMOS_GRID_CARTESC_REAL_TOTVOL       )
       call STATISTICS_total( KA, KS, KE, IA, IS, IE, JA, JS, JE, &
                              RHOT_t_CP(:,:,:), 'RHOT_t_CP',       &
                              ATMOS_GRID_CARTESC_REAL_VOL(:,:,:),  &
                              ATMOS_GRID_CARTESC_REAL_TOTVOL       )

       do iq = QS_MP, QE_MP
          call STATISTICS_total( KA, KS, KE, IA, IS, IE, JA, JS, JE, &
                                 RHOQ_t_CP(:,:,:,iq), trim(TRACER_NAME(iq))//'_t_CP', &
                                 ATMOS_GRID_CARTESC_REAL_VOL(:,:,:),                  &
                                 ATMOS_GRID_CARTESC_REAL_TOTVOL                       )
       enddo
    endif

    return
  end subroutine ATMOS_PHY_CP_driver_calc_tendency

end module mod_atmos_phy_cp_driver
