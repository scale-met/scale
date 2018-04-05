!-------------------------------------------------------------------------------
!> module atmosphere / physics / cloud microphysics
!!
!! @par Description
!!          Cloud Microphysics driver
!!
!! @author Team SCALE
!!
!<
!-------------------------------------------------------------------------------
#include "scalelib.h"
module mod_atmos_phy_mp_driver
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_stdio
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
  public :: ATMOS_PHY_MP_driver_tracer_setup
  public :: ATMOS_PHY_MP_driver_setup
  public :: ATMOS_PHY_MP_driver_calc_tendency
  public :: ATMOS_PHY_MP_driver_adjustment
  public :: ATMOS_PHY_MP_driver_qhyd2qtrc

  interface abstract
     subroutine qhyd2qtrc( &
          KA, KS, KE, IA, IS, IE, JA, JS, JE, &
          QV, QHYD, &
          QTRC,     &
          QNUM      )
       use scale_precision
       use scale_atmos_hydrometeor, only: N_HYD
       use mod_atmos_phy_mp_vars, only: QA_MP
       integer, intent(in) :: KA, KS, KE
       integer, intent(in) :: IA, IS, IE
       integer, intent(in) :: JA, JS, JE
       real(RP), intent(in) :: QV(KA,IA,JA)
       real(RP), intent(in) :: QHYD(KA,IA,JA,N_HYD)
       real(RP), intent(out) :: QTRC(KA,IA,JA,QA_MP)
       real(RP), intent(in), optional :: QNUM(KA,IA,JA,N_HYD)
     end subroutine qhyd2qtrc
  end interface abstract
  procedure(qhyd2qtrc), pointer :: ATMOS_PHY_MP_USER_qhyd2qtrc => NULL()
  public :: ATMOS_PHY_MP_USER_qhyd2qtrc

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
  logical,  private :: MP_do_precipitation   = .true.  !> apply sedimentation (precipitation)?
  logical,  private :: MP_do_negative_fixer  = .true.  !> apply negative fixer?
  real(RP), private :: MP_limit_negative    = 1.0_RP   !> Abort if abs(fixed negative vaue) > abs(MP_limit_negative)
  integer,  private :: MP_ntmax_sedimentation = 1      !> number of time step for sedimentation
  real(RP), private :: MP_max_term_vel = 10.0_RP       !> terminal velocity for calculate dt of sedimentation
  real(RP), private :: MP_cldfrac_thleshold            !> thleshold for cloud fraction
  integer,  private :: MP_NSTEP_SEDIMENTATION
  real(RP), private :: MP_RNSTEP_SEDIMENTATION
  real(DP), private :: MP_DTSEC_SEDIMENTATION

  integer, private, allocatable :: hist_vterm_id(:)
  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Config
  subroutine ATMOS_PHY_MP_driver_tracer_setup
    use scale_prc, only: &
       PRC_abort
    use scale_tracer, only: &
       TRACER_regist
    use scale_atmos_hydrometeor, only: &
       ATMOS_HYDROMETEOR_regist
    use scale_atmos_phy_mp_kessler, only: &
       ATMOS_PHY_MP_KESSLER_ntracers,            &
       ATMOS_PHY_MP_KESSLER_nwaters,             &
       ATMOS_PHY_MP_KESSLER_nices,               &
       ATMOS_PHY_MP_KESSLER_tracer_names,        &
       ATMOS_PHY_MP_KESSLER_tracer_descriptions, &
       ATMOS_PHY_MP_KESSLER_tracer_units
    use scale_atmos_phy_mp_tomita08, only: &
       ATMOS_PHY_MP_TOMITA08_ntracers,            &
       ATMOS_PHY_MP_TOMITA08_nwaters,             &
       ATMOS_PHY_MP_TOMITA08_nices,               &
       ATMOS_PHY_MP_TOMITA08_tracer_names,        &
       ATMOS_PHY_MP_TOMITA08_tracer_descriptions, &
       ATMOS_PHY_MP_TOMITA08_tracer_units
    use scale_atmos_phy_mp_sn14, only: &
       ATMOS_PHY_MP_SN14_ntracers,            &
       ATMOS_PHY_MP_SN14_nwaters,             &
       ATMOS_PHY_MP_SN14_nices,               &
       ATMOS_PHY_MP_SN14_tracer_names,        &
       ATMOS_PHY_MP_SN14_tracer_descriptions, &
       ATMOS_PHY_MP_SN14_tracer_units    
    use scale_atmos_phy_mp_suzuki10, only: &
       ATMOS_PHY_MP_suzuki10_tracer_setup,        &
       ATMOS_PHY_MP_suzuki10_ntracers,            &
       ATMOS_PHY_MP_suzuki10_nwaters,             &
       ATMOS_PHY_MP_suzuki10_nices,               &
       ATMOS_PHY_MP_suzuki10_nccn,                &
       ATMOS_PHY_MP_suzuki10_tracer_names,        &
       ATMOS_PHY_MP_suzuki10_tracer_descriptions, &
       ATMOS_PHY_MP_suzuki10_tracer_units
    use mod_atmos_admin, only: &
       ATMOS_PHY_MP_TYPE, &
       ATMOS_sw_phy_mp
    use mod_atmos_phy_mp_vars, only: &
       QA_MP, &
       QS_MP, &
       QE_MP
    implicit none

    integer :: QS2
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_PROGRESS(*) 'Module[Tracer Setup] / Categ[ATMOS PHY_MP] / Origin[SCALE-RM]'

    if ( ATMOS_sw_phy_mp ) then
       select case ( ATMOS_PHY_MP_TYPE )
       case ( 'KESSLER' )
          call ATMOS_HYDROMETEOR_regist( &
               ATMOS_PHY_MP_KESSLER_nwaters,                & ! [IN]
               ATMOS_PHY_MP_KESSLER_nices,                  & ! [IN]
               ATMOS_PHY_MP_KESSLER_tracer_names(:),        & ! [IN]
               ATMOS_PHY_MP_KESSLER_tracer_descriptions(:), & ! [IN]
               ATMOS_PHY_MP_KESSLER_tracer_units(:),        & ! [IN]
               QS_MP                                        ) ! [OUT]
          QA_MP = ATMOS_PHY_MP_KESSLER_ntracers
       case ( 'TOMITA08' )
          call ATMOS_HYDROMETEOR_regist( &
               ATMOS_PHY_MP_TOMITA08_nwaters,                & ! [IN]
               ATMOS_PHY_MP_TOMITA08_nices,                  & ! [IN]
               ATMOS_PHY_MP_TOMITA08_tracer_names(:),        & ! [IN]
               ATMOS_PHY_MP_TOMITA08_tracer_descriptions(:), & ! [IN]
               ATMOS_PHY_MP_TOMITA08_tracer_units(:),        & ! [IN]
               QS_MP                                         ) ! [OUT]
          QA_MP = ATMOS_PHY_MP_TOMITA08_ntracers
       case( 'SN14' )
          call ATMOS_HYDROMETEOR_regist( &
               ATMOS_PHY_MP_SN14_nwaters,                   & ! [IN]
               ATMOS_PHY_MP_SN14_nices,                     & ! [IN]
               ATMOS_PHY_MP_SN14_tracer_names(1:6),         & ! [IN]
               ATMOS_PHY_MP_SN14_tracer_descriptions(1:6),  & ! [IN]
               ATMOS_PHY_MP_SN14_tracer_units(1:6),         & ! [IN]
               QS_MP                                        ) ! [OUT]
          call TRACER_regist( QS2,                          & ! [OUT]
               5,                                           & ! [IN]
               ATMOS_PHY_MP_SN14_tracer_names(7:11),        & ! [IN]
               ATMOS_PHY_MP_SN14_tracer_descriptions(7:11), & ! [IN]
               ATMOS_PHY_MP_SN14_tracer_units(7:11)         ) ! [IN]
          QA_MP = ATMOS_PHY_MP_SN14_ntracers
       case ( 'SUZUKI10' )
          call ATMOS_PHY_MP_suzuki10_tracer_setup
          call ATMOS_HYDROMETEOR_regist( &
               ATMOS_PHY_MP_suzuki10_nwaters,                & ! [IN]
               ATMOS_PHY_MP_suzuki10_nices,                  & ! [IN]
               ATMOS_PHY_MP_suzuki10_tracer_names(:),        & ! [IN]
               ATMOS_PHY_MP_suzuki10_tracer_descriptions(:), & ! [IN]
               ATMOS_PHY_MP_suzuki10_tracer_units(:),        & ! [IN]
               QS_MP                                         ) ! [OUT]
          if( ATMOS_PHY_MP_suzuki10_nccn > 0 ) then
             call TRACER_regist( QS2,                                                                     & ! [OUT]
                                 ATMOS_PHY_MP_suzuki10_nccn,                                              & ! [IN]
                                 ATMOS_PHY_MP_suzuki10_tracer_names       ( ATMOS_PHY_MP_suzuki10_nwaters &
                                                                          + ATMOS_PHY_MP_suzuki10_nices   &
                                                                          + 2 : ),                        & ! [IN]
                                 ATMOS_PHY_MP_suzuki10_tracer_descriptions( ATMOS_PHY_MP_suzuki10_nwaters &
                                                                          + ATMOS_PHY_MP_suzuki10_nices   &
                                                                          + 2 : ),                        & ! [IN]
                                 ATMOS_PHY_MP_suzuki10_tracer_units       ( ATMOS_PHY_MP_suzuki10_nwaters &
                                                                          + ATMOS_PHY_MP_suzuki10_nices   &
                                                                          + 2 : )                         ) ! [IN]
          end if
          QA_MP = ATMOS_PHY_MP_suzuki10_ntracers
       case default
          LOG_ERROR("ATMOS_PHY_MP_driver_tracer_setup",*) 'ATMOS_PHY_MP_TYPE is invalud: ', trim(ATMOS_PHY_MP_TYPE)
          call PRC_abort
       end select

       QE_MP = QS_MP + QA_MP - 1

    else
       QA_MP = 0
       QS_MP = -1
       QE_MP = -1
    end if

    return
  end subroutine ATMOS_PHY_MP_driver_tracer_setup

  !-----------------------------------------------------------------------------
  !> Setup
  subroutine ATMOS_PHY_MP_driver_setup
    use scale_prc, only: &
       PRC_abort
    use scale_atmos_grid_cartesC, only: &
       CDZ => ATMOS_GRID_CARTESC_CDZ
    use scale_const, only: &
       EPS => CONST_EPS
    use scale_time, only: &
       TIME_DTSEC_ATMOS_PHY_MP
    use scale_atmos_phy_mp_kessler, only: &
       ATMOS_PHY_MP_KESSLER_setup
    use scale_atmos_phy_mp_tomita08, only: &
       ATMOS_PHY_MP_TOMITA08_setup
    use scale_atmos_phy_mp_sn14, only: &
       ATMOS_PHY_MP_SN14_setup
    use scale_atmos_phy_mp_suzuki10, only: &
       ATMOS_PHY_MP_suzuki10_setup
    use scale_file_history, only: &
       FILE_HISTORY_reg
    use mod_atmos_admin, only: &
       ATMOS_PHY_MP_TYPE, &
       ATMOS_sw_phy_mp
    use mod_atmos_phy_mp_vars, only: &
       ATMOS_PHY_MP_cldfrac_thleshold, &
       SFLX_rain => ATMOS_PHY_MP_SFLX_rain, &
       SFLX_snow => ATMOS_PHY_MP_SFLX_snow
    use mod_atmos_phy_mp_vars, only: &
       QS_MP, &
       QE_MP
    implicit none
    
    NAMELIST / PARAM_ATMOS_PHY_MP / &
       MP_do_precipitation,     &
       MP_do_negative_fixer,    &
       MP_limit_negative,       &
       MP_ntmax_sedimentation,  &
       MP_max_term_vel,         &
       MP_cldfrac_thleshold

    integer :: nstep_max

    integer :: iq
    integer :: ierr
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_PROGRESS(*) 'Module[DRIVER] / Categ[ATMOS PHY_MP] / Origin[SCALE-RM]'

    if ( ATMOS_sw_phy_mp ) then

       MP_cldfrac_thleshold = EPS

       !--- read namelist
       rewind(IO_FID_CONF)
       read(IO_FID_CONF,nml=PARAM_ATMOS_PHY_MP,iostat=ierr)
       if( ierr < 0 ) then !--- missing
          LOG_INFO("ATMOS_PHY_MP_driver_setup",*) 'Not found namelist. Default used.'
       elseif( ierr > 0 ) then !--- fatal error
          LOG_ERROR("ATMOS_PHY_MP_driver_setup",*) 'Not appropriate names in namelist PARAM_ATMOS_PHY_MP. Check!'
          call PRC_abort
       endif
       if( IO_NML ) write(IO_FID_NML,nml=PARAM_ATMOS_PHY_MP)

       nstep_max = ceiling( ( TIME_DTSEC_ATMOS_PHY_MP * MP_max_term_vel ) / minval( CDZ(:) ) )
       MP_ntmax_sedimentation = max( MP_ntmax_sedimentation, nstep_max )

       MP_NSTEP_SEDIMENTATION  = MP_ntmax_sedimentation
       MP_RNSTEP_SEDIMENTATION = 1.0_RP / real(MP_ntmax_sedimentation,kind=RP)
       MP_DTSEC_SEDIMENTATION  = TIME_DTSEC_ATMOS_PHY_MP * MP_RNSTEP_SEDIMENTATION

       ATMOS_PHY_MP_cldfrac_thleshold = MP_cldfrac_thleshold

       LOG_NEWLINE
       LOG_INFO("ATMOS_PHY_MP_driver_setup",*) 'Enable negative fixer?                    : ', MP_do_negative_fixer
       LOG_INFO("ATMOS_PHY_MP_driver_setup",*) 'Value limit of negative fixer (abs)       : ', abs(MP_limit_negative)
       LOG_INFO("ATMOS_PHY_MP_driver_setup",*) 'Enable sedimentation (precipitation)?     : ', MP_do_precipitation
       LOG_INFO("ATMOS_PHY_MP_driver_setup",*) 'Timestep of sedimentation is divided into : ', MP_ntmax_sedimentation, 'step'
       LOG_INFO("ATMOS_PHY_MP_driver_setup",*) 'DT of sedimentation                       : ', MP_DTSEC_SEDIMENTATION, '[s]'

       select case ( ATMOS_PHY_MP_TYPE )
       case ( 'KESSLER' )
          call ATMOS_PHY_MP_kessler_setup
       case ( 'TOMITA08' )
          call ATMOS_PHY_MP_tomita08_setup( &
               KA, KS, KE, IA, ISB, IEB, JA, JSB, JEB )
       case ( 'SN14' )
          call ATMOS_PHY_MP_sn14_setup( &
               KA, IA, JA )
       case ( 'SUZUKI10' )
          call ATMOS_PHY_MP_suzuki10_setup( &
               KA, IA, JA )
       end select

    else

       LOG_INFO("ATMOS_PHY_MP_driver_setup",*) 'this component is never called.'
       LOG_INFO("ATMOS_PHY_MP_driver_setup",*) 'SFLX_rain and SFLX_snow is set to zero.'
       SFLX_rain(:,:) = 0.0_RP
       SFLX_snow(:,:) = 0.0_RP

    endif


    ! history putput
    if ( MP_do_precipitation ) then
       allocate( hist_vterm_id(QS_MP+1:QE_MP) )
       do iq = QS_MP+1, QE_MP
          call FILE_HISTORY_reg( 'Vterm_'//trim(TRACER_NAME(iq)), 'terminal velocity of '//trim(TRACER_NAME(iq)), 'm/s', hist_vterm_id(iq) )
       end do
    end if


    return
  end subroutine ATMOS_PHY_MP_driver_setup

  !-----------------------------------------------------------------------------
  !> adjustment
  subroutine ATMOS_PHY_MP_driver_adjustment
    use scale_atmos_hydrometeor, only: &
       ATMOS_HYDROMETEOR_dry, &
       I_QV, &
       QHA, &
       QHS, &
       QHE
    use scale_atmos_phy_mp_common, only: &
       ATMOS_PHY_MP_negative_fixer
    use mod_atmos_vars, only: &
       DENS, &
       RHOT, &
       QTRC
    use mod_atmos_phy_mp_vars, only: &
       QA_MP

    real(RP) :: DENS0(KA,IA,JA)

    integer :: k, i, j, iq

    if ( MP_do_negative_fixer .and. (.not. ATMOS_HYDROMETEOR_dry) ) then

!OCL XFILL
       DENS0(:,:,:) = DENS(:,:,:)

       call ATMOS_PHY_MP_negative_fixer( &
            KA, KS, KE, IA, 1, IA, JA, 1, JA, QHA, &
            MP_limit_negative,                    & ! [IN]
            DENS(:,:,:),                          & ! [INOUT]
            QTRC(:,:,:,I_QV), QTRC(:,:,:,QHS:QHE) ) ! [INOUT]

       ! for non-mass tracers, such as number density
       do iq = QHE+1, QA_MP
       do j = 1, JA
       do i = 1, IA
       do k = KS, KE
          QTRC(k,i,j,iq) = max( QTRC(k,i,j,iq), 0.0_RP )
       end do
       end do
       end do
       end do

       do j = 1, JA
       do i = 1, IA
       do k = KS, KE
          RHOT(k,i,j) = RHOT(k,i,j) * DENS(k,i,j) / DENS0(k,i,j)
       end do
       end do
       end do

    end if

    return
  end subroutine ATMOS_PHY_MP_driver_adjustment

  !-----------------------------------------------------------------------------
  !> calculate tendency
  subroutine ATMOS_PHY_MP_driver_calc_tendency( update_flag )
    use scale_const, only: &
       PRE00 => CONST_PRE00
    use scale_time, only: &
       dt_MP => TIME_DTSEC_ATMOS_PHY_MP
    use scale_atmos_grid_cartesC, only: &
       Z  => ATMOS_GRID_CARTESC_CZ, &
       DZ => ATMOS_GRID_CARTESC_CDZ
    use scale_atmos_grid_cartesC_real, only: &
       REAL_CZ => ATMOS_GRID_CARTESC_REAL_CZ, &
       REAL_FZ => ATMOS_GRID_CARTESC_REAL_FZ, &
       ATMOS_GRID_CARTESC_REAL_VOL,       &
       ATMOS_GRID_CARTESC_REAL_TOTVOL,    &
       ATMOS_GRID_CARTESC_REAL_VOLWXY,    &
       ATMOS_GRID_CARTESC_REAL_TOTVOLWXY, &
       ATMOS_GRID_CARTESC_REAL_VOLZUY,    &
       ATMOS_GRID_CARTESC_REAL_TOTVOLZUY, &
       ATMOS_GRID_CARTESC_REAL_VOLZXV,    &
       ATMOS_GRID_CARTESC_REAL_TOTVOLZXV
    use scale_statistics, only: &
       STATISTICS_checktotal, &
       STATISTICS_total
    use scale_file_history, only: &
       FILE_HISTORY_in
    use scale_atmos_hydrometeor, only: &
       I_QV,  &
       N_HYD, &
       QHA,   &
       QHS,   &
       QHE,   &
       QLA,   &
       QIA
    use scale_atmos_phy_mp_common, only: &
       ATMOS_PHY_MP_precipitation, &
       ATMOS_PHY_MP_precipitation_momentum
    use scale_atmos_refstate, only: &
       REFSTATE_dens => ATMOS_REFSTATE_dens
    use scale_atmos_phy_mp_kessler, only: &
       ATMOS_PHY_MP_KESSLER_adjustment, &
       ATMOS_PHY_MP_KESSLER_terminal_velocity
    use scale_atmos_phy_mp_tomita08, only: &
       ATMOS_PHY_MP_TOMITA08_adjustment, &
       ATMOS_PHY_MP_TOMITA08_terminal_velocity
    use scale_atmos_phy_mp_sn14, only: &
       ATMOS_PHY_MP_sn14_tendency, &
       ATMOS_PHY_MP_sn14_terminal_velocity
    use scale_atmos_phy_mp_suzuki10, only: &
       ATMOS_PHY_MP_suzuki10_tendency, &
       ATMOS_PHY_MP_suzuki10_terminal_velocity
    use scale_file_history, only: &
       FILE_HISTORY_query, &
       FILE_HISTORY_put
    use mod_atmos_admin, only: &
       ATMOS_PHY_MP_TYPE
    use mod_atmos_vars, only: &
       DENS   => DENS_av, &
       MOMZ   => MOMZ_av, &
       U, &
       V, &
       W, &
       POTT, &
       QTRC   => QTRC_av, &
       DENS_t => DENS_tp, &
       MOMZ_t => MOMZ_tp, &
       RHOU_t => RHOU_tp, &
       RHOV_t => RHOV_tp, &
!       RHOT_t => RHOT_tp, &
       RHOQ_t => RHOQ_tp, &
       RHOH   => RHOH_p, &
       TEMP, &
       PRES, &
       Qdry, &
       CVtot, &
       CPtot, &
       EXNER, &
       RHOT   => RHOT_av
    use mod_atmos_phy_mp_vars, only: &
       QA_MP, &
       QS_MP, &
       QE_MP, &
       DENS_t_MP => ATMOS_PHY_MP_DENS_t,    &
       MOMZ_t_MP => ATMOS_PHY_MP_MOMZ_t,    &
!       RHOT_t_MP => ATMOS_PHY_MP_RHOT_t,    &
       RHOU_t_MP => ATMOS_PHY_MP_RHOU_t,    &
       RHOV_t_MP => ATMOS_PHY_MP_RHOV_t,    &
       RHOQ_t_MP => ATMOS_PHY_MP_RHOQ_t,    &
       RHOH_MP   => ATMOS_PHY_MP_RHOH,      &
       EVAPORATE => ATMOS_PHY_MP_EVAPORATE, &
       SFLX_rain => ATMOS_PHY_MP_SFLX_rain, &
       SFLX_snow => ATMOS_PHY_MP_SFLX_snow
    use mod_atmos_phy_ae_vars, only: &
       CCN_t => ATMOS_PHY_AE_CCN_t
    implicit none

    logical, intent(in) :: update_flag

    real(RP) :: RHOE_t(KA,IA,JA)
    real(RP) :: TEMP1 (KA,IA,JA)
    real(RP) :: PRES1 (KA,IA,JA)
    real(RP) :: CPtot1(KA,IA,JA)
    real(RP) :: CVtot1(KA,IA,JA)
    real(RP) :: CCN   (KA,IA,JA)
    real(RP) :: QHYD  (KA,IA,JA,6)
    real(RP) :: vterm (KA,QS_MP+1:QE_MP)
    real(RP), target :: QTRC1(KA,IA,JA,QS_MP:QE_MP)

    real(RP) :: FLX_hydro(KA)
    real(RP) :: DENS2    (KA)
    real(RP) :: TEMP2    (KA)
    real(RP) :: PRES2    (KA)
    real(RP) :: CPtot2   (KA)
    real(RP) :: CVtot2   (KA)
    real(RP) :: RHOE     (KA)
    real(RP) :: RHOE2    (KA)
    real(RP) :: RHOQ2    (KA,QS_MP+1:QE_MP)
    real(RP) :: mflux    (KA)
    real(RP) :: sflux    (2)  !> 1: rain, 2: snow

    real(RP) :: FDZ (KA)
    real(RP) :: RFDZ(KA)
    real(RP) :: RCDZ(KA)

    real(RP) :: CPtot_t(KA,IA,JA), CVtot_t(KA,IA,JA)
    real(RP) :: CP_t, CV_t

    real(RP) :: precip   (IA,JA)

    ! for history output
    real(RP), allocatable :: vterm_hist(:,:,:,:)
    integer :: hist_vterm_idx(QS_MP+1:QE_MP)
    logical :: flag
    integer :: ih

    integer :: k, i, j, iq
    integer :: step
    !---------------------------------------------------------------------------

    if ( update_flag ) then

       CCN(:,:,:) = CCN_t(:,:,:) * dt_MP

       select case ( ATMOS_PHY_MP_TYPE )
       case ( 'KESSLER' )
!OCL XFILL
          TEMP1(:,:,:) = TEMP(:,:,:)
!OCL XFILL
          QTRC1(:,:,:,QS_MP:QE_MP) = QTRC(:,:,:,QS_MP:QE_MP)
!OCL XFILL
          CVtot1(:,:,:) = CVtot(:,:,:)
!OCL XFILL
          CPtot1(:,:,:) = CPtot(:,:,:)

          call ATMOS_PHY_MP_kessler_adjustment( &
               KA, KS, KE, IA, ISB, IEB, JA, JSB, JEB, &
               DENS(:,:,:), PRES(:,:,:), dt_MP,                                      & ! [IN]
               TEMP1(:,:,:), QTRC1(:,:,:,QS_MP:QE_MP), CPtot1(:,:,:), CVtot1(:,:,:), & ! [INOUT]
               RHOE_t(:,:,:), EVAPORATE(:,:,:)                                       ) ! [OUT]

          do iq = QS_MP, QE_MP
          do j = JSB, JEB
          do i = ISB, IEB
          do k = KS, KE
             RHOQ_t_MP(k,i,j,iq) = ( QTRC1(k,i,j,iq) - QTRC(k,i,j,iq) ) * DENS(k,i,j) / dt_MP
          enddo
          enddo
          enddo
          enddo

          do j = JSB, JEB
          do i = ISB, IEB
          do k = KS, KE
             CPtot_t(k,i,j) = ( CPtot1(k,i,j) - CPtot(k,i,j) ) / dt_MP
             CVtot_t(k,i,j) = ( CVtot1(k,i,j) - CVtot(k,i,j) ) / dt_MP
          end do
          end do
          end do

       case ( 'TOMITA08' )
!OCL XFILL
          do j = JSB, JEB
          do i = ISB, IEB
          do k = KS, KE
             TEMP1(k,i,j) = TEMP(k,i,j)
          end do
          end do
          end do
!OCL XFILL
          do iq = QS_MP, QE_MP
          do j = JSB, JEB
          do i = ISB, IEB
          do k = KS, KE
             QTRC1(k,i,j,iq) = QTRC(k,i,j,iq)
          end do
          end do
          end do
          end do
!OCL XFILL
          do j = JSB, JEB
          do i = ISB, IEB
          do k = KS, KE
             CVtot1(k,i,j) = CVtot(k,i,j)
          end do
          end do
          end do
!OCL XFILL
          do j = JSB, JEB
          do i = ISB, IEB
          do k = KS, KE
             CPtot1(k,i,j) = CPtot(k,i,j)
          end do
          end do
          end do

          call ATMOS_PHY_MP_tomita08_adjustment( &
               KA, KS, KE, IA, ISB, IEB, JA, JSB, JEB, &
               DENS(:,:,:), PRES(:,:,:), CCN(:,:,:), dt_MP,                          & ! [IN]
               TEMP1(:,:,:), QTRC1(:,:,:,QS_MP:QE_MP), CPtot1(:,:,:), CVtot1(:,:,:), & ! [INOUT]
               RHOE_t(:,:,:), EVAPORATE(:,:,:)                                       ) ! [OUT]

          do iq = QS_MP, QE_MP
          do j = JSB, JEB
          do i = ISB, IEB
          do k = KS, KE
             RHOQ_t_MP(k,i,j,iq) = ( QTRC1(k,i,j,iq) - QTRC(k,i,j,iq) ) * DENS(k,i,j) / dt_MP
          enddo
          enddo
          enddo
          enddo

          do j = JSB, JEB
          do i = ISB, IEB
          do k = KS, KE
             CPtot_t(k,i,j) = ( CPtot1(k,i,j) - CPtot(k,i,j) ) / dt_MP
             CVtot_t(k,i,j) = ( CVtot1(k,i,j) - CVtot(k,i,j) ) / dt_MP
          end do
          end do
          end do

       case ( 'SN14' )

          call ATMOS_PHY_MP_sn14_tendency( &
               KA, KS, KE, IA, IS, IE, JA, JS, JE, &
               DENS(:,:,:), MOMZ(:,:,:), QTRC(:,:,:,QS_MP:QE_MP), PRES(:,:,:), TEMP(:,:,:),               & ! [IN]
               Qdry(:,:,:), CPtot(:,:,:), CVtot(:,:,:), CCN(:,:,:), dt_MP, Z(:), DZ(:),                   & ! [IN]
               RHOQ_t_MP(:,:,:,QS_MP:QE_MP), RHOE_t(:,:,:), CPtot_t(:,:,:), CVtot_t(:,:,:), EVAPORATE(:,:,:) ) ! [OUT]

       case ( 'SUZUKI10' )

          call ATMOS_PHY_MP_suzuki10_tendency( KA, KS,  KE, IA, ISB, IEB, JA, JSB, JEB, KIJMAX, &
                                               dt_MP,                                  & ! [IN]
                                               DENS(:,:,:),  PRES(:,:,:), TEMP(:,:,:), & ! [IN]
                                               QTRC(:,:,:,QS_MP:QE_MP), QDRY(:,:,:),   & ! [IN]
                                               CPtot(:,:,:), CVtot(:,:,:),             & ! [IN]
                                               CCN(:,:,:),                             & ! [IN] 
                                               RHOQ_t_MP(:,:,:,QS_MP:QE_MP),           & ! [OUT]
                                               RHOE_t(:,:,:),                          & ! [OUT]
                                               CPtot_t(:,:,:), CVtot_t(:,:,:),         & ! [OUT]
                                               EVAPORATE(:,:,:)                        ) ! [OUT]

       end select


       do j = JSB, JEB
       do i = ISB, IEB
       do k = KS, KE
          RHOH_MP(k,i,j) = RHOE_t(k,i,j) &
                  - ( CPtot_t(k,i,j) + log( PRES(k,i,j) / PRE00 ) * ( CVtot(k,i,j) / CPtot(k,i,j) * CPtot_t(k,i,j) - CVtot_t(k,i,j) ) ) &
                  * DENS(k,i,j) * TEMP(k,i,j)
!          RHOT_t_MP(k,i,j) = RHOE_t(k,i,j) / ( EXNER(k,i,j) * CPtot(k,i,j) ) &
!                  - RHOT(k,i,j) * CPtot_t(k,i,j) / CPtot(k,i,j) &
!                  - RHOT(k,i,j) * CVtot(k,i,j) / ( CPtot(k,i,j) - CVtot(k,i,j) ) &
!                  * log( EXNER(k,i,j) ) * ( CPtot_t(k,i,j) / CPtot(k,i,j) - CVtot_t(k,i,j) / CVtot(k,i,j) )
       end do
       end do
       end do

       if ( MP_do_precipitation ) then

          call PROF_rapstart('MP_Precipitation', 2)

          ! prepare for history output
          hist_vterm_idx(:) = -1
          ih = 0
          do iq = QS_MP+1, QE_MP
             call FILE_HISTORY_query( hist_vterm_id(iq), flag )
             if ( flag ) then
                ih = ih + 1
                hist_vterm_idx(iq) = ih
             end if
          end do
          if ( ih > 0 ) then
             allocate( vterm_hist(KA,IA,JA,ih) )
             vterm_hist(:,:,:,:) = 0.0_RP
          end if

          !$omp parallel do default(none) OMP_SCHEDULE_ collapse(2) &
          !$omp shared (KA,KS,KE,ISB,IEB,JSB,JEB,QS_MP,QE_MP,QHA,QLA,QIA, &
          !$omp         PRE00, &
          !$omp         ATMOS_PHY_MP_TYPE, &
          !$omp         dt_MP,MP_NSTEP_SEDIMENTATION,MP_DTSEC_SEDIMENTATION,MP_RNSTEP_SEDIMENTATION, &
          !$omp         REAL_CZ,REAL_FZ, &
          !$omp         DENS,MOMZ,U,V,RHOT,TEMP,PRES,QTRC,CPtot,CVtot,EXNER, &
          !$omp         DENS_t_MP,MOMZ_t_MP,RHOU_t_MP,RHOV_t_MP,RHOQ_t_MP,RHOH_MP, &
          !$omp         SFLX_rain,SFLX_snow, &
          !$omp         REFSTATE_dens, &
          !$omp         vterm_hist,hist_vterm_idx) &
          !$omp private(i,j,k,iq,step, &
          !$omp         FDZ,RFDZ,RCDZ, &
          !$omp         DENS2,TEMP2,PRES2,CPtot2,CVtot2,RHOE,RHOE2,RHOQ2, &
          !$omp         vterm,mflux,sflux,FLX_hydro,CP_t,CV_t)
          do j = JSB, JEB
          do i = ISB, IEB

             FDZ(KS-1) = REAL_CZ(KS,i,j) - REAL_FZ(KS-1,i,j)
             RFDZ(KS-1) = 1.0_RP / FDZ(KS-1)
             do k = KS, KE
                FDZ(k) = REAL_CZ(k+1,i,j) - REAL_CZ(k  ,i,j)
                RFDZ(k) = 1.0_RP / FDZ(k)
                RCDZ(k) = 1.0_RP / ( REAL_FZ(k  ,i,j) - REAL_FZ(k-1,i,j) )
             enddo

             do k = KS, KE
                DENS2(k)  = DENS(k,i,j)
                TEMP2(k)  = TEMP(k,i,j)
                PRES2(k)  = PRES(k,i,j)
                CPtot2(k) = CPtot(k,i,j)
                CVtot2(k) = CVtot(k,i,j)
                RHOE(k)   = TEMP(k,i,j) * CVtot(k,i,j) * DENS2(k)
                RHOE2(k)  = RHOE(k)
             end do
             do iq = QS_MP+1, QE_MP
             do k = KS, KE
                RHOQ2(k,iq) = DENS2(k) * QTRC(k,i,j,iq)
             end do
             end do

             SFLX_rain(i,j) = 0.0_RP
             SFLX_snow(i,j) = 0.0_RP
             FLX_hydro(:) = 0.0_RP
             do step = 1, MP_NSTEP_SEDIMENTATION

                select case ( ATMOS_PHY_MP_TYPE )
                case ( 'KESSLER' )
                   call ATMOS_PHY_MP_kessler_terminal_velocity( &
                        KA, KS, KE, &
                        DENS2(:), RHOQ2(:,:), & ! [IN]
                        REFSTATE_dens(:,i,j), & ! [IN]
                        vterm(:,:)            ) ! [OUT]
                case ( 'TOMITA08' )
                   call ATMOS_PHY_MP_tomita08_terminal_velocity( &
                        KA, KS, KE, &
                        DENS2(:), TEMP2(:), RHOQ2(:,:), & ! [IN]
                        vterm(:,:)                      ) ! [OUT]
                case ( 'SN14' )
                   call ATMOS_PHY_MP_sn14_terminal_velocity( &
                        KA, KS, KE, &
                        DENS2(:), TEMP2(:), RHOQ2(:,:), PRES2(:), & ! [IN]
                        vterm(:,:)                                ) ! [OUT]
                case ( 'SUZUKI10' )
                   call ATMOS_PHY_MP_suzuki10_terminal_velocity( &
                        KA,        & ! [IN]
                        vterm(:,:) ) ! [OUT]
                case default
                   vterm(:,:) = 0.0_RP ! tentative
                end select

                ! store to history output
                do iq = QS_MP+1, QE_MP
                   if ( hist_vterm_idx(iq) > 0 ) then
                      do k = KS, KE
                         vterm_hist(k,i,j,hist_vterm_idx(iq)) = vterm_hist(k,i,j,hist_vterm_idx(iq)) &
                                                              + vterm(k,iq) * MP_RNSTEP_SEDIMENTATION
                      end do
                   end if
                end do

                call ATMOS_PHY_MP_precipitation( &
                     KA, KS, KE, QE_MP-QS_MP, QLA, QIA, &
                     TEMP2(:), vterm(:,:),   & ! [IN]
                     FDZ(:), RCDZ(:),        & ! [IN]
                     MP_DTSEC_SEDIMENTATION, & ! [IN]
                     i, j,                   & ! [IN]
                     DENS2(:), RHOQ2(:,:),   & ! [INOUT]
                     CPtot2(:), CVtot2(:),   & ! [INOUT]
                     RHOE2(:),               & ! [INOUT]
                     mflux(:), sflux(:)      ) ! [OUT]

                do k = KS, KE
                   TEMP2(k) = RHOE2(k) / ( DENS2(k) * CVtot2(k) )
                end do

                do k = KS-1, KE-1
                   FLX_hydro(k) = FLX_hydro(k) + mflux(k) * MP_RNSTEP_SEDIMENTATION
                enddo

                SFLX_rain(i,j) = SFLX_rain(i,j) - sflux(1) * MP_RNSTEP_SEDIMENTATION
                SFLX_snow(i,j) = SFLX_snow(i,j) - sflux(2) * MP_RNSTEP_SEDIMENTATION

             enddo

!OCL XFILL
             do k = KS, KE
                DENS_t_MP(k,i,j) = ( DENS2(k) - DENS(k,i,j) ) / dt_MP
             end do

             do k = KS, KE
                CP_t = ( CPtot2(k) - CPtot(k,i,j) ) / dt_MP
                CV_t = ( CVtot2(k) - CVtot(k,i,j) ) / dt_MP
                RHOH_MP(k,i,j) = RHOH_MP(k,i,j) &
                     + ( RHOE2(k) - RHOE(k) ) / dt_MP &
                     - ( CP_t + log( PRES(k,i,j) / PRE00 ) * ( CVtot(k,i,j) / CPtot(k,i,j) * CP_t - CV_t ) ) &
                     * DENS(k,i,j) * TEMP(k,i,j)
!                RHOT_t_MP(k,i,j) = RHOT_t_MP(k,i,j) &
!                     + ( RHOE2(k) - RHOE(k) ) / ( dt_MP * EXNER(k,i,j) * CPtot(k,i,j) ) &
!                     - RHOT(k,i,j) * CP_t / CPtot(k,i,j) &
!                     - RHOT(k,i,j) * CVtot(k,i,j) / ( CPtot(k,i,j) - CVtot(k,i,j) ) &
!                     * log( EXNER(k,i,j) ) * ( CP_t / CPtot(k,i,j) - CV_t / CVtot(k,i,j) )
             end do

             do iq = QS_MP+1, QE_MP
             do k  = KS, KE
                RHOQ_t_MP(k,i,j,iq) = RHOQ_t_MP(k,i,j,iq) &
                     + ( RHOQ2(k,iq) - DENS(k,i,j) * QTRC(k,i,j,iq) ) / dt_MP
             enddo
             enddo

             call ATMOS_PHY_MP_precipitation_momentum( &
                     KA, KS, KE, &
                     DENS(:,i,j), MOMZ(:,i,j), U(:,i,j), V(:,i,j),        & ! [IN]
                     FLX_hydro(:),                                        & ! [IN]
                     RCDZ(:), RFDZ(:),                                    & ! [IN]
                     MOMZ_t_MP(:,i,j), RHOU_t_MP(:,i,j), RHOV_t_MP(:,i,j) ) ! [OUT]

          enddo
          enddo

          ! history output
          do iq = QS_MP+1, QE_MP
             if ( hist_vterm_idx(iq) > 0 ) &
                call FILE_HISTORY_put( hist_vterm_id(iq), vterm_hist(:,:,:,hist_vterm_idx(iq)) )
          end do
          if ( allocated( vterm_hist ) ) deallocate( vterm_hist )

          call PROF_rapend  ('MP_Precipitation', 2)

       end if

!OCL XFILL
       do j = JS, JE
       do i = IS, IE
          precip(i,j) = SFLX_rain(i,j) + SFLX_snow(i,j)
       end do
       end do

       call FILE_HISTORY_in( SFLX_rain(:,:),   'RAIN_MP',   'surface rain rate by MP',          'kg/m2/s',  fill_halo=.true. )
       call FILE_HISTORY_in( SFLX_snow(:,:),   'SNOW_MP',   'surface snow rate by MP',          'kg/m2/s',  fill_halo=.true. )
       call FILE_HISTORY_in( precip   (:,:),   'PREC_MP',   'surface precipitation rate by MP', 'kg/m2/s',  fill_halo=.true. )
       call FILE_HISTORY_in( EVAPORATE(:,:,:), 'EVAPORATE', 'evaporated cloud number',          'num/m3/s', fill_halo=.true. )

       call FILE_HISTORY_in( DENS_t_MP(:,:,:), 'DENS_t_MP', 'tendency DENS in MP',         'kg/m3/s'  , fill_halo=.true. )
       call FILE_HISTORY_in( MOMZ_t_MP(:,:,:), 'MOMZ_t_MP', 'tendency MOMZ in MP',         'kg/m2/s2' , fill_halo=.true. )
       call FILE_HISTORY_in( RHOU_t_MP(:,:,:), 'RHOU_t_MP', 'tendency RHOU in MP',         'kg/m2/s2' , fill_halo=.true. )
       call FILE_HISTORY_in( RHOV_t_MP(:,:,:), 'RHOV_t_MP', 'tendency RHOV in MP',         'kg/m2/s2' , fill_halo=.true. )
!       call FILE_HISTORY_in( RHOT_t_MP(:,:,:), 'RHOT_t_MP', 'tendency RHOT in MP',         'K*kg/m3/s', fill_halo=.true. )
       call FILE_HISTORY_in( RHOH_MP  (:,:,:), 'RHOH_MP',   'diabatic heating rate in MP', 'J/m3/s',    fill_halo=.true. )
       do iq = QS_MP, QE_MP
          call FILE_HISTORY_in( RHOQ_t_MP(:,:,:,iq), trim(TRACER_NAME(iq))//'_t_MP', &
                        'tendency rho*'//trim(TRACER_NAME(iq))//' in MP', 'kg/m3/s', fill_halo=.true. )
       enddo

    endif

    !$omp parallel do default(none) private(i,j,k) OMP_SCHEDULE_ collapse(2) &
    !$omp shared(KS,KE,ISB,IEB,JSB,JEB, &
    !$omp        DENS_t_MP,MOMZ_t_MP,RHOU_t_MP,RHOV_t_MP,RHOH_MP, &
    !$omp        DENS_t,MOMZ_t,RHOU_t,RHOV_t,RHOH)
    do j = JSB, JEB
    do i = ISB, IEB
    do k = KS, KE
       DENS_t(k,i,j) = DENS_t(k,i,j) + DENS_t_MP(k,i,j)
       MOMZ_t(k,i,j) = MOMZ_t(k,i,j) + MOMZ_t_MP(k,i,j)
       RHOU_t(k,i,j) = RHOU_t(k,i,j) + RHOU_t_MP(k,i,j)
       RHOV_t(k,i,j) = RHOV_t(k,i,j) + RHOV_t_MP(k,i,j)
!       RHOT_t(k,i,j) = RHOT_t(k,i,j) + RHOT_t_MP(k,i,j)
       RHOH  (k,i,j) = RHOH  (k,i,j) + RHOH_MP  (k,i,j)
    enddo
    enddo
    enddo

    do iq = QS_MP, QE_MP
    !$omp parallel do default(none) private(i,j,k) OMP_SCHEDULE_ &
    !$omp shared(JSB,JEB,ISB,IEB,KS,KE,RHOQ_t,iq,RHOQ_t_MP)
    do j = JSB, JEB
    do i = ISB, IEB
    do k = KS, KE
       RHOQ_t(k,i,j,iq) = RHOQ_t(k,i,j,iq) + RHOQ_t_MP(k,i,j,iq)
    enddo
    enddo
    enddo
    enddo

    if ( STATISTICS_checktotal ) then
       call STATISTICS_total( KA, KS, KE, IA, IS, IE, JA, JS, JE, &
                              DENS_t_MP(:,:,:), 'DENS_t_MP',       &
                              ATMOS_GRID_CARTESC_REAL_VOL(:,:,:),  &
                              ATMOS_GRID_CARTESC_REAL_TOTVOL       )
       call STATISTICS_total( KA, KS, KE, IA, IS, IE, JA, JS, JE, &
                              MOMZ_t_MP(:,:,:), 'MOMZ_t_MP',         &
                              ATMOS_GRID_CARTESC_REAL_VOLWXY(:,:,:), &
                              ATMOS_GRID_CARTESC_REAL_TOTVOLWXY      )
       call STATISTICS_total( KA, KS, KE, IA, IS, IE, JA, JS, JE, &
                              RHOH_MP(:,:,:),   'RHOH_MP',      &
                              ATMOS_GRID_CARTESC_REAL_VOL(:,:,:), &
                              ATMOS_GRID_CARTESC_REAL_TOTVOL      )
!       call STATISTICS_total( KA, KS, KE, IA, IS, IE, JA, JS, JE, &
!                              RHOT_t_MP(:,:,:), 'RHOT_t_MP',      &
!                              ATMOS_GRID_CARTESC_REAL_VOL(:,:,:), &
!                              ATMOS_GRID_CARTESC_REAL_TOTVOL      )

       do iq = QS_MP, QE_MP
          call STATISTICS_total( KA, KS, KE, IA, IS, IE, JA, JS, JE, &
                                 RHOQ_t_MP(:,:,:,iq), trim(TRACER_NAME(iq))//'_t_MP', &
                                 ATMOS_GRID_CARTESC_REAL_VOL(:,:,:),                  &
                                 ATMOS_GRID_CARTESC_REAL_TOTVOL                       )
       enddo
    endif

    return
  end subroutine ATMOS_PHY_MP_driver_calc_tendency

  subroutine ATMOS_PHY_MP_driver_qhyd2qtrc( &
       KA, KS, KE, IA, IS, IE, JA, JS, JE, &
       QV, QHYD, &
       QTRC,     &
       QNUM      )
    use scale_prc, only: &
       PRC_abort
    use scale_atmos_hydrometeor, only: &
       N_HYD
    use mod_atmos_phy_mp_vars, only: &
       QA_MP
    use mod_atmos_admin, only: &
       ATMOS_PHY_MP_TYPE
    use scale_atmos_phy_mp_kessler, only: &
       ATMOS_PHY_MP_KESSLER_qhyd2qtrc
    use scale_atmos_phy_mp_tomita08, only: &
       ATMOS_PHY_MP_TOMITA08_qhyd2qtrc
    use scale_atmos_phy_mp_sn14, only: &
       ATMOS_PHY_MP_SN14_qhyd2qtrc
    use scale_atmos_phy_mp_suzuki10, only: &
       ATMOS_PHY_MP_SUZUKI10_qhyd2qtrc
    integer, intent(in) :: KA, KS, KE
    integer, intent(in) :: IA, IS, IE
    integer, intent(in) :: JA, JS, JE

    real(RP), intent(in) :: QV  (KA,IA,JA)
    real(RP), intent(in) :: QHYD(KA,IA,JA,N_HYD)

    real(RP), intent(out) :: QTRC(KA,IA,JA,QA_MP)

    real(RP), intent(in), optional :: QNUM(KA,IA,JA,N_HYD)

    integer :: k, i, j

    select case( ATMOS_PHY_MP_TYPE )
    case ( "NONE" )
       if ( associated( ATMOS_PHY_MP_USER_qhyd2qtrc ) ) then
          call ATMOS_PHY_MP_USER_qhyd2qtrc( KA, KS, KE, IA, IS, IE, JA, JS, JE, &
                                            QV(:,:,:), QHYD(:,:,:,:), & ! [IN]
                                            QTRC(:,:,:,:),            & ! [OUT]
                                            QNUM=QNUM                 ) ! [IN]
       end if
    case ( "KESSLER" )
       !$omp parallel do OMP_SCHEDULE_
       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          QTRC(k,i,j,1) = QV(k,i,j)
       end do
       end do
       end do
       call ATMOS_PHY_MP_KESSLER_qhyd2qtrc( KA, KS, KE, IA, IS, IE, JA, JS, JE, &
                                            QHYD(:,:,:,:),  & ! [IN]
                                            QTRC(:,:,:,2:)  ) ! [OUT]
    case ( "TOMITA08" )
       !$omp parallel do OMP_SCHEDULE_
       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          QTRC(k,i,j,1) = QV(k,i,j)
       end do
       end do
       end do
       call ATMOS_PHY_MP_TOMITA08_qhyd2qtrc( KA, KS, KE, IA, IS, IE, JA, JS, JE, &
                                             QHYD(:,:,:,:),  & ! [IN]
                                             QTRC(:,:,:,2:)  ) ! [OUT]
    case ( "SN14" )
       !$omp parallel do OMP_SCHEDULE_
       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          QTRC(k,i,j,1) = QV(k,i,j)
       end do
       end do
       end do
       call ATMOS_PHY_MP_SN14_qhyd2qtrc( KA, KS, KE, IA, IS, IE, JA, JS, JE, &
                                         QHYD(:,:,:,:),  & ! [IN]
                                         QTRC(:,:,:,2:), & ! [OUT]
                                         QNUM=QNUM       ) ! [IN]
    case ( "SUZUKI10" )
       !$omp parallel do OMP_SCHEDULE_
       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          QTRC(k,i,j,1) = QV(k,i,j)
       end do
       end do
       end do
       call ATMOS_PHY_MP_SUZUKI10_qhyd2qtrc( KA, KS, KE, IA, IS, IE, JA, JS, JE, &
                                             QHYD(:,:,:,:),  & ! [IN]
                                             QTRC(:,:,:,2:), & ! [OUT]
                                             QNUM=QNUM       ) ! [IN]
    case default
       LOG_ERROR("ATMOS_PHY_MP_driver_qhyd2qtrc",*) 'ATMOS_PHY_MP_TYPE (', trim(ATMOS_PHY_MP_TYPE), ') is not supported'
       call PRC_abort
    end select

    return
  end subroutine ATMOS_PHY_MP_driver_qhyd2qtrc

end module mod_atmos_phy_mp_driver
