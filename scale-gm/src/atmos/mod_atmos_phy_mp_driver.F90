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
  public :: ATMOS_PHY_MP_driver_tracer_setup
  public :: ATMOS_PHY_MP_driver_setup
  public :: ATMOS_PHY_MP_driver_step
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
    LOG_INFO("ATMOS_PHY_MP_driver_tracer_setup",*) 'Setup'

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
    use mod_grd, only: &
       CDZ => GRD_dgz
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
    
    namelist / PARAM_ATMOS_PHY_MP / &
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
    LOG_INFO("ATMOS_PHY_MP_driver_setup",*) 'Setup'

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
       LOG_NML(PARAM_ATMOS_PHY_MP)

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
               KA, KS, KE, IA, IS, IE, JA, JS, JE )
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
       SFLX_rain(:,:,:) = 0.0_RP
       SFLX_snow(:,:,:) = 0.0_RP

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
       QLA, &
       QIA, &
       QHS, &
       QHE
    use scale_atmos_phy_mp_common, only: &
       ATMOS_PHY_MP_negative_fixer
    use mod_atmos_vars, only: &
       TEMP,  &
       CVtot, &
       CPtot, &
       DENS,  &
       RHOE,  &
       QTRC
    use mod_atmos_phy_mp_vars, only: &
       QA_MP

    real(RP) :: DENS0(KA,IA,JA,ADM_lall)

    integer :: k, i, j, iq, l

    if ( MP_do_negative_fixer .and. (.not. ATMOS_HYDROMETEOR_dry) ) then

       do l = 1, ADM_lall

          call ATMOS_PHY_MP_negative_fixer( &
               KA, KS, KE, IA, 1, IA, JA, 1, JA, QLA, QIA, &
               MP_limit_negative,                        & ! [IN]
               DENS(:,:,:,l), TEMP(:,:,:,l),             & ! [INOUT]
               CVtot(:,:,:,l), CPtot(:,:,:,l),           & ! [INOUT]
               QTRC(:,:,:,I_QV,l), QTRC(:,:,:,QHS:QHE,l) ) ! [INOUT]

          do j = 1, JA
          do i = 1, IA
          do k = KS, KE
             RHOE(k,i,j,l) = TEMP(k,i,j,l) * CVtot(k,i,j,l)
          end do
          end do
          end do

          ! for non-mass tracers, such as number density
          do iq = QHE+1, QA_MP
          do j = 1, JA
          do i = 1, IA
          do k = KS, KE
             QTRC(k,i,j,iq,l) = max( QTRC(k,i,j,iq,l), 0.0_RP )
          end do
          end do
          end do
          end do

       end do

    end if

    return
  end subroutine ATMOS_PHY_MP_driver_adjustment

  !-----------------------------------------------------------------------------
  !> time step
  subroutine ATMOS_PHY_MP_driver_step
    use scale_time, only: &
       dt_MP => TIME_DTSEC_ATMOS_PHY_MP
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
       ATMOS_PHY_MP_precipitation_upwind, &
       ATMOS_PHY_MP_precipitation_momentum
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
       ATMOS_vars_calc_diagnostics, &
       DENS, &
       RHOU, &
       RHOV, &
       MOMZ, &
       RHOE, &
       RHOQ, &
       QTRC, &
       U, &
       V, &
       W, &
       POTT, &
       TEMP, &
       PRES, &
       Qdry, &
       Rtot,  &
       CVtot, &
       CPtot, &
       EXNER, &
       CZ, &
       FZ
    use mod_atmos_phy_mp_vars, only: &
       QA_MP, &
       QS_MP, &
       QE_MP, &
       EVAPORATE => ATMOS_PHY_MP_EVAPORATE, &
       SFLX_rain => ATMOS_PHY_MP_SFLX_rain, &
       SFLX_snow => ATMOS_PHY_MP_SFLX_snow
    use mod_atmos_phy_ae_vars, only: &
       CCN => ATMOS_PHY_AE_CCN
    use mod_bsstate, only: &
       rho_bs
    implicit none

    real(RP) :: vterm  (KA,QS_MP+1:QE_MP)
    real(RP) :: MOMZ_t (KA)
    real(RP) :: RHOU_t (KA)
    real(RP) :: RHOV_t (KA)
    real(RP) :: RHOE_t (KA,IA,JA,ADM_lall)
    real(RP) :: RHOQ_t (KA,IA,JA,QS_MP:QE_MP)
    real(RP) :: CPtot_t(KA,IA,JA)
    real(RP) :: CVtot_t(KA,IA,JA)
    real(RP) :: mflux (KA)
    real(RP) :: sflux (2)  !> 1: rain, 2: snow
    real(RP) :: eflux
    real(RP) :: FLX_hydro(KA)
    real(RP) :: REFSTATE_dens(KA)


    real(RP) :: FDZ (KA)
    real(RP) :: RFDZ(KA)
    real(RP) :: RCDZ(KA)

    real(RP) :: precip   (IA,JA,ADM_lall)

    ! for history output
    real(RP), allocatable :: vterm_hist(:,:,:,:,:)
    integer :: hist_vterm_idx(QS_MP+1:QE_MP)
    logical :: flag
    integer :: ih

    integer :: KIJMAX
    integer :: k, i, j, iq, l, ij
    integer :: step
    !---------------------------------------------------------------------------

    select case ( ATMOS_PHY_MP_TYPE )
    case ( 'KESSLER' )

       do l = 1, ADM_lall

          call ATMOS_PHY_MP_kessler_adjustment( &
               KA, KS, KE, IA, IS, IE, JA, JS, JE, &
               DENS(:,:,:,l), PRES(:,:,:,l), dt_MP,                                      & ! [IN]
               TEMP(:,:,:,l), QTRC(:,:,:,QS_MP:QE_MP,l), CPtot(:,:,:,l), CVtot(:,:,:,l), & ! [INOUT]
               RHOE_t(:,:,:,l), EVAPORATE(:,:,:,l)                                         ) ! [OUT]
       end do

    case ( 'TOMITA08' )

       do l = 1, ADM_lall

          call ATMOS_PHY_MP_tomita08_adjustment( &
               KA, KS, KE, IA, IS, IE, JA, JS, JE, &
               DENS(:,:,:,l), PRES(:,:,:,l), CCN(:,:,:,l), dt_MP,                        & ! [IN]
               TEMP(:,:,:,l), QTRC(:,:,:,QS_MP:QE_MP,l), CPtot(:,:,:,l), CVtot(:,:,:,l), & ! [INOUT]
               RHOE_t(:,:,:,l), EVAPORATE(:,:,:,l)                                       ) ! [OUT]

       end do

    case ( 'SN14' )

       do l = 1, ADM_lall
          call ATMOS_PHY_MP_sn14_tendency( &
               KA, KS, KE, IA, IS, IE, JA, JS, JE, &
               DENS(:,:,:,l), W(:,:,:,l), QTRC(:,:,:,QS_MP:QE_MP,l), PRES(:,:,:,l), TEMP(:,:,:,l),            & ! [IN]
               Qdry(:,:,:,l), CPtot(:,:,:,l), CVtot(:,:,:,l), CCN(:,:,:,l), dt_MP, CZ(:,:,:,l), FZ(:,:,:,l),  & ! [IN]
               RHOQ_t(:,:,:,QS_MP:QE_MP), RHOE_t(:,:,:,l), CPtot_t(:,:,:), CVtot_t(:,:,:), EVAPORATE(:,:,:,l) ) ! [OUT]

          do j = JS, JE
          do i = IS, IE
          do k = KS, KE
             do iq = QS_MP, QE_MP
                RHOQ(k,i,j,iq,l) = RHOQ(k,i,j,iq,l) + RHOQ_t(k,i,j,iq) * dt_MP
                QTRC(k,i,j,iq,l) = RHOQ(k,i,j,iq,l) / DENS(k,i,j,l)
             end do
             CPtot(k,i,j,l) = CPtot(k,i,j,l) + CPtot_t(k,i,j) * dt_MP
             CVtot(k,i,j,l) = CVtot(k,i,j,l) + CVtot_t(k,i,j) * dt_MP
             Rtot(k,i,j,l) = CPtot(k,i,j,l) - CVtot_t(k,i,j) ! OK?
          end do
          end do
          end do

       end do

    case ( 'SUZUKI10' )

       KIJMAX = (KE - KE + 1) * IA * JA
       do l = 1, ADM_lall
          call ATMOS_PHY_MP_suzuki10_tendency( KA, KS, KE, IA, IS, IE, JA, JS, JE, KIJMAX, &
                                               dt_MP,                                  & ! [IN]
                                               DENS(:,:,:,l),  PRES(:,:,:,l), TEMP(:,:,:,l), & ! [IN]
                                               QTRC(:,:,:,QS_MP:QE_MP,l), QDRY(:,:,:,l),   & ! [IN]
                                               CPtot(:,:,:,l), CVtot(:,:,:,l),             & ! [IN]
                                               CCN(:,:,:,l),                             & ! [IN] 
                                               RHOQ_t(:,:,:,QS_MP:QE_MP),           & ! [OUT]
                                               RHOE_t(:,:,:,l),                          & ! [OUT]
                                               CPtot_t(:,:,:), CVtot_t(:,:,:),         & ! [OUT]
                                               EVAPORATE(:,:,:,l)                        ) ! [OUT]
          do j = JS, JE
          do i = IS, IE
          do k = KS, KE
             do iq = QS_MP, QE_MP
                RHOQ(k,i,j,iq,l) = RHOQ(k,i,j,iq,l) + RHOQ_t(k,i,j,iq) * dt_MP
                QTRC(k,i,j,iq,l) = RHOQ(k,i,j,iq,l) / DENS(k,i,j,l)
             end do
             CPtot(k,i,j,l) = CPtot(k,i,j,l) + CPtot_t(k,i,j) * dt_MP
             CVtot(k,i,j,l) = CVtot(k,i,j,l) + CVtot_t(k,i,j) * dt_MP
             Rtot(k,i,j,l) = CPtot(k,i,j,l) - CVtot_t(k,i,j)
          end do
          end do
          end do

       end do

    end select


    do l = 1, ADM_lall
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       RHOE(k,i,j,l) = RHOE(k,i,j,l) + RHOE_t(k,i,j,l) * dt_MP
       TEMP(k,i,j,l) = RHOE(k,i,j,l) / ( CVtot(k,i,j,l) * DENS(k,i,j,l) )
       PRES(k,i,j,l) = DENS(k,i,j,l) * Rtot(k,i,j,l) * TEMP(k,i,j,l)
    end do
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
          allocate( vterm_hist(KA,IA,JA,ADM_lall,ih) )
          vterm_hist(:,:,:,:,:) = 0.0_RP
       end if

       do l = 1, ADM_lall

          !$omp parallel do default(none) OMP_SCHEDULE_ collapse(2) &
          !$omp shared (KA,KS,KE,IE,JE,QS_MP,QE_MP,QHA,QLA,QIA,l,ADM_imin,ADM_iall, &
          !$omp         ATMOS_PHY_MP_TYPE, &
          !$omp         dt_MP,MP_NSTEP_SEDIMENTATION,MP_DTSEC_SEDIMENTATION,MP_RNSTEP_SEDIMENTATION, &
          !$omp         CZ,FZ, &
          !$omp         DENS,RHOE,RHOQ,U,V,W,TEMP,PRES,QTRC,CPtot,CVtot,EXNER,rho_bs, &
          !$omp         RHOU,RHOV,MOMZ, &
          !$omp         SFLX_rain,SFLX_snow, &
          !$omp         vterm_hist,hist_vterm_idx) &
          !$omp private(i,j,k,iq,step,ij, &
          !$omp         FDZ,RFDZ,RCDZ, &
          !$omp         RHOU_t,RHOV_t,MOMZ_t, &
          !$omp         REFSTATE_dens, &
          !$omp         vterm,mflux,sflux,eflux,FLX_hydro)
          do j = JS, JE
          do i = IS, IE

             FDZ(KS-1) = CZ(KS,i,j,l) - FZ(KS-1,i,j,l)
             RFDZ(KS-1) = 1.0_RP / FDZ(KS-1)
             do k = KS, KE
                FDZ(k) = CZ(k+1,i,j,l) - CZ(k  ,i,j,l)
                RFDZ(k) = 1.0_RP / FDZ(k)
                RCDZ(k) = 1.0_RP / ( FZ(k  ,i,j,l) - FZ(k-1,i,j,l) )
             enddo

             SFLX_rain(i,j,l) = 0.0_RP
             SFLX_snow(i,j,l) = 0.0_RP
             FLX_hydro(:) = 0.0_RP
             do step = 1, MP_NSTEP_SEDIMENTATION

                select case ( ATMOS_PHY_MP_TYPE )
                case ( 'KESSLER' )
                   do k = KS, KE
                      ij = i + ADM_imin - 1 + ( j - 1 ) * ADM_iall
                      REFSTATE_dens(k) = rho_bs(ij,k,l)
                   end do
                   call ATMOS_PHY_MP_kessler_terminal_velocity( &
                        KA, KS, KE, &
                        DENS(:,i,j,l), RHOQ(:,i,j,:,l), & ! [IN]
                        REFSTATE_dens(:),               & ! [IN]
                        vterm(:,:)                      ) ! [OUT]
                case ( 'TOMITA08' )
                   call ATMOS_PHY_MP_tomita08_terminal_velocity( &
                        KA, KS, KE, &
                        DENS(:,i,j,l), TEMP(:,i,j,l), RHOQ(:,i,j,:,l), & ! [IN]
                        vterm(:,:)                                     ) ! [OUT]
                case ( 'SN14' )
                   call ATMOS_PHY_MP_sn14_terminal_velocity( &
                        KA, KS, KE, &
                        DENS(:,i,j,l), TEMP(:,i,j,l), RHOQ(:,i,j,:,l), PRES(:,i,j,l), & ! [IN]
                        vterm(:,:)                                                    ) ! [OUT]
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
                         vterm_hist(k,i,j,hist_vterm_idx(iq),l) = vterm_hist(k,i,j,hist_vterm_idx(iq),l) &
                                                                + vterm(k,iq) * MP_RNSTEP_SEDIMENTATION
                      end do
                   end if
                end do

                call ATMOS_PHY_MP_precipitation_upwind( &
                     KA, KS, KE, QE_MP-QS_MP, QLA, QIA, &
                     TEMP(:,i,j,l), vterm(:,:),      & ! [IN]
                     FDZ(:), RCDZ(:),                & ! [IN]
                     MP_DTSEC_SEDIMENTATION,         & ! [IN]
                     i, j,                           & ! [IN]
                     DENS(:,i,j,l), RHOQ(:,i,j,:,l), & ! [INOUT]
                     CPtot(:,i,j,l), CVtot(:,i,j,l), & ! [INOUT]
                     RHOE(:,i,j,l),                  & ! [INOUT]
                     mflux(:), sflux(:), eflux       ) ! [OUT]

                do k = KS, KE
                   TEMP(k,i,j,l) = RHOE(k,i,j,l) / ( DENS(k,i,j,l) * CVtot(k,i,j,l) )
                end do

                do k = KS-1, KE-1
                   FLX_hydro(k) = FLX_hydro(k) + mflux(k) * MP_RNSTEP_SEDIMENTATION
                enddo

                SFLX_rain(i,j,l) = SFLX_rain(i,j,l) - sflux(1) * MP_RNSTEP_SEDIMENTATION
                SFLX_snow(i,j,l) = SFLX_snow(i,j,l) - sflux(2) * MP_RNSTEP_SEDIMENTATION

             enddo


             call ATMOS_PHY_MP_precipitation_momentum( &
                  KA, KS, KE, &
                  DENS(:,i,j,l), MOMZ(:,i,j,l), U(:,i,j,l), V(:,i,j,l), & ! [IN]
                  FLX_hydro(:),                                         & ! [IN]
                  RCDZ(:), RFDZ(:),                                     & ! [IN]
                  MOMZ_t(:), RHOU_t(:), RHOV_t(:)                       ) ! [OUT]

             do k = KS, KE
                RHOU(k,i,j,l) = RHOU(k,i,j,l) + RHOU_t(k) * dt_MP
                RHOV(k,i,j,l) = RHOV(k,i,j,l) + RHOV_t(k) * dt_MP
                MOMZ(k,i,j,l) = MOMZ(k,i,j,l) + MOMZ_t(k) * dt_MP
             end do

          enddo
          enddo

       end do

       ! history output
       do iq = QS_MP+1, QE_MP
          if ( hist_vterm_idx(iq) > 0 ) &
               call FILE_HISTORY_put( hist_vterm_id(iq), vterm_hist(:,:,:,:,hist_vterm_idx(iq)) )
       end do
       if ( allocated( vterm_hist ) ) deallocate( vterm_hist )

       call PROF_rapend  ('MP_Precipitation', 2)

    end if

!OCL XFILL
    do l = 1, ADM_lall
    do j = JS, JE
    do i = IS, IE
       precip(i,j,l) = SFLX_rain(i,j,l) + SFLX_snow(i,j,l)
    end do
    end do
    end do

    call FILE_HISTORY_in( SFLX_rain(:,:,:),   'RAIN_MP',   'surface rain rate by MP',          'kg/m2/s',  fill_halo=.true. )
    call FILE_HISTORY_in( SFLX_snow(:,:,:),   'SNOW_MP',   'surface snow rate by MP',          'kg/m2/s',  fill_halo=.true. )
    call FILE_HISTORY_in( precip   (:,:,:),   'PREC_MP',   'surface precipitation rate by MP', 'kg/m2/s',  fill_halo=.true. )
    call FILE_HISTORY_in( EVAPORATE(:,:,:,:), 'EVAPORATE', 'evaporated cloud number',          'num/m3/s', fill_halo=.true. )


    call ATMOS_vars_calc_diagnostics

    return
  end subroutine ATMOS_PHY_MP_driver_step

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
