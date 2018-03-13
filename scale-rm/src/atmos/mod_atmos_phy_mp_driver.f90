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
#include "inc_openmp.h"
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
  public :: ATMOS_PHY_MP_driver_resume
  public :: ATMOS_PHY_MP_driver_calc_tendency
  public :: ATMOS_PHY_MP_driver_adjustment

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
    use scale_process, only: &
       PRC_abort
    use scale_tracer, only: &
       TRACER_regist
    use scale_atmos_hydrometeor, only: &
       ATMOS_HYDROMETEOR_regist, &
       I_QC, &
       I_QR, &
       I_QI, &
       I_QS, &
       I_QG
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
    ! obsolute
    use scale_atmos_phy_mp, only: &
       ATMOS_PHY_MP_config, &
       QA_MP_obsolute => QA_MP, &
       QS_MP_obsolute => QS_MP, &
       QE_MP_obsolute => QE_MP
    implicit none

    integer :: QS2
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '++++++ Module[Tracer Setup] / Categ[ATMOS PHY_MP] / Origin[SCALE-RM]'

    if ( ATMOS_sw_phy_mp ) then
       select case ( ATMOS_PHY_MP_TYPE )
       case ( 'KESSLER' )
          call ATMOS_HYDROMETEOR_regist( &
               QS_MP,                                       & ! [OUT]
               ATMOS_PHY_MP_KESSLER_nwaters,                & ! [IN]
               ATMOS_PHY_MP_KESSLER_nices,                  & ! [IN]
               ATMOS_PHY_MP_KESSLER_tracer_names(:),        & ! [IN]
               ATMOS_PHY_MP_KESSLER_tracer_descriptions(:), & ! [IN]
               ATMOS_PHY_MP_KESSLER_tracer_units(:)         ) ! [IN]
          QA_MP = ATMOS_PHY_MP_KESSLER_ntracers
          I_QC = QS_MP+1
          I_QR = QS_MP+2
       case ( 'TOMITA08' )
          call ATMOS_HYDROMETEOR_regist( &
               QS_MP,                                        & ! [OUT]
               ATMOS_PHY_MP_TOMITA08_nwaters,                & ! [IN]
               ATMOS_PHY_MP_TOMITA08_nices,                  & ! [IN]
               ATMOS_PHY_MP_TOMITA08_tracer_names(:),        & ! [IN]
               ATMOS_PHY_MP_TOMITA08_tracer_descriptions(:), & ! [IN]
               ATMOS_PHY_MP_TOMITA08_tracer_units(:)         ) ! [IN]
          QA_MP = ATMOS_PHY_MP_TOMITA08_ntracers
          I_QC = QS_MP+1
          I_QR = QS_MP+2
          I_QI = QS_MP+3
          I_QS = QS_MP+4
          I_QG = QS_MP+5
       case ( 'SUZUKI10' )
          call ATMOS_PHY_MP_suzuki10_tracer_setup

          call ATMOS_HYDROMETEOR_regist( &
               QS_MP,                                        & ! [OUT]
               ATMOS_PHY_MP_suzuki10_nwaters,                & ! [IN]
               ATMOS_PHY_MP_suzuki10_nices,                  & ! [IN]
               ATMOS_PHY_MP_suzuki10_tracer_names(:),        & ! [IN]
               ATMOS_PHY_MP_suzuki10_tracer_descriptions(:), & ! [IN]
               ATMOS_PHY_MP_suzuki10_tracer_units(:)         ) ! [IN]

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
          call ATMOS_PHY_MP_config( ATMOS_PHY_MP_TYPE )
          QA_MP = QA_MP_obsolute
          QS_MP = QS_MP_obsolute
          ! if ( IO_L ) write(IO_FID_LOG,*) '+++ ATMOS_PHY_MP_TYPE is invalud: ', trim(ATMOS_PHY_MP_TYPE)
          ! call PRC_abort
       end select

       QE_MP = QS_MP + QA_MP - 1

    else
       QA_MP = 0
       QS_MP = -1
       QE_MP = -1
    end if

    ! tentative
    QA_MP_obsolute = QA_MP
    QS_MP_obsolute = QS_MP
    QE_MP_obsolute = QE_MP

    return
  end subroutine ATMOS_PHY_MP_driver_tracer_setup

  !-----------------------------------------------------------------------------
  !> Setup
  subroutine ATMOS_PHY_MP_driver_setup
    use scale_process, only: &
       PRC_abort
    use scale_atmos_grid_cartesC, only: &
       CDZ => ATMOS_GRID_CARTESC_CDZ
    use scale_const, only: &
       EPS => CONST_EPS
    use scale_time, only: &
       TIME_DTSEC_ATMOS_PHY_MP
    use scale_atmos_phy_mp, only: &
       ATMOS_PHY_MP_setup
    use scale_atmos_phy_mp_kessler, only: &
       ATMOS_PHY_MP_KESSLER_setup
    use scale_atmos_phy_mp_tomita08, only: &
       ATMOS_PHY_MP_TOMITA08_setup
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

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '++++++ Module[DRIVER] / Categ[ATMOS PHY_MP] / Origin[SCALE-RM]'

    if ( ATMOS_sw_phy_mp ) then

       MP_cldfrac_thleshold = EPS

       !--- read namelist
       rewind(IO_FID_CONF)
       read(IO_FID_CONF,nml=PARAM_ATMOS_PHY_MP,iostat=ierr)
       if( ierr < 0 ) then !--- missing
          if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
       elseif( ierr > 0 ) then !--- fatal error
          write(*,*) 'xxx Not appropriate names in namelist PARAM_ATMOS_PHY_MP. Check!'
          call PRC_abort
       endif
       if( IO_NML ) write(IO_FID_NML,nml=PARAM_ATMOS_PHY_MP)

       nstep_max = ceiling( ( TIME_DTSEC_ATMOS_PHY_MP * MP_max_term_vel ) / minval( CDZ(:) ) )
       MP_ntmax_sedimentation = max( MP_ntmax_sedimentation, nstep_max )

       MP_NSTEP_SEDIMENTATION  = MP_ntmax_sedimentation
       MP_RNSTEP_SEDIMENTATION = 1.0_RP / real(MP_ntmax_sedimentation,kind=RP)
       MP_DTSEC_SEDIMENTATION  = TIME_DTSEC_ATMOS_PHY_MP * MP_RNSTEP_SEDIMENTATION

       ATMOS_PHY_MP_cldfrac_thleshold = MP_cldfrac_thleshold

       if( IO_L ) write(IO_FID_LOG,*)
       if( IO_L ) write(IO_FID_LOG,*) '*** Enable negative fixer?                    : ', MP_do_negative_fixer
       if( IO_L ) write(IO_FID_LOG,*) '*** Value limit of negative fixer (abs)       : ', abs(MP_limit_negative)
       if( IO_L ) write(IO_FID_LOG,*) '*** Enable sedimentation (precipitation)?     : ', MP_do_precipitation
       if( IO_L ) write(IO_FID_LOG,*) '*** Timestep of sedimentation is divided into : ', MP_ntmax_sedimentation, 'step'
       if( IO_L ) write(IO_FID_LOG,*) '*** DT of sedimentation                       : ', MP_DTSEC_SEDIMENTATION, '[s]'

       select case ( ATMOS_PHY_MP_TYPE )
       case ( 'KESSLER' )
          call ATMOS_PHY_MP_kessler_setup
       case ( 'TOMITA08' )
          call ATMOS_PHY_MP_tomita08_setup( &
               KA, KS, KE, IA, ISB, IEB, JA, JSB, JEB )
       case ( 'SUZUKI10' )
          call ATMOS_PHY_MP_suzuki10_setup( &
               KA, IA, JA )
       case default
          ! setup library component
          call ATMOS_PHY_MP_setup
       end select

    else

       if( IO_L ) write(IO_FID_LOG,*) '*** this component is never called.'
       if( IO_L ) write(IO_FID_LOG,*) '*** SFLX_rain and SFLX_snow is set to zero.'
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
  !> resume
  subroutine ATMOS_PHY_MP_driver_resume
    use mod_atmos_admin, only: &
       ATMOS_sw_phy_mp
    implicit none

    if ( ATMOS_sw_phy_mp ) then

       ! run once (only for the diagnostic value)
       call PROF_rapstart('ATM_Microphysics', 1)
       call ATMOS_PHY_MP_driver_calc_tendency( update_flag = .true. )
       call PROF_rapend  ('ATM_Microphysics', 1)

    end if

    return
  end subroutine ATMOS_PHY_MP_driver_resume


  !-----------------------------------------------------------------------------
  !> adjustment
  subroutine ATMOS_PHY_MP_driver_adjustment
    use scale_atmos_hydrometeor, only: &
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

    if ( MP_do_negative_fixer .and. I_QV > 0) then

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
    use scale_atmos_thermodyn, only: &
       ATMOS_THERMODYN_qd,         &
       ATMOS_THERMODYN_rhoe,       &
       ATMOS_THERMODYN_rhot,       &
       ATMOS_THERMODYN_pott,       &
       ATMOS_THERMODYN_temp_pres,  &
       ATMOS_THERMODYN_temp_pres_E
    use scale_atmos_saturation, only: &
       ATMOS_SATURATION_pres2qsat_liq, &
       ATMOS_SATURATION_pres2qsat_ice
    use scale_atmos_phy_mp, only: &
       ATMOS_PHY_MP
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
    use scale_atmos_phy_mp_suzuki10, only: &
       ATMOS_PHY_MP_suzuki10_nbin, &
       ATMOS_PHY_MP_suzuki10_nbnd, &
       ATMOS_PHY_MP_suzuki10_nspc, &
       ATMOS_PHY_MP_suzuki10_ic,   &
       ATMOS_PHY_MP_suzuki10_ip,   &
       ATMOS_PHY_MP_suzuki10_id,   &
       ATMOS_PHY_MP_suzuki10_iss,  &
       ATMOS_PHY_MP_suzuki10_ig,   &
       ATMOS_PHY_MP_suzuki10_ih,   &
       ATMOS_PHY_MP_suzuki10_adjustment, &
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
       RHOT_t => RHOT_tp, &
       RHOQ_t => RHOQ_tp, &
       RHOH   => RHOH_p, &
       TEMP, &
       PRES, &
       CVtot, &
       CPtot, &
       EXNER, &
       MOMX   => MOMX_av, &
       MOMY   => MOMY_av, &
       RHOT   => RHOT_av, &
       MOMX_t => MOMX_tp, &
       MOMY_t => MOMY_tp
    use mod_atmos_phy_mp_vars, only: &
       QA_MP, &
       QS_MP, &
       QE_MP, &
       DENS_t_MP => ATMOS_PHY_MP_DENS_t,    &
       MOMZ_t_MP => ATMOS_PHY_MP_MOMZ_t,    &
       RHOT_t_MP => ATMOS_PHY_MP_RHOT_t,    &
       RHOU_t_MP => ATMOS_PHY_MP_RHOU_t,    &
       RHOV_t_MP => ATMOS_PHY_MP_RHOV_t,    &
       RHOQ_t_MP => ATMOS_PHY_MP_RHOQ_t,    &
       RHOH_MP   => ATMOS_PHY_MP_RHOH,      &
       EVAPORATE => ATMOS_PHY_MP_EVAPORATE, &
       SFLX_rain => ATMOS_PHY_MP_SFLX_rain, &
       SFLX_snow => ATMOS_PHY_MP_SFLX_snow, &
       MOMX_t_MP => ATMOS_PHY_MP_MOMX_t,    &
       MOMY_t_MP => ATMOS_PHY_MP_MOMY_t
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
    real(RP) :: QDRY  (KA,IA,JA)
    real(RP) :: QSAT_L(KA,IA,JA)
    real(RP) :: QSAT_I(KA,IA,JA)
    real(RP) :: QHYD  (KA,IA,JA,6)
    real(RP) :: vterm (KA,QS_MP+1:QE_MP)
!    real(RP), target :: QTRC1(KA,IA,JA,QS_MP:QA_MP)

    ! obsolute
    real(RP) :: DENS1(KA,IA,JA)
    real(RP) :: MOMZ1(KA,IA,JA)
    real(RP) :: MOMX1(KA,IA,JA)
    real(RP) :: MOMY1(KA,IA,JA)
    real(RP) :: RHOT1(KA,IA,JA)
    real(RP) :: POTT1(KA,IA,JA)
    real(RP), target :: QTRC1(KA,IA,JA,QA)
    logical :: integ_precip = .true.

    real(RP) :: FLX_hydro(KA)
    real(RP) :: DENS2    (KA)
    real(RP) :: TEMP2    (KA)
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

    real(RP) :: CP_t, CV_t

    real(RP) :: precip   (IA,JA)

    ! for history output
    real(RP), allocatable :: vterm_hist(:,:,:,:)
    integer :: hist_vterm_idx(QS_MP+1:QE_MP)
    logical :: flag
    integer :: ih

    integer :: k, i, j, iq
    integer :: step

    ! tentative for suzuki10
    real(RP) :: FLX_hydro_3D(KA,IA,JA,QA-1)
    real(RP) :: pflux       (KA,IA,JA,QA-1)
    real(RP) :: RHOE_3D     (KA,IA,JA)
    real(RP) :: RHOT2       (KA,IA,JA)
    real(RP) :: QTRC2       (KA,IA,JA,QA)
    real(RP) :: vterm_3D    (KA,IA,JA,QS_MP+1:QE_MP)
    integer  :: m, n
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

          do j = JSB, JEB
          do i = ISB, IEB
          do k = KS, KE
             CP_t = ( CPtot1(k,i,j) - CPtot(k,i,j) ) / dt_MP
             CV_t = ( CVtot1(k,i,j) - CVtot(k,i,j) ) / dt_MP
             RHOH_MP(k,i,j) = RHOE_t(k,i,j) &
                  - ( CP_t + log( PRES(k,i,j) / PRE00 ) * ( CVtot(k,i,j) / CPtot(k,i,j) * CP_t - CV_t ) ) &
                  * DENS(k,i,j) * TEMP(k,i,j)
          end do
          end do
          end do

          do iq = QS_MP, QE_MP
          do j = JSB, JEB
          do i = ISB, IEB
          do k = KS, KE
             RHOQ_t_MP(k,i,j,iq) = ( QTRC1(k,i,j,iq) - QTRC(k,i,j,iq) ) * DENS(k,i,j) / dt_MP
          enddo
          enddo
          enddo
          enddo

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

          do j = JSB, JEB
          do i = ISB, IEB
          do k = KS, KE
             CP_t = ( CPtot1(k,i,j) - CPtot(k,i,j) ) / dt_MP
             CV_t = ( CVtot1(k,i,j) - CVtot(k,i,j) ) / dt_MP
             RHOH_MP(k,i,j) = RHOE_t(k,i,j) &
                  - ( CP_t + log( PRES(k,i,j) / PRE00 ) * ( CVtot(k,i,j) / CPtot(k,i,j) * CP_t - CV_t ) ) &
                  * DENS(k,i,j) * TEMP(k,i,j)
!             RHOT_t_MP(k,i,j) = RHOE_t(k,i,j) / ( EXNER(k,i,j) * CPtot(k,i,j) ) &
!                  - RHOT(k,i,j) * CP_t / CPtot(k,i,j) &
!                  - RHOT(k,i,j) * CVtot(k,i,j) / ( CPtot(k,i,j) - CVtot(k,i,j) ) &
!                  * log( EXNER(k,i,j) ) * ( CP_t / CPtot(k,i,j) - CV_t / CVtot(k,i,j) )
          end do
          end do
          end do

          do iq = QS_MP, QE_MP
          do j = JSB, JEB
          do i = ISB, IEB
          do k = KS, KE
             RHOQ_t_MP(k,i,j,iq) = ( QTRC1(k,i,j,iq) - QTRC(k,i,j,iq) ) * DENS(k,i,j) / dt_MP
          enddo
          enddo
          enddo
          enddo

       case ( 'SUZUKI10' )

!OCL XFILL
          QTRC1(:,:,:,:) = QTRC(:,:,:,:) ! save

          call ATMOS_THERMODYN_qd( QDRY(:,:,:),                   & ! [OUT]
                                   QTRC1(:,:,:,:), TRACER_MASS(:) ) ! [IN]

          call ATMOS_THERMODYN_temp_pres( TEMP1(:,:,:),  & ! [OUT]
                                          PRES1(:,:,:),  & ! [OUT]
                                          DENS(:,:,:),   & ! [IN]
                                          RHOT(:,:,:),   & ! [IN]
                                          QTRC1(:,:,:,:),& ! [IN]
                                          TRACER_CV(:),  & ! [IN]
                                          TRACER_R(:),   & ! [IN]
                                          TRACER_MASS(:) ) ! [IN]

          call ATMOS_SATURATION_pres2qsat_liq( KA, KS, KE, & ! [IN]
                                               IA, IS, IE, & ! [IN]
                                               JA, JS, JE, & ! [IN]
                                               TEMP1(:,:,:), PRES1(:,:,:), QDRY(:,:,:), & ! [IN]
                                               QSAT_L(:,:,:)                            ) ! [OUT]

          call ATMOS_SATURATION_pres2qsat_ice( KA, KS, KE, & ! [IN]
                                               IA, IS, IE, & ! [IN]
                                               JA, JS, JE, & ! [IN]
                                               TEMP1(:,:,:), PRES1(:,:,:), QDRY(:,:,:), & ! [IN]
                                               QSAT_I(:,:,:)                            ) ! [OUT]

          call ATMOS_PHY_MP_suzuki10_adjustment( KA, KS,  KE,        & ! [IN]
                                                 IA, ISB, IEB,       & ! [IN]
                                                 JA, JSB, JEB,       & ! [IN]
                                                 KIJMAX,             & ! [IN]
                                                 dt_MP,              & ! [IN]
                                                 DENS     (:,:,:),   & ! [IN]
                                                 PRES1    (:,:,:),   & ! [IN]
                                                 QDRY     (:,:,:),   & ! [IN]
                                                 QSAT_L   (:,:,:),   & ! [IN]
                                                 QSAT_I   (:,:,:),   & ! [IN]
                                                 CCN      (:,:,:),   & ! [IN] 
                                                 TEMP1    (:,:,:),   & ! [INOUT]
                                                 QTRC1    (:,:,:,:), & ! [INOUT]
                                                 EVAPORATE(:,:,:)    ) ! [OUT]

          call ATMOS_THERMODYN_pott( POTT1(:,:,:),   & ! [OUT]
                                     TEMP1(:,:,:),   & ! [IN]
                                     PRES1(:,:,:),   & ! [IN]
                                     QTRC1(:,:,:,:), & ! [IN]
                                     TRACER_CV(:),   & ! [IN]
                                     TRACER_R(:),    & ! [IN]
                                     TRACER_MASS(:)  ) ! [IN]

          do j = JS, JE
          do i = IS, IE
          do k = KS, KE
             RHOT1(k,i,j) = POTT1(k,i,j) * DENS(k,i,j)
          enddo
          enddo
          enddo

!OCL XFILL
          !!$omp parallel do default(none) private(i,j,k) OMP_SCHEDULE_ collapse(2) &
          !!$omp shared(JSB,JEB,ISB,IEB,KS,KE,DENS_t_MP,DENS1,DENS,MOMZ_t_MP,MOMZ1,MOMZ,MOMX_t_MP,MOMX1) &
          !!$omp shared(MOMX,MOMY_t_MP,MOMY1,MOMY,RHOT_t_MP,RHOT1,RHOT,dt_MP)
          !do j = JSB, JEB
          !do i = ISB, IEB
          !do k = KS, KE
          !   RHOT_t_MP(k,i,j) = ( RHOT1(k,i,j) - RHOT(k,i,j) ) / dt_MP
          !enddo
          !enddo
          !enddo
!OCL XFILL!
          !do iq = QS_MP, QE_MP
          !!$omp parallel do default(none) &
          !!$omp shared(JSB,JEB,ISB,IEB,KS,KE,RHOQ_t_MP,iq,QTRC1,QTRC,DENS1,DENS,dt_MP) &
          !!$omp private(i,j,k) OMP_SCHEDULE_ collapse(2)
          !do j = JSB, JEB
          !do i = ISB, IEB
          !do k = KS, KE
          !   RHOQ_t_MP(k,i,j,iq) = ( QTRC1(k,i,j,iq) - QTRC(k,i,j,iq) ) * DENS(k,i,j) / dt_MP
          !enddo
          !enddo
          !enddo
          !enddo

          integ_precip = .false.

       case default

!OCL XFILL
          DENS1(:,:,:) = DENS(:,:,:) ! save
!OCL XFILL
          MOMZ1(:,:,:) = MOMZ(:,:,:) ! save
!OCL XFILL
          MOMX1(:,:,:) = MOMX(:,:,:) ! save
!OCL XFILL
          MOMY1(:,:,:) = MOMY(:,:,:) ! save
!OCL XFILL
          RHOT1(:,:,:) = RHOT(:,:,:) ! save
!OCL XFILL
          QTRC1(:,:,:,:) = QTRC(:,:,:,:) ! save

          call ATMOS_PHY_MP( DENS1    (:,:,:),   & ! [INOUT]
                             MOMZ1    (:,:,:),   & ! [INOUT]
                             MOMX1    (:,:,:),   & ! [INOUT]
                             MOMY1    (:,:,:),   & ! [INOUT]
                             RHOT1    (:,:,:),   & ! [INOUT]
                             QTRC1    (:,:,:,:), & ! [INOUT]
                             CCN      (:,:,:),   & ! [IN]
                             EVAPORATE(:,:,:),   & ! [OUT]
                             SFLX_rain(:,:),     & ! [OUT]
                             SFLX_snow(:,:)      ) ! [OUT]
!OCL XFILL
          !$omp parallel do default(none) private(i,j,k) OMP_SCHEDULE_ collapse(2) &
          !$omp shared(JSB,JEB,ISB,IEB,KS,KE,DENS_t_MP,DENS1,DENS,MOMZ_t_MP,MOMZ1,MOMZ,MOMX_t_MP,MOMX1) &
          !$omp shared(MOMX,MOMY_t_MP,MOMY1,MOMY,RHOT_t_MP,RHOT1,RHOT,dt_MP)
          do j = JSB, JEB
          do i = ISB, IEB
          do k = KS, KE
             DENS_t_MP(k,i,j) = ( DENS1(k,i,j) - DENS(k,i,j) ) / dt_MP
             MOMZ_t_MP(k,i,j) = ( MOMZ1(k,i,j) - MOMZ(k,i,j) ) / dt_MP
             MOMX_t_MP(k,i,j) = ( MOMX1(k,i,j) - MOMX(k,i,j) ) / dt_MP
             MOMY_t_MP(k,i,j) = ( MOMY1(k,i,j) - MOMY(k,i,j) ) / dt_MP
             RHOT_t_MP(k,i,j) = ( RHOT1(k,i,j) - RHOT(k,i,j) ) / dt_MP
          enddo
          enddo
          enddo

!OCL XFILL
          do iq = QS_MP, QE_MP
          !$omp parallel do default(none) &
          !$omp shared(JSB,JEB,ISB,IEB,KS,KE,RHOQ_t_MP,iq,QTRC1,QTRC,DENS1,DENS,dt_MP) &
          !$omp private(i,j,k) OMP_SCHEDULE_ collapse(2)
          do j = JSB, JEB
          do i = ISB, IEB
          do k = KS, KE
             RHOQ_t_MP(k,i,j,iq) = ( QTRC1(k,i,j,iq) * DENS1(k,i,j) &
                                   - QTRC (k,i,j,iq) * DENS (k,i,j) ) / dt_MP
          enddo
          enddo
          enddo
          enddo

          integ_precip = .false.

       end select

       if ( MP_do_precipitation ) then

          call PROF_rapstart('MP_Precipitation', 2)

          if( ATMOS_PHY_MP_TYPE == 'SUZUKI10' ) then
!OCL XFILL
             DENS1(:,:,:) = DENS(:,:,:) ! save
!OCL XFILL
             MOMZ1(:,:,:) = MOMZ(:,:,:) ! save
!OCL XFILL
             MOMX1(:,:,:) = MOMX(:,:,:) ! save
!OCL XFILL
             MOMY1(:,:,:) = MOMY(:,:,:) ! save
!OCL XFILL
             RHOT2(:,:,:) = RHOT1(:,:,:) ! save
!OCL XFILL
             QTRC2(:,:,:,:) = QTRC1(:,:,:,:) ! save

             call ATMOS_PHY_MP_suzuki10_terminal_velocity( &
                  KA,        & ! [IN]
                  vterm(:,:) ) ! [OUT]

             do j = 1, JA
             do i = 1, IA
                vterm_3D(:,i,j,:) = vterm(:,:)
             enddo
             enddo
             FLX_hydro_3D(:,:,:,:) = 0.0_RP

             call ATMOS_THERMODYN_rhoe( RHOE_3D(:,:,:), & ! [OUT]
                                        RHOT2(:,:,:),   & ! [IN]
                                        QTRC2(:,:,:,:), & ! [IN]
                                        TRACER_CV(:),   & ! [IN]
                                        TRACER_R(:),    & ! [IN]
                                        TRACER_MASS(:)  ) ! [IN]

             do step = 1, MP_NSTEP_SEDIMENTATION

                call ATMOS_THERMODYN_temp_pres_E( temp(:,:,:),    & ! [OUT]
                                                  pres(:,:,:),    & ! [OUT]
                                                  DENS(:,:,:),    & ! [IN]
                                                  RHOE_3D(:,:,:), & ! [IN]
                                                  QTRC2(:,:,:,:), & ! [IN]
                                                  TRACER_CV(:),   & ! [IN]
                                                  TRACER_R(:),    & ! [IN]
                                                  TRACER_MASS(:)  ) ! [IN]

                call ATMOS_PHY_MP_precipitation( QA_MP,                 & ! [IN]
                                                 QS_MP,                 & ! [IN]
                                                 pflux    (:,:,:,:),    & ! [OUT]
                                                 vterm_3D (:,:,:,:),    & ! [INOUT]
                                                 DENS1     (:,:,:),     & ! [INOUT]
                                                 MOMZ1     (:,:,:),     & ! [INOUT]
                                                 MOMX1     (:,:,:),     & ! [INOUT]
                                                 MOMY1     (:,:,:),     & ! [INOUT]
                                                 RHOE_3D  (:,:,:),      & ! [INOUT]
                                                 QTRC2    (:,:,:,:),    & ! [INOUT]
                                                 temp     (:,:,:),      & ! [IN]
                                                 TRACER_CV(:),          & ! [IN]
                                                 MP_DTSEC_SEDIMENTATION ) ! [IN]

                do iq = 1, QA-1
                do j  = JS, JE
                do i  = IS, IE
                do k  = KS-1, KE-1
                   FLX_hydro_3D(k,i,j,iq) = FLX_hydro_3D(k,i,j,iq) + pflux(k,i,j,iq) * MP_RNSTEP_SEDIMENTATION
                enddo
                enddo
                enddo
                enddo

             enddo

             call ATMOS_THERMODYN_rhot( RHOT2(:,:,:),   & ! [OUT]
                                        RHOE_3D(:,:,:), & ! [IN]
                                        QTRC2(:,:,:,:), & ! [IN]
                                        TRACER_CV(:),   & ! [IN]
                                        TRACER_R(:),    & ! [IN]
                                        TRACER_MASS(:)  ) ! [IN]
             do j = JSB, JEB
             do i = ISB, IEB
             do k = KS, KE
                DENS_t_MP(k,i,j) = ( DENS1(k,i,j) - DENS(k,i,j) ) / dt_MP
                MOMZ_t_MP(k,i,j) = ( MOMZ1(k,i,j) - MOMZ(k,i,j) ) / dt_MP
                MOMX_t_MP(k,i,j) = ( MOMX1(k,i,j) - MOMX(k,i,j) ) / dt_MP
                MOMY_t_MP(k,i,j) = ( MOMY1(k,i,j) - MOMY(k,i,j) ) / dt_MP
                RHOT_t_MP(k,i,j) = ( RHOT2(k,i,j) - RHOT(k,i,j) ) / dt_MP
                do iq = QS_MP, QE_MP
                   RHOQ_t_MP(k,i,j,iq) = ( QTRC2(k,i,j,iq) * DENS1(k,i,j) &
                                         - QTRC (k,i,j,iq) * DENS (k,i,j) ) / dt_MP
                enddo
             enddo
             enddo
             enddo

             SFLX_rain(:,:) = 0.0_RP
             SFLX_snow(:,:) = 0.0_RP

             !--- lowermost flux is saved for land process
             do n = 1, ATMOS_PHY_MP_suzuki10_nbin
                iq = n

                do j  = JS, JE
                do i  = IS, IE
                   SFLX_rain(i,j) = SFLX_rain(i,j) - FLX_hydro_3D(KS-1,i,j,iq)
                enddo
                enddo
             enddo

             if ( ATMOS_PHY_MP_suzuki10_nspc > 1 ) then
                do m = ATMOS_PHY_MP_suzuki10_ic, ATMOS_PHY_MP_suzuki10_ih
                do n = 1, ATMOS_PHY_MP_suzuki10_nbin
                   iq = (m-1)*ATMOS_PHY_MP_suzuki10_nbin + n

                   do j  = JS, JE
                   do i  = IS, IE
                      SFLX_snow(i,j) = SFLX_snow(i,j) - FLX_hydro_3D(KS-1,i,j,iq)
                   enddo
                   enddo
                enddo
                enddo
             endif

          endif

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
          !$omp         integ_precip, &
          !$omp         vterm_hist,hist_vterm_idx) &
          !$omp private(i,j,k,iq,step, &
          !$omp         FDZ,RFDZ,RCDZ, &
          !$omp         DENS2,TEMP2,CPtot2,CVtot2,RHOE,RHOE2,RHOQ2, &
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
                !case ( 'SUZUKI10' )
                !   call ATMOS_PHY_MP_suzuki10_terminal_velocity( &
                !        KA,        & ! [IN]
                !        vterm(:,:) ) ! [OUT]
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
                     KA, KS, KE, QHA, QLA, QIA, &
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

             if ( integ_precip ) then ! tentative

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

             end if

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

       ! history output for QHYD when using SUZUKI10
       if( ATMOS_PHY_MP_TYPE == 'SUZUKI10' ) then
          QHYD(:,:,:,:) = 0.0_RP

          do iq = 1, ATMOS_PHY_MP_suzuki10_nbnd
             do j = JSB, JEB
             do i = ISB, IEB
             do k = KS, KE
                QHYD(k,i,j,1) = QHYD(k,i,j,1) + QTRC1(k,i,j,iq+1)
             enddo
             enddo
             enddo
          enddo

          do iq = ATMOS_PHY_MP_suzuki10_nbnd+1, ATMOS_PHY_MP_suzuki10_nbin
             do j = JSB, JEB
             do i = ISB, IEB
             do k = KS, KE
                QHYD(k,i,j,2) = QHYD(k,i,j,2) + QTRC1(k,i,j,iq+1)
             enddo
             enddo
             enddo
          enddo

          call FILE_HISTORY_in( QHYD(:,:,:,1), 'QC', 'Mixing ratio of QC', 'kg/kg' )
          call FILE_HISTORY_in( QHYD(:,:,:,2), 'QR', 'Mixing ratio of QR', 'kg/kg' )

          if ( ATMOS_PHY_MP_suzuki10_nspc > 1 ) then
             do iq = 1, ATMOS_PHY_MP_suzuki10_nbin
                do j = JSB, JEB
                do i = ISB, IEB
                do k = KS, KE
                   ! columnar,plate,dendrite = ice
                   QHYD(k,i,j,3) = QHYD(k,i,j,3) + QTRC1(k,i,j,1+(ATMOS_PHY_MP_suzuki10_ic-1 )*ATMOS_PHY_MP_suzuki10_nbin+iq)
                   QHYD(k,i,j,3) = QHYD(k,i,j,3) + QTRC1(k,i,j,1+(ATMOS_PHY_MP_suzuki10_ip-1 )*ATMOS_PHY_MP_suzuki10_nbin+iq)
                   QHYD(k,i,j,3) = QHYD(k,i,j,3) + QTRC1(k,i,j,1+(ATMOS_PHY_MP_suzuki10_id-1 )*ATMOS_PHY_MP_suzuki10_nbin+iq)
                   QHYD(k,i,j,4) = QHYD(k,i,j,4) + QTRC1(k,i,j,1+(ATMOS_PHY_MP_suzuki10_iss-1)*ATMOS_PHY_MP_suzuki10_nbin+iq)
                   QHYD(k,i,j,5) = QHYD(k,i,j,5) + QTRC1(k,i,j,1+(ATMOS_PHY_MP_suzuki10_ig-1 )*ATMOS_PHY_MP_suzuki10_nbin+iq)
                   QHYD(k,i,j,6) = QHYD(k,i,j,6) + QTRC1(k,i,j,1+(ATMOS_PHY_MP_suzuki10_ih-1 )*ATMOS_PHY_MP_suzuki10_nbin+iq)
                enddo
                enddo
                enddo
             enddo

             call FILE_HISTORY_in( QHYD(:,:,:,3), 'QI', 'Mixing ratio of QI', 'kg/kg' )
             call FILE_HISTORY_in( QHYD(:,:,:,4), 'QS', 'Mixing ratio of QS', 'kg/kg' )
             call FILE_HISTORY_in( QHYD(:,:,:,5), 'QG', 'Mixing ratio of QG', 'kg/kg' )
             call FILE_HISTORY_in( QHYD(:,:,:,6), 'QH', 'Mixing ratio of QH', 'kg/kg' )
          endif

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
       call FILE_HISTORY_in( RHOT_t_MP(:,:,:), 'RHOT_t_MP', 'tendency RHOT in MP',         'K*kg/m3/s', fill_halo=.true. )
       call FILE_HISTORY_in( RHOH_MP  (:,:,:), 'RHOH_MP',   'diabatic heating rate in MP', 'J/m3/s',    fill_halo=.true. )
       call FILE_HISTORY_in( MOMX_t_MP(:,:,:), 'MOMX_t_MP', 'tendency MOMX in MP',         'kg/m2/s2' , fill_halo=.true. )
       call FILE_HISTORY_in( MOMY_t_MP(:,:,:), 'MOMY_t_MP', 'tendency MOMY in MP',         'kg/m2/s2' , fill_halo=.true. )

       do iq = QS_MP, QE_MP
          call FILE_HISTORY_in( RHOQ_t_MP(:,:,:,iq), trim(TRACER_NAME(iq))//'_t_MP', &
                        'tendency rho*'//trim(TRACER_NAME(iq))//' in MP', 'kg/m3/s', fill_halo=.true. )
       enddo

    endif

    !$omp parallel do default(none) private(i,j,k) OMP_SCHEDULE_ collapse(2) &
    !$omp shared(KS,KE,ISB,IEB,JSB,JEB, &
    !$omp        DENS_t_MP,MOMZ_t_MP,RHOU_t_MP,RHOV_t_MP,RHOT_t_MP,RHOH_MP,MOMX_t_MP,MOMY_t_MP, &
    !$omp        DENS_t,MOMZ_t,RHOU_t,RHOV_t,RHOT_t,RHOH,MOMX_t,MOMY_t)
    do j = JSB, JEB
    do i = ISB, IEB
    do k = KS, KE
       DENS_t(k,i,j) = DENS_t(k,i,j) + DENS_t_MP(k,i,j)
       MOMZ_t(k,i,j) = MOMZ_t(k,i,j) + MOMZ_t_MP(k,i,j)
       RHOU_t(k,i,j) = RHOU_t(k,i,j) + RHOU_t_MP(k,i,j)
       RHOV_t(k,i,j) = RHOV_t(k,i,j) + RHOV_t_MP(k,i,j)
       RHOT_t(k,i,j) = RHOT_t(k,i,j) + RHOT_t_MP(k,i,j)
       RHOH  (k,i,j) = RHOH  (k,i,j) + RHOH_MP  (k,i,j)
       MOMX_t(k,i,j) = MOMX_t(k,i,j) + MOMX_t_MP(k,i,j)
       MOMY_t(k,i,j) = MOMY_t(k,i,j) + MOMY_t_MP(k,i,j)
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
                              MOMX_t_MP(:,:,:), 'MOMX_t_MP',         &
                              ATMOS_GRID_CARTESC_REAL_VOLZUY(:,:,:), &
                              ATMOS_GRID_CARTESC_REAL_TOTVOLZUY      )
       call STATISTICS_total( KA, KS, KE, IA, IS, IE, JA, JS, JE, &
                              MOMY_t_MP(:,:,:), 'MOMY_t_MP',         &
                              ATMOS_GRID_CARTESC_REAL_VOLZXV(:,:,:), &
                              ATMOS_GRID_CARTESC_REAL_TOTVOLZXV      )
       call STATISTICS_total( KA, KS, KE, IA, IS, IE, JA, JS, JE, &
                              RHOT_t_MP(:,:,:), 'RHOT_t_MP',      &
                              ATMOS_GRID_CARTESC_REAL_VOL(:,:,:), &
                              ATMOS_GRID_CARTESC_REAL_TOTVOL      )

       do iq = QS_MP, QE_MP
          call STATISTICS_total( KA, KS, KE, IA, IS, IE, JA, JS, JE, &
                                 RHOQ_t_MP(:,:,:,iq), trim(TRACER_NAME(iq))//'_t_MP', &
                                 ATMOS_GRID_CARTESC_REAL_VOL(:,:,:),                  &
                                 ATMOS_GRID_CARTESC_REAL_TOTVOL                       )
       enddo
    endif

    return
  end subroutine ATMOS_PHY_MP_driver_calc_tendency

end module mod_atmos_phy_mp_driver
