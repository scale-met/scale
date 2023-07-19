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
  public :: ATMOS_PHY_MP_driver_finalize
  public :: ATMOS_PHY_MP_driver_calc_tendency
  public :: ATMOS_PHY_MP_driver_adjustment
  public :: ATMOS_PHY_MP_driver_qhyd2qtrc
  public :: ATMOS_PHY_MP_driver_qhyd2qtrc_onlyqv

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
  logical,  private :: MP_do_precipitation    !> apply sedimentation (precipitation)?
  logical,  private :: MP_do_negative_fixer   !> apply negative fixer?
  real(RP), private :: MP_limit_negative      !> Abort if abs(fixed negative vaue) > abs(MP_limit_negative)
  integer,  private :: MP_ntmax_sedimentation !> number of time step for sedimentation
  real(RP), private :: MP_max_term_vel        !> terminal velocity for calculate dt of sedimentation
  real(RP), private :: MP_cldfrac_thleshold   !> thleshold for cloud fraction
  integer,  private :: MP_NSTEP_SEDIMENTATION
  real(RP), private :: MP_RNSTEP_SEDIMENTATION
  real(RP), private :: MP_DTSEC_SEDIMENTATION

  integer,  private, allocatable :: HIST_hyd_id(:)
  integer,  private, allocatable :: HIST_crg_id(:)

  integer, private, allocatable :: hist_vterm_id(:)
  integer, private              :: hist_nf_rhoh_id
  integer, private              :: hist_nf_dens_id
  integer, private              :: hist_nf_engi_id
  integer, private              :: monit_nf_mass_id
  integer, private              :: monit_nf_engi_id
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
       ATMOS_sw_phy_mp, &
       ATMOS_sw_phy_lt
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
          if ( ATMOS_sw_phy_lt ) then
             LOG_ERROR("ATMOS_PHY_MP_driver_tracer_setup",*) 'ATMOS_PHY_MP_TYPE is invalud for Lightning: ', trim(ATMOS_PHY_MP_TYPE)
             LOG_ERROR("ATMOS_PHY_MP_driver_tracer_setup",*) 'ATMOS_PHY_MP_TYPE is should be TOMITA08 or SN14 or SUZUKI10 for Lightning: '
             call PRC_abort
          endif
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
       QE_MP = -2
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
    use scale_atmos_hydrometeor, only: &
       N_HYD
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
    use scale_monitor, only: &
       MONITOR_reg, &
       MONITOR_put
    use mod_atmos_admin, only: &
       ATMOS_PHY_MP_TYPE, &
       ATMOS_sw_phy_mp
    use mod_atmos_phy_mp_vars, only: &
       ATMOS_PHY_MP_cldfrac_thleshold, &
       SFLX_rain => ATMOS_PHY_MP_SFLX_rain, &
       SFLX_snow => ATMOS_PHY_MP_SFLX_snow, &
       SFLX_ENGI => ATMOS_PHY_MP_SFLX_ENGI
    use mod_atmos_phy_mp_vars, only: &
       QS_MP, &
       QE_MP
    use mod_atmos_phy_lt_vars, only: &
       flg_lt
    implicit none

    namelist / PARAM_ATMOS_PHY_MP / &
       MP_do_precipitation,     &
       MP_do_negative_fixer,    &
       MP_limit_negative,       &
       MP_ntmax_sedimentation,  &
       MP_max_term_vel,         &
       MP_cldfrac_thleshold

    character(len=H_SHORT) :: w_name(N_HYD)
    character(len=H_MID)   :: w_longname(N_HYD)
    character(len=H_SHORT) :: w_unit(N_HYD)
    data w_name / 'QC', &
                  'QR', &
                  'QI', &
                  'QS', &
                  'QG', &
                  'QH'  /
    data w_longname / &
                  'Ratio of Cloud Water mass to total mass', &
                  'Ratio of Rain Water mass to total mass', &
                  'Ratio of Ice Water mass to total mass', &
                  'Ratio of Snow Water mass to total mass', &
                  'Ratio of Graupel Water mass to total mass', &
                  'Ratio of Hail Water mass to total mass'  /
    data w_unit / &
                  'kg/kg', &
                  'kg/kg', &
                  'kg/kg', &
                  'kg/kg', &
                  'kg/kg', &
                  'kg/kg'  /

    character(len=H_SHORT) :: w_crg_name(N_HYD)
    character(len=H_MID)   :: w_crg_longname(N_HYD)
    character(len=H_SHORT) :: w_crg_unit(N_HYD)
    data w_crg_name / 'QCRG_C', &
                      'QCRG_R', &
                      'QCRG_I', &
                      'QCRG_S', &
                      'QCRG_G', &
                      'QCRG_H'  /
    data w_crg_longname / &
                  'Ratio of charge density of Cloud Water', &
                  'Ratio of charge density of Rain Water', &
                  'Ratio of charge density of Ice', &
                  'Ratio of charge density of Snow', &
                  'Ratio of charge density of Graupel', &
                  'Ratio of charge density of Hail'  /
    data w_crg_unit / &
                  'fC/kg', &
                  'fC/kg', &
                  'fC/kg', &
                  'fC/kg', &
                  'fC/kg', &
                  'fC/kg'  /

    real(RP) :: ZERO(KA,IA,JA)

    integer :: nstep_max

    integer :: iq
    integer :: ierr
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("ATMOS_PHY_MP_driver_setup",*) 'Setup'

    if ( ATMOS_sw_phy_mp ) then

       MP_cldfrac_thleshold = EPS

       MP_do_precipitation    = .true.
       MP_do_negative_fixer   = .true.
       MP_limit_negative      = 0.1_RP
       MP_ntmax_sedimentation = 1
       MP_max_term_vel        = 10.0_RP

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
                 KA, KS, KE, IA, IS, IE, JA, JS, JE, &
                 flg_lt )
       case ( 'SN14' )
          call ATMOS_PHY_MP_sn14_setup( &
               KA, IA, JA )
       case ( 'SUZUKI10' )
          call ATMOS_PHY_MP_suzuki10_setup( &
                  KA, IA, JA, &
                  flg_lt )
          allocate( HIST_hyd_id(N_HYD) )
          do iq = 1, N_HYD
             call FILE_HISTORY_reg( w_name(iq), w_longname(iq), w_unit(iq), & ! [IN]
                                    HIST_hyd_id(iq)                         ) ! [OUT]
          enddo
          if( flg_lt ) then
             allocate( HIST_crg_id(N_HYD) )
             do iq = 1, N_HYD
                call FILE_HISTORY_reg( w_crg_name(iq), w_crg_longname(iq), w_crg_unit(iq), & ! [IN]
                                       HIST_crg_id(iq)                                     ) ! [OUT]
             enddo
          endif
       end select

       ! history putput
       if ( MP_do_precipitation ) then

          allocate( hist_vterm_id(QS_MP+1:QE_MP) )
          do iq = QS_MP+1, QE_MP
             call FILE_HISTORY_reg( 'Vterm_'//trim(TRACER_NAME(iq)), 'terminal velocity of '//trim(TRACER_NAME(iq)), 'm/s', hist_vterm_id(iq) )
          end do

       else

          !$acc kernels
          SFLX_rain(:,:) = 0.0_RP
          !$acc end kernels
          !$acc kernels
          SFLX_snow(:,:) = 0.0_RP
          !$acc end kernels
          !$acc kernels
          SFLX_ENGI(:,:) = 0.0_RP
          !$acc end kernels

       end if

       ! monitor
       if ( MP_do_negative_fixer ) then
          call FILE_HISTORY_reg( "RHOH_MP_NF",   "sensible heat by the negative fixer",          "J/m3/s", & ! [IN]
                                 hist_nf_rhoh_id                                                            ) ! [OUT]
          call FILE_HISTORY_reg( "DENS_t_MP_NF", "vapor supply by the negative fixer",           "kg/m3/s", & ! [IN]
                                 hist_nf_dens_id                                                            ) ! [OUT]
          call FILE_HISTORY_reg( "ENGI_t_MP_NF", "internal energy supply by the negative fixer", "J/m3/s",  & ! [IN]
                                 hist_nf_engi_id                                                            ) ! [OUT]
          call MONITOR_reg( "QTOTTND_NF", "vapor supply by the negative fixer", "kg",           & ! [IN]
                            monit_nf_mass_id,                                                   & ! [OUT]
                            is_tendency=.true.                                                  ) ! [IN]
          call MONITOR_reg( "ENGITND_NF", "internal energy supply by the negative fixer", "J",  & ! [IN]
                            monit_nf_engi_id,                                                   & ! [OUT]
                            is_tendency=.true.                                                  ) ! [IN]
          ZERO(:,:,:) = 0.0_RP
          call MONITOR_put( MONIT_nf_mass_id, ZERO(:,:,:) )
          call MONITOR_put( MONIT_nf_engi_id, ZERO(:,:,:) )
       end if

    else

       LOG_INFO("ATMOS_PHY_MP_driver_setup",*) 'this component is never called.'
       LOG_INFO("ATMOS_PHY_MP_driver_setup",*) 'SFLX_rain and SFLX_snow is set to zero.'
       !$acc kernels
       SFLX_rain(:,:) = 0.0_RP
       !$acc end kernels
       !$acc kernels
       SFLX_snow(:,:) = 0.0_RP
       !$acc end kernels
       !$acc kernels
       SFLX_ENGI(:,:) = 0.0_RP
       !$acc end kernels

    endif

    return
  end subroutine ATMOS_PHY_MP_driver_setup

  !-----------------------------------------------------------------------------
  !> finalize
  subroutine ATMOS_PHY_MP_driver_finalize
    use mod_atmos_admin, only: &
       ATMOS_PHY_MP_TYPE, &
       ATMOS_sw_phy_mp
    use scale_atmos_phy_mp_tomita08, only: &
       ATMOS_PHY_MP_Tomita08_finalize
    use scale_atmos_phy_mp_sn14, only: &
       ATMOS_PHY_MP_sn14_finalize
    use scale_atmos_phy_mp_suzuki10, only: &
       ATMOS_PHY_MP_suzuki10_finalize
    implicit none
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("ATMOS_PHY_MP_driver_finalize",*) 'Finalize'

    if ( ATMOS_sw_phy_mp ) then
       select case ( ATMOS_PHY_MP_TYPE )
       case ( 'KESSLER' )
       case ( 'TOMITA08' )
          call ATMOS_PHY_MP_Tomita08_finalize
       case ( 'SN14' )
          call ATMOS_PHY_MP_sn14_finalize
       case ( 'SUZUKI10' )
          call ATMOS_PHY_MP_suzuki10_finalize
       end select
    end if

    if ( allocated(HIST_hyd_id) ) deallocate( HIST_hyd_id )
    if ( allocated(HIST_crg_id) ) deallocate( HIST_crg_id )
    if ( allocated(hist_vterm_id) ) deallocate( hist_vterm_id )

    return
  end subroutine ATMOS_PHY_MP_driver_finalize
  !-----------------------------------------------------------------------------
  !> adjustment
  subroutine ATMOS_PHY_MP_driver_adjustment
    use scale_const, only: &
       PRE00 => CONST_PRE00
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
       RHOT,  &
       QTRC
    use mod_atmos_phy_mp_vars, only: &
       QE_MP
    use scale_time, only: &
       dt => TIME_DTSEC
    use scale_file_history, only: &
       FILE_HISTORY_query, &
       FILE_HISTORY_put
    use scale_monitor, only: &
       MONITOR_put

    real(RP) :: RHOH  (KA,IA,JA)
    real(RP) :: DENS_d(KA,IA,JA)
    real(RP) :: ENGI_d(KA,IA,JA)
    real(RP) :: Rtot

    logical :: do_put_rhoh
    logical :: do_put_dens
    logical :: do_put_engi

    integer :: k, i, j, iq

    if ( MP_do_negative_fixer .and. (.not. ATMOS_HYDROMETEOR_dry) ) then

       !$acc data copy(DENS, TEMP, CVtot, CPtot, QTRC(:,:,:,I_QV:QE_MP)) create(RHOH, DENS_d, ENGI_d) copyout(RHOT)

       call FILE_HISTORY_query( hist_nf_rhoh_id, do_put_rhoh )
       call FILE_HISTORY_query( hist_nf_dens_id, do_put_dens )
       call FILE_HISTORY_query( hist_nf_engi_id, do_put_engi )

       if ( monit_nf_mass_id > 0 .or. monit_nf_engi_id > 0 .or. &
            do_put_rhoh .or. do_put_dens .or. do_put_engi ) then
          call ATMOS_PHY_MP_negative_fixer( &
               KA, KS, KE, IA, IS, IE, JA, JS, JE, QLA, QIA, &
               MP_limit_negative,                     & ! [IN]
               DENS(:,:,:), TEMP(:,:,:),              & ! [INOUT]
               CVtot(:,:,:), CPtot(:,:,:),            & ! [INOUT]
               QTRC(:,:,:,I_QV), QTRC(:,:,:,QHS:QHE), & ! [INOUT]
               RHOH = RHOH,                           & ! [OUT, optional]
               DENS_diff = DENS_d, ENGI_diff = ENGI_d ) ! [OUT, optional]
       else
          call ATMOS_PHY_MP_negative_fixer( &
               KA, KS, KE, IA, IS, IE, JA, JS, JE, QLA, QIA, &
               MP_limit_negative,                    & ! [IN]
               DENS(:,:,:), TEMP(:,:,:),             & ! [INOUT]
               CVtot(:,:,:), CPtot(:,:,:),           & ! [INOUT]
               QTRC(:,:,:,I_QV), QTRC(:,:,:,QHS:QHE) ) ! [INOUT]
       end if

       !$omp parallel private(Rtot)

       !$omp do
       !$acc kernels async
       !$acc loop collapse(3) independent
       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          Rtot = CPtot(k,i,j) - CVtot(k,i,j)
          RHOT(k,i,j) = PRE00 / Rtot * ( DENS(k,i,j) * TEMP(k,i,j) * Rtot / PRE00 )**( CVtot(k,i,j) / CPtot(k,i,j) )
       end do
       end do
       end do
       !$acc end kernels
       !$omp end do nowait

       ! for non-mass tracers, such as number density
       !$acc kernels async
       do iq = QHE+1, QE_MP
       !$omp do
       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          QTRC(k,i,j,iq) = max( QTRC(k,i,j,iq), 0.0_RP )
       end do
       end do
       end do
       end do
       !$acc end kernels

       !$omp end parallel

       if ( do_put_rhoh ) then
          !$omp parallel do
          !$acc kernels async
          do j = JS, JE
          do i = IS, IE
          do k = KS, KE
             RHOH(k,i,j) = RHOH(k,i,j) / dt
          end do
          end do
          end do
          !$acc end kernels
          !$acc update host(RHOH)
          call FILE_HISTORY_put( hist_nf_rhoh_id, RHOH(:,:,:) )
       end if
       if ( monit_nf_mass_id > 0 .or. do_put_dens ) then
          !$omp parallel do
          !$acc kernels async
          do j = JS, JE
          do i = IS, IE
          do k = KS, KE
             DENS_d(k,i,j) = DENS_d(k,i,j) / dt
          end do
          end do
          end do
          !$acc end kernels
          !$acc update host(DENS_d)
          call FILE_HISTORY_put( hist_nf_dens_id, DENS_d(:,:,:) )
          call MONITOR_put( monit_nf_mass_id, DENS_d(:,:,:) )
       end if
       if ( monit_nf_engi_id > 0 .or. do_put_engi ) then
          !$omp parallel do
          !$acc kernels async
          do j = JS, JE
          do i = IS, IE
          do k = KS, KE
             ENGI_d(k,i,j) = ENGI_d(k,i,j) / dt
          end do
          end do
          end do
          !$acc end kernels
          !$acc update host(ENGI_d)
          call FILE_HISTORY_put( hist_nf_engi_id, ENGI_d(:,:,:) )
          call MONITOR_put( monit_nf_engi_id, ENGI_d(:,:,:) )
       end if

       !$acc wait

       !$acc end data

    end if

    return
  end subroutine ATMOS_PHY_MP_driver_adjustment

  !-----------------------------------------------------------------------------
  !> calculate tendency
  subroutine ATMOS_PHY_MP_driver_calc_tendency( update_flag )
    use scale_const, only: &
       PRE00 => CONST_PRE00
    use scale_time, only: &
       dt_MP => TIME_DTSEC_ATMOS_PHY_MP, &
       dt_LT => TIME_DTSEC_ATMOS_PHY_LT
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
       FILE_HISTORY_in, &
       FILE_HISTORY_query, &
       FILE_HISTORY_put
    use scale_atmos_hydrometeor, only: &
       LHF,   &
       I_QV,  &
       N_HYD, &
       QHA,   &
       QHS,   &
       QHE,   &
       QLA,   &
       QLS,   &
       QLE,   &
       QIA
    use scale_atmos_phy_mp_common, only: &
       ATMOS_PHY_MP_precipitation_upwind,  &
       ATMOS_PHY_MP_precipitation_semilag, &
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
       ATMOS_PHY_MP_suzuki10_terminal_velocity, &
       ATMOS_PHY_MP_suzuki10_qtrc2qhyd, &
       ATMOS_PHY_MP_suzuki10_crg_qtrc2qhyd
    use scale_atmos_phy_lt_sato2019, only: &
       ATMOS_PHY_LT_sato2019_select_dQCRG_from_LUT
    use scale_file_history, only: &
       FILE_HISTORY_query, &
       FILE_HISTORY_put
    use mod_atmos_admin, only: &
       ATMOS_PHY_MP_TYPE,    &
       ATMOS_PHY_PRECIP_TYPE, &
       ATMOS_sw_phy_lt
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
       QS_MP, &
       QE_MP, &
       DENS_t_MP => ATMOS_PHY_MP_DENS_t,    &
       MOMZ_t_MP => ATMOS_PHY_MP_MOMZ_t,    &
!       RHOT_t_MP => ATMOS_PHY_MP_RHOT_t,    &
       RHOU_t_MP => ATMOS_PHY_MP_RHOU_t,    &
       RHOV_t_MP => ATMOS_PHY_MP_RHOV_t,    &
       RHOQ_t_MP => ATMOS_PHY_MP_RHOQ_t,    &
       RHOC_t_MP => ATMOS_PHY_MP_RHOC_t,    &
       RHOH_MP   => ATMOS_PHY_MP_RHOH,      &
       EVAPORATE => ATMOS_PHY_MP_EVAPORATE, &
       SFLX_rain => ATMOS_PHY_MP_SFLX_rain, &
       SFLX_snow => ATMOS_PHY_MP_SFLX_snow, &
       SFLX_ENGI => ATMOS_PHY_MP_SFLX_ENGI
    use mod_atmos_phy_ae_vars, only: &
       CCN_t => ATMOS_PHY_AE_CCN_t
    use mod_atmos_phy_lt_vars, only: &
       QA_LT, &
       QS_LT, &
       QE_LT, &
       d0_crg, &
       v0_crg, &
       flg_lt, &
       Sarea => ATMOS_PHY_LT_Sarea
    implicit none

    logical, intent(in) :: update_flag

    real(RP) :: RHOE_t(KA,IA,JA)
    real(RP) :: TEMP1 (KA,IA,JA)
    real(RP) :: CPtot1(KA,IA,JA)
    real(RP) :: CVtot1(KA,IA,JA)
    real(RP) :: CCN   (KA,IA,JA)
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
    real(RP) :: RHOQ     (KA,QS_MP+1:QE_MP)
    real(RP) :: RHOQ2    (KA,QS_MP+1:QE_MP)
    real(RP) :: mflux    (KA)
    real(RP) :: sflux    (2)  !> 1: rain, 2: snow
    real(RP) :: eflux

    real(RP) :: FZ  (KA)
    real(RP) :: FDZ (KA)
    real(RP) :: RFDZ(KA)
    real(RP) :: RCDZ(KA)

    real(RP) :: CPtot_t(KA,IA,JA), CVtot_t(KA,IA,JA)
    real(RP) :: CP_t, CV_t

    real(RP) :: precip   (IA,JA)

    real(RP) :: Qe(KA,IA,JA,N_HYD)
    logical  :: HIST_sw(N_HYD)
    real(RP) :: Qecrg(KA,IA,JA,N_HYD)
    logical  :: HIST_crg_sw(N_HYD)

    ! for history output
    real(RP), allocatable :: vterm_hist(:,:,:,:)
    integer :: hist_vterm_idx(QS_MP+1:QE_MP)
    logical :: flag
    integer :: ih

    integer :: k, i, j, iq
    integer :: step

    real(RP) :: QTRC1_crg(KA,IA,JA,QS_LT:QE_LT)
    real(RP) :: RHOQ2_crg(KA,QS_LT:QE_LT)
    real(RP) :: mflux_crg(KA), sflux_crg(2), eflux_crg
    real(RP) :: QSPLT_in(KA,IA,JA,3)
    real(RP) :: dqcrg(KA,IA,JA), beta_crg(KA,IA,JA)
    !---------------------------------------------------------------------------

    !$acc data copy(DENS_t, MOMZ_t, RHOU_t, RHOV_t, RHOQ_t(:,:,:,:), RHOH, &
    !$acc           DENS_t_MP, MOMZ_t_MP, RHOU_t_MP, RHOV_t_MP, RHOQ_t_MP, RHOC_t_MP, RHOH_MP, &
    !$acc           EVAPORATE, SFLX_rain, SFLX_snow, SFLX_ENGI, &
    !$acc           Sarea) &
    !$acc      copyin(DENS(:,:,:), MOMZ(:,:,:), U, V, TEMP, QTRC(:,:,:,:), PRES, CVtot, CPtot, &
    !$acc             REAL_CZ, REAL_FZ, &
    !$acc             ATMOS_PHY_MP_TYPE, ATMOS_PHY_PRECIP_TYPE) &
    !$acc      create(RHOE_t, TEMP1, CPtot1, CVtot1, CCN, QTRC1, CPtot_t, CVtot_t, precip, &
    !$acc             QTRC1_crg, QSPLT_in, dqcrg, beta_crg, &
    !$acc             hist_vterm_idx, vterm_hist)


    if ( update_flag ) then

       !$acc kernels
       CCN(:,:,:) = CCN_t(:,:,:) * dt_MP
       !$acc end kernels

       select case ( ATMOS_PHY_MP_TYPE )
       case ( 'KESSLER' )
          !$omp workshare
          !$acc kernels
!OCL XFILL
          TEMP1(:,:,:) = TEMP(:,:,:)
!OCL XFILL
          QTRC1(:,:,:,QS_MP:QE_MP) = QTRC(:,:,:,QS_MP:QE_MP)
!OCL XFILL
          CVtot1(:,:,:) = CVtot(:,:,:)
!OCL XFILL
          CPtot1(:,:,:) = CPtot(:,:,:)
          !$acc end kernels
          !$omp end workshare

          !$acc update host(DENS,PRES,TEMP1,QTRC1(:,:,:,QS_MP:QE_MP),CPtot1,CVtot1)
          call ATMOS_PHY_MP_kessler_adjustment( &
               KA, KS, KE, IA, IS, IE, JA, JS, JE, &
               DENS(:,:,:), PRES(:,:,:), dt_MP,                                      & ! [IN]
               TEMP1(:,:,:), QTRC1(:,:,:,QS_MP:QE_MP), CPtot1(:,:,:), CVtot1(:,:,:), & ! [INOUT]
               RHOE_t(:,:,:), EVAPORATE(:,:,:)                                       ) ! [OUT]
          !acc update device(TEMP1,QTRC1(:,:,:,QS_MP:QE_MP),CPtot1,CVtot1,RHOE_t,EVAPORATE)

          !$omp parallel do collapse(3)
          !$acc kernels
          do iq = QS_MP, QE_MP
          do j = JS, JE
          do i = IS, IE
          do k = KS, KE
             RHOQ_t_MP(k,i,j,iq) = ( QTRC1(k,i,j,iq) - QTRC(k,i,j,iq) ) * DENS(k,i,j) / dt_MP
          enddo
          enddo
          enddo
          enddo
          !$acc end kernels

          !$omp parallel do collapse(2)
          !$acc kernels
          do j = JS, JE
          do i = IS, IE
          do k = KS, KE
             CPtot_t(k,i,j) = ( CPtot1(k,i,j) - CPtot(k,i,j) ) / dt_MP
             CVtot_t(k,i,j) = ( CVtot1(k,i,j) - CVtot(k,i,j) ) / dt_MP
          end do
          end do
          end do
          !$acc end kernels

       case ( 'TOMITA08' )
!OCL XFILL
          !$omp parallel do collapse(2)
          !$acc parallel async
          !$acc loop collapse(2)
          do j = JS, JE
          do i = IS, IE
          do k = KS, KE
             TEMP1(k,i,j) = TEMP(k,i,j)
          end do
          end do
          end do
          !$acc end parallel
!OCL XFILL
          !$omp parallel do collapse(3)
          !$acc parallel async
          !$acc loop collapse(3)
          do iq = QS_MP, QE_MP
          do j = JS, JE
          do i = IS, IE
          do k = KS, KE
             QTRC1(k,i,j,iq) = QTRC(k,i,j,iq)
          end do
          end do
          end do
          end do
          !$acc end parallel
!OCL XFILL
          !$omp parallel do collapse(2)
          !$acc parallel async
          !$acc loop collapse(2)
          do j = JS, JE
          do i = IS, IE
          do k = KS, KE
             CVtot1(k,i,j) = CVtot(k,i,j)
          end do
          end do
          end do
          !$acc end parallel
!OCL XFILL
          !$acc parallel async
          !$acc loop collapse(2)
          do j = JS, JE
          do i = IS, IE
          do k = KS, KE
             CPtot1(k,i,j) = CPtot(k,i,j)
          end do
          end do
          end do
          !$acc end parallel

          !$acc wait

          if( flg_lt ) then
!OCL XFILL
             !$omp parallel do collapse(3)
             !$acc parallel
             !$acc loop collapse(3)
             do iq = QS_LT, QE_LT
             do j = JS, JE
             do i = IS, IE
             do k = KS, KE
                QTRC1_crg(k,i,j,iq) = QTRC(k,i,j,iq)
             end do
             end do
             end do
             end do
             !$acc end parallel

             call ATMOS_PHY_LT_sato2019_select_dQCRG_from_LUT( &
                  KA, KS, KE, IA, IS, IE, JA, JS, JE, & ! [IN]
                  QLA,                                & ! [IN]
                  TEMP1(:,:,:), DENS(:,:,:),          & ! [IN]
                  QTRC1(:,:,:,QLS:QLE),               & ! [IN]
                  dqcrg(:,:,:), beta_crg(:,:,:)       ) ! [OUT]

             call ATMOS_PHY_MP_tomita08_adjustment( &
                  KA, KS, KE, IA, IS, IE, JA, JS, JE, &
                  DENS(:,:,:), PRES(:,:,:), CCN(:,:,:), dt_MP,                          & ! [IN]
                  TEMP1(:,:,:), QTRC1(:,:,:,QS_MP:QE_MP), CPtot1(:,:,:), CVtot1(:,:,:), & ! [INOUT]
                  RHOE_t(:,:,:), EVAPORATE(:,:,:),                                      & ! [OUT]
                  flg_lt, d0_crg, v0_crg, dqcrg(:,:,:), beta_crg(:,:,:),                & ! [IN:optional]
                  QSPLT_in(:,:,:,:), Sarea(:,:,:,:),                                    & ! [OUT:optional]
                  QTRC1_crg(:,:,:,QS_LT:QE_LT)                                          ) ! [INOUT:optional]
          else
             call ATMOS_PHY_MP_tomita08_adjustment( &
                  KA, KS, KE, IA, IS, IE, JA, JS, JE, &
                  DENS(:,:,:), PRES(:,:,:), CCN(:,:,:), dt_MP,                          & ! [IN]
                  TEMP1(:,:,:), QTRC1(:,:,:,QS_MP:QE_MP), CPtot1(:,:,:), CVtot1(:,:,:), & ! [INOUT]
                  RHOE_t(:,:,:), EVAPORATE(:,:,:)                                       ) ! [OUT]
          endif

          !$omp parallel do collapse(3)
          !$acc parallel async
          !$acc loop collapse(3)
          do iq = QS_MP, QE_MP
          do j = JS, JE
          do i = IS, IE
          do k = KS, KE
             RHOQ_t_MP(k,i,j,iq) = ( QTRC1(k,i,j,iq) - QTRC(k,i,j,iq) ) * DENS(k,i,j) / dt_MP
          enddo
          enddo
          enddo
          enddo
          !$acc end parallel

          !$omp parallel do collapse(2)
          !$acc parallel async
          !$acc loop collapse(2)
          do j = JS, JE
          do i = IS, IE
          do k = KS, KE
             CPtot_t(k,i,j) = ( CPtot1(k,i,j) - CPtot(k,i,j) ) / dt_MP
             CVtot_t(k,i,j) = ( CVtot1(k,i,j) - CVtot(k,i,j) ) / dt_MP
          end do
          end do
          end do
          !$acc end parallel

          if( flg_lt ) then
             !$omp parallel do collapse(3)
             !$acc parallel async
             !$acc loop collapse(3)
             do iq = QS_LT, QE_LT
             do j = JS, JE
             do i = IS, IE
             do k = KS, KE
                RHOC_t_MP(k,i,j,iq) = ( QTRC1_crg(k,i,j,iq) - QTRC(k,i,j,iq) ) * DENS(k,i,j) / dt_MP
             enddo
             enddo
             enddo
             enddo
             !$acc end parallel
          endif

          !$acc wait

       case ( 'SN14' )

          !$acc update host(DENS,W,QTRC(:,:,:,QS_MP:QE_MP),PRES,TEMP,Qdry,CPtot,CVtot,CCN)
          if( flg_lt ) then
             call ATMOS_PHY_LT_sato2019_select_dQCRG_from_LUT( &
                  KA, KS, KE, IA, IS, IE, JA, JS, JE, & ! [IN]
                  QLA,                                & ! [IN]
                  TEMP(:,:,:), DENS(:,:,:),           & ! [IN]
                  QTRC(:,:,:,QLS:QLE),                & ! [IN]
                  dqcrg(:,:,:), beta_crg(:,:,:)       ) ! [OUT]
             !$acc update host(dqcrg,beta_crg,QTRC(:,:,:,QS_LT:QE_LT))
             call ATMOS_PHY_MP_sn14_tendency( &
                  KA, KS, KE, IA, IS, IE, JA, JS, JE, &
                  DENS(:,:,:), W(:,:,:), QTRC(:,:,:,QS_MP:QE_MP), PRES(:,:,:), TEMP(:,:,:),                     & ! [IN]
                  Qdry(:,:,:), CPtot(:,:,:), CVtot(:,:,:), CCN(:,:,:), dt_MP, REAL_CZ(:,:,:), REAL_FZ(:,:,:),   & ! [IN]
                  RHOQ_t_MP(:,:,:,QS_MP:QE_MP), RHOE_t(:,:,:), CPtot_t(:,:,:), CVtot_t(:,:,:), EVAPORATE(:,:,:),& ! [OUT]
                  flg_lt, d0_crg, v0_crg, dqcrg(:,:,:), beta_crg(:,:,:),                                        & ! [IN:optional]
                  QTRC(:,:,:,QS_LT:QE_LT),                                                                      & ! [IN:optional]
                  QSPLT_in(:,:,:,:), Sarea(:,:,:,:),                                                            & ! [OUT:optional]
                  RHOC_t_MP(:,:,:,QS_LT:QE_LT)                                                                  ) ! [OUT:optional]
             !$acc update device(QSPLT_in,Sarea,RHOC_t_MP(:,:,:,QS_LT:QE_LT))
          else
             call ATMOS_PHY_MP_sn14_tendency( &
                  KA, KS, KE, IA, IS, IE, JA, JS, JE, &
                  DENS(:,:,:), W(:,:,:), QTRC(:,:,:,QS_MP:QE_MP), PRES(:,:,:), TEMP(:,:,:),                     & ! [IN]
                  Qdry(:,:,:), CPtot(:,:,:), CVtot(:,:,:), CCN(:,:,:), dt_MP, REAL_CZ(:,:,:), REAL_FZ(:,:,:),   & ! [IN]
                  RHOQ_t_MP(:,:,:,QS_MP:QE_MP), RHOE_t(:,:,:), CPtot_t(:,:,:), CVtot_t(:,:,:), EVAPORATE(:,:,:) ) ! [OUT]
          endif
          !$acc update device(RHOQ_t_MP(:,:,:,QS_MP:QE_MP),RHOE_t,CPtot_t,CVtot_t,EVAPORATE)

       case ( 'SUZUKI10' )

          !$acc update host(DENS,PRES,TEMP,QTRC(:,:,:,QS_MP:QE_MP),QDRY,CPtot,CVtot,CCN)
          if( flg_lt ) then
             call ATMOS_PHY_LT_sato2019_select_dQCRG_from_LUT( &
                  KA, KS, KE, IA, IS, IE, JA, JS, JE, & ! [IN]
                  QLA,                                & ! [IN]
                  TEMP(:,:,:), DENS(:,:,:),           & ! [IN]
                  QTRC(:,:,:,QLS:QLE),                & ! [IN]
                  dqcrg(:,:,:), beta_crg(:,:,:)       ) ! [OUT]
             !$acc update host(dqcrg,beta_crg,QTRC(:,:,:,QS_LT:QE_LT))
             call ATMOS_PHY_MP_suzuki10_tendency( KA, KS,  KE, IA, IS, IE, JA, JS, JE, KIJMAX, &
                                                  dt_MP,                                  & ! [IN]
                                                  DENS(:,:,:),  PRES(:,:,:), TEMP(:,:,:), & ! [IN]
                                                  QTRC(:,:,:,QS_MP:QE_MP), QDRY(:,:,:),   & ! [IN]
                                                  CPtot(:,:,:), CVtot(:,:,:),             & ! [IN]
                                                  CCN(:,:,:),                             & ! [IN]
                                                  RHOQ_t_MP(:,:,:,QS_MP:QE_MP),           & ! [OUT]
                                                  RHOE_t(:,:,:),                          & ! [OUT]
                                                  CPtot_t(:,:,:), CVtot_t(:,:,:),         & ! [OUT]
                                                  EVAPORATE(:,:,:),                       & ! [OUT]
                                                  flg_lt, d0_crg, v0_crg,                 & ! [IN:optional]
                                                  dqcrg(:,:,:), beta_crg(:,:,:),          & ! [IN:optional]
                                                  QTRC(:,:,:,QS_LT:QE_LT),                & ! [IN:optional]
                                                  QSPLT_in(:,:,:,:), Sarea(:,:,:,:),      & ! [OUT:optional]
                                                  RHOC_t_MP(:,:,:,QS_LT:QE_LT)            ) ! [OUT:optional]
             !$acc update device(QSPLT_in,Sarea,RHOC_t_MP(:,:,:,QS_LT:QE_LT))
          else
             call ATMOS_PHY_MP_suzuki10_tendency( KA, KS,  KE, IA, IS, IE, JA, JS, JE, KIJMAX, &
                                                  dt_MP,                                  & ! [IN]
                                                  DENS(:,:,:),  PRES(:,:,:), TEMP(:,:,:), & ! [IN]
                                                  QTRC(:,:,:,QS_MP:QE_MP), QDRY(:,:,:),   & ! [IN]
                                                  CPtot(:,:,:), CVtot(:,:,:),             & ! [IN]
                                                  CCN(:,:,:),                             & ! [IN]
                                                  RHOQ_t_MP(:,:,:,QS_MP:QE_MP),           & ! [OUT]
                                                  RHOE_t(:,:,:),                          & ! [OUT]
                                                  CPtot_t(:,:,:), CVtot_t(:,:,:),         & ! [OUT]
                                                  EVAPORATE(:,:,:)                        ) ! [OUT]
          endif
          !$acc update device(RHOQ_t_MP(:,:,:,QS_MP:QE_mp),RHOE_t,CPtot_t,CVtot_t,EVAPORATE)

          call ATMOS_PHY_MP_suzuki10_qtrc2qhyd( KA, KS, KE, IA, IS, IE, JA, JS, JE, &  ! [IN]
                                                QTRC(:,:,:,QS_MP+1:QE_MP),          &  ! [IN]
                                                Qe(:,:,:,:)                         )  ! [OUT]
          !$acc update device(Qe)

          do iq = 1, N_HYD
             call FILE_HISTORY_query( HIST_hyd_id(iq), HIST_sw(iq) )
          enddo
          do iq = 1, N_HYD
             if( HIST_sw(iq) ) then
                call FILE_HISTORY_put( HIST_hyd_id(iq), Qe(:,:,:,iq) )
             endif
          enddo

          if( flg_lt ) then

             call ATMOS_PHY_MP_suzuki10_crg_qtrc2qhyd( KA, KS, KE, IA, IS, IE, JA, JS, JE, &  ! [IN]
                                                       QTRC(:,:,:,QS_LT:QE_LT),            &  ! [IN]
                                                       Qecrg(:,:,:,:)                      )  ! [OUT]
             !$acc update device(Qecrg)

             do iq = 1, N_HYD
                call FILE_HISTORY_query( HIST_crg_id(iq), HIST_crg_sw(iq) )
             enddo
             do iq = 1, N_HYD
                if( HIST_crg_sw(iq) ) then
                   call FILE_HISTORY_put( HIST_crg_id(iq), Qecrg(:,:,:,iq) )
                endif
             enddo

          endif

       end select


       !$omp parallel do collapse(2)
       !$acc parallel
       !$acc loop collapse(2)
       do j = JS, JE
       do i = IS, IE
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
       !$acc end parallel

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
          !$acc update device(hist_vterm_idx)
          if ( ih > 0 ) then
             allocate( vterm_hist(KA,IA,JA,ih) )
             !$acc kernels
             vterm_hist(:,:,:,:) = 0.0_RP
             !$acc end kernels
          end if

          !$omp parallel do default(none) OMP_SCHEDULE_ collapse(2) &
          !$omp shared (KA,KS,KE,IS,IE,JS,JE,QS_MP,QE_MP,QHA,QHS,QHE,QLA,QIA, &
          !$omp         QA_LT,QS_LT,QE_LT, &
          !$omp         PRE00,LHF, &
          !$omp         ATMOS_PHY_MP_TYPE, ATMOS_PHY_PRECIP_TYPE, &
          !$omp         dt_MP,MP_NSTEP_SEDIMENTATION,MP_RNSTEP_SEDIMENTATION, &
          !$omp         REAL_CZ,REAL_FZ, &
          !$omp         DENS,MOMZ,U,V,TEMP,PRES,QTRC,CPtot,CVtot, &
          !$omp         DENS_t_MP,MOMZ_t_MP,RHOU_t_MP,RHOV_t_MP,RHOQ_t_MP,RHOH_MP, &
          !$omp         SFLX_rain,SFLX_snow,SFLX_ENGI, &
          !$omp         REFSTATE_dens, &
          !$omp         flg_lt,RHOC_t_MP, &
          !$omp         MP_DTSEC_SEDIMENTATION, &
          !$omp         vterm_hist,hist_vterm_idx) &
          !$omp private(i,j,k,iq,step, &
          !$omp         FZ,FDZ,RFDZ,RCDZ, &
          !$omp         DENS2,TEMP2,PRES2,CPtot2,CVtot2,RHOE,RHOE2,RHOQ,RHOQ2, &
          !$omp         RHOQ2_crg,mflux_crg,sflux_crg,eflux_crg, &
          !$omp         vterm,mflux,sflux,eflux,FLX_hydro,CP_t,CV_t)
          !$acc parallel
          !$acc loop collapse(2) gang &
          !$acc private(vterm,FLX_hydro,DENS2,TEMP2,PRES2,CPtot2,CVtot2,RHOE,RHOE2,RHOQ,RHOQ2,mflux,sflux,eflux,FZ,FDZ,RFDZ,RCDZ, &
          !$acc         RHOQ2_crg,mflux_crg,sflux_crg,eflux_crg)
          do j = JS, JE
          do i = IS, IE

             FZ(1:KA) = REAL_FZ(1:KA,i,j)

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
             !$acc loop collapse(2)
             do iq = QS_MP+1, QE_MP
             do k = KS, KE
                RHOQ (k,iq) = DENS2(k) * QTRC(k,i,j,iq) + RHOQ_t_MP(k,i,j,iq) * dt_MP
                RHOQ2(k,iq) = RHOQ(k,iq)
             end do
             end do

             if( flg_lt ) then
                !$acc loop collapse(2)
                do iq = QS_LT, QE_LT
                do k = KS, KE
                   RHOQ2_crg(k,iq) = DENS2(k) * QTRC(k,i,j,iq)
                end do
                end do
             endif

             SFLX_rain(i,j) = 0.0_RP
             SFLX_snow(i,j) = 0.0_RP
             SFLX_ENGI(i,j) = 0.0_RP
             FLX_hydro(:) = 0.0_RP

             !$acc loop seq
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
                   !$acc loop seq
                   do iq = QS_MP+1, QE_MP
                      do k = KS, KE
                         vterm(k,iq) = 0.0_RP ! tentative
                      end do
                   end do
                end select

                ! store to history output
                !$acc loop seq
                do iq = QS_MP+1, QE_MP
                   if ( hist_vterm_idx(iq) > 0 ) then
                      !$acc loop independent
                      do k = KS, KE
                         vterm_hist(k,i,j,hist_vterm_idx(iq)) = vterm_hist(k,i,j,hist_vterm_idx(iq)) &
                                                              + vterm(k,iq) * MP_RNSTEP_SEDIMENTATION
                      end do
                   end if
                end do

                select case ( ATMOS_PHY_PRECIP_TYPE )
                case ( 'Upwind-Euler' )
                   call ATMOS_PHY_MP_precipitation_upwind( &
                        KA, KS, KE, QE_MP-QS_MP, QLA, QIA, &
                        TEMP2(:), vterm(:,:),   & ! [IN]
                        FDZ(:), RCDZ(:),        & ! [IN]
                        MP_DTSEC_SEDIMENTATION, & ! [IN]
                        i, j,                   & ! [IN]
                        DENS2(:), RHOQ2(:,:),   & ! [INOUT]
                        CPtot2(:), CVtot2(:),   & ! [INOUT]
                        RHOE2(:),               & ! [INOUT]
                        mflux(:), sflux(:),     & ! [OUT]
                        eflux                   ) ! [OUT]
                case ( 'Semilag' )
                   call ATMOS_PHY_MP_precipitation_semilag( &
                        KA, KS, KE, QE_MP-QS_MP, QLA, QIA, &
                        TEMP2(:), vterm(:,:),   & ! [IN]
                        FZ(:), FDZ(:), RCDZ(:), & ! [IN]
                        MP_DTSEC_SEDIMENTATION, & ! [IN]
                        i, j,                   & ! [IN]
                        DENS2(:), RHOQ2(:,:),   & ! [INOUT]
                        CPtot2(:), CVtot2(:),   & ! [INOUT]
                        RHOE2(:),               & ! [INOUT]
                        mflux(:), sflux(:),     & ! [OUT]
                        eflux                   ) ! [OUT]
                case default
                   mflux(:) = 0.0_RP
                   sflux(:) = 0.0_RP
                   eflux    = 0.0_RP
                end select

                if( flg_lt ) then

                   select case ( ATMOS_PHY_PRECIP_TYPE )
                   case ( 'Upwind-Euler' )
                      call ATMOS_PHY_MP_precipitation_upwind( &
                           KA, KS, KE, QA_LT, 0, 0,    & ! no mass tracer for charge density
                           TEMP2(:), vterm(:,QHS:QHE), & ! [IN]
                           FDZ(:), RCDZ(:),            & ! [IN]
                           MP_DTSEC_SEDIMENTATION,     & ! [IN]
                           i, j,                       & ! [IN]
                           DENS2(:), RHOQ2_crg(:,:),   & ! [INOUT]
                           CPtot2(:), CVtot2(:),       & ! [INOUT]
                           RHOE2(:),                   & ! [INOUT]
                           mflux_crg(:), sflux_crg(:), & ! [OUT] dummy
                           eflux_crg                   ) ! [OUT] dummy
                   case ( 'Semilag' )
                      call ATMOS_PHY_MP_precipitation_semilag( &
                           KA, KS, KE, QA_LT, 0, 0,    & ! no mass tracer for charge density
                           TEMP2(:), vterm(:,QHS:QHE), & ! [IN]
                           FZ(:), FDZ(:), RCDZ(:),     & ! [IN]
                           MP_DTSEC_SEDIMENTATION,     & ! [IN]
                           i, j,                       & ! [IN]
                           DENS2(:), RHOQ2_crg(:,:),   & ! [INOUT]
                           CPtot2(:), CVtot2(:),       & ! [INOUT]
                           RHOE2(:),                   & ! [INOUT]
                           mflux_crg(:), sflux_crg(:), & ! [OUT] dummy
                           eflux_crg                   ) ! [OUT] dummy
                   case default
                      mflux_crg(:) = 0.0_RP
                      sflux_crg(:) = 0.0_RP
                      eflux_crg    = 0.0_RP
                   end select

                endif

                !$acc loop independent
                do k = KS, KE
                   TEMP2(k) = RHOE2(k) / ( DENS2(k) * CVtot2(k) )
                end do

                !$acc loop independent
                do k = KS-1, KE-1
                   FLX_hydro(k) = FLX_hydro(k) + mflux(k) * MP_RNSTEP_SEDIMENTATION
                enddo

                SFLX_rain(i,j) = SFLX_rain(i,j) - sflux(1) * MP_RNSTEP_SEDIMENTATION
                SFLX_snow(i,j) = SFLX_snow(i,j) - sflux(2) * MP_RNSTEP_SEDIMENTATION
                SFLX_ENGI(i,j) = SFLX_ENGI(i,j) - eflux    * MP_RNSTEP_SEDIMENTATION

             enddo

             SFLX_ENGI(i,j) = SFLX_ENGI(i,j) - SFLX_snow(i,j) * LHF ! moist internal energy

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

             !$acc loop collapse(2)
             do iq = QS_MP+1, QE_MP
             do k  = KS, KE
                RHOQ_t_MP(k,i,j,iq) = RHOQ_t_MP(k,i,j,iq) &
                     + ( RHOQ2(k,iq) - RHOQ(k,iq) ) / dt_MP
             enddo
             enddo

             if( flg_lt ) then
                !$acc loop collapse(2)
                do iq = QS_LT, QE_LT
                do k  = KS, KE
                   RHOC_t_MP(k,i,j,iq) = RHOC_t_MP(k,i,j,iq) &
                        + ( RHOQ2_crg(k,iq) - DENS(k,i,j) * QTRC(k,i,j,iq) ) / dt_MP
                enddo
                enddo
             endif

             call ATMOS_PHY_MP_precipitation_momentum( &
                     KA, KS, KE, &
                     DENS(:,i,j), MOMZ(:,i,j), U(:,i,j), V(:,i,j),        & ! [IN]
                     FLX_hydro(:),                                        & ! [IN]
                     RCDZ(:), RFDZ(:),                                    & ! [IN]
                     MOMZ_t_MP(:,i,j), RHOU_t_MP(:,i,j), RHOV_t_MP(:,i,j) ) ! [OUT]

          enddo
          enddo
          !$acc end parallel

          ! history output
          do iq = QS_MP+1, QE_MP
             if ( hist_vterm_idx(iq) > 0 ) then
                call FILE_HISTORY_put( hist_vterm_id(iq), vterm_hist(:,:,:,hist_vterm_idx(iq)) )
             end if
          end do
          if ( allocated( vterm_hist ) ) then
             deallocate( vterm_hist )
          end if

          call PROF_rapend  ('MP_Precipitation', 2)

       end if

!OCL XFILL
       !$acc kernels
       do j = JS, JE
       do i = IS, IE
          precip(i,j) = SFLX_rain(i,j) + SFLX_snow(i,j)
       end do
       end do
       !$acc end kernels

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

       if ( flg_lt ) then
          call FILE_HISTORY_in( QSPLT_in(:,:,:,1), 'QSPLT_G', 'Charge split of QG by Non-inductive process', 'fC/m3/s', fill_halo=.true. )
          call FILE_HISTORY_in( QSPLT_in(:,:,:,2), 'QSPLT_I', 'Charge split of QI by Non-inductive process', 'fC/m3/s', fill_halo=.true. )
          call FILE_HISTORY_in( QSPLT_in(:,:,:,3), 'QSPLT_S', 'Charge split of QS by Non-inductive process', 'fC/m3/s', fill_halo=.true. )

          do iq = QS_LT, QE_LT
             call FILE_HISTORY_in( RHOC_t_MP(:,:,:,iq), trim(TRACER_NAME(iq))//'_t_MP', &
                                  'tendency rho*'//trim(TRACER_NAME(iq))//' in MP', 'fC/m3/s', fill_halo=.true. )
          enddo

       end if

    endif

    !$omp parallel do default(none) private(i,j,k) OMP_SCHEDULE_ collapse(2) &
    !$omp shared(KS,KE,IS,IE,JS,JE, &
    !$omp        DENS_t_MP,MOMZ_t_MP,RHOU_t_MP,RHOV_t_MP,RHOH_MP, &
    !$omp        DENS_t,MOMZ_t,RHOU_t,RHOV_t,RHOH)
    !$acc parallel
    !$acc loop collapse(2)
    do j = JS, JE
    do i = IS, IE
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
    !$acc end parallel

    !$omp parallel do default(none) OMP_SCHEDULE_ collapse(3) &
    !$omp private(iq,i,j,k) &
    !$omp shared(QS_MP,QE_MP,JS,JE,IS,IE,KS,KE,RHOQ_t,RHOQ_t_MP)
    !$acc parallel
    !$acc loop collapse(3)
    do iq = QS_MP, QE_MP
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       RHOQ_t(k,i,j,iq) = RHOQ_t(k,i,j,iq) + RHOQ_t_MP(k,i,j,iq)
    enddo
    enddo
    enddo
    enddo
    !$acc end parallel

    if( flg_lt ) then
       !$omp parallel do collapse(3)
       !$acc parallel
       !$acc loop collapse(3)
       do iq = QS_LT, QE_LT
       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          RHOQ_t(k,i,j,iq) = RHOQ_t(k,i,j,iq) &
                           + RHOC_t_MP(k,i,j,iq)
       enddo
       enddo
       enddo
       enddo
       !$acc end parallel
    endif

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

    !$acc end data

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

    !$acc data copyin(QV, QHYD) copyout(QTRC)
    !$acc data copyin(QNUM) if(present(QNUM))

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
       !$acc update host(QHYD)
       call ATMOS_PHY_MP_KESSLER_qhyd2qtrc( KA, KS, KE, IA, IS, IE, JA, JS, JE, &
                                            QHYD(:,:,:,:),  & ! [IN]
                                            QTRC(:,:,:,2:)  ) ! [OUT]
       !$acc update device(QTRC(:,:,:,2:))
    case ( "TOMITA08" )
       !$omp parallel do OMP_SCHEDULE_
       !$acc kernels
       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          QTRC(k,i,j,1) = QV(k,i,j)
       end do
       end do
       end do
       !$acc end kernels
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
       !$acc update host(QHYD)
       call ATMOS_PHY_MP_SN14_qhyd2qtrc( KA, KS, KE, IA, IS, IE, JA, JS, JE, &
                                         QHYD(:,:,:,:),  & ! [IN]
                                         QTRC(:,:,:,2:), & ! [OUT]
                                         QNUM=QNUM       ) ! [IN]
       !$acc update device(QTRC(:,:,:,2:))
       !$acc update device(QNUM) if(present(QNUM))
    case ( "SUZUKI10" )
       !$omp parallel do OMP_SCHEDULE_
       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          QTRC(k,i,j,1) = QV(k,i,j)
       end do
       end do
       end do
       !$acc update host(QHYD)
       call ATMOS_PHY_MP_SUZUKI10_qhyd2qtrc( KA, KS, KE, IA, IS, IE, JA, JS, JE, &
                                             QHYD(:,:,:,:),  & ! [IN]
                                             QTRC(:,:,:,2:), & ! [OUT]
                                             QNUM=QNUM       ) ! [IN]
       !$acc update device(QTRC(:,:,:,2:))
       !$acc update device(QNUM) if(present(QNUM))
    case default
       LOG_ERROR("ATMOS_PHY_MP_driver_qhyd2qtrc",*) 'ATMOS_PHY_MP_TYPE (', trim(ATMOS_PHY_MP_TYPE), ') is not supported'
       call PRC_abort
    end select

    !$acc end data
    !$acc end data

    return
  end subroutine ATMOS_PHY_MP_driver_qhyd2qtrc

  subroutine ATMOS_PHY_MP_driver_qhyd2qtrc_onlyqv( &
       KA, KS, KE, IA, IS, IE, JA, JS, JE, &
       QV, QHYD, &
       QTRC,     &
       QNUM      )
    use scale_atmos_hydrometeor, only: &
       N_HYD
    use mod_atmos_phy_mp_vars, only: &
       QA_MP
    integer, intent(in) :: KA, KS, KE
    integer, intent(in) :: IA, IS, IE
    integer, intent(in) :: JA, JS, JE

    real(RP), intent(in) :: QV  (KA,IA,JA)
    real(RP), intent(in) :: QHYD(KA,IA,JA,N_HYD)

    real(RP), intent(out) :: QTRC(KA,IA,JA,QA_MP)

    real(RP), intent(in), optional :: QNUM(KA,IA,JA,N_HYD)

    integer :: k, i, j

    !$omp parallel do
    !$acc kernels copyin(QV) copyout(QTRC(:,:,:,1))
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       QTRC(k,i,j,1) = QV(k,i,j)
    end do
    end do
    end do
    !$acc end kernels

    return
  end subroutine ATMOS_PHY_MP_driver_qhyd2qtrc_onlyqv
end module mod_atmos_phy_mp_driver
