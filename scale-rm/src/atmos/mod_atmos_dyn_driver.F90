!-------------------------------------------------------------------------------
!> module Atmosphere / Dynamics
!!
!! @par Description
!!          Dynamical step driver
!!
!! @author Team SCALE
!!
!! @note The coding to call DYN2 routines is a temporary measure.
!!       After improving layering of directories and generalizing API
!!       for each modules in dynamical core, we should remove the temporary codes.
!!
!<
#include "scalelib.h"
module mod_atmos_dyn_driver
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_io
  use scale_prof
  use scale_atmos_grid_cartesC_index
  use scale_index
  use scale_tracer
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: ATMOS_DYN_driver_setup
  public :: ATMOS_DYN_driver

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  character(len=H_SHORT), public :: ATMOS_DYN_TSTEP_LARGE_TYPE     = 'FVM-HEVE'
  character(len=H_SHORT), public :: ATMOS_DYN_TSTEP_TRACER_TYPE    = 'FVM-HEVE'

  character(len=H_SHORT), public :: ATMOS_DYN_TINTEG_LARGE_TYPE    = 'EULER'        ! Type of time integration
                                                                   ! 'RK3'
  character(len=H_SHORT), public :: ATMOS_DYN_TINTEG_SHORT_TYPE    = 'RK4'
                                                                   ! 'RK3WS2002'
                                                                   ! 'RK3'
  character(len=H_SHORT), public :: ATMOS_DYN_TINTEG_TRACER_TYPE   = 'RK3WS2002'
                                                                   ! 'EULER'

  character(len=H_SHORT), public :: ATMOS_DYN_FVM_FLUX_TYPE        = 'CD4'          ! Type of advective flux scheme (FVM)
  character(len=H_SHORT), public :: ATMOS_DYN_FVM_FLUX_TRACER_TYPE = 'UD3KOREN1993'
                                                                   ! 'CD2'
                                                                   ! 'UD3'
                                                                   ! 'CD4'
                                                                   ! 'UD5'
                                                                   ! 'CD6'

  ! Coriolis force
  !> If ATMOS_DYN_coriolis_type=='PLANE', then f = ATMOS_DYN_coriolis_f0 + ATMOS_DYN_coriolis_beta * ( CY - ATMOS_DYN_coriolis_y0 )
  !> If ATMOS_DYN_coriolis_type=='SPHERE', then f = 2 * CONST_OHM * sin( lat )
  character(len=H_SHORT), public :: ATMOS_DYN_coriolis_type = 'PLANE'   ! type of coriolis force: 'PLANE', 'SPHERE'
  real(RP), public :: ATMOS_DYN_coriolis_f0                 = 0.0_RP
  real(RP), public :: ATMOS_DYN_coriolis_beta               = 0.0_RP
  real(RP), public :: ATMOS_DYN_coriolis_y0                             ! default is domain center


  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  ! Numerical filter
  integer,  private :: ATMOS_DYN_NUMERICAL_DIFF_order        = 1
  real(RP), private :: ATMOS_DYN_NUMERICAL_DIFF_coef         = 1.0E-4_RP ! nondimensional numerical diffusion
  real(RP), private :: ATMOS_DYN_NUMERICAL_DIFF_coef_TRACER  = 0.0_RP    ! nondimensional numerical diffusion for tracer
  real(RP), private :: ATMOS_DYN_NUMERICAL_DIFF_sfc_fact     = 1.0_RP
  logical , private :: ATMOS_DYN_NUMERICAL_DIFF_use_refstate = .true.

  real(RP), private :: ATMOS_DYN_wdamp_tau                   = -1.0_RP   ! maximum tau for Rayleigh damping of w [s]
  real(RP), private :: ATMOS_DYN_wdamp_height                = -1.0_RP   ! height       to start apply Rayleigh damping [m]
  integer,  private :: ATMOS_DYN_wdamp_layer                 = -1        ! layer number to start apply Rayleigh damping [num]

  ! Divergence damping
  real(RP), private :: ATMOS_DYN_divdmp_coef                 = 0.0_RP    ! Divergence dumping coef

  ! Flux-Corrected Transport limiter
  logical,  private :: ATMOS_DYN_FLAG_TRACER_SPLIT_TEND      = .false.
  logical,  private :: ATMOS_DYN_FLAG_FCT_momentum           = .false.
  logical,  private :: ATMOS_DYN_FLAG_FCT_T                  = .false.
  logical,  private :: ATMOS_DYN_FLAG_FCT_TRACER             = .false.
  logical,  private :: ATMOS_DYN_FLAG_FCT_along_stream       = .true.

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine ATMOS_DYN_driver_setup
    use scale_prc, only: &
       PRC_abort
    use scale_atmos_grid_cartesC, only: &
       DOMAIN_CENTER_Y => ATMOS_GRID_CARTESC_DOMAIN_CENTER_Y, &
       CY  => ATMOS_GRID_CARTESC_CY,  &
       FZ  => ATMOS_GRID_CARTESC_FZ,  &
       CDZ => ATMOS_GRID_CARTESC_CDZ, &
       CDX => ATMOS_GRID_CARTESC_CDX, &
       CDY => ATMOS_GRID_CARTESC_CDY, &
       FDZ => ATMOS_GRID_CARTESC_FDZ, &
       FDX => ATMOS_GRID_CARTESC_FDX, &
       FDY => ATMOS_GRID_CARTESC_FDY
    use scale_atmos_grid_cartesC_real, only: &
       REAL_LAT => ATMOS_GRID_CARTESC_REAL_LAT
    use scale_time, only: &
       TIME_DTSEC_ATMOS_DYN
    use mod_atmos_admin, only: &
       ATMOS_sw_dyn,   &
       ATMOS_DYN_TYPE
    use mod_atmos_vars, only: &
       DENS, &
       MOMZ, &
       MOMX, &
       MOMY, &
       RHOT, &
       QTRC
    use mod_atmos_dyn_vars, only: &
       PROG
    use scale_atmos_dyn, only: &
       ATMOS_DYN_setup
    implicit none

    namelist / PARAM_ATMOS_DYN / &
       ATMOS_DYN_TINTEG_SHORT_TYPE,           &
       ATMOS_DYN_TINTEG_TRACER_TYPE,          &
       ATMOS_DYN_TINTEG_LARGE_TYPE,           &
       ATMOS_DYN_FVM_FLUX_TYPE,               &
       ATMOS_DYN_FVM_FLUX_TRACER_TYPE,        &
       ATMOS_DYN_NUMERICAL_DIFF_order,        &
       ATMOS_DYN_NUMERICAL_DIFF_COEF,         &
       ATMOS_DYN_NUMERICAL_DIFF_COEF_TRACER,  &
       ATMOS_DYN_NUMERICAL_DIFF_sfc_fact,     &
       ATMOS_DYN_NUMERICAL_DIFF_use_refstate, &
       ATMOS_DYN_wdamp_tau,                   &
       ATMOS_DYN_wdamp_height,                &
       ATMOS_DYN_wdamp_layer,                 &
       ATMOS_DYN_coriolis_type,               &
       ATMOS_DYN_coriolis_f0,                 &
       ATMOS_DYN_coriolis_beta,               &
       ATMOS_DYN_coriolis_y0,                 &
       ATMOS_DYN_divdmp_coef,                 &
       ATMOS_DYN_FLAG_TRACER_SPLIT_TEND,      &
       ATMOS_DYN_FLAG_FCT_momentum,           &
       ATMOS_DYN_FLAG_FCT_T,                  &
       ATMOS_DYN_FLAG_FCT_TRACER,             &
       ATMOS_DYN_FLAG_FCT_along_stream

    real(RP) :: DT
    integer  :: ierr
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("ATMOS_DYN_driver_setup",*) 'Setup'

    if ( ATMOS_sw_dyn ) then

       ATMOS_DYN_coriolis_y0 = DOMAIN_CENTER_Y
       !--- read namelist
       rewind(IO_FID_CONF)
       read(IO_FID_CONF,nml=PARAM_ATMOS_DYN,iostat=ierr)
       if( ierr < 0 ) then !--- missing
          LOG_INFO("ATMOS_DYN_driver_setup",*) 'Not found namelist. Default used.'
       elseif( ierr > 0 ) then !--- fatal error
          LOG_ERROR("ATMOS_DYN_driver_setup",*) 'Not appropriate names in namelist PARAM_ATMOS_DYN. Check!'
          call PRC_abort
       endif
       LOG_NML(PARAM_ATMOS_DYN)

       DT = real(TIME_DTSEC_ATMOS_DYN,kind=RP)

       if ( ATMOS_DYN_wdamp_layer > KMAX ) then
          LOG_ERROR("ATMOS_DYN_driver_setup",*) 'ATMOS_DYN_wdamp_layer should be less than total number of vertical layer(KA). Check!'
          call PRC_abort
       elseif( ATMOS_DYN_wdamp_layer > 0 ) then
          ATMOS_DYN_wdamp_height = FZ(ATMOS_DYN_wdamp_layer+KS-1)
       endif

       if ( ATMOS_DYN_wdamp_tau < 0.0_RP ) then
          ATMOS_DYN_wdamp_tau = DT * 10.0_RP
       elseif ( ATMOS_DYN_wdamp_tau < DT ) then
          LOG_ERROR("ATMOS_DYN_driver_setup",*) 'ATMOS_DYN_wdamp_tau should be larger than TIME_DT_ATMOS_DYN. Check!'
          call PRC_abort
       end if

       if ( ATMOS_sw_dyn ) then
          LOG_NEWLINE
          LOG_INFO("ATMOS_DYN_driver_setup",*) 'Scheme for Large time step  : ', trim(ATMOS_DYN_TINTEG_LARGE_TYPE)
          LOG_INFO("ATMOS_DYN_driver_setup",*) 'Scheme for Short time step  : ', trim(ATMOS_DYN_TINTEG_SHORT_TYPE)
          LOG_INFO("ATMOS_DYN_driver_setup",*) 'Scheme for Tracer advection : ', trim(ATMOS_DYN_TINTEG_TRACER_TYPE)
       endif

       call ATMOS_DYN_setup( ATMOS_DYN_TINTEG_SHORT_TYPE,        & ! [IN]
                             ATMOS_DYN_TINTEG_TRACER_TYPE,       & ! [IN]
                             ATMOS_DYN_TINTEG_LARGE_TYPE,        & ! [IN]
                             ATMOS_DYN_TSTEP_TRACER_TYPE,        & ! [IN]
                             ATMOS_DYN_TSTEP_LARGE_TYPE,         & ! [IN]
                             ATMOS_DYN_FVM_FLUX_TYPE,            & ! [IN]
                             ATMOS_DYN_FVM_FLUX_TRACER_TYPE,     & ! [IN]
                             DENS, MOMZ, MOMX, MOMY, RHOT, QTRC, & ! [IN]
                             PROG,                               & ! [IN]
                             CDZ, CDX, CDY, FDZ, FDX, FDY,       & ! [IN]
                             ATMOS_DYN_wdamp_tau,                & ! [IN]
                             ATMOS_DYN_wdamp_height,             & ! [IN]
                             FZ,                            & ! [IN]
                             ATMOS_DYN_coriolis_type,            & ! [IN]
                             ATMOS_DYN_coriolis_f0,              & ! [IN]
                             ATMOS_DYN_coriolis_beta,            & ! [IN]
                             ATMOS_DYN_coriolis_y0,              & ! [IN]
                             CY,                            & ! [IN]
                             REAL_LAT,                           & ! [IN]
                             none = ATMOS_DYN_TYPE=='NONE'       ) ! [IN]
    endif

    return
  end subroutine ATMOS_DYN_driver_setup

  !-----------------------------------------------------------------------------
  !> Dynamical Process (Wrapper)
  subroutine ATMOS_DYN_driver( do_flag )
    use scale_atmos_grid_cartesC, only: &
       CDZ  => ATMOS_GRID_CARTESC_CDZ,  &
       CDX  => ATMOS_GRID_CARTESC_CDX,  &
       CDY  => ATMOS_GRID_CARTESC_CDY,  &
       FDZ  => ATMOS_GRID_CARTESC_FDZ,  &
       FDX  => ATMOS_GRID_CARTESC_FDX,  &
       FDY  => ATMOS_GRID_CARTESC_FDY,  &
       RCDZ => ATMOS_GRID_CARTESC_RCDZ, &
       RCDX => ATMOS_GRID_CARTESC_RCDX, &
       RCDY => ATMOS_GRID_CARTESC_RCDY, &
       RFDZ => ATMOS_GRID_CARTESC_RFDZ, &
       RFDX => ATMOS_GRID_CARTESC_RFDX, &
       RFDY => ATMOS_GRID_CARTESC_RFDY
    use scale_atmos_grid_cartesC_real, only: &
       REAL_PHI => ATMOS_GRID_CARTESC_REAL_PHI
    use scale_atmos_grid_cartesC_metric, only: &
       GSQRT => ATMOS_GRID_CARTESC_METRIC_GSQRT, &
       J13G  => ATMOS_GRID_CARTESC_METRIC_J13G,  &
       J23G  => ATMOS_GRID_CARTESC_METRIC_J23G,  &
       J33G  => ATMOS_GRID_CARTESC_METRIC_J33G,  &
       MAPF  => ATMOS_GRID_CARTESC_METRIC_MAPF
    use scale_time, only: &
       TIME_DTSEC,           &
       TIME_DTSEC_ATMOS_DYN
    use mod_atmos_admin, only: &
       ATMOS_USE_AVERAGE
    use mod_atmos_vars, only: &
       ATMOS_vars_total,  &
       I_QV,    &
       DENS,    &
       MOMZ,    &
       MOMX,    &
       MOMY,    &
       RHOT,    &
       QTRC,    &
       DENS_av, &
       MOMZ_av, &
       MOMX_av, &
       MOMY_av, &
       RHOT_av, &
       QTRC_av, &
       DENS_tp, &
       RHOU_tp, &
       RHOV_tp, &
       RHOT_tp, &
       RHOQ_tp, &
       RHOH_p,  &
       MOMZ_tp, &
       MOMX_tp, &
       MOMY_tp, &
       CPtot,   &
       EXNER
    use mod_atmos_dyn_vars, only: &
       PROG
    use scale_atmos_refstate, only: &
       ATMOS_REFSTATE_dens, &
       ATMOS_REFSTATE_pott, &
       ATMOS_REFSTATE_qv,   &
       ATMOS_REFSTATE_pres
    use mod_atmos_bnd_driver, only: &
       ATMOS_BOUNDARY_DENS,       &
       ATMOS_BOUNDARY_VELZ,       &
       ATMOS_BOUNDARY_VELX,       &
       ATMOS_BOUNDARY_VELY,       &
       ATMOS_BOUNDARY_POTT,       &
       ATMOS_BOUNDARY_QTRC,       &
       ATMOS_BOUNDARY_alpha_DENS, &
       ATMOS_BOUNDARY_alpha_VELZ, &
       ATMOS_BOUNDARY_alpha_VELX, &
       ATMOS_BOUNDARY_alpha_VELY, &
       ATMOS_BOUNDARY_alpha_POTT, &
       ATMOS_BOUNDARY_alpha_QTRC, &
       BND_QA, &
       ATMOS_BOUNDARY_SMOOTHER_FACT
    use scale_atmos_dyn, only: &
       ATMOS_DYN
    use scale_comm_cartesC, only: &
       COMM_vars8, &
       COMM_wait
    implicit none

    logical, intent(in) :: do_flag

    integer :: k, i, j, iq
    !---------------------------------------------------------------------------

#if defined DEBUG || defined QUICKDEBUG
    RHOT_tp(   1:KS-1,:,:) = 0.0_RP
    RHOT_tp(KE+1:KA,  :,:) = 0.0_RP
    MOMX_tp(   1:KS-1,:,:) = 0.0_RP
    MOMX_tp(KE+1:KA,  :,:) = 0.0_RP
    MOMY_tp(   1:KS-1,:,:) = 0.0_RP
    MOMY_tp(KE+1:KA,  :,:) = 0.0_RP
#endif

    if ( do_flag ) then

       call COMM_vars8( RHOU_tp, 1 )
       call COMM_vars8( RHOV_tp, 2 )
       call COMM_vars8( MOMZ_tp, 3 )
       do iq = 1, QA
          call COMM_vars8( RHOQ_tp(:,:,:,iq), 4+iq)
       end do

       !$omp parallel do default(none) OMP_SCHEDULE_ collapse(2) &
       !$omp private(k,i,j) &
       !$omp shared (KA,KS,KE,ISB,IEB,JSB,JEB, &
       !$omp         RHOT_tp,RHOH_p,CPtot,EXNER)
       do j = JSB, JEB
       do i = ISB, IEB
       do k = KS, KE
          RHOT_tp(k,i,j) = RHOT_tp(k,i,j) &
                         + RHOH_p (k,i,j) / ( CPtot(k,i,j) * EXNER(k,i,j) )
       end do
       end do
       end do
       call COMM_vars8( RHOT_tp, 4 )

       call COMM_wait ( RHOU_tp, 1 )
       !$omp parallel do default(none) OMP_SCHEDULE_ collapse(2) &
       !$omp private(k,i,j) &
       !$omp shared (KA,KS,KE,IS,IE,JS,JE, &
       !$omp         MOMX_tp,RHOU_tp)
       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          MOMX_tp(k,i,j) = MOMX_tp(k,i,j) &
                         + 0.5_RP * ( RHOU_tp(k,i,j) + RHOU_tp(k,i+1,j) )
       end do
       end do
       end do
       call COMM_vars8( MOMX_tp, 1 )

       call COMM_wait ( RHOV_tp, 2 )
       !$omp parallel do default(none) OMP_SCHEDULE_ collapse(2) &
       !$omp private(k,i,j) &
       !$omp shared (KA,KS,KE,IS,IE,JS,JE, &
       !$omp         MOMY_tp,RHOV_tp)
       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          MOMY_tp(k,i,j) = MOMY_tp(k,i,j) &
                         + 0.5_RP * ( RHOV_tp(k,i,j) + RHOV_tp(k,i,j+1) )
       end do
       end do
       end do
       call COMM_vars8( MOMY_tp, 2 )

       call COMM_wait ( MOMZ_tp, 3 )
       call COMM_wait ( RHOT_tp, 4, .false. )
       do iq = 1, QA
          call COMM_wait ( RHOQ_tp(:,:,:,iq), 4+iq, .false. )
       end do
       call COMM_wait ( MOMX_tp, 1 )
       call COMM_wait ( MOMY_tp, 2 )


       call ATMOS_DYN( DENS, MOMZ, MOMX, MOMY, RHOT, QTRC,                   & ! [INOUT]
                       PROG,                                                 & ! [IN]
                       DENS_av, MOMZ_av, MOMX_av, MOMY_av, RHOT_av, QTRC_av, & ! [INOUT]
                       DENS_tp, MOMZ_tp, MOMX_tp, MOMY_tp, RHOT_tp, RHOQ_tp, & ! [IN]
                       CDZ,  CDX,  CDY,  FDZ,  FDX,  FDY,                    & ! [IN]
                       RCDZ, RCDX, RCDY, RFDZ, RFDX, RFDY,                   & ! [IN]
                       REAL_PHI,                                             & ! [IN]
                       GSQRT, J13G, J23G, J33G, MAPF,                        & ! [IN]
                       TRACER_R, TRACER_CV, TRACER_CP, TRACER_MASS,          & ! [IN]
                       ATMOS_REFSTATE_dens,                                  & ! [IN]
                       ATMOS_REFSTATE_pott,                                  & ! [IN]
                       ATMOS_REFSTATE_qv,                                    & ! [IN]
                       ATMOS_REFSTATE_pres,                                  & ! [IN]
                       ATMOS_DYN_NUMERICAL_DIFF_coef,                        & ! [IN]
                       ATMOS_DYN_NUMERICAL_DIFF_COEF_TRACER,                 & ! [IN]
                       ATMOS_DYN_NUMERICAL_DIFF_order,                       & ! [IN]
                       ATMOS_DYN_NUMERICAL_DIFF_sfc_fact,                    & ! [IN]
                       ATMOS_DYN_NUMERICAL_DIFF_use_refstate,                & ! [IN]
                       BND_QA, ATMOS_BOUNDARY_SMOOTHER_FACT,                 & ! [IN]
                       ATMOS_BOUNDARY_DENS,                                  & ! [IN]
                       ATMOS_BOUNDARY_VELZ,                                  & ! [IN]
                       ATMOS_BOUNDARY_VELX,                                  & ! [IN]
                       ATMOS_BOUNDARY_VELY,                                  & ! [IN]
                       ATMOS_BOUNDARY_POTT,                                  & ! [IN]
                       ATMOS_BOUNDARY_QTRC,                                  & ! [IN]
                       ATMOS_BOUNDARY_alpha_DENS,                            & ! [IN]
                       ATMOS_BOUNDARY_alpha_VELZ,                            & ! [IN]
                       ATMOS_BOUNDARY_alpha_VELX,                            & ! [IN]
                       ATMOS_BOUNDARY_alpha_VELY,                            & ! [IN]
                       ATMOS_BOUNDARY_alpha_POTT,                            & ! [IN]
                       ATMOS_BOUNDARY_alpha_QTRC,                            & ! [IN]
                       ATMOS_DYN_divdmp_coef,                                & ! [IN]
                       ATMOS_DYN_FLAG_TRACER_SPLIT_TEND,                     & ! [IN]
                       ATMOS_DYN_FLAG_FCT_momentum,                          & ! [IN]
                       ATMOS_DYN_FLAG_FCT_T,                                 & ! [IN]
                       ATMOS_DYN_FLAG_FCT_TRACER,                            & ! [IN]
                       ATMOS_DYN_FLAG_FCT_along_stream,                      & ! [IN]
                       ATMOS_USE_AVERAGE,                                    & ! [IN]
                       I_QV,                                                 & ! [IN]
                       TIME_DTSEC,                                           & ! [IN]
                       TIME_DTSEC_ATMOS_DYN                                  ) ! [IN]

       call ATMOS_vars_total
    endif

    return
  end subroutine ATMOS_DYN_driver

end module mod_atmos_dyn_driver
