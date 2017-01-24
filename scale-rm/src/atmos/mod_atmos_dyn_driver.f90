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
!! @par History
!! @li      2013-12-04 (S.Nishizawa)  [mod] splited from scale_atmos_dyn.f90
!<
module mod_atmos_dyn_driver
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_stdio
  use scale_prof
  use scale_grid_index
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

  ! Coriolis force
  logical,  private :: ATMOS_DYN_enable_coriolis             = .false.   ! enable coriolis force?

  ! Divergence damping
  real(RP), private :: ATMOS_DYN_divdmp_coef                 = 0.0_RP    ! Divergence dumping coef

  ! Flux-Corrected Transport limiter
  logical,  private :: ATMOS_DYN_FLAG_FCT_momentum           = .false.
  logical,  private :: ATMOS_DYN_FLAG_FCT_T                  = .false.
  logical,  private :: ATMOS_DYN_FLAG_FCT_TRACER             = .false.
  logical,  private :: ATMOS_DYN_FLAG_FCT_along_stream       = .true.

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine ATMOS_DYN_driver_setup
    use scale_process, only: &
       PRC_MPIstop
    use scale_grid, only: &
       GRID_FZ,  &
       GRID_CDZ, &
       GRID_CDX, &
       GRID_CDY, &
       GRID_FDZ, &
       GRID_FDX, &
       GRID_FDY
    use scale_grid_real, only: &
       REAL_LAT
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

    NAMELIST / PARAM_ATMOS_DYN / &
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
       ATMOS_DYN_enable_coriolis,             &
       ATMOS_DYN_divdmp_coef,                 &
       ATMOS_DYN_FLAG_FCT_momentum,           &
       ATMOS_DYN_FLAG_FCT_T,                  &
       ATMOS_DYN_FLAG_FCT_TRACER,             &
       ATMOS_DYN_FLAG_FCT_along_stream

    real(RP) :: DT
    integer  :: ierr
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '++++++ Module[DRIVER] / Categ[ATMOS DYN] / Origin[SCALE-RM]'

    if ( ATMOS_sw_dyn ) then

       !--- read namelist
       rewind(IO_FID_CONF)
       read(IO_FID_CONF,nml=PARAM_ATMOS_DYN,iostat=ierr)
       if( ierr < 0 ) then !--- missing
          if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
       elseif( ierr > 0 ) then !--- fatal error
          write(*,*) 'xxx Not appropriate names in namelist PARAM_ATMOS_DYN. Check!'
          call PRC_MPIstop
       endif
       if( IO_NML ) write(IO_FID_NML,nml=PARAM_ATMOS_DYN)

       DT = real(TIME_DTSEC_ATMOS_DYN,kind=RP)

       if ( ATMOS_DYN_wdamp_layer > KMAX ) then
          write(*,*) 'xxx ATMOS_DYN_wdamp_layer should be less than total number of vertical layer(KA). Check!'
          call PRC_MPIstop
       elseif( ATMOS_DYN_wdamp_layer > 0 ) then
          ATMOS_DYN_wdamp_height = GRID_FZ(ATMOS_DYN_wdamp_layer+KS-1)
       endif

       if ( ATMOS_DYN_wdamp_tau < 0.0_RP ) then
          ATMOS_DYN_wdamp_tau = DT * 10.0_RP
       elseif ( ATMOS_DYN_wdamp_tau < DT ) then
          write(*,*) 'xxx ATMOS_DYN_wdamp_tau should be larger than TIME_DT_ATMOS_DYN. Check!'
          call PRC_MPIstop
       end if

       if ( ATMOS_sw_dyn ) then
          if( IO_L ) write(IO_FID_LOG,*)
          if( IO_L ) write(IO_FID_LOG,*) '*** Scheme for Large time step  : ', trim(ATMOS_DYN_TINTEG_LARGE_TYPE)
          if( IO_L ) write(IO_FID_LOG,*) '*** Scheme for Short time step  : ', trim(ATMOS_DYN_TINTEG_SHORT_TYPE)
          if( IO_L ) write(IO_FID_LOG,*) '*** Scheme for Tracer advection : ', trim(ATMOS_DYN_TINTEG_TRACER_TYPE)
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
                             GRID_CDZ, GRID_CDX, GRID_CDY,       & ! [IN]
                             GRID_FDZ, GRID_FDX, GRID_FDY,       & ! [IN]
                             ATMOS_DYN_wdamp_tau,                & ! [IN]
                             ATMOS_DYN_wdamp_height,             & ! [IN]
                             GRID_FZ,                            & ! [IN]
                             ATMOS_DYN_enable_coriolis,          & ! [IN]
                             REAL_LAT,                           & ! [IN]
                             none = ATMOS_DYN_TYPE=='NONE'       ) ! [IN]
    endif

    return
  end subroutine ATMOS_DYN_driver_setup

  !-----------------------------------------------------------------------------
  !> Dynamical Process (Wrapper)
  subroutine ATMOS_DYN_driver( do_flag )
    use scale_grid, only: &
       GRID_CDZ,  &
       GRID_CDX,  &
       GRID_CDY,  &
       GRID_FDZ,  &
       GRID_FDX,  &
       GRID_FDY,  &
       GRID_RCDZ, &
       GRID_RCDX, &
       GRID_RCDY, &
       GRID_RFDZ, &
       GRID_RFDX, &
       GRID_RFDY
    use scale_grid_real, only: &
       REAL_PHI
    use scale_gridtrans, only: &
       GTRANS_GSQRT, &
       GTRANS_J13G,  &
       GTRANS_J23G,  &
       GTRANS_J33G,  &
       GTRANS_MAPF
    use scale_time, only: &
       TIME_DTSEC,           &
       TIME_DTSEC_ATMOS_DYN
    use mod_atmos_admin, only: &
       ATMOS_USE_AVERAGE
    use mod_atmos_vars, only: &
       ATMOS_vars_total,  &
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
       MOMZ_tp, &
       MOMX_tp, &
       MOMY_tp, &
       RHOT_tp, &
       RHOQ_tp
    use mod_atmos_dyn_vars, only: &
       PROG
    use scale_atmos_refstate, only: &
       ATMOS_REFSTATE_dens, &
       ATMOS_REFSTATE_pott, &
       ATMOS_REFSTATE_qv,   &
       ATMOS_REFSTATE_pres
    use scale_atmos_boundary, only: &
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
       ATMOS_BOUNDARY_alpha_QTRC
#ifndef DYN2
    use scale_atmos_dyn, only: &
       ATMOS_DYN
#else
    use scale_atmos_dyn2, only: &
       ATMOS_DYN
#endif
    implicit none

    logical, intent(in) :: do_flag
    !---------------------------------------------------------------------------

    if ( do_flag ) then
       call ATMOS_DYN( DENS, MOMZ, MOMX, MOMY, RHOT, QTRC,                   & ! [INOUT]
                       PROG,                                                 & ! [IN]
                       DENS_av, MOMZ_av, MOMX_av, MOMY_av, RHOT_av, QTRC_av, & ! [INOUT]
                       DENS_tp, MOMZ_tp, MOMX_tp, MOMY_tp, RHOT_tp, RHOQ_tp, & ! [IN]
                       GRID_CDZ,  GRID_CDX,  GRID_CDY,                       & ! [IN]
                       GRID_FDZ,  GRID_FDX,  GRID_FDY,                       & ! [IN]
                       GRID_RCDZ, GRID_RCDX, GRID_RCDY,                      & ! [IN]
                       GRID_RFDZ, GRID_RFDX, GRID_RFDY,                      & ! [IN]
                       REAL_PHI,                                             & ! [IN]
                       GTRANS_GSQRT,                                         & ! [IN]
                       GTRANS_J13G, GTRANS_J23G, GTRANS_J33G, GTRANS_MAPF,   & ! [IN]
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
                       ATMOS_DYN_FLAG_FCT_momentum,                          & ! [IN]
                       ATMOS_DYN_FLAG_FCT_T,                                 & ! [IN]
                       ATMOS_DYN_FLAG_FCT_TRACER,                            & ! [IN]
                       ATMOS_DYN_FLAG_FCT_along_stream,                      & ! [IN]
                       ATMOS_USE_AVERAGE,                                    & ! [IN]
                       TIME_DTSEC,                                           & ! [IN]
                       TIME_DTSEC_ATMOS_DYN                                  ) ! [IN]

       call ATMOS_vars_total
    endif

    return
  end subroutine ATMOS_DYN_driver

end module mod_atmos_dyn_driver
