!-------------------------------------------------------------------------------
!> module Atmosphere / Dynamics FENT + FCT
!!
!! @par Description
!!          Dynamical core for Atmospheric process
!!          Full explicit, no terrain + tracer FCT limiter
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2013-12-04 (S.Nishizawa) [mod] splited from mod_atmos_dyn.f90
!!
!<
!-------------------------------------------------------------------------------
#include "inc_openmp.h"
module mod_atmos_dyn_wrap
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use mod_stdio, only: &
     IO_FID_LOG,  &
     IO_L
  use mod_time, only: &
     TIME_rapstart, &
     TIME_rapend
#ifdef DEBUG
  use mod_debug, only: &
     CHECK
  use mod_const, only: &
     UNDEF  => CONST_UNDEF, &
     IUNDEF => CONST_UNDEF2
#endif
  use mod_atmos_vars, only: &
     ZDIR,   &
     XDIR,   &
     YDIR,   &
     I_DENS, &
     I_MOMZ, &
     I_MOMX, &
     I_MOMY, &
     I_RHOT
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: ATMOS_DYN_wrap_setup
  public :: ATMOS_DYN_wrap

  !-----------------------------------------------------------------------------
  !
  !++ included parameters
  !
  include 'inc_precision.h'
  include 'inc_index.h'
  include 'inc_tracer.h'

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
  ! time integration scheme
  integer,  private, parameter :: RK = 3             ! order of Runge-Kutta scheme

  ! numerical filter
  integer,  private, save      :: ATMOS_DYN_numerical_diff_order        = 1
  real(RP), private, save      :: ATMOS_DYN_numerical_diff_coef         = 1.0E-4_RP ! nondimensional numerical diffusion
  real(RP), private, save      :: ATMOS_DYN_numerical_diff_sfc_fact     = 1.0_RP
  logical , private, save      :: ATMOS_DYN_numerical_diff_use_refstate = .true.
  real(RP), private, save      :: ATMOS_DYN_DIFF4                                   ! for numerical filter
  real(RP), private, save      :: ATMOS_DYN_CNZ3(3,KA,2)
  real(RP), private, save      :: ATMOS_DYN_CNX3(3,IA,2)
  real(RP), private, save      :: ATMOS_DYN_CNY3(3,JA,2)
  real(RP), private, save      :: ATMOS_DYN_CNZ4(5,KA,2)
  real(RP), private, save      :: ATMOS_DYN_CNX4(5,IA,2)
  real(RP), private, save      :: ATMOS_DYN_CNY4(5,JA,2)

  ! Coriolis force
  logical,  private, save      :: ATMOS_DYN_enable_coriolis = .false. ! enable coriolis force?
  real(RP), private, save      :: ATMOS_DYN_CORIOLI(IA,JA)            ! coriolis term

  ! divergence damping
  real(RP), private, save      :: ATMOS_DYN_divdmp_coef = 0.0_RP      ! Divergence dumping coef

  ! fct
  logical,  private, save      :: ATMOS_DYN_FLAG_FCT_rho      = .false.
  logical,  private, save      :: ATMOS_DYN_FLAG_FCT_momentum = .false.
  logical,  private, save      :: ATMOS_DYN_FLAG_FCT_T        = .false.

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine ATMOS_DYN_wrap_setup
    use mod_stdio, only: &
       IO_FID_CONF
    use mod_process, only: &
       PRC_MPIstop
    use mod_time, only: &
       TIME_DTSEC_ATMOS_DYN
    use mod_grid, only: &
       GRID_CDZ, &
       GRID_CDX, &
       GRID_CDY, &
       GRID_FDZ, &
       GRID_FDX, &
       GRID_FDY
    use mod_grid_real, only: &
       REAL_LAT
    use mod_atmos_dyn, only: &
       ATMOS_DYN_setup
    implicit none

    NAMELIST / PARAM_ATMOS_DYN / &
       ATMOS_DYN_numerical_diff_order,        &
       ATMOS_DYN_numerical_diff_coef,         &
       ATMOS_DYN_numerical_diff_sfc_fact,     &
       ATMOS_DYN_numerical_diff_use_refstate, &
       ATMOS_DYN_enable_coriolis,             &
       ATMOS_DYN_divdmp_coef,                 &
       ATMOS_DYN_FLAG_FCT_rho,                &
       ATMOS_DYN_FLAG_FCT_momentum,           &
       ATMOS_DYN_FLAG_FCT_T

    real(RP) :: DT

    integer :: ierr
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[Dynamics]/Categ[ATMOS]'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_ATMOS_DYN,iostat=ierr)

    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_ATMOS_DYN. Check!'
       call PRC_MPIstop
    endif
    if( IO_L ) write(IO_FID_LOG,nml=PARAM_ATMOS_DYN)

    !-- Block size must be divisible
    if    ( mod(IMAX,IBLOCK) > 0 ) then
       if( IO_L ) write(IO_FID_LOG,*) 'xxx number of grid size IMAX must be divisible by IBLOCK! ', IMAX, IBLOCK
       call PRC_MPIstop
    elseif( mod(JMAX,JBLOCK) > 0 ) then
       if( IO_L ) write(IO_FID_LOG,*) 'xxx number of grid size JMAX must be divisible by JBLOCK! ', JMAX, JBLOCK
       call PRC_MPIstop
    endif

    DT = real(TIME_DTSEC_ATMOS_DYN,kind=RP)

    call ATMOS_DYN_init( ATMOS_DYN_DIFF4,                & ! [OUT]
                         ATMOS_DYN_CNZ3,                 & ! [OUT]
                         ATMOS_DYN_CNX3,                 & ! [OUT]
                         ATMOS_DYN_CNY3,                 & ! [OUT]
                         ATMOS_DYN_CNZ4,                 & ! [OUT]
                         ATMOS_DYN_CNX4,                 & ! [OUT]
                         ATMOS_DYN_CNY4,                 & ! [OUT]
                         ATMOS_DYN_CORIOLI,              & ! [OUT]
                         GRID_CDZ, GRID_CDX, GRID_CDY,   & ! [IN]
                         GRID_FDZ, GRID_FDX, GRID_FDY,   & ! [IN]
                         ATMOS_DYN_numerical_diff_order, & ! [IN]
                         ATMOS_DYN_numerical_diff_coef,  & ! [IN]
                         DT,                             & ! [IN]
                         ATMOS_DYN_enable_coriolis,      & ! [IN]
                         REAL_LAT                        ) ! [IN]

    return
  end subroutine ATMOS_DYN_setup

  !-----------------------------------------------------------------------------
  !> Dynamical Process (Wrapper)
  subroutine ATMOS_DYN_wrap
    use mod_time, only: &
       TIME_DTSEC,           &
       TIME_DTSEC_ATMOS_DYN, &
       TIME_NSTEP_ATMOS_DYN
    use mod_grid, only: &
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
    use mod_grid_real, only: &
       REAL_PHI
    use mod_gridtrans, only: &
       GTRANS_GSQRT, &
       GTRANS_J13G,  &
       GTRANS_J23G,  &
       GTRANS_J33G
    use mod_history, only: &
       HIST_in
    use mod_atmos_vars, only: &
       ATMOS_vars_total,  &
       ATMOS_USE_AVERAGE, &
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
       QTRC_tp
    use mod_atmos_thermodyn, only: &
       AQ_CV
    use mod_atmos_refstate, only: &
       ATMOS_REFSTATE_dens, &
       ATMOS_REFSTATE_pott, &
       ATMOS_REFSTATE_qv,   &
       ATMOS_REFSTATE_pres
    use mod_atmos_boundary, only: &
       ATMOS_BOUNDARY_var,   &
       ATMOS_BOUNDARY_alpha
    use mod_atmos_dyn, only: &
       ATMOS_DYN
    implicit none
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*) '*** Dynamics step'

    call ATMOS_DYN_main( DENS, MOMZ, MOMX, MOMY, RHOT, QTRC,                   & ! [INOUT]
                         DENS_av, MOMZ_av, MOMX_av, MOMY_av, RHOT_av, QTRC_av, & ! [INOUT]
                         DENS_tp, MOMZ_tp, MOMX_tp, MOMY_tp, RHOT_tp, QTRC_tp, & ! [IN]
                         ATMOS_DYN_CNZ3,                                       & ! [IN]
                         ATMOS_DYN_CNX3,                                       & ! [IN]
                         ATMOS_DYN_CNY3,                                       & ! [IN]
                         ATMOS_DYN_CNZ4,                                       & ! [IN]
                         ATMOS_DYN_CNX4,                                       & ! [IN]
                         ATMOS_DYN_CNY4,                                       & ! [IN]
                         GRID_CDZ,  GRID_CDX,  GRID_CDY,                       & ! [IN]
                         GRID_FDZ,  GRID_FDX,  GRID_FDY,                       & ! [IN]
                         GRID_RCDZ, GRID_RCDX, GRID_RCDY,                      & ! [IN]
                         GRID_RFDZ, GRID_RFDX, GRID_RFDY,                      & ! [IN]
                         REAL_PHI,                                             & ! [IN]
                         GTRANS_GSQRT,                                         & ! [IN]
                         GTRANS_J13G, GTRANS_J23G, GTRANS_J33G,                & ! [IN]
                         AQ_CV,                                                & ! [IN]
                         ATMOS_REFSTATE_dens,                                  & ! [IN]
                         ATMOS_REFSTATE_pott,                                  & ! [IN]
                         ATMOS_REFSTATE_qv,                                    & ! [IN]
                         ATMOS_REFSTATE_pres,                                  & ! [IN]
                         ATMOS_DYN_DIFF4,                                      & ! [IN]
                         ATMOS_DYN_numerical_diff_order,                       & ! [IN]
                         ATMOS_DYN_numerical_diff_sfc_fact,                    & ! [IN]
                         ATMOS_DYN_numerical_diff_use_refstate,                & ! [IN]
                         ATMOS_DYN_CORIOLI,                                    & ! [IN]
                         ATMOS_BOUNDARY_var,                                   & ! [IN]
                         ATMOS_BOUNDARY_alpha,                                 & ! [IN]
                         ATMOS_DYN_divdmp_coef,                                & ! [IN]
                         ATMOS_DYN_FLAG_FCT_rho,                               & ! [IN]
                         ATMOS_DYN_FLAG_FCT_momentum,                          & ! [IN]
                         ATMOS_DYN_FLAG_FCT_T,                                 & ! [IN]
                         ATMOS_USE_AVERAGE,                                    & ! [IN]
                         TIME_DTSEC,                                           & ! [IN]
                         TIME_DTSEC_ATMOS_DYN,                                 & ! [IN]
                         TIME_NSTEP_ATMOS_DYN                                  ) ! [IN]

    call ATMOS_vars_total

    return
  end subroutine ATMOS_DYN


end module mod_atmos_dyn
