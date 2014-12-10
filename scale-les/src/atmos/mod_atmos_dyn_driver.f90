!-------------------------------------------------------------------------------
!> module Atmosphere / Dynamics
!!
!! @par Description
!!          Dynamical step driver
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2013-12-04 (S.Nishizawa)  [mod] splited from scale_atmos_dyn.f90
!<
!-------------------------------------------------------------------------------
module mod_atmos_dyn_driver
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_stdio
  use scale_prof
  use scale_grid_index
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
  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  ! time integration scheme
  integer,  private, parameter :: RK = 3                                 ! order of Runge-Kutta scheme

  ! numerical filter
  integer,  private :: ATMOS_DYN_numerical_diff_order        = 1
  real(RP), private :: ATMOS_DYN_numerical_diff_coef         = 1.0E-4_RP ! nondimensional numerical diffusion
  real(RP), private :: ATMOS_DYN_numerical_diff_coef_q       = 1.0E-4_RP ! nondimensional numerical diffusion for tracer
  real(RP), private :: ATMOS_DYN_numerical_diff_sfc_fact     = 1.0_RP
  logical , private :: ATMOS_DYN_numerical_diff_use_refstate = .true.

  ! Coriolis force
  logical,  private :: ATMOS_DYN_enable_coriolis = .false.               ! enable coriolis force?

  ! divergence damping
  real(RP), private :: ATMOS_DYN_divdmp_coef = 0.0_RP                    ! Divergence dumping coef

  ! fct
  logical,  private :: ATMOS_DYN_FLAG_FCT_rho      = .false.
  logical,  private :: ATMOS_DYN_FLAG_FCT_momentum = .false.
  logical,  private :: ATMOS_DYN_FLAG_FCT_T        = .false.
  logical,  private :: ATMOS_DYN_FLAG_FCT_along_stream = .true.

  ! lateral boundary flux adjustment
  integer,  private :: ATMOS_DYN_adjust_flux_cell  = 1                   ! number of cells adjusting lateral boundary flux

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine ATMOS_DYN_driver_setup
    use scale_process, only: &
       PRC_MPIstop
    use scale_grid, only: &
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
    use scale_atmos_dyn, only: &
       ATMOS_DYN_setup
    use mod_atmos_admin, only: &
       ATMOS_DYN_TYPE, &
       ATMOS_sw_dyn
    use mod_atmos_vars, only: &
       DENS, &
       MOMZ, &
       MOMX, &
       MOMY, &
       RHOT, &
       QTRC
    implicit none

    NAMELIST / PARAM_ATMOS_DYN / &
       ATMOS_DYN_numerical_diff_order,        &
       ATMOS_DYN_numerical_diff_coef,         &
       ATMOS_DYN_numerical_diff_coef_q,       &
       ATMOS_DYN_numerical_diff_sfc_fact,     &
       ATMOS_DYN_numerical_diff_use_refstate, &
       ATMOS_DYN_enable_coriolis,             &
       ATMOS_DYN_divdmp_coef,                 &
       ATMOS_DYN_FLAG_FCT_rho,                &
       ATMOS_DYN_FLAG_FCT_momentum,           &
       ATMOS_DYN_FLAG_FCT_T,                  &
       ATMOS_DYN_FLAG_FCT_along_stream,       &
       ATMOS_DYN_adjust_flux_cell

    real(RP) :: DT
    integer  :: ierr
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[Dynamics]/Categ[ATMOS]'

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
       if( IO_L ) write(IO_FID_LOG,nml=PARAM_ATMOS_DYN)

       DT = real(TIME_DTSEC_ATMOS_DYN,kind=RP)

       call ATMOS_DYN_setup( ATMOS_DYN_TYPE,                 & ! [IN]
                             DENS, MOMZ, MOMX, MOMY, RHOT, QTRC, & ! [IN]
                             GRID_CDZ, GRID_CDX, GRID_CDY,   & ! [IN]
                             GRID_FDZ, GRID_FDX, GRID_FDY,   & ! [IN]
                             ATMOS_DYN_enable_coriolis,      & ! [IN]
                             REAL_LAT                        ) ! [IN]
    endif

    return
  end subroutine ATMOS_DYN_driver_setup

  !-----------------------------------------------------------------------------
  !> Dynamical Process (Wrapper)
  subroutine ATMOS_DYN_driver( do_flag )
    use scale_time, only: &
       TIME_DTSEC,           &
       TIME_DTSEC_ATMOS_DYN, &
       TIME_NSTEP_ATMOS_DYN
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
    use scale_atmos_thermodyn, only: &
       AQ_CV
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
    use scale_atmos_dyn, only: &
       ATMOS_DYN
    implicit none
    logical, intent(in) :: do_flag
    !---------------------------------------------------------------------------

    if ( do_flag ) then

       call ATMOS_DYN( &
         DENS, MOMZ, MOMX, MOMY, RHOT, QTRC,                   & ! [INOUT]
         DENS_av, MOMZ_av, MOMX_av, MOMY_av, RHOT_av, QTRC_av, & ! [INOUT]
         DENS_tp, MOMZ_tp, MOMX_tp, MOMY_tp, RHOT_tp, RHOQ_tp, & ! [IN]
         GRID_CDZ,  GRID_CDX,  GRID_CDY,                       & ! [IN]
         GRID_FDZ,  GRID_FDX,  GRID_FDY,                       & ! [IN]
         GRID_RCDZ, GRID_RCDX, GRID_RCDY,                      & ! [IN]
         GRID_RFDZ, GRID_RFDX, GRID_RFDY,                      & ! [IN]
         REAL_PHI,                                             & ! [IN]
         GTRANS_GSQRT,                                         & ! [IN]
         GTRANS_J13G, GTRANS_J23G, GTRANS_J33G, GTRANS_MAPF,   & ! [IN]
         AQ_CV,                                                & ! [IN]
         ATMOS_REFSTATE_dens,                                  & ! [IN]
         ATMOS_REFSTATE_pott,                                  & ! [IN]
         ATMOS_REFSTATE_qv,                                    & ! [IN]
         ATMOS_REFSTATE_pres,                                  & ! [IN]
         ATMOS_DYN_numerical_diff_coef,                        & ! [IN]
         ATMOS_DYN_numerical_diff_coef_q,                      & ! [IN]
         ATMOS_DYN_numerical_diff_order,                       & ! [IN]
         ATMOS_DYN_numerical_diff_sfc_fact,                    & ! [IN]
         ATMOS_DYN_numerical_diff_use_refstate,                & ! [IN]
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
         ATMOS_DYN_FLAG_FCT_rho,                               & ! [IN]
         ATMOS_DYN_FLAG_FCT_momentum,                          & ! [IN]
         ATMOS_DYN_FLAG_FCT_T,                                 & ! [IN]
         ATMOS_DYN_FLAG_FCT_along_stream,                      & ! [IN]
         ATMOS_USE_AVERAGE,                                    & ! [IN]
         TIME_DTSEC,                                           & ! [IN]
         TIME_DTSEC_ATMOS_DYN,                                 & ! [IN]
         TIME_NSTEP_ATMOS_DYN                                  ) ! [IN]

       call ATMOS_vars_total

    endif

    return
  end subroutine ATMOS_DYN_driver

end module mod_atmos_dyn_driver
