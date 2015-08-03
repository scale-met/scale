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
!! @li      2011-11-11 (H.Yashiro)   [new] Imported from SCALE-LES ver.2
!! @li      2011-11-11 (H.Yashiro)   [mod] Merged with Y.Miyamoto's
!! @li      2011-12-11 (H.Yashiro)   [mod] Use reference state
!! @li      2011-12-26 (Y.Miyamoto)  [mod] Add numerical diffusion into mass flux calc
!! @li      2012-01-04 (H.Yashiro)   [mod] Nonblocking communication (Y.Ohno)
!! @li      2012-01-25 (H.Yashiro)   [fix] Bugfix (Y.Miyamoto)
!! @li      2011-01-25 (H.Yashiro)   [mod] sprit as "full" FCT (Y.Miyamoto)
!! @li      2012-02-14 (H.Yashiro)   [mod] Cache tiling
!! @li      2012-03-14 (H.Yashiro)   [mod] Bugfix (Y.Miyamoto)
!! @li      2012-03-23 (H.Yashiro)   [mod] Explicit index parameter inclusion
!! @li      2012-04-09 (H.Yashiro)   [mod] Integrate RDMA communication
!! @li      2012-06-10 (Y.Miyamoto)  [mod] large-scale divergence (from H.Yashiro's)
!! @li      2012-07-13 (H.Yashiro)   [mod] prevent conditional branching in FCT
!! @li      2012-07-27 (Y.Miyamoto)  [mod] divegence damping option
!! @li      2012-08-16 (S.Nishizawa) [mod] use FCT for momentum and temperature
!! @li      2012-09-21 (Y.Sato)      [mod] merge DYCOMS-II experimental set
!! @li      2013-03-26 (Y.Sato)      [mod] modify Large scale forcing and corioli forcing
!! @li      2013-04-04 (Y.Sato)      [mod] modify Large scale forcing
!! @li      2013-06-14 (S.Nishizawa) [mod] enable to change order of numerical diffusion
!! @li      2013-06-18 (S.Nishizawa) [mod] split part of RK to other files
!! @li      2013-06-20 (S.Nishizawa) [mod] split large scale sining to other file
!!
!<
!-------------------------------------------------------------------------------
#include "inc_openmp.h"
module scale_atmos_dyn
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
#ifdef CHECK_MASS
  use mpi
#endif
  use scale_precision
  use scale_stdio
  use scale_prof
  use scale_grid_index
  use scale_index
  use scale_tracer

#ifdef DEBUG
  use scale_debug, only: &
     CHECK
  use scale_const, only: &
     UNDEF  => CONST_UNDEF, &
     IUNDEF => CONST_UNDEF2
#endif
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: ATMOS_DYN_setup
  public :: ATMOS_DYN

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
  real(RP), private, allocatable :: DENS_RK1(:,:,:) ! prognostic variables (+1/3 step)
  real(RP), private, allocatable :: MOMZ_RK1(:,:,:)
  real(RP), private, allocatable :: MOMX_RK1(:,:,:)
  real(RP), private, allocatable :: MOMY_RK1(:,:,:)
  real(RP), private, allocatable :: RHOT_RK1(:,:,:)
  real(RP), private, allocatable :: DENS_RK2(:,:,:) ! prognostic variables (+2/3 step)
  real(RP), private, allocatable :: MOMZ_RK2(:,:,:)
  real(RP), private, allocatable :: MOMX_RK2(:,:,:)
  real(RP), private, allocatable :: MOMY_RK2(:,:,:)
  real(RP), private, allocatable :: RHOT_RK2(:,:,:)

  ! tendency
  real(RP), private, allocatable :: DENS_t(:,:,:)
  real(RP), private, allocatable :: MOMZ_t(:,:,:)
  real(RP), private, allocatable :: MOMX_t(:,:,:)
  real(RP), private, allocatable :: MOMY_t(:,:,:)
  real(RP), private, allocatable :: RHOT_t(:,:,:)
  real(RP), private, allocatable :: RHOQ_t(:,:,:,:)

  real(RP), private, allocatable :: CORIOLI(:,:)            ! coriolis term
  real(RP), private, allocatable :: mflx_hi(:,:,:,:)        ! rho * vel(x,y,z) @ (u,v,w)-face high order

  real(RP), private, allocatable :: num_diff  (:,:,:,:,:)
  real(RP), private, allocatable :: num_diff_q(:,:,:,:)

  logical, private :: BND_W
  logical, private :: BND_E
  logical, private :: BND_S
  logical, private :: BND_N

  ! for communication
  integer :: I_COMM_DENS = 1
  integer :: I_COMM_MOMZ = 2
  integer :: I_COMM_MOMX = 3
  integer :: I_COMM_MOMY = 4
  integer :: I_COMM_RHOT = 5

  integer :: I_COMM_DENS_t = 1
  integer :: I_COMM_MOMZ_t = 2
  integer :: I_COMM_MOMX_t = 3
  integer :: I_COMM_MOMY_t = 4
  integer :: I_COMM_RHOT_t = 5

  integer :: I_COMM_DENS_RK1 = 1
  integer :: I_COMM_MOMZ_RK1 = 2
  integer :: I_COMM_MOMX_RK1 = 3
  integer :: I_COMM_MOMY_RK1 = 4
  integer :: I_COMM_RHOT_RK1 = 5

  integer :: I_COMM_DENS_RK2 = 1
  integer :: I_COMM_MOMZ_RK2 = 2
  integer :: I_COMM_MOMX_RK2 = 3
  integer :: I_COMM_MOMY_RK2 = 4
  integer :: I_COMM_RHOT_RK2 = 5

  integer, allocatable :: I_COMM_RHOQ_t(:)
  integer, allocatable :: I_COMM_QTRC(:)

  integer :: I_COMM_mflx_z = 1
  integer :: I_COMM_mflx_x = 2
  integer :: I_COMM_mflx_y = 3

  !-----------------------------------------------------------------------------
contains

  !-----------------------------------------------------------------------------
  !> Setup
  subroutine ATMOS_DYN_setup( &
       DYN_TYPE,         &
       DENS, MOMZ, MOMX, MOMY, RHOT, QTRC, &
       CDZ, CDX, CDY,    &
       FDZ, FDX, FDY,    &
       enable_coriolis,  &
       lat               )
    use scale_process, only: &
       PRC_MPIstop, &
       PRC_HAS_E, &
       PRC_HAS_W, &
       PRC_HAS_N, &
       PRC_HAS_S
    use scale_const, only: &
       OHM => CONST_OHM, &
       UNDEF => CONST_UNDEF
    use scale_comm, only: &
       COMM_vars8_init
    use scale_atmos_dyn_common, only: &
       ATMOS_DYN_filter_setup
    use scale_atmos_dyn_rk, only: &
       ATMOS_DYN_rk_setup
    implicit none

    character(len=H_SHORT), intent(in) :: DYN_TYPE

    ! MPI_RECV_INIT requires intent(inout)
    real(RP),               intent(inout) :: DENS(KA,IA,JA)
    real(RP),               intent(inout) :: MOMZ(KA,IA,JA)
    real(RP),               intent(inout) :: MOMX(KA,IA,JA)
    real(RP),               intent(inout) :: MOMY(KA,IA,JA)
    real(RP),               intent(inout) :: RHOT(KA,IA,JA)
    real(RP),               intent(inout) :: QTRC(KA,IA,JA,QA)

    real(RP),               intent(in) :: CDZ(KA)
    real(RP),               intent(in) :: CDX(IA)
    real(RP),               intent(in) :: CDY(JA)
    real(RP),               intent(in) :: FDZ(KA-1)
    real(RP),               intent(in) :: FDX(IA-1)
    real(RP),               intent(in) :: FDY(JA-1)
    logical,                intent(in) :: enable_coriolis
    real(RP),               intent(in) :: lat(IA,JA)

    integer :: iq
    !---------------------------------------------------------------------------

    allocate( DENS_RK1(KA,IA,JA) )
    allocate( MOMZ_RK1(KA,IA,JA) )
    allocate( MOMX_RK1(KA,IA,JA) )
    allocate( MOMY_RK1(KA,IA,JA) )
    allocate( RHOT_RK1(KA,IA,JA) )
    allocate( DENS_RK2(KA,IA,JA) )
    allocate( MOMZ_RK2(KA,IA,JA) )
    allocate( MOMX_RK2(KA,IA,JA) )
    allocate( MOMY_RK2(KA,IA,JA) )
    allocate( RHOT_RK2(KA,IA,JA) )

    allocate( DENS_t(KA,IA,JA) )
    allocate( MOMZ_t(KA,IA,JA) )
    allocate( MOMX_t(KA,IA,JA) )
    allocate( MOMY_t(KA,IA,JA) )
    allocate( RHOT_t(KA,IA,JA) )
    allocate( RHOQ_t(KA,IA,JA,QA) )

    allocate( CORIOLI(IA,JA) )
    allocate( mflx_hi(KA,IA,JA,3) )

    allocate( num_diff  (KA,IA,JA,5,3) )
    allocate( num_diff_q(KA,IA,JA,3) )

    allocate( I_COMM_RHOQ_t(QA) )
    allocate( I_COMM_QTRC(QA) )

    ! numerical diffusion
    call ATMOS_DYN_filter_setup( &
         num_diff, num_diff_q, & ! (inout)
         CDZ, CDX, CDY, FDZ, FDX, FDY ) ! (in)

    ! coriolis parameter
    if ( enable_coriolis ) then
       CORIOLI(:,:) = 2.0_RP * OHM * sin( lat(:,:) )
    else
       CORIOLI(:,:) = 0.0_RP
    endif

    BND_W = .not. PRC_HAS_W
    BND_E = .not. PRC_HAS_E
    BND_S = .not. PRC_HAS_S
    BND_N = .not. PRC_HAS_N

    call ATMOS_DYN_rk_setup( &
         DYN_TYPE, &
         BND_W, BND_E, BND_S, BND_N )

    call COMM_vars8_init( DENS, I_COMM_DENS )
    call COMM_vars8_init( MOMZ, I_COMM_MOMZ )
    call COMM_vars8_init( MOMX, I_COMM_MOMX )
    call COMM_vars8_init( MOMY, I_COMM_MOMY )
    call COMM_vars8_init( RHOT, I_COMM_RHOT )

    call COMM_vars8_init( DENS_t, I_COMM_DENS_t )
    call COMM_vars8_init( MOMZ_t, I_COMM_MOMZ_t )
    call COMM_vars8_init( MOMX_t, I_COMM_MOMX_t )
    call COMM_vars8_init( MOMY_t, I_COMM_MOMY_t )
    call COMM_vars8_init( RHOT_t, I_COMM_RHOT_t )

    call COMM_vars8_init( DENS_RK1, I_COMM_DENS_RK1 )
    call COMM_vars8_init( MOMZ_RK1, I_COMM_MOMZ_RK1 )
    call COMM_vars8_init( MOMX_RK1, I_COMM_MOMX_RK1 )
    call COMM_vars8_init( MOMY_RK1, I_COMM_MOMY_RK1 )
    call COMM_vars8_init( RHOT_RK1, I_COMM_RHOT_RK1 )

    call COMM_vars8_init( DENS_RK2, I_COMM_DENS_RK2 )
    call COMM_vars8_init( MOMZ_RK2, I_COMM_MOMZ_RK2 )
    call COMM_vars8_init( MOMX_RK2, I_COMM_MOMX_RK2 )
    call COMM_vars8_init( MOMY_RK2, I_COMM_MOMY_RK2 )
    call COMM_vars8_init( RHOT_RK2, I_COMM_RHOT_RK2 )

    do iq = 1, QA
       I_COMM_RHOQ_t(iq) = 5 + iq
       I_COMM_QTRC(iq) = 5 + iq

       call COMM_vars8_init( RHOQ_t(:,:,:,iq), I_COMM_RHOQ_t(iq) )
       call COMM_vars8_init( QTRC  (:,:,:,iq), I_COMM_QTRC(iq) )
    end do

    call COMM_vars8_init( mflx_hi(:,:,:,ZDIR), I_COMM_mflx_z )
    call COMM_vars8_init( mflx_hi(:,:,:,XDIR), I_COMM_mflx_x )
    call COMM_vars8_init( mflx_hi(:,:,:,YDIR), I_COMM_mflx_y )

    DENS_RK1(:,:,:) = UNDEF
    MOMZ_RK1(:,:,:) = UNDEF
    MOMX_RK1(:,:,:) = UNDEF
    MOMY_RK1(:,:,:) = UNDEF
    RHOT_RK1(:,:,:) = UNDEF
    DENS_RK2(:,:,:) = UNDEF
    MOMZ_RK2(:,:,:) = UNDEF
    MOMX_RK2(:,:,:) = UNDEF
    MOMY_RK2(:,:,:) = UNDEF
    RHOT_RK2(:,:,:) = UNDEF

    mflx_hi(:,:,:,:) = UNDEF

    return
  end subroutine ATMOS_DYN_setup

  !-----------------------------------------------------------------------------
  !> Dynamical Process
  subroutine ATMOS_DYN( &
       DENS,    MOMZ,    MOMX,    MOMY,    RHOT,    QTRC,                                                    &
       DENS_av, MOMZ_av, MOMX_av, MOMY_av, RHOT_av, QTRC_av,                                                 &
       DENS_tp, MOMZ_tp, MOMX_tp, MOMY_tp, RHOT_tp, RHOQ_tp,                                                 &
       CDZ, CDX, CDY, FDZ, FDX, FDY,                                                                         &
       RCDZ, RCDX, RCDY, RFDZ, RFDX, RFDY,                                                                   &
       PHI, GSQRT,                                                                                           &
       J13G, J23G, J33G, MAPF,                                                                               &
       AQ_CV,                                                                                                &
       REF_dens, REF_pott, REF_qv, REF_pres,                                                                 &
       ND_COEF, ND_COEF_Q, ND_ORDER, ND_SFC_FACT, ND_USE_RS,                                                 &
       DAMP_DENS,       DAMP_VELZ,       DAMP_VELX,       DAMP_VELY,       DAMP_POTT,       DAMP_QTRC,       &
       DAMP_alpha_DENS, DAMP_alpha_VELZ, DAMP_alpha_VELX, DAMP_alpha_VELY, DAMP_alpha_POTT, DAMP_alpha_QTRC, &
       divdmp_coef,                                                                                          &
       FLAG_FCT_RHO, FLAG_FCT_MOMENTUM, FLAG_FCT_T,                                                          &
       FLAG_FCT_ALONG_STREAM,                                                                                &
       USE_AVERAGE,                                                                                          &
       DTSEC, DTSEC_ATMOS_DYN, NSTEP_ATMOS_DYN                                                               )
    use scale_const, only: &
       Rdry   => CONST_Rdry, &
       Rvap   => CONST_Rvap, &
       CVdry  => CONST_CVdry
    use scale_process, only: &
       PRC_HAS_W, &
       PRC_HAS_E, &
       PRC_HAS_S, &
       PRC_HAS_N
    use scale_comm, only: &
#ifdef CHECK_MASS
       COMM_datatype, &
       COMM_world,    &
#endif
       COMM_vars8, &
       COMM_wait
    use scale_gridtrans, only: &
       I_XYZ, &
       I_XYW, &
       I_UYZ, &
       I_XVZ, &
       I_XY,  &
       I_UY,  &
       I_XV
#ifdef CHECK_MASS
    use scale_grid_real, only: &
       vol => REAL_VOL
#endif
    use scale_atmos_dyn_common, only: &
       FACT_N,                     &
       FACT_F,                     &
       ATMOS_DYN_numfilter_coef,   &
       ATMOS_DYN_numfilter_coef_q, &
       ATMOS_DYN_fct
    use scale_atmos_dyn_rk, only: &
       ATMOS_DYN_rk
    use scale_atmos_boundary, only: &
       BND_QA, &
       BND_SMOOTHER_FACT => ATMOS_BOUNDARY_SMOOTHER_FACT
#if defined( HIST_TEND ) || defined( CHECK_MASS )
    use scale_history, only: &
       HIST_in
#endif
    implicit none

    real(RP), intent(inout) :: DENS(KA,IA,JA)
    real(RP), intent(inout) :: MOMZ(KA,IA,JA)
    real(RP), intent(inout) :: MOMX(KA,IA,JA)
    real(RP), intent(inout) :: MOMY(KA,IA,JA)
    real(RP), intent(inout) :: RHOT(KA,IA,JA)
    real(RP), intent(inout) :: QTRC(KA,IA,JA,QA)

    real(RP), intent(inout) :: DENS_av(KA,IA,JA)
    real(RP), intent(inout) :: MOMZ_av(KA,IA,JA)
    real(RP), intent(inout) :: MOMX_av(KA,IA,JA)
    real(RP), intent(inout) :: MOMY_av(KA,IA,JA)
    real(RP), intent(inout) :: RHOT_av(KA,IA,JA)
    real(RP), intent(inout) :: QTRC_av(KA,IA,JA,QA)

    real(RP), intent(in)    :: DENS_tp(KA,IA,JA)
    real(RP), intent(in)    :: MOMZ_tp(KA,IA,JA)
    real(RP), intent(in)    :: MOMX_tp(KA,IA,JA)
    real(RP), intent(in)    :: MOMY_tp(KA,IA,JA)
    real(RP), intent(in)    :: RHOT_tp(KA,IA,JA)
    real(RP), intent(in)    :: RHOQ_tp(KA,IA,JA,QA)

    real(RP), intent(in)    :: CDZ (KA)
    real(RP), intent(in)    :: CDX (IA)
    real(RP), intent(in)    :: CDY (JA)
    real(RP), intent(in)    :: FDZ (KA-1)
    real(RP), intent(in)    :: FDX (IA-1)
    real(RP), intent(in)    :: FDY (JA-1)
    real(RP), intent(in)    :: RCDZ(KA)
    real(RP), intent(in)    :: RCDX(IA)
    real(RP), intent(in)    :: RCDY(JA)
    real(RP), intent(in)    :: RFDZ(KA-1)
    real(RP), intent(in)    :: RFDX(IA-1)
    real(RP), intent(in)    :: RFDY(JA-1)

    real(RP), intent(in)    :: PHI  (KA,IA,JA)   !< geopotential
    real(RP), intent(in)    :: GSQRT(KA,IA,JA,7) !< vertical metrics {G}^1/2
    real(RP), intent(in)    :: J13G (KA,IA,JA,7) !< (1,3) element of Jacobian matrix
    real(RP), intent(in)    :: J23G (KA,IA,JA,7) !< (2,3) element of Jacobian matrix
    real(RP), intent(in)    :: J33G              !< (3,3) element of Jacobian matrix
    real(RP), intent(in)    :: MAPF (IA,JA,2,4)  !< map factor

    real(RP), intent(in)    :: AQ_CV(QQA)

    real(RP), intent(in)    :: REF_dens(KA,IA,JA)
    real(RP), intent(in)    :: REF_pott(KA,IA,JA)
    real(RP), intent(in)    :: REF_qv  (KA,IA,JA)
    real(RP), intent(in)    :: REF_pres(KA,IA,JA)   !< reference pressure
    real(RP), intent(in)    :: ND_COEF
    real(RP), intent(in)    :: ND_COEF_Q
    integer,  intent(in)    :: ND_ORDER
    real(RP), intent(in)    :: ND_SFC_FACT
    logical,  intent(in)    :: ND_USE_RS

    real(RP), intent(in)    :: DAMP_DENS(KA,IA,JA)
    real(RP), intent(in)    :: DAMP_VELZ(KA,IA,JA)
    real(RP), intent(in)    :: DAMP_VELX(KA,IA,JA)
    real(RP), intent(in)    :: DAMP_VELY(KA,IA,JA)
    real(RP), intent(in)    :: DAMP_POTT(KA,IA,JA)
    real(RP), intent(in)    :: DAMP_QTRC(KA,IA,JA,BND_QA)

    real(RP), intent(in)    :: DAMP_alpha_DENS(KA,IA,JA)
    real(RP), intent(in)    :: DAMP_alpha_VELZ(KA,IA,JA)
    real(RP), intent(in)    :: DAMP_alpha_VELX(KA,IA,JA)
    real(RP), intent(in)    :: DAMP_alpha_VELY(KA,IA,JA)
    real(RP), intent(in)    :: DAMP_alpha_POTT(KA,IA,JA)
    real(RP), intent(in)    :: DAMP_alpha_QTRC(KA,IA,JA,BND_QA)

    real(RP), intent(in)    :: divdmp_coef

    logical,  intent(in)    :: FLAG_FCT_RHO
    logical,  intent(in)    :: FLAG_FCT_MOMENTUM
    logical,  intent(in)    :: FLAG_FCT_T
    logical,  intent(in)    :: FLAG_FCT_ALONG_STREAM

    logical,  intent(in)    :: USE_AVERAGE

    real(DP), intent(in)    :: DTSEC
    real(DP), intent(in)    :: DTSEC_ATMOS_DYN
    integer , intent(in)    :: NSTEP_ATMOS_DYN

    ! for time integartion
    real(RP) :: DENS00  (KA,IA,JA) ! saved density before small step loop
    real(RP) :: DENS0   (KA,IA,JA) ! prognostic variables (+0   step)
    real(RP) :: MOMZ0   (KA,IA,JA) !
    real(RP) :: MOMX0   (KA,IA,JA) !
    real(RP) :: MOMY0   (KA,IA,JA) !
    real(RP) :: RHOT0   (KA,IA,JA) !

    ! diagnostic variables
    real(RP) :: QDRY (KA,IA,JA) ! dry air
    real(RP) :: Rtot (KA,IA,JA) ! total R
    real(RP) :: CVtot(KA,IA,JA) ! total CV

    real(RP) :: DENS_tq(KA,IA,JA)
    real(RP) :: diff(KA,IA,JA)
    real(RP) :: damp
#ifdef HIST_TEND
    real(RP) :: damp_t(KA,IA,JA)
#endif

    ! For tracer advection
    real(RP) :: mflx_av  (KA,IA,JA,3)     ! rho * vel(x,y,z) @ (u,v,w)-face average
    real(RP) :: tflx_hi  (KA,IA,JA,3)     ! rho * theta * vel(x,y,z) @ (u,v,w)-face high order
    real(RP) :: qflx_hi  (KA,IA,JA,3)  ! rho * vel(x,y,z) * phi @ (u,v,w)-face high order
    real(RP) :: qflx_lo  (KA,IA,JA,3)  ! rho * vel(x,y,z) * phi,  monotone flux
    real(RP) :: qflx_anti(KA,IA,JA,3)  ! anti-diffusive flux

    ! lateral boundary flux
#ifdef CHECK_MASS
    real(RP) :: mflx_lb_horizontal(KA)
    real(RP) :: allmflx_lb_horizontal(KA)
    real(RP) :: mflx_lb_total
    real(RP) :: mass_total
    real(RP) :: mass_total2
    real(RP) :: allmflx_lb_total
    real(RP) :: allmass_total
    real(RP) :: allmass_total2

    integer :: ierr
#endif

    real(RP) :: dt
    real(RP) :: dtrk

    integer  :: IIS, IIE
    integer  :: JJS, JJE
    integer  :: i, j, k, iq, step
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*) '*** Dynamics step: ', NSTEP_ATMOS_DYN, ' small steps'

    call PROF_rapstart("DYN_Preparation", 2)

#ifdef DEBUG
    DENS00  (:,:,:) = UNDEF
    DENS0   (:,:,:) = UNDEF
    MOMZ0   (:,:,:) = UNDEF
    MOMX0   (:,:,:) = UNDEF
    MOMY0   (:,:,:) = UNDEF
    RHOT0   (:,:,:) = UNDEF
    DENS_RK1(:,:,:) = UNDEF
    MOMZ_RK1(:,:,:) = UNDEF
    MOMX_RK1(:,:,:) = UNDEF
    MOMY_RK1(:,:,:) = UNDEF
    RHOT_RK1(:,:,:) = UNDEF
    DENS_RK2(:,:,:) = UNDEF
    MOMZ_RK2(:,:,:) = UNDEF
    MOMX_RK2(:,:,:) = UNDEF
    MOMY_RK2(:,:,:) = UNDEF
    RHOT_RK2(:,:,:) = UNDEF

    num_diff (:,:,:,:,:) = UNDEF

    mflx_hi(:,:,:,:) = UNDEF
    tflx_hi(:,:,:,:) = UNDEF

    qflx_hi(:,:,:,:) = UNDEF
    qflx_lo(:,:,:,:) = UNDEF
#endif

!OCL XFILL
    DENS00(:,:,:) = DENS(:,:,:)

    if ( USE_AVERAGE ) then
!OCL XFILL
       DENS_av(:,:,:) = 0.0_RP
!OCL XFILL
       MOMZ_av(:,:,:) = 0.0_RP
!OCL XFILL
       MOMX_av(:,:,:) = 0.0_RP
!OCL XFILL
       MOMY_av(:,:,:) = 0.0_RP
!OCL XFILL
       RHOT_av(:,:,:) = 0.0_RP
    endif

#ifndef DRY
!OCL XFILL
    mflx_av(:,:,:,:) = 0.0_RP

!OCL XFILL
    CVtot(:,:,:) = 0.0_RP
!OCL XFILL
    QDRY (:,:,:) = 1.0_RP
    do iq = QQS, QQE
       CVtot(:,:,:) = CVtot(:,:,:) + AQ_CV(iq) * QTRC(:,:,:,iq)
       QDRY (:,:,:) = QDRY (:,:,:) - QTRC(:,:,:,iq)
    enddo
    CVtot(:,:,:) = CVdry * QDRY(:,:,:) + CVtot(:,:,:)
    Rtot (:,:,:) = Rdry  * QDRY(:,:,:) + Rvap * QTRC(:,:,:,I_QV)
#endif

#ifdef HIST_TEND
    damp_t(:,:,:) = 0.0_RP
#endif

    call PROF_rapend  ("DYN_Preparation", 2)

    !###########################################################################
    ! Update DENS,MONZ,MOMX,MOMY,MOMZ,RHOT
    !###########################################################################

    call PROF_rapstart("DYN_Tendency", 2)

!OCL XFILL
    DENS_tq(:,:,:) = 0.0_RP

    do iq = 1, BND_QA

       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
!OCL XFILL
       do j = JS-1, JE+2
       do i = IS-1, IE+2
       do k = KS, KE
          diff(k,i,j) = QTRC(k,i,j,iq) - DAMP_QTRC(k,i,j,iq)
       enddo
       enddo
       enddo
       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
!OCL XFILL
       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          damp = - DAMP_alpha_QTRC(k,i,j,iq) &
               * ( diff(k,i,j) & ! rayleigh damping
                 - ( diff(k,i-1,j) + diff(k,i+1,j) + diff(k,i,j-1) + diff(k,i,j+1) - diff(k,i,j)*4.0_RP ) &
                   * 0.125_RP * BND_SMOOTHER_FACT ) ! horizontal smoother
          RHOQ_t(k,i,j,iq) = damp * DENS00(k,i,j)
#ifdef HIST_TEND
          damp_t(k,i,j) = damp
#endif
       enddo
       enddo
       enddo
#ifdef HIST_TEND
       call HIST_in(RHOQ_tp(:,:,:,iq), trim(AQ_NAME(iq))//'_t_phys',                         &
                    'tendency of '//trim(AQ_NAME(iq))//' due to physics',          'kg/kg/s' )
       call HIST_in(damp_t,            trim(AQ_NAME(iq))//'_t_damp',                         &
                    'tendency of '//trim(AQ_NAME(iq))//' due to rayleigh damping', 'kg/kg/s' )
#endif
!OCL XFILL
       do j = JS, JE
       do i = IS, IE
          RHOQ_t(   1:KS-1,i,j,iq) = 0.0_RP
          RHOQ_t(KE+1:KA  ,i,j,iq) = 0.0_RP
       enddo
       enddo

       call COMM_vars8( RHOQ_t(:,:,:,iq), I_COMM_RHOQ_t(iq) )

       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          DENS_tq(k,i,j) = DENS_tq(k,i,j) + RHOQ_t(k,i,j,iq)
       end do
       end do
       end do

       call COMM_wait ( RHOQ_t(:,:,:,iq), I_COMM_RHOQ_t(iq), .false. )

    end do

    !$omp parallel do private(i,j,k,iq) OMP_SCHEDULE_ collapse(3)
!OCL XFILL
    do iq = BND_QA+1, QA
       do j = 1, JA
       do i = 1, IA
       do k = 1, KA
          RHOQ_t(k,i,j,iq) = 0.0_RP
       enddo
       enddo
       enddo
    end do

    call PROF_rapend  ("DYN_Tendency", 2)

    call PROF_rapstart("DYN_Boundary", 2)

    if ( BND_W ) then
       !$omp parallel do private(j,k) OMP_SCHEDULE_ collapse(2)
       do j = JS, JE
       do k = KS, KE
          mflx_hi(k,IS-1,j,XDIR) = GSQRT(k,IS-1,j,I_UYZ) * MOMX(k,IS-1,j) / MAPF(IS-1,j,2,I_UY)
          tflx_hi(k,IS-1,j,XDIR) = mflx_hi(k,IS-1,j,XDIR) &
                                 * ( FACT_N * ( DAMP_POTT(k,IS  ,j)+DAMP_POTT(k,IS-1,j) ) &
                                   + FACT_F * ( DAMP_POTT(k,IS+1,j)+DAMP_POTT(k,IS-2,j) ) )
       enddo
       enddo
    end if
    if ( BND_E ) then
       !$omp parallel do private(j,k) OMP_SCHEDULE_ collapse(2)
       do j = JS, JE
       do k = KS, KE
          mflx_hi(k,IE,j,XDIR) = GSQRT(k,IE,j,I_UYZ) * MOMX(k,IE,j) / MAPF(IE,j,2,I_UY)
          tflx_hi(k,IE,j,XDIR) = mflx_hi(k,IE,j,XDIR) &
                               * ( FACT_N * ( DAMP_POTT(k,IE+1,j)+DAMP_POTT(k,IE  ,j) ) &
                                 + FACT_F * ( DAMP_POTT(k,IE+2,j)+DAMP_POTT(k,IE-1,j) ) )
       enddo
       enddo
    end if
    if ( BND_S ) then
       !$omp parallel do private(i,k) OMP_SCHEDULE_ collapse(2)
       do i = IS, IE
       do k = KS, KE
          mflx_hi(k,i,JS-1,YDIR) = GSQRT(k,i,JS-1,I_XVZ) * MOMY(k,i,JS-1) / MAPF(i,JS-1,1,I_XV)
          tflx_hi(k,i,JS-1,YDIR) = mflx_hi(k,i,JS-1,YDIR) &
                                 * ( FACT_N * ( DAMP_POTT(k,i,JS  )+DAMP_POTT(k,i,JS-1) ) &
                                   + FACT_F * ( DAMP_POTT(k,i,JS+1)+DAMP_POTT(k,i,JS-2) ) )
       enddo
       enddo
    end if
    if ( BND_N ) then
       !$omp parallel do private(i,k) OMP_SCHEDULE_ collapse(2)
       do i = IS, IE
       do k = KS, KE
          mflx_hi(k,i,JE,YDIR) = GSQRT(k,i,JE,I_XVZ) * MOMY(k,i,JE) / MAPF(i,JE,1,I_XV)
          tflx_hi(k,i,JE,YDIR) = mflx_hi(k,i,JE,YDIR) &
                               * ( FACT_N * ( DAMP_POTT(k,i,JE+1)+DAMP_POTT(k,i,JE  ) ) &
                                 + FACT_F * ( DAMP_POTT(k,i,JE+2)+DAMP_POTT(k,i,JE-1) ) )
       enddo
       enddo
    end if

    call PROF_rapend  ("DYN_Boundary", 2)


    do step = 1, NSTEP_ATMOS_DYN

       !-----< prepare tendency >-----

       call PROF_rapstart("DYN_Tendency", 2)

       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
!OCL XFILL
       do j = JS-1, JE+2
       do i = IS-1, IE+2
       do k = KS, KE
          diff(k,i,j) = DENS(k,i,j) - DAMP_DENS(k,i,j)
       enddo
       enddo
       enddo
       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
!OCL XFILL
       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          damp = - DAMP_alpha_DENS(k,i,j) &
               * ( diff(k,i,j) & ! rayleigh damping
                 - ( diff(k,i-1,j) + diff(k,i+1,j) + diff(k,i,j-1) + diff(k,i,j+1) - diff(k,i,j)*4.0_RP ) &
                 * 0.125_RP * BND_SMOOTHER_FACT ) & ! horizontal smoother
                 + DENS_tq(k,i,j) ! dencity change due to rayleigh damping for tracers
          DENS_t(k,i,j) = DENS_tp(k,i,j) & ! tendency from physical step
                        + damp
#ifdef HIST_TEND
          damp_t(k,i,j) = damp
#endif
       enddo
       enddo
       enddo
!OCL XFILL
       do j = JS, JE
       do i = IS, IE
          DENS_t(   1:KS-1,i,j) = 0.0_RP
          DENS_t(KE+1:KA  ,i,j) = 0.0_RP
       enddo
       enddo
       call COMM_vars8( DENS_t(:,:,:), I_COMM_DENS_t )
#ifdef HIST_TEND
       call HIST_in(DENS_tp, 'DENS_t_phys', 'tendency of dencity due to physics',          'kg/m3/s' )
       call HIST_in(damp_t,  'DENS_t_damp', 'tendency of dencity due to rayleigh damping', 'kg/m3/s' )
#endif

       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
!OCL XFILL
       do j = JS-1, JE+1
       do i = IS-1, IE+1
       do k = KS, KE-1
          diff(k,i,j) = MOMZ(k,i,j) - DAMP_VELZ(k,i,j) * ( DENS(k,i,j)+DENS(k+1,i,j) ) * 0.5_RP
       enddo
       enddo
       enddo
       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
!OCL XFILL
       do j = JS, JE
       do i = IS, IE
       do k = KS, KE-1
          damp = - DAMP_alpha_VELZ(k,i,j) &
               * ( diff(k,i,j) & ! rayleigh damping
                 - ( diff(k,i-1,j) + diff(k,i+1,j) + diff(k,i,j-1) + diff(k,i,j+1) - diff(k,i,j)*4.0_RP ) &
                 * 0.125_RP * BND_SMOOTHER_FACT ) ! horizontal smoother
          MOMZ_t(k,i,j) = MOMZ_tp(k,i,j) & ! tendency from physical step
                        + damp
#ifdef HIST_TEND
          damp_t(k,i,j) = damp
#endif
       enddo
       enddo
       enddo
!OCL XFILL
       do j = JS, JE
       do i = IS, IE
          MOMZ_t( 1:KS-1,i,j) = 0.0_RP
          MOMZ_t(KE:KA  ,i,j) = 0.0_RP
       enddo
       enddo
       call COMM_vars8( MOMZ_t(:,:,:), I_COMM_MOMZ_t )
#ifdef HIST_TEND
       call HIST_in(MOMZ_tp, 'MOMZ_t_phys', 'tendency of momentum z due to physics',          'kg/m2/s2', zdim='half' )
       call HIST_in(damp_t,  'MOMZ_t_damp', 'tendency of momentum z due to rayleigh damping', 'kg/m2/s2', zdim='half' )
#endif

       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
!OCL XFILL
       do j = JS-1, JE+2
       do i = IS-2, IE+1
       do k = KS, KE
          diff(k,i,j) = MOMX(k,i,j) - DAMP_VELX(k,i,j) * ( DENS(k,i,j)+DENS(k,i+1,j) ) * 0.5_RP
       enddo
       enddo
       enddo
!OCL XFILL
       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          damp = - DAMP_alpha_VELX(k,i,j) &
                 * ( diff(k,i,j) & ! rayleigh damping
                   - ( diff(k,i-1,j) + diff(k,i+1,j) + diff(k,i,j-1) + diff(k,i,j+1) - diff(k,i,j)*4.0_RP ) &
                   * 0.125_RP * BND_SMOOTHER_FACT ) ! horizontal smoother
          MOMX_t(k,i,j) = MOMX_tp(k,i,j) & ! tendency from physical step
                        + damp
#ifdef HIST_TEND
          damp_t(k,i,j) = damp
#endif
       enddo
       enddo
       enddo
!OCL XFILL
       do j = JS, JE
       do i = IS, IE
          MOMX_t(   1:KS-1,i,j) = 0.0_RP
          MOMX_t(KE+1:KA  ,i,j) = 0.0_RP
       enddo
       enddo
       call COMM_vars8( MOMX_t(:,:,:), I_COMM_MOMX_t )
#ifdef HIST_TEND
       call HIST_in(MOMX_tp, 'MOMX_t_phys', 'tendency of momentum x due to physics',          'kg/m2/s2', xdim='half' )
       call HIST_in(damp_t,  'MOMX_t_damp', 'tendency of momentum x due to rayleigh damping', 'kg/m2/s2', xdim='half' )
#endif

       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
!OCL XFILL
       do j = JS-2, JE+1
       do i = IS-1, IE+2
       do k = KS, KE
          diff(k,i,j) = MOMY(k,i,j) - DAMP_VELY(k,i,j) * ( DENS(k,i,j)+DENS(k,i,j+1) ) * 0.5_RP
       enddo
       enddo
       enddo
!OCL XFILL
       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          damp = - DAMP_alpha_VELY(k,i,j) &
                 * ( diff(k,i,j) & ! rayleigh damping
                   - ( diff(k,i-1,j) + diff(k,i+1,j) + diff(k,i,j-1) + diff(k,i,j+1) - diff(k,i,j)*4.0_RP ) &
                   * 0.125_RP * BND_SMOOTHER_FACT ) ! horizontal smoother
          MOMY_t(k,i,j) = MOMY_tp(k,i,j) & ! tendency from physical step
                        + damp
#ifdef HIST_TEND
          damp_t(k,i,j) = damp
#endif
       enddo
       enddo
       enddo
!OCL XFILL
       do j = JS, JE
       do i = IS, IE
          MOMY_t(   1:KS-1,i,j) = 0.0_RP
          MOMY_t(KE+1:KA  ,i,j) = 0.0_RP
       enddo
       enddo
       call COMM_vars8( MOMY_t(:,:,:), I_COMM_MOMY_t )
#ifdef HIST_TEND
       call HIST_in(MOMY_tp, 'MOMY_t_phys', 'tendency of momentum y due to physics',          'kg/m2/s2', ydim='half' )
       call HIST_in(damp_t,  'MOMY_t_damp', 'tendency of momentum y due to rayleigh damping', 'kg/m2/s2', ydim='half' )
#endif

       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
!OCL XFILL
       do j = JS-1, JE+2
       do i = IS-1, IE+2
       do k = KS, KE
          diff(k,i,j) = RHOT(k,i,j) - DAMP_POTT(k,i,j) * DENS(k,i,j)
       enddo
       enddo
       enddo
!OCL XFILL
       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          damp = - DAMP_alpha_POTT(k,i,j) &
                 * ( diff(k,i,j) & ! rayleigh damping
                   - ( diff(k,i-1,j) + diff(k,i+1,j) + diff(k,i,j-1) + diff(k,i,j+1) - diff(k,i,j)*4.0_RP ) &
                   * 0.125_RP * BND_SMOOTHER_FACT ) ! horizontal smoother
          RHOT_t(k,i,j) = RHOT_tp(k,i,j) & ! tendency from physical step
                        + damp
#ifdef HIST_TEND
          damp_t(k,i,j) = damp
#endif
       enddo
       enddo
       enddo
!OCL XFILL
       do j = JS, JE
       do i = IS, IE
          RHOT_t(   1:KS-1,i,j) = 0.0_RP
          RHOT_t(KE+1:KA  ,i,j) = 0.0_RP
       enddo
       enddo
       call COMM_vars8( RHOT_t(:,:,:), I_COMM_RHOT_t )
#ifdef HIST_TEND
       call HIST_in(RHOT_tp, 'RHOT_t_phys', 'tendency of rho*theta temperature due to physics',          'K kg/m3/s' )
       call HIST_in(damp_t,  'RHOT_t_damp', 'tendency of rho*theta temperature due to rayleigh damping', 'K kg/m3/s' )
#endif

       call COMM_wait ( DENS_t(:,:,:), I_COMM_DENS_t, .false. )
       call COMM_wait ( MOMZ_t(:,:,:), I_COMM_MOMZ_t, .false. )
       call COMM_wait ( MOMX_t(:,:,:), I_COMM_MOMX_t, .false. )
       call COMM_wait ( MOMY_t(:,:,:), I_COMM_MOMY_t, .false. )
       call COMM_wait ( RHOT_t(:,:,:), I_COMM_RHOT_t, .false. )

       call PROF_rapend  ("DYN_Tendency", 2)

       call PROF_rapstart("DYN_Numfilter", 2)

       dt = real(DTSEC_ATMOS_DYN,kind=RP)

       !-----< prepare numerical diffusion coefficient >-----

       if ( ND_COEF == 0.0_RP ) then
!OCL XFILL
          num_diff(:,:,:,:,:) = 0.0_RP
       else
          call ATMOS_DYN_numfilter_coef( num_diff(:,:,:,:,:),                    & ! [OUT]
                                         DENS, MOMZ, MOMX, MOMY, RHOT,           & ! [IN]
                                         CDZ, CDX, CDY, FDZ, FDX, FDY, dt,       & ! [IN]
                                         REF_dens, REF_pott,                     & ! [IN]
                                         ND_COEF, ND_ORDER, ND_SFC_FACT, ND_USE_RS ) ! [IN]
       endif

       call PROF_rapend  ("DYN_Numfilter", 2)

       !------------------------------------------------------------------------
       ! Start RK
       !------------------------------------------------------------------------

       !##### SAVE #####
!OCL XFILL
       DENS0(:,:,:) = DENS(:,:,:)
!OCL XFILL
       MOMZ0(:,:,:) = MOMZ(:,:,:)
!OCL XFILL
       MOMX0(:,:,:) = MOMX(:,:,:)
!OCL XFILL
       MOMY0(:,:,:) = MOMY(:,:,:)
!OCL XFILL
       RHOT0(:,:,:) = RHOT(:,:,:)

       call PROF_rapstart("DYN_RK", 2)

       !##### RK1 : PROG0,PROG->PROG_RK1 #####

       dtrk = real(DTSEC_ATMOS_DYN,kind=RP) / 3.0_RP

       call ATMOS_DYN_rk( DENS_RK1, MOMZ_RK1, MOMX_RK1, MOMY_RK1, RHOT_RK1, & ! (out)
                          mflx_hi,  tflx_hi,                                & ! (out)
                          DENS0,    MOMZ0,    MOMX0,    MOMY0,    RHOT0,    & ! (in)
                          DENS,     MOMZ,     MOMX,     MOMY,     RHOT,     & ! (in)
                          DENS_t,   MOMZ_t,   MOMX_t,   MOMY_t,   RHOT_t,   & ! (in)
                          Rtot, CVtot, CORIOLI,                             & ! (in)
                          num_diff, divdmp_coef,                            & ! (in)
                          FLAG_FCT_RHO, FLAG_FCT_MOMENTUM, FLAG_FCT_T,      & ! (in)
                          FLAG_FCT_ALONG_STREAM,                            & ! (in)
                          CDZ, FDZ, FDX, FDY,                               & ! (in)
                          RCDZ, RCDX, RCDY, RFDZ, RFDX, RFDY,               & ! (in)
                          PHI, GSQRT, J13G, J23G, J33G, MAPF,               & ! (in)
                          REF_pres, REF_dens,                               & ! (in)
                          BND_W, BND_E, BND_S, BND_N,                       & ! (in)
                          dtrk, dt                                          ) ! (in)
       call PROF_rapend  ("DYN_RK", 2)

       call PROF_rapstart("DYN_Boundary", 2)

       if ( BND_W ) then
          !$omp parallel do private(j,k) OMP_SCHEDULE_ collapse(2)
!OCL XFILL
          do j = 1, JA
          do i = 1, IS-1
          do k = KS, KE
             DENS_RK1(k,i,j) = DENS0(k,i,j)
             MOMZ_RK1(k,i,j) = MOMZ0(k,i,j)
             MOMX_RK1(k,i,j) = MOMX0(k,i,j)
             MOMY_RK1(k,i,j) = MOMY0(k,i,j)
             RHOT_RK1(k,i,j) = RHOT0(k,i,j)
          enddo
          enddo
          enddo
       end if
       if ( BND_E ) then
          !$omp parallel do private(j,k) OMP_SCHEDULE_ collapse(2)
!OCL XFILL
          do j = 1, JA
          do i = IE+1, IA
          do k = KS, KE
             DENS_RK1(k,i,j) = DENS0(k,i,j)
             MOMZ_RK1(k,i,j) = MOMZ0(k,i,j)
             MOMX_RK1(k,i,j) = MOMX0(k,i,j)
             MOMY_RK1(k,i,j) = MOMY0(k,i,j)
             RHOT_RK1(k,i,j) = RHOT0(k,i,j)
          enddo
          enddo
          enddo
          !$omp parallel do private(j,k) OMP_SCHEDULE_ collapse(2)
!OCL XFILL
          do j = 1, JA
          do k = KS, KE
             MOMX_RK1(k,IE,j) = MOMX0(k,IE,j)
          enddo
          enddo
       end if
       if ( BND_S ) then
          !$omp parallel do private(j,k) OMP_SCHEDULE_ collapse(2)
!OCL XFILL
          do j = 1, JS-1
          do i = 1, IA
          do k = KS, KE
             DENS_RK1(k,i,j) = DENS0(k,i,j)
             MOMZ_RK1(k,i,j) = MOMZ0(k,i,j)
             MOMX_RK1(k,i,j) = MOMX0(k,i,j)
             MOMY_RK1(k,i,j) = MOMY0(k,i,j)
             RHOT_RK1(k,i,j) = RHOT0(k,i,j)
          enddo
          enddo
          enddo
       end if
       if ( BND_N ) then
          !$omp parallel do private(j,k) OMP_SCHEDULE_ collapse(2)
!OCL XFILL
          do j = JE+1, JA
          do i = 1, IA
          do k = KS, KE
             DENS_RK1(k,i,j) = DENS0(k,i,j)
             MOMZ_RK1(k,i,j) = MOMZ0(k,i,j)
             MOMX_RK1(k,i,j) = MOMX0(k,i,j)
             MOMY_RK1(k,i,j) = MOMY0(k,i,j)
             RHOT_RK1(k,i,j) = RHOT0(k,i,j)
          enddo
          enddo
          enddo
          !$omp parallel do private(i,k) OMP_SCHEDULE_ collapse(2)
!OCL XFILL
          do i = 1, IA
          do k = KS, KE
             MOMY_RK1(k,i,JE) = MOMY0(k,i,JE)
          enddo
          enddo
       end if

       call PROF_rapend  ("DYN_Boundary", 2)

       call COMM_vars8( DENS_RK1(:,:,:), I_COMM_DENS_RK1 )
       call COMM_vars8( MOMZ_RK1(:,:,:), I_COMM_MOMZ_RK1 )
       call COMM_vars8( MOMX_RK1(:,:,:), I_COMM_MOMX_RK1 )
       call COMM_vars8( MOMY_RK1(:,:,:), I_COMM_MOMY_RK1 )
       call COMM_vars8( RHOT_RK1(:,:,:), I_COMM_RHOT_RK1 )
       call COMM_wait ( DENS_RK1(:,:,:), I_COMM_DENS_RK1, .false. )
       call COMM_wait ( MOMZ_RK1(:,:,:), I_COMM_MOMZ_RK1, .false. )
       call COMM_wait ( MOMX_RK1(:,:,:), I_COMM_MOMX_RK1, .false. )
       call COMM_wait ( MOMY_RK1(:,:,:), I_COMM_MOMY_RK1, .false. )
       call COMM_wait ( RHOT_RK1(:,:,:), I_COMM_RHOT_RK1, .false. )

       !##### RK2 : PROG0,PROG_RK2->PROG_RK3 #####

       call PROF_rapstart("DYN_RK", 2)

       dtrk = real(DTSEC_ATMOS_DYN,kind=RP) / 2.0_RP

       call ATMOS_DYN_rk( DENS_RK2, MOMZ_RK2, MOMX_RK2, MOMY_RK2, RHOT_RK2, & ! (out)
                          mflx_hi,  tflx_hi,                                & ! (out)
                          DENS0,    MOMZ0,    MOMX0,    MOMY0,    RHOT0,    & ! (in)
                          DENS_RK1, MOMZ_RK1, MOMX_RK1, MOMY_RK1, RHOT_RK1, & ! (in)
                          DENS_t,   MOMZ_t,   MOMX_t,   MOMY_t,   RHOT_t,   & ! (in)
                          Rtot, CVtot, CORIOLI,                             & ! (in)
                          num_diff, divdmp_coef,                            & ! (in)
                          FLAG_FCT_RHO, FLAG_FCT_MOMENTUM, FLAG_FCT_T,      & ! (in)
                          FLAG_FCT_ALONG_STREAM,                            & ! (in)
                          CDZ, FDZ, FDX, FDY,                               & ! (in)
                          RCDZ, RCDX, RCDY, RFDZ, RFDX, RFDY,               & ! (in)
                          PHI, GSQRT, J13G, J23G, J33G, MAPF,               & ! (in)
                          REF_pres, REF_dens,                               & ! (in)
                          BND_W, BND_E, BND_S, BND_N,                       & ! (in)
                          dtrk, dt                                          ) ! (in)

       call PROF_rapend  ("DYN_RK", 2)

       call PROF_rapstart("DYN_Boundary", 2)

       if ( BND_W ) then
          !$omp parallel do private(j,k) OMP_SCHEDULE_ collapse(2)
!OCL XFILL
          do j = 1, JA
          do i = 1, IS-1
          do k = KS, KE
             DENS_RK2(k,i,j) = DENS0(k,i,j)
             MOMZ_RK2(k,i,j) = MOMZ0(k,i,j)
             MOMX_RK2(k,i,j) = MOMX0(k,i,j)
             MOMY_RK2(k,i,j) = MOMY0(k,i,j)
             RHOT_RK2(k,i,j) = RHOT0(k,i,j)
          enddo
          enddo
          enddo
       end if
       if ( BND_E ) then
          !$omp parallel do private(j,k) OMP_SCHEDULE_ collapse(2)
!OCL XFILL
          do j = 1, JA
          do i = IE+1, IA
          do k = KS, KE
             DENS_RK2(k,i,j) = DENS0(k,i,j)
             MOMZ_RK2(k,i,j) = MOMZ0(k,i,j)
             MOMX_RK2(k,i,j) = MOMX0(k,i,j)
             MOMY_RK2(k,i,j) = MOMY0(k,i,j)
             RHOT_RK2(k,i,j) = RHOT0(k,i,j)
          enddo
          enddo
          enddo
          !$omp parallel do private(j,k) OMP_SCHEDULE_ collapse(2)
!OCL XFILL
          do j = 1, JA
          do k = KS, KE
             MOMX_RK2(k,IE,j) = MOMX0(k,IE,j)
          enddo
          enddo
       end if
       if ( BND_S ) then
          !$omp parallel do private(j,k) OMP_SCHEDULE_ collapse(2)
!OCL XFILL
          do j = 1, JS-1
          do i = 1, IA
          do k = KS, KE
             DENS_RK2(k,i,j) = DENS0(k,i,j)
             MOMZ_RK2(k,i,j) = MOMZ0(k,i,j)
             MOMX_RK2(k,i,j) = MOMX0(k,i,j)
             MOMY_RK2(k,i,j) = MOMY0(k,i,j)
             RHOT_RK2(k,i,j) = RHOT0(k,i,j)
          enddo
          enddo
          enddo
       end if
       if ( BND_N ) then
          !$omp parallel do private(j,k) OMP_SCHEDULE_ collapse(2)
!OCL XFILL
          do j = JE+1, JA
          do i = 1, IA
          do k = KS, KE
             DENS_RK2(k,i,j) = DENS0(k,i,j)
             MOMZ_RK2(k,i,j) = MOMZ0(k,i,j)
             MOMX_RK2(k,i,j) = MOMX0(k,i,j)
             MOMY_RK2(k,i,j) = MOMY0(k,i,j)
             RHOT_RK2(k,i,j) = RHOT0(k,i,j)
          enddo
          enddo
          enddo
          !$omp parallel do private(i,k) OMP_SCHEDULE_ collapse(2)
!OCL XFILL
          do i = 1, IA
          do k = KS, KE
             MOMY_RK2(k,i,JE) = MOMY0(k,i,JE)
          enddo
          enddo
       end if

       call PROF_rapend  ("DYN_Boundary", 2)

       call COMM_vars8( DENS_RK2(:,:,:), I_COMM_DENS_RK2 )
       call COMM_vars8( MOMZ_RK2(:,:,:), I_COMM_MOMZ_RK2 )
       call COMM_vars8( MOMX_RK2(:,:,:), I_COMM_MOMX_RK2 )
       call COMM_vars8( MOMY_RK2(:,:,:), I_COMM_MOMY_RK2 )
       call COMM_vars8( RHOT_RK2(:,:,:), I_COMM_RHOT_RK2 )
       call COMM_wait ( DENS_RK2(:,:,:), I_COMM_DENS_RK2, .false. )
       call COMM_wait ( MOMZ_RK2(:,:,:), I_COMM_MOMZ_RK2, .false. )
       call COMM_wait ( MOMX_RK2(:,:,:), I_COMM_MOMX_RK2, .false. )
       call COMM_wait ( MOMY_RK2(:,:,:), I_COMM_MOMY_RK2, .false. )
       call COMM_wait ( RHOT_RK2(:,:,:), I_COMM_RHOT_RK2, .false. )

       !##### RK3 : PROG0,PROG_RK3->PROG #####

       call PROF_rapstart("DYN_RK", 2)

       dtrk = real(DTSEC_ATMOS_DYN,kind=RP)

       call ATMOS_DYN_rk( DENS,     MOMZ,     MOMX,     MOMY,     RHOT,     & ! (out)
                          mflx_hi,  tflx_hi,                                & ! (out)
                          DENS0,    MOMZ0,    MOMX0,    MOMY0,    RHOT0,    & ! (in)
                          DENS_RK2, MOMZ_RK2, MOMX_RK2, MOMY_RK2, RHOT_RK2, & ! (in)
                          DENS_t,   MOMZ_t,   MOMX_t,   MOMY_t,   RHOT_t,   & ! (in)
                          Rtot, CVtot, CORIOLI,                             & ! (in)
                          num_diff, divdmp_coef,                            & ! (in)
                          FLAG_FCT_RHO, FLAG_FCT_MOMENTUM, FLAG_FCT_T,      & ! (in)
                          FLAG_FCT_ALONG_STREAM,                            & ! (in)
                          CDZ, FDZ, FDX, FDY,                               & ! (in)
                          RCDZ, RCDX, RCDY, RFDZ, RFDX, RFDY,               & ! (in)
                          PHI, GSQRT, J13G, J23G, J33G, MAPF,               & ! (in)
                          REF_pres, REF_dens,                               & ! (in)
                          BND_W, BND_E, BND_S, BND_N,                       & ! (in)
                          dtrk, dt                                          ) ! (in)

       call PROF_rapend  ("DYN_RK", 2)

#ifdef CHECK_MASS
       call HIST_in(mflx_hi(:,:,:,ZDIR), 'MFLXZ', 'momentum flux of z-direction', 'kg/m2/s', zdim='half' )
       call HIST_in(mflx_hi(:,:,:,XDIR), 'MFLXX', 'momentum flux of x-direction', 'kg/m2/s', xdim='half' )
       call HIST_in(mflx_hi(:,:,:,YDIR), 'MFLXY', 'momentum flux of y-direction', 'kg/m2/s', ydim='half' )

       call HIST_in(tflx_hi(:,:,:,ZDIR), 'TFLXZ', 'potential temperature flux of z-direction', 'K*kg/m2/s', zdim='half' )
       call HIST_in(tflx_hi(:,:,:,XDIR), 'TFLXX', 'potential temperature flux of x-direction', 'K*kg/m2/s', xdim='half' )
       call HIST_in(tflx_hi(:,:,:,YDIR), 'TFLXY', 'potential temperature flux of y-direction', 'K*kg/m2/s', ydim='half' )

       mflx_lb_total            = 0.0_RP
       mflx_lb_horizontal(:)    = 0.0_RP
       allmflx_lb_horizontal(:) = 0.0_RP

       if ( BND_W ) then ! for western boundary
          i = IS
          do j = JS, JE
          do k = KS, KE
            mflx_lb_total = mflx_lb_total + mflx_hi(k,i-1,j,XDIR) * RCDX(i) * vol(k,i,j) &
                                          * MAPF(i,j,1,I_XY) * MAPF(i,j,2,I_XY) / GSQRT(k,i,j,I_XYZ) * dt
            mflx_lb_horizontal(k) = mflx_lb_horizontal(k) + mflx_hi(k,i-1,j,XDIR) * RCDX(i) * vol(k,i,j) &
                                                          * MAPF(i,j,1,I_XY) * MAPF(i,j,2,I_XY) / GSQRT(k,i,j,I_XYZ) * dt

          end do
          end do
       end if
       if ( BND_E ) then ! for eastern boundary
          i = IE
          do j = JS, JE
          do k = KS, KE
            mflx_lb_total = mflx_lb_total - mflx_hi(k,i,j,XDIR) * RCDX(i) * vol(k,i,j) &
                                          * MAPF(i,j,1,I_XY) * MAPF(i,j,2,I_XY) / GSQRT(k,i,j,I_XYZ) * dt
            mflx_lb_horizontal(k) = mflx_lb_horizontal(k) - mflx_hi(k,i,j,XDIR) * RCDX(i) * vol(k,i,j) &
                                                          * MAPF(i,j,1,I_XY) * MAPF(i,j,2,I_XY) / GSQRT(k,i,j,I_XYZ) * dt
          end do
          end do
       end if
       if ( BND_S ) then ! for sourthern boundary
          j = JS
          do i = IS, IE
          do k = KS, KE
            mflx_lb_total = mflx_lb_total + mflx_hi(k,i,j-1,YDIR) * RCDY(j) * vol(k,i,j) &
                                          * MAPF(i,j,1,I_XY) * MAPF(i,j,2,I_XY) / GSQRT(k,i,j,I_XYZ) * dt
            mflx_lb_horizontal(k) = mflx_lb_horizontal(k) + mflx_hi(k,i,j-1,YDIR) * RCDY(j) * vol(k,i,j) &
                                                          * MAPF(i,j,1,I_XY) * MAPF(i,j,2,I_XY) / GSQRT(k,i,j,I_XYZ) * dt
          end do
          end do
       end if
       if ( BND_N ) then ! for northern boundary
          j = JE
          do i = IS, IE
          do k = KS, KE
            mflx_lb_total = mflx_lb_total - mflx_hi(k,i,j,YDIR) * RCDY(j) * vol(k,i,j) &
                                          * MAPF(i,j,1,I_XY) * MAPF(i,j,2,I_XY) / GSQRT(k,i,j,I_XYZ) * dt
            mflx_lb_horizontal(k) = mflx_lb_horizontal(k) - mflx_hi(k,i,j,YDIR) * RCDY(j) * vol(k,i,j) &
                                                          * MAPF(i,j,1,I_XY) * MAPF(i,j,2,I_XY) / GSQRT(k,i,j,I_XYZ) * dt
          end do
          end do
       end if

       mass_total  = 0.0_RP
       mass_total2 = 0.0_RP

       ! check total mass in the inner region
       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          mass_total  = mass_total  + DENS     (k,i,j) * vol(k,i,j)
          mass_total2 = mass_total2 + DAMP_DENS(k,i,j) * vol(k,i,j)
       end do
       end do
       end do

       call MPI_Allreduce( mflx_lb_total,        &
                           allmflx_lb_total,     &
                           1,                    &
                           COMM_datatype,        &
                           MPI_SUM,              &
                           COMM_world,           &
                           ierr                  )

       if( IO_L ) write(IO_FID_LOG,'(A,1x,i1,1x,2PE24.17)') 'total mflx_lb:', step, allmflx_lb_total

       call MPI_Allreduce( mass_total,           &
                           allmass_total,        &
                           1,                    &
                           COMM_datatype,        &
                           MPI_SUM,              &
                           COMM_world,           &
                           ierr                  )

       if( IO_L ) write(IO_FID_LOG,'(A,1x,i1,1x,2PE24.17)') 'total mass   :', step, allmass_total

       call MPI_Allreduce( mass_total2,          &
                           allmass_total2,       &
                           1,                    &
                           COMM_datatype,        &
                           MPI_SUM,              &
                           COMM_world,           &
                           ierr                  )

       if( IO_L ) write(IO_FID_LOG,'(A,1x,i1,1x,2PE24.17)') 'total mass2  :', step, allmass_total2

       call MPI_Allreduce( mflx_lb_horizontal(KS:KE),    &
                           allmflx_lb_horizontal(KS:KE), &
                           KMAX,                         &
                           COMM_datatype,                &
                           MPI_SUM,                      &
                           COMM_world,                   &
                           ierr                          )

       call HIST_in(allmflx_lb_horizontal(:), 'ALLMOM_lb_hz',                           &
                    'horizontally total momentum flux from lateral boundary', 'kg/m2/s' )
#endif

       do j  = JS, JE
       do i  = IS, IE
          DENS(   1:KS-1,i,j) = DENS(KS,i,j)
          MOMZ(   1:KS-1,i,j) = MOMZ(KS,i,j)
          MOMX(   1:KS-1,i,j) = MOMX(KS,i,j)
          MOMY(   1:KS-1,i,j) = MOMY(KS,i,j)
          RHOT(   1:KS-1,i,j) = RHOT(KS,i,j)
          DENS(KE+1:KA,  i,j) = DENS(KE,i,j)
          MOMZ(KE+1:KA,  i,j) = MOMZ(KE,i,j)
          MOMX(KE+1:KA,  i,j) = MOMX(KE,i,j)
          MOMY(KE+1:KA,  i,j) = MOMY(KE,i,j)
          RHOT(KE+1:KA,  i,j) = RHOT(KE,i,j)
       enddo
       enddo

       call COMM_vars8( DENS(:,:,:), I_COMM_DENS )
       call COMM_vars8( MOMZ(:,:,:), I_COMM_MOMZ )
       call COMM_vars8( MOMX(:,:,:), I_COMM_MOMX )
       call COMM_vars8( MOMY(:,:,:), I_COMM_MOMY )
       call COMM_vars8( RHOT(:,:,:), I_COMM_RHOT )
       call COMM_wait ( DENS(:,:,:), I_COMM_DENS, .false. )
       call COMM_wait ( MOMZ(:,:,:), I_COMM_MOMZ, .false. )
       call COMM_wait ( MOMX(:,:,:), I_COMM_MOMX, .false. )
       call COMM_wait ( MOMY(:,:,:), I_COMM_MOMY, .false. )
       call COMM_wait ( RHOT(:,:,:), I_COMM_RHOT, .false. )

       if ( USE_AVERAGE ) then
          DENS_av(:,:,:) = DENS_av(:,:,:) + DENS(:,:,:)
          MOMZ_av(:,:,:) = MOMZ_av(:,:,:) + MOMZ(:,:,:)
          MOMX_av(:,:,:) = MOMX_av(:,:,:) + MOMX(:,:,:)
          MOMY_av(:,:,:) = MOMY_av(:,:,:) + MOMY(:,:,:)
          RHOT_av(:,:,:) = RHOT_av(:,:,:) + RHOT(:,:,:)
       endif

#ifndef DRY
       mflx_av(:,:,:,:) = mflx_av(:,:,:,:) + mflx_hi(:,:,:,:)
#endif

       do iq = 1, QA

          do j = JS, JE
          do i = IS, IE
          do k = KS, KE
             QTRC(k,i,j,iq) = max( QTRC(k,i,j,iq) + RHOQ_tp(k,i,j,iq) * dt / DENS00(k,i,j), 0.0_RP )
          end do
          end do
          end do

          call COMM_vars8( QTRC(:,:,:,iq), I_COMM_QTRC(iq) )
          call COMM_wait ( QTRC(:,:,:,iq), I_COMM_QTRC(iq), .false. )

          ! TODO: mass and energy conservation should be considered
       end do

    enddo ! dynamical steps

    if ( USE_AVERAGE ) then
       DENS_av(:,:,:) = DENS_av(:,:,:) / real(NSTEP_ATMOS_DYN,kind=RP)
       MOMZ_av(:,:,:) = MOMZ_av(:,:,:) / real(NSTEP_ATMOS_DYN,kind=RP)
       MOMX_av(:,:,:) = MOMX_av(:,:,:) / real(NSTEP_ATMOS_DYN,kind=RP)
       MOMY_av(:,:,:) = MOMY_av(:,:,:) / real(NSTEP_ATMOS_DYN,kind=RP)
       RHOT_av(:,:,:) = RHOT_av(:,:,:) / real(NSTEP_ATMOS_DYN,kind=RP)
    endif

#ifndef DRY
    !###########################################################################
    ! Update Tracers
    !###########################################################################

    dt = real(DTSEC,kind=RP)

!OCL XFILL
    mflx_hi(:,:,:,:) = mflx_av(:,:,:,:) / real(NSTEP_ATMOS_DYN,kind=RP)

    call COMM_vars8( mflx_hi(:,:,:,ZDIR), I_COMM_mflx_z )
    call COMM_vars8( mflx_hi(:,:,:,XDIR), I_COMM_mflx_x )
    call COMM_vars8( mflx_hi(:,:,:,YDIR), I_COMM_mflx_y )
    call COMM_wait ( mflx_hi(:,:,:,ZDIR), I_COMM_mflx_z, .false. )
    call COMM_wait ( mflx_hi(:,:,:,XDIR), I_COMM_mflx_x, .false. )
    call COMM_wait ( mflx_hi(:,:,:,YDIR), I_COMM_mflx_y, .false. )

    if ( USE_AVERAGE ) then
!OCL XFILL
       QTRC_av(:,:,:,:) = 0.0_RP
    endif

    !------------------------------------------------------------------------
    ! Update each tracer
    !------------------------------------------------------------------------
    do iq = 1, QA

       call PROF_rapstart("DYN_Numfilter", 2)

       if ( ND_COEF_Q == 0.0_RP ) then
!OCL XFILL
          num_diff_q(:,:,:,:) = 0.0_RP
       else
          call ATMOS_DYN_numfilter_coef_q( num_diff_q(:,:,:,:),                    & ! [OUT]
                                           DENS00, QTRC(:,:,:,iq),                 & ! [IN]
                                           CDZ, CDX, CDY, DT,                      & ! [IN]
                                           REF_qv, iq,                             & ! [IN]
                                           ND_COEF_Q, ND_ORDER, ND_SFC_FACT, ND_USE_RS ) ! [IN]
       endif

       call PROF_rapend  ("DYN_Numfilter", 2)

       call PROF_rapstart("DYN_Tracer", 2)

       do JJS = JS, JE, JBLOCK
       JJE = JJS+JBLOCK-1
       do IIS = IS, IE, IBLOCK
       IIE = IIS+IBLOCK-1

          !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
          do j = JJS-1, JJE+1
          do i = IIS-1, IIE+1
          do k = KS+1, KE-2
             qflx_lo(k,i,j,ZDIR) = 0.5_RP * (     mflx_hi(k,i,j,ZDIR)  * ( QTRC(k+1,i,j,iq)+QTRC(k,i,j,iq) ) &
                                            - abs(mflx_hi(k,i,j,ZDIR)) * ( QTRC(k+1,i,j,iq)-QTRC(k,i,j,iq) ) )

             qflx_hi(k,i,j,ZDIR) = mflx_hi(k,i,j,ZDIR) * ( FACT_N * ( QTRC(k+1,i,j,iq)+QTRC(k  ,i,j,iq) ) &
                                                         + FACT_F * ( QTRC(k+2,i,j,iq)+QTRC(k-1,i,j,iq) ) ) &
                                 + GSQRT(k,i,j,I_XYW) * num_diff_q(k,i,j,ZDIR)
          enddo
          enddo
          enddo

          !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
          do j = JJS-1, JJE+1
          do i = IIS-1, IIE+1
             qflx_lo(KS-1,i,j,ZDIR) = 0.0_RP
             qflx_lo(KS  ,i,j,ZDIR) = 0.5_RP * (     mflx_hi(KS  ,i,j,ZDIR)  * ( QTRC(KS+1,i,j,iq)+QTRC(KS,i,j,iq) ) &
                                               - abs(mflx_hi(KS  ,i,j,ZDIR)) * ( QTRC(KS+1,i,j,iq)-QTRC(KS,i,j,iq) ) )

             qflx_hi(KS-1,i,j,ZDIR) = 0.0_RP
             qflx_hi(KS  ,i,j,ZDIR) = 0.5_RP * mflx_hi(KS  ,i,j,ZDIR) * ( QTRC(KS+1,i,j,iq)+QTRC(KS,i,j,iq) ) &
                                    + GSQRT(k,i,j,I_XYW) * num_diff_q(KS,i,j,ZDIR)
          enddo
          enddo

          !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
          do j = JJS-1, JJE+1
          do i = IIS-1, IIE+1
             qflx_lo(KE-1,i,j,ZDIR) = 0.5_RP * (     mflx_hi(KE-1,i,j,ZDIR)  * ( QTRC(KE,i,j,iq)+QTRC(KE-1,i,j,iq) ) &
                                               - abs(mflx_hi(KE-1,i,j,ZDIR)) * ( QTRC(KE,i,j,iq)-QTRC(KE-1,i,j,iq) ) )
             qflx_lo(KE  ,i,j,ZDIR) = 0.0_RP

             qflx_hi(KE-1,i,j,ZDIR) = 0.5_RP * mflx_hi(KE-1,i,j,ZDIR) * ( QTRC(KE,i,j,iq)+QTRC(KE-1,i,j,iq) ) &
                                    + GSQRT(k,i,j,I_XYW) * num_diff_q(KE-1,i,j,ZDIR)
             qflx_hi(KE  ,i,j,ZDIR) = 0.0_RP
          enddo
          enddo

          !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
          do j = JJS-1, JJE+1
          do i = IIS-2, IIE+1
          do k = KS, KE
             qflx_lo(k,i,j,XDIR) = 0.5_RP * (     mflx_hi(k,i,j,XDIR)  * ( QTRC(k,i+1,j,iq)+QTRC(k,i,j,iq) ) &
                                            - abs(mflx_hi(k,i,j,XDIR)) * ( QTRC(k,i+1,j,iq)-QTRC(k,i,j,iq) ) )
          enddo
          enddo
          enddo

          !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
          do j = JJS,   JJE
          do i = IIS-1, IIE
          do k = KS, KE
             qflx_hi(k,i,j,XDIR) = mflx_hi(k,i,j,XDIR) * ( FACT_N * ( QTRC(k,i+1,j,iq)+QTRC(k,i  ,j,iq) ) &
                                                         + FACT_F * ( QTRC(k,i+2,j,iq)+QTRC(k,i-1,j,iq) ) ) &
                                 + GSQRT(k,i,j,I_UYZ) * num_diff_q(k,i,j,XDIR)
          enddo
          enddo
          enddo


          !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
          do j = JJS-2, JJE+1
          do i = IIS-1, IIE+1
          do k = KS, KE
             qflx_lo(k,i,j,YDIR) = 0.5_RP * (     mflx_hi(k,i,j,YDIR)  * ( QTRC(k,i,j+1,iq)+QTRC(k,i,j,iq) ) &
                                            - abs(mflx_hi(k,i,j,YDIR)) * ( QTRC(k,i,j+1,iq)-QTRC(k,i,j,iq) ) )
          enddo
          enddo
          enddo

          !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
          do j = JJS-1, JJE
          do i = IIS,   IIE
          do k = KS, KE
             qflx_hi(k,i,j,YDIR) = mflx_hi(k,i,j,YDIR) * ( FACT_N * ( QTRC(k,i,j+1,iq)+QTRC(k,i,j  ,iq) ) &
                                                         + FACT_F * ( QTRC(k,i,j+2,iq)+QTRC(k,i,j-1,iq) ) ) &
                                 + GSQRT(k,i,j,I_XVZ) * num_diff_q(k,i,j,YDIR)
          enddo
          enddo
          enddo

       enddo
       enddo

       call PROF_rapend ("DYN_Tracer", 2)

       call PROF_rapstart("DYN_Boundary", 2)

       if ( BND_W ) then
          if ( iq <= BND_QA ) then
             !$omp parallel do private(j,k) OMP_SCHEDULE_ collapse(2)
             do j = JS-1, JE+1
             do k = KS, KE
                qflx_hi(k,IS-1,j,XDIR) = mflx_hi(k,IS-1,j,XDIR) &
                                       * ( FACT_N * ( DAMP_QTRC(k,IS  ,j,iq)+DAMP_QTRC(k,IS-1,j,iq) ) &
                                         + FACT_F * ( DAMP_QTRC(k,IS+1,j,iq)+DAMP_QTRC(k,IS-2,j,iq) ) )
                qflx_lo(k,IS-1,j,XDIR) = 0.5_RP * ( &
                           mflx_hi(k,IS-1,j,XDIR)  * ( DAMP_QTRC(k,IS,j,iq)+DAMP_QTRC(k,IS-1,j,iq) ) &
                     - abs(mflx_hi(k,IS-1,j,XDIR)) * ( DAMP_QTRC(k,IS,j,iq)-DAMP_QTRC(k,IS-1,j,iq) ) )
             enddo
             enddo
          else
             !$omp parallel do private(j,k) OMP_SCHEDULE_ collapse(2)
             do j = JS-1, JE+1
             do k = KS, KE
                qflx_hi(k,IS-1,j,XDIR) = mflx_hi(k,IS-1,j,XDIR) &
                                       * ( FACT_N * ( QTRC(k,IS  ,j,iq)+QTRC(k,IS-1,j,iq) ) &
                                         + FACT_F * ( QTRC(k,IS+1,j,iq)+QTRC(k,IS-2,j,iq) ) )
                qflx_lo(k,IS-1,j,XDIR) = 0.5_RP * ( &
                           mflx_hi(k,IS-1,j,XDIR)  * ( QTRC(k,IS,j,iq)+QTRC(k,IS-1,j,iq) ) &
                     - abs(mflx_hi(k,IS-1,j,XDIR)) * ( QTRC(k,IS,j,iq)-QTRC(k,IS-1,j,iq) ) )
             enddo
             enddo
          end if
       end if
       if ( BND_E ) then
          if ( iq <= BND_QA ) then
             !$omp parallel do private(j,k) OMP_SCHEDULE_ collapse(2)
             do j = JS-1, JE+1
             do k = KS, KE
                qflx_hi(k,IE,j,XDIR) = mflx_hi(k,IE,j,XDIR) &
                                     * ( FACT_N * ( DAMP_QTRC(k,IE+1,j,iq)+DAMP_QTRC(k,IE  ,j,iq) ) &
                                       + FACT_F * ( DAMP_QTRC(k,IE+2,j,iq)+DAMP_QTRC(k,IE-1,j,iq) ) )
                qflx_lo(k,IE,j,XDIR) = 0.5_RP * ( &
                           mflx_hi(k,IE,j,XDIR)  * ( DAMP_QTRC(k,IE+1,j,iq)+DAMP_QTRC(k,IE,j,iq) ) &
                     - abs(mflx_hi(k,IE,j,XDIR)) * ( DAMP_QTRC(k,IE+1,j,iq)-DAMP_QTRC(k,IE,j,iq) ) )
             enddo
             enddo
          else
             !$omp parallel do private(j,k) OMP_SCHEDULE_ collapse(2)
             do j = JS-1, JE+1
             do k = KS, KE
                qflx_hi(k,IE,j,XDIR) = mflx_hi(k,IE,j,XDIR) &
                                     * ( FACT_N * ( QTRC(k,IE+1,j,iq)+QTRC(k,IE  ,j,iq) ) &
                                       + FACT_F * ( QTRC(k,IE+2,j,iq)+QTRC(k,IE-1,j,iq) ) )
                qflx_lo(k,IE,j,XDIR) = 0.5_RP * ( &
                           mflx_hi(k,IE,j,XDIR)  * ( QTRC(k,IE+1,j,iq)+QTRC(k,IE,j,iq) ) &
                     - abs(mflx_hi(k,IE,j,XDIR)) * ( QTRC(k,IE+1,j,iq)-QTRC(k,IE,j,iq) ) )
             enddo
             enddo
          end if
       end if
       if ( BND_S ) then
          if ( iq <= BND_QA ) then
             !$omp parallel do private(i,k) OMP_SCHEDULE_ collapse(2)
             do i = IS-1, IE+1
             do k = KS, KE
                qflx_hi(k,i,JS-1,YDIR) = mflx_hi(k,i,JS-1,YDIR) &
                                       * ( FACT_N * ( DAMP_QTRC(k,i,JS  ,iq)+DAMP_QTRC(k,i,JS-1,iq) ) &
                                         + FACT_F * ( DAMP_QTRC(k,i,JS+1,iq)+DAMP_QTRC(k,i,JS-2,iq) ) )
                qflx_lo(k,i,JS-1,YDIR) = 0.5_RP * ( &
                           mflx_hi(k,i,JS-1,YDIR)  * ( DAMP_QTRC(k,i,JS,iq)+DAMP_QTRC(k,i,JS-1,iq) ) &
                     - abs(mflx_hi(k,i,JS-1,YDIR)) * ( DAMP_QTRC(k,i,JS,iq)-DAMP_QTRC(k,i,JS-1,iq) ) )
             enddo
             enddo
          else
             !$omp parallel do private(i,k) OMP_SCHEDULE_ collapse(2)
             do i = IS-1, IE+1
             do k = KS, KE
                qflx_hi(k,i,JS-1,YDIR) = mflx_hi(k,i,JS-1,YDIR) &
                                       * ( FACT_N * ( QTRC(k,i,JS  ,iq)+QTRC(k,i,JS-1,iq) ) &
                                         + FACT_F * ( QTRC(k,i,JS+1,iq)+QTRC(k,i,JS-2,iq) ) )
                qflx_lo(k,i,JS-1,YDIR) = 0.5_RP * ( &
                           mflx_hi(k,i,JS-1,YDIR)  * ( QTRC(k,i,JS,iq)+QTRC(k,i,JS-1,iq) ) &
                     - abs(mflx_hi(k,i,JS-1,YDIR)) * ( QTRC(k,i,JS,iq)-QTRC(k,i,JS-1,iq) ) )
             enddo
             enddo
          end if
       end if
       if ( BND_N ) then
          if ( iq <= BND_QA ) then
             !$omp parallel do private(i,k) OMP_SCHEDULE_ collapse(2)
             do i = IS-1, IE+1
             do k = KS, KE
                qflx_hi(k,i,JE,YDIR) = mflx_hi(k,i,JE,YDIR) &
                                     * ( FACT_N * ( DAMP_QTRC(k,i,JE+1,iq)+DAMP_QTRC(k,i,JE  ,iq) ) &
                                       + FACT_F * ( DAMP_QTRC(k,i,JE+2,iq)+DAMP_QTRC(k,i,JE-1,iq) ) )
                qflx_lo(k,i,JE,YDIR) = 0.5_RP * ( &
                           mflx_hi(k,i,JE,YDIR)  * ( DAMP_QTRC(k,i,JE+1,iq)+DAMP_QTRC(k,i,JE,iq) ) &
                     - abs(mflx_hi(k,i,JE,YDIR)) * ( DAMP_QTRC(k,i,JE+1,iq)-DAMP_QTRC(k,i,JE,iq) ) )
             enddo
             enddo
          else
             !$omp parallel do private(i,k) OMP_SCHEDULE_ collapse(2)
             do i = IS-1, IE+1
             do k = KS, KE
                qflx_hi(k,i,JE,YDIR) = mflx_hi(k,i,JE,YDIR) &
                                     * ( FACT_N * ( QTRC(k,i,JE+1,iq)+QTRC(k,i,JE  ,iq) ) &
                                       + FACT_F * ( QTRC(k,i,JE+2,iq)+QTRC(k,i,JE-1,iq) ) )
                qflx_lo(k,i,JE,YDIR) = 0.5_RP * ( &
                           mflx_hi(k,i,JE,YDIR)  * ( QTRC(k,i,JE+1,iq)+QTRC(k,i,JE,iq) ) &
                     - abs(mflx_hi(k,i,JE,YDIR)) * ( QTRC(k,i,JE+1,iq)-QTRC(k,i,JE,iq) ) )
             enddo
             enddo
          end if
       end if

       call PROF_rapend  ("DYN_Boundary", 2)

       call PROF_rapstart("DYN_Tracer", 2)

       call ATMOS_DYN_fct( qflx_anti,                    & ! (out)
                           QTRC(:,:,:,iq), DENS00, DENS, & ! (in)
                           qflx_hi, qflx_lo,             & ! (in)
                           mflx_hi,                      & ! (in)
                           RCDZ, RCDX, RCDY,             & ! (in)
                           GSQRT(:,:,:,I_XYZ),           & ! (in)
                           MAPF(:,:,:,I_XY), dt,         & ! (in)
                           FLAG_FCT_ALONG_STREAM         ) ! (in)

       do JJS = JS, JE, JBLOCK
       JJE = JJS+JBLOCK-1
       do IIS = IS, IE, IBLOCK
       IIE = IIS+IBLOCK-1

          !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
          do j = JJS, JJE
          do i = IIS, IIE
          do k = KS, KE
             QTRC(k,i,j,iq) = ( QTRC(k,i,j,iq) * DENS00(k,i,j) &
                            + dt * ( - ( ( qflx_hi(k  ,i  ,j  ,ZDIR) - qflx_anti(k  ,i  ,j  ,ZDIR) &
                                         - qflx_hi(k-1,i  ,j  ,ZDIR) + qflx_anti(k-1,i  ,j  ,ZDIR) ) * RCDZ(k) &
                                       + ( qflx_hi(k  ,i  ,j  ,XDIR) - qflx_anti(k  ,i  ,j  ,XDIR) &
                                         - qflx_hi(k  ,i-1,j  ,XDIR) + qflx_anti(k  ,i-1,j  ,XDIR) ) * RCDX(i) &
                                       + ( qflx_hi(k  ,i  ,j  ,YDIR) - qflx_anti(k  ,i  ,j  ,YDIR) &
                                         - qflx_hi(k  ,i  ,j-1,YDIR) + qflx_anti(k  ,i  ,j-1,YDIR) ) * RCDY(j) &
                                       ) * MAPF(i,j,1,I_XY) * MAPF(i,j,2,I_XY) / GSQRT(k,i,j,I_XYZ) &
                              + RHOQ_t(k,i,j,iq) ) ) / DENS(k,i,j)
          enddo
          enddo
          enddo

       enddo
       enddo

       call PROF_rapend  ("DYN_Tracer", 2)

       do j  = JS, JE
       do i  = IS, IE
          QTRC(   1:KS-1,i,j,iq) = QTRC(KS,i,j,iq)
          QTRC(KE+1:KA,  i,j,iq) = QTRC(KE,i,j,iq)
       enddo
       enddo


       if ( USE_AVERAGE ) then
          QTRC_av(:,:,:,iq) = QTRC(:,:,:,iq)
       endif

    enddo ! scalar quantities loop

    do iq = 1, QA
       call COMM_vars8( QTRC(:,:,:,iq), I_COMM_QTRC(iq) )
    enddo
    do iq = 1, QA
       call COMM_wait ( QTRC(:,:,:,iq), I_COMM_QTRC(iq), .false. )
    enddo
#endif

    return
  end subroutine ATMOS_DYN

end module scale_atmos_dyn
