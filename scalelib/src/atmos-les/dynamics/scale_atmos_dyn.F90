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
  real(RP), private, allocatable :: CORIOLI(:,:)            ! coriolis term
  real(RP), private, allocatable :: mflx_hi(:,:,:,:)        ! rho * vel(x,y,z) @ (u,v,w)-face high order

  real(RP), private, allocatable :: num_diff  (:,:,:,:,:)
  real(RP), private, allocatable :: num_diff_q(:,:,:,:)

  logical, private :: BND_W
  logical, private :: BND_E
  logical, private :: BND_S
  logical, private :: BND_N

  logical, private :: DYN_NONE = .false.

  ! for communication
  integer :: I_COMM_DENS = 1
  integer :: I_COMM_MOMZ = 2
  integer :: I_COMM_MOMX = 3
  integer :: I_COMM_MOMY = 4
  integer :: I_COMM_RHOT = 5
  integer, allocatable :: I_COMM_PROG(:)
  integer, allocatable :: I_COMM_QTRC(:)
  !-----------------------------------------------------------------------------
contains

  !-----------------------------------------------------------------------------
  !> Setup
  subroutine ATMOS_DYN_setup( &
       DYN_Tinteg_Short_TYPE, &
       DYN_Tinteg_Large_TYPE, &
       DYN_FVM_FLUX_SCHEME, &
       DENS, MOMZ, MOMX, MOMY, RHOT, QTRC, &
       PROG,             &
       CDZ, CDX, CDY,    &
       FDZ, FDX, FDY,    &
       enable_coriolis,  &
       lat,              &
       none              )
    use scale_process, only: &
       PRC_MPIstop
    use scale_les_process, only: &
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
    use scale_atmos_dyn_tinteg_short, only: &
       ATMOS_DYN_tinteg_short_setup
    use scale_atmos_dyn_tinteg_large, only: &
       ATMOS_DYN_tinteg_large_setup
    use scale_atmos_dyn_tstep_short, only: &
       ATMOS_DYN_tstep_short_setup
    use scale_atmos_dyn_tstep_large, only: &
       ATMOS_DYN_tstep_large_setup
    use scale_atmos_dyn_fvm_flux, only: &
       ATMOS_DYN_FVM_flux_setup
    implicit none

    character(len=*), intent(in) :: DYN_Tinteg_Short_TYPE
    character(len=*), intent(in) :: DYN_Tinteg_Large_TYPE
    character(len=*), intent(in) :: DYN_FVM_FLUX_SCHEME

    ! MPI_RECV_INIT requires intent(inout)
    real(RP),               intent(inout) :: DENS(KA,IA,JA)
    real(RP),               intent(inout) :: MOMZ(KA,IA,JA)
    real(RP),               intent(inout) :: MOMX(KA,IA,JA)
    real(RP),               intent(inout) :: MOMY(KA,IA,JA)
    real(RP),               intent(inout) :: RHOT(KA,IA,JA)
    real(RP),               intent(inout) :: QTRC(KA,IA,JA,QA)

    real(RP),               intent(inout) :: PROG(KA,IA,JA,VA)

    real(RP),               intent(in) :: CDZ(KA)
    real(RP),               intent(in) :: CDX(IA)
    real(RP),               intent(in) :: CDY(JA)
    real(RP),               intent(in) :: FDZ(KA-1)
    real(RP),               intent(in) :: FDX(IA-1)
    real(RP),               intent(in) :: FDY(JA-1)
    logical,                intent(in) :: enable_coriolis
    real(RP),               intent(in) :: lat(IA,JA)
    logical, optional,      intent(in) :: none

    integer :: iv, iq
    !---------------------------------------------------------------------------

    allocate( CORIOLI(IA,JA) )
    allocate( mflx_hi(KA,IA,JA,3) )

    allocate( num_diff  (KA,IA,JA,5,3) )
    allocate( num_diff_q(KA,IA,JA,3) )

    allocate( I_COMM_PROG(max(VA,1)) )
    allocate( I_COMM_QTRC(QA) )

    call ATMOS_DYN_FVM_flux_setup( DYN_FVM_FLUX_SCHEME )

    call ATMOS_DYN_tstep_short_setup

    call ATMOS_DYN_tstep_large_setup( &
         DENS, MOMZ, MOMX, MOMY, RHOT, QTRC, PROG, &
         mflx_hi )

    call ATMOS_DYN_Tinteg_short_setup( DYN_Tinteg_Short_TYPE )

    call ATMOS_DYN_Tinteg_large_setup( DYN_Tinteg_Large_TYPE )

    if ( present(none) ) then
       DYN_NONE = none
    end if

    if ( .not. DYN_NONE ) then

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

    else

       call COMM_vars8_init( DENS, I_COMM_DENS )
       call COMM_vars8_init( MOMZ, I_COMM_MOMZ )
       call COMM_vars8_init( MOMX, I_COMM_MOMX )
       call COMM_vars8_init( MOMY, I_COMM_MOMY )
       call COMM_vars8_init( RHOT, I_COMM_RHOT )
       do iv = 1, VA
          I_COMM_PROG(iv) = 5 + iv
          call COMM_vars8_init( PROG(:,:,:,iv), I_COMM_PROG(iv) )
       end do

       do iq = 1, QA
          I_COMM_QTRC(iq) = 5 + VA + iq
          call COMM_vars8_init( QTRC  (:,:,:,iq), I_COMM_QTRC(iq) )
       end do

    end if

    BND_W = .not. PRC_HAS_W
    BND_E = .not. PRC_HAS_E
    BND_S = .not. PRC_HAS_S
    BND_N = .not. PRC_HAS_N


    mflx_hi(:,:,:,:) = UNDEF

    return
  end subroutine ATMOS_DYN_setup

  !-----------------------------------------------------------------------------
  !> Dynamical Process
  subroutine ATMOS_DYN( &
       DENS,    MOMZ,    MOMX,    MOMY,    RHOT,    QTRC,    &
       PROG,                                                 &
       DENS_av, MOMZ_av, MOMX_av, MOMY_av, RHOT_av, QTRC_av, &
       DENS_tp, MOMZ_tp, MOMX_tp, MOMY_tp, RHOT_tp, RHOQ_tp, &
       CDZ, CDX, CDY, FDZ, FDX, FDY,                         &
       RCDZ, RCDX, RCDY, RFDZ, RFDX, RFDY,                   &
       PHI, GSQRT,                                           &
       J13G, J23G, J33G, MAPF,                               &
       AQ_CV,                                                &
       REF_dens, REF_pott, REF_qv, REF_pres,                 &
       ND_COEF, ND_COEF_Q, ND_ORDER, ND_SFC_FACT, ND_USE_RS, &
       DAMP_DENS,       DAMP_VELZ,       DAMP_VELX,          &
       DAMP_VELY,       DAMP_POTT,       DAMP_QTRC,          &
       DAMP_alpha_DENS, DAMP_alpha_VELZ, DAMP_alpha_VELX,    &
       DAMP_alpha_VELY, DAMP_alpha_POTT, DAMP_alpha_QTRC,    &
       divdmp_coef,                                          &
       FLAG_FCT_MOMENTUM, FLAG_FCT_T, FLAG_FCT_TRACER,       &
       FLAG_FCT_ALONG_STREAM,                                &
       USE_AVERAGE,                                          &
       DTSEC, NSTEP_ATMOS_DYN                                )
    use scale_comm, only: &
       COMM_vars8, &
       COMM_wait
    use scale_atmos_boundary, only: &
       BND_QA
    use scale_atmos_dyn_tinteg_large, only: &
       ATMOS_DYN_tinteg_large
    implicit none

    real(RP), intent(inout) :: DENS(KA,IA,JA)
    real(RP), intent(inout) :: MOMZ(KA,IA,JA)
    real(RP), intent(inout) :: MOMX(KA,IA,JA)
    real(RP), intent(inout) :: MOMY(KA,IA,JA)
    real(RP), intent(inout) :: RHOT(KA,IA,JA)
    real(RP), intent(inout) :: QTRC(KA,IA,JA,QA)
    real(RP), intent(inout) :: PROG(KA,IA,JA,VA)

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

    logical,  intent(in)    :: FLAG_FCT_MOMENTUM
    logical,  intent(in)    :: FLAG_FCT_T
    logical,  intent(in)    :: FLAG_FCT_TRACER
    logical,  intent(in)    :: FLAG_FCT_ALONG_STREAM

    logical,  intent(in)    :: USE_AVERAGE

    real(DP), intent(in)    :: DTSEC
    integer , intent(in)    :: NSTEP_ATMOS_DYN

    ! for time integartion
    real(RP) :: DENS00  (KA,IA,JA) ! saved density before small step loop
    real(RP) :: DENS0   (KA,IA,JA) ! prognostic variables (+0   step)
    real(RP) :: MOMZ0   (KA,IA,JA) !
    real(RP) :: MOMX0   (KA,IA,JA) !
    real(RP) :: MOMY0   (KA,IA,JA) !
    real(RP) :: RHOT0   (KA,IA,JA) !
    real(RP) :: PROG0   (KA,IA,JA,VA) !

    ! diagnostic variables
    real(RP) :: QDRY (KA,IA,JA) ! dry air
    real(RP) :: Rtot (KA,IA,JA) ! total R
    real(RP) :: CVtot(KA,IA,JA) ! total CV
    real(RP) :: DDIV (KA,IA,JA) ! 3 dimensional divergence

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

    real(RP) :: dt

    integer  :: IIS, IIE
    integer  :: JJS, JJE
    integer  :: i, j, k, iq, step
    integer  :: iv

    real(RP) :: diff_coef
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*) '*** Dynamics step'


    if ( DYN_NONE ) then

       call PROF_rapstart("DYN_Tinteg", 2)

       dt = real(DTSEC, kind=RP)

       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          DENS(k,i,j) = DENS(k,i,j) + DENS_tp(k,i,j) * dt
       end do
       end do
       end do

       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          MOMZ(k,i,j) = MOMZ(k,i,j) + MOMZ_tp(k,i,j) * dt
       end do
       end do
       end do

       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          MOMX(k,i,j) = MOMX(k,i,j) + MOMX_tp(k,i,j) * dt
       end do
       end do
       end do

       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          MOMY(k,i,j) = MOMY(k,i,j) + MOMY_tp(k,i,j) * dt
       end do
       end do
       end do

       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          RHOT(k,i,j) = RHOT(k,i,j) + RHOT_tp(k,i,j) * dt
       end do
       end do
       end do

       do iq = 1, QA
       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          QTRC(k,i,j,iq) = max( 0.0_RP, &
                                QTRC(k,i,j,iq) * DENS00(k,i,j) + RHOQ_tp(k,i,j,iq) * dt &
                              ) / DENS(k,i,j)
       end do
       end do
       end do
       end do

       call COMM_vars8( DENS(:,:,:), I_COMM_DENS )
       call COMM_vars8( MOMZ(:,:,:), I_COMM_MOMZ )
       call COMM_vars8( MOMX(:,:,:), I_COMM_MOMX )
       call COMM_vars8( MOMY(:,:,:), I_COMM_MOMY )
       call COMM_vars8( RHOT(:,:,:), I_COMM_RHOT )
       do iq = 1, QA
          call COMM_vars8( QTRC(:,:,:,iq), I_COMM_QTRC(iq) )
       enddo
       call COMM_wait ( DENS(:,:,:), I_COMM_DENS, .false. )
       call COMM_wait ( MOMZ(:,:,:), I_COMM_MOMZ, .false. )
       call COMM_wait ( MOMX(:,:,:), I_COMM_MOMX, .false. )
       call COMM_wait ( MOMY(:,:,:), I_COMM_MOMY, .false. )
       call COMM_wait ( RHOT(:,:,:), I_COMM_RHOT, .false. )
       do iq = 1, QA
          call COMM_wait ( QTRC(:,:,:,iq), I_COMM_QTRC(iq), .false. )
       enddo

       if ( USE_AVERAGE ) then
          DENS_av(:,:,:) = DENS(:,:,:)
          MOMZ_av(:,:,:) = MOMZ(:,:,:)
          MOMX_av(:,:,:) = MOMX(:,:,:)
          MOMY_av(:,:,:) = MOMY(:,:,:)
          RHOT_av(:,:,:) = RHOT(:,:,:)
          QTRC_av(:,:,:,:) = QTRC(:,:,:,:)
       end if

       call PROF_rapend("DYN_Tinteg", 2)

       return
    end if ! if DYN_NONE == .true.

    call PROF_rapstart("DYN_Tinteg", 2)

    call ATMOS_DYN_tinteg_large( &
       DENS,    MOMZ,    MOMX,    MOMY,    RHOT,    QTRC,    & ! (inout)
       PROG,                                                 & ! (inout)
       DENS_av, MOMZ_av, MOMX_av, MOMY_av, RHOT_av, QTRC_av, & ! (inout)
       mflx_hi, tflx_hi,                                     & ! (out)
       num_diff, num_diff_q,                                 & ! (out; work)
       DENS_tp, MOMZ_tp, MOMX_tp, MOMY_tp, RHOT_tp, RHOQ_tp, & ! (in)
       CORIOLI,                                              & ! (in)
       CDZ, CDX, CDY, FDZ, FDX, FDY,                         & ! (in)
       RCDZ, RCDX, RCDY, RFDZ, RFDX, RFDY,                   & ! (in)
       PHI, GSQRT,                                           & ! (in)
       J13G, J23G, J33G, MAPF,                               & ! (in)
       AQ_CV,                                                & ! (in)
       REF_dens, REF_pott, REF_qv, REF_pres,                 & ! (in)
       BND_W, BND_E, BND_S, BND_N,                           & ! (in)
       ND_COEF, ND_COEF_Q, ND_ORDER, ND_SFC_FACT, ND_USE_RS, & ! (in)
       DAMP_DENS,       DAMP_VELZ,       DAMP_VELX,          & ! (in)
       DAMP_VELY,       DAMP_POTT,       DAMP_QTRC,          & ! (in)
       DAMP_alpha_DENS, DAMP_alpha_VELZ, DAMP_alpha_VELX,    & ! (in)
       DAMP_alpha_VELY, DAMP_alpha_POTT, DAMP_alpha_QTRC,    & ! (in)
       divdmp_coef,                                          & ! (in)
       FLAG_FCT_MOMENTUM, FLAG_FCT_T, FLAG_FCT_TRACER,       & ! (in)
       FLAG_FCT_ALONG_STREAM,                                & ! (in)
       USE_AVERAGE,                                          & ! (in)
       DTSEC, NSTEP_ATMOS_DYN                                ) ! (in)

    call PROF_rapend  ("DYN_Tinteg", 2)


#ifdef CHECK_MASS
       call check_mass( &
            DENS, DAMP_DENS, &
            mflx_hi, tflx_hi, &
            GSQRT, MAPF, &
            RCDX, RCDY, &
            dt, &
            BND_W, BND_E, BND_S, BND_N )
#endif

    return
  end subroutine ATMOS_DYN

#ifdef CHECK_MASS
  subroutine check_mass( &
       DENS, DAMP_DENS, &
       mflx_hi, tflx_hi, &
       GSQRT, MAPF, &
       RCDX, RCDY, &
       dt, &
       BND_W, BND_E, BND_S, BND_N &
       )
    use mpi
    use scale_grid_real, only: &
       vol => REAL_VOL
    use scale_comm, only: &
       COMM_datatype, &
       COMM_world
    use scale_history, only: &
       HIST_in
    use scale_gridtrans, only: &
       I_XYZ, &
       I_XY
    implicit none
    real(RP), intent(in) :: DENS     (KA,IA,JA)
    real(RP), intent(in) :: DAMP_DENS(KA,IA,JA)
    real(RP), intent(in) :: mflx_hi  (KA,IA,JA,3)
    real(RP), intent(in) :: tflx_hi  (KA,IA,JA,3)
    real(RP), intent(in) :: GSQRT    (KA,IA,JA,7)
    real(RP), intent(in) :: MAPF     (   IA,JA,2,7)
    real(RP), intent(in) :: RCDX(IA)
    real(RP), intent(in) :: RCDY(JA)
    real(RP), intent(in) :: dt
    logical,  intent(in) :: BND_W
    logical,  intent(in) :: BND_E
    logical,  intent(in) :: BND_S
    logical,  intent(in) :: BND_N

    ! lateral boundary flux
    real(RP) :: mflx_lb_horizontal(KA)
    real(RP) :: allmflx_lb_horizontal(KA)
    real(RP) :: mflx_lb_total
    real(RP) :: mass_total
    real(RP) :: mass_total2
    real(RP) :: allmflx_lb_total
    real(RP) :: allmass_total
    real(RP) :: allmass_total2

    integer :: k, i, j
    integer :: ierr


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

    if( IO_L ) write(IO_FID_LOG,'(A,1x,ES24.17)') 'total mflx_lb:', allmflx_lb_total

    call MPI_Allreduce( mass_total,           &
                        allmass_total,        &
                        1,                    &
                        COMM_datatype,        &
                        MPI_SUM,              &
                        COMM_world,           &
                        ierr                  )

    if( IO_L ) write(IO_FID_LOG,'(A,1x,ES24.17)') 'total mass   :', allmass_total

    call MPI_Allreduce( mass_total2,          &
                        allmass_total2,       &
                        1,                    &
                        COMM_datatype,        &
                        MPI_SUM,              &
                        COMM_world,           &
                        ierr                  )

    if( IO_L ) write(IO_FID_LOG,'(A,1x,ES24.17)') 'total mass2  :', allmass_total2

    call MPI_Allreduce( mflx_lb_horizontal(KS:KE),    &
                        allmflx_lb_horizontal(KS:KE), &
                        KMAX,                         &
                        COMM_datatype,                &
                        MPI_SUM,                      &
                        COMM_world,                   &
                        ierr                          )

    call HIST_in(allmflx_lb_horizontal(:), 'ALLMOM_lb_hz',                           &
                    'horizontally total momentum flux from lateral boundary', 'kg/m2/s' )

    return
  end subroutine check_mass
#endif

end module scale_atmos_dyn
