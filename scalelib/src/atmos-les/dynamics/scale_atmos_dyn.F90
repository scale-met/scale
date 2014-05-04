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
  real(RP), private, allocatable :: CNZ3(:,:,:)
  real(RP), private, allocatable :: CNX3(:,:,:)
  real(RP), private, allocatable :: CNY3(:,:,:)
  real(RP), private, allocatable :: CNZ4(:,:,:)
  real(RP), private, allocatable :: CNX4(:,:,:)
  real(RP), private, allocatable :: CNY4(:,:,:)

  real(RP), private, allocatable :: CORIOLI(:,:)            ! coriolis term

  !-----------------------------------------------------------------------------
contains

  !-----------------------------------------------------------------------------
  !> Setup
  subroutine ATMOS_DYN_setup( &
       DIFF4,            &
       DYN_TYPE,         &
       CDZ, CDX, CDY,    &
       FDZ, FDX, FDY,    &
       ND_ORDER,         &
       ND_COEF,          &
       DT,               &
       enable_coriolis,  &
       lat               )
    use scale_process, only: &
       PRC_MPIstop
    use scale_const, only: &
       OHM => CONST_OHM
    use scale_atmos_dyn_common, only: &
       ATMOS_DYN_numfilter_setup
    use scale_atmos_dyn_rk, only: &
       ATMOS_DYN_rk_setup
    implicit none

    real(RP),               intent(out) :: DIFF4
    character(len=H_SHORT), intent(in)  :: DYN_TYPE
    real(RP),               intent(in)  :: CDZ(KA)
    real(RP),               intent(in)  :: CDX(IA)
    real(RP),               intent(in)  :: CDY(JA)
    real(RP),               intent(in)  :: FDZ(KA-1)
    real(RP),               intent(in)  :: FDX(IA-1)
    real(RP),               intent(in)  :: FDY(JA-1)
    integer,                intent(in)  :: ND_ORDER
    real(RP),               intent(in)  :: ND_COEF
    real(RP),               intent(in)  :: DT
    logical ,               intent(in)  :: enable_coriolis
    real(RP),               intent(in)  :: lat(IA,JA)
    !---------------------------------------------------------------------------

    allocate( CNZ3(3,KA,2) )
    allocate( CNX3(3,IA,2) )
    allocate( CNY3(3,JA,2) )
    allocate( CNZ4(5,KA,2) )
    allocate( CNX4(5,IA,2) )
    allocate( CNY4(5,JA,2) )

    allocate( CORIOLI(IA,JA) )

    ! numerical diffusion
    call ATMOS_DYN_numfilter_setup( DIFF4,                              & ! [OUT]
                                    CNZ3, CNX3, CNY3, CNZ4, CNX4, CNY4, & ! [OUT]
                                    CDZ, CDX, CDY, FDZ, FDX, FDY,       & ! [IN]
                                    ND_ORDER, ND_COEF,                  & ! [IN]
                                    DT                                  ) ! [IN]

    ! coriolis parameter
    if ( enable_coriolis ) then
       CORIOLI(:,:) = 2.0_RP * OHM * sin( lat(:,:) )
    else
       CORIOLI(:,:) = 0.0_RP
    endif

    call ATMOS_DYN_rk_setup( DYN_TYPE )

    return
  end subroutine ATMOS_DYN_setup

  !-----------------------------------------------------------------------------
  !> Dynamical Process
  subroutine ATMOS_DYN( &
       DENS,    MOMZ,    MOMX,    MOMY,    RHOT,    QTRC,    &
       DENS_av, MOMZ_av, MOMX_av, MOMY_av, RHOT_av, QTRC_av, &
       DENS_tp, MOMZ_tp, MOMX_tp, MOMY_tp, RHOT_tp, QTRC_tp, &
       CDZ, CDX, CDY, FDZ, FDX, FDY,                         &
       RCDZ, RCDX, RCDY, RFDZ, RFDX, RFDY,                   &
       PHI, GSQRT,                                           &
       J13G, J23G, J33G,                                     &
       AQ_CV,                                                &
       REF_dens, REF_pott, REF_qv, REF_pres,                 &
       DIFF4, ND_ORDER, ND_SFC_FACT, ND_USE_RS,              &
       DAMP_var, DAMP_alpha,                                 &
       divdmp_coef,                                          &
       FLAG_FCT_RHO, FLAG_FCT_MOMENTUM, FLAG_FCT_T,          &
       USE_AVERAGE,                                          &
       DTSEC, DTSEC_ATMOS_DYN, NSTEP_ATMOS_DYN               )
    use scale_const, only: &
       Rdry   => CONST_Rdry,  &
       Rvap   => CONST_Rvap,  &
       CVdry  => CONST_CVdry
    use scale_comm, only: &
       COMM_vars8, &
       COMM_wait
    use scale_atmos_dyn_common, only: &
       FACT_N,                     &
       FACT_F,                     &
       ATMOS_DYN_numfilter_coef,   &
       ATMOS_DYN_numfilter_coef_q, &
       ATMOS_DYN_fct
    use scale_atmos_dyn_rk, only: &
       ATMOS_DYN_rk
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
    real(RP), intent(in)    :: QTRC_tp(KA,IA,JA,QA)

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
    real(RP), intent(in)    :: J13G (KA,IA,JA,4) !< (1,3) element of Jacobian matrix
    real(RP), intent(in)    :: J23G (KA,IA,JA,4) !< (2,3) element of Jacobian matrix
    real(RP), intent(in)    :: J33G              !< (3,3) element of Jacobian matrix

    real(RP), intent(in)    :: AQ_CV(QQA)

    real(RP), intent(in)    :: REF_dens(KA,IA,JA)
    real(RP), intent(in)    :: REF_pott(KA,IA,JA)
    real(RP), intent(in)    :: REF_qv  (KA,IA,JA)
    real(RP), intent(in)    :: REF_pres(KA,IA,JA)   !< reference pressure
    real(RP), intent(in)    :: DIFF4
    integer,  intent(in)    :: ND_ORDER
    real(RP), intent(in)    :: ND_SFC_FACT
    logical,  intent(in)    :: ND_USE_RS
    real(RP), intent(in)    :: DAMP_var  (KA,IA,JA,5)
    real(RP), intent(in)    :: DAMP_alpha(KA,IA,JA,5)
    real(RP), intent(in)    :: divdmp_coef

    logical,  intent(in)    :: FLAG_FCT_RHO
    logical,  intent(in)    :: FLAG_FCT_MOMENTUM
    logical,  intent(in)    :: FLAG_FCT_T

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
    real(RP) :: DENS_RK1(KA,IA,JA) ! prognostic variables (+1/3 step)
    real(RP) :: MOMZ_RK1(KA,IA,JA) !
    real(RP) :: MOMX_RK1(KA,IA,JA) !
    real(RP) :: MOMY_RK1(KA,IA,JA) !
    real(RP) :: RHOT_RK1(KA,IA,JA) !
    real(RP) :: DENS_RK2(KA,IA,JA) ! prognostic variables (+2/3 step)
    real(RP) :: MOMZ_RK2(KA,IA,JA) !
    real(RP) :: MOMX_RK2(KA,IA,JA) !
    real(RP) :: MOMY_RK2(KA,IA,JA) !
    real(RP) :: RHOT_RK2(KA,IA,JA) !

    ! tendency
    real(RP) :: DENS_t(KA,IA,JA)
    real(RP) :: MOMZ_t(KA,IA,JA)
    real(RP) :: MOMX_t(KA,IA,JA)
    real(RP) :: MOMY_t(KA,IA,JA)
    real(RP) :: RHOT_t(KA,IA,JA)

    ! diagnostic variables
    real(RP) :: QDRY (KA,IA,JA) ! dry air
    real(RP) :: Rtot (KA,IA,JA) ! total R
    real(RP) :: CVtot(KA,IA,JA) ! total CV

    real(RP) :: num_diff  (KA,IA,JA,5,3)
    real(RP) :: num_diff_q(KA,IA,JA,3)

    ! For tracer advection
    real(RP) :: RHOQ     (KA,IA,JA)    ! rho(previous) * phi(previous)
    real(RP) :: RHOQ_t   (KA,IA,JA)    ! tendency
    real(RP) :: mflx_hi(KA,IA,JA,3)  ! rho * vel(x,y,z) @ (u,v,w)-face high order
    real(RP) :: mflx_av(KA,IA,JA,3)  ! rho * vel(x,y,z) @ (u,v,w)-face average
    real(RP) :: qflx_hi  (KA,IA,JA,3)  ! rho * vel(x,y,z) * phi @ (u,v,w)-face high order
    real(RP) :: qflx_lo  (KA,IA,JA,3)  ! rho * vel(x,y,z) * phi,  monotone flux
    real(RP) :: qflx_anti(KA,IA,JA,3)  ! anti-diffusive flux

    real(RP) :: dt

    integer  :: IIS, IIE
    integer  :: JJS, JJE
    integer  :: i, j, k, iq, step
    !---------------------------------------------------------------------------

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

    RHOQ    (:,:,:) = UNDEF

    num_diff (:,:,:,:,:) = UNDEF

    mflx_hi(:,:,:,:) = UNDEF
    qflx_hi(:,:,:,:) = UNDEF
    qflx_lo(:,:,:,:) = UNDEF
#endif

    DENS00(:,:,:) = DENS(:,:,:)

    if ( USE_AVERAGE ) then
       DENS_av(:,:,:) = 0.0_RP
       MOMZ_av(:,:,:) = 0.0_RP
       MOMX_av(:,:,:) = 0.0_RP
       MOMY_av(:,:,:) = 0.0_RP
       RHOT_av(:,:,:) = 0.0_RP
    endif

#ifndef DRY
    mflx_av(:,:,:,:) = 0.0_RP

    CVtot(:,:,:) = 0.0_RP
    QDRY (:,:,:) = 1.0_RP
    do iq = QQS, QQE
       CVtot(:,:,:) = CVtot(:,:,:) + AQ_CV(iq) * QTRC(:,:,:,iq)
       QDRY (:,:,:) = QDRY (:,:,:) - QTRC(:,:,:,iq)
    enddo
    CVtot(:,:,:) = CVdry * QDRY(:,:,:) + CVtot(:,:,:)
    Rtot (:,:,:) = Rdry  * QDRY(:,:,:) + Rvap * QTRC(:,:,:,I_QV)
#endif

    !###########################################################################
    ! Update DENS,MONZ,MOMX,MOMY,MOMZ,RHOT
    !###########################################################################

    do step = 1, NSTEP_ATMOS_DYN

       !-----< prepare tendency >-----

       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j = JS, JE
       do i = IS, IE
          do k = KS, KE
             DENS_t(k,i,j) = DENS_tp(k,i,j) ! tendency from physical step
          enddo

          do k = KS, KE-1
             MOMZ_t(k,i,j) = MOMZ_tp(k,i,j) & ! tendency from physical step
                           - DAMP_alpha(k,i,j,I_BND_VELZ) &
                           * ( MOMZ(k,i,j) - DAMP_var(k,i,j,I_BND_VELZ) * 0.5_RP * ( DENS(k+1,i,j)+DENS(k,i,j) ) ) ! rayleigh damping
          enddo
          MOMZ_t(KE,i,j) = 0.0_RP

          do k = KS, KE
             MOMX_t(k,i,j) = MOMX_tp(k,i,j) & ! tendency from physical step
                           - DAMP_alpha(k,i,j,I_BND_VELX) &
                           * ( MOMX(k,i,j) - DAMP_var(k,i,j,I_BND_VELX) * 0.5_RP * ( DENS(k,i+1,j)+DENS(k,i,j) ) ) ! rayleigh damping
          enddo

          do k = KS, KE
             MOMY_t(k,i,j) = MOMY_tp(k,i,j) & ! tendency from physical step
                           - DAMP_alpha(k,i,j,I_BND_VELY) &
                           * ( MOMY(k,i,j) - DAMP_var(k,i,j,I_BND_VELY) * 0.5_RP * ( DENS(k,i,j+1)+DENS(k,i,j) ) ) ! rayleigh damping
          enddo

          do k = KS, KE
             RHOT_t(k,i,j) = RHOT_tp(k,i,j) & ! tendency from physical step
                           - DAMP_alpha(k,i,j,I_BND_POTT) &
                           * ( RHOT(k,i,j) - DAMP_var(k,i,j,I_BND_POTT) * DENS(k,i,j) )                            ! rayleigh damping
          enddo

          DENS_t(   1:KS-1,i,j) = 0.0_RP
          MOMZ_t(   1:KS-1,i,j) = 0.0_RP
          MOMX_t(   1:KS-1,i,j) = 0.0_RP
          MOMY_t(   1:KS-1,i,j) = 0.0_RP
          RHOT_t(   1:KS-1,i,j) = 0.0_RP
          DENS_t(KE+1:KA  ,i,j) = 0.0_RP
          MOMZ_t(KE+1:KA  ,i,j) = 0.0_RP
          MOMX_t(KE+1:KA  ,i,j) = 0.0_RP
          MOMY_t(KE+1:KA  ,i,j) = 0.0_RP
          RHOT_t(KE+1:KA  ,i,j) = 0.0_RP
       enddo
       enddo

       call COMM_vars8( DENS_t(:,:,:), 1 )
       call COMM_vars8( MOMZ_t(:,:,:), 2 )
       call COMM_vars8( MOMX_t(:,:,:), 3 )
       call COMM_vars8( MOMY_t(:,:,:), 4 )
       call COMM_vars8( RHOT_t(:,:,:), 5 )
       call COMM_wait ( DENS_t(:,:,:), 1 )
       call COMM_wait ( MOMZ_t(:,:,:), 2 )
       call COMM_wait ( MOMX_t(:,:,:), 3 )
       call COMM_wait ( MOMY_t(:,:,:), 4 )
       call COMM_wait ( RHOT_t(:,:,:), 5 )

       !-----< prepare numerical diffusion coefficient >-----

       if ( DIFF4 == 0.0_RP ) then
          num_diff(:,:,:,:,:) = 0.0_RP
       else
          call ATMOS_DYN_numfilter_coef( num_diff(:,:,:,:,:),                    & ! [OUT]
                                         DENS, MOMZ, MOMX, MOMY, RHOT,           & ! [IN]
                                         CNZ3, CNX3, CNY3, CNZ4, CNX4, CNY4,     & ! [IN]
                                         CDZ, CDX, CDY, FDZ, FDX, FDY,           & ! [IN]
                                         REF_dens, REF_pott,                     & ! [IN]
                                         DIFF4, ND_ORDER, ND_SFC_FACT, ND_USE_RS ) ! [IN]
       endif

       !------------------------------------------------------------------------
       ! Start RK
       !------------------------------------------------------------------------

       !##### SAVE #####
       DENS0(:,:,:) = DENS(:,:,:)
       MOMZ0(:,:,:) = MOMZ(:,:,:)
       MOMX0(:,:,:) = MOMX(:,:,:)
       MOMY0(:,:,:) = MOMY(:,:,:)
       RHOT0(:,:,:) = RHOT(:,:,:)

       !##### RK1 : PROG0,PROG->PROG_RK1 #####

       dt = real(DTSEC_ATMOS_DYN,kind=RP) / 3.0_RP

       call ATMOS_DYN_rk( DENS_RK1, MOMZ_RK1, MOMX_RK1, MOMY_RK1, RHOT_RK1, & ! (out)
                          mflx_hi,                                          & ! (out)
                          DENS0,    MOMZ0,    MOMX0,    MOMY0,    RHOT0,    & ! (in)
                          DENS,     MOMZ,     MOMX,     MOMY,     RHOT,     & ! (in)
                          DENS_t,   MOMZ_t,   MOMX_t,   MOMY_t,   RHOT_t,   & ! (in)
                          Rtot, CVtot, CORIOLI,                             & ! (in)
                          num_diff, divdmp_coef,                            & ! (in)
                          FLAG_FCT_RHO, FLAG_FCT_MOMENTUM, FLAG_FCT_T,      & ! (in)
                          CDZ, FDZ, FDX, FDY,                               & ! (in)
                          RCDZ, RCDX, RCDY, RFDZ, RFDX, RFDY,               & ! (in)
                          PHI, GSQRT, J13G, J23G, J33G,                     & ! (in)
                          REF_pres, REF_dens,                               & ! (in)
                          dt                                                ) ! (in)

       call COMM_vars8( DENS_RK1(:,:,:), 1 )
       call COMM_vars8( MOMZ_RK1(:,:,:), 2 )
       call COMM_vars8( MOMX_RK1(:,:,:), 3 )
       call COMM_vars8( MOMY_RK1(:,:,:), 4 )
       call COMM_vars8( RHOT_RK1(:,:,:), 5 )
       call COMM_wait ( DENS_RK1(:,:,:), 1 )
       call COMM_wait ( MOMZ_RK1(:,:,:), 2 )
       call COMM_wait ( MOMX_RK1(:,:,:), 3 )
       call COMM_wait ( MOMY_RK1(:,:,:), 4 )
       call COMM_wait ( RHOT_RK1(:,:,:), 5 )

       !##### RK2 : PROG0,PROG_RK2->PROG_RK3 #####

       dt = real(DTSEC_ATMOS_DYN,kind=RP) / 2.0_RP

       call ATMOS_DYN_rk( DENS_RK2, MOMZ_RK2, MOMX_RK2, MOMY_RK2, RHOT_RK2, & ! (out)
                          mflx_hi,                                          & ! (out)
                          DENS0,    MOMZ0,    MOMX0,    MOMY0,    RHOT0,    & ! (in)
                          DENS_RK1, MOMZ_RK1, MOMX_RK1, MOMY_RK1, RHOT_RK1, & ! (in)
                          DENS_t,   MOMZ_t,   MOMX_t,   MOMY_t,   RHOT_t,   & ! (in)
                          Rtot, CVtot, CORIOLI,                             & ! (in)
                          num_diff, divdmp_coef,                            & ! (in)
                          FLAG_FCT_RHO, FLAG_FCT_MOMENTUM, FLAG_FCT_T,      & ! (in)
                          CDZ, FDZ, FDX, FDY,                               & ! (in)
                          RCDZ, RCDX, RCDY, RFDZ, RFDX, RFDY,               & ! (in)
                          PHI, GSQRT, J13G, J23G, J33G,                     & ! (in)
                          REF_pres, REF_dens,                               & ! (in)
                          dt                                                ) ! (in)

       call COMM_vars8( DENS_RK2(:,:,:), 1 )
       call COMM_vars8( MOMZ_RK2(:,:,:), 2 )
       call COMM_vars8( MOMX_RK2(:,:,:), 3 )
       call COMM_vars8( MOMY_RK2(:,:,:), 4 )
       call COMM_vars8( RHOT_RK2(:,:,:), 5 )
       call COMM_wait ( DENS_RK2(:,:,:), 1 )
       call COMM_wait ( MOMZ_RK2(:,:,:), 2 )
       call COMM_wait ( MOMX_RK2(:,:,:), 3 )
       call COMM_wait ( MOMY_RK2(:,:,:), 4 )
       call COMM_wait ( RHOT_RK2(:,:,:), 5 )

       !##### RK3 : PROG0,PROG_RK3->PROG #####

       dt = real(DTSEC_ATMOS_DYN,kind=RP)

       call ATMOS_DYN_rk( DENS,     MOMZ,     MOMX,     MOMY,     RHOT,     & ! (out)
                          mflx_hi,                                          & ! (out)
                          DENS0,    MOMZ0,    MOMX0,    MOMY0,    RHOT0,    & ! (in)
                          DENS_RK2, MOMZ_RK2, MOMX_RK2, MOMY_RK2, RHOT_RK2, & ! (in)
                          DENS_t,   MOMZ_t,   MOMX_t,   MOMY_t,   RHOT_t,   & ! (in)
                          Rtot, CVtot, CORIOLI,                             & ! (in)
                          num_diff, divdmp_coef,                            & ! (in)
                          FLAG_FCT_RHO, FLAG_FCT_MOMENTUM, FLAG_FCT_T,      & ! (in)
                          CDZ, FDZ, FDX, FDY,                               & ! (in)
                          RCDZ, RCDX, RCDY, RFDZ, RFDX, RFDY,               & ! (in)
                          PHI, GSQRT, J13G, J23G, J33G,                     & ! (in)
                          REF_pres, REF_dens,                               & ! (in)
                          dt                                                ) ! (in)

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

       call COMM_vars8( DENS(:,:,:), 1 )
       call COMM_vars8( MOMZ(:,:,:), 2 )
       call COMM_vars8( MOMX(:,:,:), 3 )
       call COMM_vars8( MOMY(:,:,:), 4 )
       call COMM_vars8( RHOT(:,:,:), 5 )
       call COMM_wait ( DENS(:,:,:), 1 )
       call COMM_wait ( MOMZ(:,:,:), 2 )
       call COMM_wait ( MOMX(:,:,:), 3 )
       call COMM_wait ( MOMY(:,:,:), 4 )
       call COMM_wait ( RHOT(:,:,:), 5 )

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

    mflx_hi(:,:,:,:) = mflx_av(:,:,:,:) / real(NSTEP_ATMOS_DYN,kind=RP)

    call COMM_vars8( mflx_hi(:,:,:,ZDIR), 1 )
    call COMM_vars8( mflx_hi(:,:,:,XDIR), 2 )
    call COMM_vars8( mflx_hi(:,:,:,YDIR), 3 )
    call COMM_wait ( mflx_hi(:,:,:,ZDIR), 1 )
    call COMM_wait ( mflx_hi(:,:,:,XDIR), 2 )
    call COMM_wait ( mflx_hi(:,:,:,YDIR), 3 )

    if ( USE_AVERAGE ) then
       QTRC_av(:,:,:,:) = 0.0_RP
    endif

    !------------------------------------------------------------------------
    ! Update each tracer
    !------------------------------------------------------------------------
#ifdef _SDM
    do iq = 1, I_QV
#else
    do iq = 1, QA
#endif

       if ( iq == I_QV ) then
          !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
          do j = JS, JE
          do i = IS, IE
          do k = KS, KE
             RHOQ_t(k,i,j) = QTRC_tp(k,i,j,I_QV) & ! tendency from physical step
                           - DAMP_alpha(k,i,j,I_BND_QV) &
                           * ( QTRC(k,i,j,iq) - DAMP_var(k,i,j,I_BND_QV) ) * DENS00(k,i,j) ! rayleigh damping
          enddo
          enddo
          enddo
       else
          !$omp parallel do private(i,j,k,iq) OMP_SCHEDULE_ collapse(2)
          do j = JS, JE
          do i = IS, IE
          do k = KS, KE
             RHOQ_t(k,i,j) = QTRC_tp(k,i,j,iq) * DENS00(k,i,j)
          enddo
          enddo
          enddo
       endif
       do j = JS, JE
       do i = IS, IE
          RHOQ_t(   1:KS-1,i,j) = 0.0_RP
          RHOQ_t(KE+1:KA  ,i,j) = 0.0_RP
       enddo
       enddo

       call COMM_vars8( RHOQ_t(:,:,:), 5+iq )
       call COMM_wait ( RHOQ_t(:,:,:), 5+iq )

       if ( DIFF4 == 0.0_RP ) then
          num_diff_q(:,:,:,:) = 0.0_RP
       else
          call ATMOS_DYN_numfilter_coef_q( num_diff_q(:,:,:,:),                    & ! [OUT]
                                           DENS, QTRC(:,:,:,iq),                   & ! [IN]
                                           CNZ3, CNX3, CNY3, CNZ4, CNX4, CNY4,     & ! [IN]
                                           CDZ, CDX, CDY,                          & ! [IN]
                                           REF_qv, iq,                             & ! [IN]
                                           DIFF4, ND_ORDER, ND_SFC_FACT, ND_USE_RS ) ! [IN]
       endif

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
                                 + num_diff_q(k,i,j,ZDIR)
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
                                    + num_diff_q(KS,i,j,ZDIR)
          enddo
          enddo

          !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
          do j = JJS-1, JJE+1
          do i = IIS-1, IIE+1
             qflx_lo(KE-1,i,j,ZDIR) = 0.5_RP * (     mflx_hi(KE-1,i,j,ZDIR)  * ( QTRC(KE,i,j,iq)+QTRC(KE-1,i,j,iq) ) &
                                               - abs(mflx_hi(KE-1,i,j,ZDIR)) * ( QTRC(KE,i,j,iq)-QTRC(KE-1,i,j,iq) ) )
             qflx_lo(KE  ,i,j,ZDIR) = 0.0_RP

             qflx_hi(KE-1,i,j,ZDIR) = 0.5_RP * mflx_hi(KE-1,i,j,ZDIR) * ( QTRC(KE,i,j,iq)+QTRC(KE-1,i,j,iq) ) &
                                    + num_diff_q(KE-1,i,j,ZDIR)
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
                                 + num_diff_q(k,i,j,XDIR)
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
                                 + num_diff_q(k,i,j,YDIR)
          enddo
          enddo
          enddo

       enddo
       enddo

       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j = JS-1, JE+1
       do i = IS-1, IE+1
       do k = KS, KE
          RHOQ(k,i,j) = QTRC(k,i,j,iq) * DENS00(k,i,j)
       enddo
       enddo
       enddo

       call ATMOS_DYN_fct( qflx_anti,              & ! (out)
                           RHOQ, qflx_hi, qflx_lo, & ! (in)
                           RCDZ, RCDX, RCDY, dt    ) ! (in)

       do JJS = JS, JE, JBLOCK
       JJE = JJS+JBLOCK-1
       do IIS = IS, IE, IBLOCK
       IIE = IIS+IBLOCK-1

          !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
          do j = JJS, JJE
          do i = IIS, IIE
          do k = KS, KE
             QTRC(k,i,j,iq) = ( RHOQ(k,i,j) &
                              + DTSEC * ( - ( ( qflx_hi(k  ,i  ,j  ,ZDIR) + qflx_anti(k  ,i  ,j  ,ZDIR) &
                                              - qflx_hi(k-1,i  ,j  ,ZDIR) - qflx_anti(k-1,i  ,j  ,ZDIR) ) * RCDZ(k) &
                                            + ( qflx_hi(k  ,i  ,j  ,XDIR) + qflx_anti(k  ,i  ,j  ,XDIR) &
                                              - qflx_hi(k  ,i-1,j  ,XDIR) - qflx_anti(k  ,i-1,j  ,XDIR) ) * RCDX(i) &
                                            + ( qflx_hi(k  ,i  ,j  ,YDIR) + qflx_anti(k  ,i  ,j  ,YDIR) &
                                              - qflx_hi(k  ,i  ,j-1,YDIR) - qflx_anti(k  ,i  ,j-1,YDIR) ) * RCDY(j) ) &
                                          + RHOQ_t(k,i,j)                                                           ) &
                              ) / DENS(k,i,j)
          enddo
          enddo
          enddo

       enddo
       enddo

#ifdef DEBUG
       qflx_lo  (:,:,:,:) = UNDEF
       qflx_hi  (:,:,:,:) = UNDEF
       qflx_anti(:,:,:,:) = UNDEF
#endif

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
       call COMM_vars8( QTRC(:,:,:,iq), iq )
    enddo
    do iq = 1, QA
       call COMM_wait ( QTRC(:,:,:,iq), iq )
    enddo
#endif

    return
  end subroutine ATMOS_DYN

end module scale_atmos_dyn
