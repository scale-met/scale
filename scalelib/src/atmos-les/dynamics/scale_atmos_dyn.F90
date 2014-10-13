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
#ifdef DEBUG
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
  real(RP), private, allocatable :: CORIOLI(:,:)            ! coriolis term

  integer,  private              :: NUM_LBFLX               ! number of cell adjusting lateral boundary flux

  !-----------------------------------------------------------------------------
contains

  !-----------------------------------------------------------------------------
  !> Setup
  subroutine ATMOS_DYN_setup( &
       DYN_TYPE,         &
       CDZ, CDX, CDY,    &
       FDZ, FDX, FDY,    &
       enable_coriolis,  &
       adjust_flux_cell, &
       lat               )
    use scale_process, only: &
       PRC_MPIstop
    use scale_const, only: &
       OHM => CONST_OHM
    use scale_atmos_dyn_common, only: &
       ATMOS_DYN_filter_setup
    use scale_atmos_dyn_rk, only: &
       ATMOS_DYN_rk_setup
    implicit none

    character(len=H_SHORT), intent(in)  :: DYN_TYPE
    real(RP),               intent(in)  :: CDZ(KA)
    real(RP),               intent(in)  :: CDX(IA)
    real(RP),               intent(in)  :: CDY(JA)
    real(RP),               intent(in)  :: FDZ(KA-1)
    real(RP),               intent(in)  :: FDX(IA-1)
    real(RP),               intent(in)  :: FDY(JA-1)
    logical,                intent(in)  :: enable_coriolis
    integer,                intent(in)  :: adjust_flux_cell
    real(RP),               intent(in)  :: lat(IA,JA)
    !---------------------------------------------------------------------------

    NUM_LBFLX = adjust_flux_cell

    allocate( CORIOLI(IA,JA) )

    ! numerical diffusion
    call ATMOS_DYN_filter_setup( CDZ, CDX, CDY, FDZ, FDX, FDY ) ! (in)

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
       Rdry   => CONST_Rdry,  &
       Rvap   => CONST_Rvap,  &
       CVdry  => CONST_CVdry
    use scale_process, only: &
       PRC_HAS_E, &
       PRC_HAS_W, &
       PRC_HAS_N, &
       PRC_HAS_S
    use scale_comm, only: &
#ifdef DEBUG
       COMM_datatype, &
#endif
       COMM_vars8, &
       COMM_wait
    use scale_gridtrans, only: &
       I_XYZ, &
       I_XYW, &
       I_UYZ, &
       I_XVZ, &
       I_XY
#ifdef DEBUG
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
#ifdef HIST_TEND
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
    real(RP), intent(in)    :: J13G (KA,IA,JA,4) !< (1,3) element of Jacobian matrix
    real(RP), intent(in)    :: J23G (KA,IA,JA,4) !< (2,3) element of Jacobian matrix
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

    real(RP) :: diff(KA,IA,JA)
    real(RP) :: damp
#ifdef HIST_TEND
    real(RP) :: damp_t(KA,IA,JA)
#endif

    real(RP) :: num_diff  (KA,IA,JA,5,3)
    real(RP) :: num_diff_q(KA,IA,JA,3)

    ! For tracer advection
    real(RP) :: RHOQ_t   (KA,IA,JA)       ! tendency
    real(RP) :: mflx_hi  (KA,IA,JA,3)     ! rho * vel(x,y,z) @ (u,v,w)-face high order
    real(RP) :: mflx_av  (KA,IA,JA,3)     ! rho * vel(x,y,z) @ (u,v,w)-face average
    real(RP) :: tflx_hi  (KA,IA,JA,3)     ! rho * theta * vel(x,y,z) @ (u,v,w)-face high order
    real(RP) :: qflx_hi  (KA,IA,JA,3,QA)  ! rho * vel(x,y,z) * phi @ (u,v,w)-face high order
    real(RP) :: qflx_lo  (KA,IA,JA,3,QA)  ! rho * vel(x,y,z) * phi,  monotone flux
    real(RP) :: qflx_anti(KA,IA,JA,3,QA)  ! anti-diffusive flux

    ! lateral boundary flux
    real(RP) :: mflx_lb(KA,IA,JA,3)
    real(RP) :: tflx_lb(KA,IA,JA,3)
    real(RP) :: qflx_lb(KA,IA,JA,3,QA)
#ifdef DEBUG
    real(RP) :: mass_tmp (KA,IA,JA)
    real(RP) :: mass_tmp2(KA,IA,JA)
    real(RP) :: mass_total
    real(RP) :: mass_total2
    real(RP) :: allmass_total
    real(RP) :: allmass_total2
#endif

    real(RP) :: dt
    real(RP) :: dtrk

    integer  :: IIS, IIE
    integer  :: JJS, JJE
    integer  :: i, j, k, iq, step
#ifdef DEBUG
    integer :: ierr
#endif
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*) '*** Dynamics step: ', NSTEP_ATMOS_DYN, ' small steps'

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

    qflx_hi(:,:,:,:,:) = UNDEF
    qflx_lo(:,:,:,:,:) = UNDEF
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
       do j = JS-1, JE+2
       do i = IS-1, IE+2
       do k = KS, KE
          diff(k,i,j) = DENS(k,i,j)
!          diff(k,i,j) = DENS(k,i,j) - DAMP_DENS(k,i,j)
       enddo
       enddo
       enddo
       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          damp = - DAMP_alpha_DENS(k,i,j) * ( &
               - ( diff(k,i-1,j) + diff(k,i+1,j) + diff(k,i,j-1) + diff(k,i,j+1) - diff(k,i,j)*4.0_RP ) &
               * 0.125_RP * BND_SMOOTHER_FACT ) ! horizontal smoother
          DENS_t(k,i,j) = DENS_tp(k,i,j) & ! tendency from physical step
                        + damp
          DENS_t(k,i,j) = 0.0_RP
#ifdef HIST_TEND
          damp_t(k,i,j) = damp
#endif
       enddo
       enddo
       enddo
#ifdef HIST_TEND
       call HIST_in(damp_t, 'DENS_t_damp', 'tendency of dencity due to rayleigh damping', &
                            'kg/m3/s', DTSEC_ATMOS_DYN)
#endif

       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j = JS-1, JE+1
       do i = IS-1, IE+1
       do k = KS, KE-1
          diff(k,i,j) = MOMZ(k,i,j) - DAMP_VELZ(k,i,j) * ( DENS(k,i,j)+DENS(k+1,i,j) ) * 0.5_RP
       enddo
       enddo
       enddo
       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
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
       do j = JS, JE
       do i = IS, IE
          MOMZ_t(KE,i,j) = 0.0_RP
       enddo
       enddo
#ifdef HIST_TEND
       call HIST_in(damp_t, 'MOMZ_t_damp', 'tendency of momentum z due to rayleigh damping', &
                            'kg/m2/s2', DTSEC_ATMOS_DYN, zdim='half')
#endif

       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j = JS-1, JE+2
       do i = IS-2, IE+1
       do k = KS, KE
          diff(k,i,j) = MOMX(k,i,j) - DAMP_VELX(k,i,j) * ( DENS(k,i,j)+DENS(k,i+1,j) ) * 0.5_RP
       enddo
       enddo
       enddo
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
#ifdef HIST_TEND
       call HIST_in(damp_t, 'MOMX_t_damp', 'tendency of momentum x due to rayleigh damping', &
                            'kg/m2/s2', DTSEC_ATMOS_DYN, xdim='half')
#endif

       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j = JS-2, JE+1
       do i = IS-1, IE+2
       do k = KS, KE
          diff(k,i,j) = MOMY(k,i,j) - DAMP_VELY(k,i,j) * ( DENS(k,i,j)+DENS(k,i,j+1) ) * 0.5_RP
       enddo
       enddo
       enddo
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
#ifdef HIST_TEND
       call HIST_in(damp_t, 'MOMY_t_damp', 'tendency of momentum y due to rayleigh damping', &
                            'kg/m2/s2', DTSEC_ATMOS_DYN, ydim='half')
#endif

       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j = JS-1, JE+2
       do i = IS-1, IE+2
       do k = KS, KE
          diff(k,i,j) = RHOT(k,i,j) - DAMP_POTT(k,i,j) * DENS(k,i,j)
       enddo
       enddo
       enddo
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
#ifdef HIST_TEND
       call HIST_in(damp_t, 'RHOT_t_damp', 'tendency of rho*theta temperature due to rayleigh damping', &
                            'K kg/m3/s', DTSEC_ATMOS_DYN)
#endif

       do j = JS, JE
       do i = IS, IE
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

       dt = real(DTSEC_ATMOS_DYN,kind=RP)

       !-----< prepare numerical diffusion coefficient >-----

       if ( ND_COEF == 0.0_RP ) then
          num_diff(:,:,:,:,:) = 0.0_RP
       else
          call ATMOS_DYN_numfilter_coef( num_diff(:,:,:,:,:),                    & ! [OUT]
                                         DENS, MOMZ, MOMX, MOMY, RHOT,           & ! [IN]
                                         CDZ, CDX, CDY, FDZ, FDX, FDY, dt,       & ! [IN]
                                         REF_dens, REF_pott,                     & ! [IN]
                                         ND_COEF, ND_ORDER, ND_SFC_FACT, ND_USE_RS ) ! [IN]
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
                          dtrk, dt                                          ) ! (in)

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
                          dtrk, dt                                          ) ! (in)

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
                          dtrk, dt                                          ) ! (in)

       ! initialize lateral boundary flux
       mflx_lb(:,:,:,:) = 0.0_RP
       tflx_lb(:,:,:,:) = 0.0_RP

       ! adjust lateral boundary flux
       if ( .not. PRC_HAS_W ) then ! for western boundary
          call adjust_boundary_flux_dyn( &
               DENS, RHOT, mflx_lb, tflx_lb, & ! (inout)
               mflx_hi, tflx_hi, GSQRT, MAPF, dt, & ! (in)
               RCDX, RCDY, & ! (in)
               DAMP_DENS, DAMP_VELX, DAMP_VELY, DAMP_POTT, & ! (in)
               -1, 0, 1, 0, & ! (in)
               IS+IHALO+NUM_LBFLX-1, IS+IHALO, JS, JE ) ! (in)
       end if
       if ( .not. PRC_HAS_E ) then ! for eastern boundary
          call adjust_boundary_flux_dyn( &
               DENS, RHOT, mflx_lb, tflx_lb, & ! (inout)
               mflx_hi, tflx_hi, GSQRT, MAPF, dt, & ! (in)
               RCDX, RCDY, & ! (in)
               DAMP_DENS, DAMP_VELX, DAMP_VELY, DAMP_POTT, & ! (in)
               0, 0, -1, 0, & ! (in)
               IE-IHALO-NUM_LBFLX+1, IE-IHALO, JS, JE ) ! (in)
       end if
       if ( .not. PRC_HAS_S ) then ! for sourthern boundary
          call adjust_boundary_flux_dyn( &
               DENS, RHOT, mflx_lb, tflx_lb, & ! (inout)
               mflx_hi, tflx_hi, GSQRT, MAPF, dt, & ! (in)
               RCDX, RCDY, & ! (in)
               DAMP_DENS, DAMP_VELX, DAMP_VELY, DAMP_POTT, & ! (in)
               0, -1, 0, 1, & ! (in)
               IS, IE, JS+JHALO+NUM_LBFLX-1, JS+JHALO ) ! (in)
       end if
       if ( .not. PRC_HAS_N ) then ! for northern boundary
          call adjust_boundary_flux_dyn( &
               DENS, RHOT, mflx_lb, tflx_lb, & ! (inout)
               mflx_hi, tflx_hi, GSQRT, MAPF, dt, & ! (in)
               RCDX, RCDY, & ! (in)
               DAMP_DENS, DAMP_VELX, DAMP_VELY, DAMP_POTT, & ! (in)
               0, 0, 0, -1, & ! (in)
               IS, IE, JE-JHALO-NUM_LBFLX+1, JE-JHALO ) ! (in)
       end if

#ifdef HIST_TEND
       call HIST_in(mflx_lb(:,:,:,ZDIR), 'MOMZ_lb', 'momentum flux of z-direction from lateral boundary', &
                                         'kg/m2/s', dt, zdim='half'                                       )
       call HIST_in(mflx_lb(:,:,:,XDIR), 'MOMX_lb', 'momentum flux of x-direction from lateral boundary', &
                                         'kg/m2/s', dt, xdim='half'                                       )
       call HIST_in(mflx_lb(:,:,:,YDIR), 'MOMY_lb', 'momentum flux of y-direction from lateral boundary', &
                                         'kg/m2/s', dt, ydim='half'                                       )

       call HIST_in(tflx_lb(:,:,:,ZDIR), 'RHTZ_lb', 'potential temperature flux of z-direction from lateral boundary', &
                                         'K*kg/m2/s', dt, zdim='half'                                                  )
       call HIST_in(tflx_lb(:,:,:,XDIR), 'RHTX_lb', 'potential temperature flux of x-direction from lateral boundary', &
                                         'K*kg/m2/s', dt, xdim='half'                                                  )
       call HIST_in(tflx_lb(:,:,:,YDIR), 'RHTY_lb', 'potential temperature flux of y-direction from lateral boundary', &
                                         'K*kg/m2/s', dt, ydim='half'                                                  )
#endif

       ! set boundary
       if ( .not. PRC_HAS_W ) then ! for western boundary
          call set_boundary_dyn( &
               DENS, MOMZ, MOMX, MOMY, RHOT, & ! (inout)
               DAMP_DENS, DAMP_VELZ, DAMP_VELX, DAMP_VELY, DAMP_POTT, & ! (in)
               DAMP_alpha_DENS, DAMP_alpha_VELZ, DAMP_alpha_VELX, DAMP_alpha_VELY, DAMP_alpha_POTT, & ! (in)
               DAMP_VELX, & ! (in)
               -1, 0, 1, 0, & ! (in)
               IS+IHALO-1, IS, JS, JE ) ! (in)
       end if
       if ( .not. PRC_HAS_E ) then ! for eastern boundary
          call set_boundary_dyn( &
               DENS, MOMZ, MOMX, MOMY, RHOT, & ! (inout)
               DAMP_DENS, DAMP_VELZ, DAMP_VELX, DAMP_VELY, DAMP_POTT, & ! (in)
               DAMP_alpha_DENS, DAMP_alpha_VELZ, DAMP_alpha_VELX, DAMP_alpha_VELY, DAMP_alpha_POTT, & ! (in)
               DAMP_VELX, & ! (in)
               0, 0, -1, 0, & ! (in)
               IE-IHALO+1, IE, JS, JE ) ! (in)
       end if
       if ( .not. PRC_HAS_S ) then ! for sourthern boundary
          call set_boundary_dyn( &
               DENS, MOMZ, MOMX, MOMY, RHOT, & ! (inout)
               DAMP_DENS, DAMP_VELZ, DAMP_VELX, DAMP_VELY, DAMP_POTT, & ! (in)
               DAMP_alpha_DENS, DAMP_alpha_VELZ, DAMP_alpha_VELX, DAMP_alpha_VELY, DAMP_alpha_POTT, & ! (in)
               DAMP_VELY, & ! (in)
               0, -1, 0, 1, & ! (in)
               IS, IE, JS+JHALO-1, JS ) ! (in)
       end if
       if ( .not. PRC_HAS_N ) then ! for northern boundary
          call set_boundary_dyn( &
               DENS, MOMZ, MOMX, MOMY, RHOT, & ! (inout)
               DAMP_DENS, DAMP_VELZ, DAMP_VELX, DAMP_VELY, DAMP_POTT, & ! (in)
               DAMP_alpha_DENS, DAMP_alpha_VELZ, DAMP_alpha_VELX, DAMP_alpha_VELY, DAMP_alpha_POTT, & ! (in)
               DAMP_VELY, & ! (in)
               0, 0, 0, -1, & ! (in)
               IS, IE, JE-JHALO+1, JE ) ! (in)
       end if

#ifdef DEBUG
       ! check total mass in the inner region
       do j = JS, JE
       do i = IS, IE
       do k = KS, KE  
          mass_tmp (k,i,j) = DENS     (k,i,j) * vol(k,i,j)
          mass_tmp2(k,i,j) = DAMP_DENS(k,i,j) * vol(k,i,j)
       end do
       end do
       end do

       if ( .not. PRC_HAS_W ) then ! for western boundary
          do j = JS, JE
          do i = IS, IS+1
          do k = KS, KE
             mass_tmp (k,i,j) = 0.0_RP
             mass_tmp2(k,i,j) = 0.0_RP
          end do
          end do
          end do
       end if
       if ( .not. PRC_HAS_E ) then ! for eastern boundary
          do j = JS, JE
          do i = IE-1, IE
          do k = KS, KE
             mass_tmp (k,i,j) = 0.0_RP
             mass_tmp2(k,i,j) = 0.0_RP
          end do
          end do
          end do
       end if
       if ( .not. PRC_HAS_S ) then ! for sourthern boundary
          do j = JS, JS+1
          do i = IS, IE
          do k = KS, KE
             mass_tmp (k,i,j) = 0.0_RP
             mass_tmp2(k,i,j) = 0.0_RP
          end do
          end do
          end do
       end if
       if ( .not. PRC_HAS_N ) then ! for northern boundary
          do j = JE-1, JE
          do i = IS, IE
          do k = KS, KE
             mass_tmp (k,i,j) = 0.0_RP
             mass_tmp2(k,i,j) = 0.0_RP
          end do
          end do
          end do
       end if

       mass_total  = 0.0_RP
       mass_total2 = 0.0_RP

       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          mass_total  = mass_total  + mass_tmp (k,i,j)
          mass_total2 = mass_total2 + mass_tmp2(k,i,j)
       end do
       end do
       end do

       call MPI_Allreduce( mass_total,           &
                           allmass_total,        &
                           1,                    &
                           COMM_datatype,        &
                           MPI_SUM,              &
                           MPI_COMM_WORLD,       &
                           ierr                  )
       if( IO_L ) write(IO_FID_LOG,'(A,1x,i1,1x,2PE24.17)') 'total mass   :', step, allmass_total

       call MPI_Allreduce( mass_total2,          &
                           allmass_total2,       &
                           1,                    &
                           COMM_datatype,        &
                           MPI_SUM,              &
                           MPI_COMM_WORLD,       &
                           ierr                  )
       if( IO_L ) write(IO_FID_LOG,'(A,1x,i1,1x,2PE24.17)') 'total mass2  :', step, allmass_total2
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

       do iq = 1, QA
          do j = JS, JE
          do i = IS, IE
          do k = KS, KE
             QTRC(k,i,j,iq) = max( QTRC(k,i,j,iq) + RHOQ_tp(k,i,j,iq) * dt / DENS00(k,i,j), 0.0_RP )
          end do
          end do
          end do

          call COMM_vars8( QTRC(:,:,:,iq), 6 )
          call COMM_wait ( QTRC(:,:,:,iq), 6 )

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
    do iq = 1, QA

       if ( iq <= BND_QA ) then
          !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
          do j = JS-1, JE+2
          do i = IS-1, IE+2
          do k = KS, KE
             diff(k,i,j) = QTRC(k,i,j,iq) - DAMP_QTRC(k,i,j,iq)
          enddo
          enddo
          enddo
          !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
          do j = JS, JE
          do i = IS, IE
          do k = KS, KE
             damp = - DAMP_alpha_QTRC(k,i,j,iq) &
                    * ( diff(k,i,j) & ! rayleigh damping
                      - ( diff(k,i-1,j) + diff(k,i+1,j) + diff(k,i,j-1) + diff(k,i,j+1) - diff(k,i,j)*4.0_RP ) &
                      * 0.125_RP * BND_SMOOTHER_FACT ) ! horizontal smoother
             RHOQ_t(k,i,j) = damp * DENS00(k,i,j)
#ifdef HIST_TEND
             damp_t(k,i,j) = damp
#endif
          enddo
          enddo
          enddo
#ifdef HIST_TEND
          call HIST_in( damp_t, trim(AQ_NAME(iq))//'_t_damp', &
                        'tendency of '//trim(AQ_NAME(iq))//' due to rayleigh damping', 'kg/kg/s', DTSEC )
#endif
       else
          !$omp parallel do private(i,j,k,iq) OMP_SCHEDULE_ collapse(2)
          do j = JS, JE
          do i = IS, IE
          do k = KS, KE
             RHOQ_t(k,i,j) = 0.0_RP
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

       if ( ND_COEF_Q == 0.0_RP ) then
          num_diff_q(:,:,:,:) = 0.0_RP
       else
          call ATMOS_DYN_numfilter_coef_q( num_diff_q(:,:,:,:),                    & ! [OUT]
                                           DENS00, QTRC(:,:,:,iq),                 & ! [IN]
                                           CDZ, CDX, CDY, DT,                      & ! [IN]
                                           REF_qv, iq,                             & ! [IN]
                                           ND_COEF_Q, ND_ORDER, ND_SFC_FACT, ND_USE_RS ) ! [IN]
       endif

       do JJS = JS, JE, JBLOCK
       JJE = JJS+JBLOCK-1
       do IIS = IS, IE, IBLOCK
       IIE = IIS+IBLOCK-1

          !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
          do j = JJS-1, JJE+1
          do i = IIS-1, IIE+1
          do k = KS+1, KE-2
             qflx_lo(k,i,j,ZDIR,iq) = 0.5_RP * (     mflx_hi(k,i,j,ZDIR)  * ( QTRC(k+1,i,j,iq)+QTRC(k,i,j,iq) ) &
                                               - abs(mflx_hi(k,i,j,ZDIR)) * ( QTRC(k+1,i,j,iq)-QTRC(k,i,j,iq) ) )

             qflx_hi(k,i,j,ZDIR,iq) = mflx_hi(k,i,j,ZDIR) * ( FACT_N * ( QTRC(k+1,i,j,iq)+QTRC(k  ,i,j,iq) ) &
                                                            + FACT_F * ( QTRC(k+2,i,j,iq)+QTRC(k-1,i,j,iq) ) ) &
                                    + GSQRT(k,i,j,I_XYW) * num_diff_q(k,i,j,ZDIR)
          enddo
          enddo
          enddo

          !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
          do j = JJS-1, JJE+1
          do i = IIS-1, IIE+1
             qflx_lo(KS-1,i,j,ZDIR,iq) = 0.0_RP
             qflx_lo(KS  ,i,j,ZDIR,iq) = 0.5_RP * (     mflx_hi(KS  ,i,j,ZDIR)  * ( QTRC(KS+1,i,j,iq)+QTRC(KS,i,j,iq) ) &
                                                  - abs(mflx_hi(KS  ,i,j,ZDIR)) * ( QTRC(KS+1,i,j,iq)-QTRC(KS,i,j,iq) ) )

             qflx_hi(KS-1,i,j,ZDIR,iq) = 0.0_RP
             qflx_hi(KS  ,i,j,ZDIR,iq) = 0.5_RP * mflx_hi(KS  ,i,j,ZDIR) * ( QTRC(KS+1,i,j,iq)+QTRC(KS,i,j,iq) ) &
                                       + GSQRT(k,i,j,I_XYW) * num_diff_q(KS,i,j,ZDIR)
          enddo
          enddo

          !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
          do j = JJS-1, JJE+1
          do i = IIS-1, IIE+1
             qflx_lo(KE-1,i,j,ZDIR,iq) = 0.5_RP * (     mflx_hi(KE-1,i,j,ZDIR)  * ( QTRC(KE,i,j,iq)+QTRC(KE-1,i,j,iq) ) &
                                                  - abs(mflx_hi(KE-1,i,j,ZDIR)) * ( QTRC(KE,i,j,iq)-QTRC(KE-1,i,j,iq) ) )
             qflx_lo(KE  ,i,j,ZDIR,iq) = 0.0_RP

             qflx_hi(KE-1,i,j,ZDIR,iq) = 0.5_RP * mflx_hi(KE-1,i,j,ZDIR) * ( QTRC(KE,i,j,iq)+QTRC(KE-1,i,j,iq) ) &
                                       + GSQRT(k,i,j,I_XYW) * num_diff_q(KE-1,i,j,ZDIR)
             qflx_hi(KE  ,i,j,ZDIR,iq) = 0.0_RP
          enddo
          enddo

          !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
          do j = JJS-1, JJE+1
          do i = IIS-2, IIE+1
          do k = KS, KE
             qflx_lo(k,i,j,XDIR,iq) = 0.5_RP * (     mflx_hi(k,i,j,XDIR)  * ( QTRC(k,i+1,j,iq)+QTRC(k,i,j,iq) ) &
                                               - abs(mflx_hi(k,i,j,XDIR)) * ( QTRC(k,i+1,j,iq)-QTRC(k,i,j,iq) ) )
          enddo
          enddo
          enddo

          !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
          do j = JJS,   JJE
          do i = IIS-1, IIE
          do k = KS, KE
             qflx_hi(k,i,j,XDIR,iq) = mflx_hi(k,i,j,XDIR) * ( FACT_N * ( QTRC(k,i+1,j,iq)+QTRC(k,i  ,j,iq) ) &
                                                            + FACT_F * ( QTRC(k,i+2,j,iq)+QTRC(k,i-1,j,iq) ) ) &
                                    + GSQRT(k,i,j,I_UYZ) * num_diff_q(k,i,j,XDIR)
          enddo
          enddo
          enddo

          !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
          do j = JJS-2, JJE+1
          do i = IIS-1, IIE+1
          do k = KS, KE
             qflx_lo(k,i,j,YDIR,iq) = 0.5_RP * (     mflx_hi(k,i,j,YDIR)  * ( QTRC(k,i,j+1,iq)+QTRC(k,i,j,iq) ) &
                                               - abs(mflx_hi(k,i,j,YDIR)) * ( QTRC(k,i,j+1,iq)-QTRC(k,i,j,iq) ) )
          enddo
          enddo
          enddo

          !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
          do j = JJS-1, JJE
          do i = IIS,   IIE
          do k = KS, KE
             qflx_hi(k,i,j,YDIR,iq) = mflx_hi(k,i,j,YDIR) * ( FACT_N * ( QTRC(k,i,j+1,iq)+QTRC(k,i,j  ,iq) ) &
                                                            + FACT_F * ( QTRC(k,i,j+2,iq)+QTRC(k,i,j-1,iq) ) ) &
                                    + GSQRT(k,i,j,I_XVZ) * num_diff_q(k,i,j,YDIR)
          enddo
          enddo
          enddo

       enddo
       enddo

       call ATMOS_DYN_fct( qflx_anti(:,:,:,:,iq),        & ! (out)
                           QTRC(:,:,:,iq), DENS00, DENS, & ! (in)
                           qflx_hi(:,:,:,:,iq),          & ! (in)
                           qflx_lo(:,:,:,:,iq),          & ! (in)
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
                            + dt * ( - ( ( qflx_hi(k  ,i  ,j  ,ZDIR,iq) - qflx_anti(k  ,i  ,j  ,ZDIR,iq) &
                                         - qflx_hi(k-1,i  ,j  ,ZDIR,iq) + qflx_anti(k-1,i  ,j  ,ZDIR,iq) ) * RCDZ(k) &
                                       + ( qflx_hi(k  ,i  ,j  ,XDIR,iq) - qflx_anti(k  ,i  ,j  ,XDIR,iq) &
                                         - qflx_hi(k  ,i-1,j  ,XDIR,iq) + qflx_anti(k  ,i-1,j  ,XDIR,iq) ) * RCDX(i) &
                                       + ( qflx_hi(k  ,i  ,j  ,YDIR,iq) - qflx_anti(k  ,i  ,j  ,YDIR,iq) &
                                         - qflx_hi(k  ,i  ,j-1,YDIR,iq) + qflx_anti(k  ,i  ,j-1,YDIR,iq) ) * RCDY(j) &
                                       ) * MAPF(i,j,1,I_XY) * MAPF(i,j,2,I_XY) / GSQRT(k,i,j,I_XYZ) &
                              + RHOQ_t(k,i,j) ) ) / DENS(k,i,j)
          enddo
          enddo
          enddo

       enddo
       enddo

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

    ! initialize lateral boundary flux
    qflx_lb(:,:,:,:,:) = 0.0_RP

    ! adjust lateral boundary flux
    if ( .not. PRC_HAS_W ) then ! for western boundary
       call adjust_boundary_flux_qtrc( &
            QTRC, qflx_lb, & ! (inout)
            DENS, MOMX, MOMY, & ! (in)
            qflx_hi, qflx_anti, GSQRT, MAPF, dt, & ! (in)
            RCDX, RCDY, & ! (in)
            DAMP_QTRC, & ! (in)
            -1, 0, 1, 0, & ! (in)
            IS+IHALO+NUM_LBFLX-1, IS+IHALO, JS, JE ) ! (in)
    end if
    if ( .not. PRC_HAS_E ) then ! for eastern boundary
       call adjust_boundary_flux_qtrc( &
            QTRC, qflx_lb, & ! (inout)
            DENS, MOMX, MOMY, & ! (in)
            qflx_hi, qflx_anti, GSQRT, MAPF, dt, & ! (in)
            RCDX, RCDY, & ! (in)
            DAMP_QTRC, & ! (in)
            0, 0, -1, 0, & ! (in)
            IE-IHALO-NUM_LBFLX+1, IE-IHALO, JS, JE ) ! (in)
    end if
    if ( .not. PRC_HAS_S ) then ! for sourthern boundary
       call adjust_boundary_flux_qtrc( &
            QTRC, qflx_lb, & ! (inout)
            DENS, MOMX, MOMY, & ! (in)
            qflx_hi, qflx_anti, GSQRT, MAPF, dt, & ! (in)
            RCDX, RCDY, & ! (in)
            DAMP_QTRC, & ! (in)
            0, -1, 0, 1, & ! (in)
            IS, IE, JS+JHALO+NUM_LBFLX-1, JS+JHALO ) ! (in)
    end if
    if ( .not. PRC_HAS_N ) then ! for northern boundary
       call adjust_boundary_flux_qtrc( &
            QTRC, qflx_lb, & ! (inout)
            DENS, MOMX, MOMY, & ! (in)
            qflx_hi, qflx_anti, GSQRT, MAPF, dt, & ! (in)
            RCDX, RCDY, & ! (in)
            DAMP_QTRC, & ! (in)
            0, 0, 0, -1, & ! (in)
            IS, IE, JE-JHALO-NUM_LBFLX+1, JE-JHALO ) ! (in)
    end if

#ifdef HIST_TEND
    do iq = 1, QA
       call HIST_in(qflx_lb(:,:,:,ZDIR,iq), trim(AQ_NAME(iq))//'_lb',                                  &
                                            trim(AQ_NAME(iq))//'of z-direction from lateral boundary', &
                                            'kg/m2/s', DTSEC, zdim='half'                              )
       call HIST_in(qflx_lb(:,:,:,XDIR,iq), trim(AQ_NAME(iq))//'_lb',                                  &
                                            trim(AQ_NAME(iq))//'of x-direction from lateral boundary', &
                                            'kg/m2/s', DTSEC, xdim='half'                              )
       call HIST_in(qflx_lb(:,:,:,YDIR,iq), trim(AQ_NAME(iq))//'_lb',                                  &
                                            trim(AQ_NAME(iq))//'of y-direction from lateral boundary', &
                                            'kg/m2/s', DTSEC, ydim='half'                              )
    end do
#endif

    if ( .not. PRC_HAS_W ) then ! for western boundary
       ! set boundary
       call set_boundary_qtrc( &
            QTRC, & ! (inout)
            DAMP_QTRC, DAMP_alpha_QTRC, & ! (in)
            MOMX, & ! (in)
            -1, 0, 1, 0, & ! (in)
            IS+IHALO-1, IS, JS, JE ) ! (in)
    end if
    if ( .not. PRC_HAS_E ) then ! for eastern boundary
       call set_boundary_qtrc( &
            QTRC, & ! (inout)
            DAMP_QTRC, DAMP_alpha_QTRC, & ! (in)
            MOMX, & ! (in)
            0, 0, -1, 0, & ! (in)
            IE-IHALO+1, IE, JS, JE ) ! (in)
    end if
    if ( .not. PRC_HAS_S ) then ! for sourthern boundary
       call set_boundary_qtrc( &
            QTRC, & ! (inout)
            DAMP_QTRC, DAMP_alpha_QTRC, & ! (in)
            MOMY, & ! (in)
            0, -1, 0, 1, & ! (in)
            IS, IE, JS+JHALO-1, JS ) ! (in)
    end if
    if ( .not. PRC_HAS_N ) then ! for northern boundary
       call set_boundary_qtrc( &
            QTRC, & ! (inout)
            DAMP_QTRC, DAMP_alpha_QTRC, & ! (in)
            MOMY, & ! (in)
            0, 0, 0, -1, & ! (in)
            IS, IE, JE-JHALO+1, JE ) ! (in)
    end if

    do iq = 1, QA
       call COMM_vars8( QTRC(:,:,:,iq), iq )
    enddo
    do iq = 1, QA
       call COMM_wait ( QTRC(:,:,:,iq), iq )
    enddo
#endif

    return
  end subroutine ATMOS_DYN

  subroutine set_boundary_dyn( &
       DENS, MOMZ, MOMX, MOMY, RHOT, &
       DAMP_DENS, DAMP_VELZ, DAMP_VELX, DAMP_VELY, DAMP_POTT, &
       DAMP_alpha_DENS, DAMP_alpha_VELZ, DAMP_alpha_VELX, DAMP_alpha_VELY, DAMP_alpha_POTT, &
       DAMP_flow, &
       ib, jb, iu, ju, &
       i0, i1, j0, j1 )
    use scale_const, only: &
       EPS => CONST_EPS
    implicit none

    real(RP), intent(inout) :: DENS(KA,IA,JA)
    real(RP), intent(inout) :: MOMZ(KA,IA,JA)
    real(RP), intent(inout) :: MOMX(KA,IA,JA)
    real(RP), intent(inout) :: MOMY(KA,IA,JA)
    real(RP), intent(inout) :: RHOT(KA,IA,JA)

    real(RP), intent(in)    :: DAMP_DENS(KA,IA,JA)
    real(RP), intent(in)    :: DAMP_VELZ(KA,IA,JA)
    real(RP), intent(in)    :: DAMP_VELX(KA,IA,JA)
    real(RP), intent(in)    :: DAMP_VELY(KA,IA,JA)
    real(RP), intent(in)    :: DAMP_POTT(KA,IA,JA)

    real(RP), intent(in)    :: DAMP_alpha_DENS(KA,IA,JA)
    real(RP), intent(in)    :: DAMP_alpha_VELZ(KA,IA,JA)
    real(RP), intent(in)    :: DAMP_alpha_VELX(KA,IA,JA)
    real(RP), intent(in)    :: DAMP_alpha_VELY(KA,IA,JA)
    real(RP), intent(in)    :: DAMP_alpha_POTT(KA,IA,JA)

    real(RP), intent(in)    :: DAMP_flow(KA,IA,JA)

    integer,  intent(in)    :: ib ! W:-1, E: 0, S: 0, N: 0
    integer,  intent(in)    :: jb ! W: 0, E: 0, S:-1, N: 0
    integer,  intent(in)    :: iu ! W: 1, E:-1, S: 0, N: 0
    integer,  intent(in)    :: ju ! W: 0, E: 0, S: 1, N:-1
    integer,  intent(in)    :: i0
    integer,  intent(in)    :: i1
    integer,  intent(in)    :: j0
    integer,  intent(in)    :: j1

    real(RP) :: dir
    real(RP) :: sw, sw2

    integer :: i, j, k
    !---------------------------------------------------------------------------

    dir = sign(1.0_RP, real(ib+jb,kind=RP)) ! dir = -1 if ib==-1 .or. jb==-1

    do j = j0, j1, -ju+abs(iu)
    do i = i0, i1, -iu+abs(ju)
    do k = KS, KE
       sw = sign(0.5_RP, DAMP_alpha_DENS(k,i,j) - EPS) + 0.5_RP
       sw2 = sign(0.5_RP, DAMP_flow(k,i,j)*dir) + 0.5_RP ! 0:inflow, 1:outflow
       DENS(k,i,j) = DENS(k,i,j) * ( 1.0_RP - sw ) &
                   + DAMP_DENS(k,i,j) * sw
       MOMZ(k,i,j) = MOMZ(k,i,j) * ( 1.0_RP - sw ) &
                   + MOMZ(k,i+iu,j+ju) * sw2 * sw
       sw = sign(0.5_RP, DAMP_alpha_POTT(k,i,j) - EPS) + 0.5_RP
       RHOT(k,i,j) = RHOT(k,i,j) * ( 1.0_RP - sw ) &
                   + DAMP_POTT(k,i,j) * DAMP_DENS(k,i,j) * sw
    end do
    end do
    end do

    do j = j0,       j1,       -ju+abs(iu)
    do i = i0+ib+iu, i1+ib+iu, -iu+abs(ju)
    do k = KS, KE
       sw = sign(0.5_RP, DAMP_alpha_VELX(k,i,j) - EPS) + 0.5_RP
       MOMX(k,i,j) = MOMX(k,i,j) * ( 1.0_RP - sw ) &
                   + DAMP_VELX(k,i,j) * ( DENS(k,i,j)+DENS(k,i+1,j) ) * 0.5_RP * sw
    end do
    end do
    end do

    do j = j0+jb+ju, j1+jb+ju, -ju+abs(iu)
    do i = i0,       i1,       -iu+abs(ju)
    do k = KS, KE
       sw = sign(0.5_RP, DAMP_alpha_VELY(k,i,j) - EPS) + 0.5_RP
       MOMY(k,i,j) = MOMY(k,i,j) * ( 1.0_RP - sw ) &
                   + DAMP_VELY(k,i,j) * ( DENS(k,i,j)+DENS(k,i,j+1) ) * 0.5_RP * sw
    end do
    end do
    end do

    if ( iu==-1 ) then ! eastern boundary
       do j = j0, j1
       do k = KS, KE
          sw = sign(0.5_RP, DAMP_alpha_VELX(k,i1,j) - EPS) + 0.5_RP
          MOMX(k,i1,j) = MOMX(k,i1,j) * ( 1.0_RP - sw ) &
                       + MOMX(k,i1+iu,j) * sw
       end do
       end do
    end if
    if ( ju==-1 ) then ! northern boundary
       do i = i0, i1
       do k = KS, KE
          sw = sign(0.5_RP, DAMP_alpha_VELY(k,i,j1) - EPS) + 0.5_RP
          MOMY(k,i,j1) = MOMY(k,i,j1) * ( 1.0_RP - sw ) &
                       + MOMY(k,i,j1+ju) * sw
       end do
       end do
    end if

    return
  end subroutine set_boundary_dyn

  subroutine set_boundary_qtrc( &
       QTRC, &
       DAMP_QTRC, DAMP_alpha_QTRC, &
       MOM, &
       ib, jb, iu, ju, &
       i0, i1, j0, j1 )
    use scale_const, only: &
       EPS => CONST_EPS
    use scale_atmos_boundary, only: &
       BND_QA
    implicit none
    real(RP), intent(inout) :: QTRC(KA,IA,JA,QA)

    real(RP), intent(in)    :: DAMP_QTRC      (KA,IA,JA,BND_QA)
    real(RP), intent(in)    :: DAMP_alpha_QTRC(KA,IA,JA,BND_QA)

    real(RP), intent(in)    :: MOM(KA,IA,JA)

    integer,  intent(in)    :: ib ! W:-1, E: 0, S: 0, N: 0
    integer,  intent(in)    :: jb ! W: 0, E: 0, S:-1, N: 0
    integer,  intent(in)    :: iu ! W: 1, E:-1, S: 0, N: 0
    integer,  intent(in)    :: ju ! W: 0, E: 0, S: 1, N:-1
    integer,  intent(in)    :: i0
    integer,  intent(in)    :: i1
    integer,  intent(in)    :: j0
    integer,  intent(in)    :: j1

    real(RP) :: dir
    real(RP) :: sw, sw2

    integer :: i, j, k, iq
    !---------------------------------------------------------------------------

    if( BND_QA == I_QV ) then
       ! QV forcing only at boundary
       dir = sign(1.0_RP, REAL(ib+jb,RP)) ! dir = -1 if ib==-1 .or. jb==-1

       do j = j0, j1, -ju+abs(iu)
       do i = i0, i1, -iu+abs(ju)
       do k = KS, KE
          sw = sign(0.5_RP, DAMP_alpha_QTRC(k,i,j,I_QV) - EPS) + 0.5_RP
          QTRC(k,i,j,I_QV) = QTRC(k,i,j,I_QV) * ( 1.0_RP - sw ) &
                           + DAMP_QTRC(k,i,j,I_QV) * sw
       end do
       end do
       end do

       do iq = 2, QA
          do j = j0, j1, -ju+abs(iu)
          do i = i0, i1, -iu+abs(ju)
          do k = KS, KE
             sw = sign(0.5_RP, DAMP_alpha_QTRC(k,i,j,I_QV) - EPS) + 0.5_RP
             sw2 = sign(0.5_RP, MOM(k,i-ib,j-jb)*dir) + 0.5_RP ! 0:inflow, 1:outflow
             QTRC(k,i,j,iq) = QTRC(k,i,j,iq) * ( 1.0_RP - sw ) &
                            + QTRC(k,i+iu,j+ju,iq) * sw * sw2
          end do
          end do
          end do
       end do
    else
       ! ALL QTRC forcing at boundary
       do iq = 1, QA
          do j = j0, j1, -ju+abs(iu)
          do i = i0, i1, -iu+abs(ju)
          do k = KS, KE
             sw = sign(0.5_RP, DAMP_alpha_QTRC(k,i,j,iq) - EPS) + 0.5_RP
             QTRC(k,i,j,iq) = QTRC(k,i,j,iq) * ( 1.0_RP - sw ) &
                            + DAMP_QTRC(k,i,j,iq) * sw
          end do
          end do
          end do
       end do
    end if

    return
  end subroutine set_boundary_qtrc

  subroutine adjust_boundary_flux_dyn( &
       DENS, RHOT, mflx_lb, tflx_lb, &
       mflx_hi, tflx_hi, GSQRT, MAPF, dt, &
       RCDX, RCDY, &
       DAMP_DENS, DAMP_VELX, DAMP_VELY, DAMP_POTT, &
       ib, jb, iu, ju, &
       i0, i1, j0, j1 )
    use scale_gridtrans, only: &
       I_XYZ, &
       I_UYZ, &
       I_XVZ, &
       I_XY,  &
       I_UY,  &
       I_XV
    implicit none

    real(RP), intent(inout) :: DENS   (KA,IA,JA)
    real(RP), intent(inout) :: RHOT   (KA,IA,JA)
    real(RP), intent(inout) :: mflx_lb(KA,IA,JA,3)
    real(RP), intent(inout) :: tflx_lb(KA,IA,JA,3)

    real(RP), intent(in)    :: mflx_hi (KA,IA,JA,3)
    real(RP), intent(in)    :: tflx_hi (KA,IA,JA,3)
    real(RP), intent(in)    :: GSQRT   (KA,IA,JA,7)
    real(RP), intent(in)    :: MAPF    (IA,JA,2,4)
    real(RP), intent(in)    :: dt

    real(RP), intent(in)    :: RCDX(IA)
    real(RP), intent(in)    :: RCDY(JA)

    real(RP), intent(in)    :: DAMP_DENS(KA,IA,JA)
    real(RP), intent(in)    :: DAMP_VELX(KA,IA,JA)
    real(RP), intent(in)    :: DAMP_VELY(KA,IA,JA)
    real(RP), intent(in)    :: DAMP_POTT(KA,IA,JA)

    integer,  intent(in)    :: ib ! W:-1, E: 0, S: 0, N: 0
    integer,  intent(in)    :: jb ! W: 0, E: 0, S:-1, N: 0
    integer,  intent(in)    :: iu ! W: 1, E:-1, S: 0, N: 0
    integer,  intent(in)    :: ju ! W: 0, E: 0, S: 1, N:-1
    integer,  intent(in)    :: i0
    integer,  intent(in)    :: i1
    integer,  intent(in)    :: j0
    integer,  intent(in)    :: j1

    real(RP) :: DAMP_RHOT(KA,IA,JA)

    integer :: i, j, k
    !---------------------------------------------------------------------------

    do j = j0, j1, -ju+abs(iu)
    do i = i0, i1, -iu+abs(ju)
    do k = KS, KE
       DAMP_RHOT(k,i1+ib   ,j       ) = DAMP_POTT(k,i1+ib   ,j       ) * DAMP_DENS(k,i1+ib   ,j       )
       DAMP_RHOT(k,i       ,j1+jb   ) = DAMP_POTT(k,i       ,j1+jb   ) * DAMP_DENS(k,i       ,j1+jb   )
       DAMP_RHOT(k,i1-ib-iu,j       ) = DAMP_POTT(k,i1-ib-iu,j       ) * DAMP_DENS(k,i1-ib-iu,j       )
       DAMP_RHOT(k,i       ,j1-jb-ju) = DAMP_POTT(k,i       ,j1-jb-ju) * DAMP_DENS(k,i       ,j1-jb-ju)

       mflx_lb(k,i1+ib,j,XDIR) = DAMP_VELX(k,i1+ib,j) * ( DAMP_DENS(k,i1+ib,j) + DAMP_DENS(k,i1-ib-iu,j) ) * 0.5_RP &
                               * GSQRT(k,i1+ib,j,I_UYZ) / MAPF(i1+ib,j,2,I_UY) * real( abs(iu), kind=RP )
       mflx_lb(k,i,j1+jb,YDIR) = DAMP_VELY(k,i,j1+jb) * ( DAMP_DENS(k,i,j1+jb) + DAMP_DENS(k,i,j1-jb-ju) ) * 0.5_RP &
                               * GSQRT(k,i,j1+jb,I_XVZ) / MAPF(i,j1+jb,1,I_XV) * real( abs(ju), kind=RP )

       tflx_lb(k,i1+ib,j,XDIR) = DAMP_VELX(k,i1+ib,j) * ( DAMP_RHOT(k,i1+ib,j) + DAMP_RHOT(k,i1-ib-iu,j) ) * 0.5_RP &
                               * GSQRT(k,i1+ib,j,I_UYZ) / MAPF(i1+ib,j,2,I_UY) * real( abs(iu), kind=RP )
       tflx_lb(k,i,j1+jb,YDIR) = DAMP_VELY(k,i,j1+jb) * ( DAMP_RHOT(k,i,j1+jb) + DAMP_RHOT(k,i,j1-jb-ju) ) * 0.5_RP &
                               * GSQRT(k,i,j1+jb,I_XVZ) / MAPF(i,j1+jb,1,I_XV) * real( abs(ju), kind=RP )

       DENS(k,i,j) = DENS(k,i,j) &
                   + dt * ( ( mflx_lb(k,i1+ib,j    ,XDIR) - mflx_hi(k,i1+ib,j    ,XDIR) ) * real(iu,kind=RP) * RCDX(i) &
                          + ( mflx_lb(k,i    ,j1+jb,YDIR) - mflx_hi(k,i    ,j1+jb,YDIR) ) * real(ju,kind=RP) * RCDY(j) ) &
                        * MAPF(i,j,1,I_XY) * MAPF(i,j,2,I_XY) / GSQRT(k,i,j,I_XYZ) / real( NUM_LBFLX, kind=RP )

       RHOT(k,i,j) = RHOT(k,i,j) &
                   + dt * ( ( tflx_lb(k,i1+ib,j    ,XDIR) - tflx_hi(k,i1+ib,j    ,XDIR) ) * real(iu,kind=RP) * RCDX(i) &
                          + ( tflx_lb(k,i    ,j1+jb,YDIR) - tflx_hi(k,i    ,j1+jb,YDIR) ) * real(ju,kind=RP) * RCDY(j) ) &
                        * MAPF(i,j,1,I_XY) * MAPF(i,j,2,I_XY) / GSQRT(k,i,j,I_XYZ) / real( NUM_LBFLX, kind=RP )
    end do
    end do
    end do

    return
  end subroutine adjust_boundary_flux_dyn

  subroutine adjust_boundary_flux_qtrc( &
       QTRC, qflx_lb, &
       DENS, MOMX, MOMY, &
       qflx_hi, qflx_anti, GSQRT, MAPF, dt, &
       RCDX, RCDY, &
       DAMP_QTRC, &
       ib, jb, iu, ju, &
       i0, i1, j0, j1 )
    use scale_gridtrans, only: &
       I_XYZ, &
       I_UYZ, &
       I_XVZ, &
       I_XY,  &
       I_UY,  &
       I_XV
    use scale_atmos_boundary, only: &
       BND_QA
    implicit none
    real(RP), intent(inout) :: QTRC   (KA,IA,JA,QA)
    real(RP), intent(inout) :: qflx_lb(KA,IA,JA,3,QA)

    real(RP), intent(in)    :: DENS     (KA,IA,JA)
    real(RP), intent(in)    :: MOMX     (KA,IA,JA)
    real(RP), intent(in)    :: MOMY     (KA,IA,JA)

    real(RP), intent(in)    :: qflx_hi  (KA,IA,JA,3,QA)
    real(RP), intent(in)    :: qflx_anti(KA,IA,JA,3,QA)
    real(RP), intent(in)    :: GSQRT    (KA,IA,JA,7)
    real(RP), intent(in)    :: MAPF     (IA,JA,2,4)
    real(RP), intent(in)    :: dt

    real(RP), intent(in)    :: RCDX(IA)
    real(RP), intent(in)    :: RCDY(JA)

    real(RP), intent(in)    :: DAMP_QTRC      (KA,IA,JA,BND_QA)

    integer,  intent(in)    :: ib ! W:-1, E: 0, S: 0, N: 0
    integer,  intent(in)    :: jb ! W: 0, E: 0, S:-1, N: 0
    integer,  intent(in)    :: iu ! W: 1, E:-1, S: 0, N: 0
    integer,  intent(in)    :: ju ! W: 0, E: 0, S: 1, N:-1
    integer,  intent(in)    :: i0
    integer,  intent(in)    :: i1
    integer,  intent(in)    :: j0
    integer,  intent(in)    :: j1

    integer :: i, j, k, iq
    !---------------------------------------------------------------------------

    if( BND_QA == I_QV ) then
       ! QV forcing only at boundary
       do j = j0, j1, -ju+abs(iu)
       do i = i0, i1, -iu+abs(ju)
       do k = KS, KE
          qflx_lb(k,i1+ib,j,XDIR,I_QV) = MOMX(k,i1+ib,j) &
                                       * ( DAMP_QTRC(k,i1+ib,j,I_QV) + DAMP_QTRC(k,i1-ib-iu,j,I_QV) ) * 0.5_RP &
                                       * GSQRT(k,i1+ib,j,I_UYZ) / MAPF(i1+ib,j,2,I_UY) * real( abs(iu), kind=RP )
          qflx_lb(k,i,j1+jb,YDIR,I_QV) = MOMY(k,i,j1+jb) &
                                       * ( DAMP_QTRC(k,i,j1+jb,I_QV) + DAMP_QTRC(k,i,j1-jb-ju,I_QV) ) * 0.5_RP &
                                       * GSQRT(k,i,j1+jb,I_XVZ) / MAPF(i,j1+jb,1,I_XV) * real( abs(ju), kind=RP )

          QTRC(k,i,j,I_QV) = QTRC(k,i,j,I_QV) &
                           + dt * ( ( qflx_lb  (k,i1+ib,j,XDIR,I_QV) &
                                    - qflx_hi  (k,i1+ib,j,XDIR,I_QV) &
                                    + qflx_anti(k,i1+ib,j,XDIR,I_QV) ) * real(iu,kind=RP) * RCDX(i) &
                                  + ( qflx_lb  (k,i,j1+jb,YDIR,I_QV) &
                                    - qflx_hi  (k,i,j1+jb,YDIR,I_QV) &
                                    + qflx_anti(k,i,j1+jb,YDIR,I_QV) ) * real(ju,kind=RP) * RCDY(j) ) &
                                * MAPF(i,j,1,I_XY) * MAPF(i,j,2,I_XY) / GSQRT(k,i,j,I_XYZ) &
                                / DENS(k,i,j) / real( NUM_LBFLX, kind=RP )
       end do
       end do
       end do
    else
       ! ALL QTRC forcing at boundary
       do iq = 1, QA
          do j = j0, j1, -ju+abs(iu)
          do i = i0, i1, -iu+abs(ju)
          do k = KS, KE
             qflx_lb(k,i1+ib,j,XDIR,iq) = MOMX(k,i1+ib,j) &
                                        * ( DAMP_QTRC(k,i1+ib,j,iq) + DAMP_QTRC(k,i1-ib-iu,j,iq) ) * 0.5_RP &
                                        * GSQRT(k,i1+ib,j,I_UYZ) / MAPF(i1+ib,j,2,I_UY) * real( abs(iu), kind=RP )
             qflx_lb(k,i,j1+jb,YDIR,iq) = MOMY(k,i,j1+jb) &
                                        * ( DAMP_QTRC(k,i,j1+jb,iq) + DAMP_QTRC(k,i,j1-jb-ju,iq) ) * 0.5_RP &
                                        * GSQRT(k,i,j1+jb,I_XVZ) / MAPF(i,j1+jb,1,I_XV) * real( abs(ju), kind=RP )

             QTRC(k,i,j,iq) = QTRC(k,i,j,iq) &
                            + dt * ( ( qflx_lb  (k,i1+ib,j,XDIR,iq) &
                                     - qflx_hi  (k,i1+ib,j,XDIR,iq) &
                                     + qflx_anti(k,i1+ib,j,XDIR,iq) ) * real(iu,kind=RP) * RCDX(i) &
                                   + ( qflx_lb  (k,i,j1+jb,YDIR,iq) &
                                     - qflx_hi  (k,i,j1+jb,YDIR,iq) &
                                     + qflx_anti(k,i,j1+jb,YDIR,iq) ) * real(ju,kind=RP) * RCDY(j) ) &
                                 * MAPF(i,j,1,I_XY) * MAPF(i,j,2,I_XY) / GSQRT(k,i,j,I_XYZ) &
                                 / DENS(k,i,j) / real( NUM_LBFLX, kind=RP )
          end do
          end do
          end do
       end do
    end if

    return
  end subroutine adjust_boundary_flux_qtrc

end module scale_atmos_dyn
