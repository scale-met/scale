!-------------------------------------------------------------------------------
!> module Atmosphere / Dynamics
!!
!! @par Description
!!          HEVE FVM scheme for large time step in Atmospheric dynamical process
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2016-04-22 (S.Nishizawa) [new] split from scale_atmos_dyn.F90
!!
!<
!-------------------------------------------------------------------------------
#include "inc_openmp.h"
module scale_atmos_dyn_tstep_large_fvm_heve
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
  public :: ATMOS_DYN_Tstep_large_fvm_heve_setup
  public :: ATMOS_DYN_Tstep_large_fvm_heve

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
  ! tendency
  real(RP), private, allocatable :: DENS_t(:,:,:)
  real(RP), private, allocatable :: MOMZ_t(:,:,:)
  real(RP), private, allocatable :: MOMX_t(:,:,:)
  real(RP), private, allocatable :: MOMY_t(:,:,:)
  real(RP), private, allocatable :: RHOT_t(:,:,:)
  real(RP), private, allocatable :: RHOQ_t(:,:,:,:)

  real(RP), private, allocatable :: mflx_hi(:,:,:,:)        ! rho * vel(x,y,z) @ (u,v,w)-face high order

  ! for communication
  integer :: I_COMM_DENS = 1
  integer :: I_COMM_MOMZ = 2
  integer :: I_COMM_MOMX = 3
  integer :: I_COMM_MOMY = 4
  integer :: I_COMM_RHOT = 5
  integer, allocatable :: I_COMM_PROG(:)

  integer :: I_COMM_DENS_t = 1
  integer :: I_COMM_MOMZ_t = 2
  integer :: I_COMM_MOMX_t = 3
  integer :: I_COMM_MOMY_t = 4
  integer :: I_COMM_RHOT_t = 5

  integer, allocatable :: I_COMM_RHOQ_t(:)
  integer, allocatable :: I_COMM_QTRC(:)

  integer :: I_COMM_mflx_z = 1
  integer :: I_COMM_mflx_x = 2
  integer :: I_COMM_mflx_y = 3

  !-----------------------------------------------------------------------------
contains

  !-----------------------------------------------------------------------------
  !> Setup
  subroutine ATMOS_DYN_Tstep_large_fvm_heve_setup( &
       DENS, MOMZ, MOMX, MOMY, RHOT, QTRC, PROG, &
       mflx_hi )
    use scale_process, only: &
       PRC_MPIstop
    use scale_rm_process, only: &
       PRC_HAS_E, &
       PRC_HAS_W, &
       PRC_HAS_N, &
       PRC_HAS_S
    use scale_const, only: &
       OHM => CONST_OHM, &
       UNDEF => CONST_UNDEF
    use scale_comm, only: &
       COMM_vars8_init
    implicit none

    ! MPI_RECV_INIT requires intent(inout)
    real(RP),               intent(inout) :: DENS(KA,IA,JA)
    real(RP),               intent(inout) :: MOMZ(KA,IA,JA)
    real(RP),               intent(inout) :: MOMX(KA,IA,JA)
    real(RP),               intent(inout) :: MOMY(KA,IA,JA)
    real(RP),               intent(inout) :: RHOT(KA,IA,JA)
    real(RP),               intent(inout) :: QTRC(KA,IA,JA,QA)
    real(RP),               intent(inout) :: PROG(KA,IA,JA,VA)
    real(RP),               intent(inout) :: mflx_hi(KA,IA,JA,3)

    integer :: iv, iq
    !---------------------------------------------------------------------------

    allocate( DENS_t(KA,IA,JA) )
    allocate( MOMZ_t(KA,IA,JA) )
    allocate( MOMX_t(KA,IA,JA) )
    allocate( MOMY_t(KA,IA,JA) )
    allocate( RHOT_t(KA,IA,JA) )
    allocate( RHOQ_t(KA,IA,JA,QA) )

    allocate( I_COMM_PROG    (max(VA,1)) )
    allocate( I_COMM_QTRC(QA) )
    allocate( I_COMM_RHOQ_t(QA) )

    call COMM_vars8_init( 'DENS', DENS, I_COMM_DENS )
    call COMM_vars8_init( 'MOMZ', MOMZ, I_COMM_MOMZ )
    call COMM_vars8_init( 'MOMX', MOMX, I_COMM_MOMX )
    call COMM_vars8_init( 'MOMY', MOMY, I_COMM_MOMY )
    call COMM_vars8_init( 'RHOT', RHOT, I_COMM_RHOT )
    do iv = 1, VA
       I_COMM_PROG(iv) = 5 + iv
       call COMM_vars8_init( 'PROG', PROG(:,:,:,iv), I_COMM_PROG(iv) )
    end do

    call COMM_vars8_init( 'DENS_t', DENS_t, I_COMM_DENS_t )
    call COMM_vars8_init( 'MOMZ_t', MOMZ_t, I_COMM_MOMZ_t )
    call COMM_vars8_init( 'MOMX_t', MOMX_t, I_COMM_MOMX_t )
    call COMM_vars8_init( 'MOMY_t', MOMY_t, I_COMM_MOMY_t )
    call COMM_vars8_init( 'RHOT_t', RHOT_t, I_COMM_RHOT_t )

    do iq = 1, QA
       I_COMM_RHOQ_t(iq) = 5 + VA + iq
       I_COMM_QTRC(iq) = 5 + VA + iq

       call COMM_vars8_init( 'RHOQ_t', RHOQ_t(:,:,:,iq), I_COMM_RHOQ_t(iq) )
       call COMM_vars8_init( 'QTRC',   QTRC  (:,:,:,iq), I_COMM_QTRC(iq) )
    end do

    call COMM_vars8_init( 'mflx_Z', mflx_hi(:,:,:,ZDIR), I_COMM_mflx_z )
    call COMM_vars8_init( 'mflx_X', mflx_hi(:,:,:,XDIR), I_COMM_mflx_x )
    call COMM_vars8_init( 'mflx_Y', mflx_hi(:,:,:,YDIR), I_COMM_mflx_y )

    mflx_hi(:,:,:,:) = UNDEF

    return
  end subroutine ATMOS_DYN_Tstep_large_fvm_heve_setup

  !-----------------------------------------------------------------------------
  !> Dynamical Process
  subroutine ATMOS_DYN_Tstep_large_fvm_heve( &
       DENS, MOMZ, MOMX, MOMY, RHOT, QTRC, PROG,             &
       DENS_av, MOMZ_av, MOMX_av, MOMY_av, RHOT_av, QTRC_av, &
       mflx_hi, tflx_hi,                                     &
       num_diff, num_diff_q,                                 &
       QTRC0,                                                &
       DENS_tp, MOMZ_tp, MOMX_tp, MOMY_tp, RHOT_tp, RHOQ_tp, &
       CORIOLI,                                              &
       CDZ, CDX, CDY, FDZ, FDX, FDY,                         &
       RCDZ, RCDX, RCDY, RFDZ, RFDX, RFDY,                   &
       PHI, GSQRT,                                           &
       J13G, J23G, J33G, MAPF,                               &
       AQ_R, AQ_CV, AQ_CP, AQ_MASS,                          &
       REF_dens, REF_pott, REF_qv, REF_pres,                 &
       BND_W, BND_E, BND_S, BND_N,                           &
       ND_COEF, ND_COEF_Q, ND_ORDER, ND_SFC_FACT, ND_USE_RS, &
       DAMP_DENS,       DAMP_VELZ,       DAMP_VELX,          &
       DAMP_VELY,       DAMP_POTT,       DAMP_QTRC,          &
       DAMP_alpha_DENS, DAMP_alpha_VELZ, DAMP_alpha_VELX,    &
       DAMP_alpha_VELY, DAMP_alpha_POTT, DAMP_alpha_QTRC,    &
       wdamp_coef,                                           &
       divdmp_coef,                                          &
       FLAG_FCT_MOMENTUM, FLAG_FCT_T, FLAG_FCT_TRACER,       &
       FLAG_FCT_ALONG_STREAM,                                &
       USE_AVERAGE,                                          &
       DTLS, DTSS, Llast                                     )
    use scale_const, only: &
       EPS    => CONST_EPS, &
       P0     => CONST_PRE00, &
       Rdry   => CONST_Rdry, &
       CVdry  => CONST_CVdry, &
       CPdry  => CONST_CPdry
    use scale_comm, only: &
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
    use scale_atmos_dyn_common, only: &
       ATMOS_DYN_divergence, &
       ATMOS_DYN_numfilter_coef,   &
       ATMOS_DYN_numfilter_coef_q, &
       ATMOS_DYN_fct
    use scale_atmos_dyn_fvm_flux_ud1, only: &
       ATMOS_DYN_FVM_fluxZ_XYZ_ud1, &
       ATMOS_DYN_FVM_fluxX_XYZ_ud1, &
       ATMOS_DYN_FVM_fluxY_XYZ_ud1
    use scale_atmos_dyn_fvm_flux, only: &
       ATMOS_DYN_FVM_fluxZ_XYZ, &
       ATMOS_DYN_FVM_fluxX_XYZ, &
       ATMOS_DYN_FVM_fluxY_XYZ
    use scale_atmos_boundary, only: &
       BND_QA, &
       BND_SMOOTHER_FACT => ATMOS_BOUNDARY_SMOOTHER_FACT
    use scale_history, only: &
#ifdef HIST_TEND
       HIST_in, &
#endif
       HIST_switch
    use scale_atmos_dyn_tinteg_short, only: &
       ATMOS_DYN_tinteg_short
    use scale_atmos_dyn_tinteg_tracer, only: &
       ATMOS_DYN_tinteg_tracer
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

    real(RP), intent(out)   :: mflx_hi(KA,IA,JA,3)
    real(RP), intent(out)   :: tflx_hi(KA,IA,JA,3)

    real(RP), intent(out)   :: num_diff(KA,IA,JA,5,3)
    real(RP), intent(out)   :: num_diff_q(KA,IA,JA,3)

    real(RP), intent(in)    :: QTRC0(KA,IA,JA,QA)

    real(RP), intent(in)    :: DENS_tp(KA,IA,JA)
    real(RP), intent(in)    :: MOMZ_tp(KA,IA,JA)
    real(RP), intent(in)    :: MOMX_tp(KA,IA,JA)
    real(RP), intent(in)    :: MOMY_tp(KA,IA,JA)
    real(RP), intent(in)    :: RHOT_tp(KA,IA,JA)
    real(RP), intent(in)    :: RHOQ_tp(KA,IA,JA,QA)

    real(RP), intent(in)    :: CORIOLI(IA,JA)

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

    real(RP), intent(in)    :: AQ_R   (QA)
    real(RP), intent(in)    :: AQ_CV  (QA)
    real(RP), intent(in)    :: AQ_CP  (QA)
    real(RP), intent(in)    :: AQ_MASS(QA)

    real(RP), intent(in)    :: REF_dens(KA,IA,JA)
    real(RP), intent(in)    :: REF_pott(KA,IA,JA)
    real(RP), intent(in)    :: REF_qv  (KA,IA,JA)
    real(RP), intent(in)    :: REF_pres(KA,IA,JA)   !< reference pressure

    logical,  intent(in)    :: BND_W
    logical,  intent(in)    :: BND_E
    logical,  intent(in)    :: BND_S
    logical,  intent(in)    :: BND_N

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

    real(RP), intent(in)    :: wdamp_coef(KA)
    real(RP), intent(in)    :: divdmp_coef

    logical,  intent(in)    :: FLAG_FCT_MOMENTUM
    logical,  intent(in)    :: FLAG_FCT_T
    logical,  intent(in)    :: FLAG_FCT_TRACER
    logical,  intent(in)    :: FLAG_FCT_ALONG_STREAM

    logical,  intent(in)    :: USE_AVERAGE

    real(DP), intent(in)    :: DTLS
    real(DP), intent(in)    :: DTSS
    logical , intent(in)    :: Llast

    ! for time integartion
    real(RP) :: DENS00  (KA,IA,JA) ! saved density before small step loop

    ! diagnostic variables
    real(RP) :: DDIV    (KA,IA,JA) ! 3 dimensional divergence
    real(RP) :: DPRES0  (KA,IA,JA) ! pressure deviation
    real(RP) :: RT2P    (KA,IA,JA) ! factor of RHOT to PRES
    real(RP) :: REF_rhot(KA,IA,JA) ! reference of RHOT

    real(RP) :: QDRY  ! dry air
    real(RP) :: Rtot  ! total R
    real(RP) :: CVtot ! total CV
    real(RP) :: CPtot ! total CP
    real(RP) :: PRES  ! pressure

    real(RP) :: DENS_tq(KA,IA,JA)
    real(RP) :: diff(KA,IA,JA)
    real(RP) :: damp
#ifdef HIST_TEND
    real(RP) :: damp_t(KA,IA,JA)
#endif

    ! For tracer advection
    real(RP) :: mflx_av  (KA,IA,JA,3)  ! rho * vel(x,y,z) @ (u,v,w)-face average
    real(RP) :: qflx_hi  (KA,IA,JA,3)  ! rho * vel(x,y,z) * phi @ (u,v,w)-face high order
    real(RP) :: qflx_lo  (KA,IA,JA,3)  ! rho * vel(x,y,z) * phi,  monotone flux
    real(RP) :: qflx_anti(KA,IA,JA,3)  ! anti-diffusive flux

    real(RP) :: dtl
    real(RP) :: dts
    integer  :: nstep

    integer  :: IIS, IIE
    integer  :: JJS, JJE
    integer  :: i, j, k, iq, step
    integer  :: iv

    real(RP) :: diff_coef
    !---------------------------------------------------------------------------

    call PROF_rapstart("DYN_Large_Preparation", 2)

    dts   = real(DTSS, kind=RP)            ! short time step
    dtl   = real(DTLS, kind=RP)            ! large time step
    nstep = ceiling( ( dtl - eps ) / dts )
    dts   = dtl / nstep                    ! dts is divisor of dtl and smaller or equal to dtss

#ifdef DEBUG
    if( IO_L ) write(IO_FID_LOG,*)                         '*** Dynamics large time step'
    if( IO_L ) write(IO_FID_LOG,'(1x,A,F0.2,A,F0.2,A,I0)') &
    '*** -> DT_large, DT_small, DT_large/DT_small : ', dtl, ', ', dts, ', ', nstep

    DENS00  (:,:,:) = UNDEF

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

#endif

    !------------------------------------------------------------------------
    ! prepare thermodynamical data
    !   specific heat
    !   pressure data ( linearization )
    !
    ! pres = P0 * ( R * rhot / P0 )**(CP/CV)
    ! d pres / d rhot ~ CP*R/CV * ( R * rhot / P0 )**(R/CV)
    !                 = CP*R/CV * ( pres / P0 )**(R/CP)
    !                 = CP*R/CV * temp / pott
    !                 = CP/CV * pres / rhot
    ! pres ~ P0 * ( R * rhot0 / P0 ) ** (CP/CV) + CV*R/CP * ( pres / P0 )**(R/CP) * rhot'
    !------------------------------------------------------------------------
!OCL XFILL
    !$omp parallel do default(none) &
    !$omp shared(JA,IA,KS,KE,P0,Rdry,RHOT,AQ_R,AQ_CV,AQ_CP,QTRC,AQ_MASS,REF_rhot,REF_pres,CPdry,CVdry,QA,RT2P,DPRES0) &
    !$omp private(i,j,k,iq,PRES,Rtot,CVtot,CPtot,QDRY) OMP_SCHEDULE_ collapse(2)
    do j = 1, JA
    do i = 1, IA
       do k = KS, KE
#ifdef DRY
          PRES = P0 * ( Rdry * RHOT(k,i,j) / P0 )**CPovCV
          RT2P(k,i,j) = CPovCV * PRES / RHOT(k,i,j)
#else
          Rtot  = 0.0_RP
          CVtot = 0.0_RP
          CPtot = 0.0_RP
          QDRY  = 1.0_RP
          do iq = 1, QA
             Rtot  = Rtot + AQ_R(iq) * QTRC(k,i,j,iq)
             CVtot = CVtot + AQ_CV(iq) * QTRC(k,i,j,iq)
             CPtot = CPtot + AQ_CP(iq) * QTRC(k,i,j,iq)
             QDRY  = QDRY  - QTRC(k,i,j,iq) * AQ_MASS(iq)
          enddo
          Rtot  = Rtot  + Rdry  * QDRY
          CVtot = CVtot + CVdry * QDRY
          CPtot = CPtot + CPdry * QDRY
          PRES = P0 * ( Rtot * RHOT(k,i,j) / P0 )**( CPtot / CVtot )
          RT2P(k,i,j) = CPtot / CVtot * PRES / RHOT(k,i,j)
#endif
          DPRES0(k,i,j) = PRES - REF_pres(k,i,j)
          REF_rhot(k,i,j) = RHOT(k,i,j)
       end do
       DPRES0(KS-1,i,j) = DPRES0(KS+1,i,j) + ( REF_pres(KS+1,i,j) - REF_pres(KS-1,i,j) )
       DPRES0(KE+1,i,j) = DPRES0(KE-1,i,j) + ( REF_pres(KE-1,i,j) - REF_pres(KE+1,i,j) )
    end do
    end do

    call PROF_rapend  ("DYN_Large_Preparation", 2)

    !###########################################################################
    ! Update DENS,MONZ,MOMX,MOMY,MOMZ,RHOT
    !###########################################################################

    call PROF_rapstart("DYN_Large_Tendency", 2)

#ifdef HIST_TEND
!OCL XFILL
    damp_t(:,:,:) = 0.0_RP
#endif

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
       !$omp parallel do default(none) private(i,j,k) OMP_SCHEDULE_ collapse(2) &
       !$omp private(damp) &
       !$omp shared(JS,JE,IS,IE,KS,KE,DAMP_alpha_QTRC,iq,diff,BND_SMOOTHER_FACT,DENS00) &
       !$omp shared(RHOQ_t,RHOQ_tp,DENS_tq)
!OCL XFILL
       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          damp = - DAMP_alpha_QTRC(k,i,j,iq) &
               * ( diff(k,i,j) & ! rayleigh damping
                 - ( diff(k,i-1,j) + diff(k,i+1,j) + diff(k,i,j-1) + diff(k,i,j+1) - diff(k,i,j)*4.0_RP ) &
                   * 0.125_RP * BND_SMOOTHER_FACT ) ! horizontal smoother
#ifdef HIST_TEND
          damp_t(k,i,j) = damp
#endif
          damp = damp * DENS00(k,i,j)
          RHOQ_t(k,i,j,iq) = RHOQ_tp(k,i,j,iq) + damp
          DENS_tq(k,i,j) = DENS_tq(k,i,j) + damp
       enddo
       enddo
       enddo
#ifdef HIST_TEND
       call HIST_in(RHOQ_tp(:,:,:,iq), trim(TRACER_NAME(iq))//'_t_phys', &
                    'tendency of '//trim(TRACER_NAME(iq))//' due to physics', 'kg/kg/s' )
       call HIST_in(damp_t,            trim(TRACER_NAME(iq))//'_t_damp', &
                    'tendency of '//trim(TRACER_NAME(iq))//' due to damping', 'kg/kg/s' )
#endif
!OCL XFILL
       do j = JS, JE
       do i = IS, IE
          RHOQ_t(   1:KS-1,i,j,iq) = 0.0_RP
          RHOQ_t(KE+1:KA  ,i,j,iq) = 0.0_RP
       enddo
       enddo

       call COMM_vars8( RHOQ_t(:,:,:,iq), I_COMM_RHOQ_t(iq) )
       call COMM_wait ( RHOQ_t(:,:,:,iq), I_COMM_RHOQ_t(iq), .false. )

    end do

    !$omp parallel do private(i,j,k,iq) OMP_SCHEDULE_ collapse(3)
!OCL XFILL
    do iq = BND_QA+1, QA
       do j = 1, JA
       do i = 1, IA
       do k = 1, KA
          RHOQ_t(k,i,j,iq) = RHOQ_tp(k,i,j,iq)
       enddo
       enddo
       enddo
    end do

    call PROF_rapend  ("DYN_Large_Tendency", 2)

    call PROF_rapstart("DYN_Large_Boundary", 2)

    if ( BND_W ) then
       !$omp parallel do private(j,k) OMP_SCHEDULE_ collapse(2)
       do j = JS, JE
       do k = KS, KE
          mflx_hi(k,IS-1,j,XDIR) = GSQRT(k,IS-1,j,I_UYZ) * MOMX(k,IS-1,j) / MAPF(IS-1,j,2,I_UY)
       enddo
       enddo
    end if
    if ( BND_E ) then
       !$omp parallel do private(j,k) OMP_SCHEDULE_ collapse(2)
       do j = JS, JE
       do k = KS, KE
          mflx_hi(k,IE,j,XDIR) = GSQRT(k,IE,j,I_UYZ) * MOMX(k,IE,j) / MAPF(IE,j,2,I_UY)
       enddo
       enddo
    end if
    if ( BND_S ) then
       !$omp parallel do private(i,k) OMP_SCHEDULE_ collapse(2)
       do i = IS, IE
       do k = KS, KE
          mflx_hi(k,i,JS-1,YDIR) = GSQRT(k,i,JS-1,I_XVZ) * MOMY(k,i,JS-1) / MAPF(i,JS-1,1,I_XV)
       enddo
       enddo
    end if
    if ( BND_N ) then
       !$omp parallel do private(i,k) OMP_SCHEDULE_ collapse(2)
       do i = IS, IE
       do k = KS, KE
          mflx_hi(k,i,JE,YDIR) = GSQRT(k,i,JE,I_XVZ) * MOMY(k,i,JE) / MAPF(i,JE,1,I_XV)
       enddo
       enddo
    end if

    call PROF_rapend  ("DYN_Large_Boundary", 2)


    do step = 1, nstep

       call HIST_switch( Llast .AND. step == nstep )

       !-----< prepare tendency >-----

       call PROF_rapstart("DYN_Large_Tendency", 2)

       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
!OCL XFILL
       do j = JS-1, JE+1
       do i = IS-1, IE+1
       do k = KS, KE
          diff(k,i,j) = DENS(k,i,j) - DAMP_DENS(k,i,j)
       enddo
       enddo
       enddo
       !$omp parallel do default(none) private(i,j,k) OMP_SCHEDULE_ collapse(2) &
       !$omp private(damp) &
       !$omp shared(JS,JE,IS,IE,KS,KE,DAMP_alpha_DENS,diff,DENS_tq,DENS_t,DENS_tp,BND_SMOOTHER_FACT)
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
       call HIST_in(DENS_tp, 'DENS_t_phys', 'tendency of dencity due to physics', 'kg/m3/s' )
       call HIST_in(damp_t,  'DENS_t_damp', 'tendency of dencity due to damping', 'kg/m3/s' )
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
       !$omp parallel do default(none) private(i,j,k) OMP_SCHEDULE_ collapse(2) &
       !$omp private(damp) &
       !$omp shared(JS,JE,IS,IE,KS,KE,DAMP_alpha_VELZ,diff,BND_SMOOTHER_FACT,MOMZ_t,MOMZ_tp)
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
       call HIST_in(MOMZ_tp, 'MOMZ_t_phys', 'tendency of momentum z due to physics',  'kg/m2/s2', zdim='half' )
       call HIST_in(damp_t,  'MOMZ_t_damp', 'tendency of momentum z due to damping', 'kg/m2/s2', zdim='half' )
#endif

       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
!OCL XFILL
       do j = JS-1, JE+1
       do i = IS-1, IE+1
       do k = KS, KE
          diff(k,i,j) = MOMX(k,i,j) - DAMP_VELX(k,i,j) * ( DENS(k,i,j)+DENS(k,i+1,j) ) * 0.5_RP
       enddo
       enddo
       enddo
!OCL XFILL
       !$omp parallel do default(none) &
       !$omp shared(JS,JE,IS,IE,KS,KE,DAMP_alpha_VELX,diff,BND_SMOOTHER_FACT,MOMX_tp,MOMX_t) &
       !$omp private(i,j,k,damp) OMP_SCHEDULE_ collapse(2)
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
       call HIST_in(MOMX_tp, 'MOMX_t_phys', 'tendency of momentum x due to physics', 'kg/m2/s2', xdim='half' )
       call HIST_in(damp_t,  'MOMX_t_damp', 'tendency of momentum x due to damping', 'kg/m2/s2', xdim='half' )
#endif

       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
!OCL XFILL
       do j = JS-1, JE+1
       do i = IS-1, IE+1
       do k = KS, KE
          diff(k,i,j) = MOMY(k,i,j) - DAMP_VELY(k,i,j) * ( DENS(k,i,j)+DENS(k,i,j+1) ) * 0.5_RP
       enddo
       enddo
       enddo
!OCL XFILL
       !$omp parallel do default(none) &
       !$omp shared(JS,JE,IS,IE,KS,KE,DAMP_alpha_VELY,diff,BND_SMOOTHER_FACT,MOMY_tp,MOMY_t) &
       !$omp private(i,j,k,damp) OMP_SCHEDULE_ collapse(2)
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
       call HIST_in(MOMY_tp, 'MOMY_t_phys', 'tendency of momentum y due to physics', 'kg/m2/s2', ydim='half' )
       call HIST_in(damp_t,  'MOMY_t_damp', 'tendency of momentum y due to damping', 'kg/m2/s2', ydim='half' )
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
       !$omp parallel do default(none) &
       !$omp shared(JS,JE,IS,IE,KS,KE,DAMP_alpha_POTT,diff,BND_SMOOTHER_FACT,RHOT_t,RHOT_tp) &
       !$omp private(i,j,k,damp) OMP_SCHEDULE_ collapse(2)
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
       call HIST_in(RHOT_tp, 'RHOT_t_phys', 'tendency of rho*theta temperature due to physics', 'K kg/m3/s' )
       call HIST_in(damp_t,  'RHOT_t_damp', 'tendency of rho*theta temperature due to damping', 'K kg/m3/s' )
#endif

       call COMM_wait ( DENS_t(:,:,:), I_COMM_DENS_t, .false. )
       call COMM_wait ( MOMZ_t(:,:,:), I_COMM_MOMZ_t, .false. )
       call COMM_wait ( MOMX_t(:,:,:), I_COMM_MOMX_t, .false. )
       call COMM_wait ( MOMY_t(:,:,:), I_COMM_MOMY_t, .false. )
       call COMM_wait ( RHOT_t(:,:,:), I_COMM_RHOT_t, .false. )

       call PROF_rapend  ("DYN_Large_Tendency", 2)

       call PROF_rapstart("DYN_Large_Numfilter", 2)

       !-----< prepare numerical diffusion coefficient >-----

       if ( ND_COEF == 0.0_RP ) then
!OCL XFILL
          num_diff(:,:,:,:,:) = 0.0_RP
       else
          call ATMOS_DYN_numfilter_coef( num_diff(:,:,:,:,:),                    & ! [OUT]
                                         DENS, MOMZ, MOMX, MOMY, RHOT,           & ! [IN]
                                         CDZ, CDX, CDY, FDZ, FDX, FDY, dts,      & ! [IN]
                                         REF_dens, REF_pott,                     & ! [IN]
                                         ND_COEF, ND_ORDER, ND_SFC_FACT, ND_USE_RS ) ! [IN]
       endif

       call PROF_rapend  ("DYN_Large_Numfilter", 2)

       if ( divdmp_coef > 0.0_RP ) then

          call ATMOS_DYN_divergence( DDIV, & ! (out)
               MOMZ, MOMX, MOMY, & ! (in)
               GSQRT, J13G, J23G, J33G, MAPF, & ! (in)
               RCDZ, RCDX, RCDY, RFDZ, FDZ ) ! (in)

       else

!XFILL
          DDIV = 0.0_RP

       end if

       !------------------------------------------------------------------------
       ! Start short time integration
       !------------------------------------------------------------------------

       call PROF_rapstart("DYN_Short_Tinteg", 2)

       call ATMOS_DYN_tinteg_short( DENS, MOMZ, MOMX, MOMY, RHOT, PROG,       & ! (inout)
                                    mflx_hi, tflx_hi,                         & ! (inout)
                                    DENS_t, MOMZ_t, MOMX_t, MOMY_t, RHOT_t,   & ! (in)
                                    DPRES0, RT2P, CORIOLI,                    & ! (in)
                                    num_diff, wdamp_coef, divdmp_coef, DDIV,  & ! (in)
                                    FLAG_FCT_MOMENTUM, FLAG_FCT_T,            & ! (in)
                                    FLAG_FCT_ALONG_STREAM,                    & ! (in)
                                    CDZ, FDZ, FDX, FDY,                       & ! (in)
                                    RCDZ, RCDX, RCDY, RFDZ, RFDX, RFDY,       & ! (in)
                                    PHI, GSQRT, J13G, J23G, J33G, MAPF,       & ! (in)
                                    REF_dens, REF_rhot,                       & ! (in)
                                    BND_W, BND_E, BND_S, BND_N,               & ! (in)
                                    dts                                       ) ! (in)

       call PROF_rapend  ("DYN_Short_Tinteg", 2)

#ifdef CHECK_MASS
       call check_mass( &
            DENS, DAMP_DENS, &
            mflx_hi, tflx_hi, &
            GSQRT, MAPF, &
            RCDX, RCDY, &
            dts, step, &
            BND_W, BND_E, BND_S, BND_N )
#endif

       !$omp parallel do default(none) private(i,j,iv) OMP_SCHEDULE_ collapse(2) &
       !$omp shared(JS,JE,IS,IE,KS,KA,DENS,MOMZ,MOMX,MOMY,RHOT,VA,PROG,KE)
       do j  = JS, JE
       do i  = IS, IE
          DENS(   1:KS-1,i,j) = DENS(KS,i,j)
          MOMZ(   1:KS-1,i,j) = MOMZ(KS,i,j)
          MOMX(   1:KS-1,i,j) = MOMX(KS,i,j)
          MOMY(   1:KS-1,i,j) = MOMY(KS,i,j)
          RHOT(   1:KS-1,i,j) = RHOT(KS,i,j)
          do iv = 1, VA
             PROG(   1:KS-1,i,j,iv) = PROG(KS,i,j,iv)
          end do
          DENS(KE+1:KA,  i,j) = DENS(KE,i,j)
          MOMZ(KE+1:KA,  i,j) = MOMZ(KE,i,j)
          MOMX(KE+1:KA,  i,j) = MOMX(KE,i,j)
          MOMY(KE+1:KA,  i,j) = MOMY(KE,i,j)
          RHOT(KE+1:KA,  i,j) = RHOT(KE,i,j)
          do iv = 1, VA
             PROG(KE+1:KA,  i,j,iv) = PROG(KE,i,j,iv)
          end do
       enddo
       enddo

       call COMM_vars8( DENS(:,:,:), I_COMM_DENS )
       call COMM_vars8( MOMZ(:,:,:), I_COMM_MOMZ )
       call COMM_vars8( MOMX(:,:,:), I_COMM_MOMX )
       call COMM_vars8( MOMY(:,:,:), I_COMM_MOMY )
       call COMM_vars8( RHOT(:,:,:), I_COMM_RHOT )
       do iv = 1, VA
          call COMM_vars8( PROG(:,:,:,iv), I_COMM_PROG(iv) )
       end do
       call COMM_wait ( DENS(:,:,:), I_COMM_DENS, .false. )
       call COMM_wait ( MOMZ(:,:,:), I_COMM_MOMZ, .false. )
       call COMM_wait ( MOMX(:,:,:), I_COMM_MOMX, .false. )
       call COMM_wait ( MOMY(:,:,:), I_COMM_MOMY, .false. )
       call COMM_wait ( RHOT(:,:,:), I_COMM_RHOT, .false. )
       do iv = 1, VA
          call COMM_wait ( PROG(:,:,:,iv), I_COMM_PROG(iv), .false. )
       end do

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
       DENS_av(:,:,:) = DENS_av(:,:,:) / nstep
       MOMZ_av(:,:,:) = MOMZ_av(:,:,:) / nstep
       MOMX_av(:,:,:) = MOMX_av(:,:,:) / nstep
       MOMY_av(:,:,:) = MOMY_av(:,:,:) / nstep
       RHOT_av(:,:,:) = RHOT_av(:,:,:) / nstep
    endif

#ifndef DRY
    !###########################################################################
    ! Update Tracers
    !###########################################################################

!OCL XFILL
    mflx_hi(:,:,:,:) = mflx_av(:,:,:,:) / nstep

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

#ifdef _SDM
    do iq = 1, I_QV
#else
    do iq = 1, QA
#endif

       if ( TRACER_ADVC(iq) ) then

          call PROF_rapstart("DYN_Large_Numfilter", 2)

          if ( ND_COEF_Q == 0.0_RP ) then
!OCL XFILL
             num_diff_q(:,:,:,:) = 0.0_RP
          else
             call ATMOS_DYN_numfilter_coef_q( num_diff_q(:,:,:,:),                    & ! [OUT]
                                              DENS00, QTRC(:,:,:,iq),                 & ! [IN]
                                              CDZ, CDX, CDY, dtl,                     & ! [IN]
                                              REF_qv, iq,                             & ! [IN]
                                              ND_COEF_Q, ND_ORDER, ND_SFC_FACT, ND_USE_RS ) ! [IN]
          endif

          call PROF_rapend  ("DYN_Large_Numfilter", 2)

          call PROF_rapstart("DYN_Tracer_Tinteg", 2)

          call ATMOS_DYN_tinteg_tracer( &
               QTRC(:,:,:,iq), & ! (inout)
               QTRC0(:,:,:,iq), RHOQ_t(:,:,:,iq), &! (in)
               DENS00, DENS, & ! (in)
               mflx_hi, num_diff_q, & ! (in)
               GSQRT, MAPF(:,:,:,I_XY), & ! (in)
               CDZ, RCDZ, RCDX, RCDY,             & ! (in)
               dtl, & ! (in)
               Llast .AND. FLAG_FCT_TRACER, & ! (in)
               FLAG_FCT_ALONG_STREAM         ) ! (in)

          call PROF_rapend  ("DYN_Tracer_Tinteg", 2)

       else

          do j = JS, JE
          do i = IS, IE
          do k = KS, KE
             QTRC(k,i,j,iq) = ( QTRC0(k,i,j,iq) * DENS00(k,i,j) &
                              + RHOQ_t(k,i,j,iq) * dtl          ) / DENS(k,i,j)
          end do
          end do
          end do

       end if

       if ( USE_AVERAGE ) then
          QTRC_av(:,:,:,iq) = QTRC(:,:,:,iq)
       endif

       call COMM_vars8( QTRC(:,:,:,iq), I_COMM_QTRC(iq) )

    enddo ! scalar quantities loop

    do iq = 1, QA
       call COMM_wait ( QTRC(:,:,:,iq), I_COMM_QTRC(iq), .false. )
    enddo
#endif

    return
  end subroutine ATMOS_DYN_Tstep_large_fvm_heve

#ifdef CHECK_MASS
  subroutine check_mass( &
       DENS, DAMP_DENS, &
       mflx_hi, tflx_hi, &
       GSQRT, MAPF, &
       RCDX, RCDY, &
       dt, step, &
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
    integer,  intent(in) :: step
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

    if( IO_L ) write(IO_FID_LOG,'(A,1x,I1,1x,ES24.17)') 'total mflx_lb:', step, allmflx_lb_total

    call MPI_Allreduce( mass_total,           &
                        allmass_total,        &
                        1,                    &
                        COMM_datatype,        &
                        MPI_SUM,              &
                        COMM_world,           &
                        ierr                  )

    if( IO_L ) write(IO_FID_LOG,'(A,1x,I1,1x,ES24.17)') 'total mass   :', step, allmass_total

    call MPI_Allreduce( mass_total2,          &
                        allmass_total2,       &
                        1,                    &
                        COMM_datatype,        &
                        MPI_SUM,              &
                        COMM_world,           &
                        ierr                  )

    if( IO_L ) write(IO_FID_LOG,'(A,1x,I1,1x,ES24.17)') 'total mass2  :', step, allmass_total2

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

end module scale_atmos_dyn_tstep_large_fvm_heve
