!-------------------------------------------------------------------------------
!> module Atmosphere / Dynamics FENT + FCT
!!
!! @par Description
!!          Dynamical core for Atmospheric process
!!          Full explicit, no terrain + tracer FCT limiter
!!
!! @author Team SCALE
!!
!<
!-------------------------------------------------------------------------------
#include "scalelib.h"
module scale_atmos_dyn
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
  ! for communication
  integer,  private              :: I_COMM_DENS = 1
  integer,  private              :: I_COMM_MOMZ = 2
  integer,  private              :: I_COMM_MOMX = 3
  integer,  private              :: I_COMM_MOMY = 4
  integer,  private              :: I_COMM_RHOT = 5
  integer,  private, allocatable :: I_COMM_PROG(:)
  integer,  private, allocatable :: I_COMM_QTRC(:)

  logical,  private              :: BND_W
  logical,  private              :: BND_E
  logical,  private              :: BND_S
  logical,  private              :: BND_N

  real(RP), private, allocatable :: num_diff  (:,:,:,:,:)
  real(RP), private, allocatable :: num_diff_q(:,:,:,:)
  real(RP), private, allocatable :: wdamp_coef(:)         ! coefficient for Rayleigh damping of w

  logical,  private              :: DYN_NONE

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine ATMOS_DYN_setup( &
       DYN_Tinteg_Short_TYPE,        &
       DYN_Tinteg_Tracer_TYPE,       &
       DYN_Tinteg_Large_TYPE,        &
       DYN_Tstep_Tracer_TYPE,        &
       DYN_Tstep_Large_TYPE,         &
       DYN_Tstep_Short_TYPE,         &
       DYN_FVM_FLUX_TYPE,            &
       DYN_FVM_FLUX_TYPE_TRACER,     &
       DENS, MOMZ, MOMX, MOMY, RHOT, &
       QTRC, PROG,                   &
       CDZ, CDX, CDY,                &
       FDZ, FDX, FDY,                &
       wdamp_tau,                    &
       wdamp_height,                 &
       FZ,                           &
       none                          )
    use scale_prc, only: &
       PRC_abort
    use scale_prc_cartesC, only: &
       PRC_TwoD,  &
       PRC_HAS_E, &
       PRC_HAS_W, &
       PRC_HAS_N, &
       PRC_HAS_S
    use scale_const, only: &
       UNDEF => CONST_UNDEF
    use scale_comm_cartesC, only: &
       COMM_vars8_init
    use scale_atmos_dyn_common, only: &
       ATMOS_DYN_filter_setup, &
       ATMOS_DYN_wdamp_setup
    use scale_atmos_dyn_tinteg_short, only: &
       ATMOS_DYN_tinteg_short_setup
    use scale_atmos_dyn_tinteg_tracer, only: &
       ATMOS_DYN_tinteg_tracer_setup
    use scale_atmos_dyn_tinteg_large, only: &
       ATMOS_DYN_tinteg_large_setup
    use scale_atmos_dyn_tstep_short, only: &
       ATMOS_DYN_tstep_short_setup
    use scale_atmos_dyn_tstep_tracer, only: &
       ATMOS_DYN_tstep_tracer_setup
    use scale_atmos_dyn_tstep_large, only: &
       ATMOS_DYN_tstep_large_setup
    use scale_atmos_dyn_fvm_flux, only: &
       ATMOS_DYN_FVM_flux_setup
    use scale_spnudge, only: &
       SPNUDGE_setup
    implicit none

    character(len=*),  intent(in)    :: DYN_Tinteg_Short_TYPE
    character(len=*),  intent(in)    :: DYN_Tinteg_Tracer_TYPE
    character(len=*),  intent(in)    :: DYN_Tinteg_Large_TYPE
    character(len=*),  intent(in)    :: DYN_Tstep_Tracer_TYPE
    character(len=*),  intent(in)    :: DYN_Tstep_Large_TYPE
    character(len=*),  intent(in)    :: DYN_Tstep_Short_TYPE
    character(len=*),  intent(in)    :: DYN_FVM_FLUX_TYPE
    character(len=*),  intent(in)    :: DYN_FVM_FLUX_TYPE_TRACER
    real(RP),          intent(inout) :: DENS(KA,IA,JA)    ! MPI_RECV_INIT requires intent(inout)
    real(RP),          intent(inout) :: MOMZ(KA,IA,JA)
    real(RP),          intent(inout) :: MOMX(KA,IA,JA)
    real(RP),          intent(inout) :: MOMY(KA,IA,JA)
    real(RP),          intent(inout) :: RHOT(KA,IA,JA)
    real(RP),          intent(inout) :: QTRC(KA,IA,JA,QA)
    real(RP),          intent(inout) :: PROG(KA,IA,JA,VA)
    real(RP),          intent(in)    :: CDZ(KA)
    real(RP),          intent(in)    :: CDX(IA)
    real(RP),          intent(in)    :: CDY(JA)
    real(RP),          intent(in)    :: FDZ(KA-1)
    real(RP),          intent(in)    :: FDX(IA-1)
    real(RP),          intent(in)    :: FDY(JA-1)
    real(RP),          intent(in)    :: wdamp_tau
    real(RP),          intent(in)    :: wdamp_height
    real(RP),          intent(in)    :: FZ(0:KA)
    logical, optional, intent(in)    :: none

    integer :: iv, iq
    !---------------------------------------------------------------------------

    DYN_NONE = .false.

    if ( present(none) ) then
       DYN_NONE = none
    endif

    allocate( I_COMM_PROG(max(VA,1)) )
    allocate( I_COMM_QTRC(QA)        )

    if ( .NOT. DYN_NONE ) then
       BND_W = ( .NOT. PRC_HAS_W ) .and. ( .NOT. PRC_TwoD )
       BND_E = ( .NOT. PRC_HAS_E ) .and. ( .NOT. PRC_TwoD )
       BND_S = .NOT. PRC_HAS_S
       BND_N = .NOT. PRC_HAS_N

       allocate( num_diff  (KA,IA,JA,5,3) )
       allocate( num_diff_q(KA,IA,JA,3)   )
       allocate( wdamp_coef(KA)           )
       num_diff  (:,:,:,:,:) = UNDEF
       num_diff_q(:,:,:,:)   = UNDEF

       call ATMOS_DYN_FVM_flux_setup     ( DYN_FVM_FLUX_TYPE,            & ! [IN]
                                           DYN_FVM_FLUX_TYPE_TRACER      ) ! [IN]

       call ATMOS_DYN_tstep_short_setup

       call ATMOS_DYN_tstep_tracer_setup ( DYN_Tstep_Tracer_TYPE         ) ! [IN]

       call ATMOS_DYN_tstep_large_setup  ( DYN_Tstep_Large_TYPE,         & ! [IN]
                                           DENS, MOMZ, MOMX, MOMY, RHOT, & ! [INOUT]
                                           QTRC, PROG                    ) ! [INOUT]

       call ATMOS_DYN_Tinteg_short_setup ( DYN_Tinteg_Short_TYPE, DYN_Tstep_Short_TYPE  ) ! [IN]

       call ATMOS_DYN_Tinteg_tracer_setup( DYN_Tinteg_Tracer_TYPE ) ! [IN]

       call ATMOS_DYN_Tinteg_large_setup ( DYN_Tinteg_Large_TYPE  ) ! [IN]

       ! numerical diffusion
       call ATMOS_DYN_filter_setup( num_diff, num_diff_q,        & ! [INOUT]
                                    CDZ, CDX, CDY, FDZ, FDX, FDY ) ! [IN]

       ! numerical diffusion
       call ATMOS_DYN_wdamp_setup( wdamp_coef(:),           & ! [INOUT]
                                   wdamp_tau, wdamp_height, & ! [IN]
                                   FZ(:)                    ) ! [IN]

       ! setup for spectral nudging
       call SPNUDGE_setup( KA, KS, KE, IA, IS, IE, JA, JS, JE )


    else

       call COMM_vars8_init( 'DENS', DENS, I_COMM_DENS )
       call COMM_vars8_init( 'MOMZ', MOMZ, I_COMM_MOMZ )
       call COMM_vars8_init( 'MOMX', MOMX, I_COMM_MOMX )
       call COMM_vars8_init( 'MOMY', MOMY, I_COMM_MOMY )
       call COMM_vars8_init( 'RHOT', RHOT, I_COMM_RHOT )

       do iv = 1, VA
          I_COMM_PROG(iv) = 5 + iv
          call COMM_vars8_init( 'PROG', PROG(:,:,:,iv), I_COMM_PROG(iv) )
       enddo

       do iq = 1, QA
          I_COMM_QTRC(iq) = 5 + VA + iq
          call COMM_vars8_init( 'QTRC', QTRC(:,:,:,iq), I_COMM_QTRC(iq) )
       enddo

    endif

    return
  end subroutine ATMOS_DYN_setup

  !-----------------------------------------------------------------------------
  !> Dynamical Process
  subroutine ATMOS_DYN( &
       DENS,    MOMZ,    MOMX,    MOMY,    RHOT,    QTRC,    &
       PROG,                                                 &
       DENS_av, MOMZ_av, MOMX_av, MOMY_av, RHOT_av, QTRC_av, &
       DENS_tp, MOMZ_tp, MOMX_tp, MOMY_tp, RHOT_tp, RHOQ_tp, &
       CORIOLIS,                                             &
       CDZ, CDX, CDY, FDZ, FDX, FDY,                         &
       RCDZ, RCDX, RCDY, RFDZ, RFDX, RFDY,                   &
       PHI, GSQRT,                                           &
       J13G, J23G, J33G, MAPF,                               &
       AQ_R, AQ_CV, AQ_CP, AQ_MASS,                          &
       REF_dens, REF_pott, REF_qv, REF_pres,                 &
       ND_COEF, ND_COEF_Q, ND_LAPLACIAN_NUM,                 &
       ND_SFC_FACT, ND_USE_RS,                               &
       BND_QA, BND_IQ, BND_SMOOTHER_FACT,                    &
       DAMP_DENS,       DAMP_VELZ,       DAMP_VELX,          &
       DAMP_VELY,       DAMP_POTT,       DAMP_QTRC,          &
       DAMP_alpha_DENS, DAMP_alpha_VELZ, DAMP_alpha_VELX,    &
       DAMP_alpha_VELY, DAMP_alpha_POTT, DAMP_alpha_QTRC,    &
       MFLUX_OFFSET_X, MFLUX_OFFSET_Y,                       &
       divdmp_coef,                                          &
       FLAG_TRACER_SPLIT_TEND,                               &
       FLAG_FCT_MOMENTUM, FLAG_FCT_T, FLAG_FCT_TRACER,       &
       FLAG_FCT_ALONG_STREAM,                                &
       USE_AVERAGE,                                          &
       I_QV,                                                 &
       DTSEC, DTSEC_DYN                                      )
    use scale_prc_cartesC, only: &
       PRC_TwoD
    use scale_comm_cartesC, only: &
       COMM_vars8, &
       COMM_wait
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

    real(RP), intent(in)    :: CORIOLIS(IA,JA)

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
    real(RP), intent(in)    :: ND_COEF
    real(RP), intent(in)    :: ND_COEF_Q
    integer,  intent(in)    :: ND_LAPLACIAN_NUM
    real(RP), intent(in)    :: ND_SFC_FACT
    logical,  intent(in)    :: ND_USE_RS

    integer,  intent(in)    :: BND_QA
    integer,  intent(in)    :: BND_IQ(QA)
    real(RP), intent(in)    :: BND_SMOOTHER_FACT

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
    real(RP), intent(in)    :: MFLUX_OFFSET_X(KA,JA,2)
    real(RP), intent(in)    :: MFLUX_OFFSET_Y(KA,JA,2)

    real(RP), intent(in)    :: divdmp_coef

    logical,  intent(in)    :: FLAG_TRACER_SPLIT_TEND
    logical,  intent(in)    :: FLAG_FCT_MOMENTUM
    logical,  intent(in)    :: FLAG_FCT_T
    logical,  intent(in)    :: FLAG_FCT_TRACER
    logical,  intent(in)    :: FLAG_FCT_ALONG_STREAM

    logical,  intent(in)    :: USE_AVERAGE

    integer,  intent(in)    :: I_QV

    real(DP), intent(in)    :: DTSEC
    real(DP), intent(in)    :: DTSEC_DYN

    ! for time integartion
    real(RP) :: DENS00 (KA,IA,JA)   ! saved density before update

    ! For tracer advection
    real(RP) :: dt

    integer  :: i, j, k, iq
    !---------------------------------------------------------------------------

    LOG_PROGRESS(*) 'atmosphere / dynamics'

    dt = real(DTSEC, kind=RP)

    if ( DYN_NONE ) then

       call PROF_rapstart("DYN_Tinteg", 2)

       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          DENS00(k,i,j) = DENS(k,i,j) ! save previous value

          DENS  (k,i,j) = DENS(k,i,j) + DENS_tp(k,i,j) * dt
       enddo
       enddo
       enddo

       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          MOMZ(k,i,j) = MOMZ(k,i,j) + MOMZ_tp(k,i,j) * dt
       enddo
       enddo
       enddo

       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          MOMX(k,i,j) = MOMX(k,i,j) + MOMX_tp(k,i,j) * dt
       enddo
       enddo
       enddo

       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          MOMY(k,i,j) = MOMY(k,i,j) + MOMY_tp(k,i,j) * dt
       enddo
       enddo
       enddo

       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          RHOT(k,i,j) = RHOT(k,i,j) + RHOT_tp(k,i,j) * dt
       enddo
       enddo
       enddo

       do iq = 1, QA
       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          QTRC(k,i,j,iq) = max( 0.0_RP, &
                                QTRC(k,i,j,iq) * DENS00(k,i,j) + RHOQ_tp(k,i,j,iq) * dt &
                              ) / DENS(k,i,j)
       enddo
       enddo
       enddo
       enddo

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
          DENS_av(:,:,:)   = DENS(:,:,:)
          MOMZ_av(:,:,:)   = MOMZ(:,:,:)
          MOMX_av(:,:,:)   = MOMX(:,:,:)
          MOMY_av(:,:,:)   = MOMY(:,:,:)
          RHOT_av(:,:,:)   = RHOT(:,:,:)
          QTRC_av(:,:,:,:) = QTRC(:,:,:,:)
       endif

       call PROF_rapend("DYN_Tinteg", 2)

       return
    endif ! if DYN_NONE == .true.

    call PROF_rapstart("DYN_Tinteg", 2)

    call ATMOS_DYN_tinteg_large( DENS,    MOMZ,    MOMX,    MOMY,    RHOT,    QTRC,    & ! [INOUT]
                                 PROG,                                                 & ! [INOUT]
                                 DENS_av, MOMZ_av, MOMX_av, MOMY_av, RHOT_av, QTRC_av, & ! [INOUT]
                                 num_diff, num_diff_q,                                 & ! [OUT;WORK]
                                 DENS_tp, MOMZ_tp, MOMX_tp, MOMY_tp, RHOT_tp, RHOQ_tp, & ! [IN]
                                 CORIOLIS,                                             & ! [IN]
                                 CDZ, CDX, CDY, FDZ, FDX, FDY,                         & ! [IN]
                                 RCDZ, RCDX, RCDY, RFDZ, RFDX, RFDY,                   & ! [IN]
                                 PHI, GSQRT,                                           & ! [IN]
                                 J13G, J23G, J33G, MAPF,                               & ! [IN]
                                 AQ_R, AQ_CV, AQ_CP, AQ_MASS,                          & ! [IN]
                                 REF_dens, REF_pott, REF_qv, REF_pres,                 & ! [IN]
                                 BND_W, BND_E, BND_S, BND_N, PRC_TwoD,                 & ! [IN]
                                 ND_COEF, ND_COEF_Q, ND_LAPLACIAN_NUM,                 & ! [IN]
                                 ND_SFC_FACT, ND_USE_RS,                               & ! [IN]
                                 BND_QA, BND_IQ, BND_SMOOTHER_FACT,                    & ! [IN]
                                 DAMP_DENS,       DAMP_VELZ,       DAMP_VELX,          & ! [IN]
                                 DAMP_VELY,       DAMP_POTT,       DAMP_QTRC,          & ! [IN]
                                 DAMP_alpha_DENS, DAMP_alpha_VELZ, DAMP_alpha_VELX,    & ! [IN]
                                 DAMP_alpha_VELY, DAMP_alpha_POTT, DAMP_alpha_QTRC,    & ! [IN]
                                 MFLUX_OFFSET_X, MFLUX_OFFSET_Y,                       & ! [IN] 
                                 wdamp_coef, divdmp_coef,                              & ! [IN]
                                 FLAG_TRACER_SPLIT_TEND,                               & ! [IN]
                                 FLAG_FCT_MOMENTUM, FLAG_FCT_T, FLAG_FCT_TRACER,       & ! [IN]
                                 FLAG_FCT_ALONG_STREAM,                                & ! [IN]
                                 USE_AVERAGE,                                          & ! [IN]
                                 I_QV,                                                 & ! [IN]
                                 DTSEC, DTSEC_DYN                                      ) ! [IN]

    call PROF_rapend  ("DYN_Tinteg", 2)

    return
  end subroutine ATMOS_DYN

end module scale_atmos_dyn
