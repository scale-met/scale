!-------------------------------------------------------------------------------
!> module Atmosphere / Dyn Tinteg
!!
!! @par Description
!!          Temporal integration in Dynamical core for Atmospheric process
!!          the Euler scheme
!!
!! @author Team SCALE
!!
!<
!-------------------------------------------------------------------------------
#include "scalelib.h"
module scale_atmos_dyn_tinteg_large_euler
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
  public :: ATMOS_DYN_Tinteg_large_euler_setup
  public :: ATMOS_DYN_Tinteg_large_euler

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
  !-----------------------------------------------------------------------------
contains

  !-----------------------------------------------------------------------------
  !> Setup
  subroutine ATMOS_DYN_Tinteg_large_euler_setup( &
       tinteg_type )
    use scale_prc, only: &
       PRC_abort
    use scale_const, only: &
       UNDEF => CONST_UNDEF
    implicit none

    character(len=*) :: tinteg_type

    integer :: iv
    !---------------------------------------------------------------------------

    if ( tinteg_type /= 'EULER' ) then
       LOG_ERROR("ATMOS_DYN_Tinteg_large_euler_setup",*) 'TINTEG_LARGE_TYPE is not EULER. Check!'
       call PRC_abort
    end if

    return
  end subroutine ATMOS_DYN_Tinteg_large_euler_setup

  !-----------------------------------------------------------------------------
  !> RK3
  subroutine ATMOS_DYN_tinteg_large_euler( &
       DENS, MOMZ, MOMX, MOMY, RHOT, QTRC, PROG,             &
       DENS_av, MOMZ_av, MOMX_av, MOMY_av, RHOT_av, QTRC_av, &
       num_diff, num_diff_q,                                 &
       DENS_tp, MOMZ_tp, MOMX_tp, MOMY_tp, RHOT_tp, RHOQ_tp, &
       CORIOLI,                                              &
       CDZ, CDX, CDY, FDZ, FDX, FDY,                         &
       RCDZ, RCDX, RCDY, RFDZ, RFDX, RFDY,                   &
       PHI, GSQRT,                                           &
       J13G, J23G, J33G, MAPF,                               &
       AQ_R, AQ_CV, AQ_CP, AQ_MASS,                          &
       REF_dens, REF_pott, REF_qv, REF_pres,                 &
       BND_W, BND_E, BND_S, BND_N, TwoD,                     &
       ND_COEF, ND_COEF_Q, ND_LAPLACIAN_NUM,                 &
       ND_SFC_FACT, ND_USE_RS,                               &
       BND_QA, BND_IQ, BND_SMOOTHER_FACT,                    &
       DAMP_DENS,       DAMP_VELZ,       DAMP_VELX,          &
       DAMP_VELY,       DAMP_POTT,       DAMP_QTRC,          &
       DAMP_alpha_DENS, DAMP_alpha_VELZ, DAMP_alpha_VELX,    &
       DAMP_alpha_VELY, DAMP_alpha_POTT, DAMP_alpha_QTRC,    &
       MFLUX_OFFSET_X, MFLUX_OFFSET_Y,                       &
       wdamp_coef, divdmp_coef,                              &
       FLAG_TRACER_SPLIT_TEND,                               &
       FLAG_FCT_MOMENTUM, FLAG_FCT_T, FLAG_FCT_TRACER,       &
       FLAG_FCT_ALONG_STREAM,                                &
       USE_AVERAGE,                                          &
       I_QV,                                                 &
       DTL, DTS                                              )
   
    use scale_atmos_dyn_tstep_large, only: &
       ATMOS_DYN_tstep_large

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

    real(RP), intent(out)   :: num_diff  (KA,IA,JA,5,3)
    real(RP), intent(out)   :: num_diff_q(KA,IA,JA,3)

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
    logical,  intent(in)    :: TwoD

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
    real(RP), intent(in)    :: MFLUX_OFFSET_Y(KA,IA,2)

    real(RP), intent(in)    :: wdamp_coef(KA)
    real(RP), intent(in)    :: divdmp_coef

    logical,  intent(in)    :: FLAG_TRACER_SPLIT_TEND
    logical,  intent(in)    :: FLAG_FCT_MOMENTUM
    logical,  intent(in)    :: FLAG_FCT_T
    logical,  intent(in)    :: FLAG_FCT_TRACER
    logical,  intent(in)    :: FLAG_FCT_ALONG_STREAM

    logical,  intent(in)    :: USE_AVERAGE

    integer,  intent(in)    :: I_QV

    real(DP), intent(in)    :: DTL
    real(DP), intent(in)    :: DTS

    real(RP) :: QTRC0(KA,IA,JA,QA)

    !$acc data copy(QTRC) create(QTRC0)

    !$omp workshare
    !$acc kernels
    QTRC0(:,:,:,:) = QTRC(:,:,:,:)
    !$acc end kernels
    !$omp end workshare


    call ATMOS_DYN_tstep_large( &
         DENS, MOMZ, MOMX, MOMY, RHOT, QTRC, PROG,             & ! (inout)
         DENS_av, MOMZ_av, MOMX_av, MOMY_av, RHOT_av, QTRC_av, & ! (inout)
         num_diff, num_diff_q,                                 & ! (out, work)
         QTRC0,                                                & ! (in)
         DENS_tp, MOMZ_tp, MOMX_tp, MOMY_tp, RHOT_tp, RHOQ_tp, & ! (in)
         CORIOLI,                                              & ! (in)
         CDZ, CDX, CDY, FDZ, FDX, FDY,                         & ! (in)
         RCDZ, RCDX, RCDY, RFDZ, RFDX, RFDY,                   & ! (in)
         PHI, GSQRT,                                           & ! (in)
         J13G, J23G, J33G, MAPF,                               & ! (in)
         AQ_R, AQ_CV, AQ_CP, AQ_MASS,                          & ! (in)
         REF_dens, REF_pott, REF_qv, REF_pres,                 & ! (in)
         BND_W, BND_E, BND_S, BND_N, TwoD,                     & ! (in)
         ND_COEF, ND_COEF_Q, ND_LAPLACIAN_NUM,                 & ! (in)
         ND_SFC_FACT, ND_USE_RS,                               & ! (in)
         BND_QA, BND_IQ, BND_SMOOTHER_FACT,                    & ! (in)
         DAMP_DENS,       DAMP_VELZ,       DAMP_VELX,          & ! (in)
         DAMP_VELY,       DAMP_POTT,       DAMP_QTRC,          & ! (in)
         DAMP_alpha_DENS, DAMP_alpha_VELZ, DAMP_alpha_VELX,    & ! (in)
         DAMP_alpha_VELY, DAMP_alpha_POTT, DAMP_alpha_QTRC,    & ! (in)
         MFLUX_OFFSET_X, MFLUX_OFFSET_Y,                       & ! (in)
         wdamp_coef, divdmp_coef,                              & ! (in)
         FLAG_TRACER_SPLIT_TEND,                               & ! (in)
         FLAG_FCT_MOMENTUM, FLAG_FCT_T, FLAG_FCT_TRACER,       & ! (in)
         FLAG_FCT_ALONG_STREAM,                                & ! (in)
         USE_AVERAGE,                                          & ! (in)
         I_QV,                                                 & ! (in)
         DTL, DTS, .true.                                      ) ! (in)

    !$acc end data

    return
  end subroutine ATMOS_DYN_tinteg_large_euler

end module scale_atmos_dyn_tinteg_large_euler
