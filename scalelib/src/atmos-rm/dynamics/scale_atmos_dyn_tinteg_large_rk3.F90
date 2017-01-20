!-------------------------------------------------------------------------------
!> module Atmosphere / Dyn Tinteg
!!
!! @par Description
!!          Temporal integration in Dynamical core for Atmospheric process
!!          three step Runge-Kutta scheme
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2016-04-18 (S.Nishizawa) [new]
!!
!<
!-------------------------------------------------------------------------------
#include "inc_openmp.h"
module scale_atmos_dyn_tinteg_large_rk3
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
  public :: ATMOS_DYN_Tinteg_large_rk3_setup
  public :: ATMOS_DYN_Tinteg_large_rk3

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
  subroutine ATMOS_DYN_Tinteg_large_rk3_setup( &
       tinteg_type )
    use scale_process, only: &
       PRC_MPIstop
    use scale_const, only: &
       UNDEF => CONST_UNDEF
    use scale_comm, only: &
       COMM_vars8_init
    implicit none

    character(len=*) :: tinteg_type

    integer :: iv
    !---------------------------------------------------------------------------

    if ( tinteg_type /= 'RK3' ) then
       write(*,*) 'xxx TINTEG_LARGE_TYPE is not RK3. Check!'
       call PRC_MPIstop
    end if

    return
  end subroutine ATMOS_DYN_Tinteg_large_rk3_setup

  !-----------------------------------------------------------------------------
  !> RK3
  subroutine ATMOS_DYN_tinteg_large_rk3( &
       DENS,    MOMZ,    MOMX,    MOMY,    RHOT,    QTRC,    &
       PROG,                                                 &
       DENS_av, MOMZ_av, MOMX_av, MOMY_av, RHOT_av, QTRC_av, &
       mflx_hi, tflx_hi,                                     &
       num_diff, num_diff_q,                                 &
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
       DTL, DTS                                              )
    use scale_const, only: &
       Rdry   => CONST_Rdry, &
       Rvap   => CONST_Rvap, &
       CVdry  => CONST_CVdry
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
    use scale_atmos_dyn_tstep_large, only: &
       ATMOS_DYN_tstep_large
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

    real(DP), intent(in)    :: DTL
    real(DP), intent(in)    :: DTS

    ! for time integartion
    real(RP) :: DENS0   (KA,IA,JA)
    real(RP) :: MOMZ0   (KA,IA,JA)
    real(RP) :: MOMX0   (KA,IA,JA)
    real(RP) :: MOMY0   (KA,IA,JA)
    real(RP) :: RHOT0   (KA,IA,JA)
    real(RP) :: QTRC0   (KA,IA,JA,QA)
    real(RP) :: PROG0   (KA,IA,JA,VA)

    real(DP) :: dtrk

    logical :: last

    integer :: n, iq

    !---------------------------------------------------------------------------

    call PROF_rapstart("DYN_Large_RK3_Prep", 3)

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
!OCL XFILL
       QTRC0(:,:,:,:) = QTRC(:,:,:,:)

       if ( VA > 0 ) then
!OCL XFILL
          PROG0(:,:,:,:) = PROG(:,:,:,:)
       end if

       call PROF_rapend  ("DYN_Large_RK3_Prep", 3)

       !------------------------------------------------------------------------
       ! Start RK
       !------------------------------------------------------------------------

       call PROF_rapstart("DYN_Large_RK3", 3)

       do n = 1, 3

          dtrk = DTL / real( 4-n, kind=DP )
          last = n == 3

          if ( n > 1 ) then
             DENS = DENS0
             MOMZ = MOMZ0
             MOMX = MOMX0
             MOMY = MOMY0
             RHOT = RHOT0
             if ( VA > 0 ) PROG = PROG0
          end if

          call ATMOS_DYN_tstep_large( &
               DENS,    MOMZ,    MOMX,    MOMY,    RHOT,    QTRC,    PROG,  & ! (inout)
               DENS_av, MOMZ_av, MOMX_av, MOMY_av, RHOT_av, QTRC_av,        & ! (inout)
               mflx_hi, tflx_hi,                                            & ! (out)
               num_diff, num_diff_q,                                        & ! (out)
               QTRC0,                                                       & ! (in)
               DENS_tp, MOMZ_tp, MOMX_tp, MOMY_tp, RHOT_tp, RHOQ_tp,        & ! (in)
               CORIOLI,                                                     & ! (in)
               CDZ, CDX, CDY, FDZ, FDX, FDY,                                & ! (in)
               RCDZ, RCDX, RCDY, RFDZ, RFDX, RFDY,                          & ! (in)
               PHI, GSQRT,                                                  & ! (in)
               J13G, J23G, J33G, MAPF,                                      & ! (in)
               AQ_R, AQ_CV, AQ_CP, AQ_MASS,                                 & ! (in)
               REF_dens, REF_pott, REF_qv, REF_pres,                        & ! (in)
               BND_W, BND_E, BND_S, BND_N,                                  & ! (in)
               ND_COEF, ND_COEF_Q, ND_ORDER, ND_SFC_FACT, ND_USE_RS,        & ! (in)
               DAMP_DENS,       DAMP_VELZ,       DAMP_VELX,                 & ! (in)
               DAMP_VELY,       DAMP_POTT,       DAMP_QTRC,                 & ! (in)
               DAMP_alpha_DENS, DAMP_alpha_VELZ, DAMP_alpha_VELX,           & ! (in)
               DAMP_alpha_VELY, DAMP_alpha_POTT, DAMP_alpha_QTRC,           & ! (in)
               wdamp_coef,                                                  & ! (in)
               divdmp_coef,                                                 & ! (in)
               FLAG_FCT_MOMENTUM, FLAG_FCT_T, FLAG_FCT_TRACER,              & ! (in)
               FLAG_FCT_ALONG_STREAM,                                       & ! (in)
               USE_AVERAGE .AND. last,                                      & ! (in)
               dtrk, dts, last                                              ) ! (in)

       end do

       call PROF_rapend  ("DYN_Large_RK3", 3)

    return
  end subroutine ATMOS_DYN_tinteg_large_rk3

end module scale_atmos_dyn_tinteg_large_rk3
