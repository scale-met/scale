!-------------------------------------------------------------------------------
!> module Atmosphere / Dynamical scheme
!!
!! @par Description
!!          Dynamical scheme selecter for dynamical large time step
!!
!! @author Team SCALE
!!
!<
!-------------------------------------------------------------------------------
#include "scalelib.h"
module scale_atmos_dyn_tstep_large
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
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: ATMOS_DYN_Tstep_large_setup

  abstract interface
     !> calculation values at next temporal step
     subroutine large( &
          DENS, MOMZ, MOMX, MOMY, RHOT, QTRC, PROG,             &
          DENS_av, MOMZ_av, MOMX_av, MOMY_av, RHOT_av, QTRC_av, &
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
          DTL, DTS, Llast                                       )
       use scale_precision
       use scale_atmos_grid_cartesC_index
       use scale_index
       use scale_tracer
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

       logical , intent(in)    :: Llast
     end subroutine large

     subroutine finalize
     end subroutine finalize
  end interface

  procedure(large), pointer :: ATMOS_DYN_Tstep_large => NULL()
  public :: ATMOS_DYN_Tstep_large
  procedure(finalize), pointer :: ATMOS_DYN_Tstep_large_finalize => NULL()
  public :: ATMOS_DYN_Tstep_large_finalize

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
  !> Register
  subroutine ATMOS_DYN_Tstep_large_setup( &
       Tstep_large_type, &
       DENS, MOMZ, MOMX, MOMY, RHOT, QTRC, PROG )
    use scale_precision
    use scale_atmos_grid_cartesC_index
    use scale_index
    use scale_prc, only: &
       PRC_abort
    use scale_atmos_dyn_tstep_large_fvm_heve, only: &
       ATMOS_DYN_Tstep_large_fvm_heve_setup,    &
       ATMOS_DYN_Tstep_large_fvm_heve_finalize, &
       ATMOS_DYN_Tstep_large_fvm_heve
    implicit none
    character(len=*), intent(in) :: Tstep_large_type
    real(RP), intent(inout) :: DENS(KA,IA,JA)
    real(RP), intent(inout) :: MOMZ(KA,IA,JA)
    real(RP), intent(inout) :: MOMX(KA,IA,JA)
    real(RP), intent(inout) :: MOMY(KA,IA,JA)
    real(RP), intent(inout) :: RHOT(KA,IA,JA)
    real(RP), intent(inout) :: QTRC(KA,IA,JA,QA)
    real(RP), intent(inout) :: PROG(KA,IA,JA,VA)

    !---------------------------------------------------------------------------

    ATMOS_DYN_Tstep_large_finalize => DYN_Tstep_large_finalize

    select case( Tstep_large_type )
    case( 'FVM-HEVE' )
       call ATMOS_DYN_Tstep_large_fvm_heve_setup( &
            DENS, MOMZ, MOMX, MOMY, RHOT, QTRC, PROG )
       ATMOS_DYN_Tstep_large => ATMOS_DYN_Tstep_large_fvm_heve
       ATMOS_DYN_Tstep_large_finalize => ATMOS_DYN_Tstep_large_fvm_heve_finalize
    case default
       LOG_ERROR("ATMOS_DYN_Tstep_large_setup",*) 'ATMOS_DYN_Tstep_large_type is invalid: ', Tstep_large_type
       call PRC_abort
    end select

    return
  end subroutine ATMOS_DYN_Tstep_large_setup

  subroutine DYN_Tstep_large_finalize

    return
  end subroutine DYN_Tstep_large_finalize

end module scale_atmos_dyn_tstep_large
