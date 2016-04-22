!-------------------------------------------------------------------------------
!> module Atmosphere / Dynamical scheme
!!
!! @par Description
!!          Dynamical scheme selecter for dynamical large time step
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2016-04-18 (S.Nishizawa) [new]
!!
!<
!-------------------------------------------------------------------------------
#include "inc_openmp.h"
module scale_atmos_dyn_tstep_large
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
          mflx_hi, tflx_hi,                                     &
          num_diff, num_diff_q,                                 &
          DENS0, MOMZ0, MOMX0, MOMY0, RHOT0, RHOQ0, PROG0,      &
          DENS_tp, MOMZ_tp, MOMX_tp, MOMY_tp, RHOT_tp, RHOQ_tp, &
          CORIOLI,                                              &
          CDZ, CDX, CDY, FDZ, FDX, FDY,                         &
          RCDZ, RCDX, RCDY, RFDZ, RFDX, RFDY,                   &
          PHI, GSQRT,                                           &
          J13G, J23G, J33G, MAPF,                               &
          AQ_CV,                                                &
          REF_dens, REF_pott, REF_qv, REF_pres,                 &
          BND_W, BND_E, BND_S, BND_N,                           & ! (in)
          ND_COEF, ND_COEF_Q, ND_ORDER, ND_SFC_FACT, ND_USE_RS, &
          DAMP_DENS,       DAMP_VELZ,       DAMP_VELX,          &
          DAMP_VELY,       DAMP_POTT,       DAMP_QTRC,          &
          DAMP_alpha_DENS, DAMP_alpha_VELZ, DAMP_alpha_VELX,    &
          DAMP_alpha_VELY, DAMP_alpha_POTT, DAMP_alpha_QTRC,    &
          divdmp_coef,                                          &
          FLAG_FCT_MOMENTUM, FLAG_FCT_T, FLAG_FCT_TRACER,       &
          FLAG_FCT_ALONG_STREAM,                                &
          USE_AVERAGE,                                          &
          DDT, NSTEP_DYN, Llast                                 )
       use scale_precision
       use scale_grid_index
       use scale_index
       use scale_tracer
       use scale_atmos_boundary, only: &
            BND_QA
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

       real(RP), intent(in)    :: DENS0(KA,IA,JA)
       real(RP), intent(in)    :: MOMZ0(KA,IA,JA)
       real(RP), intent(in)    :: MOMX0(KA,IA,JA)
       real(RP), intent(in)    :: MOMY0(KA,IA,JA)
       real(RP), intent(in)    :: RHOT0(KA,IA,JA)
       real(RP), intent(in)    :: RHOQ0(KA,IA,JA,QA)
       real(RP), intent(in)    :: PROG0(KA,IA,JA,VA)

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

       real(RP), intent(in)    :: AQ_CV(QQA)

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

       real(RP), intent(in)    :: divdmp_coef

       logical,  intent(in)    :: FLAG_FCT_MOMENTUM
       logical,  intent(in)    :: FLAG_FCT_T
       logical,  intent(in)    :: FLAG_FCT_TRACER
       logical,  intent(in)    :: FLAG_FCT_ALONG_STREAM

       logical,  intent(in)    :: USE_AVERAGE

       real(DP), intent(in)    :: DDT
       integer , intent(in)    :: NSTEP_DYN
       logical , intent(in)    :: Llast
     end subroutine large

  end interface

  procedure(large), pointer :: ATMOS_DYN_Tstep_large => NULL()
  public :: ATMOS_DYN_Tstep_large

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
       DENS, MOMZ, MOMX, MOMY, RHOT, QTRC, PROG, &
       mflx_hi )
    use scale_precision
    use scale_grid_index
    use scale_index
    use scale_process, only: &
       PRC_MPIstop
    use scale_atmos_dyn_tstep_large_fvm_heve, only: &
       ATMOS_DYN_Tstep_large_fvm_heve_setup, &
       ATMOS_DYN_Tstep_large_fvm_heve
    implicit none
    real(RP), intent(inout) :: DENS(KA,IA,JA)
    real(RP), intent(inout) :: MOMZ(KA,IA,JA)
    real(RP), intent(inout) :: MOMX(KA,IA,JA)
    real(RP), intent(inout) :: MOMY(KA,IA,JA)
    real(RP), intent(inout) :: RHOT(KA,IA,JA)
    real(RP), intent(inout) :: QTRC(KA,IA,JA,QA)
    real(RP), intent(inout) :: PROG(KA,IA,JA,VA)
    real(RP), intent(inout) :: mflx_hi(KA,IA,JA,3)

    character(len=H_SHORT) :: Tstep_large_type = 'FVM-HEVE'

    !---------------------------------------------------------------------------

    select case ( Tstep_large_type )
    case ( 'FVM-HEVE' )
       call ATMOS_DYN_Tstep_large_fvm_heve_setup( &
            DENS, MOMZ, MOMX, MOMY, RHOT, QTRC, PROG, &
            mflx_hi )
       ATMOS_DYN_Tstep_large => ATMOS_DYN_Tstep_large_fvm_heve
    case default
       write(*,*) 'xxx ATMOS_DYN_Tstep_large_type is invalid: ', Tstep_large_type
       call PRC_MPIstop
    end select

    return
  end subroutine ATMOS_DYN_Tstep_large_setup

end module scale_atmos_dyn_tstep_large
