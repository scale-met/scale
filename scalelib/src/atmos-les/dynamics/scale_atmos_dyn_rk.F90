!-------------------------------------------------------------------------------
!> module Atmosphere / Dynamics RK
!!
!! @par Description
!!          Runge-Kutta for Atmospheric dynamical process
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2013-12-04 (S.Nishizawa) [new]
!!
!<
!-------------------------------------------------------------------------------
#include "inc_openmp.h"
module scale_atmos_dyn_rk
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_stdio
  use scale_prof
  use scale_grid_index
  use scale_tracer
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: ATMOS_DYN_rk_setup

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  !-----------------------------------------------------------------------------
  !> Runge-Kutta loop
  abstract interface
     subroutine rk( &
          DENS_RK, MOMZ_RK, MOMX_RK, MOMY_RK, RHOT_RK, &
          mflx_hi,                                     &
          DENS0,   MOMZ0,   MOMX0,   MOMY0,   RHOT0,   &
          DENS,    MOMZ,    MOMX,    MOMY,    RHOT,    &
          DENS_t,  MOMZ_t,  MOMX_t,  MOMY_t,  RHOT_t,  &
          Rtot, CVtot, CORIOLI,                        &
          num_diff, divdmp_coef,                       &
          FLAG_FCT_RHO, FLAG_FCT_MOMENTUM, FLAG_FCT_T, &
          CDZ, FDZ, FDX, FDY,                          &
          RCDZ, RCDX, RCDY, RFDZ, RFDX, RFDY,          &
          PHI, GSQRT, J13G, J23G, J33G,                &
          REF_pres, REF_dens,                          &
          dtrk                                         )
       use scale_precision
       use scale_grid_index
       real(RP), intent(out) :: DENS_RK(KA,IA,JA)   ! prognostic variables
       real(RP), intent(out) :: MOMZ_RK(KA,IA,JA)   !
       real(RP), intent(out) :: MOMX_RK(KA,IA,JA)   !
       real(RP), intent(out) :: MOMY_RK(KA,IA,JA)   !
       real(RP), intent(out) :: RHOT_RK(KA,IA,JA)   !

       real(RP), intent(out) :: mflx_hi(KA,IA,JA,3) ! mass flux

       real(RP), intent(in),target :: DENS0(KA,IA,JA) ! prognostic variables at previous dynamical time step
       real(RP), intent(in),target :: MOMZ0(KA,IA,JA) !
       real(RP), intent(in),target :: MOMX0(KA,IA,JA) !
       real(RP), intent(in),target :: MOMY0(KA,IA,JA) !
       real(RP), intent(in),target :: RHOT0(KA,IA,JA) !

       real(RP), intent(in)  :: DENS(KA,IA,JA)      ! prognostic variables at previous RK step
       real(RP), intent(in)  :: MOMZ(KA,IA,JA)      !
       real(RP), intent(in)  :: MOMX(KA,IA,JA)      !
       real(RP), intent(in)  :: MOMY(KA,IA,JA)      !
       real(RP), intent(in)  :: RHOT(KA,IA,JA)      !

       real(RP), intent(in)  :: DENS_t(KA,IA,JA)    ! tendency
       real(RP), intent(in)  :: MOMZ_t(KA,IA,JA)    !
       real(RP), intent(in)  :: MOMX_t(KA,IA,JA)    !
       real(RP), intent(in)  :: MOMY_t(KA,IA,JA)    !
       real(RP), intent(in)  :: RHOT_t(KA,IA,JA)    !

       real(RP), intent(in)  :: Rtot    (KA,IA,JA)  ! total R
       real(RP), intent(in)  :: CVtot   (KA,IA,JA)  ! total CV
       real(RP), intent(in)  :: CORIOLI (1, IA,JA)
       real(RP), intent(in)  :: num_diff(KA,IA,JA,5,3)
       real(RP), intent(in)  :: divdmp_coef

       logical,  intent(in)  :: FLAG_FCT_RHO
       logical,  intent(in)  :: FLAG_FCT_MOMENTUM
       logical,  intent(in)  :: FLAG_FCT_T

       real(RP), intent(in)  :: CDZ (KA)
       real(RP), intent(in)  :: FDZ (KA-1)
       real(RP), intent(in)  :: FDX (IA-1)
       real(RP), intent(in)  :: FDY (JA-1)
       real(RP), intent(in)  :: RCDZ(KA)
       real(RP), intent(in)  :: RCDX(IA)
       real(RP), intent(in)  :: RCDY(JA)
       real(RP), intent(in)  :: RFDZ(KA-1)
       real(RP), intent(in)  :: RFDX(IA-1)
       real(RP), intent(in)  :: RFDY(JA-1)

       real(RP), intent(in)  :: PHI     (KA,IA,JA)   !< geopotential
       real(RP), intent(in)  :: GSQRT   (KA,IA,JA,7) !< vertical metrics {G}^1/2
       real(RP), intent(in)  :: J13G    (KA,IA,JA,7) !< (1,3) element of Jacobian matrix
       real(RP), intent(in)  :: J23G    (KA,IA,JA,7) !< (2,3) element of Jacobian matrix
       real(RP), intent(in)  :: J33G                 !< (3,3) element of Jacobian matrix
       real(RP), intent(in)  :: REF_pres(KA,IA,JA)   !< reference pressure
       real(RP), intent(in)  :: REF_dens(KA,IA,JA)   !< reference density

       real(RP), intent(in)  :: dtrk
     end subroutine rk
  end interface
  procedure(rk), pointer :: ATMOS_DYN_rk => NULL()
  public :: ATMOS_DYN_rk


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
  subroutine ATMOS_DYN_rk_setup( ATMOS_TYPE_DYN )
    use scale_process, only: &
       PRC_MPIstop
    use scale_stdio, only: &
       IO_FID_LOG
#define EXTM(pre, name, post) pre ## name ## post
#define NAME(pre, name, post) EXTM(pre, name, post)
#ifdef DYNAMICS
    use NAME(scale_atmos_dyn_rk_, DYNAMICS,), only: &
       NAME(ATMOS_DYN_rk_, DYNAMICS, _setup), &
       NAME(ATMOS_DYN_rk_, DYNAMICS,)
#else
    use scale_atmos_dyn_rk_heve, only: &
       ATMOS_DYN_rk_heve_setup, &
       ATMOS_DYN_rk_heve
    use scale_atmos_dyn_rk_hevi, only: &
       ATMOS_DYN_rk_hevi_setup, &
       ATMOS_DYN_rk_hevi
    use scale_atmos_dyn_rk_hivi, only: &
       ATMOS_DYN_rk_hivi_setup, &
       ATMOS_DYN_rk_hivi
#endif
    implicit none

    character(len=H_SHORT), intent(in) :: ATMOS_TYPE_DYN
    !---------------------------------------------------------------------------

    select case ( ATMOS_TYPE_DYN )
    case ( 'HEVE' )
       call ATMOS_DYN_rk_heve_setup( ATMOS_TYPE_DYN )
       ATMOS_DYN_rk => ATMOS_DYN_rk_heve
    case ( 'HEVI' )
       call ATMOS_DYN_rk_hevi_setup( ATMOS_TYPE_DYN )
       ATMOS_DYN_rk => ATMOS_DYN_rk_hevi
    case ( 'HIVI' )
       call ATMOS_DYN_rk_hivi_setup( ATMOS_TYPE_DYN )
       ATMOS_DYN_rk => ATMOS_DYN_rk_hivi
    case default
       if( IO_L ) write(IO_FID_LOG,*) 'xxx ATMOS_TYPE_DYN is invalid'
       call PRC_MPIstop
    end select

    return
  end subroutine ATMOS_DYN_rk_setup

end module scale_atmos_dyn_rk
