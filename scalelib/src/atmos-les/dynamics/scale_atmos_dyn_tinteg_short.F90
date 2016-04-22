!-------------------------------------------------------------------------------
!> module Atmosphere / Dynamics Temporal integration
!!
!! @par Description
!!          Temporal integration scheme selecter for dynamical short time step
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2016-04-18 (S.Nishizawa) [new]
!!
!<
!-------------------------------------------------------------------------------
#include "inc_openmp.h"
module scale_atmos_dyn_tinteg_short
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
  public :: ATMOS_DYN_Tinteg_short_setup

  !> Runge-Kutta loop
  abstract interface
     subroutine short( &
          DENS, MOMZ, MOMX, MOMY, RHOT, PROG,       & ! (inout)
          mflx_hi, tflx_hi,                         & ! (inout)
          DENS0, MOMZ0, MOMX0, MOMY0, RHOT0, PROG0, & ! (in)
          DENS_t, MOMZ_t, MOMX_t, MOMY_t, RHOT_t,   & ! (in)
          Rtot, CVtot, CORIOLI,                     & ! (in)
          num_diff, divdmp_coef, DDIV,              & ! (in)
          FLAG_FCT_MOMENTUM, FLAG_FCT_T,            & ! (in)
          FLAG_FCT_ALONG_STREAM,                    & ! (in)
          CDZ, FDZ, FDX, FDY,                       & ! (in)
          RCDZ, RCDX, RCDY, RFDZ, RFDX, RFDY,       & ! (in)
          PHI, GSQRT, J13G, J23G, J33G, MAPF,       & ! (in)
          REF_pres, REF_dens,                       & ! (in)
          BND_W, BND_E, BND_S, BND_N,               & ! (in)
          dt                                        ) ! (in)
       use scale_precision
       use scale_grid_index
       use scale_index
       real(RP), intent(inout) :: DENS(KA,IA,JA)
       real(RP), intent(inout) :: MOMZ(KA,IA,JA)
       real(RP), intent(inout) :: MOMX(KA,IA,JA)
       real(RP), intent(inout) :: MOMY(KA,IA,JA)
       real(RP), intent(inout) :: RHOT(KA,IA,JA)
       real(RP), intent(inout) :: PROG(KA,IA,JA,VA)

       real(RP), intent(inout) :: mflx_hi(KA,IA,JA,3)
       real(RP), intent(inout) :: tflx_hi(KA,IA,JA,3)

       real(RP), intent(in)    :: DENS0(KA,IA,JA)
       real(RP), intent(in)    :: MOMZ0(KA,IA,JA)
       real(RP), intent(in)    :: MOMX0(KA,IA,JA)
       real(RP), intent(in)    :: MOMY0(KA,IA,JA)
       real(RP), intent(in)    :: RHOT0(KA,IA,JA)
       real(RP), intent(in)    :: PROG0(KA,IA,JA,VA)

       real(RP), intent(in)    :: DENS_t(KA,IA,JA)
       real(RP), intent(in)    :: MOMZ_t(KA,IA,JA)
       real(RP), intent(in)    :: MOMX_t(KA,IA,JA)
       real(RP), intent(in)    :: MOMY_t(KA,IA,JA)
       real(RP), intent(in)    :: RHOT_t(KA,IA,JA)

       real(RP), intent(in)    :: Rtot(KA,IA,JA)
       real(RP), intent(in)    :: CVtot(KA,IA,JA)
       real(RP), intent(in)    :: CORIOLI(IA,JA)

       real(RP), intent(in)    :: num_diff(KA,IA,JA,5,3)
       real(RP), intent(in)    :: divdmp_coef
       real(RP), intent(in)    :: DDIV(KA,IA,JA)

       logical,  intent(in)    :: FLAG_FCT_MOMENTUM
       logical,  intent(in)    :: FLAG_FCT_T
       logical,  intent(in)    :: FLAG_FCT_ALONG_STREAM

       real(RP), intent(in)    :: CDZ (KA)
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

       real(RP), intent(in)    :: REF_pres(KA,IA,JA)   !< reference pressure
       real(RP), intent(in)    :: REF_dens(KA,IA,JA)

       logical,  intent(in)    :: BND_W
       logical,  intent(in)    :: BND_E
       logical,  intent(in)    :: BND_S
       logical,  intent(in)    :: BND_N

       real(RP), intent(in)    :: dt
     end subroutine short
  end interface
  procedure(short), pointer :: ATMOS_DYN_Tinteg_short => NULL()
  public :: ATMOS_DYN_Tinteg_short

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
  subroutine ATMOS_DYN_Tinteg_short_setup( &
       ATMOS_DYN_Tinteg_short_TYPE )

    use scale_precision
    use scale_grid_index
    use scale_index
    use scale_process, only: &
       PRC_MPIstop
#define EXTM(pre, name, post) pre ## name ## post
#define NAME(pre, name, post) EXTM(pre, name, post)
#ifdef TINTEG_SHORT
    use NAME(scale_atmos_dyn_tinteg_short_, TINTEG_SHORT,), only: &
       NAME(ATMOS_DYN_rk_tinteg_short_, TINTEG_SHORT, _setup), &
       NAME(ATMOS_DYN_rk_tinteg_short_, TINTEG_SHORT,)
#else
    use scale_atmos_dyn_tinteg_short_rk3, only: &
       ATMOS_DYN_Tinteg_short_rk3_setup, &
       ATMOS_DYN_Tinteg_short_rk3
#endif
    implicit none
    character(len=*), intent(in)  :: ATMOS_DYN_Tinteg_short_TYPE
    !---------------------------------------------------------------------------

#ifdef TINTEG_SHORT
    NAME(ATMOS_DYN_Tinteg_short_, TINTEG_SHORT, _setup)( &
            ATMOS_DYN_Tinteg_short_TYPE )
    ATMOS_DYN_Tinteg_short => NAME(ATMOS_DYN_Tingeg_short_, TINTEG_SHORT,)
#else
    select case ( ATMOS_DYN_Tinteg_short_TYPE )
    case ( 'RK3' )
       call ATMOS_DYN_Tinteg_short_rk3_setup( &
            ATMOS_DYN_Tinteg_short_TYPE )
       ATMOS_DYN_Tinteg_short => ATMOS_DYN_Tinteg_short_rk3
    case ( 'OFF', 'NONE' )
       ! do nothing
    case default
       write(*,*) 'xxx ATMOS_DYN_TINTEG_SHORT_TYPE is invalid: ', ATMOS_DYN_Tinteg_short_TYPE
       call PRC_MPIstop
    end select
#endif

    return
  end subroutine ATMOS_DYN_Tinteg_short_setup

end module scale_atmos_dyn_tinteg_short
