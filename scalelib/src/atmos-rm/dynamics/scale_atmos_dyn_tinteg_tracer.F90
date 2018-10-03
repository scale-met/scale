!-------------------------------------------------------------------------------
!> module Atmosphere / Dynamics Temporal integration
!!
!! @par Description
!!          Temporal integration scheme selecter for dynamical tracer advection
!!
!! @author Team SCALE
!!
!<
!-------------------------------------------------------------------------------
#include "scalelib.h"
module scale_atmos_dyn_tinteg_tracer
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
  public :: ATMOS_DYN_Tinteg_tracer_setup

  !> Runge-Kutta loop
  abstract interface
     subroutine tinteg( &
          QTRC, & ! (out)
          QTRC0, RHOQ_t, &! (in)
          DENS0, DENS, & ! (in)
          mflx_hi, num_diff, & ! (in)
          GSQRT, MAPF, & ! (in)
          CDZ, RCDZ, RCDX, RCDY, & ! (in)
          BND_W, BND_E, BND_S, BND_N, & ! (in)
          dtl, & ! (in)
          FLAG_FCT_TRACER, & ! (in)
          FLAG_FCT_ALONG_STREAM ) ! (in)
       use scale_precision
       use scale_atmos_grid_cartesC_index
       use scale_index
       real(RP), intent(inout) :: QTRC    (KA,IA,JA)
       real(RP), intent(in)    :: QTRC0   (KA,IA,JA)
       real(RP), intent(in)    :: RHOQ_t  (KA,IA,JA)
       real(RP), intent(in)    :: DENS0   (KA,IA,JA)
       real(RP), intent(in)    :: DENS    (KA,IA,JA)
       real(RP), intent(in)    :: mflx_hi (KA,IA,JA,3)
       real(RP), intent(in)    :: num_diff(KA,IA,JA,3)
       real(RP), intent(in)    :: GSQRT   (KA,IA,JA,7)
       real(RP), intent(in)    :: MAPF    (IA,JA)
       real(RP), intent(in)    :: CDZ(KA)
       real(RP), intent(in)    :: RCDZ(KA)
       real(RP), intent(in)    :: RCDX(IA)
       real(RP), intent(in)    :: RCDY(JA)
       logical,  intent(in)    :: BND_W
       logical,  intent(in)    :: BND_E
       logical,  intent(in)    :: BND_S
       logical,  intent(in)    :: BND_N
       real(RP), intent(in)    :: dtl
       logical,  intent(in)    :: FLAG_FCT_TRACER
       logical,  intent(in)    :: FLAG_FCT_ALONG_STREAM
     end subroutine tinteg
  end interface
  procedure(tinteg), pointer :: ATMOS_DYN_Tinteg_tracer => NULL()
  public :: ATMOS_DYN_Tinteg_tracer

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
  subroutine ATMOS_DYN_Tinteg_tracer_setup( &
       ATMOS_DYN_Tinteg_tracer_TYPE )

    use scale_precision
    use scale_atmos_grid_cartesC_index
    use scale_index
    use scale_prc, only: &
       PRC_abort
    use scale_atmos_dyn_tinteg_tracer_euler, only: &
       ATMOS_DYN_Tinteg_tracer_euler_setup, &
       ATMOS_DYN_Tinteg_tracer_euler
    use scale_atmos_dyn_tinteg_tracer_rk3, only: &
       ATMOS_DYN_Tinteg_tracer_rk3_setup, &
       ATMOS_DYN_Tinteg_tracer_rk3
    implicit none
    character(len=*), intent(in)  :: ATMOS_DYN_Tinteg_tracer_TYPE
    !---------------------------------------------------------------------------

    select case( ATMOS_DYN_Tinteg_tracer_TYPE )
    case( 'EULER' )
       call ATMOS_DYN_Tinteg_tracer_euler_setup( &
            ATMOS_DYN_Tinteg_tracer_TYPE )
       ATMOS_DYN_Tinteg_tracer => ATMOS_DYN_Tinteg_tracer_euler
    case( 'RK3WS2002' )
       call ATMOS_DYN_Tinteg_tracer_rk3_setup( &
            ATMOS_DYN_Tinteg_tracer_TYPE )
       ATMOS_DYN_Tinteg_tracer => ATMOS_DYN_Tinteg_tracer_rk3
    case( 'OFF', 'NONE' )
       ! do nothing
    case default
       LOG_ERROR("ATMOS_DYN_Tinteg_tracer_setup",*) 'ATMOS_DYN_TINTEG_TRACER_TYPE is invalid: ', ATMOS_DYN_Tinteg_tracer_TYPE
       call PRC_abort
    end select

    return
  end subroutine ATMOS_DYN_Tinteg_tracer_setup

end module scale_atmos_dyn_tinteg_tracer
