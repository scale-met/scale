!-------------------------------------------------------------------------------
!> module Atmosphere / Dyn Tinteg
!!
!! @par Description
!!          Temporal integration for tracer advection for Atmospheric process
!!          Euler scheme
!!
!! @author Team SCALE
!!
!<
!-------------------------------------------------------------------------------
#include "scalelib.h"
module scale_atmos_dyn_tinteg_tracer_euler
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
  public :: ATMOS_DYN_Tinteg_tracer_euler_setup
  public :: ATMOS_DYN_Tinteg_tracer_euler

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
  subroutine ATMOS_DYN_Tinteg_tracer_euler_setup( &
       tinteg_type )
    use scale_prc, only: &
       PRC_abort
    implicit none

    character(len=*) :: tinteg_type

    integer :: iv
    !---------------------------------------------------------------------------

    if ( tinteg_type /= 'EULER' ) then
       LOG_ERROR("ATMOS_DYN_Tinteg_tracer_euler_setup",*) 'TINTEG_TRACER_TYPE is not EULER. Check!'
       call PRC_abort
    end if

    return
  end subroutine ATMOS_DYN_Tinteg_tracer_euler_setup

  !-----------------------------------------------------------------------------
  !> EULER
  subroutine ATMOS_DYN_tinteg_tracer_euler( &
       QTRC, & ! (out)
       qflx, & ! (out)
       QTRC0, RHOQ_t, &! (in)
       DENS0, DENS, & ! (in)
       mflx_hi, num_diff, & ! (in)
       GSQRT, MAPF, & ! (in)
       CDZ, RCDZ, RCDX, RCDY, & ! (in)
       BND_W, BND_E, BND_S, BND_N, & ! (in)
       TwoD, & ! (in)
       dtl, & ! (in)
       FLAG_FCT_TRACER, & ! (in)
       FLAG_FCT_ALONG_STREAM ) ! (in)
    use scale_comm_cartesC, only: &
       COMM_vars8, &
       COMM_wait
    use scale_atmos_dyn_tstep_tracer, only: &
       ATMOS_DYN_tstep_tracer
    implicit none
    real(RP), intent(inout) :: QTRC    (KA,IA,JA)
    real(RP), intent(out)   :: qflx    (KA,IA,JA,3)
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
    logical,  intent(in)    :: TwoD
    real(RP), intent(in)    :: dtl
    logical,  intent(in)    :: FLAG_FCT_TRACER
    logical,  intent(in)    :: FLAG_FCT_ALONG_STREAM

    !------------------------------------------------------------------------
    ! Start RK
    !------------------------------------------------------------------------

    call ATMOS_DYN_tstep_tracer( &
         QTRC, & ! (out)
         qflx, & ! (out)
         QTRC, QTRC0, RHOQ_t, &! (in)
         DENS0, DENS, & ! (in)
         mflx_hi, num_diff, & ! (in)
         GSQRT, MAPF, & ! (in)
         CDZ, RCDZ, RCDX, RCDY, & ! (in)
         TwoD, & ! (in)
         dtl, & ! (in)
         FLAG_FCT_TRACER, FLAG_FCT_ALONG_STREAM ) ! (in)

    return
  end subroutine ATMOS_DYN_tinteg_tracer_euler

end module scale_atmos_dyn_tinteg_tracer_euler
