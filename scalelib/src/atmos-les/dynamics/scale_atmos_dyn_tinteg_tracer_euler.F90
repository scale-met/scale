!-------------------------------------------------------------------------------
!> module Atmosphere / Dyn Tinteg
!!
!! @par Description
!!          Temporal integration for tracer advection for Atmospheric process
!!          Euler scheme
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2016-05-17 (S.Nishizawa) [new]
!!
!<
!-------------------------------------------------------------------------------
#include "inc_openmp.h"
module scale_atmos_dyn_tinteg_tracer_euler
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

    if ( tinteg_type .ne. 'EULER' ) then
       write(*,*) 'xxx TINTEG_LARGE_TYPE is not EULER. Check!'
       call PRC_MPIstop
    end if

    return
  end subroutine ATMOS_DYN_Tinteg_tracer_euler_setup

  !-----------------------------------------------------------------------------
  !> EULER
  subroutine ATMOS_DYN_tinteg_tracer_euler( &
       QTRC, & ! (out)
       QTRC0, RHOQ_t, &! (in)
       DENS0, DENS, & ! (in)
       mflx_hi, num_diff, & ! (in)
       GSQRT, MAPF, & ! (in)
       CDZ, RCDZ, RCDX, RCDY, & ! (in)
       dtl, & ! (in)
       FLAG_FCT_TRACER, & ! (in)
       FLAG_FCT_ALONG_STREAM ) ! (in)
    use scale_comm, only: &
       COMM_vars8, &
       COMM_wait
    use scale_atmos_dyn_tstep_tracer, only: &
       ATMOS_DYN_tstep_tracer
    implicit none
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
    real(RP), intent(in)    :: dtl
    logical,  intent(in)    :: FLAG_FCT_TRACER
    logical,  intent(in)    :: FLAG_FCT_ALONG_STREAM

    !------------------------------------------------------------------------
    ! Start RK
    !------------------------------------------------------------------------

    call ATMOS_DYN_tstep_tracer( &
         QTRC, & ! (out)
         QTRC, QTRC0, RHOQ_t, &! (in)
         DENS0, DENS, & ! (in)
         mflx_hi, num_diff, & ! (in)
         GSQRT, MAPF, & ! (in)
         CDZ, RCDZ, RCDX, RCDY, & ! (in)
         dtl, & ! (in)
         FLAG_FCT_TRACER, FLAG_FCT_ALONG_STREAM ) ! (in)

    return
  end subroutine ATMOS_DYN_tinteg_tracer_euler

end module scale_atmos_dyn_tinteg_tracer_euler
