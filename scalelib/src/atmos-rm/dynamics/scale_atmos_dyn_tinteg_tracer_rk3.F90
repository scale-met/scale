!-------------------------------------------------------------------------------
!> module Atmosphere / Dyn Tinteg
!!
!! @par Description
!!          Temporal integration for tracer advection for Atmospheric process
!!          three stage Runge-Kutta scheme
!!
!! @author Team SCALE
!!
!<
!-------------------------------------------------------------------------------
#include "scalelib.h"
module scale_atmos_dyn_tinteg_tracer_rk3
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
  public :: ATMOS_DYN_Tinteg_tracer_rk3_setup
  public :: ATMOS_DYN_Tinteg_tracer_rk3

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
  real(RP), allocatable :: QTRC_RK1(:,:,:)
  real(RP), allocatable :: QTRC_RK2(:,:,:)
  integer :: I_COMM_RK1 = 1
  integer :: I_COMM_RK2 = 1
  !-----------------------------------------------------------------------------
contains

  !-----------------------------------------------------------------------------
  !> Setup
  subroutine ATMOS_DYN_Tinteg_tracer_rk3_setup( &
       tinteg_type )
    use scale_prc, only: &
       PRC_abort
    use scale_comm_cartesC, only: &
       COMM_vars8_init
    implicit none

    character(len=*) :: tinteg_type

    integer :: iv
    !---------------------------------------------------------------------------

    if ( tinteg_type /= 'RK3WS2002' ) then
       LOG_ERROR("ATMOS_DYN_Tinteg_tracer_rk3_setup",*) 'TINTEG_LARGE_TYPE is not RK3WS2002. Check!'
       call PRC_abort
    end if

    allocate( QTRC_RK1(KA,IA,JA) )
    allocate( QTRC_RK2(KA,IA,JA) )

    call COMM_vars8_init( 'QTRC_RK1', QTRC_RK1, I_COMM_RK1 )
    call COMM_vars8_init( 'QTRC_RK2', QTRC_RK2, I_COMM_RK2 )

    return
  end subroutine ATMOS_DYN_Tinteg_tracer_rk3_setup

  !-----------------------------------------------------------------------------
  !> RK3
  subroutine ATMOS_DYN_tinteg_tracer_rk3( &
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
    use scale_atmos_dyn_common, only: &
       ATMOS_DYN_Copy_Boundary_tracer
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

    real(RP) :: DENS_RK(KA,IA,JA)
    real(RP) :: dtrk
    integer :: k, i, j

    !------------------------------------------------------------------------
    ! Start RK
    !------------------------------------------------------------------------

    do j = JS-1, JE+1
    do i = max(IS-1,1), min(IE+1,IA)
    do k = KS, KE
       DENS_RK(k,i,j) = DENS0(k,i,j) &
                      + ( DENS(k,i,j) - DENS0(k,i,j) ) / 3.0_RP
    end do
    end do
    end do

    dtrk = DTL / 3.0_RP
    call ATMOS_DYN_tstep_tracer( &
         QTRC_RK1, & ! (out)
         qflx, & ! (out)
         QTRC, QTRC0, RHOQ_t, &! (in)
         DENS0, DENS_RK, & ! (in)
         mflx_hi, num_diff, & ! (in)
         GSQRT, MAPF, & ! (in)
         CDZ, RCDZ, RCDX, RCDY, & ! (in)
         TwoD, dtrk, & ! (in)
         .false., FLAG_FCT_ALONG_STREAM ) ! (in)

    call ATMOS_DYN_Copy_boundary_tracer( QTRC_RK1,                   & ! [INOUT]
                                         QTRC0,                      & ! [IN]
                                         BND_W, BND_E, BND_S, BND_N, & ! [IN]
                                         TwoD                        ) ! [IN]

    call COMM_vars8( QTRC_RK1(:,:,:), I_COMM_RK1 )
    call COMM_wait ( QTRC_RK1(:,:,:), I_COMM_RK1, .false. )


    do j = JS-1, JE+1
    do i = max(IS-1,1), min(IE+1,IA)
    do k = KS, KE
       DENS_RK(k,i,j) = DENS0(k,i,j) &
                      + ( DENS(k,i,j) - DENS0(k,i,j) ) * 0.5_RP
    end do
    end do
    end do

    dtrk = DTL / 2.0_RP
    call ATMOS_DYN_tstep_tracer( &
         QTRC_RK2, & ! (out)
         qflx, & ! (out)
         QTRC_RK1, QTRC0, RHOQ_t, &! (in)
         DENS0, DENS_RK, & ! (in)
         mflx_hi, num_diff, & ! (in)
         GSQRT, MAPF, & ! (in)
         CDZ, RCDZ, RCDX, RCDY, & ! (in)
         TwoD, dtrk, & ! (in)
         .false., FLAG_FCT_ALONG_STREAM ) ! (in)

    call ATMOS_DYN_Copy_boundary_tracer( QTRC_RK2,                   & ! [INOUT]
                                         QTRC0,                      & ! [IN]
                                         BND_W, BND_E, BND_S, BND_N, & ! [IN]
                                         TwoD                        ) ! [IN]

    call COMM_vars8( QTRC_RK2(:,:,:), I_COMM_RK2 )
    call COMM_wait ( QTRC_RK2(:,:,:), I_COMM_RK2, .false. )


    dtrk = DTL
    call ATMOS_DYN_tstep_tracer( &
         QTRC, & ! (out)
         qflx, & ! (out)
         QTRC_RK2, QTRC0, RHOQ_t, &! (in)
         DENS0, DENS, & ! (in)
         mflx_hi, num_diff, & ! (in)
         GSQRT, MAPF, & ! (in)
         CDZ, RCDZ, RCDX, RCDY, & ! (in)
         TwoD, & ! (in)
         dtrk, & ! (in)
         FLAG_FCT_TRACER, FLAG_FCT_ALONG_STREAM ) ! (in)

    return
  end subroutine ATMOS_DYN_tinteg_tracer_rk3

end module scale_atmos_dyn_tinteg_tracer_rk3
