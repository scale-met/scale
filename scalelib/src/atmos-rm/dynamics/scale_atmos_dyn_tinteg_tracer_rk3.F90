!-------------------------------------------------------------------------------
!> module Atmosphere / Dyn Tinteg
!!
!! @par Description
!!          Temporal integration for tracer advection for Atmospheric process
!!          three step Runge-Kutta scheme
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2016-05-17 (S.Nishizawa) [new]
!!
!<
!-------------------------------------------------------------------------------
#include "inc_openmp.h"
module scale_atmos_dyn_tinteg_tracer_rk3
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
    use scale_process, only: &
       PRC_MPIstop
    use scale_comm, only: &
       COMM_vars8_init
    implicit none

    character(len=*) :: tinteg_type

    integer :: iv
    !---------------------------------------------------------------------------

    if ( tinteg_type /= 'RK3WS2002' ) then
       write(*,*) 'xxx TINTEG_LARGE_TYPE is not RK3WS2002. Check!'
       call PRC_MPIstop
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

    real(RP) :: DENS_RK(KA,IA,JA)
    real(RP) :: dtrk
    integer :: k, i, j

    !------------------------------------------------------------------------
    ! Start RK
    !------------------------------------------------------------------------

    do j = JS-1, JE+1
    do i = IS-1, IE+1
    do k = KS, KE
       DENS_RK(k,i,j) = DENS0(k,i,j) &
                      + ( DENS(k,i,j) - DENS0(k,i,j) ) / 3.0_RP
    end do
    end do
    end do

    dtrk = DTL / 3.0_RP
    call ATMOS_DYN_tstep_tracer( &
         QTRC_RK1, & ! (out)
         QTRC, QTRC0, RHOQ_t, &! (in)
         DENS0, DENS_RK, & ! (in)
         mflx_hi, num_diff, & ! (in)
         GSQRT, MAPF, & ! (in)
         CDZ, RCDZ, RCDX, RCDY, & ! (in)
         dtrk, & ! (in)
         .false., FLAG_FCT_ALONG_STREAM ) ! (in)

    call COMM_vars8( QTRC_RK1(:,:,:), I_COMM_RK1 )
    call COMM_wait ( QTRC_RK1(:,:,:), I_COMM_RK1, .false. )


    do j = JS-1, JE+1
    do i = IS-1, IE+1
    do k = KS, KE
       DENS_RK(k,i,j) = DENS0(k,i,j) &
                      + ( DENS(k,i,j) - DENS0(k,i,j) ) * 0.5_RP
    end do
    end do
    end do

    dtrk = DTL / 2.0_RP
    call ATMOS_DYN_tstep_tracer( &
         QTRC_RK2, & ! (out)
         QTRC_RK1, QTRC0, RHOQ_t, &! (in)
         DENS0, DENS_RK, & ! (in)
         mflx_hi, num_diff, & ! (in)
         GSQRT, MAPF, & ! (in)
         CDZ, RCDZ, RCDX, RCDY, & ! (in)
         dtrk, & ! (in)
         .false., FLAG_FCT_ALONG_STREAM ) ! (in)

    call COMM_vars8( QTRC_RK2(:,:,:), I_COMM_RK2 )
    call COMM_wait ( QTRC_RK2(:,:,:), I_COMM_RK2, .false. )


    dtrk = DTL
    call ATMOS_DYN_tstep_tracer( &
         QTRC, & ! (out)
         QTRC_RK2, QTRC0, RHOQ_t, &! (in)
         DENS0, DENS, & ! (in)
         mflx_hi, num_diff, & ! (in)
         GSQRT, MAPF, & ! (in)
         CDZ, RCDZ, RCDX, RCDY, & ! (in)
         dtrk, & ! (in)
         FLAG_FCT_TRACER, FLAG_FCT_ALONG_STREAM ) ! (in)

    return
  end subroutine ATMOS_DYN_tinteg_tracer_rk3

end module scale_atmos_dyn_tinteg_tracer_rk3
