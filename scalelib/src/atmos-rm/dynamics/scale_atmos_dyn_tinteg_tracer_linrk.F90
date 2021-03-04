!-------------------------------------------------------------------------------
!> module Atmosphere / Dyn Tinteg
!!
!! @par Description
!!          Temporal integration for tracer advection for Atmospheric process
!!         'Linear case' Runge-Kutta scheme with N stage
!!          
!! For a linear, homogeneous and time-indepent ODE system defined by equation (9) in Baldauf (2008, JCP), 
!! the 'Linear case' Runge-Kutta scheme with N stage has Nth order accuracy.
!!  
!! The 'Linear case' Runge-Kutta scheme with 3 stage corresponds to the 3-stage Runge-Kutta scheme in Wicker and Skamarock (2002). 
!! @author Team SCALE
!!
!<
!-------------------------------------------------------------------------------
#include "scalelib.h"
module scale_atmos_dyn_tinteg_tracer_linrk
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
  public :: ATMOS_DYN_Tinteg_tracer_linrk_setup
  public :: ATMOS_DYN_Tinteg_tracer_linrk_finalize
  public :: ATMOS_DYN_Tinteg_tracer_linrk

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
  real(RP), allocatable :: QTRC_RK_list(:,:,:,:)
  integer, private  :: I_COMM_RK_list(2) 
  integer, private  :: nstage = 3

  !-----------------------------------------------------------------------------
contains

  !-----------------------------------------------------------------------------
  !> Setup
  subroutine ATMOS_DYN_Tinteg_tracer_linrk_setup( &
       tinteg_type )
    use scale_prc, only: &
       PRC_abort
    use scale_comm_cartesC, only: &
       COMM_vars8_init
    implicit none

    character(len=*) :: tinteg_type

    !---------------------------------------------------------------------------

    if ( tinteg_type(1:5) /= 'LINRK' .or. tinteg_type(7:7) /= 's' ) then
       LOG_ERROR("ATMOS_DYN_Tinteg_tracer_linrkNs_setup",*) 'TINTEG_TRACER_TYPE is invalid. Check!'
       call PRC_abort
    end if
    read(tinteg_type(6:6),*) nstage

    allocate( QTRC_RK_list(KA,IA,JA,2) )
    I_COMM_RK_list(:) = (/ 1, 2 /)
    call COMM_vars8_init( 'QTRC_RK1', QTRC_RK_list(:,:,:,1), I_COMM_RK_list(1) )
    call COMM_vars8_init( 'QTRC_RK2', QTRC_RK_list(:,:,:,2), I_COMM_RK_list(2) )

    return
  end subroutine ATMOS_DYN_Tinteg_tracer_linrk_setup

  !-----------------------------------------------------------------------------
  !> finalize
  subroutine ATMOS_DYN_Tinteg_tracer_linrk_finalize

    deallocate( QTRC_RK_list )

    return
  end subroutine ATMOS_DYN_Tinteg_tracer_linrk_finalize

  !-----------------------------------------------------------------------------
  !> linear case RK
  subroutine ATMOS_DYN_tinteg_tracer_linrk(   &
       QTRC, qflx,                            & ! (out)
       QTRC0, RHOQ_t, DENS0, DENS,            & ! (in)
       mflx_hi, num_diff,                     & ! (in)
       GSQRT, MAPF,                           & ! (in)
       CDZ, RCDZ, RCDX, RCDY,                 & ! (in)
       BND_W, BND_E, BND_S, BND_N,            & ! (in)
       TwoD,                                  & ! (in)
       dtl,                                   & ! (in)
       FLAG_FCT_TRACER, FLAG_FCT_ALONG_STREAM ) ! (in)

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
    integer :: k, i, j
    integer :: nowstage
    real(RP) :: linrk_coef
    integer :: i_in, i_out, i_tmp

    !------------------------------------------------------------------------
    ! Start RK
    !------------------------------------------------------------------------

    i_in = 1; i_out = 2

    !$omp parallel do collapse(2) private(k)
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
      QTRC_RK_list(k,i,j,i_in) = QTRC0(k,i,j)
    end do
    end do
    end do    

    do nowstage=1, nstage
      linrk_coef = 1.0_RP/real(nstage - nowstage + 1, kind=RP)
    
      call ATMOS_DYN_Copy_boundary_tracer( QTRC_RK_list(:,:,:,i_in),  & ! [INOUT]
         QTRC0,                                                       & ! [IN]
         BND_W, BND_E, BND_S, BND_N,                                  & ! [IN]
         TwoD                                                         ) ! [IN]

      call COMM_vars8( QTRC_RK_list(:,:,:,i_in), I_COMM_RK_list(i_in) )
      call COMM_wait ( QTRC_RK_list(:,:,:,i_in), I_COMM_RK_list(i_in), .false. )
            
      if (nowstage < nstage) then
         !$omp parallel do collapse(2) private(k)
         do j = JS-1, JE+1
         do i = max(IS-1,1), min(IE+1,IA)
         do k = KS, KE
            DENS_RK(k,i,j) = DENS0(k,i,j) &
                           + ( DENS(k,i,j) - DENS0(k,i,j) ) * linrk_coef
         end do
         end do
         end do

         call ATMOS_DYN_tstep_tracer( &
               QTRC_RK_list(:,:,:,i_out), qflx,          & ! (out)
               QTRC_RK_list(:,:,:,i_in), QTRC0, RHOQ_t,  & ! (in)
               DENS0, DENS_RK,                           & ! (in)
               mflx_hi, num_diff,                        & ! (in)
               GSQRT, MAPF,                              & ! (in)
               CDZ, RCDZ, RCDX, RCDY,                    & ! (in)
               TwoD, DTL*linrk_coef,                     & ! (in)
               .false., FLAG_FCT_ALONG_STREAM            ) ! (in)
                           
      else ! final stage
         call ATMOS_DYN_tstep_tracer( &
               QTRC, qflx,                              & ! (out)
               QTRC_RK_list(:,:,:,i_in), QTRC0, RHOQ_t, & ! (in)
               DENS0, DENS,                             & ! (in)
               mflx_hi, num_diff,                       & ! (in)
               GSQRT, MAPF,                             & ! (in)
               CDZ, RCDZ, RCDX, RCDY,                   & ! (in)
               TwoD, DTL,                               & ! (in)
               FLAG_FCT_TRACER, FLAG_FCT_ALONG_STREAM   ) ! (in)    
           
      end if

      i_tmp = i_out
      i_out = i_in; i_in = i_tmp
    end do
    
    return
  end subroutine ATMOS_DYN_tinteg_tracer_linrk

end module scale_atmos_dyn_tinteg_tracer_linrk
