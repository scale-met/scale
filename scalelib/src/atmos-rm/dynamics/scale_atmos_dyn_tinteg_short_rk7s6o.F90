!-------------------------------------------------------------------------------
!> module Atmosphere / Dyn Tinteg
!!
!! @par Description
!!          Temporal integration in Dynamical core for Atmospheric process
!!          7 stage Runge-Kutta scheme with 6th order accuracy
!!
!! @author Team SCALE
!!
!! This module provides a 7 stage and 6th order runge=kutta method with extended region of stability proposed by Lawson (1967)
!! See scale_atmos_dyn_tinteg_rkutil.F90 for the detail of RK coffecients. 
!<
!-------------------------------------------------------------------------------
#include "scalelib.h"
module scale_atmos_dyn_tinteg_short_rk7s6o
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
#if defined DEBUG || defined QUICKDEBUG
  use scale_debug, only: &
     CHECK
  use scale_const, only: &
     UNDEF  => CONST_UNDEF, &
     IUNDEF => CONST_UNDEF2
#endif

  use scale_atmos_dyn_tinteg_rkutil, only: &
   RKUtil,                                                          &
   RKUtil_setup => ATMOS_DYN_Tinteg_RKUtil_setup,                   &
   RKUtil_rkwork_alloc   => ATMOS_DYN_Tinteg_RKUtil_rkwork_alloc,   &
   RKUtil_rkwork_dealloc => ATMOS_DYN_Tinteg_RKUtil_rkwork_dealloc, &
   RKUtil_comm  => ATMOS_DYN_Tinteg_RKUtil_comm,                    &
   RKUtil_comm_wait => ATMOS_DYN_Tinteg_RKUtil_comm_wait,           &
   RKUtil_nextstage => ATMOS_DYN_Tinteg_RKUtil_7s6o_nextstage,      &
   RKUtil_updateVar => ATMOS_DYN_Tinteg_RKUtil_7s6o_updateVar,      &
   RKUtil_updateFlux => ATMOS_DYN_Tinteg_RKUtil_7s6o_updateFlux


  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: ATMOS_DYN_Tinteg_short_rk7s6o_setup
  public :: ATMOS_DYN_Tinteg_short_rk7s6o

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
  integer, private, parameter :: RK_nstage    = 7
  integer, private, parameter :: RK_nregister = 7

  type(RKUtil), private :: rk_dynvar
  integer, private, parameter :: I_RK_DENS = 1
  integer, private, parameter :: I_RK_MOMZ = 2
  integer, private, parameter :: I_RK_MOMX = 3
  integer, private, parameter :: I_RK_MOMY = 4
  integer, private, parameter :: I_RK_RHOT = 5

  type(RKUtil), private :: rk_prgvar

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine ATMOS_DYN_Tinteg_short_rk7s6o_setup( &
       tinteg_type )
    use scale_prc, only: &
       PRC_abort
    use scale_const, only: &
       UNDEF => CONST_UNDEF
    use scale_comm_cartesC, only: &
       COMM_vars8_init
    implicit none

    character(len=*) :: tinteg_type

    integer :: iv
    character(H_SHORT) :: dynvar_name_list(5)
    character(H_MID) :: prgvar_name_list(VA)
    !---------------------------------------------------------------------------

    if ( tinteg_type /= 'RK7s6o' ) then
       LOG_ERROR("ATMOS_DYN_Tinteg_short_rk7s6o_setup",*) 'TINTEG_TYPE is not RK7s6o. Check!'
       call PRC_abort
    end if

    dynvar_name_list(1) = 'DENS'
    dynvar_name_list(2) = 'MOMZ'
    dynvar_name_list(3) = 'MOMX'
    dynvar_name_list(4) = 'MOMY'
    dynvar_name_list(5) = 'RHOT'
    call RKUtil_setup( rk_dynvar, RK_nstage, RK_nregister, dynvar_name_list, 0, .false. )

    do iv = 1, VA
      prgvar_name_list(iv) = 'PROG'
    end do
    call RKUtil_setup( rk_prgvar, RK_nstage, RK_nregister, prgvar_name_list, 5, .false. )

    !----------------------------------------------------

    return
  end subroutine ATMOS_DYN_Tinteg_short_rk7s6o_setup

  !-----------------------------------------------------------------------------
  !> RK7s6o
  subroutine ATMOS_DYN_tinteg_short_rk7s6o( &
       DENS, MOMZ, MOMX, MOMY, RHOT, PROG,      &
       mflx_hi,  tflx_hi,                       &
       DENS_t, MOMZ_t, MOMX_t, MOMY_t, RHOT_t,  &
       DPRES0, CVtot, CORIOLI,                  &
       num_diff, wdamp_coef, divdmp_coef, DDIV, &
       FLAG_FCT_MOMENTUM, FLAG_FCT_T,           &
       FLAG_FCT_ALONG_STREAM,                   &
       CDZ, FDZ, FDX, FDY,                      &
       RCDZ, RCDX, RCDY, RFDZ, RFDX, RFDY,      &
       PHI, GSQRT, J13G, J23G, J33G, MAPF,      &
       REF_pres, REF_dens,                      &
       BND_W, BND_E, BND_S, BND_N, TwoD,        &
       dt                                       )
    use scale_comm_cartesC, only: &
       COMM_vars8, &
       COMM_wait
    use scale_atmos_dyn_tstep_short, only: &
       ATMOS_DYN_tstep => ATMOS_DYN_Tstep_short
    use scale_atmos_dyn_common, only: &
       ATMOS_DYN_Copy_boundary
    implicit none

    real(RP), intent(inout) :: DENS(KA,IA,JA)
    real(RP), intent(inout) :: MOMZ(KA,IA,JA)
    real(RP), intent(inout) :: MOMX(KA,IA,JA)
    real(RP), intent(inout) :: MOMY(KA,IA,JA)
    real(RP), intent(inout) :: RHOT(KA,IA,JA)
    real(RP), intent(inout) :: PROG(KA,IA,JA,VA)

    real(RP), intent(inout) :: mflx_hi(KA,IA,JA,3)
    real(RP), intent(inout) :: tflx_hi(KA,IA,JA,3)

    real(RP), intent(in)    :: DENS_t(KA,IA,JA)
    real(RP), intent(in)    :: MOMZ_t(KA,IA,JA)
    real(RP), intent(in)    :: MOMX_t(KA,IA,JA)
    real(RP), intent(in)    :: MOMY_t(KA,IA,JA)
    real(RP), intent(in)    :: RHOT_t(KA,IA,JA)

    real(RP), intent(in)    :: DPRES0(KA,IA,JA)
    real(RP), intent(in)    :: CVtot(KA,IA,JA)
    real(RP), intent(in)    :: CORIOLI(IA,JA)
    real(RP), intent(in)    :: num_diff(KA,IA,JA,5,3)
    real(RP), intent(in)    :: wdamp_coef(KA)
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
    logical,  intent(in)    :: TwoD

    real(RP), intent(in)    :: dt

    real(RP) :: mflx_hi_RK(KA,IA,JA,3,RK_nstage)
    real(RP) :: tflx_hi_RK(KA,IA,JA,3,RK_nstage)

    integer  :: i, j, k, iv, n, s
    integer :: stage

    integer, parameter :: ko_dynvar(5) = (/ 0, 1, 0, 0, 0 /)
    integer, parameter :: io_dynvar(5) = (/ 0, 0, 1, 0, 0 /)
    integer, parameter :: jo_dynvar(5) = (/ 0, 0, 0, 1, 0 /)
    integer :: ko_prgvar(VA)
    integer :: io_prgvar(VA)
    integer :: jo_prgvar(VA)

    !---------------------------------------------------------------------------

    call PROF_rapstart("DYN_RK7s6o_Prep",3)

#ifdef DEBUG
    !$omp parallel workshare
    mflx_hi_RK(:,:,:,:,:) = UNDEF
    tflx_hi_RK(:,:,:,:,:) = UNDEF
    !$omp end parallel workshare
#endif

#ifdef QUICKDEBUG
    mflx_hi(   1:KS-1,:,:,:) = UNDEF
    mflx_hi(KE+1:KA  ,:,:,:) = UNDEF
#endif

    call RKUtil_rkwork_alloc( rk_dynvar )
    call RKUtil_rkwork_alloc( rk_prgvar )

    !$omp parallel 
   
    !$omp workshare
!OCL XFILL
    rk_dynvar%work0(:,:,:,I_RK_DENS) = DENS(:,:,:)  
!OCL XFILL
    rk_dynvar%work0(:,:,:,I_RK_MOMZ) = MOMZ(:,:,:) 
!OCL XFILL
    rk_dynvar%work0(:,:,:,I_RK_MOMX) = MOMX(:,:,:) 
!OCL XFILL
    rk_dynvar%work0(:,:,:,I_RK_MOMY) = MOMY(:,:,:) 
!OCL XFILL
    rk_dynvar%work0(:,:,:,I_RK_RHOT) = RHOT(:,:,:)
!OCL XFILL
    rk_dynvar%buf(:,:,:,:) = rk_dynvar%work0(:,:,:,:)
    !$omp end workshare

    if (VA > 0) then
      !$omp workshare
!OCL XFILL      
      rk_prgvar%work0(:,:,:,:) = PROG
!OCL XFILL
      rk_prgvar%buf(:,:,:,:) = rk_prgvar%work0(:,:,:,:)
      !$omp end workshare
    end if
    !$omp end parallel

    if ( BND_W ) then
       do j = JS, JE
       do k = KS, KE
          mflx_hi_RK(k,IS-1,j,2,:) = mflx_hi(k,IS-1,j,2)
       end do
       end do
    end if
    if ( BND_E ) then
       do j = JS, JE
       do k = KS, KE
          mflx_hi_RK(k,IE,j,2,:) = mflx_hi(k,IE,j,2)
       end do
       end do
    end if
    if ( BND_S ) then
       do i = IS, IE
       do k = KS, KE
          mflx_hi_RK(k,i,JS-1,3,:) = mflx_hi(k,i,JS-1,3)
       end do
       end do
    end if
    if ( BND_N ) then
       do i = IS, IE
       do k = KS, KE
          mflx_hi_RK(k,i,JE,3,:) = mflx_hi(k,i,JE,3)
       end do
       end do
    end if

    io_prgvar(:) = 0; jo_prgvar(:) = 0; ko_prgvar(:) = 0

    call PROF_rapend  ("DYN_RK7s6o_Prep",3)

    !------------------------------------------------------------------------
    ! Start RK
    !------------------------------------------------------------------------

    do stage = 1, RK_nstage

      if ( stage > 1) then
         call PROF_rapstart("DYN_RK7s6o_BND",3)
         call ATMOS_DYN_Copy_boundary( rk_dynvar%buf(:,:,:,I_RK_DENS),      & ! [INOUT] 
                                       rk_dynvar%buf(:,:,:,I_RK_MOMZ),      & ! [INOUT]
                                       rk_dynvar%buf(:,:,:,I_RK_MOMX),      & ! [INOUT]
                                       rk_dynvar%buf(:,:,:,I_RK_MOMY),      & ! [INOUT]
                                       rk_dynvar%buf(:,:,:,I_RK_RHOT),      & ! [INOUT]
                                       rk_prgvar%buf(:,:,:,:        ),      & ! [INOUT]
                                       rk_dynvar%work0(:,:,:,I_RK_DENS),    & ! [IN] 
                                       rk_dynvar%work0(:,:,:,I_RK_MOMZ),    & ! [IN] 
                                       rk_dynvar%work0(:,:,:,I_RK_MOMX),    & ! [IN] 
                                       rk_dynvar%work0(:,:,:,I_RK_MOMY),    & ! [IN] 
                                       rk_dynvar%work0(:,:,:,I_RK_RHOT),    & ! [IN] 
                                       rk_prgvar%work0(:,:,:,:),            & ! [IN]
                                       BND_W, BND_E, BND_S, BND_N, TwoD     ) ! [IN]
         call PROF_rapend  ("DYN_RK7s6o_BND",3)

         call RKUtil_comm( rk_dynvar )
         call RKUtil_comm( rk_prgvar )
         call RKUtil_comm_wait( rk_dynvar )
         call RKUtil_comm_wait( rk_prgvar )
      end if
      !--

      call PROF_rapstart("DYN_RK7s6o",3)

      call ATMOS_DYN_tstep( &
         rk_dynvar%work(:,:,:,I_RK_DENS,stage),                & ! [OUT] 
         rk_dynvar%work(:,:,:,I_RK_MOMZ,stage),                & ! [OUT]
         rk_dynvar%work(:,:,:,I_RK_MOMX,stage),                & ! [OUT]
         rk_dynvar%work(:,:,:,I_RK_MOMY,stage),                & ! [OUT]
         rk_dynvar%work(:,:,:,I_RK_RHOT,stage),                & ! [OUT]
         rk_prgvar%work(:,:,:,:        ,stage),                & ! [OUT]
         mflx_hi_RK(:,:,:,:,stage), tflx_hi_RK(:,:,:,:,stage), & ! [INOUT,OUT]
         rk_dynvar%work0(:,:,:,I_RK_DENS),                     & ! [IN]
         rk_dynvar%work0(:,:,:,I_RK_MOMZ),                     & ! [IN]
         rk_dynvar%work0(:,:,:,I_RK_MOMX),                     & ! [IN]
         rk_dynvar%work0(:,:,:,I_RK_MOMY),                     & ! [IN]
         rk_dynvar%work0(:,:,:,I_RK_RHOT),                     & ! [IN]
         rk_dynvar%buf(:,:,:,I_RK_DENS),                       & ! [IN] 
         rk_dynvar%buf(:,:,:,I_RK_MOMZ),                       & ! [IN]
         rk_dynvar%buf(:,:,:,I_RK_MOMX),                       & ! [IN]
         rk_dynvar%buf(:,:,:,I_RK_MOMY),                       & ! [IN]
         rk_dynvar%buf(:,:,:,I_RK_RHOT),                       & ! [IN]
         DENS_t,   MOMZ_t,   MOMX_t,   MOMY_t,   RHOT_t,       & ! [IN]
         rk_prgvar%work0, rk_prgvar%buf(:,:,:,:),              & ! [IN]
         DPRES0, CVtot, CORIOLI,                               & ! [IN]
         num_diff, wdamp_coef, divdmp_coef, DDIV,              & ! [IN]
         FLAG_FCT_MOMENTUM, FLAG_FCT_T,                        & ! [IN]
         FLAG_FCT_ALONG_STREAM,                                & ! [IN]
         CDZ, FDZ, FDX, FDY,                                   & ! [IN]
         RCDZ, RCDX, RCDY, RFDZ, RFDX, RFDY,                   & ! [IN]
         PHI, GSQRT, J13G, J23G, J33G, MAPF,                   & ! [IN]
         REF_pres, REF_dens,                                   & ! [IN]
         BND_W, BND_E, BND_S, BND_N, TwoD,                     & ! [IN]
         1.0_RP, .false.                                       ) ! [IN]
   
      if ( stage < RK_nstage) then
         call RKUtil_nextstage( rk_dynvar, stage, io_dynvar, jo_dynvar, ko_dynvar, dt )
         call RKUtil_nextstage( rk_prgvar, stage, io_prgvar, jo_prgvar, ko_prgvar, dt )
      end if

      call PROF_rapend  ("DYN_RK7s6o",3)
    end do

    call PROF_rapstart("DYN_RK7s6o",3)
    call RKUtil_updateVar( DENS, rk_dynvar, io_dynvar, jo_dynvar, ko_dynvar, I_RK_DENS, I_RK_DENS, dt )
    call RKUtil_updateVar( MOMZ, rk_dynvar, io_dynvar, jo_dynvar, ko_dynvar, I_RK_MOMZ, I_RK_MOMZ, dt )
    call RKUtil_updateVar( MOMX, rk_dynvar, io_dynvar, jo_dynvar, ko_dynvar, I_RK_MOMX, I_RK_MOMX, dt )
    call RKUtil_updateVar( MOMY, rk_dynvar, io_dynvar, jo_dynvar, ko_dynvar, I_RK_MOMY, I_RK_MOMY, dt )
    call RKUtil_updateVar( RHOT, rk_dynvar, io_dynvar, jo_dynvar, ko_dynvar, I_RK_RHOT, I_RK_RHOT, dt )
    call RKUtil_updateVar( PROG, rk_prgvar, io_prgvar, jo_prgvar, ko_prgvar, 1, VA, dt )
    call RKUtil_updateFlux( mflx_hi, mflx_hi_RK, 0, 0, 0, 3 ) 
    call RKUtil_updateFlux( tflx_hi, tflx_hi_RK, 0, 0, 0, 3 )
    call PROF_rapend("DYN_RK7s6o",3)

    call RKUtil_rkwork_dealloc( rk_dynvar )
    call RKUtil_rkwork_dealloc( rk_prgvar )

    return
  end subroutine ATMOS_DYN_tinteg_short_rk7s6o

end module scale_atmos_dyn_tinteg_short_rk7s6o
