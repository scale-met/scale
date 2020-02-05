!-------------------------------------------------------------------------------
!> module Atmosphere / Dyn Tinteg
!!
!! @par Description
!!          Temporal integration in Dynamical core for Atmospheric process
!!          7 stage and 6th order Runge-Kutta scheme
!!
!! @author Team SCALE
!!
!! This module provides a 6th order and 7 stage runge=kutta method proposed by Butcher: 
!!   y_n+1 = y_n + (11*k1 + 81*k3 + 81*k4 - 32*k5 -32*k6 + 11*k7)/120
!!  where
!!   k1 = h f(xn,yn), 
!!   k2 = h f(xn +  h/3, yn + k1/3), 
!!   k3 = h f(xn + 2h/3, yn + 2k2/3), 
!!   k4 = h f(xn +  h/3, yn +   k1/12 + k2/3  − k3/12), 
!!   k5 = h f(xn +  h/2, yn −   k1/16 + 9k2/8 − 3k3/16 − 3k4/8), 
!!   k6 = h f(xn +  h/2, yn + 9k2/8 − 3k3/8 − 3k4/4  + k5/2), and 
!!   k7 = h f(xn +  h  , yn + 9k1/44 −9k2/11 + 63k3/44 + 18k4/11 − 16k6/11). 
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
   RKUtil,                                              &
   RKUtil_setup => ATMOS_DYN_Tinteg_RKUtil_setup,       &
   RKUtil_comm  => ATMOS_DYN_Tinteg_RKUtil_comm,        &
   RKUtil_comm_wait => ATMOS_DYN_Tinteg_RKUtil_comm_wait


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
  integer, private, parameter :: RK_nregister = 8

  type(RKUtil) :: rk_dynvar
  integer, private, parameter :: I_RK_DENS = 1
  integer, private, parameter :: I_RK_MOMZ = 2
  integer, private, parameter :: I_RK_MOMX = 3
  integer, private, parameter :: I_RK_MOMY = 4
  integer, private, parameter :: I_RK_RHOT = 5

  type(RKUtil) :: rk_prgvar

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
    integer :: stage
    character(len=H_SHORT) :: stage_s

    character(H_MID) :: dynvar_name_list(5)
    character(H_MID) :: prgvar_name_list(VA)
    !---------------------------------------------------------------------------

    if ( tinteg_type /= 'RK7s6o' ) then
       LOG_ERROR("ATMOS_DYN_Tinteg_short_rk6_setup",*) 'TINTEG_TYPE is not RK4. Check!'
       call PRC_abort
    end if

    dynvar_name_list(1) = 'DENS'
    dynvar_name_list(2) = 'MOMZ'
    dynvar_name_list(3) = 'MOMX'
    dynvar_name_list(4) = 'MOMY'
    dynvar_name_list(5) = 'RHOT'
    call RKUtil_setup( rk_dynvar, RK_nregister, dynvar_name_list, 0 )

    do iv = 1, VA
      prgvar_name_list(iv) = 'PROG'
    end do
    call RKUtil_setup( rk_prgvar, RK_nregister, prgvar_name_list, 5 )
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
       BND_W, BND_E, BND_S, BND_N,              &
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

    real(RP), intent(in)    :: dt

    real(RP) :: mflx_hi_RK(KA,IA,JA,3,RK_nstage)
    real(RP) :: tflx_hi_RK(KA,IA,JA,3,RK_nstage)

    real(RP) :: dtrk(RK_nstage)

    integer  :: i, j, k, iv, n, s
    integer :: stage
    integer :: astage

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
    rk_dynvar %tmp(:,:,:,:,:) = UNDEF
    if (VA > 0) rk_prgvar%tmp(:,:,:,:,:) = UNDEF
    mflx_hi_RK(:,:,:,:,:) = UNDEF
    tflx_hi_RK(:,:,:,:,:) = UNDEF
    !$omp end parallel workshare
#endif

#ifdef QUICKDEBUG
    mflx_hi(   1:KS-1,:,:,:) = UNDEF
    mflx_hi(KE+1:KA  ,:,:,:) = UNDEF
#endif

!OCL XFILL
    rk_dynvar%tmp(:,:,:,I_RK_DENS,1) = DENS
!OCL XFILL
    rk_dynvar%tmp(:,:,:,I_RK_MOMZ,1) = MOMZ
!OCL XFILL
    rk_dynvar%tmp(:,:,:,I_RK_MOMX,1) = MOMX
!OCL XFILL
    rk_dynvar%tmp(:,:,:,I_RK_MOMY,1) = MOMY
!OCL XFILL
    rk_dynvar%tmp(:,:,:,I_RK_RHOT,1) = RHOT
!OCL XFILL
    if (VA > 0) rk_prgvar%tmp(:,:,:,:,1) = PROG

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

    dtrk(:) = (/ dt/3.0_RP, dt * 2.0_RP / 3.0_RP, 1.0_RP, 1.0_RP, 1.0_RP, 1.0_RP, 1.0_RP /)
    io_prgvar(:) = 0; jo_prgvar(:) = 0; ko_prgvar(:) = 0

    call PROF_rapend  ("DYN_RK7s6o_Prep",3)

    !------------------------------------------------------------------------
    ! Start RK
    !------------------------------------------------------------------------

    call PROF_rapstart("DYN_RK7s6o",3)
    call ATMOS_DYN_tstep( rk_dynvar%tmp(:,:,:,I_RK_DENS,2),                 & ! [OUT] 
                          rk_dynvar%tmp(:,:,:,I_RK_MOMZ,2),                 & ! [OUT]
                          rk_dynvar%tmp(:,:,:,I_RK_MOMX,2),                 & ! [OUT]
                          rk_dynvar%tmp(:,:,:,I_RK_MOMY,2),                 & ! [OUT]
                          rk_dynvar%tmp(:,:,:,I_RK_RHOT,2),                 & ! [OUT]
                          rk_prgvar%tmp(:,:,:,:       ,2),                  & ! [OUT]
                          mflx_hi_RK(:,:,:,:,1), tflx_hi_RK(:,:,:,:,1),     & ! [INOUT,OUT]
                          rk_dynvar%tmp(:,:,:,I_RK_DENS,1),                 & ! [OUT] 
                          rk_dynvar%tmp(:,:,:,I_RK_MOMZ,1),                 & ! [OUT]
                          rk_dynvar%tmp(:,:,:,I_RK_MOMX,1),                 & ! [OUT]
                          rk_dynvar%tmp(:,:,:,I_RK_MOMY,1),                 & ! [OUT]
                          rk_dynvar%tmp(:,:,:,I_RK_RHOT,1),                 & ! [OUT]
                          DENS,     MOMZ,     MOMX,     MOMY,     RHOT,     & ! [IN]
                          DENS_t,   MOMZ_t,   MOMX_t,   MOMY_t,   RHOT_t,   & ! [IN]
                          rk_prgvar%tmp(:,:,:,:,1), PROG,                   & ! [IN]
                          DPRES0, CVtot, CORIOLI,                           & ! [IN]
                          num_diff, wdamp_coef, divdmp_coef, DDIV,          & ! [IN]
                          FLAG_FCT_MOMENTUM, FLAG_FCT_T,                    & ! [IN]
                          FLAG_FCT_ALONG_STREAM,                            & ! [IN]
                          CDZ, FDZ, FDX, FDY,                               & ! [IN]
                          RCDZ, RCDX, RCDY, RFDZ, RFDX, RFDY,               & ! [IN]
                          PHI, GSQRT, J13G, J23G, J33G, MAPF,               & ! [IN]
                          REF_pres, REF_dens,                               & ! [IN]
                          BND_W, BND_E, BND_S, BND_N,                       & ! [IN]
                          dtrk(1), .false.                                  ) ! [IN]
    call PROF_rapend  ("DYN_RK7s6o",3)

    !-

    do stage = 2, RK_nstage

      astage = stage + 1

      call PROF_rapstart("DYN_RK7s6o_BND",3)
      call ATMOS_DYN_Copy_boundary(  rk_dynvar%tmp(:,:,:,I_RK_DENS,stage),                 & ! [INOUT] 
                                     rk_dynvar%tmp(:,:,:,I_RK_MOMZ,stage),                 & ! [INOUT]
                                     rk_dynvar%tmp(:,:,:,I_RK_MOMX,stage),                 & ! [INOUT]
                                     rk_dynvar%tmp(:,:,:,I_RK_MOMY,stage),                 & ! [INOUT]
                                     rk_dynvar%tmp(:,:,:,I_RK_RHOT,stage),                 & ! [INOUT]
                                     rk_prgvar%tmp(:,:,:,:        ,stage),                 & ! [INOUT]
                                     rk_dynvar%tmp(:,:,:,I_RK_DENS,    1),                 & ! [IN] 
                                     rk_dynvar%tmp(:,:,:,I_RK_MOMZ,    1),                 & ! [IN]
                                     rk_dynvar%tmp(:,:,:,I_RK_MOMX,    1),                 & ! [IN]
                                     rk_dynvar%tmp(:,:,:,I_RK_MOMY,    1),                 & ! [IN]
                                     rk_dynvar%tmp(:,:,:,I_RK_RHOT,    1),                 & ! [IN]
                                     rk_prgvar%tmp(:,:,:,:        ,    1),                 & ! [IN]
                                     BND_W, BND_E, BND_S, BND_N                            ) ! [IN]
  
      call PROF_rapend  ("DYN_RK7s6o_BND",3)

      call RKUtil_comm( rk_dynvar, stage )
      call RKUtil_comm( rk_prgvar, stage )
      call RKUtil_comm_wait( rk_dynvar, stage )
      call RKUtil_comm_wait( rk_prgvar, stage )

      !--
      call PROF_rapstart("DYN_RK7s6o",3)

      call ATMOS_DYN_tstep( &
         rk_dynvar%tmp(:,:,:,I_RK_DENS,astage),                 & ! [OUT] 
         rk_dynvar%tmp(:,:,:,I_RK_MOMZ,astage),                 & ! [OUT]
         rk_dynvar%tmp(:,:,:,I_RK_MOMX,astage),                 & ! [OUT]
         rk_dynvar%tmp(:,:,:,I_RK_MOMY,astage),                 & ! [OUT]
         rk_dynvar%tmp(:,:,:,I_RK_RHOT,astage),                 & ! [OUT]
         rk_prgvar%tmp(:,:,:,:       ,astage),                  & ! [OUT]
         mflx_hi_RK(:,:,:,:,stage), tflx_hi_RK(:,:,:,:,stage),  & ! [INOUT,OUT]
         rk_dynvar%tmp(:,:,:,I_RK_DENS,1),                      & ! [IN] 
         rk_dynvar%tmp(:,:,:,I_RK_MOMZ,1),                      & ! [IN]
         rk_dynvar%tmp(:,:,:,I_RK_MOMX,1),                      & ! [IN]
         rk_dynvar%tmp(:,:,:,I_RK_MOMY,1),                      & ! [IN]
         rk_dynvar%tmp(:,:,:,I_RK_RHOT,1),                      & ! [IN
         rk_dynvar%tmp(:,:,:,I_RK_DENS,stage),                  & ! [IN] 
         rk_dynvar%tmp(:,:,:,I_RK_MOMZ,stage),                  & ! [IN]
         rk_dynvar%tmp(:,:,:,I_RK_MOMX,stage),                  & ! [IN]
         rk_dynvar%tmp(:,:,:,I_RK_MOMY,stage),                  & ! [IN]
         rk_dynvar%tmp(:,:,:,I_RK_RHOT,stage),                  & ! [IN]
         DENS_t,   MOMZ_t,   MOMX_t,   MOMY_t,   RHOT_t,        & ! [IN]
         rk_prgvar%tmp(:,:,:,:,1),                              & ! [IN]
         rk_prgvar%tmp(:,:,:,:,stage),                          & ! [IN]
         DPRES0, CVtot, CORIOLI,                                & ! [IN]
         num_diff, wdamp_coef, divdmp_coef, DDIV,               & ! [IN]
         FLAG_FCT_MOMENTUM, FLAG_FCT_T,                         & ! [IN]
         FLAG_FCT_ALONG_STREAM,                                 & ! [IN]
         CDZ, FDZ, FDX, FDY,                                    & ! [IN]
         RCDZ, RCDX, RCDY, RFDZ, RFDX, RFDY,                    & ! [IN]
         PHI, GSQRT, J13G, J23G, J33G, MAPF,                    & ! [IN]
         REF_pres, REF_dens,                                    & ! [IN]
         BND_W, BND_E, BND_S, BND_N,                            & ! [IN]
         dtrk(stage), .false.                                   ) ! [IN]
   
      
      call calc_var_nextstage( stage, rk_dynvar%tmp(:,:,:,:,:), io_dynvar(:), jo_dynvar(:), ko_dynvar(:),  5, dt )
      call calc_var_nextstage( stage, rk_prgvar%tmp(:,:,:,:,:), io_prgvar(:), jo_prgvar(:), ko_prgvar(:), VA, dt )
      call PROF_rapend  ("DYN_RK7s6o",3)    

    end do

    call PROF_rapstart("DYN_RK7s6o",3)

    call calc_var_nextstep( DENS, rk_dynvar%tmp(:,:,:,I_RK_DENS,:), io_dynvar(I_RK_DENS), jo_dynvar(I_RK_DENS), ko_dynvar(I_RK_DENS),  1, dt )
    call calc_var_nextstep( MOMZ, rk_dynvar%tmp(:,:,:,I_RK_MOMZ,:), io_dynvar(I_RK_MOMZ), jo_dynvar(I_RK_MOMZ), ko_dynvar(I_RK_MOMZ),  2, dt )
    call calc_var_nextstep( MOMX, rk_dynvar%tmp(:,:,:,I_RK_MOMX,:), io_dynvar(I_RK_MOMX), jo_dynvar(I_RK_MOMX), ko_dynvar(I_RK_MOMX),  3, dt )
    call calc_var_nextstep( MOMY, rk_dynvar%tmp(:,:,:,I_RK_MOMY,:), io_dynvar(I_RK_MOMY), jo_dynvar(I_RK_MOMY), ko_dynvar(I_RK_MOMY),  4, dt )
    call calc_var_nextstep( RHOT, rk_dynvar%tmp(:,:,:,I_RK_RHOT,:), io_dynvar(I_RK_RHOT), jo_dynvar(I_RK_RHOT), ko_dynvar(I_RK_RHOT),  5, dt )
    call calc_var_nextstep( PROG(:,:,:,:), rk_prgvar%tmp(:,:,:,:,:), 0, 0, 0, VA, dt )
    call calc_flux_nextstep( mflx_hi, mflx_hi_RK ) 
    call calc_flux_nextstep( tflx_hi, tflx_hi_RK ) 

    call PROF_rapend("DYN_RK7s6o",3)

    return
  end subroutine ATMOS_DYN_tinteg_short_rk7s6o

  !--- private ---------------------------------------------------------------------------

  subroutine calc_var_nextstage( nowstage, var_rkwork, io, jo, ko, va_, dt )

   implicit none
   integer, intent(in) :: va_
   integer, intent(in) :: nowstage
   real(RP), intent(inout) :: var_rkwork(KA,IA,JA,va_,RK_nregister)
   integer, intent(in) :: io(va_), jo(va_), ko(va_)
   real(RP), intent(in) :: dt

   integer :: i, j, k, iv
   real(RP) :: var0
   !--------------------------------------

   select case(nowstage)
   case(3)
      !$omp parallel do private(iv,k,j,i,var0) collapse(3)
      do iv=1, va_
      do k=KS-ko(iv), KE
      do j=JS-jo(iv), JE
      do i=IS-io(iv), IE
         var0 = var_rkwork(i,j,k,iv,1)
         var_rkwork(i,j,k,iv,2) = (var_rkwork(i,j,k,iv,2) - var0) * 3.0_RP / dt
         var_rkwork(i,j,k,iv,3) = (var_rkwork(i,j,k,iv,3) - var0) * 3.0_RP / (2.0_RP * dt)
         var_rkwork(i,j,k,iv,4) =  var_rkwork(i,j,k,iv,4) - var0
   
         var_rkwork(i,j,k,iv,5) = var0 + dt * (             &
                    1.0_RP/12.0_RP * var_rkwork(i,j,k,iv,2) &
                  + 1.0_RP/ 3.0_RP * var_rkwork(i,j,k,iv,3) &
                  - 1.0_RP/12.0_RP * var_rkwork(i,j,k,iv,4) )
      end do
      end do
      end do
      end do
   case(4)
      !$omp parallel do private(iv,k,j,i,var0) collapse(3)
      do iv=1, va_
      do k=KS-ko(iv), KE
      do j=JS-jo(iv), JE
      do i=IS-io(iv), IE
         var0 = var_rkwork(i,j,k,iv,1)
         var_rkwork(i,j,k,iv,5) =  var_rkwork(i,j,k,iv,5) - var0

         var_rkwork(i,j,k,iv,6) = var0 + dt * (             &
                  - 1.0_RP/16.0_RP * var_rkwork(i,j,k,iv,2) &
                  + 9.0_RP/ 8.0_RP * var_rkwork(i,j,k,iv,3) &
                  - 3.0_RP/16.0_RP * var_rkwork(i,j,k,iv,4) &
                  - 3.0_RP/ 8.0_RP * var_rkwork(i,j,k,iv,5) )
      end do
      end do
      end do
      end do
   case(5)
      !$omp parallel do private(iv,k,j,i,var0) collapse(3)
      do iv=1, va_
      do k=KS-ko(iv), KE
      do j=JS-jo(iv), JE
      do i=IS-io(iv), IE
         var0 = var_rkwork(i,j,k,iv,1)
         var_rkwork(i,j,k,iv,6) =  var_rkwork(i,j,k,iv,6) - var0
   
         var_rkwork(i,j,k,iv,7) = var0 + dt * (            &
                  + 9.0_RP/8.0_RP * var_rkwork(i,j,k,iv,3) &
                  - 3.0_RP/8.0_RP * var_rkwork(i,j,k,iv,4) &
                  - 3.0_RP/4.0_RP * var_rkwork(i,j,k,iv,5) &
                  + 1.0_RP/2.0_RP * var_rkwork(i,j,k,iv,6) )
      end do
      end do
      end do
      end do
   case(6)
      !$omp parallel do private(iv,k,j,i,var0) collapse(3)
      do iv=1, va_
      do k=KS-ko(iv), KE
      do j=JS-jo(iv), JE
      do i=IS-io(iv), IE
         var0 = var_rkwork(i,j,k,iv,1)
         var_rkwork(i,j,k,iv,7) =  var_rkwork(i,j,k,iv,7) - var0
   
         var_rkwork(i,j,k,iv,8) = var0 + dt * (              &
                     9.0_RP/44.0_RP * var_rkwork(i,j,k,iv,2) &
                  -  9.0_RP/11.0_RP * var_rkwork(i,j,k,iv,3) &
                  + 63.0_RP/44.0_RP * var_rkwork(i,j,k,iv,4) &
                  - 18.0_RP/11.0_RP * var_rkwork(i,j,k,iv,5) &
                  - 16.0_RP/11.0_RP * var_rkwork(i,j,k,iv,7) )
      end do
      end do
      end do
      end do
   end select

   return
  end subroutine calc_var_nextstage

 subroutine calc_Var_nextstep( var, var_rkwork, io, jo, ko, va_, dt ) 
   implicit none
   integer, intent(in) :: va_
   real(RP), intent(inout) :: var(KA,IA,JA, va_)
   real(RP), intent(in) :: var_rkwork(KA,IA,JA,va_,RK_nregister)
   integer, intent(in) :: io, jo, ko
   real(RP), intent(in) :: dt

   integer :: i, j, k, iv
   real(RP) :: var0
   !--------------------------------------

   !$omp parallel do private(iv,k,j,i,var0) collapse(3)
   do iv=1, va_
   do k=KS-ko, KE
   do j=JS-jo, JE
   do i=IS-io, IE
     var0 = var_rkwork(i,j,k,iv,1)
     var(i,j,k,iv) = var0 + dt/120.0_RP * (               &
                11.0_RP *  var_rkwork(i,j,k,iv,2)         &
              + 81.0_RP *  var_rkwork(i,j,k,iv,4)         &
              + 81.0_RP *  var_rkwork(i,j,k,iv,5)         &
              - 32.0_RP *  var_rkwork(i,j,k,iv,6)         &
              - 32.0_RP *  var_rkwork(i,j,k,iv,7)         &
              + 11.0_RP * (var_rkwork(i,j,k,iv,8) - var0) )
   end do
   end do
   end do
   end do

   return
 end subroutine calc_Var_nextstep 

 subroutine calc_flux_nextstep( flux, flux_rkwork ) 
   implicit none
   real(RP), intent(inout) :: flux(KA,IA,JA,3)
   real(RP), intent(in) :: flux_rkwork(KA,IA,JA,3,RK_nstage)

   integer :: i, j, k, d

   !--------------------------------------

   !$omp parallel do private(d,k,j,i) collapse(3)
   do d=1, 3
   do k=KS, KE
   do j=JS, JE
   do i=IS, IE
     flux(i,j,k,d) = flux(i,j,k,d) + 1.0_RP/120.0_RP * (            &
                11.0_RP *  flux_rkwork(i,j,k,d,1)                   &
              + 81.0_RP *  flux_rkwork(i,j,k,d,3)                   &
              + 81.0_RP *  flux_rkwork(i,j,k,d,4)                   &
              - 32.0_RP *  flux_rkwork(i,j,k,d,5)                   &
              - 32.0_RP *  flux_rkwork(i,j,k,d,6)                   &
              + 11.0_RP * (flux_rkwork(i,j,k,d,7) - flux(i,j,k,d))  )
   end do
   end do
   end do
   end do

   return
 end subroutine calc_flux_nextstep 

end module scale_atmos_dyn_tinteg_short_rk7s6o
