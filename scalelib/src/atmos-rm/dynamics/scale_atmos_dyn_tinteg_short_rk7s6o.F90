!-------------------------------------------------------------------------------
!> module Atmosphere / Dyn Tinteg
!!
!! @par Description
!!          Temporal integration in Dynamical core for Atmospheric process
!!          7 stage and 6th order Runge-Kutta scheme
!!
!! @author Team SCALE
!!
!! This module provides a 7 stage  and 6th order runge=kutta method proposed by Butcher (1964): 
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
  integer, private, parameter :: RK_nregister = 7

  type(RKUtil) :: rk_dynvar
  integer, private, parameter :: I_RK_DENS = 1
  integer, private, parameter :: I_RK_MOMZ = 2
  integer, private, parameter :: I_RK_MOMX = 3
  integer, private, parameter :: I_RK_MOMY = 4
  integer, private, parameter :: I_RK_RHOT = 5

  type(RKUtil) :: rk_prgvar

  real(RP), parameter :: RKCoef_a(RK_nstage,RK_nstage) = reshape( &
   (/ 0.0_RP, 1.0_RP/3.0_RP,        0.0_RP,  1.0_RP/12.0_RP, -1.0_RP/16.0_RP,         0.0_RP,   9.0_RP/44.0_RP,    &
      0.0_RP,        0.0_RP, 2.0_RP/3.0_RP,   1.0_RP/3.0_RP,   9.0_RP/8.0_RP,  9.0_RP/8.0_RP,  -9.0_RP/11.0_RP,    &
      0.0_RP,        0.0_RP,        0.0_RP, -1.0_RP/12.0_RP, -3.0_RP/16.0_RP,   -3_RP/8.0_RP,  63.0_RP/44.0_RP,    &
      0.0_RP,        0.0_RP,        0.0_RP,          0.0_RP,  -3.0_RP/8.0_RP, -3.0_RP/4.0_RP,  18.0_RP/11.0_RP,    &
      0.0_RP,        0.0_RP,        0.0_RP,          0.0_RP,          0.0_RP,  1.0_RP/2.0_RP,           0.0_RP,    &
      0.0_RP,        0.0_RP,        0.0_RP,          0.0_RP,          0.0_RP,         0.0_RP, -16.0_RP/11.0_RP,    &
      0.0_RP,        0.0_RP,        0.0_RP,          0.0_RP,          0.0_RP,         0.0_RP,          0.0_RP /),  &
   (/ RK_nstage, RK_nstage /) )

  real(RP), parameter :: RKCoef_b(RK_nstage) = &
      1.0_RP/120.0_RP * (/ 11.0_RP, 0.0_RP, 81.0_RP, 81.0_RP, -32.0_RP, -32.0_RP, 11.0_RP /)

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

    integer  :: i, j, k, iv, n, s
    integer :: stage

    integer, parameter :: ko_dynvar(5) = (/ 0, 1, 0, 0, 0 /)
    integer, parameter :: io_dynvar(5) = (/ 0, 0, 1, 0, 0 /)
    integer, parameter :: jo_dynvar(5) = (/ 0, 0, 0, 1, 0 /)
    integer :: ko_prgvar(VA)
    integer :: io_prgvar(VA)
    integer :: jo_prgvar(VA)

    real(RP) :: dynvar0(KA,IA,JA,5)
    real(RP) :: prgvar0(KA,IA,JA,VA)

    !---------------------------------------------------------------------------

    call PROF_rapstart("DYN_RK7s6o_Prep",3)

#ifdef DEBUG
    !$omp parallel 
    !$omp workshare
    rk_dynvar%rkwork(:,:,:,:,:) = UNDEF
    mflx_hi_RK(:,:,:,:,:) = UNDEF
    tflx_hi_RK(:,:,:,:,:) = UNDEF
    !$omp end workshare
    if (VA > 0) then
      !$omp workshare
      rk_prgvar%rkwork(:,:,:,:,:) = UNDEF
      !$omp end workshare
    end if
    !$omp end parallel
#endif

#ifdef QUICKDEBUG
    mflx_hi(   1:KS-1,:,:,:) = UNDEF
    mflx_hi(KE+1:KA  ,:,:,:) = UNDEF
#endif

!OCL XFILL
    dynvar0(:,:,:,I_RK_DENS) = DENS    
!OCL XFILL
    dynvar0(:,:,:,I_RK_MOMZ) = MOMZ
!OCL XFILL
    dynvar0(:,:,:,I_RK_MOMX) = MOMX
!OCL XFILL
    dynvar0(:,:,:,I_RK_MOMY) = MOMY
!OCL XFILL
    dynvar0(:,:,:,I_RK_RHOT) = RHOT
!OCL XFILL
    rk_dynvar%buf(:,:,:,:) = dynvar0(:,:,:,:)

    if (VA > 0) then
      !OCL XFILL
      prgvar0(:,:,:,:) = PROG
      !OCL XFILL
      rk_prgvar%buf(:,:,:,:) = prgvar0(:,:,:,:)
    end if

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
         call ATMOS_DYN_Copy_boundary( rk_dynvar%buf(:,:,:,I_RK_DENS),           & ! [INOUT] 
                                       rk_dynvar%buf(:,:,:,I_RK_MOMZ),           & ! [INOUT]
                                       rk_dynvar%buf(:,:,:,I_RK_MOMX),           & ! [INOUT]
                                       rk_dynvar%buf(:,:,:,I_RK_MOMY),           & ! [INOUT]
                                       rk_dynvar%buf(:,:,:,I_RK_RHOT),           & ! [INOUT]
                                       rk_prgvar%buf(:,:,:,:        ),           & ! [INOUT]
                                       dynvar0(:,:,:,I_RK_DENS),                 & ! [IN] 
                                       dynvar0(:,:,:,I_RK_MOMZ),                 & ! [IN] 
                                       dynvar0(:,:,:,I_RK_MOMX),                 & ! [IN] 
                                       dynvar0(:,:,:,I_RK_MOMY),                 & ! [IN] 
                                       dynvar0(:,:,:,I_RK_RHOT),                 & ! [IN] 
                                       prgvar0(:,:,:,:),                         & ! [IN]
                                       BND_W, BND_E, BND_S, BND_N                ) ! [IN]
         call PROF_rapend  ("DYN_RK7s6o_BND",3)

         call RKUtil_comm( rk_dynvar )
         call RKUtil_comm( rk_prgvar )
         call RKUtil_comm_wait( rk_dynvar )
         call RKUtil_comm_wait( rk_prgvar )
      end if
      !--

      call PROF_rapstart("DYN_RK7s6o",3)

      call ATMOS_DYN_tstep( &
         rk_dynvar%rkwork(:,:,:,I_RK_DENS,stage),               & ! [OUT] 
         rk_dynvar%rkwork(:,:,:,I_RK_MOMZ,stage),               & ! [OUT]
         rk_dynvar%rkwork(:,:,:,I_RK_MOMX,stage),               & ! [OUT]
         rk_dynvar%rkwork(:,:,:,I_RK_MOMY,stage),               & ! [OUT]
         rk_dynvar%rkwork(:,:,:,I_RK_RHOT,stage),               & ! [OUT]
         rk_prgvar%rkwork(:,:,:,:        ,stage),               & ! [OUT]
         mflx_hi_RK(:,:,:,:,stage), tflx_hi_RK(:,:,:,:,stage),  & ! [INOUT,OUT]
         dynvar0(:,:,:,I_RK_DENS),                              & ! [IN]
         dynvar0(:,:,:,I_RK_MOMZ),                              & ! [IN]
         dynvar0(:,:,:,I_RK_MOMX),                              & ! [IN]
         dynvar0(:,:,:,I_RK_MOMY),                              & ! [IN]
         dynvar0(:,:,:,I_RK_RHOT),                              & ! [IN]
         rk_dynvar%buf(:,:,:,I_RK_DENS),                        & ! [IN] 
         rk_dynvar%buf(:,:,:,I_RK_MOMZ),                        & ! [IN]
         rk_dynvar%buf(:,:,:,I_RK_MOMX),                        & ! [IN]
         rk_dynvar%buf(:,:,:,I_RK_MOMY),                        & ! [IN]
         rk_dynvar%buf(:,:,:,I_RK_RHOT),                        & ! [IN]
         DENS_t,   MOMZ_t,   MOMX_t,   MOMY_t,   RHOT_t,        & ! [IN]
         prgvar0, rk_prgvar%buf(:,:,:,:),                       & ! [IN]
         DPRES0, CVtot, CORIOLI,                                & ! [IN]
         num_diff, wdamp_coef, divdmp_coef, DDIV,               & ! [IN]
         FLAG_FCT_MOMENTUM, FLAG_FCT_T,                         & ! [IN]
         FLAG_FCT_ALONG_STREAM,                                 & ! [IN]
         CDZ, FDZ, FDX, FDY,                                    & ! [IN]
         RCDZ, RCDX, RCDY, RFDZ, RFDX, RFDY,                    & ! [IN]
         PHI, GSQRT, J13G, J23G, J33G, MAPF,                    & ! [IN]
         REF_pres, REF_dens,                                    & ! [IN]
         BND_W, BND_E, BND_S, BND_N,                            & ! [IN]
         1.0_RP, .false.                                        ) ! [IN]
   
      if ( stage < RK_nstage) then
         call calc_var_nextstage( rk_dynvar%buf,                                      & ! (out)
            stage, dynvar0, rk_dynvar%rkwork, io_dynvar, jo_dynvar, ko_dynvar,  5, dt ) ! (in)
         
         call calc_var_nextstage( rk_prgvar%buf,                                      & ! (out)
            stage, prgvar0, rk_prgvar%rkwork, io_prgvar, jo_prgvar, ko_prgvar, VA, dt ) ! (in)
      end if

      call PROF_rapend  ("DYN_RK7s6o",3)
    end do

    call PROF_rapstart("DYN_RK7s6o",3)
    call calc_var_nextstep( DENS, dynvar0(:,:,:,I_RK_DENS), rk_dynvar%rkwork(:,:,:,I_RK_DENS,:), io_dynvar(I_RK_DENS), jo_dynvar(I_RK_DENS), ko_dynvar(I_RK_DENS), 1, dt )
    call calc_var_nextstep( MOMZ, dynvar0(:,:,:,I_RK_MOMZ), rk_dynvar%rkwork(:,:,:,I_RK_MOMZ,:), io_dynvar(I_RK_MOMZ), jo_dynvar(I_RK_MOMZ), ko_dynvar(I_RK_MOMZ), 1, dt )
    call calc_var_nextstep( MOMX, dynvar0(:,:,:,I_RK_MOMX), rk_dynvar%rkwork(:,:,:,I_RK_MOMX,:), io_dynvar(I_RK_MOMX), jo_dynvar(I_RK_MOMX), ko_dynvar(I_RK_MOMX), 1, dt )
    call calc_var_nextstep( MOMY, dynvar0(:,:,:,I_RK_MOMY), rk_dynvar%rkwork(:,:,:,I_RK_MOMY,:), io_dynvar(I_RK_MOMY), jo_dynvar(I_RK_MOMY), ko_dynvar(I_RK_MOMY), 1, dt )
    call calc_var_nextstep( RHOT, dynvar0(:,:,:,I_RK_RHOT), rk_dynvar%rkwork(:,:,:,I_RK_RHOT,:), io_dynvar(I_RK_RHOT), jo_dynvar(I_RK_RHOT), ko_dynvar(I_RK_RHOT), 1, dt )
    call calc_var_nextstep( PROG, prgvar0(:,:,:,:), rk_prgvar%rkwork(:,:,:,:,:), 0, 0, 0, VA, dt )
    call calc_flux_nextstep( mflx_hi, mflx_hi_RK, 0, 0, 0, 3 ) 
    call calc_flux_nextstep( tflx_hi, tflx_hi_RK, 0, 0, 0, 3 )
    call PROF_rapend("DYN_RK7s6o",3)

    return
  end subroutine ATMOS_DYN_tinteg_short_rk7s6o

  !--- private ---------------------------------------------------------------------------

  subroutine calc_var_nextstage( rk_next, nowstage, rkwork0, rkwork, io, jo, ko, va_, dt )

   implicit none
   integer, intent(in) :: va_
   real(RP), intent(out) :: rk_next(KA,IA,JA,va_)
   integer, intent(in) :: nowstage
   real(RP), intent(in) :: rkwork0(KA,IA,JA,va_)
   real(RP), intent(inout) :: rkwork(KA,IA,JA,va_,RK_nregister)
   integer, intent(in) :: io(va_), jo(va_), ko(va_)
   real(RP), intent(in) :: dt

   integer :: i, j, k, iv
   real(RP) :: var0

   real(RP) :: a_(RK_nstage)
   !--------------------------------------

   a_(:) = dt * RKCoef_a(nowstage+1,:)

   select case(nowstage)
   case(1,2)
      !$omp parallel do private(iv,k,j,i,var0) collapse(3)
      do iv=1, va_
      do k=KS-ko(iv), KE
      do j=JS-jo(iv), JE
      do i=IS-io(iv), IE
         var0 = rkwork0(k,i,j,iv)
         rkwork(k,i,j,iv,nowstage) = rkwork(k,i,j,iv,nowstage) - var0

         rk_next(k,i,j,iv) = var0 +                  &
            a_(nowstage) * rkwork(k,i,j,iv,nowstage)
      end do
      end do
      end do
      end do 
   case(3)
      !$omp parallel do private(iv,k,j,i,var0) collapse(3)
      do iv=1, va_
      do k=KS-ko(iv), KE
      do j=JS-jo(iv), JE
      do i=IS-io(iv), IE
         var0 = rkwork0(k,i,j,iv)
         rkwork(k,i,j,iv,3) =  rkwork(k,i,j,iv,3) - var0
   
         rk_next(k,i,j,iv) = var0 + &
            + a_(1) * rkwork(k,i,j,iv,1) + a_(2) * rkwork(k,i,j,iv,2) &
            + a_(3) * rkwork(k,i,j,iv,3) 
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
         var0 = rkwork0(k,i,j,iv)
         rkwork(k,i,j,iv,4) = rkwork(k,i,j,iv,4) - var0

         rk_next(k,i,j,iv) = var0 &
            + a_(1) * rkwork(k,i,j,iv,1) + a_(2) * rkwork(k,i,j,iv,2) &
            + a_(3) * rkwork(k,i,j,iv,3) + a_(4) * rkwork(k,i,j,iv,4) 
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
         var0 = rkwork0(k,i,j,iv)
         rkwork(k,i,j,iv,5) =  rkwork(k,i,j,iv,5) - var0
   
         rk_next(k,i,j,iv) = var0 &
            + a_(2) * rkwork(k,i,j,iv,2) + a_(3) * rkwork(k,i,j,iv,3) &
            + a_(4) * rkwork(k,i,j,iv,4) + a_(5) * rkwork(k,i,j,iv,5) 
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
         var0 = rkwork0(k,i,j,iv)
         rkwork(k,i,j,iv,6) = rkwork(k,i,j,iv,6) - var0
   
         rk_next(k,i,j,iv) = var0 + &
            + a_(1) * rkwork(k,i,j,iv,1) + a_(2) * rkwork(k,i,j,iv,2) &
            + a_(3) * rkwork(k,i,j,iv,3) + a_(4) * rkwork(k,i,j,iv,4) &
            + a_(6) * rkwork(k,i,j,iv,6) 
      end do
      end do
      end do
      end do
   end select

   return
  end subroutine calc_var_nextstage

 subroutine calc_Var_nextstep( var, rkwork0, rkwork, io, jo, ko, va_, dt ) 
   implicit none

   integer, intent(in) :: va_
   real(RP), intent(inout) :: var(KA,IA,JA, va_)
   real(RP), intent(in) :: rkwork0(KA,IA,JA,va_)
   real(RP), intent(in) :: rkwork (KA,IA,JA,va_,RK_nregister)
   integer, intent(in) :: io, jo, ko
   real(RP), intent(in) :: dt

   integer :: i, j, k, iv
   real(RP) :: var0

   real(RP) :: b_(RK_nstage)
   real(RP) :: a77
   !--------------------------------------

   b_(:) = dt * RKCoef_b(:)

   !$omp parallel do private(iv,k,j,i,var0) collapse(3)
   do iv=1, va_
   do k=KS-ko, KE
   do j=JS-jo, JE
   do i=IS-io, IE
     var0 = rkwork0(k,i,j,iv)
     var(k,i,j,iv) = var0  &
              + b_(1) *  rkwork(k,i,j,iv,1) + b_(3) *  rkwork(k,i,j,iv,3)       &
              + b_(4) *  rkwork(k,i,j,iv,4) + b_(5) *  rkwork(k,i,j,iv,5)       &
              + b_(6) *  rkwork(k,i,j,iv,6) + b_(7) * (rkwork(k,i,j,iv,7) - var0) 
   end do
   end do
   end do
   end do

   return
 end subroutine calc_Var_nextstep

 subroutine calc_flux_nextstep( flux, rkwork, io, jo, ko, va_ ) 
   implicit none

   integer, intent(in) :: va_
   real(RP), intent(inout) :: flux(KA,IA,JA, va_)
   real(RP), intent(in) :: rkwork (KA,IA,JA,va_,RK_nregister)
   integer, intent(in) :: io, jo, ko

   integer :: i, j, k, iv
   !--------------------------------------

   !$omp parallel do private(iv,k,j,i) collapse(3)
   do iv=1, va_
   do k=KS-ko, KE
   do j=JS-jo, JE
   do i=IS-io, IE
     flux(k,i,j,iv) = &
        RKCoef_b(1) * rkwork(k,i,j,iv,1) + RKCoef_b(3) * rkwork(k,i,j,iv,3)  &
      + RKCoef_b(4) * rkwork(k,i,j,iv,4) + RKCoef_b(5) * rkwork(k,i,j,iv,5)  &
      + RKCoef_b(6) * rkwork(k,i,j,iv,6) + RKCoef_b(7) * rkwork(k,i,j,iv,7)
   end do
   end do
   end do
   end do

   return
 end subroutine calc_flux_nextstep

end module scale_atmos_dyn_tinteg_short_rk7s6o
