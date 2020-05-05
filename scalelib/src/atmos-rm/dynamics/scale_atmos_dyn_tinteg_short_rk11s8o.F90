!-------------------------------------------------------------------------------
!> module Atmosphere / Dyn Tinteg
!!
!! @par Description
!!          Temporal integration in Dynamical core for Atmospheric process
!!          11 stage Runge-Kutta scheme with 8th order accuracy
!!
!! @author Team SCALE
!!
!! This module provides two type of 11 stage Runge-Kutta scheme with 8th order accuracy proposed by Cooper and Verner (1972). 
!! See this source file for the detail of RK coffecients. 
!<
!-------------------------------------------------------------------------------
#include "scalelib.h"
module scale_atmos_dyn_tinteg_short_rk11s8o
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

  use scale_atmos_dyn_tinteg_rkcommon, only: &
   RKInfo,                                                              &
   RKCommon_setup => ATMOS_DYN_Tinteg_RKCommon_setup,                   &
   RKCommon_rkwork_alloc   => ATMOS_DYN_Tinteg_RKCommon_rkwork_alloc,   &
   RKCommon_rkwork_dealloc => ATMOS_DYN_Tinteg_RKCommon_rkwork_dealloc, &
   RKCommon_comm  => ATMOS_DYN_Tinteg_RKCommon_comm,                    &
   RKCommon_comm_wait => ATMOS_DYN_Tinteg_RKCommon_comm_wait,           &
   RKCommon_nextstage => ATMOS_DYN_Tinteg_RKCommon_nextstage,           &
   RKCommon_updateVar => ATMOS_DYN_Tinteg_RKCommon_updateVar,           &
   RKCommon_updateFlux => ATMOS_DYN_Tinteg_RKCommon_updateFlux


  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: ATMOS_DYN_Tinteg_short_rk11s8o_setup
  public :: ATMOS_DYN_Tinteg_short_rk11s8o

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  !-----------------------------------------------------------------------------

  !- coeffecients for 11 stage RK with 8th order accuracy by Cooper and Verner (1972) --------------
  real(RP), parameter :: sq = 21.0_RP**0.5
  real(RP), parameter :: RC41 = 1.0_RP/7.0_RP
  real(RP), parameter :: RC42 = (-7.0_RP+3.0_RP*sq)/98.0_RP
  real(RP), parameter :: RC43 = (21.0_RP-5.0_RP*sq)/49.0_RP
  real(RP), parameter :: RC51 = (11.0_RP-sq)/84.0_RP
  real(RP), parameter :: RC53 = (18.0_RP-4.0_RP*sq)/63.0_RP
  real(RP), parameter :: RC54 = (21.0_RP+sq)/252.0_RP
  real(RP), parameter :: RC61 = (5.0_RP-sq)/48.0_RP
  real(RP), parameter :: RC63 = (9.0_RP-sq)/36.0_RP
  real(RP), parameter :: RC64 = (-231.0_RP-14.0_RP*sq)/360.0_RP
  real(RP), parameter :: RC65 = (63.0_RP+7.0_RP*sq)/80.0_RP
  real(RP), parameter :: RC71 = (10.0_RP+sq)/42.0_RP
  real(RP), parameter :: RC73 = (-432.0_RP-92.0_RP*sq)/315.0_RP
  real(RP), parameter :: RC74 = (633.0_RP+145.0_RP*sq)/90.0_RP
  real(RP), parameter :: RC75 = (-504.0_RP-115.0_RP*sq)/70.0_RP
  real(RP), parameter :: RC76 = (63.0_RP+13.0_RP*sq)/35.0_RP
  real(RP), parameter :: RC81 = 1.0_RP/14.0_RP
  real(RP), parameter :: RC85 = (14.0_RP+3.0_RP*sq)/126.0_RP
  real(RP), parameter :: RC86 = (13.0_RP+3.0_RP*sq)/63.0_RP
  real(RP), parameter :: RC87 = 1.0_RP/9.0_RP
  real(RP), parameter :: RC91 = 1.0_RP/32.0_RP
  real(RP), parameter :: RC95 = (91.0_RP+21.0_RP*sq)/576.0_RP
  real(RP), parameter :: RC96 = 11.0_RP/72.0_RP
  real(RP), parameter :: RC97 = (-385.0_RP+75.0_RP*sq)/1152.0_RP
  real(RP), parameter :: RC98 = (63.0_RP-13.0_RP*sq)/128.0_RP
  real(RP), parameter :: RC101 = 1.0_RP/14.0_RP
  real(RP), parameter :: RC105 = 1.0_RP/9.0_RP
  real(RP), parameter :: RC106 = (-733.0_RP+147.0_RP*sq)/2205.0_RP
  real(RP), parameter :: RC107 = (515.0_RP-111.0_RP*sq)/504.0_RP
  real(RP), parameter :: RC108 = (-51.0_RP+11.0_RP*sq)/56.0_RP
  real(RP), parameter :: RC109 = (132.0_RP-28.0_RP*sq)/245.0_RP
  real(RP), parameter :: RC115 = (-42.0_RP-7.0_RP*sq)/18.0_RP
  real(RP), parameter :: RC116 = (-18.0_RP-28.0_RP*sq)/45.0_RP
  real(RP), parameter :: RC117 = (-273.0_RP+53.0_RP*sq)/72.0_RP
  real(RP), parameter :: RC118 = (301.0_RP-53.0_RP*sq)/72.0_RP
  real(RP), parameter :: RC119 = (28.0_RP+28.0_RP*sq)/45.0_RP
  real(RP), parameter :: RC1110 = (49.0_RP+7.0_RP*sq)/18.0_RP
  real(RP), parameter :: ZERO = 0.0_RP

  real(RP), parameter, public :: RKCoef_a_11s8o_CooperVerner1972(11,11) = reshape( &
  (/ ZERO, 0.5_RP, 0.25_RP, RC41, RC51, RC61, RC71, RC81, RC91,  RC101,    ZERO,    &
     ZERO,   ZERO, 0.25_RP, RC42, ZERO, ZERO, ZERO, ZERO, ZERO,   ZERO,    ZERO,    &
     ZERO,   ZERO,    ZERO, RC43, RC53, RC63, RC73, ZERO, ZERO,   ZERO,    ZERO,    &
     ZERO,   ZERO,    ZERO, ZERO, RC54, RC64, RC74, ZERO, ZERO,   ZERO,    ZERO,    &
     ZERO,   ZERO,    ZERO, ZERO, ZERO, RC65, RC75, RC85, RC95,  RC105,   RC115,    &
     ZERO,   ZERO,    ZERO, ZERO, ZERO, ZERO, RC76, RC86, RC96,  RC106,   RC116,    &
     ZERO,   ZERO,    ZERO, ZERO, ZERO, ZERO, ZERO, RC87, RC97,  RC107,   RC117,    &
     ZERO,   ZERO,    ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, RC98,  RC108,   RC118,    &
     ZERO,   ZERO,    ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO,  RC109,   RC119,    &
     ZERO,   ZERO,    ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO,   ZERO,  RC1110,    &
     ZERO,   ZERO,    ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO,   ZERO,    ZERO /), &
     shape(RKCoef_a_11s8o_CooperVerner1972) )

  real(RP), parameter, public :: RKCoef_b_11s8o_CooperVerner1972(11) = &
    (/ 0.05_RP, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, 49.0_RP/180.0_RP, 16.0_RP/45.0_RP, &
                                                    49.0_RP/180.0_RP,         0.05_RP /)
 
  !-----------------------------------------------------------------------------  
  !
  !++ Private procedure
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  integer, private, parameter :: RK_nstage    = 11
  integer, private, parameter :: RK_nregister = 11

  type(RKInfo), private :: rk_dynvar
  integer, private, parameter :: I_RK_DENS = 1
  integer, private, parameter :: I_RK_MOMZ = 2
  integer, private, parameter :: I_RK_MOMX = 3
  integer, private, parameter :: I_RK_MOMY = 4
  integer, private, parameter :: I_RK_RHOT = 5

  type(RKInfo), private :: rk_prgvar

  type(RKInfo), private :: rk_mflx_hi
  type(RKInfo), private :: rk_tflx_hi

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine ATMOS_DYN_Tinteg_short_rk11s8o_setup( &
       tinteg_type )
    use scale_prc, only: &
       PRC_abort
    use scale_const, only: &
       UNDEF => CONST_UNDEF
    use scale_comm_cartesC, only: &
       COMM_vars8_init
    implicit none

    character(len=*) :: tinteg_type

    integer :: iv, d
    character(H_SHORT) :: dynvar_name_list(5)
    character(H_SHORT) :: prgvar_name_list(VA)
    character(H_SHORT) :: flux_name_list(3)

    real(RP) :: RKCoef_a(RK_nstage,RK_nstage)
    real(RP) :: RKCoef_b(RK_nstage)

    !---------------------------------------------------------------------------

    select case( trim(tinteg_type) )
    case ('RK11s8o', 'RK11s8oCooperVerner1972')
      RKCoef_a(:,:) = RKCoef_a_11s8o_CooperVerner1972
      RKCoef_b(:) = RKCoef_b_11s8o_CooperVerner1972 
    case default
       LOG_ERROR("ATMOS_DYN_Tinteg_short_rk11s8o_setup",*) 'The specified TINTEG_TYPE is invalid. Check!', tinteg_type
       call PRC_abort
    end select

    dynvar_name_list(1) = 'DENS'
    dynvar_name_list(2) = 'MOMZ'
    dynvar_name_list(3) = 'MOMX'
    dynvar_name_list(4) = 'MOMY'
    dynvar_name_list(5) = 'RHOT'
    call RKCommon_setup( rk_dynvar, RK_nstage, RK_nregister, RKCoef_a, RKCoef_b, dynvar_name_list )

    do iv = 1, VA
      flux_name_list(iv) = 'PROG'
    end do
    call RKCommon_setup(  rk_prgvar, RK_nstage, RK_nregister, RKCoef_a, RKCoef_b,prgvar_name_list, comm_id_offset=5 )

    do d = 1, 3
      flux_name_list(iv) = 'mflx_hi'
    end do
    call RKCommon_setup(  rk_mflx_hi, RK_nstage, RK_nregister, RKCoef_a, RKCoef_b, flux_name_list, is_type_flux=.true. )


    do d = 1, 3
      flux_name_list(iv) = 'tflx_hi'
    end do
    call RKCommon_setup(  rk_tflx_hi, RK_nstage, RK_nregister, RKCoef_a, RKCoef_b, flux_name_list, is_type_flux=.true. )
    
    !----------------------------------------------------

    return
  end subroutine ATMOS_DYN_Tinteg_short_rk11s8o_setup

  !-----------------------------------------------------------------------------
  !> RK7s6o
  subroutine ATMOS_DYN_tinteg_short_rk11s8o( &
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
    real(RP), intent(out) :: tflx_hi(KA,IA,JA,3)

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

    integer  :: i, j, k, iv, n, s
    integer :: stage

    integer, parameter :: ko_dynvar(5) = (/ 0, 1, 0, 0, 0 /)
    integer, parameter :: io_dynvar(5) = (/ 0, 0, 1, 0, 0 /)
    integer, parameter :: jo_dynvar(5) = (/ 0, 0, 0, 1, 0 /)
    integer :: ko_prgvar(VA)
    integer :: io_prgvar(VA)
    integer :: jo_prgvar(VA)

    !---------------------------------------------------------------------------

    call PROF_rapstart("DYN_RK11s8o_Prep",3)

    call RKCommon_rkwork_alloc( rk_dynvar )
    call RKCommon_rkwork_alloc( rk_prgvar )

#ifdef QUICKDEBUG
    rk_mflx_hi%buf(   1:KS-1,:,:,:) = UNDEF
    rk_mflx_hi%buf(KE+1:KA  ,:,:,:) = UNDEF
#endif

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
         rk_mflx_hi%buf(k,IS-1,j,2) = mflx_hi(k,IS-1,j,2)
       end do
       end do
    end if
    if ( BND_E ) then
       do j = JS, JE
       do k = KS, KE
         rk_mflx_hi%buf(k,IE,j,2) = mflx_hi(k,IE,j,2)
       end do
       end do
    end if
    if ( BND_S ) then
       do i = IS, IE
       do k = KS, KE
         rk_mflx_hi%buf(k,i,JS-1,3) = mflx_hi(k,i,JS-1,3)
       end do
       end do
    end if
    if ( BND_N ) then
       do i = IS, IE
       do k = KS, KE
         rk_mflx_hi%buf(k,i,JE,3) = mflx_hi(k,i,JE,3)
       end do
       end do
    end if

    io_prgvar(:) = 0; jo_prgvar(:) = 0; ko_prgvar(:) = 0

    call PROF_rapend  ("DYN_RK11s8o_Prep",3)

    !------------------------------------------------------------------------
    ! Start RK
    !------------------------------------------------------------------------

    do stage = 1, RK_nstage

      if ( stage > 1) then
         call PROF_rapstart("DYN_RK11s8o_BND",3)
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
         call PROF_rapend  ("DYN_RK11s8o_BND",3)

         call RKCommon_comm( rk_dynvar )
         call RKCommon_comm( rk_prgvar )
         call RKCommon_comm_wait( rk_dynvar )
         call RKCommon_comm_wait( rk_prgvar )
      end if
      !--

      call PROF_rapstart("DYN_RK11s8o",3)

      call ATMOS_DYN_tstep( &
         rk_dynvar%work(:,:,:,I_RK_DENS,stage),                & ! [OUT] 
         rk_dynvar%work(:,:,:,I_RK_MOMZ,stage),                & ! [OUT]
         rk_dynvar%work(:,:,:,I_RK_MOMX,stage),                & ! [OUT]
         rk_dynvar%work(:,:,:,I_RK_MOMY,stage),                & ! [OUT]
         rk_dynvar%work(:,:,:,I_RK_RHOT,stage),                & ! [OUT]
         rk_prgvar%work(:,:,:,:        ,stage),                & ! [OUT]
         rk_mflx_hi%buf(:,:,:,:), rk_tflx_hi%buf(:,:,:,:),     & ! [INOUT,OUT]
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
         call RKCommon_nextstage( rk_dynvar, stage, io_dynvar, jo_dynvar, ko_dynvar, dt )
         call RKCommon_nextstage( rk_prgvar, stage, io_prgvar, jo_prgvar, ko_prgvar, dt )
      else
         call RKCommon_updateVar( rk_dynvar, io_dynvar, jo_dynvar, ko_dynvar, I_RK_DENS, I_RK_DENS, dt, DENS )
         call RKCommon_updateVar( rk_dynvar, io_dynvar, jo_dynvar, ko_dynvar, I_RK_MOMZ, I_RK_MOMZ, dt, MOMZ )
         call RKCommon_updateVar( rk_dynvar, io_dynvar, jo_dynvar, ko_dynvar, I_RK_MOMX, I_RK_MOMX, dt, MOMX )
         call RKCommon_updateVar( rk_dynvar, io_dynvar, jo_dynvar, ko_dynvar, I_RK_MOMY, I_RK_MOMY, dt, MOMY )
         call RKCommon_updateVar( rk_dynvar, io_dynvar, jo_dynvar, ko_dynvar, I_RK_RHOT, I_RK_RHOT, dt, RHOT )
         call RKCommon_updateVar( rk_prgvar, io_prgvar, jo_prgvar, ko_prgvar, 1, VA, dt, PROG )     
      end if
      call RKCommon_updateFlux( rk_mflx_hi, stage, 0, 0, 0, 3, mflx_hi ) 
      call RKCommon_updateFlux( rk_tflx_hi, stage, 0, 0, 0, 3, tflx_hi )

      call PROF_rapend  ("DYN_RK11s8o",3)
    end do

    call RKCommon_rkwork_dealloc( rk_dynvar )
    call RKCommon_rkwork_dealloc( rk_prgvar )

    return
  end subroutine ATMOS_DYN_tinteg_short_rk11s8o

end module scale_atmos_dyn_tinteg_short_rk11s8o
