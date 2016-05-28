!-------------------------------------------------------------------------------
!> module Atmosphere / Dyn Tinteg
!!
!! @par Description
!!          Temporal integration in Dynamical core for Atmospheric process
!!          four step Runge-Kutta scheme
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2016-04-30 (S.Nishizawa) [new]
!!
!<
!-------------------------------------------------------------------------------
#include "inc_openmp.h"
module scale_atmos_dyn_tinteg_short_rk4
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
  public :: ATMOS_DYN_Tinteg_short_rk4_setup
  public :: ATMOS_DYN_Tinteg_short_rk4

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
  real(RP), private, allocatable :: DENS_RK1(:,:,:) ! prognostic variables (+1/4 step)
  real(RP), private, allocatable :: MOMZ_RK1(:,:,:)
  real(RP), private, allocatable :: MOMX_RK1(:,:,:)
  real(RP), private, allocatable :: MOMY_RK1(:,:,:)
  real(RP), private, allocatable :: RHOT_RK1(:,:,:)
  real(RP), private, allocatable :: PROG_RK1(:,:,:,:)
  real(RP), private, allocatable :: DENS_RK2(:,:,:) ! prognostic variables (+2/4 step)
  real(RP), private, allocatable :: MOMZ_RK2(:,:,:)
  real(RP), private, allocatable :: MOMX_RK2(:,:,:)
  real(RP), private, allocatable :: MOMY_RK2(:,:,:)
  real(RP), private, allocatable :: RHOT_RK2(:,:,:)
  real(RP), private, allocatable :: PROG_RK2(:,:,:,:)
  real(RP), private, allocatable :: DENS_RK3(:,:,:) ! prognostic variables (+3/4 step)
  real(RP), private, allocatable :: MOMZ_RK3(:,:,:)
  real(RP), private, allocatable :: MOMX_RK3(:,:,:)
  real(RP), private, allocatable :: MOMY_RK3(:,:,:)
  real(RP), private, allocatable :: RHOT_RK3(:,:,:)
  real(RP), private, allocatable :: PROG_RK3(:,:,:,:)

  ! for communication
  integer :: I_COMM_DENS_RK1 = 1
  integer :: I_COMM_MOMZ_RK1 = 2
  integer :: I_COMM_MOMX_RK1 = 3
  integer :: I_COMM_MOMY_RK1 = 4
  integer :: I_COMM_RHOT_RK1 = 5
  integer, allocatable :: I_COMM_PROG_RK1(:)

  integer :: I_COMM_DENS_RK2 = 1
  integer :: I_COMM_MOMZ_RK2 = 2
  integer :: I_COMM_MOMX_RK2 = 3
  integer :: I_COMM_MOMY_RK2 = 4
  integer :: I_COMM_RHOT_RK2 = 5
  integer, allocatable :: I_COMM_PROG_RK2(:)

  integer :: I_COMM_DENS_RK3 = 1
  integer :: I_COMM_MOMZ_RK3 = 2
  integer :: I_COMM_MOMX_RK3 = 3
  integer :: I_COMM_MOMY_RK3 = 4
  integer :: I_COMM_RHOT_RK3 = 5
  integer, allocatable :: I_COMM_PROG_RK3(:)

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine ATMOS_DYN_Tinteg_short_rk4_setup( &
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

    if ( tinteg_type /= 'RK4' ) then
       write(*,*) 'xxx TINTEG_TYPE is not RK4. Check!'
       call PRC_MPIstop
    end if

    allocate( DENS_RK1(KA,IA,JA) )
    allocate( MOMZ_RK1(KA,IA,JA) )
    allocate( MOMX_RK1(KA,IA,JA) )
    allocate( MOMY_RK1(KA,IA,JA) )
    allocate( RHOT_RK1(KA,IA,JA) )

    allocate( DENS_RK2(KA,IA,JA) )
    allocate( MOMZ_RK2(KA,IA,JA) )
    allocate( MOMX_RK2(KA,IA,JA) )
    allocate( MOMY_RK2(KA,IA,JA) )
    allocate( RHOT_RK2(KA,IA,JA) )

    allocate( DENS_RK3(KA,IA,JA) )
    allocate( MOMZ_RK3(KA,IA,JA) )
    allocate( MOMX_RK3(KA,IA,JA) )
    allocate( MOMY_RK3(KA,IA,JA) )
    allocate( RHOT_RK3(KA,IA,JA) )

    allocate( PROG_RK1(KA,IA,JA,max(VA,1)) )
    allocate( PROG_RK2(KA,IA,JA,max(VA,1)) )
    allocate( PROG_RK3(KA,IA,JA,max(VA,1)) )
    allocate( I_COMM_PROG_RK1(max(VA,1)) )
    allocate( I_COMM_PROG_RK2(max(VA,1)) )
    allocate( I_COMM_PROG_RK3(max(VA,1)) )

    call COMM_vars8_init( DENS_RK1, I_COMM_DENS_RK1 )
    call COMM_vars8_init( MOMZ_RK1, I_COMM_MOMZ_RK1 )
    call COMM_vars8_init( MOMX_RK1, I_COMM_MOMX_RK1 )
    call COMM_vars8_init( MOMY_RK1, I_COMM_MOMY_RK1 )
    call COMM_vars8_init( RHOT_RK1, I_COMM_RHOT_RK1 )
    do iv = 1, VA
       I_COMM_PROG_RK1(iv) = 5 + iv
       call COMM_vars8_init( PROG_RK1(:,:,:,iv), I_COMM_PROG_RK1(iv) )
    enddo

    call COMM_vars8_init( DENS_RK2, I_COMM_DENS_RK2 )
    call COMM_vars8_init( MOMZ_RK2, I_COMM_MOMZ_RK2 )
    call COMM_vars8_init( MOMX_RK2, I_COMM_MOMX_RK2 )
    call COMM_vars8_init( MOMY_RK2, I_COMM_MOMY_RK2 )
    call COMM_vars8_init( RHOT_RK2, I_COMM_RHOT_RK2 )
    do iv = 1, VA
       I_COMM_PROG_RK2(iv) = 5 + iv
       call COMM_vars8_init( PROG_RK2(:,:,:,iv), I_COMM_PROG_RK2(iv) )
    enddo

    call COMM_vars8_init( DENS_RK3, I_COMM_DENS_RK3 )
    call COMM_vars8_init( MOMZ_RK3, I_COMM_MOMZ_RK3 )
    call COMM_vars8_init( MOMX_RK3, I_COMM_MOMX_RK3 )
    call COMM_vars8_init( MOMY_RK3, I_COMM_MOMY_RK3 )
    call COMM_vars8_init( RHOT_RK3, I_COMM_RHOT_RK3 )
    do iv = 1, VA
       I_COMM_PROG_RK3(iv) = 5 + iv
       call COMM_vars8_init( PROG_RK3(:,:,:,iv), I_COMM_PROG_RK3(iv) )
    enddo

    DENS_RK1(:,:,:) = UNDEF
    MOMZ_RK1(:,:,:) = UNDEF
    MOMX_RK1(:,:,:) = UNDEF
    MOMY_RK1(:,:,:) = UNDEF
    RHOT_RK1(:,:,:) = UNDEF
    if ( VA > 0 ) PROG_RK1(:,:,:,:) = UNDEF

    DENS_RK2(:,:,:) = UNDEF
    MOMZ_RK2(:,:,:) = UNDEF
    MOMX_RK2(:,:,:) = UNDEF
    MOMY_RK2(:,:,:) = UNDEF
    RHOT_RK2(:,:,:) = UNDEF
    if ( VA > 0 ) PROG_RK2(:,:,:,:) = UNDEF

    DENS_RK3(:,:,:) = UNDEF
    MOMZ_RK3(:,:,:) = UNDEF
    MOMX_RK3(:,:,:) = UNDEF
    MOMY_RK3(:,:,:) = UNDEF
    RHOT_RK3(:,:,:) = UNDEF
    if ( VA > 0 ) PROG_RK3(:,:,:,:) = UNDEF

    return
  end subroutine ATMOS_DYN_Tinteg_short_rk4_setup

  !-----------------------------------------------------------------------------
  !> RK3
  subroutine ATMOS_DYN_tinteg_short_rk4( &
       DENS, MOMZ, MOMX, MOMY, RHOT, PROG,     &
       mflx_hi,  tflx_hi,                      &
       DENS_t, MOMZ_t, MOMX_t, MOMY_t, RHOT_t, &
       Rtot, CVtot, CORIOLI,                   &
       num_diff, divdmp_coef, DDIV,            &
       FLAG_FCT_MOMENTUM, FLAG_FCT_T,          &
       FLAG_FCT_ALONG_STREAM,                  &
       CDZ, FDZ, FDX, FDY,                     &
       RCDZ, RCDX, RCDY, RFDZ, RFDX, RFDY,     &
       PHI, GSQRT, J13G, J23G, J33G, MAPF,     &
       REF_pres, REF_dens,                     &
       BND_W, BND_E, BND_S, BND_N,             &
       dt                                      )
    use scale_comm, only: &
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

    real(RP), intent(in)    :: Rtot(KA,IA,JA)
    real(RP), intent(in)    :: CVtot(KA,IA,JA)
    real(RP), intent(in)    :: CORIOLI(IA,JA)
    real(RP), intent(in)    :: num_diff(KA,IA,JA,5,3)
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

    real(RP) :: DENS0(KA,IA,JA)
    real(RP) :: MOMZ0(KA,IA,JA)
    real(RP) :: MOMX0(KA,IA,JA)
    real(RP) :: MOMY0(KA,IA,JA)
    real(RP) :: RHOT0(KA,IA,JA)
    real(RP) :: PROG0(KA,IA,JA,VA)

    real(RP) :: dtrk

    integer  :: i, j, k, iv
    !---------------------------------------------------------------------------

    call PROF_rapstart("DYN_RK4_Prep",3)

#ifdef DEBUG
    DENS_RK1(:,:,:) = UNDEF
    MOMZ_RK1(:,:,:) = UNDEF
    MOMX_RK1(:,:,:) = UNDEF
    MOMY_RK1(:,:,:) = UNDEF
    RHOT_RK1(:,:,:) = UNDEF
    if ( VA > 0 ) PROG_RK1(:,:,:,:) = UNDEF

    DENS_RK2(:,:,:) = UNDEF
    MOMZ_RK2(:,:,:) = UNDEF
    MOMX_RK2(:,:,:) = UNDEF
    MOMY_RK2(:,:,:) = UNDEF
    RHOT_RK2(:,:,:) = UNDEF
    if ( VA > 0 ) PROG_RK2(:,:,:,:) = UNDEF

    DENS_RK3(:,:,:) = UNDEF
    MOMZ_RK3(:,:,:) = UNDEF
    MOMX_RK3(:,:,:) = UNDEF
    MOMY_RK3(:,:,:) = UNDEF
    RHOT_RK3(:,:,:) = UNDEF
    if ( VA > 0 ) PROG_RK3(:,:,:,:) = UNDEF

    mflx_hi(:,:,:,:) = UNDEF
    tflx_hi(:,:,:,:) = UNDEF
#endif

!OCL XFILL
    DENS0 = DENS
!OCL XFILL
    MOMZ0 = MOMZ
!OCL XFILL
    MOMX0 = MOMX
!OCL XFILL
    MOMY0 = MOMY
!OCL XFILL
    RHOT0 = RHOT
!OCL XFILL
    if ( VA > 0 ) PROG0 = PROG

    call PROF_rapend  ("DYN_RK4_Prep",3)

    !------------------------------------------------------------------------
    ! Start RK
    !------------------------------------------------------------------------

    !##### RK1 : PROG0,PROG->PROG_RK1 #####

    call PROF_rapstart("DYN_RK4",3)

    dtrk = dt / 2.0_RP

    call ATMOS_DYN_tstep( DENS_RK1, MOMZ_RK1, MOMX_RK1, MOMY_RK1, RHOT_RK1, & ! [OUT]
                          PROG_RK1,                                         & ! [OUT]
                          mflx_hi,  tflx_hi,                                & ! [OUT]
                          DENS0,    MOMZ0,    MOMX0,    MOMY0,    RHOT0,    & ! [IN]
                          DENS,     MOMZ,     MOMX,     MOMY,     RHOT,     & ! [IN]
                          DENS_t,   MOMZ_t,   MOMX_t,   MOMY_t,   RHOT_t,   & ! [IN]
                          PROG0, PROG,                                      & ! [IN]
                          Rtot, CVtot, CORIOLI,                             & ! [IN]
                          num_diff, divdmp_coef, DDIV,                      & ! [IN]
                          FLAG_FCT_MOMENTUM, FLAG_FCT_T,                    & ! [IN]
                          FLAG_FCT_ALONG_STREAM,                            & ! [IN]
                          CDZ, FDZ, FDX, FDY,                               & ! [IN]
                          RCDZ, RCDX, RCDY, RFDZ, RFDX, RFDY,               & ! [IN]
                          PHI, GSQRT, J13G, J23G, J33G, MAPF,               & ! [IN]
                          REF_pres, REF_dens,                               & ! [IN]
                          BND_W, BND_E, BND_S, BND_N,                       & ! [IN]
                          dtrk, dt                                          ) ! [IN]

    call PROF_rapend  ("DYN_RK4",3)
    call PROF_rapstart("DYN_RK4_BND",3)

    call ATMOS_DYN_Copy_boundary( DENS_RK1, MOMZ_RK1, MOMX_RK1, MOMY_RK1, RHOT_RK1, & ! [INOUT]
                                  PROG_RK1,                                         & ! [INOUT]
                                  DENS0,    MOMZ0,    MOMX0,    MOMY0,    RHOT0,    & ! [IN]
                                  PROG0,                                            & ! [IN]
                                  BND_W, BND_E, BND_S, BND_N                        ) ! [IN]

    call PROF_rapend  ("DYN_RK4_BND",3)

    call COMM_vars8( DENS_RK1(:,:,:), I_COMM_DENS_RK1 )
    call COMM_vars8( MOMZ_RK1(:,:,:), I_COMM_MOMZ_RK1 )
    call COMM_vars8( MOMX_RK1(:,:,:), I_COMM_MOMX_RK1 )
    call COMM_vars8( MOMY_RK1(:,:,:), I_COMM_MOMY_RK1 )
    call COMM_vars8( RHOT_RK1(:,:,:), I_COMM_RHOT_RK1 )
    do iv = 1, VA
       call COMM_vars8( PROG_RK1(:,:,:,iv), I_COMM_PROG_RK1(iv) )
    enddo

    call COMM_wait ( DENS_RK1(:,:,:), I_COMM_DENS_RK1, .false. )
    call COMM_wait ( MOMZ_RK1(:,:,:), I_COMM_MOMZ_RK1, .false. )
    call COMM_wait ( MOMX_RK1(:,:,:), I_COMM_MOMX_RK1, .false. )
    call COMM_wait ( MOMY_RK1(:,:,:), I_COMM_MOMY_RK1, .false. )
    call COMM_wait ( RHOT_RK1(:,:,:), I_COMM_RHOT_RK1, .false. )
    do iv = 1, VA
       call COMM_wait ( PROG_RK1(:,:,:,iv), I_COMM_PROG_RK1(iv), .false. )
    enddo

    !##### RK2 : PROG0,PROG_RK1->PROG_RK2 #####

    call PROF_rapstart("DYN_RK4",3)

    dtrk = dt / 2.0_RP

    call ATMOS_DYN_tstep( DENS_RK2, MOMZ_RK2, MOMX_RK2, MOMY_RK2, RHOT_RK2, & ! [OUT]
                          PROG_RK2,                                         & ! [OUT]
                          mflx_hi,  tflx_hi,                                & ! [OUT]
                          DENS0,    MOMZ0,    MOMX0,    MOMY0,    RHOT0,    & ! [IN]
                          DENS_RK1, MOMZ_RK1, MOMX_RK1, MOMY_RK1, RHOT_RK1, & ! [IN]
                          DENS_t,   MOMZ_t,   MOMX_t,   MOMY_t,   RHOT_t,   & ! [IN]
                          PROG0, PROG_RK1,                                  & ! [IN]
                          Rtot, CVtot, CORIOLI,                             & ! [IN]
                          num_diff, divdmp_coef, DDIV,                      & ! [IN]
                          FLAG_FCT_MOMENTUM, FLAG_FCT_T,                    & ! [IN]
                          FLAG_FCT_ALONG_STREAM,                            & ! [IN]
                          CDZ, FDZ, FDX, FDY,                               & ! [IN]
                          RCDZ, RCDX, RCDY, RFDZ, RFDX, RFDY,               & ! [IN]
                          PHI, GSQRT, J13G, J23G, J33G, MAPF,               & ! [IN]
                          REF_pres, REF_dens,                               & ! [IN]
                          BND_W, BND_E, BND_S, BND_N,                       & ! [IN]
                          dtrk, dt                                          ) ! [IN]

    call PROF_rapend  ("DYN_RK4",3)
    call PROF_rapstart("DYN_RK4_BND",3)

    call ATMOS_DYN_Copy_boundary( DENS_RK2, MOMZ_RK2, MOMX_RK2, MOMY_RK2, RHOT_RK2, & ! [INOUT]
                                  PROG_RK2,                                         & ! [INOUT]
                                  DENS0,    MOMZ0,    MOMX0,    MOMY0,    RHOT0,    & ! [IN]
                                  PROG0,                                            & ! [IN]
                                  BND_W, BND_E, BND_S, BND_N                        ) ! [IN]

    call PROF_rapend  ("DYN_RK4_BND",3)

    call COMM_vars8( DENS_RK2(:,:,:), I_COMM_DENS_RK2 )
    call COMM_vars8( MOMZ_RK2(:,:,:), I_COMM_MOMZ_RK2 )
    call COMM_vars8( MOMX_RK2(:,:,:), I_COMM_MOMX_RK2 )
    call COMM_vars8( MOMY_RK2(:,:,:), I_COMM_MOMY_RK2 )
    call COMM_vars8( RHOT_RK2(:,:,:), I_COMM_RHOT_RK2 )
    do iv = 1, VA
       call COMM_vars8( PROG_RK2(:,:,:,iv), I_COMM_PROG_RK2(iv) )
    enddo

    call COMM_wait ( DENS_RK2(:,:,:), I_COMM_DENS_RK2, .false. )
    call COMM_wait ( MOMZ_RK2(:,:,:), I_COMM_MOMZ_RK2, .false. )
    call COMM_wait ( MOMX_RK2(:,:,:), I_COMM_MOMX_RK2, .false. )
    call COMM_wait ( MOMY_RK2(:,:,:), I_COMM_MOMY_RK2, .false. )
    call COMM_wait ( RHOT_RK2(:,:,:), I_COMM_RHOT_RK2, .false. )
    do iv = 1, VA
       call COMM_wait ( PROG_RK2(:,:,:,iv), I_COMM_PROG_RK2(iv), .false. )
    enddo

    !##### RK3 : PROG0,PROG_RK2->PROG_RK3 #####

    call PROF_rapstart("DYN_RK4",3)

    dtrk = dt

    call ATMOS_DYN_tstep( DENS_RK3, MOMZ_RK3, MOMX_RK3, MOMY_RK3, RHOT_RK3, & ! [OUT]
                          PROG_RK3,                                         & ! [OUT]
                          mflx_hi,  tflx_hi,                                & ! [OUT]
                          DENS0,    MOMZ0,    MOMX0,    MOMY0,    RHOT0,    & ! [IN]
                          DENS_RK2, MOMZ_RK2, MOMX_RK2, MOMY_RK2, RHOT_RK2, & ! [IN]
                          DENS_t,   MOMZ_t,   MOMX_t,   MOMY_t,   RHOT_t,   & ! [IN]
                          PROG0, PROG_RK2,                                  & ! [IN]
                          Rtot, CVtot, CORIOLI,                             & ! [IN]
                          num_diff, divdmp_coef, DDIV,                      & ! [IN]
                          FLAG_FCT_MOMENTUM, FLAG_FCT_T,                    & ! [IN]
                          FLAG_FCT_ALONG_STREAM,                            & ! [IN]
                          CDZ, FDZ, FDX, FDY,                               & ! [IN]
                          RCDZ, RCDX, RCDY, RFDZ, RFDX, RFDY,               & ! [IN]
                          PHI, GSQRT, J13G, J23G, J33G, MAPF,               & ! [IN]
                          REF_pres, REF_dens,                               & ! [IN]
                          BND_W, BND_E, BND_S, BND_N,                       & ! [IN]
                          dtrk, dt                                          ) ! [IN]

    call PROF_rapend  ("DYN_RK4",3)
    call PROF_rapstart("DYN_RK4_BND",3)

    call ATMOS_DYN_Copy_boundary( DENS_RK3, MOMZ_RK3, MOMX_RK3, MOMY_RK3, RHOT_RK3, & ! [INOUT]
                                  PROG_RK3,                                         & ! [INOUT]
                                  DENS0,    MOMZ0,    MOMX0,    MOMY0,    RHOT0,    & ! [IN]
                                  PROG0,                                            & ! [IN]
                                  BND_W, BND_E, BND_S, BND_N                        ) ! [IN]

    call PROF_rapend  ("DYN_RK4_BND",3)

    call COMM_vars8( DENS_RK3(:,:,:), I_COMM_DENS_RK3 )
    call COMM_vars8( MOMZ_RK3(:,:,:), I_COMM_MOMZ_RK3 )
    call COMM_vars8( MOMX_RK3(:,:,:), I_COMM_MOMX_RK3 )
    call COMM_vars8( MOMY_RK3(:,:,:), I_COMM_MOMY_RK3 )
    call COMM_vars8( RHOT_RK3(:,:,:), I_COMM_RHOT_RK3 )
    do iv = 1, VA
       call COMM_vars8( PROG_RK3(:,:,:,iv), I_COMM_PROG_RK3(iv) )
    enddo

    call COMM_wait ( DENS_RK3(:,:,:), I_COMM_DENS_RK3, .false. )
    call COMM_wait ( MOMZ_RK3(:,:,:), I_COMM_MOMZ_RK3, .false. )
    call COMM_wait ( MOMX_RK3(:,:,:), I_COMM_MOMX_RK3, .false. )
    call COMM_wait ( MOMY_RK3(:,:,:), I_COMM_MOMY_RK3, .false. )
    call COMM_wait ( RHOT_RK3(:,:,:), I_COMM_RHOT_RK3, .false. )
    do iv = 1, VA
       call COMM_wait ( PROG_RK3(:,:,:,iv), I_COMM_PROG_RK3(iv), .false. )
    enddo

    !##### RK4 : PROG0,PROG_RK3->PROG #####

    call PROF_rapstart("DYN_RK4",3)

    dtrk = dt

    call ATMOS_DYN_tstep( DENS,     MOMZ,     MOMX,     MOMY,     RHOT,     & ! [OUT]
                          PROG,                                             & ! [OUT]
                          mflx_hi,  tflx_hi,                                & ! [OUT]
                          DENS0,    MOMZ0,    MOMX0,    MOMY0,    RHOT0,    & ! [IN]
                          DENS_RK3, MOMZ_RK3, MOMX_RK3, MOMY_RK3, RHOT_RK3, & ! [IN]
                          DENS_t,   MOMZ_t,   MOMX_t,   MOMY_t,   RHOT_t,   & ! [IN]
                          PROG0, PROG_RK3,                                  & ! [IN]
                          Rtot, CVtot, CORIOLI,                             & ! [IN]
                          num_diff, divdmp_coef, DDIV,                      & ! [IN]
                          FLAG_FCT_MOMENTUM, FLAG_FCT_T,                    & ! [IN]
                          FLAG_FCT_ALONG_STREAM,                            & ! [IN]
                          CDZ, FDZ, FDX, FDY,                               & ! [IN]
                          RCDZ, RCDX, RCDY, RFDZ, RFDX, RFDY,               & ! [IN]
                          PHI, GSQRT, J13G, J23G, J33G, MAPF,               & ! [IN]
                          REF_pres, REF_dens,                               & ! [IN]
                          BND_W, BND_E, BND_S, BND_N,                       & ! [IN]
                          dtrk, dt                                          ) ! [IN]

    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       DENS(k,i,j) = ( DENS_RK1(k,i,j) * 2.0_RP &
                     + DENS_RK2(k,i,j) * 4.0_RP &
                     + DENS_RK3(k,i,j) * 2.0_RP &
                     + DENS    (k,i,j)          &
                     - DENS0   (k,i,j) * 3.0_RP ) / 6.0_RP
    enddo
    enddo
    enddo

    do j = JS, JE
    do i = IS, IE
    do k = KS, KE-1
       MOMZ(k,i,j) = ( MOMZ_RK1(k,i,j) * 2.0_RP &
                     + MOMZ_RK2(k,i,j) * 4.0_RP &
                     + MOMZ_RK3(k,i,j) * 2.0_RP &
                     + MOMZ    (k,i,j)          &
                     - MOMZ0   (k,i,j) * 3.0_RP ) / 6.0_RP
    enddo
    enddo
    enddo

    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       MOMX(k,i,j) = ( MOMX_RK1(k,i,j) * 2.0_RP &
                     + MOMX_RK2(k,i,j) * 4.0_RP &
                     + MOMX_RK3(k,i,j) * 2.0_RP &
                     + MOMX    (k,i,j)          &
                     - MOMX0   (k,i,j) * 3.0_RP ) / 6.0_RP
    enddo
    enddo
    enddo

    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       MOMY(k,i,j) = ( MOMY_RK1(k,i,j) * 2.0_RP &
                     + MOMY_RK2(k,i,j) * 4.0_RP &
                     + MOMY_RK3(k,i,j) * 2.0_RP &
                     + MOMY    (k,i,j)          &
                     - MOMY0   (k,i,j) * 3.0_RP ) / 6.0_RP
    enddo
    enddo
    enddo

    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       RHOT(k,i,j) = ( RHOT_RK1(k,i,j) * 2.0_RP &
                     + RHOT_RK2(k,i,j) * 4.0_RP &
                     + RHOT_RK3(k,i,j) * 2.0_RP &
                     + RHOT    (k,i,j)          &
                     - RHOT0   (k,i,j) * 3.0_RP ) / 6.0_RP
    enddo
    enddo
    enddo

    do iv = 1, VA
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       PROG(k,i,j,iv) = ( PROG_RK1(k,i,j,iv) * 2.0_RP &
                        + PROG_RK2(k,i,j,iv) * 4.0_RP &
                        + PROG_RK3(k,i,j,iv) * 2.0_RP &
                        + PROG    (k,i,j,iv)          &
                        - PROG0   (k,i,j,iv) * 3.0_RP ) / 6.0_RP
    enddo
    enddo
    enddo
    enddo

    call PROF_rapend  ("DYN_RK4",3)

    return
  end subroutine ATMOS_DYN_tinteg_short_rk4

end module scale_atmos_dyn_tinteg_short_rk4
