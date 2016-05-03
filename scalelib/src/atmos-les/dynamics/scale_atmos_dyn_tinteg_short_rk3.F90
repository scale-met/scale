!-------------------------------------------------------------------------------
!> module Atmosphere / Dyn Tinteg
!!
!! @par Description
!!          Temporal integration in Dynamical core for Atmospheric process
!!          three step Runge-Kutta scheme
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2016-04-18 (S.Nishizawa) [mod] split from scale_atmos_dyn.F90
!!
!<
!-------------------------------------------------------------------------------
#include "inc_openmp.h"
module scale_atmos_dyn_tinteg_short_rk3
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
  public :: ATMOS_DYN_Tinteg_short_rk3_setup
  public :: ATMOS_DYN_Tinteg_short_rk3

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
  real(RP), private, allocatable :: DENS_RK1(:,:,:) ! prognostic variables (+1/3 step)
  real(RP), private, allocatable :: MOMZ_RK1(:,:,:)
  real(RP), private, allocatable :: MOMX_RK1(:,:,:)
  real(RP), private, allocatable :: MOMY_RK1(:,:,:)
  real(RP), private, allocatable :: RHOT_RK1(:,:,:)
  real(RP), private, allocatable :: PROG_RK1(:,:,:,:)
  real(RP), private, allocatable :: DENS_RK2(:,:,:) ! prognostic variables (+2/3 step)
  real(RP), private, allocatable :: MOMZ_RK2(:,:,:)
  real(RP), private, allocatable :: MOMX_RK2(:,:,:)
  real(RP), private, allocatable :: MOMY_RK2(:,:,:)
  real(RP), private, allocatable :: RHOT_RK2(:,:,:)
  real(RP), private, allocatable :: PROG_RK2(:,:,:,:)

  ! for communication
  integer :: I_COMM_DENS = 1
  integer :: I_COMM_MOMZ = 2
  integer :: I_COMM_MOMX = 3
  integer :: I_COMM_MOMY = 4
  integer :: I_COMM_RHOT = 5
  integer, allocatable :: I_COMM_PROG(:)

  integer :: I_COMM_DENS_t = 1
  integer :: I_COMM_MOMZ_t = 2
  integer :: I_COMM_MOMX_t = 3
  integer :: I_COMM_MOMY_t = 4
  integer :: I_COMM_RHOT_t = 5

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

  integer :: I_COMM_mflx_z = 1
  integer :: I_COMM_mflx_x = 2
  integer :: I_COMM_mflx_y = 3

  logical :: FLAG_WS2002
  real(RP) :: fact_dt1
  real(RP) :: fact_dt2

  !-----------------------------------------------------------------------------
contains

  !-----------------------------------------------------------------------------
  !> Setup
  subroutine ATMOS_DYN_Tinteg_short_rk3_setup( &
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

    select case ( tinteg_type )
    case ( 'RK3' )
       if( IO_L ) write(IO_FID_LOG,*) "*** RK3: Heun's method is used"
       ! Heun's method
       ! k1 = f(\phi_n); r1 = \phi_n + k1 * dt / 3
       ! k2 = f(r1);     r2 = \phi_n + k2 * dt * 2 / 3
       ! k3 = f(r2);     r3 = \phi_n + k3 * dt
       ! \phi_{n+1} = \phi_n + ( k1 + 3 * k3 ) dt / 4
       !            = \phi_n + ( (r1-\phi_n) * 3 + (r3-\phi_n) * 3 ) / 4
       !            = ( r1 * 3 + r3 * 3 - \phi_n * 2 ) / 4
       FLAG_WS2002 = .false.
       fact_dt1 = 1.0_RP / 3.0_RP
       fact_dt2 = 2.0_RP / 3.0_RP
    case ( 'RK3WS2002' )
       if( IO_L ) write(IO_FID_LOG,*) "*** RK3: Wichere and Skamarock (2002) is used"
       ! Wicher and Skamarock (2002) RK3 scheme
       ! k1 = f(\phi_n); r1 = \phi_n + k1 * dt / 3
       ! k2 = f(r1);     r2 = \phi_n + k2 * dt / 2
       ! k3 = f(r2);     r3 = \phi_n + k3 * dt
       ! \phi_{n+1} = r3
       FLAG_WS2002 = .true.
       fact_dt1 = 1.0_RP / 3.0_RP
       fact_dt2 = 1.0_RP / 2.0_RP
    case default
       write(*,*) 'xxx TINTEG_TYPE is not RK3. Check!'
       call PRC_MPIstop
    end select

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

    allocate( PROG_RK1(KA,IA,JA,max(VA,1)) )
    allocate( PROG_RK2(KA,IA,JA,max(VA,1)) )
    allocate( I_COMM_PROG_RK1(max(VA,1)) )
    allocate( I_COMM_PROG_RK2(max(VA,1)) )

    call COMM_vars8_init( DENS_RK1, I_COMM_DENS_RK1 )
    call COMM_vars8_init( MOMZ_RK1, I_COMM_MOMZ_RK1 )
    call COMM_vars8_init( MOMX_RK1, I_COMM_MOMX_RK1 )
    call COMM_vars8_init( MOMY_RK1, I_COMM_MOMY_RK1 )
    call COMM_vars8_init( RHOT_RK1, I_COMM_RHOT_RK1 )
    do iv = 1, VA
       I_COMM_PROG_RK1(iv) = 5 + iv
       call COMM_vars8_init( PROG_RK1(:,:,:,iv), I_COMM_PROG_RK1(iv) )
    end do

    call COMM_vars8_init( DENS_RK2, I_COMM_DENS_RK2 )
    call COMM_vars8_init( MOMZ_RK2, I_COMM_MOMZ_RK2 )
    call COMM_vars8_init( MOMX_RK2, I_COMM_MOMX_RK2 )
    call COMM_vars8_init( MOMY_RK2, I_COMM_MOMY_RK2 )
    call COMM_vars8_init( RHOT_RK2, I_COMM_RHOT_RK2 )
    do iv = 1, VA
       I_COMM_PROG_RK2(iv) = 5 + iv
       call COMM_vars8_init( PROG_RK2(:,:,:,iv), I_COMM_PROG_RK2(iv) )
    end do

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

    return
  end subroutine ATMOS_DYN_Tinteg_short_rk3_setup

  !-----------------------------------------------------------------------------
  !> RK3
  subroutine ATMOS_DYN_tinteg_short_rk3( &
       DENS, MOMZ, MOMX, MOMY, RHOT, PROG,       & ! (inout)
       mflx_hi,  tflx_hi,                        & ! (out)
       DENS_t, MOMZ_t, MOMX_t, MOMY_t, RHOT_t,   & ! (in)
       Rtot, CVtot, CORIOLI,                     & ! (in)
       num_diff, divdmp_coef, DDIV,              & ! (in)
       FLAG_FCT_MOMENTUM, FLAG_FCT_T,            & ! (in)
       FLAG_FCT_ALONG_STREAM,                    & ! (in)
       CDZ, FDZ, FDX, FDY,                       & ! (in)
       RCDZ, RCDX, RCDY, RFDZ, RFDX, RFDY,       & ! (in)
       PHI, GSQRT, J13G, J23G, J33G, MAPF,       & ! (in)
       REF_pres, REF_dens,                       & ! (in)
       BND_W, BND_E, BND_S, BND_N,               & ! (in)
       dt                                        ) ! (in)
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

    integer  :: IIS, IIE
    integer  :: JJS, JJE
    integer  :: i, j, k
    integer  :: iv

    !---------------------------------------------------------------------------

    call PROF_rapstart("DYN_RK3_Prep", 3)

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

    call PROF_rapend  ("DYN_RK3_Prep", 3)


    !------------------------------------------------------------------------
    ! Start RK
    !------------------------------------------------------------------------

    call PROF_rapstart("DYN_RK3", 3)

    !##### RK1 : PROG0,PROG->PROG_RK1 #####

    dtrk = dt * fact_dt1

    call ATMOS_DYN_tstep( DENS_RK1, MOMZ_RK1, MOMX_RK1, MOMY_RK1, RHOT_RK1, & ! (out)
                          PROG_RK1,                                         & ! (out)
                          mflx_hi,  tflx_hi,                                & ! (out)
                          DENS0,    MOMZ0,    MOMX0,    MOMY0,    RHOT0,    & ! (in)
                          DENS,     MOMZ,     MOMX,     MOMY,     RHOT,     & ! (in)
                          DENS_t,   MOMZ_t,   MOMX_t,   MOMY_t,   RHOT_t,   & ! (in)
                          PROG0, PROG,                                      & ! (in)
                          Rtot, CVtot, CORIOLI,                             & ! (in)
                          num_diff, divdmp_coef, DDIV,                      & ! (in)
                          FLAG_FCT_MOMENTUM, FLAG_FCT_T,                    & ! (in)
                          FLAG_FCT_ALONG_STREAM,                            & ! (in)
                          CDZ, FDZ, FDX, FDY,                               & ! (in)
                          RCDZ, RCDX, RCDY, RFDZ, RFDX, RFDY,               & ! (in)
                          PHI, GSQRT, J13G, J23G, J33G, MAPF,               & ! (in)
                          REF_pres, REF_dens,                               & ! (in)
                          BND_W, BND_E, BND_S, BND_N,                       & ! (in)
                          dtrk, dt                                          ) ! (in)

    call PROF_rapend  ("DYN_RK3", 3)

    call PROF_rapstart("DYN_RK3_BND", 3)

    call ATMOS_DYN_Copy_boundary( &
         DENS_RK1, MOMZ_RK1, MOMX_RK1, MOMY_RK1, RHOT_RK1, PROG_RK1, & ! (inout)
         DENS0,    MOMZ0,    MOMX0,    MOMY0,    RHOT0,    PROG0,    & ! (in)
         BND_W, BND_E, BND_S, BND_N ) ! (in)

    call PROF_rapend  ("DYN_RK3_BND", 3)

    call COMM_vars8( DENS_RK1(:,:,:), I_COMM_DENS_RK1 )
    call COMM_vars8( MOMZ_RK1(:,:,:), I_COMM_MOMZ_RK1 )
    call COMM_vars8( MOMX_RK1(:,:,:), I_COMM_MOMX_RK1 )
    call COMM_vars8( MOMY_RK1(:,:,:), I_COMM_MOMY_RK1 )
    call COMM_vars8( RHOT_RK1(:,:,:), I_COMM_RHOT_RK1 )
    do iv = 1, VA
       call COMM_vars8( PROG_RK1(:,:,:,iv), I_COMM_PROG_RK1(iv) )
    end do
    call COMM_wait ( DENS_RK1(:,:,:), I_COMM_DENS_RK1, .false. )
    call COMM_wait ( MOMZ_RK1(:,:,:), I_COMM_MOMZ_RK1, .false. )
    call COMM_wait ( MOMX_RK1(:,:,:), I_COMM_MOMX_RK1, .false. )
    call COMM_wait ( MOMY_RK1(:,:,:), I_COMM_MOMY_RK1, .false. )
    call COMM_wait ( RHOT_RK1(:,:,:), I_COMM_RHOT_RK1, .false. )
    do iv = 1, VA
       call COMM_wait ( PROG_RK1(:,:,:,iv), I_COMM_PROG_RK1(iv), .false. )
    end do

    !##### RK2 : PROG0,PROG_RK1->PROG_RK2 #####

    call PROF_rapstart("DYN_RK3", 3)

    dtrk = dt * fact_dt2

    call ATMOS_DYN_tstep( DENS_RK2, MOMZ_RK2, MOMX_RK2, MOMY_RK2, RHOT_RK2, & ! (out)
                          PROG_RK2,                                         & ! (out)
                          mflx_hi,  tflx_hi,                                & ! (out)
                          DENS0,    MOMZ0,    MOMX0,    MOMY0,    RHOT0,    & ! (in)
                          DENS_RK1, MOMZ_RK1, MOMX_RK1, MOMY_RK1, RHOT_RK1, & ! (in)
                          DENS_t,   MOMZ_t,   MOMX_t,   MOMY_t,   RHOT_t,   & ! (in)
                          PROG0, PROG_RK1,                                  & ! (in)
                          Rtot, CVtot, CORIOLI,                             & ! (in)
                          num_diff, divdmp_coef, DDIV,                      & ! (in)
                          FLAG_FCT_MOMENTUM, FLAG_FCT_T,                    & ! (in)
                          FLAG_FCT_ALONG_STREAM,                            & ! (in)
                          CDZ, FDZ, FDX, FDY,                               & ! (in)
                          RCDZ, RCDX, RCDY, RFDZ, RFDX, RFDY,               & ! (in)
                          PHI, GSQRT, J13G, J23G, J33G, MAPF,               & ! (in)
                          REF_pres, REF_dens,                               & ! (in)
                          BND_W, BND_E, BND_S, BND_N,                       & ! (in)
                          dtrk, dt                                          ) ! (in)

    call PROF_rapend  ("DYN_RK3", 3)

    call PROF_rapstart("DYN_RK3_BND", 3)

    call ATMOS_DYN_Copy_boundary( &
         DENS_RK2, MOMZ_RK2, MOMX_RK2, MOMY_RK2, RHOT_RK2, PROG_RK2, & ! (inout)
         DENS0,    MOMZ0,    MOMX0,    MOMY0,    RHOT0,    PROG0,    & ! (in)
         BND_W, BND_E, BND_S, BND_N ) ! (in)

    call PROF_rapend  ("DYN_RK3_BND", 3)

    call COMM_vars8( DENS_RK2(:,:,:), I_COMM_DENS_RK2 )
    call COMM_vars8( MOMZ_RK2(:,:,:), I_COMM_MOMZ_RK2 )
    call COMM_vars8( MOMX_RK2(:,:,:), I_COMM_MOMX_RK2 )
    call COMM_vars8( MOMY_RK2(:,:,:), I_COMM_MOMY_RK2 )
    call COMM_vars8( RHOT_RK2(:,:,:), I_COMM_RHOT_RK2 )
    do iv = 1, VA
       call COMM_vars8( PROG_RK2(:,:,:,iv), I_COMM_PROG_RK2(iv) )
    end do
    call COMM_wait ( DENS_RK2(:,:,:), I_COMM_DENS_RK2, .false. )
    call COMM_wait ( MOMZ_RK2(:,:,:), I_COMM_MOMZ_RK2, .false. )
    call COMM_wait ( MOMX_RK2(:,:,:), I_COMM_MOMX_RK2, .false. )
    call COMM_wait ( MOMY_RK2(:,:,:), I_COMM_MOMY_RK2, .false. )
    call COMM_wait ( RHOT_RK2(:,:,:), I_COMM_RHOT_RK2, .false. )
    do iv = 1, VA
       call COMM_wait ( PROG_RK2(:,:,:,iv), I_COMM_PROG_RK2(iv), .false. )
    end do

    !##### RK3 : PROG0,PROG_RK2->PROG #####

    call PROF_rapstart("DYN_RK3", 3)

    dtrk = dt

    call ATMOS_DYN_tstep( DENS,     MOMZ,     MOMX,     MOMY,     RHOT,     & ! (out)
                          PROG,                                             & ! (out)
                          mflx_hi,  tflx_hi,                                & ! (out)
                          DENS0,    MOMZ0,    MOMX0,    MOMY0,    RHOT0,    & ! (in)
                          DENS_RK2, MOMZ_RK2, MOMX_RK2, MOMY_RK2, RHOT_RK2, & ! (in)
                          DENS_t,   MOMZ_t,   MOMX_t,   MOMY_t,   RHOT_t,   & ! (in)
                          PROG0, PROG_RK2,                                  & ! (in)
                          Rtot, CVtot, CORIOLI,                             & ! (in)
                          num_diff, divdmp_coef, DDIV,                      & ! (in)
                          FLAG_FCT_MOMENTUM, FLAG_FCT_T,                    & ! (in)
                          FLAG_FCT_ALONG_STREAM,                            & ! (in)
                          CDZ, FDZ, FDX, FDY,                               & ! (in)
                          RCDZ, RCDX, RCDY, RFDZ, RFDX, RFDY,               & ! (in)
                          PHI, GSQRT, J13G, J23G, J33G, MAPF,               & ! (in)
                          REF_pres, REF_dens,                               & ! (in)
                          BND_W, BND_E, BND_S, BND_N,                       & ! (in)
                          dtrk, dt                                          ) ! (in)

    if ( .not. FLAG_WS2002 ) then
       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          DENS(k,i,j) = ( DENS_RK1(k,i,j) * 3.0_RP + DENS(k,i,j) * 3.0_RP - DENS0(k,i,j) * 2.0_RP ) / 4.0_RP
       end do
       end do
       end do
       do j = JS, JE
       do i = IS, IE
       do k = KS, KE-1
          MOMZ(k,i,j) = ( MOMZ_RK1(k,i,j) * 3.0_RP + MOMZ(k,i,j) * 3.0_RP - MOMZ0(k,i,j) * 2.0_RP ) / 4.0_RP
       end do
       end do
       end do
       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          MOMX(k,i,j) = ( MOMX_RK1(k,i,j) * 3.0_RP + MOMX(k,i,j) * 3.0_RP - MOMX0(k,i,j) * 2.0_RP ) / 4.0_RP
       end do
       end do
       end do
       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          MOMY(k,i,j) = ( MOMY_RK1(k,i,j) * 3.0_RP + MOMY(k,i,j) * 3.0_RP - MOMY0(k,i,j) * 2.0_RP ) / 4.0_RP
       end do
       end do
       end do
       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          RHOT(k,i,j) = ( RHOT_RK1(k,i,j) * 3.0_RP + RHOT(k,i,j) * 3.0_RP - RHOT0(k,i,j) * 2.0_RP ) / 4.0_RP
       end do
       end do
       end do
       do iv = 1, VA
       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          PROG(k,i,j,iv) = ( PROG_RK1(k,i,j,iv) * 3.0_RP + PROG(k,i,j,iv) * 3.0_RP - PROG0(k,i,j,iv) * 2.0_RP ) / 4.0_RP
       end do
       end do
       end do
       end do
    end if

    call PROF_rapend  ("DYN_RK3", 3)

    return
  end subroutine ATMOS_DYN_tinteg_short_rk3

end module scale_atmos_dyn_tinteg_short_rk3
