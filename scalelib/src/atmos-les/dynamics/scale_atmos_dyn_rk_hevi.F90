!-------------------------------------------------------------------------------
!> module Atmosphere / Dynamics RK
!!
!! @par Description
!!          Runge-Kutta for Atmospheric dynamical process
!!          HEVI, no terrain + tracer FCT limiter
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2013-06-18 (S.Nishizawa) [new] newly impremented
!! @li      2014-04-04 (S.Nishizawa) [mod] support terrain-following coordinate
!!
!<
!-------------------------------------------------------------------------------
#include "inc_openmp.h"

#ifdef PROFILE_FAPP
#define PROFILE_START(name) call fapp_start(name, 1, 1)
#define PROFILE_STOP(name)  call fapp_stop (name, 1, 1)
#elif defined(PROFILE_FINEPA)
#define PROFILE_START(name) call start_collection(name)
#define PROFILE_STOP(name)  call stop_collection (name)
#else
#define PROFILE_START(name)
#define PROFILE_STOP(name)
#endif

module scale_atmos_dyn_rk_hevi
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
#endif
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: ATMOS_DYN_rk_hevi_setup
  public :: ATMOS_DYN_rk_hevi

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
  integer, private, parameter :: NB = 1
  real(RP), private :: kappa

  integer :: IFS_OFF
  integer :: JFS_OFF
#ifdef HIST_TEND
  real(RP), allocatable :: advch_t(:,:,:,:)
  real(RP), allocatable :: advcv_t(:,:,:,:)
  real(RP), allocatable :: ddiv_t(:,:,:,:)
  real(RP), allocatable :: pg_t(:,:,:,:)
  real(RP), allocatable :: cf_t(:,:,:,:)
#endif
  !-----------------------------------------------------------------------------


contains

  subroutine ATMOS_DYN_rk_hevi_setup( &
       ATMOS_DYN_TYPE, &
       BND_W, &
       BND_E, &
       BND_S, &
       BND_N  )
    use scale_process, only: &
       PRC_MPIstop
#ifdef DRY
    use scale_const, only: &
       CVdry  => CONST_CVdry,  &
       CPdry  => CONST_CPdry
#endif
    implicit none

    character(len=*), intent(in) :: ATMOS_DYN_TYPE
    logical, intent(in) :: BND_W
    logical, intent(in) :: BND_E
    logical, intent(in) :: BND_S
    logical, intent(in) :: BND_N
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*) '*** HEVI'

    if ( ATMOS_DYN_TYPE .ne. 'HEVI' ) then
       write(*,*) 'xxx ATMOS_DYN_TYPE is not HEVI. Check!'
       call PRC_MPIstop
    end if

#ifdef HEVI_BICGSTAB
    if ( IO_L ) write(IO_FID_LOG,*) '*** USING Bi-CGSTAB'
#elif defined(HEVI_LAPACK)
    if ( IO_L ) write(IO_FID_LOG,*) '*** USING LAPACK'
#else
    if ( IO_L ) write(IO_FID_LOG,*) '*** USING DIRECT'
#endif

#ifdef DRY
    kappa = CPdry / CVdry
#endif

    IFS_OFF = 1
    JFS_OFF = 1
    if ( BND_W ) IFS_OFF = 0
    if ( BND_S ) JFS_OFF = 0


#ifdef HIST_TEND
    allocate( advch_t(KA,IA,JA,5) )
    allocate( advcv_t(KA,IA,JA,5) )
    allocate( ddiv_t(KA,IA,JA,3) )
    allocate( pg_t(KA,IA,JA,3) )
    allocate( cf_t(KA,IA,JA,2) )

    advch_t = 0.0_RP
    advcv_t = 0.0_RP
    ddiv_t = 0.0_RP
    pg_t = 0.0_RP
    cf_t = 0.0_RP
#endif


    return
  end subroutine ATMOS_DYN_rk_hevi_setup


  subroutine ATMOS_DYN_rk_hevi( &
    DENS_RK, MOMZ_RK, MOMX_RK, MOMY_RK, RHOT_RK, &
    mflx_hi, tflx_hi,                            &
    DENS0,   MOMZ0,   MOMX0,   MOMY0,   RHOT0,   &
    DENS,    MOMZ,    MOMX,    MOMY,    RHOT,    &
    DENS_t,  MOMZ_t,  MOMX_t,  MOMY_t,  RHOT_t,  &
    Rtot, CVtot, CORIOLI,                        &
    num_diff, divdmp_coef,                       &
    FLAG_FCT_RHO, FLAG_FCT_MOMENTUM, FLAG_FCT_T, &
    FLAG_FCT_ALONG_STREAM,                       &
    CDZ, FDZ, FDX, FDY,                          &
    RCDZ, RCDX, RCDY, RFDZ, RFDX, RFDY,          &
    PHI, GSQRT, J13G, J23G, J33G, MAPF,          &
    REF_pres, REF_dens,                          &
    BND_W, BND_E, BND_S, BND_N,                  &
    dtrk, dt                                     )
    use scale_grid_index
    use scale_const, only: &
#ifdef DRY
       Rdry   => CONST_Rdry,   &
       CVdry  => CONST_CVdry,  &
       CPdry  => CONST_CPdry,  &
#endif
       UNDEF  => CONST_UNDEF,  &
       IUNDEF => CONST_UNDEF2, &
       GRAV   => CONST_GRAV,   &
       P00    => CONST_PRE00
    use scale_atmos_dyn_common, only: &
       FACT_N, &
       FACT_F, &
       ATMOS_DYN_fct
    use scale_gridtrans, only: &
       I_XYZ, &
       I_XYW, &
       I_UYW, &
       I_XVW, &
       I_UYZ, &
       I_XVZ, &
       I_UVZ, &
       I_XY , &
       I_UY , &
       I_XV , &
       I_UV
#ifdef HIST_TEND
    use scale_history, only: &
       HIST_in
#endif
    implicit none

    real(RP), intent(out) :: DENS_RK(KA,IA,JA)   ! prognostic variables
    real(RP), intent(out) :: MOMZ_RK(KA,IA,JA)   !
    real(RP), intent(out) :: MOMX_RK(KA,IA,JA)   !
    real(RP), intent(out) :: MOMY_RK(KA,IA,JA)   !
    real(RP), intent(out) :: RHOT_RK(KA,IA,JA)   !

    real(RP), intent(inout) :: mflx_hi(KA,IA,JA,3) ! rho * vel(x,y,z)
    real(RP), intent(inout) :: tflx_hi(KA,IA,JA,3) ! rho * theta * vel(x,y,z)

    real(RP), intent(in),target :: DENS0(KA,IA,JA) ! prognostic variables
    real(RP), intent(in),target :: MOMZ0(KA,IA,JA) ! at previous dynamical time step
    real(RP), intent(in),target :: MOMX0(KA,IA,JA) !
    real(RP), intent(in),target :: MOMY0(KA,IA,JA) !
    real(RP), intent(in),target :: RHOT0(KA,IA,JA) !

    real(RP), intent(in) :: DENS(KA,IA,JA)   ! prognostic variables
    real(RP), intent(in) :: MOMZ(KA,IA,JA)   ! at previous RK step
    real(RP), intent(in) :: MOMX(KA,IA,JA)   !
    real(RP), intent(in) :: MOMY(KA,IA,JA)   !
    real(RP), intent(in) :: RHOT(KA,IA,JA)   !

    real(RP), intent(in) :: DENS_t(KA,IA,JA)
    real(RP), intent(in) :: MOMZ_t(KA,IA,JA)
    real(RP), intent(in) :: MOMX_t(KA,IA,JA)
    real(RP), intent(in) :: MOMY_t(KA,IA,JA)
    real(RP), intent(in) :: RHOT_t(KA,IA,JA)

    real(RP), intent(in) :: Rtot(KA,IA,JA) ! R for dry air + vapor
    real(RP), intent(in) :: CVtot(KA,IA,JA) ! CV
    real(RP), intent(in) :: CORIOLI(1,IA,JA)
    real(RP), intent(in) :: num_diff(KA,IA,JA,5,3)
    real(RP), intent(in) :: divdmp_coef

    logical,  intent(in) :: FLAG_FCT_RHO
    logical,  intent(in) :: FLAG_FCT_MOMENTUM
    logical,  intent(in) :: FLAG_FCT_T
    logical,  intent(in) :: FLAG_FCT_ALONG_STREAM

    real(RP), intent(in) :: CDZ(KA)
    real(RP), intent(in) :: FDZ(KA-1)
    real(RP), intent(in) :: FDX(IA-1)
    real(RP), intent(in) :: FDY(JA-1)
    real(RP), intent(in) :: RCDZ(KA)
    real(RP), intent(in) :: RCDX(IA)
    real(RP), intent(in) :: RCDY(JA)
    real(RP), intent(in) :: RFDZ(KA-1)
    real(RP), intent(in) :: RFDX(IA-1)
    real(RP), intent(in) :: RFDY(JA-1)

    real(RP), intent(in)  :: PHI     (KA,IA,JA)   !< geopotential
    real(RP), intent(in)  :: GSQRT   (KA,IA,JA,7) !< vertical metrics {G}^1/2
    real(RP), intent(in)  :: J13G    (KA,IA,JA,7) !< (1,3) element of Jacobian matrix
    real(RP), intent(in)  :: J23G    (KA,IA,JA,7) !< (2,3) element of Jacobian matrix
    real(RP), intent(in)  :: J33G                 !< (3,3) element of Jacobian matrix
    real(RP), intent(in)  :: MAPF    (IA,JA,2,4)  !< map factor
    real(RP), intent(in)  :: REF_pres(KA,IA,JA)   !< reference pressure
    real(RP), intent(in)  :: REF_dens(KA,IA,JA)   !< reference density

    logical,  intent(in)  :: BND_W
    logical,  intent(in)  :: BND_E
    logical,  intent(in)  :: BND_S
    logical,  intent(in)  :: BND_N

    real(RP), intent(in)  :: dtrk
    real(RP), intent(in)  :: dt


    ! diagnostic variables (work space)
    real(RP) :: PRES(KA,IA,JA) ! pressure [Pa]
    real(RP) :: VELZ(KA,IA,JA) ! velocity w [m/s]
    real(RP) :: VELX(KA,IA,JA) ! velocity u [m/s]
    real(RP) :: VELY(KA,IA,JA) ! velocity v [m/s]
    real(RP) :: POTT(KA,IA,JA) ! potential temperature [K]
    real(RP) :: DDIV(KA,IA,JA) ! divergence
    real(RP) :: DPRES(KA,IA,JA) ! pressure deviation from reference pressure

    real(RP) :: qflx_hi(KA,IA,JA,3)
    real(RP) :: qflx_J (KA,IA,JA)

    real(RP) :: advch ! horizontal advection
    real(RP) :: advcv ! vertical advection
    real(RP) :: div  ! divergence damping
    real(RP) :: pg   ! pressure gradient force
    real(RP) :: cf   ! colioris force
#ifdef HIST_TEND
    real(RP) :: advch_t(KA,IA,JA,5)
    real(RP) :: advcv_t(KA,IA,JA,5)
    real(RP) :: ddiv_t(KA,IA,JA,3)
    real(RP) :: pg_t(KA,IA,JA,3)
    real(RP) :: cf_t(KA,IA,JA,3)
    logical  :: lhist
#endif

    ! for implicit solver
    real(RP) :: A(KA,IA,JA)
    real(RP) :: B
    real(RP) :: Sr(KA,IA,JA)
    real(RP) :: Sw(KA,IA,JA)
    real(RP) :: St(KA,IA,JA)
    real(RP) :: PT(KA,IA,JA)
    real(RP) :: CPRES(KA,IA,JA) ! kappa * PRES / RHOT
    real(RP) :: CPtot(KA,IA,JA)
    real(RP) :: C(KMAX-1,IA,JA)

    real(RP) :: F1(KA,IA,JA)
    real(RP) :: F2(KA,IA,JA)
    real(RP) :: F3(KA,IA,JA)

    integer :: IIS, IIE, JJS, JJE
    integer :: k, i, j
    integer :: iss, iee

#ifdef DEBUG
    PRES(:,:,:) = UNDEF
    VELZ(:,:,:) = UNDEF
    VELX(:,:,:) = UNDEF
    VELY(:,:,:) = UNDEF
    POTT(:,:,:) = UNDEF
    DPRES(:,:,:) = UNDEF

    PT(:,:,:) = UNDEF
    CPRES(:,:,:) = UNDEF


    qflx_hi(:,:,:,:) = UNDEF
    qflx_J (:,:,:)   = UNDEF
#endif

#ifdef HIST_TEND
    lhist = dt .eq. dtrk
#endif

#ifdef PROFILE_FIPP
    call fipp_start()
#endif

    if ( .not. divdmp_coef > 0.0_RP ) then
!OCL XFILL
       DDIV(:,:,:) = 0.0_RP
    end if

    do JJS = JS, JE, JBLOCK
    JJE = JJS+JBLOCK-1
    do IIS = IS, IE, IBLOCK
    IIE = IIS+IBLOCK-1

       ! momentum -> velocity
       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j = JJS-1, JJE+2
       do i = IIS-1, IIE+2
       do k = KS, KE-1
#ifdef DEBUG
          call CHECK( __LINE__, MOMZ(k,i,j) )
          call CHECK( __LINE__, DENS(k  ,i,j) )
          call CHECK( __LINE__, DENS(k+1,i,j) )
#endif
          VELZ(k,i,j) = 2.0_RP * MOMZ(k,i,j) / ( DENS(k+1,i,j)+DENS(k,i,j) )
       enddo
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif

       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j = JJS-1, JJE+2
       do i = IIS-2, IIE+1
       do k = KS, KE
#ifdef DEBUG
          call CHECK( __LINE__, MOMX(k,i,j) )
          call CHECK( __LINE__, DENS(k,i  ,j) )
          call CHECK( __LINE__, DENS(k,i+1,j) )
#endif
          VELX(k,i,j) = 2.0_RP * MOMX(k,i,j) &
                      / ( ( DENS(k,i+1,j) + DENS(k,i  ,j) ) * MAPF(i,j,2,I_UY) )
       enddo
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif

       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j = JJS-2, JJE+1
       do i = IIS-1, IIE+2
       do k = KS, KE
#ifdef DEBUG
          call CHECK( __LINE__, MOMY(k,i,j) )
          call CHECK( __LINE__, DENS(k,i,j  ) )
          call CHECK( __LINE__, DENS(k,i,j+1) )
#endif
          VELY(k,i,j) = 2.0_RP * MOMY(k,i,j) &
                      / ( ( DENS(k,i,j+1)+DENS(k,i,j) ) * MAPF(i,j,1,I_XV) )
       enddo
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif

       PROFILE_START("hevi_pres")
       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j = JJS, JJE
       do i = IIS, IIE+1
          do k = KS, KE
#ifdef DEBUG
             call CHECK( __LINE__, RHOT(k,i,j) )
             call CHECK( __LINE__, Rtot(k,i,j) )
             call CHECK( __LINE__, CVtot(k,i,j) )
#endif
#ifdef DRY
             PRES(k,i,j) = P00 * ( RHOT(k,i,j) * Rdry / P00 )**kappa
#else
             CPtot(k,i,j) = CVtot(k,i,j) + Rtot(k,i,j)
             PRES(k,i,j) = P00 * ( RHOT(k,i,j) * Rtot(k,i,j) / P00 )**(CPtot(k,i,j)/CVtot(k,i,j))
#endif
             DPRES(k,i,j) = PRES(k,i,j) - REF_pres(k,i,j)
          enddo
          PRES (KS-1,i,j) = PRES(KS+1,i,j) - DENS(KS,i,j) * ( PHI(KS-1,i,j) - PHI(KS+1,i,j) )
          DPRES(KS-1,i,j) = PRES(KS-1,i,j) - REF_pres(KS-1,i,j)
          PRES (KE+1,i,j) = PRES(KE-1,i,j) - DENS(KE,i,j) * ( PHI(KE+1,i,j) - PHI(KE-1,i,j) )
          DPRES(KE+1,i,j) = PRES(KE+1,i,j) - REF_pres(KE+1,i,j)
       enddo
       enddo
       ! j = JJE+1
       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(1)
       do i = IIS, IIE
          do k = KS, KE
#ifdef DEBUG
             call CHECK( __LINE__, RHOT(k,i,JJE+1) )
             call CHECK( __LINE__, Rtot(k,i,JJE+1) )
             call CHECK( __LINE__, CVtot(k,i,JJE+1) )
#endif
#ifdef DRY
             PRES(k,i,JJE+1) = P00 * ( RHOT(k,i,JJE+1) * Rdry / P00 )**kappa
#else
             CPtot(k,i,JJE+1) = CVtot(k,i,JJE+1) + Rtot(k,i,JJE+1)
             PRES(k,i,JJE+1) = P00 * ( RHOT(k,i,JJE+1) * Rtot(k,i,JJE+1) / P00 )**(CPtot(k,i,JJE+1)/CVtot(k,i,JJE+1))
#endif
             DPRES(k,i,JJE+1) = PRES(k,i,JJE+1) - REF_pres(k,i,JJE+1)
          enddo
          PRES (KS-1,i,JJE+1) = PRES(KS+1,i,JJE+1) - DENS(KS,i,JJE+1) * ( PHI(KS-1,i,JJE+1) - PHI(KS+1,i,JJE+1) )
          DPRES(KS-1,i,JJE+1) = PRES(KS-1,i,JJE+1) - REF_pres(KS-1,i,JJE+1)
          PRES (KE+1,i,JJE+1) = PRES(KE-1,i,JJE+1) - DENS(KE,i,JJE+1) * ( PHI(KE+1,i,JJE+1) - PHI(KE-1,i,JJE+1) )
          DPRES(KE+1,i,JJE+1) = PRES(KE+1,i,JJE+1) - REF_pres(KE+1,i,JJE+1)
       enddo
       PROFILE_STOP("hevi_pres")

       if ( divdmp_coef > 0.0_RP ) then
          ! 3D divergence for damping
          !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
          do j = JJS, JJE+1
          do i = IIS, IIE+1
#ifdef DEBUG
             call CHECK( __LINE__, MOMZ(KS  ,i  ,j  ) )
             call CHECK( __LINE__, MOMX(KS  ,i  ,j  ) )
             call CHECK( __LINE__, MOMX(KS  ,i-1,j  ) )
             call CHECK( __LINE__, MOMY(KS  ,i  ,j  ) )
             call CHECK( __LINE__, MOMY(KS  ,i  ,j-1) )
             call CHECK( __LINE__, MOMZ(KE-1,i  ,j  ) )
             call CHECK( __LINE__, MOMX(KE  ,i  ,j  ) )
             call CHECK( __LINE__, MOMX(KE  ,i-1,j  ) )
             call CHECK( __LINE__, MOMY(KE  ,i  ,j  ) )
             call CHECK( __LINE__, MOMY(KE  ,i  ,j-1) )
#endif
             DDIV(KS,i,j) = ( MOMZ(KS,i,j)                     ) * RCDZ(KS) &
                          + MAPF(i,j,1,I_XY) * MAPF(i,j,2,I_XY) &
                          * ( ( MOMX(KS,i,j)/MAPF(i,j,2,I_UY) - MOMX(KS,i-1,j)/MAPF(i-1,j,2,I_UY) ) * RCDX(i) &
                            + ( MOMY(KS,i,j)/MAPF(i,j,1,I_XV) - MOMY(KS,i,j-1)/MAPF(i,j-1,1,I_XV) ) * RCDY(j) )
             do k = KS+1, KE-1
#ifdef DEBUG
                call CHECK( __LINE__, MOMZ(k  ,i  ,j  ) )
                call CHECK( __LINE__, MOMZ(k-1,i  ,j  ) )
                call CHECK( __LINE__, MOMX(k  ,i  ,j  ) )
                call CHECK( __LINE__, MOMX(k  ,i-1,j  ) )
                call CHECK( __LINE__, MOMY(k  ,i  ,j  ) )
                call CHECK( __LINE__, MOMY(k  ,i  ,j-1) )
#endif
                DDIV(k,i,j) = ( MOMZ(k,i,j) - MOMZ(k-1,i  ,j  ) ) * RCDZ(k) &
                            + MAPF(i,j,1,I_XY) * MAPF(i,j,2,I_XY) &
                            * ( ( MOMX(k,i,j)/MAPF(i,j,2,I_UY) - MOMX(k,i-1,j)/MAPF(i-1,j,2,I_UY) ) * RCDX(i) &
                              + ( MOMY(k,i,j)/MAPF(i,j,1,I_XV) - MOMY(k,i,j-1)/MAPF(i,j-1,1,I_XV) ) * RCDY(j) )
             enddo
             DDIV(KE,i,j) = (             - MOMZ(KE-1,i  ,j  ) ) * RCDZ(KE) &
                          + MAPF(i,j,1,I_XY) * MAPF(i,j,2,I_XY) &
                          * ( ( MOMX(KE,i,j)/MAPF(i,j,2,I_UY) - MOMX(KE,i-1,j)/MAPF(i-1,j,2,I_UY) )* RCDX(i) &
                            + ( MOMY(KE,i,j)/MAPF(i,j,1,I_XV) - MOMY(KE,i,j-1)/MAPF(i,j-1,1,I_XV) ) * RCDY(j) )
          enddo
          enddo
#ifdef DEBUG
          k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif
       end if

       !##### continuity equation #####

       PROFILE_START("hevi_mflx_z")
       ! at (x, y, w)
       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j = JJS, JJE
       do i = IIS-1, IIE
          mflx_hi(KS-1,i,j,ZDIR) = 0.0_RP
          do k = KS, KE-1
#ifdef DEBUG
             call CHECK( __LINE__, MOMX(k+1,i  ,j) )
             call CHECK( __LINE__, MOMX(k+1,i-1,j) )
             call CHECK( __LINE__, MOMX(k  ,i  ,j) )
             call CHECK( __LINE__, MOMX(k  ,i+1,j) )
             call CHECK( __LINE__, MOMY(k+1,i,j) )
             call CHECK( __LINE__, MOMY(k+1,i,j-1) )
             call CHECK( __LINE__, MOMY(k  ,i,j) )
             call CHECK( __LINE__, MOMY(k  ,i,j-1) )
#endif
             mflx_hi(k,i,j,ZDIR) = J13G(k,i,j,I_XYW) * 0.25_RP / MAPF(i,j,2,I_XY) &
                                 * ( MOMX(k+1,i,j)+MOMX(k+1,i-1,j) &
                                   + MOMX(k  ,i,j)+MOMX(k  ,i-1,j) ) &
                                 + J23G(k,i,j,I_XYW) * 0.25_RP / MAPF(i,j,1,I_XY) &
                                 * ( MOMY(k+1,i,j)+MOMY(k+1,i,j-1) &
                                   + MOMY(k  ,i,j)+MOMY(k  ,i,j-1) ) &
                                 + GSQRT(k,i,j,I_XYW) / ( MAPF(i,j,1,I_XY)*MAPF(i,j,2,I_XY) ) * num_diff(k,i,j,I_DENS,ZDIR)
          enddo
          mflx_hi(KE  ,i,j,ZDIR) = 0.0_RP
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif
       PROFILE_STOP("hevi_mflx_z")

       PROFILE_START("hevi_mflx_x")
       iss = max(IIS-1,IS-IFS_OFF)
       iee = min(IIE,IEH)
       ! at (u, y, z)
       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j = JJS, JJE
       do i = iss, iee
       do k = KS, KE
#ifdef DEBUG
          call CHECK( __LINE__, GSQRT(k,i,j,I_UYZ) )
          call CHECK( __LINE__, MOMX(k,i,j) )
          call CHECK( __LINE__, num_diff(k,i,j,I_DENS,XDIR) )
#endif
          mflx_hi(k,i,j,XDIR) = GSQRT(k,i,j,I_UYZ) / MAPF(i,j,2,I_UY) &
               * ( MOMX(k,i,j) + num_diff(k,i,j,I_DENS,XDIR) )
       enddo
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif
       PROFILE_STOP("hevi_mflx_x")

       ! at (x, v, z)
       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j = max(JJS-1,JS-JFS_OFF), min(JJE,JEH)
       do i = IIS, IIE
       do k = KS, KE
#ifdef DEBUG
          call CHECK( __LINE__, GSQRT(k,i,j,I_XVZ) )
          call CHECK( __LINE__, MOMY(k,i,j) )
          call CHECK( __LINE__, num_diff(k,i,j,I_DENS,YDIR) )
#endif
          mflx_hi(k,i,j,YDIR) = GSQRT(k,i,j,I_XVZ) / MAPF(i,j,1,I_XV) &
               * ( MOMY(k,i,j) + num_diff(k,i,j,I_DENS,YDIR) )
       enddo
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif

       !--- update density
       !$omp parallel do private(i,j,k,advch) OMP_SCHEDULE_ collapse(2)
       PROFILE_START("hevi_sr")
       do j = JJS, JJE
       do i = IIS, IIE
       do k = KS, KE
#ifdef DEBUG
          call CHECK( __LINE__, DENS0(k,i,j) )
          call CHECK( __LINE__, mflx_hi(k  ,i  ,j  ,XDIR) )
          call CHECK( __LINE__, mflx_hi(k  ,i-1,j  ,XDIR) )
          call CHECK( __LINE__, mflx_hi(k  ,i  ,j  ,YDIR) )
          call CHECK( __LINE__, mflx_hi(k  ,i  ,j-1,YDIR) )
          call CHECK( __LINE__, DENS_t(k,i,j) )
#endif
          advch = - ( ( mflx_hi(k,i,j,ZDIR)-mflx_hi(k-1,i  ,j,  ZDIR) ) * RCDZ(k) &
                    + ( mflx_hi(k,i,j,XDIR)-mflx_hi(k  ,i-1,j,  XDIR) ) * RCDX(i) &
                    + ( mflx_hi(k,i,j,YDIR)-mflx_hi(k  ,i,  j-1,YDIR) ) * RCDY(j) ) &
                * MAPF(i,j,1,I_XY) * MAPF(i,j,2,I_XY) / GSQRT(k,i,j,I_XYZ)
          Sr(k,i,j) =  advch + DENS_t(k,i,j)
#ifdef HIST_TEND
          if ( lhist ) advch_t(k,i,j,I_DENS) = advch
#endif
       enddo
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif
       PROFILE_STOP("hevi_sr")

       !##### momentum equation (z) #####

       PROFILE_START("hevi_dens_qflxhi_z")
       ! at (x, y, z)
       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j = JJS, JJE
       do i = IIS, IIE
#ifdef DEBUG
          call CHECK( __LINE__, VELZ(KS  ,i,j) )
          call CHECK( __LINE__, VELZ(KS+1,i,j) )
          call CHECK( __LINE__, MOMZ(KS  ,i,j) )
          call CHECK( __LINE__, MOMZ(KS+1,i,j) )
          call CHECK( __LINE__, num_diff(KS+1,i,j,I_MOMZ,ZDIR) )
          call CHECK( __LINE__, VELZ(KE-2,i,j) )
          call CHECK( __LINE__, VELZ(KE-1,i,j) )
          call CHECK( __LINE__, MOMZ(KE-2,i,j) )
          call CHECK( __LINE__, MOMZ(KE-1,i,j) )
          call CHECK( __LINE__, num_diff(KE-1,i,j,I_MOMZ,ZDIR) )
#endif
          qflx_hi(KS,i,j,ZDIR) = 0.0_RP
          qflx_hi(KS+1,i,j,ZDIR) = J33G * 0.5_RP * ( VELZ(KS+1,i,j)+VELZ(KS,i,j) ) &
                                        * ( FACT_N * ( MOMZ(KS+1,i,j)+MOMZ(KS,i,j) ) &
                                          + FACT_F * ( MOMZ(KS+2,i,j)            ) ) &
                                 + GSQRT(KS+1,i,j,I_XYZ) * num_diff(KS+1,i,j,I_MOMZ,ZDIR)
          do k = KS+2, KE-2
#ifdef DEBUG
             call CHECK( __LINE__, VELZ(k  ,i,j) )
             call CHECK( __LINE__, VELZ(k-1,i,j) )
             call CHECK( __LINE__, MOMZ(k-2,i,j) )
             call CHECK( __LINE__, MOMZ(k-1,i,j) )
             call CHECK( __LINE__, MOMZ(k  ,i,j) )
             call CHECK( __LINE__, MOMZ(k+1,i,j) )
             call CHECK( __LINE__, num_diff(k,i,j,I_MOMZ,ZDIR) )
#endif
             qflx_hi(k,i,j,ZDIR) = J33G * 0.5_RP * ( VELZ(k,i,j)+VELZ(k-1,i,j) ) &
                                 * ( FACT_N * ( MOMZ(k  ,i,j)+MOMZ(k-1,i,j) ) &
                                   + FACT_F * ( MOMZ(k+1,i,j)+MOMZ(k-2,i,j) ) ) &
                                 + GSQRT(k,i,j,I_XYZ) * num_diff(k,i,j,I_MOMZ,ZDIR)
          enddo
          qflx_hi(KE-1,i,j,ZDIR) = J33G * 0.5_RP * ( VELZ(KE-1,i,j)+VELZ(KE-2,i,j) ) &
                                        * ( FACT_N * ( MOMZ(KE-1,i,j)+MOMZ(KE-2,i,j) ) &
                                          + FACT_F * (                MOMZ(KE-3,i,j) ) ) &
                                 + GSQRT(KE-1,i,j,I_XYZ) * num_diff(KE-1,i,j,I_MOMZ,ZDIR)
          qflx_hi(KE,i,j,ZDIR) = 0.0_RP
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif
       PROFILE_STOP("hevi_dens_qflxhi_z")

       PROFILE_START("hevi_dens_qflxj")
       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j = JJS, JJE
       do i = IIS, IIE
#ifdef DEBUG
          call CHECK( __LINE__, VELZ(KS  ,i,j) )
          call CHECK( __LINE__, VELZ(KS+1,i,j) )
          call CHECK( __LINE__, MOMZ(KS  ,i,j) )
          call CHECK( __LINE__, MOMZ(KS+1,i,j) )
          call CHECK( __LINE__, VELZ(KE-2,i,j) )
          call CHECK( __LINE__, VELZ(KE-1,i,j) )
          call CHECK( __LINE__, MOMZ(KE-2,i,j) )
          call CHECK( __LINE__, MOMZ(KE-1,i,j) )
#endif
          qflx_J(KS,i,j) = 0.0_RP
          qflx_J(KS+1,i,j) = 0.5_RP &
                           * ( J13G(KS+1,i,j,I_XYZ) * ( VELX(KS+1,i,j)+VELX(KS+1,i-1,j) ) &
                             + J23G(KS+1,i,j,I_XYZ) * ( VELY(KS+1,i,j)+VELY(KS+1,i,j-1) ) ) &
                           * ( FACT_N * ( MOMZ(KS+1,i,j)+MOMZ(KS  ,i,j) ) &
                             + FACT_F * ( MOMZ(KS+2,i,j)                ) )
          do k = KS+2, KE-2
#ifdef DEBUG
             call CHECK( __LINE__, VELZ(k  ,i,j) )
             call CHECK( __LINE__, VELZ(k-1,i,j) )
             call CHECK( __LINE__, MOMZ(k-2,i,j) )
             call CHECK( __LINE__, MOMZ(k-1,i,j) )
             call CHECK( __LINE__, MOMZ(k  ,i,j) )
             call CHECK( __LINE__, MOMZ(k+1,i,j) )
#endif
             qflx_J(k,i,j) = 0.5_RP &
                           * ( J13G(k,i,j,I_XYZ) * ( VELX(k,i,j)+VELX(k,i-1,j) ) &
                             + J23G(k,i,j,I_XYZ) * ( VELY(k,i,j)+VELY(k,i,j-1) ) ) &
                           * ( FACT_N * ( MOMZ(k  ,i,j)+MOMZ(k-1,i,j) ) &
                             + FACT_F * ( MOMZ(k+1,i,j)+MOMZ(k-2,i,j) ) )
          enddo
          qflx_J(KE-1,i,j) = 0.5_RP &
                           * ( J13G(KE-1,i,j,I_XYZ) * ( VELX(KE-1,i,j)+VELX(KE-1,i-1,j) ) &
                             + J23G(KE-1,i,j,I_XYZ) * ( VELY(KE-1,i,j)+VELY(KE-1,i,j-1) ) ) &
                           * ( FACT_N * ( MOMZ(KE-1,i,j)+MOMZ(KE-2,i,j) ) &
                             + FACT_F * (                MOMZ(KE-3,i,j) ) )
          qflx_J(KE,i,j) = 0.0_RP
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif
       PROFILE_STOP("hevi_dens_qflxj")

       ! at (u, y, w)
       PROFILE_START("hevi_dens_qflxhi_x")
       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j = JJS,   JJE
       do i = IIS-1, IIE
       do k = KS, KE-1
#ifdef DEBUG
          call CHECK( __LINE__, VELX(k  ,i,j) )
          call CHECK( __LINE__, VELX(k+1,i,j) )
          call CHECK( __LINE__, MOMZ(k,i-1,j) )
          call CHECK( __LINE__, MOMZ(k,i  ,j) )
          call CHECK( __LINE__, MOMZ(k,i+1,j) )
          call CHECK( __LINE__, MOMZ(k,i+2,j) )
          call CHECK( __LINE__, num_diff(k,i,j,I_MOMZ,XDIR) )
#endif
          qflx_hi(k,i,j,XDIR) = GSQRT(k,i,j,I_UYW) &
                              * ( 0.5_RP * ( VELX(k+1,i,j)+VELX(k,i,j) ) &
                                * ( FACT_N * ( MOMZ(k,i+1,j)+MOMZ(k,i  ,j) ) &
                                  + FACT_F * ( MOMZ(k,i+2,j)+MOMZ(k,i-1,j) ) ) &
                                + num_diff(k,i,j,I_MOMZ,XDIR) / MAPF(i,j,2,I_UY) )
       enddo
       enddo
       enddo
       PROFILE_STOP("hevi_dens_qflxhi_x")
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif

       ! at (x, v, w)
       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j = JJS-1, JJE
       do i = IIS,   IIE
       do k = KS, KE-1
#ifdef DEBUG
          call CHECK( __LINE__, VELY(k  ,i,j) )
          call CHECK( __LINE__, VELY(k+1,i,j) )
          call CHECK( __LINE__, MOMZ(k,i,j-1) )
          call CHECK( __LINE__, MOMZ(k,i,j  ) )
          call CHECK( __LINE__, MOMZ(k,i,j+1) )
          call CHECK( __LINE__, MOMZ(k,i,j+2) )
          call CHECK( __LINE__, num_diff(k,i,j,I_MOMZ,YDIR) )
#endif
          qflx_hi(k,i,j,YDIR) = GSQRT(k,i,j,I_XVW) &
                              * ( 0.5_RP * ( VELY(k+1,i,j)+VELY(k,i,j) ) &
                                * ( FACT_N * ( MOMZ(k,i,j+1)+MOMZ(k,i,j  ) ) &
                                  + FACT_F * ( MOMZ(k,i,j+2)+MOMZ(k,i,j-1) ) ) &
                                + num_diff(k,i,j,I_MOMZ,YDIR) / MAPF(i,j,1,I_XV) )
       enddo
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif

       !--- update momentum(z)
       PROFILE_START("hevi_sw")
       !$omp parallel do private(i,j,k,advcv,advch,cf,div) OMP_SCHEDULE_ collapse(2)
       do j = JJS, JJE
       do i = IIS, IIE
       do k = KS, KE-1
#ifdef DEBUG
          call CHECK( __LINE__, qflx_hi(k+1,i  ,j  ,ZDIR) )
          call CHECK( __LINE__, qflx_hi(k  ,i  ,j  ,ZDIR) )
          call CHECK( __LINE__, qflx_J (k+1,i  ,j)        )
          call CHECK( __LINE__, qflx_J (k  ,i  ,j)        )
          call CHECK( __LINE__, qflx_hi(k  ,i  ,j  ,XDIR) )
          call CHECK( __LINE__, qflx_hi(k  ,i-1,j  ,XDIR) )
          call CHECK( __LINE__, qflx_hi(k  ,i  ,j  ,YDIR) )
          call CHECK( __LINE__, qflx_hi(k  ,i  ,j-1,YDIR) )
          call CHECK( __LINE__, DDIV(k  ,i,j) )
          call CHECK( __LINE__, DDIV(k+1,i,j) )
          call CHECK( __LINE__, MOMZ0(k,i,j) )
          call CHECK( __LINE__, MOMZ_t(k,i,j) )
#endif
          advcv = -   ( qflx_hi(k+1,i,j,ZDIR) - qflx_hi(k,i  ,j  ,ZDIR) ) * RFDZ(k)
          advch = - ( ( qflx_J (k+1,i,j)      - qflx_J (k,i  ,j  )      ) * RFDZ(k) &
                    + ( qflx_hi(k  ,i,j,XDIR) - qflx_hi(k,i-1,j  ,XDIR) ) * RCDX(i) &
                    + ( qflx_hi(k  ,i,j,YDIR) - qflx_hi(k,i  ,j-1,YDIR) ) * RCDY(j) ) &
                * MAPF(i,j,1,I_XY) * MAPF(i,j,2,I_XY)
          cf = 0.0_RP
          div = divdmp_coef * dtrk * ( DDIV(k+1,i,j)-DDIV(k,i,j) ) * FDZ(k) ! divergence damping
          Sw(k,i,j) = ( advcv + advch ) / GSQRT(k,i,j,I_XYW) &
                    + div + MOMZ_t(k,i,j)
#ifdef HIST_TEND
          if ( lhist ) then
             advcv_t(k,i,j,I_MOMZ) = advcv / GSQRT(k,i,j,I_XYW)
             advch_t(k,i,j,I_MOMZ) = advch / GSQRT(k,i,j,I_XYW)
             cf_t(k,i,j,1) = cf
             ddiv_t(k,i,j,1) = div
          end if
#endif
       enddo
       enddo
       enddo
       PROFILE_STOP("hevi_sw")
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif


       !##### Thermodynamic Equation #####

       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j = JJS-2, JJE+2
       do i = IIS-2, IIE+2
       do k = KS, KE
#ifdef DEBUG
          call CHECK( __LINE__, RHOT(k,i,j) )
          call CHECK( __LINE__, DENS(k,i,j) )
#endif
          POTT(k,i,j) = RHOT(k,i,j) / DENS(k,i,j)
       enddo
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif

       ! at (x, y, w)
       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j = JJS, JJE
       do i = IIS, IIE
          tflx_hi(KS-1,i,j,ZDIR) = 0.0_RP
          tflx_hi(KS  ,i,j,ZDIR) = mflx_hi(KS  ,i,j,ZDIR) * 0.5_RP * ( POTT(KS+1,i,j) + POTT(KS  ,i,j) ) &
                                 + GSQRT(KS,i,j,I_XYW) / (MAPF(i,j,1,I_XY)*MAPF(i,j,2,I_XY)) * num_diff(KS  ,i,j,I_RHOT,ZDIR)
          do k = KS+1, KE-2
#ifdef DEBUG
             call CHECK( __LINE__, mflx_hi(k,i,j,ZDIR) )
             call CHECK( __LINE__, POTT(k,i-1,j) )
             call CHECK( __LINE__, POTT(k,i  ,j) )
             call CHECK( __LINE__, POTT(k,i+1,j) )
             call CHECK( __LINE__, POTT(k,i+1,j) )
             call CHECK( __LINE__, num_diff(k,i,j,I_RHOT,XDIR) )
#endif
             tflx_hi(k,i,j,ZDIR) = mflx_hi(k,i,j,ZDIR) &
                                 * ( FACT_N * ( POTT(k+1,i,j) + POTT(k  ,i,j) ) &
                                   + FACT_F * ( POTT(k+2,i,j) + POTT(k-1,i,j) ) ) &
                                 + GSQRT(k,i,j,I_XYW) / (MAPF(i,j,1,I_XY)*MAPF(i,j,2,I_XY)) * num_diff(k,i,j,I_RHOT,ZDIR)
          enddo
          tflx_hi(KE-1,i,j,ZDIR) = mflx_hi(KE-1,i,j,ZDIR) * 0.5_RP * ( POTT(KE  ,i,j) + POTT(KE-1,i,j) ) &
                                 + GSQRT(KE-1,i,j,ZDIR) / (MAPF(i,j,1,I_XY)*MAPF(i,j,2,I_XY)) * num_diff(KE-1,i,j,I_RHOT,ZDIR)
          tflx_hi(KE  ,i,j,ZDIR) = 0.0_RP
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif

       ! at (u, y, z)
       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       iss = max(IIS-1,IS-IFS_OFF)
       iee = min(IIE,IEH)
       do j = JJS,   JJE
       do i = iss, iee
       do k = KS, KE
#ifdef DEBUG
          call CHECK( __LINE__, mflx_hi(k,i,j,XDIR) )
          call CHECK( __LINE__, POTT(k,i-1,j) )
          call CHECK( __LINE__, POTT(k,i  ,j) )
          call CHECK( __LINE__, POTT(k,i+1,j) )
          call CHECK( __LINE__, POTT(k,i+1,j) )
          call CHECK( __LINE__, num_diff(k,i,j,I_RHOT,XDIR) )
#endif
          tflx_hi(k,i,j,XDIR) = mflx_hi(k,i,j,XDIR) &
                                * ( FACT_N * ( POTT(k,i+1,j)+POTT(k,i  ,j) ) &
                                  + FACT_F * ( POTT(k,i+2,j)+POTT(k,i-1,j) ) ) &
                                + GSQRT(k,i,j,I_UYZ) / MAPF(i,j,2,I_UY) * num_diff(k,i,j,I_RHOT,XDIR)
       enddo
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif
       ! at (x, v, z)
       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j = max(JJS-1,JS-JFS_OFF), min(JJE,JEH)
       do i = IIS,   IIE
       do k = KS, KE
#ifdef DEBUG
          call CHECK( __LINE__, mflx_hi(k,i,j,YDIR) )
          call CHECK( __LINE__, POTT(k,i,j-1) )
          call CHECK( __LINE__, POTT(k,i,j  ) )
          call CHECK( __LINE__, POTT(k,i,j+1) )
          call CHECK( __LINE__, POTT(k,i,j+2) )
          call CHECK( __LINE__, num_diff(k,i,j,I_RHOT,YDIR) )
#endif
          tflx_hi(k,i,j,YDIR) = mflx_hi(k,i,j,YDIR) &
                                * ( FACT_N * ( POTT(k,i,j+1)+POTT(k,i,j  ) ) &
                                  + FACT_F * ( POTT(k,i,j+2)+POTT(k,i,j-1) ) ) &
                                + GSQRT(k,i,j,I_XVZ) / MAPF(i,j,1,I_XV) * num_diff(k,i,j,I_RHOT,YDIR)
       enddo
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif

       PROFILE_START("hevi_st")
       !$omp parallel do private(i,j,k,advch) OMP_SCHEDULE_ collapse(2)
       do j = JJS, JJE
       do i = IIS, IIE
       do k = KS, KE
#ifdef DEBUG
          call CHECK( __LINE__, tflx_hi(k  ,i  ,j  ,ZDIR) )
          call CHECK( __LINE__, tflx_hi(k-1,i  ,j  ,ZDIR) )
          call CHECK( __LINE__, tflx_hi(k  ,i  ,j  ,XDIR) )
          call CHECK( __LINE__, tflx_hi(k  ,i-1,j  ,XDIR) )
          call CHECK( __LINE__, tflx_hi(k  ,i  ,j  ,YDIR) )
          call CHECK( __LINE__, tflx_hi(k  ,i  ,j-1,YDIR) )
          call CHECK( __LINE__, RHOT_t(k,i,j) )
#endif
          advch = - ( ( tflx_hi(k,i,j,ZDIR) - tflx_hi(k-1,i  ,j  ,ZDIR) ) * RCDZ(k) &
                    + ( tflx_hi(k,i,j,XDIR) - tflx_hi(k  ,i-1,j  ,XDIR) ) * RCDX(i) &
                    + ( tflx_hi(k,i,j,YDIR) - tflx_hi(k  ,i  ,j-1,YDIR) ) * RCDY(j) ) &
                  * MAPF(i,j,1,I_XY) * MAPF(i,j,2,I_XY) / GSQRT(k,i,j,I_XYZ)
          St(k,i,j) = advch + RHOT_t(k,i,j)
#ifdef HIST_TEND
          if ( lhist ) then
             advch_t(k,i,j,I_RHOT) = advch
          end if
#endif
       enddo
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif
       PROFILE_STOP("hevi_st")

       ! implicit solver

       PROFILE_START("hevi_solver")

!OCL INDEPENDENT(solve_*)
       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j = JJS, JJE
       do i = IIS, IIE
          PT(KS  ,i,j) = ( POTT(KS+1,i,j) + POTT(KS  ,i,j) ) * 0.5_RP
          do k = KS+1, KE-2
             PT(k,i,j) = FACT_N * ( POTT(k+1,i,j) + POTT(k  ,i,j) ) &
                       + FACT_F * ( POTT(k+2,i,j) + POTT(k-1,i,j) )
          end do
          PT(KE-1,i,j) = ( POTT(KE  ,i,j) + POTT(KE-1,i,j) ) * 0.5_RP
          do k = KS, KE
             CPRES(k,i,j) = &
#ifdef DRY
                  kappa        * PRES(k,i,j) /                  RHOT(k,i,j)
#else
                  CPtot(k,i,j) * PRES(k,i,j) / ( CVtot(k,i,j) * RHOT(k,i,j) )
#endif
             A(k,i,j) = dtrk**2 * J33G * RCDZ(k) * CPRES(k,i,j) * J33G / GSQRT(k,i,j,I_XYZ)
          end do
          B = GRAV * dtrk**2 * J33G / ( CDZ(KS+1) + CDZ(KS) )
          F1(KS,i,j) =        - ( PT(KS+1,i,j) * RFDZ(KS) *   A(KS+1,i,j)             + B ) / GSQRT(KS,i,j,I_XYW)
          F2(KS,i,j) = 1.0_RP + ( PT(KS  ,i,j) * RFDZ(KS) * ( A(KS+1,i,j)+A(KS,i,j) )     ) / GSQRT(KS,i,j,I_XYW)
          do k = KS+1, KE-2
             B = GRAV * dtrk**2 * J33G / ( CDZ(k+1) + CDZ(k) )
             F1(k,i,j) =        - ( PT(k+1,i,j) * RFDZ(k) *   A(k+1,i,j)            + B ) / GSQRT(k,i,j,I_XYW)
             F2(k,i,j) = 1.0_RP + ( PT(k  ,i,j) * RFDZ(k) * ( A(k+1,i,j)+A(k,i,j) )     ) / GSQRT(k,i,j,I_XYW)
             F3(k,i,j) =        - ( PT(k-1,i,j) * RFDZ(k) *              A(k,i,j)   - B ) / GSQRT(k,i,j,I_XYW)
          end do
          B = GRAV * dtrk**2 * J33G / ( CDZ(KE) + CDZ(KE-1) )
          F2(KE-1,i,j) = 1.0_RP + ( PT(KE-1,i,j) * RFDZ(KE-1) * ( A(KE,i,j)+A(KE-1,i,j) )    ) / GSQRT(KE-1,i,j,I_XYW)
          F3(KE-1,i,j) =        - ( PT(KE-2,i,j) * RFDZ(KE-1) *             A(KE-1,i,j)  - B ) / GSQRT(KE-1,i,j,I_XYW)
          do k = KS, KE-1
             pg = - ( DPRES(k+1,i,j) + CPRES(k+1,i,j)*dtrk*St(k+1,i,j) &
                    - DPRES(k  ,i,j) - CPRES(k  ,i,j)*dtrk*St(k  ,i,j) ) &
                    * RFDZ(k) * J33G / GSQRT(k,i,j,I_XYW) &
                  - 0.5_RP * GRAV &
                    * ( ( DENS(k+1,i,j) - REF_dens(k+1,i,j) ) &
                      + ( DENS(k  ,i,j) - REF_dens(k  ,i,j) ) &
                      + dtrk * ( Sr(k+1,i,j) + Sr(k,i,j) ) )
             C(k-KS+1,i,j) = MOMZ(k,i,j) + dtrk * ( pg + Sw(k,i,j) )
#if HIST_TEND
             if ( lhist ) pg_t(k,i,j,1) = pg
#endif
          end do

#ifdef HEVI_BICGSTAB
          call solve_bicgstab( &
               C(:,i,j),         & ! (inout)
               F1(:,i,j), F2(:,i,j), F3(:,i,j) ) ! (in)
#elif defined(HEVI_LAPACK)
          call solve_lapack( &
               C(:,i,j),        & ! (inout)
#ifdef DEBUG
               i, j,      & ! (in)
#endif
               F1(:,i,j), F2(:,i,j), F3(:,i,j) ) ! (in)
#else
          call solve_direct( &
               C(:,i,j),        & ! (inout)
               F1(:,i,j), F2(:,i,j), F3(:,i,j) ) ! (in)
#endif

          do k = KS, KE-1
#ifdef DEBUG_HEVI2HEVE
          ! for debug (change to explicit integration)
             C(k-KS+1,i,j) = MOMZ(k,i,j)
             mflx_hi(k,i,j,ZDIR) = mflx_hi(k,i,j,ZDIR) &
                                 + J33G * MOMZ(k,i,j) / ( MAPF(i,j,1,I_XY) * MAPF(i,j,2,I_XY) )
             MOMZ_RK(k,i,j) = MOMZ0(k,i,j) &
                  + dtrk*( &
                  - J33G * ( DPRES(k+1,i,j)-DPRES(k,i,j) ) * RFDZ(k) / GSQRT(k,i,j,i_XYW) &
                  - 0.5_RP * GRAV * ( (DENS(k,i,j)-REF_dens(k,i,j))+(DENS(k+1,i,j)-REF_dens(k+1,i,j)) ) &
                  + Sw(k,i,j) )
#else
             ! z-momentum flux
             mflx_hi(k,i,j,ZDIR) = mflx_hi(k,i,j,ZDIR) &
                                 + J33G * C(k-KS+1,i,j) / ( MAPF(i,j,1,I_XY) * MAPF(i,j,2,I_XY) )
             ! z-momentum
             MOMZ_RK(k,i,j) = MOMZ0(k,i,j) &
                            + ( C(k-KS+1,i,j) - MOMZ(k,i,j) )
#endif
          end do

          ! density and rho*theta
          advcv = - C(1,i,j)            * J33G * RCDZ(KS) / GSQRT(KS,i,j,I_XYZ) ! C(0) = 0
          DENS_RK(KS,i,j) = DENS0(KS,i,j) + dtrk * ( advcv + Sr(KS,i,j) )
#ifdef HIST_TEND
          if ( lhist ) advcv_t(KS,i,j,I_DENS) = advcv
#endif
          advcv = - C(1,i,j)*PT(KS,i,j) * J33G * RCDZ(KS) / GSQRT(KS,i,j,I_XYZ) ! C(0) = 0
          RHOT_RK(KS,i,j) = RHOT0(KS,i,j) + dtrk * ( advcv + St(KS,i,j) )
#ifdef HIST_TEND
          if ( lhist ) advcv_t(KS,i,j,I_RHOT) = advcv
#endif
          do k = KS+1, KE-1
             advcv = - ( C(k-KS+1,i,j)           - C(k-KS,i,j) ) &
                   * J33G * RCDZ(k) / GSQRT(k,i,j,I_XYZ)
             DENS_RK(k,i,j) = DENS0(k,i,j) + dtrk * ( advcv + Sr(k,i,j) )
#ifdef HIST_TEND
             if ( lhist ) advcv_t(k,i,j,I_DENS) = advcv
#endif
             advcv = - ( C(k-KS+1,i,j)*PT(k,i,j) - C(k-KS,i,j)*PT(k-1,i,j) ) &
                   * J33G * RCDZ(k) / GSQRT(k,i,j,I_XYZ)
             RHOT_RK(k,i,j) = RHOT0(k,i,j) + dtrk * ( advcv + St(k,i,j) )
#ifdef HIST_TEND
             if ( lhist ) advcv_t(k,i,j,I_RHOT) = advcv
#endif
          end do
          advcv = C(KE-KS,i,j)                * J33G * RCDZ(KE) / GSQRT(KE,i,j,I_XYZ) ! C(KE-KS+1) = 0
          DENS_RK(KE,i,j) = DENS0(KE,i,j) + dtrk * ( advcv + Sr(KE,i,j) )
#ifdef HIST_TEND
          if ( lhist ) advcv_t(KE,i,j,I_DENS) = advcv
#endif
          advcv = C(KE-KS,i,j) * PT(KE-1,i,j) * J33G * RCDZ(KE) / GSQRT(KE,i,j,I_XYZ) ! C(KE-KS+1) = 0
          RHOT_RK(KE,i,j) = RHOT0(KE,i,j) + dtrk * ( advcv + St(KE,i,j) )
#ifdef HIST_TEND
          if ( lhist ) advcv_t(KE,i,j,I_RHOT) = advcv
#endif

#ifdef DEBUG
          call check_equation( &
               C(:,i,j), &
               DENS(:,i,j), MOMZ(:,i,j), RHOT(:,i,j), PRES(:,i,j), &
               Sr(:,i,j), Sw(:,i,j), St(:,i,j), &
               J33G, GSQRT(:,i,j,:), &
#ifdef DRY
               kappa, &
#else
               CPtot(:,i,j), CVtot(:,i,j), &
#endif
               dtrk, i, j )
#endif

       end do
       end do
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif

       PROFILE_STOP("hevi_solver")


       !##### momentum equation (x) #####

       ! at (u, y, w)
       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j = JJS, JJE
       do i = IIS, IIE
#ifdef DEBUG
          call CHECK( __LINE__, VELZ(KS,i  ,j) )
          call CHECK( __LINE__, VELZ(KS,i+1,j) )
          call CHECK( __LINE__, DENS(KS+1,i+1,j) )
          call CHECK( __LINE__, DENS(KS+1,i  ,j) )
          call CHECK( __LINE__, DENS(KS  ,i+1,j) )
          call CHECK( __LINE__, DENS(KS  ,i  ,j) )
          call CHECK( __LINE__, MOMX(KS+1,i,j) )
          call CHECK( __LINE__, MOMX(KS  ,i,j) )
          call CHECK( __LINE__, num_diff(KS,i,j,I_MOMX,ZDIR) )
          call CHECK( __LINE__, VELZ(KE-1,i  ,j) )
          call CHECK( __LINE__, VELZ(KE-1,i+1,j) )
          call CHECK( __LINE__, DENS(KE  ,i+1,j) )
          call CHECK( __LINE__, DENS(KE  ,i  ,j) )
          call CHECK( __LINE__, DENS(KE-1,i+1,j) )
          call CHECK( __LINE__, DENS(KE-1,i  ,j) )
          call CHECK( __LINE__, MOMX(KE-1,i,j) )
          call CHECK( __LINE__, MOMX(KE  ,i,j) )
          call CHECK( __LINE__, num_diff(KE-1,i,j,I_MOMX,ZDIR) )
#endif
          qflx_hi(KS-1,i,j,ZDIR) = 0.0_RP
          qflx_hi(KS  ,i,j,ZDIR) = J33G * 0.25_RP * ( VELZ(KS,i+1,j) + VELZ(KS,i,j) ) &
                                 * ( MOMX(KS+1,i,j)+MOMX(KS,i,j) ) &
                                 + GSQRT(KS,i,j,I_UYW) * num_diff(KS  ,i,j,I_MOMX,ZDIR)
          do k = KS+1, KE-2
#ifdef DEBUG
             call CHECK( __LINE__, VELZ(k,i  ,j) )
             call CHECK( __LINE__, VELZ(k,i+1,j) )
             call CHECK( __LINE__, DENS(k+1,i+1,j) )
             call CHECK( __LINE__, DENS(k+1,i  ,j) )
             call CHECK( __LINE__, DENS(k  ,i+1,j) )
             call CHECK( __LINE__, DENS(k  ,i  ,j) )
             call CHECK( __LINE__, MOMX(k-1,i,j) )
             call CHECK( __LINE__, MOMX(k  ,i,j) )
             call CHECK( __LINE__, MOMX(k+1,i,j) )
             call CHECK( __LINE__, MOMX(k+2,i,j) )
             call CHECK( __LINE__, num_diff(k,i,j,I_MOMX,ZDIR) )
#endif
             qflx_hi(k,i,j,ZDIR) = J33G * 0.5_RP * ( VELZ(k,i+1,j) + VELZ(k,i,j) ) &
                                 * ( FACT_N * ( MOMX(k+1,i,j)+MOMX(k  ,i,j) ) &
                                   + FACT_F * ( MOMX(k+2,i,j)+MOMX(k-1,i,j) ) ) &
                                 + GSQRT(k,i,j,I_UYW) * num_diff(k,i,j,I_MOMX,ZDIR)
          enddo
          qflx_hi(KE-1,i,j,ZDIR) = J33G * 0.25_RP * ( VELZ(KE-1,i+1,j) + VELZ(KE-1,i,j) ) &
                                 * ( MOMX(KE,i,j)+MOMX(KE-1,i,j) ) &
                                 + GSQRT(KE-1,i,j,I_UYW) * num_diff(KE-1,i,j,I_MOMX,ZDIR)
          qflx_hi(KE,i,j,ZDIR) = 0.0_RP
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif

       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j = JJS, JJE
       do i = IIS, IIE
          qflx_J(KS-1,i,j) = 0.0_RP
          qflx_J(KS  ,i,j) = &
               ( J13G(KS,i,j,I_UYW) * 0.5_RP   * ( VELX(KS+1,i,j)+VELX(KS,i,j) ) &
               + J23G(KS,i,j,I_UYW) * 0.125_RP * ( VELY(KS+1,i+1,j  )+VELY(KS,i+1,j  ) &
                                                 + VELY(KS+1,i  ,j  )+VELY(KS,i  ,j  ) &
                                                 + VELY(KS+1,i+1,j-1)+VELY(KS,i+1,j-1) &
                                                 + VELY(KS+1,i  ,j-1)+VELY(KS,i  ,j-1) ) ) &
               * 0.5_RP * ( MOMX(KS+1,i,j)+MOMX(KS,i,j) )
          do k = KS+1, KE-2
             qflx_J(k,i,j) = &
                  ( J13G(k,i,j,I_UYW) * 0.5_RP   * ( VELX(k+1,i,j)+VELX(k,i,j) ) &
                  + J23G(k,i,j,I_UYW) * 0.125_RP * ( VELY(k+1,i+1,j  )+VELY(k,i+1,j  ) &
                                                   + VELY(k+1,i  ,j  )+VELY(k,i  ,j  ) &
                                                   + VELY(k+1,i+1,j-1)+VELY(k,i+1,j-1) &
                                                   + VELY(k+1,i  ,j-1)+VELY(k,i  ,j-1) ) ) &
                  * ( FACT_N * ( MOMX(k+1,i,j)+MOMX(k  ,i,j) ) &
                    + FACT_F * ( MOMX(k+2,i,j)+MOMX(k-1,i,j) ) )
          enddo
          qflx_J(KE-1,i,j) = &
               ( J13G(KE-1,i,j,I_UYW) * 0.5_RP   * ( VELX(KE,i,j)+VELX(KE-1,i,j) ) &
               + J23G(KE-1,i,j,I_UYW) * 0.125_RP * ( VELY(KE,i+1,j  )+VELY(KE-1,i+1,j  ) &
                                                   + VELY(KE,i  ,j  )+VELY(KE-1,i  ,j  ) &
                                                   + VELY(KE,i+1,j-1)+VELY(KE-1,i+1,j-1) &
                                                   + VELY(KE,i  ,j-1)+VELY(KE-1,i  ,j-1) ) ) &
               * 0.5_RP * ( MOMX(KE,i,j)+MOMX(KE-1,i,j) )
          qflx_J(KE  ,i,j) = 0.0_RP
       enddo
       enddo

       ! at (x, y, z)
       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j = JJS, JJE
       do i = IIS, IIE+1
       do k = KS, KE
#ifdef DEBUG
          call CHECK( __LINE__, VELX(k,i  ,j) )
          call CHECK( __LINE__, VELX(k,i-1,j) )
          call CHECK( __LINE__, MOMX(k,i-2,j) )
          call CHECK( __LINE__, MOMX(k,i-1,j) )
          call CHECK( __LINE__, MOMX(k,i  ,j) )
          call CHECK( __LINE__, MOMX(k,i+1,j) )
          call CHECK( __LINE__, num_diff(k,i,j,I_MOMX,XDIR) )
#endif
          qflx_hi(k,i,j,XDIR) = GSQRT(k,i,j,I_XYZ) &
                              * ( 0.5_RP * ( VELX(k,i,j)+VELX(k,i-1,j) ) &
                                * ( FACT_N * ( MOMX(k,i  ,j)+MOMX(k,i-1,j) ) &
                                  + FACT_F * ( MOMX(k,i+1,j)+MOMX(k,i-2,j) ) ) &
                                + num_diff(k,i,j,I_MOMX,XDIR) / MAPF(i,j,2,I_XY) )
       enddo
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif
       ! at (u, v, z)
       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j = JJS-1, JJE
       do i = IIS,   IIE
       do k = KS, KE
#ifdef DEBUG
          call CHECK( __LINE__, VELY(k,i+1,j) )
          call CHECK( __LINE__, VELY(k,i  ,j) )
          call CHECK( __LINE__, MOMX(k,i,j-1) )
          call CHECK( __LINE__, MOMX(k,i,j  ) )
          call CHECK( __LINE__, MOMX(k,i,j+1) )
          call CHECK( __LINE__, MOMX(k,i,j+2) )
          call CHECK( __LINE__, num_diff(k,i,j,I_MOMX,YDIR) )
#endif
          qflx_hi(k,i,j,YDIR) = GSQRT(k,i,j,I_UVZ) &
                              * ( 0.5_RP * ( VELY(k,i+1,j)+VELY(k,i,j) ) &
                                * ( FACT_N * ( MOMX(k,i,j+1)+MOMX(k,i,j  ) ) &
                                  + FACT_F * ( MOMX(k,i,j+2)+MOMX(k,i,j-1) ) ) &
                                + num_diff(k,i,j,I_MOMX,YDIR) / MAPF(i,j,1,I_XY) )
       enddo
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif

       !--- update momentum(x)
       PROFILE_START("hevi_momx")
       iee = min(IIE,IEH)
       !$omp parallel do private(i,j,k,advch,advcv,pg,cf,div) OMP_SCHEDULE_ collapse(2)
       do j = JJS, JJE
       do i = IIS, iee
       do k = KS, KE
#ifdef DEBUG
          call CHECK( __LINE__, qflx_hi(k  ,i  ,j  ,ZDIR) )
          call CHECK( __LINE__, qflx_hi(k-1,i  ,j  ,ZDIR) )
          call CHECK( __LINE__, qflx_hi(k  ,i+1,j  ,XDIR) )
          call CHECK( __LINE__, qflx_hi(k  ,i  ,j  ,XDIR) )
          call CHECK( __LINE__, qflx_hi(k  ,i  ,j  ,YDIR) )
          call CHECK( __LINE__, qflx_hi(k  ,i  ,j-1,YDIR) )
          call CHECK( __LINE__, PRES(k,i+1,j) )
          call CHECK( __LINE__, PRES(k,i  ,j) )
          call CHECK( __LINE__, CORIOLI(1,i  ,j) )
          call CHECK( __LINE__, CORIOLI(1,i+1,j) )
          call CHECK( __LINE__, VELY(k,i  ,j  ) )
          call CHECK( __LINE__, VELY(k,i+1,j  ) )
          call CHECK( __LINE__, VELY(k,i  ,j-1) )
          call CHECK( __LINE__, VELY(k,i+1,j-1) )
          call CHECK( __LINE__, DDIV(k,i+1,j) )
          call CHECK( __LINE__, DDIV(k,i  ,j) )
          call CHECK( __LINE__, MOMX0(k,i,j) )
#endif
          advcv = -   ( qflx_hi(k,i  ,j,ZDIR) - qflx_hi(k-1,i,j  ,ZDIR) ) * RCDZ(k)
          advch = - ( ( qflx_J (k,i  ,j)      - qflx_J (k-1,i,j)        ) * RCDZ(k) &
                    + ( qflx_hi(k,i+1,j,XDIR) - qflx_hi(k  ,i,j  ,XDIR) ) * RFDX(i) &
                    + ( qflx_hi(k,i  ,j,YDIR) - qflx_hi(k  ,i,j-1,YDIR) ) * RCDY(j) ) &
                * MAPF(i,j,1,I_UY) * MAPF(i,j,2,I_UY)
          pg = - ( ( GSQRT(k,i+1,j,I_XYZ) * DPRES(k,i+1,j) &
                   - GSQRT(k,i  ,j,I_XYZ) * DPRES(k,i  ,j) ) * RFDX(i) &
                 + ( J13G(k+1,i,j,I_UYZ) * ( DPRES(k+1,i+1,j)+DPRES(k+1,i,j) ) &
                   - J13G(k-1,i,j,I_UYZ) * ( DPRES(k-1,i+1,j)+DPRES(k-1,i,j) ) ) &
                   * 0.5_RP / ( FDZ(k+1)+FDZ(k) ) ) * MAPF(i,j,1,I_UY) ! pressure gradient force
          cf = 0.125_RP * ( CORIOLI(1,i,j)+CORIOLI(1,i+1,j) ) &
                        * ( MOMY(k,i,j)+MOMY(k,i+1,j)+MOMY(k,i,j-1)+MOMY(k,i+1,j-1) ) & ! coriolis force
             + 0.25_RP * MAPF(i,j,1,I_UY) * MAPF(i,j,2,I_UY) &
             * ( MOMY(k,i,j) + MOMY(k,i,j-1) + MOMY(k,i+1,j) + MOMY(k,i+1,j-1) ) &
             * ( ( MOMY(k,i,j) + MOMY(k,i,j-1) + MOMY(k,i+1,j) + MOMY(k,i+1,j-1) ) * 0.25_RP &
                 * ( 1.0_RP/MAPF(i+1,j,2,I_XY) - 1.0_RP/MAPF(i,j,2,I_XY) ) * RCDX(i) &
               - MOMX(k,i,j) &
                 * ( 1.0_RP/MAPF(i,j,1,I_UV) - 1.0_RP/MAPF(i,j-1,1,I_UV) ) * RFDY(j) ) &
             * 2.0_RP / ( DENS(k,i+1,j) + DENS(k,i,j) ) ! metric term
          div = divdmp_coef * dtrk * ( DDIV(k,i+1,j)/MAPF(i+1,j,2,I_XY) - DDIV(k,i,j)/MAPF(i,j,1,I_XY) ) &
              * MAPF(i,j,1,I_UY) * MAPF(i,j,2,I_UY) * FDX(i) ! divergence damping
          MOMX_RK(k,i,j) = MOMX0(k,i,j) &
               + dtrk * ( ( advcv + advch + pg ) / GSQRT(k,i,j,I_UYZ) + cf + div + MOMX_t(k,i,j) )
#ifdef HIST_TEND
          if ( lhist ) then
             advcv_t(k,i,j,I_MOMX) = advcv / GSQRT(k,i,j,I_UYZ)
             advch_t(k,i,j,I_MOMX) = advch / GSQRT(k,i,j,I_UYZ)
             pg_t(k,i,j,2) = pg / GSQRT(k,i,j,I_UYZ)
             cf_t(k,i,j,2) = cf
             ddiv_t(k,i,j,2) = div
          end if
#endif
       enddo
       enddo
       enddo
       PROFILE_STOP("hevi_momx")
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif

       !##### momentum equation (y) #####
       ! at (x, v, w)
       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j = JJS, JJE
       do i = IIS, IIE
#ifdef DEBUG
          call CHECK( __LINE__, VELZ(KS  ,i,j+1) )
          call CHECK( __LINE__, VELZ(KS  ,i,j  ) )
          call CHECK( __LINE__, DENS(KS+1,i+1,j) )
          call CHECK( __LINE__, DENS(KS+1,i  ,j) )
          call CHECK( __LINE__, DENS(KS  ,i+1,j) )
          call CHECK( __LINE__, DENS(KS  ,i  ,j) )
          call CHECK( __LINE__, MOMY(KS+1,i,j) )
          call CHECK( __LINE__, MOMY(KS  ,i,j) )
          call CHECK( __LINE__, num_diff(KS,i,j,I_MOMY,ZDIR) )
          call CHECK( __LINE__, VELZ(KE-1,i,j+1) )
          call CHECK( __LINE__, VELZ(KE-1,i,j  ) )
          call CHECK( __LINE__, DENS(KE  ,i+1,j) )
          call CHECK( __LINE__, DENS(KE  ,i  ,j) )
          call CHECK( __LINE__, DENS(KE-1,i+1,j) )
          call CHECK( __LINE__, DENS(KE-1,i  ,j) )
          call CHECK( __LINE__, MOMY(KE-1,i,j) )
          call CHECK( __LINE__, MOMY(KE  ,i,j) )
          call CHECK( __LINE__, num_diff(KE-1,i,j,I_MOMY,ZDIR) )
#endif
          qflx_hi(KS-1,i,j,ZDIR) = 0.0_RP
          qflx_hi(KS  ,i,j,ZDIR) = J33G * 0.25_RP &
                                 * ( VELZ(KS,i,j+1) + VELZ(KS,i,j) ) &
                                 * ( MOMY(KS+1,i,j) + MOMY(KS,i,j) ) &
                                 + GSQRT(KS,i,j,I_XVW) * num_diff(KS  ,i,j,I_MOMY,ZDIR)
          do k = KS+1, KE-2
#ifdef DEBUG
             call CHECK( __LINE__, VELZ(k,i,j+1) )
             call CHECK( __LINE__, VELZ(k,i,j  ) )
             call CHECK( __LINE__, DENS(k+1,i+1,j) )
             call CHECK( __LINE__, DENS(k+1,i  ,j) )
             call CHECK( __LINE__, DENS(k  ,i+1,j) )
             call CHECK( __LINE__, DENS(k  ,i  ,j) )
             call CHECK( __LINE__, MOMY(k-1,i,j) )
             call CHECK( __LINE__, MOMY(k  ,i,j) )
             call CHECK( __LINE__, MOMY(k+1,i,j) )
             call CHECK( __LINE__, MOMY(k+2,i,j) )
             call CHECK( __LINE__, num_diff(k,i,j,I_MOMY,ZDIR) )
#endif
             qflx_hi(k,i,j,ZDIR) = J33G * 0.5_RP * ( VELZ(k,i,j+1) + VELZ(k,i,j) ) &
                                 * ( FACT_N * ( MOMY(k+1,i,j) + MOMY(k  ,i,j) ) &
                                   + FACT_F * ( MOMY(k+2,i,j) + MOMY(k-1,i,j) ) ) &
                                 + GSQRT(k,i,j,I_XVW) * num_diff(k,i,j,I_MOMY,ZDIR)
          enddo
          qflx_hi(KE-1,i,j,ZDIR) = J33G * 0.25_RP &
                                 * ( VELZ(KE-1,i,j+1) + VELZ(KE-1,i,j) ) &
                                 * ( MOMY(KE,i,j) + MOMY(KE-1,i,j) ) &
                                 + GSQRT(KE-1,i,j,I_XVW) * num_diff(KE-1,i,j,I_MOMY,ZDIR)
          qflx_hi(KE  ,i,j,ZDIR) = 0.0_RP
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif

       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j = JJS, JJE
       do i = IIS, IIE
          qflx_J(KS-1,i,j) = 0.0_RP
          qflx_J(KS  ,i,j) = &
               ( J13G(KS,i,j,I_XVW) * 0.125_RP * ( VELX(KS+1,i  ,j+1)+VELX(KS,i  ,j+1) &
                                                 + VELX(KS+1,i-1,j+1)+VELX(KS,i-1,j+1) &
                                                 + VELX(KS+1,i  ,j  )+VELX(KS,i  ,j  ) &
                                                 + VELX(KS+1,i-1,j  )+VELX(KS,i-1,j  ) ) &
               + J23G(KS,i,j,I_XVW) * 0.5_RP   * ( VELY(KS+1,i,j)+VELY(KS,i,j) ) ) &
               * 0.5_RP * ( MOMY(KS+1,i,j)+MOMY(KS,i,j) )
          do k = KS+1, KE-2
             qflx_J(k,i,j) = &
                  ( J13G(k,i,j,I_XVW) * 0.125_RP * ( VELX(k+1,i  ,j+1)+VELX(k,i  ,j+1) &
                                                   + VELX(k+1,i-1,j+1)+VELX(k,i-1,j+1) &
                                                   + VELX(k+1,i  ,j  )+VELX(k,i  ,j  ) &
                                                   + VELX(k+1,i-1,j  )+VELX(k,i-1,j  ) ) &
                  + J23G(k,i,j,I_XVW) * 0.5_RP   * ( VELY(k+1,i,j)+VELY(k,i,j) ) ) &
                  * ( FACT_N * ( MOMY(k+1,i,j)+MOMY(k  ,i,j) ) &
                    + FACT_F * ( MOMY(k+2,i,j)+MOMY(k-1,i,j) ) )
          enddo
          qflx_J(KE-1,i,j) = &
               ( J23G(KE-1,i,j,I_XVW) * 0.5_RP   * ( VELY(KE,i,j)+VELY(KE-1,i,j) ) &
               + J13G(KE-1,i,j,I_XVW) * 0.125_RP * ( VELX(KE,i  ,j+1)+VELX(KE-1,i  ,j+1) &
                                                   + VELX(KE,i-1,j+1)+VELX(KE-1,i-1,j+1) &
                                                   + VELX(KE,i  ,j  )+VELX(KE-1,i  ,j  ) &
                                                   + VELX(KE,i-1,j  )+VELX(KE-1,i-1,j  ) ) ) &
               * 0.5_RP * ( MOMY(KE,i,j)+MOMY(KE-1,i,j) )
          qflx_J(KE  ,i,j) = 0.0_RP
       enddo
       enddo

       ! at (u, v, z)
       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j = JJS,   JJE
       do i = IIS-1, IIE
       do k = KS, KE
#ifdef DEBUG
          call CHECK( __LINE__, VELX(k,i,j+1) )
          call CHECK( __LINE__, VELX(k,i,j  ) )
          call CHECK( __LINE__, MOMY(k,i-1,j) )
          call CHECK( __LINE__, MOMY(k,i  ,j) )
          call CHECK( __LINE__, MOMY(k,i+1,j) )
          call CHECK( __LINE__, MOMY(k,i+2,j) )
          call CHECK( __LINE__, num_diff(k,i,j,I_MOMY,XDIR) )
#endif
          qflx_hi(k,i,j,XDIR) = GSQRT(k,i,j,I_UVZ) &
                              * ( 0.5_RP * ( VELX(k,i,j+1)+VELX(k,i,j) ) &
                                * ( FACT_N * ( MOMY(k,i+1,j)+MOMY(k,i  ,j) ) &
                                  + FACT_F * ( MOMY(k,i+2,j)+MOMY(k,i-1,j) ) ) &
                                + num_diff(k,i,j,I_MOMY,XDIR) / MAPF(i,j,2,I_UV) )
       enddo
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif
       ! at (x, y, z)
       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j = JJS, JJE+1
       do i = IIS, IIE
       do k = KS, KE
#ifdef DEBUG
          call CHECK( __LINE__, VELY(k,i,j  ) )
          call CHECK( __LINE__, VELY(k,i,j-1) )
          call CHECK( __LINE__, MOMY(k,i,j-2) )
          call CHECK( __LINE__, MOMY(k,i,j-1) )
          call CHECK( __LINE__, MOMY(k,i,j  ) )
          call CHECK( __LINE__, MOMY(k,i,j+1) )
          call CHECK( __LINE__, num_diff(k,i,j,I_MOMY,YDIR) )
#endif
          qflx_hi(k,i,j,YDIR) = GSQRT(k,i,j,I_XYZ) &
                              * ( 0.5_RP * ( VELY(k,i,j)+VELY(k,i,j-1) ) &
                                * ( FACT_N * ( MOMY(k,i,j  )+MOMY(k,i,j-1) ) &
                                  + FACT_F * ( MOMY(k,i,j+1)+MOMY(k,i,j-2) ) ) &
                                + num_diff(k,i,j,I_MOMY,YDIR) / MAPF(i,j,1,I_XY) )
       enddo
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif

       !--- update momentum(y)
       PROFILE_START("hevi_momy")
       !$omp parallel do private(i,j,k,advch,advcv,pg,cf,div) OMP_SCHEDULE_ collapse(2)
       do j = JJS, min(JJE,JEH)
       do i = IIS, IIE
       do k = KS, KE
#ifdef DEBUG
          call CHECK( __LINE__, qflx_hi(k  ,i  ,j  ,ZDIR) )
          call CHECK( __LINE__, qflx_hi(k-1,i  ,j  ,ZDIR) )
          call CHECK( __LINE__, qflx_hi(k  ,i  ,j  ,XDIR) )
          call CHECK( __LINE__, qflx_hi(k  ,i-1,j  ,XDIR) )
          call CHECK( __LINE__, qflx_hi(k  ,i  ,j+1,YDIR) )
          call CHECK( __LINE__, qflx_hi(k  ,i  ,j,YDIR) )
          call CHECK( __LINE__, PRES(k,i,j  ) )
          call CHECK( __LINE__, PRES(k,i,j+1) )
          call CHECK( __LINE__, CORIOLI(1,i,j  ) )
          call CHECK( __LINE__, CORIOLI(1,i,j+1) )
          call CHECK( __LINE__, VELX(k,i  ,j  ) )
          call CHECK( __LINE__, VELX(k,i  ,j+1) )
          call CHECK( __LINE__, VELX(k,i-1,j  ) )
          call CHECK( __LINE__, VELX(k,i-1,j+1) )
          call CHECK( __LINE__, DDIV(k,i,j+1) )
          call CHECK( __LINE__, DDIV(k,i,j  ) )
          call CHECK( __LINE__, MOMY_t(k,i,j) )
          call CHECK( __LINE__, MOMY0(k,i,j) )
#endif
          advcv = -   ( qflx_hi(k,i,j  ,ZDIR) - qflx_hi(k-1,i  ,j,ZDIR) ) * RCDZ(k)
          advch = - ( ( qflx_J (k,i,j  )      - qflx_J (k-1,i  ,j)      ) * RCDZ(k) &
                    + ( qflx_hi(k,i,j  ,XDIR) - qflx_hi(k  ,i-1,j,XDIR) ) * RCDX(i) &
                    + ( qflx_hi(k,i,j+1,YDIR) - qflx_hi(k  ,i  ,j,YDIR) ) * RFDY(j) ) &
                  * MAPF(i,j,1,I_XV) * MAPF(i,j,2,I_XV)
          pg = - ( ( GSQRT(k,i,j+1,I_XYZ) * DPRES(k,i,j+1) &
                   - GSQRT(k,i,j  ,I_XYZ) * DPRES(k,i,j  ) ) * RFDY(j) &
                 + ( J23G(k+1,i,j,I_XVZ) * ( DPRES(k+1,i,j+1)+DPRES(k+1,i,j) ) &
                   - J23G(k-1,i,j,I_XVZ) * ( DPRES(k-1,i,j+1)+DPRES(k-1,i,j) ) ) &
                   * 0.5_RP / ( FDZ(k+1)+FDZ(k) ) ) &
               * MAPF(i,j,2,I_XV) ! pressure gradient force
          cf = - 0.125_RP * ( CORIOLI(1,i  ,j+1)+CORIOLI(1,i  ,j) ) &
                          * ( MOMX(k,i,j+1)+MOMX(k,i,j)+MOMX(k,i-1,j+1)+MOMX(k,i-1,j) ) & ! coriolis force
             - 0.25_RP * MAPF(i,j,1,I_XV) * MAPF(i,j,2,I_XV) * GSQRT(k,i,j,I_XVZ) &
             * ( MOMX(k,i,j) + MOMX(k,i-1,j) + MOMX(k,i,j+1) + MOMX(k,i-1,j+1) )&
             * ( MOMY(k,i,j) &
                 * ( 1.0_RP/MAPF(i,j,2,I_UV) - 1.0_RP/MAPF(i-1,j,2,I_UV) ) * RCDX(i) &
               - 0.25_RP * ( MOMX(k,i,j)+MOMX(k,i-1,j)+MOMX(k,i,j+1)+MOMX(k,i-1,j+1) ) &
                 * ( 1.0_RP/MAPF(i,j+1,1,I_XY) - 1.0_RP/MAPF(i,j,1,I_XY) ) * RFDY(j) ) &
             * 2.0_RP / ( DENS(k,i+1,j) + DENS(k,i,j) ) ! metoric term
          div = divdmp_coef * dtrk * ( DDIV(k,i,j+1)/MAPF(i,j+1,1,I_XY) - DDIV(k,i,j)/MAPF(i,j,1,I_XY) ) &
              * MAPF(i,j,1,I_XV) * MAPF(i,j,2,I_XV) * FDY(j) ! divergence damping
          MOMY_RK(k,i,j) = MOMY0(k,i,j) &
                         + dtrk * ( ( advcv + advch + pg ) / GSQRT(k,i,j,I_XVZ) + cf + div + MOMY_t(k,i,j) )
#ifdef HIST_TEND
          if ( lhist ) then
             advcv_t(k,i,j,I_MOMY) = advcv / GSQRT(k,i,j,I_XVZ)
             advch_t(k,i,j,I_MOMY) = advch / GSQRT(k,i,j,I_XVZ)
             pg_t(k,i,j,3) = pg / GSQRT(k,i,j,I_XVZ)
             cf_t(k,i,j,3) = cf
             ddiv_t(k,i,j,3) = div
          end if
#endif
       enddo
       enddo
       enddo
       PROFILE_STOP("hevi_momy")
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif

    enddo
    enddo

#ifdef PROFILE_FIPP
       call fipp_stop()
#endif

#ifdef HIST_TEND
    if ( lhist ) then
       call HIST_in(advcv_t(:,:,:,I_DENS), 'DENS_t_advcv', 'tendency of density    (vert. advection)',    'kg/m3/s'   )
       call HIST_in(advcv_t(:,:,:,I_MOMZ), 'MOMZ_t_advcv', 'tendency of momentum z (vert. advection)',    'kg/m2/s2', zdim='half')
       call HIST_in(advcv_t(:,:,:,I_MOMX), 'MOMX_t_advcv', 'tendency of momentum x (vert. advection)',    'kg/m2/s2', xdim='half')
       call HIST_in(advcv_t(:,:,:,I_MOMY), 'MOMY_t_advcv', 'tendency of momentum y (vert. advection)',    'kg/m2/s2', ydim='half')
       call HIST_in(advcv_t(:,:,:,I_RHOT), 'RHOT_t_advcv', 'tendency of rho*theta  (vert. advection)',    'K kg/m3/s' )

       call HIST_in(advch_t(:,:,:,I_DENS), 'DENS_t_advch', 'tendency of density    (horiz. advection)',   'kg/m3/s'   )
       call HIST_in(advch_t(:,:,:,I_MOMZ), 'MOMZ_t_advch', 'tendency of momentum z (horiz. advection)',   'kg/m2/s2', zdim='half')
       call HIST_in(advch_t(:,:,:,I_MOMX), 'MOMX_t_advch', 'tendency of momentum x (horiz. advection)',   'kg/m2/s2', xdim='half')
       call HIST_in(advch_t(:,:,:,I_MOMY), 'MOMY_t_advch', 'tendency of momentum y (horiz. advection)',   'kg/m2/s2', ydim='half')
       call HIST_in(advch_t(:,:,:,I_RHOT), 'RHOT_t_advch', 'tendency of rho*theta  (horiz. advection)',   'K kg/m3/s' )

       call HIST_in(pg_t   (:,:,:,1),      'MOMZ_t_pg',    'tendency of momentum z (pressure gradient)',  'kg/m2/s2', zdim='half')
       call HIST_in(pg_t   (:,:,:,2),      'MOMX_t_pg',    'tendency of momentum x (pressure gradient)',  'kg/m2/s2', xdim='half')
       call HIST_in(pg_t   (:,:,:,3),      'MOMY_t_pg',    'tendency of momentum y (pressure gradient)',  'kg/m2/s2', ydim='half')

       call HIST_in(ddiv_t (:,:,:,1),      'MOMZ_t_ddiv',  'tendency of momentum z (divergence damping)', 'kg/m2/s2', zdim='half')
       call HIST_in(ddiv_t (:,:,:,2),      'MOMX_t_ddiv',  'tendency of momentum x (divergence damping)', 'kg/m2/s2', xdim='half')
       call HIST_in(ddiv_t (:,:,:,3),      'MOMY_t_ddiv',  'tendency of momentum y (divergence damping)', 'kg/m2/s2', ydim='half')

       call HIST_in(cf_t   (:,:,:,1),      'MOMX_t_cf',    'tendency of momentum x (coliolis force)',     'kg/m2/s2', xdim='half')
       call HIST_in(cf_t   (:,:,:,2),      'MOMY_t_cf',    'tendency of momentum y (coliolis force)',     'kg/m2/s2', ydim='half')
    endif
#endif

  end subroutine ATMOS_DYN_rk_hevi

#ifdef HEVI_BICGSTAB
!OCL SERIAL
  subroutine solve_bicgstab( &
       C,        & ! (inout)
       F1, F2, F3 ) ! (in)

    use scale_process, only: &
       PRC_MPIstop
    implicit none
    real(RP), intent(inout) :: C(KMAX-1)
    real(RP), intent(in)    :: F1(KA)
    real(RP), intent(in)    :: F2(KA)
    real(RP), intent(in)    :: F3(KA)

    real(RP) :: r0(KMAX-1)

    real(RP) :: M(3,KMAX-1)
    real(RP) :: p(KMAX-1)
    real(RP) :: ap(KMAX-1)
    real(RP) :: s(KMAX-1)
    real(RP) :: as(KMAX-1)
    real(RP) :: al, be, w

    real(RP), pointer :: r(:)
    real(RP), pointer :: rn(:)
    real(RP), pointer :: swap(:)
    real(RP), target :: v0(KMAX-1)
    real(RP), target :: v1(KMAX-1)
    real(RP) :: r0r
    real(RP) :: norm, error, epsilon

    integer :: k, iter

    epsilon = 0.1_RP**(RP-1)

    M(3,1) = F1(KS)
    M(2,1) = F2(KS)
    do k = KS+1, KE-2
       M(3,k-KS+1) = F1(k)
       M(2,k-KS+1) = F2(k)
       M(1,k-KS+1) = F3(k)
    enddo
    M(2,KE-KS) = F2(KE-1)
    M(1,KE-KS) = F3(KE-1)

    norm = 0.0_RP
    do k = 1, KMAX-1
       norm = norm + C(k)**2
    end do

    r  => v0
    rn => v1

    call mul_matrix( v1, M, C )

    do k = 1, KMAX-1
       r(k) = C(k) - v1(k)
       r0(k) = r(k)
       p(k) = r(k)
    end do

    r0r = r0(1) * r(1)
    do k = 2, KMAX-1
       r0r = r0r + r0(k)*r(k)
    enddo
    do iter = 1, KMAX-1
       error = 0.0_RP
       do k = 1, KMAX-1
          error = error + r(k)**2
       end do

!       if ( error < epsilon .or. error / norm < epsilon ) then
       if ( error/norm < epsilon ) then
#ifdef DEBUG
!          write(*,*) "Bi-CGSTAB converged:", iter
#endif
          exit
       end if

       call mul_matrix( ap, M, p )
       al = r0(1) * ap(1)
       do k = 2, KMAX-1
          al = al + r0(k)*ap(k)
       end do
       al = r0r / al ! (r0,r) / (r0,Mp)
       s(:) = r(:) - al*ap(:)
       call mul_matrix( as, M, s )
       be = as(1) * s(1)  ! be is used as just work variable here
       w =  as(1) * as(1)
       do k = 2, KMAX-1
          be = be + as(k)*s(k)
          w  = w  + as(k)*as(k)
       end do
       w = be / w ! (as,s) / (as,as)

       c(:) = c(:) + al*p(:) + w*s(:)
       rn(:) = s(:) - w*as(:)
       be = al/w / r0r
       r0r = r0(1) * rn(1)
       do k = 2, KMAX-1
          r0r = r0r + r0(k)*rn(k)
       end do
       be = be * r0r ! al/w * (r0,rn)/(r0,r)
       p(:) = rn(:) + be * ( p(:) - w*ap(:) )

       swap => rn
       rn => r
       r => swap
    end do

    if ( iter >= KMAX-1 ) then
       write(*,*) 'xxx [atmos_dyn_hevi] Bi-CGSTAB'
       write(*,*) 'xxx not converged', error, norm
       call PRC_MPIstop
    end if

    return
  end subroutine solve_bicgstab


  subroutine mul_matrix(V, M, C)
    implicit none
    real(RP), intent(out) :: V(KMAX-1)
    real(RP), intent(in)  :: M(5, KMAX-1)
    real(RP), intent(in)  :: C(KMAX-1)

    integer :: k

    V(1) = M(5,1)*C(3) + M(4,1)*C(2) + M(3,1)*C(1)
    V(2) = M(5,2)*C(4) + M(4,2)*C(3) + M(3,2)*C(2) + M(2,2)*C(1)
    do k = 3, KMAX-3
       V(k) = M(5,k)*C(k+2) + M(4,k)*C(k+1) + M(3,k)*C(k) + M(2,k)*C(k-1) + M(1,k)*C(k-2)
    enddo
    ! k = KE-2
    V(KMAX-2) = M(4,KMAX-2)*C(KMAX-1) + M(3,KMAX-2)*C(KMAX-2) + M(2,KMAX-2)*C(KMAX-3) + M(1,KMAX-2)*C(KMAX-4)
    ! k = KE-1
    V(KMAX-1) = M(3,KMAX-1)*C(KMAX-1) + M(2,KMAX-1)*C(KMAX-2) + M(1,KMAX-1)*C(KMAX-3)

    return
  end subroutine mul_matrix


#elif defined(HEVI_LAPACK)
!OCL SERIAL
  subroutine solve_lapack( &
       C,   & ! (inout)
#ifdef DEBUG
       i, j, &
#endif
       F1, F2, F3 ) ! (in)
    use scale_process, only: &
         PRC_MPIstop
    implicit none

    real(RP), intent(inout) :: C(KMAX-1)
    real(RP), intent(in)    :: F1(KA)
    real(RP), intent(in)    :: F2(KA)
    real(RP), intent(in)    :: F3(KA)
#ifdef DEBUG
    integer , intent(in)    :: i
    integer , intent(in)    :: j
#endif

    real(RP) :: M(NB*3+1,KMAX-1)
    integer  :: IPIV(KMAX-1)
    integer  :: INFO

    integer :: k

#ifdef DEBUG
    real(RP) :: M2(KMAX-1,KMAX-1)
    real(RP) :: C2(KMAX-1)
    real(RP) :: sum
#endif


    ! band matrix
    do k = KS+1, KE-2
       ! k-1, +1
       M(NB+1,k-KS+1) = F1(k-1)
       ! k, 0
       M(NB+2,k-KS+1) = F2(k)
       ! k+1, -1
       M(NB+3,k-KS+1) = F3(k+1)
    end do
    ! KS, 0
    M(NB+2,1    ) = F2(KS)
    ! KS+1, -1
    M(NB+3,1    ) = F3(KS+1)
    ! KE-2, +1
    M(NB+1,KE-KS) = F1(KE-2)
    ! KE-1, 0
    M(NB+2,KE-KS) = F2(KE-1)

#ifdef DEBUG
          M2(:,:) = 0.0_RP
          do k = 1, KMAX-1
             if (k>1) M2(k-1,k) = M(NB+3,k-1)
             M2(k,k) = M(NB+2,k)
             if (k<kmax-1) M2(k+1,k) = M(NB+1,k+1)
             c2(k) = c(k)
          end do
#endif
          if ( RP == DP ) then
             call DGBSV( KMAX-1, NB, NB, 1, M, NB*3+1, IPIV, C, KMAX-1, INFO)
          else
             call SGBSV( KMAX-1, NB, NB, 1, M, NB*3+1, IPIV, C, KMAX-1, INFO)
          end if
          ! C is (\rho w)^{n+1}
#ifdef DEBUG
          if ( INFO .ne. 0 ) then
             write(*,*) "DGBSV was failed", info
             call PRC_MPIstop
          end if

          do k = 1, KMAX-1
             sum = 0.0_RP
             if (k>1) sum = sum + M2(k-1,k)*c(k-1)
             sum = sum + M2(k,k)*c(k)
             if (k<kmax-1) sum = sum + M2(k+1,k)*c(k+1)
             if ( abs(sum-c2(k)) > 1E-10_RP ) then
                write(*,*) "sum is different"
                write(*,*) k+2, i, j, sum, c2(k)
!                write(*,*) M2(k-1:k+1,k), c(k-1:k+1)
                call PRC_MPIstop
             end if
          end do
#endif
  end subroutine solve_lapack
#else
!OCL SERIAL
  subroutine solve_direct( &
       C,         & ! (inout)
       F1, F2, F3 ) ! (in)

    use scale_process, only: &
       PRC_MPIstop
    implicit none
    real(RP), intent(inout) :: C(KMAX-1)
    real(RP), intent(in)    :: F1(KA)
    real(RP), intent(in)    :: F2(KA)
    real(RP), intent(in)    :: F3(KA)

    real(RP) :: e(KMAX-2)
    real(RP) :: f(KMAX-2)

    real(RP) :: rdenom

    integer :: k

    rdenom = 1.0_RP / F2(KS)
    e(1) = - F1(KS) * rdenom
    f(1) = C(1) * rdenom
    do k = 2, KMAX-2
       rdenom = 1.0_RP / ( F2(k+KS-1) + F3(k+KS-1) * e(k-1) )
       e(k) = - F1(k+KS-1) * rdenom
       f(k) = ( C(k) - F3(k+KS-1) * f(k-1) ) * rdenom
    end do

    ! C = \rho w
    C(KMAX-1) = ( C(KMAX-1) - F3(KE-1) * f(KMAX-2) ) &
              / ( F2(KE-1) + F3(KE-1) * e(KMAX-2) ) ! C(KMAX-1) = f(KMAX-1)
    do k = KMAX-2, 1, -1
       C(k) = e(k) * C(k+1) + f(k)
    end do

    return
  end subroutine solve_direct

#endif

#ifdef DEBUG
  subroutine check_equation( &
       VECT, &
       DENS, MOMZ, RHOT, PRES, &
       Sr, Sw, St, &
       J33G, G, &
#ifdef DRY
       kappa, &
#else
       CPtot, CVtot, &
#endif
       dt, i, j )
    use scale_const, only: &
         EPS => CONST_EPS, &
         GRAV => CONST_GRAV
    use scale_process, only: &
         PRC_MPIstop
    use scale_grid, only: &
         RCDZ => GRID_RCDZ, &
         RFDZ => GRID_RFDZ
    use scale_gridtrans, only: &
       I_XYZ, &
       I_XYW
    implicit none
    real(RP), intent(in) :: VECT(KMAX-1)
    real(RP), intent(in) :: DENS(KA)
    real(RP), intent(in) :: MOMZ(KA)
    real(RP), intent(in) :: RHOT(KA)
    real(RP), intent(in) :: PRES(KA)
    real(RP), intent(in) :: Sr(KA)
    real(RP), intent(in) :: Sw(KA)
    real(RP), intent(in) :: St(KA)
    real(RP), intent(in) :: J33G
    real(RP), intent(in) :: G(KA,8)
#ifdef DRY
    real(RP), intent(in) :: kappa
#else
    real(RP), intent(in) :: CPtot(KA)
    real(RP), intent(in) :: CVtot(KA)
#endif
    real(RP), intent(in) :: dt
    integer , intent(in) :: i
    integer , intent(in) :: j

    real(RP), parameter :: small = 1e-6_RP

    real(RP) :: MOMZ_N(KA)
    real(RP) :: DENS_N(KA)
    real(RP) :: RHOT_N(KA)
    real(RP) :: PRES_N(KA)

    real(RP) :: POTT(KA)
    real(RP) :: PT(KA)

#ifndef DRY
    real(RP) :: kappa
#endif

    real(RP) :: error, lhs, rhs
    integer :: k


    do k = KS, KE-1
       MOMZ_N(k) = VECT(k-KS+1)
    end do
    MOMZ_N(:KS-1) = 0.0_RP
    MOMZ_N(KE:) = 0.0_RP

    ! density
    do k = KS+1, KE-1
       DENS_N(k) = DENS(k) &
            + dt * ( - J33G * ( MOMZ_N(k) - MOMZ_N(k-1) ) * RCDZ(k) / G(k,I_XYZ) + Sr(k) )
    end do
    DENS_N(KS) = DENS(KS) &
         + dt * ( - J33G * MOMZ_N(KS) * RCDZ(KS) / G(KS,I_XYZ) + Sr(KS) )
    DENS_N(KE) = DENS(KE) &
         + dt * ( J33G * MOMZ_N(KE-1) * RCDZ(KE) / G(KE,I_XYZ) + Sr(KE) )

    ! rho*theta
    do k = KS, KE
       POTT(k) = RHOT(k) / DENS(k)
    end do
    do k = KS+1, KE-2
       PT(k) = ( 7.0_RP * ( POTT(k+1) + POTT(k  ) ) &
                 -        ( POTT(k+2) + POTT(k-1) ) ) / 12.0_RP
    end do
    PT(KS-1) = 0.0_RP
    PT(KS  ) = ( POTT(KS+1) + POTT(KS  ) ) * 0.5_RP
    PT(KE-1) = ( POTT(KE  ) + POTT(KE-1) ) * 0.5_RP
    PT(KE  ) = 0.0_RP
    do k = KS+1, KE-1
       RHOT_N(k) = RHOT(k) &
            + dt * ( - J33G * ( MOMZ_N(k)*PT(k) - MOMZ_N(k-1)*PT(k-1) ) * RCDZ(k) / G(k,I_XYZ) &
                     + St(k) )
    end do
    RHOT_N(KS) = RHOT(KS) &
         + dt * ( - J33G * MOMZ_N(KS)*PT(KS) * RCDZ(KS) / G(KS,I_XYZ) + St(KS) )
    RHOT_N(KE) = RHOT(KE) &
         + dt * ( J33G * MOMZ_N(KE-1)*PT(KE-1) * RCDZ(KE) / G(KE-1,I_XYZ) + St(KE) )


    do k = KS, KE
#ifndef DRY
       kappa = CPtot(k) / CVtot(k)
#endif
       PRES_N(k) = PRES(k) * ( 1.0_RP + kappa * ( RHOT_N(k) - RHOT(k) ) / RHOT(k) )
    end do

    do k = KS, KE
       lhs = ( DENS_N(k) - DENS(k) ) / dt
       rhs = - J33G * ( MOMZ_N(k) - MOMZ_N(k-1) ) * RCDZ(k) / G(k,I_XYZ) + Sr(k)
       if ( abs(lhs) < small ) then
          error = rhs
       else
          error = ( lhs - rhs ) / lhs
       end if
       if ( abs(error) > small ) then
          write(*,*)"HEVI: DENS error", k, i, j, error, lhs, rhs
          write(*,*)eps
          call PRC_MPIstop
       end if
    end do

    do k = KS, KE-1
       lhs = ( MOMZ_N(k) - MOMZ(k) ) / dt
       rhs = - J33G * ( PRES_N(k+1) - PRES_N(k) ) * RFDZ(k) / G(k,I_XYW) &
             - GRAV * ( DENS_N(k+1) + DENS_N(k) ) * 0.5_RP &
             + Sw(k)
       if ( abs(lhs) < small ) then
          error = rhs
       else
          error = ( lhs - rhs ) / lhs
       end if
       if ( abs(error) > small ) then
          write(*,*)"HEVI: MOMZ error", k, i, j, error, lhs, rhs
          write(*,*) MOMZ_N(k), MOMZ(k), dt
          write(*,*) - J33G * ( PRES(k+1) - PRES(k) ) * RFDZ(k) / G(k,I_XYW) &
             - GRAV * ( DENS(k+1) + DENS(k) ) * 0.5_RP &
             + Sw(k)
          call PRC_MPIstop
       end if
    end do

    do k = KS, KE
       lhs = ( RHOT_N(k) - RHOT(k) ) / dt
       rhs = - J33G * ( MOMZ_N(k)*PT(k) - MOMZ_N(k-1)*PT(k-1) ) * RCDZ(k) / G(k,I_XYZ) + St(k)
       if ( abs(lhs) < small ) then
          error = rhs
       else
          error = ( lhs - rhs ) / lhs
       end if
       if ( abs(error) > small ) then
          write(*,*)"HEVI: RHOT error", k, i, j, error, lhs, rhs
          call PRC_MPIstop
       end if
    end do

    return
  end subroutine check_equation
#endif

end module scale_atmos_dyn_rk_hevi
