!-------------------------------------------------------------------------------
!> module atmosphere / physics / pbl / mynn
!!
!! @par Description
!!          Boundary layer turbulence model
!!          Mellor-Yamada Nakanishi-Niino model
!!
!! @author Team SCALE
!!
!! @par Reference
!! @li Nakanishi and Niino, 2009:
!!     Development of an improved turbulence closure model for the atmospheric boundary layer.
!!     J. Meteorol. Soc. Japan, 87, 895-912
!!
!<
!-------------------------------------------------------------------------------
#include "scalelib.h"
module scale_atmos_phy_bl_mynn
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_io
  use scale_prof

#if defined DEBUG || defined QUICKDEBUG
  use scale_debug, only: &
     CHECK
#endif
  use scale_const, only: &
     UNDEF  => CONST_UNDEF, &
     IUNDEF => CONST_UNDEF2
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: ATMOS_PHY_BL_mynn_tracer_setup
  public :: ATMOS_PHY_BL_mynn_setup
  public :: ATMOS_PHY_BL_mynn_finalize
  public :: ATMOS_PHY_BL_mynn_tendency

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  integer,                public              :: ATMOS_PHY_BL_MYNN_NTRACER
  character(len=H_SHORT), public, allocatable :: ATMOS_PHY_BL_MYNN_NAME(:)
  character(len=H_LONG),  public, allocatable :: ATMOS_PHY_BL_MYNN_DESC(:)
  character(len=H_SHORT), public, allocatable :: ATMOS_PHY_BL_MYNN_UNITS(:)

  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  integer,  private, parameter :: I_TKE = 1
  integer,  private, parameter :: I_TSQ = 2
  integer,  private, parameter :: I_QSQ = 3
  integer,  private, parameter :: I_COV = 4

  real(RP), private, parameter :: OneOverThree = 1.0_RP / 3.0_RP
  real(RP), private, parameter :: LT_min       = 1.E-6_RP
  real(RP), private, parameter :: FLX_LIM_FACT = 0.5_RP

  real(RP), private            :: A1
  real(RP), private            :: A2
  real(RP), private, parameter :: B1 = 24.0_RP
  real(RP), private, parameter :: B2 = 15.0_RP
  real(RP), private            :: C1
  real(RP), private, parameter :: C2 = 0.75_RP
  real(RP), private, parameter :: C3 = 0.352_RP
!  real(RP), private, parameter :: C2 = 0.70_RP  ! MYNN2004
!  real(RP), private, parameter :: C3 = 0.323_RP ! MYNN2004
  real(RP), private, parameter :: C5 = 0.2_RP
  real(RP), private, parameter :: G1 = 0.235_RP
  real(RP), private            :: G2
  real(RP), private            :: F1
  real(RP), private            :: F2
  real(RP), private            :: Rf1
  real(RP), private            :: Rf2
  real(RP), private            :: Rfc
  real(RP), private            :: AF12 !> A1 F1 / A2 F2
  real(RP), private, parameter :: PrN = 0.74_RP

  real(RP), private, parameter :: zeta_min = -5.0_RP
  real(RP), private, parameter :: zeta_max =  2.0_RP

  real(RP), private            :: SQRT_2PI
  real(RP), private            :: RSQRT_2PI
  real(RP), private            :: RSQRT_2

  logical,  private            :: initialize

  real(RP), private            :: ATMOS_PHY_BL_MYNN_PBL_MAX    = 1.E+30_RP !> maximum height of the PBL
  real(RP), private            :: ATMOS_PHY_BL_MYNN_TKE_MIN    =  1.E-20_RP
  real(RP), private            :: ATMOS_PHY_BL_MYNN_N2_MAX     =  1.E1_RP
  real(RP), private            :: ATMOS_PHY_BL_MYNN_NU_MIN     = -1.E1_RP
  real(RP), private            :: ATMOS_PHY_BL_MYNN_NU_MAX     =  1.E4_RP
  real(RP), private            :: ATMOS_PHY_BL_MYNN_KH_MIN     = -1.E1_RP
  real(RP), private            :: ATMOS_PHY_BL_MYNN_KH_MAX     =  1.E4_RP
  real(RP), private            :: ATMOS_PHY_BL_MYNN_Lt_MAX     =  2000.0_RP
  real(RP), private            :: ATMOS_PHY_BL_MYNN_Sq_fact    = 3.0_RP
  logical,  private            :: ATMOS_PHY_BL_MYNN_init_TKE   = .false.
  logical,  private            :: ATMOS_PHY_BL_MYNN_similarity = .true.

  character(len=H_SHORT), private  :: ATMOS_PHY_BL_MYNN_LEVEL = "3" ! "2.5" or "3"

  namelist / PARAM_ATMOS_PHY_BL_MYNN / &
       ATMOS_PHY_BL_MYNN_PBL_MAX,  &
       ATMOS_PHY_BL_MYNN_N2_MAX,   &
       ATMOS_PHY_BL_MYNN_NU_MIN,   &
       ATMOS_PHY_BL_MYNN_NU_MAX,   &
       ATMOS_PHY_BL_MYNN_KH_MIN,   &
       ATMOS_PHY_BL_MYNN_KH_MAX,   &
       ATMOS_PHY_BL_MYNN_Lt_MAX,   &
       ATMOS_PHY_BL_MYNN_LEVEL,    &
       ATMOS_PHY_BL_MYNN_Sq_fact,  &
       ATMOS_PHY_BL_MYNN_init_TKE, &
       ATMOS_PHY_BL_MYNN_similarity


  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> ATMOS_PHY_BL_MYNN_tracer_setup
  !! Tracer Setup
  !<
  subroutine ATMOS_PHY_BL_MYNN_tracer_setup
    use scale_prc, only: &
       PRC_abort
    implicit none

    integer :: ierr

    LOG_NEWLINE
    LOG_INFO("ATMOS_PHY_BL_MYNN_tracer_setup",*) 'Tracer Setup'
    LOG_INFO("ATMOS_PHY_BL_MYNN_tracer_setup",*) 'Mellor-Yamada Nakanishi-Niino scheme'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_ATMOS_PHY_BL_MYNN,iostat=ierr)
    if( ierr > 0 ) then !--- fatal error
       LOG_ERROR("ATMOS_PHY_BL_MYNN_tracer_setup",*) 'Not appropriate names in namelist PARAM_ATMOS_PHY_BL_MYNN. Check!'
       call PRC_abort
    endif

    select case ( ATMOS_PHY_BL_MYNN_LEVEL )
    case ( "2.5" )
       ATMOS_PHY_BL_MYNN_NTRACER = 1
       allocate( ATMOS_PHY_BL_MYNN_NAME(1), ATMOS_PHY_BL_MYNN_DESC(1), ATMOS_PHY_BL_MYNN_UNITS(1) )
       ATMOS_PHY_BL_MYNN_NAME(:) = (/ 'TKE_MYNN' /)
       ATMOS_PHY_BL_MYNN_DESC(:) = (/ 'turbulent kinetic energy (MYNN)' /)
       ATMOS_PHY_BL_MYNN_UNITS(:) = (/ 'm2/s2' /)
    case ( "3" )
       ATMOS_PHY_BL_MYNN_NTRACER = 4
       allocate( ATMOS_PHY_BL_MYNN_NAME(4), ATMOS_PHY_BL_MYNN_DESC(4), ATMOS_PHY_BL_MYNN_UNITS(4) )
       ATMOS_PHY_BL_MYNN_NAME(:) = (/ 'TKE_MYNN', &
                                      'TSQ_MYNN', &
                                      'QSQ_MYNN', &
                                      'COV_MYNN' /)
       ATMOS_PHY_BL_MYNN_DESC(:) = (/ 'turbulent kinetic energy (MYNN)                                                         ', &
                                      'sub-grid variance of liquid water potential temperature (MYNN)                          ', &
                                      'sub-grid variance of total water content (MYNN)                                         ', &
                                      'sub-grid covariance of liquid water potential temperature and total water content (MYNN)' /)
       ATMOS_PHY_BL_MYNN_UNITS(:) = (/ 'm2/s2  ', &
                                       'K2     ', &
                                       'kg2/kg2', &
                                       'K kg   '  /)
    case default
       LOG_ERROR("ATMOS_PHY_BL_MYNN_tracer_setup",*) 'only level 2.5 and 3 are supported at this moment'
       call PRC_abort
    end select

    return
  end subroutine ATMOS_PHY_BL_MYNN_tracer_setup

  !-----------------------------------------------------------------------------
  !> ATMOS_PHY_BL_MYNN_setup
  !! Setup
  !<
  subroutine ATMOS_PHY_BL_MYNN_setup( &
       KA, KS, KE, IA, IS, IE, JA, JS, JE, &
       BULKFLUX_type, &
       TKE_MIN, PBL_MAX )
    use scale_prc, only: &
       PRC_abort
    use scale_const, only: &
       PI => CONST_PI
    implicit none

    integer,  intent(in) :: KA, KS, KE
    integer,  intent(in) :: IA, IS, IE
    integer,  intent(in) :: JA, JS, JE

    character(len=*), intent(in) :: BULKFLUX_type

    real(RP), intent(in), optional :: TKE_MIN
    real(RP), intent(in), optional :: PBL_MAX

    integer :: ierr
    integer :: k, i, j
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("ATMOS_PHY_BL_MYNN_setup",*) 'Setup'
    LOG_INFO("ATMOS_PHY_BL_MYNN_setup",*) 'Mellor-Yamada Nakanishi-Niino scheme'

    if ( present(TKE_MIN) ) ATMOS_PHY_BL_MYNN_TKE_MIN = TKE_MIN
    if ( present(PBL_MAX) ) ATMOS_PHY_BL_MYNN_PBL_MAX = PBL_MAX

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_ATMOS_PHY_BL_MYNN,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("ATMOS_PHY_BL_MYNN_setup",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("ATMOS_PHY_BL_MYNN_setup",*) 'Not appropriate names in namelist PARAM_ATMOS_PHY_BL_MYNN. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_ATMOS_PHY_BL_MYNN)

    A1        = B1 * (1.0_RP - 3.0_RP * G1) / 6.0_RP
    A2        = 1.0_RP / (3.0_RP * G1 * B1**(1.0_RP/3.0_RP) * PrN )
    C1        = G1 - 1.0_RP / ( 3.0_RP * A1 * B1**(1.0_RP/3.0_RP) )
    G2        = ( 2.0_RP * A1 * (3.0_RP - 2.0_RP * C2) + B2 * (1.0_RP - C3) ) / B1
    F1        = B1 * (G1 - C1) + 2.0_RP * A1 * (3.0_RP - 2.0_RP * C2) + 3.0_RP * A2 * (1.0_RP - C2) * (1.0_RP - C5)
    F2        = B1 * (G1 + G2) - 3.0_RP * A1 * (1.0_RP - C2)

    Rf1       = B1 * (G1 - C1) / F1
    Rf2       = B1 * G1 / F2
    Rfc       = G1 / (G1 + G2)

    AF12      = A1 * F1 / ( A2 * F2 )

    SQRT_2PI  = sqrt( 2.0_RP * PI )
    RSQRT_2PI = 1.0_RP / SQRT_2PI
    RSQRT_2   = 1.0_RP / sqrt( 2.0_RP )

    initialize = ATMOS_PHY_BL_MYNN_init_TKE

    select case ( BULKFLUX_type )
    case ( "B91", "B91W01" )
       ! do nothing
    case default
       ATMOS_PHY_BL_MYNN_similarity = .false.
    end select

    return
  end subroutine ATMOS_PHY_BL_MYNN_setup

  !-----------------------------------------------------------------------------
  !! Finalize
  !<
  subroutine ATMOS_PHY_BL_MYNN_finalize
    implicit none

    deallocate( ATMOS_PHY_BL_MYNN_NAME, ATMOS_PHY_BL_MYNN_DESC, ATMOS_PHY_BL_MYNN_UNITS )

    return
  end subroutine ATMOS_PHY_BL_MYNN_finalize

  !-----------------------------------------------------------------------------
  !> ATMOS_PHY_BL_MYNN_tendency
  !! calculate tendency by the virtical eddy viscosity
  !<
  subroutine ATMOS_PHY_BL_MYNN_tendency( &
       KA, KS, KE, IA, IS, IE, JA, JS, JE, &
       DENS, U, V, POTT, PROG,             &
       PRES, EXNER, N2,                    &
       QDRY, QV, Qw, POTL, POTV, SFC_DENS, &
       SFLX_MU, SFLX_MV, SFLX_SH, SFLX_QV, &
       us, ts, qs, RLmo,                   &
       CZ, FZ, F2H, dt_DP,                 &
       BULKFLUX_type,                      &
       RHOU_t, RHOV_t, RHOT_t, RHOQV_t,    &
       RPROG_t,                            &
       Nu, Kh, Qlp, cldfrac,               &
       Zi, SFLX_BUOY                       )
    use scale_const, only: &
       EPS     => CONST_EPS,    &
       GRAV    => CONST_GRAV,   &
       KARMAN  => CONST_KARMAN, &
       CPdry   => CONST_CPdry,  &
       EPSTvap => CONST_EPSTvap
    use scale_prc, only: &
       PRC_abort
    use scale_atmos_hydrometeor, only: &
       HYDROMETEOR_LHV => ATMOS_HYDROMETEOR_LHV, &
       CP_VAPOR
    use scale_matrix, only: &
       MATRIX_SOLVER_tridiagonal
    use scale_file_history, only: &
       FILE_HISTORY_in
    implicit none

    integer, intent(in) :: KA, KS, KE
    integer, intent(in) :: IA, IS, IE
    integer, intent(in) :: JA, JS, JE

    real(RP), intent(in) :: DENS    (KA,IA,JA) !> density
    real(RP), intent(in) :: U       (KA,IA,JA) !> zonal wind
    real(RP), intent(in) :: V       (KA,IA,JA) !> meridional wind
    real(RP), intent(in) :: POTT    (KA,IA,JA) !> potential temperature
    real(RP), intent(in) :: PROG    (KA,IA,JA,ATMOS_PHY_BL_MYNN_ntracer) !> prognostic variables (TKE, TSQ, QSQ, COV)
    real(RP), intent(in) :: PRES    (KA,IA,JA) !> pressure
    real(RP), intent(in) :: EXNER   (KA,IA,JA) !> Exner function
    real(RP), intent(in) :: N2      (KA,IA,JA) !> squared Brunt-Vaisala frequency
    real(RP), intent(in) :: QDRY    (KA,IA,JA) !> dry air
    real(RP), intent(in) :: QV      (KA,IA,JA) !> vapor
    real(RP), intent(in) :: Qw      (KA,IA,JA) !> total water content
    real(RP), intent(in) :: POTL    (KA,IA,JA) !> liquid water potential temp.
    real(RP), intent(in) :: POTV    (KA,IA,JA) !> virtual potential temp.
    real(RP), intent(in) :: SFC_DENS(   IA,JA) !> surface density
    real(RP), intent(in) :: SFLX_MU (   IA,JA) !> surface flux of zonal wind
    real(RP), intent(in) :: SFLX_MV (   IA,JA) !> surface flux of meridional wind
    real(RP), intent(in) :: SFLX_SH (   IA,JA) !> surface sensible heat flux
    real(RP), intent(in) :: SFLX_QV (   IA,JA) !> surface sensible QV flux
    real(RP), intent(in) :: us      (   IA,JA) !> friction velocity
    real(RP), intent(in) :: ts      (   IA,JA) !> temperature scale
    real(RP), intent(in) :: qs      (   IA,JA) !> moisture scale
    real(RP), intent(in) :: RLmo    (   IA,JA) !> inverse of Monin-Obukhov length

    real(RP), intent(in)  :: CZ(  KA,IA,JA)
    real(RP), intent(in)  :: FZ(0:KA,IA,JA)
    real(RP), intent(in)  :: F2H(KA,2,IA,JA)
    real(DP), intent(in)  :: dt_DP

    character(len=*), intent(in) :: BULKFLUX_type

    real(RP), intent(out) :: RHOU_t (KA,IA,JA) !> tendency of dens * u
    real(RP), intent(out) :: RHOV_t (KA,IA,JA) !> tendency of dens * v
    real(RP), intent(out) :: RHOT_t (KA,IA,JA) !> tendency of dens * pt
    real(RP), intent(out) :: RHOQV_t(KA,IA,JA) !> tendency of dens * qv
    real(RP), intent(out) :: RPROG_t(KA,IA,JA,ATMOS_PHY_BL_MYNN_ntracer) !> tenddency of dens * prognostic variables (TKE, TSQ, QSQ, COV)
    real(RP), intent(out) :: Nu     (KA,IA,JA) !> eddy viscosity coefficient @ half-level
    real(RP), intent(out) :: Kh     (KA,IA,JA) !> eddy diffusion coefficient @ half-level
    real(RP), intent(out) :: Qlp    (KA,IA,JA) !> liquid-water content in partial condensation
    real(RP), intent(out) :: cldfrac(KA,IA,JA) !> cloud fraction in partial condensation
    real(RP), intent(out) :: Zi        (IA,JA) !> depth of the boundary layer (not exact)
    real(RP), intent(out) :: SFLX_BUOY (IA,JA) !> surface flux of buoyancy: g / \Theta_0 <w \theta_v> @ surface

    real(RP) :: Ri   (KA,IA,JA) !> Richardson number
    real(RP) :: Pr   (KA,IA,JA) !> Plandtle number
    real(RP) :: prod (KA,IA,JA) !> TKE production term
    real(RP) :: diss (KA,IA,JA) !> TKE dissipation term
    real(RP) :: dudz2(KA,IA,JA) !> (du/dz)^2 + (dv/dz)^2
    real(RP) :: l    (KA,IA,JA) !> length scale L
    real(RP) :: rho_h(KA)       !> dens at the half level

    real(RP) :: flxU(0:KA,IA,JA) !> dens * w * u
    real(RP) :: flxV(0:KA,IA,JA) !> dens * w * v
    real(RP) :: flxT(0:KA,IA,JA) !> dens * w * pt
    real(RP) :: flxQ(0:KA,IA,JA) !> dens * w * qv

    real(RP) :: RHO   (KA) !> dens after updated
    real(RP) :: RHONu (KA) !> dens * Nu at the half level for level 2.5
    real(RP) :: RHOKh (KA) !> dens * Kh at the half level for level 2.5
    real(RP) :: N2_new(KA) !> squared Brunt-Baisala frequency
    real(RP) :: SFLX_PT    !> surface potential temperature flux * density
    real(RP) :: sm25  (KA) !> stability function for velocity for level 2.5
    real(RP) :: sh25  (KA) !> stability function for scalars for level 2.5
    real(RP) :: Nu_f  (KA) !> Nu at the full level
    real(RP) :: Kh_f  (KA) !> Kh at the full level
    real(RP) :: q     (KA) !> q
    real(RP) :: q2_2  (KA) !> q^2 for level 2
    real(RP) :: ac    (KA) !> \alpha_c
    real(RP) :: lq    (KA) !> L * q

    ! for level 3
    real(RP) :: tsq   (KA)
    real(RP) :: qsq   (KA)
    real(RP) :: cov   (KA)
    real(RP) :: dtsq  (KA)
    real(RP) :: dqsq  (KA)
    real(RP) :: dcov  (KA)
    real(RP) :: prod_t(KA)
    real(RP) :: prod_q(KA)
    real(RP) :: prod_c(KA)
    real(RP) :: diss_p(KA)
    real(RP) :: dtldz(KA)
    real(RP) :: dqwdz(KA)
    real(RP) :: betat(KA)
    real(RP) :: betaq(KA)
    real(RP) :: smp    (KA) !> stability function for velocity for the countergradient
    real(RP) :: f_smp  (KA) !> stability function for velocity for the countergradient (factor)
    real(RP) :: shpgh  (KA) !> stability function for scalars for the countergradient
    real(RP) :: f_shpgh(KA) !> stability function for scalars for the countergradient (factor)
    real(RP) :: tltv
    real(RP) :: qwtv
    real(RP) :: tvsq
    real(RP) :: tltv25 (KA)
    real(RP) :: qwtv25 (KA)
    real(RP) :: tvsq25 (KA)
    real(RP) :: tvsq_up(KA) !> upper limit of <\theta_v^2> - <\theta_v^2>2.5
    real(RP) :: tvsq_lo(KA) !> lower limit
    real(RP) :: wtl
    real(RP) :: wqw
    real(RP) :: gammat (KA)
    real(RP) :: gammaq (KA)
    real(RP) :: f_gamma(KA) !> - E_H / q^2 * GRAV / Theta_0
    real(RP) :: rlqsm_h(KA) !> DENS * L * q * SM' @ half level

    real(RP) :: flx(0:KA)
    real(RP) :: a(KA)
    real(RP) :: b(KA)
    real(RP) :: c(KA)
    real(RP) :: d(KA)
    real(RP) :: ap
    real(RP) :: phi_N(KA)
    real(RP) :: tke_P(KA)
    real(RP) :: dummy(KA)

    real(RP) :: CPtot

    real(RP) :: sf_t
    real(RP) :: us3
    real(RP) :: zeta
    real(RP) :: phi_m, phi_h

    real(RP) :: FDZ(KA)
    real(RP) :: CDZ(KA)
    real(RP) :: z1

    logical :: mynn_level3

    real(RP) :: dt

    real(RP) :: fmin
    integer  :: kmin

    real(RP) :: tmp, sw

    integer :: KE_PBL
    integer :: k, i, j
    integer :: nit, it
    !---------------------------------------------------------------------------

    dt = real(dt_DP, RP)

    LOG_PROGRESS(*) "atmosphere / physics / pbl / MYNN"

    mynn_level3 = ( ATMOS_PHY_BL_MYNN_LEVEL == "3" )


!OCL INDEPENDENT
    !$omp parallel do default(none) &
    !$omp OMP_SCHEDULE_ collapse(2) &
    !$omp shared(KA,KS,KE,IS,IE,JS,JE, &
    !$omp        EPS,GRAV,CPdry,CP_VAPOR,EPSTvap,UNDEF,RSQRT_2,SQRT_2PI,RSQRT_2PI, &
    !$omp        ATMOS_PHY_BL_MYNN_N2_MAX,ATMOS_PHY_BL_MYNN_TKE_MIN, &
    !$omp        ATMOS_PHY_BL_MYNN_NU_MIN,ATMOS_PHY_BL_MYNN_NU_MAX, &
    !$omp        ATMOS_PHY_BL_MYNN_KH_MIN,ATMOS_PHY_BL_MYNN_KH_MAX, &
    !$omp        ATMOS_PHY_BL_MYNN_Sq_fact,ATMOS_PHY_BL_MYNN_similarity, &
    !$omp        ATMOS_PHY_BL_MYNN_PBL_MAX, &
    !$omp        RHOU_t,RHOV_t,RHOT_t,RHOQV_t,RPROG_t,Nu,Kh,Qlp,cldfrac,Zi,SFLX_BUOY, &
    !$omp        DENS,PROG,U,V,POTT,PRES,QDRY,QV,Qw,POTV,POTL,EXNER,N2, &
    !$omp        SFC_DENS,SFLX_MU,SFLX_MV,SFLX_SH,SFLX_QV,us,ts,qs,RLmo, &
    !$omp        mynn_level3,initialize, &
    !$omp        CZ,FZ,F2H,dt, &
    !$omp        BULKFLUX_type, &
    !$omp        Ri,Pr,prod,diss,dudz2,l,flxU,flxV,flxT,flxQ) &
    !$omp private(N2_new,lq,sm25,sh25,rlqsm_h,Nu_f,Kh_f,q,q2_2,ac, &
    !$omp         SFLX_PT,CPtot,RHO,RHONu,RHOKh, &
    !$omp         smp,f_smp,shpgh,f_shpgh,tltv,qwtv,tvsq,tltv25,qwtv25,tvsq25,tvsq_up,tvsq_lo,wtl,wqw, &
    !$omp         dtldz,dqwdz,betat,betaq,gammat,gammaq,f_gamma, &
    !$omp         flx,a,b,c,d,ap,rho_h,phi_n,tke_P,sf_t,zeta,phi_m,phi_h,us3,CDZ,FDZ,z1, &
    !$omp         dummy, &
    !$omp         tsq,qsq,cov,dtsq,dqsq,dcov, &
    !$omp         prod_t,prod_q,prod_c,diss_p, &
    !$omp         fmin,kmin, &
    !$omp         sw,tmp, &
    !$omp         KE_PBL,k,i,j,it,nit)
    do j = JS, JE
    do i = IS, IE

       KE_PBL = KS+1
       do k = KS+2, KE-1
          if ( ATMOS_PHY_BL_MYNN_PBL_MAX >= CZ(k,i,j) - FZ(KS-1,i,j) ) then
             KE_PBL = k
          else
             exit
          end if
       end do

       if ( ATMOS_PHY_BL_MYNN_similarity ) then

          z1 = CZ(KS,i,j) - FZ(KS-1,i,j)

          zeta = min( max( z1 * RLmo(i,j), zeta_min ), zeta_max )

          select case ( BULKFLUX_type )
!!$          case ( 'B71' )
!!$             ! Businger et al. (1971)
!!$             if ( zeta >= 0 ) then
!!$                phi_m = 4.7_RP * zeta + 1.0_RP
!!$                phi_h = 4.7_RP * zeta + 0.74_RP
!!$             else
!!$                phi_m = 1.0_RP / sqrt(sqrt( 1.0_RP - 15.0_RP * zeta ))
!!$                phi_h = 0.47_RP / sqrt( 1.0_RP - 9.0_RP * zeta )
!!$             end if
          case ( 'B91', 'B91W01' )
             ! Beljaars and Holtslag (1991)
             if ( zeta >= 0 ) then
                tmp = - 2.0_RP / 3.0_RP * ( 0.35_RP * zeta - 6.0_RP ) * exp(-0.35_RP*zeta) * zeta
                phi_m = tmp + zeta + 1.0_RP
                phi_h = tmp + zeta * sqrt( 1.0_RP + 2.0_RP * zeta / 3.0_RP ) + 1.0_RP
             else
                if ( BULKFLUX_type == 'B91W01' ) then
                   ! Wilson (2001)
                   !tmp = (-zeta)**(2.0_RP/3.0_RP)
                   tmp = abs(zeta)**(2.0_RP/3.0_RP)
                   phi_m = 1.0_RP / sqrt( 1.0_RP + 3.6_RP * tmp )
                   phi_h = 0.95_RP / sqrt( 1.0_RP + 7.9_RP * tmp )
                else
                   !                      tmp = sqrt( 1.0_RP - 16.0_RP * zeta )
                   tmp = sqrt( 1.0_RP + 16.0_RP * abs(zeta) )
                   phi_m = 1.0_RP / sqrt(tmp)
                   phi_h = 1.0_RP / tmp
                end if
             end if
          end select

       end if

       do k = KS, KE_PBL
          FDZ(k) = CZ(k+1,i,j) - CZ(k  ,i,j)
       end do
       do k = KS, KE_PBL+1
          CDZ(k) = FZ(k  ,i,j) - FZ(k-1,i,j)
       end do

       call calc_vertical_differece( KA, KS, KE_PBL, &
                                     i, j, &
                                     U(:,i,j), V(:,i,j), POTL(:,i,j), & ! (in)
                                     Qw(:,i,j), QDRY(:,i,j),          & ! (in)
                                     CDZ(:), FDZ(:), F2H(:,:,i,j),    & ! (in)
                                     dudz2(:,i,j), dtldz(:), dqwdz(:) ) ! (out)

       us3 = us(i,j)**3
       CPtot = CPdry + SFLX_QV(i,j) * ( CP_VAPOR - CPdry )
       SFLX_PT = SFLX_SH(i,j) / ( CPtot * EXNER(KS,i,j) )

       if ( initialize ) then
          do k = KS, KE_PBL
             q(k) = 1e-10_RP
          end do
       else
          do k = KS, KE_PBL
             q(k) = sqrt( max( PROG(k,i,j,I_TKE), ATMOS_PHY_BL_MYNN_TKE_MIN ) * 2.0_RP )
          end do
       end if

       if ( initialize .or. (.not. mynn_level3) ) then
          ! estimate tsq, qsq, and cov

          do k = KS, KE_PBL
             n2_new(k) = min( max( N2(k,i,j), - ATMOS_PHY_BL_MYNN_N2_MAX ), ATMOS_PHY_BL_MYNN_N2_MAX )
             !n2_new(k) = GRAV * POTV(k,i,j) * dtldz(k)
             Ri(k,i,j) = n2_new(k) / dudz2(k,i,j)
          end do

          SFLX_BUOY(i,j) = - us3 * RLmo(i,j) / KARMAN

          ! length
          call get_length( &
               KA, KS, KE_PBL, &
               i, j,           &
               q(:), n2_new(:),           & ! (in)
               SFLX_BUOY(i,j), RLmo(i,j), & ! (in)
               CZ(:,i,j), FZ(:,i,j),      & ! (in)
               l(:,i,j)                   ) ! (out)

          call get_q2_level2( &
               KA, KS, KE_PBL, &
               dudz2(:,i,j), Ri(:,i,j), l(:,i,j), & ! (in)
               q2_2(:)                            ) ! (out)

          if ( initialize ) then
             do k = KS, KE_PBL
                q(k) = sqrt( q2_2(k) )
             end do
          end if

          do k = KS, KE_PBL
             ac(k) = min( q(k) / sqrt( q2_2(k) + 1e-20_RP ), 1.0_RP )
          end do

          call get_smsh( &
               KA, KS, KE_PBL,            & ! (in)
               i, j,                      & ! (in)
               q(:), ac(:),               & ! (in)
               l(:,i,j), n2_new(:),       & ! (in)
               POTV(:,i,j), dudz2(:,i,j), & ! (in)
               dtldz(:), dqwdz(:),        & ! (in)
               betat(:), betaq(:),        & ! (in) ! dummy
               .false., .false.,          & ! (in)
               tsq(:), qsq(:), cov(:),    & ! (inout)
               sm25(:), f_smp(:),         & ! (out) ! dummy
               sh25(:), f_shpgh(:),       & ! (out) ! dymmy
               f_gamma(:),                & ! (out) ! dummy
               tltv25(:), qwtv25(:),      & ! (out) ! dummy
               tvsq25(:),                 & ! (out) ! dummy
               tvsq_up(:), tvsq_lo(:)     ) ! (out) ! dummy

       else

          do k = KS, KE_PBL
             tsq(k) = max( PROG(k,i,j,I_TSQ), 0.0_RP )
             qsq(k) = max( PROG(k,i,j,I_QSQ), 0.0_RP )
             cov(k) = PROG(k,i,j,I_COV)
             cov(k) = sign( min( abs(cov(k)), sqrt(tsq(k) * qsq(k))), cov(k) )
          end do

       end if


       flx(KS-1  ) = 0.0_RP
       flx(KE_PBL) = 0.0_RP

       if ( initialize ) then
          nit = KE_PBL - 1
       else
          nit = 1
       end if

       do it = 1, nit

          call partial_condensation( KA, KS, KE_PBL, &
                                     PRES(:,i,j), POTT(:,i,j),  & ! (in)
                                     POTL(:,i,j), Qw(:,i,j),    & ! (in)
                                     QDRY(:,i,j), EXNER(:,i,j), & ! (in)
                                     tsq(:), qsq(:), cov(:),    & ! (in)
                                     betat(:), betaq(:),        & ! (out)
                                     Qlp(:,i,j), cldfrac(:,i,j) ) ! (out)

          ! update N2
          do k = KS, KE_PBL
             n2_new(k) = min(ATMOS_PHY_BL_MYNN_N2_MAX, &
                             GRAV * ( dtldz(k) * betat(k) + dqwdz(k) * betaq(k) ) / POTV(k,i,j) )
          end do

          do k = KS, KE_PBL
             Ri(k,i,j) = n2_new(k) / dudz2(k,i,j)
          end do

          SFLX_BUOY(i,j) = GRAV / POTV(KS,i,j) * ( betat(KS) * SFLX_PT + betaq(KS) * SFLX_QV(i,j) ) / SFC_DENS(i,j)

          ! length
          call get_length( &
               KA, KS, KE_PBL, &
               i, j,           &
               q(:), n2_new(:),           & ! (in)
               SFLX_BUOY(i,j), RLmo(i,j), & ! (in)
               CZ(:,i,j), FZ(:,i,j),      & ! (in)
               l(:,i,j)                   ) ! (out)

          call get_q2_level2( &
               KA, KS, KE_PBL, &
               dudz2(:,i,j), Ri(:,i,j), l(:,i,j), & ! (in)
               q2_2(:)                            ) ! (out)

          do k = KS, KE_PBL
             ac(k) = min( q(k) / sqrt( q2_2(k) + 1e-20_RP ), 1.0_RP )
          end do

          call get_smsh( &
               KA, KS, KE_PBL,            & ! (in)
               i, j,                      & ! (in)
               q(:), ac(:),               & ! (in)
               l(:,i,j), n2_new(:),       & ! (in)
               POTV(:,i,j), dudz2(:,i,j), & ! (in)
               dtldz(:), dqwdz(:),        & ! (in)
               betat(:), betaq(:),        & ! (in)
               mynn_level3,               & ! (in)
               initialize .and. it==1,    & ! (in)
               tsq(:), qsq(:), cov(:),    & ! (inout)
               sm25(:), f_smp(:),         & ! (out)
               sh25(:), f_shpgh(:),       & ! (out)
               f_gamma(:),                & ! (out)
               tltv25(:), qwtv25(:),      & ! (out)
               tvsq25(:),                 & ! (out)
               tvsq_up(:), tvsq_lo(:)     ) ! (out)

          do k = KS, KE_PBL
             lq(k) = l(k,i,j) * q(k)
             Nu_f(k) = lq(k) * sm25(k)
             Kh_f(k) = lq(k) * sh25(k)
          end do
!          if ( ATMOS_PHY_BL_MYNN_similarity ) then
!             Nu_f(KS) = KARMAN * z1 * us(i,j) / phi_m
!             Kh_f(KS) = KARMAN * z1 * us(i,j) / phi_h
!          end if

          do k = KS, KE_PBL-1
             Nu(k,i,j) = min( F2H(k,1,i,j) * Nu_f(k+1) + F2H(k,2,i,j) * Nu_f(k), &
                              ATMOS_PHY_BL_MYNN_NU_MAX )
             Kh(k,i,j) = min( F2H(k,1,i,j) * Kh_f(k+1) + F2H(k,2,i,j) * Kh_f(k), &
                              ATMOS_PHY_BL_MYNN_KH_MAX )
          end do

          do k = KS, KE_PBL-1
             sw = 0.5_RP - sign(0.5_RP, abs(Kh(k,i,j)) - EPS)
             Pr(k,i,j) = Nu(k,i,j) / ( Kh(k,i,j) + sw ) * ( 1.0_RP - sw ) &
                       + 1.0_RP * sw
          end do

          RHO(KS) = DENS(KS,i,j) + dt * SFLX_QV(i,j) / CDZ(KS)
          do k = KS+1, KE_PBL
             RHO(k) = DENS(k,i,j)
          end do

          do k = KS, KE_PBL-1
             rho_h(k) = F2H(k,1,i,j) * RHO(k+1) + F2H(k,2,i,j) * RHO(k)
             RHONu(k) = Nu(k,i,j) * rho_h(k)
             RHOKh(k) = Kh(k,i,j) * rho_h(k)
!             RHONu(k) = F2H(k,1,i,j) * RHO(k+1) * Nu_f(k+1) + F2H(k,2,i,j) * RHO(k) * Nu_f(k)
!             RHOKh(k) = F2H(k,1,i,j) * RHO(k+1) * Kh_f(k+1) + F2H(k,2,i,j) * RHO(k) * Kh_f(k)
          end do


          if ( mynn_level3 ) then

             ! production term calculated by explicit scheme
             do k = KS, KE_PBL
                tltv = betat(k) * tsq(k) + betaq(k) * cov(k)
                qwtv = betat(k) * cov(k) + betaq(k) * qsq(k)
                gammat(k) = f_gamma(k) * ( tltv - tltv25(k) )
                gammaq(k) = f_gamma(k) * ( qwtv - qwtv25(k) )

                wtl = - lq(k) * ( sh25(k) * dtldz(k) + gammat(k) )
                wqw = - lq(k) * ( sh25(k) * dqwdz(k) + gammaq(k) )
                prod_t(k) = - 2.0_RP * wtl * dtldz(k)
                prod_q(k) = - 2.0_RP * wqw * dqwdz(k)
                prod_c(k) = - wtl * dqwdz(k) - wqw * dtldz(k)
             end do

!!$          if ( ATMOS_PHY_BL_MYNN_similarity ) then
!!$             tmp = 2.0_RP * us(i,j) * phi_h / ( KARMAN * z1 )
!!$             tmp = tmp * ( zeta / ( z1 * RLmo(i,j) ) )**2 ! correspoindint to the limitter for zeta
!!$             ! TSQ
!!$             prod_t(KS) = tmp * ts(i,j)**2
!!$             ! QSQ
!!$             prod_q(KS) = tmp * ts(i,j) * qs(i,j)
!!$             ! COV
!!$             prod_c(KS) = tmp * qs(i,j)**2
!!$          end if

             ! diffusion (explicit)
             flx(KS-1)   = 0.0_RP
             flx(KE_PBL) = 0.0_RP
             do k = KS, KE_PBL-1
                flx(k) = RHONu(k) * ( tsq(k+1) - tsq(k) ) / FDZ(k)
             end do
             do k = KS, KE_PBL
                prod_t(k) = prod_t(k) + ( flx(k) - flx(k-1) ) / CDZ(k)
             end do
             do k = KS, KE_PBL-1
                flx(k) = RHONu(k) * ( qsq(k+1) - qsq(k) ) / FDZ(k)
             end do
             do k = KS, KE_PBL
                prod_q(k) = prod_q(k) + ( flx(k) - flx(k-1) ) / CDZ(k)
             end do
             do k = KS, KE_PBL-1
                flx(k) = RHONu(k) * ( cov(k+1) - cov(k) ) / FDZ(k)
             end do
             do k = KS, KE_PBL
                prod_c(k) = prod_c(k) + ( flx(k) - flx(k-1) ) / CDZ(k)
             end do

             call get_gamma_implicit( &
                  KA, KS, KE_PBL, &
                  i, j,       &
                  tsq(:), qsq(:), cov(:),          & ! (in)
                  dtldz(:), dqwdz(:), POTV(:,i,j), & ! (in)
                  prod_t(:), prod_q(:), prod_c(:), & ! (in)
                  betat(:), betaq(:),              & ! (in)
                  f_gamma(:), l(:,i,j), q(:),      & ! (in)
                  dt,                              & ! (in)
                  dtsq(:), dqsq(:), dcov(:)        ) ! (out)

             ! update
             do k = KS, KE_PBL
                tltv = betat(k) * ( tsq(k) + dtsq(k) ) + betaq(k) * ( cov(k) + dcov(k) )
                qwtv = betat(k) * ( cov(k) + dcov(k) ) + betaq(k) * ( qsq(k) + dqsq(k) )
                gammat(k) = f_gamma(k) * ( tltv - tltv25(k) )
                gammaq(k) = f_gamma(k) * ( qwtv - qwtv25(k) )

                tvsq = max( betat(k) * tltv + betaq(k) * qwtv, 0.0_RP )
                tvsq = tvsq - tvsq25(k)
                tvsq = min( max( tvsq, tvsq_lo(k) ), tvsq_up(k) )
                smp  (k) = f_smp  (k) * tvsq
                shpgh(k) = f_shpgh(k) * tvsq

                wtl = - lq(k) * ( sh25(k) * dtldz(k) + gammat(k) )
                wqw = - lq(k) * ( sh25(k) * dqwdz(k) + gammaq(k) )
                prod_t(k) = - 2.0_RP * wtl * dtldz(k)
                prod_q(k) = - 2.0_RP * wqw * dqwdz(k)
                prod_c(k) = - wtl * dqwdz(k) - wqw * dtldz(k)
             end do

          else

             do k = KS, KE_PBL
                smp   (k) = 0.0_RP
                shpgh (k) = 0.0_RP
                gammat(k) = 0.0_RP
                gammaq(k) = 0.0_RP
             end do

          end if



          ! time integration

          if ( it == nit ) then

             do k = KS, KE_PBL-1
                rlqsm_h(k) = rho_h(k) * ( F2H(k,1,i,j) * lq(k+1) * smp(k+1) + F2H(k,2,i,j) * lq(k) * smp(k) )
             end do

             ! dens * u

             ! countergradient flux
             do k = KS, KE_PBL-1
                flx(k) = - rlqsm_h(k) * ( U(k+1,i,j) - U(k,i,j) ) / FDZ(k)
             end do

             sf_t = SFLX_MU(i,j) / CDZ(KS)
             d(KS) = ( U(KS,i,j) * DENS(KS,i,j) + dt * ( sf_t - flx(KS) / CDZ(KS) ) ) / RHO(KS)
             do k = KS+1, KE_PBL
                d(k) = U(k,i,j) - dt * ( flx(k) - flx(k-1) ) / ( CDZ(k) * RHO(k) )
             end do
             c(KS) = 0.0_RP
             do k = KS, KE_PBL-1
                ap = - dt * RHONu(k) / FDZ(k)
                a(k) = ap / ( RHO(k) * CDZ(k) )
                b(k) = - a(k) - c(k) + 1.0_RP
                c(k+1) = ap / ( RHO(k+1) * CDZ(k+1) )
             end do
             a(KE_PBL) = 0.0_RP
             b(KE_PBL) = - c(KE_PBL) + 1.0_RP

             call MATRIX_SOLVER_tridiagonal( &
                  KA, KS, KE_PBL, &
                  a(:), b(:), c(:), d(:), & ! (in)
                  dummy(:)                ) ! (out)
!                  phi_n(:)                ) ! (out)

             phi_n(KS:KE_PBL) = dummy(KS:KE_PBL)
             RHOU_t(KS,i,j) = ( phi_n(KS) * RHO(KS) - U(KS,i,j) * DENS(KS,i,j) ) / dt - sf_t
             do k = KS+1, KE_PBL
                RHOU_t(k,i,j) = ( phi_n(k) - U(k,i,j) ) * RHO(k) / dt
             end do
             do k = KE_PBL+1, KE
                RHOU_t(k,i,j) = 0.0_RP
             end do
             flxU(KS-1,i,j) = 0.0_RP
             do k = KS, KE_PBL-1
                flxU(k,i,j) = flx(k) &
                            - RHONu(k) * ( phi_n(k+1) - phi_n(k) ) / FDZ(k)
             end do
             do k = KE_PBL, KE
                flxU(k,i,j) = 0.0_RP
             end do


             ! dens * v

             ! countergradient flux
             do k = KS, KE_PBL-1
                flx(k) = - rlqsm_h(k) * ( V(k+1,i,j) - V(k,i,j) ) / FDZ(k)
             end do

             sf_t = SFLX_MV(i,j) / CDZ(KS)
             d(KS) = ( V(KS,i,j) * DENS(KS,i,j) + dt * ( sf_t - flx(KS) / CDZ(KS) ) ) / RHO(KS)
             do k = KS+1, KE_PBL
                d(k) = V(k,i,j) - dt * ( flx(k) - flx(k-1) ) / ( CDZ(k) * RHO(k) )
             end do
             ! a,b,c is the same as those for the u

             call MATRIX_SOLVER_tridiagonal( &
                  KA, KS, KE_PBL, &
                  a(:), b(:), c(:), d(:), & ! (in)
                  phi_n(:)                ) ! (out)

             RHOV_t(KS,i,j) = ( phi_n(KS) * RHO(KS) - V(KS,i,j) * DENS(KS,i,j) ) / dt - sf_t
             do k = KS+1, KE_PBL
                RHOV_t(k,i,j) = ( phi_n(k) - V(k,i,j) ) * RHO(k) / dt
             end do
             do k = KE_PBL+1, KE
                RHOV_t(k,i,j) = 0.0_RP
             end do
             flxV(KS-1,i,j) = 0.0_RP
             do k = KS, KE_PBL-1
                flxV(k,i,j) = flx(k) &
                            - RHONu(k) * ( phi_n(k+1) - phi_n(k) ) / FDZ(k)
             end do
             do k = KE_PBL, KE
                flxV(k,i,j) = 0.0_RP
             end do


             ! dens * pott

             ! countergradient flux
             do k = KS, KE_PBL-1
                flx(k) = - ( F2H(k,1,i,j) * lq(k+1) * gammat(k+1) + F2H(k,2,i,j) * lq(k) * gammat(k) ) &
                       * rho_h(k)
             end do

             sf_t = SFLX_PT / CDZ(KS)
             ! assume that induced vapor has the same PT
             d(KS) = POTT(KS,i,j) + dt * ( sf_t - flx(KS) / CDZ(KS) ) / RHO(KS)
             do k = KS+1, KE_PBL
                d(k) = POTT(k,i,j) - dt * ( flx(k) - flx(k-1) ) / ( CDZ(k) * RHO(k) )
             end do

             c(KS) = 0.0_RP
             do k = KS, KE_PBL-1
                ap = - dt * RHOKh(k) / FDZ(k)
                a(k) = ap / ( RHO(k) * CDZ(k) )
                b(k) = - a(k) - c(k) + 1.0_RP
                c(k+1) = ap / ( RHO(k+1) * CDZ(k+1) )
             end do
             a(KE_PBL) = 0.0_RP
             b(KE_PBL) = - c(KE_PBL) + 1.0_RP

             call MATRIX_SOLVER_tridiagonal( &
                  KA, KS, KE_PBL, &
                  a(:), b(:), c(:), d(:), & ! (in)
                  phi_n(:)                ) ! (out)

             RHOT_t(KS,i,j) = ( phi_n(KS) - POTT(KS,i,j) ) * RHO(KS) / dt - sf_t
             do k = KS+1, KE_PBL
                RHOT_t(k,i,j) = ( phi_n(k) - POTT(k,i,j) ) * RHO(k) / dt
             end do
             do k = KE_PBL+1, KE
                RHOT_t(k,i,j) = 0.0_RP
             end do
             flxT(KS-1,i,j) = 0.0_RP
             do k = KS, KE_PBL-1
                flxT(k,i,j) = flx(k) &
                            - RHOKh(k) * ( phi_n(k+1) - phi_n(k) ) / FDZ(k)
             end do
             do k = KE_PBL, KE
                flxT(k,i,j) = 0.0_RP
             end do

             kmin = KS-1
             if ( flxT(KS,i,j) > 1E-4_RP ) then
                fmin = flxT(KS,i,j) / rho_h(KS)
                do k = KS+1, KE_PBL-2
                   tmp = ( flxT(k-1,i,j) + flxT(k,i,j) + flxT(k+1,i,j) ) &
                        / ( rho_h(k-1) + rho_h(k) + rho_h(k+1) ) ! running mean
                   if ( fmin < 0.0_RP .and. tmp > fmin ) exit
                   if ( tmp < fmin ) then
                      fmin = tmp
                      kmin = k
                   end if
                end do
             end if
             Zi(i,j) = FZ(kmin,i,j) - FZ(KS-1,i,j)

             ! dens * qv

             ! countergradient flux
             do k = KS, KE_PBL-1
                flx(k) = - ( F2H(k,1,i,j) * lq(k+1) * gammaq(k+1) + F2H(k,2,i,j) * lq(k) * gammaq(k) ) &
                       * rho_h(k)
             end do

             sf_t = SFLX_QV(i,j) / CDZ(KS)
             d(KS) = ( QV(KS,i,j) * DENS(KS,i,j) + dt * ( sf_t - flx(KS) / CDZ(KS) ) ) / RHO(KS)
             do k = KS+1, KE_PBL
                d(k) = QV(k,i,j) - dt * ( flx(k) - flx(k-1) ) / ( CDZ(k) * RHO(k) )
             end do

             call MATRIX_SOLVER_tridiagonal( &
                  KA, KS, KE_PBL, &
                  a(:), b(:), c(:), d(:), & ! (in)
                  phi_n(:)                ) ! (out)

             RHOQV_t(KS,i,j) = ( phi_n(KS) * RHO(KS) - QV(KS,i,j) * DENS(KS,i,j) ) / dt - sf_t
             do k = KS+1, KE_PBL
                RHOQV_t(k,i,j) = ( phi_n(k) - QV(k,i,j) ) * RHO(k) / dt
             end do
             do k = KE_PBL+1, KE
                RHOQV_t(k,i,j) = 0.0_RP
             end do
             flxQ(KS-1,i,j) = 0.0_RP
             do k = KS, KE_PBL-1
                flxQ(k,i,j) = flx(k) &
                            - RHOKh(k) * ( phi_n(k+1) - phi_n(k) ) / FDZ(k)
             end do
             do k = KE_PBL, KE
                flxQ(k,i,j) = 0.0_RP
             end do

          end if


          ! dens * TKE

          ! production
          do k = KS, KE_PBL
             prod(k,i,j) = lq(k) * ( ( sm25(k) + smp(k) ) * dudz2(k,i,j) &
                                   - ( sh25(k) * n2_new(k) - shpgh(k) ) )
          end do
          if ( ATMOS_PHY_BL_MYNN_similarity ) then
             prod(KS,i,j) = us3 / ( KARMAN * z1 ) * ( phi_m - zeta ) &
                          + lq(KS) * ( smp(KS) * dudz2(KS,i,j) + shpgh(KS) )
          end if

          do k = KS, KE_PBL
             tke_p(k) = q(k)**2 * 0.5_RP
             diss(k,i,j) = - 2.0_RP * q(k) / ( B1 * l(k,i,j) )
!             prod(k,i,j) = max( prod(k,i,j), - tke_p(k) / dt - diss(k,i,j) * tke_p(k) )
          end do
          do k = KE_PBL+1, KE
             diss(k,i,j) = 0.0_RP
             prod(k,i,j) = 0.0_RP
          end do

          do k = KS, KE_PBL-1
             d(k) = tke_p(k) * DENS(k,i,j) / RHO(k) + dt * prod(k,i,j)
          end do
          d(KE_PBL) = 0.0_RP

          c(KS) = 0.0_RP
          do k = KS, KE_PBL-1
             ap = - dt * ATMOS_PHY_BL_MYNN_Sq_fact * RHONu(k) / FDZ(k)
             a(k) = ap / ( RHO(k) * CDZ(k) )
             b(k) = - a(k) - c(k) + 1.0_RP - diss(k,i,j) * dt
             c(k+1) = ap / ( RHO(k+1) * CDZ(k+1) )
          end do
          a(KE_PBL) = 0.0_RP
          b(KE_PBL) = - c(KE_PBL) + 1.0_RP - diss(KE_PBL,i,j) * dt

          call MATRIX_SOLVER_tridiagonal( &
               KA, KS, KE_PBL, &
               a(:), b(:), c(:), d(:), & ! (in)
               phi_n(:)                ) ! (out)

          do k = KS, KE_PBL-1
             phi_n(k) = max( phi_n(k), ATMOS_PHY_BL_MYNN_TKE_MIN )
          end do
          phi_n(KE_PBL) = 0.0_RP


          if ( it == nit ) then
             do k = KS, KE_PBL
                diss(k,i,j) = diss(k,i,j) * phi_n(k)
                RPROG_t(k,i,j,I_TKE) = ( phi_n(k) * RHO(k) - PROG(k,i,j,I_TKE) * DENS(k,i,j) ) / dt
             end do
             do k = KE_PBL+1, KE
                RPROG_t(k,i,j,I_TKE) = 0.0_RP
             end do
          else
             do k = KS, KE_PBL
                q(k) = sqrt( phi_n(k) * 2.0_RP )
             end do
          end if


          if ( .not. mynn_level3 ) cycle

          ! dens * tsq

          do k = KS, KE_PBL
             diss_p(k) = dt * 2.0_RP * q(k) / ( B2 * l(k,i,j) )
          end do
          do k = KS, KE_PBL-1
             d(k) = tsq(k) * DENS(k,i,j) / RHO(k) + dt * prod_t(k)
          end do
          d(KE_PBL) = 0.0_RP
          c(KS) = 0.0_RP
          do k = KS, KE_PBL-1
             ap = - dt * RHONu(k) / FDZ(k)
             a(k) = ap / ( RHO(k) * CDZ(k) )
             b(k) = - a(k) - c(k) + 1.0_RP + diss_p(k)
             c(k+1) = ap / ( RHO(k+1) * CDZ(k+1) )
          end do
          a(KE_PBL) = 0.0_RP
          b(KE_PBL) = - c(KE_PBL) + 1.0_RP + diss_p(KE_PBL)

          call MATRIX_SOLVER_tridiagonal( &
               KA, KS, KE_PBL, &
               a(:), b(:), c(:), d(:), & ! (in)
               tsq(:)                  ) ! (out)

          do k = KS, KE_PBL-1
             tsq(k) = max( tsq(k), 0.0_RP )
          end do
          tsq(KE_PBL) = 0.0_RP

          if ( it == nit ) then
             do k = KS, KE_PBL
                RPROG_t(k,i,j,I_TSQ) = ( tsq(k) * RHO(k) - PROG(k,i,j,I_TSQ) * DENS(k,i,j) ) / dt
             end do
             do k = KE_PBL+1, KE
                RPROG_t(k,i,j,I_TSQ) = 0.0_RP
             end do
          end if


          ! dens * qsq

          do k = KS, KE_PBL-1
             d(k) = qsq(k) * DENS(k,i,j) / RHO(k) + dt * prod_q(k)
          end do
          d(KE_PBL) = 0.0_RP
          ! a, b, c are same as those for tsq

          call MATRIX_SOLVER_tridiagonal( &
               KA, KS, KE_PBL, &
               a(:), b(:), c(:), d(:), & ! (in)
               qsq(:)                  ) ! (out)

          do k = KS, KE_PBL-1
             qsq(k) = max( qsq(k), 0.0_RP )
          end do
          qsq(KE_PBL) = 0.0_RP

          if ( it == nit ) then
             do k = KS, KE_PBL
                RPROG_t(k,i,j,I_QSQ) = ( qsq(k) * RHO(k) - PROG(k,i,j,I_QSQ) * DENS(k,i,j) ) / dt
             end do
             do k = KE_PBL+1, KE
                RPROG_t(k,i,j,I_QSQ) = 0.0_RP
             end do
          end if


          ! dens * cov

          do k = KS, KE_PBL-1
             d(k) = cov(k) * DENS(k,i,j) / RHO(k) + dt * prod_c(k)
          end do
          d(KE_PBL) = 0.0_RP
          ! a, b, c are same as those for tsq

          call MATRIX_SOLVER_tridiagonal( &
               KA, KS, KE_PBL, &
               a(:), b(:), c(:), d(:), & ! (in)
               cov(:)                  ) ! (out)

          do k = KS, KE_PBL-1
             cov(k) = sign( min( abs(cov(k)), sqrt(tsq(k)*qsq(k)) ), cov(k) )
          end do
          cov(KE_PBL) = 0.0_RP

          if ( it == nit ) then
             do k = KS, KE_PBL
                RPROG_t(k,i,j,I_COV) = ( cov(k) * RHO(k) - PROG(k,i,j,I_COV) * DENS(k,i,j) ) / dt
             end do
             do k = KE_PBL+1, KE
                RPROG_t(k,i,j,I_COV) = 0.0_RP
             end do
          end if

       end do

       Nu(KS-1,i,j) = 0.0_RP
       Kh(KS-1,i,j) = 0.0_RP
       do k = KE_PBL, KE
          Nu   (k,i,j) = 0.0_RP
          Kh   (k,i,j) = 0.0_RP
          Pr   (k,i,j) = 1.0_RP
       end do
       do k = KE_PBL+1, KE
          Ri     (k,i,j) = UNDEF
          prod   (k,i,j) = UNDEF
          diss   (k,i,j) = UNDEF
          dudz2  (k,i,j) = UNDEF
          l      (k,i,j) = UNDEF
          Qlp    (k,i,j) = UNDEF
          cldfrac(k,i,j) = UNDEF
       end do

    end do
    end do


    call FILE_HISTORY_in(Ri(:,:,:), 'Ri_MYNN', 'Richardson number', '1',     fill_halo=.true. )
    call FILE_HISTORY_in(Pr(:,:,:), 'Pr_MYNN', 'Prandtl number',    '1',     fill_halo=.true., dim_type="ZHXY" )
    call FILE_HISTORY_in(prod(:,:,:), 'TKE_prod_MYNN', 'TKE production',  'm2/s3', fill_halo=.true.)
    call FILE_HISTORY_in(diss(:,:,:), 'TKE_diss_MYNN', 'TKE dissipation', 'm2/s3', fill_halo=.true.)
    call FILE_HISTORY_in(dudz2(:,:,:), 'dUdZ2_MYNN', 'dudz2', 'm2/s2', fill_halo=.true.)
    call FILE_HISTORY_in(l(:,:,:), 'L_mix_MYNN', 'minxing length', 'm', fill_halo=.true.)

    call FILE_HISTORY_in(flxU(:,:,:), 'ZFLX_RHOU_BL', 'Z FLUX of RHOU (MYNN)', 'kg/m/s2', fill_halo=.true., dim_type="ZHXY" )
    call FILE_HISTORY_in(flxV(:,:,:), 'ZFLX_RHOV_BL', 'Z FLUX of RHOV (MYNN)', 'kg/m/s2', fill_halo=.true., dim_type="ZHXY" )
    call FILE_HISTORY_in(flxT(:,:,:), 'ZFLX_RHOT_BL', 'Z FLUX of RHOT (MYNN)', 'K kg/m2/s', fill_halo=.true., dim_type="ZHXY" )
    call FILE_HISTORY_in(flxQ(:,:,:), 'ZFLX_QV_BL',   'Z FLUX of RHOQV (MYNN)', 'kg/m2/s', fill_halo=.true., dim_type="ZHXY" )


    initialize = .false.

    return
  end subroutine ATMOS_PHY_BL_MYNN_tendency

  !-----------------------------------------------------------------------------
  ! private routines
  !-----------------------------------------------------------------------------

!OCL SERIAL
  subroutine get_length( &
       KA, KS, KE_PBL, &
       i, j,           &
       q, n2,           &
       SFLX_BUOY, RLmo, &
       CZ, FZ,          &
       l                )
    use scale_const, only: &
       GRAV   => CONST_GRAV, &
       KARMAN => CONST_KARMAN, &
       EPS    => CONST_EPS
    implicit none
    integer,  intent(in) :: KA, KS, KE_PBL
    integer,  intent(in) :: i, j ! for debug

    real(RP), intent(in) :: q(KA)
    real(RP), intent(in) :: n2(KA)
    real(RP), intent(in) :: SFLX_BUOY !> g/T0 <w Tv> @ surface
    real(RP), intent(in) :: RLmo      !> inverse of Obukhov length
    real(RP), intent(in) :: CZ(KA)
    real(RP), intent(in) :: FZ(0:KA)

    real(RP), intent(out) :: l(KA)

    real(RP), parameter :: ls_fact_max = 2.0_RP

    real(RP) :: ls     !> L_S
    real(RP) :: lb     !> L_B
    real(RP) :: lt     !> L_T
    real(RP) :: rlt    !> 1/L_T

    real(RP) :: qc     !> q_c
    real(RP) :: int_q  !> \int q dz
    real(RP) :: int_qz !> \int qz dz
    real(RP) :: rn2sr  !> 1/N
    real(RP) :: zeta   !> height normalized by the Obukhov length

    real(RP) :: z
    real(RP) :: qdz

    real(RP) :: sw
    integer :: k

    int_qz = 0.0_RP
    int_q = 0.0_RP
    do k = KS, KE_PBL
       qdz = q(k) * ( FZ(k) - FZ(k-1) )
       z = CZ(k) - FZ(KS-1)
       int_qz = int_qz + z * qdz
       int_q  = int_q + qdz
    end do
    ! LT
    lt = min( max(0.23_RP * int_qz / (int_q + 1e-20_RP), &
                  LT_min), &
                  ATMOS_PHY_BL_MYNN_Lt_MAX )
    rlt = 1.0_RP / lt

    qc = ( lt * max(SFLX_BUOY,0.0_RP) )**OneOverThree ! qc=0 if SFLX_BUOY<0

    do k = KS, KE_PBL
       z = CZ(k) - FZ(KS-1)
       zeta = min( max( z * RLmo, zeta_min ), zeta_max )

       ! LS
       sw = sign(0.5_RP, zeta) + 0.5_RP ! 1 for zeta >= 0, 0 for zeta < 0
       ls = KARMAN * z &
          * ( sw / (1.0_RP + 2.7_RP*zeta*sw ) &
            + min( ( (1.0_RP - 100.0_RP*zeta)*(1.0_RP-sw) )**0.2_RP, ls_fact_max) )

       ! LB
       sw  = sign(0.5_RP, n2(k)-EPS) + 0.5_RP ! 1 for dptdz >0, 0 for dptdz <= 0
       rn2sr = 1.0_RP / ( sqrt(n2(k)*sw) + 1.0_RP-sw)
       lb = (1.0_RP + 5.0_RP * sqrt(qc*rn2sr*rlt)) * q(k) * rn2sr * sw & ! qc=0 when RLmo > 0
           +  1.E10_RP * (1.0_RP-sw)

       ! L
       l(k) = 1.0_RP / ( 1.0_RP/ls + rlt + 1.0_RP/(lb+1E-20_RP) )
    end do

    return
  end subroutine get_length

!OCL SERIAL
  subroutine get_q2_level2( &
       KA, KS, KE_PBL, &
       dudz2, Ri, l, &
       q2_2          )
    implicit none
    integer,  intent(in)  :: KA, KS, KE_PBL

    real(RP), intent(in)  :: dudz2(KA)
    real(RP), intent(in)  :: Ri(KA)
    real(RP), intent(in)  :: l(KA)

    real(RP), intent(out) :: q2_2(KA)

    real(RP) :: rf   !> Rf
    real(RP) :: sm_2 !> sm for level 2
    real(RP) :: sh_2 !> sh for level 2

    integer :: k

    do k = KS, KE_PBL
       rf = min(0.5_RP / AF12 * ( Ri(k) &
                              + AF12*Rf1 &
                              - sqrt(Ri(k)**2 + 2.0_RP*AF12*(Rf1-2.0_RP*Rf2)*Ri(k) + (AF12*Rf1)**2) ), &
                Rfc)
       sh_2 = 3.0_RP * A2 * (G1+G2) * (Rfc-rf) / (1.0_RP-rf)
       sm_2 = sh_2 * AF12 * (Rf1-rf) / (Rf2-rf)
       q2_2(k) = B1 * l(k)**2 * sm_2 * (1.0_RP-rf) * dudz2(k)
    end do

    return
  end subroutine get_q2_level2

!OCL SERIAL
  subroutine get_smsh( &
       KA, KS, KE_PBL, &
       i, j,            &
       q, ac,           &
       l, n2,           &
       potv, dudz2,     &
       dtldz, dqwdz,    &
       betat, betaq,    &
       mynn_level3,     &
       initialize,      &
       tsq, qsq, cov,   &
       sm25, f_smp,     &
       sh25, f_shpgh,   &
       f_gamma,         &
       tltv25, qwtv25,  &
       tvsq25,          &
       tvsq_up, tvsq_lo )
    use scale_const, only: &
       EPS  => CONST_EPS, &
       HUGE => CONST_HUGE, &
       GRAV => CONST_GRAV
    implicit none
    integer,  intent(in)  :: KA, KS, KE_PBL
    integer,  intent(in)  :: i, j ! for debug

    real(RP), intent(in)  :: q(KA)
    real(RP), intent(in)  :: ac(KA)
    real(RP), intent(in)  :: l(KA)
    real(RP), intent(in)  :: n2(KA)
    real(RP), intent(in)  :: potv(KA)
    real(RP), intent(in)  :: dudz2(KA)
    real(RP), intent(in)  :: dtldz(KA)
    real(RP), intent(in)  :: dqwdz(KA)
    real(RP), intent(in)  :: betat(KA)
    real(RP), intent(in)  :: betaq(KA)
    logical,  intent(in)  :: mynn_level3
    logical,  intent(in)  :: initialize

    real(RP), intent(inout)  :: tsq(KA)
    real(RP), intent(inout)  :: qsq(KA)
    real(RP), intent(inout)  :: cov(KA)

    real(RP), intent(out) :: sm25   (KA) ! S_M2.5
    real(RP), intent(out) :: f_smp  (KA) ! S_M' / ( <\theta_v^2> - <\theta_v^2>2.5 )
    real(RP), intent(out) :: sh25   (KA) ! S_H2.5
    real(RP), intent(out) :: f_shpgh(KA) ! S_H' G_H * (q/L)^2 / ( <\theta_v^2> - <\theta_v^2>2.5 )
    real(RP), intent(out) :: f_gamma(KA) ! - E_H / q^2 * GRAV / Theta_0
    real(RP), intent(out) :: tltv25 (KA) !> <\theta_l \theta_v> for level 2.5
    real(RP), intent(out) :: qwtv25 (KA) !> <q_w \theta_v> for level 2.5
    real(RP), intent(out) :: tvsq25 (KA) !> <\theta_v^2> for level 2.5
    real(RP), intent(out) :: tvsq_up(KA)
    real(RP), intent(out) :: tvsq_lo(KA)

    real(RP) :: l2     !> L^2
    real(RP) :: q2     !> q^2
    real(RP) :: l2q2   !> l^2 / q^2
    real(RP) :: ac2    !> \alpha_c^2
    real(RP) :: p1q2   !> \Phi_1 * q^2
    real(RP) :: p2q2   !> \Phi_2 * q^2
    real(RP) :: p3q2   !> \Phi_3 * q^2
    real(RP) :: p4q2   !> \Phi_4 * q^2
    real(RP) :: p5q2   !> \Phi_5 * q^2
    real(RP) :: rd25q2 !> q2 / D_2.5
    real(RP) :: ghq2   !> G_H * q^2
    real(RP) :: gmq2   !> G_M * q^2
    real(RP) :: f1, f2, f3, f4

    ! for level 3
    real(RP) :: tvsq   !> <\theta_v^2>
    real(RP) :: tltv
    real(RP) :: qwtv
    real(RP) :: tsq25
    real(RP) :: qsq25
    real(RP) :: cov25
    real(RP) :: emq2   !> E_M * G_H / q^2
    real(RP) :: eh     !> E_H
    real(RP) :: ew     !> E_w
    real(RP) :: rdpq2  !> q2 / D'
    real(RP) :: cw25   !> Cw2.5
    real(RP) :: fact
    real(RP) :: tmp1, tmp2

    integer :: k

    ! level 2.5
    do k = KS, KE_PBL

       ac2 = ac(k)**2
       l2 = l(k)**2
       q2 = q(k)**2

       f1 = -  3.0_RP * ac2 * A2 * B2 * ( 1.0_RP - C3 )
       f2 = -  9.0_RP * ac2 * A1 * A2 * ( 1.0_RP - C2 )
       f3 =    9.0_RP * ac2 * A2**2   * ( 1.0_RP - C2 ) * ( 1.0_RP - C5 )
       f4 = - 12.0_RP * ac2 * A1 * A2 * ( 1.0_RP - C2 )

       if ( mynn_level3 ) then ! level 3
          ghq2 = max( - n2(k) * l2, -q2 ) ! L/q <= 1/N for N^2>0
       else
          ghq2 = - n2(k) * l2
       end if
       gmq2 = dudz2(k) * l2

       p1q2   = q2 + f1 * ghq2
       p2q2   = q2 + f2 * ghq2
       p3q2   = q2 + ( f1 + f3 ) * ghq2
       p4q2   = q2 + ( f1 + f4 ) * ghq2
       p5q2   = 6.0_RP * ac2 * A1**2 * gmq2
       rd25q2 = q2 / max( p2q2 * p4q2 + p5q2 * p3q2, 1e-20_RP )

       sm25(k) = ac(k) * A1 * ( p3q2 - 3.0_RP * C1 * p4q2 ) * rd25q2
       sh25(k) = ac(k) * A2 * ( p2q2 + 3.0_RP * C1 * p5q2 ) * rd25q2

       fact = ac(k) * B2 * l2 * sh25(k)
       tsq25 = fact * dtldz(k)**2
       qsq25 = fact * dqwdz(k)**2
       cov25 = fact * dtldz(k) * dqwdz(k)

       if ( initialize .or. (.not. mynn_level3 ) ) then
          tsq(k) = tsq25
          qsq(k) = qsq25
          cov(k) = cov25
       end if

       if ( mynn_level3 ) then ! level 3

          if ( q2 <= 1e-10_RP ) then
             tltv25 (k) = 0.0_RP
             qwtv25 (k) = 0.0_RP
             tvsq25 (k) = 0.0_RP
             tvsq_up(k) = 0.0_RP
             tvsq_lo(k) = 0.0_RP
             f_smp  (k) = 0.0_RP
             f_shpgh(k) = 0.0_RP
             f_gamma(k) = 0.0_RP
          else
             tltv25(k) = betat(k) * tsq25 + betaq(k) * cov25
             qwtv25(k) = betat(k) * cov25 + betaq(k) * qsq25
             tvsq25(k) = max(betat(k) * tltv25(k) + betaq(k) * qwtv25(k), 0.0_RP)
             cw25 = p1q2 * ( p2q2 + 3.0_RP * C1 * p5q2 ) * rd25q2 / ( 3.0_RP * q2 )

             rdpq2  = q2 / max( p2q2 * ( f4 * ghq2 + q2 ) + p5q2 * ( f3 * ghq2 + q2 ), 1e-20_RP )

             l2q2 = l2 / q2
!             l2q2 = min( l2q2, 1.0_RP / max(n2(k), EPS) )

             tltv = betat(k) * tsq(k) + betaq(k) * cov(k)
             qwtv = betat(k) * cov(k) + betaq(k) * qsq(k)
             tvsq = max(betat(k) * tltv + betaq(k) * qwtv, 0.0_RP)

             ew = ( 1.0_RP - C3 ) * ( - p2q2 * f4 - p5q2 * f3 ) * rdpq2 ! (p1-p4)/gh = -f4, (p1-p3)/gh = -f3
             if ( abs(ew) > EPS ) then
                fact = q2 * POTV(k)**2 / ( ew * l2q2 * GRAV**2 )
                if ( fact > 0.0_RP ) then
                   tmp1 = 0.76_RP
                   tmp2 = 0.12_RP
                else
                   tmp1 = 0.12_RP
                   tmp2 = 0.76_RP
                end if
                tvsq_up(k) = ( tmp1 - cw25 ) * fact
                tvsq_lo(k) = ( tmp2 - cw25 ) * fact
             else
                tvsq_up(k) =  HUGE
                tvsq_lo(k) = -HUGE
             end if

             emq2 = 3.0_RP * ac(k) * A1 * ( 1.0_RP - C3 ) * ( f3 - f4 ) * rdpq2 ! (p3-p4)/gh = (f3-f4)
             eh   = 3.0_RP * ac(k) * A2 * ( 1.0_RP - C3 ) * ( p2q2 + p5q2 ) * rdpq2

!             q2 = l2 / l2q2
             fact = GRAV / POTV(k)
             f_smp  (k) = emq2 * fact**2 * l2q2
             f_shpgh(k) = eh   * fact**2 / q2
             f_gamma(k) = - eh * fact / q2
          end if
       else ! level 2.5
          tvsq_up(k) = 0.0_RP
          tvsq_lo(k) = 0.0_RP
          f_smp  (k) = 0.0_RP
          f_shpgh(k) = 0.0_RP
          f_gamma(k) = 0.0_RP
       end if

    end do

    return
  end subroutine get_smsh

!OCL SERIAL
  subroutine get_gamma_implicit( &
       KA, KS, KE, &
       i, j,       &
       tsq, qsq, cov,          &
       dtldz, dqwdz, POTV,     &
       prod_t, prod_q, prod_c, &
       betat, betaq,           &
       f_gamma, l, q,          &
       dt,                     &
       dtsq, dqsq, dcov        )
    use scale_const, only: &
       GRAV => CONST_GRAV
    implicit none
    integer, intent(in) :: KA, KS, KE
    integer, intent(in) :: i, j ! for debug

    real(RP), intent(in) :: tsq    (KA)
    real(RP), intent(in) :: qsq    (KA)
    real(RP), intent(in) :: cov    (KA)
    real(RP), intent(in) :: dtldz  (KA)
    real(RP), intent(in) :: dqwdz  (KA)
    real(RP), intent(in) :: POTV   (KA)
    real(RP), intent(in) :: prod_t (KA)
    real(RP), intent(in) :: prod_q (KA)
    real(RP), intent(in) :: prod_c (KA)
    real(RP), intent(in) :: betat  (KA)
    real(RP), intent(in) :: betaq  (KA)
    real(RP), intent(in) :: f_gamma(KA)
    real(RP), intent(in) :: l      (KA)
    real(RP), intent(in) :: q      (KA)
    real(RP), intent(in) :: dt

    real(RP), intent(out) :: dtsq(KA)
    real(RP), intent(out) :: dqsq(KA)
    real(RP), intent(out) :: dcov(KA)

    real(RP) :: a11, a12, a21, a22, a23, a32, a33
    real(RP) :: v1, v2, v3
    real(RP) :: f1, f2
    real(RP) :: a13, a31, det

    integer :: k

    ! calculate gamma by implicit scheme

    do k = KS, KE

       ! matrix coefficient
       f1 = f_gamma(k) * l(k) * q(k)
       f2 = - 2.0_RP * q(k) / ( B2 * l(k) )

       v1 = prod_t(k) + f2 * tsq(k)
       v2 = prod_c(k) + f2 * cov(k)
       v3 = prod_q(k) + f2 * qsq(k)

       a11 = - f1 * dtldz(k) * betat(k) * 2.0_RP + 1.0_RP / dt - f2
       a12 = - f1 * dtldz(k) * betaq(k) * 2.0_RP
       a21 = - f1 * dqwdz(k) * betat(k)
       a22 = - f1 * ( dtldz(k) * betat(k) + dqwdz(k) * betaq(k) ) + 1.0_RP / dt - f2
       a23 = - f1 * dtldz(k) * betaq(k)
       a32 = - f1 * dqwdz(k) * betat(k) * 2.0_RP
       a33 = - f1 * dqwdz(k) * betaq(k) * 2.0_RP + 1.0_RP / dt - f2

       if ( q(k) < 1e-20_RP .or. abs(f_gamma(k)) * dt >= q(k) ) then
          ! solve matrix
          f1 = a21 / a11
          f2 = a23 / a33
          dcov(k) = ( v2 - f1 * v1 - f2 * v3 ) / ( a22 - f1 * a12 - f2 * a32 )
          dtsq(k) = ( v1 - a12 * dcov(k) ) / a11
          dqsq(k) = ( v3 - a32 * dcov(k) ) / a33
       else
          ! consider change in q
          ! dq = - L * G / Theta0 * ( betat * dgammat + betaq * dgammaq )
          f1 = - l(k) * GRAV / POTV(k) * f_gamma(k) / q(k)
          f2 = f1 * betat(k)**2
          a11 = a11 - f2 * v1
          a21 = a21 - f2 * v2
          a31 =     - f2 * v3
          f2 = f1 * betat(k) * betaq(k) * 2.0_RP
          a12 = a12 - f2 * v1
          a22 = a22 - f2 * v2
          a32 = a32 - f2 * v3
          f2 = f1 * betaq(k)**2
          a13 =     - f2 * v1
          a23 = a23 - f2 * v2
          a33 = a33 - f2 * v3
          ! solve matrix by Cramer's fomula
          det = a11 * ( a22 * a33 - a23 * a32 ) + a12 * ( a23 * a31 - a21 * a33 ) + a13 * ( a21 * a32 - a22 * a31 )
          dtsq(k) = ( v1 * ( a22 * a33 - a23 * a32 ) + a12 * ( a23 * v3 - v2 * a33 ) + a13 * ( v2 * a32 - a22 * v3 ) ) / det
          dcov(k) = ( a11 * ( v2 * a33 - a23 * v3 ) + v1 * ( a23 * a31 - a21 * a33 ) + a13 * ( a21 * v3 - v2 * a31 ) ) / det
          dqsq(k) = ( a11 * ( a22 * v3 - v2 * a32 ) + a12 * ( v2 * a31 - a21 * v3 ) + v3 * ( a21 * a32 - a22 * a31 ) ) / det
       end if

    end do

    return
  end subroutine get_gamma_implicit

!OCL SERIAL
  subroutine calc_vertical_differece( &
       KA, KS, KE, &
       i, j,       &
       U, V, POTL, Qw, QDRY, &
       CDZ, FDZ, F2H,        &
       dudz2, dtldz, dqwdz   )
    use scale_const, only: &
       EPS => CONST_EPS
    integer,  intent(in)  :: KA, KS, KE
    integer,  intent(in)  :: i, j  ! for debug

    real(RP), intent(in) :: U   (KA)
    real(RP), intent(in) :: V   (KA)
    real(RP), intent(in) :: POTL(KA)
    real(RP), intent(in) :: Qw  (KA)
    real(RP), intent(in) :: QDRY(KA)
    real(RP), intent(in) :: CDZ (KA)
    real(RP), intent(in) :: FDZ (KA)
    real(RP), intent(in) :: F2H (KA,2)

    real(RP), intent(out) :: dudz2(KA)
    real(RP), intent(out) :: dtldz(KA)
    real(RP), intent(out) :: dqwdz(KA)

    real(RP) :: Uh(KA), Vh(KA)
    real(RP) :: qw2(KA)
    real(RP) :: qh(KA)

    integer :: k

    do k = KS, KE
       Uh(k) = f2h(k,1) * U(k+1) + f2h(k,2) * U(k)
       Vh(k) = f2h(k,1) * V(k+1) + f2h(k,2) * V(k)
    end do

    dudz2(KS) = ( ( Uh(KS) - U(KS) )**2 + ( Vh(KS) - V(KS) )**2 ) / ( CDZ(KS) * 0.5_RP )**2
!    dudz2(KS) = ( ( Uh(KS) )**2 + ( Vh(KS) )**2 ) / CDZ(KS)**2
    dudz2(KS) = max( dudz2(KS), 1e-20_RP )
    do k = KS+1, KE
       dudz2(k) = ( ( Uh(k) - Uh(k-1) )**2 + ( Vh(k) - Vh(k-1) )**2 ) / CDZ(k)**2
!       dudz2(k) = ( &
!              ( U(k+1) * FDZ(k-1)**2 - U(k-1) * FDZ(k)**2 + U(k) * ( FDZ(k)**2 - FDZ(k-1)**2 ) )**2 &
!            + ( V(k+1) * FDZ(k-1)**2 - V(k-1) * FDZ(k)**2 + V(k) * ( FDZ(k)**2 - FDZ(k-1)**2 ) )**2 &
!            ) / ( FDZ(k-1) * FDZ(k) * ( FDZ(k-1) + FDZ(k) ) )**2
       dudz2(k) = max( dudz2(k), 1e-20_RP )
    end do


    do k = KS, KE
       qh(k) = f2h(k,1) * POTL(k+1) + f2h(k,2) * POTL(k)
    end do
    dtldz(KS) = ( qh(KS) - POTL(KS) ) / ( CDZ(KS) * 0.5_RP )
!    dtldz(KS) = ( POTL(KS+1) - POTL(KS) ) / FDZ(KS)
    do k = KS+1, KE
       dtldz(k) = ( qh(k) - qh(k-1) ) / CDZ(k)
!       dtldz(k) = ( POTL(k+1) * FDZ(k-1)**2 - POTL(k-1) * FDZ(k)**2 + POTL(k) * ( FDZ(k)**2 - FDZ(k-1)**2 ) ) &
!            / ( FDZ(k-1) * FDZ(k) * ( FDZ(k-1) + FDZ(k) ) )
    end do

    do k = KS, KE+1
       qw2(k) = Qw(k) / ( QDRY(k) + Qw(k) )
    end do
    do k = KS, KE
!       qh(k) = f2h(k,1) * Qw(k+1) + f2h(k,2) * Qw(k)
       qh(k) = f2h(k,1) * qw2(k+1) + f2h(k,2) * qw2(k)
    end do
    dqwdz(KS) = ( qh(KS) - qw2(KS) ) / ( CDZ(KS) * 0.5_RP )
!    dqwdz(KS) = ( qh(KS) - Qw(KS) ) / ( CDZ(KS) * 0.5_RP )
!    dqwdz(KS) = ( Qw(KS+1) - Qw(KS) ) / FDZ(KS)
    do k = KS+1, KE
       dqwdz(k) = ( qh(k) - qh(k-1) ) / CDZ(k)
!       dqwdz(k) = ( Qw(k+1) * FDZ(k-1)**2 - Qw(k-1) * FDZ(k)**2 + Qw(k) * ( FDZ(k)**2 - FDZ(k-1)**2 ) ) &
!            / ( FDZ(k-1) * FDZ(k) * ( FDZ(k-1) + FDZ(k) ) )
    end do

    return
  end subroutine calc_vertical_differece

!OCL SERIAL
  subroutine partial_condensation( &
       KA, KS, KE, &
       PRES, POTT, POTL, Qw, Qdry, EXNER, &
       tsq, qsq, cov,                     &
       betat, betaq, Qlp, cldfrac         )
    use scale_const, only: &
       CPdry   => CONST_CPdry,  &
       Rvap    => CONST_Rvap,   &
       EPSvap  => CONST_EPSvap, &
       EPSTvap => CONST_EPSTvap
    use scale_atmos_saturation, only: &
       ATMOS_SATURATION_psat => ATMOS_SATURATION_psat_liq
!       ATMOS_SATURATION_pres2qsat => ATMOS_SATURATION_pres2qsat_liq
    use scale_atmos_hydrometeor, only: &
       HYDROMETEOR_LHV => ATMOS_HYDROMETEOR_LHV, &
       LHV,      &
       CP_VAPOR, &
       CP_WATER
    integer,  intent(in)  :: KA, KS, KE
    real(RP), intent(in)  :: PRES (KA)
    real(RP), intent(in)  :: POTT (KA)
    real(RP), intent(in)  :: POTL (KA)
    real(RP), intent(in)  :: Qw   (KA)
    real(RP), intent(in)  :: Qdry (KA)
    real(RP), intent(in)  :: EXNER(KA)
    real(RP), intent(in)  :: tsq  (KA)
    real(RP), intent(in)  :: qsq  (KA)
    real(RP), intent(in)  :: cov  (KA)

    real(RP), intent(out) :: betat  (KA)
    real(RP), intent(out) :: betaq  (KA)
    real(RP), intent(out) :: Qlp    (KA)
    real(RP), intent(out) :: cldfrac(KA)

    real(RP) :: TEML(KA)
    real(RP) :: LHVL(KA)
    real(RP) :: psat(KA)
    real(RP) :: Qsl
    real(RP) :: dQsl
    real(RP) :: CPtot
    real(RP) :: aa, bb, cc
    real(RP) :: sigma_s
    real(RP) :: Q1
    real(RP) :: Rt

    integer :: k

    do k = KS, KE
       TEML(k) = POTL(k) * EXNER(k)
    end do

!!$    call ATMOS_SATURATION_pres2qsat( &
!!$         KA, KS, KE, &
!!$         TEML(:), PRES(:), Qdry(:), & ! (in)
!!$         Qsl(:)                     ) ! (out)
    call ATMOS_SATURATION_psat( &
         KA, KS, KE, &
         TEML(:), & ! (in)
         psat(:)  ) ! (out)

    call HYDROMETEOR_LHV( &
         KA, KS, KE, & ! (in)
         TEML(:), & ! (in)
         LHVL(:)  ) ! (out)

    do k = KS, KE

       Qsl = EPSvap * psat(k) / ( PRES(k) - ( 1.0_RP - EPSvap ) * psat(k) )
       dQsl = PRES(k) * Qsl**2 * LHVL(k) / ( EPSvap * psat(k) * Rvap * TEML(k)**2 )

       CPtot = ( 1.0_RP - Qsl ) * CPdry + Qsl * CP_VAPOR
       aa = 1.0_RP / ( 1.0_RP + LHVL(k)/CPtot * dQsl )
       bb = EXNER(k) * dQsl

       sigma_s = min( max( &
            0.5_RP * aa * sqrt( max( qsq(k) - 2.0_RP * bb * cov(k) + bb**2 * tsq(k), 1.0e-20_RP ) ), &
            aa * Qsl * 0.09_RP), aa * Qsl )
       Q1 = aa * ( Qw(k) - Qsl ) * 0.5_RP / sigma_s
       cldfrac(k) = min( max( 0.5_RP * ( 1.0_RP + erf(Q1*rsqrt_2) ), 0.0_RP ), 1.0_RP )
       Qlp(k) = min( max( 2.0_RP * sigma_s * ( cldfrac(k) * Q1 + rsqrt_2pi &
#if defined(NVIDIA) || defined(SX)
               * exp( -min( 0.5_RP*Q1**2, 1.E+3_RP ) ) & ! apply exp limiter
#else
               * exp(-0.5_RP*Q1**2) &
#endif
               ), 0.0_RP ), Qw(k) * 0.5_RP )
       cc = ( 1.0_RP + EPSTvap * Qw(k) - Qlp(k) / EPSvap ) / EXNER(k) * LHVL(k) / CPtot &
             - POTT(k) / EPSvap
       Rt = cldfrac(k) - Qlp(k) / (2.0_RP*sigma_s*sqrt_2pi) &
#if defined(NVIDIA) || defined(SX)
               * exp( -min( Q1**2 * 0.5_RP, 1.E+3_RP ) ) ! apply exp limiter
#else
               * exp(-Q1**2 * 0.5_RP)
#endif
       betat(k) = 1.0_RP + EPSTvap * Qw(k) - Qlp(k) / EPSvap - Rt * aa * bb * cc
       betaq(k) = EPSTvap * POTT(k) + Rt * aa * cc

    end do

    return
  end subroutine partial_condensation

end module scale_atmos_phy_bl_mynn
