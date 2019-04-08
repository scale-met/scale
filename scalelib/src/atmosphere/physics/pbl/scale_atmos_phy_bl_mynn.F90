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
  public :: ATMOS_PHY_BL_mynn_tendency
  public :: ATMOS_PHY_BL_mynn_tendency_tracer

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

  real(RP), private            :: SQRT_2PI
  real(RP), private            :: RSQRT_2PI
  real(RP), private            :: RSQRT_2

  integer,  private            :: KE_PBL

  logical,  private            :: initialize

  real(RP), private            :: ATMOS_PHY_BL_MYNN_PBL_MAX  = 1.E+10_RP !> maximum height of the PBL
  real(RP), private            :: ATMOS_PHY_BL_MYNN_TKE_MIN   =  1.E-20_RP
  real(RP), private            :: ATMOS_PHY_BL_MYNN_N2_MAX    =  1.E-2_RP
  real(RP), private            :: ATMOS_PHY_BL_MYNN_NU_MIN    = -1.E1_RP
  real(RP), private            :: ATMOS_PHY_BL_MYNN_NU_MAX    =  1.E4_RP
  real(RP), private            :: ATMOS_PHY_BL_MYNN_KH_MIN    = -1.E1_RP
  real(RP), private            :: ATMOS_PHY_BL_MYNN_KH_MAX    =  1.E4_RP
  real(RP), private            :: ATMOS_PHY_BL_MYNN_Lt_MAX    =   700.0_RP ! ~ 0.23 * 3 km
  real(RP), private            :: ATMOS_PHY_BL_MYNN_Sq_fact   = 3.0_RP
  logical,  private            :: ATMOS_PHY_BL_MYNN_init_TKE  = .false.
  logical,  private            :: ATMOS_PHY_BL_MYNN_prod_surf = .false.

  character(len=H_SHORT), private  :: ATMOS_PHY_BL_MYNN_LEVEL = "2.5" ! "2.5" or "3", level 3 is under experimental yet.

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
       ATMOS_PHY_BL_MYNN_prod_surf


  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> ATMOS_PHY_BL_MYNN_tracer_setup
  !! Tracer Setup
  !<
  subroutine ATMOS_PHY_BL_MYNN_tracer_setup( )
    use scale_prc, only: &
       PRC_abort

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
       CZ, &
       TKE_MIN, PBL_MAX )
    use scale_prc, only: &
       PRC_abort
    use scale_const, only: &
       PI => CONST_PI
    implicit none

    integer,  intent(in) :: KA, KS, KE
    integer,  intent(in) :: IA, IS, IE
    integer,  intent(in) :: JA, JS, JE
    real(RP), intent(in) :: CZ (KA,IA,JA)

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

    do k = KS, KE-1
       do j = JS, JE
       do i = IS, IE
          if ( ATMOS_PHY_BL_MYNN_PBL_MAX >= CZ(k,i,j) ) then
             KE_PBL = k
          end if
       end do
       end do
    end do

    if ( ATMOS_PHY_BL_MYNN_LEVEL == "3" ) then
       LOG_WARN("ATMOS_PHY_BL_MYNN_setup", *) "At this moment, level 3 is still experimental"
    end if

    initialize = ATMOS_PHY_BL_MYNN_init_TKE

    return
  end subroutine ATMOS_PHY_BL_MYNN_setup

  !-----------------------------------------------------------------------------
  !> ATMOS_PHY_BL_MYNN_tendency
  !! calculate tendency by the virtical eddy viscosity
  !<
  subroutine ATMOS_PHY_BL_MYNN_tendency( &
       KA, KS, KE, IA, IS, IE, JA, JS, JE, &
       DENS, U, V, POTT, PROG,             &
       PRES, EXNER, N2,                    &
       QDRY, QV, Qw, POTL, POTV,           &
       SFLX_MU, SFLX_MV, SFLX_SH, SFLX_QV, &
       l_mo,                               &
       CZ, FZ, dt_DP,                      &
       RHOU_t, RHOV_t, RHOT_t, RPROG_t,    &
       Nu, Nu_cg, Kh, Kh_cg                )
    use scale_const, only: &
       EPS     => CONST_EPS,    &
       GRAV    => CONST_GRAV,   &
       KARMAN  => CONST_KARMAN, &
       CPdry   => CONST_CPdry,  &
       EPSTvap => CONST_EPSTvap
    use scale_prc, only: &
       PRC_abort
    use scale_atmos_hydrometeor, only: &
       HYDROMETEOR_LHV => ATMOS_HYDROMETEOR_LHV
    use scale_matrix, only: &
       MATRIX_SOLVER_tridiagonal
    use scale_file_history, only: &
       FILE_HISTORY_in
    implicit none

    integer, intent(in) :: KA, KS, KE
    integer, intent(in) :: IA, IS, IE
    integer, intent(in) :: JA, JS, JE

    real(RP), intent(in) :: DENS   (KA,IA,JA) !> density
    real(RP), intent(in) :: U      (KA,IA,JA) !> zonal wind
    real(RP), intent(in) :: V      (KA,IA,JA) !> meridional wind
    real(RP), intent(in) :: POTT   (KA,IA,JA) !> potential temperature
    real(RP), intent(in) :: PROG   (KA,IA,JA,ATMOS_PHY_BL_MYNN_ntracer) !> prognostic variables (TKE, TSQ, QSQ, COV)
    real(RP), intent(in) :: PRES   (KA,IA,JA) !> pressure
    real(RP), intent(in) :: EXNER  (KA,IA,JA) !> Exner function
    real(RP), intent(in) :: N2     (KA,IA,JA) !> squared Brunt-Vaisala frequency
    real(RP), intent(in) :: QDRY   (KA,IA,JA) !> dry air
    real(RP), intent(in) :: QV     (KA,IA,JA) !> vapor
    real(RP), intent(in) :: Qw     (KA,IA,JA) !> total water content
    real(RP), intent(in) :: POTL   (KA,IA,JA) !> liquid water potential temp.
    real(RP), intent(in) :: POTV   (KA,IA,JA) !> virtual potential temp.
    real(RP), intent(in) :: SFLX_MU(   IA,JA) !> surface flux of zonal wind
    real(RP), intent(in) :: SFLX_MV(   IA,JA) !> surface flux of meridional wind
    real(RP), intent(in) :: SFLX_SH(   IA,JA) !> surface sensible heat flux
    real(RP), intent(in) :: SFLX_QV(   IA,JA) !> surface sensible QV flux
    real(RP), intent(in) :: l_mo   (   IA,JA) !> Monin-Obukhov length

    real(RP), intent(in)  :: CZ(  KA,IA,JA)
    real(RP), intent(in)  :: FZ(0:KA,IA,JA)
    real(DP), intent(in)  :: dt_DP

    real(RP), intent(out) :: RHOU_t(KA,IA,JA) !> tendency of dens * u
    real(RP), intent(out) :: RHOV_t(KA,IA,JA) !> tendency of dens * v
    real(RP), intent(out) :: RHOT_t(KA,IA,JA) !> tendency of dens * pt
    real(RP), intent(out) :: RPROG_t(KA,IA,JA,ATMOS_PHY_BL_MYNN_ntracer) !> tenddency of dens * prognostic variables (TKE, TSQ, QSQ, COV)
    real(RP), intent(out) :: Nu    (KA,IA,JA) !> eddy viscosity coefficient @ half-level
    real(RP), intent(out) :: Nu_cg (KA,IA,JA) !> countergradient eddy viscosity coefficient @ half-level
    real(RP), intent(out) :: Kh    (KA,IA,JA) !> eddy diffusion coefficient @ half-level
    real(RP), intent(out) :: Kh_cg (KA,IA,JA) !> countergradient eddy diffusion coefficient @ half-level

    real(RP) :: Ri   (KA,IA,JA) !> Richardson number
    real(RP) :: Pr   (KA,IA,JA) !> Plandtle number
    real(RP) :: prod (KA,IA,JA) !> TKE production term
    real(RP) :: diss (KA,IA,JA) !> TKE dissipation term
    real(RP) :: dudz2(KA,IA,JA) !> (du/dz)^2 + (dv/dz)^2
    real(RP) :: l    (KA,IA,JA) !> length scale L
    real(RP) :: rho_h      !> dens at the half level

    real(RP) :: flxU(KA,IA,JA) !> dens * w * u
    real(RP) :: flxV(KA,IA,JA) !> dens * w * v
    real(RP) :: flxT(KA,IA,JA) !> dens * w * pt

    real(RP) :: RHONu (KA) !> dens * Nu at the half level for level 2.5
    real(RP) :: RHONuc(KA) !> dens * Nu at the half level for the countergradient
    real(RP) :: RHOKh (KA) !> dens * Kh at the half level for level 2.5
    real(RP) :: RHOKhc(KA) !> dens * Kh at the half level for the countergradient
    real(RP) :: N2_new(KA) !> squared Brunt-Baisala frequency
    real(RP) :: SFLX_PT    !> surface potential temperature flux
    real(RP) :: SFLX_PTV   !> surface virtual potential temperature flux
    real(RP) :: sm25  (KA) !> stability function for velocity for level 2.5
    real(RP) :: smp   (KA) !> stability function for velocity for the countergradient
    real(RP) :: sh25  (KA) !> stability function for scalars for level 2.5
    real(RP) :: shpgh (KA) !> stability function for scalars for the countergradient
    real(RP) :: Nu_f  (KA) !> Nu at the full level
    real(RP) :: Kh_f  (KA) !> Kh at the full level
    real(RP) :: q     (KA) !> q
    real(RP) :: q2_2  (KA) !> q^2 for level 2
    real(RP) :: ac    (KA) !> !> \alpha_c

    ! for level 3
    real(RP) :: tvsq  (KA) !> <\theta_v^2>
    real(RP) :: tvsq25(KA) !> <\theta_v^2> for level 2.5
    real(RP) :: tltv  (KA)
    real(RP) :: qwtv  (KA)
    real(RP) :: tltv25(KA)
    real(RP) :: qwtv25(KA)
    real(RP) :: tsq   (KA)
    real(RP) :: qsq   (KA)
    real(RP) :: cov   (KA)
    real(RP) :: wtl   (KA)
    real(RP) :: wqw   (KA)
    real(RP) :: tsq25
    real(RP) :: qsq25
    real(RP) :: cov25
    real(RP) :: prod_t1
    real(RP) :: prod_q1
    real(RP) :: prod_c1

    real(RP) :: dtldz(KA)
    real(RP) :: dqwdz(KA)
    real(RP) :: betat(KA)
    real(RP) :: betaq(KA)

    real(RP) :: gammat(KA)
    real(RP) :: gammaq(KA)

    real(RP) :: flx(KA)
    real(RP) :: a(KA)
    real(RP) :: b(KA)
    real(RP) :: c(KA)
    real(RP) :: d(KA)
    real(RP) :: ap
    real(RP) :: phi_N(KA)
    real(RP) :: tke_P(KA)

    real(RP) :: sf_t
    real(RP) :: us, us3
    real(RP) :: zeta
    real(RP) :: phi_m, phi_h
    real(RP) :: rl_mo

    real(RP) :: FDZ(KA)
    real(RP) :: CDZ(KA)
    real(RP) :: f2h(KA,2)
    real(RP) :: z1

    logical :: mynn_level3

    real(RP) :: dt

    real(RP) :: tmp, sw

    integer :: k, i, j
    integer :: nit, it
    !---------------------------------------------------------------------------

    dt = real(dt_DP, RP)

    LOG_PROGRESS(*) "atmosphere / physics / pbl / MYNN"

    mynn_level3 = ( ATMOS_PHY_BL_MYNN_LEVEL == "3" )

    if ( initialize ) then
       nit = KE_PBL - 1
    else
       nit = 1
    end if

!OCL INDEPENDENT
    !$omp parallel do default(none) &
    !$omp OMP_SCHEDULE_ collapse(2) &
    !$omp shared(KA,KS,KE_PBL,KE,IS,IE,JS,JE, &
    !$omp        EPS,GRAV,CPdry,EPSTvap,UNDEF,RSQRT_2,SQRT_2PI,RSQRT_2PI, &
    !$omp        ATMOS_PHY_BL_MYNN_N2_MAX,ATMOS_PHY_BL_MYNN_TKE_MIN, &
    !$omp        ATMOS_PHY_BL_MYNN_NU_MIN,ATMOS_PHY_BL_MYNN_NU_MAX, &
    !$omp        ATMOS_PHY_BL_MYNN_KH_MIN,ATMOS_PHY_BL_MYNN_KH_MAX, &
    !$omp        ATMOS_PHY_BL_MYNN_Sq_fact,ATMOS_PHY_BL_MYNN_prod_surf, &
    !$omp        RHOU_t,RHOV_t,RHOT_t,RPROG_t,Nu,Nu_cg,Kh,Kh_cg, &
    !$omp        DENS,PROG,U,V,POTT,PRES,QDRY,QV,Qw,POTV,POTL,EXNER,N2,SFLX_MU,SFLX_MV,SFLX_SH,SFLX_QV,l_mo, &
    !$omp        mynn_level3,initialize,nit, &
    !$omp        CZ,FZ,dt, &
    !$omp        Ri,Pr,prod,diss,dudz2,l,flxU,flxV,flxT) &
    !$omp private(N2_new,sm25,smp,sh25,shpgh,Nu_f,Kh_f,q,q2_2,ac,SFLX_PT,SFLX_PTV,RHONu,RHONuc,RHOKh,RHOKhc, &
    !$omp         dtldz,dqwdz,betat,betaq,gammat,gammaq,wtl,wqw, &
    !$omp         flx,a,b,c,d,ap,rho_h,phi_n,tke_P,sf_t,rl_mo,zeta,phi_m,phi_h,us,us3,CDZ,FDZ,f2h,z1, &
    !$omp         tvsq,tsq,qsq,cov,tvsq25,tsq25,qsq25,cov25,tltv,qwtv,tltv25,qwtv25,prod_t1,prod_q1,prod_c1, &
    !$omp         sw,tmp, &
    !$omp         k,i,j,it)
    do j = JS, JE
    do i = IS, IE

       z1 = CZ(KS,i,j) - FZ(KS-1,i,j)

       SFLX_PT = SFLX_SH(i,j) / ( CPdry * DENS(KS,i,j) * EXNER(KS,i,j) )

       do k = KS, KE_PBL
          FDZ(k) = CZ(k+1,i,j) - CZ(k  ,i,j)
       end do
       do k = KS, KE_PBL+1
          CDZ(k) = FZ(k  ,i,j) - FZ(k-1,i,j)
       end do

       call get_f2h( &
            KA, KS, KE_PBL, &
            CDZ(:),  & ! (in)
            f2h(:,:) ) ! (out)

       call calc_vertical_differece( KA, KS, KE_PBL, &
                                     U(:,i,j), V(:,i,j), POTL(:,i,j), Qw(:,i,j), & ! (in)
                                     CDZ(:), FDZ(:), f2h(:,:),                   & ! (in)
                                     dudz2(:,i,j), dtldz(:), dqwdz(:)            ) ! (out)

       do k = KS, KE_PBL
          n2_new(k) = min( max( N2(k,i,j), - ATMOS_PHY_BL_MYNN_N2_MAX ), - ATMOS_PHY_BL_MYNN_N2_MAX )
          Ri(k,i,j) = n2_new(k) / max(dudz2(k,i,j), 1E-10_RP)
       end do

       do k = KS, KE_PBL
          q(k) = sqrt( max(PROG(k,i,j,I_TKE), ATMOS_PHY_BL_MYNN_TKE_MIN)*2.0_RP )
       end do

       if ( mynn_level3 ) then
          do k = KS, KE_PBL
             tsq(k) = max( PROG(k,i,j,I_TSQ), 0.0_RP )
             qsq(k) = max( PROG(k,i,j,I_QSQ), 0.0_RP )
             cov(k) = PROG(k,i,j,I_COV)
             cov(k) = sign( min( abs(cov(k)), sqrt(tsq(k) * qsq(k))), cov(k) )
          end do
       end if


       flx(KS-1  ) = 0.0_RP
       flx(KE_PBL) = 0.0_RP


       do it = 1, nit

          ! length
          SFLX_PTV = SFLX_PT * ( 1.0_RP + EPSTvap * Qw(KS,i,j) ) + SFLX_QV(i,j) * EPSTvap * POTT(KS,i,j)
          call get_length( &
               KA, KS, KE_PBL, &
               POTV(KS,i,j), q(:), n2_new(:), & ! (in)
               SFLX_PTV, l_mo(i,j),           & ! (in)
               FZ(:,i,j),                     & ! (in)
               l(:,i,j)                       ) ! (out)

          call get_q2_level2( &
               KA, KS, KE_PBL, &
               dudz2(:,i,j), Ri(:,i,j), l(:,i,j), & ! (in)
               q2_2(:)                            ) ! (out)

          if ( initialize .and. it==1 ) then
             do k = KS, KE_PBL
                q(k) = sqrt( q2_2(k) )
             end do
          end if

          do k = KS, KE_PBL
             ac(k) = min( q(k) / sqrt(q2_2(k)), 1.0_RP )
          end do

          call get_smsh( &
               KA, KS, KE_PBL,            & ! (in)
               q(:), ac(:),               & ! (in)
               l(:,i,j), n2_new(:),       & ! (in)
               POTV(:,i,j), dudz2(:,i,j), & ! (in)
               tvsq(:), tvsq25(:),        & ! (in) ! dummy
               tltv(:), tltv25(:),        & ! (in) ! dummy
               qwtv(:), qwtv25(:),        & ! (in) ! dummy
               .false.,                   & ! (in)
               sm25(:), smp(:),           & ! (out)
               sh25(:), shpgh(:),         & ! (out) ! dummy
               gammat(:), gammaq(:)       ) ! (out) ! dummy


          call partial_condensation( KA, KS, KE_PBL, &
                                     PRES(:,i,j), POTT(:,i,j),  & ! (in)
                                     POTL(:,i,j), Qw(:,i,j),    & ! (in)
                                     Qdry(:,i,j), EXNER(:,i,j), & ! (in)
                                     dtldz(:), dqwdz(:),        & ! (in)
                                     l(:,i,j), sh25(:), ac(:),  & ! (in)
                                     tsq(:), qsq(:), cov(:),    & ! (in)
                                     mynn_level3,               & ! (in)
                                     betat(:), betaq(:)         ) ! (out)

          if ( ATMOS_PHY_BL_MYNN_prod_surf ) then
             us3 = - l_mo(i,j) * KARMAN * GRAV * SFLX_PTV / POTV(KS,i,j) ! u_*^3
             us = max(us3, 0.0_RP)**(1.0_RP/3.0_RP)
             !us = sqrt( sqrt( SFLX_MU(i,j)**2 + SFLX_MV(i,j)**2 ) )
             rl_mo = 1.0_RP / sign( max(abs(l_mo(i,j)), EPS), l_mo(i,j) )
             zeta = z1 * rl_mo

!!$             ! Businger et al. (1971)
!!$             if ( zeta >= 0 ) then
!!$                phi_m = 4.7_RP * zeta + 1.0_RP
!!$                phi_h = 4.7_RP * zeta + 0.74_RP
!!$             else
!!$                phi_m = 1.0_RP / sqrt(sqrt( 1.0_RP - 15.0_RP * zeta ))
!!$                phi_h = 0.47_RP / sqrt( 1.0_RP - 9.0_RP * zeta )
!!$             end if

             ! Beljaars and Holtslag (1991)
             if ( zeta >= 0 ) then
                tmp = - 2.0_RP / 3.0_RP * ( 0.35_RP * zeta - 6.0_RP ) * exp(-0.35_RP*zeta)
                phi_m = tmp * zeta + zeta + 1.0_RP
                phi_h = tmp + zeta * sqrt( 1.0_RP + 2.0_RP * zeta / 3.0_RP ) + 1.0_RP
!!$             else
!!$                tmp = sqrt( 1.0_RP - 16.0_RP * zeta )
!!$                phi_m = 1.0_RP / sqrt(tmp)
!!$                phi_h = 1.0_RP / tmp
             end if

             ! Wilson (2001)
             if ( zeta < 0 ) then
                tmp = (-zeta)**(2.0_RP/3.0_RP)
                phi_m = 1.0_RP / sqrt( 1.0_RP + 3.6_RP * tmp )
                phi_h = 0.95_RP / sqrt( 1.0_RP + 7.9_RP * tmp )
             end if

          end if

          if ( mynn_level3 ) then

             if ( ATMOS_PHY_BL_MYNN_prod_surf ) then
                ! TSQ
                prod_t1 = 1.0_RP / us * phi_h / ( KARMAN * z1 ) * SFLX_PT**2
                ! QSQ
                prod_q1 = 1.0_RP / us * phi_h / ( KARMAN * z1 ) * ( SFLX_QV(i,j) / DENS(KS,i,j) )**2
                ! COV
                prod_c1 = 1.0_RP / us * phi_h * ( KARMAN * z1 ) * SFLX_PT * SFLX_QV(i,j) / DENS(KS,i,j)
             end if

             do k = KS, KE_PBL

                if ( ATMOS_PHY_BL_MYNN_prod_surf .and. k == KS ) then
                   tsq25 = prod_t1 * B2 * l(k,i,j) / q(k) * 0.5_RP
                   qsq25 = prod_q1 * B2 * l(k,i,j) / q(k) * 0.5_RP
                   cov25 = prod_c1 * B2 * l(k,i,j) / q(k) * 0.5_RP
                else
                   tsq25 = B2 * l(k,i,j)**2 * sh25(k) * dtldz(k)**2
                   qsq25 = B2 * l(k,i,j)**2 * sh25(k) * dqwdz(k)**2
                   cov25 = B2 * l(k,i,j)**2 * sh25(k) * dtldz(k) * dqwdz(k)
                end if

                tltv25(k) = betat(k) * tsq25 + betaq(k) * cov25
                qwtv25(k) = betat(k) * cov25 + betaq(k) * qsq25
                tvsq25(k) = max(betat(k) * tltv25(k) + betaq(k) * qwtv25(k), 0.0_RP)

                if ( initialize .and. it==1 ) then
                   ! initialize teq, qsq, and cov
                   tsq(k) = tsq25
                   qsq(k) = qsq25
                   cov(k) = cov25
                end if
                tltv(k) = betat(k) * tsq(k) + betaq(k) * cov(k)
                qwtv(k) = betat(k) * cov(k) + betaq(k) * qsq(k)
                tvsq(k) = max(betat(k) * tltv(k) + betaq(k) * qwtv(k), 0.0_RP)
             end do

          end if

          ! update N2
          do k = KS, KE_PBL
             n2_new(k) = min(ATMOS_PHY_BL_MYNN_N2_MAX, &
                             GRAV * ( dtldz(k) * betat(k) + dqwdz(k) * betaq(k) ) / POTV(k,i,j) )
          end do

          do k = KS, KE_PBL
             Ri(k,i,j) = n2_new(k) / max(dudz2(k,i,j), 1E-10_RP)
          end do

          ! length
          SFLX_PTV = SFLX_PT * betat(KS) + SFLX_QV(i,j) * betaq(KS)
          call get_length( &
               KA, KS, KE_PBL, &
               POTV(KS,i,j), q(:), n2_new(:), & ! (in)
               SFLX_PTV, l_mo(i,j),           & ! (in)
               FZ(:,i,j),                     & ! (in)
               l(:,i,j)                       ) ! (out)

          call get_q2_level2( &
               KA, KS, KE_PBL, &
               dudz2(:,i,j), Ri(:,i,j), l(:,i,j), & ! (in)
               q2_2(:)                            ) ! (out)

          do k = KS, KE_PBL
             ac(k) = min( q(k) / sqrt(q2_2(k)), 1.0_RP )
          end do

          call get_smsh( &
               KA, KS, KE_PBL,            & ! (in)
               q(:), ac(:),               & ! (in)
               l(:,i,j), n2_new(:),       & ! (in)
               POTV(:,i,j), dudz2(:,i,j), & ! (in)
               tvsq(:), tvsq25(:),        & ! (in)
               tltv(:), tltv25(:),        & ! (in)
               qwtv(:), qwtv25(:),        & ! (in)
               mynn_level3,               & ! (in)
               sm25(:), smp(:),           & ! (out)
               sh25(:), shpgh(:),         & ! (out)
               gammat(:), gammaq(:)       ) ! (out)


          do k = KS, KE_PBL
             Nu_f(k) = l(k,i,j) * q(k) * sm25(k)
             Kh_f(k) = l(k,i,j) * q(k) * sh25(k)
          end do
          do k = KS, KE_PBL-1
             Nu(k,i,j) = min( f2h(k,1) * Nu_f(k+1) + f2h(k,2) * Nu_f(k), &
                              ATMOS_PHY_BL_MYNN_NU_MAX )
             Kh(k,i,j) = min( f2h(k,1) * Kh_f(k+1) + f2h(k,2) * Kh_f(k), &
                              ATMOS_PHY_BL_MYNN_KH_MAX )
          end do

          do k = KS, KE_PBL
             Nu_f(k) = l(k,i,j) * q(k) * smp(k)
             Kh_f(k) = l(k,i,j) * q(k) * shpgh(k)
          end do
          do k = KS, KE_PBL-1
             tmp = f2h(k,1) * Nu_f(k+1) + f2h(k,2) * Nu_f(k)
             if ( tmp > 0.0_RP ) then
                Nu   (k,i,j) = min( Nu(k,i,j) + tmp, ATMOS_PHY_BL_MYNN_NU_MAX )
                Nu_cg(k,i,j) = 0.0_RP
             else
                Nu   (k,i,j) = max( Nu(k,i,j), - tmp * 1E-4_RP ) ! for numerical stability
                Nu_cg(k,i,j) = max( tmp, ATMOS_PHY_BL_MYNN_NU_MIN )
             end if

             tmp = - ( f2h(k,1) * n2_new(k+1) + f2h(k,2) * n2_new(k) ) ! gh * (q/L)^2
             tmp = ( f2h(k,1) * Kh_f(k+1) + f2h(k,2) * Kh_f(k) ) / sign( max(abs(tmp),EPS), tmp )
             if ( tmp > 0.0_RP ) then
                Kh   (k,i,j) = min( Kh(k,i,j) + tmp, ATMOS_PHY_BL_MYNN_KH_MAX )
                Kh_cg(k,i,j) = 0.0_RP
             else
                Kh   (k,i,j) = max( Kh(k,i,j), - tmp * 1E-4_RP ) ! for numerical stability
                Kh_cg(k,i,j) = max( tmp, ATMOS_PHY_BL_MYNN_KH_MIN )
             end if
          end do

          if ( ATMOS_PHY_BL_MYNN_prod_surf ) then
             Nu(KS,i,j) = KARMAN * CDZ(KS) * us / phi_m
             Kh(KS,i,j) = KARMAN * CDZ(KS) * us / phi_h
             Nu_cg(KS,i,j) = 0.0_RP
             Kh_cg(KS,i,j) = 0.0_RP
          end if

          do k = KS, KE_PBL-1
             sw = 0.5_RP - sign(0.5_RP, abs(Kh(k,i,j)) - EPS)
             Pr(k,i,j) = Nu(k,i,j) / ( Kh(k,i,j) + sw ) * ( 1.0_RP - sw ) &
                       + 1.0_RP * sw
          end do

          do k = KS, KE_PBL-1
             tmp = f2h(k,1) * DENS(k+1,i,j) + f2h(k,2) * DENS(k,i,j)
             RHONu (k) = tmp * Nu   (k,i,j)
             RHONuc(k) = tmp * Nu_cg(k,i,j)
             RHOKh (k) = tmp * Kh   (k,i,j)
             RHOKhc(k) = tmp * Kh_cg(k,i,j)
          end do

          ! time integration

          if ( it == nit ) then
             ! dens * u

             ! countergradient flux
             do k = KS, KE_PBL-1
                flx(k) = min( max( &
                     - RHONuc(k) * ( U(k+1,i,j) - U(k,i,j) ) / FDZ(k), &
                     - abs(U(k+1,i,j)) * DENS(k+1,i,j) * CDZ(k+1) / dt * FLX_LIM_FACT ), &
                       abs(U(k  ,i,j)) * DENS(k  ,i,j) * CDZ(k  ) / dt * FLX_LIM_FACT )
             end do

             sf_t = SFLX_MU(i,j) / CDZ(KS)
             d(KS) = U(KS,i,j) + dt * ( sf_t - flx(KS) / CDZ(KS) ) / DENS(KS,i,j)
             do k = KS+1, KE_PBL
                d(k) = U(k,i,j) - dt * ( flx(k) - flx(k-1) ) / ( CDZ(k) * DENS(k,i,j) )
             end do
             c(KS) = 0.0_RP
             do k = KS, KE_PBL-1
                ap = - dt * RHONu(k) / FDZ(k)
                a(k) = ap / ( DENS(k,i,j) * CDZ(k) )
                b(k) = - a(k) - c(k) + 1.0_RP
                c(k+1) = ap / ( DENS(k+1,i,j) * CDZ(k+1) )
             end do
             a(KE_PBL) = 0.0_RP
             b(KE_PBL) = - c(KE_PBL) + 1.0_RP

             call MATRIX_SOLVER_tridiagonal( &
                  KA, KS, KE_PBL, &
                  a(:), b(:), c(:), d(:), & ! (in)
                  phi_n(:)                ) ! (out)

             RHOU_t(KS,i,j) = ( phi_n(KS) - U(KS,i,j) ) * DENS(KS,i,j) / dt - sf_t
             do k = KS+1, KE_PBL
                RHOU_t(k,i,j) = ( phi_n(k) - U(k,i,j) ) * DENS(k,i,j) / dt
             end do
             do k = KE_PBL+1, KE
                RHOU_t(k,i,j) = 0.0_RP
             end do
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
                flx(k) = min( max( &
                     - RHONuc(k) * ( V(k+1,i,j) - V(k,i,j) ) / FDZ(k), &
                     - abs(V(k+1,i,j)) * DENS(k+1,i,j) * CDZ(k+1) / dt * FLX_LIM_FACT ), &
                       abs(V(k  ,i,j)) * DENS(k  ,i,j) * CDZ(k  ) / dt * FLX_LIM_FACT )
             end do

             sf_t = SFLX_MV(i,j) / CDZ(KS)
             d(KS) = V(KS,i,j) + dt * ( sf_t - flx(KS) / CDZ(KS) ) / DENS(KS,i,j)
             do k = KS+1, KE_PBL
                d(k) = V(k,i,j) - dt * ( flx(k) - flx(k-1) ) / ( CDZ(k) * DENS(k,i,j) )
             end do
             ! a,b,c is the same as those for the u

             call MATRIX_SOLVER_tridiagonal( &
                  KA, KS, KE_PBL, &
                  a(:), b(:), c(:), d(:), & ! (in)
                  phi_n(:)                ) ! (out)

             RHOV_t(KS,i,j) = ( phi_n(KS) - V(KS,i,j) ) * DENS(KS,i,j) / dt - sf_t
             do k = KS+1, KE_PBL
                RHOV_t(k,i,j) = ( phi_n(k) - V(k,i,j) ) * DENS(k,i,j) / dt
             end do
             do k = KE_PBL+1, KE
                RHOV_t(k,i,j) = 0.0_RP
             end do
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
                flx(k) = min( max( &
                     - RHOKhc(k) * ( POTT(k+1,i,j) - POTT(k,i,j) ) / FDZ(k), &
                     - POTT(k+1,i,j) * DENS(k+1,i,j) * CDZ(k+1) / dt * FLX_LIM_FACT ), &
                       POTT(k  ,i,j) * DENS(k  ,i,j) * CDZ(k  ) / dt * FLX_LIM_FACT )
             end do

             sf_t = SFLX_PT / CDZ(KS)
             d(KS) = POTT(KS,i,j) + dt * ( sf_t - flx(KS) / ( CDZ(KS) * DENS(KS,i,j) ) )
             do k = KS+1, KE_PBL
                d(k) = POTT(k,i,j) - dt * ( flx(k) - flx(k-1) ) / ( CDZ(k) * DENS(k,i,j) )
             end do

             c(KS) = 0.0_RP
             do k = KS, KE_PBL-1
                ap = - dt * RHOKh(k) / FDZ(k)
                a(k) = ap / ( DENS(k,i,j) * CDZ(k) )
                b(k) = - a(k) - c(k) + 1.0_RP
                c(k+1) = ap / ( DENS(k+1,i,j) * CDZ(k+1) )
             end do
             a(KE_PBL) = 0.0_RP
             b(KE_PBL) = - c(KE_PBL) + 1.0_RP

             call MATRIX_SOLVER_tridiagonal( &
                  KA, KS, KE_PBL, &
                  a(:), b(:), c(:), d(:), & ! (in)
                  phi_n(:)                ) ! (out)

             RHOT_t(KS,i,j) = ( ( phi_n(KS) - POTT(KS,i,j) ) / dt - sf_t ) * DENS(KS,i,j)
             do k = KS+1, KE_PBL
                RHOT_t(k,i,j) = ( phi_n(k) - POTT(k,i,j) ) * DENS(k,i,j) / dt
             end do
             do k = KE_PBL+1, KE
                RHOT_t(k,i,j) = 0.0_RP
             end do
             do k = KS, KE_PBL-1
                flxT(k,i,j) = flx(k) &
                            - RHOKh(k) * ( phi_n(k+1) - phi_n(k) ) / CDZ(k)
             end do
             do k = KE_PBL, KE
                flxT(k,i,j) = 0.0_RP
             end do

          end if


          ! dens * TKE

          ! production
          Nu  (KS-1,i,j) = 0.0_RP
          Nu_cg(KS-1,i,j) = 0.0_RP
          Kh  (KS-1,i,j) = 0.0_RP
          Kh_cg(KS-1,i,j) = 0.0_RP
          do k = KS, KE_PBL
             prod(k,i,j) = l(k,i,j) * q(k) * ( ( sm25(k) + smp(k) ) * dudz2(k,i,j) &
                                             - ( sh25(k) * n2_new(k) - shpgh(k) ) )
          end do
          if ( ATMOS_PHY_BL_MYNN_prod_surf ) then
             prod(KS,i,j) = us3 * phi_m / ( KARMAN * z1 )
          end if

          do k = KS, KE_PBL
             tke_p(k) = q(k)**2 * 0.5_RP
             prod(k,i,j) = max( prod(k,i,j), - tke_p(k) / dt )
          end do
          do k = KE_PBL+1, KE
             prod(k,i,j) = 0.0_RP
          end do

          ! countergradient flux
          do k = KS, KE_PBL-1
             flx(k) = min( max( &
                  - ATMOS_PHY_BL_MYNN_Sq_fact * RHONuc(k) * ( tke_p(k+1) - tke_p(k) ) / FDZ(k), &
                  - tke_p(k+1) * DENS(k+1,i,j) * CDZ(k+1) / dt * FLX_LIM_FACT ), &
                    tke_p(k  ) * DENS(k  ,i,j) * CDZ(k  ) / dt * FLX_LIM_FACT )
          end do

          do k = KS, KE_PBL
             d(k) = tke_p(k) + dt * ( prod(k,i,j) - ( flx(k) - flx(k-1) ) / ( CDZ(k) * DENS(k,i,j) ))
          end do

          c(KS) = 0.0_RP
          do k = KS, KE_PBL-1
             ap = - dt * ATMOS_PHY_BL_MYNN_Sq_fact * RHONu(k) / FDZ(k)
             a(k) = ap / ( DENS(k,i,j) * CDZ(k) )
             diss(k,i,j) = - 2.0_RP * q(k) / ( B1 * l(k,i,j) )
             b(k) = - a(k) - c(k) + 1.0_RP - diss(k,i,j) * dt
             c(k+1) = ap / ( DENS(k+1,i,j) * CDZ(k+1) )
          end do
          a(KE_PBL) = 0.0_RP
          diss(KE_PBL,i,j) = - 2.0_RP * q(KE_PBL) / ( B1 * l(KE_PBL,i,j) )
          b(KE_PBL) = - c(KE_PBL) + 1.0_RP - diss(KE_PBL,i,j) * dt
          do k = KE_PBL+1, KE
             diss(k,i,j) = 0.0_RP
          end do

          call MATRIX_SOLVER_tridiagonal( &
               KA, KS, KE_PBL, &
               a(:), b(:), c(:), d(:), & ! (in)
               phi_n(:)                ) ! (out)

          do k = KS, KE_PBL
             phi_n(k) = max( phi_n(k), ATMOS_PHY_BL_MYNN_TKE_MIN )
          end do

          if ( it == nit ) then
             do k = KS, KE_PBL
                diss(k,i,j) = diss(k,i,j) * phi_n(k)
                RPROG_t(k,i,j,I_TKE) = ( phi_n(k) - PROG(k,i,j,I_TKE) ) * DENS(k,i,j) / dt
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

          do k = KS, KE_PBL
             wtl(k) = - l(k,i,j) * q(k) * ( sh25(k) * dtldz(k) + gammat(k) )
             wqw(k) = - l(k,i,j) * q(k) * ( sh25(k) * dqwdz(k) + gammaq(k) )
          end do

          ! dens * tsq

          ! countergradient flux
          do k = KS, KE_PBL-1
             flx(k) = min( max( &
                  - RHONuc(k) * ( tsq(k+1) - tsq(k) ) / FDZ(k), &
                  - tsq(k+1) * DENS(k+1,i,j) * CDZ(k+1) / dt * FLX_LIM_FACT ), &
                    tsq(k  ) * DENS(k  ,i,j) * CDZ(k  ) / dt * FLX_LIM_FACT )
          end do

          do k = KS, KE_PBL
             d(k) = max( tsq(k) &
                       + dt * ( - 2.0_RP * wtl(k) * dtldz(k) &
                                - ( flx(k) - flx(k-1) ) / ( CDZ(k) * DENS(k,i,j) ) ), &
                       0.0_RP)
          end do
          if ( ATMOS_PHY_BL_MYNN_prod_surf ) then
             d(KS) = max( tsq(KS) + dt * ( prod_t1 - flx(KS) / CDZ(KS) ), &
                          0.0_RP)
          end if
          c(KS) = 0.0_RP
          do k = KS, KE_PBL-1
             ap = - dt * RHONu(k) / FDZ(k)
             a(k) = ap / ( DENS(k,i,j) * CDZ(k) )
             b(k) = - a(k) - c(k) + 1.0_RP + dt * 2.0_RP * q(k) / ( B2 * l(k,i,j) )
             c(k+1) = ap / ( DENS(k+1,i,j) * CDZ(k+1) )
          end do
          a(KE_PBL) = 0.0_RP
          b(KE_PBL) = - c(KE_PBL) + 1.0_RP + dt * 2.0_RP * q(KE_PBL) / ( B2 * l(KE_PBL,i,j) )

          call MATRIX_SOLVER_tridiagonal( &
               KA, KS, KE_PBL, &
               a(:), b(:), c(:), d(:), & ! (in)
               tsq(:)                  ) ! (out)

          do k = KS, KE_PBL
             tsq(k) = max( tsq(k), 0.0_RP )
          end do

          if ( it == nit ) then
             do k = KS, KE_PBL
                RPROG_t(k,i,j,I_TSQ) = ( tsq(k) - PROG(k,i,j,I_TSQ) ) * DENS(k,i,j) / dt
             end do
             do k = KE_PBL+1, KE
                RPROG_t(k,i,j,I_TSQ) = 0.0_RP
             end do
          end if


          ! dens * qsq

          ! countergradient flux
          do k = KS, KE_PBL-1
             flx(k) = min( max( &
                  - RHONuc(k) * ( qsq(k+1) - qsq(k) ) / FDZ(k), &
                  - qsq(k+1) * DENS(k+1,i,j) * CDZ(k+1) / dt * FLX_LIM_FACT ), &
                    qsq(k  ) * DENS(k  ,i,j) * CDZ(k  ) / dt * FLX_LIM_FACT )
          end do

          do k = KS, KE_PBL
             d(k) = max( qsq(k) &
                       + dt * ( - 2.0_RP * wqw(k) * dqwdz(k) &
                              - ( flx(k) - flx(k-1) ) / ( CDZ(k) * DENS(k,i,j) ) ), &
                        0.0_RP)
          end do
          if ( ATMOS_PHY_BL_MYNN_prod_surf ) then
             d(KS) = max( qsq(KS) + dt * ( prod_q1 - flx(KS) / CDZ(KS) ), &
                          0.0_RP)
          end if
          ! a, b, c are same as those for tsq

          call MATRIX_SOLVER_tridiagonal( &
               KA, KS, KE_PBL, &
               a(:), b(:), c(:), d(:), & ! (in)
               qsq(:)                  ) ! (out)

          do k = KS, KE_PBL
             qsq(k) = max( qsq(k), 0.0_RP )
          end do

          if ( it == nit ) then
             do k = KS, KE_PBL
                RPROG_t(k,i,j,I_QSQ) = ( qsq(k) - PROG(k,i,j,I_QSQ) ) * DENS(k,i,j) / dt
             end do
             do k = KE_PBL+1, KE
                RPROG_t(k,i,j,I_QSQ) = 0.0_RP
             end do
          end if


          ! dens * cov

          ! countergradient flux
          do k = KS, KE_PBL-1
             flx(k) = min( max( &
                  - RHONuc(k) * ( cov(k+1) - cov(k) ) / FDZ(k), &
                  - abs(cov(k+1)) * DENS(k+1,i,j) * CDZ(k+1) / dt * FLX_LIM_FACT ), &
                    abs(cov(k  )) * DENS(k  ,i,j) * CDZ(k  ) / dt * FLX_LIM_FACT )
          end do

          do k = KS, KE_PBL
             d(k) = cov(k) &
                  + dt * ( - wqw(k) * dtldz(k) - wtl(k) * dqwdz(k) &
                         - ( flx(k) - flx(k-1) ) / ( CDZ(k) * DENS(k,i,j) ) )
          end do
          if ( ATMOS_PHY_BL_MYNN_prod_surf ) then
             d(KS) = cov(KS) + dt * ( prod_c1 - flx(KS) / CDZ(KS) )
          end if
          ! a, b, c are same as those for tsq

          call MATRIX_SOLVER_tridiagonal( &
               KA, KS, KE_PBL, &
               a(:), b(:), c(:), d(:), & ! (in)
               cov(:)                  ) ! (out)

          do k = KS, KE_PBL
             cov(k) = sign( min( abs(cov(k)), sqrt(tsq(k)*qsq(k)) ), cov(k) )
          end do

          if ( it == nit ) then
             do k = KS, KE_PBL
                RPROG_t(k,i,j,I_COV) = ( cov(k) - PROG(k,i,j,I_COV) ) * DENS(k,i,j) / dt
             end do
             do k = KE_PBL+1, KE
                RPROG_t(k,i,j,I_COV) = 0.0_RP
             end do
          end if

       end do

    end do
    end do


    do j = JS, JE
    do i = IS, IE
       do k = KE_PBL, KE
          Nu   (k,i,j) = 0.0_RP
          Nu_cg(k,i,j) = 0.0_RP
          Kh   (k,i,j) = 0.0_RP
          Kh_cg(k,i,j) = 0.0_RP
          Pr   (k,i,j) = 1.0_RP
       end do
       do k = KE_PBL+1, KE
          Ri   (k,i,j) = UNDEF
          prod (k,i,j) = UNDEF
          diss (k,i,j) = UNDEF
          dudz2(k,i,j) = UNDEF
          l    (k,i,j) = UNDEF
       end do
    end do
    end do


    call FILE_HISTORY_in(Ri(:,:,:), 'Ri_MYNN', 'Richardson number', '1',     fill_halo=.true. )
    call FILE_HISTORY_in(Pr(:,:,:), 'Pr_MYNN', 'Prandtl number',    '1',     fill_halo=.true., dim_type="ZHXY" )
    call FILE_HISTORY_in(prod(:,:,:), 'TKE_prod_MYNN', 'TKE production',  'm2/s3', fill_halo=.true.)
    call FILE_HISTORY_in(diss(:,:,:), 'TKE_diss_MYNN', 'TKE dissipation', 'm2/s3', fill_halo=.true.)
    call FILE_HISTORY_in(dudz2(:,:,:), 'dUdZ2_MYNN', 'dudz2', 'm2/s2', fill_halo=.true.)
    call FILE_HISTORY_in(l(:,:,:), 'L_mix_MYNN', 'minxing length', 'm', fill_halo=.true.)

    call FILE_HISTORY_in(flxU(:,:,:), 'ZFLX_RHOU_MYNN', 'Z FLUX of RHOU (MYNN)', 'kg/m/s2', fill_halo=.true., dim_type="ZHXY" )
    call FILE_HISTORY_in(flxV(:,:,:), 'ZFLX_RHOV_MYNN', 'Z FLUX of RHOV (MYNN)', 'kg/m/s2', fill_halo=.true., dim_type="ZHXY" )
    call FILE_HISTORY_in(flxT(:,:,:), 'ZFLX_RHOT_MYNN', 'Z FLUX of RHOT (MYNN)', 'K kg/m2/s', fill_halo=.true., dim_type="ZHXY" )


    initialize = .false.

    return
  end subroutine ATMOS_PHY_BL_MYNN_tendency

  !-----------------------------------------------------------------------------
  !> ATMOS_PHY_BL_MYNN_tendency_tracer
  !! calculate tendency of tracers by the eddy viscosity
  !<
  subroutine ATMOS_PHY_BL_MYNN_tendency_tracer( &
       KA, KS, KE, IA, IS, IE, JA, JS, JE, &
       DENS, QTRC, SFLX_Q, &
       Kh, Kh_cg, MASS,    &
       CZ, FZ, DDT,        &
       TRACER_NAME,        &
       RHOQ_t              )
    use scale_matrix, only: &
       MATRIX_SOLVER_tridiagonal
    use scale_file_history, only: &
       FILE_HISTORY_in
    integer, intent(in) :: KA, KS, KE
    integer, intent(in) :: IA, IS, IE
    integer, intent(in) :: JA, JS, JE

    real(RP),         intent(in) :: DENS  (KA,IA,JA) !> density
    real(RP),         intent(in) :: QTRC  (KA,IA,JA) !> tracers
    real(RP),         intent(in) :: SFLX_Q(   IA,JA) !> surface flux
    real(RP),         intent(in) :: Kh    (KA,IA,JA) !> eddy diffusion coefficient @ half-level
    real(RP),         intent(in) :: Kh_cg (KA,IA,JA) !> eddy diffusion coefficient for the countergradient @ half-level
    real(RP),         intent(in) :: MASS             !> mass
    real(RP),         intent(in) :: CZ (  KA,IA,JA)  !> z at the full level
    real(RP),         intent(in) :: FZ (0:KA,IA,JA)  !> z at the half level
    real(DP),         intent(in) :: DDT              !> time step
    character(len=*), intent(in) :: TRACER_NAME      !> name of tracer (for history output)

    real(RP), intent(out) :: RHOQ_t(KA,IA,JA) !> tendency of tracers

    real(RP) :: QTRC_n(KA) !> value at the next time step
    real(RP) :: RHO   (KA)
    real(RP) :: RHOKh (KA)
    real(RP) :: RHOKhc(KA)
    real(RP) :: a(KA)
    real(RP) :: b(KA)
    real(RP) :: c(KA)
    real(RP) :: d(KA)
    real(RP) :: rho_h
    real(RP) :: ap
    real(RP) :: sf_t

    real(RP) :: flx(KA,IA,JA)

    real(RP) :: CDZ(KA)
    real(RP) :: FDZ(KA)
    real(RP) :: f2h(KA,2) !> coefficient to convert from full to half level

    real(RP) :: dt

    integer :: k, i, j

    dt = real( DDT, kind=RP )

!OCL INDEPENDENT
    !$omp parallel do default(none) OMP_SCHEDULE_ collapse(2) &
    !$omp shared(KA,KS,KE_PBL,KE,IS,IE,JS,JE) &
    !$omp shared(RHOQ_t,DENS,QTRC,SFLX_Q,Kh,Kh_cg,MASS,CZ,FZ,DT,flx) &
    !$omp private(QTRC_n,RHO,RHOKh,RHOKhc,rho_h,a,b,c,d,ap,sf_t,CDZ,FDZ,f2h) &
    !$omp private(k,i,j)
    do j = JS, JE
    do i = IS, IE

       do k = KS, KE_PBL
          CDZ(k) = FZ(k  ,i,j) - FZ(k-1,i,j)
          FDZ(k) = CZ(k+1,i,j) - CZ(k  ,i,j)
       end do

       call get_f2h( &
            KA, KS, KE_PBL-1, &
            CDZ(:),  & ! (in)
            f2h(:,:) ) ! (out)

       sf_t = SFLX_Q(i,j) / CDZ(KS)

       RHO(KS) = DENS(KS,i,j) + dt * sf_t * MASS
       do k = KS+1, KE_PBL
          RHO(k) = DENS(k,i,j)
       end do

       ! dens * coefficient at the half level
       do k = KS, KE_PBL-1
          rho_h = f2h(k,1) * DENS(k+1,i,j) + f2h(k,2) * DENS(k,i,j)
          RHOKh (k) = rho_h * Kh   (k,i,j)
          RHOKhc(k) = rho_h * Kh_cg(k,i,j)
       end do

       ! countergradient flux
       do k = KS, KE_PBL-1
          flx(k,i,j) = min( max( &
               - RHOKhc(k) * ( QTRC(k+1,i,j) - QTRC(k,i,j) ) / FDZ(k), &
               - abs(QTRC(k+1,i,j)) * DENS(k+1,i,j) * CDZ(k+1) / dt * FLX_LIM_FACT ), &
                 abs(QTRC(k  ,i,j)) * DENS(k  ,i,j) * CDZ(k  ) / dt * FLX_LIM_FACT )
       end do
       do k = KE_PBL, KE
          flx(k,i,j) = 0.0_RP
       end do

       d(KS) = ( QTRC(KS,i,j) * DENS(KS,i,j) &
               + dt * ( sf_t - flx(KS,i,j) / ( CDZ(KS) * RHO(KS) ) ) &
               ) / RHO(KS)
       do k = KS+1, KE_PBL
          d(k) = QTRC(k,i,j) - dt * ( flx(k,i,j) - flx(k-1,i,j) ) / ( RHO(k) * CDZ(k) )
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
               QTRC_n(:)               ) ! (out)

       RHOQ_t(KS,i,j) = ( QTRC_n(KS) * RHO(KS) - QTRC(KS,i,j) * DENS(KS,i,j) ) / dt - sf_t
       do k = KS+1, KE_PBL
          RHOQ_t(k,i,j) = ( QTRC_n(k) - QTRC(k,i,j) ) * RHO(k) / dt
       end do
       do k = KE_PBL+1, KE
          RHOQ_t(k,i,j) = 0.0_RP
       end do

       do k = KS, KE_PBL-1
          flx(k,i,j) = flx(k,i,j) &
                     - RHOKh(k) * ( QTRC_n(k+1) - QTRC_n(k) ) / FDZ(k)
       end do

    end do
    end do

    call FILE_HISTORY_in(flx(:,:,:), 'ZFLX_'//trim(TRACER_NAME)//'_MYNN', 'Z FLUX of DENS * '//trim(TRACER_NAME)//' (MYNN)', 'kg/m2/s', fill_halo=.true.)

    return
  end subroutine ATMOS_PHY_BL_MYNN_tendency_tracer

  !-----------------------------------------------------------------------------
  ! private routines
  !-----------------------------------------------------------------------------

!OCL SERIAL
  subroutine get_length( &
       KA, KS, KE_PBL, &
       PT0, q, n2,     &
       SFLX_PTV, l_mo, &
       FZ,             &
       l               )
    use scale_const, only: &
       GRAV   => CONST_GRAV, &
       KARMAN => CONST_KARMAN, &
       EPS    => CONST_EPS
    implicit none
    integer,  intent(in) :: KA, KS, KE_PBL

    real(RP), intent(in) :: PT0
    real(RP), intent(in) :: q(KA)
    real(RP), intent(in) :: n2(KA)
    real(RP), intent(in) :: SFLX_PTV  !> surface virtual temperature flux
    real(RP), intent(in) :: l_mo     !> Monin-Obukhov length
    real(RP), intent(in) :: FZ(0:KA)

    real(RP), intent(out) :: l(KA)

    real(RP) :: ls     !> L_S
    real(RP) :: lt     !> L_T
    real(RP) :: lb     !> L_B
    real(RP) :: rlm    !> 1/L_M
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
       z = ( FZ(k)+FZ(k-1) )*0.5_RP - FZ(KS-1)
       int_qz = int_qz + z * qdz
       int_q  = int_q + qdz
    end do
    ! LT
    lt = min( max(0.23_RP * int_qz / (int_q + EPS), &
                  LT_min), &
                  ATMOS_PHY_BL_MYNN_Lt_MAX )
    rlt = 1.0_RP / lt

    rlm = 1.0_RP / sign( max(abs(l_mo),EPS), l_mo )

    qc = ( GRAV / PT0 * max(SFLX_PTV,0.0_RP) * lt )**OneOverThree

    do k = KS, KE_PBL
       z = ( FZ(k)+FZ(k-1) )*0.5_RP - FZ(KS-1)
       zeta = z * rlm

       ! LS
       sw = sign(0.5_RP, zeta) + 0.5_RP ! 1 for zeta >= 0, 0 for zeta < 0
       ls = KARMAN * z &
          * ( sw / (1.0_RP + 2.7_RP*min(zeta,1.0_RP)*sw ) &
            + ( (1.0_RP - 100.0_RP*zeta)*(1.0_RP-sw) )**0.2_RP )

       ! LB
       sw  = sign(0.5_RP, n2(k)-EPS) + 0.5_RP ! 1 for dptdz >0, 0 for dptdz <= 0
       rn2sr = 1.0_RP / ( sqrt(n2(k)*sw) + 1.0_RP-sw)
       lb = (1.0_RP + 5.0_RP * sqrt(qc*rn2sr/lt)) * q(k) * rn2sr * sw & ! qc=0 when l_mo > 0
           +  1.E10_RP * (1.0_RP-sw)

       ! L
       l(k) = 1.0_RP / ( 1.0_RP/ls + rlt + 1.0_RP/lb )
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
    real(RP) :: q2

    integer :: k

    do k = KS, KE_PBL
       rf = min(0.5_RP / AF12 * ( Ri(k) &
                              + AF12*Rf1 &
                              - sqrt(Ri(k)**2 + 2.0_RP*AF12*(Rf1-2.0_RP*Rf2)*Ri(k) + (AF12*Rf1)**2) ), &
                Rfc)
       sh_2 = 3.0_RP * A2 * (G1+G2) * (Rfc-rf) / (1.0_RP-rf)
       sm_2 = sh_2 * AF12 * (Rf1-rf) / (Rf2-rf)
       q2 = B1 * l(k)**2 * sm_2 * (1.0_RP-rf) * dudz2(k)
       q2_2(k) = max( q2, 1.E-10_RP )
    end do

    return
  end subroutine get_q2_level2

!OCL SERIAL
  subroutine get_smsh( &
       KA, KS, KE_PBL, &
       q, ac,         &
       l, n2,         &
       potv, dudz2,   &
       tvsq, tvsq25,  &
       tltv, tltv25,  &
       qwtv, qwtv25,  &
       mynn_level3,   &
       sm25, smp,     &
       sh25, shpgh,   &
       gammat, gammaq )
    use scale_const, only: &
       EPS  => CONST_EPS, &
       GRAV => CONST_GRAV
    implicit none
    integer,  intent(in)  :: KA, KS, KE_PBL

    real(RP), intent(in)  :: q(KA)
    real(RP), intent(in)  :: ac(KA)
    real(RP), intent(in)  :: l(KA)
    real(RP), intent(in)  :: n2(KA)
    real(RP), intent(in)  :: potv(KA)
    real(RP), intent(in)  :: dudz2(KA)
    real(RP), intent(in)  :: tvsq(KA)
    real(RP), intent(in)  :: tvsq25(KA)
    real(RP), intent(in)  :: tltv(KA)
    real(RP), intent(in)  :: tltv25(KA)
    real(RP), intent(in)  :: qwtv(KA)
    real(RP), intent(in)  :: qwtv25(KA)
    logical,  intent(in)  :: mynn_level3

    real(RP), intent(out) :: sm25  (KA) ! S_M2.5
    real(RP), intent(out) :: smp   (KA) ! S_M'
    real(RP), intent(out) :: sh25  (KA) ! S_H2.5
    real(RP), intent(out) :: shpgh (KA) ! S_H' G_H * (q/L)^2
    real(RP), intent(out) :: gammat(KA) ! \Gamma_\theta
    real(RP), intent(out) :: gammaq(KA) ! \Gamma_q


    real(RP) :: l2q2 !> L^2/q^2
    real(RP) :: ac2  !> \alpha_c^2
    real(RP) :: p1   !> \Phi_1
    real(RP) :: p2   !> \Phi_2
    real(RP) :: p3   !> \Phi_3
    real(RP) :: p4   !> \Phi_4
    real(RP) :: p5   !> \Phi_5
    real(RP) :: rd25 !> 1/D_2.5
    real(RP) :: gh   !> G_H
    real(RP) :: gm   !> G_M
    real(RP) :: f1, f2, f3, f4

    ! for level 3
    real(RP) :: em   !> E_M * G_H
    real(RP) :: eh   !> E_H
    real(RP) :: ewgh !> E_w * G_H
    real(RP) :: rdp  !> 1/D'
    real(RP) :: cw25 !> Cw2.5
    real(RP) :: fact

    real(RP) :: lb, ub

    integer :: k

    ! level 2.5
    do k = KS, KE_PBL

       ac2 = ac(k)**2
       l2q2 = ( l(k) / max(q(k),EPS) )**2
       ! restriction: L/q <= 1/N
       if ( n2(k) > 0 ) l2q2 = min( l2q2, 1.0_RP/n2(k) )

       f1 = -  3.0_RP * ac2 * A2 * B2 * ( 1.0_RP - C3 )
       f2 = -  9.0_RP * ac2 * A1 * A2 * ( 1.0_RP - C2 )
       f3 =    9.0_RP * ac2 * A2**2   * ( 1.0_RP - C2 ) * ( 1.0_RP - C5 )
       f4 = - 12.0_RP * ac2 * A1 * A2 * ( 1.0_RP - C2 )
       f3 = f1 + f3
       f4 = f1 + f4

       gh = - n2(k) * l2q2
       gh = min( gh, - 1.0_RP / f4 ) ! for p1, p2, and p4 >= 0

       gm = max( dudz2(k) * l2q2, EPS )

       p1 = 1.0_RP + f1 * gh
       p2 = 1.0_RP + f2 * gh
       p3 = 1.0_RP + f3 * gh
       p4 = 1.0_RP + f4 * gh
       p5 = 6.0_RP * ac2 * A1**2 * gm
       rd25 = 1.0_RP / max( p2 * p4 + p5 * p3, EPS )

       sm25(k) = max( min( ac(k) * A1 * ( p3 - 3.0_RP * C1 * p4 ) * rd25, 0.44_RP / sqrt(gm) ), 0.0_RP )
       sh25(k) = max( min( ac(k) * A2 * ( p2 + 3.0_RP * C1 * p5 ) * rd25, 0.76_RP * B2 ), 0.0_RP )


       if ( mynn_level3 ) then ! level 3

          if ( sm25(k) < EPS ) then
             ! p3 - 3 * C1 * p4 = 0
             gh = ( 1.0_RP - 3.0_RP * C1 ) / ( 3.0_RP * C1 * f3 - f3 )
          else if ( sh25(k) < EPS ) then
              ! p2 + 3 * C1 * p5 = 0
             gh = - ( 3.0_RP * C1 * p5 ) / f2
          else
             ! sm25 / sh25 = A1 ( p3 - 3 * C1 * p4 ) / A2 / ( p2 + 3 * C1 * p5 )
             gh = ( A1 * sh25(k) * ( 1.0_RP - 3.0_RP * C1 ) - A2 * sm25(k) * ( 1.0_RP + 3.0_RP * C1 * p5 ) ) &
                / ( A2 * sm25(k) * f2 - A1 * sh25(k) * ( f3 - 3.0_RP * C1 * f4 ) )
          end if
          gh = sign( max( abs(gh), EPS ), gh )

          p1 = 1.0_RP + f1 * gh
          p2 = 1.0_RP + f2 * gh
          p3 = 1.0_RP + f3 * gh
          p4 = 1.0_RP + f4 * gh
          rd25 = 1.0_RP / max( p2 * p4 + p5 * p3, EPS )
          rdp  = 1.0_RP / max( p2 * ( p4 - p1 + 1.0_RP ) + p5 * ( p3 - p1 + 1.0_RP ), EPS )

          cw25 = p1 * ( p2 + 3.0_RP * C1 * p5 ) * rd25 / 3.0_RP

          ewgh = ( 1.0_RP - C3 ) * ( p2 * ( p1 - p4 ) + p5 * ( p1 - p3 ) ) * rdp
          ewgh  = sign( max(abs(ewgh),EPS), ewgh )
          fact = ( l2q2 * GRAV / ( l(k) * POTV(k) ) )**2 * ( tvsq(k) - tvsq25(k) )
          fact = fact * ewgh
          if ( gh >= 0.0_RP ) then
             lb = 0.12_RP - cw25
             ub = 0.76_RP - cw25
          else
             lb = 0.76_RP - cw25
             ub = 0.12_RP - cw25
          end if
          fact = min( max( fact, lb * gh ), ub * gh )
          fact = fact / ewgh

          em = 3.0_RP * ac(k) * A1 * ( 1.0_RP - C3 ) * ( p3 - p4 ) * rdp / gh
          eh = 3.0_RP * ac(k) * A2 * ( 1.0_RP - C3 ) * ( p2 + p5 ) * rdp

          smp  (k) = em * fact
          shpgh(k) = eh * fact / l2q2

          fact = tvsq(k) - tvsq25(k)
          if ( abs(fact) < EPS ) then
             fact = - eh * GRAV / ( q(k)**2 * POTV(k) )
          else
             fact = - shpgh(k) * POTV(k) / ( GRAV * fact )
          end if
          gammat(k) = ( tltv(k) - tltv25(k) ) * fact
          gammaq(k) = ( qwtv(k) - qwtv25(k) ) * fact

       else ! level 2.5
          smp  (k) = 0.0_RP
          shpgh(k) = 0.0_RP
       end if

    end do

    return
  end subroutine get_smsh

!OCL SERIAL
  subroutine calc_vertical_differece( &
       KA, KS, KE, &
       U, V, POTL, Qw,     &
       CDZ, FDZ, F2H,      &
       dudz2, dtldz, dqwdz )
    integer,  intent(in)  :: KA, KS, KE

    real(RP), intent(in) :: U   (KA)
    real(RP), intent(in) :: V   (KA)
    real(RP), intent(in) :: POTL(KA)
    real(RP), intent(in) :: Qw  (KA)
    real(RP), intent(in) :: CDZ (KA)
    real(RP), intent(in) :: FDZ (KA)
    real(RP), intent(in) :: F2H (KA,2)

    real(RP), intent(out) :: dudz2(KA)
    real(RP), intent(out) :: dtldz(KA)
    real(RP), intent(out) :: dqwdz(KA)

    real(RP) :: Uh(KA), Vh(KA)
    real(RP) :: qh(KA)

    integer :: k

    do k = KS, KE
       Uh(k) = f2h(k,1) * U(k+1) + f2h(k,2) * U(k)
       Vh(k) = f2h(k,1) * V(k+1) + f2h(k,2) * V(k)
    end do

    dudz2(KS) = ( ( Uh(KS) - U(KS) )**2 + ( Vh(KS) - V(KS) )**2 ) / ( CDZ(KS) * 0.5_RP )**2
!    dudz2(KS) = ( ( Uh(KS) )**2 + ( Vh(KS) )**2 ) / CDZ(KS)**2 &
    do k = KS+1, KE
       dudz2(k) = ( ( Uh(k) - Uh(k-1) )**2 + ( Vh(k) - Vh(k-1) )**2 ) / CDZ(k)**2
    end do


    do k = KS, KE
       qh(k) = f2h(k,1) * POTL(k+1) + f2h(k,2) * POTL(k)
    end do
    dtldz(KS) = ( qh(KS) - POTL(KS) ) / ( CDZ(KS) * 0.5_RP )
!    dtldz(KS) = ( POTL(KS+1) - POTL(KS) ) / FDZ(KS)
    do k = KS+1, KE
       dtldz(k) = ( qh(k) - qh(k-1) ) / CDZ(k)
    end do

    do k = KS, KE
       qh(k) = f2h(k,1) * Qw(k+1) + f2h(k,2) * Qw(k)
    end do
    dqwdz(KS) = ( qh(KS) - Qw(KS) ) / ( CDZ(KS) * 0.5_RP )
!    dqwdz(KS) = ( Qw(KS+1) - Qw(KS) ) / FDZ(KS)
    do k = KS+1, KE
       dqwdz(k) = ( qh(k) - qh(k-1) ) / CDZ(k)
    end do

    return
  end subroutine calc_vertical_differece

!OCL SERIAL
  subroutine partial_condensation( &
       KA, KS, KE, &
       PRES, POTT, POTL, Qw, Qdry, EXNER, &
       dtldz, dqwdz, l, sh25, ac,         &
       tsq, qsq, cov,                     &
       mynn_level3,                       &
       betat, betaq                       )
    use scale_const, only: &
       CPdry   => CONST_CPdry,  &
       Rvap    => CONST_Rvap,   &
       CPvap   => CONST_CPvap,  &
       EPSvap  => CONST_EPSvap, &
       EPSTvap => CONST_EPSTvap
    use scale_atmos_saturation, only: &
       ATMOS_SATURATION_pres2qsat => ATMOS_SATURATION_pres2qsat_liq
    use scale_atmos_hydrometeor, only: &
       HYDROMETEOR_LHV => ATMOS_HYDROMETEOR_LHV
    integer,  intent(in)  :: KA, KS, KE
    real(RP), intent(in)  :: PRES (KA)
    real(RP), intent(in)  :: POTT (KA)
    real(RP), intent(in)  :: POTL (KA)
    real(RP), intent(in)  :: Qw   (KA)
    real(RP), intent(in)  :: Qdry (KA)
    real(RP), intent(in)  :: EXNER(KA)
    real(RP), intent(in)  :: dtldz(KA)
    real(RP), intent(in)  :: dqwdz(KA)
    real(RP), intent(in)  :: l    (KA)
    real(RP), intent(in)  :: sh25 (KA)
    real(RP), intent(in)  :: ac   (KA)
    real(RP), intent(in)  :: tsq  (KA)
    real(RP), intent(in)  :: qsq  (KA)
    real(RP), intent(in)  :: cov  (KA)
    logical,  intent(in)  :: mynn_level3

    real(RP), intent(out) :: betat(KA)
    real(RP), intent(out) :: betaq(KA)

    real(RP) :: TEML(KA)
    real(RP) :: Qsl (KA)
    real(RP) :: LHVL(KA)
    real(RP) :: dQsl
    real(RP) :: CPtot
    real(RP) :: aa, bb, cc
    real(RP) :: sigma_s
    real(RP) :: Q1, Qlp
    real(RP) :: RR, Rt

    integer :: k

    do k = KS, KE
       TEML(k) = POTL(k) * EXNER(k)
    end do

    call ATMOS_SATURATION_pres2qsat( &
         KA, KS, KE, &
         TEML(:), PRES(:), Qdry(:), & ! (in)
         Qsl(:)                     ) ! (out)

    call HYDROMETEOR_LHV( &
         KA, KS, KE, & ! (in)
         TEML(:), & ! (in)
         LHVL(:)  ) ! (out)

    do k = KS, KE

       dQsl = ( 1.0_RP + Qsl(k) / ( Qdry(k) * EPSvap ) ) * Qsl(k) * LHVL(k) / ( Rvap * TEML(k)**2 )
       CPtot = Qdry(k) * CPdry + Qsl(k) * CPvap
       aa = 1.0_RP / ( 1.0_RP + LHVL(k)/CPtot * dQsl )
       bb = EXNER(k) * dQsl

       if ( mynn_level3 ) then
          sigma_s = 0.5_RP * aa * sqrt( max( qsq(k) - 2.0_RP * bb * cov(k) + bb**2 * tsq(k), 1.0e-20_RP ) )
       else
          ! level 2.5
          sigma_s = max( 0.5_RP * aa * l(k) * sqrt( ac(k) * B2 * sh25(k) ) * abs( dqwdz(k) - bb * dtldz(k) ), &
                         1.0e-10_RP )
       end if

       Q1 = aa * ( Qw(k) - Qsl(k) ) * 0.5_RP / sigma_s
       RR = min( max( 0.5_RP * ( 1.0_RP + erf(Q1*rsqrt_2) ), 0.0_RP ), 1.0_RP )
       Qlp = min( max( 2.0_RP * sigma_s * ( RR * Q1 + rsqrt_2pi &
#if defined(PGI) || defined(SX)
               * exp( -min( 0.5_RP*Q1**2, 1.E+3_RP ) ) & ! apply exp limiter
#else
               * exp(-0.5_RP*Q1**2) &
#endif
               ), 0.0_RP ), Qw(k) * 0.5_RP )
       cc = ( 1.0_RP + EPSTvap * Qw(k) - Qlp / EPSvap ) / EXNER(k) * LHVL(k) / CPtot &
             - POTT(k) / EPSvap
       Rt = min( max( RR - Qlp / (2.0_RP*sigma_s*sqrt_2pi) &
#if defined(PGI) || defined(SX)
               * exp( -min( 0.5_RP*Q1**2, 1.E+3_RP ) ) & ! apply exp limiter
#else
               * exp(-Q1**2 * 0.5_RP) &
#endif
               , 0.0_RP ), 1.0_RP )
       betat(k) = 1.0_RP + EPSTvap * Qw(k) - Qlp / EPSvap - Rt * aa * bb * cc
       betaq(k) = EPSTvap * POTT(k) + Rt * aa * cc

    end do

    return
  end subroutine partial_condensation

!OCL SERIAL
  subroutine get_f2h( &
       KA, KS, KE, &
       CDZ, &
       f2h )
    integer,  intent(in)  :: KA, KS, KE
    real(RP), intent(in)  :: CDZ(KA)
    real(RP), intent(out) :: f2h(KA,2)

    real(RP) :: dz1, dz2
    integer :: k

    do k = KS, KE
       dz1 = CDZ(k+1)
       dz2 = CDZ(k)
       f2h(k,1) = dz2 / ( dz1 + dz2 )
       f2h(k,2) = dz1 / ( dz1 + dz2 )
    end do

    return
  end subroutine get_f2h

end module scale_atmos_phy_bl_mynn
