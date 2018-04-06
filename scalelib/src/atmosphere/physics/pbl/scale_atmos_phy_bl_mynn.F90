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
  use scale_stdio
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
  public :: ATMOS_PHY_BL_mynn_setup
  public :: ATMOS_PHY_BL_mynn_tendency
  public :: ATMOS_PHY_BL_mynn_tendency_tracer

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  integer,                parameter, public :: ATMOS_PHY_BL_MYNN_NTRACER = 1
  character(len=H_SHORT), parameter, public :: ATMOS_PHY_BL_MYNN_NAME(1) = &
       (/ 'TKE_MYNN' /)
  character(len=H_LONG),  parameter, public :: ATMOS_PHY_BL_MYNN_DESC(1) = &
       (/ 'turbulent kinetic energy (MYNN)' /)
  character(len=H_SHORT), parameter, public :: ATMOS_PHY_BL_MYNN_UNITS(1) = &
       (/ 'm2/s2/kg3' /)

  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  real(RP), private, parameter :: OneOverThree = 1.0_RP / 3.0_RP
  real(RP), private, parameter :: LT_min       = 1.E-6_RP

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

  real(RP), private            :: ATMOS_PHY_BL_MYNN_TKE_MIN  =  1.E-10_RP
  real(RP), private            :: ATMOS_PHY_BL_MYNN_N2_MAX   =   1.E+3_RP
  real(RP), private            :: ATMOS_PHY_BL_MYNN_NU_MIN   =   1.E-6_RP
  real(RP), private            :: ATMOS_PHY_BL_MYNN_NU_MAX   = 10000.0_RP
  real(RP), private            :: ATMOS_PHY_BL_MYNN_KH_MIN   =   1.E-6_RP
  real(RP), private            :: ATMOS_PHY_BL_MYNN_KH_MAX   = 10000.0_RP
  real(RP), private            :: ATMOS_PHY_BL_MYNN_Lt_MAX   =   700.0_RP ! ~ 0.23 * 3 km

  real(RP), private            :: ATMOS_PHY_BL_MYNN_LEVEL    = 2.5_RP

  !-----------------------------------------------------------------------------
contains
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

    real(RP) :: ATMOS_PHY_BL_MYNN_PBL_MAX = 1.E+10_RP !> maximum height of the PBL

    NAMELIST / PARAM_ATMOS_PHY_BL_MYNN / &
       ATMOS_PHY_BL_MYNN_PBL_MAX,  &
       ATMOS_PHY_BL_MYNN_N2_MAX,   &
       ATMOS_PHY_BL_MYNN_NU_MIN,   &
       ATMOS_PHY_BL_MYNN_NU_MAX,   &
       ATMOS_PHY_BL_MYNN_KH_MIN,   &
       ATMOS_PHY_BL_MYNN_KH_MAX,   &
       ATMOS_PHY_BL_MYNN_Lt_MAX,   &
       ATMOS_PHY_BL_MYNN_LEVEL

    integer :: ierr
    integer :: k, i, j
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_PROGRESS(*) 'Module[pbl mynn] / Categ[atmosphere physics] / Origin[SCALE lib]'
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
    if( IO_NML ) write(IO_FID_NML,nml=PARAM_ATMOS_PHY_BL_MYNN)

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

    return
  end subroutine ATMOS_PHY_BL_MYNN_setup

  !-----------------------------------------------------------------------------
  !> ATMOS_PHY_BL_MYNN_tendency
  !! calculate tendency by the virtical eddy viscosity
  !<
  subroutine ATMOS_PHY_BL_MYNN_tendency( &
       KA, KS, KE, IA, IS, IE, JA, JS, JE, &
       DENS, U, V, POTT, TKE,           &
       PRES, EXNER, N2,                 &
       QDRY, QV, Qw, POTL, POTV,        &
       SFLX_MU, SFLX_MV, SFLX_SH, l_mo, &
       CZ, FZ, dt,                      &
       RHOU_t, RHOV_t, RHOT_t, RTKE_t,  &
       Nu, Kh                           )
    use scale_const, only: &
       GRAV    => CONST_GRAV,  &
       Rdry    => CONST_Rdry,  &
       Rvap    => CONST_Rvap,  &
       CPdry   => CONST_CPdry, &
       CPvap   => CONST_CPvap, &
       EPSTvap => CONST_EPSTvap
    use scale_prc, only: &
       PRC_abort
    use scale_atmos_hydrometeor, only: & !! TODO
       HYDROMETEOR_LHV => ATMOS_HYDROMETEOR_LHV, &
       LHV
    use scale_matrix, only: &
       MATRIX_SOLVER_tridiagonal
    use scale_atmos_saturation, only: &
       ATMOS_SATURATION_dens2qsat => ATMOS_SATURATION_dens2qsat_all
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
    real(RP), intent(in) :: TKE    (KA,IA,JA) !> sub-grid turbulent kinetic energy
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
    real(RP), intent(in) :: l_mo   (   IA,JA) !> Monin-Obukhov length

    real(RP), intent(in)  :: CZ(  KA,IA,JA)
    real(RP), intent(in)  :: FZ(0:KA,IA,JA)
    real(DP), intent(in)  :: dt

    real(RP), intent(out) :: RHOU_t(KA,IA,JA) !> tendency of dens * u
    real(RP), intent(out) :: RHOV_t(KA,IA,JA) !> tendency of dens * v
    real(RP), intent(out) :: RHOT_t(KA,IA,JA) !> tendency of dens * pt
    real(RP), intent(out) :: RTKE_t(KA,IA,JA) !> tenddency of dens * tke
    real(RP), intent(out) :: Nu    (KA,IA,JA) !> eddy viscosity coefficient
    real(RP), intent(out) :: Kh    (KA,IA,JA) !> eddy diffusion coefficient

    real(RP) :: Ri   (KA,IA,JA) !> Richardson number
    real(RP) :: Pr   (KA,IA,JA) !> Plandtle number
    real(RP) :: prod (KA,IA,JA) !> TKE production term
    real(RP) :: diss (KA,IA,JA) !> TKE dissipation term
    real(RP) :: dudz2(KA,IA,JA) !> (du/dz)^2 + (dv/dz)^2
    real(RP) :: l    (KA,IA,JA) !> length scale L

    real(RP) :: flxU(KA,IA,JA) !> dens * w * u
    real(RP) :: flxV(KA,IA,JA) !> dens * w * v
    real(RP) :: flxT(KA,IA,JA) !> dens * w * pt

    real(RP) :: TEML  (KA) !> liquid water temperature
    real(RP) :: RHONu (KA) !> dens * Nu at the half level
    real(RP) :: RHOKh (KA) !> dens * Kh at the half level
    real(RP) :: LHVL  (KA) !> latent heat
    real(RP) :: CPtot      !> specific heat
    real(RP) :: N2_new(KA) !> squared Brunt-Baisala frequency
    real(RP) :: SFLX_PT    !> surface potential temperature flux
    real(RP) :: sm    (KA) !> stability function for velocity
    real(RP) :: sh    (KA) !> stability function for scalars
    real(RP) :: q     (KA) !> q
    real(RP) :: q2_2  (KA) !> q^2 for level 2

    real(RP) :: qlp        !> liquid water
    real(RP) :: ac         !> \alpha_c
    real(RP) :: Q1
    real(RP) :: Qsl(KA)
    real(RP) :: dQsl
    real(RP) :: sigma_s
    real(RP) :: RR
    real(RP) :: Rt
    real(RP) :: betat
    real(RP) :: betaq
    real(RP) :: aa
    real(RP) :: bb
    real(RP) :: cc

    real(RP) :: a(KA)
    real(RP) :: b(KA)
    real(RP) :: c(KA)
    real(RP) :: d(KA)
    real(RP) :: ap
    real(RP) :: phi_N(KA)
    real(RP) :: tke_P

    real(RP) :: sf_t

    real(RP) :: f2h(KA,2)

    real(RP) :: mynn_level

    integer :: k, i, j
    !---------------------------------------------------------------------------


    LOG_PROGRESS(*) "atmosphere / physics / pbl / MYNN"

    mynn_level = ATMOS_PHY_BL_MYNN_LEVEL
    if ( mynn_level .ne. 2.5_RP ) then
       LOG_ERROR("ATMOS_PHY_BL_MYNN_tendency",*) 'only level 2.5 is supported at this moment'
       call PRC_abort
    end if

!OCL INDEPENDENT
    !$omp parallel do default(none) &
    !$omp OMP_SCHEDULE_ collapse(2) &
    !$omp shared(KA,KS,KE_PBL,KE,IS,IE,JS,JE, &
    !$omp        GRAV,CPdry,LHV,EPSTvap,UNDEF,RSQRT_2,SQRT_2PI,RSQRT_2PI, &
    !$omp        ATMOS_PHY_BL_MYNN_N2_MAX,ATMOS_PHY_BL_MYNN_TKE_MIN, &
    !$omp        ATMOS_PHY_BL_MYNN_NU_MIN,ATMOS_PHY_BL_MYNN_NU_MAX, &
    !$omp        ATMOS_PHY_BL_MYNN_KH_MIN,ATMOS_PHY_BL_MYNN_KH_MAX, &
    !$omp        RHOU_t,RHOV_t,RHOT_t,RTKE_t,Nu,Kh, &
    !$omp        DENS,TKE,U,V,POTT,PRES,QDRY,QV,Qw,POTV,POTL,EXNER,N2,SFLX_MU,SFLX_MV,SFLX_SH,l_mo, &
    !$omp        CZ,FZ,dt, &
    !$omp        Ri,Pr,prod,diss,dudz2,l,flxU,flxV,flxT) &
    !$omp private(N2_new,sm,sh,q,q2_2,SFLX_PT,TEML,RHONu,RHOKh,LHVL,CPtot,qlp,ac, &
    !$omp         Q1,Qsl,dQsl,sigma_s,RR,Rt,betat,betaq,aa,bb,cc, &
    !$omp         a,b,c,d,ap,phi_n,tke_P,sf_t,f2h, &
    !$omp         k,i,j)
    do j = JS, JE
    do i = IS, IE
       SFLX_PT = SFLX_SH(i,j) / ( CPdry * DENS(KS,i,j) * EXNER(KS,i,j) )

       dudz2(KS,i,j) = ( ( U(KS+1,i,j) - U(KS,i,j) )**2 + ( V(KS+1,i,j) - V(KS,i,j) )**2 ) &
                     / ( CZ(KS+1,i,j) - CZ(KS,i,j) )**2
       do k = KS+1, KE_PBL
          dudz2(k,i,j) = ( ( U(k+1,i,j) - U(k-1,i,j) )**2 + ( V(k+1,i,j) - V(k-1,i,j) )**2 ) &
                       / ( CZ(k+1,i,j) - CZ(k-1,i,j) )**2
       end do
       do k = KE_PBL+1, KE
          dudz2(k,i,j) = UNDEF
       end do

       do k = KS, KE_PBL
          n2_new(k) = min( ATMOS_PHY_BL_MYNN_N2_MAX, N2(k,i,j) )
          Ri(k,i,j) = n2_new(k) / max(dudz2(k,i,j), 1E-10_RP)
       end do
       do k = KE_PBL+1, KE
          Ri(k,i,j) = UNDEF
       end do

       do k = KS, KE_PBL
          q(k) = sqrt( max(TKE(k,i,j), ATMOS_PHY_BL_MYNN_TKE_MIN)*2.0_RP )
       end do


       call get_f2h( &
            KA, KS, KE, &
            FZ(:,i,j), & ! (in)
            f2h(:,:)   ) ! (out)

       ! length
       call get_length( &
            KA, KS, KE_PBL, &
            POTT(:,i,j), q(:), n2_new(:), & ! (in)
            SFLX_PT, l_mo(i,j),           & ! (in)
            FZ(:,i,j),                    & ! (in)
            l(:,i,j)                      ) ! (out)

       call get_q2_level2( &
            KA, KS, KE_PBL,                   & ! (in)
            q2_2(:),                          & ! (out)
            dudz2(:,i,j), Ri(:,i,j), l(:,i,j) ) ! (in)

       call get_smsh( &
            KA, KS, KE_PBL,                   & ! (in)
            sm(:), sh(:),                     & ! (out)
            q(:), q2_2(:),                    & ! (in)
            l(:,i,j), n2_new(:), dudz2(:,i,j) ) ! (in)


       ! liquid water temperature
       do k = KS, KE_PBL+1
          TEML(k) = POTL(k,i,j) * EXNER(k,i,j)
       end do

       call ATMOS_SATURATION_dens2qsat( &
            KA, KS, KE_PBL, &
            TEML(:), DENS(:,i,j), & ! (in)
            Qsl(:)                ) ! (out)

       call HYDROMETEOR_LHV( &
            KA, KS, KE_PBL, & ! (in)
            TEML(:), & ! (in)
            LHVL(:)  ) ! (out)

       do k = KS+1, KE_PBL

          dQsl = ( Qsl(k) * LHVL(k) / ( Rvap * TEML(k) ) - Qsl(k) ) / TEML(k)
          CPtot = Qdry(k,i,j) * CPdry + Qsl(k) * CPvap
          aa = 1.0_RP / ( 1.0_RP + LHVL(k)/CPtot * dQsl )
          bb = EXNER(k,i,j) * dQsl
          ac = min( q(k)/sqrt(q2_2(k)), 1.0_RP )
          sigma_s = max( sqrt( 0.25_RP * aa**2 * l(k,i,j)**2 * ac * B2 * sh(k) ) &
                       * abs( Qw(k+1,i,j) - Qw(k-1,i,j) - bb * ( POTL(k+1,i,j)-POTL(k-1,i,j) ) ) &
                       / ( CZ(k+1,i,j) - CZ(k-1,i,j) ), &
                       1.0e-10_RP )
          Q1 = aa * ( Qw(k,i,j) - Qsl(k) ) * 0.5_RP / sigma_s
          RR = min( max( 0.5_RP * ( 1.0_RP + erf(Q1*rsqrt_2) ), 0.0_RP ), 1.0_RP )
          Qlp = min( max( 2.0_RP * sigma_s * ( RR * Q1 + rsqrt_2pi &
#if defined(PGI) || defined(SX)
               * exp( -min( 0.5_RP*Q1**2, 1.E+3_RP ) ) & ! apply exp limiter
#else
               * exp(-0.5_RP*Q1**2) &
#endif
               ), 0.0_RP ), Qw(k,i,j) * 0.5_RP )
          cc = ( 1.0_RP + EPSTvap * Qw(k,i,j) - (1.0_RP+EPSTvap) * Qlp ) / EXNER(k,i,j) * LHVL(k) / CPtot &
             - (1.0_RP+EPSTvap) * POTT(k,i,j)
          Rt = min( max( RR - Qlp / (2.0_RP*sigma_s*sqrt_2pi) &
#if defined(PGI) || defined(SX)
               * exp( -min( 0.5_RP*Q1**2, 1.E+3_RP ) ) & ! apply exp limiter
#else
               * exp(-Q1**2 * 0.5_RP) &
#endif
               , 0.0_RP ), 1.0_RP )
          betat = 1.0_RP + EPSTvap * Qw(k,i,j) - (1.0_RP+EPSTvap) * Qlp - Rt * aa * bb * cc
          betaq = EPSTvap * POTT(k,i,j) + Rt * aa * cc
          n2_new(k) = min(ATMOS_PHY_BL_MYNN_N2_MAX, &
                          GRAV * ( ( POTL(k+1,i,j) - POTL(k-1,i,j) ) * betat &
                                 + ( Qw  (k+1,i,j) - Qw  (k-1,i,j) ) * betaq ) &
                               / ( ( CZ(k+1,i,j) - CZ(k-1,i,j) ) * POTV(k,i,j) ) )
       end do
       n2_new(KS) = n2_new(KS+1)

       do k = KS, KE_PBL
          Ri(k,i,j) = n2_new(k) / max(dudz2(k,i,j), 1E-10_RP)
       end do
       do k = KE_PBL+1, KE
          Ri    (k,i,j) = 0.0_RP
          n2_new(k    ) = 0.0_RP
       end do


       ! length
       call get_length( &
            KA, KS, KE_PBL, &
            POTT(:,i,j), q(:), n2_new(:), & ! (in)
            SFLX_PT, l_mo(i,j),           & ! (in)
            FZ(:,i,j),                    & ! (in)
            l(:,i,j)                      ) ! (out)

       call get_q2_level2( &
            KA, KS, KE_PBL,                   & ! (in)
            q2_2(:),                          & ! (out)
            dudz2(:,i,j), Ri(:,i,j), l(:,i,j) ) ! (in)

       call get_smsh( &
            KA, KS, KE_PBL,                   & ! (in)
            sm(:), sh(:),                     & ! (out)
            q(:), q2_2(:),                    & ! (in)
            l(:,i,j), n2_new(:), dudz2(:,i,j) ) ! (in)


       do k = KS, KE_PBL
          Nu(k,i,j) = max( min( l(k,i,j) * q(k) * sm(k), &
                           ATMOS_PHY_BL_MYNN_NU_MAX ), &
                           ATMOS_PHY_BL_MYNN_NU_MIN )
          Kh(k,i,j) = MAX( min( l(k,i,j) * q(k) * sh(k), &
                           ATMOS_PHY_BL_MYNN_KH_MAX ), &
                           ATMOS_PHY_BL_MYNN_KH_MIN )
          Pr(k,i,j) = Nu(k,i,j) / Kh(k,i,j)
       end do
       do k = KE_PBL+1, KE
          Nu(k,i,j) = 0.0_RP
          Kh(k,i,j) = 0.0_RP
          Pr(k,i,j) = 1.0_RP
       end do

       ! dens * coefficient at the half level
       do k = KS, KE_PBL-1
          RHONu(k) = f2h(k,1) * DENS(k+1,i,j) * Nu(k+1,i,j) &
                   + f2h(k,2) * DENS(k  ,i,j) * Nu(k  ,i,j)
          RHOKh(k) = f2h(k,1) * DENS(k+1,i,j) * Kh(k+1,i,j) &
                   + f2h(k,2) * DENS(k  ,i,j) * Kh(k  ,i,j)
       end do

       ! time integration

       ! dens * u
       sf_t = SFLX_MU(i,j) / ( FZ(KS,i,j) - FZ(KS-1,i,j) )
       d(KS) = U(KS,i,j) + dt * sf_t / DENS(KS,i,j)
       do k = KS+1, KE_PBL
          d(k) = U(k,i,j)
       end do
       c(KS) = 0.0_RP
       do k = KS, KE_PBL-1
          ap = - dt * RHONu(k) / ( CZ(k+1,i,j) - CZ(k,i,j) )
          a(k) = ap / ( DENS(k,i,j) * ( FZ(k,i,j) - FZ(k-1,i,j) ) )
          b(k) = - a(k) - c(k) + 1.0_RP
          c(k+1) = ap / ( DENS(k+1,i,j) * ( FZ(k+1,i,j) - FZ(k,i,j) ) )
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
          flxU(k,i,j) = - RHONu(k) * ( phi_n(k+1) - phi_n(k) ) / ( CZ(k+1,i,j) - CZ(k,i,j) )
       end do
       do k = KE_PBL, KE
          flxU(k,i,j) = 0.0_RP
       end do


       ! dens * v
       sf_t = SFLX_MV(i,j) / ( FZ(KS,i,j) - FZ(KS-1,i,j) )
       d(KS) = V(KS,i,j) + dt * sf_t / DENS(KS,i,j)
       do k = KS+1, KE_PBL
          d(k) = V(k,i,j)
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
          flxV(k,i,j) = - RHONu(k) * ( phi_n(k+1) - phi_n(k) ) / ( CZ(k+1,i,j) - CZ(k,i,j) )
       end do
       do k = KE_PBL, KE
          flxV(k,i,j) = 0.0_RP
       end do


       ! dens * pott
       sf_t = SFLX_PT / ( FZ(KS,i,j) - FZ(KS-1,i,j) )
       d(KS) = POTT(KS,i,j) + dt * sf_t
       do k = KS+1, KE_PBL
          d(k) = POTT(k,i,j)
       end do

       c(KS) = 0.0_RP
       do k = KS, KE_PBL-1
          ap = - dt * RHOKh(k) / ( CZ(k+1,i,j) - CZ(k,i,j) )
          a(k) = ap / ( DENS(k,i,j) * ( FZ(k,i,j) - FZ(k-1,i,j) ) )
          b(k) = - a(k) - c(k) + 1.0_RP
          c(k+1) = ap / ( DENS(k+1,i,j) * ( FZ(k+1,i,j) - FZ(k,i,j) ) )
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
          flxT(k,i,j) = - RHOKh(k) * ( phi_n(k+1) - phi_n(k) ) / ( CZ(k+1,i,j) - CZ(k,i,j) )
       end do
       do k = KE_PBL, KE
          flxT(k,i,j) = 0.0_RP
       end do


       ! dens * TKE
       do k = KS, KE_PBL
          prod(k,i,j) = Nu(k,i,j) * dudz2(k,i,j) - Kh(k,i,j) * n2_new(k)
          tke_p = q(k)**2 * 0.5_RP
          d(k) = tke_p + dt * prod(k,i,j)
       end do
       do k = KE_PBL+1, KE
          prod(k,i,j) = 0.0_RP
       end do
       c(KS) = 0.0_RP
       do k = KS, KE_PBL-1
          ap = - dt * 3.0_RP * RHONu(k) / ( CZ(k+1,i,j) - CZ(k,i,j) )
          a(k) = ap / ( DENS(k,i,j) * ( FZ(k,i,j) - FZ(k-1,i,j) ) )
          diss(k,i,j) = 2.0_RP * q(k) / ( B1 * l(k,i,j) )
          b(k) = - a(k) - c(k) + 1.0_RP + diss(k,i,j) * dt
          c(k+1) = ap / ( DENS(k+1,i,j) * ( FZ(k+1,i,j) - FZ(k,i,j) ) )
       end do
       a(KE_PBL) = 0.0_RP
       diss(KE_PBL,i,j) = 2.0_RP * q(KE_PBL) / ( B1 * l(KE_PBL,i,j) )
       b(KE_PBL) = - c(KE_PBL) + 1.0_RP + diss(KE_PBL,i,j) * dt
       do k = KE_PBL+1, KE
          diss(k,i,j) = 0.0_RP
       end do

       call MATRIX_SOLVER_tridiagonal( &
            KA, KS, KE_PBL, &
            a(:), b(:), c(:), d(:), & ! (in)
            phi_n(:)                ) ! (out)

       do k = KS, KE_PBL
          RTKE_t(k,i,j) = ( max(phi_n(k), ATMOS_PHY_BL_MYNN_TKE_MIN) - TKE(k,i,j) ) * DENS(k,i,j) / dt
       end do
       do k = KE_PBL+1, KE
          RTKE_t(k,i,j) = 0.0_RP
       end do

    end do
    end do

    l(KE_PBL+1:KE,:,:) = 0.0_RP
    call FILE_HISTORY_in(Ri(:,:,:), 'Ri_MYNN', 'Richardson number', '1',     fill_halo=.true. )
    call FILE_HISTORY_in(Pr(:,:,:), 'Pr_MYNN', 'Prandtl number',    '1',     fill_halo=.true. )
    call FILE_HISTORY_in(prod(:,:,:), 'TKE_prod_MYNN', 'TKE production',  'm2/s3', fill_halo=.true.)
    call FILE_HISTORY_in(diss(:,:,:), 'TKE_diss_MYNN', 'TKE dissipation', 'm2/s3', fill_halo=.true.)
    call FILE_HISTORY_in(dudz2(:,:,:), 'dUdZ2_MYNN', 'dudz2', 'm2/s2', fill_halo=.true.)
    call FILE_HISTORY_in(l(:,:,:), 'L_mix_MYNN', 'minxing length', 'm', fill_halo=.true.)

    call FILE_HISTORY_in(flxU(:,:,:), 'ZFLX_RHOU_MYNN', 'Z FLUX of RHOU (MYNN)', 'kg/m/s2', fill_halo=.true.)
    call FILE_HISTORY_in(flxV(:,:,:), 'ZFLX_RHOV_MYNN', 'Z FLUX of RHOV (MYNN)', 'kg/m/s2', fill_halo=.true.)
    call FILE_HISTORY_in(flxT(:,:,:), 'ZFLX_RHOT_MYNN', 'Z FLUX of RHOT (MYNN)', 'K kg/m2/s', fill_halo=.true.)

    return
  end subroutine ATMOS_PHY_BL_MYNN_tendency

  !-----------------------------------------------------------------------------
  !> ATMOS_PHY_BL_MYNN_tendency_tracer
  !! calculate tendency of tracers by the eddy viscosity
  !<
  subroutine ATMOS_PHY_BL_MYNN_tendency_tracer( &
       KA, KS, KE, IA, IS, IE, JA, JS, JE, &
       DENS, QTRC, SFLX_Q, Kh, &
       CZ, FZ, DT, name,       &
       RHOQ_t                  )
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
    real(RP),         intent(in) :: Kh    (KA,IA,JA) !> eddy diffusion coefficient
    real(RP),         intent(in) :: CZ (  KA,IA,JA) !> z at the full level
    real(RP),         intent(in) :: FZ (0:KA,IA,JA) !> z at the half level
    real(DP),         intent(in) :: DT              !> time step
    character(len=*), intent(in) :: name            !> name of tracer (for history output)

    real(RP), intent(out) :: RHOQ_t(KA,IA,JA) !> tendency of tracers

    real(RP) :: QTRC_n(KA) !> value at the next time step
    real(RP) :: RHOKh(KA)
    real(RP) :: a(KA)
    real(RP) :: b(KA)
    real(RP) :: c(KA)
    real(RP) :: d(KA)
    real(RP) :: ap
    real(RP) :: sf_t

    real(RP) :: flx(KA,IA,JA)

    real(RP) :: f2h(KA,2) !> coefficient to convert from full to half level

    integer :: k, i, j

!OCL INDEPENDENT
    !$omp parallel do default(none) OMP_SCHEDULE_ collapse(2) &
    !$omp shared(KA,KS,KE_PBL,KE,IS,IE,JS,JE) &
    !$omp shared(RHOQ_t,DENS,QTRC,SFLX_Q,Kh,CZ,FZ,F2H,DT,flx) &
    !$omp private(QTRC_n,RHOKh,a,b,c,d,ap,sf_t) &
    !$omp private(k,i,j)
    do j = JS, JE
    do i = IS, IE

       call get_f2h( &
            KA, KS, KE, &
            FZ(:,i,j), & ! (in)
            f2h(:,:)   ) ! (out)

       ! dens * coefficient at the half level
       do k = KS, KE_PBL-1
          RHOKh(k) = f2h(k,1) * DENS(k+1,i,j) * Kh(k+1,i,j) &
                   + f2h(k,2) * DENS(k  ,i,j) * Kh(k  ,i,j)
       end do

       sf_t = SFLX_Q(i,j) / ( FZ(KS,i,j) - FZ(KS-1,i,j) )
       d(KS) = QTRC(KS,i,j) + dt * sf_t / DENS(KS,i,j)
       do k = KS+1, KE_PBL
          d(k) = QTRC(k,i,j)
       end do

       c(KS) = 0.0_RP
       do k = KS, KE_PBL-1
          ap = - dt * RHOKh(k) / ( CZ(k+1,i,j) - CZ(k,i,j) )
          a(k) = ap / ( DENS(k,i,j) * ( FZ(k,i,j) - FZ(k-1,i,j) ) )
          b(k) = - a(k) - c(k) + 1.0_RP
          c(k+1) = ap / ( DENS(k+1,i,j) * ( FZ(k+1,i,j) - FZ(k,i,j) ) )
       end do
       a(KE_PBL) = 0.0_RP
       b(KE_PBL) = - c(KE_PBL) + 1.0_RP

       call MATRIX_SOLVER_tridiagonal( &
               KA, KS, KE_PBL, &
               a(:), b(:), c(:), d(:), & ! (in)
               QTRC_n(:)               ) ! (out)

       RHOQ_t(KS,i,j) = ( QTRC_n(KS) - QTRC(KS,i,j) ) * DENS(KS,i,j) / dt - sf_t
       do k = KS+1, KE_PBL
          RHOQ_t(k,i,j) = ( QTRC_n(k) - QTRC(k,i,j) ) * DENS(k,i,j) / dt
       end do
       do k = KE_PBL+1, KE
          RHOQ_t(k,i,j) = 0.0_RP
       end do

       do k = KS, KE_PBL-1
          flx(k,i,j) = - RHOKh(k) * ( QTRC_n(k+1) - QTRC_n(k) ) / ( CZ(k+1,i,j) - CZ(k,i,j) )
       end do
       do k = KE_PBL, KE
          flx(k,i,j) = 0.0_RP
       end do

    end do
    end do

    call FILE_HISTORY_in(flx(:,:,:), 'ZFLX_'//trim(name)//'_MYNN', 'Z FLUX of DENS * '//trim(name)//' (MYNN)', 'kg/m2/s', fill_halo=.true.)

    return
  end subroutine ATMOS_PHY_BL_MYNN_tendency_tracer

  !-----------------------------------------------------------------------------
  ! private routines
  !-----------------------------------------------------------------------------

  subroutine get_length( &
       KA, KS, KE_PBL, &
       PT0, q, n2,    &
       SFLX_PT, l_mo, &
       FZ,            &
       l              )
    use scale_const, only: &
       GRAV   => CONST_GRAV, &
       KARMAN => CONST_KARMAN, &
       CP     => CONST_CPdry, &
       EPS    => CONST_EPS
    implicit none
    integer,  intent(in) :: KA, KS, KE_PBL

    real(RP), intent(in) :: PT0(KA)
    real(RP), intent(in) :: q(KA)
    real(RP), intent(in) :: n2(KA)
    real(RP), intent(in) :: SFLX_PT  !> surface temerture flux
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
       int_qz = int_qz + ((FZ(k)+FZ(k-1))*0.5_RP-FZ(KS-1)) * qdz
       int_q  = int_q + qdz
    end do
    ! LT
    lt = min( max(0.23_RP * int_qz / (int_q + EPS), &
                  LT_min), &
                  ATMOS_PHY_BL_MYNN_Lt_MAX )
    rlt = 1.0_RP / lt

    rlm = 1.0_RP / l_mo

    qc = ( GRAV / PT0(KS) * max(SFLX_PT,0.0_RP) * lt )**OneOverThree

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
           +  999.E10_RP * (1.0_RP-sw)

       ! L
       l(k) = 1.0_RP / ( 1.0_RP/ls + rlt + 1.0_RP/lb )
    end do

    return
  end subroutine get_length

  subroutine get_q2_level2( &
       KA, KS, KE_PBL, &
       q2_2,           &
       dudz2, Ri, l    )
    implicit none
    integer,  intent(in)  :: KA, KS, KE_PBL

    real(RP), intent(out) :: q2_2(KA)

    real(RP), intent(in)  :: dudz2(KA)
    real(RP), intent(in)  :: Ri(KA)
    real(RP), intent(in)  :: l(KA)

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

  subroutine get_smsh( &
       KA, KS, KE_PBL, &
       sm, sh,         &
       q, q2_2,        &
       l, n2, dudz2    )
    use scale_const, only: &
         EPS => CONST_EPS
    implicit none
    integer,  intent(in)  :: KA, KS, KE_PBL

    real(RP), intent(out) :: sm(KA)
    real(RP), intent(out) :: sh(KA)

    real(RP), intent(in)  :: q(KA)
    real(RP), intent(in)  :: q2_2(KA)
    real(RP), intent(in)  :: l(KA)
    real(RP), intent(in)  :: n2(KA)
    real(RP), intent(in)  :: dudz2(KA)

    real(RP) :: l2q2 !> L^2/q^2
    real(RP) :: ac   !> \alpha_c
    real(RP) :: ac2  !> \alpha_c^2
    real(RP) :: p1   !> \Phi_1
    real(RP) :: p2   !> \Phi_2
    real(RP) :: p3   !> \Phi_3
    real(RP) :: p4   !> \Phi_4
    real(RP) :: p5   !> \Phi_5
    real(RP) :: rd25 !> 1/D_2.5
    real(RP) :: gh   !> G_H

    integer :: k

    do k = KS, KE_PBL

       ! level 2.5
       ac = min(q(k)/sqrt(q2_2(k)), 1.0_RP)
       ac2 = ac**2
       l2q2 = ( l(k) / max(q(k),EPS) )**2
       gh = - n2(k) * l2q2

       p1 = 1.0_RP - 3.0_RP * ac2 * A2 * B2 * (1.0_RP-C3) * gh
       p2 = 1.0_RP - 9.0_RP * ac2 * A1 * A2 * (1.0_RP-C2) * gh
       p3 = p1 + 9.0_RP * ac2 * A2**2 * (1.0_RP-C2) * (1.0_RP-C5) * gh
       p4 = p1 - 12.0_RP * ac2 * A1 * A2 * (1.0_RP-C2) * gh
       p5 = 6.0_RP * ac2 * A1**2 * dudz2(k) * l2q2

       rd25 = 1.0_RP / max(p2 * p4 + p5 * p3, 1.E-20_RP)
       sm(k) = max( ac * A1 * (p3 - 3.0_RP * C1 * p4) * rd25, 0.0_RP )
       sh(k) = max( ac * A2 * (p2 + 3.0_RP * C1 * p5) * rd25, 0.0_RP )

    end do

    return
  end subroutine get_smsh

  subroutine get_f2h( &
       KA, KS, KE, &
       FZ, &
       f2h )
    integer,  intent(in)  :: KA, KS, KE
    real(RP), intent(in)  :: FZ(0:KA)
    real(RP), intent(out) :: f2h(KA,2)

    real(RP) :: dz1, dz2
    integer :: k

    do k = KS, KE-1
       dz1 = FZ(k+1) - FZ(k  )
       dz2 = FZ(k)   - FZ(k-1)
       f2h(k,1) = dz2 / ( dz1 + dz2 )
       f2h(k,2) = dz1 / ( dz1 + dz2 )
    end do

    return
  end subroutine get_f2h

end module scale_atmos_phy_bl_mynn
