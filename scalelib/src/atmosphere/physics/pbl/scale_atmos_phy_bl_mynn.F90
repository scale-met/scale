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
!! @li Kitamura, 2010:
!!     Modifications to the Mellor-Yamada-Nakanishi-Niino (MYNN) model for the stable stratification case.
!!     J. Meteorol. Soc. Japan, 88, 857-864
!! @li Olson et al., 2019:
!!     A description of the MYNN-EDMF scheme and the coupling to other components in WRF-ARW.
!!     NOAA Technical Memorandum OAR GSD-61
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
     PI     => CONST_PI,    &
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
  public :: ATMOS_PHY_BL_mynn_mkinit
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

  integer,  private, parameter :: I_B71    = 1
  integer,  private, parameter :: I_B91    = 2
  integer,  private, parameter :: I_B91W01 = 3

  real(RP), private, parameter :: OneOverThree = 1.0_RP / 3.0_RP
  real(RP), private, parameter :: LT_min       = 1.E-6_RP
  real(RP), private, parameter :: FLX_LIM_FACT = 0.5_RP

  real(RP), private            :: A1
  real(RP), private            :: A2
  real(RP), private, parameter :: B1 = 24.0_RP
  real(RP), private, parameter :: B2 = 15.0_RP
  real(RP), private            :: C1
  real(RP), private            :: C2 = -1.0_RP ! N01: 0.65,  N04: 0.70,  N09: 0.75,  O19: 0.729
  real(RP), private            :: C3 = -1.0_RP ! N01: 0.294, N04: 0.323, N09: 0.352, O19: 0.34
  real(RP), private, parameter :: C5 = 0.2_RP
  real(RP), private, parameter :: G1 = 0.235_RP
  real(RP), private            :: G2
  real(RP), private            :: F2
  real(RP), private            :: Rf2
  real(RP), private            :: Rfc
  real(RP), private, parameter :: PrN = 0.74_RP
  !$acc declare create(A1, A2, C1, C2, C3, G2, F2, Rf2, Rfc)

  real(RP), private, parameter :: SQRT_2PI  = sqrt( 2.0_RP * PI )
  real(RP), private, parameter :: RSQRT_2PI = 1.0_RP / SQRT_2PI
  real(RP), private, parameter :: RSQRT_2   = 1.0_RP / sqrt( 2.0_RP )

  ! for O2019
  integer, parameter :: max_plume = 10
  integer :: nplume
  real(RP) :: dplume(max_plume)
  real(RP) :: pw(max_plume)
  logical  :: ATMOS_PHY_BL_MYNN_MF
  !$acc declare create(nplume,dplume,pw,ATMOS_PHY_BL_MYNN_MF)

  ! history
  integer,  private :: HIST_Ri
  integer,  private :: HIST_Pr
  integer,  private :: HIST_TKE_pr
  integer,  private :: HIST_TKE_di
  integer,  private :: HIST_dudz2
  integer,  private :: HIST_Lmix
  integer,  private :: HIST_flxU
  integer,  private :: HIST_flxV
  integer,  private :: HIST_flxT
  integer,  private :: HIST_flxQ
  integer,  private :: HIST_flxU2
  integer,  private :: HIST_flxV2
  integer,  private :: HIST_flxT2
  integer,  private :: HIST_flxQ2

  ! namelist
  logical, private  :: ATMOS_PHY_BL_MYNN_K2010      = .false.   !> Kitamura (2010)
  logical, private  :: ATMOS_PHY_BL_MYNN_O2019      = .false.   !> Olson et al. (2019)
  real(RP), private :: ATMOS_PHY_BL_MYNN_PBL_MAX    = 3000.0_RP !> maximum height of the PBL
  real(RP), private :: ATMOS_PHY_BL_MYNN_TKE_MIN    =  1.E-20_RP
  real(RP), private :: ATMOS_PHY_BL_MYNN_N2_MAX     =  1.E1_RP
  real(RP), private :: ATMOS_PHY_BL_MYNN_NU_MAX     =  1.E4_RP
  real(RP), private :: ATMOS_PHY_BL_MYNN_KH_MAX     =  1.E4_RP
  real(RP), private :: ATMOS_PHY_BL_MYNN_Lt_MAX     =  2000.0_RP
  logical,  private :: ATMOS_PHY_BL_MYNN_use_Zi     = .true.
  real(RP), private :: ATMOS_PHY_BL_MYNN_zeta_lim   =  10.0_RP
  real(RP), private :: ATMOS_PHY_BL_MYNN_Sq_fact    = 3.0_RP
  real(RP), private :: ATMOS_PHY_BL_MYNN_cns        = -1.0_RP ! N01: 3.5,   O19: 10.0
  real(RP), private :: ATMOS_PHY_BL_MYNN_alpha2     = -1.0_RP ! N01: 1.0,   O19: 0.3
  real(RP), private :: ATMOS_PHY_BL_MYNN_alpha4     = -1.0_RP ! N01: 100.0, O19: 10.0
  real(RP), private :: ATMOS_PHY_BL_MYNN_DUMP_coef  = 0.5_RP
  logical,  private :: ATMOS_PHY_BL_MYNN_dz_sim     = .true.
  logical,  private :: ATMOS_PHY_BL_MYNN_similarity = .true.
  !$acc declare create(ATMOS_PHY_BL_MYNN_K2010,ATMOS_PHY_BL_MYNN_O2019,ATMOS_PHY_BL_MYNN_PBL_MAX,ATMOS_PHY_BL_MYNN_TKE_MIN,ATMOS_PHY_BL_MYNN_N2_MAX,ATMOS_PHY_BL_MYNN_NU_MAX,ATMOS_PHY_BL_MYNN_KH_MAX,ATMOS_PHY_BL_MYNN_Lt_MAX,ATMOS_PHY_BL_MYNN_use_Zi,ATMOS_PHY_BL_MYNN_zeta_lim,ATMOS_PHY_BL_MYNN_Sq_fact,ATMOS_PHY_BL_MYNN_cns,ATMOS_PHY_BL_MYNN_alpha2,ATMOS_PHY_BL_MYNN_alpha4,ATMOS_PHY_BL_MYNN_DUMP_coef,ATMOS_PHY_BL_MYNN_dz_sim,ATMOS_PHY_BL_MYNN_similarity)

  character(len=H_SHORT), private  :: ATMOS_PHY_BL_MYNN_LEVEL = "2.5" ! "2.5" or "3"

  namelist / PARAM_ATMOS_PHY_BL_MYNN / &
       ATMOS_PHY_BL_MYNN_K2010,    &
       ATMOS_PHY_BL_MYNN_O2019,    &
       ATMOS_PHY_BL_MYNN_PBL_MAX,  &
       ATMOS_PHY_BL_MYNN_N2_MAX,   &
       ATMOS_PHY_BL_MYNN_NU_MAX,   &
       ATMOS_PHY_BL_MYNN_KH_MAX,   &
       ATMOS_PHY_BL_MYNN_Lt_MAX,   &
       ATMOS_PHY_BL_MYNN_use_Zi,   &
       ATMOS_PHY_BL_MYNN_zeta_lim, &
       ATMOS_PHY_BL_MYNN_LEVEL,    &
       ATMOS_PHY_BL_MYNN_Sq_fact,  &
       ATMOS_PHY_BL_MYNN_cns,      &
       ATMOS_PHY_BL_MYNN_alpha2,   &
       ATMOS_PHY_BL_MYNN_alpha4,   &
       ATMOS_PHY_BL_MYNN_DUMP_coef,&
       ATMOS_PHY_BL_MYNN_dz_sim,   &
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
       BULKFLUX_type, &
       dx, &
       TKE_MIN, PBL_MAX )
    use scale_prc, only: &
       PRC_abort
    use scale_file_history, only: &
       FILE_HISTORY_reg
    implicit none
    character(len=*), intent(in) :: BULKFLUX_type

    real(RP), intent(in), optional :: dx
    real(RP), intent(in), optional :: TKE_MIN
    real(RP), intent(in), optional :: PBL_MAX

    integer :: ierr
    integer :: n
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

    ATMOS_PHY_BL_MYNN_MF = .false.
    if ( ATMOS_PHY_BL_MYNN_O2019 ) then
       if ( ATMOS_PHY_BL_MYNN_cns < 0.0_RP ) ATMOS_PHY_BL_MYNN_cns = 3.5_RP
       if ( ATMOS_PHY_BL_MYNN_alpha2 < 0.0_RP ) ATMOS_PHY_BL_MYNN_alpha2 = 0.3_RP
       if ( ATMOS_PHY_BL_MYNN_alpha4 < 0.0_RP ) ATMOS_PHY_BL_MYNN_alpha4 = 10.0_RP
       if ( C2 < 0.0_RP ) C2 = 0.729_RP
       if ( C3 < 0.0_RP ) C3 = 0.34_RP
       ATMOS_PHY_BL_MYNN_K2010 = .true.
       ATMOS_PHY_BL_MYNN_use_Zi = .true.
       if ( ATMOS_PHY_BL_MYNN_LEVEL .ne. "3" ) then
          ATMOS_PHY_BL_MYNN_MF = .true.
       end if
       if ( .not. present(dx) ) then
          LOG_ERROR("ATMOS_PHY_BL_MYNN_setup",*) 'dx must be set with O2019'
          call PRC_abort
       end if
       if ( ATMOS_PHY_BL_MYNN_MF ) then
          nplume = max_plume
          do n = 1, max_plume
             dplume(n) = 100.0_RP * n
             pw(n) = 0.1_RP + ( 0.5_RP - 0.1_RP ) * ( n - 1 ) / ( max_plume - 1 ) ! not described in O2019
             if ( dplume(n) > dx ) then
                nplume = n - 1
                exit
             end if
          end do
       end if

    else
       if ( ATMOS_PHY_BL_MYNN_cns < 0.0_RP ) ATMOS_PHY_BL_MYNN_cns = 2.7_RP
       if ( ATMOS_PHY_BL_MYNN_alpha2 < 0.0_RP ) ATMOS_PHY_BL_MYNN_alpha2 = 1.0_RP
       if ( ATMOS_PHY_BL_MYNN_alpha4 < 0.0_RP ) ATMOS_PHY_BL_MYNN_alpha4 = 100.0_RP
       if ( C2 < 0.0_RP ) C2 = 0.75_RP
       if ( C3 < 0.0_RP ) C3 = 0.352_RP
    end if

    !$acc update device(nplume,dplume,pw,ATMOS_PHY_BL_MYNN_MF)


    A1  = B1 * (1.0_RP - 3.0_RP * G1) / 6.0_RP
    A2  = 1.0_RP / (3.0_RP * G1 * B1**(1.0_RP/3.0_RP) * PrN )
    C1  = G1 - 1.0_RP / ( 3.0_RP * A1 * B1**(1.0_RP/3.0_RP) )
    G2  = ( 2.0_RP * A1 * (3.0_RP - 2.0_RP * C2) + B2 * (1.0_RP - C3) ) / B1
    F2  = B1 * (G1 + G2) - 3.0_RP * A1 * (1.0_RP - C2)

    Rf2 = B1 * G1 / F2
    Rfc = G1 / (G1 + G2)

    !$acc update device(A1, A2, C1, C2, C3, G2, F2, Rf2, Rfc)

    select case ( BULKFLUX_type )
    case ( "B91", "B91W01" )
       ! do nothing
    case default
       ATMOS_PHY_BL_MYNN_similarity = .false.
    end select

    !$acc update device(ATMOS_PHY_BL_MYNN_K2010,ATMOS_PHY_BL_MYNN_O2019,ATMOS_PHY_BL_MYNN_PBL_MAX,ATMOS_PHY_BL_MYNN_TKE_MIN,ATMOS_PHY_BL_MYNN_N2_MAX,ATMOS_PHY_BL_MYNN_NU_MAX,ATMOS_PHY_BL_MYNN_KH_MAX,ATMOS_PHY_BL_MYNN_Lt_MAX,ATMOS_PHY_BL_MYNN_zeta_lim,ATMOS_PHY_BL_MYNN_Sq_fact,ATMOS_PHY_BL_MYNN_cns,ATMOS_PHY_BL_MYNN_alpha2,ATMOS_PHY_BL_MYNN_alpha4,ATMOS_PHY_BL_MYNN_DUMP_coef,ATMOS_PHY_BL_MYNN_dz_sim,ATMOS_PHY_BL_MYNN_similarity)


    ! history
    call FILE_HISTORY_reg('Ri_MYNN',       'Richardson number', '1',     HIST_Ri,     fill_halo=.true. )
    call FILE_HISTORY_reg('Pr_MYNN',       'Prandtl number',    '1',     HIST_Pr,     fill_halo=.true., dim_type="ZHXY" )
    call FILE_HISTORY_reg('TKE_prod_MYNN', 'TKE production',    'm2/s3', HIST_TKE_pr, fill_halo=.true.)
    call FILE_HISTORY_reg('TKE_diss_MYNN', 'TKE dissipation',   'm2/s3', HIST_TKE_di, fill_halo=.true.)
    call FILE_HISTORY_reg('dUdZ2_MYNN',    'dudz2',             'm2/s2', HIST_dudz2,  fill_halo=.true.)
    call FILE_HISTORY_reg('L_mix_MYNN',    'minxing length',    'm',     HIST_Lmix,   fill_halo=.true.)

    call FILE_HISTORY_reg('ZFLX_RHOU_BL', 'Z FLUX of RHOU (MYNN)',  'kg/m/s2',   HIST_flxU, fill_halo=.true., dim_type="ZHXY" )
    call FILE_HISTORY_reg('ZFLX_RHOV_BL', 'Z FLUX of RHOV (MYNN)',  'kg/m/s2',   HIST_flxV, fill_halo=.true., dim_type="ZHXY" )
    call FILE_HISTORY_reg('ZFLX_RHOT_BL', 'Z FLUX of RHOT (MYNN)',  'K kg/m2/s', HIST_flxT, fill_halo=.true., dim_type="ZHXY" )
    call FILE_HISTORY_reg('ZFLX_QV_BL',   'Z FLUX of RHOQV (MYNN)', 'kg/m2/s',   HIST_flxQ, fill_halo=.true., dim_type="ZHXY" )

    call FILE_HISTORY_reg('ZFLX_RHOU2_BL', 'Z FLUX of RHOU (MYNN)',  'kg/m/s2',   HIST_flxU2, fill_halo=.true., dim_type="ZHXY" )
    call FILE_HISTORY_reg('ZFLX_RHOV2_BL', 'Z FLUX of RHOV (MYNN)',  'kg/m/s2',   HIST_flxV2, fill_halo=.true., dim_type="ZHXY" )
    call FILE_HISTORY_reg('ZFLX_RHOT2_BL', 'Z FLUX of RHOT (MYNN)',  'K kg/m2/s', HIST_flxT2, fill_halo=.true., dim_type="ZHXY" )
    call FILE_HISTORY_reg('ZFLX_QV2_BL',   'Z FLUX of RHOQV (MYNN)', 'kg/m2/s',   HIST_flxQ2, fill_halo=.true., dim_type="ZHXY" )

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
  !> ATMOS_PHY_BL_MYNN_mkinit
  !! initialize TKE
  !<
  subroutine ATMOS_PHY_BL_MYNN_mkinit( &
       KA, KS, KE, IA, IS, IE, JA, JS, JE, &
       PROG,                               &
       DENS, U, V, W, POTT,                &
       PRES, EXNER, N2,                    &
       QDRY, QV, Qw, POTL, POTV, SFC_DENS, &
       SFLX_MU, SFLX_MV, SFLX_SH, SFLX_QV, &
       us, ts, qs, RLmo,                   &
       frac_land,                          &
       CZ, FZ, F2H,                        &
       BULKFLUX_type                       )
    use scale_const, only: &
       CPdry   => CONST_CPdry
    use scale_prc, only: &
       PRC_abort
    use scale_atmos_hydrometeor, only: &
       HYDROMETEOR_LHV => ATMOS_HYDROMETEOR_LHV, &
       CP_VAPOR
    use scale_matrix, only: &
       MATRIX_SOLVER_tridiagonal
    use scale_file_history, only: &
       FILE_HISTORY_query, &
       FILE_HISTORY_put
    implicit none

    integer, intent(in) :: KA, KS, KE
    integer, intent(in) :: IA, IS, IE
    integer, intent(in) :: JA, JS, JE

    real(RP), intent(out) :: PROG(KA,IA,JA,ATMOS_PHY_BL_MYNN_ntracer) !> prognostic variables (TKE, TSQ, QSQ, COV)

    real(RP), intent(in) :: DENS    (KA,IA,JA) !> density
    real(RP), intent(in) :: U       (KA,IA,JA) !> zonal wind
    real(RP), intent(in) :: V       (KA,IA,JA) !> meridional wind
    real(RP), intent(in) :: W       (KA,IA,JA) !> vertical wind
    real(RP), intent(in) :: POTT    (KA,IA,JA) !> potential temperature
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

    real(RP), intent(in) :: frac_land(IA,JA)

    real(RP), intent(in)  :: CZ(  KA,IA,JA)
    real(RP), intent(in)  :: FZ(0:KA,IA,JA)
    real(RP), intent(in)  :: F2H(KA,2,IA,JA)

    character(len=*), intent(in) :: BULKFLUX_type

    ! bulkflux
    integer :: I_B_TYPE

    real(RP) :: Nu     (KA) !> eddy viscosity coefficient @ half-level
    real(RP) :: Kh     (KA) !> eddy diffusion coefficient @ half-level
    real(RP) :: Qlp    (KA) !> liquid-water content in partial condensation
    real(RP) :: cldfrac(KA) !> cloud fraction in partial condensation
    real(RP) :: Zi          !> depth of the boundary layer (not exact)
    real(RP) :: SFLX_BUOY   !> surface flux of buoyancy: g / \Theta_0 <w \theta_v> @ surface

    real(RP) :: Ri   (KA) !> Richardson number
    real(RP) :: Pr   (KA) !> Plandtle number
    real(RP) :: prod (KA) !> TKE production term
    real(RP) :: diss (KA) !> TKE dissipation term
    real(RP) :: dudz2(KA) !> (du/dz)^2 + (dv/dz)^2
    real(RP) :: l    (KA) !> length scale L
    real(RP) :: rho_h(KA)       !> dens at the half level

    real(RP) :: RHO   (KA) !> dens after updated
    real(RP) :: RHONu (KA) !> dens * Nu at the half level for level 2.5
    real(RP) :: RHOKh (KA) !> dens * Kh at the half level for level 2.5
    real(RP) :: N2_new(KA) !> squared Brunt-Baisala frequency
    real(RP) :: SFLX_PT    !> surface potential temperature flux * density
    real(RP) :: sm25  (KA) !> stability function for velocity for level 2.5
    real(RP) :: sh25  (KA) !> stability function for scalars for level 2.5
    real(RP) :: q     (KA) !> q
    real(RP) :: lq    (KA) !> L * q

    real(RP) :: tke(KA)

    ! for level 3
    real(RP) :: tsq   (KA)
    real(RP) :: qsq   (KA)
    real(RP) :: cov   (KA)
    real(RP) :: prod_t(KA)
    real(RP) :: prod_q(KA)
    real(RP) :: prod_c(KA)
    real(RP) :: diss_p(KA)
    real(RP) :: dtldz(KA)
    real(RP) :: dqwdz(KA)
    real(RP) :: smp    (KA) !> stability function for velocity for the countergradient
    real(RP) :: shpgh(KA)   !> stability function for scalars for the countergradient
    real(RP) :: gammat (KA)
    real(RP) :: gammaq (KA)

    real(RP) :: CPtot

    real(RP) :: FDZ(KA)
    real(RP) :: CDZ(KA)
    real(RP) :: z  (KA)

    logical :: mynn_level3

    real(RP), parameter :: dt = 1.0_RP

    ! for O2019
    real(RP) :: tflux(0:KA)
    real(RP) :: qflux(0:KA)
    real(RP) :: uflux(0:KA)
    real(RP) :: vflux(0:KA)
    real(RP) :: eflux(0:KA)

    ! work
    real(RP) :: PTLV (KA)
    real(RP) :: Nu_f (KA) !> Nu at the full level
    real(RP) :: Kh_f (KA) !> Kh at the full level
    real(RP) :: q2_2 (KA) !> q^2 for level 2
    real(RP) :: ac   (KA) !> \alpha_c
    real(RP) :: betat(KA) !> \beta_t
    real(RP) :: betaq(KA) !> \beta_q

    real(RP) :: flx(0:KA)
    real(RP) :: a(KA)
    real(RP) :: b(KA)
    real(RP) :: c(KA)
    real(RP) :: d(KA)

    real(RP) :: f_smp  (KA) !> stability function for velocity for the countergradient (factor)
    real(RP) :: f_shpgh(KA) !> stability function for scalars for the countergradient (factor)
    real(RP) :: f_gamma(KA) !> - E_H / q^2 * GRAV / Theta_0
    real(RP) :: tltv25 (KA)
    real(RP) :: qwtv25 (KA)
    real(RP) :: tvsq25 (KA)
    real(RP) :: tvsq_up(KA) !> upper limit of <\theta_v^2> - <\theta_v^2>2.5
    real(RP) :: tvsq_lo(KA) !> lower limit
    real(RP) :: dtsq   (KA)
    real(RP) :: dqsq   (KA)
    real(RP) :: dcov   (KA)

    real(RP) :: mflux(0:KA)

    real(RP) :: Uh(KA), Vh(KA), Wh(KA)
    real(RP) :: qw2(KA), qh(KA)
    real(RP) :: tlh(KA), tvh(KA)
    real(RP) :: eh(KA), dh(KA)
    real(RP) :: TEML(KA), LHVL(KA), psat(KA)

#ifdef _OPENACC
    real(RP) :: work(KA,4)
#endif

    integer :: KE_PBL
    integer :: k, i, j

    !---------------------------------------------------------------------------

    mynn_level3 = ( ATMOS_PHY_BL_MYNN_LEVEL == "3" )

    select case ( BULKFLUX_type )
    case ( 'B71' )
       I_B_TYPE = I_B71
    case ( 'B91' )
       I_B_TYPE = I_B91
    case ( 'B91W01' )
       I_B_TYPE = I_B91W01
    case default
       LOG_ERROR("ATMOS_PHY_BL_MYNN_init_TKE",*) "BULKFLUX_type is invalid: ", trim(BULKFLUX_type)
       call PRC_abort
    end select

    !$acc data copyout(PROG) copyin(DENS,U,V,POTT,PRES,EXNER,N2,QDRY,QV,Qw,POTL,POTV,SFC_DENS,SFLX_MU,SFLX_MV,SFLX_SH,SFLX_QV,us,ts,qs,RLmo,CZ,FZ,F2H)


!OCL INDEPENDENT
    !$omp parallel do default(none) &
    !$omp OMP_SCHEDULE_ collapse(2) &
    !$omp shared(KA,KS,KE,IS,IE,JS,JE, &
    !$omp        CPdry,CP_VAPOR,UNDEF, &
    !$omp        ATMOS_PHY_BL_MYNN_N2_MAX,ATMOS_PHY_BL_MYNN_TKE_MIN, &
    !$omp        ATMOS_PHY_BL_MYNN_NU_MAX,ATMOS_PHY_BL_MYNN_KH_MAX, &
    !$omp        ATMOS_PHY_BL_MYNN_Sq_fact, ATMOS_PHY_BL_MYNN_zeta_lim, &
    !$omp        ATMOS_PHY_BL_MYNN_similarity,ATMOS_PHY_BL_MYNN_dz_sim, &
    !$omp        ATMOS_PHY_BL_MYNN_DUMP_coef, &
    !$omp        ATMOS_PHY_BL_MYNN_PBL_MAX, &
    !$omp        ATMOS_PHY_BL_MYNN_K2010,ATMOS_PHY_BL_MYNN_O2019,ATMOS_PHY_BL_MYNN_MF, &
    !$omp        DENS,PROG,U,V,W,POTT,PRES,QDRY,QV,Qw,POTV,POTL,EXNER,N2, &
    !$omp        SFC_DENS,SFLX_MU,SFLX_MV,SFLX_SH,SFLX_QV,us,ts,qs,RLmo, &
    !$omp        mynn_level3, &
    !$omp        frac_land, &
    !$omp        CZ,FZ,F2H, &
    !$omp        I_B_TYPE) &
    !$omp private(N2_new,lq,sm25,sh25,q, &
    !$omp         SFLX_PT,CPtot,RHO,RHONu,RHOKh, &
    !$omp         Nu,Kh,Qlp,cldfrac,Zi,SFLX_BUOY, &
    !$omp         Ri,Pr,prod,diss,dudz2,l, &
    !$omp         smp,shpgh, &
    !$omp         dtldz,dqwdz,gammat,gammaq, &
    !$omp         rho_h,tke,CDZ,FDZ,z, &
    !$omp         tsq,qsq,cov, &
    !$omp         prod_t,prod_q,prod_c,diss_p, &
    !$omp         tflux,qflux,uflux,vflux,eflux, &
    !$omp         PTLV, Nu_f, Kh_f, q2_2, ac, betat, betaq, &
    !$omp         flx, a, b, c, d, &
    !$omp         f_smp, f_shpgh, f_gamma, tltv25, qwtv25, tvsq25, &
    !$omp         tvsq_up, tvsq_lo, dtsq, dqsq, dcov, &
    !$omp         mflux, &
    !$omp         Uh, Vh, Wh, qw2, qh, tlh, tvh, eh, dh, &
    !$omp         TEML, LHVL, psat, &
    !$omp         KE_PBL,k,i,j)
    !$acc kernels
    !$acc loop independent collapse(2) &
    !$acc private(N2_new,lq,sm25,sh25,q, &
    !$acc         SFLX_PT,CPtot,RHO,RHONu,RHOKh, &
    !$acc         Nu,Kh,Qlp,cldfrac,Zi,SFLX_BUOY, &
    !$acc         Ri,Pr,prod,diss,dudz2,l, &
    !$acc         smp,shpgh, &
    !$acc         dtldz,dqwdz,gammat,gammaq, &
    !$acc         rho_h,tke,CDZ,FDZ,z, &
    !$acc         tsq,qsq,cov, &
    !$acc         prod_t,prod_q,prod_c,diss_p, &
    !$acc         tflux,qflux,uflux,vflux,eflux, &
    !$acc         work, &
    !$acc         PTLV, Nu_f, Kh_f, q2_2, ac, betat, betaq, &
    !$acc         flx, a, b, c, d, &
    !$acc         f_smp, f_shpgh, f_gamma, tltv25, qwtv25, tvsq25, &
    !$acc         tvsq_up, tvsq_lo, dtsq, dqsq, dcov, &
    !$acc         mflux, &
    !$acc         Uh, Vh, Wh, qw2, qh, tlh, tvh, eh, dh, &
    !$acc         TEML, LHVL, psat, &
    !$acc         KE_PBL,k,i,j)
    do j = JS, JE
    do i = IS, IE

       do k = KS, KE-1
          z(k) = CZ(k,i,j) - FZ(KS-1,i,j)
       end do

       KE_PBL = KS+1
       !$acc loop seq
       do k = KS+2, KE-1
          if ( ATMOS_PHY_BL_MYNN_PBL_MAX >= z(k) ) then
             KE_PBL = k
          else
             exit
          end if
       end do

       do k = KS, KE_PBL+1
          CDZ(k) = FZ(k,i,j) - FZ(k-1,i,j)
       end do
       do k = KS, KE_PBL
          FDZ(k) = CZ(k+1,i,j) - CZ(k,i,j)
       end do

       do k = KS, KE_PBL
          tke(k) = 0.01_RP
       end do

       CPtot = CPdry + SFLX_QV(i,j) * ( CP_VAPOR - CPdry )
       SFLX_PT = SFLX_SH(i,j) / ( CPtot * EXNER(KS,i,j) )

       call MYNN_main( KA, KS, KE_PBL, &
                       i, j, &
                       tke(:), tsq(:), qsq(:), cov(:),                     & ! (inout)
                       q(:), l(:), lq(:),                                  & ! (out)
                       Nu(:), RHONu(:), Kh(:), RHOKh(:), Pr(:),            & ! (out)
                       prod(:), prod_t(:), prod_q(:), prod_c(:),           & ! (out)
                       diss(:), diss_p(:),                                 & ! (out)
                       sm25(:), smp(:), sh25(:), shpgh(:),                 & ! (out)
                       gammat(:), gammaq(:),                               & ! (out)
                       dudz2(:), n2_new(:), Ri(:),                         & ! (out)
                       dtldz(:), dqwdz(:),                                 & ! (out)
                       RHO(:), rho_h(:),                                   & ! (out)
                       uflux(:), vflux(:), tflux(:), qflux(:), eflux(:),   & ! (out)
                       Qlp(:), cldfrac(:), Zi, SFLX_BUOY,                  & ! (out)
                       U(:,i,j), V(:,i,j), W(:,i,j),                       & ! (in)
                       DENS(:,i,j), PRES(:,i,j),                           & ! (in)
                       POTT(:,i,j), POTL(:,i,j), POTV(:,i,j),              & ! (in)
                       Qw(:,i,j), N2(:,i,j),                               & ! (in)
                       EXNER(:,i,j), QDRY(:,i,j),                          & ! (in)
                       SFLX_PT, SFLX_SH(i,j), SFLX_QV(i,j), SFC_DENS(i,j), & ! (in)
                       RLmo(i,j), us(i,j), ts(i,j), qs(i,j),               & ! (in)
                       z(:), CDZ(:), FDZ(:), F2H(:,:,i,j),                 & ! (in)
                       frac_land(i,j),                                     & ! (in)
                       dt,                                                 & ! (in)
                       I_B_TYPE,                                           & ! (in)
                       mynn_level3, .true.,                                & ! (in)
#ifdef _OPENACC
                       work(:,:),                                          & ! (work)
#endif
                       PTLV(:), Nu_f(:), Kh_f(:), q2_2(:), ac(:),          & ! (work)
                       betat(:), betaq(:),                                 & ! (work)
                       flx(:), a(:), b(:), c(:), d(:),                     & ! (work)
                       f_smp(:), f_shpgh(:), f_gamma(:),                   & ! (work)
                       tltv25(:), qwtv25(:), tvsq25(:),                    & ! (work)
                       tvsq_up(:), tvsq_lo(:), dtsq(:), dqsq(:), dcov(:),  & ! (work)
                       mflux(:),                                           & ! (work)
                       Uh(:), Vh(:), Wh(:), qw2(:), qh(:),                 & ! (work)
                       tlh(:), tvh(:), eh(:), dh(:),                       & ! (work)
                       TEML(:), LHVL(:), psat(:)                           ) ! (work)

       do k = KS, KE_PBL
          PROG(k,i,j,I_TKE) = tke(k)
       end do
       do k = KE_PBL+1, KE
          PROG(k,i,j,I_TKE) = 0.0_RP
       end do

       if ( mynn_level3 ) then
          do k = KS, KE_PBL
             PROG(k,i,j,I_TSQ) = tsq(k)
             PROG(k,i,j,I_QSQ) = qsq(k)
             PROG(k,i,j,I_COV) = cov(k)
          end do
          do k = KE_PBL+1, KE
             PROG(k,i,j,I_TSQ) = 0.0_RP
             PROG(k,i,j,I_QSQ) = 0.0_RP
             PROG(k,i,j,I_COV) = 0.0_RP
          end do
       end if

    end do
    end do
    !$acc end kernels

    !$acc end data

    return
  end subroutine ATMOS_PHY_BL_MYNN_mkinit

  !-----------------------------------------------------------------------------
  !> ATMOS_PHY_BL_MYNN_tendency
  !! calculate tendency by the vertical eddy viscosity
  !<
  subroutine ATMOS_PHY_BL_MYNN_tendency( &
       KA, KS, KE, IA, IS, IE, JA, JS, JE, &
       DENS, U, V, W, POTT, PROG,          &
       PRES, EXNER, N2,                    &
       QDRY, QV, Qw, POTL, POTV, SFC_DENS, &
       SFLX_MU, SFLX_MV, SFLX_SH, SFLX_QV, &
       us, ts, qs, RLmo,                   &
       frac_land,                          &
       CZ, FZ, F2H, dt_DP,                 &
       BULKFLUX_type,                      &
       RHOU_t, RHOV_t, RHOT_t, RHOQV_t,    &
       RPROG_t,                            &
       Nu, Kh, Qlp, cldfrac,               &
       Zi, SFLX_BUOY                       )
    use scale_const, only: &
       CPdry   => CONST_CPdry
    use scale_prc, only: &
       PRC_abort
    use scale_atmos_hydrometeor, only: &
       HYDROMETEOR_LHV => ATMOS_HYDROMETEOR_LHV, &
       CP_VAPOR
    use scale_matrix, only: &
       MATRIX_SOLVER_tridiagonal
    use scale_file_history, only: &
       FILE_HISTORY_query, &
       FILE_HISTORY_put
    implicit none

    integer, intent(in) :: KA, KS, KE
    integer, intent(in) :: IA, IS, IE
    integer, intent(in) :: JA, JS, JE

    real(RP), intent(in) :: DENS    (KA,IA,JA) !> density
    real(RP), intent(in) :: U       (KA,IA,JA) !> zonal wind
    real(RP), intent(in) :: V       (KA,IA,JA) !> meridional wind
    real(RP), intent(in) :: W       (KA,IA,JA) !> vertical wind
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

    real(RP), intent(in) :: frac_land(IA,JA)

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

    ! bulkflux
    integer :: I_B_TYPE

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

    real(RP) :: flxU2(0:KA,IA,JA) !> dens * w * u (counter gradient or mass flux)
    real(RP) :: flxV2(0:KA,IA,JA) !> dens * w * v (counter gradient or mass flux)
    real(RP) :: flxT2(0:KA,IA,JA) !> dens * w * pt (counter gradient or mass flux)
    real(RP) :: flxQ2(0:KA,IA,JA) !> dens * w * qv (counter gradient or mass flux)

    real(RP) :: RHO   (KA) !> dens after updated
    real(RP) :: RHONu (KA) !> dens * Nu at the half level for level 2.5
    real(RP) :: RHOKh (KA) !> dens * Kh at the half level for level 2.5
    real(RP) :: N2_new(KA) !> squared Brunt-Baisala frequency
    real(RP) :: SFLX_PT    !> surface potential temperature flux * density
    real(RP) :: sm25  (KA) !> stability function for velocity for level 2.5
    real(RP) :: sh25  (KA) !> stability function for scalars for level 2.5
    real(RP) :: q     (KA) !> q
    real(RP) :: lq    (KA) !> L * q

    real(RP) :: tke(KA)

    ! for level 3
    real(RP) :: tsq   (KA)
    real(RP) :: qsq   (KA)
    real(RP) :: cov   (KA)
    real(RP) :: prod_t(KA)
    real(RP) :: prod_q(KA)
    real(RP) :: prod_c(KA)
    real(RP) :: diss_p(KA)
    real(RP) :: dtldz(KA)
    real(RP) :: dqwdz(KA)
    real(RP) :: smp    (KA) !> stability function for velocity for the countergradient
    real(RP) :: shpgh(KA)   !> stability function for scalars for the countergradient
    real(RP) :: gammat (KA)
    real(RP) :: gammaq (KA)
    real(RP) :: rlqsm_h(KA) !> DENS * L * q * SM' @ half level

    real(RP) :: flx(0:KA)
    real(RP) :: a(KA)
    real(RP) :: b(KA)
    real(RP) :: c(KA)
    real(RP) :: d(KA)
    real(RP) :: ap
    real(RP) :: phi_N(KA)
    real(RP) :: dummy(KA)

    real(RP) :: CPtot

    real(RP) :: sf_t

    real(RP) :: FDZ(KA)
    real(RP) :: CDZ(KA)
    real(RP) :: z  (KA)

    logical :: mynn_level3

    real(RP) :: dt

    real(RP) :: fmin
    integer  :: kmin

    ! for O2019
    real(RP) :: tflux(0:KA)
    real(RP) :: qflux(0:KA)
    real(RP) :: uflux(0:KA)
    real(RP) :: vflux(0:KA)
    real(RP) :: eflux(0:KA)

    real(RP) :: tmp

    logical :: do_put

    integer :: KE_PBL
    integer :: k, i, j


    ! work
    real(RP) :: PTLV (KA)
    real(RP) :: Nu_f (KA) !> Nu at the full level
    real(RP) :: Kh_f (KA) !> Kh at the full level
    real(RP) :: q2_2 (KA) !> q^2 for level 2
    real(RP) :: ac   (KA) !> \alpha_c
    real(RP) :: betat(KA) !> \beta_t
    real(RP) :: betaq(KA) !> \beta_q

    real(RP) :: f_smp  (KA) !> stability function for velocity for the countergradient (factor)
    real(RP) :: f_shpgh(KA) !> stability function for scalars for the countergradient (factor)
    real(RP) :: f_gamma(KA) !> - E_H / q^2 * GRAV / Theta_0
    real(RP) :: tltv25 (KA)
    real(RP) :: qwtv25 (KA)
    real(RP) :: tvsq25 (KA)
    real(RP) :: tvsq_up(KA) !> upper limit of <\theta_v^2> - <\theta_v^2>2.5
    real(RP) :: tvsq_lo(KA) !> lower limit
    real(RP) :: dtsq   (KA)
    real(RP) :: dqsq   (KA)
    real(RP) :: dcov   (KA)

    real(RP) :: mflux(0:KA)

    real(RP) :: Uh(KA), Vh(KA), Wh(KA)
    real(RP) :: qw2(KA), qh(KA)
    real(RP) :: tlh(KA), tvh(KA)
    real(RP) :: eh(KA), dh(KA)
    real(RP) :: TEML(KA), LHVL(KA), psat(KA)

#ifdef _OPENACC
    real(RP) :: work(KA,4)
#endif

    !---------------------------------------------------------------------------

    dt = real(dt_DP, RP)

    LOG_PROGRESS(*) "atmosphere / physics / pbl / MYNN"

    mynn_level3 = ( ATMOS_PHY_BL_MYNN_LEVEL == "3" )

    select case ( BULKFLUX_type )
    case ( 'B71' )
       I_B_TYPE = I_B71
    case ( 'B91' )
       I_B_TYPE = I_B91
    case ( 'B91W01' )
       I_B_TYPE = I_B91W01
    case default
       LOG_ERROR("ATMOS_PHY_BL_MYNN_tendency",*) "BULKFLUX_type is invalid: ", trim(BULKFLUX_type)
       call PRC_abort
    end select

    !$acc data copyin(DENS,U,V,POTT,PROG,PRES,EXNER,N2,QDRY,QV,Qw,POTL,POTV,SFC_DENS,SFLX_MU,SFLX_MV,SFLX_SH,SFLX_QV,us,ts,qs,RLmo,CZ,FZ,F2H) &
    !$acc      copyout(RHOU_t,RHOV_t,RHOT_t,RHOQV_t,RPROG_t,Nu,Kh,Qlp,cldfrac,Zi,SFLX_BUOY), &
    !$acc      create(Ri,Pr,prod,diss,dudz2,l,flxU,flxV,flxT,flxQ)



!OCL INDEPENDENT
    !$omp parallel do default(none) &
    !$omp OMP_SCHEDULE_ collapse(2) &
    !$omp shared(KA,KS,KE,IS,IE,JS,JE, &
    !$omp        CPdry,CP_VAPOR,UNDEF, &
    !$omp        ATMOS_PHY_BL_MYNN_N2_MAX,ATMOS_PHY_BL_MYNN_TKE_MIN, &
    !$omp        ATMOS_PHY_BL_MYNN_NU_MAX,ATMOS_PHY_BL_MYNN_KH_MAX, &
    !$omp        ATMOS_PHY_BL_MYNN_Sq_fact, ATMOS_PHY_BL_MYNN_zeta_lim, &
    !$omp        ATMOS_PHY_BL_MYNN_similarity,ATMOS_PHY_BL_MYNN_dz_sim, &
    !$omp        ATMOS_PHY_BL_MYNN_DUMP_coef, &
    !$omp        ATMOS_PHY_BL_MYNN_PBL_MAX, &
    !$omp        ATMOS_PHY_BL_MYNN_K2010,ATMOS_PHY_BL_MYNN_O2019,ATMOS_PHY_BL_MYNN_MF, &
    !$omp        RHOU_t,RHOV_t,RHOT_t,RHOQV_t,RPROG_t,Nu,Kh,Qlp,cldfrac,Zi,SFLX_BUOY, &
    !$omp        DENS,PROG,U,V,W,POTT,PRES,QDRY,QV,Qw,POTV,POTL,EXNER,N2, &
    !$omp        SFC_DENS,SFLX_MU,SFLX_MV,SFLX_SH,SFLX_QV,us,ts,qs,RLmo, &
    !$omp        mynn_level3, &
    !$omp        frac_land, &
    !$omp        CZ,FZ,F2H,dt, &
    !$omp        I_B_TYPE, &
    !$omp        Ri,Pr,prod,diss,dudz2,l,flxU,flxV,flxT,flxQ,flxU2,flxV2,flxT2,flxQ2) &
    !$omp private(N2_new,lq,sm25,sh25,rlqsm_h,q, &
    !$omp         SFLX_PT,CPtot,RHO,RHONu,RHOKh, &
    !$omp         smp,shpgh, &
    !$omp         dtldz,dqwdz,gammat,gammaq, &
    !$omp         flx,a,b,c,d,ap,rho_h,phi_n,tke,sf_t,CDZ,FDZ,z, &
    !$omp         dummy, &
    !$omp         tsq,qsq,cov, &
    !$omp         prod_t,prod_q,prod_c,diss_p, &
    !$omp         fmin,kmin, &
    !$omp         tflux,qflux,uflux,vflux,eflux, &
    !$omp         tmp, &
    !$omp         PTLV, Nu_f, Kh_f, q2_2, ac, betat, betaq, &
    !$omp         f_smp, f_shpgh, f_gamma, tltv25, qwtv25, tvsq25, &
    !$omp         tvsq_up, tvsq_lo, dtsq, dqsq, dcov, &
    !$omp         mflux, &
    !$omp         Uh, Vh, Wh, qw2, qh, tlh, tvh, eh, dh, &
    !$omp         TEML, LHVL, psat, &
    !$omp         KE_PBL,k,i,j)
    !$acc kernels
    !$acc loop independent collapse(2) &
    !$acc private(N2_new,lq,sm25,sh25,rlqsm_h,q, &
    !$acc         SFLX_PT,CPtot,RHO,RHONu,RHOKh, &
    !$acc         smp,shpgh, &
    !$acc         dtldz,dqwdz,gammat,gammaq, &
    !$acc         flx,a,b,c,d,ap,rho_h,phi_n,tke,sf_t,CDZ,FDZ,z, &
    !$acc         dummy, &
    !$acc         tsq,qsq,cov, &
    !$acc         prod_t,prod_q,prod_c,diss_p, &
    !$acc         fmin,kmin, &
    !$acc         tflux,qflux,uflux,vflux,eflux, &
    !$acc         tmp, &
    !$acc         work, &
    !$acc         PTLV, Nu_f, Kh_f, q2_2, ac, betat, betaq, &
    !$acc         f_smp, f_shpgh, f_gamma, tltv25, qwtv25, tvsq25, &
    !$acc         tvsq_up, tvsq_lo, dtsq, dqsq, dcov, &
    !$acc         mflux, &
    !$acc         Uh, Vh, Wh, qw2, qh, tlh, tvh, eh, dh, &
    !$acc         TEML, LHVL, psat, &
    !$acc         KE_PBL,k,i,j)

    do j = JS, JE
    do i = IS, IE

       do k = KS, KE-1
          z(k) = CZ(k,i,j) - FZ(KS-1,i,j)
       end do

       KE_PBL = KS+1
       !$acc loop seq
       do k = KS+2, KE-1
          if ( ATMOS_PHY_BL_MYNN_PBL_MAX >= z(k) ) then
             KE_PBL = k
          else
             exit
          end if
       end do

       do k = KS, KE_PBL+1
          CDZ(k) = FZ(k,i,j) - FZ(k-1,i,j)
       end do
       do k = KS, KE_PBL
          FDZ(k) = CZ(k+1,i,j) - CZ(k,i,j)
       end do

       do k = KS, KE_PBL
          tke(k) = max( PROG(k,i,j,I_TKE), ATMOS_PHY_BL_MYNN_TKE_MIN )
       end do

       if ( mynn_level3 ) then
          do k = KS, KE_PBL
             tsq(k) = max( PROG(k,i,j,I_TSQ), 0.0_RP )
             qsq(k) = max( PROG(k,i,j,I_QSQ), 0.0_RP )
             cov(k) = PROG(k,i,j,I_COV)
             cov(k) = sign( min( abs(cov(k)), sqrt(tsq(k) * qsq(k))), cov(k) )
          end do
       end if

       CPtot = CPdry + SFLX_QV(i,j) * ( CP_VAPOR - CPdry )
       SFLX_PT = SFLX_SH(i,j) / ( CPtot * EXNER(KS,i,j) )

       call MYNN_main( KA, KS, KE_PBL, &
                       i, j, &
                       tke(:), tsq(:), qsq(:), cov(:),                      & ! (inout)
                       q(:), l(:,i,j), lq(:),                               & ! (out)
                       Nu(:,i,j), RHONu(:), Kh(:,i,j), RHOKh(:), Pr(:,i,j), & ! (out)
                       prod(:,i,j), prod_t(:), prod_q(:), prod_c(:),        & ! (out)
                       diss(:,i,j), diss_p(:),                              & ! (out)
                       sm25(:), smp(:), sh25(:), shpgh(:),                  & ! (out)
                       gammat(:), gammaq(:),                                & ! (out)
                       dudz2(:,i,j), n2_new(:), Ri(:,i,j),                  & ! (out)
                       dtldz(:), dqwdz(:),                                  & ! (out)
                       RHO(:), rho_h(:),                                    & ! (out)
                       uflux(:), vflux(:), tflux(:), qflux(:), eflux(:),    & ! (out)
                       Qlp(:,i,j), cldfrac(:,i,j), Zi(i,j), SFLX_BUOY(i,j), & ! (out)
                       U(:,i,j), V(:,i,j), W(:,i,j),                        & ! (in)
                       DENS(:,i,j), PRES(:,i,j),                            & ! (in)
                       POTT(:,i,j), POTL(:,i,j), POTV(:,i,j),               & ! (in)
                       Qw(:,i,j), N2(:,i,j),                                & ! (in)
                       EXNER(:,i,j), QDRY(:,i,j),                           & ! (in)
                       SFLX_PT, SFLX_SH(i,j), SFLX_QV(i,j), SFC_DENS(i,j),  & ! (in)
                       RLmo(i,j), us(i,j), ts(i,j), qs(i,j),                & ! (in)
                       z(:), CDZ(:), FDZ(:), F2H(:,:,i,j),                  & ! (in)
                       frac_land(i,j),                                      & ! (in)
                       dt,                                                  & ! (in)
                       I_B_TYPE,                                            & ! (in)
                       mynn_level3, .false.,                                & ! (in)
#ifdef _OPENACC
                       work(:,:),                                           & ! (work)
#endif
                       PTLV(:), Nu_f(:), Kh_f(:), q2_2(:), ac(:),           & ! (work)
                       betat(:), betaq(:),                                  & ! (work)
                       flx(:), a(:), b(:), c(:), d(:),                      & ! (work)
                       f_smp(:), f_shpgh(:), f_gamma(:),                    & ! (work)
                       tltv25(:), qwtv25(:), tvsq25(:),                     & ! (work)
                       tvsq_up(:), tvsq_lo(:), dtsq(:), dqsq(:), dcov(:),   & ! (work)
                       mflux(:),                                            & ! (work)
                       Uh(:), Vh(:), Wh(:), qw2(:), qh(:),                  & ! (work)
                       tlh(:), tvh(:), eh(:), dh(:),                        & ! (work)
                       TEML(:), LHVL(:), psat(:)                            ) ! (work)


       ! time integration

       flx(KS-1  ) = 0.0_RP
       flx(KE_PBL) = 0.0_RP

       do k = KS, KE_PBL-1
          rlqsm_h(k) = rho_h(k) * ( F2H(k,1,i,j) * lq(k+1) * smp(k+1) + F2H(k,2,i,j) * lq(k) * smp(k) )
       end do

       ! dens * u

       if ( ATMOS_PHY_BL_MYNN_MF ) then
          do k = KS, KE_PBL-1
             flx(k) = uflux(k)
          end do
       else
          ! countergradient flux
          do k = KS, KE_PBL-1
             flx(k) = - rlqsm_h(k) * ( U(k+1,i,j) - U(k,i,j) ) / FDZ(k)
          end do
       end if

       sf_t = SFLX_MU(i,j) / CDZ(KS)
       d(KS) = ( U(KS,i,j) * DENS(KS,i,j) + dt * ( sf_t - flx(KS) / CDZ(KS) ) ) / RHO(KS)
       do k = KS+1, KE_PBL
          d(k) = U(k,i,j) - dt * ( flx(k) - flx(k-1) ) / ( CDZ(k) * RHO(k) )
       end do
       c(KS) = 0.0_RP
       do k = KS, KE_PBL-1
          ap = - dt * RHONu(k) / FDZ(k)
          a(k) = ap / ( RHO(k) * CDZ(k) )
#ifdef _OPENACC
          if ( k==KS ) then
             b(k) = - a(k) + 1.0_RP
          else
             b(k) = - a(k) + dt * RHONu(k-1) / ( FDZ(k-1) * RHO(k) * CDZ(k) ) + 1.0_RP
          end if
#else
          b(k) = - a(k) - c(k) + 1.0_RP
#endif
          c(k+1) = ap / ( RHO(k+1) * CDZ(k+1) )
       end do
       a(KE_PBL) = 0.0_RP
       b(KE_PBL) = - c(KE_PBL) + 1.0_RP

       call MATRIX_SOLVER_tridiagonal( &
            KA, KS, KE_PBL, &
#ifdef _OPENACC
            work(:,:),              & ! (work)
#endif
            a(:), b(:), c(:), d(:), & ! (in)
            dummy(:)                ) ! (out)
!            phi_n(:)                ) ! (out)

       phi_n(KS:KE_PBL) = dummy(KS:KE_PBL)
       RHOU_t(KS,i,j) = ( phi_n(KS) * RHO(KS) - U(KS,i,j) * DENS(KS,i,j) ) / dt - sf_t
       do k = KS+1, KE_PBL
          RHOU_t(k,i,j) = ( phi_n(k) - U(k,i,j) ) * RHO(k) / dt
       end do
       do k = KE_PBL+1, KE
          RHOU_t(k,i,j) = 0.0_RP
       end do
       flxU(KS-1,i,j) = 0.0_RP
       flxU2(KS-1,i,j) = 0.0_RP
       do k = KS, KE_PBL-1
          flxU(k,i,j) = flx(k) &
                      - RHONu(k) * ( phi_n(k+1) - phi_n(k) ) / FDZ(k)
          flxU2(k,i,j) = flx(k)
       end do
       do k = KE_PBL, KE
          flxU(k,i,j) = 0.0_RP
          flxU2(k,i,j) = 0.0_RP
       end do


       ! dens * v

       if ( ATMOS_PHY_BL_MYNN_MF ) then
          do k = KS, KE_PBL-1
             flx(k) = vflux(k)
          end do
       else
          ! countergradient flux
          do k = KS, KE_PBL-1
             flx(k) = - rlqsm_h(k) * ( V(k+1,i,j) - V(k,i,j) ) / FDZ(k)
          end do
       end if

       sf_t = SFLX_MV(i,j) / CDZ(KS)
       d(KS) = ( V(KS,i,j) * DENS(KS,i,j) + dt * ( sf_t - flx(KS) / CDZ(KS) ) ) / RHO(KS)
       do k = KS+1, KE_PBL
          d(k) = V(k,i,j) - dt * ( flx(k) - flx(k-1) ) / ( CDZ(k) * RHO(k) )
       end do
       ! a,b,c is the same as those for the u

       call MATRIX_SOLVER_tridiagonal( &
            KA, KS, KE_PBL, &
#ifdef _OPENACC
            work(:,:),              & ! (work)
#endif
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
       flxV2(KS-1,i,j) = 0.0_RP
       do k = KS, KE_PBL-1
          flxV(k,i,j) = flx(k) &
                      - RHONu(k) * ( phi_n(k+1) - phi_n(k) ) / FDZ(k)
          flxV2(k,i,j) = flx(k)
       end do
       do k = KE_PBL, KE
          flxV(k,i,j) = 0.0_RP
          flxV2(k,i,j) = 0.0_RP
       end do


       ! dens * pott

       if ( ATMOS_PHY_BL_MYNN_MF ) then
          do k = KS, KE_PBL-1
             flx(k) = tflux(k)
          end do
       else
          ! countergradient flux
          do k = KS, KE_PBL-1
             flx(k) = - ( F2H(k,1,i,j) * lq(k+1) * gammat(k+1) + F2H(k,2,i,j) * lq(k) * gammat(k) ) &
                    * rho_h(k)
          end do
       end if

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
#ifdef _OPENACC
          if ( k==KS ) then
             b(k) = - a(k) + 1.0_RP
          else
             b(k) = - a(k) + dt * RHOKh(k-1) / ( FDZ(k-1) * RHO(k) * CDZ(k) ) + 1.0_RP
          end if
#else
          b(k) = - a(k) - c(k) + 1.0_RP
#endif
          c(k+1) = ap / ( RHO(k+1) * CDZ(k+1) )
       end do
       a(KE_PBL) = 0.0_RP
       b(KE_PBL) = - c(KE_PBL) + 1.0_RP

       call MATRIX_SOLVER_tridiagonal( &
            KA, KS, KE_PBL, &
#ifdef _OPENACC
            work(:,:),              & ! (work)
#endif
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
       flxT2(KS-1,i,j) = 0.0_RP
       do k = KS, KE_PBL-1
          flxT(k,i,j) = flx(k) &
                      - RHOKh(k) * ( phi_n(k+1) - phi_n(k) ) / FDZ(k)
          flxT2(k,i,j) = flx(k)
       end do
       do k = KE_PBL, KE
          flxT(k,i,j) = 0.0_RP
          flxT2(k,i,j) = 0.0_RP
       end do

       kmin = KS-1
       if ( flxT(KS,i,j) > 1E-4_RP ) then
          fmin = flxT(KS,i,j) / rho_h(KS)
          !$acc loop seq
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

       if ( ATMOS_PHY_BL_MYNN_MF ) then
          do k = KS, KE_PBL-1
             flx(k) = qflux(k)
          end do
       else
          ! countergradient flux
          do k = KS, KE_PBL-1
             flx(k) = - ( F2H(k,1,i,j) * lq(k+1) * gammaq(k+1) + F2H(k,2,i,j) * lq(k) * gammaq(k) ) &
                    * rho_h(k)
          end do
       end if

       sf_t = SFLX_QV(i,j) / CDZ(KS)
       d(KS) = ( QV(KS,i,j) * DENS(KS,i,j) + dt * ( sf_t - flx(KS) / CDZ(KS) ) ) / RHO(KS)
       do k = KS+1, KE_PBL
          d(k) = QV(k,i,j) - dt * ( flx(k) - flx(k-1) ) / ( CDZ(k) * RHO(k) )
       end do

       call MATRIX_SOLVER_tridiagonal( &
            KA, KS, KE_PBL, &
#ifdef _OPENACC
            work(:,:),              & ! (work)
#endif
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
       flxQ2(KS-1,i,j) = 0.0_RP
       do k = KS, KE_PBL-1
          flxQ(k,i,j) = flx(k) &
                      - RHOKh(k) * ( phi_n(k+1) - phi_n(k) ) / FDZ(k)
          flxQ2(k,i,j) = flx(k)
       end do
       do k = KE_PBL, KE
          flxQ(k,i,j) = 0.0_RP
          flxQ2(k,i,j) = 0.0_RP
       end do

       ! dens * TKE

!!$    if ( ATMOS_PHY_BL_MYNN_MF ) then
!!$       do k = KS, KE_PBL-1
!!$          flx(k) = eflux(k)
!!$       end do
!!$    else
       do k = KS, KE_PBL-1
          flx(k) = 0.0_RP
       end do
!!$    end if

       do k = KS, KE_PBL-1
          d(k) = tke(k) * DENS(k,i,j) / RHO(k) + dt * ( - ( flx(k) - flx(k-1) ) / ( CDZ(k) * RHO(k) ) + prod(k,i,j) )
       end do
       d(KE_PBL) = 0.0_RP

       c(KS) = 0.0_RP
       do k = KS, KE_PBL-1
          ap = - dt * ATMOS_PHY_BL_MYNN_Sq_fact * RHONu(k) / FDZ(k)
          a(k) = ap / ( RHO(k) * CDZ(k) )
#ifdef _OPENACC
          if ( k==KS ) then
             b(k) = - a(k) + 1.0_RP - diss(k,i,j) * dt
          else
             b(k) = - a(k) + dt * ATMOS_PHY_BL_MYNN_Sq_fact * RHONu(k-1) / ( FDZ(k-1) * RHO(k) * CDZ(k) ) + 1.0_RP - diss(k,i,j) * dt
          end if
#else
          b(k) = - a(k) - c(k) + 1.0_RP - diss(k,i,j) * dt
#endif
          c(k+1) = ap / ( RHO(k+1) * CDZ(k+1) )
       end do
       a(KE_PBL) = 0.0_RP

       call MATRIX_SOLVER_tridiagonal( &
            KA, KS, KE_PBL, &
#ifdef _OPENACC
            work(:,:),              & ! (work)
#endif
            a(:), b(:), c(:), d(:), & ! (in)
            phi_n(:)                ) ! (out)

       do k = KS, KE_PBL
          phi_n(k) = max( phi_n(k), ATMOS_PHY_BL_MYNN_TKE_MIN )
       end do

       do k = KS, KE_PBL
          diss(k,i,j) = diss(k,i,j) * phi_n(k)
          RPROG_t(k,i,j,I_TKE) = ( phi_n(k) * RHO(k) - PROG(k,i,j,I_TKE) * DENS(k,i,j) ) / dt
       end do
       !$acc loop independent
       do k = KE_PBL+1, KE
          RPROG_t(k,i,j,I_TKE) = - ATMOS_PHY_BL_MYNN_DUMP_COEF * PROG(k,i,j,I_TKE) * DENS(k,i,j) / dt
          diss(k,i,j) = 0.0_RP
          prod(k,i,j) = 0.0_RP
       end do


       if ( mynn_level3 ) then

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
#ifdef _OPENACC
             if ( k==KS ) then
                b(k) = - a(k) + diss_p(k)
             else
                b(k) = - a(k) + dt * RHONu(k-1) / ( FDZ(k-1) * RHO(k) * CDZ(k) ) + diss_p(k)
             end if
#else
             b(k) = - a(k) - c(k) + 1.0_RP + diss_p(k)
#endif
             c(k+1) = ap / ( RHO(k+1) * CDZ(k+1) )
          end do
          a(KE_PBL) = 0.0_RP
          b(KE_PBL) = - c(KE_PBL) + 1.0_RP + diss_p(KE_PBL)

          call MATRIX_SOLVER_tridiagonal( &
               KA, KS, KE_PBL, &
#ifdef _OPENACC
               work(:,:),              & ! (work)
#endif
               a(:), b(:), c(:), d(:), & ! (in)
               tsq(:)                  ) ! (out)

          do k = KS, KE_PBL
             tsq(k) = max( tsq(k), 0.0_RP )
          end do

          do k = KS, KE_PBL
             RPROG_t(k,i,j,I_TSQ) = ( tsq(k) * RHO(k) - PROG(k,i,j,I_TSQ) * DENS(k,i,j) ) / dt
          end do
          !$acc loop independent
          do k = KE_PBL+1, KE
             RPROG_t(k,i,j,I_TSQ) = - ATMOS_PHY_BL_MYNN_DUMP_COEF * PROG(k,i,j,I_TSQ) * DENS(k,i,j) / dt
          end do


          ! dens * qsq

          do k = KS, KE_PBL
             d(k) = qsq(k) * DENS(k,i,j) / RHO(k) + dt * prod_q(k)
          end do
          ! a, b, c are same as those for tsq

          call MATRIX_SOLVER_tridiagonal( &
               KA, KS, KE_PBL, &
#ifdef _OPENACC
               work(:,:),              & ! (work)
#endif
               a(:), b(:), c(:), d(:), & ! (in)
               qsq(:)                  ) ! (out)

          do k = KS, KE_PBL
             qsq(k) = max( qsq(k), 0.0_RP )
          end do

          do k = KS, KE_PBL
             RPROG_t(k,i,j,I_QSQ) = ( qsq(k) * RHO(k) - PROG(k,i,j,I_QSQ) * DENS(k,i,j) ) / dt
          end do
          !$acc loop independent
          do k = KE_PBL+1, KE
             RPROG_t(k,i,j,I_QSQ) = - ATMOS_PHY_BL_MYNN_DUMP_COEF * PROG(k,i,j,I_QSQ) * DENS(k,i,j) / dt
          end do


          ! dens * cov

          do k = KS, KE_PBL-1
             d(k) = cov(k) * DENS(k,i,j) / RHO(k) + dt * prod_c(k)
          end do
          d(KE_PBL) = 0.0_RP
          ! a, b, c are same as those for tsq

          call MATRIX_SOLVER_tridiagonal( &
               KA, KS, KE_PBL, &
#ifdef _OPENACC
               work(:,:),              & ! (work)
#endif
               a(:), b(:), c(:), d(:), & ! (in)
               cov(:)                  ) ! (out)

          do k = KS, KE_PBL
             cov(k) = sign( min( abs(cov(k)), sqrt(tsq(k)*qsq(k)) ), cov(k) )
          end do

          do k = KS, KE_PBL
             RPROG_t(k,i,j,I_COV) = ( cov(k) * RHO(k) - PROG(k,i,j,I_COV) * DENS(k,i,j) ) / dt
          end do
          !$acc loop independent
          do k = KE_PBL+1, KE
             RPROG_t(k,i,j,I_COV) = - ATMOS_PHY_BL_MYNN_DUMP_COEF * PROG(k,i,j,I_COV) * DENS(k,i,j) / dt
          end do

       end if

       Nu(KS-1,i,j) = 0.0_RP
       Kh(KS-1,i,j) = 0.0_RP
       do k = KE_PBL, KE
          Nu  (k,i,j) = 0.0_RP
          Kh  (k,i,j) = 0.0_RP
          Pr  (k,i,j) = 1.0_RP
          prod(k,i,j) = 0.0_RP
          diss(k,i,j) = 0.0_RP
       end do
       !$acc loop independent
       do k = KE_PBL+1, KE
          Ri     (k,i,j) = UNDEF
          dudz2  (k,i,j) = UNDEF
          l      (k,i,j) = UNDEF
          Qlp    (k,i,j) = UNDEF
          cldfrac(k,i,j) = UNDEF
       end do


    end do
    end do
    !$acc end kernels


    call FILE_HISTORY_query(HIST_Ri,     do_put)
    if ( do_put ) call FILE_HISTORY_put(HIST_Ri,     Ri(:,:,:))
    call FILE_HISTORY_query(HIST_Pr,     do_put)
    if ( do_put ) call FILE_HISTORY_put(HIST_Pr,     Pr(:,:,:))
    call FILE_HISTORY_query(HIST_TKE_pr, do_put)
    if ( do_put ) call FILE_HISTORY_put(HIST_TKE_pr, prod(:,:,:))
    call FILE_HISTORY_query(HIST_TKE_di, do_put)
    if ( do_put ) call FILE_HISTORY_put(HIST_TKE_di, diss(:,:,:))
    call FILE_HISTORY_query(HIST_dudz2,  do_put)
    if ( do_put ) call FILE_HISTORY_put(HIST_dudz2,  dudz2(:,:,:))
    call FILE_HISTORY_query(HIST_Lmix,   do_put)
    if ( do_put ) call FILE_HISTORY_put(HIST_Lmix,   l(:,:,:))

    call FILE_HISTORY_query(HIST_flxU, do_put)
    if ( do_put ) call FILE_HISTORY_put(HIST_flxU, flxU(1:,:,:))
    call FILE_HISTORY_query(HIST_flxV, do_put)
    if ( do_put ) call FILE_HISTORY_put(HIST_flxV, flxV(1:,:,:))
    call FILE_HISTORY_query(HIST_flxT, do_put)
    if ( do_put ) call FILE_HISTORY_put(HIST_flxT, flxT(1:,:,:))
    call FILE_HISTORY_query(HIST_flxQ, do_put)
    if ( do_put ) call FILE_HISTORY_put(HIST_flxQ, flxQ(1:,:,:))

    call FILE_HISTORY_query(HIST_flxU2, do_put)
    if ( do_put ) call FILE_HISTORY_put(HIST_flxU2, flxU2(1:,:,:))
    call FILE_HISTORY_query(HIST_flxV2, do_put)
    if ( do_put ) call FILE_HISTORY_put(HIST_flxV2, flxV2(1:,:,:))
    call FILE_HISTORY_query(HIST_flxT2, do_put)
    if ( do_put ) call FILE_HISTORY_put(HIST_flxT2, flxT2(1:,:,:))
    call FILE_HISTORY_query(HIST_flxQ2, do_put)
    if ( do_put ) call FILE_HISTORY_put(HIST_flxQ2, flxQ2(1:,:,:))

    !$acc end data

    return
  end subroutine ATMOS_PHY_BL_MYNN_tendency

  !-----------------------------------------------------------------------------
  ! private routines
  !-----------------------------------------------------------------------------
!OCL SERIAL
  subroutine MYNN_main( &
       KA, KS, KE_PBL, &
       i, j, &
       tke, tsq, qsq, cov, &
       q, l, lq, &
       Nu, RHONu, Kh, RHOKh, Pr, &
       prod, prod_t, prod_q, prod_c, &
       diss, diss_p, &
       sm25, smp, sh25, shpgh, &
       gammat, gammaq, &
       dudz2, n2_new, Ri, &
       dtldz, dqwdz, &
       RHO, RHO_h, &
       uflux, vflux, tflux, qflux, eflux, &
       Qlp, cldfrac, Zi, SFLX_BUOY, &
       U, V, W, &
       DENS, PRES, &
       POTT, POTL, POTV, &
       Qw, N2, &
       EXNER, QDRY, &
       SFLX_PT, SFLX_SH, SFLX_QV, SFC_DENS, &
       RLmo, us, ts, qs, &
       z, CDZ, FDZ, F2H, &
       frac_land, &
       dt, &
       I_B_TYPE, &
       mynn_level3, initialize, &
#ifdef _OPENACC
       work, &
#endif
       PTLV, Nu_f, Kh_f, q2_2, ac, betat, betaq, &
       flx, a, b, c, d, &
       f_smp, f_shpgh, f_gamma, tltv25, qwtv25, tvsq25, &
       tvsq_up, tvsq_lo, dtsq, dqsq, dcov, &
       mflux, &
       Uh, Vh, Wh, qw2, qh, tlh, tvh, eh, dh, &
       TEML, LHVL, psat )
    !$acc routine vector
    use scale_const, only: &
       EPS     => CONST_EPS,    &
       KARMAN  => CONST_KARMAN, &
       GRAV    => CONST_GRAV,   &
       EPSTvap => CONST_EPSTvap
    use scale_matrix, only: &
       MATRIX_SOLVER_tridiagonal
    integer, intent(in) :: KA, KS, KE_PBL
    integer, intent(in) :: i, j
    real(RP), intent(inout) :: tke(KA)
    real(RP), intent(inout) :: tsq(KA)
    real(RP), intent(inout) :: qsq(KA)
    real(RP), intent(inout) :: cov(KA)
    real(RP), intent(out) :: q(KA)
    real(RP), intent(out) :: l(KA)
    real(RP), intent(out) :: lq(KA)
    real(RP), intent(out) :: Nu(KA)
    real(RP), intent(out) :: RHONu(KA)
    real(RP), intent(out) :: Kh(KA)
    real(RP), intent(out) :: RHOKh(KA)
    real(RP), intent(out) :: Pr(KA)
    real(RP), intent(out) :: prod(KA)
    real(RP), intent(out) :: prod_t(KA)
    real(RP), intent(out) :: prod_q(KA)
    real(RP), intent(out) :: prod_c(KA)
    real(RP), intent(out) :: diss(KA)
    real(RP), intent(out) :: diss_p(KA)
    real(RP), intent(out) :: sm25(KA)
    real(RP), intent(out) :: smp(KA)
    real(RP), intent(out) :: sh25(KA)
    real(RP), intent(out) :: shpgh(KA)
    real(RP), intent(out) :: gammat(KA)
    real(RP), intent(out) :: gammaq(KA)
    real(RP), intent(out) :: dudz2(KA)
    real(RP), intent(out) :: n2_new(KA)
    real(RP), intent(out) :: Ri(KA)
    real(RP), intent(out) :: dtldz(KA)
    real(RP), intent(out) :: dqwdz(KA)
    real(RP), intent(out) :: RHO(KA)
    real(RP), intent(out) :: RHO_h(KA)
    real(RP), intent(out) :: uflux(KA)
    real(RP), intent(out) :: vflux(KA)
    real(RP), intent(out) :: tflux(KA)
    real(RP), intent(out) :: qflux(KA)
    real(RP), intent(out) :: eflux(KA)
    real(RP), intent(out) :: Qlp(KA)
    real(RP), intent(out) :: cldfrac(KA)
    real(RP), intent(out) :: Zi
    real(RP), intent(out) :: SFLX_BUOY
    real(RP), intent(in) :: U(KA)
    real(RP), intent(in) :: V(KA)
    real(RP), intent(in) :: W(KA)
    real(RP), intent(in) :: DENS(KA)
    real(RP), intent(in) :: PRES(KA)
    real(RP), intent(in) :: POTT(KA)
    real(RP), intent(in) :: POTL(KA)
    real(RP), intent(in) :: POTV(KA)
    real(RP), intent(in) :: Qw(KA)
    real(RP), intent(in) :: N2(KA)
    real(RP), intent(in) :: EXNER(KA)
    real(RP), intent(in) :: QDRY(KA)
    real(RP), intent(in) :: SFLX_PT
    real(RP), intent(in) :: SFLX_SH
    real(RP), intent(in) :: SFLX_QV
    real(RP), intent(in) :: SFC_DENS
    real(RP), intent(in) :: RLmo
    real(RP), intent(in) :: us
    real(RP), intent(in) :: ts
    real(RP), intent(in) :: qs
    real(RP), intent(in) :: z(KA)
    real(RP), intent(in) :: CDZ(KA)
    real(RP), intent(in) :: FDZ(KA)
    real(RP), intent(in) :: F2H(KA,2)
    real(RP), intent(in) :: frac_land
    real(RP), intent(in) :: dt
    integer,  intent(in) :: I_B_TYPE
    logical,  intent(in) :: mynn_level3
    logical,  intent(in) :: initialize

    ! work
    real(RP), intent(out) :: PTLV (KA)
    real(RP), intent(out) :: Nu_f (KA) !> Nu at the full level
    real(RP), intent(out) :: Kh_f (KA) !> Kh at the full level
    real(RP), intent(out) :: q2_2 (KA) !> q^2 for level 2
    real(RP), intent(out) :: ac   (KA) !> \alpha_c
    real(RP), intent(out) :: betat(KA) !> \beta_t
    real(RP), intent(out) :: betaq(KA) !> \beta_q

    real(RP), intent(out) :: flx(0:KA)
    real(RP), intent(out) :: a(KA)
    real(RP), intent(out) :: b(KA)
    real(RP), intent(out) :: c(KA)
    real(RP), intent(out) :: d(KA)

    real(RP) :: ap

    real(RP) :: zeta
    real(RP) :: phi_m, phi_h
    real(RP) :: us3

    ! for level 3 (work)
    real(RP), intent(out) :: f_smp  (KA) !> stability function for velocity for the countergradient (factor)
    real(RP), intent(out) :: f_shpgh(KA) !> stability function for scalars for the countergradient (factor)
    real(RP), intent(out) :: f_gamma(KA) !> - E_H / q^2 * GRAV / Theta_0
    real(RP), intent(out) :: tltv25 (KA)
    real(RP), intent(out) :: qwtv25 (KA)
    real(RP), intent(out) :: tvsq25 (KA)
    real(RP), intent(out) :: tvsq_up(KA) !> upper limit of <\theta_v^2> - <\theta_v^2>2.5
    real(RP), intent(out) :: tvsq_lo(KA) !> lower limit
    real(RP), intent(out) :: dtsq   (KA)
    real(RP), intent(out) :: dqsq   (KA)
    real(RP), intent(out) :: dcov   (KA)
    real(RP) :: tvsq
    real(RP) :: tltv
    real(RP) :: qwtv
    real(RP) :: wtl
    real(RP) :: wqw

    ! for O2019 (work)
    real(RP), intent(out) :: mflux(0:KA)

    ! work
    real(RP), intent(out) :: Uh(KA), Vh(KA), Wh(KA)
    real(RP), intent(out) :: qw2(KA), qh(KA)
    real(RP), intent(out) :: tlh(KA), tvh(KA)
    real(RP), intent(out) :: eh(KA), dh(KA)
    real(RP), intent(out) :: TEML(KA), LHVL(KA), psat(KA)

#ifdef _OPENACC
    real(RP), intent(out) :: work(KA,4)
#endif

    real(RP) :: sw, tmp

    integer :: k, it, nit

            
    if ( ATMOS_PHY_BL_MYNN_similarity .or. initialize ) then
       call get_phi( zeta, phi_m, phi_h, & ! (out)
                     z(KS), RLmo, I_B_TYPE ) ! (in)
    end if

    call calc_vertical_differece( KA, KS, KE_PBL, &
                                  i, j, &
                                  dudz2(:), dtldz(:), dqwdz(:), & ! (out)
                                  U(:), V(:), POTL(:),          & ! (in)
                                  Qw(:), QDRY(:),               & ! (in)
                                  CDZ(:), FDZ(:), F2H(:,:),     & ! (in)
                                  us, ts, qs,                   & ! (in)
                                  phi_m, phi_h, z(KS),          & ! (in)
                                  Uh(:), Vh(:), qw2(:), qh(:)   ) ! (work)

    us3 = us**3

    do k = KS, KE_PBL
       q(k) = sqrt( max( tke(k), ATMOS_PHY_BL_MYNN_TKE_MIN ) * 2.0_RP )
    end do

    if ( ATMOS_PHY_BL_MYNN_use_Zi ) then
       do k = KS, KE_PBL
          PTLV(k) = POTL(k) * ( 1.0_RP + EPSTvap * Qw(k) )
       end do
       call calc_zi_o2019( KA, KS, KE_PBL, &
                           Zi,              & ! (out)
                           PTLV(:), tke(:), & ! (in)
                           z(:), frac_land, & ! (in)
                           initialize )
    end if

    if ( initialize .or. (.not. mynn_level3) ) then
       ! estimate tsq, qsq, and cov

       do k = KS, KE_PBL
          n2_new(k) = min( max( N2(k), - ATMOS_PHY_BL_MYNN_N2_MAX ), ATMOS_PHY_BL_MYNN_N2_MAX )
          Ri(k) = n2_new(k) / dudz2(k)
       end do
       if ( initialize ) then
          dudz2(KS) = max( dudz2(KS), 1E-4_RP )
          n2_new(KS) = min( n2_new(KS), 0.0_RP )
          Ri(KS) = n2_new(KS) / dudz2(KS)
       end if

       SFLX_BUOY = - us3 * RLmo / KARMAN

       if ( ATMOS_PHY_BL_MYNN_MF ) then
          do k = KS, KE_PBL
             cldfrac(k) = 0.0_RP
          end do
          call calc_mflux( &
               KA, KS, KE_PBL, &
               DENS(:), POTV(:), POTL(:), Qw(:), & ! (in)
               U(:), V(:), W(:), tke(:),         & ! (in)
               cldfrac(:),                       & ! (in)
               EXNER(:),                         & ! (in)
               SFLX_SH, SFLX_BUOY,               & ! (in)
               SFLX_PT, SFLX_QV,                 & ! (in)
               SFC_DENS,                         & ! (in)
               Zi, Z(:), CDZ(:), F2H(:,:),       & ! (in)
               mflux(:),                         & ! (out)
               tflux(:), qflux(:),               & ! (out)
               uflux(:), vflux(:), eflux(:),     & ! (out)
               tlh(:), tvh(:), qh(:),            & ! (work)
               Uh(:), Vh(:), Wh(:), eh(:), dh(:) ) ! (work)
       end if

       ! length
       call get_length( &
            KA, KS, KE_PBL, &
            q(:), n2_new(:),      & ! (in)
            POTV(:), DENS(:),     & ! (in)
            mflux(:),             & ! (in)
            Zi,                   & ! (in)
            SFLX_BUOY, RLmo,      & ! (in)
            z(:), CDZ(:), FDZ(:), & ! (in)
            initialize,           & ! (in)
            l(:)                  ) ! (out)

       call get_q2_level2( &
            KA, KS, KE_PBL, &
            dudz2(:), Ri(:), l(:), & ! (in)
            q2_2(:)                ) ! (out)

       if ( initialize ) then
          do k = KS, KE_PBL
             q(k) = sqrt( q2_2(k) )
          end do
       end if

       do k = KS, KE_PBL
          ac(k) = min( q(k) / sqrt( q2_2(k) + 1e-20_RP ), 1.0_RP )
       end do

       call get_smsh( &
            KA, KS, KE_PBL, &
            i, j, &
            q(:), ac(:),            & ! (in)
            l(:), n2_new(:),        & ! (in)
            Ri(:),                  & ! (in)
            POTV(:), dudz2(:),      & ! (in)
            dtldz(:), dqwdz(:),     & ! (in)
            betat(:), betaq(:),     & ! (in) ! dummy
            .false., .false.,       & ! (in)
            tsq(:), qsq(:), cov(:), & ! (inout)
            sm25(:), f_smp(:),      & ! (out) ! dummy
            sh25(:), f_shpgh(:),    & ! (out) ! dymmy
            f_gamma(:),             & ! (out) ! dummy
            tltv25(:), qwtv25(:),   & ! (out) ! dummy
            tvsq25(:),              & ! (out) ! dummy
            tvsq_up(:), tvsq_lo(:)  ) ! (out) ! dummy

    end if

    flx(KS-1  ) = 0.0_RP
    flx(KE_PBL) = 0.0_RP

    if ( initialize ) then
       nit = KE_PBL - 1
    else
       nit = 1
    end if

    !$acc loop seq
    do it = 1, nit

       call partial_condensation( KA, KS, KE_PBL, &
                                  betat(:), betaq(:),       & ! (out)
                                  Qlp(:), cldfrac(:),       & ! (out)
                                  PRES(:), POTT(:),         & ! (in)
                                  POTL(:), Qw(:),           & ! (in)
                                  EXNER(:),                 & ! (in)
                                  tsq(:), qsq(:), cov(:),   & ! (in)
                                  TEML(:), LHVL(:), psat(:) ) ! (work)

       ! update N2
       do k = KS, KE_PBL
          n2_new(k) = min(ATMOS_PHY_BL_MYNN_N2_MAX, &
                          GRAV * ( dtldz(k) * betat(k) + dqwdz(k) * betaq(k) ) / POTV(k) )
          Ri(k) = n2_new(k) / dudz2(k)
       end do

       SFLX_BUOY = GRAV / POTV(KS) * ( betat(KS) * SFLX_PT + betaq(KS) * SFLX_QV ) / SFC_DENS


       if ( ATMOS_PHY_BL_MYNN_use_Zi ) then
          do k = KS, KE_PBL
             PTLV(k) = POTL(k) * betat(k) + Qw(k) * betaq(k)
          end do
          call calc_zi_o2019( KA, KS, KE_PBL, &
                              Zi,              & ! (out)
                              PTLV(:), tke(:), & ! (in)
                              z(:), frac_land, & ! (in)
                              .false.          ) ! (in)
       end if

       if ( ATMOS_PHY_BL_MYNN_MF ) then
          call calc_mflux( &
               KA, KS, KE_PBL, &
               DENS(:), POTV(:), POTL(:), Qw(:), & ! (in)
               U(:), V(:), W(:), tke(:),         & ! (in)
               cldfrac(:),                       & ! (in)
               EXNER(:),                         & ! (in)
               SFLX_SH, SFLX_BUOY,               & ! (in)
               SFLX_PT, SFLX_QV,                 & ! (in)
               SFC_DENS,                         & ! (in)
               Zi, Z(:), CDZ(:), F2H(:,:),       & ! (in)
               mflux(:),                         & ! (out)
               tflux(:), qflux(:),               & ! (out)
               uflux(:), vflux(:), eflux(:),     & ! (out)
               tlh(:), tvh(:), qh(:),            & ! (work)
               Uh(:), Vh(:), Wh(:), eh(:), dh(:) ) ! (work)
       end if


       ! length
       call get_length( &
            KA, KS, KE_PBL, &
            q(:), n2_new(:),      & ! (in)
            POTV(:), DENS(:),     & ! (in)
            mflux(:),             & ! (in)
            Zi,                   & ! (in)
            SFLX_BUOY, RLmo,      & ! (in)
            z(:), CDZ(:), FDZ(:), & ! (in)
            .false.,              & ! (in)
            l(:)                  ) ! (out)

       call get_q2_level2( &
            KA, KS, KE_PBL, &
            dudz2(:), Ri(:), l(:), & ! (in)
            q2_2(:)                ) ! (out)

       do k = KS, KE_PBL
          ac(k) = min( q(k) / sqrt( q2_2(k) + 1e-20_RP ), 1.0_RP )
       end do

       call get_smsh( &
            KA, KS, KE_PBL, & 
            i, j,                   & ! (in)
            q(:), ac(:),            & ! (in)
            l(:), n2_new(:),        & ! (in)
            Ri(:),                  & ! (in)
            POTV(:), dudz2(:),      & ! (in)
            dtldz(:), dqwdz(:),     & ! (in)
            betat(:), betaq(:),     & ! (in)
            mynn_level3,            & ! (in)
            initialize .and. it==1, & ! (in)
            tsq(:), qsq(:), cov(:), & ! (inout)
            sm25(:), f_smp(:),      & ! (out)
            sh25(:), f_shpgh(:),    & ! (out)
            f_gamma(:),             & ! (out)
            tltv25(:), qwtv25(:),   & ! (out)
            tvsq25(:),              & ! (out)
            tvsq_up(:), tvsq_lo(:)  ) ! (out)

       do k = KS, KE_PBL
          lq(k) = l(k) * q(k)
          Nu_f(k) = lq(k) * sm25(k)
          Kh_f(k) = lq(k) * sh25(k)
       end do
       if ( ATMOS_PHY_BL_MYNN_similarity ) then
          Nu_f(KS) = KARMAN * z(KS) * us / phi_m
          Kh_f(KS) = KARMAN * z(KS) * us / phi_h
       end if

       do k = KS, KE_PBL-1
          Nu(k) = min( F2H(k,1) * Nu_f(k+1) + F2H(k,2) * Nu_f(k), &
                       ATMOS_PHY_BL_MYNN_NU_MAX )
          Kh(k) = min( F2H(k,1) * Kh_f(k+1) + F2H(k,2) * Kh_f(k), &
                        ATMOS_PHY_BL_MYNN_KH_MAX )
       end do

       do k = KS, KE_PBL-1
          sw = 0.5_RP - sign(0.5_RP, abs(Kh(k)) - EPS)
          Pr(k) = Nu(k) / ( Kh(k) + sw ) * ( 1.0_RP - sw ) &
                + 1.0_RP * sw
       end do

       RHO(KS) = DENS(KS) + dt * SFLX_QV / CDZ(KS)
       do k = KS+1, KE_PBL
          RHO(k) = DENS(k)
       end do

       do k = KS, KE_PBL-1
          rho_h(k) = F2H(k,1) * RHO(k+1) + F2H(k,2) * RHO(k)
          RHONu(k) = max( Nu(k) * rho_h(k), 1e-20_RP )
          RHOKh(k) = max( Kh(k) * rho_h(k), 1e-20_RP )
!          RHONu(k) = F2H(k,1) * RHO(k+1) * Nu_f(k+1) + F2H(k,2) * RHO(k) * Nu_f(k)
!          RHOKh(k) = F2H(k,1) * RHO(k+1) * Kh_f(k+1) + F2H(k,2) * RHO(k) * Kh_f(k)
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

          if ( ATMOS_PHY_BL_MYNN_similarity .or. initialize ) then
             tmp = 2.0_RP * us * phi_h / ( KARMAN * z(KS) )
             tmp = tmp * ( zeta / ( z(KS) * RLmo ) )**2 ! correspoinding to the limitter for zeta
             ! TSQ
             prod_t(KS) = tmp * ts**2
             ! QSQ
             prod_q(KS) = tmp * ts * qs
             ! COV
             prod_c(KS) = tmp * qs**2
          end if

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
               dtldz(:), dqwdz(:), POTV(:),     & ! (in)
               prod_t(:), prod_q(:), prod_c(:), & ! (in)
               betat(:), betaq(:),              & ! (in)
               f_gamma(:), l(:), q(:),          & ! (in)
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

          if ( ATMOS_PHY_BL_MYNN_similarity .or. initialize ) then
             smp(KS)    = 0.0_RP
             shpgh(KS)  = 0.0_RP
             gammat(KS) = 0.0_RP
             gammaq(KS) = 0.0_RP

             tmp = 2.0_RP * us * phi_h / ( KARMAN * z(KS) )
             tmp = tmp * ( zeta / ( z(KS) * RLmo ) )**2 ! correspoinding to the limitter for zeta
             ! TSQ
             prod_t(KS) = tmp * ts**2
             ! QSQ
             prod_q(KS) = tmp * ts * qs
             ! COV
             prod_c(KS) = tmp * qs**2
          end if

       else

          do k = KS, KE_PBL
             smp   (k) = 0.0_RP
             shpgh (k) = 0.0_RP
             gammat(k) = 0.0_RP
             gammaq(k) = 0.0_RP
          end do

       end if


       ! production of TKE
       do k = KS, KE_PBL
          prod(k) = lq(k) * ( ( sm25(k) + smp(k) ) * dudz2(k) &
                            - ( sh25(k) * n2_new(k) - shpgh(k) ) )
       end do
       if ( ATMOS_PHY_BL_MYNN_similarity .or. initialize ) then
          prod(KS) = us3 / ( KARMAN * z(KS) ) * ( phi_m - zeta )
       end if

       do k = KS, KE_PBL
          diss(k) = - 2.0_RP * q(k) / ( B1 * l(k) )
!          prod(k) = max( prod(k), - tke(k) / dt - diss(k) * tke(k) )
          diss_p(k) = dt * 2.0_RP * q(k) / ( B2 * l(k) )
       end do


       if ( .not. initialize ) exit

       ! dens * TKE

!!$    if ( ATMOS_PHY_BL_MYNN_MF ) then
!!$       do k = KS, KE_PBL-1
!!$          flx(k) = eflux(k)
!!$       end do
!!$    else
       do k = KS, KE_PBL-1
          flx(k) = 0.0_RP
       end do
!!$    end if

       do k = KS, KE_PBL-1
          d(k) = dt * ( - ( flx(k) - flx(k-1) ) / ( CDZ(k) * RHO(k) ) + prod(k) )
       end do
       d(KE_PBL) = 0.0_RP

       c(KS) = 0.0_RP
       do k = KS, KE_PBL-1
          ap = - dt * ATMOS_PHY_BL_MYNN_Sq_fact * RHONu(k) / FDZ(k)
          a(k) = ap / ( RHO(k) * CDZ(k) )
#ifdef _OPENACC
          if ( k==KS ) then
             b(k) = - a(k) - diss(k) * dt
          else
             b(k) = - a(k) + dt * ATMOS_PHY_BL_MYNN_Sq_fact * RHONu(k-1) / ( FDZ(k-1) * RHO(k) * CDZ(k) ) - diss(k) * dt
          end if
#else
          b(k) = - a(k) - c(k) - diss(k) * dt
#endif
          c(k+1) = ap / ( RHO(k+1) * CDZ(k+1) )
       end do
       a(KE_PBL) = 0.0_RP
       b(KE_PBL) = - c(KE_PBL) - diss(KE_PBL) * dt

       call MATRIX_SOLVER_tridiagonal( &
            KA, KS, KE_PBL, &
#ifdef _OPENACC
            work(:,:),              & ! (work)
#endif
            a(:), b(:), c(:), d(:), & ! (in)
            tke(:)                ) ! (out)

       do k = KS, KE_PBL-1
          tke(k) = min( max( tke(k), ATMOS_PHY_BL_MYNN_TKE_MIN ), 100.0_RP )
       end do
       tke(KE_PBL) = 0.0_RP


       do k = KS, KE_PBL
          q(k) = ( q(k) + sqrt( tke(k) * 2.0_RP ) ) * 0.5_RP ! to avoid oscillation
       end do


       if ( .not. mynn_level3 ) cycle

       ! dens * tsq

       do k = KS, KE_PBL-1
          d(k) = dt * prod_t(k)
       end do
       d(KE_PBL) = 0.0_RP
       c(KS) = 0.0_RP
       do k = KS, KE_PBL-1
          ap = - dt * RHONu(k) / FDZ(k)
          a(k) = ap / ( RHO(k) * CDZ(k) )
#ifdef _OPENACC
          if ( k==KS ) then
             b(k) = - a(k) + diss_p(k)
          else
             b(k) = - a(k) + dt * RHONu(k-1) / ( FDZ(k-1) * RHO(k) * CDZ(k) ) + diss_p(k)
          end if
#else
          b(k) = - a(k) - c(k) + 1.0_RP + diss_p(k)
#endif
          c(k+1) = ap / ( RHO(k+1) * CDZ(k+1) )
       end do
       a(KE_PBL) = 0.0_RP
       b(KE_PBL) = - c(KE_PBL) + diss_p(KE_PBL)

       call MATRIX_SOLVER_tridiagonal( &
            KA, KS, KE_PBL, &
#ifdef _OPENACC
            work(:,:),              & ! (work)
#endif
            a(:), b(:), c(:), d(:), & ! (in)
            tsq(:)                  ) ! (out)

       do k = KS, KE_PBL-1
          tsq(k) = max( tsq(k), 0.0_RP )
       end do
       tsq(KE_PBL) = 0.0_RP


       ! dens * qsq

       do k = KS, KE_PBL-1
          d(k) = dt * prod_q(k)
       end do
       d(KE_PBL) = 0.0_RP
       ! a, b, c are same as those for tsq

       call MATRIX_SOLVER_tridiagonal( &
            KA, KS, KE_PBL, &
#ifdef _OPENACC
            work(:,:),              & ! (work)
#endif
            a(:), b(:), c(:), d(:), & ! (in)
            qsq(:)                  ) ! (out)

       do k = KS, KE_PBL-1
          qsq(k) = max( qsq(k), 0.0_RP )
       end do
       qsq(KE_PBL) = 0.0_RP


       ! dens * cov

       do k = KS, KE_PBL-1
          d(k) = dt * prod_c(k)
       end do
       d(KE_PBL) = 0.0_RP
       ! a, b, c are same as those for tsq

       call MATRIX_SOLVER_tridiagonal( &
            KA, KS, KE_PBL, &
#ifdef _OPENACC
            work(:,:),              & ! (work)
#endif
            a(:), b(:), c(:), d(:), & ! (in)
            cov(:)                  ) ! (out)

       do k = KS, KE_PBL-1
          cov(k) = sign( min( abs(cov(k)), sqrt(tsq(k)*qsq(k)) ), cov(k) )
       end do
       cov(KE_PBL) = 0.0_RP

    end do

    return
  end subroutine MYNN_main

!OCL SERIAL
  subroutine get_length( &
       KA, KS, KE_PBL, &
       q, n2,           &
       POTV, DENS,      &
       mflux_gl,        &
       Zi,              &
       SFLX_BUOY, RLmo, &
       z, CDZ, FDZ,     &
       initialize,      &
       l                )
    !$acc routine vector
    use scale_const, only: &
!       GRAV   => CONST_GRAV, &
       KARMAN => CONST_KARMAN, &
       EPS    => CONST_EPS
    implicit none
    integer,  intent(in), value :: KA, KS, KE_PBL

    real(RP), intent(in) :: q(KA)
    real(RP), intent(in) :: n2(KA)
    real(RP), intent(in) :: POTV(KA)
    real(RP), intent(in) :: DENS(KA)
    real(RP), intent(in) :: mflux_gl(0:KA)
    real(RP), intent(in) :: Zi
    real(RP), intent(in) :: SFLX_BUOY !> g/T0 <w Tv> @ surface
    real(RP), intent(in) :: RLmo      !> inverse of Obukhov length
    real(RP), intent(in) :: z(KA)
    real(RP), intent(in) :: CDZ(KA)
    real(RP), intent(in) :: FDZ(KA)
    logical,  intent(in) :: initialize

    real(RP), intent(out) :: l(KA)

    real(RP), parameter :: sqrt05 = sqrt( 0.5_RP )

    real(RP) :: ls     !> L_S
    real(RP) :: lb     !> L_B
    real(RP) :: lt     !> L_T
    real(RP) :: rlt    !> 1/L_T

    real(RP) :: qc     !> q_c
    real(RP) :: int_q  !> \int q dz
    real(RP) :: int_qz !> \int qz dz
    real(RP) :: rn2sr  !> 1/N
    real(RP) :: zeta   !> height normalized by the Obukhov length

    real(RP) :: qdz

    ! for Olson et al. (2019)
!    real(RP) :: lbl, lup, ldown
    real(RP) :: tau

    real(RP) :: sw
    integer :: kmax
    integer :: k

    if ( ATMOS_PHY_BL_MYNN_use_Zi ) then
       kmax = KE_PBL
       !$acc loop seq
       do k = KS, KE_PBL
          if ( z(k) > Zi * 1.3_RP ) then
             kmax = k - 1
             exit
          end if
       end do
    else
       kmax = KE_PBL
    end if

    int_qz = 0.0_RP
    int_q = 0.0_RP
    !$acc loop private(qdz) reduction(+:int_qz,int_q)
    do k = KS, kmax
       qdz = q(k) * CDZ(k)
       int_qz = int_qz + z(k) * qdz
       int_q  = int_q + qdz
    end do
    ! LT
    lt = min( max(0.23_RP * int_qz / (int_q + 1e-20_RP), &
                  LT_min), &
                  ATMOS_PHY_BL_MYNN_Lt_MAX )
    rlt = 1.0_RP / lt

    qc = ( lt * max(SFLX_BUOY,0.0_RP) )**OneOverThree ! qc=0 if SFLX_BUOY<0
    if ( ATMOS_PHY_BL_MYNN_O2019 ) then
       tau = 0.5_RP * zi / max(SFLX_BUOY,1.0E-10_RP)**OneOverThree
    end if

    !$acc loop private(zeta,sw,ls,rn2sr,tau,lb)
    do k = KS, KE_PBL
       zeta = z(k) * RLmo

       ! LS
       sw = sign(0.5_RP, zeta) + 0.5_RP ! 1 for zeta >= 0, 0 for zeta < 0
       ls = KARMAN * z(k) &
          * ( sw / (1.0_RP + ATMOS_PHY_BL_MYNN_cns*zeta*sw ) &
            + ( (1.0_RP - ATMOS_PHY_BL_MYNN_alpha4*zeta)*(1.0_RP-sw) )**0.2_RP )

       ! LB
       sw  = sign(0.5_RP, n2(k)-EPS) + 0.5_RP ! 1 for dptdz >0, 0 for dptdz <= 0
       rn2sr = 1.0_RP / ( sqrt(n2(k)*sw) + 1.0_RP-sw)
       if ( ATMOS_PHY_BL_MYNN_O2019 ) then
          if ( z(k) > zi ) tau = 50.0_RP
!!$          qtmp = q(k)**2 * 0.5_RP * POTV(k) / GRAV
!!$          lup = z(KE_PBL) - z(k)
!!$          do kk = k+1, KE_PBL
!!$             qtmp = qtmp - ( POTV(kk) - POTV(k) ) * FDZ(kk-1)
!!$             if ( qtmp < 0.0_RP ) then
!!$                lup = z(kk) - z(k)
!!$                exit
!!$             end if
!!$          end do
!!$          qtmp = q(k)**2 * 0.5_RP * POTV(k) / GRAV
!!$          ldown = z(k)
!!$          do kk = k-1, KS
!!$             qtmp = qtmp - ( POTV(k) - POTV(kk) ) * FDZ(kk)
!!$             if ( qtmp < 0.0_RP ) then
!!$                ldown = z(k) - z(kk)
!!$                exit
!!$             end if
!!$          end do
!!$          lbl = sqrt( lup**2 + ldown**2 )
!!$          we = 0.5_RP * tanh( ( z(k) - zi * 1.3_RP ) /  ( zi * 0.15_RP ) ) + 0.5_RP
!!$          lb = lb * ( 1.0_RP - we ) + lbl * we
          lb = ATMOS_PHY_BL_MYNN_alpha2 * max( q(k), ( mflux_gl(k)+mflux_gl(k-1))/DENS(k) ) * rn2sr * sw &
             + tau * q(k) * sqrt05 * (1.0_RP-sw)
       else
          lb = ATMOS_PHY_BL_MYNN_alpha2 * (1.0_RP + 5.0_RP * sqrt(qc*rn2sr*rlt)) * q(k) * rn2sr * sw & ! qc=0 when RLmo > 0
             + 1.E10_RP * (1.0_RP-sw)
       end if

       ! L
       if ( initialize ) then
          l(k) = min( ls, lb )
       else if ( ATMOS_PHY_BL_MYNN_O2019 ) then
          l(k) = min( 1.0_RP / ( 1.0_RP/ls + rlt ), lb )
       else
          l(k) = 1.0_RP / ( 1.0_RP/ls + rlt + 1.0_RP/(lb+1E-20_RP) )
       end if
    end do

    return
  end subroutine get_length

!OCL SERIAL
  subroutine get_q2_level2( &
       KA, KS, KE_PBL, &
       dudz2, Ri, l, &
       q2_2          )
    !$acc routine vector
    implicit none
    integer,  intent(in), value  :: KA, KS, KE_PBL

    real(RP), intent(in)  :: dudz2(KA)
    real(RP), intent(in)  :: Ri(KA)
    real(RP), intent(in)  :: l(KA)

    real(RP), intent(out) :: q2_2(KA)

    real(RP) :: rf   !> Rf
    real(RP) :: sm_2 !> sm for level 2
    real(RP) :: sh_2 !> sh for level 2

    ! for K2010
    real(RP) :: A2k, F1, Rf1, AF12

    integer :: k

    do k = KS, KE_PBL
       if ( ATMOS_PHY_BL_MYNN_K2010 .and. Ri(k) > 0.0_RP ) then
          A2k = A2 / ( 1.0_RP + Ri(k) )
       else
          A2k = A2
       end if
       F1 = B1 * ( G1 - C1 ) + 2.0 * A1 * ( 3.0_RP - 2.0_RP * C2 ) + 3.0 * A2k * ( 1.0_RP - C2 ) * ( 1.0_RP - C5 )
       Rf1 = B1 * ( G1 - C1 ) / F1
       AF12 = A1 * F1 / ( A2k * F2 )
       rf = min(0.5_RP / AF12 * ( Ri(k) &
                                + AF12*Rf1 &
                                - sqrt(Ri(k)**2 + 2.0_RP*AF12*(Rf1-2.0_RP*Rf2)*Ri(k) + (AF12*Rf1)**2) ), &
                Rfc)
       sh_2 = 3.0_RP * A2k * (G1+G2) * (Rfc-rf) / (1.0_RP-rf)
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
       l, n2, Ri,       &
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
       tvsq25, tvsq_up, &
       tvsq_lo          )
    !$acc routine vector
    use scale_const, only: &
       EPS  => CONST_EPS, &
       HUGE => CONST_HUGE, &
       GRAV => CONST_GRAV
    implicit none
    integer,  intent(in), value  :: KA, KS, KE_PBL
    integer,  intent(in), value  :: i, j ! for debug

    real(RP), intent(in)  :: q(KA)
    real(RP), intent(in)  :: ac(KA)
    real(RP), intent(in)  :: l(KA)
    real(RP), intent(in)  :: n2(KA)
    real(RP), intent(in)  :: Ri(KA)
    real(RP), intent(in)  :: potv(KA)
    real(RP), intent(in)  :: dudz2(KA)
    real(RP), intent(in)  :: dtldz(KA)
    real(RP), intent(in)  :: dqwdz(KA)
    real(RP), intent(in)  :: betat(KA)
    real(RP), intent(in)  :: betaq(KA)
    logical,  intent(in), value :: mynn_level3
    logical,  intent(in), value :: initialize

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

    ! for K2010
    real(RP) :: A2k

    integer :: k

    ! level 2.5
    do k = KS, KE_PBL

       if ( ATMOS_PHY_BL_MYNN_K2010 .and. Ri(k) > 0.0_RP ) then
          A2k = A2 / ( 1.0_RP + Ri(k) )
       else
          A2k = A2
       end if

       ac2 = ac(k)**2
       l2 = l(k)**2
       q2 = q(k)**2

       f1 = -  3.0_RP * ac2 * A2k * B2  * ( 1.0_RP - C3 )
       f2 = -  9.0_RP * ac2 * A1  * A2k * ( 1.0_RP - C2 )
       f3 =    9.0_RP * ac2 * A2k**2    * ( 1.0_RP - C2 ) * ( 1.0_RP - C5 )
       f4 = - 12.0_RP * ac2 * A1  * A2k * ( 1.0_RP - C2 )

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

       sm25(k) = ac(k) * A1  * ( p3q2 - 3.0_RP * C1 * p4q2 ) * rd25q2
       sh25(k) = ac(k) * A2k * ( p2q2 + 3.0_RP * C1 * p5q2 ) * rd25q2

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

             l2q2 = l2 / q2
             l2q2 = min( l2q2, 1.0_RP / max(n2(k), EPS) )
             q2 = l2 / l2q2

             rdpq2  = q2 / max( p2q2 * ( f4 * ghq2 + q2 ) + p5q2 * ( f3 * ghq2 + q2 ), 1e-20_RP )

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

             emq2 = 3.0_RP * ac(k) * A1  * ( 1.0_RP - C3 ) * ( f3 - f4 ) * rdpq2 ! (p3-p4)/gh = (f3-f4)
             eh   = 3.0_RP * ac(k) * A2k * ( 1.0_RP - C3 ) * ( p2q2 + p5q2 ) * rdpq2

!             q2 = l2 / l2q2
             fact = GRAV / POTV(k)
             f_smp  (k) = emq2 * fact**2 * l2q2
             f_shpgh(k) = eh   * fact**2 / q2
             f_gamma(k) = - eh * fact / q2
          end if
       else ! level 2.5
          tvsq_up(k) = 0.0_RP
          tvsq_lo(k) = 0.0_RP
          f_smp(k)   = 0.0_RP
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
       dtsq, dqsq, dcov )
    !$acc routine vector
    use scale_const, only: &
       GRAV => CONST_GRAV
    implicit none
    integer, intent(in), value :: KA, KS, KE
    integer, intent(in), value :: i, j ! for debug

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
    real(RP), intent(in), value :: dt

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
       dudz2, dtldz, dqwdz,  &
       U, V, POTL, Qw, QDRY, &
       CDZ, FDZ, F2H,        &
       us, ts, qs,           &
       phi_m, phi_h, z1,     &
       Uh, Vh, qw2, qh       )
    !$acc routine vector
    use scale_const, only: &
       KARMAN => CONST_KARMAN
    integer,  intent(in), value  :: KA, KS, KE
    integer,  intent(in), value  :: i, j  ! for debug

    real(RP), intent(out) :: dudz2(KA)
    real(RP), intent(out) :: dtldz(KA)
    real(RP), intent(out) :: dqwdz(KA)

    real(RP), intent(in) :: U   (KA)
    real(RP), intent(in) :: V   (KA)
    real(RP), intent(in) :: POTL(KA)
    real(RP), intent(in) :: Qw  (KA)
    real(RP), intent(in) :: QDRY(KA)
    real(RP), intent(in) :: CDZ (KA)
    real(RP), intent(in) :: FDZ (KA)
    real(RP), intent(in) :: F2H (KA,2)
    real(RP), intent(in) :: us
    real(RP), intent(in) :: ts
    real(RP), intent(in) :: qs
    real(RP), intent(in) :: phi_m
    real(RP), intent(in) :: phi_h
    real(RP), intent(in) :: z1

    real(RP), intent(out) :: Uh (KA) ! work
    real(RP), intent(out) :: Vh (KA) ! work
    real(RP), intent(out) :: qw2(KA) ! work
    real(RP), intent(out) :: qh (KA) ! work

    integer :: k

    do k = KS, KE
       Uh(k) = f2h(k,1) * U(k+1) + f2h(k,2) * U(k)
       Vh(k) = f2h(k,1) * V(k+1) + f2h(k,2) * V(k)
    end do

    if ( ATMOS_PHY_BL_MYNN_dz_sim ) then
       dudz2(KS) = ( us * phi_m / ( KARMAN * z1 ) )**2
    else
       dudz2(KS) = ( ( Uh(KS) - U(KS) )**2 + ( Vh(KS) - V(KS) )**2 ) / ( CDZ(KS) * 0.5_RP )**2
!       dudz2(KS) = ( ( Uh(KS) )**2 + ( Vh(KS) )**2 ) / CDZ(KS)**2
    end if
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
    if ( ATMOS_PHY_BL_MYNN_dz_sim ) then
       dtldz(KS) = ts * phi_h / ( KARMAN * z1 )
    else
       dtldz(KS) = ( qh(KS) - POTL(KS) ) / ( CDZ(KS) * 0.5_RP )
!       dtldz(KS) = ( POTL(KS+1) - POTL(KS) ) / FDZ(KS)
    end if
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
    if ( ATMOS_PHY_BL_MYNN_dz_sim ) then
       dqwdz(KS) = qs * phi_h / ( KARMAN * z1 )
    else
       dqwdz(KS) = ( qh(KS) - qw2(KS) ) / ( CDZ(KS) * 0.5_RP )
!       dqwdz(KS) = ( qh(KS) - Qw(KS) ) / ( CDZ(KS) * 0.5_RP )
!       dqwdz(KS) = ( Qw(KS+1) - Qw(KS) ) / FDZ(KS)
    end if
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
       betat, betaq, Qlp, cldfrac,  &
       PRES, POTT, POTL, Qw, EXNER, &
       tsq, qsq, cov,               &
       TEML, LHVL, psat             )
    use scale_const, only: &
       CPdry   => CONST_CPdry,  &
       Rvap    => CONST_Rvap,   &
       EPSvap  => CONST_EPSvap, &
       EPSTvap => CONST_EPSTvap
    !$acc routine vector
    use scale_atmos_saturation, only: &
       ATMOS_SATURATION_psat => ATMOS_SATURATION_psat_liq
!       ATMOS_SATURATION_pres2qsat => ATMOS_SATURATION_pres2qsat_liq
    use scale_atmos_hydrometeor, only: &
       HYDROMETEOR_LHV => ATMOS_HYDROMETEOR_LHV, &
       CP_VAPOR
    integer,  intent(in), value  :: KA, KS, KE
    real(RP), intent(in)  :: PRES (KA)
    real(RP), intent(in)  :: POTT (KA)
    real(RP), intent(in)  :: POTL (KA)
    real(RP), intent(in)  :: Qw   (KA)
    real(RP), intent(in)  :: EXNER(KA)
    real(RP), intent(in)  :: tsq  (KA)
    real(RP), intent(in)  :: qsq  (KA)
    real(RP), intent(in)  :: cov  (KA)

    real(RP), intent(out) :: betat  (KA)
    real(RP), intent(out) :: betaq  (KA)
    real(RP), intent(out) :: Qlp    (KA)
    real(RP), intent(out) :: cldfrac(KA)

    real(RP), intent(out) :: TEML(KA)
    real(RP), intent(out) :: LHVL(KA)
    real(RP), intent(out) :: psat(KA)

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

  subroutine get_phi( &
       zeta, phi_m, phi_h, &
       z1, RLmo, I_B_TYPE )
    !$acc routine seq
    real(RP), intent(out) :: zeta
    real(RP), intent(out) :: phi_m
    real(RP), intent(out) :: phi_h
    real(RP), intent(in)  :: z1
    real(RP), intent(in)  :: RLmo
    integer,  intent(in)  :: I_B_TYPE

    real(RP) :: tmp

    zeta = min( max( z1 * RLmo, -ATMOS_PHY_BL_MYNN_zeta_lim ), ATMOS_PHY_BL_MYNN_zeta_lim )

    select case ( I_B_TYPE )
    case ( I_B71 )
       ! Businger et al. (1971)
       if ( zeta >= 0 ) then
          phi_m = 4.7_RP * zeta + 1.0_RP
          phi_h = 4.7_RP * zeta + 0.74_RP
       else
          phi_m = 1.0_RP / sqrt(sqrt( 1.0_RP - 15.0_RP * zeta ))
          phi_h = 0.47_RP / sqrt( 1.0_RP - 9.0_RP * zeta )
       end if
    case ( I_B91, I_B91W01 )
       ! Beljaars and Holtslag (1991)
       if ( zeta >= 0 ) then
          tmp = - 2.0_RP / 3.0_RP * ( 0.35_RP * zeta - 6.0_RP ) * exp(-0.35_RP*zeta) * zeta
          phi_m = tmp + zeta + 1.0_RP
          phi_h = tmp + zeta * sqrt( 1.0_RP + 2.0_RP * zeta / 3.0_RP ) + 1.0_RP
       else
          if ( I_B_TYPE == I_B91W01 ) then
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

    return
  end subroutine get_phi

!OCL SERIAL
  subroutine calc_zi_o2019( &
       KA, KS, KE, &
       Zi, &
       PTLV, &
       tke, &
       z, frac_land, &
       initialize )
    !$acc routine seq
    integer, intent(in) :: KA, KS, KE
    real(RP), intent(out) :: Zi
    real(RP), intent(in)  :: PTLV(KA)
    real(RP), intent(in)  :: tke(KA)
    real(RP), intent(in)  :: z(KA)
    real(RP), intent(in)  :: frac_land
    logical,  intent(in)  :: initialize

    real(RP) :: zit, zie
    real(RP) :: tmin, dtv
    real(RP) :: tke_min
    real(RP) :: we
    integer :: k, k2

    tmin = 999e10_RP
    k2 = KE
    do k = KS, KE
       tmin = min( tmin, PTLV(k) )
       if ( z(k) > 200.0_RP ) then
          k2 = k
          exit
       end if
    end do
    if ( frac_land >= 0.5_RP ) then ! over land
       dtv = 1.25_RP
    else ! over water
       dtv = 0.75_RP
    end if
    zit = z(KE)
    do k = k2+1, KE
       if ( tmin + dtv <= PTLV(k) ) then
          zit = z(k)
          exit
       end if
    end do
    zie = 1000.0_RP
    if ( .not. initialize ) then
       tke_min = max( tke(KS), 0.02_RP )
       do k = KS+1, KE
          if ( tke(k) <= tke_min * 0.05_RP ) then
             zie = z(k)
             exit
          end if
       end do
    end if
    we = 0.5_RP * tanh( ( zit - 200.0_RP ) / 400.0_RP ) + 0.5_RP
    Zi = zit * ( 1.0_RP - we ) + zie * we

    return
  end subroutine calc_zi_o2019

!OCL SERIAL
  subroutine calc_mflux( &
       KA, KS, KE, &
       DENS, &
       POTV, POTL, Qw, &
       U, V, W, &
       tke, &
       cldfrac, &
       EXNER, &
       SHFLX, SFLX_BUOY, SFLX_PT, SFLX_QV, &
       SFC_DENS, &
       Zi, &
       Z, CDZ, F2H, &
       mflux, tflux, qflux, uflux, vflux, eflux, &
       tlh, tvh, qh, uh, vh, wh, eh, dh )

    !$acc routine vector
    use scale_const, only: &
       GRAV    => CONST_GRAV, &
       CPdry => CONST_CPdry, &
       CVdry => CONST_CVdry, &
       EPSTvap => CONST_EPSTvap
    use scale_atmos_hydrometeor, only: &
       CP_VAPOR, &
       CV_VAPOR, &
       LHV
    use scale_atmos_saturation, only: &
       ATMOS_SATURATION_moist_conversion_dens_liq
    integer, intent(in) :: KA, KS, KE

    real(RP), intent(in) :: DENS(KA)
    real(RP), intent(in) :: POTV(KA)
    real(RP), intent(in) :: POTL(KA)
    real(RP), intent(in) :: Qw(KA)
    real(RP), intent(in) :: U(KA)
    real(RP), intent(in) :: V(KA)
    real(RP), intent(in) :: W(KA)
    real(RP), intent(in) :: tke(KA)
    real(RP), intent(in) :: cldfrac(KA)
    real(RP), intent(in) :: EXNER(KA)
    real(RP), intent(in) :: SHFLX
    real(RP), intent(in) :: SFLX_BUOY
    real(RP), intent(in) :: SFLX_PT
    real(RP), intent(in) :: SFLX_QV
    real(RP), intent(in) :: SFC_DENS
    real(RP), intent(in) :: Zi
    real(RP), intent(in) :: Z(KA)
    real(RP), intent(in) :: CDZ(KA)
    real(RP), intent(in) :: F2H(KA,2)

    real(RP), intent(out) :: mflux(0:KA)
    real(RP), intent(out) :: tflux(0:KA)
    real(RP), intent(out) :: qflux(0:KA)
    real(RP), intent(out) :: uflux(0:KA)
    real(RP), intent(out) :: vflux(0:KA)
    real(RP), intent(out) :: eflux(0:KA)

    ! work
    real(RP), intent(out) :: tlh(KA)
    real(RP), intent(out) :: tvh(KA)
    real(RP), intent(out) :: qh(KA)
    real(RP), intent(out) :: uh(KA)
    real(RP), intent(out) :: vh(KA)
    real(RP), intent(out) :: wh(KA)
    real(RP), intent(out) :: eh(KA)
    real(RP), intent(out) :: dh(KA)

    real(RP), parameter :: amax = 0.1_RP
    real(RP), parameter :: Ce = 0.35_RP
    real(RP), parameter :: x = -1.9_RP
    real(RP), parameter :: Cwt = 0.58_RP
    real(RP), parameter :: Cs = 1.34_RP
    real(RP), parameter :: zs = 50.0_RP
    real(RP), parameter :: a = 2.0_RP

    real(RP) :: au ! total area
    real(RP) :: ap(nplume) ! area
    real(RP) :: wp, tp, qp, up, vp, ep
    real(RP) :: zmax
    real(RP) :: er ! entrainment rate
    real(RP) :: sw, sq, st ! sigma_w,q,t
    real(RP) :: ws, ts, qs ! scales
    real(RP) :: buoy ! buoyancy
    real(RP) :: mf
    real(RP) :: b
    real(RP) :: dd
    real(RP) :: temp, qc
    real(RP) :: tps, qps
    real(RP) :: Emoist
    real(RP) :: cp, cv
    real(RP) :: tmp, fact
    logical  :: converged
    integer :: np
    integer :: k, n

    do k = KS-1, KE
       mflux(k) = 0.0_RP
       tflux(k) = 0.0_RP
       qflux(k) = 0.0_RP
       uflux(k) = 0.0_RP
       vflux(k) = 0.0_RP
       eflux(k) = 0.0_RP
    end do

    ! must be a positive surface buoyancy flux
    if ( SFLX_BUOY <= 0.0_RP ) then
       return
    end if

    ! must be superadiabatic in the lowest 50 m
    !$acc loop seq
    do k = KS, KE-1
       if ( z(k) > 50.0_RP ) exit
       if ( POTV(k+1) >= POTV(k) ) then
          return
       end if
    end do

    zmax = zi
    !$acc loop seq
    do k = KS, KE
       if ( cldfrac(k) > 0.5_RP .and. z(k) <= 2000.0_RP ) then
          zmax = min( zmax, z(k) * 0.5_RP ) ! zc
          exit
       end if
    end do
    np = nplume
    !$acc loop seq
    do n = 1, nplume
       if ( dplume(n) > zmax ) then
          np = n - 1
          exit
       end if
    end do
    if ( np == 0 ) return

    ! updraft area
    au = amax * ( tanh( ( SHFLX - 20.0_RP ) / 90.0_RP ) * 0.5_RP + 0.5_RP )
    fact = 0.0_RP
    !$acc loop reduction(+:fact)
    do n = 1, np
       if ( n == 1 ) then
          dd = dplume(2) - dplume(1)
       else if ( n == np ) then
          dd = dplume(np) - dplume(np-1)
       else
          dd = ( dplume(n+1) - dplume(n-1) ) * 0.5_RP
       end if
       ap(n) = dplume(n)**(x + 2.0_RP) * dd
       fact = fact + ap(n)
    end do
    fact = au / fact
    do n = 1, np
       ap(n) = ap(n) * fact
    end do

    ! @half levels
    do k = KS, KE-1
       tlh(k) = F2H(k,1) * POTL(k) + F2H(k,2) * POTL(k+1)
       tvh(k) = F2H(k,1) * POTV(k) + F2H(k,2) * POTV(k+1)
       qh(k) = F2H(k,1) * Qw  (k) + F2H(k,2) * Qw  (k+1)
       uh(k) = F2H(k,1) * U   (k) + F2H(k,2) * U   (k+1)
       vh(k) = F2H(k,1) * V   (k) + F2H(k,2) * V   (k+1)
       wh(k) = F2H(k,1) * W   (k) + F2H(k,2) * W   (k+1)
       eh(k) = F2H(k,1) * tke (k) + F2H(k,2) * tke (k+1)
       dh(k) = F2H(k,1) * DENS(k) + F2H(k,2) * DENS(k+1)
    end do

    ws = max( ( zi * SFLX_BUOY )**OneOverThree, 1.0E-10_RP )
    ts = SFLX_PT / ( SFC_DENS * ws )
    qs = max( SFLX_QV, 0.0_RP ) / ( SFC_DENS * ws )
    tmp = ( zs / zi )**OneOverThree
    sw = Cs * ws * tmp * ( 1.0_RP - 0.8_RP*zs/zi )
    sq = Cs * qs * tmp
    st = Cs * ts * tmp
    !$acc loop seq
    do n = 1, np
       ! initial properties of the plume
       wp = min( pw(n) * sw, 0.5_RP )
       tp = tlh(KS) + wp * Cwt * st / sw
       qp = qh(KS) + wp * Cwt * sq / sw
       up = uh(KS)
       vp = vh(KS)
       ep = eh(KS)
       !$acc loop seq
       do k = KS, KE-1
          mf = ap(n) * ( wp - wh(k) ) * dh(k)
          mflux(k) = mflux(k) + mf
          tflux(k) = tflux(k) + mf * ( tp - tlh(k) )
          qflux(k) = qflux(k) + mf * ( qp - qh(k) )
          uflux(k) = uflux(k) + mf * ( up - uh(k) )
          vflux(k) = vflux(k) + mf * ( vp - vh(k) )
          eflux(k) = eflux(k) + mf * ( ep - eh(k) )

          er = Ce / ( wp * dplume(n) )

          ! d/dz phi = - er ( phi - phi_env )
          fact = exp( - er * CDZ(k+1) )
          tp = POTL(k+1) + ( tp - POTL(k+1) ) * fact
          qp = Qw(k+1) + ( qp - Qw(k+1) ) * fact
          up = U(k+1) + ( up - U(k+1) ) * fact
          vp = V(k+1) + ( vp - V(k+1) ) * fact
          ep = tke(k+1) + ( ep - tke(k+1) ) * fact

          qps = qp
          qc = 0.0_RP
          temp = tp * EXNER(k+1)
          cp = CPdry * ( 1.0_RP - qp ) + CP_VAPOR * qp
          cv = CVdry * ( 1.0_RP - qp ) + CV_VAPOR * qp
          Emoist = temp * cv + qp * LHV
          call ATMOS_SATURATION_moist_conversion_dens_liq( DENS(k+1), Emoist, & ! (in)
                                                           temp, qps, qc,     & ! (inout)
                                                           cp, cv,            & ! (inout)
                                                           converged          ) ! (out)
          if ( converged ) then
             tps = temp / EXNER(k+1)
          else
             LOG_WARN("calc_mflux",*) "moist_conversion did not converged"
             ! tentative
             qps = qp
             tps = tp
          end if
          buoy = GRAV * ( tps * ( 1.0_RP + EPSTvap * qps ) / POTV(k+1) - 1.0_RP )
          if ( buoy > 0.0_RP ) then
             b = 0.15_RP
          else
             b = 0.2_RP
          end if
          ! d/dz w = - er a w + b buoy / w
          wp = sqrt( max( 0.0_RP, ( er * a * wp**2 - b * buoy ) * exp( - 2.0_RP * er * a * CDZ(k+1) ) + b * buoy ) / ( er * a ) )
          if ( wp <= 0.0_RP ) then
             exit
          end if
       end do
    end do


    return
  end subroutine calc_mflux

end module scale_atmos_phy_bl_mynn
