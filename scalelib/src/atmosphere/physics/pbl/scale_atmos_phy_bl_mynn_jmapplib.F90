!-------------------------------------------------------------------------------
!> module atmosphere / physics / pbl / mynn-jmapplib
!!
!! @par Description
!!          Boundary layer turbulence model
!!          Mellor-Yamada Nakanishi-Niino model implemented in the JMA Physics Process Library
!!
!! @author Team SCALE
!!
!<
!-------------------------------------------------------------------------------
#include "scalelib.h"
module scale_atmos_phy_bl_mynn_jmapplib
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_io
  use scale_prof

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
  public :: ATMOS_PHY_BL_MYNN_JMAPPLIB_setup
  public :: ATMOS_PHY_BL_MYNN_JMAPPLIB_finalize
  public :: ATMOS_PHY_BL_MYNN_JMAPPLIB_tendency

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
#ifdef JMAPPLIB
  integer,                public, parameter :: ATMOS_PHY_BL_MYNN_JMAPPLIB_NTRACER = 4
  character(len=H_SHORT), public :: ATMOS_PHY_BL_MYNN_JMAPPLIB_NAME(4) = &
       (/ 'TKE_MYNN', 'TSQ_MYNN', 'QSQ_MYNN', 'COV_MYNN' /)
  character(len=H_LONG),  public :: ATMOS_PHY_BL_MYNN_JMAPPLIB_DESC(4) = &
       (/ 'doubled turbulent kinetic energy (MYNN)                                                 ', &
          'sub-grid variance of liquid water potential temperature (MYNN)                          ', &
          'sub-grid variance of total water content (MYNN)                                         ', &
          'sub-grid covariance of liquid water potential temperature and total water content (MYNN)' /)
  character(len=H_SHORT), public :: ATMOS_PHY_BL_MYNN_JMAPPLIB_UNITS(4) = &
       (/ 'm2/s2  ', 'K2     ', 'kg2/kg2', 'K kg   ' /)
#else
  integer,                public, parameter :: ATMOS_PHY_BL_MYNN_JMAPPLIB_NTRACER = 0
  character(len=H_SHORT), public :: ATMOS_PHY_BL_MYNN_JMAPPLIB_NAME(1) = (/ '' /)
  character(len=H_LONG),  public :: ATMOS_PHY_BL_MYNN_JMAPPLIB_DESC(1) = (/ '' /)
  character(len=H_SHORT), public :: ATMOS_PHY_BL_MYNN_JMAPPLIB_UNITS(1) = (/ '' /)
#endif

  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
#ifdef JMAPPLIB
  integer,  private, parameter :: I_TKE = 1
  integer,  private, parameter :: I_TSQ = 2
  integer,  private, parameter :: I_QSQ = 3
  integer,  private, parameter :: I_COV = 4

  integer,  private :: KE_PBL

  character(len=3), private :: ATMOS_PHY_BL_MYNN_JMAPPLIB_LEVEL           = "3"
  logical,          private :: ATMOS_PHY_BL_MYNN_JMAPPLIB_SHCU_BUOY_FLAG  = .true.
  real(RP),         private :: ATMOS_PHY_BL_MYNN_JMAPPLIB_PBL_MAX         = 1.E+99_RP !> maximum height of the PBL
  real(RP),         private :: ATMOS_PHY_BL_MYNN_JMAPPLIB_SHCU_MAX        = 1.E+99_RP !> maximum height of the shallow cumulus
  real(RP),         private :: ATMOS_PHY_BL_MYNN_JMAPPLIB_SGM_MIN_FCT     = 0.09_RP

  namelist / PARAM_ATMOS_PHY_BL_MYNN_JMAPPLIB / &
       ATMOS_PHY_BL_MYNN_JMAPPLIB_LEVEL, &
       ATMOS_PHY_BL_MYNN_JMAPPLIB_SHCU_BUOY_FLAG, &
       ATMOS_PHY_BL_MYNN_JMAPPLIB_PBL_MAX, &
       ATMOS_PHY_BL_MYNN_JMAPPLIB_SHCU_MAX, &
       ATMOS_PHY_BL_MYNN_JMAPPLIB_SGM_MIN_FCT
#endif

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> ATMOS_PHY_BL_MYNN_JMAPPLIB_setup
  !! Setup
  !<
  subroutine ATMOS_PHY_BL_MYNN_JMAPPLIB_setup( &
       KA, KS, KE, &
       CZ, &
       dt, PBL_MAX, SHCU_MAX )
    use scale_prc, only: &
       PRC_abort
    use scale_const, only: &
       PRE00 => CONST_PRE00
#ifdef JMAPPLIB
    use pbl_const, only: &
       pbl_const_ini
    use pbl_parm, only: &
       pbl_parm_ini
    use pbl_grid, only: &
       pbl_grid_ini
    use pbl_mym_option, only: &
       pbl_mym_option_ini
    use pbl_mym_option_symbol, only: &
       mymodel3, &
       mymodel25
    use pbl_mym_parm, only: &
       pbl_mym_parm_ini
    use pbl_mym_const, only: &
       pbl_mym_const_ini
#endif
    implicit none

    integer,  intent(in) :: KA, KS, KE

    real(RP), intent(in) :: CZ(KA)

    real(DP), intent(in), optional :: dt
    real(RP), intent(in), optional :: PBL_MAX
    real(RP), intent(in), optional :: SHCU_MAX

    real(RP) :: e_emit

    integer :: KE_shcu
    integer :: nz_pbl, nz_shcu

    integer :: ierr
    integer :: k
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("ATMOS_PHY_BL_MYNN_JMAPPLIB_setup",*) 'Setup'
    LOG_INFO("ATMOS_PHY_BL_MYNN_JMAPPLIB_setup",*) 'Mellor-Yamada Nakanishi-Niino scheme implemented in the JMA Physics Process Library'

#ifdef JMAPPLIB
    if ( present(PBL_MAX) ) ATMOS_PHY_BL_MYNN_JMAPPLIB_PBL_MAX = PBL_MAX
    if ( present(SHCU_MAX) ) ATMOS_PHY_BL_MYNN_JMAPPLIB_SHCU_MAX = SHCU_MAX

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_ATMOS_PHY_BL_MYNN_JMAPPLIB,iostat=ierr)
    if( ierr > 0 ) then !--- fatal error
       LOG_ERROR("ATMOS_PHY_BL_MYNN_JMAPPLIB_setup",*) 'Not appropriate names in namelist PARAM_ATMOS_PHY_BL_MYNN_JMAPPLIB. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_ATMOS_PHY_BL_MYNN_JMAPPLIB)


    KE_PBL  = KS+1
    KE_SHCU = KS+1
    do k = KS+2, KE-1
       if ( ATMOS_PHY_BL_MYNN_JMAPPLIB_PBL_MAX >= CZ(k) ) then
          KE_PBL = k
       end if
       if ( ATMOS_PHY_BL_MYNN_JMAPPLIB_SHCU_MAX >= CZ(k) ) then
          KE_SHCU = k
       end if
    end do
    nz_pbl  = KE_PBL - KS + 1
    nz_shcu = min(KE_SHCU - KS + 1, nz_pbl)

    e_emit = 1.0_RP
    call pbl_const_ini(pref_in = PRE00, timestep_in = real(dt,RP), e_emit_in = e_emit)
    call pbl_parm_ini
    call pbl_grid_ini(nz_pbl, shcu_levels_in = nz_shcu)

    select case ( ATMOS_PHY_BL_MYNN_JMAPPLIB_LEVEL )
    case ( "3" )
       call pbl_mym_option_ini(levflag_in = mymodel3, l_shcu_buoy_in = ATMOS_PHY_BL_MYNN_JMAPPLIB_SHCU_BUOY_FLAG)

    case ( "2.5" )
       call pbl_mym_option_ini(levflag_in = mymodel25, l_shcu_buoy_in = ATMOS_PHY_BL_MYNN_JMAPPLIB_SHCU_BUOY_FLAG)
    case default
       LOG_ERROR("ATMOS_PHY_BL_MYNN_JMAPPLIB_setup",*) 'only level 2.5 and 3 are supported at this moment'
       call PRC_abort
    end select

    call pbl_mym_parm_ini(my_sgm_min_fct_in = ATMOS_PHY_BL_MYNN_JMAPPLIB_SGM_MIN_FCT)
    call pbl_mym_const_ini

#else

    LOG_ERROR("ATMOS_PHY_BL_MYNN_JMAPPLIB_setup",*) 'To use "MYNN-JMAPPLIB", compile SCALE with "SCALE_ENABLE_JMAPPLIB=T" option.'
    call PRC_abort

#endif

    return
  end subroutine ATMOS_PHY_BL_MYNN_JMAPPLIB_setup

  !-----------------------------------------------------------------------------
  !! Finalize
  !<
  subroutine ATMOS_PHY_BL_MYNN_JMAPPLIB_finalize
    implicit none

    return
  end subroutine ATMOS_PHY_BL_MYNN_JMAPPLIB_finalize

  !-----------------------------------------------------------------------------
  !> ATMOS_PHY_BL_MYNN_JMAPPLIB_tendency
  !! calculate tendency by the virtical eddy viscosity
  !<
  subroutine ATMOS_PHY_BL_MYNN_JMAPPLIB_tendency( &
       KA, KS, KE, IA, IS, IE, JA, JS, JE, &
       DENS, U, V, POTT, PROG,             &
       PRES, EXNER,                        &
       QDRY, QV, QC, QI,                   &
       SFC_DENS, SFC_PRES,                 &
       SFLX_MU, SFLX_MV, SFLX_SH, SFLX_QV, &
       us, RLmo,                           &
       CZ, FZ, F2H, dt,                    &
       RHOU_t, RHOV_t, RHOT_t, RHOQV_t,    &
       RPROG_t,                            &
       Nu, Kh,                             &
       Zi, SFLX_BUOY                       )
    use scale_const, only: &
       CPdry => CONST_CPdry
    use scale_atmos_hydrometeor, only: &
       CP_VAPOR
    use scale_file_history, only: &
       FILE_HISTORY_in
#ifdef JMAPPLIB
    use pbl_mym_main, only: &
       pbl_mym_main_level3, &
       pbl_mym_main_level25
    use pbl_coupler, only: &
       pbl_coupler_flx_force_tend_run
    use pbl_diag, only: &
       pbl_diag_pbl_height
#endif
    implicit none

    integer, intent(in) :: KA, KS, KE
    integer, intent(in) :: IA, IS, IE
    integer, intent(in) :: JA, JS, JE

    real(RP), intent(in) :: DENS    (KA,IA,JA) !> density
    real(RP), intent(in) :: U       (KA,IA,JA) !> zonal wind
    real(RP), intent(in) :: V       (KA,IA,JA) !> meridional wind
    real(RP), intent(in) :: POTT    (KA,IA,JA) !> potential temperature
    real(RP), intent(in) :: PROG    (KA,IA,JA,ATMOS_PHY_BL_MYNN_JMAPPLIB_ntracer) !> prognostic variables (TKE, TSQ, QSQ, COV)
    real(RP), intent(in) :: PRES    (KA,IA,JA) !> pressure
    real(RP), intent(in) :: EXNER   (KA,IA,JA) !> Exner function
    real(RP), intent(in) :: QDRY    (KA,IA,JA) !> dry air
    real(RP), intent(in) :: QV      (KA,IA,JA) !> vapor
    real(RP), intent(in) :: QC      (KA,IA,JA) !> cloud water
    real(RP), intent(in) :: QI      (KA,IA,JA) !> cloud ice water
    real(RP), intent(in) :: SFC_DENS(   IA,JA) !> surface density
    real(RP), intent(in) :: SFC_PRES(   IA,JA) !> surface pressure
    real(RP), intent(in) :: SFLX_MU (   IA,JA) !> surface flux of zonal wind
    real(RP), intent(in) :: SFLX_MV (   IA,JA) !> surface flux of meridional wind
    real(RP), intent(in) :: SFLX_SH (   IA,JA) !> surface sensible heat flux
    real(RP), intent(in) :: SFLX_QV (   IA,JA) !> surface sensible QV flux
    real(RP), intent(in) :: us      (   IA,JA) !> friction velocity
    real(RP), intent(in) :: RLmo    (   IA,JA) !> inverse of Monin-Obukhov length

    real(RP), intent(in) :: CZ(  KA,IA,JA)
    real(RP), intent(in) :: FZ(0:KA,IA,JA)
    real(RP), intent(in) :: F2H(KA,2,IA,JA)
    real(DP), intent(in) :: dt

    real(RP), intent(out) :: RHOU_t (KA,IA,JA) !> tendency of dens * u
    real(RP), intent(out) :: RHOV_t (KA,IA,JA) !> tendency of dens * v
    real(RP), intent(out) :: RHOT_t (KA,IA,JA) !> tendency of dens * pt
    real(RP), intent(out) :: RHOQV_t(KA,IA,JA) !> tendency of dens * qv
    real(RP), intent(out) :: RPROG_t(KA,IA,JA,ATMOS_PHY_BL_MYNN_JMAPPLIB_ntracer) !> tendency of dens * prognostic variables (TKE, TSQ, QSQ, COV)
    real(RP), intent(out) :: Nu     (KA,IA,JA) !> eddy viscosity coefficient @ half-level
    real(RP), intent(out) :: Kh     (KA,IA,JA) !> eddy diffusion coefficient @ half-level
    real(RP), intent(out) :: Zi        (IA,JA) !> PBL height
    real(RP), intent(out) :: SFLX_BUOY (IA,JA)

#ifdef JMAPPLIB
    real(RP) :: l(KA,IA,JA)

    real(RP) :: qke(KS:KE_PBL)
    real(RP) :: qv_lc(KS:KE_PBL)
    real(RP) :: qc_lc(KS:KE_PBL)
    real(RP) :: qi_lc(KS:KE_PBL)
    real(RP) :: dens_lc(KS:KE_PBL)
    real(RP) :: dfm(KS:KE_PBL)
    real(RP) :: dfh(KS:KE_PBL)
    real(RP) :: taux_ex(KS:KE_PBL)
    real(RP) :: tauy_ex(KS:KE_PBL)
    real(RP) :: ftl_ex(KS:KE_PBL)
    real(RP) :: fqw_ex(KS:KE_PBL)
    real(RP) :: tend_qke(KS:KE_PBL)
    real(RP) :: tend_tsq(KS:KE_PBL)
    real(RP) :: tend_qsq(KS:KE_PBL)
    real(RP) :: tend_cov(KS:KE_PBL)
    real(RP) :: tend_u(KS:KE_PBL)
    real(RP) :: tend_v(KS:KE_PBL)
    real(RP) :: tend_pt(KS:KE_PBL)
    real(RP) :: tend_qv(KS:KE_PBL)
    real(RP) :: z_f(KS:KE_PBL)
    real(RP) :: dz_f(KS:KE_PBL)
    real(RP) :: rdz_f(KS:KE_PBL)
    real(RP) :: rdz_h(KS:KE_PBL)
    real(RP) :: h2f_m(KS:KE_PBL)
    real(RP) :: h2f_p(KS:KE_PBL)
    real(RP) :: fb_surf
    real(RP) :: rho_ov_rhoa
    real(RP) :: SFLX_U
    real(RP) :: SFLX_V
    real(RP) :: SFLX_PT
    real(RP) :: SFLX_Q
    real(RP) :: CPtot

    real(RP) :: diss(KA,IA,JA) !> TKE dissipation term
    real(RP) :: dummy(KA,IA,JA)

    integer :: k, i, j
    !---------------------------------------------------------------------------

    LOG_PROGRESS(*) "atmosphere / physics / pbl / MYNN-JMAPPLIB"

    h2f_m(:) = 0.5_RP
    h2f_p(:) = 0.5_RP

    !$omp parallel do &
    !$omp private(qke,qv_lc,qc_lc,qi_lc,dens_lc,taux_ex,tauy_ex,ftl_ex,fqw_ex,&
    !$omp         tend_qke,tend_tsq,tend_qsq,tend_cov,tend_u,tend_v,tend_pt,tend_qv, &
    !$omp         z_f,dz_f,rdz_f,rdz_h, &
    !$omp         rho_ov_rhoa,SFLX_U,SFLX_V,SFLX_PT,SFLX_Q,CPtot)
    do j = JS, JE
    do i = IS, IE

       do k = KS, KE_PBL
          qke(k) = PROG(k,i,j,I_TKE) * 2.0_RP
          rho_ov_rhoa = 1.0_RP / ( QDRY(k,i,j) + QV(k,i,j) )
          qv_lc(k) = QV(k,i,j) * rho_ov_rhoa
          qc_lc(k) = QC(k,i,j) * rho_ov_rhoa
          qi_lc(k) = QI(k,i,j) * rho_ov_rhoa
          z_f(k) = CZ(k,i,j) - FZ(KS-1,i,j)
          dz_f(k) = FZ(k,i,j) - FZ(k-1,i,j)
          rdz_f(k) = 1.0_RP / dz_f(k)
          rdz_h(k) = 1.0_RP / ( CZ(k+1,i,j) - CZ(k,i,j) )
       end do

       CPtot = CPdry + SFLX_QV(i,j) * ( CP_VAPOR - CPdry )

       SFLX_U  = SFLX_MU(i,j) / SFC_DENS(i,j)
       SFLX_V  = SFLX_MV(i,j) / SFC_DENS(i,j)
       SFLX_PT = SFLX_SH(i,j) / ( CPtot * EXNER(KS,i,j) * SFC_DENS(i,j) )
       SFLX_Q  = SFLX_QV(i,j) / SFC_DENS(i,j)

       select case ( ATMOS_PHY_BL_MYNN_JMAPPLIB_LEVEL )
       case ( "3" )
          call pbl_mym_main_level3( &
               SFLX_U, SFLX_V, SFLX_PT, SFLX_Q, & ! (in)
               us(i,j), RLmo(i,j), SFC_PRES(i,j), & ! (in)
               U(KS:KE_PBL,i,j), V(KS:KE_PBL,i,j), POTT(KS:KE_PBL,i,j), & ! (in)
               qv_lc, qc_lc, qi_lc, PRES(KS:KE_PBL,i,j), EXNER(KS:KE_PBL,i,j), & ! (in)
               qke(:), PROG(KS:KE_PBL,i,j,I_TSQ), PROG(KS:KE_PBL,i,j,I_QSQ), PROG(KS:KE_PBL,i,j,I_COV), & ! (in)
               z_f, dz_f, rdz_f, rdz_h, F2H(KS:KE_PBL,2,i,j), F2H(KS:KE_PBL,1,i,j), h2f_m, h2f_p, & ! (in)
               SFLX_BUOY(i,j), & ! (out)
               dfm(:), dfh(:), tend_qke(:), tend_tsq(:), tend_qsq(:), tend_cov(:), & ! (out)
               taux_ex, tauy_ex, ftl_ex, fqw_ex, & ! (out)
               l_e = l(KS:KE_PBL,i,j), tke_dis = diss(KS:KE_PBL,i,j), tsq_dis=dummy(KS:KE_PBL,i,j) ) ! (out)
       case ( "2.5" )
          call pbl_mym_main_level25( &
               SFLX_U, SFLX_V, SFLX_PT, SFLX_Q, & ! (in)
               us(i,j), RLmo(i,j), SFC_PRES(i,j), & ! (in)
               U(KS:KE_PBL,i,j), V(KS:KE_PBL,i,j), POTT(KS:KE_PBL,i,j), & ! (in)
               qv_lc, qc_lc, qi_lc, PRES(KS:KE_PBL,i,j), EXNER(KS:KE_PBL,i,j), & ! (in)
               qke(:), PROG(KS:KE_PBL,i,j,I_TSQ), PROG(KS:KE_PBL,i,j,I_QSQ), PROG(KS:KE_PBL,i,j,I_COV), & ! (in)
               z_f, dz_f, rdz_f, rdz_h, F2H(KS:KE_PBL,2,i,j), F2H(KS:KE_PBL,1,i,j), h2f_m, h2f_p, & ! (in)
               SFLX_BUOY(i,j), & ! (out)
               dfm(:), dfh(:), tend_qke(:), tend_tsq(:), tend_qsq(:), tend_cov(:), & ! (out)
               taux_ex, tauy_ex, ftl_ex, fqw_ex, & ! (out)
               l_e = l(KS:KE_PBL,i,j), tke_dis = diss(KS:KE_PBL,i,j), tsq_dis=dummy(KS:KE_PBL,i,j) ) ! (out)
          do k = KS, KE_PBL
             tend_tsq(k) = ( tend_tsq(k) - PROG(k,i,j,I_TSQ) ) / dt
             tend_qsq(k) = ( tend_qsq(k) - PROG(k,i,j,I_QSQ) ) / dt
             tend_cov(k) = ( tend_cov(k) - PROG(k,i,j,I_COV) ) / dt
          end do
       end select

       do k = KS, KE_PBL
          dens_lc(k) = DENS(k,i,j) * ( QDRY(k,i,j) + QV(k,i,j) )
       end do

       call pbl_coupler_flx_force_tend_run( &
            SFC_DENS(i,j), SFLX_U, SFLX_V, SFLX_PT, SFLX_Q, & ! (in)
            dfm(:), dfh(:), dens_lc(:),       & ! (in)
            taux_ex, tauy_ex, ftl_ex, fqw_ex, & ! (in)
            F2H(KS:KE_PBL,2,i,j), F2H(KS:KE_PBL,1,i,j), rdz_f, rdz_h, & ! (in)
            tend_u, tend_v, tend_pt, tend_qv ) ! (out)

       do k = KS, KE_PBL
          RHOU_t(k,i,j) = tend_u(k) * DENS(k,i,j)
          RHOV_t(k,i,j) = tend_v(k) * DENS(k,i,j)
          RHOT_t(k,i,j) = tend_pt(k) * DENS(k,i,j)
          RHOQV_t(k,i,j) = tend_qv(k) * DENS(k,i,j) * ( QDRY(k,i,j) + QV(k,i,j) )
          RPROG_t(k,i,j,I_TKE) = tend_qke(k) * DENS(k,i,j) * 0.5_RP
          RPROG_t(k,i,j,I_TSQ) = tend_tsq(k) * DENS(k,i,j)
          RPROG_t(k,i,j,I_QSQ) = tend_qsq(k) * DENS(k,i,j)
          RPROG_t(k,i,j,I_COV) = tend_cov(k) * DENS(k,i,j)
       end do
       RHOU_t (KS,i,j) = RHOU_t (KS,i,j) - SFLX_MU(i,j) * rdz_f(KS)
       RHOV_t (KS,i,j) = RHOV_t (KS,i,j) - SFLX_MV(i,j) * rdz_f(KS)
       RHOT_t (KS,i,j) = RHOT_t (KS,i,j) - SFLX_PT * SFC_DENS(i,j) * rdz_f(KS)
       RHOQV_t(KS,i,j) = RHOQV_t(KS,i,j) - SFLX_QV(i,j) * rdz_f(KS)
       do k = KE_PBL+1, KE
          RHOU_t(k,i,j) = 0.0_RP
          RHOV_t(k,i,j) = 0.0_RP
          RHOT_t(k,i,j) = 0.0_RP
          RHOQV_t(k,i,j) = 0.0_RP
          RPROG_t(k,i,j,I_TKE) = 0.0_RP
          RPROG_t(k,i,j,I_TSQ) = 0.0_RP
          RPROG_t(k,i,j,I_QSQ) = 0.0_RP
          RPROG_t(k,i,j,I_COV) = 0.0_RP
       end do

       Nu(KS-1,i,j) = 0.0_RP
       Kh(KS-1,i,j) = 0.0_RP
       do k = KS, KE_PBL-1
          Nu(k,i,j) = dfm(k+1) * F2H(k,1,i,j) + dfm(k) * F2H(k,2,i,j)
          Kh(k,i,j) = dfh(k+1) * F2H(k,1,i,j) + dfh(k) * F2H(k,2,i,j)
       end do
       do k = KE_PBL, KE
          Nu(k,i,j) = 0.0_RP
          Kh(k,i,j) = 0.0_RP
       end do

       call pbl_diag_pbl_height( &
            SFLX_PT, us(i,j), RLmo(i,j), & ! (in)
            U(KS:KE_PBL,i,j), V(KS:KE_PBL,i,j), POTT(KS:KE_PBL,i,j), qv_lc, z_f, & ! (in)
            Zi(i,j) ) ! (out)

    end do
    end do


    l(KE_PBL+1:KE,:,:) = UNDEF
    call FILE_HISTORY_in(l(:,:,:), 'L_mix_MYNN', 'minxing length', 'm', fill_halo=.true.)

    diss(KS:KE_PBL,:,:) = - diss(KS:KE_PBL,:,:)
    diss(KE_PBL+1:KE,:,:) = UNDEF
    call FILE_HISTORY_in(diss(:,:,:), 'TKE_diss_MYNN', 'TKE dissipation', 'm2/s3', fill_halo=.true.)
#endif

    return
  end subroutine ATMOS_PHY_BL_MYNN_JMAPPLIB_tendency

end module scale_atmos_phy_bl_mynn_jmapplib
