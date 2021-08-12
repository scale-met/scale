!-------------------------------------------------------------------------------
!> module atmosphere / physics / cumulus / kf-jmapplib
!!
!! @par Description
!!          Cumulus parameterization
!!          Kain-Fritsch scheme implemented in JMA Physics Process library
!!
!! @author Team SCALE
!!
!<
!-------------------------------------------------------------------------------
#include "scalelib.h"
module scale_atmos_phy_cp_kf_jmapplib
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
  public :: ATMOS_PHY_CP_KF_JMAPPLIB_setup
  public :: ATMOS_PHY_CP_KF_JMAPPLIB_finalize
  public :: ATMOS_PHY_CP_KF_JMAPPLIB_tendency

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
#ifdef JMAPPLIB
  character(len=6), private :: ATMOS_PHY_CP_KF_JMAPPLIB_type = "KF1701"  !> KF1701 for MSM; KF for LFM
  real(RP), allocatable, private  :: ishall_counter(:,:)
  real(RP), allocatable, private  :: ideep_counter (:,:)
#endif
  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> ATMOS_PHY_CP_KF_JMAPPLIB_setup
  !! Setup
  !<
  subroutine ATMOS_PHY_CP_KF_JMAPPLIB_setup( &
       KA, KS, KE, IA, JA, &
       dx, dt )
    use scale_prc, only: &
       PRC_abort
    use scale_const, only: &
       PRE00 => CONST_PRE00
#ifdef JMAPPLIB
    use conv_grid, only: &
       conv_grid_ini
    use conv_kf_const, only: &
       conv_kf_const_ini
    use conv_kf_parm, only: &
       conv_kf_parm_ini, &
       conv_kf_parm_get_stepcu
    use conv_kf_parm1701, only: &
       conv_kf_parm_ini1701, &
       conv_kf_parm_get_stepcu1701
    use conv_kf_lut, only: &
       conv_kf_lut_ini
#endif
    implicit none

    integer,  intent(in) :: KA, KS, KE
    integer,  intent(in) :: IA, JA

    real(RP), intent(in) :: dx
    real(DP), intent(in) :: dt

    integer  :: ATMOS_PHY_CP_KF_JMAPPLIB_radius_type   = 3         !> 3 for KF1701
                                                                   !> 2 for KF Narita and Moriyasu (2010)
                                                                   !> 1 for Narita (2008)
                                                                   !> 0 for Kain (2004)
    real(RP) :: ATMOS_PHY_CP_KF_JMAPPLIB_dtlcl_fact    = 0.5_RP    !> ignored for KF1701
    real(RP) :: ATMOS_PHY_CP_KF_JMAPPLIB_dpmin         = 5000.0_RP !> 5000 for KF1701; 1000 for KF
    real(RP) :: ATMOS_PHY_CP_KF_JMAPPLIB_dlifetime     = 600.0_RP  !> 600 for KF1701; 3600 for KF
    real(RP) :: ATMOS_PHY_CP_KF_JMAPPLIB_slifetime     = 600.0_RP  !> 600 for KF1701; 3600 for KF
    real(RP) :: ATMOS_PHY_CP_KF_JMAPPLIB_detlq_qr_rate = 0.0_RP    !> 0.0 for KF1701; 0.0 for KF
    real(RP) :: ATMOS_PHY_CP_KF_JMAPPLIB_detic_qs_rate = 1.0_RP    !> 1.0 for KF1701; 0.0 for KF

    namelist / PARAM_ATMOS_PHY_CP_KF_JMAPPLIB / &
         ATMOS_PHY_CP_KF_JMAPPLIB_type, &
         ATMOS_PHY_CP_KF_JMAPPLIB_radius_type, &
         ATMOS_PHY_CP_KF_JMAPPLIB_dtlcl_fact, &
         ATMOS_PHY_CP_KF_JMAPPLIB_dpmin, &
         ATMOS_PHY_CP_KF_JMAPPLIB_dlifetime, &
         ATMOS_PHY_CP_KF_JMAPPLIB_slifetime, &
         ATMOS_PHY_CP_KF_JMAPPLIB_detlq_qr_rate, &
         ATMOS_PHY_CP_KF_JMAPPLIB_detic_qs_rate

    integer :: KMAX

    integer :: ierr
    integer :: k
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("ATMOS_PHY_CP_KF_JMAPPLIB_setup",*) 'Setup'
    LOG_INFO("ATMOS_PHY_CP_KF_JMAPPLIB_setup",*) 'KF scheme implemented in the JMA Physics Process Library'

#ifdef JMAPPLIB

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_ATMOS_PHY_CP_KF_JMAPPLIB,iostat=ierr)
    if( ierr > 0 ) then !--- fatal error
       LOG_ERROR("ATMOS_PHY_CP_KF_JMAPPLIB_setup",*) 'Not appropriate names in namelist PARAM_ATMOS_PHY_CP_KF_JMAPPLIB. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_ATMOS_PHY_CP_KF_JMAPPLIB)


    KMAX = KE - KS + 1
    call conv_grid_ini( KMAX )
    call conv_kf_const_ini( dx, PRE00, dt )

    select case ( ATMOS_PHY_CP_KF_JMAPPLIB_type )
    case ( "KF" )

       call conv_kf_parm_ini( &
            dt,                                                           &
            kfrad_var_in        = ATMOS_PHY_CP_KF_JMAPPLIB_radius_type,   &
            cudt_in             = -1.0_RP,                                &
            dtlcl_fct_in        = ATMOS_PHY_CP_KF_JMAPPLIB_dtlcl_fact,    &
            dpmin_in            = ATMOS_PHY_CP_KF_JMAPPLIB_dpmin,         &
            cu_lifetime_min_in  = ATMOS_PHY_CP_KF_JMAPPLIB_dlifetime,     &
            shallow_lifetime_in = ATMOS_PHY_CP_KF_JMAPPLIB_slifetime,     &
            detlq_qr_r_in       = ATMOS_PHY_CP_KF_JMAPPLIB_detlq_qr_rate, &
            detic_qs_r_in       = ATMOS_PHY_CP_KF_JMAPPLIB_detic_qs_rate  )

!       call conv_kf_parm_get_stepcu(stepcu) ! (out)

    case ( "KF1701" )

       call conv_kf_parm_ini1701( &
            dt,                                                           &
            kfrad_var_in        = ATMOS_PHY_CP_KF_JMAPPLIB_radius_type,   &
            kfrad_in            = 1000.0_RP,                              &
            cudt_in             = -1.0_RP,                                &
            dpmin_in            = ATMOS_PHY_CP_KF_JMAPPLIB_dpmin,         &
            cu_lifetime_min_in  = ATMOS_PHY_CP_KF_JMAPPLIB_dlifetime,     &
            shallow_lifetime_in = ATMOS_PHY_CP_KF_JMAPPLIB_slifetime,     &
            detlq_qr_r_in       = ATMOS_PHY_CP_KF_JMAPPLIB_detlq_qr_rate, &
            detic_qs_r_in       = ATMOS_PHY_CP_KF_JMAPPLIB_detic_qs_rate, &
            tkemax_in           = 1.0_RP                                  )
            
!       call conv_kf_parm_get_stepcu1701(stepcu) ! (out)

    case default
       LOG_ERROR("ATMOS_PHY_CP_KF_JMAPPLIB_setup",*) 'ATMOS_PHY_CP_KF_TYPE must be either "KF" or "KF1701"'
       call PRC_abort
    end select

    call conv_kf_lut_ini

#else

    LOG_ERROR("ATMOS_PHY_CP_KF_JMAPPLIB_setup",*) 'To use "KF-JMAPPLIB", compile SCALE with "SCALE_ENABLE_JMAPPLIB=T" option.'
    call PRC_abort

#endif

    allocate( ishall_counter(IA,JA), ideep_counter(IA,JA) )
    ishall_counter(:,:) = 0.0_RP
    ideep_counter (:,:) = 0.0_RP


    return
  end subroutine ATMOS_PHY_CP_KF_JMAPPLIB_setup

  !-----------------------------------------------------------------------------
  !! Finalize
  !<
  subroutine ATMOS_PHY_CP_KF_JMAPPLIB_finalize
    implicit none

    deallocate( ishall_counter, ideep_counter )

    return
  end subroutine ATMOS_PHY_CP_KF_JMAPPLIB_finalize

  !-----------------------------------------------------------------------------
  !> ATMOS_PHY_CP_KF_JMAPPLIB_tendency
  !! calculate tendency by the virtical eddy viscosity
  !<
  subroutine ATMOS_PHY_CP_KF_JMAPPLIB_tendency( &
       KA, KS, KE, IA, IS, IE, JA, JS, JE, &
       DENS, U, V, W, TEMP, POTT,       &
       PRES, EXNER,                     &
       QDRY, QV, QC, QI,                &
       us, PBLH, SFLX_BUOY,             &
       CZ, FZ,                          &
       DENS_t, RHOT_t,                  &
       RHOQV_t, RHOQ_t,                 &
       w0avg, nca,                      &
       SFLX_rain, SFLX_snow, SFLX_prec, &
       cloudtop, cloudbase              )
    use scale_atmos_hydrometeor, only: &
       N_HYD, &
       I_HC, &
       I_HR, &
       I_HI, &
       I_HS, &
       I_HG
#ifdef JMAPPLIB
    use conv_kf_main, only: &
       conv_kf_run
    use conv_kf_main1701, only: &
       conv_kf_run1701
    use pbl_diag, only: &
       pbl_diag_buoyancy_excess
#endif
    implicit none

    integer, intent(in) :: KA, KS, KE
    integer, intent(in) :: IA, IS, IE
    integer, intent(in) :: JA, JS, JE

    real(RP), intent(in) :: DENS (KA,IA,JA)
    real(RP), intent(in) :: U    (KA,IA,JA)
    real(RP), intent(in) :: V    (KA,IA,JA)
    real(RP), intent(in) :: W    (KA,IA,JA)
    real(RP), intent(in) :: TEMP (KA,IA,JA)
    real(RP), intent(in) :: POTT (KA,IA,JA)
    real(RP), intent(in) :: PRES (KA,IA,JA)
    real(RP), intent(in) :: EXNER(KA,IA,JA)
    real(RP), intent(in) :: QDRY (KA,IA,JA)
    real(RP), intent(in) :: QV   (KA,IA,JA)
    real(RP), intent(in) :: QC   (KA,IA,JA)
    real(RP), intent(in) :: QI   (KA,IA,JA)
    real(RP), intent(in) :: us       (IA,JA)
    real(RP), intent(in) :: PBLH     (IA,JA)
    real(RP), intent(in) :: SFLX_BUOY(IA,JA)
    real(RP), intent(in) :: CZ   (KA,IA,JA)
    real(RP), intent(in) :: FZ (0:KA,IA,JA)

    real(RP), intent(out)   :: DENS_t (KA,IA,JA)
    real(RP), intent(inout) :: RHOT_t (KA,IA,JA)
    real(RP), intent(inout) :: RHOQV_t(KA,IA,JA)
    real(RP), intent(inout) :: RHOQ_t (KA,IA,JA,N_HYD)
    real(RP), intent(inout) :: w0avg  (KA,IA,JA)
    real(RP), intent(inout) :: nca      (IA,JA)
    real(RP), intent(inout) :: SFLX_rain(IA,JA)
    real(RP), intent(inout) :: SFLX_snow(IA,JA)
    real(RP), intent(inout) :: SFLX_prec(IA,JA)
    real(RP), intent(inout) :: cloudtop (IA,JA)
    real(RP), intent(inout) :: cloudbase(IA,JA)

#ifdef JMAPPLIB

    integer, parameter :: iconv_run = 1

    real(RP) :: tend_pt(KS:KE)
    real(RP) :: tend_qv(KS:KE)
    real(RP) :: tend_qc(KS:KE)
    real(RP) :: tend_qr(KS:KE)
    real(RP) :: tend_qi(KS:KE)
    real(RP) :: tend_qs(KS:KE)
    real(RP) :: tend_qg(KS:KE)
    real(RP) :: dens_lc(KS:KE)
    real(RP) :: qv_lc  (KS:KE)
    real(RP) :: qc_lc
    real(RP) :: qi_lc
    real(RP) :: z_f    (KS:KE)
    real(RP) :: dz_f   (KS:KE)

    real(RP) :: lcl_temp(IA,JA)
    real(RP) :: abe_out
    real(RP) :: tv_ex
    real(RP) :: ptv_ex

    integer :: k, i, j
    !---------------------------------------------------------------------------

    LOG_PROGRESS(*) "atmosphere / physics / cumulus / KF-JMAPPLIB"

    !$omp parallel do &
    !$omp private(tend_pt,tend_qv,tend_qc,tend_qr,tend_qi,tend_qs,tend_qg, &
    !$omp         dens_lc,qv_lc,qc_lc,qi_lc,z_f,dz_f, &
    !$omp         ptv_ex,tv_ex,abe_out)
    do j = JS, JE
    do i = IS, IE

       do k = KS, KE
          dens_lc(k) = DENS(k,i,j) * QDRY(k,i,j)
          qv_lc  (k) = QV(k,i,j) / QDRY(k,i,j)
          tend_pt(k) = RHOT_t(k,i,j) / DENS(k,i,j)
          tend_qv(k) = RHOQV_t(k,i,j) / dens_lc(k)
          tend_qc(k) = RHOQ_t(k,i,j,I_HC) / dens_lc(k)
          tend_qr(k) = RHOQ_t(k,i,j,I_HR) / dens_lc(k)
          tend_qi(k) = RHOQ_t(k,i,j,I_HI) / dens_lc(k)
          tend_qs(k) = RHOQ_t(k,i,j,I_HS) / dens_lc(k)
          tend_qg(k) = RHOQ_t(k,i,j,I_HG) / dens_lc(k)

          z_f (k) = CZ(k,i,j) - FZ(KS-1,i,j)
          dz_f(k) = FZ(k,i,j) - FZ(k-1,i,j)
       end do

       select case ( ATMOS_PHY_CP_KF_JMAPPLIB_TYPE )
       case ( "KF" )

          call conv_kf_run( &
               iconv_run,                                                           & ! (in)
               U(KS:KE,i,j), V(KS:KE,i,j), W(KS:KE,i,j), TEMP(KS:KE,i,j), qv_lc(:), & ! (in)
               PRES(KS:KE,i,j), EXNER(KS:KE,i,j), dens_lc(:), z_f(:), dz_f(:),      & ! (in)
               nca(i,j),                                                            & ! (inout)
               SFLX_prec(i,j), SFLX_rain(i,j), SFLX_snow(i,j),                      & ! (inout)
               cloudtop(i,j), cloudbase(i,j), abe_out,                              & ! (inout)
               ishall_counter(i,j), ideep_counter(i,j),                             & ! (inout)
               lcl_temp(i,j),                                                       & ! (inout)
               w0avg(KS:KE,i,j),                                                    & ! (inout)
               tend_pt(:), tend_qv(:), tend_qc(:),                                  & ! (inout)
               tend_qi(:), tend_qr(:), tend_qs(:)                                   ) ! (inout)

          tend_qg(KS:KE) = 0.0_RP

       case ( "KF1701" )

          qc_lc = QC(KS,i,j) / QDRY(KS,i,j)
          qi_lc = QI(KS,i,j) / QDRY(KS,i,j)
          call pbl_diag_buoyancy_excess( &
               POTT(KS,i,j), qv_lc(KS),qc_lc, qi_lc,              & ! (in)
               PBLH(i,j), us(i,j), SFLX_BUOY(i,j), EXNER(KS,i,j), & ! (in)
               ptv_ex, tv_ex                                      ) ! (out)

          call conv_kf_run1701( &
               iconv_run,                                                           & ! (in)
               U(KS:KE,i,j), V(KS:KE,i,j), W(KS:KE,i,j), TEMP(KS:KE,i,j), qv_lc(:), & ! (in)
               PRES(KS:KE,i,j), EXNER(KS:KE,i,j), dens_lc(:), z_f(:), dz_f(:),      & ! (in)
               PBLH(i,j), tv_ex,                                                    & ! (in)
               nca(i,j),                                                            & ! (inout)
               SFLX_prec(i,j), SFLX_rain(i,j), SFLX_snow(i,j),                      & ! (inout)
               cloudtop(i,j), cloudbase(i,j), abe_out,                              & ! (inout)
               ishall_counter(i,j), ideep_counter(i,j),                             & ! (inout)
               lcl_temp(i,j),                                                       & ! (inout)
               w0avg(KS:KE,i,j),                                                    & ! (inout)
               tend_pt(:), tend_qv(:), tend_qc(:),                                  & ! (inout)
               tend_qi(:), tend_qr(:),                                              & ! (inout)
               tend_qs(:), tend_qg(:)                                               ) ! (inout)

       end select

       if ( nca(i,j) < -1 ) nca(i,j) = -1

       do k = KS, KE
          RHOT_t (k,i,j)      = tend_pt(k) * DENS(k,i,j)
          RHOQV_t(k,i,j)      = tend_qv(k) * dens_lc(k)
          RHOQ_t (k,i,j,I_HC) = tend_qc(k) * dens_lc(k)
          RHOQ_t (k,i,j,I_HR) = tend_qr(k) * dens_lc(k)
          RHOQ_t (k,i,j,I_HI) = tend_qi(k) * dens_lc(k)
          RHOQ_t (k,i,j,I_HS) = tend_qs(k) * dens_lc(k)
          RHOQ_t (k,i,j,I_HG) = tend_qg(k) * dens_lc(k)
          DENS_t (k,i,j)      = RHOQV_t(k,i,j) &
                              + RHOQ_t(k,i,j,I_HC) + RHOQ_t(k,i,j,I_HR) &
                              + RHOQ_t(k,i,j,I_HI) + RHOQ_t(k,i,j,I_HS) + RHOQ_t(k,i,j,I_HG)
       end do

    end do
    end do

#endif

    return
  end subroutine ATMOS_PHY_CP_KF_JMAPPLIB_tendency

end module scale_atmos_phy_cp_kf_jmapplib
