!-------------------------------------------------------------------------------
!> module atmosphere / physics / PBL
!!
!! @par Description
!!          Planetary boundary layer turbulence
!!
!! @author Team SCALE
!!
!<
!-------------------------------------------------------------------------------
#include "scalelib.h"
module mod_atmos_phy_bl_driver
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_io
  use scale_prof
  use scale_atmos_grid_icoA_index
  use scale_tracer
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: ATMOS_PHY_BL_driver_tracer_setup
  public :: ATMOS_PHY_BL_driver_setup
  public :: ATMOS_PHY_BL_driver_step

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
  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Config
  subroutine ATMOS_PHY_BL_driver_tracer_setup
    use scale_prc, only: &
       PRC_abort
    use scale_tracer, only: &
       TRACER_regist
    use scale_atmos_phy_bl_mynn, only: &
       ATMOS_PHY_BL_MYNN_tracer_setup, &
       ATMOS_PHY_BL_MYNN_NTRACER, &
       ATMOS_PHY_BL_MYNN_NAME, &
       ATMOS_PHY_BL_MYNN_DESC, &
       ATMOS_PHY_BL_MYNN_UNITS
    use mod_atmos_admin, only: &
       ATMOS_PHY_BL_TYPE, &
       ATMOS_sw_phy_bl
    use mod_atmos_phy_bl_vars, only: &
       QS, QE
    implicit none
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("ATMOS_PHY_BL_driver_tracer_setup",*) 'Setup'

    if ( ATMOS_sw_phy_bl ) then
       select case ( ATMOS_PHY_BL_TYPE )
       case ( 'MYNN' )
          call ATMOS_PHY_BL_MYNN_tracer_setup
          call TRACER_regist( &
               QS, &
               ATMOS_PHY_BL_MYNN_NTRACER, &
               ATMOS_PHY_BL_MYNN_NAME,    &
               ATMOS_PHY_BL_MYNN_DESC,    &
               ATMOS_PHY_BL_MYNN_UNITS    )
          QE = QS + ATMOS_PHY_BL_MYNN_NTRACER - 1
       case default
          LOG_ERROR("ATMOS_PHY_BL_driver_tracer_setup",*) 'ATMOS_PHY_BL_TYPE is invalid: ', trim(ATMOS_PHY_BL_TYPE)
          call PRC_abort
       end select
    end if

    return
  end subroutine ATMOS_PHY_BL_driver_tracer_setup

  !-----------------------------------------------------------------------------
  !> Setup
  subroutine ATMOS_PHY_BL_driver_setup
    use mod_atmos_vars, only: &
         CZ
    use scale_atmos_phy_bl_mynn, only: &
       ATMOS_PHY_BL_MYNN_setup
    use mod_atmos_admin, only: &
       ATMOS_PHY_BL_TYPE, &
       ATMOS_sw_phy_bl
    implicit none

    real(RP) :: CZ2(KA,IA,JA)

    integer :: k, i, j, l
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("ATMOS_PHY_BL_driver_setup",*) 'Setup'

    if ( ATMOS_sw_phy_bl ) then
       select case ( ATMOS_PHY_BL_TYPE )
       case ( 'MYNN' )
          do j = JS, JE
          do i = IS, IE
          do k = KS, KE
             CZ2(k,i,j) = minval( CZ(k,i,j,:) )
          end do
          end do
          end do
          call ATMOS_PHY_BL_MYNN_setup( &
               KA, KS, KE, IA, IS, IE, JA, JS, JE, &
               CZ2(:,:,:) ) ! (in)
       end select
    else
       LOG_INFO("ATMOS_PHY_BL_driver_setup",*) 'this component is never called.'
    endif

    return
  end subroutine ATMOS_PHY_BL_driver_setup

  !-----------------------------------------------------------------------------
  !> time step
  subroutine ATMOS_PHY_BL_driver_step
    use scale_const, only: &
       PRE00 => CONST_PRE00
    use scale_file_history, only: &
       FILE_HISTORY_in
    use scale_time, only: &
       dt_BL => TIME_DTSEC_ATMOS_PHY_BL
    use scale_atmos_phy_bl_mynn, only: &
       ATMOS_PHY_BL_MYNN_tendency, &
       ATMOS_PHY_BL_MYNN_tendency_tracer
    use mod_atmos_admin, only: &
       ATMOS_PHY_BL_TYPE, &
       ATMOS_sw_phy_bl
    use mod_atmos_vars, only: &
       ATMOS_vars_calc_diagnostics, &
       DENS, &
       RHOU, &
       RHOV, &
       RHOE, &
       RHOQ, &
       QTRC, &
       U,     &
       V,     &
       TEMP,  &
       POTT,  &
       PRES,  &
       EXNER, &
       QDRY,  &
       Rtot,  &
       CVtot, &
       CPtot, &
       QV,    &
       QC,    &
       QI,    &
       ATMOS_vars_get_diagnostic
    use mod_atmos_phy_bl_vars, only: &
       QS, QE
    use mod_atmos_phy_sf_vars, only: &
       SFLX_MU => ATMOS_PHY_SF_SFLX_MU, &
       SFLX_MV => ATMOS_PHY_SF_SFLX_MV, &
       SFLX_SH => ATMOS_PHY_SF_SFLX_SH, &
       SFLX_Q  => ATMOS_PHY_SF_SFLX_QTRC, &
       SFLX_QV => ATMOS_PHY_SF_SFLX_QV, &
       Ustar   => ATMOS_PHY_SF_Ustar, &
       RLmo    => ATMOS_PHY_SF_RLmo
    use mod_atmos_vars, only: &
       CZ, &
       FZ
    implicit none

    real(RP) :: Nu   (KA,IA,JA,ADM_lall) !> eddy viscosity
    real(RP) :: Nu_cg(KA,IA,JA,ADM_lall) !> eddy viscosity for the counter gradient
    real(RP) :: Kh   (KA,IA,JA,ADM_lall) !> eddy diffution
    real(RP) :: Kh_cg(KA,IA,JA,ADM_lall) !> eddy diffution for the counter gradient
    real(RP) :: QW(KA,IA,JA)          !> total water

    real(RP) :: N2  (KA,IA,JA,ADM_lall) !> static stability
    real(RP) :: POTL(KA,IA,JA,ADM_lall) !> liquid water potential temperature
    real(RP) :: POTV(KA,IA,JA,ADM_lall) !> virtual potential temperature

    real(RP) :: RHOU_t(KA,IA,JA)
    real(RP) :: RHOV_t(KA,IA,JA)
    real(RP) :: RHOT_t(KA,IA,JA)
    real(RP) :: RHOQ_t(KA,IA,JA,QA)

    integer  :: k, i, j, iq, l
    !---------------------------------------------------------------------------

    select case ( ATMOS_PHY_BL_TYPE )
    case ( 'MYNN' )
       call ATMOS_vars_get_diagnostic( "N2",   N2   )
       call ATMOS_vars_get_diagnostic( "POTL", POTL )
       call ATMOS_vars_get_diagnostic( "POTV", POTV )
       do l = 1, ADM_lall

          do j = JS, JE
          do i = IS, IE
          do k = KS, KE
             QW(k,i,j) = QV(k,i,j,l) + QC(k,i,j,l) + QI(k,i,j,l)
          end do
          end do
          end do
          call ATMOS_PHY_BL_MYNN_tendency( &
               KA, KS, KE, IA, IS, IE, JA, JS, JE, &
               DENS(:,:,:,l), U(:,:,:,l), V(:,:,:,l),                          & ! (in)
               POTT(:,:,:,l), QTRC(:,:,:,QS:QE,l),                             & ! (in)
               PRES(:,:,:,l), EXNER(:,:,:,l), N2(:,:,:,l),                     & ! (in)
               QDRY(:,:,:,l), QV(:,:,:,l), QW(:,:,:),                          & ! (in)
               POTL(:,:,:,l), POTV(:,:,:,l),                                   & ! (in)
               SFLX_MU(:,:,l), SFLX_MV(:,:,l), SFLX_SH(:,:,l), SFLX_QV(:,:,l), & ! (in)
               Ustar(:,:,l), RLmo(:,:,l),                                      & ! (in)
               CZ(:,:,:,l), FZ(:,:,:,l), dt_BL,                                & ! (in)
               RHOU_t(:,:,:), RHOV_t(:,:,:),                                   & ! (out)
               RHOT_t(:,:,:), RHOQ_t(:,:,:,QS:QE),                             & ! (out)
               Nu(:,:,:,l), Nu_cg(:,:,:,l), Kh(:,:,:,l), Kh_cg(:,:,:,l)        ) ! (out)
          do j = JS, JE
          do i = IS, IE
          do k = KS, KE
             RHOU(k,i,j,l) = RHOU(k,i,j,l) + RHOU_t(k,i,j) * dt_BL
             RHOV(k,i,j,l) = RHOV(k,i,j,l) + RHOV_t(k,i,j) * dt_BL
             RHOE(k,i,j,l) = RHOE(k,i,j,l) + RHOT_t(k,i,j) * CPtot(k,i,j,l) * EXNER(k,i,j,l) * dt_BL
             RHOQ(k,i,j,QS:QE,l) = RHOQ(k,i,j,QS:QE,l) + RHOQ_t(k,i,j,QS:QE) * dt_BL
          end do
          end do
          end do

          do iq = 1, QA
             if ( ( .not. TRACER_ADVC(iq) ) .or. (iq>=QS .and. iq<=QE) ) cycle
             call ATMOS_PHY_BL_MYNN_tendency_tracer( &
                  KA, KS, KE, IA, IS, IE, JA, JS, JE, &
                  DENS(:,:,:,l), QTRC(:,:,:,iq,l), SFLX_Q(:,:,iq,l), &
                  Kh(:,:,:,l), Kh_cg(:,:,:,l), TRACER_MASS(iq),      & ! (in)
                  CZ(:,:,:,l), FZ(:,:,:,l),                          & ! (in)
                  dt_BL, TRACER_NAME(iq),                            & ! (in)
                  RHOQ_t(:,:,:,iq)                                   ) ! (out)
             do j = JS, JE
             do i = IS, IE
             do k = KS, KE
                RHOQ(k,i,j,iq,l) = RHOQ(k,i,j,iq,l) + RHOQ_t(k,i,j,iq) * dt_BL
             end do
             end do
             end do
          end do

       end do
    end select

    call FILE_HISTORY_in( Nu(:,:,:,:),        'Nu_BL',     'eddy viscosity',     'm2/s',      fill_halo=.true. )
    call FILE_HISTORY_in( Kh(:,:,:,:),        'Kh_BL',     'eddy diffusion',     'm2/s',      fill_halo=.true. )

    call ATMOS_vars_calc_diagnostics

    return
  end subroutine ATMOS_PHY_BL_driver_step

end module mod_atmos_phy_bl_driver
