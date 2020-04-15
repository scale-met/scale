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
  use scale_atmos_grid_cartesC_index
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
  public :: ATMOS_PHY_BL_driver_calc_tendency

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
    use scale_atmos_grid_cartesC, only: &
       CDZ => ATMOS_GRID_CARTESC_CDZ, &
       CDX => ATMOS_GRID_CARTESC_CDX, &
       CDY => ATMOS_GRID_CARTESC_CDY
    use scale_atmos_grid_cartesC_real, only: &
       REAL_CZ  => ATMOS_GRID_CARTESC_REAL_CZ
    use scale_atmos_phy_bl_mynn, only: &
       ATMOS_PHY_BL_MYNN_setup
    use mod_atmos_admin, only: &
       ATMOS_PHY_BL_TYPE, &
       ATMOS_sw_phy_bl
    implicit none
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("ATMOS_PHY_BL_driver_setup",*) 'Setup'

    if ( ATMOS_sw_phy_bl ) then
       select case ( ATMOS_PHY_BL_TYPE )
       case ( 'MYNN' )
          call ATMOS_PHY_BL_MYNN_setup( &
               KA, KS, KE, IA, IS, IE, JA, JS, JE, &
               REAL_CZ ) ! (in)
       end select
    else
       LOG_INFO("ATMOS_PHY_BL_driver_setup",*) 'this component is never called.'
    endif

    return
  end subroutine ATMOS_PHY_BL_driver_setup

  !-----------------------------------------------------------------------------
  !> calculate tendency
  subroutine ATMOS_PHY_BL_driver_calc_tendency( update_flag )
    use scale_statistics, only: &
       STATISTICS_checktotal, &
       STATISTICS_total
    use scale_file_history, only: &
       FILE_HISTORY_in
    use scale_time, only: &
       dt_BL => TIME_DTSEC_ATMOS_PHY_BL
    use scale_atmos_phy_bl_mynn, only: &
       ATMOS_PHY_BL_MYNN_tendency, &
       ATMOS_PHY_BL_MYNN_tendency_tracer
    use scale_atmos_grid_cartesC_real, only: &
       CZ => ATMOS_GRID_CARTESC_REAL_CZ, &
       FZ => ATMOS_GRID_CARTESC_REAL_FZ, &
       ATMOS_GRID_CARTESC_REAL_VOL, &
       ATMOS_GRID_CARTESC_REAL_TOTVOL
    use scale_atmos_hydrometeor, only: &
       I_QV
    use mod_atmos_admin, only: &
       ATMOS_PHY_BL_TYPE, &
       ATMOS_sw_phy_bl
    use mod_atmos_vars, only: &
       DENS => DENS_av, &
       QTRC => QTRC_av, &
       U,     &
       V,     &
       POTT,  &
       PRES,  &
       EXNER, &
       QDRY,  &
       QV,    &
       QC,    &
       QI,    &
       RHOU_t => RHOU_tp, &
       RHOV_t => RHOV_tp, &
       RHOT_t => RHOT_tp, &
       RHOQ_t => RHOQ_tp, &
       ATMOS_vars_get_diagnostic
    use mod_atmos_phy_bl_vars, only: &
       QS, QE, &
       RHOU_t_BL => ATMOS_PHY_BL_RHOU_t, &
       RHOV_t_BL => ATMOS_PHY_BL_RHOV_t, &
       RHOT_t_BL => ATMOS_PHY_BL_RHOT_t, &
       RHOQ_t_BL => ATMOS_PHY_BL_RHOQ_t
    use mod_atmos_phy_sf_vars, only: &
       SFC_DENS => ATMOS_PHY_SF_SFC_DENS,  &
       SFLX_MU  => ATMOS_PHY_SF_SFLX_MU,   &
       SFLX_MV  => ATMOS_PHY_SF_SFLX_MV,   &
       SFLX_SH  => ATMOS_PHY_SF_SFLX_SH,   &
       SFLX_Q   => ATMOS_PHY_SF_SFLX_QTRC, &
       SFLX_QV  => ATMOS_PHY_SF_SFLX_QV,   &
       Ustar    => ATMOS_PHY_SF_Ustar,     &
       RLmo     => ATMOS_PHY_SF_RLmo
    implicit none

    logical, intent(in) :: update_flag

    real(RP) :: Nu(KA,IA,JA) !> eddy viscosity
    real(RP) :: Kh(KA,IA,JA) !> eddy diffution

    real(RP) :: QW(KA,IA,JA) !> total water

    real(RP) :: N2  (KA,IA,JA) !> static stability
    real(RP) :: POTL(KA,IA,JA) !> liquid water potential temperature
    real(RP) :: POTV(KA,IA,JA) !> virtual potential temperature

    real(RP), pointer :: RHOQV_t(:,:,:)

    integer  :: k, i, j, iq
    !---------------------------------------------------------------------------

    if ( update_flag ) then

       RHOQ_t_BL(:,:,:,:) = 0.0_RP

       select case ( ATMOS_PHY_BL_TYPE )
       case ( 'MYNN' )
          call ATMOS_vars_get_diagnostic( "N2",   N2   )
          call ATMOS_vars_get_diagnostic( "POTL", POTL )
          call ATMOS_vars_get_diagnostic( "POTV", POTV )
          do j = JSB, JEB
          do i = ISB, IEB
          do k = KS, KE
             QW(k,i,j) = QV(k,i,j) + QC(k,i,j) + QI(k,i,j)
          end do
          end do
          end do
          if ( I_QV > 0 ) then
             RHOQV_t => RHOQ_t_BL(:,:,:,I_QV)
          else
             allocate( RHOQV_T(KA,IA,JA) )
          end if
          call ATMOS_PHY_BL_MYNN_tendency( &
               KA, KS, KE, IA, IS, IE, JA, JS, JE, &
               DENS(:,:,:), U(:,:,:), V(:,:,:),                        & ! (in)
               POTT(:,:,:), QTRC(:,:,:,QS:QE),                         & ! (in)
               PRES(:,:,:), EXNER(:,:,:), N2(:,:,:),                   & ! (in)
               QDRY(:,:,:), QV(:,:,:), QW(:,:,:),                      & ! (in)
               POTL(:,:,:), POTV(:,:,:),                               & ! (in)
               SFC_DENS(:,:),                                          & ! (in)
               SFLX_MU(:,:), SFLX_MV(:,:), SFLX_SH(:,:), SFLX_QV(:,:), & ! (in)
               Ustar(:,:), RLmo(:,:),                                  & ! (in)
               CZ(:,:,:), FZ(:,:,:), dt_BL,                            & ! (in)
               RHOU_t_BL(:,:,:), RHOV_t_BL(:,:,:), RHOT_t_BL(:,:,:),   & ! (out)
               RHOQV_t(:,:,:), RHOQ_t_BL(:,:,:,QS:QE),                 & ! (out)
               Nu(:,:,:), Kh(:,:,:)                                    ) ! (out)
          if ( I_QV <= 0 ) deallocate( RHOQV_T )
          do iq = 1, QA
             if ( ( .not. TRACER_ADVC(iq) ) .or. iq==I_QV .or. (iq>=QS .and. iq<=QE) ) cycle
             call ATMOS_PHY_BL_MYNN_tendency_tracer( &
                  KA, KS, KE, IA, IS, IE, JA, JS, JE, &
                  DENS(:,:,:), QTRC(:,:,:,iq), & ! (in)
                  SFLX_Q(:,:,iq),              & ! (in)
                  Kh(:,:,:),                   & ! (in)
                  TRACER_MASS(iq),             & ! (in)
                  CZ(:,:,:), FZ(:,:,:),        & ! (in)
                  dt_BL, TRACER_NAME(iq),      & ! (in)
                  RHOQ_t_BL(:,:,:,iq)          ) ! (out)
          end do
       end select

       call FILE_HISTORY_in( Nu   (:,:,:),     'Nu_BL',     'eddy viscosity',                         'm2/s',      fill_halo=.true., dim_type="ZHXY" )
       call FILE_HISTORY_in( Kh   (:,:,:),     'Kh_BL',     'eddy diffusion',                         'm2/s',      fill_halo=.true., dim_type="ZHXY" )

       call FILE_HISTORY_in( RHOU_t_BL(:,:,:), 'RHOU_t_BL', 'MOMX tendency (BL)',                     'kg/m2/s2',  fill_halo=.true. )
       call FILE_HISTORY_in( RHOV_t_BL(:,:,:), 'RHOV_t_BL', 'MOMY tendency (BL)',                     'kg/m2/s2',  fill_halo=.true. )
       call FILE_HISTORY_in( RHOT_t_BL(:,:,:), 'RHOT_t_BL', 'RHOT tendency (BL)',                     'K.kg/m3/s', fill_halo=.true. )

       do iq = 1, QA
          if ( .not. TRACER_ADVC(iq) ) cycle
          call FILE_HISTORY_in( RHOQ_t_BL(:,:,:,iq), trim(TRACER_NAME(iq))//'_t_BL',                      &
                        'RHO*'//trim(TRACER_NAME(iq))//' tendency (BL)', 'kg/m3/s', fill_halo=.true. )
       enddo

       if ( STATISTICS_checktotal ) then
          call STATISTICS_total( KA, KS, KE, IA, IS, IE, JA, JS, JE, &
                                 RHOU_t_BL(:,:,:), 'RHOU_t_BL',      &
                                 ATMOS_GRID_CARTESC_REAL_VOL(:,:,:), &
                                 ATMOS_GRID_CARTESC_REAL_TOTVOL      )
          call STATISTICS_total( KA, KS, KE, IA, IS, IE, JA, JS, JE, &
                                 RHOV_t_BL(:,:,:), 'RHOV_t_BL',      &
                                 ATMOS_GRID_CARTESC_REAL_VOL(:,:,:), &
                                 ATMOS_GRID_CARTESC_REAL_TOTVOL      )
          call STATISTICS_total( KA, KS, KE, IA, IS, IE, JA, JS, JE, &
                                 RHOT_t_BL(:,:,:), 'RHOT_t_BL',      &
                                 ATMOS_GRID_CARTESC_REAL_VOL(:,:,:), &
                                 ATMOS_GRID_CARTESC_REAL_TOTVOL      )
          call STATISTICS_total( KA, KS, KE, IA, IS, IE, JA, JS, JE, &
                                 Nu(:,:,:),        'Nu_BL',          &
                                 ATMOS_GRID_CARTESC_REAL_VOL(:,:,:), &
                                 ATMOS_GRID_CARTESC_REAL_TOTVOL      )
          call STATISTICS_total( KA, KS, KE, IA, IS, IE, JA, JS, JE, &
                                 Kh(:,:,:),         'Kh_BL',         &
                                 ATMOS_GRID_CARTESC_REAL_VOL(:,:,:), &
                                 ATMOS_GRID_CARTESC_REAL_TOTVOL      )

          do iq = 1, QA
             if ( .not. TRACER_ADVC(iq) ) cycle
             call STATISTICS_total( KA, KS, KE, IA, IS, IE, JA, JS, JE, &
                                    RHOQ_t_BL(:,:,:,iq), trim(TRACER_NAME(iq))//'_t_BL', &
                                    ATMOS_GRID_CARTESC_REAL_VOL(:,:,:),                  &
                                    ATMOS_GRID_CARTESC_REAL_TOTVOL                       )
          enddo
       endif

    endif

    !$omp parallel do default(none) private(i,j,k) OMP_SCHEDULE_ collapse(2) &
    !$omp shared(JS,JE,IS,IE,KS,KE,RHOU_t,RHOU_t_BL,RHOV_t,RHOV_t_BL,RHOT_t,RHOT_t_BL)
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       RHOU_t(k,i,j) = RHOU_t(k,i,j) + RHOU_t_BL(k,i,j)
       RHOV_t(k,i,j) = RHOV_t(k,i,j) + RHOV_t_BL(k,i,j)
       RHOT_t(k,i,j) = RHOT_t(k,i,j) + RHOT_t_BL(k,i,j)
    enddo
    enddo
    enddo

    do iq = 1,  QA
       if ( .not. TRACER_ADVC(iq) ) cycle
       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(3)
       do j  = JS, JE
       do i  = IS, IE
       do k  = KS, KE
          RHOQ_t(k,i,j,iq) = RHOQ_t(k,i,j,iq) + RHOQ_t_BL(k,i,j,iq)
       enddo
       enddo
       enddo
    enddo

    return
  end subroutine ATMOS_PHY_BL_driver_calc_tendency

end module mod_atmos_phy_bl_driver
