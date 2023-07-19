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
  public :: ATMOS_PHY_BL_driver_finalize
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
    use scale_atmos_phy_bl_mynn_jmapplib, only: &
       ATMOS_PHY_BL_MYNN_JMAPPLIB_NTRACER, &
       ATMOS_PHY_BL_MYNN_JMAPPLIB_NAME, &
       ATMOS_PHY_BL_MYNN_JMAPPLIB_DESC, &
       ATMOS_PHY_BL_MYNN_JMAPPLIB_UNITS
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
       case ( 'MYNN-JMAPPLIB' )
          call TRACER_regist( &
               QS, &
               ATMOS_PHY_BL_MYNN_JMAPPLIB_NTRACER, &
               ATMOS_PHY_BL_MYNN_JMAPPLIB_NAME,    &
               ATMOS_PHY_BL_MYNN_JMAPPLIB_DESC,    &
               ATMOS_PHY_BL_MYNN_JMAPPLIB_UNITS    )
          QE = QS + ATMOS_PHY_BL_MYNN_JMAPPLIB_NTRACER - 1
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
    use scale_atmos_phy_bl_mynn, only: &
       ATMOS_PHY_BL_MYNN_setup
    use scale_atmos_phy_bl_mynn_jmapplib, only: &
       ATMOS_PHY_BL_MYNN_JMAPPLIB_setup
    use scale_atmos_grid_cartesC, only: &
       CZ => ATMOS_GRID_CARTESC_CZ
    use mod_atmos_admin, only: &
       ATMOS_PHY_BL_TYPE, &
       ATMOS_sw_phy_bl
    use mod_atmos_phy_bl_vars, only: &
       ATMOS_PHY_BL_Zi, &
       ATMOS_PHY_BL_SFLX_BUOY
    use scale_bulkflux, only: &
       BULKFLUX_type
    use scale_time, only: &
       dt_BL => TIME_DTSEC_ATMOS_PHY_BL
    implicit none
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("ATMOS_PHY_BL_driver_setup",*) 'Setup'

    if ( ATMOS_sw_phy_bl ) then
       select case ( ATMOS_PHY_BL_TYPE )
       case ( 'MYNN' )
          call ATMOS_PHY_BL_MYNN_setup( &
               KA, KS, KE, IA, IS, IE, JA, JS, JE, &
               BULKFLUX_type ) ! (in)
       case ( 'MYNN-JMAPPLIB' )
          call ATMOS_PHY_BL_MYNN_JMAPPLIB_setup( &
               KA, KS, KE, &
               CZ(:), &
               dt_BL )
       end select
    else
       LOG_INFO("ATMOS_PHY_BL_driver_setup",*) 'this component is never called.'
       !$acc kernels
       ATMOS_PHY_BL_Zi       (:,:) = 0.0_RP
       ATMOS_PHY_BL_SFLX_BUOY(:,:) = 0.0_RP
       !$acc end kernels
    endif

    return
  end subroutine ATMOS_PHY_BL_driver_setup

  !-----------------------------------------------------------------------------
  !> Finalize
  subroutine ATMOS_PHY_BL_driver_finalize
    use scale_prc, only: &
       PRC_abort
    use scale_atmos_phy_bl_mynn, only: &
       ATMOS_PHY_BL_MYNN_finalize
    use scale_atmos_phy_bl_mynn_jmapplib, only: &
       ATMOS_PHY_BL_MYNN_JMAPPLIB_finalize
    use mod_atmos_admin, only: &
       ATMOS_PHY_BL_TYPE, &
       ATMOS_sw_phy_bl
    implicit none
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("ATMOS_PHY_BL_driver_finalize",*) 'Finalize'

    if ( ATMOS_sw_phy_bl ) then
       select case ( ATMOS_PHY_BL_TYPE )
       case ( 'MYNN' )
          call ATMOS_PHY_BL_MYNN_finalize
       case ( 'MYNN-JMAPPLIB' )
          call ATMOS_PHY_BL_MYNN_JMAPPLIB_finalize
       case default
          LOG_ERROR("ATMOS_PHY_BL_driver_finalize",*) 'ATMOS_PHY_BL_TYPE is invalid: ', trim(ATMOS_PHY_BL_TYPE)
          call PRC_abort
       end select
    end if

    return
  end subroutine ATMOS_PHY_BL_driver_finalize

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
    use scale_atmos_phy_bl_common, only: &
       ATMOS_PHY_BL_tendency_tracer
    use scale_atmos_phy_bl_mynn, only: &
       ATMOS_PHY_BL_MYNN_tendency
    use scale_atmos_phy_bl_mynn_jmapplib, only: &
       ATMOS_PHY_BL_MYNN_JMAPPLIB_tendency
    use scale_atmos_grid_cartesC_real, only: &
       CZ => ATMOS_GRID_CARTESC_REAL_CZ, &
       FZ => ATMOS_GRID_CARTESC_REAL_FZ, &
       F2H => ATMOS_GRID_CARTESC_REAL_F2H, &
       ATMOS_GRID_CARTESC_REAL_VOL, &
       ATMOS_GRID_CARTESC_REAL_TOTVOL
    use scale_atmos_hydrometeor, only: &
       I_QV
    use scale_bulkflux, only: &
       BULKFLUX_type
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
       ATMOS_PHY_BL_MIX_TRACERS, &
       RHOU_t_BL => ATMOS_PHY_BL_RHOU_t,    &
       RHOV_t_BL => ATMOS_PHY_BL_RHOV_t,    &
       RHOT_t_BL => ATMOS_PHY_BL_RHOT_t,    &
       RHOQ_t_BL => ATMOS_PHY_BL_RHOQ_t,    &
       Zi        => ATMOS_PHY_BL_Zi,        &
       SFLX_BUOY => ATMOS_PHY_BL_SFLX_BUOY, &
       QL        => ATMOS_PHY_BL_QL,        &
       cldfrac   => ATMOS_PHY_BL_cldfrac
    use mod_atmos_phy_sf_vars, only: &
       SFC_DENS => ATMOS_PHY_SF_SFC_DENS,  &
       SFC_PRES => ATMOS_PHY_SF_SFC_PRES,  &
       SFLX_MU  => ATMOS_PHY_SF_SFLX_MU,   &
       SFLX_MV  => ATMOS_PHY_SF_SFLX_MV,   &
       SFLX_SH  => ATMOS_PHY_SF_SFLX_SH,   &
       SFLX_Q   => ATMOS_PHY_SF_SFLX_QTRC, &
       SFLX_QV  => ATMOS_PHY_SF_SFLX_QV,   &
       Ustar    => ATMOS_PHY_SF_Ustar,     &
       Tstar    => ATMOS_PHY_SF_Tstar,     &
       Qstar    => ATMOS_PHY_SF_Qstar,     &
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

       !$acc data create(Nu,Kh,QW,N2,POTL,POTV)

       !$acc kernels
       RHOQ_t_BL(:,:,:,:) = 0.0_RP
       !$acc end kernels

       select case ( ATMOS_PHY_BL_TYPE )
       case ( 'MYNN' )
          call ATMOS_vars_get_diagnostic( "N2",   N2   )
          call ATMOS_vars_get_diagnostic( "POTL", POTL )
          call ATMOS_vars_get_diagnostic( "POTV", POTV )
          !$acc kernels
          do j = JSB, JEB
          do i = ISB, IEB
          do k = KS, KE
             QW(k,i,j) = QV(k,i,j) + QC(k,i,j) + QI(k,i,j)
          end do
          end do
          end do
          !$acc end kernels
          if ( I_QV > 0 ) then
             RHOQV_t => RHOQ_t_BL(:,:,:,I_QV)
          else
             allocate( RHOQV_t(KA,IA,JA) )
             !$acc enter data create(RHOQV_t)
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
               Ustar(:,:), Tstar(:,:), Qstar(:,:), RLmo(:,:),          & ! (in)
               CZ(:,:,:), FZ(:,:,:), F2H(:,:,:,:), dt_BL,              & ! (in)
               BULKFLUX_type,                                          & ! (in)
               RHOU_t_BL(:,:,:), RHOV_t_BL(:,:,:), RHOT_t_BL(:,:,:),   & ! (out)
               RHOQV_t(:,:,:), RHOQ_t_BL(:,:,:,QS:QE),                 & ! (out)
               Nu(:,:,:), Kh(:,:,:),                                   & ! (out)
               QL(:,:,:), cldfrac(:,:,:), Zi(:,:), SFLX_BUOY(:,:)      ) ! (out)
          if ( I_QV <= 0 ) then
             !$acc exit data delete(RHOQV_t)
             deallocate( RHOQV_t )
          end if
       case ( 'MYNN-JMAPPLIB' )
          if ( I_QV > 0 ) then
             RHOQV_t => RHOQ_t_BL(:,:,:,I_QV)
          else
             allocate( RHOQV_t(KA,IA,JA) )
             !$acc enter data create(RHOQV_t)
          end if
          !$acc update host(DENS,U,V,POTT,QTRC(:,:,:,QS:QE),PRES,EXNER,QDRY,QV,QC,QI,SFC_DENS,SFC_PRES,SFLX_MU,SFLX_MV,SFLX_SH,SFLX_QV,Ustar,RLmo)
          call ATMOS_PHY_BL_MYNN_JMAPPLIB_tendency( &
               KA, KS, KE, IA, IS, IE, JA, JS, JE, &
               DENS(:,:,:), U(:,:,:), V(:,:,:),                        & ! (in)
               POTT(:,:,:), QTRC(:,:,:,QS:QE),                         & ! (in)
               PRES(:,:,:), EXNER(:,:,:),                              & ! (in)
               QDRY(:,:,:), QV(:,:,:), QC(:,:,:), QI(:,:,:),           & ! (in)
               SFC_DENS(:,:), SFC_PRES(:,:),                           & ! (in)
               SFLX_MU(:,:), SFLX_MV(:,:), SFLX_SH(:,:), SFLX_QV(:,:), & ! (in)
               Ustar(:,:), RLmo(:,:),                                  & ! (in)
               CZ(:,:,:), FZ(:,:,:), F2H(:,:,:,:), dt_BL,              & ! (in)
               RHOU_t_BL(:,:,:), RHOV_t_BL(:,:,:), RHOT_t_BL(:,:,:),   & ! (out)
               RHOQV_t(:,:,:), RHOQ_t_BL(:,:,:,QS:QE),                 & ! (out)
               Nu(:,:,:), Kh(:,:,:),                                   & ! (out)
               Zi = Zi(:,:), SFLX_BUOY = SFLX_BUOY(:,:)                ) ! (out)
          !$acc update device(RHOU_t_BL,RHOV_t_BL,RHOT_t_BL,RHOQV_t,RHOQ_t_BL,Nu,Kh,Zi,SFLX_BUOY)
          if ( I_QV <= 0 ) then
             !$acc exit data delete(RHOQV_t)
             deallocate( RHOQV_t )
          end if
       end select

       if ( ATMOS_PHY_BL_MIX_TRACERS ) then
          do iq = 1, QA
             if ( ( .not. TRACER_ADVC(iq) ) .or. iq==I_QV .or. (iq>=QS .and. iq<=QE) ) cycle
             call ATMOS_PHY_BL_tendency_tracer( &
                  KA, KS, KE, IA, IS, IE, JA, JS, JE, &
                  DENS(:,:,:), QTRC(:,:,:,iq),        & ! (in)
                  SFLX_Q(:,:,iq),                     & ! (in)
                  Kh(:,:,:),                          & ! (in)
                  TRACER_MASS(iq),                    & ! (in)
                  CZ(:,:,:), FZ(:,:,:), F2H(:,:,:,:), & ! (in)
                  dt_BL, TRACER_NAME(iq),             & ! (in)
                  RHOQ_t_BL(:,:,:,iq)                 ) ! (out)
          end do
       end if

       call FILE_HISTORY_in( Nu   (:,:,:),     'Nu_BL',     'eddy viscosity',                         'm2/s',      fill_halo=.true., dim_type="ZHXY" )
       call FILE_HISTORY_in( Kh   (:,:,:),     'Kh_BL',     'eddy diffusion',                         'm2/s',      fill_halo=.true., dim_type="ZHXY" )

       call FILE_HISTORY_in( QL    (:,:,:),    'QL_BL',      'liquid water content in partial condensation', 'kg/kg', fill_halo=.true. )
       call FILE_HISTORY_in( cldfrac(:,:,:),   'cldfrac_BL', 'cloud fraction in partial condensation',       '1',     fill_halo=.true. )

       call FILE_HISTORY_in( Zi     (:,:),     'Zi_BL',      'depth of the boundary layer', 'm', fill_halo=.true. )

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

       !$acc end data

    endif

    !$omp parallel do default(none) private(i,j,k) OMP_SCHEDULE_ collapse(2) &
    !$omp shared(JS,JE,IS,IE,KS,KE,RHOU_t,RHOU_t_BL,RHOV_t,RHOV_t_BL,RHOT_t,RHOT_t_BL)
    !$acc kernels
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       RHOU_t(k,i,j) = RHOU_t(k,i,j) + RHOU_t_BL(k,i,j)
       RHOV_t(k,i,j) = RHOV_t(k,i,j) + RHOV_t_BL(k,i,j)
       RHOT_t(k,i,j) = RHOT_t(k,i,j) + RHOT_t_BL(k,i,j)
    enddo
    enddo
    enddo
    !$acc end kernels

    !$omp parallel private(iq,i,j,k)
    !$acc kernels
    do iq = 1,  QA
#ifndef _OPENACC
       if ( .not. TRACER_ADVC(iq) ) cycle
#endif
       !$omp do OMP_SCHEDULE_ collapse(2)
       do j  = JS, JE
       do i  = IS, IE
       do k  = KS, KE
#ifdef _OPENACC
          if ( TRACER_ADVC(iq) ) &
#endif
          RHOQ_t(k,i,j,iq) = RHOQ_t(k,i,j,iq) + RHOQ_t_BL(k,i,j,iq)
       enddo
       enddo
       enddo
       !$omp end do nowait
    enddo
    !$acc end kernels
    !$omp end parallel

    return
  end subroutine ATMOS_PHY_BL_driver_calc_tendency

end module mod_atmos_phy_bl_driver
