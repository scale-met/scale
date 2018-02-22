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
#include "inc_openmp.h"
module mod_atmos_phy_bl_driver
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_stdio
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
  public :: ATMOS_PHY_BL_driver_resume
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
    use scale_process, only: &
       PRC_abort
    use scale_tracer, only: &
       TRACER_regist
    use scale_atmos_phy_bl_mynn, only: &
       ATMOS_PHY_BL_MYNN_NTRACER, &
       ATMOS_PHY_BL_MYNN_NAME, &
       ATMOS_PHY_BL_MYNN_DESC, &
       ATMOS_PHY_BL_MYNN_UNITS
    use mod_atmos_admin, only: &
       ATMOS_PHY_BL_TYPE, &
       ATMOS_sw_phy_bl
    use mod_atmos_phy_bl_vars, only: &
       I_TKE
    implicit none
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '++++++ Module[Tracer Setup] / Categ[ATMOS PHY_BL] / Origin[SCALE-RM]'

    if ( ATMOS_sw_phy_bl ) then
       select case ( ATMOS_PHY_BL_TYPE )
       case ( 'MYNN' )
          call TRACER_regist( &
               I_TKE, &
               ATMOS_PHY_BL_MYNN_NTRACER, &
               ATMOS_PHY_BL_MYNN_NAME,    &
               ATMOS_PHY_BL_MYNN_DESC,    &
               ATMOS_PHY_BL_MYNN_UNITS    )
       case default
          if ( IO_L ) write(IO_FID_LOG,*) '+++ ATMOS_PHY_BL_TYPE is invalid: ', trim(ATMOS_PHY_BL_TYPE)
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

    integer :: i, j
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '++++++ Module[DRIVER] / Categ[ATMOS PHY_BL] / Origin[SCALE-RM]'

    if ( ATMOS_sw_phy_bl ) then
       select case ( ATMOS_PHY_BL_TYPE )
       case ( 'MYNN' )
          call ATMOS_PHY_BL_MYNN_setup( &
               KA, KS, KE, IA, IS, IE, JA, JS, JE, &
               REAL_CZ ) ! (in)
       end select
    else
       if( IO_L ) write(IO_FID_LOG,*) '*** this component is never called.'
    endif

    return
  end subroutine ATMOS_PHY_BL_driver_setup

  !-----------------------------------------------------------------------------
  !> Resume
  subroutine ATMOS_PHY_BL_driver_resume
    use mod_atmos_admin, only: &
       ATMOS_PHY_BL_TYPE, &
       ATMOS_sw_phy_bl
    implicit none

    if ( ATMOS_sw_phy_bl ) then

       ! run once (only for the diagnostic value)
       call PROF_rapstart('ATM_Turbulence', 1)
       call ATMOS_PHY_BL_driver_calc_tendency( update_flag = .true. )
       call PROF_rapend  ('ATM_Turbulence', 1)

    end if

    return
  end subroutine ATMOS_PHY_BL_driver_resume

  !-----------------------------------------------------------------------------
  !> calculate tendency
  subroutine ATMOS_PHY_BL_driver_calc_tendency( update_flag )
    use scale_rm_statistics, only: &
       STATISTICS_checktotal, &
       STAT_total
    use scale_file_history, only: &
       FILE_HISTORY_in
    use scale_time, only: &
       dt_BL => TIME_DTSEC_ATMOS_PHY_BL
    use scale_atmos_phy_bl_mynn, only: &
       ATMOS_PHY_BL_MYNN_tendency, &
       ATMOS_PHY_BL_MYNN_tendency_tracer
    use scale_atmos_hydrometeor, only: &
       I_QV
    use scale_atmos_grid_cartesC_real, only: &
       CZ => ATMOS_GRID_CARTESC_REAL_CZ, &
       FZ => ATMOS_GRID_CARTESC_REAL_FZ
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
       I_TKE, &
       RHOU_t_BL => ATMOS_PHY_BL_RHOU_t, &
       RHOV_t_BL => ATMOS_PHY_BL_RHOV_t, &
       RHOT_t_BL => ATMOS_PHY_BL_RHOT_t, &
       RHOQ_t_BL => ATMOS_PHY_BL_RHOQ_t
    use mod_atmos_phy_sf_vars, only: &
       SFLX_MU => ATMOS_PHY_SF_SFLX_MU, &
       SFLX_MV => ATMOS_PHY_SF_SFLX_MV, &
       SFLX_SH => ATMOS_PHY_SF_SFLX_SH, &
       SFLX_Q  => ATMOS_PHY_SF_SFLX_QTRC, &
       l_mo    => ATMOS_PHY_SF_l_mo
    implicit none

    logical, intent(in) :: update_flag

    real(RP) :: Nu(KA,IA,JA) !> eddy viscosity
    real(RP) :: Kh(KA,IA,JA) !> eddy diffution
    real(RP) :: QW(KA,IA,JA) !> total water

    real(RP) :: N2  (KA,IA,JA) !> static stability
    real(RP) :: POTL(KA,IA,JA) !> liquid water potential temperature
    real(RP) :: POTV(KA,IA,JA) !> virtual potential temperature

    real(RP) :: total ! dummy

    integer  :: k, i, j, iq
    !---------------------------------------------------------------------------

    if ( update_flag ) then

       RHOQ_t_BL = 0.0_RP

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
          call ATMOS_PHY_BL_MYNN_tendency( &
               KA, KS, KE, IA, IS, IE, JA, JS, JE, &
               DENS(:,:,:), U(:,:,:), V(:,:,:),                      & ! (in)
               POTT(:,:,:), QTRC(:,:,:,I_TKE),                       & ! (in)
               PRES(:,:,:), EXNER(:,:,:), N2(:,:,:),                 & ! (in)
               QDRY(:,:,:), QV(:,:,:), QW(:,:,:),                    & ! (in)
               POTL(:,:,:), POTV(:,:,:),                             & ! (in)
               SFLX_MU(:,:), SFLX_MV(:,:), SFLX_SH(:,:), l_mo(:,:),  & ! (in)
               CZ(:,:,:), FZ(:,:,:), dt_BL,                          & ! (in)
               RHOU_t_BL(:,:,:), RHOV_t_BL(:,:,:),                   & ! (out)
               RHOT_t_BL(:,:,:), RHOQ_t_BL(:,:,:,I_TKE),             & ! (out)
               Nu(:,:,:), Kh(:,:,:)                                  ) ! (out)
          do iq = 1, QA
             if ( ( .not. TRACER_ADVC(iq) ) .or. iq==I_TKE ) cycle
             call ATMOS_PHY_BL_MYNN_tendency_tracer( &
                  KA, KS, KE, IA, IS, IE, JA, JS, JE, &
                  DENS(:,:,:), QTRC(:,:,:,iq), & ! (in)
                  SFLX_Q(:,:,iq), Kh(:,:,:),   & ! (in)
                  CZ(:,:,:), FZ(:,:,:),        & ! (in)
                  dt_BL, TRACER_NAME(iq),      & ! (in)
                  RHOQ_t_BL(:,:,:,iq)          ) ! (out)
          end do
       end select

       call FILE_HISTORY_in( Nu(:,:,:),        'Nu_BL',     'eddy viscosity',     'm2/s',      fill_halo=.true. )
       call FILE_HISTORY_in( Kh(:,:,:),        'Ku_BL',     'eddy diffusion',     'm2/s',      fill_halo=.true. )

       call FILE_HISTORY_in( RHOU_t_BL(:,:,:), 'RHOU_t_BL', 'MOMX tendency (BL)', 'kg/m2/s2',  fill_halo=.true. )
       call FILE_HISTORY_in( RHOV_t_BL(:,:,:), 'RHOV_t_BL', 'MOMY tendency (BL)', 'kg/m2/s2',  fill_halo=.true. )
       call FILE_HISTORY_in( RHOT_t_BL(:,:,:), 'RHOT_t_BL', 'RHOT tendency (BL)', 'K.kg/m3/s', fill_halo=.true. )

       do iq = 1, QA
          if ( .not. TRACER_ADVC(iq) ) cycle
          call FILE_HISTORY_in( RHOQ_t_BL(:,:,:,iq), trim(TRACER_NAME(iq))//'_t_BL',                      &
                        'RHO*'//trim(TRACER_NAME(iq))//' tendency (BL)', 'kg/m3/s', fill_halo=.true. )
       enddo

       if ( STATISTICS_checktotal ) then
          call STAT_total( total, RHOU_t_BL(:,:,:), 'RHOU_t_BL' )
          call STAT_total( total, RHOV_t_BL(:,:,:), 'RHOV_t_BL' )
          call STAT_total( total, RHOT_t_BL(:,:,:), 'RHOT_t_BL' )
          call STAT_total( total, Nu(:,:,:), 'Nu_BL' )
          call STAT_total( total, Kh(:,:,:), 'Kh_BL' )

          do iq = 1, QA
             if ( .not. TRACER_ADVC(iq) ) cycle
             call STAT_total( total, RHOQ_t_BL(:,:,:,iq), trim(TRACER_NAME(iq))//'_t_BL' )
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
