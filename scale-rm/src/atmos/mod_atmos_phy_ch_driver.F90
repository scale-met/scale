!-------------------------------------------------------------------------------
!> module ATMOSPHERE / Physics Chemistry
!!
!! @par Description
!!          Atmospheric chemistry driver
!!
!! @author Team SCALE
!!
!<
!-------------------------------------------------------------------------------
#include "scalelib.h"
module mod_atmos_phy_ch_driver
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_stdio
  use scale_prof
  use scale_atmos_grid_cartesC_index
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: ATMOS_PHY_CH_driver_tracer_setup
  public :: ATMOS_PHY_CH_driver_setup
  public :: ATMOS_PHY_CH_driver_calc_tendency

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
  subroutine ATMOS_PHY_CH_driver_tracer_setup
    use mod_atmos_admin, only: &
       ATMOS_PHY_CH_TYPE, &
       ATMOS_sw_phy_ch
    use scale_tracer, only: &
       TRACER_regist
    use scale_atmos_phy_ch_rn222, only: &
       ATMOS_PHY_CH_rn222_ntracers, &
       ATMOS_PHY_CH_rn222_NAME, &
       ATMOS_PHY_CH_rn222_DESC, &
       ATMOS_PHY_CH_rn222_UNIT
    use mod_atmos_phy_ch_vars, only: &
       QA_CH, &
       QS_CH, &
       QE_CH
    use scale_prc, only: &
       PRC_abort
    implicit none
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_PROGRESS(*) 'Module[Tracer Setup] / Categ[ATMOS PHY_CH] / Origin[SCALE-RM]'

    if ( ATMOS_sw_phy_ch ) then

       select case ( ATMOS_PHY_CH_TYPE )
       case ( 'OFF', 'NONE' )
          LOG_INFO("ATMOS_PHY_CH_driver_tracer_setup",*) 'this component is never called.'
       case ( 'RN222' )

          call TRACER_regist( QS_CH,                        & ! [OUT]
                              ATMOS_PHY_CH_rn222_ntracers,  & ! [IN]
                              ATMOS_PHY_CH_rn222_NAME(:),   & ! [IN]
                              ATMOS_PHY_CH_rn222_DESC(:),   & ! [IN]
                              ATMOS_PHY_CH_rn222_UNIT(:)    ) ! [IN]
          QA_CH = ATMOS_PHY_CH_rn222_ntracers

       case default
          LOG_ERROR("ATMOS_PHY_CH_driver_tracer_setup",*) 'invalid chemistry type(', ATMOS_PHY_CH_TYPE, '). CHECK!'
          call PRC_abort
       end select

       QE_CH = QS_CH + QA_CH - 1

    else
       QA_CH = 0
       QS_CH = -1
       QE_CH = -1
    endif

    return
  end subroutine ATMOS_PHY_CH_driver_tracer_setup

  !-----------------------------------------------------------------------------
  !> Setup
  subroutine ATMOS_PHY_CH_driver_setup
    use mod_atmos_admin, only: &
       ATMOS_PHY_CH_TYPE, &
       ATMOS_sw_phy_ch
    use scale_atmos_phy_ch_rn222, only: &
       ATMOS_PHY_CH_rn222_setup
    implicit none
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_PROGRESS(*) 'Module[DRIVER] / Categ[ATMOS PHY_CH] / Origin[SCALE-RM]'

    if ( ATMOS_sw_phy_ch ) then

       call ATMOS_PHY_CH_rn222_setup

    else
       LOG_INFO("ATMOS_PHY_CH_driver_setup",*) 'this component is never called.'
    endif

    return
  end subroutine ATMOS_PHY_CH_driver_setup

  !-----------------------------------------------------------------------------
  !> Driver
  subroutine ATMOS_PHY_CH_driver_calc_tendency( update_flag )
    use scale_tracer, only: &
       TRACER_NAME
    use scale_time, only: &
       dt_CH => TIME_DTSEC_ATMOS_PHY_CH
    use scale_statistics, only: &
       STATISTICS_checktotal, &
       STATISTICS_total
    use scale_atmos_grid_cartesC_real, only: &
       ATMOS_GRID_CARTESC_REAL_VOL, &
       ATMOS_GRID_CARTESC_REAL_TOTVOL
    use scale_file_history, only: &
       FILE_HISTORY_in
    use mod_atmos_vars, only: &
       DENS   => DENS_av, &
       QTRC   => QTRC_av, &
       RHOQ_t => RHOQ_tp
    use mod_atmos_phy_ch_vars, only: &
       RHOQ_t_CH => ATMOS_PHY_CH_RHOQ_t, &
       QA_CH,        &
       QS_CH,        &
       QE_CH
!       O3        => ATMOS_PHY_CH_O3
    use scale_landuse, only: &
       LANDUSE_fact_land
    use scale_atmos_grid_cartesC_real, only: &
       ATMOS_GRID_CARTESC_REAL_FZ
    use mod_atmos_admin, only: &
       ATMOS_PHY_CH_TYPE
    use scale_atmos_grid_cartesC_index
    use scale_atmos_phy_ch_rn222, only: &
       ATMOS_PHY_CH_rn222_tendency
    use scale_prc, only: &
       PRC_abort
    implicit none

    logical, intent(in) :: update_flag

    integer  :: k, i, j, iq
    !---------------------------------------------------------------------------

    if ( update_flag ) then

!OCL XFILL
       RHOQ_t_CH(:,:,:,:) = 0.0_RP

       select case ( ATMOS_PHY_CH_TYPE )
       case ( 'RN222' )
          call ATMOS_PHY_CH_RN222_TENDENCY( KA, KS, KE, IA, IS, IE, JA, JS, JE, QA_CH, &
                                            DENS                      (:,:,:),             & ! [IN]
                                            QTRC                      (:,:,:,QS_CH:QE_CH), & ! [IN]
                                            ATMOS_GRID_CARTESC_REAL_FZ(:,:,:),             & ! [IN]
                                            LANDUSE_fact_land         (:,:),               & ! [IN]
                                            RHOQ_t_CH                 (:,:,:,QS_CH:QE_CH)  ) ! [INOUT]
       end select

       do iq = QS_CH, QE_CH
          call FILE_HISTORY_in( RHOQ_t_CH(:,:,:,iq), trim(TRACER_NAME(iq))//'_t_CH', &
                        'tendency rho*'//trim(TRACER_NAME(iq))//' in CH', 'kg/m3/s', fill_halo=.true. )
       enddo
    endif

    do iq = QS_CH, QE_CH
       !$omp parallel do private(i,j,k) OMP_SCHEDULE_
       do j  = JS, JE
       do i  = IS, IE
       do k  = KS, KE
          RHOQ_t(k,i,j,iq) = RHOQ_t(k,i,j,iq) + RHOQ_t_CH(k,i,j,iq)
       enddo
       enddo
       enddo
    enddo

    if ( STATISTICS_checktotal ) then
       do iq = QS_CH, QE_CH
          call STATISTICS_total( KA, KS, KE, IA, IS, IE, JA, JS, JE, &
                                 RHOQ_t_CH(:,:,:,iq), trim(TRACER_NAME(iq))//'_t_CH', &
                                 ATMOS_GRID_CARTESC_REAL_VOL(:,:,:),                  &
                                 ATMOS_GRID_CARTESC_REAL_TOTVOL                       )
       enddo
    endif

    return
  end subroutine ATMOS_PHY_CH_driver_calc_tendency

end module mod_atmos_phy_ch_driver
