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
  public :: ATMOS_PHY_CH_driver_tracer_setup
  public :: ATMOS_PHY_CH_driver_setup
  public :: ATMOS_PHY_CH_driver_calc_tendency
  public :: ATMOS_PHY_CH_driver_OCEAN_flux
  public :: ATMOS_PHY_CH_driver_LAND_flux
  public :: ATMOS_PHY_CH_driver_URBAN_flux

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
    LOG_INFO("ATMOS_PHY_CH_driver_tracer_setup",*) 'Setup'

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
       QE_CH = -2
    endif

    return
  end subroutine ATMOS_PHY_CH_driver_tracer_setup

  !-----------------------------------------------------------------------------
  !> Setup
  subroutine ATMOS_PHY_CH_driver_setup
    use scale_atmos_grid_cartesC_real, only: &
       REAL_LON => ATMOS_GRID_CARTESC_REAL_LON, &
       REAL_LAT => ATMOS_GRID_CARTESC_REAL_LAT
    use scale_atmos_phy_ch_rn222, only: &
       ATMOS_PHY_CH_rn222_setup
    use scale_atmos_sfc_ch_rn222, only: &
       ATMOS_SFC_CH_rn222_setup
    use mod_atmos_admin, only: &
       ATMOS_PHY_CH_TYPE, &
       ATMOS_sw_phy_ch
    implicit none
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("ATMOS_PHY_CH_driver_setup",*) 'Setup'

    if ( ATMOS_sw_phy_ch ) then

       call ATMOS_PHY_CH_rn222_setup
       call ATMOS_SFC_CH_rn222_setup( IA, JA,            & ! [IN]
                                      REAL_LON, REAL_LAT ) ! [IN]

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
    use mod_atmos_phy_sf_vars, only: &
       SFLX_QTRC => ATMOS_PHY_SF_SFLX_QTRC
    use mod_atmos_phy_ch_vars, only: &
       RHOQ_t_CH => ATMOS_PHY_CH_RHOQ_t, &
       QA_CH,        &
       QS_CH,        &
       QE_CH
!       O3        => ATMOS_PHY_CH_O3
    use mod_atmos_admin, only: &
       ATMOS_PHY_CH_TYPE
    use scale_atmos_phy_ch_rn222, only: &
       ATMOS_PHY_CH_rn222_tendency
    implicit none

    logical, intent(in) :: update_flag

    integer  :: k, i, j, iq
    !---------------------------------------------------------------------------

    if ( update_flag ) then

!OCL XFILL
       RHOQ_t_CH(:,:,:,:) = 0.0_RP

       select case( ATMOS_PHY_CH_TYPE )
       case( 'RN222' )
          call ATMOS_PHY_CH_rn222_TENDENCY( KA, KS, KE,                   & ! [IN]
                                            IA, IS, IE,                   & ! [IN]
                                            JA, JS, JE,                   & ! [IN]
                                            QA_CH,                        & ! [IN]
                                            DENS     (:,:,:),             & ! [IN]
                                            QTRC     (:,:,:,QS_CH:QE_CH), & ! [IN]
                                            RHOQ_t_CH(:,:,:,QS_CH:QE_CH)  ) ! [INOUT]
       end select

       do iq = QS_CH, QE_CH
          call FILE_HISTORY_in( RHOQ_t_CH(:,:,:,iq), trim(TRACER_NAME(iq))//'_t_CH', &
                                'tendency rho*'//trim(TRACER_NAME(iq))//' in CH',    &
                                'kg/m3/s', fill_halo=.true.                          )
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

  !-----------------------------------------------------------------------------
  !> Driver
  subroutine ATMOS_PHY_CH_driver_OCEAN_flux( &
       SFLX_QTRC )
    use scale_ocean_grid_cartesC_index
    use scale_atmos_sfc_ch_rn222, only: &
       ATMOS_SFC_CH_rn222_OCEAN_flux
    use mod_atmos_admin, only: &
       ATMOS_PHY_CH_TYPE
    use mod_atmos_phy_ch_vars, only: &
       QA_CH, &
       QS_CH, &
       QE_CH
    implicit none

    real(RP), intent(inout) :: SFLX_QTRC(OIA,OJA,QA)

    !---------------------------------------------------------------------------

    select case( ATMOS_PHY_CH_TYPE )
    case( 'RN222' )
       call ATMOS_SFC_CH_rn222_OCEAN_flux( OIA, OIS, OIE,             & ! [IN]
                                           OJA, OJS, OJE,             & ! [IN]
                                           QA_CH,                     & ! [IN]
                                           SFLX_QTRC(:,:,QS_CH:QE_CH) ) ! [INOUT]
    end select

    return
  end subroutine ATMOS_PHY_CH_driver_OCEAN_flux

  !-----------------------------------------------------------------------------
  !> Driver
  subroutine ATMOS_PHY_CH_driver_LAND_flux( &
       SFLX_QTRC   )
    use scale_land_grid_cartesC_index
    use scale_time, only: &
       TIME_NOWDATE
    use scale_atmos_sfc_ch_rn222, only: &
       ATMOS_SFC_CH_rn222_LAND_flux
    use mod_atmos_admin, only: &
       ATMOS_PHY_CH_TYPE
    use mod_atmos_phy_ch_vars, only: &
       QA_CH, &
       QS_CH, &
       QE_CH
    implicit none

    real(RP), intent(inout) :: SFLX_QTRC(LIA,LJA,QA)

    !---------------------------------------------------------------------------

    select case( ATMOS_PHY_CH_TYPE )
    case( 'RN222' )
       call ATMOS_SFC_CH_rn222_LAND_flux( LIA, LIS, LIE,             & ! [IN]
                                          LJA, LJS, LJE,             & ! [IN]
                                          QA_CH,                     & ! [IN]
                                          TIME_NOWDATE(:),           & ! [IN]
                                          SFLX_QTRC(:,:,QS_CH:QE_CH) ) ! [INOUT]
    end select

    return
  end subroutine ATMOS_PHY_CH_driver_LAND_flux

  !-----------------------------------------------------------------------------
  !> Driver
  subroutine ATMOS_PHY_CH_driver_URBAN_flux( &
       SFLX_QTRC )
    use scale_urban_grid_cartesC_index
    use scale_time, only: &
       TIME_NOWDATE
    use scale_atmos_sfc_ch_rn222, only: &
       ATMOS_SFC_CH_rn222_LAND_flux
    use mod_atmos_admin, only: &
       ATMOS_PHY_CH_TYPE
    use mod_atmos_phy_ch_vars, only: &
       QA_CH, &
       QS_CH, &
       QE_CH
    implicit none

    real(RP), intent(inout) :: SFLX_QTRC(UIA,UJA,QA)

    !---------------------------------------------------------------------------

    select case( ATMOS_PHY_CH_TYPE )
    case( 'RN222' )
       call ATMOS_SFC_CH_rn222_LAND_flux( UIA, UIS, UIE,             & ! [IN]
                                          UJA, UJS, UJE,             & ! [IN]
                                          QA_CH,                     & ! [IN]
                                          TIME_NOWDATE(:),           & ! [IN]
                                          SFLX_QTRC(:,:,QS_CH:QE_CH) ) ! [INOUT]
    end select

    return
  end subroutine ATMOS_PHY_CH_driver_URBAN_flux

end module mod_atmos_phy_ch_driver
