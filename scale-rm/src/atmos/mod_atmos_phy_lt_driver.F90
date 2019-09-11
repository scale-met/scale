!-------------------------------------------------------------------------------
!> module ATMOSPHERE / Physics Chemistry
!!
!! @par Description
!!          Atmospheric lightning driver
!!
!! @author Team SCALE
!!
!<
!-------------------------------------------------------------------------------
#include "scalelib.h"
module mod_atmos_phy_lt_driver
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
  public :: ATMOS_PHY_LT_driver_tracer_setup
  public :: ATMOS_PHY_LT_driver_setup
  public :: ATMOS_PHY_LT_driver_calc_tendency

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
  subroutine ATMOS_PHY_LT_driver_tracer_setup
    use mod_atmos_admin, only: &
       ATMOS_PHY_MP_TYPE, &
       ATMOS_PHY_LT_TYPE, &
       ATMOS_sw_phy_lt
    use scale_tracer, only: &
       TRACER_regist
    use scale_atmos_phy_mp_suzuki10, only: &
       ATMOS_PHY_MP_suzuki10_nwaters, &
       ATMOS_PHY_MP_suzuki10_nices, &
       ATMOS_PHY_MP_suzuki10_nccn
    use mod_atmos_phy_lt_vars, only: &
       QA_LT, &
       QS_LT, &
       QE_LT
    use scale_atmos_hydrometeor, only: &
       QHA, &
       QHS
    use scale_prc, only: &
       PRC_abort
    implicit none

    character(len=H_SHORT), allocatable :: NAME(:)
    character(len=H_MID  ), allocatable :: DESC(:)
    character(len=H_SHORT), allocatable :: UNIT(:)
    integer :: iq
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("ATMOS_PHY_LT_driver_tracer_setup",*) 'Setup'

    if ( ATMOS_sw_phy_lt ) then

       select case ( ATMOS_PHY_LT_TYPE )
       case ( 'OFF', 'NONE' )
          LOG_INFO("ATMOS_PHY_LT_driver_tracer_setup",*) 'this component is never called.'
       case ( 'SATO2019' )
          select case ( ATMOS_PHY_MP_TYPE )
          case ( 'TOMITA08', 'SN14' )
             ! do nothing
          case ( 'SUZUKI10' )
             if( ATMOS_PHY_MP_suzuki10_nccn /= 0 ) then
                LOG_ERROR("ATMOS_PHY_LT_driver_tracer_setup",*) 'nccn in SUZUKI10 should be 0 for lithgning component(', ATMOS_PHY_MP_suzuki10_nccn, '). CHECK!'
                call PRC_abort
             endif
             if ( ATMOS_PHY_MP_suzuki10_nices == 0 ) then
                LOG_ERROR("ATMOS_PHY_LT_driver_tracer_setup",*) 'ICEFLG in SUZUKI10 should be 1 for lithgning component. CHECK!'
                call PRC_abort
             endif
          case ( 'KESSLER' )
             LOG_ERROR("ATMOS_PHY_LT_driver_tracer_setup",*) 'ATMOS_PHY_MP_TYPE should be TOMITA08, or SN14, or SUZUKI10 (', ATMOS_PHY_MP_TYPE, '). CHECK!'
             call PRC_abort
          end select
       case default
          LOG_ERROR("ATMOS_PHY_LT_driver_tracer_setup",*) 'invalid lithgning type(', ATMOS_PHY_LT_TYPE, '). CHECK!'
          call PRC_abort
       end select

       QA_LT = QHA

       allocate( NAME(QA_LT), DESC(QA_LT), UNIT(QA_LT) )
       do iq = 1, QA_LT
          NAME(iq) = 'QCRG_'//trim(TRACER_NAME(QHS+iq-1)(2:))
          DESC(iq) = 'Ratio of charge density of '//trim(TRACER_NAME(QHS+iq-1))
          UNIT(iq) = 'fC/kg'
       end do
       call TRACER_regist( QS_LT,   & ! [OUT]
                           QA_LT,   & ! [IN]
                           NAME(:), & ! [IN]
                           DESC(:), & ! [IN]
                           UNIT(:)  ) ! [IN]
       deallocate( NAME, DESC, UNIT )

       QE_LT = QS_LT - 1 + QA_LT

    else
       QA_LT = 0
       QS_LT = -1
       QE_LT = -2
    endif

    return
  end subroutine ATMOS_PHY_LT_driver_tracer_setup

  !-----------------------------------------------------------------------------
  !> Setup
  subroutine ATMOS_PHY_LT_driver_setup
    use scale_time, only: &
       dt_LT => TIME_DTSEC_ATMOS_PHY_LT, &
       dt_MP => TIME_DTSEC_ATMOS_PHY_MP
    use scale_atmos_grid_cartesC_real, only: &
       REAL_LON => ATMOS_GRID_CARTESC_REAL_LON, &
       REAL_LAT => ATMOS_GRID_CARTESC_REAL_LAT
    use scale_atmos_grid_cartesC, only: &
       CDX => ATMOS_GRID_CARTESC_CDX, &
       CDY => ATMOS_GRID_CARTESC_CDY
    use scale_atmos_phy_lt_sato2019, only: &
       ATMOS_PHY_LT_sato2019_setup
    use mod_atmos_phy_lt_vars, only: &
       flg_lt
    use mod_atmos_admin, only: &
       ATMOS_PHY_MP_TYPE, &
       ATMOS_PHY_LT_TYPE, &
       ATMOS_sw_phy_lt
    use scale_prc, only: &
       PRC_abort
    use scale_atmos_hydrometeor, only: &
       QLA, &
       QIA
    implicit none
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("ATMOS_PHY_LT_driver_setup",*) 'Setup'

    if ( ATMOS_sw_phy_lt ) then

       if( dt_LT /= dt_MP ) then
          LOG_ERROR("ATMOS_PHY_LT_driver_setup",*) 'dt_LT should be dt_MP. CHECK!'
          call PRC_abort
       endif

       flg_lt = .true.
       call ATMOS_PHY_LT_sato2019_setup( KA, KS, KE, & ! [IN]
                                         IA, IS, IE, & ! [IN]
                                         JA, JS, JE, & ! [IN]
                                         IMAXG,      & ! [IN]
                                         JMAXG,      & ! [IN]
                                         KMAX,       & ! [IN]
                                         ATMOS_PHY_MP_TYPE, & ! [IN]
                                         QLA, QIA,          & ! [IN]
                                         CDX, CDY           ) ! [IN]

    else

       flg_lt = .false.
       LOG_INFO("ATMOS_PHY_LT_driver_setup",*) 'This component is never called.'

    endif

    return
  end subroutine ATMOS_PHY_LT_driver_setup

  !-----------------------------------------------------------------------------
  !> Driver
  subroutine ATMOS_PHY_LT_driver_calc_tendency( update_flag )
    use scale_tracer, only: &
       TRACER_NAME
    use scale_time, only: &
       dt_LT => TIME_DTSEC_ATMOS_PHY_LT
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
    use mod_atmos_phy_lt_vars, only: &
       RHOQ_t_LT => ATMOS_PHY_LT_RHOQ_t, &
       RHOQ_t_LT_mp => ATMOS_PHY_LT_RHOQ_mp_t, &
       QA_LT,        &
       QS_LT,        &
       QE_LT
    implicit none

    logical, intent(in) :: update_flag

    integer  :: k, i, j, iq
    !---------------------------------------------------------------------------

    if ( update_flag ) then

       do iq = QS_LT, QE_LT
          call FILE_HISTORY_in( RHOQ_t_LT(:,:,:,iq), trim(TRACER_NAME(iq))//'_t_LT', &
                                'tendency rho*'//trim(TRACER_NAME(iq))//' in LT',    &
                                'kg/m3/s', fill_halo=.true.                          )
       enddo
    endif

    do iq = QS_LT, QE_LT
       !$omp parallel do private(i,j,k) OMP_SCHEDULE_
       do j  = JS, JE
       do i  = IS, IE
       do k  = KS, KE
          RHOQ_t(k,i,j,iq) = RHOQ_t(k,i,j,iq) &
                           + RHOQ_t_LT(k,i,j,iq) &
                           + RHOQ_t_LT_mp(k,i,j,iq)
       enddo
       enddo
       enddo
    enddo

    if ( STATISTICS_checktotal ) then
       do iq = QS_LT, QE_LT
          call STATISTICS_total( KA, KS, KE, IA, IS, IE, JA, JS, JE, &
                                 RHOQ_t_LT(:,:,:,iq), trim(TRACER_NAME(iq))//'_t_LT', &
                                 ATMOS_GRID_CARTESC_REAL_VOL(:,:,:),                  &
                                 ATMOS_GRID_CARTESC_REAL_TOTVOL                       )
       enddo
    endif

    return
  end subroutine ATMOS_PHY_LT_driver_calc_tendency

end module mod_atmos_phy_lt_driver
