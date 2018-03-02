!-------------------------------------------------------------------------------
!> module ATMOSPHERE / Physics Chemistry
!!
!! @par Description
!!          Atmospheric chemistry driver
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2014-05-04 (H.Yashiro)    [new]
!<
!-------------------------------------------------------------------------------
#include "inc_openmp.h"
module mod_atmos_phy_ch_driver
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
  public :: ATMOS_PHY_CH_driver_tracer_setup
  public :: ATMOS_PHY_CH_driver_setup
  public :: ATMOS_PHY_CH_driver_resume
  public :: ATMOS_PHY_CH_driver_tendency

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
    use scale_atmos_phy_ch, only: &
       ATMOS_PHY_CH_setup
    use mod_atmos_admin, only: &
       ATMOS_PHY_CH_TYPE, &
       ATMOS_sw_phy_ch
    use scale_atmos_phy_ch_rn222, only: &
       QS_CH, &
       QE_CH, &
       QA_CH, &
       ATMOS_PHY_CH_rn222_NAME, &
       ATMOS_PHY_CH_rn222_DESC, &
       ATMOS_PHY_CH_rn222_UNIT
    use scale_atmos_phy_ch, only: &
       QA => QA_CH, &
       QS => QS_CH, &
       QE => QE_CH
    use scale_process, only: &
       PRC_MPIstop
    implicit none
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '++++++ Module[Tracer Setup] / Categ[ATMOS PHY_CH] / Origin[SCALE-RM]'

    if ( ATMOS_sw_phy_ch ) then

       select case ( ATMOS_PHY_CH_TYPE )
       case ( 'RN222' )
          call TRACER_regist( QS_CH,                   & ! [OUT]
                              QA_CH,                   & ! [IN]
                              ATMOS_PHY_CH_rn222_NAME, & ! [IN]
                              ATMOS_PHY_CH_rn222_DESC, & ! [IN]
                              ATMOS_PHY_CH_rn222_UNIT  ) ! [IN]

          QA    = QA_CH
          QS    = QS_CH
          QE_CH = QS_CH + QA_CH - 1
          QE    = QE_CH
       case default
          write(*,*) 'xxx invalid chemistry type(', ATMOS_PHY_CH_TYPE, '). CHECK!'
          call PRC_MPIstop
       end select

    else
       if( IO_L ) write(IO_FID_LOG,*) '*** this component is never called.'
    endif

    return
  end subroutine ATMOS_PHY_CH_driver_tracer_setup

  !-----------------------------------------------------------------------------
  !> Setup
  subroutine ATMOS_PHY_CH_driver_setup
    use scale_atmos_phy_ch, only: &
       ATMOS_PHY_CH_setup
    use mod_atmos_admin, only: &
       ATMOS_PHY_CH_TYPE, &
       ATMOS_sw_phy_ch
    use scale_atmos_phy_ch_rn222, only: &
       ATMOS_PHY_CH_rn222_setup
    implicit none
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '++++++ Module[DRIVER] / Categ[ATMOS PHY_CH] / Origin[SCALE-RM]'

    if ( ATMOS_sw_phy_ch ) then

       call ATMOS_PHY_CH_rn222_setup

    else
       if( IO_L ) write(IO_FID_LOG,*) '*** this component is never called.'
    endif

    return
  end subroutine ATMOS_PHY_CH_driver_setup

  !-----------------------------------------------------------------------------
  !> Resume
  subroutine ATMOS_PHY_CH_driver_resume
    use mod_atmos_admin, only: &
       ATMOS_sw_phy_ch
    implicit none

    if ( ATMOS_sw_phy_ch ) then

       ! run once (only for the diagnostic value)
       call PROF_rapstart('ATM_Chemistry', 1)
       call ATMOS_PHY_CH_driver_tendency( update_flag = .true. )
       call PROF_rapend  ('ATM_Chemistry', 1)

    end if

    return
  end subroutine ATMOS_PHY_CH_driver_resume

  !-----------------------------------------------------------------------------
  !> Driver
  subroutine ATMOS_PHY_CH_driver_tendency( update_flag )
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
    use scale_atmos_phy_ch, only: &
       ATMOS_PHY_CH, &
       QA_CH,        &
       QS_CH,        &
       QE_CH
    use mod_atmos_vars, only: &
       DENS   => DENS_av, &
       QTRC   => QTRC_av, &
       RHOQ_t => RHOQ_tp
    use mod_atmos_phy_ch_vars, only: &
       RHOQ_t_CH => ATMOS_PHY_CH_RHOQ_t
!       O3        => ATMOS_PHY_CH_O3
    use scale_landuse, only: &
       LANDUSE_fact_land
    use scale_atmos_grid_cartesC_real, only: &
       ATMOS_GRID_CARTESC_REAL_FZ
    use mod_atmos_admin, only: &
       ATMOS_PHY_CH_TYPE
    use scale_atmos_grid_cartesC_index
    use scale_atmos_phy_ch_rn222, only: &
       ATMOS_PHY_CH_rn222_tendency, &
       I_ch_rn222
    use scale_process, only: &
       PRC_MPIstop
    implicit none

    logical, intent(in) :: update_flag

    integer  :: k, i, j, iq
    !---------------------------------------------------------------------------

    iq = QS_CH - 1 + I_ch_rn222

    if ( update_flag ) then

!OCL XFILL
       RHOQ_t_CH(:,:,:,:) = 0.0_RP

       select case ( ATMOS_PHY_CH_TYPE )
       case ( 'RN222' )
          call ATMOS_PHY_CH_RN222_TENDENCY( KA, KS, KE,                                    & ! [IN]
                                            IA, IS, IE,                                    & ! [IN]
                                            JA, JS, JE,                                    & ! [IN]
                                            QA_CH,                                         & ! [IN]
                                            DENS                      (:,:,:),             & ! [IN]
                                            QTRC                      (:,:,:,iq), & ! [IN]
                                            ATMOS_GRID_CARTESC_REAL_FZ(:,:,:),             & ! [IN]
                                            LANDUSE_fact_land         (:,:),               & ! [IN]
                                            RHOQ_t_CH                 (:,:,:,:)            ) ! [INOUT]
       case default
          write(*,*) 'xxx invalid chemistry type(', ATMOS_PHY_CH_TYPE, '). CHECK!'
          call PRC_MPIstop
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
  end subroutine ATMOS_PHY_CH_driver_tendency

end module mod_atmos_phy_ch_driver
