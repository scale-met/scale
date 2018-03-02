!-------------------------------------------------------------------------------
!> module ATMOSPHERE / Physics Aerosol Microphysics
!!
!! @par Description
!!          Aerosol Microphysics driver
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2013-12-06 (S.Nishizawa)  [new]
!<
!-------------------------------------------------------------------------------
#include "inc_openmp.h"
module mod_atmos_phy_ae_driver
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
  public :: ATMOS_PHY_AE_driver_tracer_setup
  public :: ATMOS_PHY_AE_driver_setup
  public :: ATMOS_PHY_AE_driver_resume
  public :: ATMOS_PHY_AE_driver_adjustment
  public :: ATMOS_PHY_AE_driver_calc_tendency

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
  !> Setup
  subroutine ATMOS_PHY_AE_driver_tracer_setup
    use mod_atmos_admin, only: &
       ATMOS_PHY_AE_TYPE, &
       ATMOS_sw_phy_ae
    use scale_atmos_phy_ae_kajino13, only: &
       ATMOS_PHY_AE_kajino13_tracer_setup, &
       ATMOS_PHY_AE_kajino13_NAME, &
       ATMOS_PHY_AE_kajino13_DESC, &
       ATMOS_PHY_AE_kajino13_UNIT
    use scale_process, only: &
       PRC_abort
    use mod_atmos_phy_ae_vars, only: &
       QA_AE, &
       QS_AE, &
       QE_AE
    implicit none

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '++++++ Module[Tracer Setup] / Categ[ATMOS PHY_AE] / Origin[SCALE-RM]'

    if ( ATMOS_sw_phy_ae ) then
       select case ( ATMOS_PHY_AE_TYPE )
       case ( 'OFF', 'NONE' )
          if( IO_L ) write(IO_FID_LOG,*) '*** this component is never called.'
       case ( 'KAJINO13' )
          call ATMOS_PHY_AE_kajino13_tracer_setup( QA_AE ) ! [OUT]

          call TRACER_regist( QS_AE,                         & ! [OUT]
                              QA_AE,                         & ! [IN]
                              ATMOS_PHY_AE_kajino13_NAME(:), & ! [IN]
                              ATMOS_PHY_AE_kajino13_DESC(:), & ! [IN]
                              ATMOS_PHY_AE_kajino13_UNIT(:)  ) ! [IN]
       case default
          write(*,*) 'xxx invalid aerosol type(', ATMOS_PHY_AE_TYPE, '). CHECK!'
          call PRC_abort
       end select

       QE_AE = QS_AE + QA_AE - 1

    else
       QA_AE = 0
       QS_AE = -1
       QE_AE = -1
    end if

    return
  end subroutine ATMOS_PHY_AE_driver_tracer_setup

  !-----------------------------------------------------------------------------
  !> Setup
  subroutine ATMOS_PHY_AE_driver_setup
    use scale_atmos_phy_ae, only: &
       ATMOS_PHY_AE_setup
    use mod_atmos_admin, only: &
       ATMOS_PHY_AE_TYPE, &
       ATMOS_sw_phy_ae
    use scale_atmos_phy_ae_kajino13, only: &
        ATMOS_PHY_AE_kajino13_setup
    use scale_process, only: &
       PRC_abort
    implicit none
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '++++++ Module[DRIVER] / Categ[ATMOS PHY_AE] / Origin[SCALE-RM]'

    if ( ATMOS_sw_phy_ae ) then

       select case ( ATMOS_PHY_AE_TYPE )
       case ( 'KAJINO13' )
          call ATMOS_PHY_AE_kajino13_setup
       case default
          write(*,*) 'xxx invalid aerosol type(', ATMOS_PHY_AE_TYPE, '). CHECK!'
          call PRC_abort
       end select

    endif

    return
  end subroutine ATMOS_PHY_AE_driver_setup


  !-----------------------------------------------------------------------------
  !> Resume
  subroutine ATMOS_PHY_AE_driver_resume
     use mod_atmos_admin, only: &
        ATMOS_sw_phy_ae
    implicit none

    if ( ATMOS_sw_phy_ae ) then

       ! run once (only for the diagnostic value)
       call PROF_rapstart('ATM_Aerosol', 1)
       call ATMOS_PHY_AE_driver_calc_tendency( update_flag = .true. )
       call PROF_rapend  ('ATM_Aerosol', 1)

    endif

    return
  end subroutine ATMOS_PHY_AE_driver_resume

  !-----------------------------------------------------------------------------
  !> adjustment
  subroutine ATMOS_PHY_AE_driver_adjustment
    use mod_atmos_vars, only: &
       QTRC
    use mod_atmos_phy_ae_vars, only: &
       QA_AE, &
       QS_AE, &
       QE_AE
    use mod_atmos_admin, only: &
       ATMOS_PHY_AE_TYPE, &
       ATMOS_sw_phy_ae
    use scale_atmos_phy_ae_kajino13, only: &
        ATMOS_PHY_AE_kajino13_negative_fixer
    implicit none

    if ( ATMOS_sw_phy_ae ) then
       
       select case ( ATMOS_PHY_AE_TYPE )
       case ( 'KAJINO13' )
          call ATMOS_PHY_AE_kajino13_negative_fixer( KA, KS, KE, IA, ISB, IEB, JA, JSB, JEB, QA_AE, &
                                                     QTRC(:,:,:,QS_AE:QE_AE) ) ! [INOUT]

       end select

    end if

    return
  end subroutine ATMOS_PHY_AE_driver_adjustment

  !-----------------------------------------------------------------------------
  !> Driver
  subroutine ATMOS_PHY_AE_driver_calc_tendency( update_flag )
    use scale_process, only: &
       PRC_abort
    use scale_time, only: &
       dt_AE => TIME_DTSEC_ATMOS_PHY_AE
    use scale_rm_statistics, only: &
       STATISTICS_checktotal, &
       STAT_total
    use scale_file_history, only: &
       FILE_HISTORY_in
    use scale_atmos_phy_ae, only: &
       ATMOS_PHY_AE
    use mod_atmos_vars, only: &
       DENS => DENS_av, &
       QTRC => QTRC_av, &
       QDRY, &
       PRES, &
       TEMP, &
       QV,   &
       RHOQ_t => RHOQ_tp
    use mod_atmos_phy_ae_vars, only: &
       QA_AE, &
       QS_AE, &
       QE_AE, &
       RHOQ_t_AE => ATMOS_PHY_AE_RHOQ_t, &
       CCN       => ATMOS_PHY_AE_CCN,   &
       CCN_t     => ATMOS_PHY_AE_CCN_t, &
       AE_EMIT   => ATMOS_PHY_AE_EMIT
    use mod_atmos_phy_mp_vars, only: &
       EVAPORATE => ATMOS_PHY_MP_EVAPORATE
    use mod_atmos_admin, only: &
       ATMOS_PHY_AE_TYPE
    use scale_atmos_phy_ae_kajino13, only: &
       ATMOS_PHY_AE_kajino13_tendency
    implicit none

    logical, intent(in) :: update_flag

    real(RP) :: CN(KA,IA,JA)
    real(RP) :: NREG(KA,IA,JA)

    real(RP) :: total ! dummy

    integer  :: k, i, j, iq
    !---------------------------------------------------------------------------

    if ( update_flag ) then

!OCL XFILL
       CCN      (:,:,:)   = 0.0_RP ! reset
!OCL XFILL
       RHOQ_t_AE(:,:,:,:) = 0.0_RP ! reset

       do j  = JS, JE
       do i  = IS, IE
       do k  = KS, KE
          NREG(k,i,j) = EVAPORATE(k,i,j) * dt_AE
       enddo
       enddo
       enddo

       select case ( ATMOS_PHY_AE_TYPE )
       case ( 'KAJINO13' )
          call ATMOS_PHY_AE_kajino13_tendency( KA, KS, KE, IA, IS, IE, JA, JS, JE, QA_AE, &
                                               TEMP     (:,:,:),             & ! [IN]
                                               PRES     (:,:,:),             & ! [IN]
                                               QDRY     (:,:,:),             & ! [IN]
                                               NREG     (:,:,:),             & ! [IN]
                                               DENS     (:,:,:),             & ! [IN]
                                               QV       (:,:,:),             & ! [IN]
                                               QTRC     (:,:,:,QS_AE:QE_AE), & ! [IN]
                                               AE_EMIT  (:,:,:,QS_AE:QE_AE), & ! [IN]
                                               RHOQ_t_AE(:,:,:,QS_AE:QE_AE), & ! [OUT]
                                               CN       (:,:,:),             & ! [OUT]
                                               CCN      (:,:,:)              ) ! [OUT]
       end select

       CCN_t(:,:,:) = CCN(:,:,:) / dt_AE

       call FILE_HISTORY_in( CN (:,:,:)*1.E-6_RP, 'CN',  'condensation nucrei',       'num/cc' )
       call FILE_HISTORY_in( CCN(:,:,:)*1.E-6_RP, 'CCN', 'cloud condensation nucrei', 'num/cc' )

    endif

    do iq = QS_AE, QE_AE
       !$omp parallel do private(i,j,k) OMP_SCHEDULE_
       do j  = JS, JE
       do i  = IS, IE
       do k  = KS, KE
          RHOQ_t(k,i,j,iq) = RHOQ_t(k,i,j,iq) + RHOQ_t_AE(k,i,j,iq)
       enddo
       enddo
       enddo
    enddo

    if ( STATISTICS_checktotal ) then
       do iq = QS_AE, QE_AE
          call STAT_total( total, RHOQ_t_AE(:,:,:,iq), trim(TRACER_NAME(iq))//'_t_AE' )
       enddo
    endif

    return
  end subroutine ATMOS_PHY_AE_driver_calc_tendency

end module mod_atmos_phy_ae_driver
