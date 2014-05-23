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
  use scale_grid_index
  use scale_tracer
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: ATMOS_PHY_AE_driver_setup
  public :: ATMOS_PHY_AE_driver

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
  subroutine ATMOS_PHY_AE_driver_setup
    use scale_atmos_phy_ae, only: &
       ATMOS_PHY_AE_setup
    use mod_atmos_admin, only: &
       ATMOS_PHY_AE_TYPE, &
       ATMOS_sw_phy_ae
    implicit none
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '++++++ Module[DRIVER] / Categ[ATMOS PHY_AE] / Origin[SCALE-LES]'

    ! note: tentatively, aerosol module should be called at all time. we need dummy subprogram.
!    if ( ATMOS_sw_phy_ae ) then

       ! setup library component
       call ATMOS_PHY_AE_setup( ATMOS_PHY_AE_TYPE )

       ! run once (only for the diagnostic value)
       call ATMOS_PHY_AE_driver( .true., .false. )

!    endif

    return
  end subroutine ATMOS_PHY_AE_driver_setup

  !-----------------------------------------------------------------------------
  !> Driver
  subroutine ATMOS_PHY_AE_driver( update_flag, history_flag )
    use scale_time, only: &
       dt_AE => TIME_DTSEC_ATMOS_PHY_AE
    use scale_stats, only: &
       STAT_checktotal, &
       STAT_total
    use scale_history, only: &
       HIST_in
    use scale_atmos_phy_ae, only: &
       ATMOS_PHY_AE
    use mod_atmos_vars, only: &
       DENS,              &
       MOMZ,              &
       MOMX,              &
       MOMY,              &
       RHOT,              &
       QTRC,              &
       QTRC_t => QTRC_tp
    use mod_atmos_phy_ae_vars, only: &
       QTRC_t_AE => ATMOS_PHY_AE_QTRC_t, &
       CCN       => ATMOS_PHY_AE_CCN
    implicit none

    logical, intent(in) :: update_flag
    logical, intent(in) :: history_flag

    real(RP) :: QTRC0(KA,IA,JA,QA)

    real(RP) :: RHOQ(KA,IA,JA)
    real(RP) :: total ! dummy

    integer  :: k, i, j, iq
    !---------------------------------------------------------------------------

    if ( update_flag ) then

       if( IO_L ) write(IO_FID_LOG,*) '*** Physics step, aerosol'

       do iq = 1, QA
       do j  = JS, JE
       do i  = IS, IE
       do k  = KS, KE
          QTRC0(k,i,j,iq) = QTRC(k,i,j,iq) ! save
       enddo
       enddo
       enddo
       enddo

       call ATMOS_PHY_AE( DENS, & ! [IN]
                          MOMZ, & ! [IN]
                          MOMX, & ! [IN]
                          MOMY, & ! [IN]
                          RHOT, & ! [IN]
                          QTRC0 ) ! [INOUT]

       do iq = 1, QA
       do j  = JS, JE
       do i  = IS, IE
       do k  = KS, KE
          QTRC_t_AE(k,i,j,iq) = QTRC0(k,i,j,iq) - QTRC(k,i,j,iq)
       enddo
       enddo
       enddo
       enddo

       CCN(:,:,:) = 0.0_RP ! tentative

       if ( history_flag ) then
          call HIST_in( CCN(:,:,:), 'CCN', 'cloud condensation nucrei', '', dt_AE )
       endif

    endif

    !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(3)
    do iq = 1, QA
    do j  = JS, JE
    do i  = IS, IE
    do k  = KS, KE
       QTRC_t(k,i,j,iq) = QTRC_t(k,i,j,iq) + QTRC_t_AE(k,i,j,iq)
    enddo
    enddo
    enddo
    enddo

    if ( STAT_checktotal ) then
       do iq = 1, QA
          RHOQ(:,:,:) = DENS(:,:,:) * QTRC_t_AE(:,:,:,iq)

          call STAT_total( total, RHOQ(:,:,:), trim(AQ_NAME(iq))//'_t_AE' )
       enddo
    endif

    return
  end subroutine ATMOS_PHY_AE_driver

end module mod_atmos_phy_ae_driver
