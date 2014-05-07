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
  use scale_grid_index
  use scale_tracer
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: ATMOS_PHY_CH_driver_setup
  public :: ATMOS_PHY_CH_driver

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
  subroutine ATMOS_PHY_CH_driver_setup( CH_TYPE )
!    use scale_atmos_phy_ch, only: &
!       ATMOS_PHY_CH_setup
    implicit none

    character(len=H_SHORT), intent(in) :: CH_TYPE
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[Physics-CH]/Categ[ATMOS]'

!    call ATMOS_PHY_CH_setup( CH_TYPE )

    call ATMOS_PHY_CH_driver( .true., .false. )

    return
  end subroutine ATMOS_PHY_CH_driver_setup

  !-----------------------------------------------------------------------------
  !> Driver
  subroutine ATMOS_PHY_CH_driver( update_flag, history_flag )
    use scale_time, only: &
       dt_CH => TIME_DTSEC_ATMOS_PHY_CH
    use scale_history, only: &
       HIST_in
!    use scale_atmos_phy_ch, only: &
!       ATMOS_PHY_CH
    use mod_atmos_vars, only: &
       DENS,              &
       MOMZ,              &
       MOMX,              &
       MOMY,              &
       RHOT,              &
       QTRC,              &
       QTRC_t => QTRC_tp
    use mod_atmos_phy_ch_vars, only: &
       QTRC_t_CH => ATMOS_PHY_CH_QTRC_t, &
       O3        => ATMOS_PHY_CH_O3
    implicit none

    logical, intent(in) :: update_flag
    logical, intent(in) :: history_flag

    integer :: k, i, j, iq
    !---------------------------------------------------------------------------

    if ( update_flag ) then
!       call ATMOS_PHY_CH( DENS,          & ! [IN]
!                          MOMZ,          & ! [IN]
!                          MOMX,          & ! [IN]
!                          MOMY,          & ! [IN]
!                          RHOT,          & ! [IN]
!                          QTRC,          & ! [IN]
!                          QTRC_t_CH,     & ! [INOUT]
!                          O3             ) ! [INOUT]

       ! tentative!
       QTRC_t_CH(:,:,:,:) = 0.0_RP
       O3       (:,:,:)   = 0.0_RP

       if ( history_flag ) then
          call HIST_in( O3(:,:,:), 'Ozone', 'Ozone', 'PPM', dt_CH )
       endif
    endif

    !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(3)
    do iq = 1, QA
    do j  = JS, JE
    do i  = IS, IE
    do k  = KS, KE
       QTRC_t(k,i,j,iq) = QTRC_t(k,i,j,iq) + QTRC_t_CH(k,i,j,iq)
    enddo
    enddo
    enddo
    enddo

    return
  end subroutine ATMOS_PHY_CH_driver

end module mod_atmos_phy_ch_driver
