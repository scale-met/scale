!-------------------------------------------------------------------------------
!> module ATMOSPHERE / Physics Cloud Microphysics
!!
!! @par Description
!!          Cloud Microphysics driver
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2013-12-06 (S.Nishizawa)  [new]
!<
!-------------------------------------------------------------------------------
#include "inc_openmp.h"
module mod_atmos_phy_mp_driver
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
  public :: ATMOS_PHY_MP_driver_setup
  public :: ATMOS_PHY_MP_driver

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
  subroutine ATMOS_PHY_MP_driver_setup
    use scale_atmos_phy_mp, only: &
       ATMOS_PHY_MP_setup
    use mod_atmos_admin, only: &
       ATMOS_PHY_MP_TYPE, &
       ATMOS_sw_phy_mp
    use mod_atmos_phy_mp_vars, only: &
       SFLX_rain => ATMOS_PHY_MP_SFLX_rain, &
       SFLX_snow => ATMOS_PHY_MP_SFLX_snow
    implicit none
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '++++++ Module[DRIVER] / Categ[ATMOS PHY_MP] / Origin[SCALE-LES]'

    if ( ATMOS_sw_phy_mp ) then

       ! setup library component
       call ATMOS_PHY_MP_setup( ATMOS_PHY_MP_TYPE )

       ! run once (only for the diagnostic value)
       call ATMOS_PHY_MP_driver( .true., .false. )

    else

       if( IO_L ) write(IO_FID_LOG,*) '*** ATMOS_PHY_MP is disabled.'
       if( IO_L ) write(IO_FID_LOG,*) '*** SFLX_rain and SFLX_snow is set to zero.'
       SFLX_rain(:,:) = 0.0_RP
       SFLX_snow(:,:) = 0.0_RP

    endif

    return
  end subroutine ATMOS_PHY_MP_driver_setup

  !-----------------------------------------------------------------------------
  !> Driver
  subroutine ATMOS_PHY_MP_driver( update_flag, history_flag )
    use scale_time, only: &
       dt_MP => TIME_DTSEC_ATMOS_PHY_MP
    use scale_stats, only: &
       STAT_checktotal, &
       STAT_total
    use scale_history, only: &
       HIST_in
    use scale_atmos_phy_mp, only: &
       ATMOS_PHY_MP
    use mod_atmos_vars, only: &
       DENS,              &
       MOMZ,              &
       MOMX,              &
       MOMY,              &
       RHOT,              &
       QTRC,              &
       DENS_t => DENS_tp, &
       MOMZ_t => MOMZ_tp, &
       MOMX_t => MOMX_tp, &
       MOMY_t => MOMY_tp, &
       RHOT_t => RHOT_tp, &
       QTRC_t => QTRC_tp
    use mod_atmos_phy_mp_vars, only: &
       DENS_t_MP => ATMOS_PHY_MP_DENS_t,    &
       MOMZ_t_MP => ATMOS_PHY_MP_MOMZ_t,    &
       MOMX_t_MP => ATMOS_PHY_MP_MOMX_t,    &
       MOMY_t_MP => ATMOS_PHY_MP_MOMY_t,    &
       RHOT_t_MP => ATMOS_PHY_MP_RHOT_t,    &
       QTRC_t_MP => ATMOS_PHY_MP_QTRC_t,    &
       SFLX_rain => ATMOS_PHY_MP_SFLX_rain, &
       SFLX_snow => ATMOS_PHY_MP_SFLX_snow
    implicit none

    logical, intent(in) :: update_flag
    logical, intent(in) :: history_flag

    real(RP) :: DENS0(KA,IA,JA)
    real(RP) :: MOMZ0(KA,IA,JA)
    real(RP) :: MOMX0(KA,IA,JA)
    real(RP) :: MOMY0(KA,IA,JA)
    real(RP) :: RHOT0(KA,IA,JA)
    real(RP) :: QTRC0(KA,IA,JA,QA)

    real(RP) :: precip(IA,JA)
    real(RP) :: RHOQ(KA,IA,JA)
    real(RP) :: total ! dummy

    integer :: k, i, j, iq
    !---------------------------------------------------------------------------

    if ( update_flag ) then

       do j  = JS, JE
       do i  = IS, IE
       do k  = KS, KE
          DENS0(k,i,j) = DENS(k,i,j) ! save
          MOMZ0(k,i,j) = MOMZ(k,i,j) ! save
          MOMX0(k,i,j) = MOMX(k,i,j) ! save
          MOMY0(k,i,j) = MOMY(k,i,j) ! save
          RHOT0(k,i,j) = RHOT(k,i,j) ! save
       enddo
       enddo
       enddo

       do iq = 1, QA
       do j  = JS, JE
       do i  = IS, IE
       do k  = KS, KE
          QTRC0(k,i,j,iq) = QTRC(k,i,j,iq) ! save
       enddo
       enddo
       enddo
       enddo

       call ATMOS_PHY_MP( DENS0    (:,:,:),   & ! [INOUT]
                          MOMZ0    (:,:,:),   & ! [INOUT]
                          MOMX0    (:,:,:),   & ! [INOUT]
                          MOMY0    (:,:,:),   & ! [INOUT]
                          RHOT0    (:,:,:),   & ! [INOUT]
                          QTRC0    (:,:,:,:), & ! [INOUT]
                          SFLX_rain(:,:),     & ! [OUT]
                          SFLX_snow(:,:)      ) ! [OUT]

       do j  = JS, JE
       do i  = IS, IE
       do k  = KS, KE
          DENS_t_MP(k,i,j) = DENS0(k,i,j) - DENS(k,i,j)
          MOMZ_t_MP(k,i,j) = MOMZ0(k,i,j) - MOMZ(k,i,j)
          MOMX_t_MP(k,i,j) = MOMX0(k,i,j) - MOMX(k,i,j)
          MOMY_t_MP(k,i,j) = MOMY0(k,i,j) - MOMY(k,i,j)
          RHOT_t_MP(k,i,j) = RHOT0(k,i,j) - RHOT(k,i,j)
       enddo
       enddo
       enddo

       do iq = 1, QA
       do j  = JS, JE
       do i  = IS, IE
       do k  = KS, KE
          QTRC_t_MP(k,i,j,iq) = QTRC0(k,i,j,iq) - QTRC(k,i,j,iq)
       enddo
       enddo
       enddo
       enddo

       precip(:,:) = SFLX_rain(:,:) + SFLX_snow(:,:)

       if ( history_flag ) then
          call HIST_in( SFLX_rain(:,:), 'RAIN', 'surface rain rate',          'kg/m2/s', dt_MP )
          call HIST_in( SFLX_snow(:,:), 'SNOW', 'surface snow rate',          'kg/m2/s', dt_MP )
          call HIST_in( precip   (:,:), 'PREC', 'surface precipitation rate', 'kg/m2/s', dt_MP )
       endif

    endif

    !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       DENS_t(k,i,j) = DENS_t(k,i,j) + DENS_t_MP(k,i,j)
       MOMZ_t(k,i,j) = MOMZ_t(k,i,j) + MOMZ_t_MP(k,i,j)
       MOMX_t(k,i,j) = MOMX_t(k,i,j) + MOMX_t_MP(k,i,j)
       MOMY_t(k,i,j) = MOMY_t(k,i,j) + MOMY_t_MP(k,i,j)
       RHOT_t(k,i,j) = RHOT_t(k,i,j) + RHOT_t_MP(k,i,j)
    enddo
    enddo
    enddo

    !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(3)
    do iq = 1,  QA
    do j  = JS, JE
    do i  = IS, IE
    do k  = KS, KE
       QTRC_t(k,i,j,iq) = QTRC_t(k,i,j,iq) + QTRC_t_MP(k,i,j,iq)
    enddo
    enddo
    enddo
    enddo

    if ( STAT_checktotal ) then
       call STAT_total( total, DENS_t_MP(:,:,:), 'DENS_t_MP' )
       call STAT_total( total, MOMZ_t_MP(:,:,:), 'MOMZ_t_MP' )
       call STAT_total( total, MOMX_t_MP(:,:,:), 'MOMX_t_MP' )
       call STAT_total( total, MOMY_t_MP(:,:,:), 'MOMY_t_MP' )
       call STAT_total( total, RHOT_t_MP(:,:,:), 'RHOT_t_MP' )

       do iq = 1, QA
          RHOQ(:,:,:) = DENS(:,:,:) * QTRC_t_MP(:,:,:,iq)

          call STAT_total( total, RHOQ(:,:,:), trim(AQ_NAME(iq))//'_t_MP' )
       enddo
    endif

    return
  end subroutine ATMOS_PHY_MP_driver

end module mod_atmos_phy_mp_driver
