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
       call PROF_rapstart('ATM_Aerosol', 1)
       call ATMOS_PHY_AE_driver( update_flag = .true. )
       call PROF_rapend  ('ATM_Aerosol', 1)

!    else
!       if( IO_L ) write(IO_FID_LOG,*) '*** this component is never called.'
!    endif

    return
  end subroutine ATMOS_PHY_AE_driver_setup

  !-----------------------------------------------------------------------------
  !> Driver
  subroutine ATMOS_PHY_AE_driver( update_flag )
    use scale_time, only: &
       dt_AE => TIME_DTSEC_ATMOS_PHY_AE
    use scale_les_statistics, only: &
       STATISTICS_checktotal, &
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
       RHOQ_t => RHOQ_tp
    use mod_atmos_phy_ae_vars, only: &
       RHOQ_t_AE => ATMOS_PHY_AE_RHOQ_t, &
       CCN       => ATMOS_PHY_AE_CCN,   &
       CCN_t     => ATMOS_PHY_AE_CCN_t, &
       AE_EMIT   => ATMOS_PHY_AE_EMIT
    use mod_atmos_phy_mp_vars, only: &
       EVAPORATE => ATMOS_PHY_MP_EVAPORATE
    implicit none

    logical, intent(in) :: update_flag

    real(RP) :: QTRC0(KA,IA,JA,QA)
    real(RP) :: CN(KA,IA,JA)
    real(RP) :: NREG(KA,IA,JA)

    real(RP) :: total ! dummy

    integer  :: k, i, j, iq
    !---------------------------------------------------------------------------

    if ( update_flag ) then

!OCL XFILL
       do iq = 1, QA
       do j  = JS, JE
       do i  = IS, IE
       do k  = KS, KE
          QTRC0(k,i,j,iq) = QTRC(k,i,j,iq) ! save
       enddo
       enddo
       enddo
       enddo

       CCN(:,:,:) = 0.0_RP ! reset
       CCN_t(:,:,:) = 0.0_RP ! reset
       RHOQ_t_AE(:,:,:,:) = 0.0_RP ! reset

       NREG(:,:,:) = EVAPORATE(:,:,:) * dt_AE 

       call ATMOS_PHY_AE( DENS, & ! [IN]
                          MOMZ, & ! [IN]
                          MOMX, & ! [IN]
                          MOMY, & ! [IN]
                          RHOT, & ! [IN]
                          AE_EMIT, & ! [IN]
                          NREG,    & ! [IN]
                          CN ,  & ! [OUT]
                          CCN,  & ! [OUT]
                          QTRC0 ) ! [INOUT]

       do iq = 1, QA
       do j  = JS, JE
       do i  = IS, IE
       do k  = KS, KE
          RHOQ_t_AE(k,i,j,iq) = ( QTRC0(k,i,j,iq) - QTRC(k,i,j,iq) ) * DENS(k,i,j) / dt_AE
       enddo
       enddo
       enddo
       enddo

       CCN_t(:,:,:) = CCN(:,:,:) / dt_AE 

!       CCN(:,:,:) = 0.0_RP ! tentative

       call HIST_in( CN(:,:,:)*1e-6_RP,  'CN',  'condensation nucrei', '/cc' )
       call HIST_in( CCN(:,:,:)*1e-6_RP, 'CCN', 'cloud condensation nucrei', '/cc' )

    endif

    !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(3)
    do iq = 1, QA
    do j  = JS, JE
    do i  = IS, IE
    do k  = KS, KE
       RHOQ_t(k,i,j,iq) = RHOQ_t(k,i,j,iq) + RHOQ_t_AE(k,i,j,iq)
    enddo
    enddo
    enddo
    enddo

    if ( STATISTICS_checktotal ) then
       do iq = 1, QA
          call STAT_total( total, RHOQ_t_AE(:,:,:,iq), trim(AQ_NAME(iq))//'_t_AE' )
       enddo
    endif

    return
  end subroutine ATMOS_PHY_AE_driver

end module mod_atmos_phy_ae_driver
