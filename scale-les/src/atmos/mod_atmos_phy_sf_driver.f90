!-------------------------------------------------------------------------------
!> module ATMOSPHERE / Physics Surface fluxes
!!
!! @par Description
!!          Flux from/to bottom boundary of atmosphere (surface)
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2013-12-05 (S.Nishizawa)  [new]
!<
!-------------------------------------------------------------------------------
#include "inc_openmp.h"
module mod_atmos_phy_sf_driver
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_stdio
  use scale_prof
  use scale_grid_index
  use scale_tracer

  use scale_const, only: &
     I_SW  => CONST_I_SW, &
     I_LW  => CONST_I_LW
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: ATMOS_PHY_SF_driver_setup
  public :: ATMOS_PHY_SF_driver

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
  subroutine ATMOS_PHY_SF_driver_setup
    use scale_atmos_phy_sf, only: &
       ATMOS_PHY_SF_setup
    use mod_atmos_admin, only: &
       ATMOS_PHY_SF_TYPE, &
       ATMOS_sw_phy_sf
    use mod_atmos_phy_sf_vars, only: &
       SFC_Z0M   => ATMOS_PHY_SF_SFC_Z0M,   &
       SFC_Z0H   => ATMOS_PHY_SF_SFC_Z0H,   &
       SFC_Z0E   => ATMOS_PHY_SF_SFC_Z0E,   &
       SFLX_MW   => ATMOS_PHY_SF_SFLX_MW,   &
       SFLX_MU   => ATMOS_PHY_SF_SFLX_MU,   &
       SFLX_MV   => ATMOS_PHY_SF_SFLX_MV,   &
       SFLX_SH   => ATMOS_PHY_SF_SFLX_SH,   &
       SFLX_LH   => ATMOS_PHY_SF_SFLX_LH,   &
       SFLX_QTRC => ATMOS_PHY_SF_SFLX_QTRC
    use mod_cpl_admin, only: &
       CPL_sw
    implicit none
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '++++++ Module[DRIVER] / Categ[ATMOS PHY_SF] / Origin[SCALE-LES]'

    if ( ATMOS_sw_phy_sf ) then

       ! setup library component
       call ATMOS_PHY_SF_setup( ATMOS_PHY_SF_TYPE )

       if ( .NOT. CPL_sw ) then
          if( IO_L ) write(IO_FID_LOG,*) '*** Coupler is disabled.'
          if( IO_L ) write(IO_FID_LOG,*) '*** SFC_Z0[MHE] is assumed to be 0.'
          SFC_Z0M(:,:) = 0.0_RP
          SFC_Z0H(:,:) = 0.0_RP
          SFC_Z0E(:,:) = 0.0_RP
       endif

       ! run once (only for the diagnostic value)
       call PROF_rapstart('ATM SurfaceFlux', 1)
       call ATMOS_PHY_SF_driver( update_flag = .true. )
       call PROF_rapend  ('ATM SurfaceFlux', 1)

    else

       if( IO_L ) write(IO_FID_LOG,*) '*** this component is never called.'
       if( IO_L ) write(IO_FID_LOG,*) '*** surface fluxes are set to zero.'
       SFLX_MW  (:,:)   = 0.0_RP
       SFLX_MU  (:,:)   = 0.0_RP
       SFLX_MV  (:,:)   = 0.0_RP
       SFLX_SH  (:,:)   = 0.0_RP
       SFLX_LH  (:,:)   = 0.0_RP
       SFLX_QTRC(:,:,:) = 0.0_RP

       if( IO_L ) write(IO_FID_LOG,*) '*** SFC_TEMP, SFC_albedo is set in ATMOS_PHY_SF_vars.'

    endif

    return
  end subroutine ATMOS_PHY_SF_driver_setup

  !-----------------------------------------------------------------------------
  !> Driver
  subroutine ATMOS_PHY_SF_driver( update_flag )
    use scale_const, only: &
       CPdry => CONST_CPdry
    use scale_comm, only: &
       COMM_vars8, &
       COMM_wait
    use scale_grid, only: &
       RCDZ => GRID_RCDZ, &
       RFDZ => GRID_RFDZ
    use scale_grid_real, only: &
       REAL_CZ, &
       REAL_Z1
    use scale_gridtrans, only: &
       GSQRT => GTRANS_GSQRT, &
       I_XYZ,                 &
       I_UYZ,                 &
       I_XVZ,                 &
       I_XYW
    use scale_topography, only: &
       TOPO_Zsfc
    use scale_time, only: &
       dt_SF => TIME_DTSEC_ATMOS_PHY_SF
    use scale_statistics, only: &
       STATISTICS_checktotal, &
       STAT_total
    use scale_history, only: &
       HIST_in
    use scale_atmos_bottom, only: &
       BOTTOM_estimate => ATMOS_BOTTOM_estimate
    use scale_atmos_hydrostatic, only: &
       barometric_law_mslp => ATMOS_HYDROSTATIC_barometric_law_mslp
    use scale_atmos_phy_sf, only: &
       ATMOS_PHY_SF
    use mod_atmos_vars, only: &
       DENS,              &
       MOMZ,              &
       MOMX,              &
       MOMY,              &
       RHOT,              &
       QTRC,              &
       TEMP,              &
       PRES,              &
       W,                 &
       U,                 &
       V,                 &
       DENS_t => DENS_tp, &
       MOMZ_t => MOMZ_tp, &
       MOMX_t => MOMX_tp, &
       MOMY_t => MOMY_tp, &
       RHOT_t => RHOT_tp, &
       RHOQ_t => RHOQ_tp
    use mod_atmos_phy_rd_vars, only: &
       SFLX_LW_dn => ATMOS_PHY_RD_SFLX_LW_dn, &
       SFLX_SW_dn => ATMOS_PHY_RD_SFLX_SW_dn
    use mod_atmos_phy_sf_vars, only: &
       DENS_t_SF  => ATMOS_PHY_SF_DENS_t,     &
       MOMZ_t_SF  => ATMOS_PHY_SF_MOMZ_t,     &
       MOMX_t_SF  => ATMOS_PHY_SF_MOMX_t,     &
       MOMY_t_SF  => ATMOS_PHY_SF_MOMY_t,     &
       RHOT_t_SF  => ATMOS_PHY_SF_RHOT_t,     &
       RHOQ_t_SF  => ATMOS_PHY_SF_RHOQ_t,     &
       SFC_DENS   => ATMOS_PHY_SF_SFC_DENS,   &
       SFC_PRES   => ATMOS_PHY_SF_SFC_PRES,   &
       SFC_TEMP   => ATMOS_PHY_SF_SFC_TEMP,   &
       SFC_albedo => ATMOS_PHY_SF_SFC_albedo, &
       SFC_Z0M    => ATMOS_PHY_SF_SFC_Z0M,    &
       SFC_Z0H    => ATMOS_PHY_SF_SFC_Z0H,    &
       SFC_Z0E    => ATMOS_PHY_SF_SFC_Z0E,    &
       SFLX_MW    => ATMOS_PHY_SF_SFLX_MW,    &
       SFLX_MU    => ATMOS_PHY_SF_SFLX_MU,    &
       SFLX_MV    => ATMOS_PHY_SF_SFLX_MV,    &
       SFLX_SH    => ATMOS_PHY_SF_SFLX_SH,    &
       SFLX_LH    => ATMOS_PHY_SF_SFLX_LH,    &
       SFLX_GH    => ATMOS_PHY_SF_SFLX_GH,    &
       SFLX_QTRC  => ATMOS_PHY_SF_SFLX_QTRC,  &
       U10        => ATMOS_PHY_SF_U10,        &
       V10        => ATMOS_PHY_SF_V10,        &
       T2         => ATMOS_PHY_SF_T2,         &
       Q2         => ATMOS_PHY_SF_Q2
    use mod_cpl_admin, only: &
       CPL_sw
    implicit none

    logical, intent(in) :: update_flag

    real(RP) :: Uabs10(IA,JA) ! 10m absolute wind [m/s]
    real(RP) :: MSLP  (IA,JA) ! mean sea-level pressure [Pa]
    real(RP) :: beta  (IA,JA)
    real(RP) :: total ! dummy
    real(RP) :: work

    integer :: i, j, iq
    !---------------------------------------------------------------------------

    if ( update_flag ) then

       ! update surface density, surface pressure
       call BOTTOM_estimate( DENS     (:,:,:), & ! [IN]
                             PRES     (:,:,:), & ! [IN]
                             REAL_CZ  (:,:,:), & ! [IN]
                             TOPO_Zsfc(:,:),   & ! [IN]
                             REAL_Z1  (:,:),   & ! [IN]
                             SFC_DENS (:,:),   & ! [OUT]
                             SFC_PRES (:,:)    ) ! [OUT]

       if ( .NOT. CPL_sw ) then
          beta(:,:) = 1.0_RP

          call ATMOS_PHY_SF( TEMP      (KS,:,:),   & ! [IN]
                             PRES      (KS,:,:),   & ! [IN]
                             W         (KS,:,:),   & ! [IN]
                             U         (KS,:,:),   & ! [IN]
                             V         (KS,:,:),   & ! [IN]
                             DENS      (KS,:,:),   & ! [IN]
                             QTRC      (KS,:,:,:), & ! [IN]
                             REAL_Z1   (:,:),      & ! [IN]
                             dt_SF,                & ! [IN]
                             SFC_DENS  (:,:),      & ! [IN]
                             SFC_PRES  (:,:),      & ! [IN]
                             SFLX_LW_dn(:,:),      & ! [IN]
                             SFLX_SW_dn(:,:),      & ! [IN]
                             SFC_TEMP  (:,:),      & ! [IN]
                             SFC_albedo(:,:,:),    & ! [IN]
                             beta      (:,:),      & ! [IN]
                             SFC_Z0M   (:,:),      & ! [INOUT]
                             SFC_Z0H   (:,:),      & ! [INOUT]
                             SFC_Z0E   (:,:),      & ! [INOUT]
                             SFLX_MW   (:,:),      & ! [OUT]
                             SFLX_MU   (:,:),      & ! [OUT]
                             SFLX_MV   (:,:),      & ! [OUT]
                             SFLX_SH   (:,:),      & ! [OUT]
                             SFLX_LH   (:,:),      & ! [OUT]
                             SFLX_QTRC (:,:,:),    & ! [OUT]
                             U10       (:,:),      & ! [OUT]
                             V10       (:,:),      & ! [OUT]
                             T2        (:,:),      & ! [OUT]
                             Q2        (:,:)       ) ! [OUT]
       endif

       do j = JS, JE
       do i = IS, IE
          Uabs10(i,j) = sqrt( U10(i,j)**2 + V10(i,j)**2 )
       end do
       end do

       call barometric_law_mslp( MSLP(:,:), SFC_PRES(:,:), T2(:,:), TOPO_Zsfc(:,:) )

       call HIST_in( SFC_DENS  (:,:),      'SFC_DENS',   'surface atmospheric density',       'kg/m3'   )
       call HIST_in( SFC_PRES  (:,:),      'SFC_PRES',   'surface atmospheric pressure',      'Pa'      )
       call HIST_in( SFC_TEMP  (:,:),      'SFC_TEMP',   'surface skin temperature (merged)', 'K'       )
       call HIST_in( SFC_albedo(:,:,I_LW), 'SFC_ALB_LW', 'surface albedo (longwave,merged)',  '0-1'     )
       call HIST_in( SFC_albedo(:,:,I_SW), 'SFC_ALB_SW', 'surface albedo (shortwave,merged)', '0-1'     )
       call HIST_in( SFC_Z0M   (:,:),      'SFC_Z0M',    'roughness length (momentum)',       'm'       )
       call HIST_in( SFC_Z0H   (:,:),      'SFC_Z0H',    'roughness length (heat)',           'm'       )
       call HIST_in( SFC_Z0E   (:,:),      'SFC_Z0E',    'roughness length (vapor)',          'm'       )
       call HIST_in( SFLX_MW   (:,:),      'MWFLX',      'w-momentum flux (merged)',          'kg/m2/s' )
       call HIST_in( SFLX_MU   (:,:),      'MUFLX',      'u-momentum flux (merged)',          'kg/m2/s' )
       call HIST_in( SFLX_MV   (:,:),      'MVFLX',      'v-momentum flux (merged)',          'kg/m2/s' )
       call HIST_in( SFLX_SH   (:,:),      'SHFLX',      'sensible heat flux (merged)',       'W/m2'    )
       call HIST_in( SFLX_LH   (:,:),      'LHFLX',      'latent heat flux (merged)',         'W/m2'    )
       call HIST_in( SFLX_GH   (:,:),      'GHFLX',      'ground heat flux (merged)',         'W/m2'    )
       call HIST_in( Uabs10    (:,:),      'Uabs10',     '10m absolute wind',                 'm/s'     )
       call HIST_in( U10       (:,:),      'U10',        '10m x-wind',                        'm/s'     )
       call HIST_in( V10       (:,:),      'V10',        '10m y-wind',                        'm/s'     )
       call HIST_in( T2        (:,:),      'T2 ',        '2m temperature',                    'K'       )
       call HIST_in( Q2        (:,:),      'Q2 ',        '2m water vapor',                    'kg/kg'   )
       call HIST_in( MSLP      (:,:),      'MSLP',       'mean sea-level pressure',           'Pa'      )

       call COMM_vars8( SFLX_MU(:,:), 1 )
       call COMM_vars8( SFLX_MV(:,:), 2 )
       call COMM_wait ( SFLX_MU(:,:), 1 )
       call COMM_wait ( SFLX_MV(:,:), 2 )

       do j = JS, JE
       do i = IS, IE
          MOMZ_t_SF(i,j) = SFLX_MW(i,j) * RFDZ(KS) / GSQRT(KS,i,j,I_XYW)
       enddo
       enddo

       do j = JS, JE
       do i = IS, IE
          MOMX_t_SF(i,j) = 0.5_RP * ( SFLX_MU(i,j) + SFLX_MU(i+1,j) ) * RCDZ(KS) / GSQRT(KS,i,j,I_UYZ)
       enddo
       enddo

       do j = JS, JE
       do i = IS, IE
          MOMY_t_SF(i,j) = 0.5_RP * ( SFLX_MV(i,j) + SFLX_MV(i,j+1) ) * RCDZ(KS) / GSQRT(KS,i,j,I_XVZ)
       enddo
       enddo

       do j = JS, JE
       do i = IS, IE
          RHOT_t_SF(i,j) = SFLX_SH(i,j) / CPdry * RCDZ(KS) / GSQRT(KS,i,j,I_XYZ)
       enddo
       enddo

       DENS_t_SF(:,:) = 0.0_RP
       do iq = I_QV, I_QV
       do j  = JS, JE
       do i  = IS, IE
          work = SFLX_QTRC(i,j,iq) * RCDZ(KS) / GSQRT(KS,i,j,I_XYZ)
          DENS_t_SF(i,j)    = DENS_t_SF(i,j) + work
          RHOQ_t_SF(i,j,iq) = work
       enddo
       enddo
       enddo

       do j  = JS, JE
       do i  = IS, IE
          RHOT_t_SF(i,j) = RHOT_t_SF(i,j) + DENS_t_SF(i,j) * RHOT(KS,i,j) / DENS(KS,i,j)
       enddo
       enddo

    endif

    do j = JS, JE
    do i = IS, IE
       DENS_t(KS,i,j) = DENS_t(KS,i,j) + DENS_t_SF(i,j)
       MOMZ_t(KS,i,j) = MOMZ_t(KS,i,j) + MOMZ_t_SF(i,j)
       MOMX_t(KS,i,j) = MOMX_t(KS,i,j) + MOMX_t_SF(i,j)
       MOMY_t(KS,i,j) = MOMY_t(KS,i,j) + MOMY_t_SF(i,j)
       RHOT_t(KS,i,j) = RHOT_t(KS,i,j) + RHOT_t_SF(i,j)
    enddo
    enddo

    do iq = I_QV, I_QV
    do j  = JS, JE
    do i  = IS, IE
       RHOQ_t(KS,i,j,iq) = RHOQ_t(KS,i,j,iq) + RHOQ_t_SF(i,j,iq)
    enddo
    enddo
    enddo

    if ( STATISTICS_checktotal ) then
       call STAT_total( total, DENS_t_SF(:,:), 'DENS_t_SF' )
       call STAT_total( total, MOMZ_t_SF(:,:), 'MOMZ_t_SF' )
       call STAT_total( total, MOMX_t_SF(:,:), 'MOMX_t_SF' )
       call STAT_total( total, MOMY_t_SF(:,:), 'MOMY_t_SF' )
       call STAT_total( total, RHOT_t_SF(:,:), 'RHOT_t_SF' )

       do iq = I_QV, I_QV
          call STAT_total( total, RHOQ_t_SF(:,:,iq), trim(AQ_NAME(iq))//'_t_SF' )
       enddo
    endif

    return
  end subroutine ATMOS_PHY_SF_driver

end module mod_atmos_phy_sf_driver
