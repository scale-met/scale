!-------------------------------------------------------------------------------
!> module ATMOSPHERE / Physics Surface fluxes
!!
!! @par Description
!!          Flux from/to bottom boundary of atmosphere (surface)
!!
!! @author Team SCALE
!!
!<
!-------------------------------------------------------------------------------
#include "scalelib.h"
module mod_atmos_phy_sf_driver
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_stdio
  use scale_prof
  use scale_atmos_grid_cartesC_index
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
  public :: ATMOS_PHY_SF_driver_calc_tendency

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
    use scale_prc, only: &
       PRC_abort
    use scale_atmos_phy_sf_bulk, only: &
       ATMOS_PHY_SF_bulk_setup
    use scale_atmos_phy_sf_const, only: &
       ATMOS_PHY_SF_const_setup
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

    LOG_NEWLINE
    LOG_PROGRESS(*) 'Module[DRIVER] / Categ[ATMOS PHY_SF] / Origin[SCALE-RM]'

    if ( ATMOS_sw_phy_sf ) then

       if ( CPL_sw ) then
          LOG_INFO("ATMOS_PHY_SF_driver_setup",*) 'Coupler is enabled.'
       else
          ! setup library component
          select case( ATMOS_PHY_SF_TYPE )
          case ( 'BULK' )
             call ATMOS_PHY_SF_bulk_setup
          case ( 'CONST' )
             call ATMOS_PHY_SF_const_setup
          case default
             LOG_ERROR("ATMOS_PHY_SF_driver_setup",*) 'invalid Surface flux type(', trim(ATMOS_PHY_SF_TYPE), '). CHECK!'
             call PRC_abort
          end select
       endif

    else

       LOG_INFO("ATMOS_PHY_SF_driver_setup",*) 'this component is never called.'
       LOG_INFO("ATMOS_PHY_SF_driver_setup",*) 'surface fluxes are set to zero.'
       SFLX_MW  (:,:)   = 0.0_RP
       SFLX_MU  (:,:)   = 0.0_RP
       SFLX_MV  (:,:)   = 0.0_RP
       SFLX_SH  (:,:)   = 0.0_RP
       SFLX_LH  (:,:)   = 0.0_RP
       LOG_INFO("ATMOS_PHY_SF_driver_setup",*) 'SFC_TEMP, SFC_albedo is set in ATMOS_PHY_SF_vars.'

    endif

    SFLX_QTRC(:,:,:) = 0.0_RP

    return
  end subroutine ATMOS_PHY_SF_driver_setup

  !-----------------------------------------------------------------------------
  !> calculation tendency
  subroutine ATMOS_PHY_SF_driver_calc_tendency( update_flag )
    use scale_const, only: &
       GRAV   => CONST_GRAV,   &
       KARMAN => CONST_KARMAN, &
       CPdry  => CONST_CPdry
    use scale_atmos_grid_cartesC_real, only: &
       CZ => ATMOS_GRID_CARTESC_REAL_CZ, &
       FZ => ATMOS_GRID_CARTESC_REAL_FZ, &
       Z1 => ATMOS_GRID_CARTESC_REAL_Z1, &
       ATMOS_GRID_CARTESC_REAL_AREA, &
       ATMOS_GRID_CARTESC_REAL_TOTAREA
    use scale_topography, only: &
       TOPO_Zsfc
    use scale_time, only: &
       dt_SF => TIME_DTSEC_ATMOS_PHY_SF
    use scale_statistics, only: &
       STATISTICS_checktotal, &
       STATISTICS_total
    use scale_atmos_bottom, only: &
       BOTTOM_estimate => ATMOS_BOTTOM_estimate
    use scale_atmos_hydrometeor, only: &
       ATMOS_HYDROMETEOR_dry, &
       I_QV
    use scale_atmos_phy_sf_bulk, only: &
       ATMOS_PHY_SF_bulk_flux
    use scale_atmos_phy_sf_const, only: &
       ATMOS_PHY_SF_const_flux
    use mod_atmos_vars, only: &
       DENS   => DENS_av, &
       RHOT   => RHOT_av, &
       POTT,              &
       TEMP,              &
       PRES,              &
       W,                 &
       U,                 &
       V,                 &
       QV,                &
       DENS_t => DENS_tp, &
       MOMZ_t => MOMZ_tp, &
       RHOU_t => RHOU_tp, &
       RHOV_t => RHOV_tp, &
       RHOH   => RHOH_p,  &
       RHOT_t => RHOT_tp, &
       RHOQ_t => RHOQ_tp
    use mod_atmos_phy_rd_vars, only: &
       SFLX_LW_dn => ATMOS_PHY_RD_SFLX_LW_dn, &
       SFLX_SW_dn => ATMOS_PHY_RD_SFLX_SW_dn
    use mod_atmos_phy_bl_vars, only: &
       PBL_Zi => ATMOS_PHY_BL_Zi
    use mod_atmos_phy_sf_vars, only: &
       DENS_t_SF  => ATMOS_PHY_SF_DENS_t,     &
       MOMZ_t_SF  => ATMOS_PHY_SF_MOMZ_t,     &
       RHOU_t_SF  => ATMOS_PHY_SF_RHOU_t,     &
       RHOV_t_SF  => ATMOS_PHY_SF_RHOV_t,     &
       RHOH_SF    => ATMOS_PHY_SF_RHOH,       &
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
       Q2         => ATMOS_PHY_SF_Q2,         &
       l_mo       => ATMOS_PHY_SF_l_mo
    use mod_cpl_admin, only: &
       CPL_sw
    use mod_atmos_admin, only: &
       ATMOS_PHY_SF_TYPE
    implicit none

    logical, intent(in) :: update_flag

    real(RP) :: ATM_W   (IA,JA)
    real(RP) :: ATM_U   (IA,JA)
    real(RP) :: ATM_V   (IA,JA)
    real(RP) :: ATM_DENS(IA,JA)
    real(RP) :: ATM_TEMP(IA,JA)
    real(RP) :: ATM_PRES(IA,JA)
    real(RP) :: ATM_QV  (IA,JA)
    real(RP) :: SFLX_QV (IA,JA)

    real(RP) :: Z0M_t(IA,JA)
    real(RP) :: Z0H_t(IA,JA)
    real(RP) :: Z0E_t(IA,JA)

    real(RP) :: q(QA)
    real(RP) :: qdry
    real(RP) :: Rtot
    real(RP) :: CPtot

    real(RP) :: us, SFLX_PT

    real(RP) :: work

    integer :: i, j, iq
    !---------------------------------------------------------------------------

    if ( update_flag ) then

       ! update surface density, surface pressure
       call BOTTOM_estimate( KA, KS, KE, IA, ISB, IEB, JA, JSB, JEB, &
                             DENS(:,:,:), PRES(:,:,:),           & ! [IN]
                             CZ(:,:,:), TOPO_Zsfc(:,:), Z1(:,:), & ! [IN]
                             SFC_DENS(:,:), SFC_PRES(:,:)        ) ! [OUT]

       if ( .NOT. CPL_sw ) then

          !omp parallel do
          do j = JSB, JEB
          do i = ISB, IEB
             ATM_W   (i,j) = W   (KS,i,j)
             ATM_U   (i,j) = U   (KS,i,j)
             ATM_V   (i,j) = V   (KS,i,j)
             ATM_DENS(i,j) = DENS(KS,i,j)
             ATM_TEMP(i,j) = TEMP(KS,i,j)
             ATM_PRES(i,j) = PRES(KS,i,j)
             ATM_QV  (i,j) = QV  (KS,i,j)
          end do
          end do

          select case ( ATMOS_PHY_SF_TYPE )
          case ( 'BULK' )
             call ATMOS_PHY_SF_bulk_flux( &
                  IA, ISB, IEB, JA, JSB, JEB, &
                  ATM_W(:,:), ATM_U(:,:), ATM_V(:,:),          & ! [IN]
                  ATM_TEMP(:,:), ATM_PRES(:,:), ATM_QV(:,:),   & ! [IN]
                  SFC_DENS(:,:), SFC_TEMP(:,:), SFC_PRES(:,:), & ! [IN]
                  SFC_Z0M(:,:), SFC_Z0H(:,:), SFC_Z0E(:,:),    & ! [IN]
                  PBL_Zi(:,:), Z1(:,:),                        & ! [IN]
                  SFLX_MW(:,:), SFLX_MU(:,:), SFLX_MV(:,:),    & ! [OUT]
                  SFLX_SH(:,:), SFLX_LH(:,:), SFLX_QV(:,:),    & ! [OUT]
                  U10(:,:), V10(:,:), T2(:,:), Q2(:,:)         ) ! [OUT]

          case ( 'CONST' )

             call ATMOS_PHY_SF_const_flux( &
                  IA, IS, IE, JA, JS, JE, &
                  ATM_W(:,:), ATM_U(:,:), ATM_V(:,:), ATM_TEMP(:,:), & ! [IN]
                  Z1(:,:), SFC_DENS(:,:),                            & ! [IN]
                  SFLX_MW(:,:), SFLX_MU(:,:), SFLX_MV(:,:),          & ! [OUT]
                  SFLX_SH(:,:), SFLX_LH(:,:), SFLX_QV(:,:),          & ! [OUT]
                  U10(:,:), V10(:,:)                                 ) ! [OUT]
             T2(:,:) = ATM_TEMP(:,:)
             Q2(:,:) = ATM_QV(:,:)

          end select

          if ( .not. ATMOS_HYDROMETEOR_dry ) then
             SFLX_QTRC(:,:,I_QV) = SFLX_QV(:,:)
          end if

       endif

       ! temtative
       !$omp parallel do
       do j = JSB, JEB
       do i = ISB, IEB
          us = max( 1.E-6_RP, &
                    sqrt( sqrt( SFLX_MU(i,j)**2 + SFLX_MV(i,j)**2 ) / DENS(KS,i,j) ) ) ! frictional velocity
          SFLX_PT = SFLX_SH(i,j) / ( CPdry * DENS(KS,i,j) ) &
                  * POTT(KS,i,j) / TEMP(KS,i,j)
          l_mo(i,j) = - us**3 * POTT(KS,i,j) / ( KARMAN * GRAV * SFLX_PT )
       end do
       end do

       call history_output

       !omp parallel do
!OCL XFILL
       do j = JSB, JEB
       do i = ISB, IEB
          MOMZ_t_SF(i,j) = SFLX_MW(i,j) / ( CZ(KS+1,i,j) - CZ(KS,i,j) )
          RHOU_t_SF(i,j) = SFLX_MU(i,j) / ( FZ(KS,i,j) - FZ(KS-1,i,j) )
          RHOV_t_SF(i,j) = SFLX_MV(i,j) / ( FZ(KS,i,j) - FZ(KS-1,i,j) )
          RHOH_SF  (i,j) = SFLX_SH(i,j) / ( FZ(KS,i,j) - FZ(KS-1,i,j) )
       enddo
       enddo

       if ( .not. ATMOS_HYDROMETEOR_dry ) then
          !omp parallel do
          do j = JSB, JEB
          do i = ISB, IEB
             work = SFLX_QTRC(i,j,I_QV) / ( FZ(KS,i,j) - FZ(KS-1,i,j) )
             DENS_t_SF(i,j)      = work
             RHOQ_t_SF(i,j,I_QV) = work
             RHOT_t_SF(i,j) = work * RHOT(KS,i,j) / DENS(KS,i,j)
          enddo
          enddo
       end if

    endif

    !omp parallel do
    do j = JSB, JEB
    do i = ISB, IEB
       MOMZ_t(KS,i,j) = MOMZ_t(KS,i,j) + MOMZ_t_SF(i,j)
       RHOU_t(KS,i,j) = RHOU_t(KS,i,j) + RHOU_t_SF(i,j)
       RHOV_t(KS,i,j) = RHOV_t(KS,i,j) + RHOV_t_SF(i,j)
       RHOH  (KS,i,j) = RHOH  (KS,i,j) + RHOH_SF  (i,j)
    enddo
    enddo

    if ( .not. ATMOS_HYDROMETEOR_dry ) then
       !omp parallel do
       do j  = JS, JE
       do i  = IS, IE
          DENS_t(KS,i,j) = DENS_t(KS,i,j) + DENS_t_SF(i,j)
          RHOQ_t(KS,i,j,I_QV) = RHOQ_t(KS,i,j,I_QV) + RHOQ_t_SF(i,j,I_QV)
          RHOT_t(KS,i,j) = RHOT_t(KS,i,j) + RHOT_t_SF(i,j)
       enddo
       enddo
    end if

    if ( STATISTICS_checktotal ) then
       call STATISTICS_total( IA, IS, IE, JA, JS, JE, &
                              DENS_t_SF(:,:), 'DENS_t_SF',       &
                              ATMOS_GRID_CARTESC_REAL_AREA(:,:), &
                              ATMOS_GRID_CARTESC_REAL_TOTAREA    )
       call STATISTICS_total( IA, IS, IE, JA, JS, JE, &
                              MOMZ_t_SF(:,:), 'MOMZ_t_SF',       &
                              ATMOS_GRID_CARTESC_REAL_AREA(:,:), &
                              ATMOS_GRID_CARTESC_REAL_TOTAREA    )
       call STATISTICS_total( IA, IS, IE, JA, JS, JE, &
                              RHOU_t_SF(:,:), 'RHOU_t_SF',       &
                              ATMOS_GRID_CARTESC_REAL_AREA(:,:), &
                              ATMOS_GRID_CARTESC_REAL_TOTAREA    )
       call STATISTICS_total( IA, IS, IE, JA, JS, JE, &
                              RHOV_t_SF(:,:), 'RHOV_t_SF',       &
                              ATMOS_GRID_CARTESC_REAL_AREA(:,:), &
                              ATMOS_GRID_CARTESC_REAL_TOTAREA    )
       call STATISTICS_total( IA, IS, IE, JA, JS, JE, &
                              RHOH_SF  (:,:), 'RHOH_SF',         &
                              ATMOS_GRID_CARTESC_REAL_AREA(:,:), &
                              ATMOS_GRID_CARTESC_REAL_TOTAREA    )

       if ( .not. ATMOS_HYDROMETEOR_dry ) then
          call STATISTICS_total( IA, IS, IE, JA, JS, JE, &
                                 RHOT_t_SF(:,:)     , 'RHOT_t_SF',                      &
                                 ATMOS_GRID_CARTESC_REAL_AREA(:,:),                     &
                                 ATMOS_GRID_CARTESC_REAL_TOTAREA                        )
          call STATISTICS_total( IA, IS, IE, JA, JS, JE, &
                                 RHOQ_t_SF(:,:,I_QV), trim(TRACER_NAME(I_QV))//'_t_SF', &
                                 ATMOS_GRID_CARTESC_REAL_AREA(:,:),                     &
                                 ATMOS_GRID_CARTESC_REAL_TOTAREA                        )
       end if
    endif

    return
  end subroutine ATMOS_PHY_SF_driver_calc_tendency

  subroutine history_output
    use scale_const, only: &
       UNDEF => CONST_UNDEF
    use scale_file_history, only: &
       FILE_HISTORY_in
    use scale_atmos_hydrostatic, only: &
       barometric_law_mslp => ATMOS_HYDROSTATIC_barometric_law_mslp
    use scale_topography, only: &
       TOPO_Zsfc
    use mod_atmos_phy_sf_vars, only: &
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
       U10        => ATMOS_PHY_SF_U10,        &
       V10        => ATMOS_PHY_SF_V10,        &
       T2         => ATMOS_PHY_SF_T2,         &
       Q2         => ATMOS_PHY_SF_Q2
    real(RP) :: MSLP  (IA,JA) ! mean sea-level pressure [Pa]

    real(RP) :: Uabs10(IA,JA) ! 10m absolute wind [m/s]

    integer :: i, j

!OCL XFILL
    !$omp parallel do
    do j = JSB, JEB
    do i = ISB, IEB
       Uabs10(i,j) = sqrt( U10(i,j)**2 + V10(i,j)**2 )
    end do
    end do

    
    call barometric_law_mslp( IA, IS, IE, JA, JS, JE, &
                              SFC_PRES(:,:), T2(:,:), TOPO_Zsfc(:,:), & ! [IN]
                              MSLP(:,:)                               ) ! [OUT]

    call FILE_HISTORY_in( SFC_DENS  (:,:),      'SFC_DENS',   'surface atmospheric density',       'kg/m3'   )
    call FILE_HISTORY_in( SFC_PRES  (:,:),      'SFC_PRES',   'surface atmospheric pressure',      'Pa'      )
    call FILE_HISTORY_in( SFC_TEMP  (:,:),      'SFC_TEMP',   'surface skin temperature (merged)', 'K'       )
    call FILE_HISTORY_in( SFC_albedo(:,:,I_LW), 'SFC_ALB_LW', 'surface albedo (longwave,merged)',  '1'       , fill_halo=.true. )
    call FILE_HISTORY_in( SFC_albedo(:,:,I_SW), 'SFC_ALB_SW', 'surface albedo (shortwave,merged)', '1'       , fill_halo=.true. )
    call FILE_HISTORY_in( SFC_Z0M   (:,:),      'SFC_Z0M',    'roughness length (momentum)',       'm'       , fill_halo=.true. )
    call FILE_HISTORY_in( SFC_Z0H   (:,:),      'SFC_Z0H',    'roughness length (heat)',           'm'       , fill_halo=.true. )
    call FILE_HISTORY_in( SFC_Z0E   (:,:),      'SFC_Z0E',    'roughness length (vapor)',          'm'       , fill_halo=.true. )
    call FILE_HISTORY_in( SFLX_MW   (:,:),      'MWFLX',      'w-momentum flux (merged)',          'kg/m/s2' )
    call FILE_HISTORY_in( SFLX_MU   (:,:),      'MUFLX',      'u-momentum flux (merged)',          'kg/m/s2' )
    call FILE_HISTORY_in( SFLX_MV   (:,:),      'MVFLX',      'v-momentum flux (merged)',          'kg/m/s2' )
    call FILE_HISTORY_in( SFLX_SH   (:,:),      'SHFLX',      'sensible heat flux (merged)',       'W/m2'    , fill_halo=.true. )
    call FILE_HISTORY_in( SFLX_LH   (:,:),      'LHFLX',      'latent heat flux (merged)',         'W/m2'    , fill_halo=.true. )
    call FILE_HISTORY_in( SFLX_GH   (:,:),      'GHFLX',      'ground heat flux (merged)',         'W/m2'    , fill_halo=.true. )
    call FILE_HISTORY_in( Uabs10    (:,:),      'Uabs10',     '10m absolute wind',                 'm/s'     , fill_halo=.true. )
    call FILE_HISTORY_in( U10       (:,:),      'U10',        '10m x-wind',                        'm/s'     , fill_halo=.true. )
    call FILE_HISTORY_in( V10       (:,:),      'V10',        '10m y-wind',                        'm/s'     , fill_halo=.true. )
    call FILE_HISTORY_in( T2        (:,:),      'T2 ',        '2m air temperature',                'K'       , fill_halo=.true. )
    call FILE_HISTORY_in( Q2        (:,:),      'Q2 ',        '2m specific humidity',              'kg/kg'   , fill_halo=.true. )
    call FILE_HISTORY_in( MSLP      (:,:),      'MSLP',       'mean sea-level pressure',           'Pa'      )

    return
  end subroutine history_output

end module mod_atmos_phy_sf_driver
