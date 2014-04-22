!-------------------------------------------------------------------------------
!> module ATMOSPHERE / Physics Radiation
!!
!! @par Description
!!          Atmospheric radiation transfer process wrapper
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2013-12-06 (S.Nishizawa)   [new]
!!
!<
!-------------------------------------------------------------------------------
module mod_atmos_phy_rd_driver
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
  public :: ATMOS_PHY_RD_driver_setup
  public :: ATMOS_PHY_RD_driver

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
  real(RP), private, allocatable :: RHOT_t(:,:,:)

  integer, private :: ktop
  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine ATMOS_PHY_RD_driver_setup( RD_TYPE )
    use scale_process, only: &
       PRC_MPIstop
    use scale_atmos_phy_rd, only: &
       ATMOS_PHY_RD_setup
    implicit none

    character(len=H_SHORT), intent(in) :: RD_TYPE
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[Physics-RD]/Categ[ATMOS]'

    allocate( RHOT_t(KA,IA,JA) )

    call ATMOS_PHY_RD_setup( RD_TYPE )

    call ATMOS_PHY_RD_driver( .true., .false. )

    return
  end subroutine ATMOS_PHY_RD_driver_setup

  !-----------------------------------------------------------------------------
  !> Radiation main
  subroutine ATMOS_PHY_RD_driver( update_flag, history_flag )
    use scale_grid_real, only: &
       REAL_CZ,  &
       REAL_FZ,  &
       REAL_LON, &
       REAL_LAT
    use scale_landuse, only: &
       LANDUSE_frac_ocean
    use scale_time, only: &
       dt_RD => TIME_DTSEC_ATMOS_PHY_RD, &
       TIME_NOWDATE
    use scale_history, only: &
       HIST_in
    use scale_atmos_thermodyn, only: &
       THERMODYN_qd        => ATMOS_THERMODYN_qd,       &
       THERMODYN_cv        => ATMOS_THERMODYN_cv,       &
       THERMODYN_rhoe      => ATMOS_THERMODYN_rhoe,     &
       THERMODYN_rhot      => ATMOS_THERMODYN_rhot,     &
       THERMODYN_temp_pres => ATMOS_THERMODYN_temp_pres
    use scale_atmos_solarins, only: &
       SOLARINS_insolation => ATMOS_SOLARINS_insolation
    use scale_atmos_phy_rd, only: &
       ATMOS_PHY_RD
    use scale_atmos_phy_rd_common, only: &
       RD_heating => ATMOS_PHY_RD_heating, &
       I_SW, &
       I_LW, &
       I_dn, &
       I_up
    use mod_atmos_vars, only: &
       ATMOS_vars_fillhalo, &
       ATMOS_vars_total, &
       DENS, &
       RHOT, &
       QTRC, &
       RHOT_tp
    use mod_atmos_vars_sf, only: &
       ATMOS_vars_sf_fillhalo, &
       SWD,                    &
       LWD
    use mod_cpl_vars, only: &
       sw_AtmLnd => CPL_sw_AtmLnd, &
       sw_AtmOcn => CPL_sw_AtmOcn, &
       SkinT,                      &
       ALBW,                       &
       ALBG
    implicit none

    logical, intent(in) :: update_flag
    logical, intent(in), optional :: history_flag

    real(RP) :: RHOE    (KA,IA,JA)
    real(RP) :: RHOE_t  (KA,IA,JA,2)
    real(RP) :: RHOT_tmp(KA,IA,JA)
    real(RP) :: TEMP_t  (KA,IA,JA)
    real(RP) :: QDRY    (KA,IA,JA)
    real(RP) :: CVtot   (KA,IA,JA)
    real(RP) :: TEMP    (KA,IA,JA)
    real(RP) :: PRES    (KA,IA,JA)

    real(RP) :: flux_rad(KA,IA,JA,2,2)
    real(RP) :: flux_net(KA,IA,JA,2)
    real(RP) :: flux_up (KA,IA,JA,2)
    real(RP) :: flux_dn (KA,IA,JA,2)

    real(RP) :: flux_rad_top(IA,JA,2)
    real(RP) :: flux_rad_sfc(IA,JA,2)

    real(RP) :: solins  (IA,JA)
    real(RP) :: cosSZA  (IA,JA)
    real(RP) :: temp_sfc(IA,JA)
    real(RP) :: albedo_land(IA,JA,2)

    integer :: k, i, j
    !---------------------------------------------------------------------------

    if ( update_flag ) then

       if ( sw_AtmLnd .or. sw_AtmOcn ) then ! tentative
          temp_sfc   (:,:) = SkinT(:,:)
       else
          call THERMODYN_temp_pres( TEMP(:,:,:),  & ! [OUT]
                                    PRES(:,:,:),  & ! [OUT]
                                    DENS(:,:,:),  & ! [IN]
                                    RHOT(:,:,:),  & ! [IN]
                                    QTRC(:,:,:,:) ) ! [IN]

          temp_sfc(:,:) = TEMP(KS,:,:)
       endif

       if ( sw_AtmLnd ) then ! tentative
          albedo_land(:,:,I_SW) = ALBG(:,:,1)
          albedo_land(:,:,I_LW) = ALBG(:,:,2)
       elseif ( sw_AtmOcn ) then ! tentative
          albedo_land(:,:,I_SW) = ALBW(:,:,1)
          albedo_land(:,:,I_LW) = ALBW(:,:,2)
       else
          albedo_land(:,:,:) = 0.5_RP
       endif

       call SOLARINS_insolation( solins  (:,:),  & ! [OUT]
                                 cosSZA  (:,:),  & ! [OUT]
                                 REAL_LON(:,:),  & ! [IN]
                                 REAL_LAT(:,:),  & ! [IN]
                                 TIME_NOWDATE(:) ) ! [IN]

       call ATMOS_PHY_RD( DENS, RHOT, QTRC,      & ! [IN]
                          REAL_CZ, REAL_FZ,      & ! [IN]
                          LANDUSE_frac_ocean,    & ! [IN]
                          temp_sfc, albedo_land, & ! [IN]
                          solins, cosSZA,        & ! [IN]
                          flux_rad,              & ! [OUT]
                          flux_rad_top           ) ! [OUT]

       call THERMODYN_rhoe( RHOE(:,:,:),  & ! [OUT]
                            RHOT(:,:,:),  & ! [IN]
                            QTRC(:,:,:,:) ) ! [IN]

       ! apply radiative flux convergence -> heating rate
       call RD_heating( flux_rad(:,:,:,:,:),  & ! [IN]
                        REAL_FZ (:,:,:),      & ! [IN]
                        dt_RD,                & ! [IN]
                        RHOE_t  (:,:,:,:),    & ! [OUT]
                        RHOE    (:,:,:)       ) ! [INOUT]

       ! update rhot
       call THERMODYN_rhot( RHOT_tmp(:,:,:),& ! [OUT]
                            RHOE(:,:,:),    & ! [IN]
                            QTRC(:,:,:,:)   ) ! [IN]

       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          RHOT_t(k,i,j) = ( RHOT_tmp(k,i,j)-RHOT(k,i,j) ) / dt_RD
       enddo
       enddo
       enddo

       if ( sw_AtmLnd .or. sw_AtmOcn ) then
          do j = JS, JE
          do i = IS, IE
             LWD(i,j) = flux_rad(KS-1,i,j,I_LW,I_dn)
             SWD(i,j) = flux_rad(KS-1,i,j,I_SW,I_dn)
          enddo
          enddo
       endif

       if ( present(history_flag) ) then
       if ( history_flag ) then

          call HIST_in( solins(:,:), 'SOLINS', 'solar insolation',        'W/m2', dt_RD )
          call HIST_in( cosSZA(:,:), 'COSZ',   'cos(solar zenith angle)', '0-1',  dt_RD )

          do j = JS, JE
          do i = IS, IE
          do k = KS, KE
             flux_net(k,i,j,I_LW) = 0.5_RP * ( ( flux_rad(k-1,i,j,I_LW,I_up) - flux_rad(k-1,i,j,I_LW,I_dn) ) &
                                             + ( flux_rad(k  ,i,j,I_LW,I_up) - flux_rad(k  ,i,j,I_LW,I_dn) ) )
             flux_net(k,i,j,I_SW) = 0.5_RP * ( ( flux_rad(k-1,i,j,I_SW,I_up) - flux_rad(k-1,i,j,I_SW,I_dn) ) &
                                             + ( flux_rad(k  ,i,j,I_SW,I_up) - flux_rad(k  ,i,j,I_SW,I_dn) ) )

             flux_up (k,i,j,I_LW) = 0.5_RP * ( flux_rad(k-1,i,j,I_LW,I_up) + flux_rad(k,i,j,I_LW,I_up) )
             flux_up (k,i,j,I_SW) = 0.5_RP * ( flux_rad(k-1,i,j,I_SW,I_up) + flux_rad(k,i,j,I_SW,I_up) )
             flux_dn (k,i,j,I_LW) = 0.5_RP * ( flux_rad(k-1,i,j,I_LW,I_dn) + flux_rad(k,i,j,I_LW,I_dn) )
             flux_dn (k,i,j,I_SW) = 0.5_RP * ( flux_rad(k-1,i,j,I_SW,I_dn) + flux_rad(k,i,j,I_SW,I_dn) )
          enddo
          enddo
          enddo
          call HIST_in( flux_net(:,:,:,I_LW), 'RADFLUX_LW',  'net radiation flux(LW)', 'W/m2', dt_RD )
          call HIST_in( flux_net(:,:,:,I_SW), 'RADFLUX_SW',  'net radiation flux(SW)', 'W/m2', dt_RD )
          call HIST_in( flux_up (:,:,:,I_LW), 'RADFLUX_LWUP', 'up radiation flux(LW)', 'W/m2', dt_RD )
          call HIST_in( flux_up (:,:,:,I_SW), 'RADFLUX_SWUP', 'up radiation flux(SW)', 'W/m2', dt_RD )
          call HIST_in( flux_dn (:,:,:,I_LW), 'RADFLUX_LWDN', 'dn radiation flux(LW)', 'W/m2', dt_RD )
          call HIST_in( flux_dn (:,:,:,I_SW), 'RADFLUX_SWDN', 'dn radiation flux(SW)', 'W/m2', dt_RD )

          call HIST_in( flux_rad_top(:,:,I_LW), 'OLR', 'TOA     longwave  radiation', 'W/m2', dt_RD )
          call HIST_in( flux_rad_top(:,:,I_SW), 'OSR', 'TOA     shortwave radiation', 'W/m2', dt_RD )
          do j = JS, JE
          do i = IS, IE
             flux_rad_sfc(i,j,I_LW) = flux_rad(KS-1,i,j,I_LW,I_up)-flux_rad(KS-1,i,j,I_LW,I_dn)
             flux_rad_sfc(i,j,I_SW) = flux_rad(KS-1,i,j,I_SW,I_up)-flux_rad(KS-1,i,j,I_SW,I_dn)
          enddo
          enddo
          call HIST_in( flux_rad_sfc(:,:,I_LW), 'SLR', 'Surface longwave  radiation', 'W/m2', dt_RD )
          call HIST_in( flux_rad_sfc(:,:,I_SW), 'SSR', 'Surface shortwave radiation', 'W/m2', dt_RD )

          call THERMODYN_qd( QDRY(:,:,:),  & ! [OUT]
                             QTRC(:,:,:,:) ) ! [IN]

          call THERMODYN_cv( CVtot(:,:,:),   & ! [OUT]
                             QTRC (:,:,:,:), & ! [IN]
                             QDRY (:,:,:)    ) ! [IN]

          TEMP_t(:,:,:) = RHOE_t(:,:,:,I_LW) / CVtot(:,:,:) * 86400.0_RP
          call HIST_in( TEMP_t, 'TEMP_t_rd_LW', 'tendency of temp in rd(LW)', 'K/day', dt_RD )
          TEMP_t(:,:,:) = RHOE_t(:,:,:,I_SW) / CVtot(:,:,:) * 86400.0_RP
          call HIST_in( TEMP_t, 'TEMP_t_rd_SW', 'tendency of temp in rd(SW)', 'K/day', dt_RD )
          TEMP_t(:,:,:) = ( RHOE_t(:,:,:,I_LW) + RHOE_t(:,:,:,I_SW) ) / CVtot(:,:,:) * 86400.0_RP
          call HIST_in( TEMP_t, 'TEMP_t_rd', 'tendency of temp in rd', 'K/day', dt_RD )
       endif
       endif
    endif

    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       RHOT_tp(k,i,j) = RHOT_tp(k,i,j) + RHOT_t(k,i,j)
    enddo
    enddo
    enddo

    return
  end subroutine ATMOS_PHY_RD_driver

end module mod_atmos_phy_rd_driver
