!-------------------------------------------------------------------------------
!> module ATMOSPHERE / Physics Radiation
!!
!! @par Description
!!          Atmospheric radiation transfer process driver
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2013-12-06 (S.Nishizawa)  [new]
!<
!-------------------------------------------------------------------------------
#include "inc_openmp.h"
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
  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine ATMOS_PHY_RD_driver_setup
    use scale_atmos_phy_rd, only: &
       ATMOS_PHY_RD_setup
    use mod_atmos_admin, only: &
       ATMOS_PHY_RD_TYPE, &
       ATMOS_sw_phy_rd
    use mod_atmos_phy_rd_vars, only: &
       SFCFLX_LW_up => ATMOS_PHY_RD_SFLX_LW_up,   &
       SFCFLX_LW_dn => ATMOS_PHY_RD_SFLX_LW_dn,   &
       SFCFLX_SW_up => ATMOS_PHY_RD_SFLX_SW_up,   &
       SFCFLX_SW_dn => ATMOS_PHY_RD_SFLX_SW_dn,   &
       TOAFLX_LW_up => ATMOS_PHY_RD_TOAFLX_LW_up, &
       TOAFLX_LW_dn => ATMOS_PHY_RD_TOAFLX_LW_dn, &
       TOAFLX_SW_up => ATMOS_PHY_RD_TOAFLX_SW_up, &
       TOAFLX_SW_dn => ATMOS_PHY_RD_TOAFLX_SW_dn
    implicit none
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '++++++ Module[DRIVER] / Categ[ATMOS PHY_RD] / Origin[SCALE-LES]'

    if ( ATMOS_sw_phy_rd ) then

       ! setup library component
       call ATMOS_PHY_RD_setup( ATMOS_PHY_RD_TYPE )

       ! run once (only for the diagnostic value)
       call ATMOS_PHY_RD_driver( .true., .false. )

    else

       if( IO_L ) write(IO_FID_LOG,*) '*** this component is never called.'
       if( IO_L ) write(IO_FID_LOG,*) '*** radiation fluxes are set to zero.'
       SFCFLX_LW_up(:,:) = 0.0_RP
       SFCFLX_LW_dn(:,:) = 0.0_RP
       SFCFLX_SW_up(:,:) = 0.0_RP
       SFCFLX_SW_dn(:,:) = 0.0_RP
       TOAFLX_LW_up(:,:) = 0.0_RP
       TOAFLX_LW_dn(:,:) = 0.0_RP
       TOAFLX_SW_up(:,:) = 0.0_RP
       TOAFLX_SW_dn(:,:) = 0.0_RP

    endif

    return
  end subroutine ATMOS_PHY_RD_driver_setup

  !-----------------------------------------------------------------------------
  !> Driver
  subroutine ATMOS_PHY_RD_driver( update_flag, history_flag )
    use scale_grid_real, only: &
       REAL_CZ,            &
       REAL_FZ,            &
       REAL_LON,           &
       REAL_LAT,           &
       REAL_BASEPOINT_LON, &
       REAL_BASEPOINT_LAT
    use scale_const, only: &
       PRE00 => CONST_PRE00, &
       Rdry  => CONST_Rdry,  &
       CPdry => CONST_CPdry
    use scale_landuse, only: &
       LANDUSE_fact_ocean, &
       LANDUSE_fact_land,  &
       LANDUSE_fact_urban
    use scale_time, only: &
       dt_RD => TIME_DTSEC_ATMOS_PHY_RD, &
       TIME_NOWDATE
    use scale_statistics, only: &
       STATISTICS_checktotal, &
       STAT_total
    use scale_history, only: &
       HIST_in
    use mod_atmos_admin, only: &
       ATMOS_PHY_RD_TYPE
    use scale_atmos_solarins, only: &
       SOLARINS_fixedlatlon => ATMOS_SOLARINS_fixedlatlon, &
       SOLARINS_insolation  => ATMOS_SOLARINS_insolation
    use scale_atmos_phy_rd, only: &
       ATMOS_PHY_RD
    use scale_atmos_phy_rd_mm5sw, only: &
       SWRAD
    use scale_atmos_phy_rd_common, only: &
       RD_heating => ATMOS_PHY_RD_heating, &
       I_SW, &
       I_LW, &
       I_dn, &
       I_up
    use mod_atmos_vars, only: &
       TEMP,              &
       PRES,              &
       DENS,              &
       RHOT,              &
       QTRC,              &
       RHOT_t => RHOT_tp
    use mod_atmos_phy_sf_vars, only: &
       SFC_TEMP   => ATMOS_PHY_SF_SFC_TEMP,  &
       SFC_albedo => ATMOS_PHY_SF_SFC_albedo
    use mod_atmos_phy_rd_vars, only: &
       RHOT_t_RD    => ATMOS_PHY_RD_RHOT_t,       &
       SFCFLX_LW_up => ATMOS_PHY_RD_SFLX_LW_up,   &
       SFCFLX_LW_dn => ATMOS_PHY_RD_SFLX_LW_dn,   &
       SFCFLX_SW_up => ATMOS_PHY_RD_SFLX_SW_up,   &
       SFCFLX_SW_dn => ATMOS_PHY_RD_SFLX_SW_dn,   &
       TOAFLX_LW_up => ATMOS_PHY_RD_TOAFLX_LW_up, &
       TOAFLX_LW_dn => ATMOS_PHY_RD_TOAFLX_LW_dn, &
       TOAFLX_SW_up => ATMOS_PHY_RD_TOAFLX_SW_up, &
       TOAFLX_SW_dn => ATMOS_PHY_RD_TOAFLX_SW_dn
    implicit none

    logical, intent(in) :: update_flag
    logical, intent(in) :: history_flag

    real(RP) :: TEMP_t      (KA,IA,JA,3)
    real(RP) :: flux_rad    (KA,IA,JA,2,2)
    real(RP) :: flux_rad_top(   IA,JA,2,2)

    real(RP) :: flux_up     (KA,IA,JA,2)
    real(RP) :: flux_dn     (KA,IA,JA,2)
    real(RP) :: flux_net    (KA,IA,JA,2)
    real(RP) :: flux_net_toa(   IA,JA,2)
    real(RP) :: flux_net_sfc(   IA,JA,2)

    real(RP) :: solins(IA,JA)
    real(RP) :: cosSZA(IA,JA)
    real(RP) :: LON   (IA,JA)
    real(RP) :: LAT   (IA,JA)

    ! for WRF radiation scheme added by Adachi; array order is (i,k,j)
    real(RP) :: RTHRATENSW(IA,KA,JA)
    real(RP) :: SDOWN3D   (IA,KA,JA)  ! downward short wave flux (W/m2)
    real(RP) :: GSW       (IA,JA)     ! net short wave flux at ground surface (W/m2)
    real(RP) :: RHO3D     (IA,KA,JA)
    real(RP) :: T3D       (IA,KA,JA)
    real(RP) :: P3D       (IA,KA,JA)
    real(RP) :: PI3D      (IA,KA,JA)
    real(RP) :: DZ8W      (IA,KA,JA)
    real(RP) :: QV3D      (IA,KA,JA)
    real(RP) :: QC3D      (IA,KA,JA)
    real(RP) :: QR3D      (IA,KA,JA)
    real(RP) :: QI3D      (IA,KA,JA)
    real(RP) :: QS3D      (IA,KA,JA)
    real(RP) :: QG3D      (IA,KA,JA)
    real(RP) :: fact1, fact2
    real(RP) :: flux_rad_org(KA,IA,JA,2,2)

    real(RP) :: total ! dummy

    integer :: k, i, j
    !---------------------------------------------------------------------------

    if ( update_flag ) then

       ! calc solar insolation
       if ( SOLARINS_fixedlatlon ) then
          LON(:,:) = REAL_BASEPOINT_LON
          LAT(:,:) = REAL_BASEPOINT_LAT
       else
          LON(:,:) = REAL_LON(:,:)
          LAT(:,:) = REAL_LAT(:,:)
       endif

       call SOLARINS_insolation( solins(:,:),    & ! [OUT]
                                 cosSZA(:,:),    & ! [OUT]
                                 LON   (:,:),    & ! [IN]
                                 LAT   (:,:),    & ! [IN]
                                 TIME_NOWDATE(:) ) ! [IN]

       call ATMOS_PHY_RD( DENS, RHOT, QTRC,   & ! [IN]
                          REAL_CZ, REAL_FZ,   & ! [IN]
                          LANDUSE_fact_ocean, & ! [IN]
                          LANDUSE_fact_land,  & ! [IN]
                          LANDUSE_fact_urban, & ! [IN]
                          SFC_TEMP,           & ! [IN]
                          SFC_albedo,         & ! [IN]
                          solins, cosSZA,     & ! [IN]
                          flux_rad,           & ! [OUT]
                          flux_rad_top        ) ! [OUT]

       ! apply radiative flux convergence -> heating rate
       call RD_heating( flux_rad (:,:,:,:,:), & ! [IN]
                        RHOT     (:,:,:),     & ! [IN]
                        QTRC     (:,:,:,:),   & ! [IN]
                        REAL_FZ  (:,:,:),     & ! [IN]
                        dt_RD,                & ! [IN]
                        TEMP_t   (:,:,:,:),   & ! [OUT]
                        RHOT_t_RD(:,:,:)      ) ! [OUT]


       do j = JS, JE
       do i = IS, IE
          SFCFLX_LW_up(i,j)      = flux_rad(KS-1,i,j,I_LW,I_up)
          SFCFLX_LW_dn(i,j)      = flux_rad(KS-1,i,j,I_LW,I_dn)
          SFCFLX_SW_up(i,j)      = flux_rad(KS-1,i,j,I_SW,I_up)
          SFCFLX_SW_dn(i,j)      = flux_rad(KS-1,i,j,I_SW,I_dn)

          flux_net_sfc(i,j,I_LW) = flux_rad(KS-1,i,j,I_LW,I_up) - flux_rad(KS-1,i,j,I_LW,I_dn)
          flux_net_sfc(i,j,I_SW) = flux_rad(KS-1,i,j,I_SW,I_up) - flux_rad(KS-1,i,j,I_SW,I_dn)
       enddo
       enddo

       do j = JS, JE
       do i = IS, IE
          TOAFLX_LW_up(i,j)      = flux_rad_top(i,j,I_LW,I_up)
          TOAFLX_LW_dn(i,j)      = flux_rad_top(i,j,I_LW,I_dn)
          TOAFLX_SW_up(i,j)      = flux_rad_top(i,j,I_SW,I_up)
          TOAFLX_SW_dn(i,j)      = flux_rad_top(i,j,I_SW,I_dn)

          flux_net_toa(i,j,I_LW) = flux_rad_top(i,j,I_LW,I_up) - flux_rad_top(i,j,I_LW,I_dn)
          flux_net_toa(i,j,I_SW) = flux_rad_top(i,j,I_SW,I_up) - flux_rad_top(i,j,I_SW,I_dn)
       enddo
       enddo

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

       if ( ATMOS_PHY_RD_TYPE == 'WRF' ) then

          if ( history_flag ) then

             flux_rad_org(:,:,:,:,:) = flux_rad(:,:,:,:,:)
             RTHRATENSW = 0.0_RP
             SDOWN3D    = 0.0_RP
             GSW        = 0.0_RP

             do j = 1, JA
             do i = 1, IA
             do k = 1, KA
                QV3D (i,k,j) = QTRC(k,i,j,I_QV)
                QC3D (i,k,j) = QTRC(k,i,j,I_QC)
                QR3D (i,k,j) = QTRC(k,i,j,I_QR)
                QI3D (i,k,j) = QTRC(k,i,j,I_QI)
                QS3D (i,k,j) = QTRC(k,i,j,I_QS)
                QG3D (i,k,j) = QTRC(k,i,j,I_QG)
                T3D  (i,k,j) = TEMP(k,i,j)                       ! temperature (K)
                RHO3D(i,k,j) = DENS(k,i,j)                       ! density (kg/m^3)
                P3D  (i,k,j) = PRES(k,i,j)                       ! pressure (Pa)
                PI3D (i,k,j) = (PRES(k,i,j)/PRE00)**(Rdry/CPdry) ! exner function (dimensionless)
                DZ8W (i,k,j) = REAL_FZ(k,i,j)-REAL_FZ(k-1,i,j)   ! dz between full levels(m)
             enddo
             enddo
             enddo

             call SWRAD( dt_RD,                & ! [IN]
                         RTHRATENSW,           & ! [INOUT]
                         SDOWN3D,              & ! [INOUT]
                         GSW,                  & ! [INOUT]
                         LAT,                  & ! [IN]
                         LON,                  & ! [IN]
                         SFC_albedo(:,:,I_SW), & ! [IN]
                         RHO3D,                & ! [IN]
                         T3D,                  & ! [IN]
                         P3D,                  & ! [IN]
                         PI3D,                 & ! [IN]
                         DZ8W,                 & ! [IN]
                         solins(:,:),          & ! [IN]
                         cosSZA(:,:),          & ! [IN]
                         QV3D,                 & ! [IN]
                         QC3D,                 & ! [IN]
                         QR3D,                 & ! [IN]
                         QI3D,                 & ! [IN]
                         QS3D,                 & ! [IN]
                         QG3D,                 & ! [IN]
                         F_QV      = .true.,   & ! [IN]
                         F_QC      = .true.,   & ! [IN]
                         F_QR      = .true.,   & ! [IN]
                         F_QI      = .true.,   & ! [IN]
                         F_QS      = .true.,   & ! [IN]
                         F_QG      = .true.,   & ! [IN]
                         icloud    = 1,        & ! [IN]
                         warm_rain = .true.    ) ! [IN]

             do j = JS, JE
             do i = IS, IE
                flux_net_sfc(i,j,I_SW) = GSW(i,j)
                do k = KS-1, KE
                   flux_rad(k,i,j,I_SW,I_up) = 0.0_RP
                   flux_rad(k,i,j,I_SW,I_dn) = SDOWN3D(i,k,j)
                enddo

                do k = 1, KS-2
                   flux_rad(k,i,j,I_SW,I_dn) = SDOWN3D(i,KS-1,j)
                enddo

                do k = KE+1, KA
                   flux_rad(k,i,j,I_SW,I_dn) = SDOWN3D(i,KE,j)
                enddo
             enddo
             enddo

             do j = JS, JE
             do i = IS, IE
                SFCFLX_SW_up(i,j) = flux_rad(KS-1,i,j,I_SW,I_dn) * SFC_albedo(i,j,I_SW)
                SFCFLX_SW_dn(i,j) = flux_rad(KS-1,i,j,I_SW,I_dn)

                TOAFLX_SW_up(i,j) = 0.0_RP
                TOAFLX_SW_dn(i,j) = flux_rad(KE,i,j,I_SW,I_dn) ! sometimes TOA altitude is very low
             enddo
             enddo

             do j = JS, JE
             do i = IS, IE
                flux_net_toa(i,j,I_SW) = 0.0_RP
                do k = KS, KE
                   flux_net(k,i,j,I_SW) = 0.0_RP
                   flux_up (k,i,j,I_SW) = 0.0_RP
                   flux_dn (k,i,j,I_SW) = 0.5_RP * ( flux_rad(k-1,i,j,I_SW,I_dn) + flux_rad(k,i,j,I_SW,I_dn) )
                enddo
             enddo
             enddo

             do j = JS, JE
             do i = IS, IE
                flux_net_sfc(i,j,I_SW) = flux_rad(KS-1,i,j,I_SW,I_dn)*SFC_albedo(i,j,I_SW)-flux_rad(KS-1,i,j,I_SW,I_dn)
             enddo
             enddo
          endif

       endif


       if ( history_flag ) then

          call HIST_in( solins(:,:), 'SOLINS', 'solar insolation',        'W/m2', dt_RD, nohalo=.true. )
          call HIST_in( cosSZA(:,:), 'COSZ',   'cos(solar zenith angle)', '0-1',  dt_RD , nohalo=.true.)

          call HIST_in( SFCFLX_LW_up(:,:), 'SFLX_LW_up',   'SFC upward   longwave  radiation flux', 'W/m2', dt_RD, nohalo=.true. )
          call HIST_in( SFCFLX_LW_dn(:,:), 'SFLX_LW_dn',   'SFC downward longwave  radiation flux', 'W/m2', dt_RD, nohalo=.true. )
          call HIST_in( SFCFLX_SW_up(:,:), 'SFLX_SW_up',   'SFC upward   shortwave radiation flux', 'W/m2', dt_RD, nohalo=.true. )
          call HIST_in( SFCFLX_SW_dn(:,:), 'SFLX_SW_dn',   'SFC downward shortwave radiation flux', 'W/m2', dt_RD, nohalo=.true. )

          call HIST_in( TOAFLX_LW_up(:,:), 'TOAFLX_LW_up', 'TOA upward   longwave  radiation flux', 'W/m2', dt_RD, nohalo=.true. )
          call HIST_in( TOAFLX_LW_dn(:,:), 'TOAFLX_LW_dn', 'TOA downward longwave  radiation flux', 'W/m2', dt_RD, nohalo=.true. )
          call HIST_in( TOAFLX_SW_up(:,:), 'TOAFLX_SW_up', 'TOA upward   shortwave radiation flux', 'W/m2', dt_RD, nohalo=.true. )
          call HIST_in( TOAFLX_SW_dn(:,:), 'TOAFLX_SW_dn', 'TOA downward shortwave radiation flux', 'W/m2', dt_RD, nohalo=.true. )

          call HIST_in( flux_net_sfc(:,:,I_LW), 'SLR', 'SFC net longwave  radiation flux', 'W/m2', dt_RD, nohalo=.true. )
          call HIST_in( flux_net_sfc(:,:,I_SW), 'SSR', 'SFC net shortwave radiation flux', 'W/m2', dt_RD, nohalo=.true. )
          call HIST_in( flux_net_toa(:,:,I_LW), 'OLR', 'TOA net longwave  radiation flux', 'W/m2', dt_RD, nohalo=.true. )
          call HIST_in( flux_net_toa(:,:,I_SW), 'OSR', 'TOA net shortwave radiation flux', 'W/m2', dt_RD, nohalo=.true. )

          call HIST_in( flux_up (:,:,:,I_LW), 'RADFLUX_LWUP', 'upward   longwave  radiation flux', 'W/m2', dt_RD, nohalo=.true. )
          call HIST_in( flux_dn (:,:,:,I_LW), 'RADFLUX_LWDN', 'downward longwave  radiation flux', 'W/m2', dt_RD, nohalo=.true. )
          call HIST_in( flux_net(:,:,:,I_LW), 'RADFLUX_LW',   'net      longwave  radiation flux', 'W/m2', dt_RD, nohalo=.true. )
          call HIST_in( flux_up (:,:,:,I_SW), 'RADFLUX_SWUP', 'upward   shortwave radiation flux', 'W/m2', dt_RD, nohalo=.true. )
          call HIST_in( flux_dn (:,:,:,I_SW), 'RADFLUX_SWDN', 'downward shortwave radiation flux', 'W/m2', dt_RD, nohalo=.true. )
          call HIST_in( flux_net(:,:,:,I_SW), 'RADFLUX_SW',   'net      shortwave radiation flux', 'W/m2', dt_RD, nohalo=.true. )

          call HIST_in( TEMP_t(:,:,:,I_LW), 'TEMP_t_rd_LW', 'tendency of temp in rd(LW)', 'K/day', dt_RD, nohalo=.true. )
          call HIST_in( TEMP_t(:,:,:,I_SW), 'TEMP_t_rd_SW', 'tendency of temp in rd(SW)', 'K/day', dt_RD, nohalo=.true. )
          call HIST_in( TEMP_t(:,:,:,3   ), 'TEMP_t_rd',    'tendency of temp in rd',     'K/day', dt_RD, nohalo=.true. )
       endif

       if ( ATMOS_PHY_RD_TYPE ==  'WRF' ) then
          ! revert all radiation flux from MM5 scheme to default
          flux_rad = flux_rad_org

          do j = JS, JE
          do i = IS, IE
             SFCFLX_SW_up(i,j) = flux_rad(KS-1,i,j,I_SW,I_up)
             SFCFLX_SW_dn(i,j) = flux_rad(KS-1,i,j,I_SW,I_dn)

             TOAFLX_SW_up(i,j) = flux_rad_top(i,j,I_SW,I_up) ! mstrnx
             TOAFLX_SW_dn(i,j) = flux_rad_top(i,j,I_SW,I_dn) ! mstrnx
          enddo
          enddo

       endif

    else
       ! to avoid adding UNDEF to RHOT_t in the next block
       RHOT_t_RD(:,:,:) = 0.0_RP

    endif

    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       RHOT_t(k,i,j) = RHOT_t(k,i,j) + RHOT_t_RD(k,i,j)
    enddo
    enddo
    enddo

    if ( STATISTICS_checktotal ) then
       call STAT_total( total, RHOT_t_RD(:,:,:), 'RHOT_t_RD' )
    endif

    return
  end subroutine ATMOS_PHY_RD_driver

end module mod_atmos_phy_rd_driver
