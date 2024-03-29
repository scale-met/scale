!-------------------------------------------------------------------------------
!> module ATMOSPHERE / Physics Radiation
!!
!! @par Description
!!          Atmospheric radiation transfer process driver
!!
!! @author Team SCALE
!!
!<
!-------------------------------------------------------------------------------
#include "scalelib.h"
module mod_atmos_phy_rd_driver
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_io
  use scale_prof
  use scale_atmos_grid_icoA_index
  use scale_tracer
  use scale_cpl_sfc_index
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: ATMOS_PHY_RD_driver_setup
  public :: ATMOS_PHY_RD_driver_step

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
    use scale_prc, only: &
       PRC_abort
    use scale_atmos_phy_rd_mstrnx, only: &
       ATMOS_PHY_RD_mstrnx_setup
    use scale_atmos_phy_rd_offline, only: &
       ATMOS_PHY_RD_offline_setup
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
       TOAFLX_SW_dn => ATMOS_PHY_RD_TOAFLX_SW_dn, &
       SFLX_rad_dn  => ATMOS_PHY_RD_SFLX_down,    &
       solins       => ATMOS_PHY_RD_solins,       &
       cosSZA       => ATMOS_PHY_RD_cosSZA
    use mod_grd, only: &
       GRD_gz, &
       GRD_gzh
    implicit none

    real(RP) :: FZ(0:KA)
    integer :: k
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("ATMOS_PHY_RD_driver_setup",*) 'Setup'

    if ( ATMOS_sw_phy_rd ) then

       select case ( ATMOS_PHY_RD_TYPE )
       case ( "MSTRNX" )
          do k = KS-1, KE
             FZ(k) = GRD_gzh(k+1)
          end do
          call ATMOS_PHY_RD_MSTRNX_setup( KA, KS, KE, GRD_gz(:), FZ(:) )
       case ( "OFFLINE" )
          call ATMOS_PHY_RD_offline_setup
       case default
          LOG_ERROR("ATMOS_PHY_RD_driver_setup",*) 'invalid Radiation type(', trim(ATMOS_PHY_RD_TYPE), '). CHECK!'
          call PRC_abort
       end select

    else

       LOG_INFO("ATMOS_PHY_RD_driver_setup",*) 'this component is never called.'
       LOG_INFO("ATMOS_PHY_RD_driver_setup",*) 'radiation fluxes are set to zero.'
       SFCFLX_LW_up(:,:,:)     = 0.0_RP
       SFCFLX_LW_dn(:,:,:)     = 0.0_RP
       SFCFLX_SW_up(:,:,:)     = 0.0_RP
       SFCFLX_SW_dn(:,:,:)     = 0.0_RP
       TOAFLX_LW_up(:,:,:)     = 0.0_RP
       TOAFLX_LW_dn(:,:,:)     = 0.0_RP
       TOAFLX_SW_up(:,:,:)     = 0.0_RP
       TOAFLX_SW_dn(:,:,:)     = 0.0_RP
       SFLX_rad_dn (:,:,:,:,:) = 0.0_RP
       solins      (:,:,:)     = 0.0_RP
       cosSZA      (:,:,:)     = 0.0_RP

    endif

    return
  end subroutine ATMOS_PHY_RD_driver_setup

  !-----------------------------------------------------------------------------
  !> time step
  subroutine ATMOS_PHY_RD_driver_step
    use scale_time, only: &
       TIME_NOWDATE,                     &
       TIME_NOWDAYSEC,                   &
       TIME_OFFSET_YEAR
   !  use scale_file_history, only: &
   !     FILE_HISTORY_in
    use mod_history, only: &
       history_in
    use scale_atmos_hydrometeor, only: &
       N_HYD
    use scale_atmos_aerosol, only: &
       N_AE
    use mod_atmos_admin, only: &
       ATMOS_PHY_RD_TYPE
    use scale_atmos_solarins, only: &
       SOLARINS_insolation => ATMOS_SOLARINS_insolation
    use scale_atmos_phy_rd_mstrnx, only: &
       ATMOS_PHY_RD_MSTRNX_flux
    use scale_atmos_phy_rd_offline, only: &
       ATMOS_PHY_RD_OFFLINE_flux
    use scale_atmos_phy_rd_common, only: &
       ATMOS_PHY_RD_calc_heating, &
       I_SW, &
       I_LW, &
       I_dn, &
       I_up
    use scale_time, only: &
       dt_RD => TIME_DTSEC_ATMOS_PHY_RD
    use mod_atmos_vars, only: &
       ATMOS_vars_calc_diagnostics, &
       DENS,  &
       RHOE,  &
       QTRC,  &
       TEMP,  &
       PRES,  &
       QV,    &
       CVtot, &
       CZ,    &
       FZ,    &
       LON,   &
       LAT,   &
       fact_ocean, &
       fact_land,  &
       fact_urban
    use mod_atmos_phy_sf_vars, only: &
       SFC_TEMP    => ATMOS_PHY_SF_SFC_TEMP,   &
       SFC_albedo  => ATMOS_PHY_SF_SFC_albedo
    use mod_atmos_phy_rd_vars, only: &
       SFCFLX_LW_up => ATMOS_PHY_RD_SFLX_LW_up,   &
       SFCFLX_LW_dn => ATMOS_PHY_RD_SFLX_LW_dn,   &
       SFCFLX_SW_up => ATMOS_PHY_RD_SFLX_SW_up,   &
       SFCFLX_SW_dn => ATMOS_PHY_RD_SFLX_SW_dn,   &
       TOAFLX_LW_up => ATMOS_PHY_RD_TOAFLX_LW_up, &
       TOAFLX_LW_dn => ATMOS_PHY_RD_TOAFLX_LW_dn, &
       TOAFLX_SW_up => ATMOS_PHY_RD_TOAFLX_SW_up, &
       TOAFLX_SW_dn => ATMOS_PHY_RD_TOAFLX_SW_dn, &
       SFLX_rad_dn  => ATMOS_PHY_RD_SFLX_down,    &
       solins       => ATMOS_PHY_RD_solins,       &
       cosSZA       => ATMOS_PHY_RD_cosSZA
    use mod_atmos_vars, only: &
       ATMOS_vars_get_diagnostic
    use mod_atmos_phy_mp_vars, only: &
       ATMOS_PHY_MP_vars_get_diagnostic
    use mod_atmos_phy_ae_vars, only: &
       ATMOS_PHY_AE_vars_get_diagnostic
    implicit none

    real(RP) :: RHOH        (KA,IA,JA,ADM_lall)
    real(RP) :: TEMP_t      (KA,IA,JA,3,ADM_lall)
    real(RP) :: flux_rad    (KA,IA,JA,2,2,2,ADM_lall)
    real(RP) :: flux_rad_top(   IA,JA,2,2,2,ADM_lall)

    real(RP) :: dtau_s      (KA,IA,JA,ADM_lall) ! 0.67 micron cloud optical depth
    real(RP) :: dem_s       (KA,IA,JA,ADM_lall) ! 10.5 micron cloud emissivity

    real(RP) :: flux_up     (KA,IA,JA,ADM_lall,2)
    real(RP) :: flux_dn     (KA,IA,JA,ADM_lall,2)
    real(RP) :: flux_net    (KA,IA,JA,ADM_lall,2)
    real(RP) :: flux_net_toa(   IA,JA,ADM_lall,2)
    real(RP) :: flux_net_sfc(   IA,JA,ADM_lall,2)

    real(RP) :: SFCFLX_LW_up_c(IA,JA,ADM_lall)
    real(RP) :: SFCFLX_LW_dn_c(IA,JA,ADM_lall)
    real(RP) :: SFCFLX_SW_up_c(IA,JA,ADM_lall)
    real(RP) :: SFCFLX_SW_dn_c(IA,JA,ADM_lall)
    real(RP) :: TOAFLX_LW_up_c(IA,JA,ADM_lall)
    real(RP) :: TOAFLX_LW_dn_c(IA,JA,ADM_lall)
    real(RP) :: TOAFLX_SW_up_c(IA,JA,ADM_lall)
    real(RP) :: TOAFLX_SW_dn_c(IA,JA,ADM_lall)

    real(RP) :: CLDFRAC(KA,IA,JA,ADM_lall)
    real(RP) :: MP_Re  (KA,IA,JA,ADM_lall,N_HYD)
    real(RP) :: MP_Qe  (KA,IA,JA,ADM_lall,N_HYD)
    real(RP) :: AE_Re  (KA,IA,JA,ADM_lall,N_AE)
    real(RP) :: AE_Qe  (KA,IA,JA,ADM_lall,N_AE)

    real(RP) :: RH(KA,IA,JA,ADM_lall)

    integer  :: k, i, j, l
    !---------------------------------------------------------------------------

    do l = 1, ADM_lall
       call SOLARINS_insolation( IA, IS, IE, JA, JS, JE, &
                                 LON(:,:,l), LAT(:,:,l),            & ! [IN]
                                 TIME_NOWDATE(:), TIME_OFFSET_YEAR, & ! [IN]
                                 solins(:,:,l), cosSZA  (:,:,l)     ) ! [OUT]
    end do

    call ATMOS_PHY_MP_vars_get_diagnostic( &
         DENS(:,:,:,:), TEMP(:,:,:,:), QTRC(:,:,:,:,:), & ! [IN]
         CLDFRAC=CLDFRAC, Re=MP_Re, Qe=MP_Qe            ) ! [IN]

    call ATMOS_vars_get_diagnostic( "RH", RH )
    call ATMOS_PHY_AE_vars_get_diagnostic( &
         QTRC(:,:,:,:,:), RH(:,:,:,:), & ! [IN]
         Re=AE_Re, Qe=AE_Qe            ) ! [IN]


    select case ( ATMOS_PHY_RD_TYPE )
    case ( "MSTRNX" )

       do l = 1, ADM_lall
          call ATMOS_PHY_RD_MSTRNX_flux( &
               KA, KS, KE, IA, IS, IE, JA, JS, JE, &
               DENS(:,:,:,l), TEMP(:,:,:,l), PRES(:,:,:,l), QV(:,:,:,l), & ! [IN]
               CZ(:,:,:,l), FZ(:,:,:,l),                   & ! [IN]
               fact_ocean(:,:,l), fact_land(:,:,l), fact_urban(:,:,l), & ! [IN]
               SFC_TEMP(:,:,l), SFC_albedo(:,:,:,:,l),               & ! [IN]
               solins(:,:,l), cosSZA(:,:,l),                         & ! [IN]
               CLDFRAC(:,:,:,l), MP_Re(:,:,:,l,:), MP_Qe(:,:,:,l,:),   & ! [IN]
               AE_Re(:,:,:,l,:), AE_Qe(:,:,:,l,:),                   & ! [IN]
               flux_rad(:,:,:,:,:,:,l),                            & ! [OUT]
               flux_rad_top(:,:,:,:,:,l), SFLX_rad_dn(:,:,:,:,l),    & ! [OUT]
               dtau_s = dtau_s(:,:,:,l), dem_s = dem_s(:,:,:,l)      ) ! [OUT]
       end do

    case ( "OFFLINE" )

       do l = 1, ADM_lall
          call ATMOS_PHY_RD_offline_flux( &
               KA, KS, KE, IA, IS, IE, JA, JS, JE, &
               TIME_NOWDAYSEC,        & ! [IN]
               flux_rad(:,:,:,:,:,2,l), & ! [OUT]
               SFLX_rad_dn(:,:,:,:,l)   ) ! [OUT]
          flux_rad(:,:,:,:,:,1,l)   = 0.0_RP ! clear sky
       end do
       flux_rad_top(:,:,:,:,:,:) = 0.0_RP
       dtau_s(:,:,:,:) = 0.0_RP
       dem_s(:,:,:,:) = 0.0_RP

    end select


    do l = 1, ADM_lall

!OCL XFILL
       do j = JS, JE
       do i = IS, IE
          ! for clear-sky
          SFCFLX_LW_up_c(i,j,l)    = flux_rad(KS-1,i,j,I_LW,I_up,1,l)
          SFCFLX_LW_dn_c(i,j,l)    = flux_rad(KS-1,i,j,I_LW,I_dn,1,l)
          SFCFLX_SW_up_c(i,j,l)    = flux_rad(KS-1,i,j,I_SW,I_up,1,l)
          SFCFLX_SW_dn_c(i,j,l)    = flux_rad(KS-1,i,j,I_SW,I_dn,1,l)
          ! for all-sky
          SFCFLX_LW_up  (i,j,l)    = flux_rad(KS-1,i,j,I_LW,I_up,2,l)
          SFCFLX_LW_dn  (i,j,l)    = flux_rad(KS-1,i,j,I_LW,I_dn,2,l)
          SFCFLX_SW_up  (i,j,l)    = flux_rad(KS-1,i,j,I_SW,I_up,2,l)
          SFCFLX_SW_dn  (i,j,l)    = flux_rad(KS-1,i,j,I_SW,I_dn,2,l)

          flux_net_sfc(i,j,l,I_LW) = SFCFLX_LW_up(i,j,l) - SFCFLX_LW_dn(i,j,l)
          flux_net_sfc(i,j,l,I_SW) = SFCFLX_SW_up(i,j,l) - SFCFLX_SW_dn(i,j,l)
       enddo
       enddo

!OCL XFILL
       do j = JS, JE
       do i = IS, IE
          ! for clear-sky
          TOAFLX_LW_up_c(i,j,l)    = flux_rad_top(i,j,I_LW,I_up,1,l)
          TOAFLX_LW_dn_c(i,j,l)    = flux_rad_top(i,j,I_LW,I_dn,1,l)
          TOAFLX_SW_up_c(i,j,l)    = flux_rad_top(i,j,I_SW,I_up,1,l)
          TOAFLX_SW_dn_c(i,j,l)    = flux_rad_top(i,j,I_SW,I_dn,1,l)
          ! for all-sky
          TOAFLX_LW_up  (i,j,l)    = flux_rad_top(i,j,I_LW,I_up,2,l)
          TOAFLX_LW_dn  (i,j,l)    = flux_rad_top(i,j,I_LW,I_dn,2,l)
          TOAFLX_SW_up  (i,j,l)    = flux_rad_top(i,j,I_SW,I_up,2,l)
          TOAFLX_SW_dn  (i,j,l)    = flux_rad_top(i,j,I_SW,I_dn,2,l)

          flux_net_toa(i,j,l,I_LW) = TOAFLX_LW_up(i,j,l) - TOAFLX_LW_dn(i,j,l)
          flux_net_toa(i,j,l,I_SW) = TOAFLX_SW_up(i,j,l) - TOAFLX_SW_dn(i,j,l)
       enddo
       enddo

!OCL XFILL
       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          flux_up (k,i,j,l,I_LW) = 0.5_RP * ( flux_rad(k-1,i,j,I_LW,I_up,2,l) + flux_rad(k,i,j,I_LW,I_up,2,l) )
          flux_dn (k,i,j,l,I_LW) = 0.5_RP * ( flux_rad(k-1,i,j,I_LW,I_dn,2,l) + flux_rad(k,i,j,I_LW,I_dn,2,l) )
          flux_up (k,i,j,l,I_SW) = 0.5_RP * ( flux_rad(k-1,i,j,I_SW,I_up,2,l) + flux_rad(k,i,j,I_SW,I_up,2,l) )
          flux_dn (k,i,j,l,I_SW) = 0.5_RP * ( flux_rad(k-1,i,j,I_SW,I_dn,2,l) + flux_rad(k,i,j,I_SW,I_dn,2,l) )

          flux_net(k,i,j,l,I_LW) = flux_up(k,i,j,l,I_LW) - flux_dn(k,i,j,l,I_LW)
          flux_net(k,i,j,l,I_SW) = flux_up(k,i,j,l,I_SW) - flux_dn(k,i,j,l,I_SW)
       enddo
       enddo
       enddo

       ! apply radiative flux convergence -> heating rate
       call ATMOS_PHY_RD_calc_heating( &
            KA, KS, KE, IA, IS, IE, JA, JS, JE, &
            flux_rad(:,:,:,:,:,2,l),      & ! [IN]
            DENS(:,:,:,l), TEMP(:,:,:,l), & ! [IN]
            CVtot(:,:,:,l),               & ! [IN]
            FZ(:,:,:,l),                  & ! [IN]
            RHOH(:,:,:,l),                & ! [OUT]
            temp_t = TEMP_t(:,:,:,:,l)    ) ! [OUT]

       !$omp parallel do
       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          RHOE(k,i,j,l) = RHOE(k,i,j,l) + RHOH(k,i,j,l) * dt_RD
       end do
       end do
       end do

    end do


   !  call FILE_HISTORY_in( solins(:,:,:), 'SOLINS', 'solar insolation',        'W/m2', fill_halo=.true. )
   !  call FILE_HISTORY_in( cosSZA(:,:,:), 'COSZ',   'cos(solar zenith angle)', '1',    fill_halo=.true. )

   !  call FILE_HISTORY_in( SFCFLX_LW_up_c(:,:,:),      'SFLX_LW_up_c',   'SFC upward   longwave  radiation flux (clr)', 'W/m2', fill_halo=.true. )
   !  call FILE_HISTORY_in( SFCFLX_LW_dn_c(:,:,:),      'SFLX_LW_dn_c',   'SFC downward longwave  radiation flux (clr)', 'W/m2', fill_halo=.true. )
   !  call FILE_HISTORY_in( SFCFLX_SW_up_c(:,:,:),      'SFLX_SW_up_c',   'SFC upward   shortwave radiation flux (clr)', 'W/m2', fill_halo=.true. )
   !  call FILE_HISTORY_in( SFCFLX_SW_dn_c(:,:,:),      'SFLX_SW_dn_c',   'SFC downward shortwave radiation flux (clr)', 'W/m2', fill_halo=.true. )

   !  call FILE_HISTORY_in( SFCFLX_LW_up  (:,:,:),      'SFLX_LW_up',     'SFC upward   longwave  radiation flux',       'W/m2', fill_halo=.true. )
   !  call FILE_HISTORY_in( SFCFLX_LW_dn  (:,:,:),      'SFLX_LW_dn',     'SFC downward longwave  radiation flux',       'W/m2', fill_halo=.true. )
   !  call FILE_HISTORY_in( SFCFLX_SW_up  (:,:,:),      'SFLX_SW_up',     'SFC upward   shortwave radiation flux',       'W/m2', fill_halo=.true. )
   !  call FILE_HISTORY_in( SFCFLX_SW_dn  (:,:,:),      'SFLX_SW_dn',     'SFC downward shortwave radiation flux',       'W/m2', fill_halo=.true. )

   !  call FILE_HISTORY_in( flux_net_sfc  (:,:,:,I_LW), 'SFLX_LW_net',    'SFC net      longwave  radiation flux',       'W/m2', fill_halo=.true. )
   !  call FILE_HISTORY_in( flux_net_sfc  (:,:,:,I_SW), 'SFLX_SW_net',    'SFC net      shortwave radiation flux',       'W/m2', fill_halo=.true. )

   !  call FILE_HISTORY_in( TOAFLX_LW_up_c(:,:,:),      'TOAFLX_LW_up_c', 'TOA upward   longwave  radiation flux (clr)', 'W/m2', fill_halo=.true. )
   !  call FILE_HISTORY_in( TOAFLX_LW_dn_c(:,:,:),      'TOAFLX_LW_dn_c', 'TOA downward longwave  radiation flux (clr)', 'W/m2', fill_halo=.true. )
   !  call FILE_HISTORY_in( TOAFLX_SW_up_c(:,:,:),      'TOAFLX_SW_up_c', 'TOA upward   shortwave radiation flux (clr)', 'W/m2', fill_halo=.true. )
   !  call FILE_HISTORY_in( TOAFLX_SW_dn_c(:,:,:),      'TOAFLX_SW_dn_c', 'TOA downward shortwave radiation flux (clr)', 'W/m2', fill_halo=.true. )

   !  call FILE_HISTORY_in( TOAFLX_LW_up  (:,:,:),      'TOAFLX_LW_up',   'TOA upward   longwave  radiation flux',       'W/m2', fill_halo=.true. )
   !  call FILE_HISTORY_in( TOAFLX_LW_dn  (:,:,:),      'TOAFLX_LW_dn',   'TOA downward longwave  radiation flux',       'W/m2', fill_halo=.true. )
   !  call FILE_HISTORY_in( TOAFLX_SW_up  (:,:,:),      'TOAFLX_SW_up',   'TOA upward   shortwave radiation flux',       'W/m2', fill_halo=.true. )
   !  call FILE_HISTORY_in( TOAFLX_SW_dn  (:,:,:),      'TOAFLX_SW_dn',   'TOA downward shortwave radiation flux',       'W/m2', fill_halo=.true. )

   !  call FILE_HISTORY_in( flux_net_toa  (:,:,:,I_LW), 'TOAFLX_LW_net',  'TOA net      longwave  radiation flux',       'W/m2', fill_halo=.true. )
   !  call FILE_HISTORY_in( flux_net_toa  (:,:,:,I_SW), 'TOAFLX_SW_net',  'TOA net      shortwave radiation flux',       'W/m2', fill_halo=.true. )

   !  call FILE_HISTORY_in( flux_net_sfc(:,:,:,I_LW),   'SLR',          'SFC net longwave  radiation flux',  'W/m2', fill_halo=.true. )
   !  call FILE_HISTORY_in( flux_net_sfc(:,:,:,I_SW),   'SSR',          'SFC net shortwave radiation flux',  'W/m2', fill_halo=.true. )
   !  call FILE_HISTORY_in( TOAFLX_LW_up(:,:,:),        'OLR',          'outgoing longwave  radiation flux', 'W/m2', fill_halo=.true. )
   !  call FILE_HISTORY_in( TOAFLX_SW_up(:,:,:),        'OSR',          'outgoing shortwave radiation flux', 'W/m2', fill_halo=.true. )

   !  call FILE_HISTORY_in( flux_up     (:,:,:,:,I_LW), 'RADFLUX_LWUP', 'upward   longwave  radiation flux', 'W/m2', fill_halo=.true. )
   !  call FILE_HISTORY_in( flux_dn     (:,:,:,:,I_LW), 'RADFLUX_LWDN', 'downward longwave  radiation flux', 'W/m2', fill_halo=.true. )
   !  call FILE_HISTORY_in( flux_net    (:,:,:,:,I_LW), 'RADFLUX_LW',   'net      longwave  radiation flux', 'W/m2', fill_halo=.true. )
   !  call FILE_HISTORY_in( flux_up     (:,:,:,:,I_SW), 'RADFLUX_SWUP', 'upward   shortwave radiation flux', 'W/m2', fill_halo=.true. )
   !  call FILE_HISTORY_in( flux_dn     (:,:,:,:,I_SW), 'RADFLUX_SWDN', 'downward shortwave radiation flux', 'W/m2', fill_halo=.true. )
   !  call FILE_HISTORY_in( flux_net    (:,:,:,:,I_SW), 'RADFLUX_SW',   'net      shortwave radiation flux', 'W/m2', fill_halo=.true. )

   !  call FILE_HISTORY_in( SFLX_rad_dn(:,:,I_R_direct ,I_R_IR ,:), 'SFLX_IR_dn_dir',  'SFC downward radiation flux (direct ,IR )', 'W/m2', fill_halo=.true. )
   !  call FILE_HISTORY_in( SFLX_rad_dn(:,:,I_R_diffuse,I_R_IR ,:), 'SFLX_IR_dn_dif',  'SFC downward radiation flux (diffuse,IR )', 'W/m2', fill_halo=.true. )
   !  call FILE_HISTORY_in( SFLX_rad_dn(:,:,I_R_direct ,I_R_NIR,:), 'SFLX_NIR_dn_dir', 'SFC downward radiation flux (direct ,NIR)', 'W/m2', fill_halo=.true. )
   !  call FILE_HISTORY_in( SFLX_rad_dn(:,:,I_R_diffuse,I_R_NIR,:), 'SFLX_NIR_dn_dif', 'SFC downward radiation flux (diffuse,NIR)', 'W/m2', fill_halo=.true. )
   !  call FILE_HISTORY_in( SFLX_rad_dn(:,:,I_R_direct ,I_R_VIS,:), 'SFLX_VIS_dn_dir', 'SFC downward radiation flux (direct ,VIS)', 'W/m2', fill_halo=.true. )
   !  call FILE_HISTORY_in( SFLX_rad_dn(:,:,I_R_diffuse,I_R_VIS,:), 'SFLX_VIS_dn_dif', 'SFC downward radiation flux (diffuse,VIS)', 'W/m2', fill_halo=.true. )

   !  call FILE_HISTORY_in( TEMP_t   (:,:,:,I_LW,:), 'TEMP_t_rd_LW', 'tendency of temp in rd(LW)',  'K/day',     fill_halo=.true. )
   !  call FILE_HISTORY_in( TEMP_t   (:,:,:,I_SW,:), 'TEMP_t_rd_SW', 'tendency of temp in rd(SW)',  'K/day',     fill_halo=.true. )
   !  call FILE_HISTORY_in( TEMP_t   (:,:,:,3   ,:), 'TEMP_t_rd',    'tendency of temp in rd',      'K/day',     fill_halo=.true. )
   !  call FILE_HISTORY_in( RHOH     (:,:,:,:),      'RHOH_RD',      'diabatic heating rate in rd', 'J/m3/s',    fill_halo=.true. )

   !  ! output of raw data, for offline output
   !  call FILE_HISTORY_in( flux_rad(:,:,:,I_LW,I_up,2,:), 'RFLX_LW_up', 'upward   longwave  radiation flux (cell face)', 'W/m2', fill_halo=.true. )
   !  call FILE_HISTORY_in( flux_rad(:,:,:,I_LW,I_dn,2,:), 'RFLX_LW_dn', 'downward longwave  radiation flux (cell face)', 'W/m2', fill_halo=.true. )
   !  call FILE_HISTORY_in( flux_rad(:,:,:,I_SW,I_up,2,:), 'RFLX_SW_up', 'upward   shortwave radiation flux (cell face)', 'W/m2', fill_halo=.true. )
   !  call FILE_HISTORY_in( flux_rad(:,:,:,I_SW,I_dn,2,:), 'RFLX_SW_dn', 'downward shortwave radiation flux (cell face)', 'W/m2', fill_halo=.true. )

   !  call FILE_HISTORY_in( dtau_s(:,:,:,:), 'dtau_s', '0.67 micron cloud optical depth', '1', fill_halo=.true. )
   !  call FILE_HISTORY_in( dem_s (:,:,:,:), 'dem_s',  '10.5 micron cloud emissivity',    '1', fill_halo=.true. )

    do l = 1, ADM_lall

       call history_in( 'SOLINS',         var_g(  solins(:,:,l)             ) )  ! 'solar insolation',        'W/m2'
       call history_in( 'COSZ'  ,         var_g(  cosSZA(:,:,l)             ) )  ! 'cos(solar zenith angle)', '1'

       call history_in( 'SFLX_LW_up_c',   var_g(  SFCFLX_LW_up_c(:,:,l)     ) )  ! 'SFC upward   longwave  radiation flux (clr)', 'W/m2'
       call history_in( 'SFLX_LW_dn_c',   var_g(  SFCFLX_LW_dn_c(:,:,l)     ) )  ! 'SFC downward longwave  radiation flux (clr)', 'W/m2'
       call history_in( 'SFLX_SW_up_c',   var_g(  SFCFLX_SW_up_c(:,:,l)     ) )  ! 'SFC upward   shortwave radiation flux (clr)', 'W/m2'
       call history_in( 'SFLX_SW_dn_c',   var_g(  SFCFLX_SW_dn_c(:,:,l)     ) )  ! 'SFC downward shortwave radiation flux (clr)', 'W/m2'

       call history_in( 'SFLX_LW_up',     var_g(  SFCFLX_LW_up  (:,:,l)     ) )  ! 'SFC upward   longwave  radiation flux',       'W/m2'
       call history_in( 'SFLX_LW_dn',     var_g(  SFCFLX_LW_dn  (:,:,l)     ) )  ! 'SFC downward longwave  radiation flux',       'W/m2'
       call history_in( 'SFLX_SW_up',     var_g(  SFCFLX_SW_up  (:,:,l)     ) )  ! 'SFC upward   shortwave radiation flux',       'W/m2'
       call history_in( 'SFLX_SW_dn',     var_g(  SFCFLX_SW_dn  (:,:,l)     ) )  ! 'SFC downward shortwave radiation flux',       'W/m2'

       call history_in( 'SFLX_LW_net',    var_g(  flux_net_sfc  (:,:,l,I_LW) ) ) ! 'SFC net      longwave  radiation flux',       'W/m2'
       call history_in( 'SFLX_SW_net',    var_g(  flux_net_sfc  (:,:,l,I_SW) ) ) ! 'SFC net      shortwave radiation flux',       'W/m2'

       call history_in( 'TOAFLX_LW_up_c', var_g(  TOAFLX_LW_up_c(:,:,l)      ) )  ! 'TOA upward   longwave  radiation flux (clr)', 'W/m2'
       call history_in( 'TOAFLX_LW_dn_c', var_g(  TOAFLX_LW_dn_c(:,:,l)      ) )  ! 'TOA downward longwave  radiation flux (clr)', 'W/m2'
       call history_in( 'TOAFLX_SW_up_c', var_g(  TOAFLX_SW_up_c(:,:,l)      ) )  ! 'TOA upward   shortwave radiation flux (clr)', 'W/m2'
       call history_in( 'TOAFLX_SW_dn_c', var_g(  TOAFLX_SW_dn_c(:,:,l)      ) )  ! 'TOA downward shortwave radiation flux (clr)', 'W/m2'

       call history_in( 'TOAFLX_LW_up',   var_g(  TOAFLX_LW_up  (:,:,l)      ) )  ! 'TOA upward   longwave  radiation flux',       'W/m2'
       call history_in( 'TOAFLX_LW_dn',   var_g(  TOAFLX_LW_dn  (:,:,l)      ) )  ! 'TOA downward longwave  radiation flux',       'W/m2'
       call history_in( 'TOAFLX_SW_up',   var_g(  TOAFLX_SW_up  (:,:,l)      ) )  ! 'TOA upward   shortwave radiation flux',       'W/m2'
       call history_in( 'TOAFLX_SW_dn',   var_g(  TOAFLX_SW_dn  (:,:,l)      ) )  ! 'TOA downward shortwave radiation flux',       'W/m2'

       call history_in( 'TOAFLX_LW_net',  var_g(  flux_net_toa  (:,:,l,I_LW) ) )  ! 'TOA net      longwave  radiation flux',       'W/m2'
       call history_in( 'TOAFLX_SW_net',  var_g(  flux_net_toa  (:,:,l,I_SW) ) )  ! 'TOA net      shortwave radiation flux',       'W/m2'

       call history_in( 'SLR',            var_g(  flux_net_sfc  (:,:,l,I_LW) ) )  ! 'SFC net longwave  radiation flux',            'W/m2'
       call history_in( 'SSR',            var_g(  flux_net_sfc  (:,:,l,I_SW) ) )  ! 'SFC net shortwave radiation flux',            'W/m2'
       call history_in( 'OLR',            var_g(  TOAFLX_LW_up  (:,:,l)      ) )  ! 'outgoing longwave  radiation flux',           'W/m2'
       call history_in( 'OSR',            var_g(  TOAFLX_SW_up  (:,:,l)      ) )  ! 'outgoing shortwave radiation flux',           'W/m2'

       call history_in( 'RADFLUX_LWUP',   var_gk( flux_up     (:,:,:,l,I_LW) ) )  ! 'upward   longwave  radiation flux',           'W/m2'
       call history_in( 'RADFLUX_LWDN',   var_gk( flux_dn     (:,:,:,l,I_LW) ) )  ! 'downward longwave  radiation flux',           'W/m2'
       call history_in( 'RADFLUX_LW',     var_gk( flux_net    (:,:,:,l,I_LW) ) )  ! 'net      longwave  radiation flux',           'W/m2'
       call history_in( 'RADFLUX_SWUP',   var_gk( flux_up     (:,:,:,l,I_SW) ) )  ! 'upward   shortwave radiation flux',           'W/m2'
       call history_in( 'RADFLUX_SWDN',   var_gk( flux_dn     (:,:,:,l,I_SW) ) )  ! 'downward shortwave radiation flux',           'W/m2'
       call history_in( 'RADFLUX_SW',     var_gk( flux_net    (:,:,:,l,I_SW) ) )  ! 'net      shortwave radiation flux',           'W/m2'

       call history_in( 'SFLX_IR_dn_dir', var_g(  SFLX_rad_dn(:,:,I_R_direct ,I_R_IR ,l) ) ) ! 'SFC downward radiation flux (direct ,IR )', 'W/m2'
       call history_in( 'SFLX_IR_dn_dif', var_g(  SFLX_rad_dn(:,:,I_R_diffuse,I_R_IR ,l) ) ) ! 'SFC downward radiation flux (diffuse,IR )', 'W/m2'
       call history_in( 'SFLX_NIR_dn_dir',var_g(  SFLX_rad_dn(:,:,I_R_direct ,I_R_NIR,l) ) ) ! 'SFC downward radiation flux (direct ,NIR)', 'W/m2'
       call history_in( 'SFLX_NIR_dn_dif',var_g(  SFLX_rad_dn(:,:,I_R_diffuse,I_R_NIR,l) ) ) ! 'SFC downward radiation flux (diffuse,NIR)', 'W/m2'
       call history_in( 'SFLX_VIS_dn_dir',var_g(  SFLX_rad_dn(:,:,I_R_direct ,I_R_VIS,l) ) ) ! 'SFC downward radiation flux (direct ,VIS)', 'W/m2'
       call history_in( 'SFLX_VIS_dn_dif',var_g(  SFLX_rad_dn(:,:,I_R_diffuse,I_R_VIS,l) ) ) ! 'SFC downward radiation flux (diffuse,VIS)', 'W/m2'

       call history_in( 'TEMP_t_rd_LW',   var_gk( TEMP_t   (:,:,:,I_LW,l)    ) ) !  'tendency of temp in rd(LW)',                    'K/day'
       call history_in( 'TEMP_t_rd_SW',   var_gk( TEMP_t   (:,:,:,I_SW,l)    ) ) !  'tendency of temp in rd(SW)',                    'K/day'
       call history_in( 'TEMP_t_rd',      var_gk( TEMP_t   (:,:,:,3   ,l)    ) ) !  'tendency of temp in rd',                        'K/day'
       call history_in( 'RHOH_RD',        var_gk( RHOH     (:,:,:,l)         ) ) !  'diabatic heating rate in rd',                   'J/m3/s'

       ! output of raw data, for offline output
       call history_in( 'RFLX_LW_up',     var_gk( flux_rad(:,:,:,I_LW,I_up,2,l) ) ) ! 'upward   longwave  radiation flux (cell face)', 'W/m2'
       call history_in( 'RFLX_LW_dn',     var_gk( flux_rad(:,:,:,I_LW,I_dn,2,l) ) ) ! 'downward longwave  radiation flux (cell face)', 'W/m2'
       call history_in( 'RFLX_SW_up',     var_gk( flux_rad(:,:,:,I_SW,I_up,2,l) ) ) ! 'upward   shortwave radiation flux (cell face)', 'W/m2'
       call history_in( 'RFLX_SW_dn',     var_gk( flux_rad(:,:,:,I_SW,I_dn,2,l) ) ) ! 'downward shortwave radiation flux (cell face)', 'W/m2'

       call history_in( 'dtau_s',         var_gk( dtau_s(:,:,:,l)            ) ) ! '0.67 micron cloud optical depth',                 '1'
       call history_in( 'dem_s',          var_gk( dem_s (:,:,:,l)            ) ) ! '10.5 micron cloud emissivity',                    '1'

    enddo


    call ATMOS_vars_calc_diagnostics

    return
  end subroutine ATMOS_PHY_RD_driver_step

  ! the following two functions are conversion to use history_in(), which will be replaced by FILE_HISTOTY_in in future.
  function var_g(var_ij)
     implicit none
     real(RP) :: var_ij (IA,JA)
     real(RP) :: var_g  (ADM_gall_in, ADM_KNONE)
     integer  :: i, j, g
     do j = 1, JA
     do i = 1, IA
        g = i + ( j - 1 ) * ADM_imax
        var_g(g,ADM_KNONE) = var_ij(i,j)
     enddo
     enddo
  end function var_g

  function var_gk(var_kij)
     implicit none
     real(RP) :: var_kij (KA, IA, JA)
     real(RP) :: var_gk  (ADM_gall_in, ADM_Kall)
     integer  :: i, j, k, g
     do j = 1, JA
     do i = 1, IA
        g = i + ( j - 1 ) * ADM_imax
        do k = 1, KA
           var_gk(g,k) = var_kij(k,i,j)
        enddo
     enddo
     enddo
  end function var_gk


end module mod_atmos_phy_rd_driver
