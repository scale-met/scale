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
  use scale_atmos_grid_cartesC_index
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
  public :: ATMOS_PHY_RD_driver_finalize
  public :: ATMOS_PHY_RD_driver_calc_tendency

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
  logical, private :: RD_use_PBL_cloud = .false.

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
    use scale_atmos_grid_cartesC, only: &
       CZ => ATMOS_GRID_CARTESC_CZ, &
       FZ => ATMOS_GRID_CARTESC_FZ
    use mod_atmos_admin, only: &
       ATMOS_PHY_RD_TYPE, &
       ATMOS_sw_phy_rd
    use mod_atmos_phy_rd_vars, only: &
       SFCFLX_LW_up => ATMOS_PHY_RD_SFLX_LW_up,   &
       SFCFLX_LW_dn => ATMOS_PHY_RD_SFLX_LW_dn,   &
       SFCFLX_SW_up => ATMOS_PHY_RD_SFLX_SW_up,   &
       SFCFLX_SW_dn => ATMOS_PHY_RD_SFLX_SW_dn,   &
       TOMFLX_LW_up => ATMOS_PHY_RD_TOMFLX_LW_up, &
       TOMFLX_LW_dn => ATMOS_PHY_RD_TOMFLX_LW_dn, &
       TOMFLX_SW_up => ATMOS_PHY_RD_TOMFLX_SW_up, &
       TOMFLX_SW_dn => ATMOS_PHY_RD_TOMFLX_SW_dn, &
       SFLX_rad_dn  => ATMOS_PHY_RD_SFLX_down,    &
       solins       => ATMOS_PHY_RD_solins,       &
       cosSZA       => ATMOS_PHY_RD_cosSZA
    implicit none
    
    namelist / PARAM_ATMOS_PHY_RD / &
         RD_use_PBL_cloud

    integer :: ierr
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("ATMOS_PHY_RD_driver_setup",*) 'Setup'

    if ( ATMOS_sw_phy_rd ) then

       !--- read namelist
       rewind(IO_FID_CONF)
       read(IO_FID_CONF,nml=PARAM_ATMOS_PHY_RD,iostat=ierr)
       if( ierr < 0 ) then !--- missing
          LOG_INFO("ATMOS_PHY_RD_driver_setup",*) 'Not found namelist. Default used.'
       elseif( ierr > 0 ) then !--- fatal error
          LOG_ERROR("ATMOS_PHY_RD_driver_setup",*) 'Not appropriate names in namelist PARAM_ATMOS_PHY_RD. Check!'
          call PRC_abort
       endif
       LOG_NML(PARAM_ATMOS_PHY_RD)

       select case ( ATMOS_PHY_RD_TYPE )
       case ( "MSTRNX" )
          call ATMOS_PHY_RD_MSTRNX_setup( KA, KS, KE, CZ(:), FZ(:) )
       case ( "OFFLINE" )
          call ATMOS_PHY_RD_offline_setup
       case default
          LOG_ERROR("ATMOS_PHY_RD_driver_setup",*) 'invalid Radiation type(', trim(ATMOS_PHY_RD_TYPE), '). CHECK!'
          call PRC_abort
       end select

    else

       LOG_INFO("ATMOS_PHY_RD_driver_setup",*) 'this component is never called.'
       LOG_INFO("ATMOS_PHY_RD_driver_setup",*) 'radiation fluxes are set to zero.'
       !$acc kernels
       SFCFLX_LW_up(:,:)     = 0.0_RP
       SFCFLX_LW_dn(:,:)     = 0.0_RP
       SFCFLX_SW_up(:,:)     = 0.0_RP
       SFCFLX_SW_dn(:,:)     = 0.0_RP
       TOMFLX_LW_up(:,:)     = 0.0_RP
       TOMFLX_LW_dn(:,:)     = 0.0_RP
       TOMFLX_SW_up(:,:)     = 0.0_RP
       TOMFLX_SW_dn(:,:)     = 0.0_RP
       SFLX_rad_dn (:,:,:,:) = 0.0_RP
       solins      (:,:)     = 0.0_RP
       cosSZA      (:,:)     = 0.0_RP
       !$acc end kernels

    endif

    return
  end subroutine ATMOS_PHY_RD_driver_setup

  !-----------------------------------------------------------------------------
  !> finalize
  subroutine ATMOS_PHY_RD_driver_finalize
    use scale_atmos_phy_rd_mstrnx, only: &
       ATMOS_PHY_RD_MSTRNX_finalize
    use scale_atmos_phy_rd_profile, only: &
       ATMOS_PHY_RD_PROFILE_finalize
    use mod_atmos_admin, only: &
       ATMOS_PHY_RD_TYPE, &
       ATMOS_sw_phy_rd
    implicit none
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("ATMOS_PHY_RD_driver_finalize",*) 'Finalize'

    if ( ATMOS_sw_phy_rd ) then
       select case ( ATMOS_PHY_RD_TYPE )
       case ( "MSTRNX" )
          call ATMOS_PHY_RD_MSTRNX_finalize
       case ( "OFFLINE" )
       end select
    end if

    call ATMOS_PHY_RD_PROFILE_finalize

    return
  end subroutine ATMOS_PHY_RD_driver_finalize

  !-----------------------------------------------------------------------------
  !> Driver
  subroutine ATMOS_PHY_RD_driver_calc_tendency( update_flag )
    use scale_const, only: &
       EPS => CONST_EPS
    use scale_atmos_grid_cartesC_real, only: &
       REAL_CZ  => ATMOS_GRID_CARTESC_REAL_CZ,  &
       REAL_FZ  => ATMOS_GRID_CARTESC_REAL_FZ,  &
       REAL_LON => ATMOS_GRID_CARTESC_REAL_LON, &
       REAL_LAT => ATMOS_GRID_CARTESC_REAL_LAT, &
       ATMOS_GRID_CARTESC_REAL_VOL, &
       ATMOS_GRID_CARTESC_REAL_TOTVOL
    use scale_landuse, only: &
       fact_ocean => LANDUSE_fact_ocean, &
       fact_land  => LANDUSE_fact_land,  &
       fact_urban => LANDUSE_fact_urban
    use scale_time, only: &
       TIME_NOWDATE,                     &
       TIME_NOWDAYSEC,                   &
       TIME_OFFSET_YEAR
    use scale_statistics, only: &
       STATISTICS_checktotal, &
       STATISTICS_total
    use scale_file_history, only: &
       FILE_HISTORY_in
    use scale_atmos_hydrometeor, only: &
       N_HYD, &
       I_HC, &
       I_HI
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
    use mod_atmos_vars, only: &
       TEMP,              &
       PRES,              &
       QV,                &
       CVtot,             &
       DENS   => DENS_av, &
       QTRC   => QTRC_av, &
       RHOH   => RHOH_p
    use mod_atmos_phy_sf_vars, only: &
       SFC_TEMP    => ATMOS_PHY_SF_SFC_TEMP,   &
       SFC_albedo  => ATMOS_PHY_SF_SFC_albedo
    use mod_atmos_phy_rd_vars, only: &
       RHOH_RD      => ATMOS_PHY_RD_RHOH,         &
       SFCFLX_LW_up => ATMOS_PHY_RD_SFLX_LW_up,   &
       SFCFLX_LW_dn => ATMOS_PHY_RD_SFLX_LW_dn,   &
       SFCFLX_SW_up => ATMOS_PHY_RD_SFLX_SW_up,   &
       SFCFLX_SW_dn => ATMOS_PHY_RD_SFLX_SW_dn,   &
       TOMFLX_LW_up => ATMOS_PHY_RD_TOMFLX_LW_up, &
       TOMFLX_LW_dn => ATMOS_PHY_RD_TOMFLX_LW_dn, &
       TOMFLX_SW_up => ATMOS_PHY_RD_TOMFLX_SW_up, &
       TOMFLX_SW_dn => ATMOS_PHY_RD_TOMFLX_SW_dn, &
       SFLX_rad_dn  => ATMOS_PHY_RD_SFLX_down,    &
       solins       => ATMOS_PHY_RD_solins,       &
       cosSZA       => ATMOS_PHY_RD_cosSZA
    use mod_atmos_phy_bl_vars, only: &
       ATMOS_PHY_BL_Zi,     &
       ATMOS_PHY_BL_QL,     &
       ATMOS_PHY_BL_cldfrac
    use mod_atmos_vars, only: &
       ATMOS_vars_get_diagnostic
    use mod_atmos_phy_mp_vars, only: &
       ATMOS_PHY_MP_vars_get_diagnostic
    use mod_atmos_phy_ae_vars, only: &
       ATMOS_PHY_AE_vars_get_diagnostic
    implicit none

    logical, intent(in) :: update_flag

    real(RP) :: TEMP_t      (KA,IA,JA,3)
    real(RP) :: flux_rad    (KA,IA,JA,2,2,2)
    real(RP) :: flux_rad_top(   IA,JA,2,2,2)

    real(RP) :: dtau_s      (KA,IA,JA) ! 0.67 micron cloud optical depth
    real(RP) :: dem_s       (KA,IA,JA) ! 10.5 micron cloud emissivity

    real(RP) :: flux_up     (KA,IA,JA,2)
    real(RP) :: flux_dn     (KA,IA,JA,2)
    real(RP) :: flux_net    (KA,IA,JA,2)
    real(RP) :: flux_net_sfc(   IA,JA,2)
    real(RP) :: flux_net_toa(   IA,JA,2)
    real(RP) :: flux_net_tom(   IA,JA,2)

    real(RP) :: TOAFLX_LW_up(IA,JA)
    real(RP) :: TOAFLX_LW_dn(IA,JA)
    real(RP) :: TOAFLX_SW_up(IA,JA)
    real(RP) :: TOAFLX_SW_dn(IA,JA)

    real(RP) :: SFCFLX_LW_up_c(IA,JA)
    real(RP) :: SFCFLX_LW_dn_c(IA,JA)
    real(RP) :: SFCFLX_SW_up_c(IA,JA)
    real(RP) :: SFCFLX_SW_dn_c(IA,JA)
    real(RP) :: TOAFLX_LW_up_c(IA,JA)
    real(RP) :: TOAFLX_LW_dn_c(IA,JA)
    real(RP) :: TOAFLX_SW_up_c(IA,JA)
    real(RP) :: TOAFLX_SW_dn_c(IA,JA)
    real(RP) :: TOMFLX_LW_up_c(IA,JA)
    real(RP) :: TOMFLX_LW_dn_c(IA,JA)
    real(RP) :: TOMFLX_SW_up_c(IA,JA)
    real(RP) :: TOMFLX_SW_dn_c(IA,JA)

    real(RP) :: CLDFRAC(KA,IA,JA)
    real(RP) :: MP_Re  (KA,IA,JA,N_HYD)
    real(RP) :: MP_Qe  (KA,IA,JA,N_HYD)
    real(RP) :: AE_Re  (KA,IA,JA,N_AE)
    real(RP) :: AE_Qe  (KA,IA,JA,N_AE)

    real(RP) :: RH(KA,IA,JA)

    integer  :: k, i, j
    !---------------------------------------------------------------------------

    if ( update_flag ) then

       !$acc data create(TEMP_t,flux_rad,flux_rad_top,dtau_s,dem_s,flux_up,flux_dn,flux_net,flux_net_sfc,flux_net_toa,flux_net_tom,TOAFLX_LW_up,TOAFLX_LW_dn,TOAFLX_SW_up,TOAFLX_SW_dn,SFCFLX_LW_up_c,SFCFLX_LW_dn_c,SFCFLX_SW_up_c,SFCFLX_SW_dn_c,TOAFLX_LW_up_c,TOAFLX_LW_dn_c,TOAFLX_SW_up_c,TOAFLX_SW_dn_c,TOMFLX_LW_up_c,TOMFLX_LW_dn_c,TOMFLX_SW_up_c,TOMFLX_SW_dn_c,CLDFRAC,MP_Re,MP_Qe,AE_Re,AE_Qe)

       call SOLARINS_insolation( IA, IS, IE, JA, JS, JE, &
                                 REAL_LON(:,:), REAL_LAT(:,:),      & ! [IN]
                                 TIME_NOWDATE(:), TIME_OFFSET_YEAR, & ! [IN]
                                 solins(:,:), cosSZA  (:,:)         ) ! [OUT]

       call ATMOS_PHY_MP_vars_get_diagnostic( &
            DENS(:,:,:), TEMP(:,:,:), QTRC(:,:,:,:), & ! [IN]
            CLDFRAC=CLDFRAC, Re=MP_Re, Qe=MP_Qe      ) ! [OUT]

       if ( RD_use_PBL_cloud ) then
          !$omp parallel do
          !$acc kernels
          do j = JS, JE
          do i = IS, IE
          do k = KS, KE
             if ( REAL_CZ(k,i,j) < ATMOS_PHY_BL_Zi(i,j) + REAL_FZ(KS-1,i,j) ) then
                CLDFRAC(k,i,j) = ATMOS_PHY_BL_cldfrac(k,i,j)
                if ( ATMOS_PHY_BL_cldfrac(k,i,j) > EPS ) then
                   MP_Qe(k,i,j,I_HC) = ATMOS_PHY_BL_QL(k,i,j) / ATMOS_PHY_BL_cldfrac(k,i,j)
                else
                   MP_Qe(k,i,j,I_HC) = 0.0_RP
                end if
                MP_Qe(k,i,j,I_HI) = 0.0_RP
                MP_Re(k,i,j,I_HC) = 8.0E-4
             end if
          end do
          end do
          end do
          !$acc end kernels
          
       end if

       call ATMOS_vars_get_diagnostic( "RH", RH )
       call ATMOS_PHY_AE_vars_get_diagnostic( &
            QTRC(:,:,:,:), RH(:,:,:), & ! [IN]
            Re=AE_Re, Qe=AE_Qe        ) ! [IN]


       select case ( ATMOS_PHY_RD_TYPE )
       case ( "MSTRNX" )

          !$acc update host(DENS,TEMP,PRES,QV,SFC_TEMP,SFC_albedo,solins,cosSZA,CLDFRAC,MP_Re,MP_Qe,AE_Re,AE_Qe)
          call ATMOS_PHY_RD_MSTRNX_flux( &
               KA, KS, KE, IA, IS, IE, JA, JS, JE, &
               DENS(:,:,:), TEMP(:,:,:), PRES(:,:,:), QV(:,:,:), & ! [IN]
               REAL_CZ(:,:,:), REAL_FZ(:,:,:),                   & ! [IN]
               fact_ocean(:,:), fact_land(:,:), fact_urban(:,:), & ! [IN]
               SFC_TEMP(:,:), SFC_albedo(:,:,:,:),               & ! [IN]
               solins(:,:), cosSZA(:,:),                         & ! [IN]
               CLDFRAC(:,:,:), MP_Re(:,:,:,:), MP_Qe(:,:,:,:),   & ! [IN]
               AE_Re(:,:,:,:), AE_Qe(:,:,:,:),                   & ! [IN]
               flux_rad(:,:,:,:,:,:),                            & ! [OUT]
               flux_rad_top(:,:,:,:,:), SFLX_rad_dn(:,:,:,:),    & ! [OUT]
               dtau_s = dtau_s(:,:,:), dem_s = dem_s(:,:,:)      ) ! [OUT]
          !$acc update device(flux_rad,flux_rad_top,SFLX_rad_dn,dtau_s,dem_s)

       case ( "OFFLINE" )

          call ATMOS_PHY_RD_offline_flux( &
               KA, KS, KE, IA, IS, IE, JA, JS, JE, &
               TIME_NOWDAYSEC,        & ! [IN]
               flux_rad(:,:,:,:,:,2), & ! [OUT]
               SFLX_rad_dn(:,:,:,:)   ) ! [OUT]
          !$acc update device(flux_rad,SFLX_rad_dn)
          !$acc kernels
          flux_rad(:,:,:,:,:,1)   = 0.0_RP ! clear sky
          flux_rad_top(:,:,:,:,:) = 0.0_RP
          dtau_s(:,:,:) = 0.0_RP
          dem_s(:,:,:) = 0.0_RP
          !$acc end kernels

       end select


       ! surface
!OCL XFILL
       !$omp parallel do
       !$acc kernels
       do j = JS, JE
       do i = IS, IE
          ! for clear-sky
          SFCFLX_LW_up_c(i,j)    = flux_rad(KS-1,i,j,I_LW,I_up,1)
          SFCFLX_LW_dn_c(i,j)    = flux_rad(KS-1,i,j,I_LW,I_dn,1)
          SFCFLX_SW_up_c(i,j)    = flux_rad(KS-1,i,j,I_SW,I_up,1)
          SFCFLX_SW_dn_c(i,j)    = flux_rad(KS-1,i,j,I_SW,I_dn,1)
          ! for all-sky
          SFCFLX_LW_up  (i,j)    = flux_rad(KS-1,i,j,I_LW,I_up,2)
          SFCFLX_LW_dn  (i,j)    = flux_rad(KS-1,i,j,I_LW,I_dn,2)
          SFCFLX_SW_up  (i,j)    = flux_rad(KS-1,i,j,I_SW,I_up,2)
          SFCFLX_SW_dn  (i,j)    = flux_rad(KS-1,i,j,I_SW,I_dn,2)

          flux_net_sfc(i,j,I_LW) = SFCFLX_LW_up(i,j) - SFCFLX_LW_dn(i,j)
          flux_net_sfc(i,j,I_SW) = SFCFLX_SW_up(i,j) - SFCFLX_SW_dn(i,j)
       enddo
       enddo
       !$acc end kernels

       ! top of the atmosphere
!OCL XFILL
       !$omp parallel do
       !$acc kernels
       do j = JS, JE
       do i = IS, IE
          ! for clear-sky
          TOAFLX_LW_up_c(i,j)    = flux_rad_top(i,j,I_LW,I_up,1)
          TOAFLX_LW_dn_c(i,j)    = flux_rad_top(i,j,I_LW,I_dn,1)
          TOAFLX_SW_up_c(i,j)    = flux_rad_top(i,j,I_SW,I_up,1)
          TOAFLX_SW_dn_c(i,j)    = flux_rad_top(i,j,I_SW,I_dn,1)
          ! for all-sky
          TOAFLX_LW_up  (i,j)    = flux_rad_top(i,j,I_LW,I_up,2)
          TOAFLX_LW_dn  (i,j)    = flux_rad_top(i,j,I_LW,I_dn,2)
          TOAFLX_SW_up  (i,j)    = flux_rad_top(i,j,I_SW,I_up,2)
          TOAFLX_SW_dn  (i,j)    = flux_rad_top(i,j,I_SW,I_dn,2)

          flux_net_toa(i,j,I_LW) = TOAFLX_LW_up(i,j) - TOAFLX_LW_dn(i,j)
          flux_net_toa(i,j,I_SW) = TOAFLX_SW_up(i,j) - TOAFLX_SW_dn(i,j)
       enddo
       enddo
       !$acc end kernels

       ! top of the model
!OCL XFILL
       !$omp parallel do
       !$acc kernels
       do j = JS, JE
       do i = IS, IE
          ! for clear-sky
          TOMFLX_LW_up_c(i,j)    = flux_rad(KE,i,j,I_LW,I_up,1)
          TOMFLX_LW_dn_c(i,j)    = flux_rad(KE,i,j,I_LW,I_dn,1)
          TOMFLX_SW_up_c(i,j)    = flux_rad(KE,i,j,I_SW,I_up,1)
          TOMFLX_SW_dn_c(i,j)    = flux_rad(KE,i,j,I_SW,I_dn,1)
          ! for all-sky
          TOMFLX_LW_up  (i,j)    = flux_rad(KE,i,j,I_LW,I_up,2)
          TOMFLX_LW_dn  (i,j)    = flux_rad(KE,i,j,I_LW,I_dn,2)
          TOMFLX_SW_up  (i,j)    = flux_rad(KE,i,j,I_SW,I_up,2)
          TOMFLX_SW_dn  (i,j)    = flux_rad(KE,i,j,I_SW,I_dn,2)

          flux_net_tom(i,j,I_LW) = TOMFLX_LW_up(i,j) - TOMFLX_LW_dn(i,j)
          flux_net_tom(i,j,I_SW) = TOMFLX_SW_up(i,j) - TOMFLX_SW_dn(i,j)
       enddo
       enddo
       !$acc end kernels

!OCL XFILL
       !$omp parallel do collapse(2)
       !$acc kernels
       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          flux_up (k,i,j,I_LW) = 0.5_RP * ( flux_rad(k-1,i,j,I_LW,I_up,2) + flux_rad(k,i,j,I_LW,I_up,2) )
          flux_dn (k,i,j,I_LW) = 0.5_RP * ( flux_rad(k-1,i,j,I_LW,I_dn,2) + flux_rad(k,i,j,I_LW,I_dn,2) )
          flux_up (k,i,j,I_SW) = 0.5_RP * ( flux_rad(k-1,i,j,I_SW,I_up,2) + flux_rad(k,i,j,I_SW,I_up,2) )
          flux_dn (k,i,j,I_SW) = 0.5_RP * ( flux_rad(k-1,i,j,I_SW,I_dn,2) + flux_rad(k,i,j,I_SW,I_dn,2) )

          flux_net(k,i,j,I_LW) = flux_up(k,i,j,I_LW) - flux_dn(k,i,j,I_LW)
          flux_net(k,i,j,I_SW) = flux_up(k,i,j,I_SW) - flux_dn(k,i,j,I_SW)
       enddo
       enddo
       enddo
       !$acc end kernels

       ! apply radiative flux convergence -> heating rate
       !$acc update host(CVtot)
       call ATMOS_PHY_RD_calc_heating( &
            KA, KS, KE, IA, IS, IE, JA, JS, JE, &
            flux_rad(:,:,:,:,:,2),     & ! [IN]
            DENS(:,:,:), TEMP(:,:,:),  & ! [IN]
            CVtot(:,:,:),              & ! [IN]
            REAL_FZ(:,:,:),            & ! [IN]
            RHOH_RD(:,:,:),            & ! [OUT]
            temp_t = TEMP_t(:,:,:,:)   ) ! [OUT]
       !$acc update device(RHOH_RD,TEMP_t)



       call FILE_HISTORY_in( solins(:,:), 'SOLINS', 'solar insolation',        'W/m2', fill_halo=.true. )
       call FILE_HISTORY_in( cosSZA(:,:), 'COSZ',   'cos(solar zenith angle)', '1',    fill_halo=.true. )

       call FILE_HISTORY_in( SFCFLX_LW_up_c(:,:),      'SFLX_LW_up_c',   'SFC upward   longwave  radiation flux (clr)', 'W/m2', fill_halo=.true. )
       call FILE_HISTORY_in( SFCFLX_LW_dn_c(:,:),      'SFLX_LW_dn_c',   'SFC downward longwave  radiation flux (clr)', 'W/m2', fill_halo=.true. )
       call FILE_HISTORY_in( SFCFLX_SW_up_c(:,:),      'SFLX_SW_up_c',   'SFC upward   shortwave radiation flux (clr)', 'W/m2', fill_halo=.true. )
       call FILE_HISTORY_in( SFCFLX_SW_dn_c(:,:),      'SFLX_SW_dn_c',   'SFC downward shortwave radiation flux (clr)', 'W/m2', fill_halo=.true. )

       call FILE_HISTORY_in( SFCFLX_LW_up  (:,:),      'SFLX_LW_up',     'SFC upward   longwave  radiation flux',       'W/m2', fill_halo=.true. )
       call FILE_HISTORY_in( SFCFLX_LW_dn  (:,:),      'SFLX_LW_dn',     'SFC downward longwave  radiation flux',       'W/m2', fill_halo=.true. )
       call FILE_HISTORY_in( SFCFLX_SW_up  (:,:),      'SFLX_SW_up',     'SFC upward   shortwave radiation flux',       'W/m2', fill_halo=.true. )
       call FILE_HISTORY_in( SFCFLX_SW_dn  (:,:),      'SFLX_SW_dn',     'SFC downward shortwave radiation flux',       'W/m2', fill_halo=.true. )

       call FILE_HISTORY_in( flux_net_sfc  (:,:,I_LW), 'SFLX_LW_net',    'SFC net      longwave  radiation flux',       'W/m2', fill_halo=.true. )
       call FILE_HISTORY_in( flux_net_sfc  (:,:,I_SW), 'SFLX_SW_net',    'SFC net      shortwave radiation flux',       'W/m2', fill_halo=.true. )

       call FILE_HISTORY_in( TOAFLX_LW_up_c(:,:),      'TOAFLX_LW_up_c', 'TOA upward   longwave  radiation flux (clr)', 'W/m2', fill_halo=.true. )
       call FILE_HISTORY_in( TOAFLX_LW_dn_c(:,:),      'TOAFLX_LW_dn_c', 'TOA downward longwave  radiation flux (clr)', 'W/m2', fill_halo=.true. )
       call FILE_HISTORY_in( TOAFLX_SW_up_c(:,:),      'TOAFLX_SW_up_c', 'TOA upward   shortwave radiation flux (clr)', 'W/m2', fill_halo=.true. )
       call FILE_HISTORY_in( TOAFLX_SW_dn_c(:,:),      'TOAFLX_SW_dn_c', 'TOA downward shortwave radiation flux (clr)', 'W/m2', fill_halo=.true. )

       call FILE_HISTORY_in( TOAFLX_LW_up  (:,:),      'TOAFLX_LW_up',   'TOA upward   longwave  radiation flux',       'W/m2', fill_halo=.true. )
       call FILE_HISTORY_in( TOAFLX_LW_dn  (:,:),      'TOAFLX_LW_dn',   'TOA downward longwave  radiation flux',       'W/m2', fill_halo=.true. )
       call FILE_HISTORY_in( TOAFLX_SW_up  (:,:),      'TOAFLX_SW_up',   'TOA upward   shortwave radiation flux',       'W/m2', fill_halo=.true. )
       call FILE_HISTORY_in( TOAFLX_SW_dn  (:,:),      'TOAFLX_SW_dn',   'TOA downward shortwave radiation flux',       'W/m2', fill_halo=.true. )

       call FILE_HISTORY_in( flux_net_toa  (:,:,I_LW), 'TOAFLX_LW_net',  'TOA net      longwave  radiation flux',       'W/m2', fill_halo=.true. )
       call FILE_HISTORY_in( flux_net_toa  (:,:,I_SW), 'TOAFLX_SW_net',  'TOA net      shortwave radiation flux',       'W/m2', fill_halo=.true. )

       call FILE_HISTORY_in( TOMFLX_LW_up_c(:,:),      'TOMFLX_LW_up_c', 'TOM upward   longwave  radiation flux (clr)', 'W/m2', fill_halo=.true. )
       call FILE_HISTORY_in( TOMFLX_LW_dn_c(:,:),      'TOMFLX_LW_dn_c', 'TOM downward longwave  radiation flux (clr)', 'W/m2', fill_halo=.true. )
       call FILE_HISTORY_in( TOMFLX_SW_up_c(:,:),      'TOMFLX_SW_up_c', 'TOM upward   shortwave radiation flux (clr)', 'W/m2', fill_halo=.true. )
       call FILE_HISTORY_in( TOMFLX_SW_dn_c(:,:),      'TOMFLX_SW_dn_c', 'TOM downward shortwave radiation flux (clr)', 'W/m2', fill_halo=.true. )

       call FILE_HISTORY_in( TOMFLX_LW_up  (:,:),      'TOMFLX_LW_up',   'TOM upward   longwave  radiation flux',       'W/m2', fill_halo=.true. )
       call FILE_HISTORY_in( TOMFLX_LW_dn  (:,:),      'TOMFLX_LW_dn',   'TOM downward longwave  radiation flux',       'W/m2', fill_halo=.true. )
       call FILE_HISTORY_in( TOMFLX_SW_up  (:,:),      'TOMFLX_SW_up',   'TOM upward   shortwave radiation flux',       'W/m2', fill_halo=.true. )
       call FILE_HISTORY_in( TOMFLX_SW_dn  (:,:),      'TOMFLX_SW_dn',   'TOM downward shortwave radiation flux',       'W/m2', fill_halo=.true. )

       call FILE_HISTORY_in( flux_net_tom  (:,:,I_LW), 'TOMFLX_LW_net',  'TOM net      longwave  radiation flux',       'W/m2', fill_halo=.true. )
       call FILE_HISTORY_in( flux_net_tom  (:,:,I_SW), 'TOMFLX_SW_net',  'TOM net      shortwave radiation flux',       'W/m2', fill_halo=.true. )

       call FILE_HISTORY_in( flux_net_sfc(:,:,I_LW),   'SLR',          'SFC net longwave  radiation flux',  'W/m2', fill_halo=.true. )
       call FILE_HISTORY_in( flux_net_sfc(:,:,I_SW),   'SSR',          'SFC net shortwave radiation flux',  'W/m2', fill_halo=.true. )
       call FILE_HISTORY_in( TOAFLX_LW_up(:,:),        'OLR',          'outgoing longwave  radiation flux', 'W/m2', fill_halo=.true. )
       call FILE_HISTORY_in( TOAFLX_SW_up(:,:),        'OSR',          'outgoing shortwave radiation flux', 'W/m2', fill_halo=.true. )

       call FILE_HISTORY_in( flux_up     (:,:,:,I_LW), 'RADFLUX_LWUP', 'upward   longwave  radiation flux', 'W/m2', fill_halo=.true. )
       call FILE_HISTORY_in( flux_dn     (:,:,:,I_LW), 'RADFLUX_LWDN', 'downward longwave  radiation flux', 'W/m2', fill_halo=.true. )
       call FILE_HISTORY_in( flux_net    (:,:,:,I_LW), 'RADFLUX_LW',   'net      longwave  radiation flux', 'W/m2', fill_halo=.true. )
       call FILE_HISTORY_in( flux_up     (:,:,:,I_SW), 'RADFLUX_SWUP', 'upward   shortwave radiation flux', 'W/m2', fill_halo=.true. )
       call FILE_HISTORY_in( flux_dn     (:,:,:,I_SW), 'RADFLUX_SWDN', 'downward shortwave radiation flux', 'W/m2', fill_halo=.true. )
       call FILE_HISTORY_in( flux_net    (:,:,:,I_SW), 'RADFLUX_SW',   'net      shortwave radiation flux', 'W/m2', fill_halo=.true. )

       call FILE_HISTORY_in( SFLX_rad_dn(:,:,I_R_direct ,I_R_IR ), 'SFLX_IR_dn_dir',  'SFC downward radiation flux (direct; IR)',   'W/m2', fill_halo=.true. )
       call FILE_HISTORY_in( SFLX_rad_dn(:,:,I_R_diffuse,I_R_IR ), 'SFLX_IR_dn_dif',  'SFC downward radiation flux (diffuse; IR)',  'W/m2', fill_halo=.true. )
       call FILE_HISTORY_in( SFLX_rad_dn(:,:,I_R_direct ,I_R_NIR), 'SFLX_NIR_dn_dir', 'SFC downward radiation flux (direct; NIR)',  'W/m2', fill_halo=.true. )
       call FILE_HISTORY_in( SFLX_rad_dn(:,:,I_R_diffuse,I_R_NIR), 'SFLX_NIR_dn_dif', 'SFC downward radiation flux (diffuse; NIR)', 'W/m2', fill_halo=.true. )
       call FILE_HISTORY_in( SFLX_rad_dn(:,:,I_R_direct ,I_R_VIS), 'SFLX_VIS_dn_dir', 'SFC downward radiation flux (direct; VIS)',  'W/m2', fill_halo=.true. )
       call FILE_HISTORY_in( SFLX_rad_dn(:,:,I_R_diffuse,I_R_VIS), 'SFLX_VIS_dn_dif', 'SFC downward radiation flux (diffuse; VIS)', 'W/m2', fill_halo=.true. )

       call FILE_HISTORY_in( TEMP_t   (:,:,:,I_LW), 'TEMP_t_rd_LW', 'tendency of temp in rd(LW)',  'K/day',     fill_halo=.true. )
       call FILE_HISTORY_in( TEMP_t   (:,:,:,I_SW), 'TEMP_t_rd_SW', 'tendency of temp in rd(SW)',  'K/day',     fill_halo=.true. )
       call FILE_HISTORY_in( TEMP_t   (:,:,:,3   ), 'TEMP_t_rd',    'tendency of temp in rd',      'K/day',     fill_halo=.true. )
       call FILE_HISTORY_in( RHOH_RD  (:,:,:),      'RHOH_RD',      'diabatic heating rate in rd', 'J/m3/s',    fill_halo=.true. )

       ! output of raw data, for offline output
       call FILE_HISTORY_in( flux_rad(:,:,:,I_LW,I_up,2), 'RFLX_LW_up', 'upward   longwave  radiation flux (cell face)', 'W/m2', fill_halo=.true. )
       call FILE_HISTORY_in( flux_rad(:,:,:,I_LW,I_dn,2), 'RFLX_LW_dn', 'downward longwave  radiation flux (cell face)', 'W/m2', fill_halo=.true. )
       call FILE_HISTORY_in( flux_rad(:,:,:,I_SW,I_up,2), 'RFLX_SW_up', 'upward   shortwave radiation flux (cell face)', 'W/m2', fill_halo=.true. )
       call FILE_HISTORY_in( flux_rad(:,:,:,I_SW,I_dn,2), 'RFLX_SW_dn', 'downward shortwave radiation flux (cell face)', 'W/m2', fill_halo=.true. )

       call FILE_HISTORY_in( dtau_s(:,:,:), 'dtau_s', '0.67 micron cloud optical depth', '1', fill_halo=.true. )
       call FILE_HISTORY_in( dem_s (:,:,:), 'dem_s',  '10.5 micron cloud emissivity',    '1', fill_halo=.true. )

       !$acc end data

    endif

    !$acc kernels
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       RHOH(k,i,j) = RHOH(k,i,j) + RHOH_RD(k,i,j)
    enddo
    enddo
    enddo
    !$acc end kernels

    if ( STATISTICS_checktotal ) then
       call STATISTICS_total( KA, KS, KE, IA, IS, IE, JA, JS, JE, &
                              RHOH_RD(:,:,:), 'RHOH_RD',          &
                              ATMOS_GRID_CARTESC_REAL_VOL(:,:,:), &
                              ATMOS_GRID_CARTESC_REAL_TOTVOL      )
    endif

    return
  end subroutine ATMOS_PHY_RD_driver_calc_tendency

end module mod_atmos_phy_rd_driver
