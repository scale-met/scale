!-------------------------------------------------------------------------------
!> module LAND driver
!!
!! @par Description
!!          Land model driver
!!
!! @author Team SCALE
!<
!-------------------------------------------------------------------------------
#include "scalelib.h"
module mod_land_driver
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_io
  use scale_prof
  use scale_land_grid_cartesC_index
  use scale_tracer
  use scale_cpl_sfc_index

  use scale_const, only: &
     I_LW => CONST_I_LW, &
     I_SW => CONST_I_SW
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: LAND_driver_setup
  public :: LAND_driver_calc_tendency
  public :: LAND_driver_update
  public :: LAND_SURFACE_GET
  public :: LAND_SURFACE_SET

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
  subroutine LAND_driver_setup
    use scale_prc, only: &
       PRC_abort
    use mod_land_admin, only: &
       LAND_do, &
       LAND_DYN_TYPE, &
       LAND_SFC_TYPE, &
       SNOW_TYPE
    use scale_land_dyn_bucket, only: &
       LAND_DYN_BUCKET_setup
    use scale_land_phy_snow_ky90, only: &
       LAND_PHY_SNOW_KY90_setup
    use scale_cpl_phy_sfc_skin, only: &
       CPL_PHY_SFC_skin_setup
    use scale_cpl_phy_sfc_fixed_temp, only: &
       CPL_PHY_SFC_fixed_temp_setup
    implicit none
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("LAND_driver_setup",*) 'Setup'

    if ( LAND_do ) then

       select case ( LAND_DYN_TYPE )
       case ( 'BUCKET' )
          call LAND_DYN_BUCKET_setup
       case ( 'INIT' )
          ! do nothing
       case default
          LOG_ERROR("LAND_driver_setup",*) 'LAND_DYN_TYPE is invalid: ', trim(LAND_DYN_TYPE)
          call PRC_abort
       end select

       select case ( LAND_SFC_TYPE )
       case ( 'SKIN' )
          call CPL_PHY_SFC_skin_setup
       case ( 'FIXED-TEMP' )
          call CPL_PHY_SFC_fixed_temp_setup
       case default
          LOG_ERROR("LAND_driver_setup",*) 'LAND_SFC_TYPE is invalid: ', trim(LAND_SFC_TYPE)
          call PRC_abort
       end select

       select case ( SNOW_TYPE )
       case ( 'NONE', 'OFF' )
       case ( 'KY90' )
          LOG_WARN("LAND_driver_setup",*) 'SNOW model is enabled'
          LOG_WARN("LAND_driver_setup",*) 'SNOW model is on experimental stage.'
          LOG_WARN("LAND_driver_setup",*) 'Use this with your own risk.'
          call LAND_PHY_SNOW_KY90_setup
       case default
          LOG_ERROR("LAND_driver_setup",*) 'SNOW_TYPE is invalid: ', trim(SNOW_TYPE)
          call PRC_abort
       end select

    end if

    return
  end subroutine LAND_driver_setup

  !-----------------------------------------------------------------------------
  !> Calculate tendency
  subroutine LAND_driver_calc_tendency( force )
    use scale_const, only: &
       TEM00 => CONST_TEM00
    use scale_time, only: &
       dt => TIME_DTSEC_LAND
    use scale_file_history, only: &
       FILE_HISTORY_in
    use scale_atmos_grid_cartesC_real, only: &
       REAL_Z1 => ATMOS_GRID_CARTESC_REAL_Z1
    use scale_topography, only: &
       TanSL_X => TOPOGRAPHY_TanSL_X, &
       TanSL_Y => TOPOGRAPHY_TanSL_Y
    use scale_atmos_hydrometeor, only: &
       HYDROMETEOR_LHV => ATMOS_HYDROMETEOR_LHV, &
       HYDROMETEOR_LHS => ATMOS_HYDROMETEOR_LHS, &
       ATMOS_HYDROMETEOR_dry,                    &
       CV_WATER,                                 &
       CV_ICE,                                   &
       LHF,                                      &
       I_QV
    use scale_land_grid_cartesC, only: &
       LCZ => LAND_GRID_CARTESC_CZ, &
       CDZ => LAND_GRID_CARTESC_CDZ
    use scale_land_phy_snow_ky90, only: &
       LAND_PHY_SNOW_KY90
    use scale_land_phy_snow_diagnos, only: &
       LAND_PHY_SNOW_DIAGS
    use scale_cpl_phy_sfc_skin, only: &
       CPL_PHY_SFC_skin
    use scale_cpl_phy_sfc_fixed_temp, only: &
       CPL_PHY_SFC_fixed_temp
    use scale_bulkflux, only: &
       BULKFLUX_diagnose
    use mod_atmos_admin, only: &
       ATMOS_sw_phy_ch
    use mod_atmos_phy_ch_driver, only: &
       ATMOS_PHY_CH_driver_LAND_flux
    use mod_land_admin, only: &
       LAND_SFC_TYPE, &
       SNOW_TYPE
    use mod_land_vars, only: &
       I_WaterLimit,      &
       I_WaterCritical,   &
       I_StomataResist,   &
       I_ThermalCond,     &
       I_HeatCapacity,    &
       I_WaterDiff,       &
       I_ALBLW,           &
       I_ALBSW,           &
       I_Z0M,             &
       I_Z0H,             &
       I_Z0E,             &
       SNOW_flag,         &
       LAND_PROPERTY,     &
       LAND_TEMP,         &
       LAND_WATER,        &
       LAND_ICE,          &
       LAND_SFC_TEMP,     &
       LAND_SFC_albedo,   &
       SNOW_SFC_TEMP,     &
       SNOW_SWE,          &
       SNOW_Depth,        &
       SNOW_Dzero,        &
       SNOW_nosnowsec,    &
       LAND_TEMP_t,       &
       LAND_WATER_t,      &
       LAND_ICE_t,        &
       LAND_SFLX_GH,      &
       LAND_SFLX_water,   &
       LAND_SFLX_ENGI,    &
       LAND_SFLX_MW,      &
       LAND_SFLX_MU,      &
       LAND_SFLX_MV,      &
       LAND_SFLX_SH,      &
       LAND_SFLX_LH,      &
       LAND_SFLX_QTRC,    &
       LAND_U10,          &
       LAND_V10,          &
       LAND_T2,           &
       LAND_Q2,           &
       LAND_Ustar,        &
       LAND_Tstar,        &
       LAND_Qstar,        &
       LAND_Wstar,        &
       LAND_RLmo,         &
       SOIL_Ustar,        &
       SOIL_Tstar,        &
       SOIL_Qstar,        &
       SOIL_Wstar,        &
       SOIL_RLmo,         &
       SNOW_Ustar,        &
       SNOW_Tstar,        &
       SNOW_Qstar,        &
       SNOW_Wstar,        &
       SNOW_RLmo,         &
       ATMOS_TEMP,        &
       ATMOS_PRES,        &
       ATMOS_U,           &
       ATMOS_V,           &
       ATMOS_DENS,        &
       ATMOS_QV,          &
       ATMOS_PBL,         &
       ATMOS_SFC_DENS,    &
       ATMOS_SFC_PRES,    &
       ATMOS_SFLX_rad_dn, &
       ATMOS_SFLX_water,  &
       ATMOS_SFLX_ENGI
    use scale_landuse, only: &
       LANDUSE_fact_land, &
       LANDUSE_exists_land
    implicit none

    logical, intent(in) :: force

    ! parameters
    real(RP), parameter :: BETA_MAX = 1.0_RP

    ! works
    real(RP) :: SNOW_QVEF (LIA,LJA)
    real(RP) :: LAND_WSTR (LIA,LJA)
    real(RP) :: LAND_QVEF (LIA,LJA)
    real(RP) :: LAND_TC_dZ(LIA,LJA)
    real(RP) :: SFLX_QV   (LIA,LJA)
    real(RP) :: SFLX_ENGI (LIA,LJA)
    real(RP) :: LH        (LIA,LJA) ! latent heat of vaporization [J/kg]
    real(RP) :: ATMOS_W   (LIA,LJA)
    real(RP) :: total

    ! for snow
    real(RP) :: SNOW_albedo         (LIA,LJA,2)
    real(RP) :: SNOW_ATMOS_SFLX_SH  (LIA,LJA)
    real(RP) :: SNOW_ATMOS_SFLX_LH  (LIA,LJA)
    real(RP) :: SNOW_ATMOS_SFLX_GH  (LIA,LJA)
    real(RP) :: SNOW_ATMOS_SFLX_QV  (LIA,LJA)
    real(RP) :: SNOW_LAND_SFLX_GH   (LIA,LJA)
    real(RP) :: SNOW_LAND_SFLX_water(LIA,LJA)
    real(RP) :: SNOW_LAND_SFLX_ENGI (LIA,LJA)
    real(RP) :: SNOW_frac           (LIA,LJA)

    real(RP) :: SNOW_ATMOS_SFLX_MW  (LIA,LJA)
    real(RP) :: SNOW_ATMOS_SFLX_MU  (LIA,LJA)
    real(RP) :: SNOW_ATMOS_SFLX_MV  (LIA,LJA)
    real(RP) :: SNOW_U10            (LIA,LJA)
    real(RP) :: SNOW_V10            (LIA,LJA)
    real(RP) :: SNOW_T2             (LIA,LJA)
    real(RP) :: SNOW_Q2             (LIA,LJA)

    ! monitor
    !real(RP) :: MONIT_WCONT0        (LIA,LJA)
    !real(RP) :: MONIT_WCONT1        (LIA,LJA)
    !real(RP) :: MONIT_ENG0          (LIA,LJA)
    !real(RP) :: MONIT_ENG1          (LIA,LJA)
    !
    !real(RP) :: MONIT_SNOW_heat     (LIA,LJA)
    !real(RP) :: MONIT_SNOW_water    (LIA,LJA)
    !real(RP) :: MONIT_LAND_heat     (LIA,LJA)
    !real(RP) :: MONIT_LAND_water    (LIA,LJA)

    integer :: k, i, j, iq, idir
    !---------------------------------------------------------------------------

    call PROF_rapstart('LND_CalcTend', 1)

    !########## Get Surface Boundary from coupler ##########
    call LAND_SURFACE_GET

    !########## reset tendencies ##########
!OCL XFILL
    !$omp parallel do
    do j = LJS, LJE
    do i = LIS, LIE
    do k = LKS, LKE
       LAND_TEMP_t (k,i,j) = 0.0_RP
       LAND_WATER_t(k,i,j) = 0.0_RP
       LAND_ICE_t  (k,i,j) = 0.0_RP
    enddo
    enddo
    enddo
!OCL XFILL
    do iq = 1, QA
    !$omp parallel do
    do j  = LJS, LJE
    do i  = LIS, LIE
       LAND_SFLX_QTRC(i,j,iq) = 0.0_RP
    enddo
    enddo
    enddo

    !$omp parallel do
    do j = LJS, LJE
    do i = LIS, LIE
       ATMOS_W(i,j) = ATMOS_U(i,j) * TanSL_X(i,j) + ATMOS_V(i,j) * TanSL_Y(i,j)
    end do
    end do

    if ( SNOW_flag ) then
       !------------------------------------------------------------------------
       !> snow area

!OCL XFILL
       !$omp parallel do
       do j = LJS, LJE
       do i = LIS, LIE
          ! This is for debug---adachi start
          !if(( int(SNOW_frac(i,j)) == 1 ).and.( abs(SNOW_SFC_TEMP(i,j)-LAND_SFC_TEMP(i,j))/=0 ))then
          !   LOG_ERROR("LAND_driver_calc_tendency",*) "Error please check SNOW_SFC_TEMP routine"
          !   call PRC_abort
          !endif
          ! This is for debug---adachi end
          SNOW_SFC_TEMP(i,j) = LAND_SFC_TEMP(i,j)
       end do
       end do

       select case ( SNOW_TYPE )
       case ( 'KY90' )
          ! accumulation and melt of snow if there is snow

          !MONIT_WCONT0 = 0.0_RP
          !call monitor_snow_water(  SNOW_Depth          (:,:),   & ! [IN]
          !                          SNOW_Dzero          (:,:),   & ! [IN]
          !                          MONIT_WCONT0        (:,:)    ) ! [OUT]

          call LAND_PHY_SNOW_KY90( LIA, LIS, LIE, LJA, LJS, LJE, &
                                   ATMOS_SFLX_water(:,:), ATMOS_SFLX_ENGI(:,:),      & ! [IN]
                                   ATMOS_PRES(:,:), ATMOS_TEMP(:,:), ATMOS_QV(:,:),  & ! [IN]
                                   ATMOS_W(:,:), ATMOS_U(:,:), ATMOS_V(:,:),         & ! [IN]
                                   ATMOS_SFC_DENS(:,:),                              & ! [IN]
                                   ATMOS_SFLX_rad_dn(:,:,:,:),                       & ! [IN]
                                   LANDUSE_fact_land(:,:), dt,                       & ! [IN]
                                   SNOW_SFC_TEMP(:,:), SNOW_SWE(:,:),                & ! [INOUT]
                                   SNOW_Depth(:,:), SNOW_Dzero(:,:),                 & ! [INOUT]
                                   SNOW_nosnowsec(:,:),                              & ! [INOUT]
                                   SNOW_albedo(:,:,:),                               & ! [OUT]
                                   SNOW_ATMOS_SFLX_SH(:,:),                          & ! [OUT]
                                   SNOW_ATMOS_SFLX_LH(:,:), SNOW_ATMOS_SFLX_QV(:,:), & ! [OUT]
                                   SFLX_ENGI(:,:),                                   & ! [OUT]
                                   SNOW_ATMOS_SFLX_GH(:,:), SNOW_LAND_SFLX_GH(:,:),  & ! [OUT]
                                   SNOW_LAND_SFLX_water(:,:),                        & ! [OUT]
                                   SNOW_frac           (:,:)                         ) ! [OUT]

!OCL XFILL
          !$omp parallel do
          do j = LJS, LJE
          do i = LIS, LIE
             SNOW_LAND_SFLX_ENGI(i,j) = ATMOS_SFLX_ENGI(i,j) & ! internal energy of precipitation
                                      - SFLX_ENGI(i,j)         ! internal energy of evapolation
          enddo
          enddo
       end select

!OCL XFILL
       !call monitor_snow_water(  SNOW_Depth          (:,:),   & ! [IN]
       !                          SNOW_Dzero          (:,:),   & ! [IN]
       !                          MONIT_WCONT1        (:,:)    ) ! [OUT]

       !call monitor_land_regidual( ATMOS_SFLX_water    (:,:),   & ! [IN] ! downward at surface
       !                            ATMOS_SFLX_ENGI     (:,:),   & ! [IN] ! downward at surface
       !                            SNOW_ATMOS_SFLX_evap(:,:),   & ! [IN] ! upward   at surface
       !                            SNOW_LAND_SFLX_water(:,:),   & ! [IN] ! downward at bottom
       !                            MONIT_WCONT0        (:,:),   & ! [IN]
       !                            MONIT_WCONT1        (:,:),   & ! [IN]
       !                            MONIT_SNOW_water    (:,:)    ) ! [OUT]

!OCL XFILL
       !$omp parallel do
       do j = LJS, LJE
       do i = LIS, LIE
          SNOW_QVEF(i,j) = 1.0_RP ! tentative
       end do
       end do

       ! momentum fluxes and diagnostic variables above snowpack
       call LAND_PHY_SNOW_DIAGS( LIA, LIS, LIE, LJA, LJS, LJE, &
                                 SNOW_frac(:,:),                                                & ! [IN]
                                 ATMOS_TEMP(:,:), ATMOS_PRES(:,:),                              & ! [IN]
                                 ATMOS_W(:,:), ATMOS_U(:,:), ATMOS_V(:,:),                      & ! [IN]
                                 ATMOS_DENS(:,:), ATMOS_QV(:,:),                                & ! [IN]
                                 REAL_Z1(:,:), ATMOS_PBL(:,:),                                  & ! [IN]
                                 ATMOS_SFC_DENS (:,:), ATMOS_SFC_PRES(:,:), SNOW_SFC_TEMP(:,:), & ! [IN]
                                 SNOW_QVEF(:,:),                                                & ! [IN]
                                 LAND_PROPERTY(:,:,I_Z0M),                                      & ! [IN]
                                 LAND_PROPERTY(:,:,I_Z0H),                                      & ! [IN]
                                 LAND_PROPERTY(:,:,I_Z0E),                                      & ! [IN]
                                 SNOW_ATMOS_SFLX_MW(:,:),                                       & ! [OUT]
                                 SNOW_ATMOS_SFLX_MU(:,:),                                       & ! [OUT]
                                 SNOW_ATMOS_SFLX_MV(:,:),                                       & ! [OUT]
                                 SNOW_Ustar(:,:), SNOW_Tstar(:,:), SNOW_Qstar(:,:),             & ! [OUT]
                                 SNOW_Wstar(:,:),                                               & ! [OUT]
                                 SNOW_RLmo(:,:),                                                & ! [OUT]
                                 SNOW_U10(:,:), SNOW_V10(:,:),                                  & ! [OUT]
                                 SNOW_T2(:,:), SNOW_Q2(:,:)                                     ) ! [OUT]

       call FILE_HISTORY_in( SNOW_frac         (:,:),      'LAND_SNOW_frac',    'Snow fraction on land subgrid',           '1'      )
       call FILE_HISTORY_in( SNOW_albedo       (:,:,I_SW), 'LAND_SNOW_ALB_SW',  'Snow surface albedo (short wave)',        '1'      )
       call FILE_HISTORY_in( SNOW_albedo       (:,:,I_LW), 'LAND_SNOW_ALB_LW',  'Snow surface albedo (long wave)',         '1'      )
       call FILE_HISTORY_in( SNOW_ATMOS_SFLX_SH(:,:),      'LAND_SNOW_SFLX_SH', 'Snow surface sensible heat flux',         'J/m2/s' )
       call FILE_HISTORY_in( SNOW_ATMOS_SFLX_LH(:,:),      'LAND_SNOW_SFLX_LH', 'Snow surface latent heat flux',           'J/m2/s' )
       call FILE_HISTORY_in( SNOW_ATMOS_SFLX_GH(:,:),      'LAND_SNOW_SFLX_GH', 'Snowpack received heat flux',             'J/m2/s' )
       call FILE_HISTORY_in( SNOW_ATMOS_SFLX_MW(:,:),      'LAND_SNOW_SFLX_MW', 'Snow surface w-momentum flux',            'J/m2/s' )
       call FILE_HISTORY_in( SNOW_ATMOS_SFLX_MU(:,:),      'LAND_SNOW_SFLX_MU', 'Snow surface u-momentum flux',            'J/m2/s' )
       call FILE_HISTORY_in( SNOW_ATMOS_SFLX_MV(:,:),      'LAND_SNOW_SFLX_MV', 'Snow surface v-momentum flux',            'J/m2/s' )
       call FILE_HISTORY_in( SNOW_U10          (:,:),      'LAND_SNOW_U10',     'Wind velocity u at 10 m on snow surface', 'm/s'    )
       call FILE_HISTORY_in( SNOW_V10          (:,:),      'LAND_SNOW_V10',     'Wind velocity v at 10 m on snow surface', 'm/s'    )
       call FILE_HISTORY_in( SNOW_T2           (:,:),      'LAND_SNOW_T2',      'Air temperature at 2m on snow surface',   'K'      )
       call FILE_HISTORY_in( SNOW_Q2           (:,:),      'LAND_SNOW_Q2',      'Specific humidity at 2m on snow surface', 'kg/kg'  )

       call FILE_HISTORY_in( SNOW_LAND_SFLX_GH   (:,:), 'LAND_SNOW_LAND_SFLX_GH',    'land surface ground heat flux under snow',     'J/m2/s'  )
       call FILE_HISTORY_in( SNOW_LAND_SFLX_water(:,:), 'LAND_SNOW_LAND_SFLX_water', 'land surface water mass flux under snow',      'kg/m2/s' )
       call FILE_HISTORY_in( SNOW_LAND_SFLX_ENGI (:,:), 'LAND_SNOW_LAND_SFLX_ENGI',  'land surface internal energy flux under snow', 'kg/m2/s' )
    endif


!OCL XFILL
    !$omp parallel do &
    !$omp private(total)
    do j = LJS, LJE
    do i = LIS, LIE
       total = LAND_WATER(LKS,i,j) + LAND_ICE(LKS,i,j)
       LAND_WSTR(i,j) = total * CDZ(LKS) &
                      + dt * ( ATMOS_SFLX_water(i,j) &
                             + max( 0.0_RP, 2.0_RP * LAND_PROPERTY(i,j,I_WaterDiff) &
                                    * ( LAND_WATER(LKS+1,i,j) - LAND_WATER(LKS,i,j) ) / ( LCZ(LKS) + LCZ(LKS+1) ) ) )
       if ( ATMOS_HYDROMETEOR_dry ) then
          LAND_QVEF(i,j) = 0.0_RP
       else
          LAND_QVEF(i,j) = min( total / LAND_PROPERTY(i,j,I_WaterCritical), BETA_MAX )
       end if

       ! eq.(12) in Merlin et al.(2011) but simplified P=0.5 used
       !sw = 0.5_RP + sign(0.5_RP,LAND_WATER(LKS,i,j)-LAND_PROPERTY(i,j,I_WaterCritical)) ! if W > Wc, sw = 1
       !LAND_QVEF(i,j) = (        sw ) * 1.0_RP &
       !               + ( 1.0_RP-sw ) * sqrt( 0.5_RP - 0.5_RP * cos( PI * LAND_WATER(LKS,i,j) / LAND_PROPERTY(i,j,I_WaterCritical) ) )

       LAND_TC_dZ(i,j) = LAND_PROPERTY(i,j,I_ThermalCond) / LCZ(LKS)
    end do
    end do


    !------------------------------------------------------------------------
    !> all land area without snow model or no snow area with snow model


    !$omp parallel do
    do j = LJS, LJE
    do i = LIS, LIE
       if ( LAND_ICE(LKS,i,j) > 0.0_RP ) then
          call HYDROMETEOR_LHS( LIA, LIS, LIE, LJA, LJS, LJE, &
                                LAND_SFC_TEMP(:,:), LH(:,:) )
       else
          call HYDROMETEOR_LHV( LIA, LIS, LIE, LJA, LJS, LJE, &
                                LAND_SFC_TEMP(:,:), LH(:,:) )
       end if
    end do
    end do


    select case ( LAND_SFC_TYPE )
    case ( 'SKIN' )
!OCL XFILL
       !$omp parallel do
       do j = LJS, LJE
       do i = LIS, LIE
          LAND_SFC_albedo(i,j,I_R_direct ,I_R_IR ) = LAND_PROPERTY(i,j,I_ALBLW)
          LAND_SFC_albedo(i,j,I_R_diffuse,I_R_IR ) = LAND_PROPERTY(i,j,I_ALBLW)
          LAND_SFC_albedo(i,j,I_R_direct ,I_R_NIR) = LAND_PROPERTY(i,j,I_ALBSW)
          LAND_SFC_albedo(i,j,I_R_diffuse,I_R_NIR) = LAND_PROPERTY(i,j,I_ALBSW)
          LAND_SFC_albedo(i,j,I_R_direct ,I_R_VIS) = LAND_PROPERTY(i,j,I_ALBSW)
          LAND_SFC_albedo(i,j,I_R_diffuse,I_R_VIS) = LAND_PROPERTY(i,j,I_ALBSW)
       end do
       end do

       call CPL_PHY_SFC_skin( LIA, LIS, LIE, LJA, LJS, LJE, &
                              ATMOS_TEMP(:,:), ATMOS_PRES(:,:),                        & ! [IN]
                              ATMOS_W(:,:), ATMOS_U(:,:), ATMOS_V(:,:),                & ! [IN]
                              ATMOS_DENS(:,:), ATMOS_QV(:,:),                          & ! [IN]
                              LH(:,:), REAL_Z1(:,:), ATMOS_PBL(:,:),                   & ! [IN]
                              ATMOS_SFC_DENS(:,:), ATMOS_SFC_PRES(:,:),                & ! [IN]
                              ATMOS_SFLX_rad_dn(:,:,:,:),                              & ! [IN]
                              LAND_TEMP(LKS,:,:), LAND_WSTR(:,:), LAND_QVEF(:,:),      & ! [IN]
                              LAND_SFC_albedo(:,:,:,:),                                & ! [IN]
                              LAND_PROPERTY(:,:,I_StomataResist),                      & ! [IN]
                              LAND_TC_dZ(:,:),                                         & ! [IN]
                              LAND_PROPERTY(:,:,I_Z0M),                                & ! [IN]
                              LAND_PROPERTY(:,:,I_Z0H),                                & ! [IN]
                              LAND_PROPERTY(:,:,I_Z0E),                                & ! [IN]
                              LANDUSE_exists_land(:,:), dt,                            & ! [IN]
                              'LAND',                                                  & ! [IN]
                              LAND_SFC_TEMP(:,:),                                      & ! [INOUT]
                              LAND_SFLX_MW(:,:), LAND_SFLX_MU(:,:), LAND_SFLX_MV(:,:), & ! [OUT]
                              LAND_SFLX_SH(:,:), LAND_SFLX_LH(:,:), SFLX_QV(:,:),      & ! [OUT]
                              LAND_SFLX_GH(:,:),                                       & ! [OUT]
                              SOIL_Ustar(:,:), SOIL_Tstar(:,:), SOIL_Qstar(:,:),       & ! [OUT]
                              SOIL_Wstar(:,:),                                         & ! [OUT]
                              SOIL_RLmo(:,:),                                          & ! [OUT]
                              LAND_U10(:,:), LAND_V10(:,:), LAND_T2(:,:), LAND_Q2(:,:) ) ! [OUT]

    case ( 'FIXED-TEMP' )
!OCL XFILL
       !$omp parallel do
       do j = LJS, LJE
       do i = LIS, LIE
          LAND_SFC_TEMP(i,j) = LAND_TEMP(LKS,i,j)
       end do
       end do
!OCL XFILL
       !$omp parallel do
       do j = LJS, LJE
       do i = LIS, LIE
          LAND_SFC_albedo(i,j,I_R_direct ,I_R_IR ) = LAND_PROPERTY(i,j,I_ALBLW)
          LAND_SFC_albedo(i,j,I_R_diffuse,I_R_IR ) = LAND_PROPERTY(i,j,I_ALBLW)
          LAND_SFC_albedo(i,j,I_R_direct ,I_R_NIR) = LAND_PROPERTY(i,j,I_ALBSW)
          LAND_SFC_albedo(i,j,I_R_diffuse,I_R_NIR) = LAND_PROPERTY(i,j,I_ALBSW)
          LAND_SFC_albedo(i,j,I_R_direct ,I_R_VIS) = LAND_PROPERTY(i,j,I_ALBSW)
          LAND_SFC_albedo(i,j,I_R_diffuse,I_R_VIS) = LAND_PROPERTY(i,j,I_ALBSW)
       end do
       end do

       call CPL_PHY_SFC_fixed_temp( LIA, LIS, LIE, LJA, LJS, LJE, &
                                    ATMOS_TEMP(:,:), ATMOS_PRES(:,:),                        & ! [IN]
                                    ATMOS_W(:,:), ATMOS_U(:,:), ATMOS_V(:,:),                & ! [IN]
                                    ATMOS_DENS(:,:), ATMOS_QV(:,:), LH(:,:),                 & ! [IN]
                                    REAL_Z1(:,:), ATMOS_PBL(:,:),                            & ! [IN]
                                    ATMOS_SFC_DENS(:,:), ATMOS_SFC_PRES(:,:),                & ! [IN]
                                    ATMOS_SFLX_rad_dn(:,:,:,:),                              & ! [IN]
                                    LAND_SFC_TEMP(:,:), LAND_WSTR(:,:), LAND_QVEF(:,:),      & ! [IN]
                                    LAND_SFC_albedo(:,:,:,:),                                & ! [IN]
                                    LAND_PROPERTY(:,:,I_StomataResist),                      & ! [IN]
                                    LAND_PROPERTY(:,:,I_Z0M),                                & ! [IN]
                                    LAND_PROPERTY(:,:,I_Z0H),                                & ! [IN]
                                    LAND_PROPERTY(:,:,I_Z0E),                                & ! [IN]
                                    LANDUSE_exists_land(:,:), dt,                            & ! [IN]
                                    LAND_SFLX_MW(:,:), LAND_SFLX_MU(:,:), LAND_SFLX_MV(:,:), & ! [OUT]
                                    LAND_SFLX_SH(:,:), LAND_SFLX_LH(:,:), SFLX_QV(:,:),      & ! [OUT]
                                    LAND_SFLX_GH(:,:),                                       & ! [OUT]
                                    SOIL_Ustar(:,:), SOIL_Tstar(:,:), SOIL_Qstar(:,:),       & ! [OUT]
                                    SOIL_Wstar(:,:),                                         & ! [OUT]
                                    SOIL_RLmo(:,:),                                          & ! [OUT]
                                    LAND_U10(:,:), LAND_V10(:,:),                            & ! [OUT]
                                    LAND_T2(:,:), LAND_Q2(:,:)                               ) ! [OUT]
    end select

    !$omp parallel do
    do j = LJS, LJE
    do i = LIS, LIE
       if ( LAND_ICE(LKS,i,j) > 0.0_RP ) then
          SFLX_ENGI(i,j) = ( CV_ICE * LAND_SFC_TEMP(i,j) - LHF ) * SFLX_QV(i,j)
       else
          SFLX_ENGI(i,j) = CV_WATER * LAND_SFC_TEMP(i,j) * SFLX_QV(i,j)
       end if
    end do
    end do

    ! LAND_SFLX_* are positive for downward
!OCL XFILL
    !$omp parallel do
    do j = LJS, LJE
    do i = LIS, LIE
       LAND_SFLX_water(i,j) = ATMOS_SFLX_water(i,j) - SFLX_QV(i,j)
       LAND_SFLX_ENGI (i,j) = ATMOS_SFLX_ENGI(i,j) & ! internal energy of precipitation
                            - SFLX_ENGI(i,j)         ! internal energy of evapolation or sublimation
    end do
    end do

    if ( SNOW_flag ) then

       call FILE_HISTORY_in( LAND_SFLX_MW(:,:), 'SOIL_SFLX_MW',  'soil surface w-momentum flux (upward)',    'kg/m2/s' )
       call FILE_HISTORY_in( LAND_SFLX_MU(:,:), 'SOIL_SFLX_MU',  'soil surface u-momentum flux (upward)',    'kg/m2/s' )
       call FILE_HISTORY_in( LAND_SFLX_MV(:,:), 'SOIL_SFLX_MV',  'soil surface v-momentum flux (upward)',    'kg/m2/s' )
       call FILE_HISTORY_in( LAND_SFLX_SH(:,:), 'SOIL_SFLX_SH',  'soil surface sensible heat flux (upward)', 'J/m2/s'  )
       call FILE_HISTORY_in( LAND_SFLX_LH(:,:), 'SOIL_SFLX_LH',  'soil surface latent heat flux (upward)',   'J/m2/s'  )
       call FILE_HISTORY_in( LAND_U10    (:,:), 'LAND_SOIL_U10', 'Wind velocity u at 10 m on soil surface',  'm/s'     )
       call FILE_HISTORY_in( LAND_V10    (:,:), 'LAND_SOIL_V10', 'Wind velocity v at 10 m on soil surface',  'm/s'     )
       call FILE_HISTORY_in( LAND_T2     (:,:), 'LAND_SOIL_T2',  'Air temperature at 2m on soil surface',    'K'       )
       call FILE_HISTORY_in( LAND_Q2     (:,:), 'LAND_SOIL_Q2',  'Specific humidity at 2m on soil surface',  'kg/kg'   )

       ! marge land surface and snow surface !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!OCL XFILL
       !$omp parallel do
       do j = LJS, LJE
       do i = LIS, LIE
          LAND_SFC_TEMP(i,j) = (        SNOW_frac(i,j) ) * SNOW_SFC_TEMP(i,j) &
                             + ( 1.0_RP-SNOW_frac(i,j) ) * LAND_SFC_TEMP(i,j)

          do idir = I_R_direct, I_R_diffuse
             LAND_SFC_albedo(i,j,idir,I_R_IR ) = (        SNOW_frac(i,j) ) * SNOW_albedo    (i,j,I_LW)         &
                                               + ( 1.0_RP-SNOW_frac(i,j) ) * LAND_SFC_albedo(i,j,idir,I_R_IR)
             LAND_SFC_albedo(i,j,idir,I_R_NIR) = (        SNOW_frac(i,j) ) * SNOW_albedo    (i,j,I_SW)         &
                                               + ( 1.0_RP-SNOW_frac(i,j) ) * LAND_SFC_albedo(i,j,idir,I_R_NIR)
             LAND_SFC_albedo(i,j,idir,I_R_VIS) = (        SNOW_frac(i,j) ) * SNOW_albedo    (i,j,I_SW)         &
                                               + ( 1.0_RP-SNOW_frac(i,j) ) * LAND_SFC_albedo(i,j,idir,I_R_VIS)
          enddo

          ! flux to the soil
          LAND_SFLX_GH   (i,j) = (        SNOW_frac(i,j) ) * SNOW_LAND_SFLX_GH   (i,j) &
                               + ( 1.0_RP-SNOW_frac(i,j) ) *      LAND_SFLX_GH   (i,j)
          LAND_SFLX_water(i,j) = (        SNOW_frac(i,j) ) * SNOW_LAND_SFLX_water(i,j) &
                               + ( 1.0_RP-SNOW_frac(i,j) ) *      LAND_SFLX_water(i,j)
          LAND_SFLX_ENGI (i,j) = (        SNOW_frac(i,j) ) * SNOW_LAND_SFLX_ENGI (i,j) &
                               + ( 1.0_RP-SNOW_frac(i,j) ) *      LAND_SFLX_ENGI (i,j)
          ! flux to the atmosphere
          LAND_SFLX_MW(i,j) = (        SNOW_frac(i,j) ) * SNOW_ATMOS_SFLX_MW(i,j) &
                            + ( 1.0_RP-SNOW_frac(i,j) ) *       LAND_SFLX_MW(i,j)
          LAND_SFLX_MU(i,j) = (        SNOW_frac(i,j) ) * SNOW_ATMOS_SFLX_MU(i,j) &
                            + ( 1.0_RP-SNOW_frac(i,j) ) *       LAND_SFLX_MU(i,j)
          LAND_SFLX_MV(i,j) = (        SNOW_frac(i,j) ) * SNOW_ATMOS_SFLX_MV(i,j) &
                            + ( 1.0_RP-SNOW_frac(i,j) ) *       LAND_SFLX_MV(i,j)
          LAND_SFLX_SH(i,j) = (        SNOW_frac(i,j) ) * SNOW_ATMOS_SFLX_SH(i,j) &
                            + ( 1.0_RP-SNOW_frac(i,j) ) *       LAND_SFLX_SH(i,j)
          LAND_SFLX_LH(i,j) = (        SNOW_frac(i,j) ) * SNOW_ATMOS_SFLX_LH(i,j) &
                            + ( 1.0_RP-SNOW_frac(i,j) ) *       LAND_SFLX_LH(i,j)
               SFLX_QV(i,j) = (        SNOW_frac(i,j) ) * SNOW_ATMOS_SFLX_QV(i,j) &
                            + ( 1.0_RP-SNOW_frac(i,j) ) *            SFLX_QV(i,j)
          ! diagnostics
          LAND_U10(i,j) = (        SNOW_frac(i,j) ) * SNOW_U10(i,j) &
                        + ( 1.0_RP-SNOW_frac(i,j) ) * LAND_U10(i,j)
          LAND_V10(i,j) = (        SNOW_frac(i,j) ) * SNOW_V10(i,j) &
                        + ( 1.0_RP-SNOW_frac(i,j) ) * LAND_V10(i,j)
          LAND_T2 (i,j) = (        SNOW_frac(i,j) ) * SNOW_T2 (i,j) &
                        + ( 1.0_RP-SNOW_frac(i,j) ) * LAND_T2 (i,j)
          LAND_Q2 (i,j) = (        SNOW_frac(i,j) ) * SNOW_Q2 (i,j) &
                        + ( 1.0_RP-SNOW_frac(i,j) ) * LAND_Q2 (i,j)

       enddo
       enddo

       call BULKFLUX_diagnose( LIA, LIS, LIE, LJA, LJS, LJE, &
                               LAND_SFLX_MW(:,:), LAND_SFLX_MU(:,:), LAND_SFLX_MV(:,:), & ! [IN]
                               LAND_SFLX_SH(:,:), SFLX_QV(:,:),                         & ! [IN]
                               ATMOS_SFC_DENS(:,:), LAND_SFC_TEMP(:,:), ATMOS_PBL(:,:), & ! [IN]
                               LAND_Ustar(:,:), LAND_Tstar(:,:), LAND_Qstar(:,:),       & ! [OUT]
                               LAND_Wstar(:,:), LAND_RLmo(:,:)                          ) ! [OUT]

    end if

    if ( .NOT. ATMOS_HYDROMETEOR_dry ) then
       !$omp parallel do
       do j = LJS, LJE
       do i = LIS, LIE
          LAND_SFLX_QTRC(i,j,I_QV) = SFLX_QV(i,j)
       enddo
       enddo
    end if


    ! Surface flux for chemical tracers
    if ( ATMOS_sw_phy_ch ) then
       call ATMOS_PHY_CH_driver_LAND_flux( LAND_SFLX_QTRC(:,:,:) ) ! [INOUT]
    endif

    !########## Set Surface Boundary to coupler ##########
    call LAND_SURFACE_SET( countup=.true. )

    call PROF_rapend  ('LND_CalcTend', 1)

    return
  end subroutine LAND_driver_calc_tendency

  !-----------------------------------------------------------------------------
  !> Land step
  subroutine LAND_driver_update
    use scale_time, only: &
       dt => TIME_DTSEC_LAND
    use mod_land_vars, only: &
       LAND_PROPERTY,     &
       I_WaterLimit,      &
       I_ThermalCond,     &
       I_HeatCapacity,    &
       I_WaterDiff,       &
       LAND_TEMP,         &
       LAND_WATER,        &
       LAND_ICE,          &
       LAND_SFLX_GH,      &
       LAND_SFLX_water,   &
       LAND_SFLX_ENGI,    &
       LAND_RUNOFF,       &
       LAND_RUNOFF_ENGI,  &
       LAND_SFC_TEMP,     &
       LAND_TEMP_t,       &
       LAND_WATER_t,      &
       LAND_ICE_t,        &
       LAND_vars_total
    use scale_land_grid_cartesC, only: &
       LCDZ => LAND_GRID_CARTESC_CDZ
    use scale_land_dyn_bucket, only: &
       LAND_DYN_bucket
    use scale_landuse, only: &
       LANDUSE_fact_land
    use scale_time, only: &
       NOWDAYSEC => TIME_NOWDAYSEC
    use scale_file_history, only: &
       FILE_HISTORY_in
    use mod_land_admin, only: &
       LAND_DYN_TYPE
    implicit none

    integer :: k, i, j
    !---------------------------------------------------------------------------

    call PROF_rapstart('LND_Update', 2)

    !########## Get Surface Boundary from coupler ##########
    call LAND_SURFACE_GET

    !########## Dynamics / Update variables ##########
    select case ( LAND_DYN_TYPE )
    case ( 'BUCKET' )
       call LAND_DYN_bucket( LKMAX, LKS, LKE, LIA, LIS, LIE, LJA, LJS, LJE, &
                             LAND_TEMP_t(:,:,:),                     & ! [IN]
                             LAND_WATER_t(:,:,:), LAND_ICE_t(:,:,:), & ! [IN]
                             LAND_PROPERTY(:,:,I_WaterLimit),        & ! [IN]
                             LAND_PROPERTY(:,:,I_ThermalCond),       & ! [IN]
                             LAND_PROPERTY(:,:,I_HeatCapacity),      & ! [IN]
                             LAND_PROPERTY(:,:,I_WaterDiff),         & ! [IN]
                             LAND_SFLX_GH(:,:),                      & ! [IN]
                             LAND_SFLX_water(:,:),                   & ! [IN]
                             LAND_SFLX_ENGI(:,:),                    & ! [IN]
                             LANDUSE_fact_land(:,:), LCDZ(:),        & ! [IN]
                             dt, NOWDAYSEC,                          & ! [IN]
                             LAND_TEMP(:,:,:),                       & ! [INOUT]
                             LAND_WATER(:,:,:), LAND_ICE(:,:,:),     & ! [INOUT]
                             LAND_RUNOFF(:,:), LAND_RUNOFF_ENGI(:,:) ) ! [OUT]
    case ( 'INIT' )
       ! Never update LAND_TEMP and LAND_WATER from initial condition
    end select

    !########## Negative Fixer ##########
    !$omp parallel do
    do j = LJS, LJE
    do i = LIS, LIE
    do k = LKS, LKE
       LAND_WATER(k,i,j) = max( LAND_WATER(k,i,j), 0.0_RP )
       LAND_ICE  (k,i,j) = max( LAND_ICE  (k,i,j), 0.0_RP )
    enddo
    enddo
    enddo

    call LAND_vars_total

    call PROF_rapend  ('LND_Update', 1)

    return
  end subroutine LAND_driver_update

  !-----------------------------------------------------------------------------
  !> Get surface boundary from other model
  subroutine LAND_SURFACE_GET
    use mod_land_admin, only: &
       LAND_do
    use mod_land_vars, only: &
       ATMOS_TEMP,        &
       ATMOS_PRES,        &
       ATMOS_W,           &
       ATMOS_U,           &
       ATMOS_V,           &
       ATMOS_DENS,        &
       ATMOS_QV,          &
       ATMOS_PBL,         &
       ATMOS_SFC_DENS,    &
       ATMOS_SFC_PRES,    &
       ATMOS_SFLX_rad_dn, &
       ATMOS_cosSZA,      &
       ATMOS_SFLX_water,  &
       ATMOS_SFLX_ENGI
    use mod_cpl_vars, only: &
       CPL_getATM_LND
    implicit none
    !---------------------------------------------------------------------------

    call PROF_rapstart('LND_SfcExch', 2)

    if ( LAND_do ) then
       call CPL_getATM_LND( ATMOS_TEMP       (:,:),     & ! [OUT]
                            ATMOS_PRES       (:,:),     & ! [OUT]
                            ATMOS_W          (:,:),     & ! [OUT]
                            ATMOS_U          (:,:),     & ! [OUT]
                            ATMOS_V          (:,:),     & ! [OUT]
                            ATMOS_DENS       (:,:),     & ! [OUT]
                            ATMOS_QV         (:,:),     & ! [OUT]
                            ATMOS_PBL        (:,:),     & ! [OUT]
                            ATMOS_SFC_DENS   (:,:),     & ! [OUT]
                            ATMOS_SFC_PRES   (:,:),     & ! [OUT]
                            ATMOS_SFLX_rad_dn(:,:,:,:), & ! [OUT]
                            ATMOS_cosSZA     (:,:),     & ! [OUT]
                            ATMOS_SFLX_water (:,:),     & ! [OUT]
                            ATMOS_SFLX_ENGI  (:,:)      ) ! [OUT]
    endif

    call PROF_rapend  ('LND_SfcExch', 2)

    return
  end subroutine LAND_SURFACE_GET

  !-----------------------------------------------------------------------------
  !> Put surface boundary to other model
  subroutine LAND_SURFACE_SET( countup )
    use mod_land_admin, only: &
       LAND_do
    use mod_land_vars, only: &
       LAND_PROPERTY,   &
       I_Z0M,           &
       I_Z0H,           &
       I_Z0E,           &
       LAND_SFC_TEMP,   &
       LAND_SFC_albedo, &
       LAND_SFLX_MW,    &
       LAND_SFLX_MU,    &
       LAND_SFLX_MV,    &
       LAND_SFLX_SH,    &
       LAND_SFLX_LH,    &
       LAND_SFLX_GH,    &
       LAND_SFLX_QTRC,  &
       LAND_U10,        &
       LAND_V10,        &
       LAND_T2,         &
       LAND_Q2
    use mod_cpl_vars, only: &
       CPL_putLND
    implicit none

    ! arguments
    logical, intent(in) :: countup
    !---------------------------------------------------------------------------

    call PROF_rapstart('LND_SfcExch', 2)

    if ( LAND_do ) then
       call CPL_putLND( LAND_SFC_TEMP  (:,:),       & ! [IN]
                        LAND_SFC_albedo(:,:,:,:),   & ! [IN]
                        LAND_PROPERTY  (:,:,I_Z0M), & ! [IN]
                        LAND_PROPERTY  (:,:,I_Z0H), & ! [IN]
                        LAND_PROPERTY  (:,:,I_Z0E), & ! [IN]
                        LAND_SFLX_MW   (:,:),       & ! [IN]
                        LAND_SFLX_MU   (:,:),       & ! [IN]
                        LAND_SFLX_MV   (:,:),       & ! [IN]
                        LAND_SFLX_SH   (:,:),       & ! [IN]
                        LAND_SFLX_LH   (:,:),       & ! [IN]
                        LAND_SFLX_GH   (:,:),       & ! [IN]
                        LAND_SFLX_QTRC (:,:,:),     & ! [IN]
                        LAND_U10       (:,:),       & ! [IN]
                        LAND_V10       (:,:),       & ! [IN]
                        LAND_T2        (:,:),       & ! [IN]
                        LAND_Q2        (:,:),       & ! [IN]
                        countup                     ) ! [IN]
    endif

    call PROF_rapend  ('LND_SfcExch', 2)

    return
  end subroutine LAND_SURFACE_SET

end module mod_land_driver
