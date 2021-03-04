!-------------------------------------------------------------------------------
!> module OCEAN driver
!!
!! @par Description
!!          Ocean model driver
!!
!! @author Team SCALE
!<
!-------------------------------------------------------------------------------
#include "scalelib.h"
module mod_ocean_driver
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_io
  use scale_prof
  use scale_ocean_grid_cartesC_index
  use scale_tracer
  use scale_cpl_sfc_index
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: OCEAN_driver_setup
  public :: OCEAN_driver_finalize
  public :: OCEAN_driver_calc_tendency
  public :: OCEAN_driver_update
  public :: OCEAN_SURFACE_GET
  public :: OCEAN_SURFACE_SET

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
  real(RP), private, allocatable :: WSTR   (:,:)
  real(RP), private, allocatable :: QVEF   (:,:)
  real(RP), private, allocatable :: SR     (:,:)
  real(RP), private, allocatable :: ATMOS_W(:,:)

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine OCEAN_driver_setup
    use scale_prc, only: &
       PRC_abort
    use scale_const, only: &
       HUGE => CONST_HUGE
    use mod_ocean_admin, only: &
       OCEAN_do,       &
       OCEAN_DYN_TYPE, &
       OCEAN_SFC_TYPE, &
       OCEAN_ICE_TYPE, &
       OCEAN_ALB_TYPE, &
       OCEAN_RGN_TYPE
    use scale_cpl_phy_sfc_fixed_temp, only: &
       CPL_PHY_SFC_fixed_temp_setup
    use scale_atmos_hydrometeor, only: &
       ATMOS_HYDROMETEOR_dry
    use scale_ocean_dyn_slab, only: &
       OCEAN_DYN_SLAB_setup
    use scale_ocean_dyn_offline, only: &
       OCEAN_DYN_OFFLINE_setup
    use scale_ocean_phy_ice_simple, only: &
       OCEAN_PHY_ICE_setup
    use scale_ocean_phy_albedo, only: &
       OCEAN_PHY_ALBEDO_const_setup, &
       OCEAN_PHY_ALBEDO_seaice_setup
    use scale_ocean_phy_albedo_nakajima00, only: &
       OCEAN_PHY_ALBEDO_nakajima00_setup
    use scale_ocean_phy_roughness, only: &
       OCEAN_PHY_ROUGHNESS_const_setup, &
       OCEAN_PHY_ROUGHNESS_seaice_setup
    use scale_ocean_phy_roughness_miller92, only: &
       OCEAN_PHY_ROUGHNESS_miller92_setup
    use scale_ocean_phy_roughness_moon07, only: &
       OCEAN_PHY_ROUGHNESS_moon07_setup
    use scale_ocean_phy_tc, only: &
       OCEAN_PHY_TC_seaice_setup
    use scale_ocean_grid_cartesC, only: &
       CDZ => OCEAN_GRID_CARTESC_CDZ
    implicit none
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("OCEAN_driver_setup",*) 'Setup'

    if ( OCEAN_do ) then

       select case ( OCEAN_DYN_TYPE )
       case ( 'SLAB' )
          call OCEAN_DYN_SLAB_setup( CDZ(OKS) )
       case ( 'OFFLINE' )
          call OCEAN_DYN_OFFLINE_setup
       case ( 'INIT' )
          ! do nothing
       case default
          LOG_ERROR("OCEAN_driver_setup",*) 'OCEAN_DYN_TYPE is invalid: ', trim(OCEAN_DYN_TYPE)
          call PRC_abort
       end select

       select case ( OCEAN_SFC_TYPE )
       case ( 'FIXED-TEMP' )
          call CPL_PHY_SFC_fixed_temp_setup
       case default
          LOG_ERROR("OCEAN_driver_setup",*) 'OCEAN_SFC_TYPE is invalid: ', trim(OCEAN_SFC_TYPE)
          call PRC_abort
       end select

       select case ( OCEAN_ICE_TYPE )
       case ( 'NONE' )
       case ( 'SIMPLE' )
          call OCEAN_PHY_ICE_setup
       case ( 'INIT' )
       case default
          LOG_ERROR("OCEAN_driver_setup",*) 'OCEAN_ICE_TYPE is invalid: ', trim(OCEAN_ICE_TYPE)
          call PRC_abort
       end select

       ! surface albedo
       select case ( OCEAN_ALB_TYPE )
       case ( 'NAKAJIMA00' )
          call OCEAN_PHY_ALBEDO_const_setup
          call OCEAN_PHY_ALBEDO_nakajima00_setup
          call OCEAN_PHY_ALBEDO_seaice_setup
       case ( 'CONST' )
          call OCEAN_PHY_ALBEDO_const_setup
       case ( 'INIT' )
          ! do nothing
       case default
          LOG_ERROR("OCEAN_driver_setup",*) 'OCEAN_ALB_TYPE is invalid: ', trim(OCEAN_ALB_TYPE)
          call PRC_abort
       end select

       ! surface roughness length
       select case ( OCEAN_RGN_TYPE )
       case ( 'MILLER92' )
          call OCEAN_PHY_ROUGHNESS_miller92_setup
          call OCEAN_PHY_ROUGHNESS_seaice_setup
       case ( 'MOON07' )
          call OCEAN_PHY_ROUGHNESS_moon07_setup
          call OCEAN_PHY_ROUGHNESS_seaice_setup
       case ( 'CONST' )
          call OCEAN_PHY_ROUGHNESS_const_setup
       case ( 'INIT' )
          ! do nothing
       case default
          LOG_ERROR("OCEAN_driver_setup",*) 'OCEAN_RGN_TYPE is invalid: ', trim(OCEAN_RGN_TYPE)
          call PRC_abort
       end select

       ! thermal conductivity
       call OCEAN_PHY_TC_seaice_setup

       allocate( WSTR   (OIA,OJA) )
       allocate( QVEF   (OIA,OJA) )
       allocate( SR     (OIA,OJA) )
       allocate( ATMOS_W(OIA,OJA) )
       WSTR(:,:) = HUGE
       if ( ATMOS_HYDROMETEOR_dry ) then
          QVEF   (:,:) = 0.0_RP
       else
          QVEF   (:,:) = 1.0_RP
       end if
       SR     (:,:) = 0.0_RP
       ATMOS_W(:,:) = 0.0_RP ! slope of the sea surface is zero

    endif

    return
  end subroutine OCEAN_driver_setup

  !-----------------------------------------------------------------------------
  !> Finalize
  subroutine OCEAN_driver_finalize
    use mod_ocean_admin, only: &
       OCEAN_do,       &
       OCEAN_DYN_TYPE, &
       OCEAN_ICE_TYPE, &
       OCEAN_ALB_TYPE, &
       OCEAN_RGN_TYPE
    implicit none
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("OCEAN_driver_finalize",*) 'Finalize'

    if ( OCEAN_do ) then

       select case ( OCEAN_DYN_TYPE )
       case ( 'SLAB' )
       case ( 'OFFLINE' )
       case ( 'INIT' )
       end select

       select case ( OCEAN_ICE_TYPE )
       case ( 'NONE' )
       case ( 'SIMPLE' )
       case ( 'INIT' )
       end select

       ! surface albedo
       select case ( OCEAN_ALB_TYPE )
       case ( 'NAKAJIMA00' )
       case ( 'CONST' )
       case ( 'INIT' )
       end select

       ! surface roughness length
       select case ( OCEAN_RGN_TYPE )
       case ( 'MILLER92' )
       case ( 'MOON07' )
       case ( 'CONST' )
       case ( 'INIT' )
       end select

       deallocate( WSTR    )
       deallocate( QVEF    )
       deallocate( SR      )
       deallocate( ATMOS_W )

    endif

    return
  end subroutine OCEAN_driver_finalize

  !-----------------------------------------------------------------------------
  !> Calculate tendency
  subroutine OCEAN_driver_calc_tendency( force )
    use scale_const, only: &
       EPS => CONST_EPS
    use scale_time, only: &
       dt => TIME_DTSEC_OCEAN
    use scale_statistics, only: &
       STATISTICS_checktotal, &
       STATISTICS_total
    use scale_landuse, only: &
       exists_ocean => LANDUSE_exists_ocean
    use scale_atmos_grid_cartesC_real, only: &
       REAL_Z1 => ATMOS_GRID_CARTESC_REAL_Z1
    use scale_atmos_hydrometeor, only: &
       HYDROMETEOR_LHV => ATMOS_HYDROMETEOR_LHV, &
       HYDROMETEOR_LHS => ATMOS_HYDROMETEOR_LHS, &
       ATMOS_HYDROMETEOR_dry, &
       CV_WATER,              &
       CV_ICE,                &
       LHF,                   &
       I_QV
    use scale_ocean_grid_cartesC_real, only: &
       OCEAN_GRID_CARTESC_REAL_VOL,    &
       OCEAN_GRID_CARTESC_REAL_TOTVOL, &
       OCEAN_GRID_CARTESC_REAL_AREA,   &
       OCEAN_GRID_CARTESC_REAL_TOTAREA
    use scale_ocean_grid_cartesC_real, only: &
       OCEAN_GRID_CARTESC_REAL_AREA, &
       OCEAN_GRID_CARTESC_REAL_TOTAREA
    use scale_cpl_phy_sfc_fixed_temp, only: &
       CPL_PHY_SFC_fixed_temp
    use scale_ocean_phy_albedo, only: &
       OCEAN_PHY_ALBEDO_const, &
       OCEAN_PHY_ALBEDO_seaice
    use scale_ocean_phy_albedo_nakajima00, only: &
       OCEAN_PHY_ALBEDO_nakajima00
    use scale_OCEAN_PHY_roughness, only: &
       OCEAN_PHY_ROUGHNESS_const, &
       OCEAN_PHY_ROUGHNESS_seaice
    use scale_ocean_phy_roughness_miller92, only: &
       OCEAN_PHY_ROUGHNESS_miller92
    use scale_ocean_phy_roughness_moon07, only: &
       OCEAN_PHY_ROUGHNESS_moon07
    use scale_ocean_phy_tc, only: &
       OCEAN_PHY_TC_seaice
    use scale_ocean_phy_ice_simple, only: &
       OCEAN_PHY_ICE_simple
    use scale_bulkflux, only: &
       BULKFLUX_diagnose_scales
    use mod_atmos_admin, only: &
       ATMOS_sw_phy_ch
    use mod_atmos_phy_ch_driver, only: &
       ATMOS_PHY_CH_driver_OCEAN_flux
    use mod_ocean_admin, only: &
       OCEAN_SFC_TYPE, &
       OCEAN_ICE_TYPE, &
       OCEAN_ALB_TYPE, &
       OCEAN_RGN_TYPE
    use mod_ocean_vars, only: &
       ICE_FLAG,          &
       OCEAN_TEMP,        &
       OCEAN_OCN_Z0M,     &
       OCEAN_ICE_TEMP,    &
       OCEAN_ICE_MASS,    &
       OCEAN_SFC_TEMP,    &
       OCEAN_SFC_albedo,  &
       OCEAN_SFC_Z0M,     &
       OCEAN_SFC_Z0H,     &
       OCEAN_SFC_Z0E,     &
       OCEAN_TEMP_t,      &
       OCEAN_SALT_t,      &
       OCEAN_UVEL_t,      &
       OCEAN_VVEL_t,      &
       OCEAN_ICE_TEMP_t,  &
       OCEAN_ICE_MASS_t,  &
       ATMOS_TEMP,        &
       ATMOS_PRES,        &
       ATMOS_U,           &
       ATMOS_V,           &
       ATMOS_DENS,        &
       ATMOS_QV,          &
       ATMOS_PBL,         &
       ATMOS_cosSZA,      &
       ATMOS_SFC_DENS,    &
       ATMOS_SFC_PRES,    &
       ATMOS_SFLX_rad_dn, &
       ATMOS_SFLX_water,  &
       ATMOS_SFLX_ENGI,   &
       OCEAN_SFLX_MW,     &
       OCEAN_SFLX_MU,     &
       OCEAN_SFLX_MV,     &
       OCEAN_SFLX_SH,     &
       OCEAN_SFLX_LH,     &
       OCEAN_SFLX_QTRC,   &
       OCEAN_U10,         &
       OCEAN_V10,         &
       OCEAN_T2,          &
       OCEAN_Q2,          &
       OCEAN_Ustar,       &
       OCEAN_Tstar,       &
       OCEAN_Qstar,       &
       OCEAN_Wstar,       &
       OCEAN_RLmo,        &
       OCEAN_OCN_Ustar,   &
       OCEAN_OCN_Tstar,   &
       OCEAN_OCN_Qstar,   &
       OCEAN_OCN_Wstar,   &
       OCEAN_OCN_RLmo,    &
       OCEAN_ICE_Ustar,   &
       OCEAN_ICE_Tstar,   &
       OCEAN_ICE_Qstar,   &
       OCEAN_ICE_Wstar,   &
       OCEAN_ICE_RLmo,    &
       OCEAN_SFLX_GH,     &
       OCEAN_SFLX_water,  &
       OCEAN_SFLX_ENGI,   &
       OCEAN_OFLX_GH,     &
       OCEAN_OFLX_water,  &
       OCEAN_OFLX_ENGI,   &
       OCEAN_ICE_FRAC
    use scale_file_history, only: &
       FILE_HISTORY_in
    implicit none

    logical, intent(in) :: force

    real(RP) :: LHV          (OIA,OJA)
    real(RP) :: LHS          (OIA,OJA)
    real(RP) :: ATMOS_Uabs   (OIA,OJA)
    real(RP) :: sfc_temp     (OIA,OJA)
    real(RP) :: sfc_albedo   (OIA,OJA,N_RAD_DIR,N_RAD_RGN)
    real(RP) :: sfc_Z0M      (OIA,OJA)
    real(RP) :: sfc_Z0H      (OIA,OJA)
    real(RP) :: sfc_Z0E      (OIA,OJA)
    real(RP) :: subsfc_temp  (OIA,OJA)
    real(RP) :: TC_dz        (OIA,OJA)
    real(RP) :: sflx_MW      (OIA,OJA)
    real(RP) :: sflx_MU      (OIA,OJA)
    real(RP) :: sflx_MV      (OIA,OJA)
    real(RP) :: sflx_SH      (OIA,OJA)
    real(RP) :: sflx_LH      (OIA,OJA)
    real(RP) :: sflx_QV      (OIA,OJA)
    real(RP) :: OCEAN_SFLX_QV(OIA,OJA)
    real(RP) :: U10          (OIA,OJA)
    real(RP) :: V10          (OIA,OJA)
    real(RP) :: T2           (OIA,OJA)
    real(RP) :: Q2           (OIA,OJA)
    real(RP) :: sflx_hbalance(OIA,OJA)
    real(RP) :: sflx_GH      (OIA,OJA)
    real(RP) :: sflx_water   (OIA,OJA)
    real(RP) :: sflx_engi    (OIA,OJA)
    real(RP) :: ice_mass     (OIA,OJA)
    logical  :: exists_ice   (OIA,OJA)
    real(RP) :: sw

    real(RP) :: sfc_frac

    integer  :: k, i, j, iq, idir, irgn
    !---------------------------------------------------------------------------

    call PROF_rapstart('OCN_CalcTend', 1)

    !########## Get Surface Boundary from coupler ##########
    call OCEAN_SURFACE_GET

    !$omp parallel do
    do j = OJS, OJE
    do i = OIS, OIE
       ATMOS_Uabs(i,j) = sqrt( ATMOS_U(i,j)**2 + ATMOS_V(i,j)**2 )
    enddo
    enddo

    !$omp parallel do
    do j = OJS, OJE
    do i = OIS, OIE
       exists_ice(i,j) = .false.
       if( exists_ocean(i,j) .AND. OCEAN_ICE_FRAC(i,j) > 0.0_RP ) exists_ice(i,j) = .true.
    enddo
    enddo

    !########## reset tendencies ##########

    !$omp parallel do
    do j = OJS, OJE
    do i = OIS, OIE
    do k = OKS, OKE
       OCEAN_TEMP_t(k,i,j) = 0.0_RP
       OCEAN_SALT_t(k,i,j) = 0.0_RP
       OCEAN_UVEL_t(k,i,j) = 0.0_RP
       OCEAN_VVEL_t(k,i,j) = 0.0_RP
    enddo
    enddo
    enddo

    do iq = 1, QA
    !$omp parallel do
    do j  = OJS, OJE
    do i  = OIS, OIE
       OCEAN_SFLX_QTRC(i,j,iq) = 0.0_RP
    enddo
    enddo
    enddo

    if ( ICE_flag ) then
       !$omp parallel do
       do j = OJS, OJE
       do i = OIS, OIE
          OCEAN_ICE_TEMP_t(i,j) = 0.0_RP
          OCEAN_ICE_MASS_t(i,j) = 0.0_RP
       enddo
       enddo
    end if



    !########## surface process (ice-free ocean) ##########

    !$omp parallel do
    do j = OJS, OJE
    do i = OIS, OIE
       OCEAN_SFC_TEMP(i,j) = OCEAN_TEMP(OKS,i,j)
    enddo
    enddo

    ! albedo
    select case ( OCEAN_ALB_TYPE )
    case ( 'NAKAJIMA00' )
       ! for Near-IR, IR
       call OCEAN_PHY_ALBEDO_const     ( OIA, OIS, OIE,        & ! [IN]
                                         OJA, OJS, OJE,        & ! [IN]
                                         sfc_albedo  (:,:,:,:) ) ! [OUT]
       ! for VIS (overwrite)
       call OCEAN_PHY_ALBEDO_nakajima00( OIA, OIS, OIE,              & ! [IN]
                                         OJA, OJS, OJE,              & ! [IN]
                                         ATMOS_cosSZA(:,:),          & ! [IN]
                                         sfc_albedo  (:,:,:,I_R_VIS) ) ! [OUT]
    case ( 'CONST' )
       call OCEAN_PHY_ALBEDO_const     ( OIA, OIS, OIE,        & ! [IN]
                                         OJA, OJS, OJE,        & ! [IN]
                                         sfc_albedo  (:,:,:,:) ) ! [OUT]
    case ( 'INIT' )
       ! Never update OCEAN_SFC_albedo from initial condition
       do irgn = I_R_IR, I_R_VIS
       do idir = I_R_direct, I_R_diffuse
       do j    = OJS, OJE
       do i    = OIS, OIE
       if ( exists_ocean(i,j) ) then
          sfc_albedo(i,j,idir,irgn) = OCEAN_SFC_albedo(i,j,idir,irgn)
       end if
       enddo
       enddo
       enddo
       enddo
    end select

    ! roughness length
    select case ( OCEAN_RGN_TYPE )
    case ( 'MILLER92' )
       call OCEAN_PHY_ROUGHNESS_miller92( OIA, OIS, OIE, & ! [IN]
                                          OJA, OJS, OJE, & ! [IN]
                                          ATMOS_Uabs(:,:),    & ! [IN]
                                          OCEAN_OCN_Z0M(:,:), & ! [OUT]
                                          OCEAN_SFC_Z0H(:,:), & ! [OUT]
                                          OCEAN_SFC_Z0E(:,:)  ) ! [OUT]
    case ( 'MOON07' )
       call OCEAN_PHY_ROUGHNESS_moon07  ( OIA, OIS, OIE, & ! [IN]
                                          OJA, OJS, OJE, & ! [IN]
                                          ATMOS_Uabs   (:,:), & ! [IN]
                                          REAL_Z1      (:,:), & ! [IN]
                                          exists_ocean (:,:), & ! [IN]
                                          OCEAN_OCN_Z0M(:,:), & ! [INOUT]
                                          OCEAN_SFC_Z0H(:,:), & ! [OUT]
                                          OCEAN_SFC_Z0E(:,:)  ) ! [OUT]
    case ( 'CONST' )
       call OCEAN_PHY_ROUGHNESS_const   ( OIA, OIS, OIE, & ! [IN]
                                          OJA, OJS, OJE, & ! [IN]
                                          OCEAN_OCN_Z0M(:,:), & ! [OUT]
                                          OCEAN_SFC_Z0H(:,:), & ! [OUT]
                                          OCEAN_SFC_Z0E(:,:)  ) ! [OUT]
    case ( 'INIT' )
       ! Never update from initial condition
    end select

    ! tendency
    select case ( OCEAN_SFC_TYPE )
    case ( 'FIXED-TEMP' )

       call HYDROMETEOR_LHV( OIA, OIS, OIE, OJA, OJS, OJE, & ! [IN]
                             OCEAN_SFC_TEMP(:,:), & ! [IN]
                             LHV(:,:)             ) ! [OUT]

       call CPL_PHY_SFC_fixed_temp( OIA, OIS, OIE,              & ! [IN]
                                    OJA, OJS, OJE,              & ! [IN]
                                    ATMOS_TEMP       (:,:),     & ! [IN]
                                    ATMOS_PRES       (:,:),     & ! [IN]
                                    ATMOS_W          (:,:),     & ! [IN]
                                    ATMOS_U          (:,:),     & ! [IN]
                                    ATMOS_V          (:,:),     & ! [IN]
                                    ATMOS_DENS       (:,:),     & ! [IN]
                                    ATMOS_QV         (:,:),     & ! [IN]
                                    LHV              (:,:),     & ! [IN]
                                    REAL_Z1          (:,:),     & ! [IN]
                                    ATMOS_PBL        (:,:),     & ! [IN]
                                    ATMOS_SFC_DENS   (:,:),     & ! [IN]
                                    ATMOS_SFC_PRES   (:,:),     & ! [IN]
                                    ATMOS_SFLX_rad_dn(:,:,:,:), & ! [IN]
                                    OCEAN_SFC_TEMP   (:,:),     & ! [IN]
                                    WSTR             (:,:),     & ! [IN]
                                    QVEF             (:,:),     & ! [IN]
                                    sfc_albedo       (:,:,:,:), & ! [IN]
                                    SR               (:,:),     & ! [IN]
                                    OCEAN_OCN_Z0M    (:,:),     & ! [IN]
                                    OCEAN_SFC_Z0H    (:,:),     & ! [IN]
                                    OCEAN_SFC_Z0E    (:,:),     & ! [IN]
                                    exists_ocean     (:,:),     & ! [IN]
                                    dt,                         & ! [IN]
                                    OCEAN_SFLX_MW    (:,:),     & ! [OUT]
                                    OCEAN_SFLX_MU    (:,:),     & ! [OUT]
                                    OCEAN_SFLX_MV    (:,:),     & ! [OUT]
                                    OCEAN_SFLX_SH    (:,:),     & ! [OUT]
                                    OCEAN_SFLX_LH    (:,:),     & ! [OUT]
                                    OCEAN_SFLX_QV    (:,:),     & ! [OUT]
                                    OCEAN_SFLX_GH    (:,:),     & ! [OUT]
                                    OCEAN_OCN_Ustar  (:,:),     & ! [OUT]
                                    OCEAN_OCN_Tstar  (:,:),     & ! [OUT]
                                    OCEAN_OCN_Qstar  (:,:),     & ! [OUT]
                                    OCEAN_OCN_Wstar  (:,:),     & ! [OUT]
                                    OCEAN_OCN_RLmo   (:,:),     & ! [OUT]
                                    OCEAN_U10        (:,:),     & ! [OUT]
                                    OCEAN_V10        (:,:),     & ! [OUT]
                                    OCEAN_T2         (:,:),     & ! [OUT]
                                    OCEAN_Q2         (:,:)      ) ! [OUT]
    end select

    !$omp parallel do
    do j = OJS, OJE
    do i = OIS, OIE
    if ( exists_ocean(i,j) ) then
       OCEAN_SFLX_water(i,j) = ATMOS_SFLX_water(i,j) - OCEAN_SFLX_QV(i,j)
       OCEAN_SFLX_ENGI (i,j) = ATMOS_SFLX_ENGI(i,j)                              & ! internal energy of precipitation
                             - OCEAN_SFLX_QV(i,j) * CV_WATER * OCEAN_SFC_TEMP(i,j) ! internal energy of evaporation water
    end if
    enddo
    enddo

    ! weighted average

    if ( ICE_flag ) then

       ! history (open ocean)
       call FILE_HISTORY_in( OCEAN_U10(:,:), 'OCEAN_OCN_U10', 'Wind velocity u at 10 m on open ocean surface', 'm/s'   )
       call FILE_HISTORY_in( OCEAN_V10(:,:), 'OCEAN_OCN_V10', 'Wind velocity v at 10 m on open ocean surface', 'm/s'   )
       call FILE_HISTORY_in( OCEAN_T2 (:,:), 'OCEAN_OCN_T2',  'Air temperature at 2m on open ocean surface',   'K'     )
       call FILE_HISTORY_in( OCEAN_Q2 (:,:), 'OCEAN_OCN_Q2',  'Specific humidity at 2m on open ocean surface', 'kg/kg' )

       !$omp parallel do &
       !$omp private(sfc_frac)
       do j = OJS, OJE
       do i = OIS, OIE
       if ( exists_ocean(i,j) ) then
          sfc_frac = 1.0_RP - OCEAN_ICE_FRAC(i,j)

          OCEAN_SFC_TEMP(i,j) = OCEAN_SFC_TEMP(i,j) * sfc_frac
          OCEAN_SFC_Z0M (i,j) = OCEAN_OCN_Z0M (i,j) * sfc_frac
          OCEAN_SFC_Z0H (i,j) = OCEAN_SFC_Z0H (i,j) * sfc_frac
          OCEAN_SFC_Z0E (i,j) = OCEAN_SFC_Z0E (i,j) * sfc_frac
          OCEAN_SFLX_MW (i,j) = OCEAN_SFLX_MW (i,j) * sfc_frac
          OCEAN_SFLX_MU (i,j) = OCEAN_SFLX_MU (i,j) * sfc_frac
          OCEAN_SFLX_MV (i,j) = OCEAN_SFLX_MV (i,j) * sfc_frac
          OCEAN_SFLX_SH (i,j) = OCEAN_SFLX_SH (i,j) * sfc_frac
          OCEAN_SFLX_LH (i,j) = OCEAN_SFLX_LH (i,j) * sfc_frac
          OCEAN_SFLX_QV (i,j) = OCEAN_SFLX_QV (i,j) * sfc_frac
          OCEAN_U10     (i,j) = OCEAN_U10     (i,j) * sfc_frac
          OCEAN_V10     (i,j) = OCEAN_V10     (i,j) * sfc_frac
          OCEAN_T2      (i,j) = OCEAN_T2      (i,j) * sfc_frac
          OCEAN_Q2      (i,j) = OCEAN_Q2      (i,j) * sfc_frac

          OCEAN_SFLX_GH   (i,j) = OCEAN_SFLX_GH   (i,j) * sfc_frac
          OCEAN_SFLX_water(i,j) = OCEAN_SFLX_water(i,j) * sfc_frac
          OCEAN_SFLX_ENGI (i,j) = OCEAN_SFLX_ENGI (i,j) * sfc_frac

          do irgn = I_R_IR, I_R_VIS
          do idir = I_R_direct, I_R_diffuse
             OCEAN_SFC_albedo(i,j,idir,irgn) = sfc_albedo(i,j,idir,irgn) * sfc_frac
          enddo
          enddo
       end if
       enddo
       enddo

    end if


    !########## surface process (ice) ##########

    if ( ICE_flag ) then

       !$omp parallel do
       do j = OJS, OJE
       do i = OIS, OIE
          subsfc_temp(i,j) = OCEAN_TEMP(OKS,i,j)
       enddo
       enddo

       ! albedo
       select case ( OCEAN_ALB_TYPE )
       case ( 'NAKAJIMA00' )
          call OCEAN_PHY_ALBEDO_seaice( OIA, OIS, OIE,      & ! [IN]
                                        OJA, OJS, OJE,      & ! [IN]
                                        sfc_albedo(:,:,:,:) ) ! [OUT]
       case ( 'CONST' )
          call OCEAN_PHY_ALBEDO_const ( OIA, OIS, OIE,      & ! [IN]
                                        OJA, OJS, OJE,      & ! [IN]
                                        sfc_albedo(:,:,:,:) ) ! [OUT]
       case ( 'INIT' )
          ! Never update OCEAN_SFC_albedo from initial condition
          do irgn = I_R_IR, I_R_VIS
          do idir = I_R_direct, I_R_diffuse
          do j    = OJS, OJE
          do i    = OIS, OIE
             sfc_albedo(i,j,idir,irgn) = OCEAN_SFC_albedo(i,j,idir,irgn)
          enddo
          enddo
          enddo
          enddo
       end select

       ! roughness length
       select case ( OCEAN_RGN_TYPE )
       case ( 'MILLER92', 'MOON07' )
          call OCEAN_PHY_ROUGHNESS_seaice( OIA, OIS, OIE, & ! [IN]
                                           OJA, OJS, OJE, & ! [IN]
                                           sfc_Z0M(:,:),  & ! [OUT]
                                           sfc_Z0H(:,:),  & ! [OUT]
                                           sfc_Z0E(:,:)   ) ! [OUT]
       case ( 'CONST' )
          call OCEAN_PHY_ROUGHNESS_const ( OIA, OIS, OIE, & ! [IN]
                                           OJA, OJS, OJE, & ! [IN]
                                           sfc_Z0M(:,:),  & ! [OUT]
                                           sfc_Z0H(:,:),  & ! [OUT]
                                           sfc_Z0E(:,:)   ) ! [OUT]
       case ( 'INIT' )
          ! Never update OCEAN_SFC_Z0M/H/E from initial condition
          !$omp parallel do
          do j = OJS, OJE
          do i = OIS, OIE
             sfc_Z0M(i,j) = OCEAN_SFC_Z0M(i,j)
             sfc_Z0H(i,j) = OCEAN_SFC_Z0H(i,j)
             sfc_Z0E(i,j) = OCEAN_SFC_Z0E(i,j)
          enddo
          enddo
       end select

       ! thermal conductivity / depth
       call OCEAN_PHY_TC_seaice( OIA, OIS, OIE,       & ! [IN]
                                 OJA, OJS, OJE,       & ! [IN]
                                 OCEAN_ICE_MASS(:,:), & ! [IN]
                                 OCEAN_ICE_FRAC(:,:), & ! [IN]
                                 exists_ice    (:,:), & ! [IN]
                                 TC_dz         (:,:)  ) ! [OUT]

       ! tendency
       select case ( OCEAN_SFC_TYPE )
       case ( 'FIXED-TEMP' )
          !$omp parallel do private(sw)
          do j = OJS, OJE
          do i = OIS, OIE
          if ( exists_ocean(i,j) ) then
             sw = 0.5_RP + sign(0.5_RP, OCEAN_ICE_FRAC(i,j)-EPS)
             ice_mass(i,j) = OCEAN_ICE_MASS(i,j) * sw / ( OCEAN_ICE_FRAC(i,j) + 1.0_RP - sw )
          end if
          end do
          end do

          call HYDROMETEOR_LHS( OIA, OIS, OIE, OJA, OJS, OJE, & ! [IN]
                                OCEAN_ICE_TEMP(:,:), & ! [IN]
                                LHS(:,:)             ) ! [OUT]

          call CPL_PHY_SFC_fixed_temp( OIA, OIS, OIE,              & ! [IN]
                                       OJA, OJS, OJE,              & ! [IN]
                                       ATMOS_TEMP       (:,:),     & ! [IN]
                                       ATMOS_PRES       (:,:),     & ! [IN]
                                       ATMOS_W          (:,:),     & ! [IN]
                                       ATMOS_U          (:,:),     & ! [IN]
                                       ATMOS_V          (:,:),     & ! [IN]
                                       ATMOS_DENS       (:,:),     & ! [IN]
                                       ATMOS_QV         (:,:),     & ! [IN]
                                       LHS              (:,:),     & ! [IN]
                                       REAL_Z1          (:,:),     & ! [IN]
                                       ATMOS_PBL        (:,:),     & ! [IN]
                                       ATMOS_SFC_DENS   (:,:),     & ! [IN]
                                       ATMOS_SFC_PRES   (:,:),     & ! [IN]
                                       ATMOS_SFLX_rad_dn(:,:,:,:), & ! [IN]
                                       OCEAN_ICE_TEMP   (:,:),     & ! [IN]
                                       ice_mass         (:,:),     & ! [IN]
                                       QVEF             (:,:),     & ! [IN]
                                       sfc_albedo       (:,:,:,:), & ! [IN]
                                       SR               (:,:),     & ! [IN]
                                       sfc_Z0M          (:,:),     & ! [IN]
                                       sfc_Z0H          (:,:),     & ! [IN]
                                       sfc_Z0E          (:,:),     & ! [IN]
                                       exists_ice       (:,:),     & ! [IN]
                                       dt,                         & ! [IN]
                                       sflx_MW          (:,:),     & ! [OUT]
                                       sflx_MU          (:,:),     & ! [OUT]
                                       sflx_MV          (:,:),     & ! [OUT]
                                       sflx_SH          (:,:),     & ! [OUT]
                                       sflx_LH          (:,:),     & ! [OUT]
                                       sflx_QV          (:,:),     & ! [OUT]
                                       sflx_GH          (:,:),     & ! [OUT]
                                       OCEAN_ICE_Ustar  (:,:),     & ! [OUT]
                                       OCEAN_ICE_Tstar  (:,:),     & ! [OUT]
                                       OCEAN_ICE_Qstar  (:,:),     & ! [OUT]
                                       OCEAN_ICE_Wstar  (:,:),     & ! [OUT]
                                       OCEAN_ICE_RLmo   (:,:),     & ! [OUT]
                                       U10              (:,:),     & ! [OUT]
                                       V10              (:,:),     & ! [OUT]
                                       T2               (:,:),     & ! [OUT]
                                       Q2               (:,:)      ) ! [OUT]
       end select

       !$omp parallel do
       do j = OJS, OJE
       do i = OIS, OIE
       if ( exists_ocean(i,j) ) then
          sflx_water(i,j) = ATMOS_SFLX_water(i,j) - sflx_QV(i,j)
          sflx_engi (i,j) = ATMOS_SFLX_ENGI(i,j)                           & ! internal energy of precipitation
                          - sflx_QV (i,j) * ( CV_ICE * OCEAN_ICE_TEMP(i,j) - LHF ) ! internal energy of evaporation water
       end if
       enddo
       enddo

       ! weighted average
       !$omp parallel do
       do j = OJS, OJE
       do i = OIS, OIE
       if ( exists_ocean(i,j) ) then
          OCEAN_SFC_TEMP(i,j) = OCEAN_SFC_TEMP(i,j) + OCEAN_ICE_TEMP(i,j) * OCEAN_ICE_FRAC(i,j)
          OCEAN_SFC_Z0M (i,j) = OCEAN_SFC_Z0M (i,j) + sfc_Z0M       (i,j) * OCEAN_ICE_FRAC(i,j)
          OCEAN_SFC_Z0H (i,j) = OCEAN_SFC_Z0H (i,j) + sfc_Z0H       (i,j) * OCEAN_ICE_FRAC(i,j)
          OCEAN_SFC_Z0E (i,j) = OCEAN_SFC_Z0E (i,j) + sfc_Z0E       (i,j) * OCEAN_ICE_FRAC(i,j)
          OCEAN_SFLX_MW (i,j) = OCEAN_SFLX_MW (i,j) + sflx_MW       (i,j) * OCEAN_ICE_FRAC(i,j)
          OCEAN_SFLX_MU (i,j) = OCEAN_SFLX_MU (i,j) + sflx_MU       (i,j) * OCEAN_ICE_FRAC(i,j)
          OCEAN_SFLX_MV (i,j) = OCEAN_SFLX_MV (i,j) + sflx_MV       (i,j) * OCEAN_ICE_FRAC(i,j)
          OCEAN_SFLX_SH (i,j) = OCEAN_SFLX_SH (i,j) + sflx_SH       (i,j) * OCEAN_ICE_FRAC(i,j)
          OCEAN_SFLX_QV (i,j) = OCEAN_SFLX_QV (i,j) + sflx_QV       (i,j) * OCEAN_ICE_FRAC(i,j)
          OCEAN_SFLX_LH (i,j) = OCEAN_SFLX_LH (i,j) + sflx_LH       (i,j) * OCEAN_ICE_FRAC(i,j)
          OCEAN_U10     (i,j) = OCEAN_U10     (i,j) + U10           (i,j) * OCEAN_ICE_FRAC(i,j)
          OCEAN_V10     (i,j) = OCEAN_V10     (i,j) + V10           (i,j) * OCEAN_ICE_FRAC(i,j)
          OCEAN_T2      (i,j) = OCEAN_T2      (i,j) + T2            (i,j) * OCEAN_ICE_FRAC(i,j)
          OCEAN_Q2      (i,j) = OCEAN_Q2      (i,j) + Q2            (i,j) * OCEAN_ICE_FRAC(i,j)

          OCEAN_OFLX_GH   (i,j) = OCEAN_SFLX_GH   (i,j)
          OCEAN_OFLX_water(i,j) = OCEAN_SFLX_water(i,j)
          OCEAN_OFLX_ENGI (i,j) = OCEAN_SFLX_ENGI (i,j)

          OCEAN_SFLX_GH   (i,j) = OCEAN_SFLX_GH   (i,j) + sflx_GH   (i,j) * OCEAN_ICE_FRAC(i,j)
          OCEAN_SFLX_water(i,j) = OCEAN_SFLX_water(i,j) + sflx_water(i,j) * OCEAN_ICE_FRAC(i,j)
          OCEAN_SFLX_ENGI (i,j) = OCEAN_SFLX_ENGI (i,j) + sflx_engi (i,j) * OCEAN_ICE_FRAC(i,j)
       end if
       end do
       end do

       !$omp parallel do
       do irgn = I_R_IR, I_R_VIS
       do idir = I_R_direct, I_R_diffuse
       do j    = OJS, OJE
       do i    = OIS, OIE
       if ( exists_ocean(i,j) ) then
          OCEAN_SFC_albedo(i,j,idir,irgn) = OCEAN_SFC_albedo(i,j,idir,irgn) + sfc_albedo(i,j,idir,irgn) * OCEAN_ICE_FRAC(i,j)
       end if
       enddo
       enddo
       enddo
       enddo


       ! seaice
       select case ( OCEAN_ICE_TYPE )
       case ( 'SIMPLE' )

          !$omp parallel do
          do j = OJS, OJE
          do i = OIS, OIE
          if ( exists_ocean(i,j) ) then
             sflx_hbalance(i,j) = sflx_GH(i,j) + sflx_engi(i,j)
          end if
          enddo
          enddo

          call OCEAN_PHY_ICE_simple( OIA, OIS, OIE,         & ! [IN]
                                     OJA, OJS, OJE,         & ! [IN]
                                     sflx_water      (:,:), & ! [IN]
                                     sflx_hbalance   (:,:), & ! [IN]
                                     subsfc_temp     (:,:), & ! [IN]
                                     TC_dz           (:,:), & ! [IN]
                                     OCEAN_ICE_TEMP  (:,:), & ! [IN]
                                     OCEAN_ICE_MASS  (:,:), & ! [IN]
                                     OCEAN_ICE_FRAC  (:,:), & ! [IN]
                                     exists_ice      (:,:), & ! [IN]
                                     dt,                    & ! [IN]
                                     OCEAN_ICE_TEMP_t(:,:), & ! [OUT]
                                     OCEAN_ICE_MASS_t(:,:), & ! [OUT]
                                     sflx_GH         (:,:), & ! [OUT]
                                     sflx_water      (:,:), & ! [OUT]
                                     sflx_engi       (:,:)  ) ! [OUT]
       case ( 'INIT' )
          !$omp parallel do
          do j = OJS, OJE
          do i = OIS, OIE
          if ( exists_ocean(i,j) ) then
             sflx_GH   (i,j) = sflx_GH(i,j) * OCEAN_ICE_FRAC(i,j)
             sflx_water(i,j) = 0.0_RP ! no flux from seaice to ocean
             sflx_engi (i,j) = 0.0_RP ! no flux from seaice to ocean
          end if
          enddo
          enddo
       end select

       ! history (sea ice)
       call FILE_HISTORY_in( U10(:,:), 'OCEAN_ICE_U10', 'Wind velocity u at 10 m on sea ice surface', 'm/s'   )
       call FILE_HISTORY_in( V10(:,:), 'OCEAN_ICE_V10', 'Wind velocity v at 10 m on sea ice surface', 'm/s'   )
       call FILE_HISTORY_in( T2 (:,:), 'OCEAN_ICE_T2',  'Air temperature at 2m on sea ice surface',   'K'     )
       call FILE_HISTORY_in( Q2 (:,:), 'OCEAN_ICE_Q2',  'Specific humidity at 2m on sea ice surface', 'kg/kg' )

       ! weighted average
       !$omp parallel do
       do j = OJS, OJE
       do i = OIS, OIE
       if ( exists_ocean(i,j) ) then
          OCEAN_OFLX_GH   (i,j) = OCEAN_OFLX_GH   (i,j) + sflx_GH   (i,j)
          OCEAN_OFLX_water(i,j) = OCEAN_OFLX_water(i,j) + sflx_water(i,j)
          OCEAN_OFLX_ENGI (i,j) = OCEAN_OFLX_ENGI (i,j) + sflx_engi (i,j)
       end if
       enddo
       enddo

       call BULKFLUX_diagnose_scales( OIA, OIS, OIE, OJA, OJS, OJE, &
                                      OCEAN_SFLX_MW(:,:), OCEAN_SFLX_MU(:,:), OCEAN_SFLX_MV(:,:), & ! [IN]
                                      OCEAN_SFLX_SH(:,:), OCEAN_SFLX_QV(:,:),                     & ! [IN]
                                      ATMOS_SFC_DENS(:,:), OCEAN_SFC_TEMP(:,:), ATMOS_PBL(:,:),   & ! [IN]
                                      OCEAN_Ustar(:,:), OCEAN_Tstar(:,:), OCEAN_Qstar(:,:),       & ! [OUT]
                                      OCEAN_Wstar(:,:), OCEAN_RLmo(:,:),                          & ! [OUT]
                                      mask = exists_ocean(:,:)                                    ) ! [IN]




    endif ! ICE process?

    if ( .NOT. ATMOS_HYDROMETEOR_dry ) then
       !$omp parallel do
       do j = OJS, OJE
       do i = OIS, OIE
       if ( exists_ocean(i,j) ) then
          OCEAN_SFLX_QTRC(i,j,I_QV) = OCEAN_SFLX_QV(i,j)
       end if
       enddo
       enddo
    endif


    ! Surface flux for chemical tracers
    if ( ATMOS_sw_phy_ch ) then
       call ATMOS_PHY_CH_driver_OCEAN_flux( OCEAN_SFLX_QTRC(:,:,:) ) ! [INOUT]
    endif

    if ( STATISTICS_checktotal ) then
       call STATISTICS_total( OKA, OKS, OKE, OIA, OIS, OIE, OJA, OJS, OJE, &
                              OCEAN_TEMP_t(:,:,:), 'OCEAN_TEMP_t',         &
                              OCEAN_GRID_CARTESC_REAL_VOL(:,:,:),          &
                              OCEAN_GRID_CARTESC_REAL_TOTVOL               )
       if ( ICE_flag ) then
          call STATISTICS_total( OIA, OIS, OIE, OJA, OJS, OJE,                &
                                 OCEAN_ICE_TEMP_t(:,:), 'OCEAN_ICE_TEMP_t',   &
                                 OCEAN_GRID_CARTESC_REAL_AREA(:,:),           &
                                 OCEAN_GRID_CARTESC_REAL_TOTAREA              )
          call STATISTICS_total( OIA, OIS, OIE, OJA, OJS, OJE,                &
                                 OCEAN_ICE_MASS_t(:,:), 'OCEAN_ICE_MASS_t',   &
                                 OCEAN_GRID_CARTESC_REAL_AREA(:,:),           &
                                 OCEAN_GRID_CARTESC_REAL_TOTAREA              )
       end if
    endif

    call PROF_rapend  ('OCN_CalcTend', 1)

    !########## Set Surface Boundary to coupler ##########
    call OCEAN_SURFACE_SET( countup=.true. )

    return
  end subroutine OCEAN_driver_calc_tendency

  !-----------------------------------------------------------------------------
  !> Ocean step
  subroutine OCEAN_driver_update
    use scale_time, only: &
       NOWDAYSEC => TIME_NOWDAYSEC,  &
       dt        => TIME_DTSEC_OCEAN
    use scale_landuse, only: &
       exists_ocean => LANDUSE_exists_ocean
    use mod_ocean_admin, only: &
       OCEAN_DYN_TYPE, &
       OCEAN_ICE_TYPE
    use mod_ocean_vars, only: &
       OCEAN_TEMP,       &
       OCEAN_SALT,       &
       OCEAN_UVEL,       &
       OCEAN_VVEL,       &
       OCEAN_ICE_TEMP,   &
       OCEAN_ICE_MASS,   &
       OCEAN_ICE_FRAC,   &
       OCEAN_TEMP_t,     &
       OCEAN_SALT_t,     &
       OCEAN_UVEL_t,     &
       OCEAN_VVEL_t,     &
       OCEAN_ICE_TEMP_t, &
       OCEAN_ICE_MASS_t, &
       OCEAN_OFLX_GH,    &
       OCEAN_OFLX_water, &
       OCEAN_OFLX_ENGI,  &
       OCEAN_MASS_SUPL,  &
       OCEAN_ENGI_SUPL,  &
       OCEAN_vars_check
    use scale_ocean_dyn_slab, only: &
       OCEAN_DYN_SLAB
    use scale_ocean_dyn_offline, only: &
       OCEAN_DYN_OFFLINE
    use scale_ocean_phy_ice_simple, only: &
       OCEAN_PHY_ICE_adjustment, &
       OCEAN_PHY_ICE_fraction
    use scale_ocean_grid_cartesC, only: &
       CDZ => OCEAN_GRID_CARTESC_CDZ
    implicit none

    real(RP) :: MASS_FLUX(OIA,OJA)
    real(RP) :: ENGI_FLUX(OIA,OJA)
    real(RP) :: MASS_SUPL(OIA,OJA)
    real(RP) :: ENGI_SUPL(OIA,OJA)

    real(RP) :: sflx_GH(OIA,OJA)

    integer :: i, j
    !---------------------------------------------------------------------------

    call PROF_rapstart('OCN_Update', 1)

    !########## Get Surface Boundary from coupler ##########
    call OCEAN_SURFACE_GET

    !########## Dynamics / Update ##########
    !$omp parallel do
    do j = OJS, OJE
    do i = OIS, OIE
       OCEAN_MASS_SUPL(i,j) = 0.0_RP
       OCEAN_ENGI_SUPL(i,j) = 0.0_RP
    end do
    end do

    select case ( OCEAN_DYN_TYPE )
    case ( 'SLAB' )

       !$omp parallel do
       do j = OJS, OJE
       do i = OIS, OIE
       if ( exists_ocean(i,j) ) then
          sflx_GH(i,j) = OCEAN_OFLX_GH(i,j) + OCEAN_OFLX_ENGI(i,j)
       end if
       end do
       end do

       call OCEAN_DYN_SLAB( OKMAX, OKS, OKE,         & ! [IN]
                            OIA,   OIS, OIE,         & ! [IN]
                            OJA,   OJS, OJE,         & ! [IN]
                            OCEAN_TEMP_t    (:,:,:), & ! [IN]
                            sflx_GH         (:,:),   & ! [IN]
                            OCEAN_OFLX_water(:,:),   & ! [IN]
                            exists_ocean    (:,:),   & ! [IN]
                            dt, NOWDAYSEC,           & ! [IN]
                            OCEAN_TEMP      (:,:,:), & ! [INOUT]
                            MASS_SUPL       (:,:),   & ! [OUT]
                            ENGI_SUPL       (:,:)    ) ! [OUT]

       !$omp parallel do
       do j = OJS, OJE
       do i = OIS, OIE
       if ( exists_ocean(i,j) ) then
          OCEAN_MASS_SUPL(i,j) = OCEAN_MASS_SUPL(i,j) + MASS_SUPL(i,j)
          OCEAN_ENGI_SUPL(i,j) = OCEAN_ENGI_SUPL(i,j) + ENGI_SUPL(i,j)
       end if
       end do
       end do

    case ( 'OFFLINE' )

       call OCEAN_DYN_OFFLINE( OKMAX, OKS, OKE,         & ! [IN]
                               OIA,   OIS, OIE,         & ! [IN]
                               OJA,   OJS, OJE,         & ! [IN]
                               exists_ocean    (:,:),   & ! [IN]
                               dt, NOWDAYSEC,           & ! [IN]
                               OCEAN_TEMP      (:,:,:)  ) ! [INOUT]

    case ( 'INIT' )
       ! Never update OCEAN_TEMP from initial condition
    end select

    !########## Ice / Update ##########
    select case ( OCEAN_ICE_TYPE )
    case ( 'SIMPLE' )

       !$omp parallel do
       do j = OJS, OJE
       do i = OIS, OIE
       if ( exists_ocean(i,j) ) then
          OCEAN_ICE_TEMP(i,j) = OCEAN_ICE_TEMP(i,j) + OCEAN_ICE_TEMP_t(i,j) * dt
          OCEAN_ICE_MASS(i,j) = OCEAN_ICE_MASS(i,j) + OCEAN_ICE_MASS_t(i,j) * dt
       end if
       enddo
       enddo

       ! ice adjustment
       call OCEAN_PHY_ICE_adjustment( OIA, OIS, OIE,           & ! [IN]
                                      OJA, OJS, OJE,           & ! [IN]
                                      exists_ocean  (:,:),     & ! [IN]
                                      CDZ(OKS),                & ! [IN]
                                      OCEAN_TEMP    (OKS,:,:), & ! [INOUT]
                                      OCEAN_ICE_TEMP(:,:),     & ! [INOUT]
                                      OCEAN_ICE_MASS(:,:),     & ! [INOUT]
                                      MASS_FLUX     (:,:),     & ! [OUT]
                                      ENGI_FLUX     (:,:),     & ! [OUT]
                                      MASS_SUPL     (:,:),     & ! [OUT]
                                      ENGI_SUPL     (:,:)      ) ! [OUT]

       !$omp parallel do
       do j = OJS, OJE
       do i = OIS, OIE
       if ( exists_ocean(i,j) ) then
          OCEAN_OFLX_water(i,j) = OCEAN_OFLX_water(i,j) - MASS_FLUX(i,j) / dt
          OCEAN_OFLX_ENGI (i,j) = OCEAN_OFLX_ENGI (i,j) - ENGI_FLUX(i,j) / dt
          OCEAN_MASS_SUPL (i,j) = OCEAN_MASS_SUPL (i,j) + MASS_SUPL(i,j) / dt
          OCEAN_ENGI_SUPL (i,j) = OCEAN_ENGI_SUPL (i,j) + ENGI_SUPL(i,j) / dt
       end if
       end do
       end do


       ! update ice fraction
       call OCEAN_PHY_ICE_fraction( OIA, OIS, OIE,       & ! [IN]
                                    OJA, OJS, OJE,       & ! [IN]
                                    OCEAN_ICE_MASS(:,:), & ! [IN]
                                    OCEAN_ICE_FRAC(:,:)  ) ! [OUT]

    case ( 'INIT' )
       ! Never update OCEAN_ICE_TEMP, OCEAN_ICE_MASS, OCEAN_ICE_FRAC from initial condition
    end select

    call OCEAN_vars_check

    call PROF_rapend  ('OCN_Update', 1)

    return
  end subroutine OCEAN_driver_update

  !-----------------------------------------------------------------------------
  !> Get surface boundary from other model
  subroutine OCEAN_SURFACE_GET
    use mod_ocean_admin, only: &
       OCEAN_do
    use mod_ocean_vars, only: &
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
       CPL_getATM_OCN
    implicit none
    !---------------------------------------------------------------------------

    call PROF_rapstart('OCN_SfcExch', 2)

    if ( OCEAN_do ) then
       call CPL_getATM_OCN( ATMOS_TEMP       (:,:),     & ! [OUT]
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

    call PROF_rapend  ('OCN_SfcExch', 2)

    return
  end subroutine OCEAN_SURFACE_GET

  !-----------------------------------------------------------------------------
  !> Put surface boundary to other model
  subroutine OCEAN_SURFACE_SET( countup )
    use mod_ocean_admin, only: &
       OCEAN_do
    use mod_ocean_vars, only: &
       OCEAN_SFC_TEMP,   &
       OCEAN_SFC_albedo, &
       OCEAN_SFC_Z0M,    &
       OCEAN_SFC_Z0H,    &
       OCEAN_SFC_Z0E,    &
       OCEAN_SFLX_MW,    &
       OCEAN_SFLX_MU,    &
       OCEAN_SFLX_MV,    &
       OCEAN_SFLX_SH,    &
       OCEAN_SFLX_LH,    &
       OCEAN_SFLX_QTRC,  &
       OCEAN_U10,        &
       OCEAN_V10,        &
       OCEAN_T2,         &
       OCEAN_Q2,         &
       OCEAN_SFLX_GH
    use mod_cpl_vars, only: &
       CPL_putOCN
    use scale_landuse, only: &
       exists_ocean => LANDUSE_exists_ocean
    implicit none

    logical, intent(in) :: countup
    !---------------------------------------------------------------------------

    call PROF_rapstart('OCN_SfcExch', 2)

    if ( OCEAN_do ) then
       call CPL_putOCN( OCEAN_SFC_TEMP  (:,:),     & ! [IN]
                        OCEAN_SFC_albedo(:,:,:,:), & ! [IN]
                        OCEAN_SFC_Z0M   (:,:),     & ! [IN]
                        OCEAN_SFC_Z0H   (:,:),     & ! [IN]
                        OCEAN_SFC_Z0E   (:,:),     & ! [IN]
                        OCEAN_SFLX_MW   (:,:),     & ! [IN]
                        OCEAN_SFLX_MU   (:,:),     & ! [IN]
                        OCEAN_SFLX_MV   (:,:),     & ! [IN]
                        OCEAN_SFLX_SH   (:,:),     & ! [IN]
                        OCEAN_SFLX_LH   (:,:),     & ! [IN]
                        OCEAN_SFLX_GH   (:,:),     & ! [IN]
                        OCEAN_SFLX_QTRC (:,:,:),   & ! [IN]
                        OCEAN_U10       (:,:),     & ! [IN]
                        OCEAN_V10       (:,:),     & ! [IN]
                        OCEAN_T2        (:,:),     & ! [IN]
                        OCEAN_Q2        (:,:),     & ! [IN]
                        exists_ocean    (:,:),     & ! [IN]
                        countup                    ) ! [IN]
    endif

    call PROF_rapend  ('OCN_SfcExch', 2)

    return
  end subroutine OCEAN_SURFACE_SET

end module mod_ocean_driver
