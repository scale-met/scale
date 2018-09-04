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
  real(RP), private, allocatable :: QVEF(:,:)
  real(RP), private, allocatable :: SR  (:,:)

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine OCEAN_driver_setup
    use scale_prc, only: &
       PRC_abort
    use mod_ocean_admin, only: &
       OCEAN_do,       &
       OCEAN_DYN_TYPE, &
       OCEAN_SFC_TYPE, &
       OCEAN_ICE_TYPE, &
       OCEAN_ALB_TYPE, &
       OCEAN_RGN_TYPE
    use scale_cpl_phy_sfc_fixed_temp, only: &
       CPL_PHY_SFC_fixed_temp_setup
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
    implicit none
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("OCEAN_driver_setup",*) 'Setup'

    if ( OCEAN_do ) then

       select case ( OCEAN_DYN_TYPE )
       case ( 'SLAB' )
          call OCEAN_DYN_SLAB_setup
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
          ! do nothing
       case ( 'SIMPLE' )
          call OCEAN_PHY_ICE_setup
       case ( 'INIT' )
          ! do nothing
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

       allocate( QVEF(OIA,OJA) )
       allocate( SR  (OIA,OJA) )
       QVEF(:,:) = 1.0_RP
       SR  (:,:) = 0.0_RP

    endif

    return
  end subroutine OCEAN_driver_setup

  !-----------------------------------------------------------------------------
  !> Calculate tendency
  subroutine OCEAN_driver_calc_tendency( force )
    use scale_time, only: &
       dt => TIME_DTSEC_OCEAN
    use scale_statistics, only: &
       STATISTICS_checktotal, &
       STATISTICS_total
    use scale_landuse, only: &
       exists_ocean => LANDUSE_exists_ocean
    use scale_file_history, only: &
       FILE_HISTORY_in
    use scale_atmos_grid_cartesC_real, only: &
       REAL_Z1 => ATMOS_GRID_CARTESC_REAL_Z1
    use scale_atmos_hydrometeor, only: &
       HYDROMETEOR_LHV => ATMOS_HYDROMETEOR_LHV, &
       HYDROMETEOR_LHS => ATMOS_HYDROMETEOR_LHS, &
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
       ATMOS_W,           &
       ATMOS_U,           &
       ATMOS_V,           &
       ATMOS_DENS,        &
       ATMOS_QV,          &
       ATMOS_PBL,         &
       ATMOS_cosSZA,      &
       ATMOS_SFC_DENS,    &
       ATMOS_SFC_PRES,    &
       ATMOS_SFLX_rad_dn, &
       ATMOS_SFLX_rain,   &
       ATMOS_SFLX_snow,   &
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
       OCEAN_SFLX_G,      &
       OCEAN_SFLX_water,  &
       OCEAN_SFLX_ice,    &
       OCEAN_ICE_FRAC, &
       OCEAN_vars_total
    implicit none

    logical, intent(in) :: force

    real(RP) :: LHV          (OIA,OJA) ! latent heat of vaporization [J/kg]
    real(RP) :: LHS          (OIA,OJA) ! latent heat of sublimation  [J/kg]
    real(RP) :: ATMOS_Uabs   (OIA,OJA)
    real(RP) :: sfc_frac     (OIA,OJA)
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
    real(RP) :: sflx_QV      (OIA,OJA)
    real(RP) :: U10          (OIA,OJA)
    real(RP) :: V10          (OIA,OJA)
    real(RP) :: T2           (OIA,OJA)
    real(RP) :: Q2           (OIA,OJA)
    real(RP) :: sflx_hbalance(OIA,OJA)
    real(RP) :: sflx_G       (OIA,OJA)
    real(RP) :: sflx_water   (OIA,OJA)
    real(RP) :: sflx_ice     (OIA,OJA)
    logical  :: exists_ice   (OIA,OJA)

    integer  :: k, i, j, iq, idir, irgn
    !---------------------------------------------------------------------------

    call PROF_rapstart('OCN_CalcTend', 1)

    !########## Get Surface Boundary from coupler ##########
    call OCEAN_SURFACE_GET

    call HYDROMETEOR_LHV( OIA, OIS, OIE, OJA, OJS, OJE, & ! [IN]
                          ATMOS_TEMP(:,:),              & ! [IN]
                          LHV       (:,:)               ) ! [OUT]
    call HYDROMETEOR_LHS( OIA, OIS, OIE, OJA, OJS, OJE, & ! [IN]
                          ATMOS_TEMP(:,:),              & ! [IN]
                          LHS       (:,:)               ) ! [OUT]

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

    !$omp parallel do
    do j = OJS, OJE
    do i = OIS, OIE
       OCEAN_ICE_TEMP_t(i,j) = 0.0_RP
       OCEAN_ICE_MASS_t(i,j) = 0.0_RP
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



    !########## surface process (ice-free ocean) ##########

    !$omp parallel do
    do j = OJS, OJE
    do i = OIS, OIE
       sfc_frac(i,j) = 1.0_RP - OCEAN_ICE_FRAC(i,j)
       sfc_temp(i,j) = OCEAN_TEMP(OKS,i,j)
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
          sfc_albedo(i,j,idir,irgn) = OCEAN_SFC_albedo(i,j,idir,irgn)
       enddo
       enddo
       enddo
       enddo
    end select

    ! roughness length
    select case ( OCEAN_RGN_TYPE )
    case ( 'MILLER92' )
       call OCEAN_PHY_ROUGHNESS_miller92( OIA, OIS, OIE,   & ! [IN]
                                          OJA, OJS, OJE,   & ! [IN]
                                          ATMOS_Uabs(:,:), & ! [IN]
                                          sfc_Z0M   (:,:), & ! [OUT]
                                          sfc_Z0H   (:,:), & ! [OUT]
                                          sfc_Z0E   (:,:)  ) ! [OUT]
    case ( 'MOON07' )
       call OCEAN_PHY_ROUGHNESS_moon07  ( OIA, OIS, OIE,      & ! [IN]
                                          OJA, OJS, OJE,      & ! [IN]
                                          ATMOS_Uabs   (:,:), & ! [IN]
                                          REAL_Z1      (:,:), & ! [IN]
                                          OCEAN_OCN_Z0M(:,:), & ! [INOUT]
                                          sfc_Z0H      (:,:), & ! [OUT]
                                          sfc_Z0E      (:,:)  ) ! [OUT]

       !$omp parallel do
       do j = OJS, OJE
       do i = OIS, OIE
          sfc_Z0M(i,j) = OCEAN_OCN_Z0M(i,j)
       enddo
       enddo
    case ( 'CONST' )
       call OCEAN_PHY_ROUGHNESS_const   ( OIA, OIS, OIE,   & ! [IN]
                                          OJA, OJS, OJE,   & ! [IN]
                                          sfc_Z0M   (:,:), & ! [OUT]
                                          sfc_Z0H   (:,:), & ! [OUT]
                                          sfc_Z0E   (:,:)  ) ! [OUT]
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

    ! tendency
    select case ( OCEAN_SFC_TYPE )
    case ( 'FIXED-TEMP' )
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
                                    sfc_temp         (:,:),     & ! [IN]
                                    QVEF             (:,:),     & ! [IN]
                                    sfc_albedo       (:,:,:,:), & ! [IN]
                                    SR               (:,:),     & ! [IN]
                                    sfc_Z0M          (:,:),     & ! [IN]
                                    sfc_Z0H          (:,:),     & ! [IN]
                                    sfc_Z0E          (:,:),     & ! [IN]
                                    exists_ocean     (:,:),     & ! [IN]
                                    dt,                         & ! [IN]
                                    sflx_MW          (:,:),     & ! [OUT]
                                    sflx_MU          (:,:),     & ! [OUT]
                                    sflx_MV          (:,:),     & ! [OUT]
                                    sflx_SH          (:,:),     & ! [OUT]
                                    sflx_QV          (:,:),     & ! [OUT]
                                    sflx_G           (:,:),     & ! [OUT]
                                    U10              (:,:),     & ! [OUT]
                                    V10              (:,:),     & ! [OUT]
                                    T2               (:,:),     & ! [OUT]
                                    Q2               (:,:)      ) ! [OUT]
    end select

    !$omp parallel do
    do j = OJS, OJE
    do i = OIS, OIE
       sflx_water(i,j) = ATMOS_SFLX_rain(i,j) - sflx_QV(i,j)
       sflx_ice  (i,j) = ATMOS_SFLX_snow(i,j)
    enddo
    enddo

    ! weighted average

    !$omp parallel do
    do j = OJS, OJE
    do i = OIS, OIE
       OCEAN_SFC_TEMP  (i,j)      = sfc_temp  (i,j) * sfc_frac(i,j)
       OCEAN_SFC_Z0M   (i,j)      = sfc_Z0M   (i,j) * sfc_frac(i,j)
       OCEAN_SFC_Z0H   (i,j)      = sfc_Z0H   (i,j) * sfc_frac(i,j)
       OCEAN_SFC_Z0E   (i,j)      = sfc_Z0E   (i,j) * sfc_frac(i,j)
       OCEAN_SFLX_MW   (i,j)      = sflx_MW   (i,j) * sfc_frac(i,j)
       OCEAN_SFLX_MU   (i,j)      = sflx_MU   (i,j) * sfc_frac(i,j)
       OCEAN_SFLX_MV   (i,j)      = sflx_MV   (i,j) * sfc_frac(i,j)
       OCEAN_SFLX_SH   (i,j)      = sflx_SH   (i,j) * sfc_frac(i,j)
       OCEAN_SFLX_QTRC (i,j,I_QV) = sflx_QV   (i,j) * sfc_frac(i,j)
       OCEAN_U10       (i,j)      = U10       (i,j) * sfc_frac(i,j)
       OCEAN_V10       (i,j)      = V10       (i,j) * sfc_frac(i,j)
       OCEAN_T2        (i,j)      = T2        (i,j) * sfc_frac(i,j)
       OCEAN_Q2        (i,j)      = Q2        (i,j) * sfc_frac(i,j)

       OCEAN_SFLX_G    (i,j)      = sflx_G    (i,j) * sfc_frac(i,j) * (-1.0_RP) ! upward to downward
       OCEAN_SFLX_water(i,j)      = sflx_water(i,j) * sfc_frac(i,j)
       OCEAN_SFLX_ice  (i,j)      = sflx_ice  (i,j) * sfc_frac(i,j)
    enddo
    enddo

    !$omp parallel do
    do irgn = I_R_IR, I_R_VIS
    do idir = I_R_direct, I_R_diffuse
    do j    = OJS, OJE
    do i    = OIS, OIE
       OCEAN_SFC_albedo(i,j,idir,irgn) = sfc_albedo(i,j,idir,irgn) * sfc_frac(i,j)
    enddo
    enddo
    enddo
    enddo



    !########## surface process (ice) ##########

    if ( OCEAN_ICE_TYPE /= 'NONE' ) then

       !$omp parallel do
       do j = OJS, OJE
       do i = OIS, OIE
          sfc_frac   (i,j) = OCEAN_ICE_FRAC(i,j)
          sfc_temp   (i,j) = OCEAN_ICE_TEMP(i,j)
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
                                 TC_dz         (:,:)  ) ! [OUT]

       ! tendency
       select case ( OCEAN_SFC_TYPE )
       case ( 'FIXED-TEMP' )
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
                                       sfc_temp         (:,:),     & ! [IN]
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
                                       sflx_QV          (:,:),     & ! [OUT]
                                       sflx_G           (:,:),     & ! [OUT]
                                       U10              (:,:),     & ! [OUT]
                                       V10              (:,:),     & ! [OUT]
                                       T2               (:,:),     & ! [OUT]
                                       Q2               (:,:)      ) ! [OUT]
       end select

       ! seaice
       select case ( OCEAN_ICE_TYPE )
       case ( 'SIMPLE' )

          !$omp parallel do
          do j = OJS, OJE
          do i = OIS, OIE
             sflx_hbalance(i,j) = - sflx_G(i,j) ! upward to downward
          enddo
          enddo

          call OCEAN_PHY_ICE_simple( OIA, OIS, OIE,         & ! [IN]
                                     OJA, OJS, OJE,         & ! [IN]
                                     sflx_QV         (:,:), & ! [IN]
                                     ATMOS_SFLX_rain (:,:), & ! [IN]
                                     ATMOS_SFLX_snow (:,:), & ! [IN]
                                     sflx_hbalance   (:,:), & ! [IN]
                                     subsfc_temp     (:,:), & ! [IN]
                                     TC_dz           (:,:), & ! [IN]
                                     OCEAN_ICE_TEMP  (:,:), & ! [IN]
                                     OCEAN_ICE_MASS  (:,:), & ! [IN]
                                     exists_ice      (:,:), & ! [IN]
                                     dt,                    & ! [IN]
                                     OCEAN_ICE_TEMP_t(:,:), & ! [OUT]
                                     OCEAN_ICE_MASS_t(:,:), & ! [OUT]
                                     sflx_G          (:,:), & ! [OUT]
                                     sflx_water      (:,:), & ! [OUT]
                                     sflx_ice        (:,:)  ) ! [OUT]
       case ( 'INIT' )
          !$omp parallel do
          do j = OJS, OJE
          do i = OIS, OIE
             sflx_G    (i,j) = - sflx_G(i,j) ! upward to downward
             sflx_water(i,j) = 0.0_RP ! no flux from seaice to ocean
             sflx_ice  (i,j) = 0.0_RP ! no flux from seaice to ocean
          enddo
          enddo
       end select

       ! weighted average

       !$omp parallel do
       do j = OJS, OJE
       do i = OIS, OIE
          OCEAN_SFC_TEMP  (i,j)      = OCEAN_SFC_TEMP  (i,j)      + sfc_temp  (i,j) * sfc_frac(i,j)
          OCEAN_SFC_Z0M   (i,j)      = OCEAN_SFC_Z0M   (i,j)      + sfc_Z0M   (i,j) * sfc_frac(i,j)
          OCEAN_SFC_Z0H   (i,j)      = OCEAN_SFC_Z0H   (i,j)      + sfc_Z0H   (i,j) * sfc_frac(i,j)
          OCEAN_SFC_Z0E   (i,j)      = OCEAN_SFC_Z0E   (i,j)      + sfc_Z0E   (i,j) * sfc_frac(i,j)
          OCEAN_SFLX_MW   (i,j)      = OCEAN_SFLX_MW   (i,j)      + sflx_MW   (i,j) * sfc_frac(i,j)
          OCEAN_SFLX_MU   (i,j)      = OCEAN_SFLX_MU   (i,j)      + sflx_MU   (i,j) * sfc_frac(i,j)
          OCEAN_SFLX_MV   (i,j)      = OCEAN_SFLX_MV   (i,j)      + sflx_MV   (i,j) * sfc_frac(i,j)
          OCEAN_SFLX_SH   (i,j)      = OCEAN_SFLX_SH   (i,j)      + sflx_SH   (i,j) * sfc_frac(i,j)
          OCEAN_SFLX_QTRC (i,j,I_QV) = OCEAN_SFLX_QTRC (i,j,I_QV) + sflx_QV   (i,j) * sfc_frac(i,j)
          OCEAN_U10       (i,j)      = OCEAN_U10       (i,j)      + U10       (i,j) * sfc_frac(i,j)
          OCEAN_V10       (i,j)      = OCEAN_V10       (i,j)      + V10       (i,j) * sfc_frac(i,j)
          OCEAN_T2        (i,j)      = OCEAN_T2        (i,j)      + T2        (i,j) * sfc_frac(i,j)
          OCEAN_Q2        (i,j)      = OCEAN_Q2        (i,j)      + Q2        (i,j) * sfc_frac(i,j)

          OCEAN_SFLX_G    (i,j)      = OCEAN_SFLX_G    (i,j)      + sflx_G    (i,j) * sfc_frac(i,j)
          OCEAN_SFLX_water(i,j)      = OCEAN_SFLX_water(i,j)      + sflx_water(i,j) * sfc_frac(i,j)
          OCEAN_SFLX_ice  (i,j)      = OCEAN_SFLX_ice  (i,j)      + sflx_ice  (i,j) * sfc_frac(i,j)
       enddo
       enddo

       !$omp parallel do
       do irgn = I_R_IR, I_R_VIS
       do idir = I_R_direct, I_R_diffuse
       do j    = OJS, OJE
       do i    = OIS, OIE
          OCEAN_SFC_albedo(i,j,idir,irgn) = OCEAN_SFC_albedo(i,j,idir,irgn) + sfc_albedo(i,j,idir,irgn) * sfc_frac(i,j)
       enddo
       enddo
       enddo
       enddo

    endif ! ICE process?

    !$omp parallel do
    do j = OJS, OJE
    do i = OIS, OIE
       OCEAN_SFLX_LH(i,j) = OCEAN_SFLX_QTRC(i,j,I_QV) * LHV(i,j) ! always LHV
    enddo
    enddo

    ! Surface flux for chemical tracers
    if ( ATMOS_sw_phy_ch ) then
       call ATMOS_PHY_CH_driver_OCEAN_flux( OCEAN_SFLX_QTRC(:,:,:) ) ! [INOUT]
    endif

    if ( STATISTICS_checktotal ) then
       call STATISTICS_total( OKA, OKS, OKE, OIA, OIS, OIE, OJA, OJS, OJE, &
                              OCEAN_TEMP_t(:,:,:), 'OCEAN_TEMP_t',         &
                              OCEAN_GRID_CARTESC_REAL_VOL(:,:,:),          &
                              OCEAN_GRID_CARTESC_REAL_TOTVOL               )
       call STATISTICS_total( OIA, OIS, OIE, OJA, OJS, OJE,                &
                              OCEAN_ICE_TEMP_t(:,:), 'OCEAN_ICE_TEMP_t',   &
                              OCEAN_GRID_CARTESC_REAL_AREA(:,:),           &
                              OCEAN_GRID_CARTESC_REAL_TOTAREA              )
       call STATISTICS_total( OIA, OIS, OIE, OJA, OJS, OJE,                &
                              OCEAN_ICE_MASS_t(:,:), 'OCEAN_ICE_MASS_t',   &
                              OCEAN_GRID_CARTESC_REAL_AREA(:,:),           &
                              OCEAN_GRID_CARTESC_REAL_TOTAREA              )
    endif

    !########## Set Surface Boundary to coupler ##########
    call OCEAN_SURFACE_SET( countup=.true. )

    call PROF_rapend  ('OCN_CalcTend', 1)

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
       OCEAN_TEMP,         &
       OCEAN_SALT,         &
       OCEAN_UVEL,         &
       OCEAN_VVEL,         &
       OCEAN_ICE_TEMP,     &
       OCEAN_ICE_MASS,     &
       OCEAN_ICE_FRAC,     &
       OCEAN_TEMP_t,       &
       OCEAN_SALT_t,       &
       OCEAN_UVEL_t,       &
       OCEAN_VVEL_t,       &
       OCEAN_ICE_TEMP_t,   &
       OCEAN_ICE_MASS_t,   &
       OCEAN_SFLX_G,       &
       OCEAN_SFLX_water,   &
       OCEAN_SFLX_ice,     &
       OCEAN_vars_total,   &
       OCEAN_vars_history
    use scale_ocean_dyn_slab, only: &
       OCEAN_DYN_SLAB, &
       OCEAN_DYN_SLAB_DEPTH
    use scale_ocean_dyn_offline, only: &
       OCEAN_DYN_OFFLINE
    use scale_ocean_phy_ice_simple, only: &
       OCEAN_PHY_ICE_adjustment, &
       OCEAN_PHY_ICE_fraction
    implicit none

    integer :: i, j
    !---------------------------------------------------------------------------

    call PROF_rapstart('OCN_Update', 2)

    !########## Get Surface Boundary from coupler ##########
    call OCEAN_SURFACE_GET

    !########## Dynamics / Update ##########
    select case ( OCEAN_DYN_TYPE )
    case ( 'SLAB' )

       call OCEAN_DYN_SLAB( OKMAX, OKS, OKE,         & ! [IN]
                            OIA,   OIS, OIE,         & ! [IN]
                            OJA,   OJS, OJE,         & ! [IN]
                            OCEAN_TEMP_t    (:,:,:), & ! [IN]
                            OCEAN_SFLX_G    (:,:),   & ! [IN]
                            OCEAN_SFLX_water(:,:),   & ! [IN]
                            OCEAN_SFLX_ice  (:,:),   & ! [IN]
                            exists_ocean    (:,:),   & ! [IN]
                            dt, NOWDAYSEC,           & ! [IN]
                            OCEAN_TEMP      (:,:,:)  ) ! [INOUT]

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
          OCEAN_ICE_TEMP(i,j) = OCEAN_ICE_TEMP(i,j) + OCEAN_ICE_TEMP_t(i,j) * dt
          OCEAN_ICE_MASS(i,j) = OCEAN_ICE_MASS(i,j) + OCEAN_ICE_MASS_t(i,j) * dt
       enddo
       enddo

       ! ice adjustment
       call OCEAN_PHY_ICE_adjustment( OIA, OIS, OIE,           & ! [IN]
                                      OJA, OJS, OJE,           & ! [IN]
                                      exists_ocean  (:,:),     & ! [IN]
                                      OCEAN_DYN_SLAB_DEPTH,    & ! [IN]
                                      OCEAN_TEMP    (OKS,:,:), & ! [INOUT]
                                      OCEAN_ICE_TEMP(:,:),     & ! [INOUT]
                                      OCEAN_ICE_MASS(:,:)      ) ! [INOUT]

       ! update ice fraction
       call OCEAN_PHY_ICE_fraction( OIA, OIS, OIE,       & ! [IN]
                                    OJA, OJS, OJE,       & ! [IN]
                                    OCEAN_ICE_MASS(:,:), & ! [IN]
                                    OCEAN_ICE_FRAC(:,:)  ) ! [OUT]

    case ( 'INIT' )
       ! Never update OCEAN_ICE_TEMP, OCEAN_ICE_MASS, OCEAN_ICE_FRAC from initial condition
    end select

    call OCEAN_vars_total

    !########## History & Monitor ##########
    call OCEAN_vars_history

    call PROF_rapend  ('OCN_Update', 2)

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
       ATMOS_SFLX_rain,   &
       ATMOS_SFLX_snow
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
                            ATMOS_SFLX_rain  (:,:),     & ! [OUT]
                            ATMOS_SFLX_snow  (:,:)      ) ! [OUT]
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
       OCEAN_SFLX_G
    use mod_cpl_vars, only: &
       CPL_putOCN
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
                        OCEAN_SFLX_G    (:,:),     & ! [IN]
                        OCEAN_SFLX_QTRC (:,:,:),   & ! [IN]
                        OCEAN_U10       (:,:),     & ! [IN]
                        OCEAN_V10       (:,:),     & ! [IN]
                        OCEAN_T2        (:,:),     & ! [IN]
                        OCEAN_Q2        (:,:),     & ! [IN]
                        countup                    ) ! [IN]
    endif

    call PROF_rapend  ('OCN_SfcExch', 2)

    return
  end subroutine OCEAN_SURFACE_SET

end module mod_ocean_driver
