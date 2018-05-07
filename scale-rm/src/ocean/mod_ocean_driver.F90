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
  real(RP), allocatable :: QVEF(:,:)
  real(RP), allocatable :: SR(:,:)
  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine OCEAN_driver_setup
    use scale_prc, only: &
       PRC_abort
    use mod_ocean_admin, only: &
       OCEAN_do, &
       OCEAN_DYN_TYPE, &
       OCEAN_SFC_TYPE, &
       OCEAN_ALB_TYPE, &
       OCEAN_RGN_TYPE
    use scale_ocean_dyn_slab, only: &
       OCEAN_DYN_SLAB_setup
    use scale_cpl_phy_sfc_fixed_temp, only: &
       CPL_PHY_SFC_fixed_temp_setup
    use scale_ocean_phy_albedo_nakajima00, only: &
       OCEAN_PHY_ALBEDO_nakajima00_setup
    use scale_ocean_phy_roughness_miller92, only: &
       OCEAN_PHY_ROUGHNESS_miller92_setup
    use scale_ocean_phy_roughness_moon07, only: &
       OCEAN_PHY_ROUGHNESS_moon07_setup
    implicit none
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("OCEAN_driver_setup",*) 'Setup'

    if ( OCEAN_do ) then

       select case ( OCEAN_DYN_TYPE )
       case ( 'SLAB' )
          call OCEAN_DYN_SLAB_setup
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

       select case ( OCEAN_ALB_TYPE )
       case ( 'NAKAJIMA00' )
          call OCEAN_PHY_ALBEDO_nakajima00_setup
       case ( 'INIT' )
          ! do nothing
       case default
          LOG_ERROR("OCEAN_driver_setup",*) 'OCEAN_ALB_TYPE is invalid: ', trim(OCEAN_ALB_TYPE)
          call PRC_abort
       end select

       select case ( OCEAN_RGN_TYPE )
       case ( 'MILLER92' )
          call OCEAN_PHY_ROUGHNESS_miller92_setup
       case ( 'MOON07' )
          call OCEAN_PHY_ROUGHNESS_moon07_setup
       case ( 'INIT' )
          ! do nothing
       case default
          LOG_ERROR("OCEAN_driver_setup",*) 'OCEAN_RGN_TYPE is invalid: ', trim(OCEAN_RGN_TYPE)
          call PRC_abort
       end select

       allocate( QVEF(OIA,OJA) )
       allocate( SR(OIA,OJA) )
       QVEF(:,:) = 1.0_RP
       SR(:,:) = 0.0_RP

    end if

    return
  end subroutine OCEAN_driver_setup

  !-----------------------------------------------------------------------------
  !> Calculate tendency
  subroutine OCEAN_driver_calc_tendency( force )
    use scale_atmos_hydrometeor, only: &
       HYDROMETEOR_LHV => ATMOS_HYDROMETEOR_LHV
    use scale_time, only: &
       dt => TIME_DTSEC_OCEAN
    use scale_statistics, only: &
       STATISTICS_checktotal, &
       STATISTICS_total
    use scale_ocean_grid_cartesC_real, only: &
       OCEAN_GRID_CARTESC_REAL_VOL,    &
       OCEAN_GRID_CARTESC_REAL_TOTVOL, &
       OCEAN_GRID_CARTESC_REAL_AREA,   &
       OCEAN_GRID_CARTESC_REAL_TOTAREA
    use scale_file_history, only: &
       FILE_HISTORY_in
    use scale_ocean_grid_cartesC_real, only: &
       OCEAN_GRID_CARTESC_REAL_AREA, &
       OCEAN_GRID_CARTESC_REAL_TOTAREA
    use scale_atmos_grid_cartesC_real, only: &
       REAL_Z1 => ATMOS_GRID_CARTESC_REAL_Z1
    use scale_ocean_phy_albedo_nakajima00, only: &
       OCEAN_PHY_albedo_nakajima00
    use scale_ocean_phy_roughness_miller92, only: &
       OCEAN_PHY_ROUGHNESS_miller92
    use scale_ocean_phy_roughness_moon07, only: &
       OCEAN_PHY_ROUGHNESS_moon07
    use scale_cpl_phy_sfc_fixed_temp, only: &
       CPL_PHY_SFC_fixed_temp
    use mod_ocean_admin, only: &
       OCEAN_SFC_TYPE, &
       OCEAN_ALB_TYPE, &
       OCEAN_RGN_TYPE
    use mod_ocean_vars, only: &
       OCEAN_TEMP,         &
!!$       OCEAN_SALT,         &
!!$       OCEAN_UVEL,         &
!!$       OCEAN_VVEL,         &
       OCEAN_SFC_TEMP,     &
       OCEAN_SFC_albedo,   &
       OCEAN_SFC_Z0M,      &
       OCEAN_SFC_Z0H,      &
       OCEAN_SFC_Z0E,      &
       OCEAN_TEMP_t,       &
       OCEAN_SALT_t,       &
       OCEAN_UVEL_t,       &
       OCEAN_VVEL_t,       &
       OCEAN_SFLX_MW,      &
       OCEAN_SFLX_MU,      &
       OCEAN_SFLX_MV,      &
       OCEAN_SFLX_SH,      &
       OCEAN_SFLX_LH,      &
       OCEAN_SFLX_evap,    &
       OCEAN_SFLX_WH,      &
       OCEAN_SFLX_water,   &
       OCEAN_SFLX_ice,     &
       OCEAN_U10,          &
       OCEAN_V10,          &
       OCEAN_T2,           &
       OCEAN_Q2,           &
       ATMOS_TEMP,         &
       ATMOS_PRES,         &
       ATMOS_W,            &
       ATMOS_U,            &
       ATMOS_V,            &
       ATMOS_DENS,         &
       ATMOS_QV,           &
       ATMOS_PBL,          &
       ATMOS_cosSZA,       &
       ATMOS_SFC_DENS,     &
       ATMOS_SFC_PRES,     &
       ATMOS_SFLX_rad_dn,  &
       ATMOS_SFLX_rain,    &
       ATMOS_SFLX_snow
    use scale_landuse, only: &
       LANDUSE_fact_ocean
    implicit none
    logical, intent(in) :: force

    real(RP) :: LHV(OIA,OJA) ! latent heat of vaporization [J/kg]
    real(RP) :: ATMOS_Uabs(OIA,OJA)

    integer :: k, i, j
    !---------------------------------------------------------------------------

    call PROF_rapstart('OCN_CalcTend', 1)

    !########## Get Surface Boundary from coupler ##########
    call OCEAN_SURFACE_GET

    !########## reset tendencies ##########
!OCL XFILL
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
       ATMOS_Uabs(i,j) = sqrt( ATMOS_U(i,j)**2 + ATMOS_V(i,j)**2 )
    end do
    end do

    !########## ROUGHNESS ##########
    select case ( OCEAN_RGN_TYPE )
    case ( 'MILLER92' )
       call OCEAN_PHY_ROUGHNESS_miller92( OIA, OIS, OIE, OJA, OJS, OJE, &
                                          ATMOS_Uabs(:,:),    & ! [IN]
                                          OCEAN_SFC_Z0M(:,:), & ! [OUT]
                                          OCEAN_SFC_Z0H(:,:), & ! [OUT]
                                          OCEAN_SFC_Z0E(:,:)  ) ! [OUT]
    case ( 'MOON07' )
       call OCEAN_PHY_ROUGHNESS_moon07( OIA, OIS, OIE, OJA, OJS, OJE, &
                                        ATMOS_Uabs(:,:), REAL_Z1(:,:),         & ! [IN]
                                        OCEAN_SFC_Z0M(:,:),                    & ! [INOUT]
                                        OCEAN_SFC_Z0H(:,:), OCEAN_SFC_Z0E(:,:) ) ! [OUT]
    case ( 'INIT' )
       ! Never update OCEAN_SFC_Z0M/H/E from initial condition
    end select

    !########## ALBEDO ##########
    select case ( OCEAN_ALB_TYPE )
    case ( 'NAKAJIMA00' )
       call OCEAN_PHY_albedo_nakajima00( OIA, OIS, OIE, OJA, OJS, OJE,            & ! [IN]
                                         ATMOS_cosSZA    (:,:),                   & ! [IN]
                                         OCEAN_SFC_albedo(:,:,I_R_direct,I_R_VIS) ) ! [OUT]


       OCEAN_SFC_albedo(:,:,I_R_direct ,I_R_NIR) = OCEAN_SFC_albedo(:,:,I_R_direct ,I_R_VIS)

       OCEAN_SFC_albedo(:,:,I_R_diffuse,I_R_NIR) = 0.06_RP
       OCEAN_SFC_albedo(:,:,I_R_diffuse,I_R_VIS) = 0.06_RP
    case ( 'INIT' )
       ! Never update OCEAN_SFC_albedo from initial condition
    end select


    !########## tendency  ##########
    call HYDROMETEOR_LHV( OIA, OIS, OIE, OJA, OJS, OJE, &
                          ATMOS_TEMP(:,:), LHV(:,:) )

    select case ( OCEAN_SFC_TYPE )
    case ( 'FIXED-TEMP' )
!OCL XFILL
       !$omp parallel do
       do j = OJS, OJE
       do i = OIS, OIE
          OCEAN_SFC_TEMP(i,j) = OCEAN_TEMP(OKS,i,j)
       end do
       end do

       call CPL_PHY_SFC_fixed_temp( OIA, OIS, OIE, OJA, OJS, OJE, &
                                    ATMOS_TEMP(:,:), ATMOS_PRES(:,:),                           & ! [IN]
                                    ATMOS_W(:,:), ATMOS_U(:,:), ATMOS_V(:,:),                   & ! [IN]
                                    ATMOS_DENS(:,:), ATMOS_QV(:,:), LHV(:,:),                   & ! [IN]
                                    REAL_Z1(:,:), ATMOS_PBL(:,:),                               & ! [IN]
                                    ATMOS_SFC_DENS(:,:), ATMOS_SFC_PRES(:,:),                   & ! [IN]
                                    ATMOS_SFLX_rad_dn(:,:,:,:),                                 & ! [IN]
                                    OCEAN_SFC_TEMP(:,:), QVEF(:,:),                             & ! [IN]
                                    OCEAN_SFC_albedo(:,:,:,:),                                  & ! [IN]
                                    SR(:,:),                                                    & ! [IN]
                                    OCEAN_SFC_Z0M(:,:), OCEAN_SFC_Z0H(:,:), OCEAN_SFC_Z0E(:,:), & ! [IN]
                                    LANDUSE_fact_ocean, dt,                                     & ! [IN]
                                    OCEAN_SFLX_MW(:,:), OCEAN_SFLX_MU(:,:), OCEAN_SFLX_MV(:,:), & ! [OUT]
                                    OCEAN_SFLX_SH(:,:), OCEAN_SFLX_evap(:,:), OCEAN_SFLX_WH(:,:), & ! [OUT]
                                    OCEAN_U10(:,:), OCEAN_V10(:,:),                             & ! [OUT]
                                    OCEAN_T2(:,:), OCEAN_Q2(:,:)                                ) ! [OUT]
!OCL XFILL
       !$omp parallel do
       do j = OJS, OJE
       do i = OIS, OIE
          OCEAN_SFLX_WH(i,j) = - OCEAN_SFLX_WH(i,j) ! upward to downward
       end do
       end do
    end select

!OCL XFILL
    !$omp parallel do
    do j = OJS, OJE
    do i = OIS, OIE
       OCEAN_SFLX_LH   (i,j) = OCEAN_SFLX_evap(i,j) * LHV(i,j)
       OCEAN_SFLX_water(i,j) = ATMOS_SFLX_rain(i,j) - OCEAN_SFLX_evap(i,j)
       OCEAN_SFLX_ice  (i,j) = ATMOS_SFLX_snow(i,j)
    end do
    end do


    if ( STATISTICS_checktotal ) then
       call STATISTICS_total( OKA, OKS, OKE, OIA, OIS, OIE, OJA, OJS, OJE, &
                              OCEAN_TEMP_t      (:,:,:),    'OCEAN_TEMP_t', &
                              OCEAN_GRID_CARTESC_REAL_VOL(:,:,:),           &
                              OCEAN_GRID_CARTESC_REAL_TOTVOL                )
    end if


    !########## Set Surface Boundary to coupler ##########
    call OCEAN_SURFACE_SET( countup=.true. )

    call PROF_rapend  ('OCN_CalcTend', 1)

    return
  end subroutine OCEAN_driver_calc_tendency

  !-----------------------------------------------------------------------------
  !> Ocean step
  subroutine OCEAN_driver_update
    use scale_time, only: &
       dt => TIME_DTSEC_OCEAN
    use mod_ocean_vars, only: &
       OCEAN_TEMP,         &
!!$       OCEAN_SALT,         &
!!$       OCEAN_UVEL,         &
!!$       OCEAN_VVEL,         &
       OCEAN_SFC_TEMP,     &
       OCEAN_SFC_Z0M,      &
       OCEAN_SFC_Z0H,      &
       OCEAN_SFC_Z0E,      &
       OCEAN_SFLX_WH,      &
       OCEAN_SFLX_water,   &
       OCEAN_SFLX_ice,     &
       OCEAN_TEMP_t,       &
       OCEAN_SALT_t,       &
       OCEAN_UVEL_t,       &
       OCEAN_VVEL_t,       &
       OCEAN_vars_total,   &
       OCEAN_vars_history
    use scale_ocean_dyn_slab, only: &
       OCEAN_DYN_slab
    use scale_landuse, only: &
       LANDUSE_fact_ocean
    use scale_time, only: &
       NOWDAYSEC => TIME_NOWDAYSEC
    use scale_file_history, only: &
       FILE_HISTORY_in
    use mod_ocean_admin, only: &
       OCEAN_DYN_TYPE
    use scale_landuse, only: &
       LANDUSE_fact_ocean
    implicit none

    integer :: k, i, j
    !---------------------------------------------------------------------------

    call PROF_rapstart('OCN_Update', 2)

    !########## Get Surface Boundary from coupler ##########
    call OCEAN_SURFACE_GET

    !########## Dynamics / Update ##########
    select case ( OCEAN_DYN_TYPE )
    case ( 'SLAB' )
       call OCEAN_DYN_slab( OKMAX, OKS, OKE, OIA, OIS, OIE, OJA, OJS, OJE, &
                            OCEAN_TEMP_t(:,:,:),                        & ! [IN]
                            OCEAN_SFLX_WH(:,:),                         & ! [IN]
                            OCEAN_SFLX_water(:,:), OCEAN_SFLX_ice(:,:), & ! [IN]
                            LANDUSE_fact_ocean(:,:),                    & ! [IN]
                            dt, NOWDAYSEC,                              & ! [IN]
                            OCEAN_TEMP(:,:,:)                           ) ! [INOUT]
    case ( 'INIT' )
       ! Never update OCEAN_TEMP from initial condition
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

    integer  :: i, j
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
       OCEAN_SFC_TEMP,      &
       OCEAN_SFC_albedo,    &
       OCEAN_SFC_Z0M,       &
       OCEAN_SFC_Z0H,       &
       OCEAN_SFC_Z0E,       &
       OCEAN_SFLX_MW,       &
       OCEAN_SFLX_MU,       &
       OCEAN_SFLX_MV,       &
       OCEAN_SFLX_SH,       &
       OCEAN_SFLX_LH,       &
       OCEAN_SFLX_WH,       &
       OCEAN_SFLX_evap,     &
       OCEAN_U10,           &
       OCEAN_V10,           &
       OCEAN_T2,            &
       OCEAN_Q2
    use mod_cpl_vars, only: &
       CPL_putOCN
    implicit none

    ! arguments
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
                        OCEAN_SFLX_WH   (:,:),     & ! [IN]
                        OCEAN_SFLX_evap (:,:),     & ! [IN]
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
