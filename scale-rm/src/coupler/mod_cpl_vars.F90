!-------------------------------------------------------------------------------
!> module COUPLER Variables
!!
!! @par Description
!!          Container for coupler variables
!!
!! @author Team SCALE
!!
!<
!-------------------------------------------------------------------------------
#include "scalelib.h"
module mod_cpl_vars
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_io
  use scale_prof
  use scale_debug
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
  public :: CPL_vars_setup
  public :: CPL_vars_finalize

  public :: CPL_putATM
  public :: CPL_putOCN
  public :: CPL_putLND
  public :: CPL_putURB
  public :: CPL_getSFC_ATM
  public :: CPL_getATM_OCN
  public :: CPL_getATM_LND
  public :: CPL_getATM_URB

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !

  ! Input from ocean model
  real(RP), public, allocatable :: OCN_SFC_TEMP  (:,:)     ! ocean surface skin temperature [K]
  real(RP), public, allocatable :: OCN_SFC_albedo(:,:,:,:) ! ocean surface albedo (direct/diffuse,IR/near-IR/VIS) (0-1)
  real(RP), public, allocatable :: OCN_SFC_Z0M   (:,:)     ! ocean surface roughness length for momemtum [m]
  real(RP), public, allocatable :: OCN_SFC_Z0H   (:,:)     ! ocean surface roughness length for heat [m]
  real(RP), public, allocatable :: OCN_SFC_Z0E   (:,:)     ! ocean surface roughness length for vapor [m]
  real(RP), public, allocatable :: OCN_SFLX_MW   (:,:)     ! ocean surface w-momentum flux [kg/m/s2]
  real(RP), public, allocatable :: OCN_SFLX_MU   (:,:)     ! ocean surface u-momentum flux [kg/m/s2]
  real(RP), public, allocatable :: OCN_SFLX_MV   (:,:)     ! ocean surface v-momentum flux [kg/m/s2]
  real(RP), public, allocatable :: OCN_SFLX_SH   (:,:)     ! ocean surface sensible heat flux [J/m/s2]
  real(RP), public, allocatable :: OCN_SFLX_LH   (:,:)     ! ocean surface latent heat flux [J/m2/s]
  real(RP), public, allocatable :: OCN_SFLX_GH   (:,:)     ! ocean surface water heat flux [J/m2/s]
  real(RP), public, allocatable :: OCN_SFLX_QTRC (:,:,:)   ! ocean surface tracer flux [kg/m2/s]
  real(RP), public, allocatable :: OCN_SFLX_ENGI (:,:)     ! ocean surface internal energy flux [J/m2/s]
  real(RP), public, allocatable :: OCN_U10       (:,:)     ! ocean velocity u at 10m [m/s]
  real(RP), public, allocatable :: OCN_V10       (:,:)     ! ocean velocity v at 10m [m/s]
  real(RP), public, allocatable :: OCN_T2        (:,:)     ! ocean temperature at 2m [K]
  real(RP), public, allocatable :: OCN_Q2        (:,:)     ! ocean water vapor at 2m [kg/kg]

  ! Input from land model
  real(RP), public, allocatable :: LND_SFC_TEMP  (:,:)     ! land surface skin temperature [K]
  real(RP), public, allocatable :: LND_SFC_albedo(:,:,:,:) ! land surface albedo (direct/diffuse,IR/near-IR/VIS) (0-1)
  real(RP), public, allocatable :: LND_SFC_Z0M   (:,:)     ! land surface roughness length for momemtum [m]
  real(RP), public, allocatable :: LND_SFC_Z0H   (:,:)     ! land surface roughness length for heat [m]
  real(RP), public, allocatable :: LND_SFC_Z0E   (:,:)     ! land surface roughness length for vapor [m]
  real(RP), public, allocatable :: LND_SFLX_MW   (:,:)     ! land surface w-momentum flux [kg/m/s2]
  real(RP), public, allocatable :: LND_SFLX_MU   (:,:)     ! land surface u-momentum flux [kg/m/s2]
  real(RP), public, allocatable :: LND_SFLX_MV   (:,:)     ! land surface v-momentum flux [kg/m/s2]
  real(RP), public, allocatable :: LND_SFLX_SH   (:,:)     ! land surface sensible heat flux [J/m2/s]
  real(RP), public, allocatable :: LND_SFLX_LH   (:,:)     ! land surface latent heat flux [J/m2/s]
  real(RP), public, allocatable :: LND_SFLX_GH   (:,:)     ! land surface ground heat flux [J/m2/s]
  real(RP), public, allocatable :: LND_SFLX_QTRC (:,:,:)   ! land surface tracer flux [kg/m2/s]
  real(RP), public, allocatable :: LND_SFLX_ENGI (:,:)     ! land surface internal energy flux [J/m2/s]
  real(RP), public, allocatable :: LND_U10       (:,:)     ! land velocity u at 10m [m/s]
  real(RP), public, allocatable :: LND_V10       (:,:)     ! land velocity v at 10m [m/s]
  real(RP), public, allocatable :: LND_T2        (:,:)     ! land temperature at 2m [K]
  real(RP), public, allocatable :: LND_Q2        (:,:)     ! land water vapor at 2m [kg/kg]

  ! Input from urban model
  real(RP), public, allocatable :: URB_SFC_TEMP  (:,:)     ! urban surface skin temperature [K]
  real(RP), public, allocatable :: URB_SFC_albedo(:,:,:,:) ! urban surface albedo (direct/diffuse,IR/near-IR/VIS) (0-1)
  real(RP), public, allocatable :: URB_SFC_Z0M   (:,:)     ! urban surface roughness length for momemtum [m]
  real(RP), public, allocatable :: URB_SFC_Z0H   (:,:)     ! urban surface roughness length for heat [m]
  real(RP), public, allocatable :: URB_SFC_Z0E   (:,:)     ! urban surface roughness length for vapor [m]
  real(RP), public, allocatable :: URB_SFLX_MW   (:,:)     ! urban surface w-momentum flux [kg/m/s2]
  real(RP), public, allocatable :: URB_SFLX_MU   (:,:)     ! urban surface u-momentum flux [kg/m/s2]
  real(RP), public, allocatable :: URB_SFLX_MV   (:,:)     ! urban surface v-momentum flux [kg/m/s2]
  real(RP), public, allocatable :: URB_SFLX_SH   (:,:)     ! urban surface sensible heat flux [J/m2/s]
  real(RP), public, allocatable :: URB_SFLX_LH   (:,:)     ! urban surface latent heat flux [J/m2/s]
  real(RP), public, allocatable :: URB_SFLX_SHEX (:,:)     ! urban surface extra sensible heat flux [J/m2/s]
  real(RP), public, allocatable :: URB_SFLX_QVEX (:,:)     ! urban surface extra latent heat flux [kg/kg/m2/s]
  real(RP), public, allocatable :: URB_SFLX_GH   (:,:)     ! urban surface ground heat flux [J/m2/s]
  real(RP), public, allocatable :: URB_SFLX_QTRC (:,:,:)   ! urban surface tracer flux [kg/m2/s]
  real(RP), public, allocatable :: URB_SFLX_ENGI (:,:)     ! urban surface internal energy flux [J/m2/s]
  real(RP), public, allocatable :: URB_U10       (:,:)     ! urban velocity u at 10m [m/s]
  real(RP), public, allocatable :: URB_V10       (:,:)     ! urban velocity v at 10m [m/s]
  real(RP), public, allocatable :: URB_T2        (:,:)     ! urban temperature at 2m [K]
  real(RP), public, allocatable :: URB_Q2        (:,:)     ! urban water vapor at 2m [kg/kg]

  ! Output to ocean model
  real(RP), public, allocatable :: OCN_ATM_TEMP       (:,:)     ! temperature at the lowermost atmosphere layer [K]
  real(RP), public, allocatable :: OCN_ATM_PRES       (:,:)     ! pressure at the lowermost atmosphere layer [Pa]
  real(RP), public, allocatable :: OCN_ATM_W          (:,:)     ! velocity w at the lowermost atmosphere layer [m/s]
  real(RP), public, allocatable :: OCN_ATM_U          (:,:)     ! velocity u at the lowermost atmosphere layer [m/s]
  real(RP), public, allocatable :: OCN_ATM_V          (:,:)     ! velocity v at the lowermost atmosphere layer [m/s]
  real(RP), public, allocatable :: OCN_ATM_DENS       (:,:)     ! density at the lowermost atmosphere layer [kg/m3]
  real(RP), public, allocatable :: OCN_ATM_QV         (:,:)     ! water vapor at the lowermost atmosphere layer [kg/kg]
  real(RP), public, allocatable :: OCN_ATM_PBL        (:,:)     ! the top of atmospheric mixing layer [m]
  real(RP), public, allocatable :: OCN_ATM_SFC_DENS   (:,:)     ! surface density  [kg/m3]
  real(RP), public, allocatable :: OCN_ATM_SFC_PRES   (:,:)     ! surface pressure [Pa]
  real(RP), public, allocatable :: OCN_ATM_SFLX_rad_dn(:,:,:,:) ! downward radiation flux (direct/diffuse,IR/near-IR/VIS) [J/m2/s]
  real(RP), public, allocatable :: OCN_ATM_cosSZA     (:,:)     ! cos(solar zenith angle) (0-1)
  real(RP), public, allocatable :: OCN_ATM_SFLX_water  (:,:)    ! precipitation flux [kg/m2/s]
  real(RP), public, allocatable :: OCN_ATM_SFLX_ENGI  (:,:)     ! internal energy flux [J/m2/s]

  ! Output to land model
  real(RP), public, allocatable :: LND_ATM_TEMP       (:,:)     ! temperature at the lowermost atmosphere layer [K]
  real(RP), public, allocatable :: LND_ATM_PRES       (:,:)     ! pressure at the lowermost atmosphere layer [Pa]
  real(RP), public, allocatable :: LND_ATM_W          (:,:)     ! velocity w at the lowermost atmosphere layer [m/s]
  real(RP), public, allocatable :: LND_ATM_U          (:,:)     ! velocity u at the lowermost atmosphere layer [m/s]
  real(RP), public, allocatable :: LND_ATM_V          (:,:)     ! velocity v at the lowermost atmosphere layer [m/s]
  real(RP), public, allocatable :: LND_ATM_DENS       (:,:)     ! density at the lowermost atmosphere layer [kg/m3]
  real(RP), public, allocatable :: LND_ATM_QV         (:,:)     ! water vapor at the lowermost atmosphere layer [kg/kg]
  real(RP), public, allocatable :: LND_ATM_PBL        (:,:)     ! the top of atmospheric mixing layer [m]
  real(RP), public, allocatable :: LND_ATM_SFC_DENS   (:,:)     ! surface density  [kg/m3]
  real(RP), public, allocatable :: LND_ATM_SFC_PRES   (:,:)     ! surface pressure [Pa]
  real(RP), public, allocatable :: LND_ATM_SFLX_rad_dn(:,:,:,:) ! downward radiation flux (direct/diffuse,IR/near-IR/VIS) [J/m2/s]
  real(RP), public, allocatable :: LND_ATM_cosSZA     (:,:)     ! cos(solar zenith angle) (0-1)
  real(RP), public, allocatable :: LND_ATM_SFLX_water (:,:)     ! precipitation flux [kg/m2/s]
  real(RP), public, allocatable :: LND_ATM_SFLX_ENGI  (:,:)     ! internal energy flux [J/m2/s]

  ! Output to urban model
  real(RP), public, allocatable :: URB_ATM_TEMP       (:,:)     ! temperature at the lowermost atmosphere layer [K]
  real(RP), public, allocatable :: URB_ATM_PRES       (:,:)     ! pressure at the lowermost atmosphere layer [Pa]
  real(RP), public, allocatable :: URB_ATM_W          (:,:)     ! velocity w at the lowermost atmosphere layer [m/s]
  real(RP), public, allocatable :: URB_ATM_U          (:,:)     ! velocity u at the lowermost atmosphere layer [m/s]
  real(RP), public, allocatable :: URB_ATM_V          (:,:)     ! velocity v at the lowermost atmosphere layer [m/s]
  real(RP), public, allocatable :: URB_ATM_DENS       (:,:)     ! density at the lowermost atmosphere layer [kg/m3]
  real(RP), public, allocatable :: URB_ATM_QV         (:,:)     ! water vapor at the lowermost atmosphere layer [kg/kg]
  real(RP), public, allocatable :: URB_ATM_PBL        (:,:)     ! the top of atmospheric mixing layer [m]
  real(RP), public, allocatable :: URB_ATM_SFC_DENS   (:,:)     ! surface density  [kg/m3]
  real(RP), public, allocatable :: URB_ATM_SFC_PRES   (:,:)     ! surface pressure [Pa]
  real(RP), public, allocatable :: URB_ATM_SFLX_rad_dn(:,:,:,:) ! downward radiation flux (direct/diffuse,IR/near-IR/VIS) [J/m2/s]
  real(RP), public, allocatable :: URB_ATM_cosSZA     (:,:)     ! cos(solar zenith angle) (0-1)
  real(RP), public, allocatable :: URB_ATM_SFLX_water (:,:)     ! precipitation flux [kg/m2/s]
  real(RP), public, allocatable :: URB_ATM_SFLX_ENGI  (:,:)     ! internal energy flux [J/m2/s]

  ! counter
  real(RP), public :: CNT_putATM_OCN ! put counter for atmos to ocean
  real(RP), public :: CNT_putATM_LND ! put counter for atmos to land
  real(RP), public :: CNT_putATM_URB ! put counter for atmos to urban
  real(RP), public :: CNT_putOCN     ! put counter for ocean
  real(RP), public :: CNT_putLND     ! put counter for land
  real(RP), public :: CNT_putURB     ! put counter for urban

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
  subroutine CPL_vars_setup
    use scale_const, only: &
       UNDEF => CONST_UNDEF
    use scale_prc, only: &
       PRC_abort
    use scale_landuse, only: &
       LANDUSE_fact_ocean, &
       LANDUSE_fact_land,  &
       LANDUSE_fact_urban, &
       LANDUSE_fact_lake
    use scale_atmos_hydrometeor, only: &
       ATMOS_HYDROMETEOR_dry
    use mod_ocean_admin, only: &
       OCEAN_do
    use mod_land_admin, only: &
       LAND_do
    use mod_urban_admin, only: &
       URBAN_do
    use mod_lake_admin, only: &
       LAKE_do
    implicit none

    real(RP) :: checkfact
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("CPL_vars_setup",*) 'Setup'
    LOG_INFO("CPL_vars_setup",*) 'No namelists.'

    ! Check consistency of OCEAN_do and LANDUSE_fact_ocean
    checkfact = maxval( LANDUSE_fact_ocean(:,:) )
    if ( .NOT. OCEAN_do .AND. checkfact > 0.0_RP ) then
       LOG_ERROR("CPL_vars_setup",*) 'Ocean fraction exists, but ocean component has not been called. Please check this inconsistency. STOP.', checkfact
       call PRC_abort
    endif

    ! Check consistency of LAND_do and LANDUSE_fact_land
    checkfact = maxval( LANDUSE_fact_land(:,:) )
    if ( .NOT. LAND_do .AND. checkfact > 0.0_RP ) then
       LOG_ERROR("CPL_vars_setup",*) 'Land fraction exists, but land component has not been called. Please check this inconsistency. STOP.', checkfact
       call PRC_abort
    endif

    ! Check consistency of URBAN_do and LANDUSE_fact_urban
    checkfact = maxval( LANDUSE_fact_urban(:,:) )
    if ( .NOT. URBAN_do .AND. checkfact > 0.0_RP ) then
       LOG_ERROR("CPL_vars_setup",*) 'Urban fraction exists, but urban component has not been called. Please check this inconsistency. STOP.', checkfact
       call PRC_abort
    endif

    ! Check consistency of LAKE_do and LANDUSE_fact_lake
    checkfact = maxval( LANDUSE_fact_lake(:,:) )
    if ( .NOT. LAKE_do .AND. checkfact > 0.0_RP ) then
       LOG_ERROR("CPL_vars_setup",*) 'Lake fraction exists, but lake component has not been called. Please check this inconsistency. STOP.', checkfact
       call PRC_abort
    endif


    allocate( OCN_SFC_TEMP  (IA,JA)                     )
    allocate( OCN_SFC_albedo(IA,JA,N_RAD_DIR,N_RAD_RGN) )
    allocate( OCN_SFC_Z0M   (IA,JA)                     )
    allocate( OCN_SFC_Z0H   (IA,JA)                     )
    allocate( OCN_SFC_Z0E   (IA,JA)                     )
    allocate( OCN_SFLX_MU   (IA,JA)                     )
    allocate( OCN_SFLX_MV   (IA,JA)                     )
    allocate( OCN_SFLX_MW   (IA,JA)                     )
    allocate( OCN_SFLX_SH   (IA,JA)                     )
    allocate( OCN_SFLX_LH   (IA,JA)                     )
    allocate( OCN_SFLX_GH   (IA,JA)                     )
    allocate( OCN_SFLX_QTRC (IA,JA,QA)                  )
    allocate( OCN_SFLX_ENGI (IA,JA)                     )
    allocate( OCN_U10       (IA,JA)                     )
    allocate( OCN_V10       (IA,JA)                     )
    allocate( OCN_T2        (IA,JA)                     )
    allocate( OCN_Q2        (IA,JA)                     )
    OCN_SFC_TEMP  (:,:)     = UNDEF
    OCN_SFC_albedo(:,:,:,:) = UNDEF
    OCN_SFC_Z0M   (:,:)     = UNDEF
    OCN_SFC_Z0H   (:,:)     = UNDEF
    OCN_SFC_Z0E   (:,:)     = UNDEF
    OCN_SFLX_MU   (:,:)     = UNDEF
    OCN_SFLX_MV   (:,:)     = UNDEF
    OCN_SFLX_MW   (:,:)     = UNDEF
    OCN_SFLX_SH   (:,:)     = UNDEF
    OCN_SFLX_LH   (:,:)     = UNDEF
    OCN_SFLX_GH   (:,:)     = UNDEF
    OCN_SFLX_QTRC (:,:,:)   = UNDEF
    OCN_SFLX_ENGI (:,:)     = UNDEF
    OCN_U10       (:,:)     = UNDEF
    OCN_V10       (:,:)     = UNDEF
    OCN_T2        (:,:)     = UNDEF
    OCN_Q2        (:,:)     = UNDEF

    allocate( LND_SFC_TEMP  (IA,JA)                     )
    allocate( LND_SFC_albedo(IA,JA,N_RAD_DIR,N_RAD_RGN) )
    allocate( LND_SFC_Z0M   (IA,JA)                     )
    allocate( LND_SFC_Z0H   (IA,JA)                     )
    allocate( LND_SFC_Z0E   (IA,JA)                     )
    allocate( LND_SFLX_MU   (IA,JA)                     )
    allocate( LND_SFLX_MV   (IA,JA)                     )
    allocate( LND_SFLX_MW   (IA,JA)                     )
    allocate( LND_SFLX_SH   (IA,JA)                     )
    allocate( LND_SFLX_LH   (IA,JA)                     )
    allocate( LND_SFLX_GH   (IA,JA)                     )
    allocate( LND_SFLX_QTRC (IA,JA,QA)                  )
    allocate( LND_SFLX_ENGI (IA,JA)                     )
    allocate( LND_U10       (IA,JA)                     )
    allocate( LND_V10       (IA,JA)                     )
    allocate( LND_T2        (IA,JA)                     )
    allocate( LND_Q2        (IA,JA)                     )
    LND_SFC_TEMP  (:,:)     = UNDEF
    LND_SFC_albedo(:,:,:,:) = UNDEF
    LND_SFC_Z0M   (:,:)     = UNDEF
    LND_SFC_Z0H   (:,:)     = UNDEF
    LND_SFC_Z0E   (:,:)     = UNDEF
    LND_SFLX_MU   (:,:)     = UNDEF
    LND_SFLX_MV   (:,:)     = UNDEF
    LND_SFLX_MW   (:,:)     = UNDEF
    LND_SFLX_SH   (:,:)     = UNDEF
    LND_SFLX_LH   (:,:)     = UNDEF
    LND_SFLX_GH   (:,:)     = UNDEF
    LND_SFLX_QTRC (:,:,:)   = UNDEF
    LND_SFLX_ENGI (:,:)     = UNDEF
    LND_U10       (:,:)     = UNDEF
    LND_V10       (:,:)     = UNDEF
    LND_T2        (:,:)     = UNDEF
    LND_Q2        (:,:)     = UNDEF

    allocate( URB_SFC_TEMP  (IA,JA)                     )
    allocate( URB_SFC_albedo(IA,JA,N_RAD_DIR,N_RAD_RGN) )
    allocate( URB_SFC_Z0M   (IA,JA)                     )
    allocate( URB_SFC_Z0H   (IA,JA)                     )
    allocate( URB_SFC_Z0E   (IA,JA)                     )
    allocate( URB_SFLX_MU   (IA,JA)                     )
    allocate( URB_SFLX_MV   (IA,JA)                     )
    allocate( URB_SFLX_MW   (IA,JA)                     )
    allocate( URB_SFLX_SH   (IA,JA)                     )
    allocate( URB_SFLX_LH   (IA,JA)                     )
    allocate( URB_SFLX_SHEX (IA,JA)                     )
    allocate( URB_SFLX_QVEX (IA,JA)                     )
    allocate( URB_SFLX_GH   (IA,JA)                     )
    allocate( URB_SFLX_QTRC (IA,JA,QA)                  )
    allocate( URB_SFLX_ENGI (IA,JA)                     )
    allocate( URB_U10       (IA,JA)                     )
    allocate( URB_V10       (IA,JA)                     )
    allocate( URB_T2        (IA,JA)                     )
    allocate( URB_Q2        (IA,JA)                     )
    URB_SFC_TEMP  (:,:)     = UNDEF
    URB_SFC_albedo(:,:,:,:) = UNDEF
    URB_SFC_Z0M   (:,:)     = UNDEF
    URB_SFC_Z0H   (:,:)     = UNDEF
    URB_SFC_Z0E   (:,:)     = UNDEF
    URB_SFLX_MU   (:,:)     = UNDEF
    URB_SFLX_MV   (:,:)     = UNDEF
    URB_SFLX_MW   (:,:)     = UNDEF
    URB_SFLX_SH   (:,:)     = UNDEF
    URB_SFLX_LH   (:,:)     = UNDEF
    URB_SFLX_SHEX (:,:)     = UNDEF
    URB_SFLX_QVEX (:,:)     = UNDEF
    URB_SFLX_GH   (:,:)     = UNDEF
    URB_SFLX_QTRC (:,:,:)   = UNDEF
    URB_SFLX_ENGI (:,:)     = UNDEF
    URB_U10       (:,:)     = UNDEF
    URB_V10       (:,:)     = UNDEF
    URB_T2        (:,:)     = UNDEF
    URB_Q2        (:,:)     = UNDEF

    allocate( OCN_ATM_TEMP       (IA,JA)                     )
    allocate( OCN_ATM_PRES       (IA,JA)                     )
    allocate( OCN_ATM_W          (IA,JA)                     )
    allocate( OCN_ATM_U          (IA,JA)                     )
    allocate( OCN_ATM_V          (IA,JA)                     )
    allocate( OCN_ATM_DENS       (IA,JA)                     )
    allocate( OCN_ATM_QV         (IA,JA)                     )
    allocate( OCN_ATM_PBL        (IA,JA)                     )
    allocate( OCN_ATM_SFC_DENS   (IA,JA)                     )
    allocate( OCN_ATM_SFC_PRES   (IA,JA)                     )
    allocate( OCN_ATM_SFLX_rad_dn(IA,JA,N_RAD_DIR,N_RAD_RGN) )
    allocate( OCN_ATM_cosSZA     (IA,JA)                     )
    allocate( OCN_ATM_SFLX_water (IA,JA)                     )
    allocate( OCN_ATM_SFLX_ENGI  (IA,JA)                     )
    OCN_ATM_TEMP       (:,:)     = UNDEF
    OCN_ATM_PRES       (:,:)     = UNDEF
    OCN_ATM_W          (:,:)     = UNDEF
    OCN_ATM_U          (:,:)     = UNDEF
    OCN_ATM_V          (:,:)     = UNDEF
    OCN_ATM_DENS       (:,:)     = UNDEF
    OCN_ATM_QV         (:,:)     = UNDEF
    OCN_ATM_PBL        (:,:)     = UNDEF
    OCN_ATM_SFC_DENS   (:,:)     = UNDEF
    OCN_ATM_SFC_PRES   (:,:)     = UNDEF
    OCN_ATM_SFLX_rad_dn(:,:,:,:) = UNDEF
    OCN_ATM_cosSZA     (:,:)     = UNDEF
    OCN_ATM_SFLX_water (:,:)     = UNDEF
    OCN_ATM_SFLX_ENGI  (:,:)     = UNDEF

    allocate( LND_ATM_TEMP       (IA,JA)                     )
    allocate( LND_ATM_PRES       (IA,JA)                     )
    allocate( LND_ATM_W          (IA,JA)                     )
    allocate( LND_ATM_U          (IA,JA)                     )
    allocate( LND_ATM_V          (IA,JA)                     )
    allocate( LND_ATM_DENS       (IA,JA)                     )
    allocate( LND_ATM_QV         (IA,JA)                     )
    allocate( LND_ATM_PBL        (IA,JA)                     )
    allocate( LND_ATM_SFC_DENS   (IA,JA)                     )
    allocate( LND_ATM_SFC_PRES   (IA,JA)                     )
    allocate( LND_ATM_SFLX_rad_dn(IA,JA,N_RAD_DIR,N_RAD_RGN) )
    allocate( LND_ATM_cosSZA     (IA,JA)                     )
    allocate( LND_ATM_SFLX_water (IA,JA)                     )
    allocate( LND_ATM_SFLX_ENGI  (IA,JA)                     )
    LND_ATM_TEMP       (:,:)     = UNDEF
    LND_ATM_PRES       (:,:)     = UNDEF
    LND_ATM_W          (:,:)     = UNDEF
    LND_ATM_U          (:,:)     = UNDEF
    LND_ATM_V          (:,:)     = UNDEF
    LND_ATM_DENS       (:,:)     = UNDEF
    LND_ATM_QV         (:,:)     = UNDEF
    LND_ATM_PBL        (:,:)     = UNDEF
    LND_ATM_SFC_DENS   (:,:)     = UNDEF
    LND_ATM_SFC_PRES   (:,:)     = UNDEF
    LND_ATM_SFLX_rad_dn(:,:,:,:) = UNDEF
    LND_ATM_cosSZA     (:,:)     = UNDEF
    LND_ATM_SFLX_water (:,:)     = UNDEF
    LND_ATM_SFLX_ENGI  (:,:)     = UNDEF

    allocate( URB_ATM_TEMP       (IA,JA)                     )
    allocate( URB_ATM_PRES       (IA,JA)                     )
    allocate( URB_ATM_W          (IA,JA)                     )
    allocate( URB_ATM_U          (IA,JA)                     )
    allocate( URB_ATM_V          (IA,JA)                     )
    allocate( URB_ATM_DENS       (IA,JA)                     )
    allocate( URB_ATM_QV         (IA,JA)                     )
    allocate( URB_ATM_PBL        (IA,JA)                     )
    allocate( URB_ATM_SFC_DENS   (IA,JA)                     )
    allocate( URB_ATM_SFC_PRES   (IA,JA)                     )
    allocate( URB_ATM_SFLX_rad_dn(IA,JA,N_RAD_DIR,N_RAD_RGN) )
    allocate( URB_ATM_cosSZA     (IA,JA)                     )
    allocate( URB_ATM_SFLX_water (IA,JA)                     )
    allocate( URB_ATM_SFLX_ENGI  (IA,JA)                     )
    URB_ATM_TEMP       (:,:)     = UNDEF
    URB_ATM_PRES       (:,:)     = UNDEF
    URB_ATM_W          (:,:)     = UNDEF
    URB_ATM_U          (:,:)     = UNDEF
    URB_ATM_V          (:,:)     = UNDEF
    URB_ATM_DENS       (:,:)     = UNDEF
    URB_ATM_QV         (:,:)     = UNDEF
    URB_ATM_PBL        (:,:)     = UNDEF
    URB_ATM_SFC_DENS   (:,:)     = UNDEF
    URB_ATM_SFC_PRES   (:,:)     = UNDEF
    URB_ATM_SFLX_rad_dn(:,:,:,:) = UNDEF
    URB_ATM_cosSZA     (:,:)     = UNDEF
    URB_ATM_SFLX_water (:,:)     = UNDEF
    URB_ATM_SFLX_ENGI  (:,:)     = UNDEF

    ! counter intialize
    CNT_putATM_OCN = 0.0_RP
    CNT_putATM_LND = 0.0_RP
    CNT_putATM_URB = 0.0_RP
    CNT_putOCN     = 0.0_RP
    CNT_putLND     = 0.0_RP
    CNT_putURB     = 0.0_RP

    if ( ATMOS_HYDROMETEOR_dry ) then
       OCN_ATM_QV = 0.0_RP
       LND_ATM_QV = 0.0_RP
       URB_ATM_QV = 0.0_RP
    endif

    return
  end subroutine CPL_vars_setup

  !-----------------------------------------------------------------------------
  !> Finalize
  subroutine CPL_vars_finalize
    implicit none
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("CPL_vars_finalize",*) 'Finalize'

    deallocate( OCN_SFC_TEMP   )
    deallocate( OCN_SFC_albedo )
    deallocate( OCN_SFC_Z0M    )
    deallocate( OCN_SFC_Z0H    )
    deallocate( OCN_SFC_Z0E    )
    deallocate( OCN_SFLX_MU    )
    deallocate( OCN_SFLX_MV    )
    deallocate( OCN_SFLX_MW    )
    deallocate( OCN_SFLX_SH    )
    deallocate( OCN_SFLX_LH    )
    deallocate( OCN_SFLX_GH    )
    deallocate( OCN_SFLX_QTRC  )
    deallocate( OCN_SFLX_ENGI  )
    deallocate( OCN_U10        )
    deallocate( OCN_V10        )
    deallocate( OCN_T2         )
    deallocate( OCN_Q2         )

    deallocate( LND_SFC_TEMP   )
    deallocate( LND_SFC_albedo )
    deallocate( LND_SFC_Z0M    )
    deallocate( LND_SFC_Z0H    )
    deallocate( LND_SFC_Z0E    )
    deallocate( LND_SFLX_MU    )
    deallocate( LND_SFLX_MV    )
    deallocate( LND_SFLX_MW    )
    deallocate( LND_SFLX_SH    )
    deallocate( LND_SFLX_LH    )
    deallocate( LND_SFLX_GH    )
    deallocate( LND_SFLX_QTRC  )
    deallocate( LND_SFLX_ENGI  )
    deallocate( LND_U10        )
    deallocate( LND_V10        )
    deallocate( LND_T2         )
    deallocate( LND_Q2         )

    deallocate( URB_SFC_TEMP   )
    deallocate( URB_SFC_albedo )
    deallocate( URB_SFC_Z0M    )
    deallocate( URB_SFC_Z0H    )
    deallocate( URB_SFC_Z0E    )
    deallocate( URB_SFLX_MU    )
    deallocate( URB_SFLX_MV    )
    deallocate( URB_SFLX_MW    )
    deallocate( URB_SFLX_SH    )
    deallocate( URB_SFLX_LH    )
    deallocate( URB_SFLX_SHEX  )
    deallocate( URB_SFLX_QVEX  )
    deallocate( URB_SFLX_GH    )
    deallocate( URB_SFLX_QTRC  )
    deallocate( URB_SFLX_ENGI  )
    deallocate( URB_U10        )
    deallocate( URB_V10        )
    deallocate( URB_T2         )
    deallocate( URB_Q2         )

    deallocate( OCN_ATM_TEMP        )
    deallocate( OCN_ATM_PRES        )
    deallocate( OCN_ATM_W           )
    deallocate( OCN_ATM_U           )
    deallocate( OCN_ATM_V           )
    deallocate( OCN_ATM_DENS        )
    deallocate( OCN_ATM_QV          )
    deallocate( OCN_ATM_PBL         )
    deallocate( OCN_ATM_SFC_DENS    )
    deallocate( OCN_ATM_SFC_PRES    )
    deallocate( OCN_ATM_SFLX_rad_dn )
    deallocate( OCN_ATM_cosSZA      )
    deallocate( OCN_ATM_SFLX_water  )
    deallocate( OCN_ATM_SFLX_ENGI   )

    deallocate( LND_ATM_TEMP        )
    deallocate( LND_ATM_PRES        )
    deallocate( LND_ATM_W           )
    deallocate( LND_ATM_U           )
    deallocate( LND_ATM_V           )
    deallocate( LND_ATM_DENS        )
    deallocate( LND_ATM_QV          )
    deallocate( LND_ATM_PBL         )
    deallocate( LND_ATM_SFC_DENS    )
    deallocate( LND_ATM_SFC_PRES    )
    deallocate( LND_ATM_SFLX_rad_dn )
    deallocate( LND_ATM_cosSZA      )
    deallocate( LND_ATM_SFLX_water  )
    deallocate( LND_ATM_SFLX_ENGI   )

    deallocate( URB_ATM_TEMP        )
    deallocate( URB_ATM_PRES        )
    deallocate( URB_ATM_W           )
    deallocate( URB_ATM_U           )
    deallocate( URB_ATM_V           )
    deallocate( URB_ATM_DENS        )
    deallocate( URB_ATM_QV          )
    deallocate( URB_ATM_PBL         )
    deallocate( URB_ATM_SFC_DENS    )
    deallocate( URB_ATM_SFC_PRES    )
    deallocate( URB_ATM_SFLX_rad_dn )
    deallocate( URB_ATM_cosSZA      )
    deallocate( URB_ATM_SFLX_water  )
    deallocate( URB_ATM_SFLX_ENGI   )

    return
  end subroutine CPL_vars_finalize

  !-----------------------------------------------------------------------------
  subroutine CPL_putATM( &
       TEMP,        &
       PRES,        &
       W,           &
       U,           &
       V,           &
       DENS,        &
       QV,          &
       PBL,         &
       SFC_DENS,    &
       SFC_PRES,    &
       SFLX_rad_dn, &
       cosSZA,      &
       SFLX_water,  &
       SFLX_ENGI,   &
       countup      )
    use scale_atmos_hydrometeor, only: &
       ATMOS_HYDROMETEOR_dry, &
       CV_WATER, &
       CV_ICE,   &
       LHF
    implicit none

    real(RP), intent(in) :: TEMP       (IA,JA)
    real(RP), intent(in) :: PRES       (IA,JA)
    real(RP), intent(in) :: W          (IA,JA)
    real(RP), intent(in) :: U          (IA,JA)
    real(RP), intent(in) :: V          (IA,JA)
    real(RP), intent(in) :: DENS       (IA,JA)
    real(RP), intent(in) :: QV         (IA,JA)
    real(RP), intent(in) :: PBL        (IA,JA)
    real(RP), intent(in) :: SFC_DENS   (IA,JA)
    real(RP), intent(in) :: SFC_PRES   (IA,JA)
    real(RP), intent(in) :: SFLX_rad_dn(IA,JA,N_RAD_DIR,N_RAD_RGN)
    real(RP), intent(in) :: cosSZA     (IA,JA)
    real(RP), intent(in) :: SFLX_water (IA,JA)
    real(RP), intent(in) :: SFLX_ENGI  (IA,JA)
    logical,  intent(in) :: countup

    integer :: i, j, idir, irgn
    !---------------------------------------------------------------------------

    !$omp parallel do default(none) private(i,j,idir,irgn) OMP_SCHEDULE_ &
    !$omp shared(JS,JE,IS,IE,OCN_ATM_TEMP,OCN_ATM_PRES,OCN_ATM_W,OCN_ATM_U) &
    !$omp shared(OCN_ATM_V,OCN_ATM_DENS,CNT_putATM_OCN,TEMP,PRES,W,U,V,DENS) &
    !$omp shared(ATMOS_HYDROMETEOR_dry,OCN_ATM_QV,OCN_ATM_PBL,OCN_ATM_SFC_DENS,OCN_ATM_SFC_PRES,OCN_ATM_SFLX_rad_dn) &
    !$omp shared(OCN_ATM_cosSZA,OCN_ATM_SFLX_water,OCN_ATM_SFLX_ENGI,QV,PBL,SFC_DENS,SFC_PRES) &
    !$omp shared(SFLX_rad_dn,cosSZA,SFLX_water,SFLX_ENGI) &
    !$omp shared(LND_ATM_TEMP,LND_ATM_PRES,LND_ATM_W,LND_ATM_U,LND_ATM_V,LND_ATM_DENS) &
    !$omp shared(LND_ATM_QV,LND_ATM_PBL,LND_ATM_SFC_DENS,LND_ATM_SFC_PRES,LND_ATM_SFLX_rad_dn,LND_ATM_cosSZA) &
    !$omp shared(LND_ATM_SFLX_water,LND_ATM_SFLX_ENGI) &
    !$omp shared(URB_ATM_TEMP,URB_ATM_PRES,URB_ATM_W,URB_ATM_U,URB_ATM_V,URB_ATM_DENS) &
    !$omp shared(URB_ATM_QV,URB_ATM_PBL,URB_ATM_SFC_DENS,URB_ATM_SFC_PRES,URB_ATM_SFLX_rad_dn,URB_ATM_cosSZA) &
    !$omp shared(URB_ATM_SFLX_water,URB_ATM_SFLX_ENGI,CNT_putATM_URB,CNT_putATM_LND)
    do j = JS, JE
    do i = IS, IE
       ! for ocean
       OCN_ATM_TEMP     (i,j) = OCN_ATM_TEMP     (i,j) * CNT_putATM_OCN + TEMP     (i,j)
       OCN_ATM_PRES     (i,j) = OCN_ATM_PRES     (i,j) * CNT_putATM_OCN + PRES     (i,j)
       OCN_ATM_W        (i,j) = OCN_ATM_W        (i,j) * CNT_putATM_OCN + W        (i,j)
       OCN_ATM_U        (i,j) = OCN_ATM_U        (i,j) * CNT_putATM_OCN + U        (i,j)
       OCN_ATM_V        (i,j) = OCN_ATM_V        (i,j) * CNT_putATM_OCN + V        (i,j)
       OCN_ATM_DENS     (i,j) = OCN_ATM_DENS     (i,j) * CNT_putATM_OCN + DENS     (i,j)
       if ( .NOT. ATMOS_HYDROMETEOR_dry ) &
       OCN_ATM_QV       (i,j) = OCN_ATM_QV       (i,j) * CNT_putATM_OCN + QV       (i,j)
       OCN_ATM_PBL      (i,j) = OCN_ATM_PBL      (i,j) * CNT_putATM_OCN + PBL      (i,j)
       OCN_ATM_SFC_DENS (i,j) = OCN_ATM_SFC_DENS (i,j) * CNT_putATM_OCN + SFC_DENS (i,j)
       OCN_ATM_SFC_PRES (i,j) = OCN_ATM_SFC_PRES (i,j) * CNT_putATM_OCN + SFC_PRES (i,j)
       OCN_ATM_cosSZA   (i,j) = OCN_ATM_cosSZA   (i,j) * CNT_putATM_OCN + cosSZA   (i,j)
       OCN_ATM_SFLX_water(i,j) = OCN_ATM_SFLX_water(i,j) * CNT_putATM_OCN + SFLX_water(i,j)
       OCN_ATM_SFLX_ENGI(i,j) = OCN_ATM_SFLX_ENGI(i,j) * CNT_putATM_OCN  + SFLX_ENGI(i,j)
       do irgn = I_R_IR, I_R_VIS
       do idir = I_R_direct, I_R_diffuse
          OCN_ATM_SFLX_rad_dn(i,j,idir,irgn) = OCN_ATM_SFLX_rad_dn(i,j,idir,irgn) * CNT_putATM_OCN + SFLX_rad_dn(i,j,idir,irgn)
       enddo
       enddo

       ! for land
       LND_ATM_TEMP     (i,j) = LND_ATM_TEMP     (i,j) * CNT_putATM_LND + TEMP     (i,j)
       LND_ATM_PRES     (i,j) = LND_ATM_PRES     (i,j) * CNT_putATM_LND + PRES     (i,j)
       LND_ATM_W        (i,j) = LND_ATM_W        (i,j) * CNT_putATM_LND + W        (i,j)
       LND_ATM_U        (i,j) = LND_ATM_U        (i,j) * CNT_putATM_LND + U        (i,j)
       LND_ATM_V        (i,j) = LND_ATM_V        (i,j) * CNT_putATM_LND + V        (i,j)
       LND_ATM_DENS     (i,j) = LND_ATM_DENS     (i,j) * CNT_putATM_LND + DENS     (i,j)
       if ( .NOT. ATMOS_HYDROMETEOR_dry ) &
       LND_ATM_QV       (i,j) = LND_ATM_QV       (i,j) * CNT_putATM_LND + QV       (i,j)
       LND_ATM_PBL      (i,j) = LND_ATM_PBL      (i,j) * CNT_putATM_LND + PBL      (i,j)
       LND_ATM_SFC_DENS (i,j) = LND_ATM_SFC_DENS (i,j) * CNT_putATM_LND + SFC_DENS (i,j)
       LND_ATM_SFC_PRES (i,j) = LND_ATM_SFC_PRES (i,j) * CNT_putATM_LND + SFC_PRES (i,j)
       LND_ATM_cosSZA   (i,j) = LND_ATM_cosSZA   (i,j) * CNT_putATM_LND + cosSZA   (i,j)
       LND_ATM_SFLX_water(i,j) = LND_ATM_SFLX_water(i,j) * CNT_putATM_LND + SFLX_water(i,j)
       LND_ATM_SFLX_ENGI(i,j) = LND_ATM_SFLX_ENGI(i,j) * CNT_putATM_LND + SFLX_ENGI(i,j)
       do irgn = I_R_IR, I_R_VIS
       do idir = I_R_direct, I_R_diffuse
          LND_ATM_SFLX_rad_dn(i,j,idir,irgn) = LND_ATM_SFLX_rad_dn(i,j,idir,irgn) * CNT_putATM_LND + SFLX_rad_dn(i,j,idir,irgn)
       enddo
       enddo

       ! for urban
       URB_ATM_TEMP     (i,j) = URB_ATM_TEMP     (i,j) * CNT_putATM_URB + TEMP     (i,j)
       URB_ATM_PRES     (i,j) = URB_ATM_PRES     (i,j) * CNT_putATM_URB + PRES     (i,j)
       URB_ATM_W        (i,j) = URB_ATM_W        (i,j) * CNT_putATM_URB + W        (i,j)
       URB_ATM_U        (i,j) = URB_ATM_U        (i,j) * CNT_putATM_URB + U        (i,j)
       URB_ATM_V        (i,j) = URB_ATM_V        (i,j) * CNT_putATM_URB + V        (i,j)
       URB_ATM_DENS     (i,j) = URB_ATM_DENS     (i,j) * CNT_putATM_URB + DENS     (i,j)
       if ( .NOT. ATMOS_HYDROMETEOR_dry ) &
       URB_ATM_QV       (i,j) = URB_ATM_QV       (i,j) * CNT_putATM_URB + QV       (i,j)
       URB_ATM_PBL      (i,j) = URB_ATM_PBL      (i,j) * CNT_putATM_URB + PBL      (i,j)
       URB_ATM_SFC_DENS (i,j) = URB_ATM_SFC_DENS (i,j) * CNT_putATM_URB + SFC_DENS (i,j)
       URB_ATM_SFC_PRES (i,j) = URB_ATM_SFC_PRES (i,j) * CNT_putATM_URB + SFC_PRES (i,j)
       URB_ATM_cosSZA   (i,j) = URB_ATM_cosSZA   (i,j) * CNT_putATM_URB + cosSZA   (i,j)
       URB_ATM_SFLX_water(i,j) = URB_ATM_SFLX_water(i,j) * CNT_putATM_URB + SFLX_water(i,j)
       URB_ATM_SFLX_ENGI(i,j) = URB_ATM_SFLX_ENGI(i,j) * CNT_putATM_URB + SFLX_ENGI(i,j)
       do irgn = I_R_IR, I_R_VIS
       do idir = I_R_direct, I_R_diffuse
          URB_ATM_SFLX_rad_dn(i,j,idir,irgn) = URB_ATM_SFLX_rad_dn(i,j,idir,irgn) * CNT_putATM_URB + SFLX_rad_dn(i,j,idir,irgn)
       enddo
       enddo


       ! for ocean
       OCN_ATM_TEMP     (i,j) = OCN_ATM_TEMP     (i,j) / ( CNT_putATM_OCN + 1.0_RP )
       OCN_ATM_PRES     (i,j) = OCN_ATM_PRES     (i,j) / ( CNT_putATM_OCN + 1.0_RP )
       OCN_ATM_W        (i,j) = OCN_ATM_W        (i,j) / ( CNT_putATM_OCN + 1.0_RP )
       OCN_ATM_U        (i,j) = OCN_ATM_U        (i,j) / ( CNT_putATM_OCN + 1.0_RP )
       OCN_ATM_V        (i,j) = OCN_ATM_V        (i,j) / ( CNT_putATM_OCN + 1.0_RP )
       OCN_ATM_DENS     (i,j) = OCN_ATM_DENS     (i,j) / ( CNT_putATM_OCN + 1.0_RP )
       OCN_ATM_QV       (i,j) = OCN_ATM_QV       (i,j) / ( CNT_putATM_OCN + 1.0_RP )
       OCN_ATM_PBL      (i,j) = OCN_ATM_PBL      (i,j) / ( CNT_putATM_OCN + 1.0_RP )
       OCN_ATM_SFC_DENS (i,j) = OCN_ATM_SFC_DENS (i,j) / ( CNT_putATM_OCN + 1.0_RP )
       OCN_ATM_SFC_PRES (i,j) = OCN_ATM_SFC_PRES (i,j) / ( CNT_putATM_OCN + 1.0_RP )
       OCN_ATM_cosSZA   (i,j) = OCN_ATM_cosSZA   (i,j) / ( CNT_putATM_OCN + 1.0_RP )
       OCN_ATM_SFLX_water(i,j) = OCN_ATM_SFLX_water(i,j) / ( CNT_putATM_OCN + 1.0_RP )
       OCN_ATM_SFLX_ENGI(i,j) = OCN_ATM_SFLX_ENGI(i,j) / ( CNT_putATM_OCN + 1.0_RP )
       do irgn = I_R_IR, I_R_VIS
       do idir = I_R_direct, I_R_diffuse
          OCN_ATM_SFLX_rad_dn(i,j,idir,irgn) = OCN_ATM_SFLX_rad_dn(i,j,idir,irgn) / ( CNT_putATM_OCN + 1.0_RP )
       enddo
       enddo

       ! for land
       LND_ATM_TEMP     (i,j) = LND_ATM_TEMP     (i,j) / ( CNT_putATM_LND + 1.0_RP )
       LND_ATM_PRES     (i,j) = LND_ATM_PRES     (i,j) / ( CNT_putATM_LND + 1.0_RP )
       LND_ATM_W        (i,j) = LND_ATM_W        (i,j) / ( CNT_putATM_LND + 1.0_RP )
       LND_ATM_U        (i,j) = LND_ATM_U        (i,j) / ( CNT_putATM_LND + 1.0_RP )
       LND_ATM_V        (i,j) = LND_ATM_V        (i,j) / ( CNT_putATM_LND + 1.0_RP )
       LND_ATM_DENS     (i,j) = LND_ATM_DENS     (i,j) / ( CNT_putATM_LND + 1.0_RP )
       LND_ATM_QV       (i,j) = LND_ATM_QV       (i,j) / ( CNT_putATM_LND + 1.0_RP )
       LND_ATM_PBL      (i,j) = LND_ATM_PBL      (i,j) / ( CNT_putATM_LND + 1.0_RP )
       LND_ATM_SFC_DENS (i,j) = LND_ATM_SFC_DENS (i,j) / ( CNT_putATM_LND + 1.0_RP )
       LND_ATM_SFC_PRES (i,j) = LND_ATM_SFC_PRES (i,j) / ( CNT_putATM_LND + 1.0_RP )
       LND_ATM_cosSZA   (i,j) = LND_ATM_cosSZA   (i,j) / ( CNT_putATM_LND + 1.0_RP )
       LND_ATM_SFLX_water(i,j) = LND_ATM_SFLX_water(i,j) / ( CNT_putATM_LND + 1.0_RP )
       LND_ATM_SFLX_ENGI(i,j) = LND_ATM_SFLX_ENGI(i,j) / ( CNT_putATM_LND + 1.0_RP )
       do irgn = I_R_IR, I_R_VIS
       do idir = I_R_direct, I_R_diffuse
          LND_ATM_SFLX_rad_dn(i,j,idir,irgn) = LND_ATM_SFLX_rad_dn(i,j,idir,irgn) / ( CNT_putATM_LND + 1.0_RP )
       enddo
       enddo

       ! for urban
       URB_ATM_TEMP     (i,j) = URB_ATM_TEMP     (i,j) / ( CNT_putATM_URB + 1.0_RP )
       URB_ATM_PRES     (i,j) = URB_ATM_PRES     (i,j) / ( CNT_putATM_URB + 1.0_RP )
       URB_ATM_W        (i,j) = URB_ATM_W        (i,j) / ( CNT_putATM_URB + 1.0_RP )
       URB_ATM_U        (i,j) = URB_ATM_U        (i,j) / ( CNT_putATM_URB + 1.0_RP )
       URB_ATM_V        (i,j) = URB_ATM_V        (i,j) / ( CNT_putATM_URB + 1.0_RP )
       URB_ATM_DENS     (i,j) = URB_ATM_DENS     (i,j) / ( CNT_putATM_URB + 1.0_RP )
       URB_ATM_QV       (i,j) = URB_ATM_QV       (i,j) / ( CNT_putATM_URB + 1.0_RP )
       URB_ATM_PBL      (i,j) = URB_ATM_PBL      (i,j) / ( CNT_putATM_URB + 1.0_RP )
       URB_ATM_SFC_DENS (i,j) = URB_ATM_SFC_DENS (i,j) / ( CNT_putATM_URB + 1.0_RP )
       URB_ATM_SFC_PRES (i,j) = URB_ATM_SFC_PRES (i,j) / ( CNT_putATM_URB + 1.0_RP )
       URB_ATM_cosSZA   (i,j) = URB_ATM_cosSZA   (i,j) / ( CNT_putATM_URB + 1.0_RP )
       URB_ATM_SFLX_water(i,j) = URB_ATM_SFLX_water(i,j) / ( CNT_putATM_URB + 1.0_RP )
       URB_ATM_SFLX_ENGI(i,j) = URB_ATM_SFLX_ENGI(i,j) / ( CNT_putATM_URB + 1.0_RP )
       do irgn = I_R_IR, I_R_VIS
       do idir = I_R_direct, I_R_diffuse
          URB_ATM_SFLX_rad_dn(i,j,idir,irgn) = URB_ATM_SFLX_rad_dn(i,j,idir,irgn) / ( CNT_putATM_URB + 1.0_RP )
       enddo
       enddo
    enddo
    enddo

    if( countup ) then
       CNT_putATM_OCN = CNT_putATM_OCN + 1.0_RP
       CNT_putATM_LND = CNT_putATM_LND + 1.0_RP
       CNT_putATM_URB = CNT_putATM_URB + 1.0_RP
    endif

    return
  end subroutine CPL_putATM

  !-----------------------------------------------------------------------------
  subroutine CPL_putOCN( &
       SFC_TEMP,   &
       SFC_albedo, &
       SFC_Z0M,    &
       SFC_Z0H,    &
       SFC_Z0E,    &
       SFLX_MW,    &
       SFLX_MU,    &
       SFLX_MV,    &
       SFLX_SH,    &
       SFLX_LH,    &
       SFLX_GH,    &
       SFLX_QTRC,  &
       U10,        &
       V10,        &
       T2,         &
       Q2,         &
       mask,       &
       countup     )
    implicit none

    real(RP), intent(in) :: SFC_TEMP  (IA,JA)
    real(RP), intent(in) :: SFC_albedo(IA,JA,N_RAD_DIR,N_RAD_RGN)
    real(RP), intent(in) :: SFC_Z0M   (IA,JA)
    real(RP), intent(in) :: SFC_Z0H   (IA,JA)
    real(RP), intent(in) :: SFC_Z0E   (IA,JA)
    real(RP), intent(in) :: SFLX_MW   (IA,JA)
    real(RP), intent(in) :: SFLX_MU   (IA,JA)
    real(RP), intent(in) :: SFLX_MV   (IA,JA)
    real(RP), intent(in) :: SFLX_SH   (IA,JA)
    real(RP), intent(in) :: SFLX_LH   (IA,JA)
    real(RP), intent(in) :: SFLX_GH   (IA,JA)
    real(RP), intent(in) :: SFLX_QTRC (IA,JA,QA)
    real(RP), intent(in) :: U10       (IA,JA)
    real(RP), intent(in) :: V10       (IA,JA)
    real(RP), intent(in) :: T2        (IA,JA)
    real(RP), intent(in) :: Q2        (IA,JA)
    logical,  intent(in) :: mask      (IA,JA)
    logical,  intent(in) :: countup

    integer :: i, j, iq, idir, irgn
    !---------------------------------------------------------------------------

    !$omp parallel do default(none) private(i,j,idir,irgn) OMP_SCHEDULE_ &
    !$omp shared(JS,JE,IS,IE,QA) &
    !$omp shared(OCN_SFC_TEMP,OCN_SFC_albedo,OCN_SFC_Z0M,OCN_SFC_Z0H,OCN_SFC_Z0E) &
    !$omp shared(OCN_SFLX_MW,OCN_SFLX_MU,OCN_SFLX_MV,OCN_SFLX_SH,OCN_SFLX_LH,OCN_SFLX_GH,OCN_SFLX_QTRC,OCN_SFLX_ENGI,OCN_U10,OCN_V10,OCN_T2,OCN_Q2) &
    !$omp shared(SFC_TEMP,SFC_albedo,SFC_Z0M,SFC_Z0H,SFC_Z0E,SFLX_MW,SFLX_MU,SFLX_MV,SFLX_SH,SFLX_LH,SFLX_GH,SFLX_QTRC,U10,V10,T2,Q2) &
    !$omp shared(mask) &
    !$omp shared(CNT_putOCN) &
    !$omp shared(TRACER_CV,TRACER_ENGI0)
    do j = JS, JE
    do i = IS, IE
    if ( mask(i,j) ) then
       OCN_SFC_TEMP (i,j)   = OCN_SFC_TEMP (i,j)   * CNT_putOCN + SFC_TEMP (i,j)
       OCN_SFC_Z0M  (i,j)   = OCN_SFC_Z0M  (i,j)   * CNT_putOCN + SFC_Z0M  (i,j)
       OCN_SFC_Z0H  (i,j)   = OCN_SFC_Z0H  (i,j)   * CNT_putOCN + SFC_Z0H  (i,j)
       OCN_SFC_Z0E  (i,j)   = OCN_SFC_Z0E  (i,j)   * CNT_putOCN + SFC_Z0E  (i,j)
       OCN_SFLX_MW  (i,j)   = OCN_SFLX_MW  (i,j)   * CNT_putOCN + SFLX_MW  (i,j)
       OCN_SFLX_MU  (i,j)   = OCN_SFLX_MU  (i,j)   * CNT_putOCN + SFLX_MU  (i,j)
       OCN_SFLX_MV  (i,j)   = OCN_SFLX_MV  (i,j)   * CNT_putOCN + SFLX_MV  (i,j)
       OCN_SFLX_SH  (i,j)   = OCN_SFLX_SH  (i,j)   * CNT_putOCN + SFLX_SH  (i,j)
       OCN_SFLX_LH  (i,j)   = OCN_SFLX_LH  (i,j)   * CNT_putOCN + SFLX_LH  (i,j)
       OCN_SFLX_GH  (i,j)   = OCN_SFLX_GH  (i,j)   * CNT_putOCN + SFLX_GH  (i,j)
       OCN_SFLX_ENGI(i,j)   = OCN_SFLX_ENGI(i,j)   * CNT_putOCN
       do iq = 1, QA
          OCN_SFLX_QTRC(i,j,iq) = OCN_SFLX_QTRC(i,j,iq) * CNT_putOCN + SFLX_QTRC(i,j,iq)
          OCN_SFLX_ENGI(i,j)    = OCN_SFLX_ENGI(i,j) + SFLX_QTRC(i,j,iq) * ( TRACER_CV(iq) * SFC_TEMP(i,j) + TRACER_ENGI0(iq) )
       end do
       OCN_U10      (i,j)   = OCN_U10      (i,j)   * CNT_putOCN + U10      (i,j)
       OCN_V10      (i,j)   = OCN_V10      (i,j)   * CNT_putOCN + V10      (i,j)
       OCN_T2       (i,j)   = OCN_T2       (i,j)   * CNT_putOCN + T2       (i,j)
       OCN_Q2       (i,j)   = OCN_Q2       (i,j)   * CNT_putOCN + Q2       (i,j)
       do irgn = I_R_IR, I_R_VIS
       do idir = I_R_direct, I_R_diffuse
          OCN_SFC_albedo(i,j,idir,irgn) = OCN_SFC_albedo(i,j,idir,irgn) * CNT_putOCN + SFC_albedo(i,j,idir,irgn)
       enddo
       enddo

       OCN_SFC_TEMP (i,j)   = OCN_SFC_TEMP (i,j)   / ( CNT_putOCN + 1.0_RP )
       OCN_SFC_Z0M  (i,j)   = OCN_SFC_Z0M  (i,j)   / ( CNT_putOCN + 1.0_RP )
       OCN_SFC_Z0H  (i,j)   = OCN_SFC_Z0H  (i,j)   / ( CNT_putOCN + 1.0_RP )
       OCN_SFC_Z0E  (i,j)   = OCN_SFC_Z0E  (i,j)   / ( CNT_putOCN + 1.0_RP )
       OCN_SFLX_MW  (i,j)   = OCN_SFLX_MW  (i,j)   / ( CNT_putOCN + 1.0_RP )
       OCN_SFLX_MU  (i,j)   = OCN_SFLX_MU  (i,j)   / ( CNT_putOCN + 1.0_RP )
       OCN_SFLX_MV  (i,j)   = OCN_SFLX_MV  (i,j)   / ( CNT_putOCN + 1.0_RP )
       OCN_SFLX_SH  (i,j)   = OCN_SFLX_SH  (i,j)   / ( CNT_putOCN + 1.0_RP )
       OCN_SFLX_LH  (i,j)   = OCN_SFLX_LH  (i,j)   / ( CNT_putOCN + 1.0_RP )
       OCN_SFLX_GH  (i,j)   = OCN_SFLX_GH  (i,j)   / ( CNT_putOCN + 1.0_RP )
       OCN_SFLX_QTRC(i,j,:) = OCN_SFLX_QTRC(i,j,:) / ( CNT_putOCN + 1.0_RP )
       OCN_SFLX_ENGI(i,j)   = OCN_SFLX_ENGI(i,j)   / ( CNT_putOCN + 1.0_RP )
       OCN_U10      (i,j)   = OCN_U10      (i,j)   / ( CNT_putOCN + 1.0_RP )
       OCN_V10      (i,j)   = OCN_V10      (i,j)   / ( CNT_putOCN + 1.0_RP )
       OCN_T2       (i,j)   = OCN_T2       (i,j)   / ( CNT_putOCN + 1.0_RP )
       OCN_Q2       (i,j)   = OCN_Q2       (i,j)   / ( CNT_putOCN + 1.0_RP )
       do irgn = I_R_IR, I_R_VIS
       do idir = I_R_direct, I_R_diffuse
          OCN_SFC_albedo(i,j,idir,irgn) = OCN_SFC_albedo(i,j,idir,irgn) / ( CNT_putOCN + 1.0_RP )
       enddo
       enddo
       end if
    enddo
    enddo

    if( countup ) then
       CNT_putOCN = CNT_putOCN + 1.0_RP
    endif

    return
  end subroutine CPL_putOCN

  !-----------------------------------------------------------------------------
  subroutine CPL_putLND( &
       SFC_TEMP,   &
       SFC_albedo, &
       SFC_Z0M,    &
       SFC_Z0H,    &
       SFC_Z0E,    &
       SFLX_MW,    &
       SFLX_MU,    &
       SFLX_MV,    &
       SFLX_SH,    &
       SFLX_LH,    &
       SFLX_GH,    &
       SFLX_QTRC,  &
       U10,        &
       V10,        &
       T2,         &
       Q2,         &
       mask,       &
       countup     )
    implicit none

    real(RP), intent(in) :: SFC_TEMP  (IA,JA)
    real(RP), intent(in) :: SFC_albedo(IA,JA,N_RAD_DIR,N_RAD_RGN)
    real(RP), intent(in) :: SFC_Z0M   (IA,JA)
    real(RP), intent(in) :: SFC_Z0H   (IA,JA)
    real(RP), intent(in) :: SFC_Z0E   (IA,JA)
    real(RP), intent(in) :: SFLX_MW   (IA,JA)
    real(RP), intent(in) :: SFLX_MU   (IA,JA)
    real(RP), intent(in) :: SFLX_MV   (IA,JA)
    real(RP), intent(in) :: SFLX_SH   (IA,JA)
    real(RP), intent(in) :: SFLX_LH   (IA,JA)
    real(RP), intent(in) :: SFLX_GH   (IA,JA)
    real(RP), intent(in) :: SFLX_QTRC (IA,JA,QA)
    real(RP), intent(in) :: U10       (IA,JA)
    real(RP), intent(in) :: V10       (IA,JA)
    real(RP), intent(in) :: T2        (IA,JA)
    real(RP), intent(in) :: Q2        (IA,JA)
    logical,  intent(in) :: mask      (IA,JA)
    logical,  intent(in) :: countup

    integer :: i, j, iq, idir, irgn
    !---------------------------------------------------------------------------

    !$omp parallel do default(none) &
    !$omp shared(JS,JE,IS,IE,QA) &
    !$omp shared(LND_SFC_TEMP,LND_SFC_albedo,LND_SFC_Z0M,LND_SFC_Z0H,LND_SFC_Z0E) &
    !$omp shared(LND_SFLX_MW,LND_SFLX_MU,LND_SFLX_MV,LND_SFLX_SH,LND_SFLX_LH,LND_SFLX_GH,LND_SFLX_QTRC,LND_SFLX_ENGI) &
    !$omp shared(LND_U10,LND_V10,LND_T2,LND_Q2,CNT_putLND,SFC_TEMP,SFC_albedo,SFC_Z0M,SFC_Z0H) &
    !$omp shared(SFC_Z0E,SFLX_MW,SFLX_MU,SFLX_MV,SFLX_SH,SFLX_LH,SFLX_GH,SFLX_QTRC,U10,V10,T2,Q2) &
    !$omp shared(mask) &
    !$omp shared(TRACER_CV,TRACER_ENGI0) &
    !$omp private(i,j,idir,irgn) OMP_SCHEDULE_
    do j = JS, JE
    do i = IS, IE
    if ( mask(i,j) ) then
       LND_SFC_TEMP (i,j)   = LND_SFC_TEMP (i,j)   * CNT_putLND + SFC_TEMP (i,j)
       LND_SFC_Z0M  (i,j)   = LND_SFC_Z0M  (i,j)   * CNT_putLND + SFC_Z0M  (i,j)
       LND_SFC_Z0H  (i,j)   = LND_SFC_Z0H  (i,j)   * CNT_putLND + SFC_Z0H  (i,j)
       LND_SFC_Z0E  (i,j)   = LND_SFC_Z0E  (i,j)   * CNT_putLND + SFC_Z0E  (i,j)
       LND_SFLX_MW  (i,j)   = LND_SFLX_MW  (i,j)   * CNT_putLND + SFLX_MW  (i,j)
       LND_SFLX_MU  (i,j)   = LND_SFLX_MU  (i,j)   * CNT_putLND + SFLX_MU  (i,j)
       LND_SFLX_MV  (i,j)   = LND_SFLX_MV  (i,j)   * CNT_putLND + SFLX_MV  (i,j)
       LND_SFLX_SH  (i,j)   = LND_SFLX_SH  (i,j)   * CNT_putLND + SFLX_SH  (i,j)
       LND_SFLX_LH  (i,j)   = LND_SFLX_LH  (i,j)   * CNT_putLND + SFLX_LH  (i,j)
       LND_SFLX_GH  (i,j)   = LND_SFLX_GH  (i,j)   * CNT_putLND + SFLX_GH  (i,j)
       LND_SFLX_ENGI(i,j)   = LND_SFLX_ENGI(i,j)   * CNT_putLND
       do iq = 1, QA
          LND_SFLX_QTRC(i,j,iq) = LND_SFLX_QTRC(i,j,iq) * CNT_putLND + SFLX_QTRC(i,j,iq)
          LND_SFLX_ENGI(i,j)    = LND_SFLX_ENGI(i,j) + SFLX_QTRC(i,j,iq) * ( TRACER_CV(iq) * SFC_TEMP(i,j) + TRACER_ENGI0(iq) )
       end do
       LND_U10      (i,j)   = LND_U10      (i,j)   * CNT_putLND + U10      (i,j)
       LND_V10      (i,j)   = LND_V10      (i,j)   * CNT_putLND + V10      (i,j)
       LND_T2       (i,j)   = LND_T2       (i,j)   * CNT_putLND + T2       (i,j)
       LND_Q2       (i,j)   = LND_Q2       (i,j)   * CNT_putLND + Q2       (i,j)
       do irgn = I_R_IR, I_R_VIS
       do idir = I_R_direct, I_R_diffuse
          LND_SFC_albedo(i,j,idir,irgn) = LND_SFC_albedo(i,j,idir,irgn) * CNT_putLND + SFC_albedo(i,j,idir,irgn)
       enddo
       enddo

       LND_SFC_TEMP (i,j)   = LND_SFC_TEMP (i,j)   / ( CNT_putLND + 1.0_RP )
       LND_SFC_Z0M  (i,j)   = LND_SFC_Z0M  (i,j)   / ( CNT_putLND + 1.0_RP )
       LND_SFC_Z0H  (i,j)   = LND_SFC_Z0H  (i,j)   / ( CNT_putLND + 1.0_RP )
       LND_SFC_Z0E  (i,j)   = LND_SFC_Z0E  (i,j)   / ( CNT_putLND + 1.0_RP )
       LND_SFLX_MW  (i,j)   = LND_SFLX_MW  (i,j)   / ( CNT_putLND + 1.0_RP )
       LND_SFLX_MU  (i,j)   = LND_SFLX_MU  (i,j)   / ( CNT_putLND + 1.0_RP )
       LND_SFLX_MV  (i,j)   = LND_SFLX_MV  (i,j)   / ( CNT_putLND + 1.0_RP )
       LND_SFLX_SH  (i,j)   = LND_SFLX_SH  (i,j)   / ( CNT_putLND + 1.0_RP )
       LND_SFLX_LH  (i,j)   = LND_SFLX_LH  (i,j)   / ( CNT_putLND + 1.0_RP )
       LND_SFLX_GH  (i,j)   = LND_SFLX_GH  (i,j)   / ( CNT_putLND + 1.0_RP )
       LND_SFLX_QTRC(i,j,:) = LND_SFLX_QTRC(i,j,:) / ( CNT_putLND + 1.0_RP )
       LND_SFLX_ENGI(i,j)   = LND_SFLX_ENGI(i,j)   / ( CNT_putLND + 1.0_RP )
       LND_U10      (i,j)   = LND_U10      (i,j)   / ( CNT_putLND + 1.0_RP )
       LND_V10      (i,j)   = LND_V10      (i,j)   / ( CNT_putLND + 1.0_RP )
       LND_T2       (i,j)   = LND_T2       (i,j)   / ( CNT_putLND + 1.0_RP )
       LND_Q2       (i,j)   = LND_Q2       (i,j)   / ( CNT_putLND + 1.0_RP )
       do irgn = I_R_IR, I_R_VIS
       do idir = I_R_direct, I_R_diffuse
          LND_SFC_albedo(i,j,idir,irgn) = LND_SFC_albedo(i,j,idir,irgn) / ( CNT_putLND + 1.0_RP )
       enddo
       enddo
       end if
    enddo
    enddo

    if( countup ) then
       CNT_putLND = CNT_putLND + 1.0_RP
    endif

    return
  end subroutine CPL_putLND

  !-----------------------------------------------------------------------------
  subroutine CPL_putURB( &
       SFC_TEMP,   &
       SFC_albedo, &
       SFC_Z0M,    &
       SFC_Z0H,    &
       SFC_Z0E,    &
       SFLX_MW,    &
       SFLX_MU,    &
       SFLX_MV,    &
       SFLX_SH,    &
       SFLX_LH,    &
       SFLX_SHEX,  &
       SFLX_QVEX,  &
       SFLX_GH,    &
       SFLX_QTRC,  &
       U10,        &
       V10,        &
       T2,         &
       Q2,         &
       mask,       &
       countup     )
    implicit none

    real(RP), intent(in) :: SFC_TEMP  (IA,JA)
    real(RP), intent(in) :: SFC_albedo(IA,JA,N_RAD_DIR,N_RAD_RGN)
    real(RP), intent(in) :: SFC_Z0M   (IA,JA)
    real(RP), intent(in) :: SFC_Z0H   (IA,JA)
    real(RP), intent(in) :: SFC_Z0E   (IA,JA)
    real(RP), intent(in) :: SFLX_MW   (IA,JA)
    real(RP), intent(in) :: SFLX_MU   (IA,JA)
    real(RP), intent(in) :: SFLX_MV   (IA,JA)
    real(RP), intent(in) :: SFLX_SH   (IA,JA)
    real(RP), intent(in) :: SFLX_LH   (IA,JA)
    real(RP), intent(in) :: SFLX_SHEX (IA,JA)
    real(RP), intent(in) :: SFLX_QVEX (IA,JA)
    real(RP), intent(in) :: SFLX_GH   (IA,JA)
    real(RP), intent(in) :: SFLX_QTRC (IA,JA,QA)
    real(RP), intent(in) :: U10       (IA,JA)
    real(RP), intent(in) :: V10       (IA,JA)
    real(RP), intent(in) :: T2        (IA,JA)
    real(RP), intent(in) :: Q2        (IA,JA)
    logical,  intent(in) :: mask      (IA,JA)
    logical,  intent(in) :: countup

    integer :: i, j, iq, idir, irgn
    !---------------------------------------------------------------------------

    !$omp parallel do default(none)  OMP_SCHEDULE_ &
    !$omp shared(JS,JE,IS,IE,QA, &
    !$omp        URB_SFC_TEMP,URB_SFC_albedo,URB_SFC_Z0M,URB_SFC_Z0H,URB_SFC_Z0E, &
    !$omp        URB_SFLX_MW,URB_SFLX_MU,URB_SFLX_MV,URB_SFLX_SH,URB_SFLX_LH,URB_SFLX_SHEX,URB_SFLX_QVEX,URB_SFLX_GH, &
    !$omp        URB_SFLX_QTRC,URB_SFLX_ENGI,URB_U10,URB_V10,URB_T2,URB_Q2,CNT_putURB, &
    !$omp        SFC_TEMP,SFC_albedo,SFC_Z0M,SFC_Z0H,SFC_Z0E,SFLX_MW,SFLX_MU,SFLX_MV,SFLX_SH,SFLX_LH,SFLX_SHEX,SFLX_QVEX,SFLX_GH,SFLX_QTRC,U10,V10,T2,Q2, &
    !$omp        mask, &
    !$omp        TRACER_CV,TRACER_ENGI0)
    do j = JS, JE
    do i = IS, IE
    if ( mask(i,j) ) then
       URB_SFC_TEMP (i,j)   = URB_SFC_TEMP (i,j)   * CNT_putURB + SFC_TEMP (i,j)
       URB_SFC_Z0M  (i,j)   = URB_SFC_Z0M  (i,j)   * CNT_putURB + SFC_Z0M  (i,j)
       URB_SFC_Z0H  (i,j)   = URB_SFC_Z0H  (i,j)   * CNT_putURB + SFC_Z0H  (i,j)
       URB_SFC_Z0E  (i,j)   = URB_SFC_Z0E  (i,j)   * CNT_putURB + SFC_Z0E  (i,j)
       URB_SFLX_MW  (i,j)   = URB_SFLX_MW  (i,j)   * CNT_putURB + SFLX_MW  (i,j)
       URB_SFLX_MU  (i,j)   = URB_SFLX_MU  (i,j)   * CNT_putURB + SFLX_MU  (i,j)
       URB_SFLX_MV  (i,j)   = URB_SFLX_MV  (i,j)   * CNT_putURB + SFLX_MV  (i,j)
       URB_SFLX_SH  (i,j)   = URB_SFLX_SH  (i,j)   * CNT_putURB + SFLX_SH  (i,j)
       URB_SFLX_LH  (i,j)   = URB_SFLX_LH  (i,j)   * CNT_putURB + SFLX_LH  (i,j)
       URB_SFLX_SHEX(i,j)   = URB_SFLX_SHEX(i,j)   * CNT_putURB + SFLX_SHEX(i,j)
       URB_SFLX_QVEX(i,j)   = URB_SFLX_QVEX(i,j)   * CNT_putURB + SFLX_QVEX(i,j)
       URB_SFLX_GH  (i,j)   = URB_SFLX_GH  (i,j)   * CNT_putURB + SFLX_GH  (i,j)
       URB_SFLX_ENGI(i,j)   = URB_SFLX_ENGI(i,j)   * CNT_putURB
       do iq = 1, QA
          URB_SFLX_QTRC(i,j,iq) = URB_SFLX_QTRC(i,j,iq) * CNT_putURB + SFLX_QTRC(i,j,iq)
          URB_SFLX_ENGI(i,j)    = URB_SFLX_ENGI(i,j) + SFLX_QTRC(i,j,iq) * ( TRACER_CV(iq) * SFC_TEMP(i,j) + TRACER_ENGI0(iq) )
       end do
       URB_U10      (i,j)   = URB_U10      (i,j)   * CNT_putURB + U10      (i,j)
       URB_V10      (i,j)   = URB_V10      (i,j)   * CNT_putURB + V10      (i,j)
       URB_T2       (i,j)   = URB_T2       (i,j)   * CNT_putURB + T2       (i,j)
       URB_Q2       (i,j)   = URB_Q2       (i,j)   * CNT_putURB + Q2       (i,j)
       do irgn = I_R_IR, I_R_VIS
       do idir = I_R_direct, I_R_diffuse
          URB_SFC_albedo(i,j,idir,irgn) = URB_SFC_albedo(i,j,idir,irgn) * CNT_putURB + SFC_albedo(i,j,idir,irgn)
       enddo
       enddo

       URB_SFC_TEMP (i,j)   = URB_SFC_TEMP (i,j)   / ( CNT_putURB + 1.0_RP )
       URB_SFC_Z0M  (i,j)   = URB_SFC_Z0M  (i,j)   / ( CNT_putURB + 1.0_RP )
       URB_SFC_Z0H  (i,j)   = URB_SFC_Z0H  (i,j)   / ( CNT_putURB + 1.0_RP )
       URB_SFC_Z0E  (i,j)   = URB_SFC_Z0E  (i,j)   / ( CNT_putURB + 1.0_RP )
       URB_SFLX_MW  (i,j)   = URB_SFLX_MW  (i,j)   / ( CNT_putURB + 1.0_RP )
       URB_SFLX_MU  (i,j)   = URB_SFLX_MU  (i,j)   / ( CNT_putURB + 1.0_RP )
       URB_SFLX_MV  (i,j)   = URB_SFLX_MV  (i,j)   / ( CNT_putURB + 1.0_RP )
       URB_SFLX_SH  (i,j)   = URB_SFLX_SH  (i,j)   / ( CNT_putURB + 1.0_RP )
       URB_SFLX_LH  (i,j)   = URB_SFLX_LH  (i,j)   / ( CNT_putURB + 1.0_RP )
       URB_SFLX_SHEX(i,j)   = URB_SFLX_SHEX(i,j)   / ( CNT_putURB + 1.0_RP )
       URB_SFLX_QVEX(i,j)   = URB_SFLX_QVEX(i,j)   / ( CNT_putURB + 1.0_RP )
       URB_SFLX_GH  (i,j)   = URB_SFLX_GH  (i,j)   / ( CNT_putURB + 1.0_RP )
       URB_SFLX_QTRC(i,j,:) = URB_SFLX_QTRC(i,j,:) / ( CNT_putURB + 1.0_RP )
       URB_SFLX_ENGI(i,j)   = URB_SFLX_ENGI(i,j)   / ( CNT_putURB + 1.0_RP )
       URB_U10      (i,j)   = URB_U10      (i,j)   / ( CNT_putURB + 1.0_RP )
       URB_V10      (i,j)   = URB_V10      (i,j)   / ( CNT_putURB + 1.0_RP )
       URB_T2       (i,j)   = URB_T2       (i,j)   / ( CNT_putURB + 1.0_RP )
       URB_Q2       (i,j)   = URB_Q2       (i,j)   / ( CNT_putURB + 1.0_RP )
       do irgn = I_R_IR, I_R_VIS
       do idir = I_R_direct, I_R_diffuse
          URB_SFC_albedo(i,j,idir,irgn) = URB_SFC_albedo(i,j,idir,irgn) / ( CNT_putURB + 1.0_RP )
       enddo
       enddo
       end if
    enddo
    enddo

    if( countup ) then
       CNT_putURB = CNT_putURB + 1.0_RP
    endif

    return
  end subroutine CPL_putURB

  !-----------------------------------------------------------------------------
  subroutine CPL_getSFC_ATM( &
       SFC_TEMP,   &
       SFC_albedo, &
       SFC_Z0M,    &
       SFC_Z0H,    &
       SFC_Z0E,    &
       SFLX_MW,    &
       SFLX_MU,    &
       SFLX_MV,    &
       SFLX_SH,    &
       SFLX_LH,    &
       SFLX_SHEX,  &
       SFLX_QVEX,  &
       SFLX_GH,    &
       SFLX_QTRC,  &
       SFLX_ENGI,  &
       U10,        &
       V10,        &
       T2,         &
       Q2          )
    use scale_landuse, only: &
       fact_ocean => LANDUSE_fact_ocean, &
       fact_land  => LANDUSE_fact_land,  &
       fact_urban => LANDUSE_fact_urban
    implicit none

    real(RP), intent(out) :: SFC_TEMP  (IA,JA)
    real(RP), intent(out) :: SFC_albedo(IA,JA,N_RAD_DIR,N_RAD_RGN)
    real(RP), intent(out) :: SFC_Z0M   (IA,JA)
    real(RP), intent(out) :: SFC_Z0H   (IA,JA)
    real(RP), intent(out) :: SFC_Z0E   (IA,JA)
    real(RP), intent(out) :: SFLX_MW   (IA,JA)
    real(RP), intent(out) :: SFLX_MU   (IA,JA)
    real(RP), intent(out) :: SFLX_MV   (IA,JA)
    real(RP), intent(out) :: SFLX_SH   (IA,JA)
    real(RP), intent(out) :: SFLX_LH   (IA,JA)
    real(RP), intent(out) :: SFLX_SHEX (IA,JA)
    real(RP), intent(out) :: SFLX_QVEX (IA,JA)
    real(RP), intent(out) :: SFLX_GH   (IA,JA)
    real(RP), intent(out) :: SFLX_QTRC (IA,JA,QA)
    real(RP), intent(out) :: SFLX_ENGI (IA,JA)
    real(RP), intent(out) :: U10       (IA,JA)
    real(RP), intent(out) :: V10       (IA,JA)
    real(RP), intent(out) :: T2        (IA,JA)
    real(RP), intent(out) :: Q2        (IA,JA)

    integer :: i, j, idir, irgn, iq
    !---------------------------------------------------------------------------

    !$omp parallel do default(none) &
    !$omp shared(JS,JE,IS,IE,QA,SFLX_QTRC,SFLX_ENGI,SFC_TEMP,SFC_albedo,SFC_Z0M,SFC_Z0H,SFC_Z0E) &
    !$omp shared(SFLX_MW,SFLX_MU,SFLX_MV,SFLX_SH,SFLX_LH,SFLX_SHEX,SFLX_QVEX,SFLX_GH,U10,V10,T2,Q2) &
    !$omp shared(fact_ocean,fact_land,fact_urban,OCN_SFC_TEMP,LND_SFC_TEMP,URB_SFC_TEMP,OCN_SFC_albedo) &
    !$omp shared(LND_SFC_albedo,URB_SFC_albedo,OCN_SFC_Z0M,LND_SFC_Z0M,URB_SFC_Z0M) &
    !$omp shared(OCN_SFC_Z0H,LND_SFC_Z0H,URB_SFC_Z0H,OCN_SFC_Z0E,LND_SFC_Z0E,URB_SFC_Z0E,OCN_SFLX_MW) &
    !$omp shared(LND_SFLX_MW,URB_SFLX_MW,OCN_SFLX_MU,LND_SFLX_MU,URB_SFLX_MU,OCN_SFLX_MV,LND_SFLX_MV) &
    !$omp shared(URB_SFLX_MV,OCN_SFLX_SH,LND_SFLX_SH,URB_SFLX_SH,URB_SFLX_SHEX,OCN_SFLX_LH,LND_SFLX_LH,URB_SFLX_LH,URB_SFLX_QVEX) &
    !$omp shared(OCN_SFLX_GH,LND_SFLX_GH,URB_SFLX_GH,OCN_SFLX_QTRC,LND_SFLX_QTRC,URB_SFLX_QTRC,OCN_SFLX_ENGI,LND_SFLX_ENGI,URB_SFLX_ENGI) &
    !$omp shared(OCN_U10,LND_U10,URB_U10,OCN_V10,LND_V10,URB_V10,OCN_T2,LND_T2,URB_T2,OCN_Q2,LND_Q2,URB_Q2) &
    !$omp private(i,j,iq) OMP_SCHEDULE_
    do j = JS, JE
    do i = IS, IE
       SFC_TEMP (i,j)      =   fact_ocean(i,j) * OCN_SFC_TEMP (i,j) &
                             + fact_land (i,j) * LND_SFC_TEMP (i,j) &
                             + fact_urban(i,j) * URB_SFC_TEMP (i,j)

       do irgn = I_R_IR, I_R_VIS
       do idir = I_R_direct, I_R_diffuse
          SFC_albedo(i,j,idir,irgn) = fact_ocean(i,j) * OCN_SFC_albedo(i,j,idir,irgn) &
                                    + fact_land (i,j) * LND_SFC_albedo(i,j,idir,irgn) &
                                    + fact_urban(i,j) * URB_SFC_albedo(i,j,idir,irgn)
       enddo
       enddo

       SFC_Z0M  (i,j)      =   fact_ocean(i,j) * OCN_SFC_Z0M  (i,j) &
                             + fact_land (i,j) * LND_SFC_Z0M  (i,j) &
                             + fact_urban(i,j) * URB_SFC_Z0M  (i,j)

       SFC_Z0H  (i,j)      =   fact_ocean(i,j) * OCN_SFC_Z0H  (i,j) &
                             + fact_land (i,j) * LND_SFC_Z0H  (i,j) &
                             + fact_urban(i,j) * URB_SFC_Z0H  (i,j)

       SFC_Z0E  (i,j)      =   fact_ocean(i,j) * OCN_SFC_Z0E  (i,j) &
                             + fact_land (i,j) * LND_SFC_Z0E  (i,j) &
                             + fact_urban(i,j) * URB_SFC_Z0E  (i,j)

       SFLX_MW  (i,j)      =   fact_ocean(i,j) * OCN_SFLX_MW  (i,j) &
                             + fact_land (i,j) * LND_SFLX_MW  (i,j) &
                             + fact_urban(i,j) * URB_SFLX_MW  (i,j)

       SFLX_MU  (i,j)      =   fact_ocean(i,j) * OCN_SFLX_MU  (i,j) &
                             + fact_land (i,j) * LND_SFLX_MU  (i,j) &
                             + fact_urban(i,j) * URB_SFLX_MU  (i,j)

       SFLX_MV  (i,j)      =   fact_ocean(i,j) * OCN_SFLX_MV  (i,j) &
                             + fact_land (i,j) * LND_SFLX_MV  (i,j) &
                             + fact_urban(i,j) * URB_SFLX_MV  (i,j)

       SFLX_SH  (i,j)      =   fact_ocean(i,j) * OCN_SFLX_SH  (i,j) &
                             + fact_land (i,j) * LND_SFLX_SH  (i,j) &
                             + fact_urban(i,j) * URB_SFLX_SH  (i,j)

       SFLX_LH  (i,j)      =   fact_ocean(i,j) * OCN_SFLX_LH  (i,j) &
                             + fact_land (i,j) * LND_SFLX_LH  (i,j) &
                             + fact_urban(i,j) * URB_SFLX_LH  (i,j)

       SFLX_SHEX(i,j)      = + fact_urban(i,j) * URB_SFLX_SHEX(i,j)
            
       SFLX_QVEX(i,j)      = + fact_urban(i,j) * URB_SFLX_QVEX(i,j)

       SFLX_GH  (i,j)      =   fact_ocean(i,j) * OCN_SFLX_GH  (i,j) &
                             + fact_land (i,j) * LND_SFLX_GH  (i,j) &
                             + fact_urban(i,j) * URB_SFLX_GH  (i,j)

       do iq = 1, QA
          SFLX_QTRC(i,j,iq) =   fact_ocean(i,j) * OCN_SFLX_QTRC(i,j,iq) &
                              + fact_land (i,j) * LND_SFLX_QTRC(i,j,iq) &
                              + fact_urban(i,j) * URB_SFLX_QTRC(i,j,iq)
       enddo

       SFLX_ENGI(i,j)      =   fact_ocean(i,j) * OCN_SFLX_ENGI(i,j) &
                             + fact_land (i,j) * LND_SFLX_ENGI(i,j) &
                             + fact_urban(i,j) * URB_SFLX_ENGI(i,j)

       U10      (i,j)      =   fact_ocean(i,j) * OCN_U10      (i,j) &
                             + fact_land (i,j) * LND_U10      (i,j) &
                             + fact_urban(i,j) * URB_U10      (i,j)

       V10      (i,j)      =   fact_ocean(i,j) * OCN_V10      (i,j) &
                             + fact_land (i,j) * LND_V10      (i,j) &
                             + fact_urban(i,j) * URB_V10      (i,j)

       T2       (i,j)      =   fact_ocean(i,j) * OCN_T2       (i,j) &
                             + fact_land (i,j) * LND_T2       (i,j) &
                             + fact_urban(i,j) * URB_T2       (i,j)

       Q2       (i,j)      =   fact_ocean(i,j) * OCN_Q2       (i,j) &
                             + fact_land (i,j) * LND_Q2       (i,j) &
                             + fact_urban(i,j) * URB_Q2       (i,j)
    enddo
    enddo

    CNT_putOCN = 0.0_RP
    CNT_putLND = 0.0_RP
    CNT_putURB = 0.0_RP

    return
  end subroutine CPL_getSFC_ATM

  !-----------------------------------------------------------------------------
  subroutine CPL_getATM_OCN( &
       TEMP,        &
       PRES,        &
       W,           &
       U,           &
       V,           &
       DENS,        &
       QV,          &
       PBL,         &
       SFC_DENS,    &
       SFC_PRES,    &
       SFLX_rad_dn, &
       cosSZA,      &
       SFLX_water,   &
       SFLX_ENGI    )
    implicit none

    real(RP), intent(out) :: TEMP       (IA,JA)
    real(RP), intent(out) :: PRES       (IA,JA)
    real(RP), intent(out) :: W          (IA,JA)
    real(RP), intent(out) :: U          (IA,JA)
    real(RP), intent(out) :: V          (IA,JA)
    real(RP), intent(out) :: DENS       (IA,JA)
    real(RP), intent(out) :: QV         (IA,JA)
    real(RP), intent(out) :: PBL        (IA,JA)
    real(RP), intent(out) :: SFC_DENS   (IA,JA)
    real(RP), intent(out) :: SFC_PRES   (IA,JA)
    real(RP), intent(out) :: SFLX_rad_dn(IA,JA,N_RAD_DIR,N_RAD_RGN)
    real(RP), intent(out) :: cosSZA     (IA,JA)
    real(RP), intent(out) :: SFLX_water (IA,JA)
    real(RP), intent(out) :: SFLX_ENGI  (IA,JA)

    integer :: i, j, idir, irgn
    !---------------------------------------------------------------------------

!OCL XFILL
    !$omp parallel do default(none) private(i,j,idir,irgn) OMP_SCHEDULE_ &
    !$omp shared(JS,JE,IS,IE,TEMP,PRES,W,U,V,DENS,QV,PBL,SFC_DENS,SFC_PRES,SFLX_rad_dn,cosSZA,SFLX_water,SFLX_ENGI) &
    !$omp shared(OCN_ATM_TEMP,OCN_ATM_PRES,OCN_ATM_W,OCN_ATM_U,OCN_ATM_V,OCN_ATM_DENS,OCN_ATM_QV) &
    !$omp shared(OCN_ATM_PBL,OCN_ATM_SFC_DENS,OCN_ATM_SFC_PRES,OCN_ATM_SFLX_rad_dn,OCN_ATM_cosSZA,OCN_ATM_SFLX_water,OCN_ATM_SFLX_ENGI)
    do j = JS, JE
    do i = IS, IE
       TEMP     (i,j) = OCN_ATM_TEMP     (i,j)
       PRES     (i,j) = OCN_ATM_PRES     (i,j)
       W        (i,j) = OCN_ATM_W        (i,j)
       U        (i,j) = OCN_ATM_U        (i,j)
       V        (i,j) = OCN_ATM_V        (i,j)
       DENS     (i,j) = OCN_ATM_DENS     (i,j)
       QV       (i,j) = OCN_ATM_QV       (i,j)
       PBL      (i,j) = OCN_ATM_PBL      (i,j)
       SFC_DENS (i,j) = OCN_ATM_SFC_DENS (i,j)
       SFC_PRES (i,j) = OCN_ATM_SFC_PRES (i,j)
       cosSZA   (i,j) = OCN_ATM_cosSZA   (i,j)
       SFLX_water(i,j) = OCN_ATM_SFLX_water(i,j)
       SFLX_ENGI(i,j) = OCN_ATM_SFLX_ENGI(i,j)
       do irgn = I_R_IR, I_R_VIS
       do idir = I_R_direct, I_R_diffuse
          SFLX_rad_dn(i,j,idir,irgn) = OCN_ATM_SFLX_rad_dn(i,j,idir,irgn)
       enddo
       enddo
    enddo
    enddo

    CNT_putATM_OCN = 0.0_RP

    return
  end subroutine CPL_getATM_OCN

  !-----------------------------------------------------------------------------
  subroutine CPL_getATM_LND( &
       TEMP,        &
       PRES,        &
       W,           &
       U,           &
       V,           &
       DENS,        &
       QV,          &
       PBL,         &
       SFC_DENS,    &
       SFC_PRES,    &
       SFLX_rad_dn, &
       cosSZA,      &
       SFLX_water,  &
       SFLX_ENGI    )
    implicit none

    real(RP), intent(out) :: TEMP       (IA,JA)
    real(RP), intent(out) :: PRES       (IA,JA)
    real(RP), intent(out) :: W          (IA,JA)
    real(RP), intent(out) :: U          (IA,JA)
    real(RP), intent(out) :: V          (IA,JA)
    real(RP), intent(out) :: DENS       (IA,JA)
    real(RP), intent(out) :: QV         (IA,JA)
    real(RP), intent(out) :: PBL        (IA,JA)
    real(RP), intent(out) :: SFC_DENS   (IA,JA)
    real(RP), intent(out) :: SFC_PRES   (IA,JA)
    real(RP), intent(out) :: SFLX_rad_dn(IA,JA,N_RAD_DIR,N_RAD_RGN)
    real(RP), intent(out) :: cosSZA     (IA,JA)
    real(RP), intent(out) :: SFLX_water (IA,JA)
    real(RP), intent(out) :: SFLX_ENGI  (IA,JA)

    integer :: i, j, idir, irgn
    !---------------------------------------------------------------------------

!OCL XFILL
    !$omp parallel do default(none) private(i,j,idir,irgn) OMP_SCHEDULE_ &
    !$omp shared(JS,JE,IS,IE,TEMP,PRES,W,U,V,DENS,QV,PBL,SFC_DENS,SFC_PRES,SFLX_rad_dn,cosSZA,SFLX_water,SFLX_ENGI) &
    !$omp shared(LND_ATM_TEMP,LND_ATM_PRES,LND_ATM_W,LND_ATM_U,LND_ATM_V,LND_ATM_DENS,LND_ATM_QV) &
    !$omp shared(LND_ATM_PBL,LND_ATM_SFC_DENS,LND_ATM_SFC_PRES,LND_ATM_SFLX_rad_dn,LND_ATM_cosSZA,LND_ATM_SFLX_water,LND_ATM_SFLX_ENGI)
    do j = JS, JE
    do i = IS, IE
       TEMP     (i,j) = LND_ATM_TEMP     (i,j)
       PRES     (i,j) = LND_ATM_PRES     (i,j)
       W        (i,j) = LND_ATM_W        (i,j)
       U        (i,j) = LND_ATM_U        (i,j)
       V        (i,j) = LND_ATM_V        (i,j)
       DENS     (i,j) = LND_ATM_DENS     (i,j)
       QV       (i,j) = LND_ATM_QV       (i,j)
       PBL      (i,j) = LND_ATM_PBL      (i,j)
       SFC_DENS (i,j) = LND_ATM_SFC_DENS (i,j)
       SFC_PRES (i,j) = LND_ATM_SFC_PRES (i,j)
       cosSZA   (i,j) = LND_ATM_cosSZA   (i,j)
       SFLX_water(i,j) = LND_ATM_SFLX_water(i,j)
       SFLX_ENGI(i,j) = LND_ATM_SFLX_ENGI(i,j)
       do irgn = I_R_IR, I_R_VIS
       do idir = I_R_direct, I_R_diffuse
          SFLX_rad_dn(i,j,idir,irgn) = LND_ATM_SFLX_rad_dn(i,j,idir,irgn)
       enddo
       enddo
    enddo
    enddo

    CNT_putATM_LND = 0.0_RP

    return
  end subroutine CPL_getATM_LND

  !-----------------------------------------------------------------------------
  subroutine CPL_getATM_URB( &
       TEMP,        &
       PRES,        &
       W,           &
       U,           &
       V,           &
       DENS,        &
       QV,          &
       PBL,         &
       SFC_DENS,    &
       SFC_PRES,    &
       SFLX_rad_dn, &
       cosSZA,      &
       SFLX_water,  &
       SFLX_ENGI    )
    implicit none

    real(RP), intent(out) :: TEMP       (IA,JA)
    real(RP), intent(out) :: PRES       (IA,JA)
    real(RP), intent(out) :: W          (IA,JA)
    real(RP), intent(out) :: U          (IA,JA)
    real(RP), intent(out) :: V          (IA,JA)
    real(RP), intent(out) :: DENS       (IA,JA)
    real(RP), intent(out) :: QV         (IA,JA)
    real(RP), intent(out) :: PBL        (IA,JA)
    real(RP), intent(out) :: SFC_DENS   (IA,JA)
    real(RP), intent(out) :: SFC_PRES   (IA,JA)
    real(RP), intent(out) :: SFLX_rad_dn(IA,JA,N_RAD_DIR,N_RAD_RGN)
    real(RP), intent(out) :: cosSZA     (IA,JA)
    real(RP), intent(out) :: SFLX_water (IA,JA)
    real(RP), intent(out) :: SFLX_ENGI  (IA,JA)

    integer :: i, j, idir, irgn
    !---------------------------------------------------------------------------

!OCL XFILL
    !$omp parallel do default(none) private(i,j,idir,irgn) OMP_SCHEDULE_ &
    !$omp shared(JS,JE,IS,IE,TEMP,PRES,W,U,V,DENS,QV,PBL,SFC_DENS,SFC_PRES,SFLX_rad_dn,cosSZA,SFLX_water,SFLX_ENGI) &
    !$omp shared(URB_ATM_TEMP,URB_ATM_PRES,URB_ATM_W,URB_ATM_U,URB_ATM_V,URB_ATM_DENS,URB_ATM_QV) &
    !$omp shared(URB_ATM_PBL,URB_ATM_SFC_DENS,URB_ATM_SFC_PRES,URB_ATM_SFLX_rad_dn,URB_ATM_cosSZA,URB_ATM_SFLX_water,URB_ATM_SFLX_ENGI)
    do j = JS, JE
    do i = IS, IE
       TEMP     (i,j) = URB_ATM_TEMP     (i,j)
       PRES     (i,j) = URB_ATM_PRES     (i,j)
       W        (i,j) = URB_ATM_W        (i,j)
       U        (i,j) = URB_ATM_U        (i,j)
       V        (i,j) = URB_ATM_V        (i,j)
       DENS     (i,j) = URB_ATM_DENS     (i,j)
       QV       (i,j) = URB_ATM_QV       (i,j)
       PBL      (i,j) = URB_ATM_PBL      (i,j)
       SFC_DENS (i,j) = URB_ATM_SFC_DENS (i,j)
       SFC_PRES (i,j) = URB_ATM_SFC_PRES (i,j)
       cosSZA   (i,j) = URB_ATM_cosSZA   (i,j)
       SFLX_water(i,j) = URB_ATM_SFLX_water(i,j)
       SFLX_ENGI(i,j) = URB_ATM_SFLX_ENGI(i,j)
       do irgn = I_R_IR, I_R_VIS
       do idir = I_R_direct, I_R_diffuse
          SFLX_rad_dn(i,j,idir,irgn) = URB_ATM_SFLX_rad_dn(i,j,idir,irgn)
       enddo
       enddo
    enddo
    enddo

    CNT_putATM_URB = 0.0_RP

    return
  end subroutine CPL_getATM_URB

end module mod_cpl_vars
