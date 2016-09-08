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
module mod_cpl_vars
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_stdio
  use scale_prof
  use scale_debug
  use scale_grid_index
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
  public :: CPL_vars_setup

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
  real(RP), public, allocatable :: OCN_SFC_TEMP  (:,:)   ! ocean surface skin temperature [K]
  real(RP), public, allocatable :: OCN_SFC_albedo(:,:,:) ! ocean surface albedo [0-1]
  real(RP), public, allocatable :: OCN_SFC_Z0M   (:,:)   ! ocean surface roughness length for momemtum [m]
  real(RP), public, allocatable :: OCN_SFC_Z0H   (:,:)   ! ocean surface roughness length for heat [m]
  real(RP), public, allocatable :: OCN_SFC_Z0E   (:,:)   ! ocean surface roughness length for vapor [m]
  real(RP), public, allocatable :: OCN_SFLX_MW   (:,:)   ! ocean surface w-momentum flux [kg/m2/s]
  real(RP), public, allocatable :: OCN_SFLX_MU   (:,:)   ! ocean surface u-momentum flux [kg/m2/s]
  real(RP), public, allocatable :: OCN_SFLX_MV   (:,:)   ! ocean surface v-momentum flux [kg/m2/s]
  real(RP), public, allocatable :: OCN_SFLX_SH   (:,:)   ! ocean surface sensible heat flux [J/m2/s]
  real(RP), public, allocatable :: OCN_SFLX_LH   (:,:)   ! ocean surface latent heat flux [J/m2/s]
  real(RP), public, allocatable :: OCN_SFLX_WH   (:,:)   ! ocean surface water heat flux [J/m2/s]
  real(RP), public, allocatable :: OCN_SFLX_evap (:,:)   ! ocean surface water vapor flux [kg/m2/s]
  real(RP), public, allocatable :: OCN_U10       (:,:)   ! ocean velocity u at 10m [m/s]
  real(RP), public, allocatable :: OCN_V10       (:,:)   ! ocean velocity v at 10m [m/s]
  real(RP), public, allocatable :: OCN_T2        (:,:)   ! ocean temperature at 2m [K]
  real(RP), public, allocatable :: OCN_Q2        (:,:)   ! ocean water vapor at 2m [kg/kg]

  ! Input from land model
  real(RP), public, allocatable :: LND_SFC_TEMP  (:,:)   ! land surface skin temperature [K]
  real(RP), public, allocatable :: LND_SFC_albedo(:,:,:) ! land surface albedo [0-1]
  real(RP), public, allocatable :: LND_SFC_Z0M   (:,:)   ! land surface roughness length for momemtum [m]
  real(RP), public, allocatable :: LND_SFC_Z0H   (:,:)   ! land surface roughness length for heat [m]
  real(RP), public, allocatable :: LND_SFC_Z0E   (:,:)   ! land surface roughness length for vapor [m]
  real(RP), public, allocatable :: LND_SFLX_MW   (:,:)   ! land surface w-momentum flux [kg/m2/s]
  real(RP), public, allocatable :: LND_SFLX_MU   (:,:)   ! land surface u-momentum flux [kg/m2/s]
  real(RP), public, allocatable :: LND_SFLX_MV   (:,:)   ! land surface v-momentum flux [kg/m2/s]
  real(RP), public, allocatable :: LND_SFLX_SH   (:,:)   ! land surface sensible heat flux [J/m2/s]
  real(RP), public, allocatable :: LND_SFLX_LH   (:,:)   ! land surface latent heat flux [J/m2/s]
  real(RP), public, allocatable :: LND_SFLX_GH   (:,:)   ! land surface ground heat flux [J/m2/s]
  real(RP), public, allocatable :: LND_SFLX_evap (:,:)   ! land surface water vapor flux [kg/m2/s]
  real(RP), public, allocatable :: LND_U10       (:,:)   ! land velocity u at 10m [m/s]
  real(RP), public, allocatable :: LND_V10       (:,:)   ! land velocity v at 10m [m/s]
  real(RP), public, allocatable :: LND_T2        (:,:)   ! land temperature at 2m [K]
  real(RP), public, allocatable :: LND_Q2        (:,:)   ! land water vapor at 2m [kg/kg]

  ! Input from urban model
  real(RP), public, allocatable :: URB_SFC_TEMP  (:,:)   ! urban surface skin temperature [K]
  real(RP), public, allocatable :: URB_SFC_albedo(:,:,:) ! urban surface albedo [0-1]
  real(RP), public, allocatable :: URB_SFC_Z0M   (:,:)   ! urban surface roughness length for momemtum [m]
  real(RP), public, allocatable :: URB_SFC_Z0H   (:,:)   ! urban surface roughness length for heat [m]
  real(RP), public, allocatable :: URB_SFC_Z0E   (:,:)   ! urban surface roughness length for vapor [m]
  real(RP), public, allocatable :: URB_SFLX_MW   (:,:)   ! urban surface w-momentum flux [kg/m2/s]
  real(RP), public, allocatable :: URB_SFLX_MU   (:,:)   ! urban surface u-momentum flux [kg/m2/s]
  real(RP), public, allocatable :: URB_SFLX_MV   (:,:)   ! urban surface v-momentum flux [kg/m2/s]
  real(RP), public, allocatable :: URB_SFLX_SH   (:,:)   ! urban surface sensible heat flux [J/m2/s]
  real(RP), public, allocatable :: URB_SFLX_LH   (:,:)   ! urban surface latent heat flux [J/m2/s]
  real(RP), public, allocatable :: URB_SFLX_GH   (:,:)   ! urban surface ground heat flux [J/m2/s]
  real(RP), public, allocatable :: URB_SFLX_evap (:,:)   ! urban surface water vapor flux [kg/m2/s]
  real(RP), public, allocatable :: URB_U10       (:,:)   ! urban velocity u at 10m [m/s]
  real(RP), public, allocatable :: URB_V10       (:,:)   ! urban velocity v at 10m [m/s]
  real(RP), public, allocatable :: URB_T2        (:,:)   ! urban temperature at 2m [K]
  real(RP), public, allocatable :: URB_Q2        (:,:)   ! urban water vapor at 2m [kg/kg]

  ! Output to ocean model
  real(RP), public, allocatable :: OCN_ATM_TEMP       (:,:)     ! temperature at the lowermost atmosphere layer [K]
  real(RP), public, allocatable :: OCN_ATM_PRES       (:,:)     ! pressure at the lowermost atmosphere layer [Pa]
  real(RP), public, allocatable :: OCN_ATM_W          (:,:)     ! velocity w at the lowermost atmosphere layer [m/s]
  real(RP), public, allocatable :: OCN_ATM_U          (:,:)     ! velocity u at the lowermost atmosphere layer [m/s]
  real(RP), public, allocatable :: OCN_ATM_V          (:,:)     ! velocity v at the lowermost atmosphere layer [m/s]
  real(RP), public, allocatable :: OCN_ATM_DENS       (:,:)     ! density at the lowermost atmosphere layer [kg/m3]
  real(RP), public, allocatable :: OCN_ATM_QV         (:,:)     ! water vapor at the lowermost atmosphere layer [kg/kg]
  real(RP), public, allocatable :: OCN_ATM_PBL        (:,:)     ! the top of atmospheric mixing layer [m]
  real(RP), public, allocatable :: OCN_ATM_SFC_PRES   (:,:)     ! surface pressure [Pa]
  real(RP), public, allocatable :: OCN_ATM_SFLX_rad_dn(:,:,:,:) ! downward radiation flux (SW/LW,direct/diffuse) [J/m2/s]
  real(RP), public, allocatable :: OCN_ATM_cosSZA     (:,:)     ! cos(solar zenith angle) [0-1]
  real(RP), public, allocatable :: OCN_ATM_SFLX_rain  (:,:)     ! liquid water flux [kg/m2/s]
  real(RP), public, allocatable :: OCN_ATM_SFLX_snow  (:,:)     ! ice    water flux [kg/m2/s]

  ! Output to land model
  real(RP), public, allocatable :: LND_ATM_TEMP       (:,:)     ! temperature at the lowermost atmosphere layer [K]
  real(RP), public, allocatable :: LND_ATM_PRES       (:,:)     ! pressure at the lowermost atmosphere layer [Pa]
  real(RP), public, allocatable :: LND_ATM_W          (:,:)     ! velocity w at the lowermost atmosphere layer [m/s]
  real(RP), public, allocatable :: LND_ATM_U          (:,:)     ! velocity u at the lowermost atmosphere layer [m/s]
  real(RP), public, allocatable :: LND_ATM_V          (:,:)     ! velocity v at the lowermost atmosphere layer [m/s]
  real(RP), public, allocatable :: LND_ATM_DENS       (:,:)     ! density at the lowermost atmosphere layer [kg/m3]
  real(RP), public, allocatable :: LND_ATM_QV         (:,:)     ! water vapor at the lowermost atmosphere layer [kg/kg]
  real(RP), public, allocatable :: LND_ATM_PBL        (:,:)     ! the top of atmospheric mixing layer [m]
  real(RP), public, allocatable :: LND_ATM_SFC_PRES   (:,:)     ! surface pressure [Pa]
  real(RP), public, allocatable :: LND_ATM_SFLX_rad_dn(:,:,:,:) ! downward radiation flux (SW/LW,direct/diffuse) [J/m2/s]
  real(RP), public, allocatable :: LND_ATM_cosSZA     (:,:)     ! cos(solar zenith angle) [0-1]
  real(RP), public, allocatable :: LND_ATM_SFLX_rain  (:,:)     ! liquid water flux [kg/m2/s]
  real(RP), public, allocatable :: LND_ATM_SFLX_snow  (:,:)     ! ice    water flux [kg/m2/s]

  ! Output to urban model
  real(RP), public, allocatable :: URB_ATM_TEMP       (:,:)     ! temperature at the lowermost atmosphere layer [K]
  real(RP), public, allocatable :: URB_ATM_PRES       (:,:)     ! pressure at the lowermost atmosphere layer [Pa]
  real(RP), public, allocatable :: URB_ATM_W          (:,:)     ! velocity w at the lowermost atmosphere layer [m/s]
  real(RP), public, allocatable :: URB_ATM_U          (:,:)     ! velocity u at the lowermost atmosphere layer [m/s]
  real(RP), public, allocatable :: URB_ATM_V          (:,:)     ! velocity v at the lowermost atmosphere layer [m/s]
  real(RP), public, allocatable :: URB_ATM_DENS       (:,:)     ! density at the lowermost atmosphere layer [kg/m3]
  real(RP), public, allocatable :: URB_ATM_QV         (:,:)     ! water vapor at the lowermost atmosphere layer [kg/kg]
  real(RP), public, allocatable :: URB_ATM_PBL        (:,:)     ! the top of atmospheric mixing layer [m]
  real(RP), public, allocatable :: URB_ATM_SFC_PRES   (:,:)     ! surface pressure [Pa]
  real(RP), public, allocatable :: URB_ATM_SFLX_rad_dn(:,:,:,:) ! downward radiation flux (SW/LW,direct/diffuse) [J/m2/s]
  real(RP), public, allocatable :: URB_ATM_cosSZA     (:,:)     ! cos(solar zenith angle) [0-1]
  real(RP), public, allocatable :: URB_ATM_SFLX_rain  (:,:)     ! liquid water flux [kg/m2/s]
  real(RP), public, allocatable :: URB_ATM_SFLX_snow  (:,:)     ! ice    water flux [kg/m2/s]

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
    implicit none
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '++++++ Module[VARS] / Categ[CPL] / Origin[SCALE-RM]'

    allocate( OCN_SFC_TEMP  (IA,JA)   )
    allocate( OCN_SFC_albedo(IA,JA,2) )
    allocate( OCN_SFC_Z0M   (IA,JA)   )
    allocate( OCN_SFC_Z0H   (IA,JA)   )
    allocate( OCN_SFC_Z0E   (IA,JA)   )
    allocate( OCN_SFLX_MU   (IA,JA)   )
    allocate( OCN_SFLX_MV   (IA,JA)   )
    allocate( OCN_SFLX_MW   (IA,JA)   )
    allocate( OCN_SFLX_SH   (IA,JA)   )
    allocate( OCN_SFLX_LH   (IA,JA)   )
    allocate( OCN_SFLX_WH   (IA,JA)   )
    allocate( OCN_SFLX_evap (IA,JA)   )
    allocate( OCN_U10       (IA,JA)   )
    allocate( OCN_V10       (IA,JA)   )
    allocate( OCN_T2        (IA,JA)   )
    allocate( OCN_Q2        (IA,JA)   )
    OCN_SFC_TEMP  (:,:)   = UNDEF
    OCN_SFC_albedo(:,:,:) = UNDEF
    OCN_SFC_Z0M   (:,:)   = UNDEF
    OCN_SFC_Z0H   (:,:)   = UNDEF
    OCN_SFC_Z0E   (:,:)   = UNDEF
    OCN_SFLX_MU   (:,:)   = UNDEF
    OCN_SFLX_MV   (:,:)   = UNDEF
    OCN_SFLX_MW   (:,:)   = UNDEF
    OCN_SFLX_SH   (:,:)   = UNDEF
    OCN_SFLX_LH   (:,:)   = UNDEF
    OCN_SFLX_WH   (:,:)   = UNDEF
    OCN_SFLX_evap (:,:)   = UNDEF
    OCN_U10       (:,:)   = UNDEF
    OCN_V10       (:,:)   = UNDEF
    OCN_T2        (:,:)   = UNDEF
    OCN_Q2        (:,:)   = UNDEF

    allocate( LND_SFC_TEMP  (IA,JA)   )
    allocate( LND_SFC_albedo(IA,JA,2) )
    allocate( LND_SFC_Z0M   (IA,JA)   )
    allocate( LND_SFC_Z0H   (IA,JA)   )
    allocate( LND_SFC_Z0E   (IA,JA)   )
    allocate( LND_SFLX_MU   (IA,JA)   )
    allocate( LND_SFLX_MV   (IA,JA)   )
    allocate( LND_SFLX_MW   (IA,JA)   )
    allocate( LND_SFLX_SH   (IA,JA)   )
    allocate( LND_SFLX_LH   (IA,JA)   )
    allocate( LND_SFLX_GH   (IA,JA)   )
    allocate( LND_SFLX_evap (IA,JA)   )
    allocate( LND_U10       (IA,JA)   )
    allocate( LND_V10       (IA,JA)   )
    allocate( LND_T2        (IA,JA)   )
    allocate( LND_Q2        (IA,JA)   )
    LND_SFC_TEMP  (:,:)   = UNDEF
    LND_SFC_albedo(:,:,:) = UNDEF
    LND_SFC_Z0M   (:,:)   = UNDEF
    LND_SFC_Z0H   (:,:)   = UNDEF
    LND_SFC_Z0E   (:,:)   = UNDEF
    LND_SFLX_MU   (:,:)   = UNDEF
    LND_SFLX_MV   (:,:)   = UNDEF
    LND_SFLX_MW   (:,:)   = UNDEF
    LND_SFLX_SH   (:,:)   = UNDEF
    LND_SFLX_LH   (:,:)   = UNDEF
    LND_SFLX_GH   (:,:)   = UNDEF
    LND_SFLX_evap (:,:)   = UNDEF
    LND_U10       (:,:)   = UNDEF
    LND_V10       (:,:)   = UNDEF
    LND_T2        (:,:)   = UNDEF
    LND_Q2        (:,:)   = UNDEF

    allocate( URB_SFC_TEMP  (IA,JA)   )
    allocate( URB_SFC_albedo(IA,JA,2) )
    allocate( URB_SFC_Z0M   (IA,JA)   )
    allocate( URB_SFC_Z0H   (IA,JA)   )
    allocate( URB_SFC_Z0E   (IA,JA)   )
    allocate( URB_SFLX_MU   (IA,JA)   )
    allocate( URB_SFLX_MV   (IA,JA)   )
    allocate( URB_SFLX_MW   (IA,JA)   )
    allocate( URB_SFLX_SH   (IA,JA)   )
    allocate( URB_SFLX_LH   (IA,JA)   )
    allocate( URB_SFLX_GH   (IA,JA)   )
    allocate( URB_SFLX_evap (IA,JA)   )
    allocate( URB_U10       (IA,JA)   )
    allocate( URB_V10       (IA,JA)   )
    allocate( URB_T2        (IA,JA)   )
    allocate( URB_Q2        (IA,JA)   )
    URB_SFC_TEMP  (:,:)   = UNDEF
    URB_SFC_albedo(:,:,:) = UNDEF
    URB_SFC_Z0M   (:,:)   = UNDEF
    URB_SFC_Z0H   (:,:)   = UNDEF
    URB_SFC_Z0E   (:,:)   = UNDEF
    URB_SFLX_MU   (:,:)   = UNDEF
    URB_SFLX_MV   (:,:)   = UNDEF
    URB_SFLX_MW   (:,:)   = UNDEF
    URB_SFLX_SH   (:,:)   = UNDEF
    URB_SFLX_LH   (:,:)   = UNDEF
    URB_SFLX_GH   (:,:)   = UNDEF
    URB_SFLX_evap (:,:)   = UNDEF
    URB_U10       (:,:)   = UNDEF
    URB_V10       (:,:)   = UNDEF
    URB_T2        (:,:)   = UNDEF
    URB_Q2        (:,:)   = UNDEF

    allocate( OCN_ATM_TEMP       (IA,JA)     )
    allocate( OCN_ATM_PRES       (IA,JA)     )
    allocate( OCN_ATM_W          (IA,JA)     )
    allocate( OCN_ATM_U          (IA,JA)     )
    allocate( OCN_ATM_V          (IA,JA)     )
    allocate( OCN_ATM_DENS       (IA,JA)     )
    allocate( OCN_ATM_QV         (IA,JA)     )
    allocate( OCN_ATM_PBL        (IA,JA)     )
    allocate( OCN_ATM_SFC_PRES   (IA,JA)     )
    allocate( OCN_ATM_SFLX_rad_dn(IA,JA,2,2) )
    allocate( OCN_ATM_cosSZA     (IA,JA)     )
    allocate( OCN_ATM_SFLX_rain  (IA,JA)     )
    allocate( OCN_ATM_SFLX_snow  (IA,JA)     )
    OCN_ATM_TEMP       (:,:)     = UNDEF
    OCN_ATM_PRES       (:,:)     = UNDEF
    OCN_ATM_W          (:,:)     = UNDEF
    OCN_ATM_U          (:,:)     = UNDEF
    OCN_ATM_V          (:,:)     = UNDEF
    OCN_ATM_DENS       (:,:)     = UNDEF
    OCN_ATM_QV         (:,:)     = UNDEF
    OCN_ATM_PBL        (:,:)     = UNDEF
    OCN_ATM_SFC_PRES   (:,:)     = UNDEF
    OCN_ATM_SFLX_rad_dn(:,:,:,:) = UNDEF
    OCN_ATM_cosSZA     (:,:)     = UNDEF
    OCN_ATM_SFLX_rain  (:,:)     = UNDEF
    OCN_ATM_SFLX_snow  (:,:)     = UNDEF

    allocate( LND_ATM_TEMP       (IA,JA)     )
    allocate( LND_ATM_PRES       (IA,JA)     )
    allocate( LND_ATM_W          (IA,JA)     )
    allocate( LND_ATM_U          (IA,JA)     )
    allocate( LND_ATM_V          (IA,JA)     )
    allocate( LND_ATM_DENS       (IA,JA)     )
    allocate( LND_ATM_QV         (IA,JA)     )
    allocate( LND_ATM_PBL        (IA,JA)     )
    allocate( LND_ATM_SFC_PRES   (IA,JA)     )
    allocate( LND_ATM_SFLX_rad_dn(IA,JA,2,2) )
    allocate( LND_ATM_cosSZA     (IA,JA)     )
    allocate( LND_ATM_SFLX_rain  (IA,JA)     )
    allocate( LND_ATM_SFLX_snow  (IA,JA)     )
    LND_ATM_TEMP       (:,:)     = UNDEF
    LND_ATM_PRES       (:,:)     = UNDEF
    LND_ATM_W          (:,:)     = UNDEF
    LND_ATM_U          (:,:)     = UNDEF
    LND_ATM_V          (:,:)     = UNDEF
    LND_ATM_DENS       (:,:)     = UNDEF
    LND_ATM_QV         (:,:)     = UNDEF
    LND_ATM_PBL        (:,:)     = UNDEF
    LND_ATM_SFC_PRES   (:,:)     = UNDEF
    LND_ATM_SFLX_rad_dn(:,:,:,:) = UNDEF
    LND_ATM_cosSZA     (:,:)     = UNDEF
    LND_ATM_SFLX_rain  (:,:)     = UNDEF
    LND_ATM_SFLX_snow  (:,:)     = UNDEF

    allocate( URB_ATM_TEMP       (IA,JA)     )
    allocate( URB_ATM_PRES       (IA,JA)     )
    allocate( URB_ATM_W          (IA,JA)     )
    allocate( URB_ATM_U          (IA,JA)     )
    allocate( URB_ATM_V          (IA,JA)     )
    allocate( URB_ATM_DENS       (IA,JA)     )
    allocate( URB_ATM_QV         (IA,JA)     )
    allocate( URB_ATM_PBL        (IA,JA)     )
    allocate( URB_ATM_SFC_PRES   (IA,JA)     )
    allocate( URB_ATM_SFLX_rad_dn(IA,JA,2,2) )
    allocate( URB_ATM_cosSZA     (IA,JA)     )
    allocate( URB_ATM_SFLX_rain  (IA,JA)     )
    allocate( URB_ATM_SFLX_snow  (IA,JA)     )
    URB_ATM_TEMP       (:,:)     = UNDEF
    URB_ATM_PRES       (:,:)     = UNDEF
    URB_ATM_W          (:,:)     = UNDEF
    URB_ATM_U          (:,:)     = UNDEF
    URB_ATM_V          (:,:)     = UNDEF
    URB_ATM_DENS       (:,:)     = UNDEF
    URB_ATM_QV         (:,:)     = UNDEF
    URB_ATM_PBL        (:,:)     = UNDEF
    URB_ATM_SFC_PRES   (:,:)     = UNDEF
    URB_ATM_SFLX_rad_dn(:,:,:,:) = UNDEF
    URB_ATM_cosSZA     (:,:)     = UNDEF
    URB_ATM_SFLX_rain  (:,:)     = UNDEF
    URB_ATM_SFLX_snow  (:,:)     = UNDEF

    ! counter intialize
    CNT_putATM_OCN = 0.0_RP
    CNT_putATM_LND = 0.0_RP
    CNT_putATM_URB = 0.0_RP
    CNT_putOCN     = 0.0_RP
    CNT_putLND     = 0.0_RP
    CNT_putURB     = 0.0_RP

    return
  end subroutine CPL_vars_setup

  !-----------------------------------------------------------------------------
  subroutine CPL_putATM( &
       TEMP,        &
       PRES,        &
       W,           &
       U,           &
       V,           &
       DENS,        &
       QTRC,        &
       PBL,         &
       SFC_PRES,    &
       SFLX_rad_dn, &
       cosSZA,      &
       SFLX_rain,   &
       SFLX_snow,   &
       countup      )
    implicit none

    ! arguments
    real(RP), intent(in) :: TEMP       (IA,JA)
    real(RP), intent(in) :: PRES       (IA,JA)
    real(RP), intent(in) :: W          (IA,JA)
    real(RP), intent(in) :: U          (IA,JA)
    real(RP), intent(in) :: V          (IA,JA)
    real(RP), intent(in) :: DENS       (IA,JA)
    real(RP), intent(in) :: QTRC       (IA,JA,QA)
    real(RP), intent(in) :: PBL        (IA,JA)
    real(RP), intent(in) :: SFC_PRES   (IA,JA)
    real(RP), intent(in) :: SFLX_rad_dn(IA,JA,2,2)
    real(RP), intent(in) :: cosSZA     (IA,JA)
    real(RP), intent(in) :: SFLX_rain  (IA,JA)
    real(RP), intent(in) :: SFLX_snow  (IA,JA)

    logical,  intent(in) :: countup

    ! works
    integer :: i, j
    !---------------------------------------------------------------------------

    do j = JS, JE
    do i = IS, IE
       ! for ocean
       OCN_ATM_TEMP       (i,j)     = OCN_ATM_TEMP       (i,j)     * CNT_putATM_OCN + TEMP       (i,j)
       OCN_ATM_PRES       (i,j)     = OCN_ATM_PRES       (i,j)     * CNT_putATM_OCN + PRES       (i,j)
       OCN_ATM_W          (i,j)     = OCN_ATM_W          (i,j)     * CNT_putATM_OCN + W          (i,j)
       OCN_ATM_U          (i,j)     = OCN_ATM_U          (i,j)     * CNT_putATM_OCN + U          (i,j)
       OCN_ATM_V          (i,j)     = OCN_ATM_V          (i,j)     * CNT_putATM_OCN + V          (i,j)
       OCN_ATM_DENS       (i,j)     = OCN_ATM_DENS       (i,j)     * CNT_putATM_OCN + DENS       (i,j)
       OCN_ATM_QV         (i,j)     = OCN_ATM_QV         (i,j)     * CNT_putATM_OCN + QTRC       (i,j,I_QV)
       OCN_ATM_PBL        (i,j)     = OCN_ATM_PBL        (i,j)     * CNT_putATM_OCN + PBL        (i,j)
       OCN_ATM_SFC_PRES   (i,j)     = OCN_ATM_SFC_PRES   (i,j)     * CNT_putATM_OCN + SFC_PRES   (i,j)
       OCN_ATM_SFLX_rad_dn(i,j,:,:) = OCN_ATM_SFLX_rad_dn(i,j,:,:) * CNT_putATM_OCN + SFLX_rad_dn(i,j,:,:)
       OCN_ATM_cosSZA     (i,j)     = OCN_ATM_cosSZA     (i,j)     * CNT_putATM_OCN + cosSZA     (i,j)
       OCN_ATM_SFLX_rain  (i,j)     = OCN_ATM_SFLX_rain  (i,j)     * CNT_putATM_OCN + SFLX_rain  (i,j)
       OCN_ATM_SFLX_snow  (i,j)     = OCN_ATM_SFLX_snow  (i,j)     * CNT_putATM_OCN + SFLX_snow  (i,j)
       ! for land
       LND_ATM_TEMP       (i,j)     = LND_ATM_TEMP       (i,j)     * CNT_putATM_LND + TEMP       (i,j)
       LND_ATM_PRES       (i,j)     = LND_ATM_PRES       (i,j)     * CNT_putATM_LND + PRES       (i,j)
       LND_ATM_W          (i,j)     = LND_ATM_W          (i,j)     * CNT_putATM_LND + W          (i,j)
       LND_ATM_U          (i,j)     = LND_ATM_U          (i,j)     * CNT_putATM_LND + U          (i,j)
       LND_ATM_V          (i,j)     = LND_ATM_V          (i,j)     * CNT_putATM_LND + V          (i,j)
       LND_ATM_DENS       (i,j)     = LND_ATM_DENS       (i,j)     * CNT_putATM_LND + DENS       (i,j)
       LND_ATM_QV         (i,j)     = LND_ATM_QV         (i,j)     * CNT_putATM_LND + QTRC       (i,j,I_QV)
       LND_ATM_PBL        (i,j)     = LND_ATM_PBL        (i,j)     * CNT_putATM_LND + PBL        (i,j)
       LND_ATM_SFC_PRES   (i,j)     = LND_ATM_SFC_PRES   (i,j)     * CNT_putATM_LND + SFC_PRES   (i,j)
       LND_ATM_SFLX_rad_dn(i,j,:,:) = LND_ATM_SFLX_rad_dn(i,j,:,:) * CNT_putATM_LND + SFLX_rad_dn(i,j,:,:)
       LND_ATM_cosSZA     (i,j)     = LND_ATM_cosSZA     (i,j)     * CNT_putATM_LND + cosSZA     (i,j)
       LND_ATM_SFLX_rain  (i,j)     = LND_ATM_SFLX_rain  (i,j)     * CNT_putATM_LND + SFLX_rain  (i,j)
       LND_ATM_SFLX_snow  (i,j)     = LND_ATM_SFLX_snow  (i,j)     * CNT_putATM_LND + SFLX_snow  (i,j)
       ! for urban
       URB_ATM_TEMP       (i,j)     = URB_ATM_TEMP       (i,j)     * CNT_putATM_URB + TEMP       (i,j)
       URB_ATM_PRES       (i,j)     = URB_ATM_PRES       (i,j)     * CNT_putATM_URB + PRES       (i,j)
       URB_ATM_W          (i,j)     = URB_ATM_W          (i,j)     * CNT_putATM_URB + W          (i,j)
       URB_ATM_U          (i,j)     = URB_ATM_U          (i,j)     * CNT_putATM_URB + U          (i,j)
       URB_ATM_V          (i,j)     = URB_ATM_V          (i,j)     * CNT_putATM_URB + V          (i,j)
       URB_ATM_DENS       (i,j)     = URB_ATM_DENS       (i,j)     * CNT_putATM_URB + DENS       (i,j)
       URB_ATM_QV         (i,j)     = URB_ATM_QV         (i,j)     * CNT_putATM_URB + QTRC       (i,j,I_QV)
       URB_ATM_PBL        (i,j)     = URB_ATM_PBL        (i,j)     * CNT_putATM_URB + PBL        (i,j)
       URB_ATM_SFC_PRES   (i,j)     = URB_ATM_SFC_PRES   (i,j)     * CNT_putATM_URB + SFC_PRES   (i,j)
       URB_ATM_SFLX_rad_dn(i,j,:,:) = URB_ATM_SFLX_rad_dn(i,j,:,:) * CNT_putATM_URB + SFLX_rad_dn(i,j,:,:)
       URB_ATM_cosSZA     (i,j)     = URB_ATM_cosSZA     (i,j)     * CNT_putATM_URB + cosSZA     (i,j)
       URB_ATM_SFLX_rain  (i,j)     = URB_ATM_SFLX_rain  (i,j)     * CNT_putATM_URB + SFLX_rain  (i,j)
       URB_ATM_SFLX_snow  (i,j)     = URB_ATM_SFLX_snow  (i,j)     * CNT_putATM_URB + SFLX_snow  (i,j)
    enddo
    enddo

    do j = JS, JE
    do i = IS, IE
       ! for ocean
       OCN_ATM_TEMP       (i,j)     = OCN_ATM_TEMP       (i,j)     / ( CNT_putATM_OCN + 1.0_RP )
       OCN_ATM_PRES       (i,j)     = OCN_ATM_PRES       (i,j)     / ( CNT_putATM_OCN + 1.0_RP )
       OCN_ATM_W          (i,j)     = OCN_ATM_W          (i,j)     / ( CNT_putATM_OCN + 1.0_RP )
       OCN_ATM_U          (i,j)     = OCN_ATM_U          (i,j)     / ( CNT_putATM_OCN + 1.0_RP )
       OCN_ATM_V          (i,j)     = OCN_ATM_V          (i,j)     / ( CNT_putATM_OCN + 1.0_RP )
       OCN_ATM_DENS       (i,j)     = OCN_ATM_DENS       (i,j)     / ( CNT_putATM_OCN + 1.0_RP )
       OCN_ATM_QV         (i,j)     = OCN_ATM_QV         (i,j)     / ( CNT_putATM_OCN + 1.0_RP )
       OCN_ATM_PBL        (i,j)     = OCN_ATM_PBL        (i,j)     / ( CNT_putATM_OCN + 1.0_RP )
       OCN_ATM_SFC_PRES   (i,j)     = OCN_ATM_SFC_PRES   (i,j)     / ( CNT_putATM_OCN + 1.0_RP )
       OCN_ATM_SFLX_rad_dn(i,j,:,:) = OCN_ATM_SFLX_rad_dn(i,j,:,:) / ( CNT_putATM_OCN + 1.0_RP )
       OCN_ATM_cosSZA     (i,j)     = OCN_ATM_cosSZA     (i,j)     / ( CNT_putATM_OCN + 1.0_RP )
       OCN_ATM_SFLX_rain  (i,j)     = OCN_ATM_SFLX_rain  (i,j)     / ( CNT_putATM_OCN + 1.0_RP )
       OCN_ATM_SFLX_snow  (i,j)     = OCN_ATM_SFLX_snow  (i,j)     / ( CNT_putATM_OCN + 1.0_RP )
       ! for ocean
       LND_ATM_TEMP       (i,j)     = LND_ATM_TEMP       (i,j)     / ( CNT_putATM_LND + 1.0_RP )
       LND_ATM_PRES       (i,j)     = LND_ATM_PRES       (i,j)     / ( CNT_putATM_LND + 1.0_RP )
       LND_ATM_W          (i,j)     = LND_ATM_W          (i,j)     / ( CNT_putATM_LND + 1.0_RP )
       LND_ATM_U          (i,j)     = LND_ATM_U          (i,j)     / ( CNT_putATM_LND + 1.0_RP )
       LND_ATM_V          (i,j)     = LND_ATM_V          (i,j)     / ( CNT_putATM_LND + 1.0_RP )
       LND_ATM_DENS       (i,j)     = LND_ATM_DENS       (i,j)     / ( CNT_putATM_LND + 1.0_RP )
       LND_ATM_QV         (i,j)     = LND_ATM_QV         (i,j)     / ( CNT_putATM_LND + 1.0_RP )
       LND_ATM_PBL        (i,j)     = LND_ATM_PBL        (i,j)     / ( CNT_putATM_LND + 1.0_RP )
       LND_ATM_SFC_PRES   (i,j)     = LND_ATM_SFC_PRES   (i,j)     / ( CNT_putATM_LND + 1.0_RP )
       LND_ATM_SFLX_rad_dn(i,j,:,:) = LND_ATM_SFLX_rad_dn(i,j,:,:) / ( CNT_putATM_LND + 1.0_RP )
       LND_ATM_cosSZA     (i,j)     = LND_ATM_cosSZA     (i,j)     / ( CNT_putATM_LND + 1.0_RP )
       LND_ATM_SFLX_rain  (i,j)     = LND_ATM_SFLX_rain  (i,j)     / ( CNT_putATM_LND + 1.0_RP )
       LND_ATM_SFLX_snow  (i,j)     = LND_ATM_SFLX_snow  (i,j)     / ( CNT_putATM_LND + 1.0_RP )
       ! for ocean
       URB_ATM_TEMP       (i,j)     = URB_ATM_TEMP       (i,j)     / ( CNT_putATM_URB + 1.0_RP )
       URB_ATM_PRES       (i,j)     = URB_ATM_PRES       (i,j)     / ( CNT_putATM_URB + 1.0_RP )
       URB_ATM_W          (i,j)     = URB_ATM_W          (i,j)     / ( CNT_putATM_URB + 1.0_RP )
       URB_ATM_U          (i,j)     = URB_ATM_U          (i,j)     / ( CNT_putATM_URB + 1.0_RP )
       URB_ATM_V          (i,j)     = URB_ATM_V          (i,j)     / ( CNT_putATM_URB + 1.0_RP )
       URB_ATM_DENS       (i,j)     = URB_ATM_DENS       (i,j)     / ( CNT_putATM_URB + 1.0_RP )
       URB_ATM_QV         (i,j)     = URB_ATM_QV         (i,j)     / ( CNT_putATM_URB + 1.0_RP )
       URB_ATM_PBL        (i,j)     = URB_ATM_PBL        (i,j)     / ( CNT_putATM_URB + 1.0_RP )
       URB_ATM_SFC_PRES   (i,j)     = URB_ATM_SFC_PRES   (i,j)     / ( CNT_putATM_URB + 1.0_RP )
       URB_ATM_SFLX_rad_dn(i,j,:,:) = URB_ATM_SFLX_rad_dn(i,j,:,:) / ( CNT_putATM_URB + 1.0_RP )
       URB_ATM_cosSZA     (i,j)     = URB_ATM_cosSZA     (i,j)     / ( CNT_putATM_URB + 1.0_RP )
       URB_ATM_SFLX_rain  (i,j)     = URB_ATM_SFLX_rain  (i,j)     / ( CNT_putATM_URB + 1.0_RP )
       URB_ATM_SFLX_snow  (i,j)     = URB_ATM_SFLX_snow  (i,j)     / ( CNT_putATM_URB + 1.0_RP )
    enddo
    enddo

    if( countup ) then
       CNT_putATM_OCN = CNT_putATM_OCN + 1.0_RP
       CNT_putATM_LND = CNT_putATM_LND + 1.0_RP
       CNT_putATM_URB = CNT_putATM_URB + 1.0_RP
    end if

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
       SFLX_WH,    &
       SFLX_evap,  &
       U10,        &
       V10,        &
       T2,         &
       Q2,         &
       countup     )
    implicit none

    ! arguments
    real(RP), intent(in) :: SFC_TEMP  (IA,JA)
    real(RP), intent(in) :: SFC_albedo(IA,JA,2)
    real(RP), intent(in) :: SFC_Z0M   (IA,JA)
    real(RP), intent(in) :: SFC_Z0H   (IA,JA)
    real(RP), intent(in) :: SFC_Z0E   (IA,JA)
    real(RP), intent(in) :: SFLX_MW   (IA,JA)
    real(RP), intent(in) :: SFLX_MU   (IA,JA)
    real(RP), intent(in) :: SFLX_MV   (IA,JA)
    real(RP), intent(in) :: SFLX_SH   (IA,JA)
    real(RP), intent(in) :: SFLX_LH   (IA,JA)
    real(RP), intent(in) :: SFLX_WH   (IA,JA)
    real(RP), intent(in) :: SFLX_evap (IA,JA)
    real(RP), intent(in) :: U10       (IA,JA)
    real(RP), intent(in) :: V10       (IA,JA)
    real(RP), intent(in) :: T2        (IA,JA)
    real(RP), intent(in) :: Q2        (IA,JA)

    logical,  intent(in) :: countup

    ! works
    integer :: i, j
    !---------------------------------------------------------------------------

    do j = JS, JE
    do i = IS, IE
       OCN_SFC_TEMP  (i,j)      = OCN_SFC_TEMP  (i,j)      * CNT_putOCN + SFC_TEMP  (i,j)
       OCN_SFC_albedo(i,j,I_LW) = OCN_SFC_albedo(i,j,I_LW) * CNT_putOCN + SFC_albedo(i,j,I_LW)
       OCN_SFC_albedo(i,j,I_SW) = OCN_SFC_albedo(i,j,I_SW) * CNT_putOCN + SFC_albedo(i,j,I_SW)
       OCN_SFC_Z0M   (i,j)      = OCN_SFC_Z0M   (i,j)      * CNT_putOCN + SFC_Z0M   (i,j)
       OCN_SFC_Z0H   (i,j)      = OCN_SFC_Z0H   (i,j)      * CNT_putOCN + SFC_Z0H   (i,j)
       OCN_SFC_Z0E   (i,j)      = OCN_SFC_Z0E   (i,j)      * CNT_putOCN + SFC_Z0E   (i,j)
       OCN_SFLX_MW   (i,j)      = OCN_SFLX_MW   (i,j)      * CNT_putOCN + SFLX_MW   (i,j)
       OCN_SFLX_MU   (i,j)      = OCN_SFLX_MU   (i,j)      * CNT_putOCN + SFLX_MU   (i,j)
       OCN_SFLX_MV   (i,j)      = OCN_SFLX_MV   (i,j)      * CNT_putOCN + SFLX_MV   (i,j)
       OCN_SFLX_SH   (i,j)      = OCN_SFLX_SH   (i,j)      * CNT_putOCN + SFLX_SH   (i,j)
       OCN_SFLX_LH   (i,j)      = OCN_SFLX_LH   (i,j)      * CNT_putOCN + SFLX_LH   (i,j)
       OCN_SFLX_WH   (i,j)      = OCN_SFLX_WH   (i,j)      * CNT_putOCN + SFLX_WH   (i,j)
       OCN_SFLX_evap (i,j)      = OCN_SFLX_evap (i,j)      * CNT_putOCN + SFLX_evap (i,j)
       OCN_U10       (i,j)      = OCN_U10       (i,j)      * CNT_putOCN + U10       (i,j)
       OCN_V10       (i,j)      = OCN_V10       (i,j)      * CNT_putOCN + V10       (i,j)
       OCN_T2        (i,j)      = OCN_T2        (i,j)      * CNT_putOCN + T2        (i,j)
       OCN_Q2        (i,j)      = OCN_Q2        (i,j)      * CNT_putOCN + Q2        (i,j)
    enddo
    enddo

    do j = JS, JE
    do i = IS, IE
       OCN_SFC_TEMP  (i,j)      = OCN_SFC_TEMP  (i,j)      / ( CNT_putOCN + 1.0_RP )
       OCN_SFC_albedo(i,j,I_LW) = OCN_SFC_albedo(i,j,I_LW) / ( CNT_putOCN + 1.0_RP )
       OCN_SFC_albedo(i,j,I_SW) = OCN_SFC_albedo(i,j,I_SW) / ( CNT_putOCN + 1.0_RP )
       OCN_SFC_Z0M   (i,j)      = OCN_SFC_Z0M   (i,j)      / ( CNT_putOCN + 1.0_RP )
       OCN_SFC_Z0H   (i,j)      = OCN_SFC_Z0H   (i,j)      / ( CNT_putOCN + 1.0_RP )
       OCN_SFC_Z0E   (i,j)      = OCN_SFC_Z0E   (i,j)      / ( CNT_putOCN + 1.0_RP )
       OCN_SFLX_MW   (i,j)      = OCN_SFLX_MW   (i,j)      / ( CNT_putOCN + 1.0_RP )
       OCN_SFLX_MU   (i,j)      = OCN_SFLX_MU   (i,j)      / ( CNT_putOCN + 1.0_RP )
       OCN_SFLX_MV   (i,j)      = OCN_SFLX_MV   (i,j)      / ( CNT_putOCN + 1.0_RP )
       OCN_SFLX_SH   (i,j)      = OCN_SFLX_SH   (i,j)      / ( CNT_putOCN + 1.0_RP )
       OCN_SFLX_LH   (i,j)      = OCN_SFLX_LH   (i,j)      / ( CNT_putOCN + 1.0_RP )
       OCN_SFLX_WH   (i,j)      = OCN_SFLX_WH   (i,j)      / ( CNT_putOCN + 1.0_RP )
       OCN_SFLX_evap (i,j)      = OCN_SFLX_evap (i,j)      / ( CNT_putOCN + 1.0_RP )
       OCN_U10       (i,j)      = OCN_U10       (i,j)      / ( CNT_putOCN + 1.0_RP )
       OCN_V10       (i,j)      = OCN_V10       (i,j)      / ( CNT_putOCN + 1.0_RP )
       OCN_T2        (i,j)      = OCN_T2        (i,j)      / ( CNT_putOCN + 1.0_RP )
       OCN_Q2        (i,j)      = OCN_Q2        (i,j)      / ( CNT_putOCN + 1.0_RP )
    enddo
    enddo

    if( countup ) then
       CNT_putOCN = CNT_putOCN + 1.0_RP
    end if

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
       SFLX_evap,  &
       U10,        &
       V10,        &
       T2,         &
       Q2,         &
       countup     )
    implicit none

    ! arguments
    real(RP), intent(in) :: SFC_TEMP  (IA,JA)
    real(RP), intent(in) :: SFC_albedo(IA,JA,2)
    real(RP), intent(in) :: SFC_Z0M   (IA,JA)
    real(RP), intent(in) :: SFC_Z0H   (IA,JA)
    real(RP), intent(in) :: SFC_Z0E   (IA,JA)
    real(RP), intent(in) :: SFLX_MW   (IA,JA)
    real(RP), intent(in) :: SFLX_MU   (IA,JA)
    real(RP), intent(in) :: SFLX_MV   (IA,JA)
    real(RP), intent(in) :: SFLX_SH   (IA,JA)
    real(RP), intent(in) :: SFLX_LH   (IA,JA)
    real(RP), intent(in) :: SFLX_GH   (IA,JA)
    real(RP), intent(in) :: SFLX_evap (IA,JA)
    real(RP), intent(in) :: U10       (IA,JA)
    real(RP), intent(in) :: V10       (IA,JA)
    real(RP), intent(in) :: T2        (IA,JA)
    real(RP), intent(in) :: Q2        (IA,JA)

    logical,  intent(in) :: countup

    ! works
    integer :: i, j
    !---------------------------------------------------------------------------

    do j = JS, JE
    do i = IS, IE
       LND_SFC_TEMP  (i,j)      = LND_SFC_TEMP  (i,j)      * CNT_putLND + SFC_TEMP  (i,j)
       LND_SFC_albedo(i,j,I_LW) = LND_SFC_albedo(i,j,I_LW) * CNT_putLND + SFC_albedo(i,j,I_LW)
       LND_SFC_albedo(i,j,I_SW) = LND_SFC_albedo(i,j,I_SW) * CNT_putLND + SFC_albedo(i,j,I_SW)
       LND_SFC_Z0M   (i,j)      = LND_SFC_Z0M   (i,j)      * CNT_putLND + SFC_Z0M   (i,j)
       LND_SFC_Z0H   (i,j)      = LND_SFC_Z0H   (i,j)      * CNT_putLND + SFC_Z0H   (i,j)
       LND_SFC_Z0E   (i,j)      = LND_SFC_Z0E   (i,j)      * CNT_putLND + SFC_Z0E   (i,j)
       LND_SFLX_MW   (i,j)      = LND_SFLX_MW   (i,j)      * CNT_putLND + SFLX_MW   (i,j)
       LND_SFLX_MU   (i,j)      = LND_SFLX_MU   (i,j)      * CNT_putLND + SFLX_MU   (i,j)
       LND_SFLX_MV   (i,j)      = LND_SFLX_MV   (i,j)      * CNT_putLND + SFLX_MV   (i,j)
       LND_SFLX_SH   (i,j)      = LND_SFLX_SH   (i,j)      * CNT_putLND + SFLX_SH   (i,j)
       LND_SFLX_LH   (i,j)      = LND_SFLX_LH   (i,j)      * CNT_putLND + SFLX_LH   (i,j)
       LND_SFLX_GH   (i,j)      = LND_SFLX_GH   (i,j)      * CNT_putLND + SFLX_GH   (i,j)
       LND_SFLX_evap (i,j)      = LND_SFLX_evap (i,j)      * CNT_putLND + SFLX_evap (i,j)
       LND_U10       (i,j)      = LND_U10       (i,j)      * CNT_putLND + U10       (i,j)
       LND_V10       (i,j)      = LND_V10       (i,j)      * CNT_putLND + V10       (i,j)
       LND_T2        (i,j)      = LND_T2        (i,j)      * CNT_putLND + T2        (i,j)
       LND_Q2        (i,j)      = LND_Q2        (i,j)      * CNT_putLND + Q2        (i,j)
    enddo
    enddo

    do j = JS, JE
    do i = IS, IE
       LND_SFC_TEMP  (i,j)      = LND_SFC_TEMP  (i,j)      / ( CNT_putLND + 1.0_RP )
       LND_SFC_albedo(i,j,I_LW) = LND_SFC_albedo(i,j,I_LW) / ( CNT_putLND + 1.0_RP )
       LND_SFC_albedo(i,j,I_SW) = LND_SFC_albedo(i,j,I_SW) / ( CNT_putLND + 1.0_RP )
       LND_SFC_Z0M   (i,j)      = LND_SFC_Z0M   (i,j)      / ( CNT_putLND + 1.0_RP )
       LND_SFC_Z0H   (i,j)      = LND_SFC_Z0H   (i,j)      / ( CNT_putLND + 1.0_RP )
       LND_SFC_Z0E   (i,j)      = LND_SFC_Z0E   (i,j)      / ( CNT_putLND + 1.0_RP )
       LND_SFLX_MW   (i,j)      = LND_SFLX_MW   (i,j)      / ( CNT_putLND + 1.0_RP )
       LND_SFLX_MU   (i,j)      = LND_SFLX_MU   (i,j)      / ( CNT_putLND + 1.0_RP )
       LND_SFLX_MV   (i,j)      = LND_SFLX_MV   (i,j)      / ( CNT_putLND + 1.0_RP )
       LND_SFLX_SH   (i,j)      = LND_SFLX_SH   (i,j)      / ( CNT_putLND + 1.0_RP )
       LND_SFLX_LH   (i,j)      = LND_SFLX_LH   (i,j)      / ( CNT_putLND + 1.0_RP )
       LND_SFLX_GH   (i,j)      = LND_SFLX_GH   (i,j)      / ( CNT_putLND + 1.0_RP )
       LND_SFLX_evap (i,j)      = LND_SFLX_evap (i,j)      / ( CNT_putLND + 1.0_RP )
       LND_U10       (i,j)      = LND_U10       (i,j)      / ( CNT_putLND + 1.0_RP )
       LND_V10       (i,j)      = LND_V10       (i,j)      / ( CNT_putLND + 1.0_RP )
       LND_T2        (i,j)      = LND_T2        (i,j)      / ( CNT_putLND + 1.0_RP )
       LND_Q2        (i,j)      = LND_Q2        (i,j)      / ( CNT_putLND + 1.0_RP )
    enddo
    enddo

    if( countup ) then
       CNT_putLND = CNT_putLND + 1.0_RP
    end if

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
       SFLX_GH,    &
       SFLX_evap,  &
       U10,        &
       V10,        &
       T2,         &
       Q2,         &
       countup     )
    implicit none

    ! arguments
    real(RP), intent(in) :: SFC_TEMP  (IA,JA)
    real(RP), intent(in) :: SFC_albedo(IA,JA,2)
    real(RP), intent(in) :: SFC_Z0M   (IA,JA)
    real(RP), intent(in) :: SFC_Z0H   (IA,JA)
    real(RP), intent(in) :: SFC_Z0E   (IA,JA)
    real(RP), intent(in) :: SFLX_MW   (IA,JA)
    real(RP), intent(in) :: SFLX_MU   (IA,JA)
    real(RP), intent(in) :: SFLX_MV   (IA,JA)
    real(RP), intent(in) :: SFLX_SH   (IA,JA)
    real(RP), intent(in) :: SFLX_LH   (IA,JA)
    real(RP), intent(in) :: SFLX_GH   (IA,JA)
    real(RP), intent(in) :: SFLX_evap (IA,JA)
    real(RP), intent(in) :: U10       (IA,JA)
    real(RP), intent(in) :: V10       (IA,JA)
    real(RP), intent(in) :: T2        (IA,JA)
    real(RP), intent(in) :: Q2        (IA,JA)

    logical,  intent(in) :: countup

    ! works
    integer :: i, j
    !---------------------------------------------------------------------------

    do j = JS, JE
    do i = IS, IE
       URB_SFC_TEMP  (i,j)      = URB_SFC_TEMP  (i,j)      * CNT_putURB + SFC_TEMP  (i,j)
       URB_SFC_albedo(i,j,I_LW) = URB_SFC_albedo(i,j,I_LW) * CNT_putURB + SFC_albedo(i,j,I_LW)
       URB_SFC_albedo(i,j,I_SW) = URB_SFC_albedo(i,j,I_SW) * CNT_putURB + SFC_albedo(i,j,I_SW)
       URB_SFC_Z0M   (i,j)      = URB_SFC_Z0M   (i,j)      * CNT_putURB + SFC_Z0M   (i,j)
       URB_SFC_Z0H   (i,j)      = URB_SFC_Z0H   (i,j)      * CNT_putURB + SFC_Z0H   (i,j)
       URB_SFC_Z0E   (i,j)      = URB_SFC_Z0E   (i,j)      * CNT_putURB + SFC_Z0E   (i,j)
       URB_SFLX_MW   (i,j)      = URB_SFLX_MW   (i,j)      * CNT_putURB + SFLX_MW   (i,j)
       URB_SFLX_MU   (i,j)      = URB_SFLX_MU   (i,j)      * CNT_putURB + SFLX_MU   (i,j)
       URB_SFLX_MV   (i,j)      = URB_SFLX_MV   (i,j)      * CNT_putURB + SFLX_MV   (i,j)
       URB_SFLX_SH   (i,j)      = URB_SFLX_SH   (i,j)      * CNT_putURB + SFLX_SH   (i,j)
       URB_SFLX_LH   (i,j)      = URB_SFLX_LH   (i,j)      * CNT_putURB + SFLX_LH   (i,j)
       URB_SFLX_GH   (i,j)      = URB_SFLX_GH   (i,j)      * CNT_putURB + SFLX_GH   (i,j)
       URB_SFLX_evap (i,j)      = URB_SFLX_evap (i,j)      * CNT_putURB + SFLX_evap (i,j)
       URB_U10       (i,j)      = URB_U10       (i,j)      * CNT_putURB + U10       (i,j)
       URB_V10       (i,j)      = URB_V10       (i,j)      * CNT_putURB + V10       (i,j)
       URB_T2        (i,j)      = URB_T2        (i,j)      * CNT_putURB + T2        (i,j)
       URB_Q2        (i,j)      = URB_Q2        (i,j)      * CNT_putURB + Q2        (i,j)
    enddo
    enddo

    do j = JS, JE
    do i = IS, IE
       URB_SFC_TEMP  (i,j)      = URB_SFC_TEMP  (i,j)      / ( CNT_putURB + 1.0_RP )
       URB_SFC_albedo(i,j,I_LW) = URB_SFC_albedo(i,j,I_LW) / ( CNT_putURB + 1.0_RP )
       URB_SFC_albedo(i,j,I_SW) = URB_SFC_albedo(i,j,I_SW) / ( CNT_putURB + 1.0_RP )
       URB_SFC_Z0M   (i,j)      = URB_SFC_Z0M   (i,j)      / ( CNT_putURB + 1.0_RP )
       URB_SFC_Z0H   (i,j)      = URB_SFC_Z0H   (i,j)      / ( CNT_putURB + 1.0_RP )
       URB_SFC_Z0E   (i,j)      = URB_SFC_Z0E   (i,j)      / ( CNT_putURB + 1.0_RP )
       URB_SFLX_MW   (i,j)      = URB_SFLX_MW   (i,j)      / ( CNT_putURB + 1.0_RP )
       URB_SFLX_MU   (i,j)      = URB_SFLX_MU   (i,j)      / ( CNT_putURB + 1.0_RP )
       URB_SFLX_MV   (i,j)      = URB_SFLX_MV   (i,j)      / ( CNT_putURB + 1.0_RP )
       URB_SFLX_SH   (i,j)      = URB_SFLX_SH   (i,j)      / ( CNT_putURB + 1.0_RP )
       URB_SFLX_LH   (i,j)      = URB_SFLX_LH   (i,j)      / ( CNT_putURB + 1.0_RP )
       URB_SFLX_GH   (i,j)      = URB_SFLX_GH   (i,j)      / ( CNT_putURB + 1.0_RP )
       URB_SFLX_evap (i,j)      = URB_SFLX_evap (i,j)      / ( CNT_putURB + 1.0_RP )
       URB_U10       (i,j)      = URB_U10       (i,j)      / ( CNT_putURB + 1.0_RP )
       URB_V10       (i,j)      = URB_V10       (i,j)      / ( CNT_putURB + 1.0_RP )
       URB_T2        (i,j)      = URB_T2        (i,j)      / ( CNT_putURB + 1.0_RP )
       URB_Q2        (i,j)      = URB_Q2        (i,j)      / ( CNT_putURB + 1.0_RP )
    enddo
    enddo

    if( countup ) then
       CNT_putURB = CNT_putURB + 1.0_RP
    end if

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
       SFLX_GH,    &
       SFLX_QTRC,  &
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
    real(RP), intent(out) :: SFC_albedo(IA,JA,2)
    real(RP), intent(out) :: SFC_Z0M   (IA,JA)
    real(RP), intent(out) :: SFC_Z0H   (IA,JA)
    real(RP), intent(out) :: SFC_Z0E   (IA,JA)
    real(RP), intent(out) :: SFLX_MW   (IA,JA)
    real(RP), intent(out) :: SFLX_MU   (IA,JA)
    real(RP), intent(out) :: SFLX_MV   (IA,JA)
    real(RP), intent(out) :: SFLX_SH   (IA,JA)
    real(RP), intent(out) :: SFLX_LH   (IA,JA)
    real(RP), intent(out) :: SFLX_GH   (IA,JA)
    real(RP), intent(out) :: SFLX_QTRC (IA,JA,QA)
    real(RP), intent(out) :: U10       (IA,JA)
    real(RP), intent(out) :: V10       (IA,JA)
    real(RP), intent(out) :: T2        (IA,JA)
    real(RP), intent(out) :: Q2        (IA,JA)

    integer :: i, j, iq
    !---------------------------------------------------------------------------

    do j = JS, JE
    do i = IS, IE
       do iq = 1, QA
          SFLX_QTRC(i,j,iq) = 0.0_RP ! tentative
       enddo

       SFC_TEMP  (i,j)      = fact_ocean(i,j) * OCN_SFC_TEMP  (i,j) &
                            + fact_land (i,j) * LND_SFC_TEMP  (i,j) &
                            + fact_urban(i,j) * URB_SFC_TEMP  (i,j)

       SFC_albedo(i,j,I_LW) = fact_ocean(i,j) * OCN_SFC_albedo(i,j,I_LW) &
                            + fact_land (i,j) * LND_SFC_albedo(i,j,I_LW) &
                            + fact_urban(i,j) * URB_SFC_albedo(i,j,I_LW)

       SFC_albedo(i,j,I_SW) = fact_ocean(i,j) * OCN_SFC_albedo(i,j,I_SW) &
                            + fact_land (i,j) * LND_SFC_albedo(i,j,I_SW) &
                            + fact_urban(i,j) * URB_SFC_albedo(i,j,I_SW)

       SFC_Z0M   (i,j)      = fact_ocean(i,j) * OCN_SFC_Z0M   (i,j) &
                            + fact_land (i,j) * LND_SFC_Z0M   (i,j) &
                            + fact_urban(i,j) * URB_SFC_Z0M   (i,j)

       SFC_Z0H   (i,j)      = fact_ocean(i,j) * OCN_SFC_Z0H   (i,j) &
                            + fact_land (i,j) * LND_SFC_Z0H   (i,j) &
                            + fact_urban(i,j) * URB_SFC_Z0H   (i,j)

       SFC_Z0E   (i,j)      = fact_ocean(i,j) * OCN_SFC_Z0E   (i,j) &
                            + fact_land (i,j) * LND_SFC_Z0E   (i,j) &
                            + fact_urban(i,j) * URB_SFC_Z0E   (i,j)

       SFLX_MW   (i,j)      = fact_ocean(i,j) * OCN_SFLX_MW   (i,j) &
                            + fact_land (i,j) * LND_SFLX_MW   (i,j) &
                            + fact_urban(i,j) * URB_SFLX_MW   (i,j)

       SFLX_MU   (i,j)      = fact_ocean(i,j) * OCN_SFLX_MU   (i,j) &
                            + fact_land (i,j) * LND_SFLX_MU   (i,j) &
                            + fact_urban(i,j) * URB_SFLX_MU   (i,j)

       SFLX_MV   (i,j)      = fact_ocean(i,j) * OCN_SFLX_MV   (i,j) &
                            + fact_land (i,j) * LND_SFLX_MV   (i,j) &
                            + fact_urban(i,j) * URB_SFLX_MV   (i,j)

       SFLX_SH   (i,j)      = fact_ocean(i,j) * OCN_SFLX_SH   (i,j) &
                            + fact_land (i,j) * LND_SFLX_SH   (i,j) &
                            + fact_urban(i,j) * URB_SFLX_SH   (i,j)

       SFLX_LH   (i,j)      = fact_ocean(i,j) * OCN_SFLX_LH   (i,j) &
                            + fact_land (i,j) * LND_SFLX_LH   (i,j) &
                            + fact_urban(i,j) * URB_SFLX_LH   (i,j)

       SFLX_GH   (i,j)      = fact_ocean(i,j) * OCN_SFLX_WH   (i,j) &
                            + fact_land (i,j) * LND_SFLX_GH   (i,j) &
                            + fact_urban(i,j) * URB_SFLX_GH   (i,j)

       SFLX_QTRC (i,j,I_QV) = fact_ocean(i,j) * OCN_SFLX_evap (i,j) &
                            + fact_land (i,j) * LND_SFLX_evap (i,j) &
                            + fact_urban(i,j) * URB_SFLX_evap (i,j)

       U10       (i,j)      = fact_ocean(i,j) * OCN_U10       (i,j) &
                            + fact_land (i,j) * LND_U10       (i,j) &
                            + fact_urban(i,j) * URB_U10       (i,j)

       V10       (i,j)      = fact_ocean(i,j) * OCN_V10       (i,j) &
                            + fact_land (i,j) * LND_V10       (i,j) &
                            + fact_urban(i,j) * URB_V10       (i,j)

       T2        (i,j)      = fact_ocean(i,j) * OCN_T2        (i,j) &
                            + fact_land (i,j) * LND_T2        (i,j) &
                            + fact_urban(i,j) * URB_T2        (i,j)

       Q2        (i,j)      = fact_ocean(i,j) * OCN_Q2        (i,j) &
                            + fact_land (i,j) * LND_Q2        (i,j) &
                            + fact_urban(i,j) * URB_Q2        (i,j)
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
       SFC_PRES,    &
       SFLX_rad_dn, &
       cosSZA,      &
       SFLX_rain,   &
       SFLX_snow    )
    implicit none

    real(RP), intent(out) :: TEMP       (IA,JA)
    real(RP), intent(out) :: PRES       (IA,JA)
    real(RP), intent(out) :: W          (IA,JA)
    real(RP), intent(out) :: U          (IA,JA)
    real(RP), intent(out) :: V          (IA,JA)
    real(RP), intent(out) :: DENS       (IA,JA)
    real(RP), intent(out) :: QV         (IA,JA)
    real(RP), intent(out) :: PBL        (IA,JA)
    real(RP), intent(out) :: SFC_PRES   (IA,JA)
    real(RP), intent(out) :: SFLX_rad_dn(IA,JA,2,2)
    real(RP), intent(out) :: cosSZA     (IA,JA)
    real(RP), intent(out) :: SFLX_rain  (IA,JA)
    real(RP), intent(out) :: SFLX_snow  (IA,JA)

    integer :: i, j
    !---------------------------------------------------------------------------

!OCL XFILL
    do j = JS, JE
    do i = IS, IE
       TEMP       (i,j)     = OCN_ATM_TEMP       (i,j)
       PRES       (i,j)     = OCN_ATM_PRES       (i,j)
       W          (i,j)     = OCN_ATM_W          (i,j)
       U          (i,j)     = OCN_ATM_U          (i,j)
       V          (i,j)     = OCN_ATM_V          (i,j)
       DENS       (i,j)     = OCN_ATM_DENS       (i,j)
       QV         (i,j)     = OCN_ATM_QV         (i,j)
       PBL        (i,j)     = OCN_ATM_PBL        (i,j)
       SFC_PRES   (i,j)     = OCN_ATM_SFC_PRES   (i,j)
       SFLX_rad_dn(i,j,:,:) = OCN_ATM_SFLX_rad_dn(i,j,:,:)
       cosSZA     (i,j)     = OCN_ATM_cosSZA     (i,j)
       SFLX_rain  (i,j)     = OCN_ATM_SFLX_rain  (i,j)
       SFLX_snow  (i,j)     = OCN_ATM_SFLX_snow  (i,j)
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
       SFC_PRES,    &
       SFLX_rad_dn, &
       cosSZA,      &
       SFLX_rain,   &
       SFLX_snow    )
    implicit none

    real(RP), intent(out) :: TEMP       (IA,JA)
    real(RP), intent(out) :: PRES       (IA,JA)
    real(RP), intent(out) :: W          (IA,JA)
    real(RP), intent(out) :: U          (IA,JA)
    real(RP), intent(out) :: V          (IA,JA)
    real(RP), intent(out) :: DENS       (IA,JA)
    real(RP), intent(out) :: QV         (IA,JA)
    real(RP), intent(out) :: PBL        (IA,JA)
    real(RP), intent(out) :: SFC_PRES   (IA,JA)
    real(RP), intent(out) :: SFLX_rad_dn(IA,JA,2,2)
    real(RP), intent(out) :: cosSZA     (IA,JA)
    real(RP), intent(out) :: SFLX_rain  (IA,JA)
    real(RP), intent(out) :: SFLX_snow  (IA,JA)

    integer :: i, j
    !---------------------------------------------------------------------------

!OCL XFILL
    do j = JS, JE
    do i = IS, IE
       TEMP       (i,j)     = LND_ATM_TEMP       (i,j)
       PRES       (i,j)     = LND_ATM_PRES       (i,j)
       W          (i,j)     = LND_ATM_W          (i,j)
       U          (i,j)     = LND_ATM_U          (i,j)
       V          (i,j)     = LND_ATM_V          (i,j)
       DENS       (i,j)     = LND_ATM_DENS       (i,j)
       QV         (i,j)     = LND_ATM_QV         (i,j)
       PBL        (i,j)     = LND_ATM_PBL        (i,j)
       SFC_PRES   (i,j)     = LND_ATM_SFC_PRES   (i,j)
       SFLX_rad_dn(i,j,:,:) = LND_ATM_SFLX_rad_dn(i,j,:,:)
       cosSZA     (i,j)     = LND_ATM_cosSZA     (i,j)
       SFLX_rain  (i,j)     = LND_ATM_SFLX_rain  (i,j)
       SFLX_snow  (i,j)     = LND_ATM_SFLX_snow  (i,j)
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
       SFC_PRES,    &
       SFLX_rad_dn, &
       cosSZA,      &
       SFLX_rain,   &
       SFLX_snow    )
    implicit none

    real(RP), intent(out) :: TEMP       (IA,JA)
    real(RP), intent(out) :: PRES       (IA,JA)
    real(RP), intent(out) :: W          (IA,JA)
    real(RP), intent(out) :: U          (IA,JA)
    real(RP), intent(out) :: V          (IA,JA)
    real(RP), intent(out) :: DENS       (IA,JA)
    real(RP), intent(out) :: QV         (IA,JA)
    real(RP), intent(out) :: PBL        (IA,JA)
    real(RP), intent(out) :: SFC_PRES   (IA,JA)
    real(RP), intent(out) :: SFLX_rad_dn(IA,JA,2,2)
    real(RP), intent(out) :: cosSZA     (IA,JA)
    real(RP), intent(out) :: SFLX_rain  (IA,JA)
    real(RP), intent(out) :: SFLX_snow  (IA,JA)

    integer :: i, j
    !---------------------------------------------------------------------------

!OCL XFILL
    do j = JS, JE
    do i = IS, IE
       TEMP       (i,j)     = URB_ATM_TEMP       (i,j)
       PRES       (i,j)     = URB_ATM_PRES       (i,j)
       W          (i,j)     = URB_ATM_W          (i,j)
       U          (i,j)     = URB_ATM_U          (i,j)
       V          (i,j)     = URB_ATM_V          (i,j)
       DENS       (i,j)     = URB_ATM_DENS       (i,j)
       QV         (i,j)     = URB_ATM_QV         (i,j)
       PBL        (i,j)     = URB_ATM_PBL        (i,j)
       SFC_PRES   (i,j)     = URB_ATM_SFC_PRES   (i,j)
       SFLX_rad_dn(i,j,:,:) = URB_ATM_SFLX_rad_dn(i,j,:,:)
       cosSZA     (i,j)     = URB_ATM_cosSZA     (i,j)
       SFLX_rain  (i,j)     = URB_ATM_SFLX_rain  (i,j)
       SFLX_snow  (i,j)     = URB_ATM_SFLX_snow  (i,j)
    enddo
    enddo

    CNT_putATM_URB = 0.0_RP

    return
  end subroutine CPL_getATM_URB

end module mod_cpl_vars
