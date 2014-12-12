!-------------------------------------------------------------------------------
!> module COUPLER Variables
!!
!! @par Description
!!          Container for coupler variables
!!
!! @author Team SCALE
!! @li      2013-08-31 (T.Yamaura)  [new]
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
  public :: CPL_vars_merge

  public :: CPL_putAtm_setup
  public :: CPL_putOcn_setup
  public :: CPL_putLnd_setup
  public :: CPL_putUrb_setup
  public :: CPL_putAtm_restart
  public :: CPL_putOcn_restart
  public :: CPL_putLnd_restart
  public :: CPL_putUrb_restart
  public :: CPL_putAtm
  public :: CPL_putOcn
  public :: CPL_putLnd
  public :: CPL_putUrb
  public :: CPL_getAtm
  public :: CPL_getAtm_SF
  public :: CPL_getOcn
  public :: CPL_getLnd
  public :: CPL_getUrb

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !

  !##### INPUT: Submodel->Coupler (flux is upward(Surface submodel->Atmosphere) positive) #####

  ! Input form atmosphere model
  real(RP), public, allocatable :: CPL_fromAtm_ATM_DENS  (:,:) ! density     at the lowermost atmosphere layer [kg/m3]
  real(RP), public, allocatable :: CPL_fromAtm_ATM_U     (:,:) ! velocity u  at the lowermost atmosphere layer [m/s]
  real(RP), public, allocatable :: CPL_fromAtm_ATM_V     (:,:) ! velocity v  at the lowermost atmosphere layer [m/s]
  real(RP), public, allocatable :: CPL_fromAtm_ATM_W     (:,:) ! velocity w  at the lowermost atmosphere layer [m/s]
  real(RP), public, allocatable :: CPL_fromAtm_ATM_TEMP  (:,:) ! temperature at the lowermost atmosphere layer [K]
  real(RP), public, allocatable :: CPL_fromAtm_ATM_PRES  (:,:) ! pressure    at the lowermost atmosphere layer [Pa]
  real(RP), public, allocatable :: CPL_fromAtm_ATM_QV    (:,:) ! water vapor at the lowermost atmosphere layer [kg/kg]
  real(RP), public, allocatable :: CPL_fromAtm_ATM_PBL   (:,:) ! the top of atmospheric mixing layer           [m]
  real(RP), public, allocatable :: CPL_fromAtm_SFC_PRES  (:,:) ! pressure    at the surface                    [Pa]
  real(RP), public, allocatable :: CPL_fromAtm_FLX_precip(:,:) ! liquid water                 flux [kg/m2/s]
  real(RP), public, allocatable :: CPL_fromAtm_FLX_SW_dn (:,:) ! downward shortwave radiation flux [J/m2/s]
  real(RP), public, allocatable :: CPL_fromAtm_FLX_LW_dn (:,:) ! downward longwave  radiation flux [J/m2/s]
  real(RP), public, allocatable :: CPL_fromAtm_ATM_Z1    (:,:) ! height of lowermost atmosphere layer (cell center) [m]

  ! Input form ocean model
  real(RP), public, allocatable :: CPL_fromOcn_SFC_TEMP  (:,:)   ! (first time only) surface skin temperature [K]
  real(RP), public, allocatable :: CPL_fromOcn_SFC_albedo(:,:,:) ! (first time only) surface albedo [0-1]
  real(RP), public, allocatable :: CPL_fromOcn_SFC_Z0M   (:,:)   ! roughness length for momemtum [m]
  real(RP), public, allocatable :: CPL_fromOcn_SFC_Z0H   (:,:)   ! roughness length for heat [m]
  real(RP), public, allocatable :: CPL_fromOcn_SFC_Z0E   (:,:)   ! roughness length for vapor [m]
  real(RP), public, allocatable :: CPL_fromOcn_OCN_TEMP  (:,:)   ! temperature at the uppermost ocean layer [K]

  ! Input form land model
  real(RP), public, allocatable :: CPL_fromLnd_SFC_TEMP  (:,:)   ! (first time only) surface skin temperature [K]
  real(RP), public, allocatable :: CPL_fromLnd_SFC_albedo(:,:,:) ! (first time only) surface albedo           [0-1]
  real(RP), public, allocatable :: CPL_fromLnd_LND_TCS   (:,:)   ! (first time only) thermal conductivity for soil [W/m/K]
  real(RP), public, allocatable :: CPL_fromLnd_LND_DZ    (:,:)   ! (first time only) soil depth [m]
  real(RP), public, allocatable :: CPL_fromLnd_SFC_Z0M   (:,:)   ! (first time only) roughness length for momemtum [m]
  real(RP), public, allocatable :: CPL_fromLnd_SFC_Z0H   (:,:)   ! (first time only) roughness length for heat     [m]
  real(RP), public, allocatable :: CPL_fromLnd_SFC_Z0E   (:,:)   ! (first time only) roughness length for vapor    [m]
  real(RP), public, allocatable :: CPL_fromLnd_LND_TEMP  (:,:)   ! temperature at the uppermost land layer [K]
  real(RP), public, allocatable :: CPL_fromLnd_LND_BETA  (:,:)   ! efficiency of evaporation [0-1]

  ! Input form urban model
  real(RP), public, allocatable :: CPL_fromUrb_SFC_Z0M   (:,:)   ! (first time only) roughness length for momemtum [m]
  real(RP), public, allocatable :: CPL_fromUrb_SFC_Z0H   (:,:)   ! (first time only) roughness length for momemtum [m]
  real(RP), public, allocatable :: CPL_fromUrb_SFC_Z0E   (:,:)   ! (first time only) roughness length for momemtum [m]
  real(RP), public, allocatable :: CPL_fromUrb_SFC_TEMP  (:,:)   ! (first time only) surface skin temperature [K]
  real(RP), public, allocatable :: CPL_fromUrb_SFC_albedo(:,:,:) ! (first time only) surface albedo           [0-1]

  !##### OUTPUT: Coupler->Submodel (flux is upward(Surface submodel->Atmosphere) positive) #####

  ! Output for atmosphere model (merged)
  real(RP), public, allocatable :: CPL_Merged_SFC_TEMP  (:,:)   ! Merged surface skin temperature [K]
  real(RP), public, allocatable :: CPL_Merged_SFC_albedo(:,:,:) ! Merged surface albedo           [0-1]
  real(RP), public, allocatable :: CPL_Merged_FLX_MU    (:,:)   ! Merged w-momentum    flux [kg/m2/s]
  real(RP), public, allocatable :: CPL_Merged_FLX_MV    (:,:)   ! Merged u-momentum    flux [kg/m2/s]
  real(RP), public, allocatable :: CPL_Merged_FLX_MW    (:,:)   ! Merged v-momentum    flux [kg/m2/s]
  real(RP), public, allocatable :: CPL_Merged_FLX_SH    (:,:)   ! Merged sensible heat flux [J/m2/s]
  real(RP), public, allocatable :: CPL_Merged_FLX_LH    (:,:)   ! Merged latent heat   flux [J/m2/s]
  real(RP), public, allocatable :: CPL_Merged_FLX_QV    (:,:)   ! Merged water vapor   flux [kg/m2/s]
  real(RP), public, allocatable :: CPL_Merged_Z0M       (:,:)   ! Merged roughness length   [m]
  real(RP), public, allocatable :: CPL_Merged_Z0H       (:,:)   ! Merged roughness length   [m]
  real(RP), public, allocatable :: CPL_Merged_Z0E       (:,:)   ! Merged roughness length   [m]
  real(RP), public, allocatable :: CPL_Merged_U10       (:,:)   ! Merged velocity u at 10m  [m/s]
  real(RP), public, allocatable :: CPL_Merged_V10       (:,:)   ! Merged velocity v at 10m  [m/s]
  real(RP), public, allocatable :: CPL_Merged_T2        (:,:)   ! Merged temperature at 2m  [K]
  real(RP), public, allocatable :: CPL_Merged_Q2        (:,:)   ! Merged water vapor at 2m  [kg/kg]
  ! Output for surface model (merged, for history)
  real(RP), public, allocatable :: CPL_Merged_FLX_heat  (:,:)   ! Merged heat flux [J/m2/s]

  ! Atmosphere-Ocean coupler: Output for atmosphere model
  real(RP), public, allocatable :: CPL_AtmOcn_ATM_FLX_MW    (:,:) ! w-momentum    flux [kg/m2/s]
  real(RP), public, allocatable :: CPL_AtmOcn_ATM_FLX_MU    (:,:) ! u-momentum    flux [kg/m2/s]
  real(RP), public, allocatable :: CPL_AtmOcn_ATM_FLX_MV    (:,:) ! v-momentum    flux [kg/m2/s]
  real(RP), public, allocatable :: CPL_AtmOcn_ATM_FLX_SH    (:,:) ! sensible heat flux [J/m2/s]
  real(RP), public, allocatable :: CPL_AtmOcn_ATM_FLX_LH    (:,:) ! latent heat   flux [J/m2/s]
  real(RP), public, allocatable :: CPL_AtmOcn_ATM_FLX_evap  (:,:) ! water vapor   flux [kg/m2/s]
  real(RP), public, allocatable :: CPL_AtmOcn_ATM_U10       (:,:) ! velocity u at 10m  [m/s]
  real(RP), public, allocatable :: CPL_AtmOcn_ATM_V10       (:,:) ! velocity v at 10m  [m/s]
  real(RP), public, allocatable :: CPL_AtmOcn_ATM_T2        (:,:) ! temperature at 2m  [K]
  real(RP), public, allocatable :: CPL_AtmOcn_ATM_Q2        (:,:) ! water vapor at 2m  [kg/kg]
  !  Atmosphere-Ocean coupler: Output for ocean model
  real(RP), public, allocatable :: CPL_AtmOcn_OCN_FLX_heat  (:,:) ! heat         flux  [J/m2/s]
  real(RP), public, allocatable :: CPL_AtmOcn_OCN_FLX_precip(:,:) ! liquid water flux  [kg/m2/s]
  real(RP), public, allocatable :: CPL_AtmOcn_OCN_FLX_evap  (:,:) ! water vapor  flux  [kg/m2/s]

  ! Atmosphere-Land coupler: Output for atmosphere model
  real(RP), public, allocatable :: CPL_AtmLnd_ATM_FLX_MW    (:,:) ! w-momentum    flux [kg/m2/s]
  real(RP), public, allocatable :: CPL_AtmLnd_ATM_FLX_MU    (:,:) ! u-momentum    flux [kg/m2/s]
  real(RP), public, allocatable :: CPL_AtmLnd_ATM_FLX_MV    (:,:) ! v-momentum    flux [kg/m2/s]
  real(RP), public, allocatable :: CPL_AtmLnd_ATM_FLX_SH    (:,:) ! sensible heat flux [J/m2/s]
  real(RP), public, allocatable :: CPL_AtmLnd_ATM_FLX_LH    (:,:) ! latent heat   flux [J/m2/s]
  real(RP), public, allocatable :: CPL_AtmLnd_ATM_FLX_evap  (:,:) ! water vapor   flux [kg/m2/s]
  real(RP), public, allocatable :: CPL_AtmLnd_ATM_U10       (:,:) ! velocity u at 10m  [m/s]
  real(RP), public, allocatable :: CPL_AtmLnd_ATM_V10       (:,:) ! velocity v at 10m  [m/s]
  real(RP), public, allocatable :: CPL_AtmLnd_ATM_T2        (:,:) ! temperature at 2m  [K]
  real(RP), public, allocatable :: CPL_AtmLnd_ATM_Q2        (:,:) ! water vapor at 2m  [kg/kg]
  ! Atmosphere-Land coupler: Output for land model
  real(RP), public, allocatable :: CPL_AtmLnd_LND_FLX_heat  (:,:) ! heat         flux  [J/m2/s]
  real(RP), public, allocatable :: CPL_AtmLnd_LND_FLX_precip(:,:) ! liquid water flux  [kg/m2/s]
  real(RP), public, allocatable :: CPL_AtmLnd_LND_FLX_evap  (:,:) ! water vapor  flux  [kg/m2/s]

  ! Atmosphere-Urban coupler: Output for atmosphere model
  real(RP), public, allocatable :: CPL_AtmUrb_ATM_FLX_MW    (:,:) ! w-momentum    flux [kg/m2/s]
  real(RP), public, allocatable :: CPL_AtmUrb_ATM_FLX_MU    (:,:) ! u-momentum    flux [kg/m2/s]
  real(RP), public, allocatable :: CPL_AtmUrb_ATM_FLX_MV    (:,:) ! v-momentum    flux [kg/m2/s]
  real(RP), public, allocatable :: CPL_AtmUrb_ATM_FLX_SH    (:,:) ! sensible heat flux [J/m2/s]
  real(RP), public, allocatable :: CPL_AtmUrb_ATM_FLX_LH    (:,:) ! latent heat   flux [J/m2/s]
  real(RP), public, allocatable :: CPL_AtmUrb_ATM_FLX_evap  (:,:) ! water vapor   flux [kg/m2/s]
  real(RP), public, allocatable :: CPL_AtmUrb_ATM_U10       (:,:) ! velocity u at 10m  [m/s]
  real(RP), public, allocatable :: CPL_AtmUrb_ATM_V10       (:,:) ! velocity v at 10m  [m/s]
  real(RP), public, allocatable :: CPL_AtmUrb_ATM_T2        (:,:) ! temperature at 2m  [K]
  real(RP), public, allocatable :: CPL_AtmUrb_ATM_Q2        (:,:) ! water vapor at 2m  [kg/kg]
  ! Atmosphere-Urban coupler: Output for urban model
  real(RP), public, allocatable :: CPL_AtmUrb_URB_FLX_heat  (:,:) ! heat         flux  [J/m2/s]
  real(RP), public, allocatable :: CPL_AtmUrb_URB_FLX_precip(:,:) ! liquid water flux  [kg/m2/s]
  real(RP), public, allocatable :: CPL_AtmUrb_URB_FLX_evap  (:,:) ! water vapor  flux  [kg/m2/s]

  ! counter
  real(RP), public :: CNT_putAtm    ! put counter for atmos flux
  real(RP), public :: CNT_putLnd    ! put counter for land flux
  real(RP), public :: CNT_putUrb    ! put counter for urban flux
  real(RP), public :: CNT_putOcn    ! put counter for ocean flux
  real(RP), public :: CNT_getAtmLnd ! get counter for atmos flux by land
  real(RP), public :: CNT_getAtmUrb ! get counter for atmos flux by urban
  real(RP), public :: CNT_getAtmOcn ! get counter for atmos flux by ocean
  real(RP), public :: CNT_getLnd    ! get counter for land flux
  real(RP), public :: CNT_getUrb    ! get counter for urban flux
  real(RP), public :: CNT_getOcn    ! get counter for ocean flux

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
    if( IO_L ) write(IO_FID_LOG,*) '++++++ Module[VARS] / Categ[CPL] / Origin[SCALE-LES]'

    allocate( CPL_fromAtm_ATM_DENS  (IA,JA) )
    allocate( CPL_fromAtm_ATM_U     (IA,JA) )
    allocate( CPL_fromAtm_ATM_V     (IA,JA) )
    allocate( CPL_fromAtm_ATM_W     (IA,JA) )
    allocate( CPL_fromAtm_ATM_TEMP  (IA,JA) )
    allocate( CPL_fromAtm_ATM_PRES  (IA,JA) )
    allocate( CPL_fromAtm_ATM_QV    (IA,JA) )
    allocate( CPL_fromAtm_ATM_PBL   (IA,JA) )
    allocate( CPL_fromAtm_SFC_PRES  (IA,JA) )
    allocate( CPL_fromAtm_FLX_precip(IA,JA) )
    allocate( CPL_fromAtm_FLX_SW_dn (IA,JA) )
    allocate( CPL_fromAtm_FLX_LW_dn (IA,JA) )
    CPL_fromAtm_ATM_DENS  (:,:) = UNDEF
    CPL_fromAtm_ATM_U     (:,:) = UNDEF
    CPL_fromAtm_ATM_V     (:,:) = UNDEF
    CPL_fromAtm_ATM_W     (:,:) = UNDEF
    CPL_fromAtm_ATM_TEMP  (:,:) = UNDEF
    CPL_fromAtm_ATM_PRES  (:,:) = UNDEF
    CPL_fromAtm_ATM_QV    (:,:) = UNDEF
    CPL_fromAtm_ATM_PBL   (:,:) = UNDEF
    CPL_fromAtm_SFC_PRES  (:,:) = UNDEF
    CPL_fromAtm_FLX_precip(:,:) = UNDEF
    CPL_fromAtm_FLX_SW_dn (:,:) = UNDEF
    CPL_fromAtm_FLX_LW_dn (:,:) = UNDEF

    allocate( CPL_fromOcn_SFC_TEMP  (IA,JA) )
    allocate( CPL_fromOcn_SFC_albedo(IA,JA,2) )
    allocate( CPL_fromOcn_SFC_Z0M   (IA,JA) )
    allocate( CPL_fromOcn_SFC_Z0H   (IA,JA) )
    allocate( CPL_fromOcn_SFC_Z0E   (IA,JA) )
    allocate( CPL_fromOcn_OCN_TEMP  (IA,JA) )
    CPL_fromOcn_SFC_TEMP  (:,:)   = UNDEF
    CPL_fromOcn_SFC_albedo(:,:,:) = UNDEF
    CPL_fromOcn_SFC_Z0M   (:,:)   = UNDEF
    CPL_fromOcn_SFC_Z0H   (:,:)   = UNDEF
    CPL_fromOcn_SFC_Z0E   (:,:)   = UNDEF
    CPL_fromOcn_OCN_TEMP  (:,:)   = UNDEF

    allocate( CPL_fromLnd_SFC_TEMP  (IA,JA) )
    allocate( CPL_fromLnd_SFC_albedo(IA,JA,2) )
    allocate( CPL_fromLnd_LND_TCS   (IA,JA) )
    allocate( CPL_fromLnd_LND_DZ    (IA,JA) )
    allocate( CPL_fromLnd_SFC_Z0M   (IA,JA) )
    allocate( CPL_fromLnd_SFC_Z0H   (IA,JA) )
    allocate( CPL_fromLnd_SFC_Z0E   (IA,JA) )
    allocate( CPL_fromLnd_LND_TEMP  (IA,JA) )
    allocate( CPL_fromLnd_LND_BETA  (IA,JA) )
    CPL_fromLnd_SFC_TEMP  (:,:)   = UNDEF
    CPL_fromLnd_SFC_albedo(:,:,:) = UNDEF
    CPL_fromLnd_LND_TCS   (:,:)   = UNDEF
    CPL_fromLnd_LND_DZ    (:,:)   = UNDEF
    CPL_fromLnd_SFC_Z0M   (:,:)   = UNDEF
    CPL_fromLnd_SFC_Z0H   (:,:)   = UNDEF
    CPL_fromLnd_SFC_Z0E   (:,:)   = UNDEF
    CPL_fromLnd_LND_TEMP  (:,:)   = UNDEF
    CPL_fromLnd_LND_BETA  (:,:)   = UNDEF

    allocate( CPL_fromUrb_SFC_Z0M   (IA,JA) )
    allocate( CPL_fromUrb_SFC_Z0H   (IA,JA) )
    allocate( CPL_fromUrb_SFC_Z0E   (IA,JA) )
    allocate( CPL_fromUrb_SFC_TEMP  (IA,JA) )
    allocate( CPL_fromUrb_SFC_albedo(IA,JA,2) )
    CPL_fromUrb_SFC_Z0M   (:,:)   = UNDEF
    CPL_fromUrb_SFC_Z0H   (:,:)   = UNDEF
    CPL_fromUrb_SFC_Z0E   (:,:)   = UNDEF
    CPL_fromUrb_SFC_TEMP  (:,:)   = UNDEF
    CPL_fromUrb_SFC_albedo(:,:,:) = UNDEF

    allocate( CPL_Merged_SFC_TEMP  (IA,JA) )
    allocate( CPL_Merged_SFC_albedo(IA,JA,2) )
    allocate( CPL_Merged_FLX_MU    (IA,JA) )
    allocate( CPL_Merged_FLX_MV    (IA,JA) )
    allocate( CPL_Merged_FLX_MW    (IA,JA) )
    allocate( CPL_Merged_FLX_SH    (IA,JA) )
    allocate( CPL_Merged_FLX_LH    (IA,JA) )
    allocate( CPL_Merged_FLX_QV    (IA,JA) )
    allocate( CPL_Merged_Z0M       (IA,JA) )
    allocate( CPL_Merged_Z0H       (IA,JA) )
    allocate( CPL_Merged_Z0E       (IA,JA) )
    allocate( CPL_Merged_U10       (IA,JA) )
    allocate( CPL_Merged_V10       (IA,JA) )
    allocate( CPL_Merged_T2        (IA,JA) )
    allocate( CPL_Merged_Q2        (IA,JA) )
    allocate( CPL_Merged_FLX_heat  (IA,JA) )
    CPL_Merged_SFC_TEMP  (:,:)   = UNDEF
    CPL_Merged_SFC_albedo(:,:,:) = UNDEF
    CPL_Merged_FLX_MU    (:,:)   = UNDEF
    CPL_Merged_FLX_MV    (:,:)   = UNDEF
    CPL_Merged_FLX_MW    (:,:)   = UNDEF
    CPL_Merged_FLX_SH    (:,:)   = UNDEF
    CPL_Merged_FLX_LH    (:,:)   = UNDEF
    CPL_Merged_FLX_QV    (:,:)   = UNDEF
    CPL_Merged_Z0M       (:,:)   = UNDEF
    CPL_Merged_Z0H       (:,:)   = UNDEF
    CPL_Merged_Z0E       (:,:)   = UNDEF
    CPL_Merged_U10       (:,:)   = UNDEF
    CPL_Merged_V10       (:,:)   = UNDEF
    CPL_Merged_T2        (:,:)   = UNDEF
    CPL_Merged_Q2        (:,:)   = UNDEF
    CPL_Merged_FLX_heat  (:,:)   = UNDEF

    allocate( CPL_AtmOcn_ATM_FLX_MU    (IA,JA) )
    allocate( CPL_AtmOcn_ATM_FLX_MV    (IA,JA) )
    allocate( CPL_AtmOcn_ATM_FLX_MW    (IA,JA) )
    allocate( CPL_AtmOcn_ATM_FLX_SH    (IA,JA) )
    allocate( CPL_AtmOcn_ATM_FLX_LH    (IA,JA) )
    allocate( CPL_AtmOcn_ATM_FLX_evap  (IA,JA) )
    allocate( CPL_AtmOcn_ATM_U10       (IA,JA) )
    allocate( CPL_AtmOcn_ATM_V10       (IA,JA) )
    allocate( CPL_AtmOcn_ATM_T2        (IA,JA) )
    allocate( CPL_AtmOcn_ATM_Q2        (IA,JA) )
    allocate( CPL_AtmOcn_OCN_FLX_heat  (IA,JA) )
    allocate( CPL_AtmOcn_OCN_FLX_precip(IA,JA) )
    allocate( CPL_AtmOcn_OCN_FLX_evap  (IA,JA) )
    CPL_AtmOcn_ATM_FLX_MU    (:,:) = UNDEF
    CPL_AtmOcn_ATM_FLX_MV    (:,:) = UNDEF
    CPL_AtmOcn_ATM_FLX_MW    (:,:) = UNDEF
    CPL_AtmOcn_ATM_FLX_SH    (:,:) = UNDEF
    CPL_AtmOcn_ATM_FLX_LH    (:,:) = UNDEF
    CPL_AtmOcn_ATM_FLX_evap  (:,:) = UNDEF
    CPL_AtmOcn_ATM_U10       (:,:) = UNDEF
    CPL_AtmOcn_ATM_V10       (:,:) = UNDEF
    CPL_AtmOcn_ATM_T2        (:,:) = UNDEF
    CPL_AtmOcn_ATM_Q2        (:,:) = UNDEF
    CPL_AtmOcn_OCN_FLX_heat  (:,:) = UNDEF
    CPL_AtmOcn_OCN_FLX_precip(:,:) = UNDEF
    CPL_AtmOcn_OCN_FLX_evap  (:,:) = UNDEF

    allocate( CPL_AtmLnd_ATM_FLX_MU    (IA,JA) )
    allocate( CPL_AtmLnd_ATM_FLX_MV    (IA,JA) )
    allocate( CPL_AtmLnd_ATM_FLX_MW    (IA,JA) )
    allocate( CPL_AtmLnd_ATM_FLX_SH    (IA,JA) )
    allocate( CPL_AtmLnd_ATM_FLX_LH    (IA,JA) )
    allocate( CPL_AtmLnd_ATM_FLX_evap  (IA,JA) )
    allocate( CPL_AtmLnd_ATM_U10       (IA,JA) )
    allocate( CPL_AtmLnd_ATM_V10       (IA,JA) )
    allocate( CPL_AtmLnd_ATM_T2        (IA,JA) )
    allocate( CPL_AtmLnd_ATM_Q2        (IA,JA) )
    allocate( CPL_AtmLnd_LND_FLX_heat  (IA,JA) )
    allocate( CPL_AtmLnd_LND_FLX_precip(IA,JA) )
    allocate( CPL_AtmLnd_LND_FLX_evap  (IA,JA) )
    CPL_AtmLnd_ATM_FLX_MU    (:,:) = UNDEF
    CPL_AtmLnd_ATM_FLX_MV    (:,:) = UNDEF
    CPL_AtmLnd_ATM_FLX_MW    (:,:) = UNDEF
    CPL_AtmLnd_ATM_FLX_SH    (:,:) = UNDEF
    CPL_AtmLnd_ATM_FLX_LH    (:,:) = UNDEF
    CPL_AtmLnd_ATM_FLX_evap  (:,:) = UNDEF
    CPL_AtmLnd_ATM_U10       (:,:) = UNDEF
    CPL_AtmLnd_ATM_V10       (:,:) = UNDEF
    CPL_AtmLnd_ATM_T2        (:,:) = UNDEF
    CPL_AtmLnd_ATM_Q2        (:,:) = UNDEF
    CPL_AtmLnd_LND_FLX_heat  (:,:) = UNDEF
    CPL_AtmLnd_LND_FLX_precip(:,:) = UNDEF
    CPL_AtmLnd_LND_FLX_evap  (:,:) = UNDEF

    allocate( CPL_AtmUrb_ATM_FLX_MU    (IA,JA) )
    allocate( CPL_AtmUrb_ATM_FLX_MV    (IA,JA) )
    allocate( CPL_AtmUrb_ATM_FLX_MW    (IA,JA) )
    allocate( CPL_AtmUrb_ATM_FLX_SH    (IA,JA) )
    allocate( CPL_AtmUrb_ATM_FLX_LH    (IA,JA) )
    allocate( CPL_AtmUrb_ATM_FLX_evap  (IA,JA) )
    allocate( CPL_AtmUrb_ATM_U10       (IA,JA) )
    allocate( CPL_AtmUrb_ATM_V10       (IA,JA) )
    allocate( CPL_AtmUrb_ATM_T2        (IA,JA) )
    allocate( CPL_AtmUrb_ATM_Q2        (IA,JA) )
    allocate( CPL_AtmUrb_URB_FLX_heat  (IA,JA) )
    allocate( CPL_AtmUrb_URB_FLX_precip(IA,JA) )
    allocate( CPL_AtmUrb_URB_FLX_evap  (IA,JA) )
    CPL_AtmUrb_ATM_FLX_MU    (:,:) = UNDEF
    CPL_AtmUrb_ATM_FLX_MV    (:,:) = UNDEF
    CPL_AtmUrb_ATM_FLX_MW    (:,:) = UNDEF
    CPL_AtmUrb_ATM_FLX_SH    (:,:) = UNDEF
    CPL_AtmUrb_ATM_FLX_LH    (:,:) = UNDEF
    CPL_AtmUrb_ATM_FLX_evap  (:,:) = UNDEF
    CPL_AtmUrb_ATM_U10       (:,:) = UNDEF
    CPL_AtmUrb_ATM_V10       (:,:) = UNDEF
    CPL_AtmUrb_ATM_T2        (:,:) = UNDEF
    CPL_AtmUrb_ATM_Q2        (:,:) = UNDEF
    CPL_AtmUrb_URB_FLX_heat  (:,:) = UNDEF
    CPL_AtmUrb_URB_FLX_precip(:,:) = UNDEF
    CPL_AtmUrb_URB_FLX_evap  (:,:) = UNDEF

    CNT_putAtm    = 0.0_RP
    CNT_putOcn    = 0.0_RP
    CNT_putLnd    = 0.0_RP
    CNT_putUrb    = 0.0_RP
    CNT_getAtmOcn = 0.0_RP
    CNT_getAtmLnd = 0.0_RP
    CNT_getAtmUrb = 0.0_RP
    CNT_getOcn    = 0.0_RP
    CNT_getLnd    = 0.0_RP
    CNT_getUrb    = 0.0_RP

    return
  end subroutine CPL_vars_setup

  !-----------------------------------------------------------------------------
  subroutine CPL_vars_merge
    use scale_landuse, only: &
       fact_ocean => LANDUSE_fact_ocean, &
       fact_land  => LANDUSE_fact_land,  &
       fact_urban => LANDUSE_fact_urban
    use scale_statistics, only: &
       STATISTICS_checktotal, &
       STAT_total
    implicit none

    real(RP) :: total ! dummy

    integer :: i, j
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*) '*** Coupler: Merge values of surface submodel'

    do j = JS, JE
    do i = IS, IE
       CPL_Merged_SFC_TEMP(i,j)     = fact_ocean(i,j) * CPL_fromOcn_SFC_TEMP  (i,j) &
                                    + fact_land (i,j) * CPL_fromLnd_SFC_TEMP  (i,j) &
                                    + fact_urban(i,j) * CPL_fromUrb_SFC_TEMP  (i,j)
    end do
    end do

    do j = JS, JE
    do i = IS, IE
       CPL_Merged_SFC_albedo(i,j,I_LW) = fact_ocean(i,j) * CPL_fromOcn_SFC_albedo(i,j,I_LW) &
                                       + fact_land (i,j) * CPL_fromLnd_SFC_albedo(i,j,I_LW) &
                                       + fact_urban(i,j) * CPL_fromUrb_SFC_albedo(i,j,I_LW)
    end do
    end do

    do j = JS, JE
    do i = IS, IE
       CPL_Merged_SFC_albedo(i,j,I_SW) = fact_ocean(i,j) * CPL_fromOcn_SFC_albedo(i,j,I_SW) &
                                       + fact_land (i,j) * CPL_fromLnd_SFC_albedo(i,j,I_SW) &
                                       + fact_urban(i,j) * CPL_fromUrb_SFC_albedo(i,j,I_SW)
    end do
    end do

    do j = JS, JE
    do i = IS, IE
       CPL_Merged_FLX_MW(i,j) = fact_ocean(i,j) * CPL_AtmOcn_ATM_FLX_MW  (i,j) &
                              + fact_land (i,j) * CPL_AtmLnd_ATM_FLX_MW  (i,j) &
                              + fact_urban(i,j) * CPL_AtmUrb_ATM_FLX_MW  (i,j)
    end do
    end do

    do j = JS, JE
    do i = IS, IE
       CPL_Merged_FLX_MU(i,j) = fact_ocean(i,j) * CPL_AtmOcn_ATM_FLX_MU  (i,j) &
                              + fact_land (i,j) * CPL_AtmLnd_ATM_FLX_MU  (i,j) &
                              + fact_urban(i,j) * CPL_AtmUrb_ATM_FLX_MU  (i,j)
    end do
    end do

    do j = JS, JE
    do i = IS, IE
       CPL_Merged_FLX_MV(i,j) = fact_ocean(i,j) * CPL_AtmOcn_ATM_FLX_MV  (i,j) &
                              + fact_land (i,j) * CPL_AtmLnd_ATM_FLX_MV  (i,j) &
                              + fact_urban(i,j) * CPL_AtmUrb_ATM_FLX_MV  (i,j)
    end do
    end do

    do j = JS, JE
    do i = IS, IE
       CPL_Merged_FLX_SH(i,j) = fact_ocean(i,j) * CPL_AtmOcn_ATM_FLX_SH  (i,j) &
                              + fact_land (i,j) * CPL_AtmLnd_ATM_FLX_SH  (i,j) &
                              + fact_urban(i,j) * CPL_AtmUrb_ATM_FLX_SH  (i,j)
    end do
    end do

    do j = JS, JE
    do i = IS, IE
       CPL_Merged_FLX_LH(i,j) = fact_ocean(i,j) * CPL_AtmOcn_ATM_FLX_LH  (i,j) &
                              + fact_land (i,j) * CPL_AtmLnd_ATM_FLX_LH  (i,j) &
                              + fact_urban(i,j) * CPL_AtmUrb_ATM_FLX_LH  (i,j)
    end do
    end do

    do j = JS, JE
    do i = IS, IE
       CPL_Merged_FLX_QV(i,j) = fact_ocean(i,j) * CPL_AtmOcn_ATM_FLX_evap(i,j) &
                              + fact_land (i,j) * CPL_AtmLnd_ATM_FLX_evap(i,j) &
                              + fact_urban(i,j) * CPL_AtmUrb_ATM_FLX_evap(i,j)
    end do
    end do

    do j = JS, JE
    do i = IS, IE
       CPL_Merged_Z0M   (i,j) = fact_ocean(i,j) * CPL_fromOcn_SFC_Z0M  (i,j) &
                              + fact_land (i,j) * CPL_fromLnd_SFC_Z0M  (i,j) &
                              + fact_urban(i,j) * CPL_fromUrb_SFC_Z0M  (i,j)
    end do
    end do

    do j = JS, JE
    do i = IS, IE
       CPL_Merged_Z0H   (i,j) = fact_ocean(i,j) * CPL_fromOcn_SFC_Z0H  (i,j) &
                              + fact_land (i,j) * CPL_fromLnd_SFC_Z0H  (i,j) &
                              + fact_urban(i,j) * CPL_fromUrb_SFC_Z0H  (i,j)
    end do
    end do

    do j = JS, JE
    do i = IS, IE
       CPL_Merged_Z0E   (i,j) = fact_ocean(i,j) * CPL_fromOcn_SFC_Z0E  (i,j) &
                              + fact_land (i,j) * CPL_fromLnd_SFC_Z0E  (i,j) &
                              + fact_urban(i,j) * CPL_fromUrb_SFC_Z0E  (i,j)
    end do
    end do

    do j = JS, JE
    do i = IS, IE
       CPL_Merged_U10   (i,j) = fact_ocean(i,j) * CPL_AtmOcn_ATM_U10     (i,j) &
                              + fact_land (i,j) * CPL_AtmLnd_ATM_U10     (i,j) &
                              + fact_urban(i,j) * CPL_AtmUrb_ATM_U10     (i,j)
    end do
    end do

    do j = JS, JE
    do i = IS, IE
       CPL_Merged_V10   (i,j) = fact_ocean(i,j) * CPL_AtmOcn_ATM_V10     (i,j) &
                              + fact_land (i,j) * CPL_AtmLnd_ATM_V10     (i,j) &
                              + fact_urban(i,j) * CPL_AtmUrb_ATM_V10     (i,j)
    end do
    end do

    do j = JS, JE
    do i = IS, IE
       CPL_Merged_T2    (i,j) = fact_ocean(i,j) * CPL_AtmOcn_ATM_T2      (i,j) &
                              + fact_land (i,j) * CPL_AtmLnd_ATM_T2      (i,j) &
                              + fact_urban(i,j) * CPL_AtmUrb_ATM_T2      (i,j)
    end do
    end do

    do j = JS, JE
    do i = IS, IE
       CPL_Merged_Q2    (i,j) = fact_ocean(i,j) * CPL_AtmOcn_ATM_Q2      (i,j) &
                              + fact_land (i,j) * CPL_AtmLnd_ATM_Q2      (i,j) &
                              + fact_urban(i,j) * CPL_AtmUrb_ATM_Q2      (i,j)
    end do
    end do

    do j = JS, JE
    do i = IS, IE
       CPL_Merged_FLX_heat(i,j) = fact_ocean(i,j) * CPL_AtmOcn_OCN_FLX_heat(i,j) &
                                + fact_land (i,j) * CPL_AtmLnd_LND_FLX_heat(i,j) &
                                + fact_urban(i,j) * CPL_AtmUrb_URB_FLX_heat(i,j)
    end do
    end do

    if ( STATISTICS_checktotal ) then
       call STAT_total( total, CPL_Merged_SFC_TEMP  (:,:),      'SFC_TEMP  ' )
       call STAT_total( total, CPL_Merged_SFC_albedo(:,:,I_LW), 'SFC_ALB_LW' )
       call STAT_total( total, CPL_Merged_SFC_albedo(:,:,I_SW), 'SFC_ALB_SW' )
       call STAT_total( total, CPL_Merged_FLX_MU    (:,:),      'FLX_MU    ' )
       call STAT_total( total, CPL_Merged_FLX_MV    (:,:),      'FLX_MV    ' )
       call STAT_total( total, CPL_Merged_FLX_MW    (:,:),      'FLX_MW    ' )
       call STAT_total( total, CPL_Merged_FLX_SH    (:,:),      'FLX_SH    ' )
       call STAT_total( total, CPL_Merged_FLX_LH    (:,:),      'FLX_LH    ' )
       call STAT_total( total, CPL_Merged_FLX_QV    (:,:),      'FLX_QV    ' )
       call STAT_total( total, CPL_Merged_U10       (:,:),      'U10       ' )
       call STAT_total( total, CPL_Merged_V10       (:,:),      'V10       ' )
       call STAT_total( total, CPL_Merged_T2        (:,:),      'T2        ' )
       call STAT_total( total, CPL_Merged_Q2        (:,:),      'Q2        ' )
       call STAT_total( total, CPL_Merged_FLX_heat  (:,:),      'FLX_heat  ' )
    endif

    return
  end subroutine CPL_vars_merge

  !-----------------------------------------------------------------------------
  subroutine CPL_putAtm_setup( &
       ATM_Z1 )
    implicit none

    real(RP), intent(in) :: ATM_Z1(IA,JA)
    integer :: i, j
    !---------------------------------------------------------------------------

    do j = JS, JE
    do i = IS, IE
       CPL_fromAtm_ATM_Z1(i,j) = ATM_Z1(i,j)
    end do
    end do

    return
  end subroutine CPL_putAtm_setup

  !-----------------------------------------------------------------------------
  subroutine CPL_putOcn_setup( &
       OCN_TEMP,   &
       SFC_TEMP,   &
       SFC_albedo, &
       SFC_Z0M,    &
       SFC_Z0H,    &
       SFC_Z0E     )
    implicit none

    real(RP), intent(in) :: OCN_TEMP  (IA,JA)
    real(RP), intent(in) :: SFC_TEMP  (IA,JA)
    real(RP), intent(in) :: SFC_albedo(IA,JA,2)
    real(RP), intent(in) :: SFC_Z0M   (IA,JA)
    real(RP), intent(in) :: SFC_Z0H   (IA,JA)
    real(RP), intent(in) :: SFC_Z0E   (IA,JA)

    integer :: i, j
    !---------------------------------------------------------------------------

    do j = JS, JE
    do i = IS, IE
       CPL_fromOcn_OCN_TEMP  (i,j)   = OCN_TEMP  (i,j)
    end do
    end do
    do j = JS, JE
    do i = IS, IE
       CPL_fromOcn_SFC_TEMP  (i,j)   = SFC_TEMP  (i,j)
    end do
    end do
    do j = JS, JE
    do i = IS, IE
       CPL_fromOcn_SFC_albedo(i,j,1) = SFC_albedo(i,j,1)
       CPL_fromOcn_SFC_albedo(i,j,2) = SFC_albedo(i,j,2)
    end do
    end do
    do j = JS, JE
    do i = IS, IE
       CPL_fromOcn_SFC_Z0M   (i,j)   = SFC_Z0M   (i,j)
    end do
    end do
    do j = JS, JE
    do i = IS, IE
       CPL_fromOcn_SFC_Z0H   (i,j)   = SFC_Z0H   (i,j)
    end do
    end do
    do j = JS, JE
    do i = IS, IE
       CPL_fromOcn_SFC_Z0E   (i,j)   = SFC_Z0E   (i,j)
    end do
    end do

    return
  end subroutine CPL_putOcn_setup

  !-----------------------------------------------------------------------------
  subroutine CPL_putLnd_setup( &
       LND_TEMP,   &
       LND_BETA,   &
       SFC_TEMP,   &
       SFC_albedo, &
       LND_TCS,    &
       LND_DZ,     &
       SFC_Z0M,    &
       SFC_Z0H,    &
       SFC_Z0E     )
    implicit none

    real(RP), intent(in) :: LND_TEMP  (IA,JA)
    real(RP), intent(in) :: LND_BETA  (IA,JA)
    real(RP), intent(in) :: SFC_TEMP  (IA,JA)
    real(RP), intent(in) :: SFC_albedo(IA,JA,2)
    real(RP), intent(in) :: LND_TCS   (IA,JA)
    real(RP), intent(in) :: LND_DZ    (IA,JA)
    real(RP), intent(in) :: SFC_Z0M   (IA,JA)
    real(RP), intent(in) :: SFC_Z0H   (IA,JA)
    real(RP), intent(in) :: SFC_Z0E   (IA,JA)

    integer :: i, j
    !---------------------------------------------------------------------------


    do j = JS, JE
    do i = IS, IE
       CPL_fromLnd_LND_TEMP  (i,j)   = LND_TEMP  (i,j)
    end do
    end do
    do j = JS, JE
    do i = IS, IE
       CPL_fromLnd_LND_BETA  (i,j)   = LND_BETA  (i,j)
    end do
    end do
    do j = JS, JE
    do i = IS, IE
       CPL_fromLnd_SFC_TEMP  (i,j)   = SFC_TEMP  (i,j)
    end do
    end do
    do j = JS, JE
    do i = IS, IE
       CPL_fromLnd_SFC_albedo(i,j,1) = SFC_albedo(i,j,1)
       CPL_fromLnd_SFC_albedo(i,j,2) = SFC_albedo(i,j,2)
    end do
    end do
    do j = JS, JE
    do i = IS, IE
       CPL_fromLnd_LND_TCS   (i,j)   = LND_TCS   (i,j)
    end do
    end do
    do j = JS, JE
    do i = IS, IE
       CPL_fromLnd_LND_DZ    (i,j)   = LND_DZ    (i,j)
    end do
    end do
    do j = JS, JE
    do i = IS, IE
       CPL_fromLnd_SFC_Z0M   (i,j)   = SFC_Z0M   (i,j)
    end do
    end do
    do j = JS, JE
    do i = IS, IE
       CPL_fromLnd_SFC_Z0H   (i,j)   = SFC_Z0H   (i,j)
    end do
    end do
    do j = JS, JE
    do i = IS, IE
       CPL_fromLnd_SFC_Z0E   (i,j)   = SFC_Z0E   (i,j)
    end do
    end do

    return
  end subroutine CPL_putLnd_setup

  !-----------------------------------------------------------------------------
  subroutine CPL_putUrb_setup( &
       SFC_TEMP,  &
       SFC_albedo )
    implicit none

    real(RP), intent(in) :: SFC_TEMP  (IA,JA)
    real(RP), intent(in) :: SFC_albedo(IA,JA,2)

    integer :: i, j
    !---------------------------------------------------------------------------

    do j = JS, JE
    do i = IS, IE
       CPL_fromUrb_SFC_TEMP  (i,j)   = SFC_TEMP  (i,j)
    end do
    end do
    do j = JS, JE
    do i = IS, IE
       CPL_fromUrb_SFC_albedo(i,j,1) = SFC_albedo(i,j,1)
       CPL_fromUrb_SFC_albedo(i,j,2) = SFC_albedo(i,j,2)
    end do
    end do

    return
  end subroutine CPL_putUrb_setup

  !-----------------------------------------------------------------------------
  subroutine CPL_putAtm_restart( &
       SFC_TEMP,   &
       SFC_albedo, &
       Z0M,        &
       Z0H,        &
       Z0E,        &
       SFLX_MW,    &
       SFLX_MU,    &
       SFLX_MV,    &
       SFLX_SH,    &
       SFLX_LH,    &
       SFLX_QTRC   )
    implicit none

    real(RP), intent(in) :: SFC_TEMP  (IA,JA)
    real(RP), intent(in) :: SFC_albedo(IA,JA,2)
    real(RP), intent(in) :: Z0M       (IA,JA)
    real(RP), intent(in) :: Z0H       (IA,JA)
    real(RP), intent(in) :: Z0E       (IA,JA)
    real(RP), intent(in) :: SFLX_MW   (IA,JA)
    real(RP), intent(in) :: SFLX_MU   (IA,JA)
    real(RP), intent(in) :: SFLX_MV   (IA,JA)
    real(RP), intent(in) :: SFLX_SH   (IA,JA)
    real(RP), intent(in) :: SFLX_LH   (IA,JA)
    real(RP), intent(in) :: SFLX_QTRC (IA,JA,QA)

    integer :: i, j
    !---------------------------------------------------------------------------

    do j = JS, JE
    do i = IS, IE
       CPL_Merged_SFC_TEMP  (i,j)   = SFC_TEMP  (i,j)  
    end do
    end do
    do j = JS, JE
    do i = IS, IE
       CPL_Merged_SFC_albedo(i,j,1) = SFC_albedo(i,j,1)
       CPL_Merged_SFC_albedo(i,j,2) = SFC_albedo(i,j,2)
    end do
    end do
    do j = JS, JE
    do i = IS, IE
       CPL_Merged_Z0M       (i,j)   = Z0M       (i,j)  
    end do
    end do
    do j = JS, JE
    do i = IS, IE
       CPL_Merged_Z0H       (i,j)   = Z0H       (i,j)  
    end do
    end do
    do j = JS, JE
    do i = IS, IE
       CPL_Merged_Z0E       (i,j)   = Z0E       (i,j)  
    end do
    end do
    do j = JS, JE
    do i = IS, IE
       CPL_Merged_FLX_MW    (i,j)   = SFLX_MW   (i,j)
    end do
    end do
    do j = JS, JE
    do i = IS, IE
       CPL_Merged_FLX_MU    (i,j)   = SFLX_MU   (i,j)
    end do
    end do
    do j = JS, JE
    do i = IS, IE
       CPL_Merged_FLX_MV    (i,j)   = SFLX_MV   (i,j)
    end do
    end do
    do j = JS, JE
    do i = IS, IE
       CPL_Merged_FLX_SH    (i,j)   = SFLX_SH   (i,j)
    end do
    end do
    do j = JS, JE
    do i = IS, IE
       CPL_Merged_FLX_LH    (i,j)   = SFLX_LH   (i,j)
    end do
    end do
    do j = JS, JE
    do i = IS, IE
       CPL_Merged_FLX_QV    (i,j)   = SFLX_QTRC (i,j,1)  ! tentative
    end do
    end do

    return
  end subroutine CPL_putAtm_restart

  !-----------------------------------------------------------------------------
  subroutine CPL_putOcn_restart( &
       FLX_heat,   & ! (in)
       FLX_precip, & ! (in)
       FLX_evap    ) ! (in)
    implicit none

    real(RP), intent(in) :: FLX_heat  (IA,JA)
    real(RP), intent(in) :: FLX_precip(IA,JA)
    real(RP), intent(in) :: FLX_evap  (IA,JA)

    integer :: i, j
    !---------------------------------------------------------------------------

    do j = JS, JE
    do i = IS, IE
       CPL_AtmOcn_OCN_FLX_heat  (i,j)   = FLX_heat  (i,j)  
    end do
    end do
    do j = JS, JE
    do i = IS, IE
       CPL_AtmOcn_OCN_FLX_precip(i,j)   = FLX_precip(i,j)  
    end do
    end do
    do j = JS, JE
    do i = IS, IE
       CPL_AtmOcn_OCN_FLX_evap  (i,j)   = FLX_evap  (i,j)  
    end do
    end do

    return
  end subroutine CPL_putOcn_restart

  !-----------------------------------------------------------------------------
  subroutine CPL_putLnd_restart( &
       FLX_heat,   & ! (in)
       FLX_precip, & ! (in)
       FLX_evap    ) ! (in)
    implicit none

    real(RP), intent(in) :: FLX_heat  (IA,JA)
    real(RP), intent(in) :: FLX_precip(IA,JA)
    real(RP), intent(in) :: FLX_evap  (IA,JA)

    integer :: i, j
    !---------------------------------------------------------------------------

    do j = JS, JE
    do i = IS, IE
       CPL_AtmLnd_LND_FLX_heat  (i,j)   = FLX_heat  (i,j)  
    end do
    end do
    do j = JS, JE
    do i = IS, IE
       CPL_AtmLnd_LND_FLX_precip(i,j)   = FLX_precip(i,j)  
    end do
    end do
    do j = JS, JE
    do i = IS, IE
       CPL_AtmLnd_LND_FLX_evap  (i,j)   = FLX_evap  (i,j)  
    end do
    end do

    return
  end subroutine CPL_putLnd_restart

  !-----------------------------------------------------------------------------
  subroutine CPL_putUrb_restart( &
       FLX_heat,   & ! (in)
       FLX_precip, & ! (in)
       FLX_evap    ) ! (in)
    implicit none

    real(RP), intent(in) :: FLX_heat  (IA,JA)
    real(RP), intent(in) :: FLX_precip(IA,JA)
    real(RP), intent(in) :: FLX_evap  (IA,JA)

    integer :: i, j
    !---------------------------------------------------------------------------

    do j = JS, JE
    do i = IS, IE
       CPL_AtmUrb_URB_FLX_heat  (i,j)   = FLX_heat  (i,j)  
    end do
    end do
    do j = JS, JE
    do i = IS, IE
       CPL_AtmUrb_URB_FLX_precip(i,j)   = FLX_precip(i,j)  
    end do
    end do
    do j = JS, JE
    do i = IS, IE
       CPL_AtmUrb_URB_FLX_evap  (i,j)   = FLX_evap  (i,j)  
    end do
    end do

    return
  end subroutine CPL_putUrb_restart

  !-----------------------------------------------------------------------------
  subroutine CPL_putAtm( &
       ATM_TEMP,   &
       ATM_PRES,   &
       ATM_W,      &
       ATM_U,      &
       ATM_V,      &
       ATM_DENS,   &
       ATM_QTRC,   &
       ATM_PBL,    &
       SFC_PRES,   &
       SFLX_LW_dn, &
       SFLX_SW_dn, &
       SFLX_rain,  &
       SFLX_snow   )
    implicit none

    real(RP), intent(in) :: ATM_TEMP  (IA,JA)
    real(RP), intent(in) :: ATM_PRES  (IA,JA)
    real(RP), intent(in) :: ATM_W     (IA,JA)
    real(RP), intent(in) :: ATM_U     (IA,JA)
    real(RP), intent(in) :: ATM_V     (IA,JA)
    real(RP), intent(in) :: ATM_DENS  (IA,JA)
    real(RP), intent(in) :: ATM_QTRC  (IA,JA,QA)
    real(RP), intent(in) :: ATM_PBL   (IA,JA)
    real(RP), intent(in) :: SFC_PRES  (IA,JA)
    real(RP), intent(in) :: SFLX_LW_dn(IA,JA)
    real(RP), intent(in) :: SFLX_SW_dn(IA,JA)
    real(RP), intent(in) :: SFLX_rain (IA,JA)
    real(RP), intent(in) :: SFLX_snow (IA,JA)

    integer :: i, j
    !---------------------------------------------------------------------------

    do j = JS, JE
    do i = IS, IE
       CPL_fromAtm_ATM_TEMP  (i,j) = ( CPL_fromAtm_ATM_TEMP  (i,j) * CNT_putAtm + ATM_TEMP  (i,j)                   ) &
                                   / ( CNT_putAtm + 1.0_RP )
    end do
    end do
    do j = JS, JE
    do i = IS, IE
       CPL_fromAtm_ATM_PRES  (i,j) = ( CPL_fromAtm_ATM_PRES  (i,j) * CNT_putAtm + ATM_PRES  (i,j)                   ) &
                                   / ( CNT_putAtm + 1.0_RP )
    end do
    end do
    do j = JS, JE
    do i = IS, IE
       CPL_fromAtm_ATM_W     (i,j) = ( CPL_fromAtm_ATM_W     (i,j) * CNT_putAtm + ATM_W     (i,j)                   ) &
                                   / ( CNT_putAtm + 1.0_RP )
    end do
    end do
    do j = JS, JE
    do i = IS, IE
       CPL_fromAtm_ATM_U     (i,j) = ( CPL_fromAtm_ATM_U     (i,j) * CNT_putAtm + ATM_U     (i,j)                   ) &
                                   / ( CNT_putAtm + 1.0_RP )
    end do
    end do
    do j = JS, JE
    do i = IS, IE
       CPL_fromAtm_ATM_V     (i,j) = ( CPL_fromAtm_ATM_V     (i,j) * CNT_putAtm + ATM_V     (i,j)                   ) &
                                   / ( CNT_putAtm + 1.0_RP )
    end do
    end do
    do j = JS, JE
    do i = IS, IE
       CPL_fromAtm_ATM_DENS  (i,j) = ( CPL_fromAtm_ATM_DENS  (i,j) * CNT_putAtm + ATM_DENS  (i,j)                   ) &
                                   / ( CNT_putAtm + 1.0_RP )
    end do
    end do
    do j = JS, JE
    do i = IS, IE
       CPL_fromAtm_ATM_QV    (i,j) = ( CPL_fromAtm_ATM_QV    (i,j) * CNT_putAtm + ATM_QTRC  (i,j,I_QV)              ) &
                                   / ( CNT_putAtm + 1.0_RP )
    end do
    end do
    do j = JS, JE
    do i = IS, IE
       CPL_fromAtm_ATM_PBL   (i,j) = ( CPL_fromAtm_ATM_PBL   (i,j) * CNT_putAtm + ATM_PBL   (i,j)                   ) &
                                   / ( CNT_putAtm + 1.0_RP )
    end do
    end do
    do j = JS, JE
    do i = IS, IE
       CPL_fromAtm_SFC_PRES  (i,j) = ( CPL_fromAtm_SFC_PRES  (i,j) * CNT_putAtm + SFC_PRES  (i,j)                   ) &
                                   / ( CNT_putAtm + 1.0_RP )
    end do
    end do
    do j = JS, JE
    do i = IS, IE
       CPL_fromAtm_FLX_LW_dn (i,j) = ( CPL_fromAtm_FLX_LW_dn (i,j) * CNT_putAtm + SFLX_LW_dn(i,j)                   ) &
                                   / ( CNT_putAtm + 1.0_RP )
    end do
    end do
    do j = JS, JE
    do i = IS, IE
       CPL_fromAtm_FLX_SW_dn (i,j) = ( CPL_fromAtm_FLX_SW_dn (i,j) * CNT_putAtm + SFLX_SW_dn(i,j)                   ) &
                                   / ( CNT_putAtm + 1.0_RP )
    end do
    end do
    do j = JS, JE
    do i = IS, IE
       CPL_fromAtm_FLX_precip(i,j) = ( CPL_fromAtm_FLX_precip(i,j) * CNT_putAtm + SFLX_rain (i,j) + SFLX_snow (i,j) ) &
                                   / ( CNT_putAtm + 1.0_RP )
    end do
    end do

    CNT_putAtm = CNT_putAtm + 1.0_RP

    return
  end subroutine CPL_putAtm

  !-----------------------------------------------------------------------------
  subroutine CPL_putOcn( &
       OCN_TEMP )
    implicit none

    real(RP), intent(in) :: OCN_TEMP(IA,JA)
    integer :: i, j
    !---------------------------------------------------------------------------

    do j = JS, JE
    do i = IS, IE
       CPL_fromOcn_OCN_TEMP(i,j) = ( CPL_fromOcn_OCN_TEMP(i,j) * CNT_putOcn + OCN_TEMP(i,j) ) &
                                 / ( CNT_putOcn + 1.0_RP )
    end do
    end do

    CNT_putOcn = CNT_putOcn + 1.0_RP

    return
  end subroutine CPL_putOcn

  !-----------------------------------------------------------------------------
  subroutine CPL_putLnd( &
       LND_TEMP, &
       LND_BETA  )
    implicit none

    real(RP), intent(in) :: LND_TEMP(IA,JA)
    real(RP), intent(in) :: LND_BETA(IA,JA)
    integer :: i, j
    !---------------------------------------------------------------------------

    do j = JS, JE
    do i = IS, IE
       CPL_fromLnd_LND_TEMP(i,j) = ( CPL_fromLnd_LND_TEMP(i,j) * CNT_putLnd + LND_TEMP(i,j) ) &
                                 / ( CNT_putLnd + 1.0_RP )
    end do
    end do
    do j = JS, JE
    do i = IS, IE
       CPL_fromLnd_LND_BETA(i,j) = ( CPL_fromLnd_LND_BETA(i,j) * CNT_putLnd + LND_BETA(i,j) ) &
                                 / ( CNT_putLnd + 1.0_RP )
    end do
    end do

    CNT_putLnd = CNT_putLnd + 1.0_RP

    return
  end subroutine CPL_putLnd

  !-----------------------------------------------------------------------------
  subroutine CPL_putUrb
    implicit none
    !---------------------------------------------------------------------------

    return
  end subroutine CPL_putUrb

  !-----------------------------------------------------------------------------
  subroutine CPL_getAtm( &
       SFC_TEMP,   &
       SFC_albedo, &
       Z0M,        &
       Z0H,        &
       Z0E,        &
       Uabs10,     &
       U10,        &
       V10,        &
       T2,         &
       Q2,         &
       FLX_heat    )
    implicit none

    real(RP), intent(out) :: SFC_TEMP  (IA,JA)
    real(RP), intent(out) :: SFC_albedo(IA,JA,2)
    real(RP), intent(out) :: Z0M       (IA,JA)
    real(RP), intent(out) :: Z0H       (IA,JA)
    real(RP), intent(out) :: Z0E       (IA,JA)
    real(RP), intent(out) :: Uabs10    (IA,JA)
    real(RP), intent(out) :: U10       (IA,JA)
    real(RP), intent(out) :: V10       (IA,JA)
    real(RP), intent(out) :: T2        (IA,JA)
    real(RP), intent(out) :: Q2        (IA,JA)
    real(RP), intent(out) :: FLX_heat  (IA,JA)

    integer :: i, j
    !---------------------------------------------------------------------------

    do j = JS, JE
    do i = IS, IE
       SFC_TEMP  (i,j)   = CPL_Merged_SFC_TEMP  (i,j)
    end do
    end do
    do j = JS, JE
    do i = IS, IE
       SFC_albedo(i,j,1) = CPL_Merged_SFC_albedo(i,j,1)
       SFC_albedo(i,j,2) = CPL_Merged_SFC_albedo(i,j,2)
    end do
    end do
    do j = JS, JE
    do i = IS, IE
       Z0M       (i,j)   = CPL_Merged_Z0M       (i,j)
    end do
    end do
    do j = JS, JE
    do i = IS, IE
       Z0H       (i,j)   = CPL_Merged_Z0H       (i,j)
    end do
    end do
    do j = JS, JE
    do i = IS, IE
       Z0E       (i,j)   = CPL_Merged_Z0E       (i,j)
    end do
    end do
    do j = JS, JE
    do i = IS, IE
       U10       (i,j)   = CPL_Merged_U10       (i,j)
    end do
    end do
    do j = JS, JE
    do i = IS, IE
       V10       (i,j)   = CPL_Merged_V10       (i,j)
    end do
    end do
    do j = JS, JE
    do i = IS, IE
       T2        (i,j)   = CPL_Merged_T2        (i,j)
    end do
    end do
    do j = JS, JE
    do i = IS, IE
       Q2        (i,j)   = CPL_Merged_Q2        (i,j)
    end do
    end do
    do j = JS, JE
    do i = IS, IE
       FLX_heat  (i,j)   = CPL_Merged_FLX_heat  (i,j)
    end do
    end do

    do j = JS, JE
    do i = IS, IE
       Uabs10(i,j) = sqrt( CPL_Merged_U10(i,j)*CPL_Merged_U10(i,j) &
                         + CPL_Merged_V10(i,j)*CPL_Merged_V10(i,j) )
    end do
    end do

    return
  end subroutine CPL_getAtm

  !-----------------------------------------------------------------------------
  subroutine CPL_getAtm_SF( &
       SFLX_MW,   &
       SFLX_MU,   &
       SFLX_MV,   &
       SFLX_SH,   &
       SFLX_LH,   &
       SFLX_QTRC  )
    implicit none

    real(RP), intent(out) :: SFLX_MW  (IA,JA)
    real(RP), intent(out) :: SFLX_MU  (IA,JA)
    real(RP), intent(out) :: SFLX_MV  (IA,JA)
    real(RP), intent(out) :: SFLX_SH  (IA,JA)
    real(RP), intent(out) :: SFLX_LH  (IA,JA)
    real(RP), intent(out) :: SFLX_QTRC(IA,JA,QA)

    integer :: i, j, iq
    !---------------------------------------------------------------------------

    do j = JS, JE
    do i = IS, IE
       SFLX_MW  (i,j)   = CPL_Merged_FLX_MW(i,j)
    end do
    end do
    do j = JS, JE
    do i = IS, IE
       SFLX_MU  (i,j)   = CPL_Merged_FLX_MU(i,j)
    end do
    end do
    do j = JS, JE
    do i = IS, IE
       SFLX_MV  (i,j)   = CPL_Merged_FLX_MV(i,j)
    end do
    end do
    do j = JS, JE
    do i = IS, IE
       SFLX_SH  (i,j)   = CPL_Merged_FLX_SH(i,j)
    end do
    end do
    do j = JS, JE
    do i = IS, IE
       SFLX_LH  (i,j)   = CPL_Merged_FLX_LH(i,j)
    end do
    end do
    do iq = 1, QA
    do j = JS, JE
    do i = IS, IE
       SFLX_QTRC(i,j,iq) = 0.0_RP                 ! tentative
    end do
    end do
    end do
    do j = JS, JE
    do i = IS, IE
       SFLX_QTRC(i,j,1) = CPL_Merged_FLX_QV(i,j) ! tentative
    end do
    end do

    CNT_getAtmLnd = 0.0_RP
    CNT_getAtmUrb = 0.0_RP
    CNT_getAtmOcn = 0.0_RP

    return
  end subroutine CPL_getAtm_SF

  !-----------------------------------------------------------------------------
  subroutine CPL_getOcn( &
       SFC_TEMP,   & ! (out)
       SFC_albedo, & ! (out)
       SFC_Z0M,    & ! (out)
       SFC_Z0H,    & ! (out)
       SFC_Z0E,    & ! (out)
       FLX_heat,   & ! (out)
       FLX_precip, & ! (out)
       FLX_evap    ) ! (out)
    implicit none

    real(RP), intent(out) :: SFC_TEMP  (IA,JA)
    real(RP), intent(out) :: SFC_albedo(IA,JA,2)
    real(RP), intent(out) :: SFC_Z0M   (IA,JA)
    real(RP), intent(out) :: SFC_Z0H   (IA,JA)
    real(RP), intent(out) :: SFC_Z0E   (IA,JA)
    real(RP), intent(out) :: FLX_heat  (IA,JA)
    real(RP), intent(out) :: FLX_precip(IA,JA)
    real(RP), intent(out) :: FLX_evap  (IA,JA)

    integer :: i, j
    !---------------------------------------------------------------------------

    do j = JS, JE
    do i = IS, IE
       SFC_TEMP  (i,j)   = CPL_fromOcn_SFC_TEMP     (i,j)
    end do
    end do
    do j = JS, JE
    do i = IS, IE
       SFC_albedo(i,j,1) = CPL_fromOcn_SFC_albedo   (i,j,1)
       SFC_albedo(i,j,2) = CPL_fromOcn_SFC_albedo   (i,j,2)
    end do
    end do
    do j = JS, JE
    do i = IS, IE
       SFC_Z0M   (i,j)   = CPL_fromOcn_SFC_Z0M      (i,j)
    end do
    end do
    do j = JS, JE
    do i = IS, IE
       SFC_Z0H   (i,j)   = CPL_fromOcn_SFC_Z0H      (i,j)
    end do
    end do
    do j = JS, JE
    do i = IS, IE
       SFC_Z0E   (i,j)   = CPL_fromOcn_SFC_Z0E      (i,j)
    end do
    end do
    do j = JS, JE
    do i = IS, IE
       FLX_heat  (i,j)   = CPL_AtmOcn_OCN_FLX_heat  (i,j)
    end do
    end do
    do j = JS, JE
    do i = IS, IE
       FLX_precip(i,j)   = CPL_AtmOcn_OCN_FLX_precip(i,j)
    end do
    end do
    do j = JS, JE
    do i = IS, IE
       FLX_evap  (i,j)   = CPL_AtmOcn_OCN_FLX_evap  (i,j)
    end do
    end do

    CNT_getOcn = 0.0_RP

    return
  end subroutine CPL_getOcn

  !-----------------------------------------------------------------------------
  subroutine CPL_getLnd( &
      SFC_TEMP,   & ! (out)
      SFC_albedo, & ! (out)
      FLX_heat,   & ! (out)
      FLX_precip, & ! (out)
      FLX_evap    ) ! (out)
    implicit none

    real(RP), intent(out) :: SFC_TEMP  (IA,JA)
    real(RP), intent(out) :: SFC_albedo(IA,JA,2)
    real(RP), intent(out) :: FLX_heat  (IA,JA)
    real(RP), intent(out) :: FLX_precip(IA,JA)
    real(RP), intent(out) :: FLX_evap  (IA,JA)

    integer :: i, j
    !---------------------------------------------------------------------------

    do j = JS, JE
    do i = IS, JE
       SFC_TEMP  (i,j)   = CPL_fromLnd_SFC_TEMP  (i,j)
    end do
    end do
    do j = JS, JE
    do i = IS, JE
       SFC_albedo(i,j,1) = CPL_fromLnd_SFC_albedo(i,j,1)
       SFC_albedo(i,j,2) = CPL_fromLnd_SFC_albedo(i,j,2)
    end do
    end do
    do j = JS, JE
    do i = IS, IE
       FLX_heat  (i,j) = CPL_AtmLnd_LND_FLX_heat  (i,j)
    end do
    end do
    do j = JS, JE
    do i = IS, IE
       FLX_precip(i,j) = CPL_AtmLnd_LND_FLX_precip(i,j)
    end do
    end do
    do j = JS, JE
    do i = IS, IE
       FLX_evap  (i,j) = CPL_AtmLnd_LND_FLX_evap  (i,j)
    end do
    end do

    CNT_getLnd = 0.0_RP

    return
  end subroutine CPL_getLnd

  !-----------------------------------------------------------------------------
  subroutine CPL_getUrb( &
      SFC_TEMP,   & ! (out)
      SFC_albedo, & ! (out)
      FLX_heat,   & ! (out)
      FLX_precip, & ! (out)
      FLX_evap    ) ! (out)
    implicit none

    real(RP), intent(out) :: SFC_TEMP  (IA,JA)
    real(RP), intent(out) :: SFC_albedo(IA,JA,2)
    real(RP), intent(out) :: FLX_heat  (IA,JA)
    real(RP), intent(out) :: FLX_precip(IA,JA)
    real(RP), intent(out) :: FLX_evap  (IA,JA)

    integer :: i, j
    !---------------------------------------------------------------------------

    do j = JS, JE
    do i = IS, JE
       SFC_TEMP      (i,j)   = CPL_fromUrb_SFC_TEMP     (i,j)
    end do
    end do
    do j = JS, JE
    do i = IS, JE
       SFC_albedo    (i,j,1) = CPL_fromUrb_SFC_albedo   (i,j,1)
       SFC_albedo    (i,j,2) = CPL_fromUrb_SFC_albedo   (i,j,2)
    end do
    end do
    do j = JS, JE
    do i = IS, JE
       FLX_heat  (i,j)   = CPL_AtmUrb_URB_FLX_heat  (i,j)
    end do
    end do
    do j = JS, JE
    do i = IS, JE
       FLX_precip(i,j)   = CPL_AtmUrb_URB_FLX_precip(i,j)
    end do
    end do
    do j = JS, JE
    do i = IS, JE
       FLX_evap  (i,j)   = CPL_AtmUrb_URB_FLX_evap  (i,j)
    end do
    end do

    CNT_getUrb = 0.0_RP

    return
  end subroutine CPL_getUrb

end module mod_CPL_vars
