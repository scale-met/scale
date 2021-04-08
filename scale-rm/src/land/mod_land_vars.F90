!-------------------------------------------------------------------------------
!> module LAND Variables
!!
!! @par Description
!!          Container for land variables
!!
!! @author Team SCALE
!<
!-------------------------------------------------------------------------------
#include "scalelib.h"
module mod_land_vars
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_io
  use scale_prof
  use scale_debug
  use scale_land_grid_cartesC_index
  use scale_tracer
  use scale_cpl_sfc_index
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: LAND_vars_setup
  public :: LAND_vars_finalize
  public :: LAND_vars_restart_read
  public :: LAND_vars_restart_write
  public :: LAND_vars_history
  public :: LAND_vars_monitor
  public :: LAND_vars_check

  public :: LAND_vars_restart_create
  public :: LAND_vars_restart_open
  public :: LAND_vars_restart_def_var
  public :: LAND_vars_restart_enddef
  public :: LAND_vars_restart_close

  public :: convert_WS2VWC

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  logical,                public :: LAND_RESTART_OUTPUT                = .false.         !< Output restart file?

  character(len=H_LONG),  public :: LAND_RESTART_IN_BASENAME           = ''              !< Basename of the input  file
  logical,                public :: LAND_RESTART_IN_AGGREGATE                            !< Switch to use aggregate file
  logical,                public :: LAND_RESTART_IN_POSTFIX_TIMELABEL  = .false.         !< Add timelabel to the basename of input  file?
  character(len=H_LONG),  public :: LAND_RESTART_OUT_BASENAME          = ''              !< Basename of the output file
  logical,                public :: LAND_RESTART_OUT_AGGREGATE                           !< Switch to use aggregate file
  logical,                public :: LAND_RESTART_OUT_POSTFIX_TIMELABEL = .true.          !< Add timelabel to the basename of output file?
  character(len=H_MID),   public :: LAND_RESTART_OUT_TITLE             = 'LAND restart' !< Title    of the output file
  character(len=H_SHORT), public :: LAND_RESTART_OUT_DTYPE             = 'DEFAULT'       !< REAL4 or REAL8

  ! prognostic variables
  real(RP), public, allocatable :: LAND_TEMP      (:,:,:)   !< temperature of each soil layer [K]
  real(RP), public, allocatable :: LAND_WATER     (:,:,:)   !< moisture of each soil layer    [m3/m3]
  real(RP), public, allocatable :: LAND_ICE       (:,:,:)   !< ice of each soil layer         [m3/m3]
  real(RP), public, allocatable :: LAND_SFC_TEMP  (:,:)     !< land surface skin temperature  [K]
  real(RP), public, allocatable :: LAND_SFC_albedo(:,:,:,:) !< land surface albedo (direct/diffuse,IR/near-IR/VIS) (0-1)

  ! for snow model
  real(RP), public, allocatable :: SNOW_SFC_TEMP (:,:) !< snow surface temperature     [K]
  real(RP), public, allocatable :: SNOW_SWE      (:,:) !< snow water equivalent        [kg/m2]
  real(RP), public, allocatable :: SNOW_Depth    (:,:) !< snow depth                   [m]
  real(RP), public, allocatable :: SNOW_Dzero    (:,:) !< snow depth at melting point  [m]
  real(RP), public, allocatable :: SNOW_nosnowsec(:,:) !< sec while no snow            [s]

  ! tendency variables
  real(RP), public, allocatable :: LAND_TEMP_t (:,:,:) !< tendency of LAND_TEMP
  real(RP), public, allocatable :: LAND_WATER_t(:,:,:) !< tendency of LAND_WATER
  real(RP), public, allocatable :: LAND_ICE_t  (:,:,:) !< tendency of LAND_ICE

  ! surface flux for land
  real(RP), public, allocatable :: LAND_SFLX_GH   (:,:) !< land surface heat flux            [J/m2/s]
  real(RP), public, allocatable :: LAND_SFLX_water(:,:) !< land surface water flux           [kg/m2/s]
  real(RP), public, allocatable :: LAND_SFLX_ENGI (:,:) !< land surface internal energy flux [J/m2/s]

  ! surface flux for atmosphere
  real(RP), public, allocatable :: LAND_SFLX_MW  (:,:)   !< land surface w-momentum flux    [kg/m2/s]
  real(RP), public, allocatable :: LAND_SFLX_MU  (:,:)   !< land surface u-momentum flux    [kg/m2/s]
  real(RP), public, allocatable :: LAND_SFLX_MV  (:,:)   !< land surface v-momentum flux    [kg/m2/s]
  real(RP), public, allocatable :: LAND_SFLX_SH  (:,:)   !< land surface sensible heat flux [J/m2/s]
  real(RP), public, allocatable :: LAND_SFLX_LH  (:,:)   !< land surface latent heat flux   [J/m2/s]
  real(RP), public, allocatable :: LAND_SFLX_QTRC(:,:,:) !< land surface tracer flux        [kg/m2/s]
  real(RP), public, allocatable :: LAND_U10      (:,:)   !< land surface velocity u at 10m [m/s]
  real(RP), public, allocatable :: LAND_V10      (:,:)   !< land surface velocity v at 10m [m/s]
  real(RP), public, allocatable :: LAND_T2       (:,:)   !< land surface temperature at 2m [K]
  real(RP), public, allocatable :: LAND_Q2       (:,:)   !< land surface water vapor at 2m [kg/kg]
  real(RP), public, allocatable, target :: LAND_Ustar(:,:) !< friction velocity         [m/s]
  real(RP), public, allocatable, target :: LAND_Tstar(:,:) !< temperature scale         [K]
  real(RP), public, allocatable, target :: LAND_Qstar(:,:) !< moisture scale            [kg/kg]
  real(RP), public, allocatable, target :: LAND_Wstar(:,:) !< convective velocity scale [m/s]
  real(RP), public, allocatable, target :: LAND_RLmo (:,:) !< inversed Obukhov length   [1/m]
  real(RP), public, pointer :: SOIL_Ustar(:,:)
  real(RP), public, pointer :: SOIL_Tstar(:,:)
  real(RP), public, pointer :: SOIL_Qstar(:,:)
  real(RP), public, pointer :: SOIL_Wstar(:,:)
  real(RP), public, pointer :: SOIL_RLmo (:,:)
  real(RP), public, allocatable :: SNOW_Ustar(:,:)
  real(RP), public, allocatable :: SNOW_Tstar(:,:)
  real(RP), public, allocatable :: SNOW_Qstar(:,:)
  real(RP), public, allocatable :: SNOW_Wstar(:,:)
  real(RP), public, allocatable :: SNOW_RLmo (:,:)

  real(RP), public, allocatable :: LAND_RUNOFF     (:,:) !< runoff of the land water      [kg/m2/s]
  real(RP), public, allocatable :: LAND_RUNOFF_ENGI(:,:) !< internal energy of the runoff [J/m2/s]


  ! recieved atmospheric variables
  real(RP), public, allocatable :: ATMOS_TEMP       (:,:)
  real(RP), public, allocatable :: ATMOS_PRES       (:,:)
  real(RP), public, allocatable :: ATMOS_W          (:,:)
  real(RP), public, allocatable :: ATMOS_U          (:,:)
  real(RP), public, allocatable :: ATMOS_V          (:,:)
  real(RP), public, allocatable :: ATMOS_DENS       (:,:)
  real(RP), public, allocatable :: ATMOS_QV         (:,:)
  real(RP), public, allocatable :: ATMOS_PBL        (:,:)
  real(RP), public, allocatable :: ATMOS_SFC_DENS   (:,:)
  real(RP), public, allocatable :: ATMOS_SFC_PRES   (:,:)
  real(RP), public, allocatable :: ATMOS_SFLX_rad_dn(:,:,:,:)
  real(RP), public, allocatable :: ATMOS_cosSZA     (:,:)
  real(RP), public, allocatable :: ATMOS_SFLX_water (:,:)
  real(RP), public, allocatable :: ATMOS_SFLX_ENGI  (:,:)


  logical, public :: SNOW_flag


  real(RP), public, allocatable :: LAND_PROPERTY(:,:,:) !< land surface property

  character(len=H_LONG), public :: LAND_PROPERTY_IN_FILENAME  = '' !< the file of land parameter table

  integer,  public, parameter   :: LAND_PROPERTY_nmax = 11
  integer,  public, parameter   :: I_WaterLimit       =  1 ! maximum  soil moisture           [m3/m3]
  integer,  public, parameter   :: I_WaterCritical    =  2 ! critical soil moisture           [m3/m3]
  integer,  public, parameter   :: I_StomataResist    =  3 ! stomata resistance               [1/s]
  integer,  public, parameter   :: I_ThermalCond      =  4 ! thermal conductivity for soil    [W/K/m]
  integer,  public, parameter   :: I_HeatCapacity     =  5 ! heat capacity for soil           [J/K/m3]
  integer,  public, parameter   :: I_WaterDiff        =  6 ! moisture diffusivity in the soil [m2/s]
  integer,  public, parameter   :: I_ALBLW            =  7 ! surface albedo for long  wave    [1]
  integer,  public, parameter   :: I_ALBSW            =  8 ! surface albedo for short wave    [1]
  integer,  public, parameter   :: I_Z0M              =  9 ! roughness length for momemtum    [m]
  integer,  public, parameter   :: I_Z0H              = 10 ! roughness length for heat        [m]
  integer,  public, parameter   :: I_Z0E              = 11 ! roughness length for vapor       [m]

  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  private :: LAND_param_read

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  logical, private :: LAND_VARS_CHECKRANGE = .false.

  ! index of the prognostic variables
  integer, private, parameter :: VMAX              = 16 !< number of the variables
  integer, private, parameter :: I_TEMP            =  1
  integer, private, parameter :: I_WATER           =  2
  integer, private, parameter :: I_ICE             =  3
  integer, private, parameter :: I_WATERDS         =  4
  integer, private, parameter :: I_SFC_TEMP        =  5
  integer, private, parameter :: I_SFC_ALB_IR_dir  =  6
  integer, private, parameter :: I_SFC_ALB_IR_dif  =  7
  integer, private, parameter :: I_SFC_ALB_NIR_dir =  8
  integer, private, parameter :: I_SFC_ALB_NIR_dif =  9
  integer, private, parameter :: I_SFC_ALB_VIS_dir = 10
  integer, private, parameter :: I_SFC_ALB_VIS_dif = 11
  integer, private, parameter :: I_SNOW_SFC_TEMP   = 12
  integer, private, parameter :: I_SNOW_SWE        = 13
  integer, private, parameter :: I_SNOW_Depth      = 14
  integer, private, parameter :: I_SNOW_Dzero      = 15
  integer, private, parameter :: I_SNOW_nosnowsec  = 16

  character(len=H_SHORT), private            :: VAR_NAME(VMAX) !< name  of the variables
  character(len=H_MID),   private            :: VAR_DESC(VMAX) !< desc. of the variables
  character(len=H_MID),   private            :: VAR_STDN(VMAX) !< standard name of the variables
  character(len=H_SHORT), private            :: VAR_UNIT(VMAX) !< unit  of the variables
  integer,                private            :: VAR_ID(VMAX)   !< ID    of the variables
  integer,                private            :: restart_fid = -1  ! file ID

  data VAR_NAME / 'LAND_TEMP',            &
                  'LAND_WATER',           &
                  'LAND_ICE',             &
                  'LAND_DSAT',            &
                  'LAND_SFC_TEMP',        &
                  'LAND_SFC_ALB_IR_dir',  &
                  'LAND_SFC_ALB_IR_dif',  &
                  'LAND_SFC_ALB_NIR_dir', &
                  'LAND_SFC_ALB_NIR_dif', &
                  'LAND_SFC_ALB_VIS_dir', &
                  'LAND_SFC_ALB_VIS_dif', &
                  'LAND_SNOW_SFC_TEMP',   &
                  'LAND_SNOW_SWE',        &
                  'LAND_SNOW_Depth',      &
                  'LAND_SNOW_Dzero',      &
                  'LAND_SNOW_nosnowsec'   /

  data VAR_DESC / 'temperature at each soil layer',          &
                  'moisture at each soil layer',             &
                  'ice at each soil layer',                  &
                  'degree of saturation at each soil layer', &
                  'land surface skin temperature',           &
                  'land surface albedo for IR (direct)',     &
                  'land surface albedo for IR (diffuse)',    &
                  'land surface albedo for NIR (direct)',    &
                  'land surface albedo for NIR (diffuse)',   &
                  'land surface albedo for VIS (direct)',    &
                  'land surface albedo for VIS (diffuse)',   &
                  'Snow surface temperature',                &
                  'Snow water equivalent',                   &
                  'Snow depth',                              &
                  'Snow depth at melting point',             &
                  'Time duration without snow'               /

  data VAR_STDN / 'soil_temperature', &
                  'volume_fraction_of_condensed_water_in_soil', &
                  '', &
                  'volume_fraction_of_condensed_water_in_soil_at_field_capacity', &
                  'surface_temperature_where_land', &
                  '', &
                  '', &
                  '', &
                  '', &
                  '', &
                  '', &
                  '', &
                  '', &
                  '', &
                  '', &
                  '' /

  data VAR_UNIT / 'K',     &
                  'm3/m3', &
                  'm3/m3', &
                  '1',     &
                  'K',     &
                  '1',     &
                  '1',     &
                  '1',     &
                  '1',     &
                  '1',     &
                  '1',     &
                  'K',     &
                  'kg/m2', &
                  'm',     &
                  'm',     &
                  's'      /

  real(RP), private, allocatable :: LAND_PROPERTY_table(:,:)

  logical,  private :: LAND_RESTART_IN_CHECK_COORDINATES = .true.

  ! for monitor
  integer, parameter :: IM_TEMP      = 1
  integer, parameter :: IM_WATER     = 2
  integer, parameter :: IM_ICE       = 3
  integer, parameter :: IM_SFC       = 4
  integer, parameter :: IM_ROFF      = 5
  integer, parameter :: IM_MASFLX    = 6
  integer, parameter :: IM_ENGI      = 7
  integer, parameter :: IM_W_ENGI    = 8
  integer, parameter :: IM_I_ENGI    = 9
  integer, parameter :: IM_ENGSFC_GH = 10
  integer, parameter :: IM_ENGSFC_EI = 11
  integer, parameter :: IM_ROFF_EI   = 12
  integer, parameter :: IM_ENGFLX    = 13
  integer, parameter :: IM_max = 13
  integer, private   :: MONIT_id(IM_max)

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine LAND_vars_setup
    use scale_prc, only: &
       PRC_abort
    use scale_const, only: &
       UNDEF => CONST_UNDEF
    use scale_comm_cartesC, only: &
       COMM_vars8, &
       COMM_wait
    use scale_landuse, only: &
       LANDUSE_index_PFT, &
       LANDUSE_PFT_nmin, &
       LANDUSE_PFT_nmax
    use mod_land_admin, only: &
       SNOW_TYPE
    use scale_monitor, only: &
       MONITOR_reg
    implicit none

    namelist / PARAM_LAND_VARS / &
       LAND_RESTART_IN_BASENAME,           &
       LAND_RESTART_IN_AGGREGATE,          &
       LAND_RESTART_IN_POSTFIX_TIMELABEL,  &
       LAND_RESTART_IN_CHECK_COORDINATES,  &
       LAND_RESTART_OUTPUT,                &
       LAND_RESTART_OUT_BASENAME,          &
       LAND_RESTART_OUT_AGGREGATE,         &
       LAND_RESTART_OUT_POSTFIX_TIMELABEL, &
       LAND_RESTART_OUT_TITLE,             &
       LAND_RESTART_OUT_DTYPE,             &
       LAND_VARS_CHECKRANGE

    integer :: ierr
    integer :: i, j, iv, p
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("LAND_vars_setup",*) 'Setup'

    select case ( SNOW_TYPE )
    case ( 'NONE', 'OFF' )
       SNOW_flag = .false.
    case default
       SNOW_flag = .true.
    end select

    allocate( LAND_TEMP      (LKMAX,LIA,LJA)               )
    allocate( LAND_WATER     (LKMAX,LIA,LJA)               )
    allocate( LAND_ICE       (LKMAX,LIA,LJA)               )
    allocate( LAND_SFC_TEMP  (LIA,LJA)                     )
    allocate( LAND_SFC_albedo(LIA,LJA,N_RAD_DIR,N_RAD_RGN) )

    LAND_TEMP      (:,:,:)   = UNDEF
    LAND_WATER     (:,:,:)   = UNDEF
    LAND_ICE       (:,:,:)   = UNDEF
    LAND_SFC_TEMP  (:,:)     = UNDEF
    LAND_SFC_albedo(:,:,:,:) = UNDEF

    if ( SNOW_flag ) then
       allocate( SNOW_SFC_TEMP (LIA,LJA) )
       allocate( SNOW_SWE      (LIA,LJA) )
       allocate( SNOW_Depth    (LIA,LJA) )
       allocate( SNOW_Dzero    (LIA,LJA) )
       allocate( SNOW_nosnowsec(LIA,LJA) )
       SNOW_SFC_TEMP (:,:) = UNDEF
       SNOW_SWE      (:,:) = UNDEF
       SNOW_Depth    (:,:) = UNDEF
       SNOW_Dzero    (:,:) = UNDEF
       SNOW_nosnowsec(:,:) = UNDEF
    end if

    allocate( LAND_TEMP_t (LKMAX,LIA,LJA) )
    allocate( LAND_WATER_t(LKMAX,LIA,LJA) )
    allocate( LAND_ICE_t  (LKMAX,LIA,LJA) )
    LAND_TEMP_t (:,:,:) = UNDEF
    LAND_WATER_t(:,:,:) = UNDEF
    LAND_ICE_t  (:,:,:) = UNDEF

    allocate( LAND_SFLX_GH   (LIA,LJA) )
    allocate( LAND_SFLX_water(LIA,LJA) )
    allocate( LAND_SFLX_ENGI (LIA,LJA) )
    LAND_SFLX_GH   (:,:) = UNDEF
    LAND_SFLX_water(:,:) = UNDEF
    LAND_SFLX_ENGI  (:,:) = UNDEF

    allocate( LAND_RUNOFF     (LIA,LJA) )
    allocate( LAND_RUNOFF_ENGI(LIA,LJA) )
    LAND_RUNOFF     (:,:) = UNDEF
    LAND_RUNOFF_ENGI(:,:) = UNDEF

    allocate( LAND_SFLX_MW  (LIA,LJA) )
    allocate( LAND_SFLX_MU  (LIA,LJA) )
    allocate( LAND_SFLX_MV  (LIA,LJA) )
    allocate( LAND_SFLX_SH  (LIA,LJA) )
    allocate( LAND_SFLX_LH  (LIA,LJA) )
    allocate( LAND_SFLX_QTRC(LIA,LJA,QA) )
    LAND_SFLX_MW  (:,:)   = UNDEF
    LAND_SFLX_MU  (:,:)   = UNDEF
    LAND_SFLX_MV  (:,:)   = UNDEF
    LAND_SFLX_SH  (:,:)   = UNDEF
    LAND_SFLX_LH  (:,:)   = UNDEF
    LAND_SFLX_QTRC(:,:,:) = UNDEF

    allocate( LAND_U10      (LIA,LJA) )
    allocate( LAND_V10      (LIA,LJA) )
    allocate( LAND_T2       (LIA,LJA) )
    allocate( LAND_Q2       (LIA,LJA) )
    LAND_U10      (:,:) = UNDEF
    LAND_V10      (:,:) = UNDEF
    LAND_T2       (:,:) = UNDEF
    LAND_Q2       (:,:) = UNDEF

    allocate( LAND_Ustar(LIA,LJA) )
    allocate( LAND_Tstar(LIA,LJA) )
    allocate( LAND_Qstar(LIA,LJA) )
    allocate( LAND_Wstar(LIA,LJA) )
    allocate( LAND_RLmo (LIA,LJA) )
    LAND_Ustar(:,:) = UNDEF
    LAND_Tstar(:,:) = UNDEF
    LAND_Qstar(:,:) = UNDEF
    LAND_Wstar(:,:) = UNDEF
    LAND_RLmo (:,:) = UNDEF
    if ( SNOW_flag ) then
       allocate( SOIL_Ustar(LIA,LJA) )
       allocate( SOIL_Tstar(LIA,LJA) )
       allocate( SOIL_Qstar(LIA,LJA) )
       allocate( SOIL_Wstar(LIA,LJA) )
       allocate( SOIL_RLmo (LIA,LJA) )
       SOIL_Ustar(:,:) = UNDEF
       SOIL_Tstar(:,:) = UNDEF
       SOIL_Qstar(:,:) = UNDEF
       SOIL_Wstar(:,:) = UNDEF
       SOIL_RLmo (:,:) = UNDEF
    else
       SOIL_Ustar => LAND_Ustar
       SOIL_Tstar => LAND_Tstar
       SOIL_Qstar => LAND_Qstar
       SOIL_Wstar => LAND_Wstar
       SOIL_RLmo  => LAND_RLmo
    end if
    if ( SNOW_flag ) then
       allocate( SNOW_Ustar(LIA,LJA) )
       allocate( SNOW_Tstar(LIA,LJA) )
       allocate( SNOW_Qstar(LIA,LJA) )
       allocate( SNOW_Wstar(LIA,LJA) )
       allocate( SNOW_RLmo (LIA,LJA) )
       SNOW_Ustar(:,:) = UNDEF
       SNOW_Tstar(:,:) = UNDEF
       SNOW_Qstar(:,:) = UNDEF
       SNOW_Wstar(:,:) = UNDEF
       SNOW_RLmo (:,:) = UNDEF
    end if

    allocate( ATMOS_TEMP       (LIA,LJA)                     )
    allocate( ATMOS_PRES       (LIA,LJA)                     )
    allocate( ATMOS_W          (LIA,LJA)                     )
    allocate( ATMOS_U          (LIA,LJA)                     )
    allocate( ATMOS_V          (LIA,LJA)                     )
    allocate( ATMOS_DENS       (LIA,LJA)                     )
    allocate( ATMOS_QV         (LIA,LJA)                     )
    allocate( ATMOS_PBL        (LIA,LJA)                     )
    allocate( ATMOS_SFC_DENS   (LIA,LJA)                     )
    allocate( ATMOS_SFC_PRES   (LIA,LJA)                     )
    allocate( ATMOS_SFLX_rad_dn(LIA,LJA,N_RAD_DIR,N_RAD_RGN) )
    allocate( ATMOS_cosSZA     (LIA,LJA)                     )
    allocate( ATMOS_SFLX_water (LIA,LJA)                     )
    allocate( ATMOS_SFLX_ENGI  (LIA,LJA)                     )
    ATMOS_TEMP       (:,:)     = UNDEF
    ATMOS_PRES       (:,:)     = UNDEF
    ATMOS_W          (:,:)     = UNDEF
    ATMOS_U          (:,:)     = UNDEF
    ATMOS_V          (:,:)     = UNDEF
    ATMOS_DENS       (:,:)     = UNDEF
    ATMOS_QV         (:,:)     = UNDEF
    ATMOS_PBL        (:,:)     = UNDEF
    ATMOS_SFC_DENS   (:,:)     = UNDEF
    ATMOS_SFC_PRES   (:,:)     = UNDEF
    ATMOS_SFLX_rad_dn(:,:,:,:) = UNDEF
    ATMOS_cosSZA     (:,:)     = UNDEF
    ATMOS_SFLX_water (:,:)     = UNDEF
    ATMOS_SFLX_ENGI  (:,:)     = UNDEF

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_LAND_VARS,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("LAND_vars_setup",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("LAND_vars_setup",*) 'Not appropriate names in namelist PARAM_LAND_VARS. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_LAND_VARS)

    LOG_NEWLINE
    LOG_INFO("LAND_vars_setup",*) 'List of prognostic variables (LAND) '
    LOG_INFO_CONT('(1x,A,A24,A,A48,A,A12,A)') &
               '      |', 'VARNAME                 ','|', &
               'DESCRIPTION                                     ', '[', 'UNIT        ', ']'
    do iv = 1, VMAX
       LOG_INFO_CONT('(1x,A,I3,A,A24,A,A48,A,A12,A)') &
                  'NO.',iv,'|',VAR_NAME(iv),'|',VAR_DESC(iv),'[',VAR_UNIT(iv),']'
    enddo

    LOG_NEWLINE
    if ( LAND_RESTART_IN_BASENAME /= '' ) then
       LOG_INFO("LAND_vars_setup",*) 'Restart input?  : YES, file = ', trim(LAND_RESTART_IN_BASENAME)
       LOG_INFO("LAND_vars_setup",*) 'Add timelabel?  : ', LAND_RESTART_IN_POSTFIX_TIMELABEL
    else
       LOG_INFO("LAND_vars_setup",*) 'Restart input?  : NO'
    endif
    if (       LAND_RESTART_OUTPUT             &
         .AND. LAND_RESTART_OUT_BASENAME /= '' ) then
       LOG_INFO("LAND_vars_setup",*) 'Restart output? : YES, file = ', trim(LAND_RESTART_OUT_BASENAME)
       LOG_INFO("LAND_vars_setup",*) 'Add timelabel?  : ', LAND_RESTART_OUT_POSTFIX_TIMELABEL
    else
       LOG_INFO("LAND_vars_setup",*) 'Restart output? : NO'
       LAND_RESTART_OUTPUT = .false.
    endif

    ! Read land property table
    allocate( LAND_PROPERTY_table(LANDUSE_PFT_nmin:LANDUSE_PFT_nmax,LAND_PROPERTY_nmax) )
    LAND_PROPERTY_table(:,:) = UNDEF

    call LAND_param_read

    ! Apply land property to 2D map
    allocate( LAND_PROPERTY(LIA,LJA,LAND_PROPERTY_nmax) )

    ! tentative, mosaic is off
    do p = 1, LAND_PROPERTY_nmax
    do j = LJS, LJE
    do i = LIS, LIE
       LAND_PROPERTY(i,j,p) = LAND_PROPERTY_table( LANDUSE_index_PFT(i,j,1), p )
    enddo
    enddo
    enddo

    do p = 1, LAND_PROPERTY_nmax
       call COMM_vars8( LAND_PROPERTY(:,:,p), p )
    enddo
    do p = 1, LAND_PROPERTY_nmax
       call COMM_wait ( LAND_PROPERTY(:,:,p), p )
    enddo

    ! monitor
    call MONITOR_reg( 'LND_TEMP',       'land temperature',                'K m3', & ! (in)
                      MONIT_id(IM_TEMP),                                           & ! (out)
                      dim_type='LXY', is_tendency=.false.                          ) ! (in)
    call MONITOR_reg( 'LND_WATER',      'land water',                      'kg',   & ! (in)
                      MONIT_id(IM_WATER),                                          & ! (out)
                      dim_type='LXY', is_tendency=.false.                          ) ! (in)
    call MONITOR_reg( 'LND_ICE',        'land ice',                        'kg',   & ! (in)
                      MONIT_id(IM_ICE),                                            & ! (out)
                      dim_type='LXY', is_tendency=.false.                          ) ! (in)
    call MONITOR_reg( 'LND_MASSFC',     'SFC water flux',                  'kg',   & ! (in)
                      MONIT_id(IM_SFC),                                            & ! (out)
                      dim_type='XY',  is_tendency=.true.                           ) ! (in)
    call MONITOR_reg( 'LND_ROFF',       'runoff water',                    'kg',   & ! (in)
                      MONIT_id(IM_ROFF),                                           & ! (out)
                      dim_type='XY',  is_tendency=.true.                           ) ! (in)
    call MONITOR_reg( 'LND_MASFLX',     'total mass change',               'kg',   & ! (in)
                      MONIT_id(IM_MASFLX),                                         & ! (out)
                      dim_type='XY',  is_tendency=.true.                           ) ! (in)
    call MONITOR_reg( 'LND_ENGI',       'total internal energy',           'J',    & ! (in)
                      MONIT_id(IM_ENGI),                                           & ! (out)
                      dim_type='LXY', is_tendency=.false.                          ) ! (in)
    call MONITOR_reg( 'LND_WTR_ENGI',   'water internal energy',           'J',    & ! (in)
                      MONIT_id(IM_W_ENGI),                                         & ! (out)
                      dim_type='LXY', is_tendency=.false.                          ) ! (in)
    call MONITOR_reg( 'LND_ICE_ENGI',   'ice internal energy',             'J',    & ! (in)
                      MONIT_id(IM_I_ENGI),                                         & ! (out)
                      dim_type='LXY', is_tendency=.false.                          ) ! (in)
    call MONITOR_reg( 'LND_ENGSFC_GH',  'SFC ground heat flux',            'J',    & ! (in)
                      MONIT_id(IM_ENGSFC_GH),                                      & ! (out)
                      dim_type='XY',  is_tendency=.true.                           ) ! (in)
    call MONITOR_reg( 'LND_ENGSFC_EI',  'SFC internal energy flux',        'J',    & ! (in)
                      MONIT_id(IM_ENGSFC_EI),                                      & ! (out)
                      dim_type='XY',  is_tendency=.true.                           ) ! (in)
    call MONITOR_reg( 'LND_ROFF_EI',    'internal energy of runoff water', 'J',    & ! (in)
                      MONIT_id(IM_ROFF_EI),                                        & ! (out)
                      dim_type='XY',  is_tendency=.true.                           ) ! (in)
    call MONITOR_reg( 'LND_ENGFLX',     'total internal energy change',    'J',    & ! (in)
                      MONIT_id(IM_ENGFLX),                                         & ! (out)
                      dim_type='XY',  is_tendency=.true.                           ) ! (in)

    return
  end subroutine LAND_vars_setup

  !-----------------------------------------------------------------------------
  !> Finalize
  subroutine LAND_vars_finalize
    use scale_landuse, only: &
       LANDUSE_index_PFT, &
       LANDUSE_PFT_nmin, &
       LANDUSE_PFT_nmax
    use mod_land_admin, only: &
       SNOW_TYPE
    implicit none
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("LAND_vars_finalize",*) 'Finalize'

    select case ( SNOW_TYPE )
    case ( 'NONE', 'OFF' )
       SNOW_flag = .false.
    case default
       SNOW_flag = .true.
    end select

    deallocate( LAND_TEMP       )
    deallocate( LAND_WATER      )
    deallocate( LAND_ICE        )
    deallocate( LAND_SFC_TEMP   )
    deallocate( LAND_SFC_albedo )

    if ( SNOW_flag ) then
       deallocate( SNOW_SFC_TEMP  )
       deallocate( SNOW_SWE       )
       deallocate( SNOW_Depth     )
       deallocate( SNOW_Dzero     )
       deallocate( SNOW_nosnowsec )
    end if

    deallocate( LAND_TEMP_t  )
    deallocate( LAND_WATER_t )
    deallocate( LAND_ICE_t   )

    deallocate( LAND_SFLX_GH    )
    deallocate( LAND_SFLX_water )
    deallocate( LAND_SFLX_ENGI  )

    deallocate( LAND_RUNOFF      )
    deallocate( LAND_RUNOFF_ENGI )

    deallocate( LAND_SFLX_MW   )
    deallocate( LAND_SFLX_MU   )
    deallocate( LAND_SFLX_MV   )
    deallocate( LAND_SFLX_SH   )
    deallocate( LAND_SFLX_LH   )
    deallocate( LAND_SFLX_QTRC )

    deallocate( LAND_U10       )
    deallocate( LAND_V10       )
    deallocate( LAND_T2        )
    deallocate( LAND_Q2        )

    deallocate( LAND_Ustar )
    deallocate( LAND_Tstar )
    deallocate( LAND_Qstar )
    deallocate( LAND_Wstar )
    deallocate( LAND_RLmo  )
    if ( SNOW_flag ) then
       deallocate( SOIL_Ustar )
       deallocate( SOIL_Tstar )
       deallocate( SOIL_Qstar )
       deallocate( SOIL_Wstar )
       deallocate( SOIL_RLmo  )
    end if
    if ( SNOW_flag ) then
       deallocate( SNOW_Ustar )
       deallocate( SNOW_Tstar )
       deallocate( SNOW_Qstar )
       deallocate( SNOW_Wstar )
       deallocate( SNOW_RLmo  )
    end if

    deallocate( ATMOS_TEMP        )
    deallocate( ATMOS_PRES        )
    deallocate( ATMOS_W           )
    deallocate( ATMOS_U           )
    deallocate( ATMOS_V           )
    deallocate( ATMOS_DENS        )
    deallocate( ATMOS_QV          )
    deallocate( ATMOS_PBL         )
    deallocate( ATMOS_SFC_DENS    )
    deallocate( ATMOS_SFC_PRES    )
    deallocate( ATMOS_SFLX_rad_dn )
    deallocate( ATMOS_cosSZA      )
    deallocate( ATMOS_SFLX_water  )
    deallocate( ATMOS_SFLX_ENGI   )

    deallocate( LAND_PROPERTY_table )

    deallocate( LAND_PROPERTY )

    return
  end subroutine LAND_vars_finalize

  !-----------------------------------------------------------------------------
  !> Open land restart file for read
  subroutine LAND_vars_restart_open
    use scale_time, only: &
       TIME_gettimelabel
    use scale_file_cartesC, only: &
       FILE_CARTESC_open, &
       FILE_CARTESC_check_coordinates
    use mod_land_admin, only: &
       LAND_do
    implicit none

    character(len=19)     :: timelabel
    character(len=H_LONG) :: basename
    !---------------------------------------------------------------------------

    call PROF_rapstart('LND_Restart', 1)

    LOG_NEWLINE
    LOG_INFO("LAND_vars_restart_open",*) 'Open restart file (LAND) '

    if ( LAND_do .and. LAND_RESTART_IN_BASENAME /= '' ) then

       if ( LAND_RESTART_IN_POSTFIX_TIMELABEL ) then
          call TIME_gettimelabel( timelabel )
          basename = trim(LAND_RESTART_IN_BASENAME)//'_'//trim(timelabel)
       else
          basename = trim(LAND_RESTART_IN_BASENAME)
       endif

       LOG_INFO("LAND_vars_restart_open",*) 'basename: ', trim(basename)

       call FILE_CARTESC_open( basename, restart_fid, aggregate=LAND_RESTART_IN_AGGREGATE )

       if ( LAND_RESTART_IN_CHECK_COORDINATES ) then
          call FILE_CARTESC_check_coordinates( restart_fid, land=.true. )
       end if

    else
       LOG_INFO("LAND_vars_restart_open",*) 'restart file for land is not specified.'
    endif

    call PROF_rapend('LND_Restart', 1)

    return
  end subroutine LAND_vars_restart_open

  !-----------------------------------------------------------------------------
  !> Read land restart
  subroutine LAND_vars_restart_read
    use scale_prc, only: &
       PRC_abort
    use scale_file, only: &
       FILE_get_aggregate
    use scale_file_cartesC, only: &
       FILE_CARTESC_read, &
       FILE_CARTESC_flush
    implicit none
    !---------------------------------------------------------------------------

    call PROF_rapstart('LND_Restart', 1)

    if ( restart_fid /= -1 ) then
       LOG_NEWLINE
       LOG_INFO("LAND_vars_restart_read",*) 'Read from restart file (LAND) '

       call FILE_CARTESC_read( restart_fid, VAR_NAME(I_TEMP),            'LXY', & ! [IN]
                               LAND_TEMP      (:,:,:)                           ) ! [OUT]
       call FILE_CARTESC_read( restart_fid, VAR_NAME(I_WATER),           'LXY', & ! [IN]
                               LAND_WATER     (:,:,:)                           ) ! [OUT]
       call FILE_CARTESC_read( restart_fid, VAR_NAME(I_ICE),             'LXY', & ! [IN]
                               LAND_ICE       (:,:,:),                          & ! [OUT]
                               allow_missing = .true.                           ) ! [IN]
       call FILE_CARTESC_read( restart_fid, VAR_NAME(I_SFC_TEMP),        'XY',  & ! [IN]
                               LAND_SFC_TEMP  (:,:)                             ) ! [OUT]
       call FILE_CARTESC_read( restart_fid, VAR_NAME(I_SFC_ALB_IR_dir),  'XY',  & ! [IN]
                               LAND_SFC_albedo(:,:,I_R_direct ,I_R_IR )         ) ! [OUT]
       call FILE_CARTESC_read( restart_fid, VAR_NAME(I_SFC_ALB_IR_dif),  'XY',  & ! [IN]
                               LAND_SFC_albedo(:,:,I_R_diffuse,I_R_IR )         ) ! [OUT]
       call FILE_CARTESC_read( restart_fid, VAR_NAME(I_SFC_ALB_NIR_dir), 'XY',  & ! [IN]
                               LAND_SFC_albedo(:,:,I_R_direct ,I_R_NIR)         ) ! [OUT]
       call FILE_CARTESC_read( restart_fid, VAR_NAME(I_SFC_ALB_NIR_dif), 'XY',  & ! [IN]
                               LAND_SFC_albedo(:,:,I_R_diffuse,I_R_NIR)         ) ! [OUT]
       call FILE_CARTESC_read( restart_fid, VAR_NAME(I_SFC_ALB_VIS_dir), 'XY',  & ! [IN]
                               LAND_SFC_albedo(:,:,I_R_direct ,I_R_VIS)         ) ! [OUT]
       call FILE_CARTESC_read( restart_fid, VAR_NAME(I_SFC_ALB_VIS_dif), 'XY',  & ! [IN]
                               LAND_SFC_albedo(:,:,I_R_diffuse,I_R_VIS)         ) ! [OUT]

       if ( SNOW_flag ) then
          call FILE_CARTESC_read( restart_fid, VAR_NAME(I_SNOW_SFC_TEMP),  'XY', & ! [OUT]
                                  SNOW_SFC_TEMP(:,:)                             ) ! [IN]
          call FILE_CARTESC_read( restart_fid, VAR_NAME(I_SNOW_SWE),       'XY', & ! [OUT]
                                  SNOW_SWE(:,:)                                  ) ! [IN]
          call FILE_CARTESC_read( restart_fid, VAR_NAME(I_SNOW_Depth),     'XY', & ! [OUT]
                                  SNOW_Depth(:,:)                                ) ! [IN]
          call FILE_CARTESC_read( restart_fid, VAR_NAME(I_SNOW_Dzero),     'XY', & ! [OUT]
                                  SNOW_Dzero(:,:)                                ) ! [IN]
          call FILE_CARTESC_read( restart_fid, VAR_NAME(I_SNOW_nosnowsec), 'XY', & ! [OUT]
                                  SNOW_nosnowsec(:,:)                            ) ! [IN]
       end if

       if( FILE_get_AGGREGATE(restart_fid) ) call FILE_CARTESC_flush( restart_fid ) ! commit all pending read requests

       call LAND_vars_check( force = .true. )
    else
       LOG_ERROR("LAND_vars_restart_read",*) 'invalid restart file ID for land.'
       call PRC_abort
    endif

    call PROF_rapend('LND_Restart', 1)

    return
  end subroutine LAND_vars_restart_read

  !-----------------------------------------------------------------------------
  !> History output set for land variables
  subroutine LAND_vars_history
    use scale_file_history, only: &
       FILE_HISTORY_in
    use scale_atmos_hydrometeor, only: &
       I_QV
    implicit none

    real(RP) :: LAND_WATERDS(LKMAX,LIA,LJA)
    integer :: k, i, j
    !---------------------------------------------------------------------------

    call PROF_rapstart('LND_History', 1)


    call FILE_HISTORY_in( LAND_TEMP (:,:,:), VAR_NAME(I_TEMP),  VAR_DESC(I_TEMP),  VAR_UNIT(I_TEMP),  dim_type='LXY', standard_name=VAR_STDN(I_TEMP) )
    call FILE_HISTORY_in( LAND_WATER(:,:,:), VAR_NAME(I_WATER), VAR_DESC(I_WATER), VAR_UNIT(I_WATER), dim_type='LXY', standard_name=VAR_STDN(I_WATER) )
    call FILE_HISTORY_in( LAND_ICE  (:,:,:), VAR_NAME(I_ICE),   VAR_DESC(I_ICE),   VAR_UNIT(I_ICE), dim_type='LXY', standard_name=VAR_STDN(I_ICE) )
    do j = LJS, LJE
    do i = LIS, LIE
    do k = 1, LKMAX
       LAND_WATERDS(k,i,j) = LAND_WATER(k,i,j) / LAND_PROPERTY(i,j,I_WaterLimit)
    end do
    end do
    end do
    call FILE_HISTORY_in( LAND_WATERDS(:,:,:), VAR_NAME(I_WATERDS), VAR_DESC(I_WATERDS), VAR_UNIT(I_WATERDS), dim_type='LXY', fill_halo=.true., standard_name=VAR_STDN(I_WATERDS) )


    call FILE_HISTORY_in( LAND_SFC_TEMP  (:,:),                     VAR_NAME(I_SFC_TEMP),                                     &
                          VAR_DESC(I_SFC_TEMP),        VAR_UNIT(I_SFC_TEMP),        standard_name=VAR_STDN(I_SFC_TEMP)        )
    call FILE_HISTORY_in( LAND_SFC_albedo(:,:,I_R_direct ,I_R_IR ), VAR_NAME(I_SFC_ALB_IR_dir),                               &
                          VAR_DESC(I_SFC_ALB_IR_dir),  VAR_UNIT(I_SFC_ALB_IR_dir),  standard_name=VAR_STDN(I_SFC_ALB_IR_dir)  )
    call FILE_HISTORY_in( LAND_SFC_albedo(:,:,I_R_diffuse,I_R_IR ), VAR_NAME(I_SFC_ALB_IR_dif),                               &
                          VAR_DESC(I_SFC_ALB_IR_dif),  VAR_UNIT(I_SFC_ALB_IR_dif),  standard_name=VAR_STDN(I_SFC_ALB_IR_dif)  )
    call FILE_HISTORY_in( LAND_SFC_albedo(:,:,I_R_direct ,I_R_NIR), VAR_NAME(I_SFC_ALB_NIR_dir),                              &
                          VAR_DESC(I_SFC_ALB_NIR_dir), VAR_UNIT(I_SFC_ALB_NIR_dir), standard_name=VAR_STDN(I_SFC_ALB_NIR_dir) )
    call FILE_HISTORY_in( LAND_SFC_albedo(:,:,I_R_diffuse,I_R_NIR), VAR_NAME(I_SFC_ALB_NIR_dif),                              &
                          VAR_DESC(I_SFC_ALB_NIR_dif), VAR_UNIT(I_SFC_ALB_NIR_dif), standard_name=VAR_STDN(I_SFC_ALB_NIR_dif) )
    call FILE_HISTORY_in( LAND_SFC_albedo(:,:,I_R_direct ,I_R_VIS), VAR_NAME(I_SFC_ALB_VIS_dir),                              &
                          VAR_DESC(I_SFC_ALB_VIS_dir), VAR_UNIT(I_SFC_ALB_VIS_dir), standard_name=VAR_STDN(I_SFC_ALB_VIS_dir) )
    call FILE_HISTORY_in( LAND_SFC_albedo(:,:,I_R_diffuse,I_R_VIS), VAR_NAME(I_SFC_ALB_VIS_dif),                              &
                          VAR_DESC(I_SFC_ALB_VIS_dif), VAR_UNIT(I_SFC_ALB_VIS_dif), standard_name=VAR_STDN(I_SFC_ALB_VIS_dif) )
    if ( SNOW_flag ) then
       ! snow model
       call FILE_HISTORY_in( SNOW_SFC_TEMP (:,:), VAR_NAME(I_SNOW_SFC_TEMP),                                               &
                             VAR_DESC(I_SNOW_SFC_TEMP),  VAR_UNIT(I_SNOW_SFC_TEMP),  standard_name=VAR_STDN(I_SNOW_SFC_TEMP)  )
       call FILE_HISTORY_in( SNOW_SWE      (:,:), VAR_NAME(I_SNOW_SWE),                                                    &
                             VAR_DESC(I_SNOW_SWE),       VAR_UNIT(I_SNOW_SWE),       standard_name=VAR_STDN(I_SNOW_SWE)       )
       call FILE_HISTORY_in( SNOW_Depth    (:,:), VAR_NAME(I_SNOW_Depth),                                                  &
                             VAR_DESC(I_SNOW_Depth),     VAR_UNIT(I_SNOW_Depth),     standard_name=VAR_STDN(I_SNOW_Depth)     )
       call FILE_HISTORY_in( SNOW_Dzero    (:,:), VAR_NAME(I_SNOW_Dzero),                                                  &
                             VAR_DESC(I_SNOW_Dzero),     VAR_UNIT(I_SNOW_Dzero),     standard_name=VAR_STDN(I_SNOW_Dzero)     )
       call FILE_HISTORY_in( SNOW_nosnowsec(:,:), VAR_NAME(I_SNOW_nosnowsec),                                                 &
                             VAR_DESC(I_SNOW_nosnowsec), VAR_UNIT(I_SNOW_nosnowsec), standard_name=VAR_STDN(I_SNOW_nosnowsec) )
    end if

    call FILE_HISTORY_in( LAND_SFLX_GH   (:,:), 'LAND_SFLX_GH',                     &
                          'land subsurface heat flux (downward)',         'J/m2/s'  )
    call FILE_HISTORY_in( LAND_SFLX_water(:,:), 'LAND_SFLX_water',                  &
                          'land surface water mass flux (downward)',      'kg/m2/s' )
    call FILE_HISTORY_in( LAND_SFLX_ENGI (:,:), 'LAND_SFLX_ENGI',                   &
                          'land surface internal energy flux (downward)', 'kg/m2/s' )

    call FILE_HISTORY_in( LAND_RUNOFF     (:,:), 'LAND_RUNOFF',        &
                          'runoff water',                    'kg/m2/s' )
    call FILE_HISTORY_in( LAND_RUNOFF_ENGI(:,:), 'LAND_RUNOFF_ENGI',   &
                          'internal energy of runoff water', 'J/m2/s'  )

    call FILE_HISTORY_in( LAND_SFLX_MW   (:,:),      'LAND_SFLX_MW',            &
                          'land surface w-momentum flux (upward)',    'kg/m2/s' )
    call FILE_HISTORY_in( LAND_SFLX_MU   (:,:),      'LAND_SFLX_MU',            &
                          'land surface u-momentum flux (upward)',    'kg/m2/s' )
    call FILE_HISTORY_in( LAND_SFLX_MV   (:,:),      'LAND_SFLX_MV',            &
                          'land surface v-momentum flux (upward)',    'kg/m2/s' )
    call FILE_HISTORY_in( LAND_SFLX_SH   (:,:),      'LAND_SFLX_SH',            &
                          'land surface sensible heat flux (upward)', 'J/m2/s'  )
    call FILE_HISTORY_in( LAND_SFLX_LH   (:,:),      'LAND_SFLX_LH',            &
                          'land surface latent heat flux (upward)',   'J/m2/s'  )
    if ( I_QV > 0 ) &
    call FILE_HISTORY_in( LAND_SFLX_QTRC (:,:,I_QV), 'LAND_SFLX_evap',          &
                          'land surface water vapor flux (upward)',   'kg/m2/s' )

    call FILE_HISTORY_in( LAND_U10       (:,:),  'LAND_U10',                    &
                          'land 10m x-wind',                          'm/s'     )
    call FILE_HISTORY_in( LAND_V10       (:,:),  'LAND_V10',                    &
                          'land 10m y-wind',                          'm/s'     )
    call FILE_HISTORY_in( LAND_T2        (:,:),  'LAND_T2',                     &
                          'land 2m temperature',                      'K'       )
    call FILE_HISTORY_in( LAND_Q2        (:,:),  'LAND_Q2',                     &
                          'land 2m specific humidity',                'kg/kg'   )

    call FILE_HISTORy_in( LAND_Ustar     (:,:), 'LAND_Ustar',                   &
                          'land friction velocity',                   'm/s'     )
    call FILE_HISTORY_in( LAND_Tstar     (:,:), 'LAND_Tstar',                   &
                          'land temperature scale',                   'K'       )
    call FILE_HISTORY_in( LAND_Qstar     (:,:), 'LAND_Qstar',                   &
                          'land moisture scale',                      'kg/kg'   )
    call FILE_HISTORY_in( LAND_Wstar     (:,:),  'LAND_Wstar',                  &
                          'land convective velocity scale',           'm/s'     )
    call FILE_HISTORY_in( LAND_RLmo      (:,:),  'LAND_RLmo',                   &
                          'land inversed Obukhov length',             '1/m'     )
    if ( SNOW_flag ) then
       ! soil
       call FILE_HISTORy_in( SOIL_Ustar     (:,:), 'SOIL_Ustar',                   &
                             'soil friction velocity',                   'm/s'     )
       call FILE_HISTORY_in( SOIL_Tstar     (:,:), 'SOIL_Tstar',                   &
                             'soil temperature scale',                   'K'       )
       call FILE_HISTORY_in( SOIL_Qstar     (:,:), 'SOIL_Qstar',                   &
                             'soil moisture scale',                      'kg/kg'   )
       call FILE_HISTORY_in( SOIL_Wstar     (:,:),  'SOIL_Wstar',                  &
                             'soil convective velocity scale',           'm/s'     )
       call FILE_HISTORY_in( SOIL_RLmo      (:,:),  'SOIL_RLmo',                   &
                             'soil inversed Obukhov length',             '1/m'     )
       ! snow pack
       call FILE_HISTORy_in( SNOW_Ustar     (:,:), 'SNOW_Ustar',                   &
                             'snow friction velocity',                   'm/s'     )
       call FILE_HISTORY_in( SNOW_Tstar     (:,:), 'SNOW_Tstar',                   &
                             'snow temperature scale',                   'K'       )
       call FILE_HISTORY_in( SNOW_Qstar     (:,:), 'SNOW_Qstar',                   &
                             'snow moisture scale',                      'kg/kg'   )
       call FILE_HISTORY_in( SNOW_Wstar     (:,:),  'SNOW_Wstar',                  &
                             'snow convective velocity scale',           'm/s'     )
       call FILE_HISTORY_in( SNOW_RLmo      (:,:),  'SNOW_RLmo',                   &
                             'snow inversed Obukhov length',             '1/m'     )
    end if

    call PROF_rapend  ('LND_History', 1)

    return
  end subroutine LAND_vars_history

  !-----------------------------------------------------------------------------
  !> monitor output
  subroutine LAND_vars_monitor
    use scale_const, only: &
       DWATR => CONST_DWATR, &
       DICE  => CONST_DICE
    use scale_atmos_hydrometeor, only: &
       CV_WATER, &
       CV_ICE,   &
       LHF
    use scale_monitor, only: &
       MONITOR_put
    implicit none

    real(RP)               :: WORK3D(LKA,LIA,LJA)
    real(RP)               :: WORK2D(LIA,LJA)

    integer  :: k, i, j
    !---------------------------------------------------------------------------

    call MONITOR_put( MONIT_id(IM_TEMP),  LAND_TEMP(:,:,:) )
    if ( MONIT_id(IM_WATER) > 0 ) then
       !$omp parallel do
       do j = LJS, LJE
       do i = LIS, LIE
       do k = LKS, LKE
          WORK3D(k,i,j) = LAND_WATER(k,i,j) * DWATR
       end do
       end do
       end do
       call MONITOR_put( MONIT_id(IM_WATER), WORK3D(:,:,:) )
    end if
    if ( MONIT_id(IM_ICE) > 0 ) then
       !$omp parallel do
       do j = LJS, LJE
       do i = LIS, LIE
       do k = LKS, LKE
          WORK3D(k,i,j) = LAND_ICE(k,i,j) * DICE
       end do
       end do
       end do
       call MONITOR_put( MONIT_id(IM_ICE),   WORK3D(:,:,:) )
    end if


    ! mass budget
    call MONITOR_put( MONIT_id(IM_SFC),  LAND_SFLX_water(:,:) )
    call MONITOR_put( MONIT_id(IM_ROFF), LAND_RUNOFF    (:,:) )
    if ( MONIT_id(IM_MASFLX) > 0 ) then
       !$omp parallel do
       do j = LJS, LJE
       do i = LIS, LIE
          WORK2D(i,j) = LAND_SFLX_water(i,j) - LAND_RUNOFF(i,j)
       end do
       end do
       call MONITOR_put( MONIT_id(IM_MASFLX), WORK2D(:,:) )
    end if

    ! energy budget
    if ( MONIT_id(IM_ENGI) > 0 ) then
       !$omp parallel do
       do j = LJS, LJE
       do i = LIS, LIE
       do k = LKS, LKE
          WORK3D(k,i,j) = ( LAND_PROPERTY(i,j,I_HeatCapacity) * ( 1.0_RP - LAND_PROPERTY(i,j,I_WaterLimit) )  & ! soil particles
                          + CV_WATER * DWATR * LAND_WATER(k,i,j) & ! land water
                          + CV_ICE   * DICE  * LAND_ICE  (k,i,j) & ! land ice
                          ) * LAND_TEMP(k,i,j) &
                        - LHF * DICE * LAND_ICE(k,i,j)
       end do
       end do
       end do
       call MONITOR_put( MONIT_id(IM_ENGI), WORK3D(:,:,:) )
    end if
    if ( MONIT_id(IM_W_ENGI) > 0 ) then
       !$omp parallel do
       do j = LJS, LJE
       do i = LIS, LIE
       do k = LKS, LKE
          WORK3D(k,i,j) = CV_WATER * DWATR * LAND_WATER(k,i,j) * LAND_TEMP(k,i,j)
       end do
       end do
       end do
       call MONITOR_put( MONIT_id(IM_W_ENGI), WORK3D(:,:,:) )
    end if
    if ( MONIT_id(IM_I_ENGI) > 0 ) then
       !$omp parallel do
       do j = LJS, LJE
       do i = LIS, LIE
       do k = LKS, LKE
          WORK3D(k,i,j) = ( CV_ICE * LAND_TEMP(k,i,j) - LHF ) * DICE * LAND_ICE(k,i,j)
       end do
       end do
       end do
       call MONITOR_put( MONIT_id(IM_I_ENGI), WORK3D(:,:,:) )
    end if


    call MONITOR_put( MONIT_id(IM_ENGSFC_GH), LAND_SFLX_GH    (:,:) )
    call MONITOR_put( MONIT_id(IM_ENGSFC_EI), LAND_SFLX_ENGI  (:,:) )
    call MONITOR_put( MONIT_id(IM_ROFF_EI),   LAND_RUNOFF_ENGI(:,:) )
    if ( MONIT_id(IM_ENGFLX) > 0 ) then
       !$omp parallel do
       do j = LJS, LJE
       do i = LIS, LIE
          WORK2D(i,j) = LAND_SFLX_GH(i,j) + LAND_SFLX_ENGI(i,j) &
                      - LAND_RUNOFF_ENGI(i,j)
       end do
       end do
       call MONITOR_put( MONIT_id(IM_ENGFLX), WORK2D(:,:) )
    end if


    return
  end subroutine LAND_vars_monitor

  !-----------------------------------------------------------------------------
  !> Budget monitor for land
  subroutine LAND_vars_check( force )
    use scale_statistics, only: &
       STATISTICS_checktotal, &
       STATISTICS_total
    use scale_land_grid_cartesC_real, only: &
       LAND_GRID_CARTESC_REAL_AREA,    &
       LAND_GRID_CARTESC_REAL_TOTAREA, &
       LAND_GRID_CARTESC_REAL_VOL,     &
       LAND_GRID_CARTESC_REAL_TOTVOL
    use scale_landuse, only: &
       LANDUSE_exists_land
    implicit none
    logical, intent(in), optional :: force
    logical :: check
    !---------------------------------------------------------------------------

    if ( present(force) ) then
       check = force
    else
       check = LAND_VARS_CHECKRANGE
    end if

    if ( check ) then
       call VALCHECK( LKA, LKS, LKE, LIA, LIS, LIE, LJA, LJS, LJE,                 &
                      LAND_TEMP      (:,:,:),                   0.0_RP, 1000.0_RP, &
                      VAR_NAME(I_TEMP),                        __FILE__, __LINE__, &
                      mask = LANDUSE_exists_land(:,:)                              )
       call VALCHECK( LKA, LKS, LKE, LIA, LIS, LIE, LJA, LJS, LJE,                 &
                      LAND_WATER     (:,:,:),                   0.0_RP,    1.0_RP, &
                      VAR_NAME(I_WATER),                       __FILE__, __LINE__, &
                      mask = LANDUSE_exists_land(:,:)                              )
       call VALCHECK( LKA, LKS, LKE, LIA, LIS, LIE, LJA, LJS, LJE,                 &
                      LAND_ICE       (:,:,:),                   0.0_RP,    1.0_RP, &
                      VAR_NAME(I_ICE),                         __FILE__, __LINE__, &
                      mask = LANDUSE_exists_land(:,:)                              )
       call VALCHECK( LIA, LIS, LIE, LJA, LJS, LJE,                                &
                      LAND_SFC_TEMP  (:,:),                     0.0_RP, 1000.0_RP, &
                      VAR_NAME(I_SFC_TEMP),                    __FILE__, __LINE__, &
                      mask = LANDUSE_exists_land(:,:)                              )
       call VALCHECK( LIA, LIS, LIE, LJA, LJS, LJE,                                &
                      LAND_SFC_albedo(:,:,I_R_direct ,I_R_IR ), 0.0_RP,    2.0_RP, &
                      VAR_NAME(I_SFC_ALB_IR_dir ),             __FILE__, __LINE__, &
                      mask = LANDUSE_exists_land(:,:)                              )
       call VALCHECK( LIA, LIS, LIE, LJA, LJS, LJE,                                &
                      LAND_SFC_albedo(:,:,I_R_diffuse,I_R_IR ), 0.0_RP,    2.0_RP, &
                      VAR_NAME(I_SFC_ALB_IR_dif ),             __FILE__, __LINE__, &
                      mask = LANDUSE_exists_land(:,:)                              )
       call VALCHECK( LIA, LIS, LIE, LJA, LJS, LJE,                                &
                      LAND_SFC_albedo(:,:,I_R_direct ,I_R_NIR), 0.0_RP,    2.0_RP, &
                      VAR_NAME(I_SFC_ALB_NIR_dir),             __FILE__, __LINE__, &
                      mask = LANDUSE_exists_land(:,:)                              )
       call VALCHECK( LIA, LIS, LIE, LJA, LJS, LJE,                                &
                      LAND_SFC_albedo(:,:,I_R_diffuse,I_R_NIR), 0.0_RP,    2.0_RP, &
                      VAR_NAME(I_SFC_ALB_NIR_dif),             __FILE__, __LINE__, &
                      mask = LANDUSE_exists_land(:,:)                              )
       call VALCHECK( LIA, LIS, LIE, LJA, LJS, LJE,                                &
                      LAND_SFC_albedo(:,:,I_R_direct ,I_R_VIS), 0.0_RP,    2.0_RP, &
                      VAR_NAME(I_SFC_ALB_VIS_dir),             __FILE__, __LINE__, &
                      mask = LANDUSE_exists_land(:,:)                              )
       call VALCHECK( LIA, LIS, LIE, LJA, LJS, LJE,                                &
                      LAND_SFC_albedo(:,:,I_R_diffuse,I_R_VIS), 0.0_RP,    2.0_RP, &
                      VAR_NAME(I_SFC_ALB_VIS_dif),             __FILE__, __LINE__, &
                      mask = LANDUSE_exists_land(:,:)                              )

       if ( SNOW_flag ) then
          call VALCHECK( LIA, LIS, LIE, LJA, LJS, LJE,                      &
                         SNOW_SFC_TEMP(:,:), 0.0_RP, 1000.0_RP,             &
                         VAR_NAME(I_SNOW_SFC_TEMP),     __FILE__, __LINE__, &
                         mask = LANDUSE_exists_land(:,:)                    )
          call VALCHECK( LIA, LIS, LIE, LJA, LJS, LJE,                      &
                         SNOW_SWE     (:,:), 0.0_RP, 1000.0_RP,             &
                         VAR_NAME(I_SNOW_SWE),          __FILE__, __LINE__, &
                         mask = LANDUSE_exists_land(:,:)                    )
          call VALCHECK( LIA, LIS, LIE, LJA, LJS, LJE,                      &
                         SNOW_Depth   (:,:), 0.0_RP, 1000.0_RP,             &
                         VAR_NAME(I_SNOW_Depth),        __FILE__, __LINE__, &
                         mask = LANDUSE_exists_land(:,:)                    )

          call VALCHECK( LIA, LIS, LIE, LJA, LJS, LJE,                      &
                         SNOW_Dzero   (:,:), 0.0_RP, 1000.0_RP,             &
                         VAR_NAME(I_SNOW_Dzero),        __FILE__, __LINE__, &
                         mask = LANDUSE_exists_land(:,:)                    )
       endif

    end if

    if ( present(force) ) then
       check = force
    else
       check = STATISTICS_checktotal
    end if

    if ( check ) then

       ! 3D
       call STATISTICS_total( LKA, LKS, LKE, LIA, LIS, LIE, LJA, LJS, LJE, &
                              LAND_TEMP (:,:,:), VAR_NAME(I_TEMP),  & ! (in)
                              LAND_GRID_CARTESC_REAL_VOL(:,:,:),    & ! (in)
                              LAND_GRID_CARTESC_REAL_TOTVOL         ) ! (in)
       call STATISTICS_total( LKA, LKS, LKE, LIA, LIS, LIE, LJA, LJS, LJE, &
                              LAND_WATER(:,:,:), VAR_NAME(I_WATER), & ! (in)
                              LAND_GRID_CARTESC_REAL_VOL(:,:,:),    & ! (in)
                              LAND_GRID_CARTESC_REAL_TOTVOL         ) ! (in)
       call STATISTICS_total( LKA, LKS, LKE, LIA, LIS, LIE, LJA, LJS, LJE, &
                              LAND_ICE(:,:,:), VAR_NAME(I_ICE),     & ! (in)
                              LAND_GRID_CARTESC_REAL_VOL(:,:,:),    & ! (in)
                              LAND_GRID_CARTESC_REAL_TOTVOL         ) ! (in)

       ! 2D
       call STATISTICS_total( LIA, LIS, LIE, LJA, LJS, LJE,                                          & ! [IN]
                              LAND_SFC_TEMP  (:,:),                     VAR_NAME(I_SFC_TEMP),        & ! [IN]
                              LAND_GRID_CARTESC_REAL_AREA(:,:),                                      & ! [IN]
                              LAND_GRID_CARTESC_REAL_TOTAREA                                         ) ! [IN]
       call STATISTICS_total( LIA, LIS, LIE, LJA, LJS, LJE,                                          & ! [IN]
                              LAND_SFC_albedo(:,:,I_R_direct ,I_R_IR ), VAR_NAME(I_SFC_ALB_IR_dir),  & ! [IN]
                              LAND_GRID_CARTESC_REAL_AREA(:,:),                                      & ! [IN]
                              LAND_GRID_CARTESC_REAL_TOTAREA                                         ) ! [IN]
       call STATISTICS_total( LIA, LIS, LIE, LJA, LJS, LJE,                                          & ! [IN]
                              LAND_SFC_albedo(:,:,I_R_diffuse,I_R_IR ), VAR_NAME(I_SFC_ALB_IR_dif),  & ! [IN]
                              LAND_GRID_CARTESC_REAL_AREA(:,:),                                      & ! [IN]
                              LAND_GRID_CARTESC_REAL_TOTAREA                                         ) ! [IN]
       call STATISTICS_total( LIA, LIS, LIE, LJA, LJS, LJE,                                          & ! [IN]
                              LAND_SFC_albedo(:,:,I_R_direct ,I_R_NIR), VAR_NAME(I_SFC_ALB_NIR_dir), & ! [IN]
                              LAND_GRID_CARTESC_REAL_AREA(:,:),                                      & ! [IN]
                              LAND_GRID_CARTESC_REAL_TOTAREA                                         ) ! [IN]
       call STATISTICS_total( LIA, LIS, LIE, LJA, LJS, LJE,                                          & ! [IN]
                              LAND_SFC_albedo(:,:,I_R_diffuse,I_R_NIR), VAR_NAME(I_SFC_ALB_NIR_dif), & ! [IN]
                              LAND_GRID_CARTESC_REAL_AREA(:,:),                                      & ! [IN]
                              LAND_GRID_CARTESC_REAL_TOTAREA                                         ) ! [IN]
       call STATISTICS_total( LIA, LIS, LIE, LJA, LJS, LJE,                                          & ! [IN]
                              LAND_SFC_albedo(:,:,I_R_direct ,I_R_VIS), VAR_NAME(I_SFC_ALB_VIS_dir), & ! [IN]
                              LAND_GRID_CARTESC_REAL_AREA(:,:),                                      & ! [IN]
                              LAND_GRID_CARTESC_REAL_TOTAREA                                         ) ! [IN]
       call STATISTICS_total( LIA, LIS, LIE, LJA, LJS, LJE,                                          & ! [IN]
                              LAND_SFC_albedo(:,:,I_R_diffuse,I_R_VIS), VAR_NAME(I_SFC_ALB_VIS_dif), & ! [IN]
                              LAND_GRID_CARTESC_REAL_AREA(:,:),                                      & ! [IN]
                              LAND_GRID_CARTESC_REAL_TOTAREA                                         ) ! [IN]

       if ( SNOW_flag ) then
          call STATISTICS_total( LIA, LIS, LIE, LJA, LJS, LJE, &
                                 SNOW_SFC_TEMP (:,:), VAR_NAME(I_SNOW_SFC_TEMP),  & ! (in)
                                 LAND_GRID_CARTESC_REAL_AREA(:,:),                & ! (in)
                                 LAND_GRID_CARTESC_REAL_TOTAREA                   ) ! (in)
          call STATISTICS_total( LIA, LIS, LIE, LJA, LJS, LJE, &
                                 SNOW_SWE      (:,:), VAR_NAME(I_SNOW_SWE),       & ! (in)
                                 LAND_GRID_CARTESC_REAL_AREA(:,:),                & ! (in)
                                 LAND_GRID_CARTESC_REAL_TOTAREA                   ) ! (in)
          call STATISTICS_total( LIA, LIS, LIE, LJA, LJS, LJE, &
                                 SNOW_Depth    (:,:), VAR_NAME(I_SNOW_Depth),     & ! (in)
                                 LAND_GRID_CARTESC_REAL_AREA(:,:),                & ! (in)
                                 LAND_GRID_CARTESC_REAL_TOTAREA                   ) ! (in)
          call STATISTICS_total( LIA, LIS, LIE, LJA, LJS, LJE, &
                                 SNOW_Dzero    (:,:), VAR_NAME(I_SNOW_Dzero),     & ! (in)
                                 LAND_GRID_CARTESC_REAL_AREA(:,:),                & ! (in)
                                 LAND_GRID_CARTESC_REAL_TOTAREA                   ) ! (in)
          call STATISTICS_total( LIA, LIS, LIE, LJA, LJS, LJE, &
                                 SNOW_nosnowsec(:,:), VAR_NAME(I_SNOW_nosnowsec), & ! (in)
                                 LAND_GRID_CARTESC_REAL_AREA(:,:),                & ! (in)
                                 LAND_GRID_CARTESC_REAL_TOTAREA                   ) ! (in)
       end if

    endif

    return
  end subroutine LAND_vars_check

  !-----------------------------------------------------------------------------
  !> Budget monitor for land
  subroutine LAND_param_read
    use scale_prc, only: &
       PRC_abort
    use scale_landuse, only: &
       LANDUSE_PFT_nmin, &
       LANDUSE_PFT_nmax
    implicit none

    integer              :: index
    character(len=H_MID) :: description
    real(RP)             :: STRGMAX ! Water Limit             [0-1]
    real(RP)             :: STRGCRT ! Water Critical          [0-1]
    real(RP)             :: RSTOMA  ! Stomata Resistance      [0-1]
    real(RP)             :: TCS     ! Thermal Conductivity    [W m-1 K-1]
    real(RP)             :: HCS     ! Dencity x Heat Capacity [J m-3 K-1]
    real(RP)             :: DFW     ! Water Diffusivity       [m2 s-1]
    real(RP)             :: ALBLW   ! Albedo Long Wave        [0-1]
    real(RP)             :: ALBSW   ! Albedo Short Wave       [0-1]
    real(RP)             :: Z0M     ! Z0 for momentum         [m]
    real(RP)             :: Z0H     ! Z0 for heat             [m]
    real(RP)             :: Z0E     ! Z0 for vapor            [m]

    namelist / PARAM_LAND_PROPERTY / &
       LAND_PROPERTY_IN_FILENAME

    namelist / PARAM_LAND_DATA / &
       index,       &
       description, &
       STRGMAX,     &
       STRGCRT,     &
       RSTOMA,      &
       TCS,         &
       HCS,         &
       DFW,         &
       ALBLW,       &
       ALBSW,       &
       Z0M,         &
       Z0H,         &
       Z0E

    integer :: n
    integer :: ierr

    character(len=H_LONG) :: fname

    integer :: IO_FID_LAND_PROPERTY
    !---------------------------------------------------------------------------

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_LAND_PROPERTY,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("LAND_param_read",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("LAND_param_read",*) 'Not appropriate names in namelist PARAM_LAND_PROPERTY. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_LAND_PROPERTY)

    if( LAND_PROPERTY_IN_FILENAME /= '' ) then
      !--- Open land parameter file
      IO_FID_LAND_PROPERTY = IO_get_available_fid()
      call IO_get_fname(fname, LAND_PROPERTY_IN_FILENAME)
      open( IO_FID_LAND_PROPERTY, &
            file   = fname,       &
            form   = 'formatted', &
            status = 'old',       &
            iostat = ierr         )

      if ( ierr /= 0 ) then
         LOG_ERROR("LAND_param_read",*) 'Failed to open land parameter file! :', trim(fname)
         call PRC_abort
      else
        LOG_NEWLINE
        LOG_INFO("LAND_param_read",*) 'Properties for each plant functional type (PFT)'
        LOG_INFO_CONT('(12(1x,A))') '                         PFT DESCRIPTION', &
                                    'Max Stg', &
                                    'CRT Stg', &
                                    'Stm.Res', &
                                    'T condu', &
                                    'H capac', &
                                    'DFC Wat', &
                                    'LW ALB',  &
                                    'SW ALB',  &
                                    ' Z0(m)',  &
                                    ' Z0(h)',  &
                                    ' Z0(e)'

        !--- read namelist
        rewind(IO_FID_LAND_PROPERTY)

        do n = LANDUSE_PFT_nmin, LANDUSE_PFT_nmax
           ! default value
           ALBSW   =  0.2_RP
           STRGMAX =  0.2_RP
           STRGCRT =  0.1_RP
           RSTOMA  = 50.0_RP
           TCS     =  1.0_RP
           HCS     =  2.E+6_RP
           DFW     =  1.E-6_RP
           ALBLW   =  0.04_RP
           ALBSW   =  0.22_RP
           Z0M     =  0.1_RP
           Z0H     = -1.0_RP
           Z0E     = -1.0_RP

           read(IO_FID_LAND_PROPERTY,nml=PARAM_LAND_DATA,iostat=ierr)
           if ( ierr < 0 ) then !--- no more data
              exit
           elseif( ierr > 0 ) then !--- fatal error
              LOG_ERROR("LAND_param_read",*) 'Not appropriate names in namelist PARAM_LAND_DATA. Check!'
              call PRC_abort
           endif

           if( Z0H < 0.0_RP ) then
             Z0H = Z0M / 7.4_RP ! defined by Garratt and Francey (1978)
           endif
           if( Z0E < 0.0_RP ) then
             Z0E = Z0M / 7.4_RP ! defined by Garratt and Francey (1978)
           endif

           LAND_PROPERTY_table(index,I_WaterLimit   ) = STRGMAX
           LAND_PROPERTY_table(index,I_WaterCritical) = STRGCRT
           LAND_PROPERTY_table(index,I_StomataResist) = RSTOMA
           LAND_PROPERTY_table(index,I_ThermalCond  ) = TCS
           LAND_PROPERTY_table(index,I_HeatCapacity ) = HCS
           LAND_PROPERTY_table(index,I_WaterDiff    ) = DFW
           LAND_PROPERTY_table(index,I_ALBLW        ) = ALBLW
           LAND_PROPERTY_table(index,I_ALBSW        ) = ALBSW
           LAND_PROPERTY_table(index,I_Z0M          ) = Z0M
           LAND_PROPERTY_table(index,I_Z0H          ) = Z0H
           LAND_PROPERTY_table(index,I_Z0E          ) = Z0E

           LOG_INFO_CONT('(1x,A4,I4.3,1x,A32,4(1x,F7.3),2(1x,ES7.1),5(1x,F6.3))') &
                                         'IDX=', index, &
                                         trim(description), &
                                         STRGMAX, &
                                         STRGCRT, &
                                         RSTOMA,  &
                                         TCS,     &
                                         HCS,     &
                                         DFW,     &
                                         ALBLW,   &
                                         ALBSW,   &
                                         Z0M,     &
                                         Z0H,     &
                                         Z0E
        enddo

      end if

      close( IO_FID_LAND_PROPERTY )

    endif

    return
  end subroutine LAND_param_read

  !-----------------------------------------------------------------------------
  !> conversion from water saturation [fraction] to volumetric water content [m3/m3]
  function convert_WS2VWC( WS, critical ) result( VWC )
    implicit none

    real(RP), intent(in) :: WS(LIA,LJA) ! water saturation [fraction]
    logical,  intent(in) :: critical  ! is I_WaterCritical used?

    real(RP) :: VWC(LIA,LJA) ! volumetric water content [m3/m3]

    ! work
    integer :: i, j, num
    !---------------------------------------------------------------------------

    if( critical ) then
      num = I_WaterCritical
    else
      num = I_WaterLimit
    end if

    do j = LJS, LJE
    do i = LIS, LIE
      VWC(i,j) = max( min( WS(i,j)*LAND_PROPERTY(i,j,num), LAND_PROPERTY(i,j,num) ), 0.0_RP )
    end do
    end do

    return
  end function convert_WS2VWC

  !-----------------------------------------------------------------------------
  !> Create land restart file
  subroutine LAND_vars_restart_create
    use scale_time, only: &
       TIME_gettimelabel
    use scale_file_cartesC, only: &
       FILE_CARTESC_create
    use mod_land_admin, only: &
       LAND_do
    implicit none

    character(len=19)     :: timelabel
    character(len=H_LONG) :: basename
    !---------------------------------------------------------------------------

    call PROF_rapstart('LND_Restart', 1)

    if ( LAND_do .and. LAND_RESTART_OUT_BASENAME /= '' ) then

       LOG_NEWLINE
       LOG_INFO("LAND_vars_restart_create",*) 'Create restart file (LAND) '

       if ( LAND_RESTART_OUT_POSTFIX_TIMELABEL ) then
          call TIME_gettimelabel( timelabel )
          basename = trim(LAND_RESTART_OUT_BASENAME)//'_'//trim(timelabel)
       else
          basename = trim(LAND_RESTART_OUT_BASENAME)
       endif

       LOG_INFO("LAND_vars_restart_create",*) 'basename: ', trim(basename)

       call FILE_CARTESC_create( &
            basename, LAND_RESTART_OUT_TITLE, LAND_RESTART_OUT_DTYPE, & ! [IN]
            restart_fid,                                              & ! [OUT]
            aggregate=LAND_RESTART_OUT_AGGREGATE                      ) ! [IN]

    endif

    call PROF_rapend('LND_Restart', 1)

    return
  end subroutine LAND_vars_restart_create

  !-----------------------------------------------------------------------------
  !> Exit netCDF define mode
  subroutine LAND_vars_restart_enddef
    use scale_file_cartesC, only: &
       FILE_CARTESC_enddef
    implicit none

    call PROF_rapstart('LND_Restart', 1)

    if ( restart_fid /= -1 ) then
       call FILE_CARTESC_enddef( restart_fid ) ! [IN]
    endif

    call PROF_rapend('LND_Restart', 1)

    return
  end subroutine LAND_vars_restart_enddef

  !-----------------------------------------------------------------------------
  !> Close restart file
  subroutine LAND_vars_restart_close
    use scale_file_cartesC, only: &
       FILE_CARTESC_close
    implicit none
    !---------------------------------------------------------------------------

    call PROF_rapstart('LND_Restart', 1)

    if ( restart_fid /= -1 ) then
       LOG_NEWLINE
       LOG_INFO("LAND_vars_restart_close",*) 'Close restart file (LAND) '

       call FILE_CARTESC_close( restart_fid ) ! [IN]

       restart_fid = -1
    endif

    call PROF_rapend('LND_Restart', 1)

    return
  end subroutine LAND_vars_restart_close

  !-----------------------------------------------------------------------------
  !> Define land variables in restart file
  subroutine LAND_vars_restart_def_var
    use scale_file_cartesC, only: &
       FILE_CARTESC_def_var
    implicit none
    integer :: i
    !---------------------------------------------------------------------------

    call PROF_rapstart('LND_Restart', 1)

    if ( restart_fid /= -1 ) then

       do i = I_TEMP, I_ICE
          call FILE_CARTESC_def_var( restart_fid,     & ! [IN]
               VAR_NAME(i), VAR_DESC(i), VAR_UNIT(i), & ! [IN]
               'LXY', LAND_RESTART_OUT_DTYPE,         & ! [IN]
               VAR_ID(i),                             & ! [OUT]
               standard_name=VAR_STDN(i)              ) ! [IN]
       end do
       do i = I_SFC_TEMP, I_SFC_ALB_VIS_dif
          call FILE_CARTESC_def_var( restart_fid,     & ! [IN]
               VAR_NAME(i), VAR_DESC(i), VAR_UNIT(i), & ! [IN]
               'XY', LAND_RESTART_OUT_DTYPE,          & ! [IN]
               VAR_ID(i),                             & ! [OUT]
               standard_name=VAR_STDN(i)              ) ! [IN]
       end do

       if ( SNOW_flag ) then
          do i = I_SNOW_SFC_TEMP, I_SNOW_nosnowsec
             call FILE_CARTESC_def_var( restart_fid,     & ! [IN]
                  VAR_NAME(i), VAR_DESC(i), VAR_UNIT(i), & ! [IN]
                  'XY', LAND_RESTART_OUT_DTYPE,          & ! [IN]
                  VAR_ID(i),                             & ! [OUT]
                  standard_name=VAR_STDN(i)              ) ! [IN]
          end do
       end if

    endif

    call PROF_rapend('LND_Restart', 1)

    return
  end subroutine LAND_vars_restart_def_var

  !-----------------------------------------------------------------------------
  !> Write land variables to restart file
  subroutine LAND_vars_restart_write
    use scale_file_cartesC, only: &
       FILE_CARTESC_write_var
    implicit none
    !---------------------------------------------------------------------------

    call PROF_rapstart('LND_Restart', 1)

    if ( restart_fid /= -1 ) then

       call LAND_vars_check( force = .true. )

       call FILE_CARTESC_write_var( restart_fid, VAR_ID(I_TEMP),                         & ! [IN]
                                    LAND_TEMP(:,:,:),                                    & ! [IN]
                                    VAR_NAME(I_TEMP),            'LXY', fill_halo=.true. ) ! [IN]
       call FILE_CARTESC_write_var( restart_fid, VAR_ID(I_WATER),                        & ! [IN]
                                    LAND_WATER(:,:,:),                                   & ! [IN]
                                    VAR_NAME(I_WATER),           'LXY', fill_halo=.true. ) ! [IN]
       call FILE_CARTESC_write_var( restart_fid, VAR_ID(I_ICE),                          & ! [IN]
                                    LAND_ICE(:,:,:),                                     & ! [IN]
                                    VAR_NAME(I_ICE),             'LXY', fill_halo=.true. ) ! [IN]
       call FILE_CARTESC_write_var( restart_fid, VAR_ID(I_SFC_TEMP),                     & ! [IN]
                                    LAND_SFC_TEMP(:,:),                                  & ! [IN]
                                    VAR_NAME(I_SFC_TEMP),        'XY',  fill_halo=.true. ) ! [IN]
       call FILE_CARTESC_write_var( restart_fid, VAR_ID(I_SFC_ALB_IR_dir),               & ! [IN]
                                    LAND_SFC_albedo(:,:,I_R_direct ,I_R_IR ),            & ! [IN]
                                    VAR_NAME(I_SFC_ALB_IR_dir),  'XY',  fill_halo=.true. ) ! [IN]
       call FILE_CARTESC_write_var( restart_fid, VAR_ID(I_SFC_ALB_IR_dif),               & ! [IN]
                                    LAND_SFC_albedo(:,:,I_R_diffuse,I_R_IR ),            & ! [IN]
                                    VAR_NAME(I_SFC_ALB_IR_dif),  'XY',  fill_halo=.true. ) ! [IN]
       call FILE_CARTESC_write_var( restart_fid, VAR_ID(I_SFC_ALB_NIR_dir),              & ! [IN]
                                    LAND_SFC_albedo(:,:,I_R_direct ,I_R_NIR),            & ! [IN]
                                    VAR_NAME(I_SFC_ALB_NIR_dir), 'XY',  fill_halo=.true. ) ! [IN]
       call FILE_CARTESC_write_var( restart_fid, VAR_ID(I_SFC_ALB_NIR_dif),              & ! [IN]
                                    LAND_SFC_albedo(:,:,I_R_diffuse,I_R_NIR),            & ! [IN]
                                    VAR_NAME(I_SFC_ALB_NIR_dif), 'XY',  fill_halo=.true. ) ! [IN]
       call FILE_CARTESC_write_var( restart_fid, VAR_ID(I_SFC_ALB_VIS_dir),              & ! [IN]
                                    LAND_SFC_albedo(:,:,I_R_direct ,I_R_VIS),            & ! [IN]
                                    VAR_NAME(I_SFC_ALB_VIS_dir), 'XY',  fill_halo=.true. ) ! [IN]
       call FILE_CARTESC_write_var( restart_fid, VAR_ID(I_SFC_ALB_VIS_dif),              & ! [IN]
                                    LAND_SFC_albedo(:,:,I_R_diffuse,I_R_VIS),            & ! [IN]
                                    VAR_NAME(I_SFC_ALB_VIS_dif), 'XY',  fill_halo=.true. ) ! [IN]

       if ( SNOW_flag ) then
          call FILE_CARTESC_write_var( restart_fid, VAR_ID(I_SNOW_SFC_TEMP), SNOW_SFC_TEMP(:,:),   &
                                       VAR_NAME(I_SNOW_SFC_TEMP),  'XY', fill_halo=.true.          )
          call FILE_CARTESC_write_var( restart_fid, VAR_ID(I_SNOW_SWE), SNOW_SWE(:,:),             &
                                       VAR_NAME(I_SNOW_SWE),       'XY', fill_halo=.true.          )
          call FILE_CARTESC_write_var( restart_fid, VAR_ID(I_SNOW_Depth), SNOW_Depth(:,:),         &
                                       VAR_NAME(I_SNOW_Depth),     'XY', fill_halo=.true.          )
          call FILE_CARTESC_write_var( restart_fid, VAR_ID(I_SNOW_Dzero), SNOW_Dzero(:,:),         &
                                       VAR_NAME(I_SNOW_Dzero),     'XY', fill_halo=.true.          )
          call FILE_CARTESC_write_var( restart_fid, VAR_ID(I_SNOW_nosnowsec), SNOW_nosnowsec(:,:), &
                                       VAR_NAME(I_SNOW_nosnowsec), 'XY', fill_halo=.true.          )
       end if

    endif

    call PROF_rapend('LND_Restart', 1)

    return
  end subroutine LAND_vars_restart_write

end module mod_land_vars
