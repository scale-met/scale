!-------------------------------------------------------------------------------
!> module OCEAN Variables
!!
!! @par Description
!!          Container for ocean variables
!!
!! @author Team SCALE
!!
!<
!-------------------------------------------------------------------------------
#include "scalelib.h"
module mod_ocean_vars
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_io
  use scale_prof
  use scale_debug
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
  public :: OCEAN_vars_setup
  public :: OCEAN_vars_restart_read
  public :: OCEAN_vars_restart_write
  public :: OCEAN_vars_history
  public :: OCEAN_vars_total

  public :: OCEAN_vars_restart_create
  public :: OCEAN_vars_restart_open
  public :: OCEAN_vars_restart_def_var
  public :: OCEAN_vars_restart_enddef
  public :: OCEAN_vars_restart_close

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  logical,                public :: OCEAN_RESTART_OUTPUT                = .false.         !< Output restart file?

  character(len=H_LONG),  public :: OCEAN_RESTART_IN_BASENAME           = ''              !< Basename of the input  file
  logical,                public :: OCEAN_RESTART_IN_AGGREGATE                            !< Switch to use aggregate file
  logical,                public :: OCEAN_RESTART_IN_POSTFIX_TIMELABEL  = .false.         !< Add timelabel to the basename of input  file?
  character(len=H_LONG),  public :: OCEAN_RESTART_OUT_BASENAME          = ''              !< Basename of the output file
  logical,                public :: OCEAN_RESTART_OUT_AGGREGATE                           !< Switch to use aggregate file
  logical,                public :: OCEAN_RESTART_OUT_POSTFIX_TIMELABEL = .true.          !< Add timelabel to the basename of output file?
  character(len=H_MID),   public :: OCEAN_RESTART_OUT_TITLE             = 'OCEAN restart' !< Title    of the output file
  character(len=H_SHORT), public :: OCEAN_RESTART_OUT_DTYPE             = 'DEFAULT'       !< REAL4 or REAL8

  ! prognostic variables
  real(RP), public, allocatable :: OCEAN_TEMP       (:,:,:)   !< ocean temperature         [K]
  real(RP), public, allocatable :: OCEAN_SALT       (:,:,:)   !< ocean salinity            [PSU]
  real(RP), public, allocatable :: OCEAN_UVEL       (:,:,:)   !< ocean zonal velocity      [m/s]
  real(RP), public, allocatable :: OCEAN_VVEL       (:,:,:)   !< ocean meridional velocity [m/s]

  real(RP), public, allocatable :: OCEAN_OCN_Z0M    (:,:)     !< surface roughness length for momentum, open ocean [m]

  real(RP), public, allocatable :: OCEAN_ICE_TEMP   (:,:)     !< sea ice temperature [K]
  real(RP), public, allocatable :: OCEAN_ICE_MASS   (:,:)     !< sea ice mass        [kg]

  real(RP), public, allocatable :: OCEAN_SFC_TEMP   (:,:)     !< ocean surface skin temperature [K]
  real(RP), public, allocatable :: OCEAN_SFC_albedo (:,:,:,:) !< ocean surface albedo (direct/diffuse,IR/near-IR/VIS) (0-1)
  real(RP), public, allocatable :: OCEAN_SFC_Z0M    (:,:)     !< ocean surface roughness length for momentum [m]
  real(RP), public, allocatable :: OCEAN_SFC_Z0H    (:,:)     !< ocean surface roughness length for heat     [m]
  real(RP), public, allocatable :: OCEAN_SFC_Z0E    (:,:)     !< ocean surface roughness length for vapor    [m]

  ! tendency variables
  real(RP), public, allocatable :: OCEAN_TEMP_t     (:,:,:)   !< tendency of OCEAN_OCN_TEMP
  real(RP), public, allocatable :: OCEAN_SALT_t     (:,:,:)   !< tendency of OCEAN_OCN_SALT
  real(RP), public, allocatable :: OCEAN_UVEL_t     (:,:,:)   !< tendency of OCEAN_OCN_UVEL
  real(RP), public, allocatable :: OCEAN_VVEL_t     (:,:,:)   !< tendency of OCEAN_OCN_VVEL

  real(RP), public, allocatable :: OCEAN_ICE_TEMP_t (:,:)     !< tendency of OCEAN_ICE_TEMP
  real(RP), public, allocatable :: OCEAN_ICE_MASS_t (:,:)     !< tendency of OCEAN_ICE_MASS

  ! recieved from lowermost atmosphere
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
  real(RP), public, allocatable :: ATMOS_SFLX_rain  (:,:)
  real(RP), public, allocatable :: ATMOS_SFLX_snow  (:,:)

  ! send to lowermost atmosphere
  real(RP), public, allocatable :: OCEAN_SFLX_MW    (:,:)     !< ocean surface w-momentum flux    [kg/m/s2]
  real(RP), public, allocatable :: OCEAN_SFLX_MU    (:,:)     !< ocean surface u-momentum flux    [kg/m/s2]
  real(RP), public, allocatable :: OCEAN_SFLX_MV    (:,:)     !< ocean surface v-momentum flux    [kg/m/s2]
  real(RP), public, allocatable :: OCEAN_SFLX_SH    (:,:)     !< ocean surface sensible heat flux [J/m2/s]
  real(RP), public, allocatable :: OCEAN_SFLX_LH    (:,:)     !< ocean surface latent heat flux   [J/m2/s]
  real(RP), public, allocatable :: OCEAN_SFLX_QTRC  (:,:,:)   !< ocean surface tracer flux        [kg/m2/s]
  real(RP), public, allocatable :: OCEAN_U10        (:,:)     !< ocean surface velocity u at 10m  [m/s]
  real(RP), public, allocatable :: OCEAN_V10        (:,:)     !< ocean surface velocity v at 10m  [m/s]
  real(RP), public, allocatable :: OCEAN_T2         (:,:)     !< ocean surface temperature at 2m  [K]
  real(RP), public, allocatable :: OCEAN_Q2         (:,:)     !< ocean surface water vapor at 2m  [kg/kg]

  ! send to uppermost ocean
  real(RP), public, allocatable :: OCEAN_SFLX_G     (:,:)     !< ocean surface water heat flux   [J/m2/s]
  real(RP), public, allocatable :: OCEAN_SFLX_water (:,:)     !< ocean surface liquid water flux [kg/m2/s]
  real(RP), public, allocatable :: OCEAN_SFLX_ice   (:,:)     !< ocean surface ice water flux    [kg/m2/s]

  ! diagnostic value
  real(RP), public, allocatable :: OCEAN_ICE_FRAC   (:,:)     !< area fraction of sea ice [1]

  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  logical,                private :: OCEAN_VARS_CHECKRANGE              = .false.
  logical,                private :: OCEAN_RESTART_IN_CHECK_COORDINATES = .true.

  integer,                private, parameter :: VMAX              = 27 !< number of the variables
  integer,                private, parameter :: I_TEMP            =  1
  integer,                private, parameter :: I_SALT            =  2
  integer,                private, parameter :: I_UVEL            =  3
  integer,                private, parameter :: I_VVEL            =  4
  integer,                private, parameter :: I_OCN_Z0M         =  5
  integer,                private, parameter :: I_ICE_TEMP        =  6
  integer,                private, parameter :: I_ICE_MASS        =  7
  integer,                private, parameter :: I_SFC_TEMP        =  8
  integer,                private, parameter :: I_SFC_ALB_IR_dir  =  9
  integer,                private, parameter :: I_SFC_ALB_IR_dif  = 10
  integer,                private, parameter :: I_SFC_ALB_NIR_dir = 11
  integer,                private, parameter :: I_SFC_ALB_NIR_dif = 12
  integer,                private, parameter :: I_SFC_ALB_VIS_dir = 13
  integer,                private, parameter :: I_SFC_ALB_VIS_dif = 14
  integer,                private, parameter :: I_SFC_Z0M         = 15
  integer,                private, parameter :: I_SFC_Z0H         = 16
  integer,                private, parameter :: I_SFC_Z0E         = 17
  integer,                private, parameter :: I_SFLX_MW         = 18
  integer,                private, parameter :: I_SFLX_MU         = 19
  integer,                private, parameter :: I_SFLX_MV         = 20
  integer,                private, parameter :: I_SFLX_SH         = 21
  integer,                private, parameter :: I_SFLX_LH         = 22
  integer,                private, parameter :: I_SFLX_evap       = 23
  integer,                private, parameter :: I_SFLX_G          = 24
  integer,                private, parameter :: I_SFLX_water      = 25
  integer,                private, parameter :: I_SFLX_ice        = 26
  integer,                private, parameter :: I_ICE_FRAC        = 27

  character(len=H_SHORT), private            :: VAR_NAME(VMAX)    !< name  of the variables
  character(len=H_MID),   private            :: VAR_DESC(VMAX)    !< desc. of the variables
  character(len=H_MID),   private            :: VAR_STDN(VMAX)    !< standard name of the variables
  character(len=H_SHORT), private            :: VAR_UNIT(VMAX)    !< unit  of the variables
  integer,                private            :: VAR_ID(VMAX)      !< ID    of the variables
  integer,                private            :: restart_fid = -1  ! file ID

  data VAR_NAME / 'OCEAN_TEMP',            &
                  'OCEAN_SALT',            &
                  'OCEAN_UVEL',            &
                  'OCEAN_VVEL',            &
                  'OCEAN_OCN_Z0M',         &
                  'OCEAN_ICE_TEMP',        &
                  'OCEAN_ICE_MASS',        &
                  'OCEAN_SFC_TEMP',        &
                  'OCEAN_SFC_ALB_IR_dir',  &
                  'OCEAN_SFC_ALB_IR_dif',  &
                  'OCEAN_SFC_ALB_NIR_dir', &
                  'OCEAN_SFC_ALB_NIR_dif', &
                  'OCEAN_SFC_ALB_VIS_dir', &
                  'OCEAN_SFC_ALB_VIS_dif', &
                  'OCEAN_SFC_Z0M',         &
                  'OCEAN_SFC_Z0H',         &
                  'OCEAN_SFC_Z0E',         &
                  'OCEAN_SFLX_MW',         &
                  'OCEAN_SFLX_MU',         &
                  'OCEAN_SFLX_MV',         &
                  'OCEAN_SFLX_SH',         &
                  'OCEAN_SFLX_LH',         &
                  'OCEAN_SFLX_evap',       &
                  'OCEAN_SFLX_G',          &
                  'OCEAN_SFLX_water',      &
                  'OCEAN_SFLX_ice',        &
                  'OCEAN_ICE_FRAC'         /

  data VAR_DESC / 'ocean temperature',                          &
                  'ocean salinity',                             &
                  'ocean u-velocity',                           &
                  'ocean v-velocity',                           &
                  'open ocean roughness length (momentum)',     &
                  'seaice temperature',                         &
                  'seaice mass',                                &
                  'ocean surface skin temperature',             &
                  'ocean surface albedo for IR (direct)',       &
                  'ocean surface albedo for IR (diffuse)',      &
                  'ocean surface albedo for NIR (direct)',      &
                  'ocean surface albedo for NIR (diffuse)',     &
                  'ocean surface albedo for VIS (direct)',      &
                  'ocean surface albedo for VIS (diffuse)',     &
                  'ocean surface roughness length (momentum)',  &
                  'ocean surface roughness length (heat)',      &
                  'ocean surface roughness length (vapor)',     &
                  'ocean surface w-momentum flux (upward)',     &
                  'ocean surface u-momentum flux (upward)',     &
                  'ocean surface v-momentum flux (upward)',     &
                  'ocean surface sensible heat flux (upward)',  &
                  'ocean surface latent heat flux (upward)',    &
                  'ocean surface water vapor flux (upward)',    &
                  'ocean subsurface heat flux (downward)',      &
                  'ocean surface liquid water flux (downward)', &
                  'ocean surface ice water flux (downward)',    &
                  'seaice fraction'                             /

  data VAR_STDN / 'sea_water_temperature',        &
                  'sea_water_salinity',           &
                  'eastward_sea_water_velocity',  &
                  'northward_sea_water_velocity', &
                  '', &
                  '', &
                  '', &
                  'sea_surface_skin_temperature', &
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
                  '', &
                  '', &
                  '', &
                  '', &
                  '', &
                  '', &
                  '', &
                  '', &
                  '' /

  data VAR_UNIT / 'K',       &
                  'PSU',     &
                  'm/s',     &
                  'm/s',     &
                  'm',       &
                  'K',       &
                  'kg/m2',   &
                  'K',       &
                  '1',       &
                  '1',       &
                  '1',       &
                  '1',       &
                  '1',       &
                  '1',       &
                  'm',       &
                  'm',       &
                  'm',       &
                  'kg/m2/s', &
                  'kg/m2/s', &
                  'kg/m2/s', &
                  'J/m2/s',  &
                  'J/m2/s',  &
                  'kg/m2/s', &
                  'J/m2/s',  &
                  'kg/m2/s', &
                  'kg/m2/s', &
                  '1'        /

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine OCEAN_vars_setup
    use scale_prc, only: &
       PRC_abort
    use scale_const, only: &
       UNDEF => CONST_UNDEF
    implicit none

    namelist / PARAM_OCEAN_VARS / &
       OCEAN_RESTART_IN_BASENAME,           &
       OCEAN_RESTART_IN_AGGREGATE,          &
       OCEAN_RESTART_IN_POSTFIX_TIMELABEL,  &
       OCEAN_RESTART_IN_CHECK_COORDINATES,  &
       OCEAN_RESTART_OUTPUT,                &
       OCEAN_RESTART_OUT_BASENAME,          &
       OCEAN_RESTART_OUT_AGGREGATE,         &
       OCEAN_RESTART_OUT_POSTFIX_TIMELABEL, &
       OCEAN_RESTART_OUT_TITLE,             &
       OCEAN_RESTART_OUT_DTYPE,             &
       OCEAN_VARS_CHECKRANGE

    integer :: ierr
    integer :: iv
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("OCEAN_vars_setup",*) 'Setup'

    allocate( OCEAN_TEMP       (OKMAX,OIA,OJA) )
    allocate( OCEAN_SALT       (OKMAX,OIA,OJA) )
    allocate( OCEAN_UVEL       (OKMAX,OIA,OJA) )
    allocate( OCEAN_VVEL       (OKMAX,OIA,OJA) )
    OCEAN_TEMP       (:,:,:)   = UNDEF
    OCEAN_SALT       (:,:,:)   = UNDEF
    OCEAN_UVEL       (:,:,:)   = UNDEF
    OCEAN_VVEL       (:,:,:)   = UNDEF

    allocate( OCEAN_OCN_Z0M    (OIA,OJA) )
    OCEAN_OCN_Z0M    (:,:)     = UNDEF

    allocate( OCEAN_ICE_TEMP   (OIA,OJA) )
    allocate( OCEAN_ICE_MASS   (OIA,OJA) )
    OCEAN_ICE_TEMP   (:,:)     = UNDEF
    OCEAN_ICE_MASS   (:,:)     = UNDEF

    allocate( OCEAN_SFC_TEMP   (OIA,OJA)                     )
    allocate( OCEAN_SFC_albedo (OIA,OJA,N_RAD_DIR,N_RAD_RGN) )
    allocate( OCEAN_SFC_Z0M    (OIA,OJA)                     )
    allocate( OCEAN_SFC_Z0H    (OIA,OJA)                     )
    allocate( OCEAN_SFC_Z0E    (OIA,OJA)                     )
    OCEAN_SFC_TEMP   (:,:)     = UNDEF
    OCEAN_SFC_albedo (:,:,:,:) = UNDEF
    OCEAN_SFC_Z0M    (:,:)     = UNDEF
    OCEAN_SFC_Z0H    (:,:)     = UNDEF
    OCEAN_SFC_Z0E    (:,:)     = UNDEF

    allocate( OCEAN_TEMP_t     (OKMAX,OIA,OJA) )
    allocate( OCEAN_SALT_t     (OKMAX,OIA,OJA) )
    allocate( OCEAN_UVEL_t     (OKMAX,OIA,OJA) )
    allocate( OCEAN_VVEL_t     (OKMAX,OIA,OJA) )
    OCEAN_TEMP_t     (:,:,:)   = UNDEF
    OCEAN_SALT_t     (:,:,:)   = UNDEF
    OCEAN_UVEL_t     (:,:,:)   = UNDEF
    OCEAN_VVEL_t     (:,:,:)   = UNDEF

    allocate( OCEAN_ICE_TEMP_t (OIA,OJA) )
    allocate( OCEAN_ICE_MASS_t (OIA,OJA) )
    OCEAN_ICE_TEMP_t (:,:)     = UNDEF
    OCEAN_ICE_MASS_t (:,:)     = UNDEF

    allocate( ATMOS_TEMP       (OIA,OJA)                     )
    allocate( ATMOS_PRES       (OIA,OJA)                     )
    allocate( ATMOS_W          (OIA,OJA)                     )
    allocate( ATMOS_U          (OIA,OJA)                     )
    allocate( ATMOS_V          (OIA,OJA)                     )
    allocate( ATMOS_DENS       (OIA,OJA)                     )
    allocate( ATMOS_QV         (OIA,OJA)                     )
    allocate( ATMOS_PBL        (OIA,OJA)                     )
    allocate( ATMOS_SFC_DENS   (OIA,OJA)                     )
    allocate( ATMOS_SFC_PRES   (OIA,OJA)                     )
    allocate( ATMOS_SFLX_rad_dn(OIA,OJA,N_RAD_DIR,N_RAD_RGN) )
    allocate( ATMOS_cosSZA     (OIA,OJA)                     )
    allocate( ATMOS_SFLX_rain  (OIA,OJA)                     )
    allocate( ATMOS_SFLX_snow  (OIA,OJA)                     )
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
    ATMOS_SFLX_rain  (:,:)     = UNDEF
    ATMOS_SFLX_snow  (:,:)     = UNDEF

    allocate( OCEAN_SFLX_MW    (OIA,OJA) )
    allocate( OCEAN_SFLX_MU    (OIA,OJA) )
    allocate( OCEAN_SFLX_MV    (OIA,OJA) )
    allocate( OCEAN_SFLX_SH    (OIA,OJA) )
    allocate( OCEAN_SFLX_LH    (OIA,OJA) )
    allocate( OCEAN_SFLX_QTRC  (OIA,OJA,QA) )
    allocate( OCEAN_U10        (OIA,OJA) )
    allocate( OCEAN_V10        (OIA,OJA) )
    allocate( OCEAN_T2         (OIA,OJA) )
    allocate( OCEAN_Q2         (OIA,OJA) )
    OCEAN_SFLX_MW    (:,:)     = UNDEF
    OCEAN_SFLX_MU    (:,:)     = UNDEF
    OCEAN_SFLX_MV    (:,:)     = UNDEF
    OCEAN_SFLX_SH    (:,:)     = UNDEF
    OCEAN_SFLX_LH    (:,:)     = UNDEF
    OCEAN_SFLX_QTRC  (:,:,:)   = UNDEF
    OCEAN_U10        (:,:)     = UNDEF
    OCEAN_V10        (:,:)     = UNDEF
    OCEAN_T2         (:,:)     = UNDEF
    OCEAN_Q2         (:,:)     = UNDEF

    allocate( OCEAN_SFLX_G     (OIA,OJA) )
    allocate( OCEAN_SFLX_water (OIA,OJA) )
    allocate( OCEAN_SFLX_ice   (OIA,OJA) )
    OCEAN_SFLX_G     (:,:)     = UNDEF
    OCEAN_SFLX_water (:,:)     = UNDEF
    OCEAN_SFLX_ice   (:,:)     = UNDEF

    allocate( OCEAN_ICE_FRAC   (OIA,OJA) )
    OCEAN_ICE_FRAC   (:,:)     = UNDEF

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_OCEAN_VARS,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("OCEAN_vars_setup",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("OCEAN_vars_setup",*) 'Not appropriate names in namelist PARAM_OCEAN_VARS. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_OCEAN_VARS)

    LOG_NEWLINE
    LOG_INFO("OCEAN_vars_setup",*) 'List of prognostic variables (OCEAN) '
    LOG_INFO_CONT('(1x,A,A24,A,A48,A,A12,A)') &
                  '      |','VARNAME                 ','|','DESCRIPTION                                     ','[','UNIT        ',']'
    do iv = 1, VMAX
       LOG_INFO_CONT('(1x,A,I3,A,A24,A,A48,A,A12,A)') &
                     'NO.',iv,'|',VAR_NAME(iv),'|',VAR_DESC(iv),'[',VAR_UNIT(iv),']'
    enddo

    LOG_NEWLINE
    if ( OCEAN_RESTART_IN_BASENAME /= '' ) then
       LOG_INFO("OCEAN_vars_setup",*) 'Restart input?  : YES, file = ', trim(OCEAN_RESTART_IN_BASENAME)
       LOG_INFO_CONT(*)               'Add timelabel?  : ', OCEAN_RESTART_IN_POSTFIX_TIMELABEL
    else
       LOG_INFO("OCEAN_vars_setup",*) 'Restart input?  : NO'
    endif
    if (       OCEAN_RESTART_OUTPUT             &
         .AND. OCEAN_RESTART_OUT_BASENAME /= '' ) then
       LOG_INFO("OCEAN_vars_setup",*) 'Restart output? : YES, file = ', trim(OCEAN_RESTART_OUT_BASENAME)
       LOG_INFO_CONT(*)               'Add timelabel?  : ', OCEAN_RESTART_OUT_POSTFIX_TIMELABEL
    else
       LOG_INFO("OCEAN_vars_setup",*) 'Restart output? : NO'
       OCEAN_RESTART_OUTPUT = .false.
    endif

    return
  end subroutine OCEAN_vars_setup

  !-----------------------------------------------------------------------------
  !> Open ocean restart file for read
  subroutine OCEAN_vars_restart_open
    use scale_time, only: &
       TIME_gettimelabel
    use scale_file_cartesC, only: &
       FILE_CARTESC_open, &
       FILE_CARTESC_check_coordinates
    use mod_ocean_admin, only: &
       OCEAN_do
    implicit none

    character(len=19)     :: timelabel
    character(len=H_LONG) :: basename
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("OCEAN_vars_restart_open",*) 'Open restart file (OCEAN) '

    if ( OCEAN_do .and. OCEAN_RESTART_IN_BASENAME /= '' ) then

       if ( OCEAN_RESTART_IN_POSTFIX_TIMELABEL ) then
          call TIME_gettimelabel( timelabel )
          basename = trim(OCEAN_RESTART_IN_BASENAME)//'_'//trim(timelabel)
       else
          basename = trim(OCEAN_RESTART_IN_BASENAME)
       endif
       LOG_INFO_CONT(*) 'basename: ', trim(basename)

       call FILE_CARTESC_open( basename, restart_fid, aggregate=OCEAN_RESTART_IN_AGGREGATE )

       if( OCEAN_RESTART_IN_CHECK_COORDINATES ) call FILE_CARTESC_check_coordinates( restart_fid )

    else
       LOG_INFO_CONT(*) 'restart file for ocean is not specified.'
    endif

    return
  end subroutine OCEAN_vars_restart_open

  !-----------------------------------------------------------------------------
  !> Read ocean restart
  subroutine OCEAN_vars_restart_read
    use scale_prc, only: &
       PRC_abort
    use scale_file, only: &
       FILE_get_aggregate
    use scale_file_cartesC, only: &
       FILE_CARTESC_read, &
       FILE_CARTESC_flush
    use scale_ocean_phy_ice_simple, only: &
       OCEAN_PHY_ICE_freezetemp, &
       OCEAN_PHY_ICE_fraction
    use mod_ocean_admin, only: &
       OCEAN_ICE_TYPE
    implicit none
    !---------------------------------------------------------------------------

    if ( restart_fid /= -1 ) then
       LOG_NEWLINE
       LOG_INFO("OCEAN_vars_restart_read",*) 'Read from restart file (OCEAN) '

       call FILE_CARTESC_read( restart_fid, VAR_NAME(I_TEMP),            'OXY', & ! [IN]
                               OCEAN_TEMP(:,:,:)                                ) ! [OUT]
!      call FILE_CARTESC_read( restart_fid, VAR_NAME(I_SALT),            'OXY', & ! [IN]
!                              OCEAN_SALT(:,:,:)                                ) ! [OUT]
!      call FILE_CARTESC_read( restart_fid, VAR_NAME(I_UVEL),            'OXY', & ! [IN]
!                              OCEAN_UVEL(:,:,:)                                ) ! [OUT]
!      call FILE_CARTESC_read( restart_fid, VAR_NAME(I_VVEL),            'OXY', & ! [IN]
!                              OCEAN_VVEL(:,:,:)                                ) ! [OUT]
       call FILE_CARTESC_read( restart_fid, VAR_NAME(I_OCN_Z0M),         'XY',  & ! [IN]
                               OCEAN_OCN_Z0M(:,:)                               ) ! [OUT]
       call FILE_CARTESC_read( restart_fid, VAR_NAME(I_ICE_TEMP),        'XY',  & ! [IN]
                               OCEAN_ICE_TEMP(:,:)                              ) ! [OUT]
       call FILE_CARTESC_read( restart_fid, VAR_NAME(I_ICE_MASS),        'XY',  & ! [IN]
                               OCEAN_ICE_MASS(:,:)                              ) ! [OUT]
       call FILE_CARTESC_read( restart_fid, VAR_NAME(I_SFC_TEMP),        'XY',  & ! [IN]
                               OCEAN_SFC_TEMP(:,:)                              ) ! [OUT]
       call FILE_CARTESC_read( restart_fid, VAR_NAME(I_SFC_ALB_IR_dir),  'XY',  & ! [IN]
                               OCEAN_SFC_albedo(:,:,I_R_direct ,I_R_IR )        ) ! [OUT]
       call FILE_CARTESC_read( restart_fid, VAR_NAME(I_SFC_ALB_IR_dif),  'XY',  & ! [IN]
                               OCEAN_SFC_albedo(:,:,I_R_diffuse,I_R_IR )        ) ! [OUT]
       call FILE_CARTESC_read( restart_fid, VAR_NAME(I_SFC_ALB_NIR_dir), 'XY',  & ! [IN]
                               OCEAN_SFC_albedo(:,:,I_R_direct ,I_R_NIR)        ) ! [OUT]
       call FILE_CARTESC_read( restart_fid, VAR_NAME(I_SFC_ALB_NIR_dif), 'XY',  & ! [IN]
                               OCEAN_SFC_albedo(:,:,I_R_diffuse,I_R_NIR)        ) ! [OUT]
       call FILE_CARTESC_read( restart_fid, VAR_NAME(I_SFC_ALB_VIS_dir), 'XY',  & ! [IN]
                               OCEAN_SFC_albedo(:,:,I_R_direct ,I_R_VIS)        ) ! [OUT]
       call FILE_CARTESC_read( restart_fid, VAR_NAME(I_SFC_ALB_VIS_dif), 'XY',  & ! [IN]
                               OCEAN_SFC_albedo(:,:,I_R_diffuse,I_R_VIS)        ) ! [OUT]
       call FILE_CARTESC_read( restart_fid, VAR_NAME(I_SFC_Z0M),         'XY',  & ! [IN]
                               OCEAN_SFC_Z0M(:,:)                               ) ! [OUT]
       call FILE_CARTESC_read( restart_fid, VAR_NAME(I_SFC_Z0H),         'XY',  & ! [IN]
                               OCEAN_SFC_Z0H(:,:)                               ) ! [OUT]
       call FILE_CARTESC_read( restart_fid, VAR_NAME(I_SFC_Z0E),         'XY',  & ! [IN]
                               OCEAN_SFC_Z0E(:,:)                               ) ! [OUT]

       if( FILE_get_AGGREGATE(restart_fid) ) call FILE_CARTESC_flush( restart_fid ) ! commit all pending read requests


       if ( OCEAN_ICE_TYPE == 'NONE' ) then
          OCEAN_ICE_TEMP(:,:) = OCEAN_PHY_ICE_freezetemp
          OCEAN_ICE_FRAC(:,:) = 0.0_RP
       else
          OCEAN_ICE_TEMP(:,:) = min( OCEAN_ICE_TEMP(:,:), OCEAN_PHY_ICE_freezetemp )
          call OCEAN_PHY_ICE_fraction( OIA, OIS, OIE,       & ! [IN]
                                       OJA, OJS, OJE,       & ! [IN]
                                       OCEAN_ICE_MASS(:,:), & ! [IN]
                                       OCEAN_ICE_FRAC(:,:)  ) ! [OUT]
       endif

       call OCEAN_vars_total
    else
       LOG_ERROR("OCEAN_vars_restart_read",*) 'invalid restart file ID for ocean.'
       call PRC_abort
    endif

    return
  end subroutine OCEAN_vars_restart_read

  !-----------------------------------------------------------------------------
  !> History output set for ocean variables
  subroutine OCEAN_vars_history
    use scale_file_history, only: &
       FILE_HISTORY_in
    use scale_atmos_hydrometeor, only: &
       I_QV
    implicit none
    !---------------------------------------------------------------------------

    call PROF_rapstart('OCN_History', 1)

    if ( OCEAN_VARS_CHECKRANGE ) then
       call VALCHECK( OCEAN_TEMP      (OKS:OKE,OIS:OIE,OJS:OJE),             0.0_RP, 1000.0_RP, &
                      VAR_NAME(I_TEMP),                                      __FILE__, __LINE__ )
!      call VALCHECK( OCEAN_SALT      (OKS:OKE,OIS:OIE,OJS:OJE),             0.0_RP, 1000.0_RP, &
!                     VAR_NAME(I_SALT),                                      __FILE__, __LINE__ )
!      call VALCHECK( OCEAN_UVEL      (OKS:OKE,OIS:OIE,OJS:OJE),             0.0_RP, 1000.0_RP, &
!                     VAR_NAME(I_UVEL),                                      __FILE__, __LINE__ )
!      call VALCHECK( OCEAN_VVEL      (OKS:OKE,OIS:OIE,OJS:OJE),             0.0_RP, 1000.0_RP, &
!                     VAR_NAME(I_VVEL),                                      __FILE__, __LINE__ )

       call VALCHECK( OCEAN_OCN_Z0M   (OIS:OIE,OJS:OJE),                     0.0_RP, 1000.0_RP, &
                      VAR_NAME(I_OCN_Z0M),                                   __FILE__, __LINE__ )
       call VALCHECK( OCEAN_ICE_TEMP  (OIS:OIE,OJS:OJE),                     0.0_RP, 1000.0_RP, &
                      VAR_NAME(I_ICE_TEMP),                                  __FILE__, __LINE__ )
       call VALCHECK( OCEAN_ICE_MASS  (OIS:OIE,OJS:OJE),                     0.0_RP,   5E+5_RP, &
                      VAR_NAME(I_ICE_MASS),                                  __FILE__, __LINE__ )
       call VALCHECK( OCEAN_SFC_TEMP  (OIS:OIE,OJS:OJE),                     0.0_RP, 1000.0_RP, &
                      VAR_NAME(I_SFC_TEMP),                                  __FILE__, __LINE__ )

       call VALCHECK( OCEAN_SFC_albedo(OIS:OIE,OJS:OJE,I_R_direct ,I_R_IR ), 0.0_RP,    2.0_RP, &
                      VAR_NAME(I_SFC_ALB_IR_dir ),                           __FILE__, __LINE__ )
       call VALCHECK( OCEAN_SFC_albedo(OIS:OIE,OJS:OJE,I_R_diffuse,I_R_IR ), 0.0_RP,    2.0_RP, &
                      VAR_NAME(I_SFC_ALB_IR_dif ),                           __FILE__, __LINE__ )
       call VALCHECK( OCEAN_SFC_albedo(OIS:OIE,OJS:OJE,I_R_direct ,I_R_NIR), 0.0_RP,    2.0_RP, &
                      VAR_NAME(I_SFC_ALB_NIR_dir),                           __FILE__, __LINE__ )
       call VALCHECK( OCEAN_SFC_albedo(OIS:OIE,OJS:OJE,I_R_diffuse,I_R_NIR), 0.0_RP,    2.0_RP, &
                      VAR_NAME(I_SFC_ALB_NIR_dif),                           __FILE__, __LINE__ )
       call VALCHECK( OCEAN_SFC_albedo(OIS:OIE,OJS:OJE,I_R_direct ,I_R_VIS), 0.0_RP,    2.0_RP, &
                      VAR_NAME(I_SFC_ALB_VIS_dir),                           __FILE__, __LINE__ )
       call VALCHECK( OCEAN_SFC_albedo(OIS:OIE,OJS:OJE,I_R_diffuse,I_R_VIS), 0.0_RP,    2.0_RP, &
                      VAR_NAME(I_SFC_ALB_VIS_dif),                           __FILE__, __LINE__ )
       call VALCHECK( OCEAN_SFC_Z0M   (OIS:OIE,OJS:OJE),                     0.0_RP, 1000.0_RP, &
                      VAR_NAME(I_SFC_Z0M),                                   __FILE__, __LINE__ )
       call VALCHECK( OCEAN_SFC_Z0H   (OIS:OIE,OJS:OJE),                     0.0_RP, 1000.0_RP, &
                      VAR_NAME(I_SFC_Z0H),                                   __FILE__, __LINE__ )
       call VALCHECK( OCEAN_SFC_Z0E   (OIS:OIE,OJS:OJE),                     0.0_RP, 1000.0_RP, &
                      VAR_NAME(I_SFC_Z0E),                                   __FILE__, __LINE__ )
    endif

    call FILE_HISTORY_in( OCEAN_TEMP      (:,:,:),                                               &
                          VAR_NAME(I_TEMP),            VAR_DESC(I_TEMP),                         &
                          VAR_UNIT(I_TEMP),            standard_name=VAR_STDN(I_TEMP),           &
                          dim_type="OXY"                                                         )
!   call FILE_HISTORY_in( OCEAN_SALT      (:,:,:),                                               &
!                         VAR_NAME(I_SALT),            VAR_DESC(I_SALT),                         &
!                         VAR_UNIT(I_SALT),            standard_name=VAR_STDN(I_SALT),           &
!                         dim_type="OXY"                                                         )
!   call FILE_HISTORY_in( OCEAN_UVEL      (:,:,:),                                               &
!                         VAR_NAME(I_UVEL),            VAR_DESC(I_UVEL),                         &
!                         VAR_UNIT(I_UVEL),            standard_name=VAR_STDN(I_UVEL),           &
!                         dim_type="OXY"                                                         )
!   call FILE_HISTORY_in( OCEAN_VVEL      (:,:,:),                                               &
!                         VAR_NAME(I_VVEL),            VAR_DESC(I_VVEL),                         &
!                         VAR_UNIT(I_VVEL),            standard_name=VAR_STDN(I_VVEL),           &
!                         dim_type="OXY"                                                         )

    call FILE_HISTORY_in( OCEAN_OCN_Z0M   (:,:),                                                 &
                          VAR_NAME(I_OCN_Z0M),         VAR_DESC(I_OCN_Z0M),                      &
                          VAR_UNIT(I_OCN_Z0M),         standard_name=VAR_STDN(I_OCN_Z0M)         )
    call FILE_HISTORY_in( OCEAN_ICE_TEMP  (:,:),                                                 &
                          VAR_NAME(I_ICE_TEMP),        VAR_DESC(I_ICE_TEMP),                     &
                          VAR_UNIT(I_ICE_TEMP),        standard_name=VAR_STDN(I_ICE_TEMP)        )
    call FILE_HISTORY_in( OCEAN_ICE_MASS  (:,:),                                                 &
                          VAR_NAME(I_ICE_MASS),        VAR_DESC(I_ICE_MASS),                     &
                          VAR_UNIT(I_ICE_MASS),        standard_name=VAR_STDN(I_ICE_MASS)        )

    call FILE_HISTORY_in( OCEAN_SFC_TEMP  (:,:),                                                 &
                          VAR_NAME(I_SFC_TEMP),        VAR_DESC(I_SFC_TEMP),                     &
                          VAR_UNIT(I_SFC_TEMP),        standard_name=VAR_STDN(I_SFC_TEMP)        )
    call FILE_HISTORY_in( OCEAN_SFC_albedo(:,:,I_R_direct ,I_R_IR ),                             &
                          VAR_NAME(I_SFC_ALB_IR_dir),  VAR_DESC(I_SFC_ALB_IR_dir),               &
                          VAR_UNIT(I_SFC_ALB_IR_dir),  standard_name=VAR_STDN(I_SFC_ALB_IR_dir)  )
    call FILE_HISTORY_in( OCEAN_SFC_albedo(:,:,I_R_diffuse,I_R_IR ),                             &
                          VAR_NAME(I_SFC_ALB_IR_dif),  VAR_DESC(I_SFC_ALB_IR_dif),               &
                          VAR_UNIT(I_SFC_ALB_IR_dif),  standard_name=VAR_STDN(I_SFC_ALB_IR_dif)  )
    call FILE_HISTORY_in( OCEAN_SFC_albedo(:,:,I_R_direct ,I_R_NIR),                             &
                          VAR_NAME(I_SFC_ALB_NIR_dir), VAR_DESC(I_SFC_ALB_NIR_dir),              &
                          VAR_UNIT(I_SFC_ALB_NIR_dir), standard_name=VAR_STDN(I_SFC_ALB_NIR_dir) )
    call FILE_HISTORY_in( OCEAN_SFC_albedo(:,:,I_R_diffuse,I_R_NIR),                             &
                          VAR_NAME(I_SFC_ALB_NIR_dif), VAR_DESC(I_SFC_ALB_NIR_dif),              &
                          VAR_UNIT(I_SFC_ALB_NIR_dif), standard_name=VAR_STDN(I_SFC_ALB_NIR_dif) )
    call FILE_HISTORY_in( OCEAN_SFC_albedo(:,:,I_R_direct ,I_R_VIS),                             &
                          VAR_NAME(I_SFC_ALB_VIS_dir), VAR_DESC(I_SFC_ALB_VIS_dir),              &
                          VAR_UNIT(I_SFC_ALB_VIS_dir), standard_name=VAR_STDN(I_SFC_ALB_VIS_dir) )
    call FILE_HISTORY_in( OCEAN_SFC_albedo(:,:,I_R_diffuse,I_R_VIS),                             &
                          VAR_NAME(I_SFC_ALB_VIS_dif), VAR_DESC(I_SFC_ALB_VIS_dif),              &
                          VAR_UNIT(I_SFC_ALB_VIS_dif), standard_name=VAR_STDN(I_SFC_ALB_VIS_dif) )
    call FILE_HISTORY_in( OCEAN_SFC_Z0M   (:,:),                                                 &
                          VAR_NAME(I_SFC_Z0M),         VAR_DESC(I_SFC_Z0M),                      &
                          VAR_UNIT(I_SFC_Z0M),         standard_name=VAR_STDN(I_SFC_Z0M)         )
    call FILE_HISTORY_in( OCEAN_SFC_Z0H   (:,:),                                                 &
                          VAR_NAME(I_SFC_Z0H),         VAR_DESC(I_SFC_Z0H),                      &
                          VAR_UNIT(I_SFC_Z0H),         standard_name=VAR_STDN(I_SFC_Z0H)         )
    call FILE_HISTORY_in( OCEAN_SFC_Z0E   (:,:),                                                 &
                          VAR_NAME(I_SFC_Z0E),         VAR_DESC(I_SFC_Z0E),                      &
                          VAR_UNIT(I_SFC_Z0E),         standard_name=VAR_STDN(I_SFC_Z0H)         )

    call FILE_HISTORY_in( OCEAN_SFLX_MW   (:,:),                                                 &
                          VAR_NAME(I_SFLX_MW),         VAR_DESC(I_SFLX_MW),                      &
                          VAR_UNIT(I_SFLX_MW),         standard_name=VAR_STDN(I_SFLX_MW)         )
    call FILE_HISTORY_in( OCEAN_SFLX_MU   (:,:),                                                 &
                          VAR_NAME(I_SFLX_MU),         VAR_DESC(I_SFLX_MU),                      &
                          VAR_UNIT(I_SFLX_MU),         standard_name=VAR_STDN(I_SFLX_MU)         )
    call FILE_HISTORY_in( OCEAN_SFLX_MV   (:,:),                                                 &
                          VAR_NAME(I_SFLX_MV),         VAR_DESC(I_SFLX_MV),                      &
                          VAR_UNIT(I_SFLX_MV),         standard_name=VAR_STDN(I_SFLX_MV)         )
    call FILE_HISTORY_in( OCEAN_SFLX_SH   (:,:),                                                 &
                          VAR_NAME(I_SFLX_SH),         VAR_DESC(I_SFLX_SH),                      &
                          VAR_UNIT(I_SFLX_SH),         standard_name=VAR_STDN(I_SFLX_SH)         )
    call FILE_HISTORY_in( OCEAN_SFLX_LH   (:,:),                                                 &
                          VAR_NAME(I_SFLX_LH),         VAR_DESC(I_SFLX_LH),                      &
                          VAR_UNIT(I_SFLX_LH),         standard_name=VAR_STDN(I_SFLX_LH)         )
    if ( I_QV > 0 ) then
    call FILE_HISTORY_in( OCEAN_SFLX_QTRC (:,:,I_QV),                                            &
                          VAR_NAME(I_SFLX_evap),       VAR_DESC(I_SFLX_evap),                    &
                          VAR_UNIT(I_SFLX_evap),       standard_name=VAR_STDN(I_SFLX_evap)       )
    endif
    call FILE_HISTORY_in( OCEAN_SFLX_G    (:,:),                                                 &
                          VAR_NAME(I_SFLX_G),          VAR_DESC(I_SFLX_G),                       &
                          VAR_UNIT(I_SFLX_G),          standard_name=VAR_STDN(I_SFLX_G)          )
    call FILE_HISTORY_in( OCEAN_SFLX_water(:,:),                                                 &
                          VAR_NAME(I_SFLX_water),      VAR_DESC(I_SFLX_water),                   &
                          VAR_UNIT(I_SFLX_water),      standard_name=VAR_STDN(I_SFLX_water)      )
    call FILE_HISTORY_in( OCEAN_SFLX_ice  (:,:),                                                 &
                          VAR_NAME(I_SFLX_ice),        VAR_DESC(I_SFLX_ice),                     &
                          VAR_UNIT(I_SFLX_ice),        standard_name=VAR_STDN(I_SFLX_ice)        )
    call FILE_HISTORY_in( OCEAN_ICE_FRAC  (:,:),                                                 &
                          VAR_NAME(I_ICE_FRAC),        VAR_DESC(I_ICE_FRAC),                     &
                          VAR_UNIT(I_ICE_FRAC),        standard_name=VAR_STDN(I_ICE_FRAC)        )

    call PROF_rapend  ('OCN_History', 1)

    return
  end subroutine OCEAN_vars_history

  !-----------------------------------------------------------------------------
  !> Budget monitor for ocean
  subroutine OCEAN_vars_total
    use scale_statistics, only: &
       STATISTICS_checktotal, &
       STATISTICS_total
    use scale_atmos_hydrometeor, only: &
       I_QV
    use scale_ocean_grid_cartesC_real, only: &
       OCEAN_GRID_CARTESC_REAL_AREA,    &
       OCEAN_GRID_CARTESC_REAL_TOTAREA, &
       OCEAN_GRID_CARTESC_REAL_VOL,     &
       OCEAN_GRID_CARTESC_REAL_TOTVOL
    implicit none
    !---------------------------------------------------------------------------

    if ( STATISTICS_checktotal ) then

       call STATISTICS_total( OKA, OKS, OKE, OIA, OIS, OIE, OJA, OJS, OJE, & ! [IN]
                              OCEAN_TEMP(:,:,:), VAR_NAME(I_TEMP),         & ! [IN]
                              OCEAN_GRID_CARTESC_REAL_VOL(:,:,:),          & ! [IN]
                              OCEAN_GRID_CARTESC_REAL_TOTVOL               ) ! [IN]
!      call STATISTICS_total( OKA, OKS, OKE, OIA, OIS, OIE, OJA, OJS, OJE, & ! [IN]
!                             OCEAN_SALT(:,:,:), VAR_NAME(I_SALT),         & ! [IN]
!                             OCEAN_GRID_CARTESC_REAL_VOL(:,:,:),          & ! [IN]
!                             OCEAN_GRID_CARTESC_REAL_TOTVOL               ) ! [IN]
!      call STATISTICS_total( OKA, OKS, OKE, OIA, OIS, OIE, OJA, OJS, OJE, & ! [IN]
!                             OCEAN_UVEL(:,:,:), VAR_NAME(I_UVEL),         & ! [IN]
!                             OCEAN_GRID_CARTESC_REAL_VOL(:,:,:),          & ! [IN]
!                             OCEAN_GRID_CARTESC_REAL_TOTVOL               ) ! [IN]
!      call STATISTICS_total( OKA, OKS, OKE, OIA, OIS, OIE, OJA, OJS, OJE, & ! [IN]
!                             OCEAN_VVEL(:,:,:), VAR_NAME(I_VVEL),         & ! [IN]
!                             OCEAN_GRID_CARTESC_REAL_VOL(:,:,:),          & ! [IN]
!                             OCEAN_GRID_CARTESC_REAL_TOTVOL               ) ! [IN]

       call STATISTICS_total( OIA, OIS, OIE, OJA, OJS, OJE,                                           & ! [IN]
                              OCEAN_OCN_Z0M   (:,:),                     VAR_NAME(I_OCN_Z0M),         & ! [IN]
                              OCEAN_GRID_CARTESC_REAL_AREA(:,:),                                      & ! [IN]
                              OCEAN_GRID_CARTESC_REAL_TOTAREA                                         ) ! [IN]
       call STATISTICS_total( OIA, OIS, OIE, OJA, OJS, OJE,                                           & ! [IN]
                              OCEAN_ICE_TEMP  (:,:),                     VAR_NAME(I_ICE_TEMP),        & ! [IN]
                              OCEAN_GRID_CARTESC_REAL_AREA(:,:),                                      & ! [IN]
                              OCEAN_GRID_CARTESC_REAL_TOTAREA                                         ) ! [IN]
       call STATISTICS_total( OIA, OIS, OIE, OJA, OJS, OJE,                                           & ! [IN]
                              OCEAN_ICE_MASS  (:,:),                     VAR_NAME(I_ICE_MASS),        & ! [IN]
                              OCEAN_GRID_CARTESC_REAL_AREA(:,:),                                      & ! [IN]
                              OCEAN_GRID_CARTESC_REAL_TOTAREA                                         ) ! [IN]

       call STATISTICS_total( OIA, OIS, OIE, OJA, OJS, OJE,                                           & ! [IN]
                              OCEAN_SFC_TEMP  (:,:),                     VAR_NAME(I_SFC_TEMP),        & ! [IN]
                              OCEAN_GRID_CARTESC_REAL_AREA(:,:),                                      & ! [IN]
                              OCEAN_GRID_CARTESC_REAL_TOTAREA                                         ) ! [IN]
       call STATISTICS_total( OIA, OIS, OIE, OJA, OJS, OJE,                                           & ! [IN]
                              OCEAN_SFC_albedo(:,:,I_R_direct ,I_R_IR ), VAR_NAME(I_SFC_ALB_IR_dir),  & ! [IN]
                              OCEAN_GRID_CARTESC_REAL_AREA(:,:),                                      & ! [IN]
                              OCEAN_GRID_CARTESC_REAL_TOTAREA                                         ) ! [IN]
       call STATISTICS_total( OIA, OIS, OIE, OJA, OJS, OJE,                                           & ! [IN]
                              OCEAN_SFC_albedo(:,:,I_R_diffuse,I_R_IR ), VAR_NAME(I_SFC_ALB_IR_dif),  & ! [IN]
                              OCEAN_GRID_CARTESC_REAL_AREA(:,:),                                      & ! [IN]
                              OCEAN_GRID_CARTESC_REAL_TOTAREA                                         ) ! [IN]
       call STATISTICS_total( OIA, OIS, OIE, OJA, OJS, OJE,                                           & ! [IN]
                              OCEAN_SFC_albedo(:,:,I_R_direct ,I_R_NIR), VAR_NAME(I_SFC_ALB_NIR_dir), & ! [IN]
                              OCEAN_GRID_CARTESC_REAL_AREA(:,:),                                      & ! [IN]
                              OCEAN_GRID_CARTESC_REAL_TOTAREA                                         ) ! [IN]
       call STATISTICS_total( OIA, OIS, OIE, OJA, OJS, OJE,                                           & ! [IN]
                              OCEAN_SFC_albedo(:,:,I_R_diffuse,I_R_NIR), VAR_NAME(I_SFC_ALB_NIR_dif), & ! [IN]
                              OCEAN_GRID_CARTESC_REAL_AREA(:,:),                                      & ! [IN]
                              OCEAN_GRID_CARTESC_REAL_TOTAREA                                         ) ! [IN]
       call STATISTICS_total( OIA, OIS, OIE, OJA, OJS, OJE,                                           & ! [IN]
                              OCEAN_SFC_albedo(:,:,I_R_direct ,I_R_VIS), VAR_NAME(I_SFC_ALB_VIS_dir), & ! [IN]
                              OCEAN_GRID_CARTESC_REAL_AREA(:,:),                                      & ! [IN]
                              OCEAN_GRID_CARTESC_REAL_TOTAREA                                         ) ! [IN]
       call STATISTICS_total( OIA, OIS, OIE, OJA, OJS, OJE,                                           & ! [IN]
                              OCEAN_SFC_albedo(:,:,I_R_diffuse,I_R_VIS), VAR_NAME(I_SFC_ALB_VIS_dif), & ! [IN]
                              OCEAN_GRID_CARTESC_REAL_AREA(:,:),                                      & ! [IN]
                              OCEAN_GRID_CARTESC_REAL_TOTAREA                                         ) ! [IN]
       call STATISTICS_total( OIA, OIS, OIE, OJA, OJS, OJE,                                           & ! [IN]
                              OCEAN_SFC_Z0M   (:,:),                     VAR_NAME(I_SFC_Z0M),         & ! [IN]
                              OCEAN_GRID_CARTESC_REAL_AREA(:,:),                                      & ! [IN]
                              OCEAN_GRID_CARTESC_REAL_TOTAREA                                         ) ! [IN]
       call STATISTICS_total( OIA, OIS, OIE, OJA, OJS, OJE,                                           & ! [IN]
                              OCEAN_SFC_Z0H   (:,:),                     VAR_NAME(I_SFC_Z0H),         & ! [IN]
                              OCEAN_GRID_CARTESC_REAL_AREA(:,:),                                      & ! [IN]
                              OCEAN_GRID_CARTESC_REAL_TOTAREA                                         ) ! [IN]
       call STATISTICS_total( OIA, OIS, OIE, OJA, OJS, OJE,                                           & ! [IN]
                              OCEAN_SFC_Z0E   (:,:),                     VAR_NAME(I_SFC_Z0E),         & ! [IN]
                              OCEAN_GRID_CARTESC_REAL_AREA(:,:),                                      & ! [IN]
                              OCEAN_GRID_CARTESC_REAL_TOTAREA                                         ) ! [IN]

       call STATISTICS_total( OIA, OIS, OIE, OJA, OJS, OJE,                  & ! [IN]
                              OCEAN_SFLX_MW   (:,:), VAR_NAME(I_SFLX_MW),    & ! [IN]
                              OCEAN_GRID_CARTESC_REAL_AREA(:,:),             & ! [IN]
                              OCEAN_GRID_CARTESC_REAL_TOTAREA                ) ! [IN]
       call STATISTICS_total( OIA, OIS, OIE, OJA, OJS, OJE,                  & ! [IN]
                              OCEAN_SFLX_MU   (:,:), VAR_NAME(I_SFLX_MU),    & ! [IN]
                              OCEAN_GRID_CARTESC_REAL_AREA(:,:),             & ! [IN]
                              OCEAN_GRID_CARTESC_REAL_TOTAREA                ) ! [IN]
       call STATISTICS_total( OIA, OIS, OIE, OJA, OJS, OJE,                  & ! [IN]
                              OCEAN_SFLX_MV   (:,:), VAR_NAME(I_SFLX_MV),    & ! [IN]
                              OCEAN_GRID_CARTESC_REAL_AREA(:,:),             & ! [IN]
                              OCEAN_GRID_CARTESC_REAL_TOTAREA                ) ! [IN]
       call STATISTICS_total( OIA, OIS, OIE, OJA, OJS, OJE,                  & ! [IN]
                              OCEAN_SFLX_SH   (:,:), VAR_NAME(I_SFLX_SH),    & ! [IN]
                              OCEAN_GRID_CARTESC_REAL_AREA(:,:),             & ! [IN]
                              OCEAN_GRID_CARTESC_REAL_TOTAREA                ) ! [IN]
       call STATISTICS_total( OIA, OIS, OIE, OJA, OJS, OJE,                  & ! [IN]
                              OCEAN_SFLX_LH   (:,:), VAR_NAME(I_SFLX_LH),    & ! [IN]
                              OCEAN_GRID_CARTESC_REAL_AREA(:,:),             & ! [IN]
                              OCEAN_GRID_CARTESC_REAL_TOTAREA                ) ! [IN]
       if ( I_QV > 0 ) then
       call STATISTICS_total( OIA, OIS, OIE, OJA, OJS, OJE,                      & ! [IN]
                              OCEAN_SFLX_QTRC (:,:,I_QV), VAR_NAME(I_SFLX_evap), & ! [IN]
                              OCEAN_GRID_CARTESC_REAL_AREA(:,:),                 & ! [IN]
                              OCEAN_GRID_CARTESC_REAL_TOTAREA                    ) ! [IN]
       endif
       call STATISTICS_total( OIA, OIS, OIE, OJA, OJS, OJE,                  & ! [IN]
                              OCEAN_SFLX_G    (:,:), VAR_NAME(I_SFLX_G),     & ! [IN]
                              OCEAN_GRID_CARTESC_REAL_AREA(:,:),             & ! [IN]
                              OCEAN_GRID_CARTESC_REAL_TOTAREA                ) ! [IN]
       call STATISTICS_total( OIA, OIS, OIE, OJA, OJS, OJE,                  & ! [IN]
                              OCEAN_SFLX_water(:,:), VAR_NAME(I_SFLX_water), & ! [IN]
                              OCEAN_GRID_CARTESC_REAL_AREA(:,:),             & ! [IN]
                              OCEAN_GRID_CARTESC_REAL_TOTAREA                ) ! [IN]
       call STATISTICS_total( OIA, OIS, OIE, OJA, OJS, OJE,                  & ! [IN]
                              OCEAN_SFLX_ice  (:,:), VAR_NAME(I_SFLX_ice),   & ! [IN]
                              OCEAN_GRID_CARTESC_REAL_AREA(:,:),             & ! [IN]
                              OCEAN_GRID_CARTESC_REAL_TOTAREA                ) ! [IN]

       call STATISTICS_total( OIA, OIS, OIE, OJA, OJS, OJE,                  & ! [IN]
                              OCEAN_ICE_FRAC  (:,:), VAR_NAME(I_ICE_FRAC),   & ! [IN]
                              OCEAN_GRID_CARTESC_REAL_AREA(:,:),             & ! [IN]
                              OCEAN_GRID_CARTESC_REAL_TOTAREA                ) ! [IN]
    endif

    return
  end subroutine OCEAN_vars_total

  !-----------------------------------------------------------------------------
  !> Create ocean restart file
  subroutine OCEAN_vars_restart_create
    use scale_time, only: &
       TIME_gettimelabel
    use scale_file_cartesC, only: &
       FILE_CARTESC_create
    use mod_ocean_admin, only: &
       OCEAN_do
    implicit none

    character(len=19)     :: timelabel
    character(len=H_LONG) :: basename
    !---------------------------------------------------------------------------

    if ( OCEAN_do .and. OCEAN_RESTART_OUT_BASENAME /= '' ) then
       LOG_NEWLINE
       LOG_INFO("OCEAN_vars_restart_create",*) 'Create restart file (OCEAN) '

       if ( OCEAN_RESTART_OUT_POSTFIX_TIMELABEL ) then
          call TIME_gettimelabel( timelabel )
          basename = trim(OCEAN_RESTART_OUT_BASENAME)//'_'//trim(timelabel)
       else
          basename = trim(OCEAN_RESTART_OUT_BASENAME)
       endif
       LOG_INFO_CONT(*) 'basename: ', trim(basename)

       call FILE_CARTESC_create( basename,                               & ! [IN]
                                 OCEAN_RESTART_OUT_TITLE,                & ! [IN]
                                 OCEAN_RESTART_OUT_DTYPE,                & ! [IN]
                                 restart_fid,                            & ! [OUT]
                                 aggregate = OCEAN_RESTART_OUT_AGGREGATE ) ! [IN]
    endif

    return
  end subroutine OCEAN_vars_restart_create

  !-----------------------------------------------------------------------------
  !> Exit netCDF define mode
  subroutine OCEAN_vars_restart_enddef
    use scale_file_cartesC, only: &
       FILE_CARTESC_enddef
    implicit none
    !---------------------------------------------------------------------------

    if ( restart_fid /= -1 ) then
       call FILE_CARTESC_enddef( restart_fid ) ! [IN]
    endif

    return
  end subroutine OCEAN_vars_restart_enddef

  !-----------------------------------------------------------------------------
  !> Close restart file
  subroutine OCEAN_vars_restart_close
    use scale_file_cartesC, only: &
       FILE_CARTESC_close
    implicit none
    !---------------------------------------------------------------------------

    if ( restart_fid /= -1 ) then
       LOG_NEWLINE
       LOG_INFO("OCEAN_vars_restart_close",*) 'Close restart file (OCEAN) '

       call FILE_CARTESC_close( restart_fid ) ! [IN]

       restart_fid = -1
    endif

    return
  end subroutine OCEAN_vars_restart_close

  !-----------------------------------------------------------------------------
  !> Define ocean variables in restart file
  subroutine OCEAN_vars_restart_def_var
    use scale_file_cartesC, only: &
       FILE_CARTESC_def_var
    implicit none

    integer :: i
    !---------------------------------------------------------------------------

    if ( restart_fid /= -1 ) then
       do i = I_TEMP, I_TEMP
          call FILE_CARTESC_def_var( restart_fid,                           & ! [IN]
                                     VAR_NAME(i), VAR_DESC(i), VAR_UNIT(i), & ! [IN]
                                     'OXY', OCEAN_RESTART_OUT_DTYPE,        & ! [IN]
                                     VAR_ID(i),                             & ! [OUT]
                                     standard_name=VAR_STDN(i)              ) ! [IN]
       enddo

       do i = I_OCN_Z0M, I_SFC_Z0E
          call FILE_CARTESC_def_var( restart_fid,                           & ! [IN]
                                     VAR_NAME(i), VAR_DESC(i), VAR_UNIT(i), & ! [IN]
                                     'XY', OCEAN_RESTART_OUT_DTYPE,         & ! [IN]
                                     VAR_ID(i),                             & ! [OUT]
                                     standard_name=VAR_STDN(i)              ) ! [IN]
       enddo
    endif

    return
  end subroutine OCEAN_vars_restart_def_var

  !-----------------------------------------------------------------------------
  !> Write ocean variables to restart file
  subroutine OCEAN_vars_restart_write
    use scale_file_cartesC, only: &
       FILE_CARTESC_write_var
    implicit none
    !---------------------------------------------------------------------------

    if ( restart_fid /= -1 ) then
       call OCEAN_vars_total

       call FILE_CARTESC_write_var( restart_fid, VAR_ID(I_TEMP),                         & ! [IN]
                                    OCEAN_TEMP(:,:,:),                                   & ! [IN]
                                    VAR_NAME(I_TEMP),            'OXY', fill_halo=.true. ) ! [IN]
!      call FILE_CARTESC_write_var( restart_fid, VAR_ID(I_SALT),                         & ! [IN]
!                                   OCEAN_SALT(:,:,:),                                   & ! [IN]
!                                   VAR_NAME(I_SALT),            'OXY', fill_halo=.true. ) ! [IN]
!      call FILE_CARTESC_write_var( restart_fid, VAR_ID(I_UVEL),                         & ! [IN]
!                                   OCEAN_UVEL(:,:,:),                                   & ! [IN]
!                                   VAR_NAME(I_UVEL),            'OXY', fill_halo=.true. ) ! [IN]
!      call FILE_CARTESC_write_var( restart_fid, VAR_ID(I_VVEL),                         & ! [IN]
!                                   OCEAN_VVEL(:,:,:),                                   & ! [IN]
!                                   VAR_NAME(I_VVEL),            'OXY', fill_halo=.true. ) ! [IN]

       call FILE_CARTESC_write_var( restart_fid, VAR_ID(I_OCN_Z0M),                      & ! [IN]
                                    OCEAN_OCN_Z0M(:,:),                                  & ! [IN]
                                    VAR_NAME(I_OCN_Z0M),         'XY',  fill_halo=.true. ) ! [IN]
       call FILE_CARTESC_write_var( restart_fid, VAR_ID(I_ICE_TEMP),                     & ! [IN]
                                    OCEAN_ICE_TEMP(:,:),                                 & ! [IN]
                                    VAR_NAME(I_ICE_TEMP),        'XY',  fill_halo=.true. ) ! [IN]
       call FILE_CARTESC_write_var( restart_fid, VAR_ID(I_ICE_MASS),                     & ! [IN]
                                    OCEAN_ICE_MASS(:,:),                                 & ! [IN]
                                    VAR_NAME(I_ICE_MASS),        'XY',  fill_halo=.true. ) ! [IN]
       call FILE_CARTESC_write_var( restart_fid, VAR_ID(I_SFC_TEMP),                     & ! [IN]
                                    OCEAN_SFC_TEMP(:,:),                                 & ! [IN]
                                    VAR_NAME(I_SFC_TEMP),        'XY',  fill_halo=.true. ) ! [IN]
       call FILE_CARTESC_write_var( restart_fid, VAR_ID(I_SFC_ALB_IR_dir),               & ! [IN]
                                    OCEAN_SFC_albedo(:,:,I_R_direct ,I_R_IR ),           & ! [IN]
                                    VAR_NAME(I_SFC_ALB_IR_dir),  'XY',  fill_halo=.true. ) ! [IN]
       call FILE_CARTESC_write_var( restart_fid, VAR_ID(I_SFC_ALB_IR_dif),               & ! [IN]
                                    OCEAN_SFC_albedo(:,:,I_R_diffuse,I_R_IR ),           & ! [IN]
                                    VAR_NAME(I_SFC_ALB_IR_dif),  'XY',  fill_halo=.true. ) ! [IN]
       call FILE_CARTESC_write_var( restart_fid, VAR_ID(I_SFC_ALB_NIR_dir),              & ! [IN]
                                    OCEAN_SFC_albedo(:,:,I_R_direct ,I_R_NIR),           & ! [IN]
                                    VAR_NAME(I_SFC_ALB_NIR_dir), 'XY',  fill_halo=.true. ) ! [IN]
       call FILE_CARTESC_write_var( restart_fid, VAR_ID(I_SFC_ALB_NIR_dif),              & ! [IN]
                                    OCEAN_SFC_albedo(:,:,I_R_diffuse,I_R_NIR),           & ! [IN]
                                    VAR_NAME(I_SFC_ALB_NIR_dif), 'XY',  fill_halo=.true. ) ! [IN]
       call FILE_CARTESC_write_var( restart_fid, VAR_ID(I_SFC_ALB_VIS_dir),              & ! [IN]
                                    OCEAN_SFC_albedo(:,:,I_R_direct ,I_R_VIS),           & ! [IN]
                                    VAR_NAME(I_SFC_ALB_VIS_dir), 'XY',  fill_halo=.true. ) ! [IN]
       call FILE_CARTESC_write_var( restart_fid, VAR_ID(I_SFC_ALB_VIS_dif),              & ! [IN]
                                    OCEAN_SFC_albedo(:,:,I_R_diffuse,I_R_VIS),           & ! [IN]
                                    VAR_NAME(I_SFC_ALB_VIS_dif), 'XY',  fill_halo=.true. ) ! [IN]
       call FILE_CARTESC_write_var( restart_fid, VAR_ID(I_SFC_Z0M),                      & ! [IN]
                                    OCEAN_SFC_Z0M(:,:),                                  & ! [IN]
                                    VAR_NAME(I_SFC_Z0M),         'XY',  fill_halo=.true. ) ! [IN]
       call FILE_CARTESC_write_var( restart_fid, VAR_ID(I_SFC_Z0H),                      & ! [IN]
                                    OCEAN_SFC_Z0H(:,:),                                  & ! [IN]
                                    VAR_NAME(I_SFC_Z0H),         'XY',  fill_halo=.true. ) ! [IN]
       call FILE_CARTESC_write_var( restart_fid, VAR_ID(I_SFC_Z0E),                      & ! [IN]
                                    OCEAN_SFC_Z0E(:,:),                                  & ! [IN]
                                    VAR_NAME(I_SFC_Z0E),         'XY',  fill_halo=.true. ) ! [IN]
    endif

    return
  end subroutine OCEAN_vars_restart_write

end module mod_ocean_vars
