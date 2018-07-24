!-------------------------------------------------------------------------------
!> module Atmosphere / Physics Radiation
!!
!! @par Description
!!          Container for mod_atmos_phy_rd
!!
!! @author Team SCALE
!!
!<
!-------------------------------------------------------------------------------
#include "scalelib.h"
module mod_atmos_phy_rd_vars
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_io
  use scale_prof
  use scale_atmos_grid_icoA_index
  use scale_tracer
  use scale_cpl_sfc_index
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: ATMOS_PHY_RD_vars_setup
  public :: ATMOS_PHY_RD_vars_fillhalo
  public :: ATMOS_PHY_RD_vars_restart_read
  public :: ATMOS_PHY_RD_vars_restart_write

  public :: ATMOS_PHY_RD_vars_restart_create
  public :: ATMOS_PHY_RD_vars_restart_open
  public :: ATMOS_PHY_RD_vars_restart_def_var
  public :: ATMOS_PHY_RD_vars_restart_enddef
  public :: ATMOS_PHY_RD_vars_restart_close

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  logical,               public :: ATMOS_PHY_RD_RESTART_OUTPUT                 = .false.                !< output restart file?

  character(len=H_LONG),  public :: ATMOS_PHY_RD_RESTART_IN_BASENAME           = ''                     !< Basename of the input  file
  logical,                public :: ATMOS_PHY_RD_RESTART_IN_AGGREGATE                                   !< Switch to use aggregate file
  logical,                public :: ATMOS_PHY_RD_RESTART_IN_POSTFIX_TIMELABEL  = .false.                !< Add timelabel to the basename of input  file?
  character(len=H_LONG),  public :: ATMOS_PHY_RD_RESTART_OUT_BASENAME          = ''                     !< Basename of the output file
  logical,                public :: ATMOS_PHY_RD_RESTART_OUT_AGGREGATE                                  !< Switch to use aggregate file
  logical,                public :: ATMOS_PHY_RD_RESTART_OUT_POSTFIX_TIMELABEL = .true.                 !< Add timelabel to the basename of output file?
  character(len=H_MID),   public :: ATMOS_PHY_RD_RESTART_OUT_TITLE             = 'ATMOS_PHY_RD restart' !< title    of the output file
  character(len=H_SHORT), public :: ATMOS_PHY_RD_RESTART_OUT_DTYPE             = 'DEFAULT'              !< REAL4 or REAL8

  real(RP), public, allocatable :: ATMOS_PHY_RD_SFLX_LW_up  (:,:,:) ! surface upward   longwave  flux [J/m2/s]
  real(RP), public, allocatable :: ATMOS_PHY_RD_SFLX_LW_dn  (:,:,:) ! surface downward longwave  flux [J/m2/s]
  real(RP), public, allocatable :: ATMOS_PHY_RD_SFLX_SW_up  (:,:,:) ! surface upward   shortwave flux [J/m2/s]
  real(RP), public, allocatable :: ATMOS_PHY_RD_SFLX_SW_dn  (:,:,:) ! surface downward shortwave flux [J/m2/s]

  real(RP), public, allocatable :: ATMOS_PHY_RD_TOAFLX_LW_up(:,:,:) ! TOA upward   longwave  flux [J/m2/s]
  real(RP), public, allocatable :: ATMOS_PHY_RD_TOAFLX_LW_dn(:,:,:) ! TOA downward longwave  flux [J/m2/s]
  real(RP), public, allocatable :: ATMOS_PHY_RD_TOAFLX_SW_up(:,:,:) ! TOA upward   shortwave flux [J/m2/s]
  real(RP), public, allocatable :: ATMOS_PHY_RD_TOAFLX_SW_dn(:,:,:) ! TOA downward shortwave flux [J/m2/s]

  real(RP), public, allocatable :: ATMOS_PHY_RD_SFLX_down   (:,:,:,:,:) ! surface downward flux (LW/SW,direct/diffuse) [J/m2/s]

  real(RP), public, allocatable :: ATMOS_PHY_RD_solins      (:,:,:) ! solar insolation flux   [J/m2/s]
  real(RP), public, allocatable :: ATMOS_PHY_RD_cosSZA      (:,:,:) ! cos(solar zenith angle) (0-1)

  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  integer,                private, parameter :: VMAX              = 14 !< number of the variables
  integer,                private, parameter :: I_SFLX_LW_up      =  1
  integer,                private, parameter :: I_SFLX_LW_dn      =  2
  integer,                private, parameter :: I_SFLX_SW_up      =  3
  integer,                private, parameter :: I_SFLX_SW_dn      =  4
  integer,                private, parameter :: I_TOAFLX_LW_up    =  5
  integer,                private, parameter :: I_TOAFLX_LW_dn    =  6
  integer,                private, parameter :: I_TOAFLX_SW_up    =  7
  integer,                private, parameter :: I_TOAFLX_SW_dn    =  8
  integer,                private, parameter :: I_SFLX_IR_dn_dir  =  9
  integer,                private, parameter :: I_SFLX_IR_dn_dif  = 10
  integer,                private, parameter :: I_SFLX_NIR_dn_dir = 11
  integer,                private, parameter :: I_SFLX_NIR_dn_dif = 12
  integer,                private, parameter :: I_SFLX_VIS_dn_dir = 13
  integer,                private, parameter :: I_SFLX_VIS_dn_dif = 14

  character(len=H_SHORT), private            :: VAR_NAME(VMAX) !< name  of the variables
  character(len=H_MID),   private            :: VAR_DESC(VMAX) !< desc. of the variables
  character(len=H_MID),   private            :: VAR_STDN(VMAX) !< standard name of the variables
  character(len=H_SHORT), private            :: VAR_UNIT(VMAX) !< unit  of the variables
  integer,                private            :: VAR_ID(VMAX)   !< ID    of the variables
  integer,                private            :: restart_fid = -1  ! file ID

  data VAR_NAME / 'SFLX_LW_up',      &
                  'SFLX_LW_dn',      &
                  'SFLX_SW_up',      &
                  'SFLX_SW_dn',      &
                  'TOAFLX_LW_up',    &
                  'TOAFLX_LW_dn',    &
                  'TOAFLX_SW_up',    &
                  'TOAFLX_SW_dn',    &
                  'SFLX_IR_dn_dir',  &
                  'SFLX_IR_dn_dif',  &
                  'SFLX_NIR_dn_dir', &
                  'SFLX_NIR_dn_dif', &
                  'SFLX_VIS_dn_dir', &
                  'SFLX_VIS_dn_dif'  /
  data VAR_DESC / 'surface upward   longwave  flux', &
                  'surface downward longwave  flux', &
                  'surface upward   shortwave flux', &
                  'surface downward shortwave flux', &
                  'TOA upward   longwave  flux',     &
                  'TOA downward longwave  flux',     &
                  'TOA upward   shortwave flux',     &
                  'TOA downward shortwave flux',     &
                  'sfc. down. flux direct  IR',      &
                  'sfc. down. flux diffuse IR',      &
                  'sfc. down. flux direct  NIR',     &
                  'sfc. down. flux diffuse NIR',     &
                  'sfc. down. flux direct  VIS',     &
                  'sfc. down. flux diffuse VIS'      /
  data VAR_STDN / 'surface_upwelling_longwave_flux_in_air',    &
                  'surface_downwelling_longwave_flux_in_air',  &
                  'surface_upwelling_shortwave_flux_in_air',   &
                  'surface_downwelling_shortwave_flux_in_air', &
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
  data VAR_UNIT / 'W/m2', &
                  'W/m2', &
                  'W/m2', &
                  'W/m2', &
                  'W/m2', &
                  'W/m2', &
                  'W/m2', &
                  'W/m2', &
                  'W/m2', &
                  'W/m2', &
                  'W/m2', &
                  'W/m2', &
                  'W/m2', &
                  'W/m2'  /

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine ATMOS_PHY_RD_vars_setup
    use scale_prc, only: &
       PRC_abort
    use scale_const, only: &
       UNDEF => CONST_UNDEF
    implicit none

    namelist / PARAM_ATMOS_PHY_RD_VARS / &
       ATMOS_PHY_RD_RESTART_IN_BASENAME,           &
       ATMOS_PHY_RD_RESTART_IN_AGGREGATE,          &
       ATMOS_PHY_RD_RESTART_IN_POSTFIX_TIMELABEL,  &
       ATMOS_PHY_RD_RESTART_OUTPUT,                &
       ATMOS_PHY_RD_RESTART_OUT_BASENAME,          &
       ATMOS_PHY_RD_RESTART_OUT_AGGREGATE,         &
       ATMOS_PHY_RD_RESTART_OUT_POSTFIX_TIMELABEL, &
       ATMOS_PHY_RD_RESTART_OUT_TITLE,             &
       ATMOS_PHY_RD_RESTART_OUT_DTYPE

    integer :: ierr
    integer :: iv
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("ATMOS_PHY_RD_vars_setup",*) 'Setup'

    allocate( ATMOS_PHY_RD_SFLX_LW_up  (IA,JA,ADM_lall) )
    allocate( ATMOS_PHY_RD_SFLX_LW_dn  (IA,JA,ADM_lall) )
    allocate( ATMOS_PHY_RD_SFLX_SW_up  (IA,JA,ADM_lall) )
    allocate( ATMOS_PHY_RD_SFLX_SW_dn  (IA,JA,ADM_lall) )
    allocate( ATMOS_PHY_RD_TOAFLX_LW_up(IA,JA,ADM_lall) )
    allocate( ATMOS_PHY_RD_TOAFLX_LW_dn(IA,JA,ADM_lall) )
    allocate( ATMOS_PHY_RD_TOAFLX_SW_up(IA,JA,ADM_lall) )
    allocate( ATMOS_PHY_RD_TOAFLX_SW_dn(IA,JA,ADM_lall) )
    ATMOS_PHY_RD_SFLX_LW_up  (:,:,:) = UNDEF
    ATMOS_PHY_RD_SFLX_LW_dn  (:,:,:) = UNDEF
    ATMOS_PHY_RD_SFLX_SW_up  (:,:,:) = UNDEF
    ATMOS_PHY_RD_SFLX_SW_dn  (:,:,:) = UNDEF
    ATMOS_PHY_RD_TOAFLX_LW_up(:,:,:) = UNDEF
    ATMOS_PHY_RD_TOAFLX_LW_dn(:,:,:) = UNDEF
    ATMOS_PHY_RD_TOAFLX_SW_up(:,:,:) = UNDEF
    ATMOS_PHY_RD_TOAFLX_SW_dn(:,:,:) = UNDEF

    allocate( ATMOS_PHY_RD_SFLX_down(IA,JA,N_RAD_DIR,N_RAD_RGN,ADM_lall) )
    ATMOS_PHY_RD_SFLX_down(:,:,:,:,:) = UNDEF

    allocate( ATMOS_PHY_RD_solins(IA,JA,ADM_lall) )
    allocate( ATMOS_PHY_RD_cosSZA(IA,JA,ADM_lall) )
    ATMOS_PHY_RD_solins(:,:,:) = UNDEF
    ATMOS_PHY_RD_cosSZA(:,:,:) = UNDEF

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_ATMOS_PHY_RD_VARS,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("ATMOS_PHY_RD_vars_setup",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("ATMOS_PHY_RD_vars_setup",*) 'Not appropriate names in namelist PARAM_ATMOS_PHY_RD_VARS. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_ATMOS_PHY_RD_VARS)

    LOG_NEWLINE
    LOG_INFO("ATMOS_PHY_RD_vars_setup",*) '[ATMOS_PHY_RD] prognostic/diagnostic variables'
    LOG_INFO_CONT('(1x,A,A24,A,A48,A,A12,A)') &
               '      |', 'VARNAME                 ','|', &
               'DESCRIPTION                                     ', '[', 'UNIT        ', ']'
    do iv = 1, VMAX
       LOG_INFO_CONT('(1x,A,I3,A,A24,A,A48,A,A12,A)') &
                  'NO.',iv,'|',VAR_NAME(iv),'|',VAR_DESC(iv),'[',VAR_UNIT(iv),']'
    enddo

    LOG_NEWLINE
    if ( ATMOS_PHY_RD_RESTART_IN_BASENAME /= '' ) then
       LOG_INFO("ATMOS_PHY_RD_vars_setup",*) 'Restart input?  : YES, file = ', trim(ATMOS_PHY_RD_RESTART_IN_BASENAME)
       LOG_INFO("ATMOS_PHY_RD_vars_setup",*) 'Add timelabel?  : ', ATMOS_PHY_RD_RESTART_IN_POSTFIX_TIMELABEL
    else
       LOG_INFO("ATMOS_PHY_RD_vars_setup",*) 'Restart input?  : NO'
    endif
    if (       ATMOS_PHY_RD_RESTART_OUTPUT             &
         .AND. ATMOS_PHY_RD_RESTART_OUT_BASENAME /= '' ) then
       LOG_INFO("ATMOS_PHY_RD_vars_setup",*) 'Restart output? : YES, file = ', trim(ATMOS_PHY_RD_RESTART_OUT_BASENAME)
       LOG_INFO("ATMOS_PHY_RD_vars_setup",*) 'Add timelabel?  : ', ATMOS_PHY_RD_RESTART_OUT_POSTFIX_TIMELABEL
    else
       LOG_INFO("ATMOS_PHY_RD_vars_setup",*) 'Restart output? : NO'
       ATMOS_PHY_RD_RESTART_OUTPUT = .false.
    endif

    return
  end subroutine ATMOS_PHY_RD_vars_setup

  !-----------------------------------------------------------------------------
  !> HALO Communication
  subroutine ATMOS_PHY_RD_vars_fillhalo

    return
  end subroutine ATMOS_PHY_RD_vars_fillhalo

  !-----------------------------------------------------------------------------
  !> Open restart file for read
  subroutine ATMOS_PHY_RD_vars_restart_open

    return
  end subroutine ATMOS_PHY_RD_vars_restart_open

  !-----------------------------------------------------------------------------
  !> Read restart
  subroutine ATMOS_PHY_RD_vars_restart_read

    return
  end subroutine ATMOS_PHY_RD_vars_restart_read

  !-----------------------------------------------------------------------------
  !> Create restart file
  subroutine ATMOS_PHY_RD_vars_restart_create

    return
  end subroutine ATMOS_PHY_RD_vars_restart_create

  !-----------------------------------------------------------------------------
  !> Exit netCDF define mode
  subroutine ATMOS_PHY_RD_vars_restart_enddef

    return
  end subroutine ATMOS_PHY_RD_vars_restart_enddef

  !-----------------------------------------------------------------------------
  !> Close restart file
  subroutine ATMOS_PHY_RD_vars_restart_close

    return
  end subroutine ATMOS_PHY_RD_vars_restart_close

  !-----------------------------------------------------------------------------
  !> Define variables in restart file
  subroutine ATMOS_PHY_RD_vars_restart_def_var

    return
  end subroutine ATMOS_PHY_RD_vars_restart_def_var

  !-----------------------------------------------------------------------------
  !> Write variables to restart file
  subroutine ATMOS_PHY_RD_vars_restart_write

    return
  end subroutine ATMOS_PHY_RD_vars_restart_write

  !-----------------------------------------------------------------------------
  subroutine ATMOS_PHY_RD_vars_checktotal

    return
  end subroutine ATMOS_PHY_RD_vars_checktotal

end module mod_atmos_phy_rd_vars
