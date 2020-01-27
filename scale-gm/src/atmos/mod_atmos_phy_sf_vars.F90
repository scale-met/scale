!-------------------------------------------------------------------------------
!> module ATMOSPHERIC Surface Variables
!!
!! @par Description
!!          Container for mod_atmos_phy_sf
!!
!! @author Team SCALE
!!
!<
!-------------------------------------------------------------------------------
#include "scalelib.h"
module mod_atmos_phy_sf_vars
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
  public :: ATMOS_PHY_SF_vars_setup
  public :: ATMOS_PHY_SF_vars_fillhalo
  public :: ATMOS_PHY_SF_vars_restart_read
  public :: ATMOS_PHY_SF_vars_restart_write

  public :: ATMOS_PHY_SF_vars_restart_create
  public :: ATMOS_PHY_SF_vars_restart_open
  public :: ATMOS_PHY_SF_vars_restart_def_var
  public :: ATMOS_PHY_SF_vars_restart_enddef
  public :: ATMOS_PHY_SF_vars_restart_close

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  logical,               public :: ATMOS_PHY_SF_RESTART_OUTPUT                 = .false.                !< output restart file?

  character(len=H_LONG),  public :: ATMOS_PHY_SF_RESTART_IN_BASENAME           = ''                     !< Basename of the input  file
  logical,                public :: ATMOS_PHY_SF_RESTART_IN_AGGREGATE                                   !< Switch to use aggregate file
  logical,                public :: ATMOS_PHY_SF_RESTART_IN_POSTFIX_TIMELABEL  = .false.                !< Add timelabel to the basename of input  file?
  character(len=H_LONG),  public :: ATMOS_PHY_SF_RESTART_OUT_BASENAME          = ''                     !< Basename of the output file
  logical,                public :: ATMOS_PHY_SF_RESTART_OUT_AGGREGATE                                  !< Switch to use aggregate file
  logical,                public :: ATMOS_PHY_SF_RESTART_OUT_POSTFIX_TIMELABEL = .true.                 !< Add timelabel to the basename of output file?
  character(len=H_MID),   public :: ATMOS_PHY_SF_RESTART_OUT_TITLE             = 'ATMOS_PHY_SF restart' !< title    of the output file
  character(len=H_SHORT), public :: ATMOS_PHY_SF_RESTART_OUT_DTYPE             = 'DEFAULT'              !< REAL4 or REAL8

  real(RP), public, allocatable :: ATMOS_PHY_SF_SFC_TEMP  (:,:,:)     ! surface skin temperature             [K]
  real(RP), public, allocatable :: ATMOS_PHY_SF_SFC_albedo(:,:,:,:,:) ! surface albedo (direct/diffuse,IR/near-IR/VIS) (0-1)
  real(RP), public, allocatable :: ATMOS_PHY_SF_SFC_Z0M   (:,:,:)     ! surface roughness length, ocean only [m]
  real(RP), public, allocatable :: ATMOS_PHY_SF_SFC_Z0H   (:,:,:)     ! surface roughness length, ocean only [m]
  real(RP), public, allocatable :: ATMOS_PHY_SF_SFC_Z0E   (:,:,:)     ! surface roughness length, ocean only [m]

  ! surface diagnostic variables
  real(RP), public, allocatable :: ATMOS_PHY_SF_SFC_DENS  (:,:,:)     ! surface atmosphere density  [kg/m3]
  real(RP), public, allocatable :: ATMOS_PHY_SF_SFC_PRES  (:,:,:)     ! surface atmosphere pressure [Pa]

  real(RP), public, allocatable :: ATMOS_PHY_SF_SFLX_MW   (:,:,:)     ! z-momentum flux (area center) [m/s*kg/m2/s]
  real(RP), public, allocatable :: ATMOS_PHY_SF_SFLX_MU   (:,:,:)     ! x-momentum flux (area center) [m/s*kg/m2/s]
  real(RP), public, allocatable :: ATMOS_PHY_SF_SFLX_MV   (:,:,:)     ! y-momentum flux (area center) [m/s*kg/m2/s]
  real(RP), public, allocatable :: ATMOS_PHY_SF_SFLX_SH   (:,:,:)     ! sensible heat flux [J/m2/s]
  real(RP), public, allocatable :: ATMOS_PHY_SF_SFLX_LH   (:,:,:)     ! latent heat flux [J/m2/s]
  real(RP), public, allocatable :: ATMOS_PHY_SF_SFLX_GH   (:,:,:)     ! ground heat flux [J/m2/s]
  real(RP), public, allocatable, target :: ATMOS_PHY_SF_SFLX_QTRC (:,:,:,:)   ! tracer mass flux [kg/m2/s]
  real(RP), public, pointer     :: ATMOS_PHY_SF_SFLX_QV   (:,:,:)

  real(RP), public, allocatable :: ATMOS_PHY_SF_U10       (:,:,:)     ! 10m x-wind [m/s]
  real(RP), public, allocatable :: ATMOS_PHY_SF_V10       (:,:,:)     ! 10m y-wind [m/s]
  real(RP), public, allocatable :: ATMOS_PHY_SF_T2        (:,:,:)     ! 2m temperature [K]
  real(RP), public, allocatable :: ATMOS_PHY_SF_Q2        (:,:,:)     ! 2m specific humidity [kg/kg]

  real(RP), public, allocatable :: ATMOS_PHY_SF_Ustar     (:,:,:)     ! friction velocity
  real(RP), public, allocatable :: ATMOS_PHY_SF_Tstar     (:,:,:)     ! temperature scale
  real(RP), public, allocatable :: ATMOS_PHY_SF_Qstar     (:,:,:)     ! moisture scale
  real(RP), public, allocatable :: ATMOS_PHY_SF_Wstar     (:,:,:)     ! convective velocity scale
  real(RP), public, allocatable :: ATMOS_PHY_SF_RLmo      (:,:,:)     ! inverced monin-obukov length

!  real(RP), public, allocatable :: ATMOS_PHY_SF_SFLX_QEMIS(:,:,:,:) ! tracer emission   flux [kg/m2/s]
!  real(RP), public, allocatable :: ATMOS_PHY_SF_SFLX_QDEP (:,:,:,:) ! tracer deposition flux [kg/m2/s]
!  real(RP), public, allocatable :: ATMOS_PHY_SF_SFLX_VDEP (:,:,:,:) ! tracer deposition velocity [m/s]

  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  integer,                private, parameter :: VMAX              = 10 !< number of the variables
  integer,                private, parameter :: I_SFC_TEMP        =  1
  integer,                private, parameter :: I_SFC_ALB_IR_dir  =  2
  integer,                private, parameter :: I_SFC_ALB_IR_dif  =  3
  integer,                private, parameter :: I_SFC_ALB_NIR_dir =  4
  integer,                private, parameter :: I_SFC_ALB_NIR_dif =  5
  integer,                private, parameter :: I_SFC_ALB_VIS_dir =  6
  integer,                private, parameter :: I_SFC_ALB_VIS_dif =  7
  integer,                private, parameter :: I_SFC_Z0M         =  8
  integer,                private, parameter :: I_SFC_Z0H         =  9
  integer,                private, parameter :: I_SFC_Z0E         = 10

  character(len=H_SHORT), private            :: VAR_NAME(VMAX) !< name  of the variables
  character(len=H_MID),   private            :: VAR_DESC(VMAX) !< desc. of the variables
  character(len=H_MID),   private            :: VAR_STDN(VMAX) !< standard name of the variables
  character(len=H_SHORT), private            :: VAR_UNIT(VMAX) !< unit  of the variables
  integer,                private            :: VAR_ID(VMAX)   !< ID    of the variables
  integer,                private            :: restart_fid = -1  ! file ID

  data VAR_NAME / 'SFC_TEMP',        &
                  'SFC_ALB_IR_dir',  &
                  'SFC_ALB_IR_dif',  &
                  'SFC_ALB_NIR_dir', &
                  'SFC_ALB_NIR_dif', &
                  'SFC_ALB_VIS_dir', &
                  'SFC_ALB_VIS_dif', &
                  'SFC_Z0M',         &
                  'SFC_Z0H',         &
                  'SFC_Z0E'          /

  data VAR_DESC / 'surface skin temperature',            &
                  'surface albedo for IR,  direct ',     &
                  'surface albedo for IR,  diffuse',     &
                  'surface albedo for NIR, direct ',     &
                  'surface albedo for NIR, diffuse',     &
                  'surface albedo for VIS, direct ',     &
                  'surface albedo for VIS, diffuse',     &
                  'surface roughness length (momentum)', &
                  'surface roughness length (heat)',     &
                  'surface roughness length (vapor)'     /

  data VAR_STDN / 'surface_temp', &
                  '', &
                  '', &
                  '', &
                  '', &
                  '', &
                  '', &
                  'surface_roughness_length_for_momentum_in_air', &
                  'surface_roughness_length_for_heat_in_air', &
                  '' /

  data VAR_UNIT / 'K', &
                  '1', &
                  '1', &
                  '1', &
                  '1', &
                  '1', &
                  '1', &
                  'm', &
                  'm', &
                  'm'  /

  real(RP), private :: ATMOS_PHY_SF_DEFAULT_SFC_TEMP       = 300.0_RP
  real(RP), private :: ATMOS_PHY_SF_DEFAULT_SFC_albedo_IR  = 0.04_RP
  real(RP), private :: ATMOS_PHY_SF_DEFAULT_SFC_albedo_NIR = 0.04_RP
  real(RP), private :: ATMOS_PHY_SF_DEFAULT_SFC_albedo_VIS = 0.10_RP
  real(RP), private :: ATMOS_PHY_SF_DEFAULT_SFC_Z0M        = 1E-4_RP
  real(RP), private :: ATMOS_PHY_SF_DEFAULT_SFC_Z0H        = 1E-4_RP
  real(RP), private :: ATMOS_PHY_SF_DEFAULT_SFC_Z0E        = 1E-4_RP

  real(RP), allocatable, target :: ZERO(:,:,:)

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine ATMOS_PHY_SF_vars_setup
    use scale_prc, only: &
       PRC_abort
    use scale_const, only: &
       UNDEF => CONST_UNDEF
    use scale_atmos_hydrometeor, only: &
       I_QV
    implicit none

    namelist / PARAM_ATMOS_PHY_SF_VARS / &
       ATMOS_PHY_SF_RESTART_IN_BASENAME,           &
       ATMOS_PHY_SF_RESTART_IN_AGGREGATE,          &
       ATMOS_PHY_SF_RESTART_IN_POSTFIX_TIMELABEL,  &
       ATMOS_PHY_SF_RESTART_OUTPUT,                &
       ATMOS_PHY_SF_RESTART_OUT_BASENAME,          &
       ATMOS_PHY_SF_RESTART_OUT_AGGREGATE,         &
       ATMOS_PHY_SF_RESTART_OUT_POSTFIX_TIMELABEL, &
       ATMOS_PHY_SF_RESTART_OUT_TITLE,             &
       ATMOS_PHY_SF_RESTART_OUT_DTYPE,             &
       ATMOS_PHY_SF_DEFAULT_SFC_TEMP,              &
       ATMOS_PHY_SF_DEFAULT_SFC_albedo_IR,         &
       ATMOS_PHY_SF_DEFAULT_SFC_albedo_NIR,        &
       ATMOS_PHY_SF_DEFAULT_SFC_albedo_VIS,        &
       ATMOS_PHY_SF_DEFAULT_SFC_Z0M,               &
       ATMOS_PHY_SF_DEFAULT_SFC_Z0H,               &
       ATMOS_PHY_SF_DEFAULT_SFC_Z0E

    integer :: ierr
    integer :: iv
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("ATMOS_PHY_SF_vars_setup",*) 'Setup'

    allocate( ATMOS_PHY_SF_SFC_TEMP  (IA,JA,ADM_lall)                     )
    allocate( ATMOS_PHY_SF_SFC_albedo(IA,JA,N_RAD_DIR,N_RAD_RGN,ADM_lall) )
    allocate( ATMOS_PHY_SF_SFC_Z0M   (IA,JA,ADM_lall)                     )
    allocate( ATMOS_PHY_SF_SFC_Z0H   (IA,JA,ADM_lall)                     )
    allocate( ATMOS_PHY_SF_SFC_Z0E   (IA,JA,ADM_lall)                     )
    ATMOS_PHY_SF_SFC_TEMP  (:,:,:)           = ATMOS_PHY_SF_DEFAULT_SFC_TEMP
    ATMOS_PHY_SF_SFC_albedo(:,:,:,I_R_IR, :) = ATMOS_PHY_SF_DEFAULT_SFC_albedo_IR
    ATMOS_PHY_SF_SFC_albedo(:,:,:,I_R_NIR,:) = ATMOS_PHY_SF_DEFAULT_SFC_albedo_NIR
    ATMOS_PHY_SF_SFC_albedo(:,:,:,I_R_VIS,:) = ATMOS_PHY_SF_DEFAULT_SFC_albedo_VIS
    ATMOS_PHY_SF_SFC_Z0M   (:,:,:)           = ATMOS_PHY_SF_DEFAULT_SFC_Z0M
    ATMOS_PHY_SF_SFC_Z0H   (:,:,:)           = ATMOS_PHY_SF_DEFAULT_SFC_Z0H
    ATMOS_PHY_SF_SFC_Z0E   (:,:,:)           = ATMOS_PHY_SF_DEFAULT_SFC_Z0E

    allocate( ATMOS_PHY_SF_SFC_DENS  (IA,JA,ADM_lall) )
    allocate( ATMOS_PHY_SF_SFC_PRES  (IA,JA,ADM_lall) )
    ATMOS_PHY_SF_SFC_DENS  (:,:,:)     = UNDEF
    ATMOS_PHY_SF_SFC_PRES  (:,:,:)     = UNDEF

    allocate( ATMOS_PHY_SF_SFLX_MW   (IA,JA,ADM_lall) )
    allocate( ATMOS_PHY_SF_SFLX_MU   (IA,JA,ADM_lall) )
    allocate( ATMOS_PHY_SF_SFLX_MV   (IA,JA,ADM_lall) )
    allocate( ATMOS_PHY_SF_SFLX_SH   (IA,JA,ADM_lall) )
    allocate( ATMOS_PHY_SF_SFLX_LH   (IA,JA,ADM_lall) )
    allocate( ATMOS_PHY_SF_SFLX_GH   (IA,JA,ADM_lall) )
    allocate( ATMOS_PHY_SF_SFLX_QTRC (IA,JA,max(QA,1),ADM_lall) )
    ATMOS_PHY_SF_SFLX_MW   (:,:,:)     = UNDEF
    ATMOS_PHY_SF_SFLX_MU   (:,:,:)     = UNDEF
    ATMOS_PHY_SF_SFLX_MV   (:,:,:)     = UNDEF
    ATMOS_PHY_SF_SFLX_SH   (:,:,:)     = UNDEF
    ATMOS_PHY_SF_SFLX_LH   (:,:,:)     = UNDEF
    ATMOS_PHY_SF_SFLX_GH   (:,:,:)     = UNDEF
    ATMOS_PHY_SF_SFLX_QTRC (:,:,:,:)   = UNDEF

    allocate( ATMOS_PHY_SF_U10       (IA,JA,ADM_lall) )
    allocate( ATMOS_PHY_SF_V10       (IA,JA,ADM_lall) )
    allocate( ATMOS_PHY_SF_T2        (IA,JA,ADM_lall) )
    allocate( ATMOS_PHY_SF_Q2        (IA,JA,ADM_lall) )
    allocate( ATMOS_PHY_SF_Ustar     (IA,JA,ADM_lall) )
    allocate( ATMOS_PHY_SF_Tstar     (IA,JA,ADM_lall) )
    allocate( ATMOS_PHY_SF_Qstar     (IA,JA,ADM_lall) )
    allocate( ATMOS_PHY_SF_Wstar     (IA,JA,ADM_lall) )
    allocate( ATMOS_PHY_SF_RLmo      (IA,JA,ADM_lall) )
    ATMOS_PHY_SF_U10       (:,:,:)     = UNDEF
    ATMOS_PHY_SF_V10       (:,:,:)     = UNDEF
    ATMOS_PHY_SF_T2        (:,:,:)     = UNDEF
    ATMOS_PHY_SF_Q2        (:,:,:)     = UNDEF
    ATMOS_PHY_SF_Ustar     (:,:,:)     = UNDEF
    ATMOS_PHY_SF_Tstar     (:,:,:)     = UNDEF
    ATMOS_PHY_SF_Qstar     (:,:,:)     = UNDEF
    ATMOS_PHY_SF_Wstar     (:,:,:)     = UNDEF
    ATMOS_PHY_SF_RLmo      (:,:,:)     = UNDEF

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_ATMOS_PHY_SF_VARS,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("ATMOS_PHY_SF_vars_setup",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("ATMOS_PHY_SF_vars_setup",*) 'Not appropriate names in namelist PARAM_ATMOS_PHY_SF_VARS. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_ATMOS_PHY_SF_VARS)

    LOG_NEWLINE
    LOG_INFO("ATMOS_PHY_SF_vars_setup",*) '[ATMOS_PHY_SF] prognostic/diagnostic variables'
    LOG_INFO_CONT('(1x,A,A24,A,A48,A,A12,A)') &
               '      |', 'VARNAME                 ','|', &
               'DESCRIPTION                                     ', '[', 'UNIT        ', ']'
    do iv = 1, VMAX
       LOG_INFO_CONT('(1x,A,I3,A,A24,A,A48,A,A12,A)') &
                  'NO.',iv,'|',VAR_NAME(iv),'|',VAR_DESC(iv),'[',VAR_UNIT(iv),']'
    enddo

    LOG_NEWLINE
    if ( ATMOS_PHY_SF_RESTART_IN_BASENAME /= '' ) then
       LOG_INFO("ATMOS_PHY_SF_vars_setup",*) 'Restart input?  : YES, file = ', trim(ATMOS_PHY_SF_RESTART_IN_BASENAME)
       LOG_INFO("ATMOS_PHY_SF_vars_setup",*) 'Add timelabel?  : ', ATMOS_PHY_SF_RESTART_IN_POSTFIX_TIMELABEL
    else
       LOG_INFO("ATMOS_PHY_SF_vars_setup",*) 'Restart input?  : NO'
    endif
    if (       ATMOS_PHY_SF_RESTART_OUTPUT             &
         .AND. ATMOS_PHY_SF_RESTART_OUT_BASENAME /= '' ) then
       LOG_INFO("ATMOS_PHY_SF_vars_setup",*) 'Restart output? : YES, file = ', trim(ATMOS_PHY_SF_RESTART_OUT_BASENAME)
       LOG_INFO("ATMOS_PHY_SF_vars_setup",*) 'Add timelabel?  : ', ATMOS_PHY_SF_RESTART_OUT_POSTFIX_TIMELABEL
    else
       LOG_INFO("ATMOS_PHY_SF_vars_setup",*) 'Restart output? : NO'
       ATMOS_PHY_SF_RESTART_OUTPUT = .false.
    endif

    if ( I_QV > 0 ) then
       ATMOS_PHY_SF_SFLX_QV => ATMOS_PHY_SF_SFLX_QTRC(:,:,:,I_QV)
    else
       allocate( ZERO(IA,JA,ADM_lall) )
       ATMOS_PHY_SF_SFLX_QV => ZERO
    end if

    return
  end subroutine ATMOS_PHY_SF_vars_setup

  !-----------------------------------------------------------------------------
  !> HALO Communication
  subroutine ATMOS_PHY_SF_vars_fillhalo

    return
  end subroutine ATMOS_PHY_SF_vars_fillhalo

  !-----------------------------------------------------------------------------
  !> Open restart file for read
  subroutine ATMOS_PHY_SF_vars_restart_open

    return
  end subroutine ATMOS_PHY_SF_vars_restart_open

  !-----------------------------------------------------------------------------
  !> Read restart
  subroutine ATMOS_PHY_SF_vars_restart_read

    return
  end subroutine ATMOS_PHY_SF_vars_restart_read

  !-----------------------------------------------------------------------------
  !> Create restart file
  subroutine ATMOS_PHY_SF_vars_restart_create

    return
  end subroutine ATMOS_PHY_SF_vars_restart_create

  !-----------------------------------------------------------------------------
  !> Exit netCDF define mode
  subroutine ATMOS_PHY_SF_vars_restart_enddef

    return
  end subroutine ATMOS_PHY_SF_vars_restart_enddef

  !-----------------------------------------------------------------------------
  !> Close restart file
  subroutine ATMOS_PHY_SF_vars_restart_close

    return
  end subroutine ATMOS_PHY_SF_vars_restart_close

  !-----------------------------------------------------------------------------
  !> Write restart
  subroutine ATMOS_PHY_SF_vars_restart_def_var

    return
  end subroutine ATMOS_PHY_SF_vars_restart_def_var

  !-----------------------------------------------------------------------------
  !> Write variables to restart file
  subroutine ATMOS_PHY_SF_vars_restart_write

    return
  end subroutine ATMOS_PHY_SF_vars_restart_write

  !-----------------------------------------------------------------------------
  subroutine ATMOS_PHY_SF_vars_checktotal

    return
  end subroutine ATMOS_PHY_SF_vars_checktotal

end module mod_atmos_phy_sf_vars
