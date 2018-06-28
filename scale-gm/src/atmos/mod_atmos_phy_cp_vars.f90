!-------------------------------------------------------------------------------
!> module Atmosphere / Physics Cumulus
!!
!! @par Description
!!          Container for mod_atmos_phy_cp
!!
!! @author Team SCALE
!!
!<
!-------------------------------------------------------------------------------
#include "scalelib.h"
module mod_atmos_phy_cp_vars
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_io
  use scale_prof
  use scale_atmos_grid_icoA_index
  use scale_tracer
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: ATMOS_PHY_CP_vars_setup
  public :: ATMOS_PHY_CP_vars_fillhalo
  public :: ATMOS_PHY_CP_vars_restart_read
  public :: ATMOS_PHY_CP_vars_restart_write

  public :: ATMOS_PHY_CP_vars_restart_create
  public :: ATMOS_PHY_CP_vars_restart_open
  public :: ATMOS_PHY_CP_vars_restart_def_var
  public :: ATMOS_PHY_CP_vars_restart_enddef
  public :: ATMOS_PHY_CP_vars_restart_close

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  logical,               public :: ATMOS_PHY_CP_RESTART_OUTPUT                 = .false.                !< output restart file?

  character(len=H_LONG),  public :: ATMOS_PHY_CP_RESTART_IN_BASENAME           = ''                     !< Basename of the input  file
  logical,                public :: ATMOS_PHY_CP_RESTART_IN_POSTFIX_TIMELABEL  = .false.                !< Add timelabel to the basename of input  file?
  logical,                public :: ATMOS_PHY_CP_RESTART_IN_AGGREGATE                                   !< Switch to use aggregate file
  character(len=H_LONG),  public :: ATMOS_PHY_CP_RESTART_OUT_BASENAME          = ''                     !< Basename of the output file
  logical,                public :: ATMOS_PHY_CP_RESTART_OUT_AGGREGATE                                  !< Switch to use aggregate file
  logical,                public :: ATMOS_PHY_CP_RESTART_OUT_POSTFIX_TIMELABEL = .true.                 !< Add timelabel to the basename of output file?
  character(len=H_MID),   public :: ATMOS_PHY_CP_RESTART_OUT_TITLE             = 'ATMOS_PHY_CP restart' !< title    of the output file
  character(len=H_SHORT), public :: ATMOS_PHY_CP_RESTART_OUT_DTYPE             = 'DEFAULT'              !< REAL4 or REAL8

  real(RP), public, allocatable :: ATMOS_PHY_CP_MFLX_cloudbase(:,:,:)   ! cloud base mass flux [kg/m2/s]
  real(RP), public, allocatable :: ATMOS_PHY_CP_SFLX_rain     (:,:,:)   ! convective rain [kg/m2/s]
  real(RP), public, allocatable :: ATMOS_PHY_CP_cloudtop      (:,:,:)   ! cloud top  height [m]
  real(RP), public, allocatable :: ATMOS_PHY_CP_cloudbase     (:,:,:)   ! cloud base height [m]
  real(RP), public, allocatable :: ATMOS_PHY_CP_cldfrac_dp    (:,:,:,:) ! cloud fraction (deep    convection) (0-1)
  real(RP), public, allocatable :: ATMOS_PHY_CP_cldfrac_sh    (:,:,:,:) ! cloud fraction (shallow convection) (0-1)
  real(RP), public, allocatable :: ATMOS_PHY_CP_w0mean        (:,:,:,:) ! running mean vertical wind velocity [m/s]
  ! only for K-F scheme
  real(RP), public, allocatable :: ATMOS_PHY_CP_kf_nca        (:,:,:)   ! advection/cumulus convection timescale/dt for KF[step]

  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  integer,                private, parameter :: VMAX = 8       !< number of the variables
  integer,                private, parameter :: I_MFLX_cloudbase = 1
  integer,                private, parameter :: I_SFLX_convrain  = 2
  integer,                private, parameter :: I_cloudtop       = 3
  integer,                private, parameter :: I_cloudbase      = 4
  integer,                private, parameter :: I_cldfrac_dp     = 5
  integer,                private, parameter :: I_cldfrac_sh     = 6
  integer,                private, parameter :: I_w0mean         = 7
  integer,                private, parameter :: I_kf_nca         = 8


  character(len=H_SHORT), private            :: VAR_NAME(VMAX) !< name  of the variables
  character(len=H_MID),   private            :: VAR_DESC(VMAX) !< desc. of the variables
  character(len=H_SHORT), private            :: VAR_UNIT(VMAX) !< unit  of the variables
  character(len=H_SHORT), private            :: VAR_DIM (VMAX) !< dimension type
  integer,                private            :: VAR_ID  (VMAX) !< ID    of the variables
  integer,                private            :: restart_fid = -1  ! file ID

  data VAR_NAME / 'MFLX_cloudbase',  &
                  'SFLX_convrain',   &
                  'cloudtop',        &
                  'cloudbase',       &
                  'cldfrac_dp',      &
                  'cldfrac_sh',      &
                  'w0mean',          &
                  'kf_nca'           /
  data VAR_DESC / 'cloud base mass flux',                             &
                  'convective rain',                                  &
                  'cloud top height',                                 &
                  'cloud base height',                                &
                  'cloud fraction (deep convection)',                 &
                  'cloud fraction (shallow convection)',              &
                  'running mean vertical wind velocity',              &
                  'advection/cumulus convection timescale/dt for KF'  /
  data VAR_UNIT / 'kg/m2/s', &
                  'kg/m2/s', &
                  'm',       &
                  'm',       &
                  '1',       &
                  '1',       &
                  'm/s',     &
                  'step'     /
  data VAR_DIM  / 'XY',  &
                  'XY',  &
                  'XY',  &
                  'XY',  &
                  'ZXY', &
                  'ZXY', &
                  'ZXY', &
                  'XY'   /

  ! tendency names
  integer,                private              :: VMAX_t       !< number of the tendency variables dens+rhot+qv+N_HYD
  integer,                private              :: I_cp_dens_t = 1
  integer,                private              :: I_cp_rhot_t = 2
  integer,                private              :: I_cp_qv_t   = 3

  character(len=H_SHORT), private, allocatable :: VAR_t_NAME(:) !< name  of the variables
  character(len=H_MID),   private, allocatable :: VAR_t_DESC(:) !< desc. of the variables
  character(len=H_SHORT), private, allocatable :: VAR_t_UNIT(:) !< unit  of the variables
  integer,                private, allocatable :: VAR_t_ID  (:) !< ID    of the variables

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine ATMOS_PHY_CP_vars_setup
    use scale_prc, only: &
       PRC_abort
    use scale_const, only: &
       UNDEF => CONST_UNDEF
    use scale_atmos_hydrometeor, only: &
       N_HYD, &
       HYD_NAME
    implicit none

    namelist / PARAM_ATMOS_PHY_CP_VARS / &
       ATMOS_PHY_CP_RESTART_IN_BASENAME,           &
       ATMOS_PHY_CP_RESTART_IN_AGGREGATE,          &
       ATMOS_PHY_CP_RESTART_IN_POSTFIX_TIMELABEL,  &
       ATMOS_PHY_CP_RESTART_OUTPUT,                &
       ATMOS_PHY_CP_RESTART_OUT_BASENAME,          &
       ATMOS_PHY_CP_RESTART_OUT_AGGREGATE,         &
       ATMOS_PHY_CP_RESTART_OUT_POSTFIX_TIMELABEL, &
       ATMOS_PHY_CP_RESTART_OUT_TITLE,             &
       ATMOS_PHY_CP_RESTART_OUT_DTYPE

    integer :: ierr
    integer :: iv
    integer :: iq
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("ATMOS_PHY_CP_vars_setup",*) 'Setup'

    allocate( ATMOS_PHY_CP_MFLX_cloudbase(IA,JA,ADM_lall)    )
    allocate( ATMOS_PHY_CP_SFLX_rain     (IA,JA,ADM_lall)    )
    allocate( ATMOS_PHY_CP_cloudtop      (IA,JA,ADM_lall)    )
    allocate( ATMOS_PHY_CP_cloudbase     (IA,JA,ADM_lall)    )
    allocate( ATMOS_PHY_CP_cldfrac_dp    (KA,IA,JA,ADM_lall) )
    allocate( ATMOS_PHY_CP_cldfrac_sh    (KA,IA,JA,ADM_lall) )
    allocate( ATMOS_PHY_CP_w0mean        (KA,IA,JA,ADM_lall) )
    allocate( ATMOS_PHY_CP_kf_nca        (IA,JA,ADM_lall)    )
    ATMOS_PHY_CP_MFLX_cloudbase(:,:,:)   =    0.0_RP
    ATMOS_PHY_CP_SFLX_rain     (:,:,:)   =    0.0_RP
    ATMOS_PHY_CP_cloudtop      (:,:,:)   =    0.0_RP
    ATMOS_PHY_CP_cloudbase     (:,:,:)   =    0.0_RP
    ATMOS_PHY_CP_cldfrac_dp    (:,:,:,:) =    0.0_RP
    ATMOS_PHY_CP_cldfrac_sh    (:,:,:,:) =    0.0_RP
    ATMOS_PHY_CP_w0mean        (:,:,:,:) =    0.0_RP
    ATMOS_PHY_CP_kf_nca        (:,:,:)   = -100.0_RP

    ! for tendency restart
    VMAX_t = 3 + N_HYD
    allocate( VAR_t_NAME(VMAX_t) )
    allocate( VAR_t_DESC(VMAX_t) )
    allocate( VAR_t_UNIT(VMAX_t) )
    allocate( VAR_t_ID  (VMAX_t) )

    VAR_t_NAME(I_cp_dens_t) = 'DENS_t_CP'
    VAR_t_DESC(I_cp_dens_t) = 'tendency DENS in CP'
    VAR_t_UNIT(I_cp_dens_t) = 'kg/m3/s'
    VAR_t_NAME(I_cp_rhot_t) = 'RHOT_t_CP'
    VAR_t_DESC(I_cp_rhot_t) = 'tendency RHOT in CP'
    VAR_t_UNIT(I_cp_rhot_t) = 'K*kg/m3/s'

    VAR_t_NAME(I_cp_qv_t) = 'QV_t_CP'
    VAR_t_DESC(I_cp_qv_t) = 'tendency rho*QV in CP'
    VAR_t_UNIT(I_cp_qv_t) = 'kg/m3/s'
    do iq = 1, N_HYD
       VAR_t_NAME(3+iq) = trim(HYD_NAME(iq))//'_t_CP'
       VAR_t_DESC(3+iq) = 'tendency rho*'//trim(HYD_NAME(iq))//' in CP'
       VAR_t_UNIT(3+iq) = 'kg/m3/s'
    enddo

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_ATMOS_PHY_CP_VARS,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("ATMOS_PHY_CP_vars_setup",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("ATMOS_PHY_CP_vars_setup",*) 'Not appropriate names in namelist PARAM_ATMOS_PHY_CP_VARS. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_ATMOS_PHY_CP_VARS)

    LOG_NEWLINE
    LOG_INFO("ATMOS_PHY_CP_vars_setup",*) '[ATMOS_PHY_CP] prognostic/diagnostic variables'
    LOG_INFO_CONT('(1x,A,A24,A,A48,A,A12,A)') &
               '      |', 'VARNAME                 ','|', &
               'DESCRIPTION                                     ', '[', 'UNIT        ', ']'
    do iv = 1, VMAX
       LOG_INFO_CONT('(1x,A,I3,A,A24,A,A48,A,A12,A)') &
                  'NO.',iv,'|',VAR_NAME(iv),'|',VAR_DESC(iv),'[',VAR_UNIT(iv),']'
    enddo

    ! tendency
    do iv = 1, VMAX_t
       LOG_INFO_CONT('(1x,A,I3,A,A24,A,A48,A,A12,A)') &
                  'NO.',iv+VMAX,'|',VAR_t_NAME(iv),'|',VAR_t_DESC(iv),'[',VAR_t_UNIT(iv),']'
    enddo

    LOG_NEWLINE
    if ( ATMOS_PHY_CP_RESTART_IN_BASENAME /= '' ) then
       LOG_INFO("ATMOS_PHY_CP_vars_setup",*) 'Restart input?  : YES, file = ', trim(ATMOS_PHY_CP_RESTART_IN_BASENAME)
       LOG_INFO("ATMOS_PHY_CP_vars_setup",*) 'Add timelabel?  : ', ATMOS_PHY_CP_RESTART_IN_POSTFIX_TIMELABEL
    else
       LOG_INFO("ATMOS_PHY_CP_vars_setup",*) 'Restart input?  : NO'
    endif
    if (       ATMOS_PHY_CP_RESTART_OUTPUT             &
         .AND. ATMOS_PHY_CP_RESTART_OUT_BASENAME /= '' ) then
       LOG_INFO("ATMOS_PHY_CP_vars_setup",*) 'Restart output? : YES, file = ', trim(ATMOS_PHY_CP_RESTART_OUT_BASENAME)
       LOG_INFO("ATMOS_PHY_CP_vars_setup",*) 'Add timelabel?  : ', ATMOS_PHY_CP_RESTART_OUT_POSTFIX_TIMELABEL
    else
       LOG_INFO("ATMOS_PHY_CP_vars_setup",*) 'Restart output? : NO'
       ATMOS_PHY_CP_RESTART_OUTPUT = .false.
    endif

    return
  end subroutine ATMOS_PHY_CP_vars_setup

  !-----------------------------------------------------------------------------
  !> HALO Communication
  subroutine ATMOS_PHY_CP_vars_fillhalo

    return
  end subroutine ATMOS_PHY_CP_vars_fillhalo

  !-----------------------------------------------------------------------------
  !> Open restart file for read
  subroutine ATMOS_PHY_CP_vars_restart_open

    return
  end subroutine ATMOS_PHY_CP_vars_restart_open

  !-----------------------------------------------------------------------------
  !> Read restart
  subroutine ATMOS_PHY_CP_vars_restart_read

    return
  end subroutine ATMOS_PHY_CP_vars_restart_read

  !-----------------------------------------------------------------------------
  !> Create restart file
  subroutine ATMOS_PHY_CP_vars_restart_create

  end subroutine ATMOS_PHY_CP_vars_restart_create

  !-----------------------------------------------------------------------------
  !> End def
  subroutine ATMOS_PHY_CP_vars_restart_enddef
    return
  end subroutine ATMOS_PHY_CP_vars_restart_enddef

  !-----------------------------------------------------------------------------
  !> Close restart file
  subroutine ATMOS_PHY_CP_vars_restart_close

    return
  end subroutine ATMOS_PHY_CP_vars_restart_close

  !-----------------------------------------------------------------------------
  !> Write restart
  subroutine ATMOS_PHY_CP_vars_restart_def_var

    return
  end subroutine ATMOS_PHY_CP_vars_restart_def_var

  !-----------------------------------------------------------------------------
  !> Write restart
  subroutine ATMOS_PHY_CP_vars_restart_write

    return
  end subroutine ATMOS_PHY_CP_vars_restart_write

  subroutine ATMOS_PHY_CP_vars_checktotal

    return
  end subroutine ATMOS_PHY_CP_vars_checktotal

end module mod_atmos_phy_cp_vars
