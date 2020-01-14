!-------------------------------------------------------------------------------
!> module Atmosphere / Physics Cloud Microphysics
!!
!! @par Description
!!          Container for mod_atmos_phy_mp
!!
!! @author Team SCALE
!!
!<
!-------------------------------------------------------------------------------
#include "scalelib.h"
module mod_atmos_phy_mp_vars
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_io
  use scale_prof
  use scale_atmos_grid_cartesC_index
  use scale_tracer
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: ATMOS_PHY_MP_vars_setup
  public :: ATMOS_PHY_MP_vars_fillhalo
  public :: ATMOS_PHY_MP_vars_restart_read
  public :: ATMOS_PHY_MP_vars_restart_write

  public :: ATMOS_PHY_MP_vars_restart_create
  public :: ATMOS_PHY_MP_vars_restart_open
  public :: ATMOS_PHY_MP_vars_restart_def_var
  public :: ATMOS_PHY_MP_vars_restart_enddef
  public :: ATMOS_PHY_MP_vars_restart_close

  public :: ATMOS_PHY_MP_vars_history

  public :: ATMOS_PHY_MP_vars_get_diagnostic
  public :: ATMOS_PHY_MP_vars_reset_diagnostics

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  logical,                public :: ATMOS_PHY_MP_RESTART_OUTPUT                = .false.                !< output restart file?

  character(len=H_LONG),  public :: ATMOS_PHY_MP_RESTART_IN_BASENAME           = ''                     !< Basename of the input  file
  logical,                public :: ATMOS_PHY_MP_RESTART_IN_AGGREGATE                                   !< Switch to use aggregate file
  logical,                public :: ATMOS_PHY_MP_RESTART_IN_POSTFIX_TIMELABEL  = .false.                !< Add timelabel to the basename of input  file?
  character(len=H_LONG),  public :: ATMOS_PHY_MP_RESTART_OUT_BASENAME          = ''                     !< Basename of the output file
  logical,                public :: ATMOS_PHY_MP_RESTART_OUT_AGGREGATE                                  !< Switch to use aggregate file
  logical,                public :: ATMOS_PHY_MP_RESTART_OUT_POSTFIX_TIMELABEL = .true.                 !< Add timelabel to the basename of output file?
  character(len=H_MID),   public :: ATMOS_PHY_MP_RESTART_OUT_TITLE             = 'ATMOS_PHY_MP restart' !< title    of the output file
  character(len=H_SHORT), public :: ATMOS_PHY_MP_RESTART_OUT_DTYPE             = 'DEFAULT'              !< REAL4 or REAL8

  real(RP), public :: ATMOS_PHY_MP_cldfrac_thleshold

  real(RP), public, allocatable :: ATMOS_PHY_MP_DENS_t(:,:,:)    ! tendency DENS [kg/m3/s]
  real(RP), public, allocatable :: ATMOS_PHY_MP_MOMZ_t(:,:,:)    ! tendency MOMZ [kg/m2/s2]
  real(RP), public, allocatable :: ATMOS_PHY_MP_RHOU_t(:,:,:)    ! tendency dens*U [kg/m2/s2]
  real(RP), public, allocatable :: ATMOS_PHY_MP_RHOV_t(:,:,:)    ! tendency dens*V [kg/m2/s2]
  real(RP), public, allocatable :: ATMOS_PHY_MP_RHOT_t(:,:,:)    ! tendency RHOT [K*kg/m3/s]
  real(RP), public, allocatable :: ATMOS_PHY_MP_RHOQ_t(:,:,:,:)  ! tendency rho*QTRC [kg/kg/s]
  real(RP), public, allocatable :: ATMOS_PHY_MP_RHOH  (:,:,:)    ! diabatic heating rate [J/kg/s]

  real(RP), public, allocatable :: ATMOS_PHY_MP_EVAPORATE(:,:,:) ! number concentration of evaporated cloud [/m3]
  real(RP), public, allocatable :: ATMOS_PHY_MP_SFLX_rain(:,:)   ! precipitation flux (liquid) [kg/m2/s]
  real(RP), public, allocatable :: ATMOS_PHY_MP_SFLX_snow(:,:)   ! precipitation flux (solid)  [kg/m2/s]
  real(RP), public, allocatable :: ATMOS_PHY_MP_SFLX_ENGI(:,:)   ! internal energy flux        [J/m2/s]

  integer, public :: QA_MP = 0
  integer, public :: QS_MP = -1
  integer, public :: QE_MP = -2

  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  integer,                private, parameter :: VMAX = 3       !< number of the variables
  integer,                private, parameter :: I_SFLX_rain = 1
  integer,                private, parameter :: I_SFLX_snow = 2
  integer,                private, parameter :: I_SFLX_ENGI = 3

  character(len=H_SHORT), private            :: VAR_NAME(VMAX) !< name  of the variables
  character(len=H_MID),   private            :: VAR_DESC(VMAX) !< desc. of the variables
  character(len=H_SHORT), private            :: VAR_UNIT(VMAX) !< unit  of the variables
  integer,                private            :: VAR_ID(VMAX)   !< ID    of the variables
  integer,                private            :: restart_fid = -1  ! file ID

  data VAR_NAME / 'SFLX_rain', &
                  'SFLX_snow', &
                  'SFLX_ENGI'  /
  data VAR_DESC / 'precipitation flux (liquid)', &
                  'precipitation flux (solid)',  &
                  'internal energy flux'         /
  data VAR_UNIT / 'kg/m2/s', &
                  'kg/m2/s', &
                  'J/m2/s'   /


  ! for diagnostics
  real(RP), private, allocatable :: ATMOS_PHY_MP_CLDFRAC(:,:,:)
  real(RP), private, allocatable :: ATMOS_PHY_MP_Re     (:,:,:,:)
  real(RP), private, allocatable :: ATMOS_PHY_MP_Qe     (:,:,:,:)
  real(RP), private, allocatable :: ATMOS_PHY_MP_Ne     (:,:,:,:)
  logical, private :: DIAG_CLDFRAC
  logical, private :: DIAG_Re
  logical, private :: DIAG_Qe
  logical, private :: DIAG_Ne

  ! for history
  integer, private              :: HIST_CLDFRAC_id
  logical, private              :: HIST_Re
  logical, private              :: HIST_Qe
  logical, private              :: HIST_Ne
  integer, private, allocatable :: HIST_Re_id(:)
  integer, private, allocatable :: HIST_Qe_id(:)
  integer, private, allocatable :: HIST_Ne_id(:)

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine ATMOS_PHY_MP_vars_setup
    use scale_prc, only: &
       PRC_abort
    use scale_const, only: &
       UNDEF => CONST_UNDEF
    use scale_atmos_hydrometeor, only: &
       N_HYD,    &
       QHA,      &
       HYD_NAME, &
       NUM_NAME, &
       HYD_DESC
    use scale_file_history, only: &
       FILE_HISTORY_reg
    implicit none

    namelist / PARAM_ATMOS_PHY_MP_VARS / &
       ATMOS_PHY_MP_RESTART_IN_BASENAME,           &
       ATMOS_PHY_MP_RESTART_IN_AGGREGATE,          &
       ATMOS_PHY_MP_RESTART_IN_POSTFIX_TIMELABEL,  &
       ATMOS_PHY_MP_RESTART_OUTPUT,                &
       ATMOS_PHY_MP_RESTART_OUT_BASENAME,          &
       ATMOS_PHY_MP_RESTART_OUT_AGGREGATE,         &
       ATMOS_PHY_MP_RESTART_OUT_POSTFIX_TIMELABEL, &
       ATMOS_PHY_MP_RESTART_OUT_TITLE,             &
       ATMOS_PHY_MP_RESTART_OUT_DTYPE

    integer :: ierr
    integer :: iv, ih
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("ATMOS_PHY_MP_vars_setup",*) 'Setup'

    allocate( ATMOS_PHY_MP_DENS_t   (KA,IA,JA)    )
    allocate( ATMOS_PHY_MP_MOMZ_t   (KA,IA,JA)    )
    allocate( ATMOS_PHY_MP_RHOU_t   (KA,IA,JA)    )
    allocate( ATMOS_PHY_MP_RHOV_t   (KA,IA,JA)    )
    allocate( ATMOS_PHY_MP_RHOT_t   (KA,IA,JA)    )
    allocate( ATMOS_PHY_MP_RHOQ_t   (KA,IA,JA,QS_MP:QE_MP) )
    allocate( ATMOS_PHY_MP_RHOH     (KA,IA,JA)    )
    allocate( ATMOS_PHY_MP_EVAPORATE(KA,IA,JA)    )
    ! tentative approach
    ATMOS_PHY_MP_DENS_t   (:,:,:)   = 0.0_RP
    ATMOS_PHY_MP_MOMZ_t   (:,:,:)   = 0.0_RP
    ATMOS_PHY_MP_RHOU_t   (:,:,:)   = 0.0_RP
    ATMOS_PHY_MP_RHOV_t   (:,:,:)   = 0.0_RP
    ATMOS_PHY_MP_RHOT_t   (:,:,:)   = 0.0_RP
    ATMOS_PHY_MP_RHOQ_t   (:,:,:,:) = 0.0_RP
    ATMOS_PHY_MP_RHOH     (:,:,:)   = 0.0_RP
    ATMOS_PHY_MP_EVAPORATE(:,:,:)   = 0.0_RP

    allocate( ATMOS_PHY_MP_SFLX_rain(IA,JA) )
    allocate( ATMOS_PHY_MP_SFLX_snow(IA,JA) )
    allocate( ATMOS_PHY_MP_SFLX_ENGI(IA,JA) )
    ATMOS_PHY_MP_SFLX_rain(:,:) = UNDEF
    ATMOS_PHY_MP_SFLX_snow(:,:) = UNDEF
    ATMOS_PHY_MP_SFLX_ENGI(:,:) = UNDEF

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_ATMOS_PHY_MP_VARS,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("ATMOS_PHY_MP_vars_setup",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("ATMOS_PHY_MP_vars_setup",*) 'Not appropriate names in namelist PARAM_ATMOS_PHY_MP_VARS. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_ATMOS_PHY_MP_VARS)

    LOG_NEWLINE
    LOG_INFO("ATMOS_PHY_MP_vars_setup",*) '[ATMOS_PHY_MP] prognostic/diagnostic variables'
    LOG_INFO_CONT('(1x,A,A24,A,A48,A,A12,A)') &
               '      |', 'VARNAME                 ','|', &
               'DESCRIPTION                                     ', '[', 'UNIT        ', ']'
    do iv = 1, VMAX
       LOG_INFO_CONT('(1x,A,I3,A,A24,A,A48,A,A12,A)') &
                  'NO.',iv,'|',VAR_NAME(iv),'|',VAR_DESC(iv),'[',VAR_UNIT(iv),']'
    enddo

    LOG_NEWLINE
    if ( ATMOS_PHY_MP_RESTART_IN_BASENAME /= '' ) then
       LOG_INFO("ATMOS_PHY_MP_vars_setup",*) 'Restart input?  : YES, file = ', trim(ATMOS_PHY_MP_RESTART_IN_BASENAME)
       LOG_INFO("ATMOS_PHY_MP_vars_setup",*) 'Add timelabel?  : ', ATMOS_PHY_MP_RESTART_IN_POSTFIX_TIMELABEL
    else
       LOG_INFO("ATMOS_PHY_MP_vars_setup",*) 'Restart input?  : NO'
    endif
    if (       ATMOS_PHY_MP_RESTART_OUTPUT             &
         .AND. ATMOS_PHY_MP_RESTART_OUT_BASENAME /= '' ) then
       LOG_INFO("ATMOS_PHY_MP_vars_setup",*) 'Restart output? : YES, file = ', trim(ATMOS_PHY_MP_RESTART_OUT_BASENAME)
       LOG_INFO("ATMOS_PHY_MP_vars_setup",*) 'Add timelabel?  : ', ATMOS_PHY_MP_RESTART_OUT_POSTFIX_TIMELABEL
    else
       LOG_INFO("ATMOS_PHY_MP_vars_setup",*) 'Restart output? : NO'
       ATMOS_PHY_MP_RESTART_OUTPUT = .false.
    endif


    ! diagnostices
    allocate( ATMOS_PHY_MP_CLDFRAC(KA,IA,JA) )
    allocate( ATMOS_PHY_MP_Re     (KA,IA,JA,N_HYD) )
    allocate( ATMOS_PHY_MP_Qe     (KA,IA,JA,N_HYD) )
    allocate( ATMOS_PHY_MP_Ne     (KA,IA,JA,N_HYD) )
!OCL XFILL
    ATMOS_PHY_MP_CLDFRAC(:,:,:) = UNDEF
!OCL XFILL
    ATMOS_PHY_MP_Re     (:,:,:,:) = UNDEF
!OCL XFILL
    ATMOS_PHY_MP_Qe     (:,:,:,:) = UNDEF
    ATMOS_PHY_MP_Ne     (:,:,:,:) = UNDEF
    DIAG_CLDFRAC = .false.
    DIAG_Re      = .false.
    DIAG_Qe      = .false.
    DIAG_Ne      = .false.

    ! history
    allocate( HIST_Re_id(N_HYD) )
    allocate( HIST_Qe_id(N_HYD) )
    allocate( HIST_Ne_id(N_HYD) )

    call FILE_HISTORY_reg( 'CLDFRAC', 'cloud fraction', '1', HIST_CLDFRAC_id, fill_halo=.true., dim_type='ZXY' )

    HIST_Re = .false.
    do ih = 1, N_HYD
       call FILE_HISTORY_reg( 'Re_'//trim(HYD_NAME(ih)), 'effective radius of '//trim(HYD_DESC(ih)), 'cm', HIST_Re_id(ih), fill_halo=.true., dim_type='ZXY' )
       if( HIST_Re_id(ih) > 0 ) HIST_Re = .true.
    enddo

    HIST_Qe = .false.
    do ih = 1, N_HYD
       call FILE_HISTORY_reg( trim(HYD_NAME(ih))//'_hyd', 'mass ratio of '//trim(HYD_DESC(ih)), 'kg/kg', HIST_Qe_id(ih), fill_halo=.true., dim_type='ZXY' )
       if( HIST_Qe_id(ih) > 0 ) HIST_Qe = .true.
    enddo

    HIST_Ne = .false.
    do ih = 1, N_HYD
       call FILE_HISTORY_reg( trim(NUM_NAME(ih))//'_hyd', 'number concentration of '//trim(HYD_DESC(ih)), '1/m3', HIST_Ne_id(ih), fill_halo=.true., dim_type='ZXY' )
       if( HIST_Ne_id(ih) > 0 ) HIST_Ne = .true.
    enddo

    return
  end subroutine ATMOS_PHY_MP_vars_setup

  !-----------------------------------------------------------------------------
  !> HALO Communication
  subroutine ATMOS_PHY_MP_vars_fillhalo
    use scale_comm_cartesC, only: &
       COMM_vars8, &
       COMM_wait
    implicit none
    !---------------------------------------------------------------------------

    call COMM_vars8( ATMOS_PHY_MP_SFLX_rain(:,:), 1 )
    call COMM_vars8( ATMOS_PHY_MP_SFLX_snow(:,:), 2 )
    call COMM_vars8( ATMOS_PHY_MP_SFLX_ENGI(:,:), 3 )
    call COMM_wait ( ATMOS_PHY_MP_SFLX_rain(:,:), 1 )
    call COMM_wait ( ATMOS_PHY_MP_SFLX_snow(:,:), 2 )
    call COMM_wait ( ATMOS_PHY_MP_SFLX_ENGI(:,:), 3 )

    return
  end subroutine ATMOS_PHY_MP_vars_fillhalo

  !-----------------------------------------------------------------------------
  !> Open restart file for read
  subroutine ATMOS_PHY_MP_vars_restart_open
    use scale_time, only: &
       TIME_gettimelabel
    use scale_file_cartesC, only: &
       FILE_CARTESC_open
    implicit none

    character(len=19)     :: timelabel
    character(len=H_LONG) :: basename
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("ATMOS_PHY_MP_vars_restart_open",*) 'Open restart file (ATMOS_PHY_MP) '

    if ( ATMOS_PHY_MP_RESTART_IN_BASENAME /= '' ) then

       if ( ATMOS_PHY_MP_RESTART_IN_POSTFIX_TIMELABEL ) then
          call TIME_gettimelabel( timelabel )
          basename = trim(ATMOS_PHY_MP_RESTART_IN_BASENAME)//'_'//trim(timelabel)
       else
          basename = trim(ATMOS_PHY_MP_RESTART_IN_BASENAME)
       endif

       LOG_INFO("ATMOS_PHY_MP_vars_restart_open",*) 'basename: ', trim(basename)

       call FILE_CARTESC_open( basename, restart_fid, aggregate=ATMOS_PHY_MP_RESTART_IN_AGGREGATE )
    else
       LOG_INFO("ATMOS_PHY_MP_vars_restart_open",*) 'restart file for ATMOS_PHY_MP is not specified.'
    endif

    return
  end subroutine ATMOS_PHY_MP_vars_restart_open

  !-----------------------------------------------------------------------------
  !> Read restart
  subroutine ATMOS_PHY_MP_vars_restart_read
    use scale_statistics, only: &
       STATISTICS_checktotal, &
       STATISTICS_total
    use scale_atmos_grid_cartesC_real, only: &
       ATMOS_GRID_CARTESC_REAL_AREA, &
       ATMOS_GRID_CARTESC_REAL_TOTAREA
    use scale_file, only: &
       FILE_get_aggregate
    use scale_file_cartesC, only: &
       FILE_CARTESC_read, &
       FILE_CARTESC_flush
    implicit none

    !---------------------------------------------------------------------------

    if ( restart_fid /= -1 ) then
       LOG_NEWLINE
       LOG_INFO("ATMOS_PHY_MP_vars_restart_read",*) 'Read from restart file (ATMOS_PHY_MP) '

       call FILE_CARTESC_read( restart_fid, VAR_NAME(1), 'XY', & ! [IN]
                               ATMOS_PHY_MP_SFLX_rain(:,:)     ) ! [OUT]
       call FILE_CARTESC_read( restart_fid, VAR_NAME(2), 'XY', & ! [IN]
                               ATMOS_PHY_MP_SFLX_snow(:,:)     ) ! [OUT]
       call FILE_CARTESC_read( restart_fid, VAR_NAME(3), 'XY', & ! [IN]
                               ATMOS_PHY_MP_SFLX_ENGI(:,:)     ) ! [OUT]

       if ( FILE_get_AGGREGATE(restart_fid) ) then
          call FILE_CARTESC_flush( restart_fid ) ! X/Y halos have been read from file
       else
          call ATMOS_PHY_MP_vars_fillhalo
       end if

       if ( STATISTICS_checktotal ) then
          call STATISTICS_total( IA, IS, IE, JA, JS, JE, &
                                 ATMOS_PHY_MP_SFLX_rain(:,:), VAR_NAME(1), & ! (in)
                                 ATMOS_GRID_CARTESC_REAL_AREA(:,:),        & ! (in)
                                 ATMOS_GRID_CARTESC_REAL_TOTAREA           ) ! (in)
          call STATISTICS_total( IA, IS, IE, JA, JS, JE, &
                                 ATMOS_PHY_MP_SFLX_snow(:,:), VAR_NAME(2), & ! (in)
                                 ATMOS_GRID_CARTESC_REAL_AREA(:,:),        & ! (in)
                                 ATMOS_GRID_CARTESC_REAL_TOTAREA           ) ! (in)
          call STATISTICS_total( IA, IS, IE, JA, JS, JE, &
                                 ATMOS_PHY_MP_SFLX_ENGI(:,:), VAR_NAME(3), & ! (in)
                                 ATMOS_GRID_CARTESC_REAL_AREA(:,:),        & ! (in)
                                 ATMOS_GRID_CARTESC_REAL_TOTAREA           ) ! (in)
       endif
    else
       LOG_INFO("ATMOS_PHY_MP_vars_restart_read",*) 'invalid restart file ID for ATMOS_PHY_MP.'
    endif

    return
  end subroutine ATMOS_PHY_MP_vars_restart_read

  !-----------------------------------------------------------------------------
  !> Create restart file
  subroutine ATMOS_PHY_MP_vars_restart_create
    use scale_time, only: &
       TIME_gettimelabel
    use scale_file_cartesC, only: &
       FILE_CARTESC_create
    implicit none

    character(len=19)     :: timelabel
    character(len=H_LONG) :: basename
    !---------------------------------------------------------------------------

    if ( ATMOS_PHY_MP_RESTART_OUT_BASENAME /= '' ) then

       LOG_NEWLINE
       LOG_INFO("ATMOS_PHY_MP_vars_restart_create",*) 'Create restart file (ATMOS_PHY_AE) '

       if ( ATMOS_PHY_MP_RESTART_OUT_POSTFIX_TIMELABEL ) then
          call TIME_gettimelabel( timelabel )
          basename = trim(ATMOS_PHY_MP_RESTART_OUT_BASENAME)//'_'//trim(timelabel)
       else
          basename = trim(ATMOS_PHY_MP_RESTART_OUT_BASENAME)
       endif

       LOG_INFO("ATMOS_PHY_MP_vars_restart_create",*) 'basename: ', trim(basename)

       call FILE_CARTESC_create( &
            basename, ATMOS_PHY_MP_RESTART_OUT_TITLE, ATMOS_PHY_MP_RESTART_OUT_DTYPE, & ! [IN]
            restart_fid,                                                              & ! [OUT]
            aggregate=ATMOS_PHY_MP_RESTART_OUT_AGGREGATE                              ) ! [IN]

    endif

    return
  end subroutine ATMOS_PHY_MP_vars_restart_create

  !-----------------------------------------------------------------------------
  !> Exit netCDF define mode
  subroutine ATMOS_PHY_MP_vars_restart_enddef
    use scale_file_cartesC, only: &
       FILE_CARTESC_enddef
    implicit none

    if ( restart_fid /= -1 ) then
       call FILE_CARTESC_enddef( restart_fid ) ! [IN]
    endif

    return
  end subroutine ATMOS_PHY_MP_vars_restart_enddef

  !-----------------------------------------------------------------------------
  !> Close restart file
  subroutine ATMOS_PHY_MP_vars_restart_close
    use scale_file_cartesC, only: &
       FILE_CARTESC_close
    implicit none
    !---------------------------------------------------------------------------

    if ( restart_fid /= -1 ) then
       LOG_NEWLINE
       LOG_INFO("ATMOS_PHY_MP_vars_restart_close",*) 'Close restart file (ATMOS_PHY_MP) '

       call FILE_CARTESC_close( restart_fid ) ! [IN]

       restart_fid = -1
    endif

    return
  end subroutine ATMOS_PHY_MP_vars_restart_close

  !-----------------------------------------------------------------------------
  !> Define variables in restart file
  subroutine ATMOS_PHY_MP_vars_restart_def_var
    use scale_file_cartesC, only: &
       FILE_CARTESC_def_var
    implicit none
    !---------------------------------------------------------------------------

    if ( restart_fid /= -1 ) then

       call FILE_CARTESC_def_var( restart_fid, VAR_NAME(1), VAR_DESC(1), VAR_UNIT(1), 'XY', ATMOS_PHY_MP_RESTART_OUT_DTYPE, &
                                  VAR_ID(1) )
       call FILE_CARTESC_def_var( restart_fid, VAR_NAME(2), VAR_DESC(2), VAR_UNIT(2), 'XY', ATMOS_PHY_MP_RESTART_OUT_DTYPE, &
                                  VAR_ID(2) )
       call FILE_CARTESC_def_var( restart_fid, VAR_NAME(3), VAR_DESC(3), VAR_UNIT(3), 'XY', ATMOS_PHY_MP_RESTART_OUT_DTYPE, &
                                  VAR_ID(3) )
    endif

    return
  end subroutine ATMOS_PHY_MP_vars_restart_def_var

  !-----------------------------------------------------------------------------
  !> Write restart
  subroutine ATMOS_PHY_MP_vars_restart_write
    use scale_statistics, only: &
       STATISTICS_checktotal, &
       STATISTICS_total
    use scale_atmos_grid_cartesC_real, only: &
       ATMOS_GRID_CARTESC_REAL_AREA, &
       ATMOS_GRID_CARTESC_REAL_TOTAREA
    use scale_file_cartesC, only: &
       FILE_CARTESC_write_var
    implicit none
    !---------------------------------------------------------------------------

    if ( restart_fid /= -1 ) then

       call ATMOS_PHY_MP_vars_fillhalo

       if ( STATISTICS_checktotal ) then
          call STATISTICS_total( IA, IS, IE, JA, JS, JE, &
                                 ATMOS_PHY_MP_SFLX_rain(:,:), VAR_NAME(1), & ! (in)
                                 ATMOS_GRID_CARTESC_REAL_AREA(:,:),        & ! (in)
                                 ATMOS_GRID_CARTESC_REAL_TOTAREA           ) ! (in)
          call STATISTICS_total( IA, IS, IE, JA, JS, JE, &
                                 ATMOS_PHY_MP_SFLX_snow(:,:), VAR_NAME(2), & ! (in)
                                 ATMOS_GRID_CARTESC_REAL_AREA(:,:),        & ! (in)
                                 ATMOS_GRID_CARTESC_REAL_TOTAREA           ) ! (in)
          call STATISTICS_total( IA, IS, IE, JA, JS, JE, &
                                 ATMOS_PHY_MP_SFLX_ENGI(:,:), VAR_NAME(3), & ! (in)
                                 ATMOS_GRID_CARTESC_REAL_AREA(:,:),        & ! (in)
                                 ATMOS_GRID_CARTESC_REAL_TOTAREA           ) ! (in)
       endif

       call FILE_CARTESC_write_var( restart_fid, VAR_ID(1), ATMOS_PHY_MP_SFLX_rain(:,:), &
                              VAR_NAME(1), 'XY' ) ! [IN]
       call FILE_CARTESC_write_var( restart_fid, VAR_ID(2), ATMOS_PHY_MP_SFLX_snow(:,:), &
                              VAR_NAME(2), 'XY' ) ! [IN]
       call FILE_CARTESC_write_var( restart_fid, VAR_ID(3), ATMOS_PHY_MP_SFLX_ENGI(:,:), &
                              VAR_NAME(3), 'XY' ) ! [IN]

    endif

    return
  end subroutine ATMOS_PHY_MP_vars_restart_write

  !-----------------------------------------------------------------------------
  subroutine ATMOS_PHY_MP_vars_history( &
       DENS, TEMP, QTRC )
    use scale_atmos_hydrometeor, only: &
       N_HYD
    use scale_file_history, only: &
       FILE_HISTORY_query, &
       FILE_HISTORY_put
    use mod_atmos_admin, only: &
       ATMOS_PHY_MP_TYPE
    implicit none

    real(RP), intent(in) :: DENS(KA,IA,JA)
    real(RP), intent(in) :: TEMP(KA,IA,JA)
    real(RP), intent(in) :: QTRC(KA,IA,JA,QA)

    real(RP) :: WORK  (KA,IA,JA,N_HYD)
    logical  :: do_put
    integer  :: ih
    !---------------------------------------------------------------------------

    if ( HIST_CLDFRAC_id > 0 ) then
       call FILE_HISTORY_query( HIST_CLDFRAC_id, do_put )

       if ( do_put ) then
          call ATMOS_PHY_MP_vars_get_diagnostic( &
               DENS(:,:,:), TEMP(:,:,:), QTRC(:,:,:,:), & ! [IN]
               CLDFRAC=WORK(:,:,:,1)                    ) ! [OUT]
          call FILE_HISTORY_put( HIST_CLDFRAC_id, WORK(:,:,:,1) )
       end if
    end if

    if ( HIST_Re ) then
       do ih = 1, N_HYD
          if ( HIST_Re_id(ih) > 0 ) then
             call FILE_HISTORY_query( HIST_Re_id(ih), do_put )
             if ( do_put ) then
                call ATMOS_PHY_MP_vars_get_diagnostic( &
                     DENS(:,:,:), TEMP(:,:,:), QTRC(:,:,:,:), & ! [IN]
                     Re=WORK(:,:,:,:)                         ) ! [OUT]
                exit
             end if
          end if
       end do
       if ( do_put ) then
          do ih = 1, N_HYD
             if ( HIST_Re_id(ih) > 0 ) then
                call FILE_HISTORY_query( HIST_Re_id(ih), do_put )
                if ( do_put ) call FILE_HISTORY_put( HIST_Re_id(ih), WORK(:,:,:,ih) )
             end if
          end do
       end if
    end if

    if ( HIST_Qe ) then
       do ih = 1, N_HYD
          if ( HIST_Qe_id(ih) > 0 ) then
             call FILE_HISTORY_query( HIST_Qe_id(ih), do_put )
             if ( do_put ) then
                call ATMOS_PHY_MP_vars_get_diagnostic( &
                     DENS(:,:,:), TEMP(:,:,:), QTRC(:,:,:,:), & ! [IN]
                     Qe=WORK(:,:,:,:)                         ) ! [OUT]
                exit
             end if
          end if
       end do
       if ( do_put ) then
          do ih = 1, N_HYD
             if ( HIST_Qe_id(ih) > 0 ) then
                call FILE_HISTORY_query( HIST_Qe_id(ih), do_put )
                if( do_put ) call FILE_HISTORY_put( HIST_Qe_id(ih), WORK(:,:,:,ih) )
             end if
          end do
       end if
    end if

    if ( HIST_Ne ) then
       do ih = 1, N_HYD
          if ( HIST_Ne_id(ih) > 0 ) then
             call FILE_HISTORY_query( HIST_Ne_id(ih), do_put )
             if ( do_put ) then
                call ATMOS_PHY_MP_vars_get_diagnostic( &
                     DENS(:,:,:), TEMP(:,:,:), QTRC(:,:,:,:), & ! [IN]
                     Ne=WORK(:,:,:,:)                         ) ! [OUT]
                exit
             end if
          end if
       end do
       if ( do_put ) then
          do ih = 1, N_HYD
             if ( HIST_Ne_id(ih) > 0 ) then
                call FILE_HISTORY_query( HIST_Ne_id(ih), do_put )
                if( do_put ) call FILE_HISTORY_put( HIST_Ne_id(ih), WORK(:,:,:,ih) )
             end if
          end do
       end if
    end if

    return
  end subroutine ATMOS_PHY_MP_vars_history

  subroutine ATMOS_PHY_MP_vars_get_diagnostic( &
       DENS, TEMP, QTRC,   &
       CLDFRAC, Re, Qe, Ne )
    use scale_atmos_hydrometeor, only: &
       N_HYD, &
       I_HC,  &
       I_HR,  &
       I_HI,  &
       I_HS,  &
       I_HG,  &
       I_HH
    use scale_atmos_phy_mp_kessler, only: &
       ATMOS_PHY_MP_KESSLER_qtrc2qhyd, &
       ATMOS_PHY_MP_KESSLER_effective_radius, &
       ATMOS_PHY_MP_KESSLER_cloud_fraction
    use scale_atmos_phy_mp_tomita08, only: &
       ATMOS_PHY_MP_TOMITA08_qtrc2qhyd, &
       ATMOS_PHY_MP_TOMITA08_effective_radius, &
       ATMOS_PHY_MP_TOMITA08_cloud_fraction
    use scale_atmos_phy_mp_sn14, only: &
       ATMOS_PHY_MP_SN14_qtrc2qhyd, &
       ATMOS_PHY_MP_SN14_qtrc2nhyd, &
       ATMOS_PHY_MP_SN14_effective_radius, &
       ATMOS_PHY_MP_SN14_cloud_fraction
    use scale_atmos_phy_mp_suzuki10, only: &
       ATMOS_PHY_MP_suzuki10_qtrc2qhyd, &
       ATMOS_PHY_MP_suzuki10_qtrc2nhyd, &
       ATMOS_PHY_MP_suzuki10_effective_radius, &
       ATMOS_PHY_MP_suzuki10_cloud_fraction
    use mod_atmos_admin, only: &
       ATMOS_PHY_MP_TYPE
    implicit none

    real(RP), intent(in) :: DENS(KA,IA,JA)
    real(RP), intent(in) :: TEMP(KA,IA,JA)
    real(RP), intent(in) :: QTRC(KA,IA,JA,QA)
    real(RP), intent(out), optional :: CLDFRAC(KA,IA,JA)       !> cloud fraction [0-1]
    real(RP), intent(out), optional :: Re     (KA,IA,JA,N_HYD) !> effective radius [cm]
    real(RP), intent(out), optional :: Qe     (KA,IA,JA,N_HYD) !> mass ratio [kg/kg]
    real(RP), intent(out), optional :: Ne     (KA,IA,JA,N_HYD) !> number concentratio [1/m3]

    integer :: k, i, j, ih

    if ( present(CLDFRAC) ) then
       if ( .not. DIAG_CLDFRAC ) then
          select case ( ATMOS_PHY_MP_TYPE )
          case ( 'KESSLER' )
             call ATMOS_PHY_MP_kessler_cloud_fraction( &
                  KA, KS, KE, IA, IS, IE, JA, JS, JE, &
                  QTRC(:,:,:,QS_MP+1:QE_MP), ATMOS_PHY_MP_cldfrac_thleshold, & ! [IN]
                  ATMOS_PHY_MP_CLDFRAC(:,:,:)                                ) ! [OUT]
          case ( 'TOMITA08' )
             call ATMOS_PHY_MP_tomita08_cloud_fraction( &
                  KA, KS, KE, IA, ISB, IEB, JA, JSB, JEB, &
                  QTRC(:,:,:,QS_MP+1:QE_MP), ATMOS_PHY_MP_cldfrac_thleshold, & ! [IN]
                  ATMOS_PHY_MP_CLDFRAC(:,:,:)                                ) ! [OUT]
          case ( 'SN14' )
             call ATMOS_PHY_MP_sn14_cloud_fraction( &
                  KA, KS, KE, IA, ISB, IEB, JA, JSB, JEB, &
                  QTRC(:,:,:,QS_MP+1:QE_MP), ATMOS_PHY_MP_cldfrac_thleshold, & ! [IN]
                  ATMOS_PHY_MP_CLDFRAC(:,:,:)                                ) ! [OUT]
          case ( 'SUZUKI10' )
             call ATMOS_PHY_MP_suzuki10_cloud_fraction( &
                  KA, KS, KE, IA, ISB, IEB, JA, JSB, JEB, &
                  QTRC(:,:,:,QS_MP+1:QE_MP), ATMOS_PHY_MP_cldfrac_thleshold, & ! [IN]
                  ATMOS_PHY_MP_CLDFRAC(:,:,:)                                ) ! [OUT]
          case default
!OCL XFILL
             ATMOS_PHY_MP_CLDFRAC(:,:,:) = 0.0_RP
          end select
          DIAG_CLDFRAC = .true.
       end if
!OCL XFILL
       do j = JSB, JEB
       do i = ISB, IEB
       do k = KS, KE
          CLDFRAC(k,i,j) = ATMOS_PHY_MP_CLDFRAC(k,i,j)
       end do
       end do
       end do
    end if

    if ( present(Re) ) then
       if ( .not. DIAG_Re ) then
          select case ( ATMOS_PHY_MP_TYPE )
          case ( 'KESSLER' )
             call ATMOS_PHY_MP_kessler_effective_radius( &
                  KA, KS, KE, IA, IS, IE, JA, JS, JE, &
                  DENS(:,:,:), TEMP(:,:,:), QTRC(:,:,:,QS_MP+1:QE_MP), & ! [IN]
                  ATMOS_PHY_MP_Re(:,:,:,:)                             ) ! [OUT]
          case ( 'TOMITA08' )
             call ATMOS_PHY_MP_tomita08_effective_radius( &
                  KA, KS, KE, IA, ISB, IEB, JA, JSB, JEB, &
                  DENS(:,:,:), TEMP(:,:,:), QTRC(:,:,:,QS_MP+1:QE_MP), & ! [IN]
                  ATMOS_PHY_MP_Re(:,:,:,:)                             ) ! [OUT]
          case ( 'SN14' )
             call ATMOS_PHY_MP_sn14_effective_radius( &
                  KA, KS, KE, IA, ISB, IEB, JA, JSB, JEB, &
                  DENS(:,:,:), TEMP(:,:,:), QTRC(:,:,:,QS_MP+1:QE_MP), & ! [IN]
                  ATMOS_PHY_MP_Re(:,:,:,:)                             ) ! [OUT]
          case ( 'SUZUKI10' )
             call ATMOS_PHY_MP_suzuki10_effective_radius( &
                  KA, KS, KE, IA, ISB, IEB, JA, JSB, JEB, &
                  DENS(:,:,:), TEMP(:,:,:), QTRC(:,:,:,QS_MP+1:QE_MP), & ! [IN]
                  ATMOS_PHY_MP_Re(:,:,:,:)                             ) ! [OUT]
          case default
!OCL XFILL
             ATMOS_PHY_MP_Re(:,:,:,:) = 0.0_RP
          end select
          DIAG_Re = .true.
       end if
!OCL XFILL
       do ih = 1, N_HYD
       do j = JSB, JEB
       do i = ISB, IEB
       do k = KS, KE
          Re(k,i,j,ih) = ATMOS_PHY_MP_Re(k,i,j,ih)
       end do
       end do
       end do
       end do
    end if

    if ( present(Qe) ) then
       if ( .not. DIAG_Qe ) then
          select case ( ATMOS_PHY_MP_TYPE )
          case ( 'KESSLER' )
             call ATMOS_PHY_MP_kessler_qtrc2qhyd( &
                  KA, KS, KE, IA, IS, IE, JA, JS, JE, &
                  QTRC(:,:,:,QS_MP+1:QE_MP), & ! [IN]
                  ATMOS_PHY_MP_Qe(:,:,:,:)   ) ! [OUT]
          case ( 'TOMITA08' )
             call ATMOS_PHY_MP_tomita08_qtrc2qhyd( &
                  KA, KS, KE, IA, ISB, IEB, JA, JSB, JEB, &
                  QTRC(:,:,:,QS_MP+1:QE_MP), & ! [IN]
                  ATMOS_PHY_MP_Qe(:,:,:,:)   ) ! [OUT]
          case ( 'SN14' )
             call ATMOS_PHY_MP_sn14_qtrc2qhyd( &
                  KA, KS, KE, IA, ISB, IEB, JA, JSB, JEB, &
                  QTRC(:,:,:,QS_MP+1:QE_MP), & ! [IN]
                  ATMOS_PHY_MP_Qe(:,:,:,:)   ) ! [OUT]
          case ( 'SUZUKI10' )
             call ATMOS_PHY_MP_suzuki10_qtrc2qhyd( &
                  KA, KS, KE, IA, ISB, IEB, JA, JSB, JEB, &
                  QTRC(:,:,:,QS_MP+1:QE_MP), & ! [IN]
                  ATMOS_PHY_MP_Qe(:,:,:,:)   ) ! [OUT]
          case default
!OCL XFILL
             ATMOS_PHY_MP_Qe(:,:,:,:) = 0.0_RP
          end select
          DIAG_Qe = .true.
       end if
!OCL XFILL
       do ih = 1, N_HYD
       do j = JSB, JEB
       do i = ISB, IEB
       do k = KS, KE
          Qe(k,i,j,ih) = ATMOS_PHY_MP_Qe(k,i,j,ih)
       end do
       end do
       end do
       end do
    end if

    if ( present(Ne) ) then
       if ( .not. DIAG_Ne ) then
          select case ( ATMOS_PHY_MP_TYPE )
          case ( 'KESSLER', 'TOMITA08' )
             ! do nothing
          case ( 'SN14' )
             call ATMOS_PHY_MP_sn14_qtrc2nhyd( &
                  KA, KS, KE, IA, ISB, IEB, JA, JSB, JEB, &
                  QTRC(:,:,:,QS_MP+1:QE_MP), & ! [IN]
                  ATMOS_PHY_MP_Ne(:,:,:,:)   ) ! [OUT]
          case ( 'SUZUKI10' )
             call ATMOS_PHY_MP_suzuki10_qtrc2nhyd( &
                  KA, KS, KE, IA, ISB, IEB, JA, JSB, JEB, &
                  DENS(:,:,:), QTRC(:,:,:,QS_MP+1:QE_MP), & ! [IN]
                  ATMOS_PHY_MP_Ne(:,:,:,:)                ) ! [OUT]
          end select
          DIAG_Ne = .true.
       end if
!OCL XFILL
       do ih = 1, N_HYD
       do j = JSB, JEB
       do i = ISB, IEB
       do k = KS, KE
          Ne(k,i,j,ih) = ATMOS_PHY_MP_Ne(k,i,j,ih)
       end do
       end do
       end do
       end do
    end if

    return
  end subroutine ATMOS_PHY_MP_vars_get_diagnostic

  subroutine ATMOS_PHY_MP_vars_reset_diagnostics
    DIAG_CLDFRAC = .false.
    DIAG_Re      = .false.
    DIAG_Qe      = .false.
    DIAG_Ne      = .false.

    return
  end subroutine ATMOS_PHY_MP_vars_reset_diagnostics

end module mod_atmos_phy_mp_vars
