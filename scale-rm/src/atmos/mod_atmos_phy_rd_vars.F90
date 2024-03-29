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
  public :: ATMOS_PHY_RD_vars_setup
  public :: ATMOS_PHY_RD_vars_finalize
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

  real(RP), public, allocatable :: ATMOS_PHY_RD_RHOH(:,:,:)   ! diavatic heating rate [J/m3/s]

  ! surface
  real(RP), public, allocatable :: ATMOS_PHY_RD_SFLX_LW_up  (:,:) ! upward   longwave  flux [J/m2/s]
  real(RP), public, allocatable :: ATMOS_PHY_RD_SFLX_LW_dn  (:,:) ! downward longwave  flux [J/m2/s]
  real(RP), public, allocatable :: ATMOS_PHY_RD_SFLX_SW_up  (:,:) ! upward   shortwave flux [J/m2/s]
  real(RP), public, allocatable :: ATMOS_PHY_RD_SFLX_SW_dn  (:,:) ! downward shortwave flux [J/m2/s]

  ! TOM (top of the model)
  real(RP), public, allocatable :: ATMOS_PHY_RD_TOMFLX_LW_up(:,:) ! upward   longwave  flux [J/m2/s]
  real(RP), public, allocatable :: ATMOS_PHY_RD_TOMFLX_LW_dn(:,:) ! downward longwave  flux [J/m2/s]
  real(RP), public, allocatable :: ATMOS_PHY_RD_TOMFLX_SW_up(:,:) ! upward   shortwave flux [J/m2/s]
  real(RP), public, allocatable :: ATMOS_PHY_RD_TOMFLX_SW_dn(:,:) ! downward shortwave flux [J/m2/s]

  real(RP), public, allocatable :: ATMOS_PHY_RD_SFLX_down   (:,:,:,:) ! surface downward flux (LW/SW,direct/diffuse) [J/m2/s]

  real(RP), public, allocatable :: ATMOS_PHY_RD_solins      (:,:) ! solar insolation flux   [J/m2/s]
  real(RP), public, allocatable :: ATMOS_PHY_RD_cosSZA      (:,:) ! cos(solar zenith angle) (0-1)

  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !

  ! restart variables
  integer,                private, parameter :: VMAX              =  4 !< number of the variables
  integer,                private, parameter :: I_SFLX_LW_up      =  1
  integer,                private, parameter :: I_SFLX_LW_dn      =  2
  integer,                private, parameter :: I_SFLX_SW_up      =  3
  integer,                private, parameter :: I_SFLX_SW_dn      =  4

  character(len=H_SHORT), private            :: VAR_NAME(VMAX) !< name  of the variables
  character(len=H_MID),   private            :: VAR_DESC(VMAX) !< desc. of the variables
  character(len=H_MID),   private            :: VAR_STDN(VMAX) !< standard name of the variables
  character(len=H_SHORT), private            :: VAR_UNIT(VMAX) !< unit  of the variables
  integer,                private            :: VAR_ID(VMAX)   !< ID    of the variables
  integer,                private            :: restart_fid = -1  ! file ID

  data VAR_NAME / 'SFLX_LW_up', &
                  'SFLX_LW_dn', &
                  'SFLX_SW_up', &
                  'SFLX_SW_dn'  /
  data VAR_DESC / 'surface upward   longwave  flux', &
                  'surface downward longwave  flux', &
                  'surface upward   shortwave flux', &
                  'surface downward shortwave flux'  /
  data VAR_STDN / 'surface_upwelling_longwave_flux_in_air',   &
                  'surface_downwelling_longwave_flux_in_air', &
                  'surface_upwelling_shortwave_flux_in_air',  &
                  'surface_downwelling_shortwave_flux_in_air' /
  data VAR_UNIT / 'W/m2', &
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

    allocate( ATMOS_PHY_RD_RHOH(KA,IA,JA) )
    ATMOS_PHY_RD_RHOH(:,:,:) = UNDEF
    !$acc enter data create(ATMOS_PHY_RD_RHOH)

    allocate( ATMOS_PHY_RD_SFLX_LW_up  (IA,JA) )
    allocate( ATMOS_PHY_RD_SFLX_LW_dn  (IA,JA) )
    allocate( ATMOS_PHY_RD_SFLX_SW_up  (IA,JA) )
    allocate( ATMOS_PHY_RD_SFLX_SW_dn  (IA,JA) )
    ATMOS_PHY_RD_SFLX_LW_up  (:,:) = UNDEF
    ATMOS_PHY_RD_SFLX_LW_dn  (:,:) = UNDEF
    ATMOS_PHY_RD_SFLX_SW_up  (:,:) = UNDEF
    ATMOS_PHY_RD_SFLX_SW_dn  (:,:) = UNDEF
    !$acc enter data create(ATMOS_PHY_RD_SFLX_LW_up,ATMOS_PHY_RD_SFLX_LW_dn,ATMOS_PHY_RD_SFLX_SW_up,ATMOS_PHY_RD_SFLX_SW_dn)

    allocate( ATMOS_PHY_RD_TOMFLX_LW_up(IA,JA) )
    allocate( ATMOS_PHY_RD_TOMFLX_LW_dn(IA,JA) )
    allocate( ATMOS_PHY_RD_TOMFLX_SW_up(IA,JA) )
    allocate( ATMOS_PHY_RD_TOMFLX_SW_dn(IA,JA) )
    ATMOS_PHY_RD_TOMFLX_LW_up(:,:) = UNDEF
    ATMOS_PHY_RD_TOMFLX_LW_dn(:,:) = UNDEF
    ATMOS_PHY_RD_TOMFLX_SW_up(:,:) = UNDEF
    ATMOS_PHY_RD_TOMFLX_SW_dn(:,:) = UNDEF
    !$acc enter data create(ATMOS_PHY_RD_TOMFLX_LW_up,ATMOS_PHY_RD_TOMFLX_LW_dn,ATMOS_PHY_RD_TOMFLX_SW_up,ATMOS_PHY_RD_TOMFLX_SW_dn)

    allocate( ATMOS_PHY_RD_SFLX_down(IA,JA,N_RAD_DIR,N_RAD_RGN) )
    ATMOS_PHY_RD_SFLX_down(:,:,:,:) = UNDEF
    !$acc enter data create(ATMOS_PHY_RD_SFLX_down)

    allocate( ATMOS_PHY_RD_solins(IA,JA) )
    allocate( ATMOS_PHY_RD_cosSZA(IA,JA) )
    ATMOS_PHY_RD_solins(:,:) = UNDEF
    ATMOS_PHY_RD_cosSZA(:,:) = UNDEF
    !$acc enter data create(ATMOS_PHY_RD_solins,ATMOS_PHY_RD_cosSZA)

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
  !> Finalize
  subroutine ATMOS_PHY_RD_vars_finalize
    implicit none
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("ATMOS_PHY_RD_vars_finalize",*) 'Finalize'

    !$acc exit data delete(ATMOS_PHY_RD_RHOH)
    deallocate( ATMOS_PHY_RD_RHOH )

    !$acc exit data delete(ATMOS_PHY_RD_SFLX_LW_up,ATMOS_PHY_RD_SFLX_LW_dn,ATMOS_PHY_RD_SFLX_SW_up,ATMOS_PHY_RD_SFLX_SW_dn)
    deallocate( ATMOS_PHY_RD_SFLX_LW_up )
    deallocate( ATMOS_PHY_RD_SFLX_LW_dn )
    deallocate( ATMOS_PHY_RD_SFLX_SW_up )
    deallocate( ATMOS_PHY_RD_SFLX_SW_dn )

    !$acc exit data delete(ATMOS_PHY_RD_TOMFLX_LW_up,ATMOS_PHY_RD_TOMFLX_LW_dn,ATMOS_PHY_RD_TOMFLX_SW_up,ATMOS_PHY_RD_TOMFLX_SW_dn)
    deallocate( ATMOS_PHY_RD_TOMFLX_LW_up )
    deallocate( ATMOS_PHY_RD_TOMFLX_LW_dn )
    deallocate( ATMOS_PHY_RD_TOMFLX_SW_up )
    deallocate( ATMOS_PHY_RD_TOMFLX_SW_dn )

    !$acc exit data delete(ATMOS_PHY_RD_SFLX_down)
    deallocate( ATMOS_PHY_RD_SFLX_down )

    !$acc exit data delete(ATMOS_PHY_RD_solins,ATMOS_PHY_RD_cosSZA)
    deallocate( ATMOS_PHY_RD_solins )
    deallocate( ATMOS_PHY_RD_cosSZA )

    return
  end subroutine ATMOS_PHY_RD_vars_finalize

  !-----------------------------------------------------------------------------
  !> HALO Communication
  subroutine ATMOS_PHY_RD_vars_fillhalo
    use scale_comm_cartesC, only: &
       COMM_vars8, &
       COMM_wait
    implicit none

    integer :: n ,idir, irgn
    !---------------------------------------------------------------------------

    call COMM_vars8( ATMOS_PHY_RD_SFLX_LW_up  (:,:),  1 )
    call COMM_vars8( ATMOS_PHY_RD_SFLX_LW_dn  (:,:),  2 )
    call COMM_vars8( ATMOS_PHY_RD_SFLX_SW_up  (:,:),  3 )
    call COMM_vars8( ATMOS_PHY_RD_SFLX_SW_dn  (:,:),  4 )
    call COMM_vars8( ATMOS_PHY_RD_TOMFLX_LW_up(:,:),  5 )
    call COMM_vars8( ATMOS_PHY_RD_TOMFLX_LW_dn(:,:),  6 )
    call COMM_vars8( ATMOS_PHY_RD_TOMFLX_SW_up(:,:),  7 )
    call COMM_vars8( ATMOS_PHY_RD_TOMFLX_SW_dn(:,:),  8 )

    n = 8
    do irgn = I_R_IR, I_R_VIS
    do idir = I_R_direct, I_R_diffuse
       n = n + 1
       call COMM_vars8( ATMOS_PHY_RD_SFLX_down(:,:,idir,irgn), n )
    enddo
    enddo

    call COMM_wait ( ATMOS_PHY_RD_SFLX_LW_up  (:,:),  1 )
    call COMM_wait ( ATMOS_PHY_RD_SFLX_LW_dn  (:,:),  2 )
    call COMM_wait ( ATMOS_PHY_RD_SFLX_SW_up  (:,:),  3 )
    call COMM_wait ( ATMOS_PHY_RD_SFLX_SW_dn  (:,:),  4 )
    call COMM_wait ( ATMOS_PHY_RD_TOMFLX_LW_up(:,:),  5 )
    call COMM_wait ( ATMOS_PHY_RD_TOMFLX_LW_dn(:,:),  6 )
    call COMM_wait ( ATMOS_PHY_RD_TOMFLX_SW_up(:,:),  7 )
    call COMM_wait ( ATMOS_PHY_RD_TOMFLX_SW_dn(:,:),  8 )

    n = 8
    do irgn = I_R_IR, I_R_VIS
    do idir = I_R_direct, I_R_diffuse
       n = n + 1
       call COMM_wait ( ATMOS_PHY_RD_SFLX_down(:,:,idir,irgn), n )
    enddo
    enddo

    return
  end subroutine ATMOS_PHY_RD_vars_fillhalo

  !-----------------------------------------------------------------------------
  !> Open restart file for read
  subroutine ATMOS_PHY_RD_vars_restart_open
    use scale_time, only: &
       TIME_gettimelabel
    use scale_file_cartesC, only: &
       FILE_CARTESC_open
    implicit none

    character(len=19)     :: timelabel
    character(len=H_LONG) :: basename
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("ATMOS_PHY_RD_vars_restart_open",*) 'Open restart file (ATMOS_PHY_RD) '

    if ( ATMOS_PHY_RD_RESTART_IN_BASENAME /= '' ) then

       if ( ATMOS_PHY_RD_RESTART_IN_POSTFIX_TIMELABEL ) then
          call TIME_gettimelabel( timelabel )
          basename = trim(ATMOS_PHY_RD_RESTART_IN_BASENAME)//'_'//trim(timelabel)
       else
          basename = trim(ATMOS_PHY_RD_RESTART_IN_BASENAME)
       endif

       LOG_INFO("ATMOS_PHY_RD_vars_restart_open",*) 'basename: ', trim(basename)

       call FILE_CARTESC_open( basename, restart_fid, aggregate=ATMOS_PHY_RD_RESTART_IN_AGGREGATE )
    else
       LOG_INFO("ATMOS_PHY_RD_vars_restart_open",*) 'restart file for ATMOS_PHY_RD is not specified.'
    endif

    return
  end subroutine ATMOS_PHY_RD_vars_restart_open

  !-----------------------------------------------------------------------------
  !> Read restart
  subroutine ATMOS_PHY_RD_vars_restart_read
    use scale_file, only: &
       FILE_get_aggregate
    use scale_file_cartesC, only: &
       FILE_CARTESC_read, &
       FILE_CARTESC_flush
    implicit none
    !---------------------------------------------------------------------------

    if ( restart_fid /= -1 ) then
       LOG_NEWLINE
       LOG_INFO("ATMOS_PHY_RD_vars_restart_read",*) 'Read from restart file (ATMOS_PHY_RD) '

       call FILE_CARTESC_read( restart_fid, VAR_NAME(I_SFLX_LW_up),      'XY', & ! [IN]
                               ATMOS_PHY_RD_SFLX_LW_up(:,:)                    ) ! [OUT]
       call FILE_CARTESC_read( restart_fid, VAR_NAME(I_SFLX_LW_dn),      'XY', & ! [IN]
                               ATMOS_PHY_RD_SFLX_LW_dn(:,:)                    ) ! [OUT]
       call FILE_CARTESC_read( restart_fid, VAR_NAME(I_SFLX_SW_up),      'XY', & ! [IN]
                               ATMOS_PHY_RD_SFLX_SW_up(:,:)                    ) ! [OUT]
       call FILE_CARTESC_read( restart_fid, VAR_NAME(I_SFLX_SW_dn),      'XY', & ! [IN]
                               ATMOS_PHY_RD_SFLX_SW_dn(:,:)                    ) ! [OUT]

       if ( FILE_get_AGGREGATE(restart_fid) ) then
          call FILE_CARTESC_flush( restart_fid ) ! X/Y halos have been read from file
          !$acc update device(ATMOS_PHY_RD_SFLX_LW_up,ATMOS_PHY_RD_SFLX_LW_dn,ATMOS_PHY_RD_SFLX_SW_up,ATMOS_PHY_RD_SFLX_SW_dn)
       else
          call ATMOS_PHY_RD_vars_fillhalo
       end if

       call ATMOS_PHY_RD_vars_check

    else
       LOG_INFO("ATMOS_PHY_RD_vars_restart_read",*) 'invalid restart file ID for ATMOS_PHY_RD.'
    endif

    return
  end subroutine ATMOS_PHY_RD_vars_restart_read

  !-----------------------------------------------------------------------------
  !> Create restart file
  subroutine ATMOS_PHY_RD_vars_restart_create
    use scale_time, only: &
       TIME_gettimelabel
    use scale_file_cartesC, only: &
       FILE_CARTESC_create
    implicit none

    character(len=19)     :: timelabel
    character(len=H_LONG) :: basename
    !---------------------------------------------------------------------------

    if ( ATMOS_PHY_RD_RESTART_OUT_BASENAME /= '' ) then

       LOG_NEWLINE
       LOG_INFO("ATMOS_PHY_RD_vars_restart_create",*) 'Create restart file (ATMOS_PHY_AE) '

       if ( ATMOS_PHY_RD_RESTART_OUT_POSTFIX_TIMELABEL ) then
          call TIME_gettimelabel( timelabel )
          basename = trim(ATMOS_PHY_RD_RESTART_OUT_BASENAME)//'_'//trim(timelabel)
       else
          basename = trim(ATMOS_PHY_RD_RESTART_OUT_BASENAME)
       endif

       LOG_INFO("ATMOS_PHY_RD_vars_restart_create",*) 'basename: ', trim(basename)

       call FILE_CARTESC_create( &
            basename, ATMOS_PHY_RD_RESTART_OUT_TITLE, ATMOS_PHY_RD_RESTART_OUT_DTYPE, & ! [IN]
            restart_fid,                                                              & ! [OUT]
            aggregate=ATMOS_PHY_RD_RESTART_OUT_AGGREGATE                              ) ! [IN]

    endif

    return
  end subroutine ATMOS_PHY_RD_vars_restart_create

  !-----------------------------------------------------------------------------
  !> Exit netCDF define mode
  subroutine ATMOS_PHY_RD_vars_restart_enddef
    use scale_file_cartesC, only: &
       FILE_CARTESC_enddef
    implicit none

    if ( restart_fid /= -1 ) then
       call FILE_CARTESC_enddef( restart_fid ) ! [IN]
    endif

    return
  end subroutine ATMOS_PHY_RD_vars_restart_enddef

  !-----------------------------------------------------------------------------
  !> Close restart file
  subroutine ATMOS_PHY_RD_vars_restart_close
    use scale_file_cartesC, only: &
       FILE_CARTESC_close
    implicit none
    !---------------------------------------------------------------------------

    if ( restart_fid /= -1 ) then
       LOG_NEWLINE
       LOG_INFO("ATMOS_PHY_RD_vars_restart_close",*) 'Close restart file (ATMOS_PHY_RD) '

       call FILE_CARTESC_close( restart_fid ) ! [IN]

       restart_fid = -1
    endif

    return
  end subroutine ATMOS_PHY_RD_vars_restart_close

  !-----------------------------------------------------------------------------
  !> Define variables in restart file
  subroutine ATMOS_PHY_RD_vars_restart_def_var
    use scale_file_cartesC, only: &
       FILE_CARTESC_def_var
    implicit none

    integer :: i
    !---------------------------------------------------------------------------

    if ( restart_fid /= -1 ) then

       do i = 1, VMAX
          call FILE_CARTESC_def_var( restart_fid,                           & ! [IN]
                                     VAR_NAME(i), VAR_DESC(i), VAR_UNIT(i), & ! [IN]
                                     'XY', ATMOS_PHY_RD_RESTART_OUT_DTYPE,  & ! [IN]
                                     VAR_ID(i),                             & ! [OUT]
                                     standard_name=VAR_STDN(i)              ) ! [IN]
       end do
    endif

    return
  end subroutine ATMOS_PHY_RD_vars_restart_def_var

  !-----------------------------------------------------------------------------
  !> Write variables to restart file
  subroutine ATMOS_PHY_RD_vars_restart_write
    use scale_file_cartesC, only: &
       FILE_CARTESC_write_var
    implicit none
    !---------------------------------------------------------------------------

    if ( restart_fid /= -1 ) then

       call ATMOS_PHY_RD_vars_fillhalo

       call ATMOS_PHY_RD_vars_check

       call FILE_CARTESC_write_var( restart_fid, VAR_ID(I_SFLX_LW_up),               & ! [IN]
                                    ATMOS_PHY_RD_SFLX_LW_up(:,:),                    & ! [IN]
                                    VAR_NAME(I_SFLX_LW_up),      'XY'                ) ! [IN]
       call FILE_CARTESC_write_var( restart_fid, VAR_ID(I_SFLX_LW_dn),               & ! [IN]
                                    ATMOS_PHY_RD_SFLX_LW_dn(:,:),                    & ! [IN]
                                    VAR_NAME(I_SFLX_LW_dn),      'XY'                ) ! [IN]
       call FILE_CARTESC_write_var( restart_fid, VAR_ID(I_SFLX_SW_up),               & ! [IN]
                                    ATMOS_PHY_RD_SFLX_SW_up(:,:),                    & ! [IN]
                                    VAR_NAME(I_SFLX_SW_up),      'XY'                ) ! [IN]
       call FILE_CARTESC_write_var( restart_fid, VAR_ID(I_SFLX_SW_dn),               & ! [IN]
                                    ATMOS_PHY_RD_SFLX_SW_dn(:,:),                    & ! [IN]
                                    VAR_NAME(I_SFLX_SW_dn),      'XY'                ) ! [IN]

    endif

    return
  end subroutine ATMOS_PHY_RD_vars_restart_write

  !-----------------------------------------------------------------------------
  subroutine ATMOS_PHY_RD_vars_check
    use scale_statistics, only: &
       STATISTICS_total
    use scale_atmos_grid_cartesC_real, only: &
       ATMOS_GRID_CARTESC_REAL_AREA,   &
       ATMOS_GRID_CARTESC_REAL_TOTAREA
    implicit none
    !---------------------------------------------------------------------------

    call VALCHECK( IA, IS, IE, JA, JS, JE, &
                   ATMOS_PHY_RD_SFLX_LW_up(:,:),             & ! (in)
                   0.0_RP, 1.0E4_RP, VAR_NAME(I_SFLX_LW_up), & ! (in)
                   __FILE__, __LINE__                        ) ! (in)
    call VALCHECK( IA, IS, IE, JA, JS, JE, &
                   ATMOS_PHY_RD_SFLX_LW_dn(:,:),             & ! (in)
                   0.0_RP, 1.0E4_RP, VAR_NAME(I_SFLX_LW_dn), & ! (in)
                   __FILE__, __LINE__                        ) ! (in)
    call VALCHECK( IA, IS, IE, JA, JS, JE, &
                   ATMOS_PHY_RD_SFLX_SW_up(:,:),             & ! (in)
                   0.0_RP, 1.0E4_RP, VAR_NAME(I_SFLX_SW_up), & ! (in)
                   __FILE__, __LINE__                        ) ! (in)
    call VALCHECK( IA, IS, IE, JA, JS, JE, &
                   ATMOS_PHY_RD_SFLX_SW_dn(:,:),             & ! (in)
                   0.0_RP, 1.0E4_RP, VAR_NAME(I_SFLX_SW_dn), & ! (in)
                   __FILE__, __LINE__                        ) ! (in)

    call STATISTICS_total( IA, IS, IE, JA, JS, JE,                          & ! [IN]
                           ATMOS_PHY_RD_SFLX_LW_up(:,:),                    & ! [IN]
                           VAR_NAME(I_SFLX_LW_up),                          & ! [IN]
                           ATMOS_GRID_CARTESC_REAL_AREA(:,:),               & ! [IN]
                           ATMOS_GRID_CARTESC_REAL_TOTAREA                  ) ! [IN]
    call STATISTICS_total( IA, IS, IE, JA, JS, JE,                          & ! [IN]
                           ATMOS_PHY_RD_SFLX_LW_dn(:,:),                    & ! [IN]
                           VAR_NAME(I_SFLX_LW_dn),                          & ! [IN]
                           ATMOS_GRID_CARTESC_REAL_AREA(:,:),               & ! [IN]
                           ATMOS_GRID_CARTESC_REAL_TOTAREA                  ) ! [IN]
    call STATISTICS_total( IA, IS, IE, JA, JS, JE,                          & ! [IN]
                           ATMOS_PHY_RD_SFLX_SW_up(:,:),                    & ! [IN]
                           VAR_NAME(I_SFLX_SW_up),                          & ! [IN]
                           ATMOS_GRID_CARTESC_REAL_AREA(:,:),               & ! [IN]
                           ATMOS_GRID_CARTESC_REAL_TOTAREA                  ) ! [IN]
    call STATISTICS_total( IA, IS, IE, JA, JS, JE,                          & ! [IN]
                           ATMOS_PHY_RD_SFLX_SW_dn(:,:),                    & ! [IN]
                           VAR_NAME(I_SFLX_SW_dn),                          & ! [IN]
                           ATMOS_GRID_CARTESC_REAL_AREA(:,:),               & ! [IN]
                           ATMOS_GRID_CARTESC_REAL_TOTAREA                  ) ! [IN]

    return
  end subroutine ATMOS_PHY_RD_vars_check

end module mod_atmos_phy_rd_vars
