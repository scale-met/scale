!-------------------------------------------------------------------------------
!> module atmosphere / physics / PBL
!!
!! @par Description
!!          Container for mod_atmos_phy_bl
!!
!! @author Team SCALE
!!
!<
!-------------------------------------------------------------------------------
#include "scalelib.h"
module mod_atmos_phy_bl_vars
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
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: ATMOS_PHY_BL_vars_setup
  public :: ATMOS_PHY_BL_vars_finalize
  public :: ATMOS_PHY_BL_vars_fillhalo
  public :: ATMOS_PHY_BL_vars_restart_read
  public :: ATMOS_PHY_BL_vars_restart_write

  public :: ATMOS_PHY_BL_vars_restart_create
  public :: ATMOS_PHY_BL_vars_restart_open
  public :: ATMOS_PHY_BL_vars_restart_def_var
  public :: ATMOS_PHY_BL_vars_restart_enddef
  public :: ATMOS_PHY_BL_vars_restart_close

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  integer,  public              :: QS, QE

  logical,               public :: ATMOS_PHY_BL_RESTART_OUTPUT                 = .false.                !< output restart file?

  character(len=H_LONG),  public :: ATMOS_PHY_BL_RESTART_IN_BASENAME           = ''                     !< Basename of the input  file
  logical,                public :: ATMOS_PHY_BL_RESTART_IN_AGGREGATE                                   !< Switch to use aggregate file
  logical,                public :: ATMOS_PHY_BL_RESTART_IN_POSTFIX_TIMELABEL  = .false.                !< Add timelabel to the basename of input  file?
  character(len=H_LONG),  public :: ATMOS_PHY_BL_RESTART_OUT_BASENAME          = ''                     !< Basename of the output file
  logical,                public :: ATMOS_PHY_BL_RESTART_OUT_AGGREGATE                                  !< Switch to use aggregate file
  logical,                public :: ATMOS_PHY_BL_RESTART_OUT_POSTFIX_TIMELABEL = .true.                 !< Add timelabel to the basename of output file?
  character(len=H_MID),   public :: ATMOS_PHY_BL_RESTART_OUT_TITLE             = 'ATMOS_PHY_BL restart' !< title    of the output file
  character(len=H_SHORT), public :: ATMOS_PHY_BL_RESTART_OUT_DTYPE             = 'DEFAULT'              !< REAL4 or REAL8

  logical,                public :: ATMOS_PHY_BL_MIX_TRACERS                   = .true.

  real(RP), public, allocatable :: ATMOS_PHY_BL_RHOU_t(:,:,:)   ! tendency RHOU [kg/m2/s2]
  real(RP), public, allocatable :: ATMOS_PHY_BL_RHOV_t(:,:,:)   ! tendency RHOV [kg/m2/s2]
  real(RP), public, allocatable :: ATMOS_PHY_BL_RHOT_t(:,:,:)   ! tendency RHOT [K*kg/m3/s]

  real(RP), public, allocatable, target :: ATMOS_PHY_BL_RHOQ_t(:,:,:,:) ! tendency rho*QTRC [kg/kg/s]

  real(RP), public, allocatable :: ATMOS_PHY_BL_Zi       (:,:)  ! depth of the PBL
  real(RP), public, allocatable :: ATMOS_PHY_BL_SFLX_BUOY(:,:)  ! surface flux of buoyancy

  real(RP), public, allocatable :: ATMOS_PHY_BL_QL    (:,:,:)   ! liquid water content in partial condensation
  real(RP), public, allocatable :: ATMOS_PHY_BL_cldfrac(:,:,:)  ! cloud fraction in partial condensation

  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  integer,                private, parameter :: VMAX = 1       !< number of the variables
  integer,                private, parameter :: I_Zi = 1

  character(len=H_SHORT), private            :: VAR_NAME(VMAX) !< name  of the variables
  character(len=H_MID),   private            :: VAR_DESC(VMAX) !< desc. of the variables
  character(len=H_SHORT), private            :: VAR_UNIT(VMAX) !< unit  of the variables
  integer,                private            :: VAR_ID(VMAX)   !< ID    of the variables
  integer,                private            :: restart_fid = -1  ! file ID

  data VAR_NAME / 'PBL_Zi' /

  data VAR_DESC / 'depth of the boundary layer' /

  data VAR_UNIT / 'm' /

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine ATMOS_PHY_BL_vars_setup
    use scale_prc, only: &
       PRC_abort
    use scale_const, only: &
       UNDEF => CONST_UNDEF
    implicit none

    namelist / PARAM_ATMOS_PHY_BL_VARS / &
       ATMOS_PHY_BL_RESTART_IN_BASENAME,           &
       ATMOS_PHY_BL_RESTART_IN_AGGREGATE,          &
       ATMOS_PHY_BL_RESTART_IN_POSTFIX_TIMELABEL,  &
       ATMOS_PHY_BL_RESTART_OUTPUT,                &
       ATMOS_PHY_BL_RESTART_OUT_BASENAME,          &
       ATMOS_PHY_BL_RESTART_OUT_AGGREGATE,         &
       ATMOS_PHY_BL_RESTART_OUT_POSTFIX_TIMELABEL, &
       ATMOS_PHY_BL_RESTART_OUT_TITLE,             &
       ATMOS_PHY_BL_RESTART_OUT_DTYPE,             &
       ATMOS_PHY_BL_MIX_TRACERS

    integer :: ierr
    integer :: iv
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("ATMOS_PHY_BL_vars_setup",*) 'Setup'

    allocate( ATMOS_PHY_BL_RHOU_t(KA,IA,JA)    )
    allocate( ATMOS_PHY_BL_RHOV_t(KA,IA,JA)    )
    allocate( ATMOS_PHY_BL_RHOT_t(KA,IA,JA)    )
    allocate( ATMOS_PHY_BL_RHOQ_t(KA,IA,JA,QA) )
    ATMOS_PHY_BL_RHOU_t(:,:,:)   = UNDEF
    ATMOS_PHY_BL_RHOV_t(:,:,:)   = UNDEF
    ATMOS_PHY_BL_RHOT_t(:,:,:)   = UNDEF
    ATMOS_PHY_BL_RHOQ_t(:,:,:,:) = UNDEF
    !$acc enter data create(ATMOS_PHY_BL_RHOU_t,ATMOS_PHY_BL_RHOV_t,ATMOS_PHY_BL_RHOT_t,ATMOS_PHY_BL_RHOQ_t)

    allocate( ATMOS_PHY_BL_Zi      (IA,JA) )
    ATMOS_PHY_BL_Zi      (:,:) = UNDEF
    !$acc enter data create(ATMOS_PHY_BL_Zi)

    allocate( ATMOS_PHY_BL_SFLX_BUOY(IA,JA) )
    ATMOS_PHY_BL_SFLX_BUOY(:,:) = UNDEF
    !$acc enter data create(ATMOS_PHY_BL_SFLX_BUOY)

    allocate( ATMOS_PHY_BL_QL     (KA,IA,JA) )
    allocate( ATMOS_PHY_BL_cldfrac(KA,IA,JA) )
    ATMOS_PHY_BL_QL     (:,:,:) = UNDEF
    ATMOS_PHY_BL_cldfrac(:,:,:) = UNDEF
    !$acc enter data create(ATMOS_PHY_BL_QL,ATMOS_PHY_BL_cldfrac)

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_ATMOS_PHY_BL_VARS,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("ATMOS_PHY_BL_vars_setup",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("ATMOS_PHY_BL_vars_setup",*) 'Not appropriate names in namelist PARAM_ATMOS_PHY_BL_VARS. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_ATMOS_PHY_BL_VARS)

    LOG_NEWLINE
    LOG_INFO("ATMOS_PHY_BL_vars_setup",*) '[ATMOS_PHY_BL] prognostic/diagnostic variables'
    LOG_INFO_CONT('(1x,A,A24,A,A48,A,A12,A)') &
               '      |', 'VARNAME                 ','|', &
               'DESCRIPTION                                     ', '[', 'UNIT        ', ']'
    do iv = 1, VMAX
       LOG_INFO_CONT('(1x,A,I3,A,A24,A,A48,A,A12,A)') &
                  'NO.',iv,'|',VAR_NAME(iv),'|',VAR_DESC(iv),'[',VAR_UNIT(iv),']'
    enddo

    LOG_NEWLINE
    if ( ATMOS_PHY_BL_RESTART_IN_BASENAME /= '' ) then
       LOG_INFO("ATMOS_PHY_BL_vars_setup",*) 'Restart input?  : YES, file = ', trim(ATMOS_PHY_BL_RESTART_IN_BASENAME)
       LOG_INFO("ATMOS_PHY_BL_vars_setup",*) 'Add timelabel?  : ', ATMOS_PHY_BL_RESTART_IN_POSTFIX_TIMELABEL
    else
       LOG_INFO("ATMOS_PHY_BL_vars_setup",*) 'Restart input?  : NO'
    endif
    if (       ATMOS_PHY_BL_RESTART_OUTPUT             &
         .AND. ATMOS_PHY_BL_RESTART_OUT_BASENAME /= '' ) then
       LOG_INFO("ATMOS_PHY_BL_vars_setup",*) 'Restart output? : YES, file = ', trim(ATMOS_PHY_BL_RESTART_OUT_BASENAME)
       LOG_INFO("ATMOS_PHY_BL_vars_setup",*) 'Add timelabel?  : ', ATMOS_PHY_BL_RESTART_OUT_POSTFIX_TIMELABEL
    else
       LOG_INFO("ATMOS_PHY_BL_vars_setup",*) 'Restart output? : NO'
       ATMOS_PHY_BL_RESTART_OUTPUT = .false.
    endif

    return
  end subroutine ATMOS_PHY_BL_vars_setup

  !-----------------------------------------------------------------------------
  !> Finalize
  subroutine ATMOS_PHY_BL_vars_finalize
    implicit none
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("ATMOS_PHY_BL_vars_finalize",*) 'Finalize'

    !$acc exit data delete(ATMOS_PHY_BL_RHOU_t,ATMOS_PHY_BL_RHOV_t,ATMOS_PHY_BL_RHOT_t,ATMOS_PHY_BL_RHOQ_t)
    deallocate( ATMOS_PHY_BL_RHOU_t )
    deallocate( ATMOS_PHY_BL_RHOV_t )
    deallocate( ATMOS_PHY_BL_RHOT_t )
    deallocate( ATMOS_PHY_BL_RHOQ_t )

    !$acc exit data delete(ATMOS_PHY_BL_Zi)
    deallocate( ATMOS_PHY_BL_Zi )

    !$acc exit data delete(ATMOS_PHY_BL_SFLX_BUOY)
    deallocate( ATMOS_PHY_BL_SFLX_BUOY )

    !$acc exit data delete(ATMOS_PHY_BL_QL,ATMOS_PHY_BL_cldfrac)
    deallocate( ATMOS_PHY_BL_QL      )
    deallocate( ATMOS_PHY_BL_cldfrac )

    return
  end subroutine ATMOS_PHY_BL_vars_finalize

  !-----------------------------------------------------------------------------
  !> HALO Communication
  subroutine ATMOS_PHY_BL_vars_fillhalo
    use scale_comm_cartesC, only: &
       COMM_vars8, &
       COMM_wait
    implicit none
    !---------------------------------------------------------------------------

    call COMM_vars8( ATMOS_PHY_BL_Zi(:,:), 1 )
    call COMM_wait ( ATMOS_PHY_BL_Zi(:,:), 1 )

    return
  end subroutine ATMOS_PHY_BL_vars_fillhalo

  !-----------------------------------------------------------------------------
  !> Open restart file for read
  subroutine ATMOS_PHY_BL_vars_restart_open
    use scale_time, only: &
       TIME_gettimelabel
    use scale_file_cartesC, only: &
       FILE_CARTESC_open
    implicit none

    character(len=19)     :: timelabel
    character(len=H_LONG) :: basename
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("ATMOS_PHY_BL_vars_restart_open",*) 'Open restart file (ATMOS_PHY_BL) '

    if ( ATMOS_PHY_BL_RESTART_IN_BASENAME /= '' ) then

       if ( ATMOS_PHY_BL_RESTART_IN_POSTFIX_TIMELABEL ) then
          call TIME_gettimelabel( timelabel )
          basename = trim(ATMOS_PHY_BL_RESTART_IN_BASENAME)//'_'//trim(timelabel)
       else
          basename = trim(ATMOS_PHY_BL_RESTART_IN_BASENAME)
       endif

       LOG_INFO("ATMOS_PHY_BL_vars_restart_open",*) 'basename: ', trim(basename)

       call FILE_CARTESC_open( basename, restart_fid, aggregate=ATMOS_PHY_BL_RESTART_IN_AGGREGATE )
    else
       LOG_INFO("ATMOS_PHY_BL_vars_restart_open",*) 'restart file for ATMOS_PHY_BL is not specified.'
    endif

    return
  end subroutine ATMOS_PHY_BL_vars_restart_open

  !-----------------------------------------------------------------------------
  !> Read restart
  subroutine ATMOS_PHY_BL_vars_restart_read
    use scale_file, only: &
       FILE_get_aggregate
    use scale_file_cartesC, only: &
       FILE_CARTESC_read, &
       FILE_CARTESC_flush
    implicit none

    integer  :: i, j
    !---------------------------------------------------------------------------

    if ( restart_fid /= -1 ) then
       LOG_NEWLINE
       LOG_INFO("ATMOS_PHY_BL_vars_restart_read",*) 'Read from restart file (ATMOS_PHY_BL) '

       call FILE_CARTESC_read( restart_fid, VAR_NAME(1), 'XY', & ! [IN]
                               ATMOS_PHY_BL_Zi(:,:)            ) ! [OUT]

       if ( FILE_get_AGGREGATE(restart_fid) ) then
          call FILE_CARTESC_flush( restart_fid ) ! X/Y halos have been read from file
          !$acc update device(ATMOS_PHY_BL_Zi)
       else
          call ATMOS_PHY_BL_vars_fillhalo
       end if

       call ATMOS_PHY_BL_vars_check

    else
       LOG_INFO("ATMOS_PHY_BL_vars_restart_read",*) 'invalid restart file ID for ATMOS_PHY_BL.'
    endif

    return
  end subroutine ATMOS_PHY_BL_vars_restart_read

  !-----------------------------------------------------------------------------
  !> Create restart file
  subroutine ATMOS_PHY_BL_vars_restart_create
    use scale_time, only: &
       TIME_gettimelabel
    use scale_file_cartesC, only: &
       FILE_CARTESC_create
    implicit none

    character(len=19)     :: timelabel
    character(len=H_LONG) :: basename
    !---------------------------------------------------------------------------

    if ( ATMOS_PHY_BL_RESTART_OUT_BASENAME /= '' ) then

       LOG_NEWLINE
       LOG_INFO("ATMOS_PHY_BL_vars_restart_create",*) 'Create restart file (ATMOS_PHY_AE) '

       if ( ATMOS_PHY_BL_RESTART_OUT_POSTFIX_TIMELABEL ) then
          call TIME_gettimelabel( timelabel )
          basename = trim(ATMOS_PHY_BL_RESTART_OUT_BASENAME)//'_'//trim(timelabel)
       else
          basename = trim(ATMOS_PHY_BL_RESTART_OUT_BASENAME)
       endif

       LOG_INFO("ATMOS_PHY_BL_vars_restart_create",*) 'basename: ', trim(basename)

    call FILE_CARTESC_create( &
         basename, ATMOS_PHY_BL_RESTART_OUT_TITLE, ATMOS_PHY_BL_RESTART_OUT_DTYPE, & ! [IN]
         restart_fid,                                                              & ! [OUT]
         aggregate=ATMOS_PHY_BL_RESTART_OUT_AGGREGATE                              ) ! [IN]
    endif

    return
  end subroutine ATMOS_PHY_BL_vars_restart_create

  !-----------------------------------------------------------------------------
  !> Exit netCDF define mode
  subroutine ATMOS_PHY_BL_vars_restart_enddef
    use scale_file_cartesC, only: &
       FILE_CARTESC_enddef
    implicit none

    if ( restart_fid /= -1 ) then
       call FILE_CARTESC_enddef( restart_fid ) ! [IN]
    endif

    return
  end subroutine ATMOS_PHY_BL_vars_restart_enddef

  !-----------------------------------------------------------------------------
  !> Close restart file
  subroutine ATMOS_PHY_BL_vars_restart_close
    use scale_file_cartesC, only: &
       FILE_CARTESC_close
    implicit none
    !---------------------------------------------------------------------------

    if ( restart_fid /= -1 ) then
       LOG_NEWLINE
       LOG_INFO("ATMOS_PHY_BL_vars_restart_close",*) 'Close restart file (ATMOS_PHY_BL) '

       call FILE_CARTESC_close( restart_fid ) ! [IN]

       restart_fid = -1
    endif

    return
  end subroutine ATMOS_PHY_BL_vars_restart_close

  !-----------------------------------------------------------------------------
  !> Write restart
  subroutine ATMOS_PHY_BL_vars_restart_def_var
    use scale_file_cartesC, only: &
       FILE_CARTESC_def_var
    implicit none
    integer :: iv
    !---------------------------------------------------------------------------

    if ( restart_fid /= -1 ) then

       do iv = 1, VMAX
          call FILE_CARTESC_def_var( restart_fid,        & ! [IN]
               VAR_NAME(iv), VAR_DESC(iv), VAR_UNIT(iv), & ! [IN]
               'XY', ATMOS_PHY_BL_RESTART_OUT_DTYPE,     & ! [IN]
               VAR_ID(iv)                                ) ! [OUT]
       end do

    endif

    return
  end subroutine ATMOS_PHY_BL_vars_restart_def_var

  !-----------------------------------------------------------------------------
  !> Write restart
  subroutine ATMOS_PHY_BL_vars_restart_write
    use scale_file_cartesC, only: &
       FILE_CARTESC_write_var
    implicit none

    !---------------------------------------------------------------------------

    if ( restart_fid /= -1 ) then

       call ATMOS_PHY_BL_vars_fillhalo

       call ATMOS_PHY_BL_vars_check

       call FILE_CARTESC_write_var( restart_fid, VAR_ID(1), ATMOS_PHY_BL_Zi(:,:), &
                                    VAR_NAME(1), 'XY' ) ! [IN]

    endif

    return
  end subroutine ATMOS_PHY_BL_vars_restart_write

  subroutine ATMOS_PHY_BL_vars_check
    use scale_statistics, only: &
       STATISTICS_total
    use scale_atmos_grid_cartesC_real, only: &
       ATMOS_GRID_CARTESC_REAL_AREA, &
       ATMOS_GRID_CARTESC_REAL_TOTAREA
    implicit none

    call VALCHECK( IA, IS, IE, JA, JS, JE, &
                   ATMOS_PHY_BL_Zi(:,:),          & ! (in)
                   0.0_RP, 1.0E5_RP, VAR_NAME(1), & ! (in)
                   __FILE__, __LINE__             ) ! (in)

    call STATISTICS_total( IA, IS, IE, JA, JS, JE, &
                           ATMOS_PHY_BL_Zi(:,:), VAR_NAME(1), &
                           ATMOS_GRID_CARTESC_REAL_AREA(:,:), & ! (in)
                           ATMOS_GRID_CARTESC_REAL_TOTAREA    ) ! (in)

    return
  end subroutine ATMOS_PHY_BL_vars_check
  
end module mod_atmos_phy_bl_vars
