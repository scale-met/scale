!-------------------------------------------------------------------------------
!> module Atmosphere / Physics Turbulence
!!
!! @par Description
!!          Container for mod_atmos_phy_tb
!!
!! @author Team SCALE
!!
!<
!-------------------------------------------------------------------------------
#include "scalelib.h"
module mod_atmos_phy_tb_vars
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
  public :: ATMOS_PHY_TB_vars_setup
  public :: ATMOS_PHY_TB_vars_finalize
  public :: ATMOS_PHY_TB_vars_fillhalo
  public :: ATMOS_PHY_TB_vars_restart_read
  public :: ATMOS_PHY_TB_vars_restart_write

  public :: ATMOS_PHY_TB_vars_restart_create
  public :: ATMOS_PHY_TB_vars_restart_open
  public :: ATMOS_PHY_TB_vars_restart_def_var
  public :: ATMOS_PHY_TB_vars_restart_enddef
  public :: ATMOS_PHY_TB_vars_restart_close

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  logical,               public :: ATMOS_PHY_TB_RESTART_OUTPUT                 = .false.                !< output restart file?

  character(len=H_LONG),  public :: ATMOS_PHY_TB_RESTART_IN_BASENAME           = ''                     !< Basename of the input  file
  logical,                public :: ATMOS_PHY_TB_RESTART_IN_AGGREGATE                                   !< Switch to use aggregate file
  logical,                public :: ATMOS_PHY_TB_RESTART_IN_POSTFIX_TIMELABEL  = .false.                !< Add timelabel to the basename of input  file?
  character(len=H_LONG),  public :: ATMOS_PHY_TB_RESTART_OUT_BASENAME          = ''                     !< Basename of the output file
  logical,                public :: ATMOS_PHY_TB_RESTART_OUT_AGGREGATE                                  !< Switch to use aggregate file
  logical,                public :: ATMOS_PHY_TB_RESTART_OUT_POSTFIX_TIMELABEL = .true.                 !< Add timelabel to the basename of output file?
  character(len=H_MID),   public :: ATMOS_PHY_TB_RESTART_OUT_TITLE             = 'ATMOS_PHY_TB restart' !< title    of the output file
  character(len=H_SHORT), public :: ATMOS_PHY_TB_RESTART_OUT_DTYPE             = 'DEFAULT'              !< REAL4 or REAL8

  real(RP), public, allocatable :: ATMOS_PHY_TB_MOMZ_t(:,:,:)   ! tendency MOMZ [kg/m2/s2]
  real(RP), public, allocatable :: ATMOS_PHY_TB_MOMX_t(:,:,:)   ! tendency MOMX [kg/m2/s2]
  real(RP), public, allocatable :: ATMOS_PHY_TB_MOMY_t(:,:,:)   ! tendency MOMY [kg/m2/s2]
  real(RP), public, allocatable :: ATMOS_PHY_TB_RHOT_t(:,:,:)   ! tendency RHOT [K*kg/m3/s]
  real(RP), public, allocatable :: ATMOS_PHY_TB_RHOQ_t(:,:,:,:) ! tendency rho*QTRC [kg/kg/s]

  integer, public :: I_TKE = -1

  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
!  integer,                private, parameter :: VMAX = 0       !< number of the variables
!  character(len=H_SHORT), private            :: VAR_NAME(VMAX) !< name  of the variables
!  character(len=H_MID),   private            :: VAR_DESC(VMAX) !< desc. of the variables
!  character(len=H_SHORT), private            :: VAR_UNIT(VMAX) !< unit  of the variables
!  integer,                private            :: VAR_ID(VMAX)   !< ID    of the variables
!  integer,                private            :: restart_fid = -1  ! file ID

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine ATMOS_PHY_TB_vars_setup
    use scale_prc, only: &
       PRC_abort
    use scale_const, only: &
       UNDEF => CONST_UNDEF
    implicit none

    namelist / PARAM_ATMOS_PHY_TB_VARS / &
       ATMOS_PHY_TB_RESTART_IN_BASENAME,           &
       ATMOS_PHY_TB_RESTART_IN_AGGREGATE,          &
       ATMOS_PHY_TB_RESTART_IN_POSTFIX_TIMELABEL,  &
       ATMOS_PHY_TB_RESTART_OUTPUT,                &
       ATMOS_PHY_TB_RESTART_OUT_BASENAME,          &
       ATMOS_PHY_TB_RESTART_OUT_AGGREGATE,         &
       ATMOS_PHY_TB_RESTART_OUT_POSTFIX_TIMELABEL, &
       ATMOS_PHY_TB_RESTART_OUT_TITLE,             &
       ATMOS_PHY_TB_RESTART_OUT_DTYPE

    integer :: ierr
    integer :: iv
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("ATMOS_PHY_TB_vars_setup",*) 'Setup'

    allocate( ATMOS_PHY_TB_MOMZ_t(KA,IA,JA)    )
    allocate( ATMOS_PHY_TB_MOMX_t(KA,IA,JA)    )
    allocate( ATMOS_PHY_TB_MOMY_t(KA,IA,JA)    )
    allocate( ATMOS_PHY_TB_RHOT_t(KA,IA,JA)    )
    allocate( ATMOS_PHY_TB_RHOQ_t(KA,IA,JA,QA) )
    ATMOS_PHY_TB_MOMZ_t(:,:,:)   = UNDEF
    ATMOS_PHY_TB_MOMX_t(:,:,:)   = UNDEF
    ATMOS_PHY_TB_MOMY_t(:,:,:)   = UNDEF
    ATMOS_PHY_TB_RHOT_t(:,:,:)   = UNDEF
    ATMOS_PHY_TB_RHOQ_t(:,:,:,:) = UNDEF
    !$acc enter data create(ATMOS_PHY_TB_MOMZ_t, ATMOS_PHY_TB_MOMX_t, ATMOS_PHY_TB_MOMY_t, ATMOS_PHY_TB_RHOT_t, ATMOS_PHY_TB_RHOQ_t)

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_ATMOS_PHY_TB_VARS,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("ATMOS_PHY_TB_vars_setup",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("ATMOS_PHY_TB_vars_setup",*) 'Not appropriate names in namelist PARAM_ATMOS_PHY_TB_VARS. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_ATMOS_PHY_TB_VARS)

!    LOG_NEWLINE
!    LOG_INFO("ATMOS_PHY_TB_vars_setup",*) '[ATMOS_PHY_TB] prognostic/diagnostic variables'
!    LOG_INFO_CONT('(1x,A,A24,A,A48,A,A12,A)') &
!               '      |', 'VARNAME                 ','|', &
!               'DESCRIPTION                                     ', '[', 'UNIT        ', ']'
!    do iv = 1, VMAX
!       LOG_INFO_CONT('(1x,A,I3,A,A24,A,A48,A,A12,A)') &
!                  'NO.',iv,'|',VAR_NAME(iv),'|',VAR_DESC(iv),'[',VAR_UNIT(iv),']'
!    enddo

!    LOG_NEWLINE
!    if ( ATMOS_PHY_TB_RESTART_IN_BASENAME /= '' ) then
!       LOG_INFO("ATMOS_PHY_TB_vars_setup",*) 'Restart input?  : YES, file = ', trim(ATMOS_PHY_TB_RESTART_IN_BASENAME)
!       LOG_INFO("ATMOS_PHY_TB_vars_setup",*) 'Add timelabel?  : ', ATMOS_PHY_TB_RESTART_IN_POSTFIX_TIMELABEL
!    else
!       LOG_INFO("ATMOS_PHY_TB_vars_setup",*) 'Restart input?  : NO'
!    endif
!    if (       ATMOS_PHY_TB_RESTART_OUTPUT             &
!         .AND. ATMOS_PHY_TB_RESTART_OUT_BASENAME /= '' ) then
!       LOG_INFO("ATMOS_PHY_TB_vars_setup",*) 'Restart output? : YES, file = ', trim(ATMOS_PHY_TB_RESTART_OUT_BASENAME)
!       LOG_INFO("ATMOS_PHY_TB_vars_setup",*) 'Add timelabel?  : ', ATMOS_PHY_TB_RESTART_OUT_POSTFIX_TIMELABEL
!    else
!       LOG_INFO("ATMOS_PHY_TB_vars_setup",*) 'Restart output? : NO'
!       ATMOS_PHY_TB_RESTART_OUTPUT = .false.
!    endif

    return
  end subroutine ATMOS_PHY_TB_vars_setup

  !-----------------------------------------------------------------------------
  !> Finalize
  subroutine ATMOS_PHY_TB_vars_finalize
    implicit none
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("ATMOS_PHY_TB_vars_finalize",*) 'Finalize'

    !$acc exit data delete(ATMOS_PHY_TB_MOMZ_t, ATMOS_PHY_TB_MOMX_t, ATMOS_PHY_TB_MOMY_t, ATMOS_PHY_TB_RHOT_t, ATMOS_PHY_TB_RHOQ_t)
    deallocate( ATMOS_PHY_TB_MOMZ_t )
    deallocate( ATMOS_PHY_TB_MOMX_t )
    deallocate( ATMOS_PHY_TB_MOMY_t )
    deallocate( ATMOS_PHY_TB_RHOT_t )
    deallocate( ATMOS_PHY_TB_RHOQ_t )

    return
  end subroutine ATMOS_PHY_TB_vars_finalize

  !-----------------------------------------------------------------------------
  !> HALO Communication
  subroutine ATMOS_PHY_TB_vars_fillhalo
    use scale_comm_cartesC, only: &
       COMM_vars8, &
       COMM_wait
    implicit none
    !---------------------------------------------------------------------------

!    do j  = 1, JA
!    do i  = 1, IA
!       ATMOS_PHY_TB_??(   1:KS-1,i,j) = ATMOS_PHY_TB_??(KS,i,j)
!       ATMOS_PHY_TB_??(KE+1:KA,  i,j) = ATMOS_PHY_TB_??(KE,i,j)
!    enddo
!    enddo

!    call COMM_vars8( ATMOS_PHY_TB_??(:,:,:), 1 )
!    call COMM_wait ( ATMOS_PHY_TB_??(:,:,:), 1 )

    return
  end subroutine ATMOS_PHY_TB_vars_fillhalo

  !-----------------------------------------------------------------------------
  !> Open restart file for read
  subroutine ATMOS_PHY_TB_vars_restart_open
    use scale_time, only: &
       TIME_gettimelabel
    use scale_file_cartesC, only: &
       FILE_CARTESC_open
    implicit none

    character(len=19)     :: timelabel
    character(len=H_LONG) :: basename
    !---------------------------------------------------------------------------

!    LOG_NEWLINE
!    LOG_INFO("ATMOS_PHY_TB_vars_restart_open",*) 'Open restart file (ATMOS_PHY_TB) '

!    if ( ATMOS_PHY_TB_RESTART_IN_BASENAME /= '' ) then

!       if ( ATMOS_PHY_TB_RESTART_IN_POSTFIX_TIMELABEL ) then
!          call TIME_gettimelabel( timelabel )
!          basename = trim(ATMOS_PHY_TB_RESTART_IN_BASENAME)//'_'//trim(timelabel)
!       else
!          basename = trim(ATMOS_PHY_TB_RESTART_IN_BASENAME)
!       endif

!       LOG_INFO("ATMOS_PHY_TB_vars_restart_open",*) 'basename: ', trim(basename)

!       call FILE_CARTESC_open( basename, restart_fid, aggregate=ATMOS_PHY_TB_RESTART_IN_AGGREGATE )
!    else
!       LOG_INFO("ATMOS_PHY_TB_vars_restart_open",*) 'restart file for ATMOS_PHY_TB is not specified.'
!    endif

    return
  end subroutine ATMOS_PHY_TB_vars_restart_open

  !-----------------------------------------------------------------------------
  !> Read restart
  subroutine ATMOS_PHY_TB_vars_restart_read
!!$    use scale_file, only: &
!!$       FILE_get_aggregate
!!$    use scale_file_cartesC, only: &
!!$       FILE_CARTESC_read, &
!!$       FILE_CARTESC_flush
!!$    implicit none

    integer  :: i, j
    !---------------------------------------------------------------------------

!    if ( restart_fid /= -1 ) then
!       LOG_NEWLINE
!       LOG_INFO("ATMOS_PHY_TB_vars_restart_read",*) 'Read from restart file (ATMOS_PHY_TB) '
!
!       call FILE_CARTESC_read( restart_fid, VAR_NAME(1), 'ZXY', & ! [IN]
!                               ATMOS_PHY_TB_??(:,:,:)           ) ! [OUT]
!
!       if ( FILE_get_AGGREGATE(restart_fid) ) then
!          call FILE_CARTESC_flush( restart_fid ) ! X/Y halos have been read from file
!
!          ! fill k halos
!          do j  = 1, JA
!          do i  = 1, IA
!             ATMOS_PHY_TB_??(   1:KS-1,i,j) = ATMOS_PHY_TB_??(KS,i,j)
!             ATMOS_PHY_TB_??(KE+1:KA,  i,j) = ATMOS_PHY_TB_??(KE,i,j)
!          enddo
!          enddo
!       else
!          call ATMOS_PHY_TB_vars_fillhalo
!       end if
!
!       call ATMOS_PHY_TB_vars_check
!
!    else
!       LOG_INFO("ATMOS_PHY_TB_vars_restart_read",*) 'invalid restart file ID for ATMOS_PHY_TB.'
!    endif

    return
  end subroutine ATMOS_PHY_TB_vars_restart_read

  !-----------------------------------------------------------------------------
  !> Create restart file
  subroutine ATMOS_PHY_TB_vars_restart_create
    use scale_time, only: &
       TIME_gettimelabel
    use scale_file_cartesC, only: &
       FILE_CARTESC_create
    implicit none

    character(len=19)     :: timelabel
    character(len=H_LONG) :: basename
    !---------------------------------------------------------------------------

!    if ( ATMOS_PHY_TB_RESTART_OUT_BASENAME /= '' ) then
!
!       LOG_NEWLINE
!       LOG_INFO("ATMOS_PHY_TB_vars_restart_create",*) 'Create restart file (ATMOS_PHY_AE) '
!
!       if ( ATMOS_PHY_TB_RESTART_OUT_POSTFIX_TIMELABEL ) then
!          call TIME_gettimelabel( timelabel )
!          basename = trim(ATMOS_PHY_TB_RESTART_OUT_BASENAME)//'_'//trim(timelabel)
!       else
!          basename = trim(ATMOS_PHY_TB_RESTART_OUT_BASENAME)
!       endif
!
!       LOG_INFO("ATMOS_PHY_TB_vars_restart_create",*) 'basename: ', trim(basename)
!
!    call FILE_CARTESC_create( &
!         basename, ATMOS_PHY_TB_RESTART_OUT_TITLE, ATMOS_PHY_TB_RESTART_OUT_DTYPE, & ! [IN]
!         restart_fid,                                                              & ! [OUT]
!         aggregate=ATMOS_PHY_TB_RESTART_OUT_AGGREGATE                              ) ! [IN]
!    endif

    return
  end subroutine ATMOS_PHY_TB_vars_restart_create

  !-----------------------------------------------------------------------------
  !> Exit netCDF define mode
  subroutine ATMOS_PHY_TB_vars_restart_enddef
    use scale_file_cartesC, only: &
       FILE_CARTESC_enddef
    implicit none

!    if ( restart_fid /= -1 ) then
!       call FILE_CARTESC_enddef( restart_fid ) ! [IN]
!    endif

    return
  end subroutine ATMOS_PHY_TB_vars_restart_enddef

  !-----------------------------------------------------------------------------
  !> Close restart file
  subroutine ATMOS_PHY_TB_vars_restart_close
    use scale_file_cartesC, only: &
       FILE_CARTESC_close
    implicit none
    !---------------------------------------------------------------------------

!    if ( restart_fid /= -1 ) then
!       LOG_NEWLINE
!       LOG_INFO("ATMOS_PHY_TB_vars_restart_close",*) 'Close restart file (ATMOS_PHY_TB) '
!
!       call FILE_CARTESC_close( restart_fid ) ! [IN]
!
!       restart_fid = -1
!    endif

    return
  end subroutine ATMOS_PHY_TB_vars_restart_close

  !-----------------------------------------------------------------------------
  !> Write restart
  subroutine ATMOS_PHY_TB_vars_restart_def_var
    use scale_file_cartesC, only: &
       FILE_CARTESC_def_var
    implicit none
    !---------------------------------------------------------------------------

!    if ( restart_fid /= -1 ) then

!       call FILE_CARTESC_def_var( restart_fid,     & ! [IN]
!            VAR_NAME(1), VAR_DESC(1), VAR_UNIT(1), & ! [IN]
!            'ZXY', ATMOS_PHY_TB_RESTART_OUT_DTYPE, & ! [IN]
!            VAR_ID(1)                              ) ! [OUT]

!    endif

    return
  end subroutine ATMOS_PHY_TB_vars_restart_def_var

  !-----------------------------------------------------------------------------
  !> Write restart
  subroutine ATMOS_PHY_TB_vars_restart_write
    use scale_file_cartesC, only: &
       FILE_CARTESC_write_var
    implicit none

    !---------------------------------------------------------------------------

!    if ( restart_fid /= -1 ) then
!
!       call ATMOS_PHY_TB_vars_fillhalo
!
!       call ATMOS_PHY_TB_vars_check
!
!       call FILE_CARTESC_write_var( restart_fid, VAR_ID(1), ATMOS_PHY_TB_??(:,:,:), &
!                              VAR_NAME(1), 'ZXY' ) ! [IN]
!
!    endif

    return
  end subroutine ATMOS_PHY_TB_vars_restart_write

  subroutine ATMOS_PHY_TB_vars_check
    use scale_statistics, only: &
       STATISTICS_total
    use scale_atmos_grid_cartesC_real, only: &
       ATMOS_GRID_CARTESC_REAL_VOL, &
       ATMOS_GRID_CARTESC_REAL_TOTVOL
    implicit none

!!$    call VALCHECK( KA, KS, KE, IA, IS, IE, JA, JS, JE, &
!!$                   ATMOS_PHY_TB_??(:,:),        & ! (in)
!!$                   0.0_RP, 0.0_RP, VAR_NAME(1), & ! (in)
!!$                   __FILE__, __LINE__           ) ! (in)
!!$
!!$    call STATISTICS_total( KA, KS, KE, IA, IS, IE, JA, JS, JE, &
!!$                           ATMOS_PHY_TB_??(:,:,:), VAR_NAME(1), &
!!$                           ATMOS_GRID_CARTESC_REAL_VOL(:,:,:),  & ! (in)
!!$                           ATMOS_GRID_CARTESC_REAL_TOTVOL       ) ! (in)

    return
  end subroutine ATMOS_PHY_TB_vars_check
  
end module mod_atmos_phy_tb_vars
