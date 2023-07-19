!-------------------------------------------------------------------------------
!> module Atmosphere / Dynamics
!!
!! @par Description
!!          Container for mod_atmos_dyn
!!
!! @author Team SCALE
!!
!<
!-------------------------------------------------------------------------------
#include "scalelib.h"
module mod_atmos_dyn_vars
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_io
  use scale_prof
  use scale_debug
  use scale_atmos_grid_cartesC_index
  use scale_index
  use scale_tracer
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: ATMOS_DYN_vars_setup
  public :: ATMOS_DYN_vars_finalize
  public :: ATMOS_DYN_vars_fillhalo
  public :: ATMOS_DYN_vars_restart_read
  public :: ATMOS_DYN_vars_restart_write

  public :: ATMOS_DYN_vars_restart_create
  public :: ATMOS_DYN_vars_restart_open
  public :: ATMOS_DYN_vars_restart_def_var
  public :: ATMOS_DYN_vars_restart_enddef
  public :: ATMOS_DYN_vars_restart_close

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  logical,               public :: ATMOS_DYN_RESTART_OUTPUT                 = .false.             !< output restart file?

  character(len=H_LONG),  public :: ATMOS_DYN_RESTART_IN_BASENAME           = ''                  !< Basename of the input  file
  logical,                public :: ATMOS_DYN_RESTART_IN_AGGREGATE                                !< Switch to use aggregate file
  logical,                public :: ATMOS_DYN_RESTART_IN_POSTFIX_TIMELABEL  = .false.             !< Add timelabel to the basename of input  file?
  character(len=H_LONG),  public :: ATMOS_DYN_RESTART_OUT_BASENAME          = ''                  !< Basename of the output file
  logical,                public :: ATMOS_DYN_RESTART_OUT_AGGREGATE                               !< Switch to use aggregate file
  logical,                public :: ATMOS_DYN_RESTART_OUT_POSTFIX_TIMELABEL = .true.              !< Add timelabel to the basename of output file?
  character(len=H_MID),   public :: ATMOS_DYN_RESTART_OUT_TITLE             = 'ATMOS_DYN restart' !< title    of the output file
  character(len=H_SHORT), public :: ATMOS_DYN_RESTART_OUT_DTYPE             = 'DEFAULT'           !< REAL4 or REAL8

  ! prognostic variables
  real(RP), public, allocatable :: PROG(:,:,:,:)

  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  integer,                private, parameter :: VMAX = 100     !< max of number of the variables
  character(len=H_SHORT), private            :: VAR_NAME(VMAX) !< name  of the variables
  character(len=H_MID),   private            :: VAR_DESC(VMAX) !< desc. of the variables
  character(len=H_SHORT), private            :: VAR_UNIT(VMAX) !< unit  of the variables
  integer,                private            :: VAR_ID(VMAX)   !< ID    of the variables
  integer,                private            :: restart_fid = -1  ! file ID

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine ATMOS_DYN_vars_setup
    use scale_prc, only: &
       PRC_abort
    use scale_const, only: &
       UNDEF => CONST_UNDEF
    use mod_atmos_admin, only: &
       ATMOS_DYN_TYPE
    use scale_atmos_dyn_tstep_short, only: &
       ATMOS_DYN_Tstep_short_regist
    implicit none

    namelist / PARAM_ATMOS_DYN_VARS / &
       ATMOS_DYN_RESTART_IN_BASENAME,           &
       ATMOS_DYN_RESTART_IN_AGGREGATE,          &
       ATMOS_DYN_RESTART_IN_POSTFIX_TIMELABEL,  &
       ATMOS_DYN_RESTART_OUTPUT,                &
       ATMOS_DYN_RESTART_OUT_BASENAME,          &
       ATMOS_DYN_RESTART_OUT_AGGREGATE,         &
       ATMOS_DYN_RESTART_OUT_POSTFIX_TIMELABEL, &
       ATMOS_DYN_RESTART_OUT_TITLE,             &
       ATMOS_DYN_RESTART_OUT_DTYPE

    integer :: ierr
    integer :: iv
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("ATMOS_DYN_vars_setup",*) 'Setup'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_ATMOS_DYN_VARS,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("ATMOS_DYN_vars_setup",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("ATMOS_DYN_vars_setup",*) 'Not appropriate names in namelist PARAM_ATMOS_DYN_VARS. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_ATMOS_DYN_VARS)

    call ATMOS_DYN_Tstep_short_regist( ATMOS_DYN_TYPE, & ! [IN]
                                       VA,             & ! [OUT]
                                       VAR_NAME,       & ! [OUT]
                                       VAR_DESC,       & ! [OUT]
                                       VAR_UNIT        ) ! [OUT]

    if ( VA > 0 ) then ! additional prognostic variables

       if ( VA > VMAX ) then
          LOG_ERROR("ATMOS_DYN_vars_setup",*) 'number of the prognostic variables is exceed the limit', VA, ' > ', VMAX
          call PRC_abort
       endif
       allocate( PROG(KA,IA,JA,VA) )
       PROG(:,:,:,:) = UNDEF

       LOG_NEWLINE
       LOG_INFO("ATMOS_DYN_vars_setup",*) '[ATMOS_DYN] prognostic/diagnostic variables'
       LOG_INFO_CONT('(1x,A,A24,A,A48,A,A12,A)') &
                  '      |', 'VARNAME                 ','|', &
                  'DESCRIPTION                                     ', '[', 'UNIT        ', ']'
       do iv = 1, VA
          LOG_INFO_CONT('(1x,A,I3,A,A24,A,A48,A,A12,A)') &
                     'NO.',iv,'|',VAR_NAME(iv),'|',VAR_DESC(iv),'[',VAR_UNIT(iv),']'
       enddo

       LOG_NEWLINE
       if ( ATMOS_DYN_RESTART_IN_BASENAME /= '' ) then
          LOG_INFO("ATMOS_DYN_vars_setup",*) 'Restart input?  : YES, file = ', trim(ATMOS_DYN_RESTART_IN_BASENAME)
          LOG_INFO("ATMOS_DYN_vars_setup",*) 'Add timelabel?  : ', ATMOS_DYN_RESTART_IN_POSTFIX_TIMELABEL
       else
          LOG_INFO("ATMOS_DYN_vars_setup",*) 'Restart input?  : NO'
       endif
       if (       ATMOS_DYN_RESTART_OUTPUT             &
            .AND. ATMOS_DYN_RESTART_OUT_BASENAME /= '' ) then
          LOG_INFO("ATMOS_DYN_vars_setup",*) 'Restart output? : YES, file = ', trim(ATMOS_DYN_RESTART_OUT_BASENAME)
          LOG_INFO("ATMOS_DYN_vars_setup",*) 'Add timelabel?  : ', ATMOS_DYN_RESTART_OUT_POSTFIX_TIMELABEL
       else
          LOG_INFO("ATMOS_DYN_vars_setup",*) 'Restart output? : NO'
          ATMOS_DYN_RESTART_OUTPUT = .false.
       endif

    else ! no additional prognostic variables

       allocate( PROG(KA,IA,JA,1) ) ! for safety
       PROG(:,:,:,:) = UNDEF
       ATMOS_DYN_RESTART_OUTPUT = .false.

    endif

    !$acc enter data create(PROG)

    return
  end subroutine ATMOS_DYN_vars_setup

  !-----------------------------------------------------------------------------
  !> Finalize
  subroutine ATMOS_DYN_vars_finalize
    implicit none
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("ATMOS_DYN_vars_finalize",*) 'Finalize'

    if ( allocated( PROG ) ) then
       !$acc exit data delete(PROG)
       deallocate( PROG )
    end if

    return
  end subroutine ATMOS_DYN_vars_finalize

  !-----------------------------------------------------------------------------
  !> HALO Communication
  subroutine ATMOS_DYN_vars_fillhalo
    use scale_comm_cartesC, only: &
       COMM_vars8, &
       COMM_wait
    implicit none

    integer :: i, j, iv
    !---------------------------------------------------------------------------

    do iv = 1, VA
       !acc kernels
       do j  = JS, JE
       do i  = IS, IE
          PROG(   1:KS-1,i,j,iv) = PROG(KS,i,j,iv)
          PROG(KE+1:KA,  i,j,iv) = PROG(KE,i,j,iv)
       enddo
       enddo
       !acc end kernels

       call COMM_vars8( PROG(:,:,:,iv), iv )
    enddo

    do iv = 1, VA
       call COMM_wait ( PROG(:,:,:,iv), iv )
    enddo

    return
  end subroutine ATMOS_DYN_vars_fillhalo

  !-----------------------------------------------------------------------------
  !> Open restart file for read
  subroutine ATMOS_DYN_vars_restart_open
    use scale_time, only: &
       TIME_gettimelabel
    use scale_file_cartesC, only: &
       FILE_CARTESC_open
    implicit none

    character(len=19)     :: timelabel
    character(len=H_LONG) :: basename
    !---------------------------------------------------------------------------

    if ( VA < 1 ) return

    LOG_NEWLINE
    LOG_INFO("ATMOS_DYN_vars_restart_open",*) 'Open restart file (ATMOS_DYN) '

    if ( ATMOS_DYN_RESTART_IN_BASENAME /= '' ) then

       if ( ATMOS_DYN_RESTART_IN_POSTFIX_TIMELABEL ) then
          call TIME_gettimelabel( timelabel )
          basename = trim(ATMOS_DYN_RESTART_IN_BASENAME)//'_'//trim(timelabel)
       else
          basename = trim(ATMOS_DYN_RESTART_IN_BASENAME)
       endif

       LOG_INFO("ATMOS_DYN_vars_restart_open",*) 'basename: ', trim(basename)

       call FILE_CARTESC_open( basename, restart_fid, aggregate=ATMOS_DYN_RESTART_IN_AGGREGATE )
    else
       LOG_INFO("ATMOS_DYN_vars_restart_open",*) 'restart file for ATMOS_DYN is not specified.'
    endif

    return
  end subroutine ATMOS_DYN_vars_restart_open

  !-----------------------------------------------------------------------------
  !> Read restart
  subroutine ATMOS_DYN_vars_restart_read
    use scale_prc, only: &
       PRC_abort
    use scale_file, only: &
       FILE_get_AGGREGATE
    use scale_file_cartesC, only: &
       FILE_CARTESC_read, &
       FILE_CARTESC_flush
    implicit none

    integer  :: i, j, iv
    !---------------------------------------------------------------------------

    if ( restart_fid /= -1 ) then
       LOG_NEWLINE
       LOG_INFO("ATMOS_DYN_vars_restart_read",*) 'Read from restart file (ATMOS_DYN) '

       do iv = 1, VA
          call FILE_CARTESC_read( restart_fid, VAR_NAME(iv), 'ZXY', & ! [IN]
                                  PROG(:,:,:,iv)                    ) ! [OUT]
       enddo

       if ( FILE_get_AGGREGATE(restart_fid) ) then
          call FILE_CARTESC_flush( restart_fid )
          !$acc update device(PROG)
          ! X/Y halos have been read from file

          ! fill K halos
          !acc kernels copy(PROG)
          do iv = 1, VA
             do j  = 1, JA
             do i  = 1, IA
                PROG(   1:KS-1,i,j,iv) = PROG(KS,i,j,iv)
                PROG(KE+1:KA,  i,j,iv) = PROG(KE,i,j,iv)
             enddo
             enddo
          enddo
          !acc end kernels
       else
          call ATMOS_DYN_vars_fillhalo
       end if

       call ATMOS_DYN_vars_check

    else
       if ( VA > 0 ) then
          LOG_ERROR("ATMOS_DYN_vars_restart_read",*) 'invalid restart file ID for ATMOS_DYN.'
          call PRC_abort
       end if
    endif

    return
  end subroutine ATMOS_DYN_vars_restart_read

  !> Create restart file
  subroutine ATMOS_DYN_vars_restart_create
    use scale_time, only: &
       TIME_gettimelabel
    use scale_file_cartesC, only: &
       FILE_CARTESC_create
    implicit none

    character(len=19)     :: timelabel
    character(len=H_LONG) :: basename
    !---------------------------------------------------------------------------

    if ( VA < 1 ) return

    if ( ATMOS_DYN_RESTART_OUT_BASENAME /= '' ) then

       LOG_NEWLINE
       LOG_INFO("ATMOS_DYN_vars_restart_create",*) 'Create restart file (ATMOS_DYN) '

       if ( ATMOS_DYN_RESTART_OUT_POSTFIX_TIMELABEL ) then
          call TIME_gettimelabel( timelabel )
          basename = trim(ATMOS_DYN_RESTART_OUT_BASENAME)//'_'//trim(timelabel)
       else
          basename = trim(ATMOS_DYN_RESTART_OUT_BASENAME)
       endif

       LOG_INFO("ATMOS_DYN_vars_restart_create",*) 'basename: ', trim(basename)

       call FILE_CARTESC_create( &
            basename, ATMOS_DYN_RESTART_OUT_TITLE, ATMOS_DYN_RESTART_OUT_DTYPE, & ! [IN]
            restart_fid,                                                        & ! [OUT]
            aggregate=ATMOS_DYN_RESTART_OUT_AGGREGATE                           ) ! [IN]

    endif

    return
  end subroutine ATMOS_DYN_vars_restart_create

  !-----------------------------------------------------------------------------
  !> Exit netCDF define mode
  subroutine ATMOS_DYN_vars_restart_enddef
    use scale_file_cartesC, only: &
       FILE_CARTESC_enddef
    implicit none

    if ( restart_fid /= -1 .and. VA > 0 ) then
       call FILE_CARTESC_enddef( restart_fid ) ! [IN]
    endif

    return
  end subroutine ATMOS_DYN_vars_restart_enddef

  !-----------------------------------------------------------------------------
  !> Close restart file
  subroutine ATMOS_DYN_vars_restart_close
    use scale_file_cartesC, only: &
       FILE_CARTESC_close
    implicit none
    !---------------------------------------------------------------------------

    if ( restart_fid /= -1 ) then
       LOG_NEWLINE
       LOG_INFO("ATMOS_DYN_vars_restart_close",*) 'Close restart file (ATMOS_DYN) '

       call FILE_CARTESC_close( restart_fid ) ! [IN]

       restart_fid = -1
    endif

    return
  end subroutine ATMOS_DYN_vars_restart_close

  !-----------------------------------------------------------------------------
  !> Define variables in restart file
  subroutine ATMOS_DYN_vars_restart_def_var
    use scale_file_cartesC, only: &
       FILE_CARTESC_def_var
    implicit none

    integer iv
    !---------------------------------------------------------------------------

    if ( restart_fid /= -1 ) then

       do iv = 1, VA
          call FILE_CARTESC_def_var( restart_fid, VAR_NAME(iv), VAR_DESC(iv), VAR_UNIT(iv), 'ZXY', ATMOS_DYN_RESTART_OUT_DTYPE, &
                                     VAR_ID(iv) )
       enddo

    endif

    return
  end subroutine ATMOS_DYN_vars_restart_def_var

  !-----------------------------------------------------------------------------
  !> Write variables to restart file
  subroutine ATMOS_DYN_vars_restart_write
    use scale_file_cartesC, only: &
       FILE_CARTESC_write_var
    implicit none

    integer  :: iv
    !---------------------------------------------------------------------------

    if ( restart_fid /= -1 ) then

       call ATMOS_DYN_vars_fillhalo

       call ATMOS_DYN_vars_check

       !$acc update host(PROG)
       do iv = 1, VA
          call FILE_CARTESC_write_var( restart_fid, VAR_ID(iv), PROG(:,:,:,iv), VAR_NAME(iv), 'ZXY' ) ! [IN]
       enddo

    endif

    return
  end subroutine ATMOS_DYN_vars_restart_write

  subroutine ATMOS_DYN_vars_check
    use scale_statistics, only: &
       STATISTICS_total
    use scale_atmos_grid_cartesC_real, only: &
       ATMOS_GRID_CARTESC_REAL_VOL, &
       ATMOS_GRID_CARTESC_REAL_TOTVOL
    implicit none
    integer :: iv

    do iv = 1, VA
       call VALCHECK( KA, KS, KE, IA, IS, IE, JA, JS, JE, &
                      PROG(:,:,:,iv),                      & ! (in)
                      -1.0E20_RP, 1.0E20_RP, VAR_NAME(iv), & ! (in)
                      __FILE__, __LINE__                   ) ! (in)
       call STATISTICS_total( KA, KS, KE, IA, IS, IE, JA, JS, JE, &
                              PROG(:,:,:,iv), VAR_NAME(iv),       & ! (in)
                              ATMOS_GRID_CARTESC_REAL_VOL(:,:,:), & ! (in)
                              ATMOS_GRID_CARTESC_REAL_TOTVOL      ) ! (in)
    enddo

    return
  end subroutine ATMOS_DYN_vars_check

end module mod_atmos_dyn_vars
