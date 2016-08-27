!-------------------------------------------------------------------------------
!> module Atmosphere / Dynamics
!!
!! @par Description
!!          Container for mod_atmos_dyn
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2014-05-04 (H.Yashiro)    [new]
!<
!-------------------------------------------------------------------------------
#include "inc_openmp.h"
module mod_atmos_dyn_vars
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_stdio
  use scale_prof
  use scale_grid_index
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
  logical,               public :: ATMOS_DYN_RESTART_OUTPUT                = .false.             !< output restart file?

  character(len=H_LONG), public :: ATMOS_DYN_RESTART_IN_BASENAME           = ''                  !< Basename of the input  file
  logical,               public :: ATMOS_DYN_RESTART_IN_POSTFIX_TIMELABEL  = .false.             !< Add timelabel to the basename of input  file?
  character(len=H_LONG), public :: ATMOS_DYN_RESTART_OUT_BASENAME          = ''                  !< Basename of the output file
  logical,               public :: ATMOS_DYN_RESTART_OUT_POSTFIX_TIMELABEL = .true.              !< Add timelabel to the basename of output file?
  character(len=H_MID),  public :: ATMOS_DYN_RESTART_OUT_TITLE             = 'ATMOS_DYN restart' !< title    of the output file
  character(len=H_MID),  public :: ATMOS_DYN_RESTART_OUT_DTYPE             = 'DEFAULT'           !< REAL4 or REAL8

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
    use scale_process, only: &
       PRC_MPIstop
    use scale_const, only: &
       UNDEF => CONST_UNDEF
    use mod_atmos_admin, only: &
       ATMOS_DYN_TYPE
    use scale_atmos_dyn_tstep_short, only: &
       ATMOS_DYN_Tstep_short_regist
    implicit none

    NAMELIST / PARAM_ATMOS_DYN_VARS / &
       ATMOS_DYN_RESTART_IN_BASENAME,           &
       ATMOS_DYN_RESTART_IN_POSTFIX_TIMELABEL,  &
       ATMOS_DYN_RESTART_OUTPUT,                &
       ATMOS_DYN_RESTART_OUT_BASENAME,          &
       ATMOS_DYN_RESTART_OUT_POSTFIX_TIMELABEL, &
       ATMOS_DYN_RESTART_OUT_TITLE,             &
       ATMOS_DYN_RESTART_OUT_DTYPE

    integer :: ierr
    integer :: iv
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '++++++ Module[VARS] / Categ[ATMOS DYN] / Origin[SCALE-RM]'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_ATMOS_DYN_VARS,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_ATMOS_DYN_VARS. Check!'
       call PRC_MPIstop
    endif
    if( IO_LNML ) write(IO_FID_LOG,nml=PARAM_ATMOS_DYN_VARS)

    call ATMOS_DYN_Tstep_short_regist( ATMOS_DYN_TYPE, & ! [IN]
                                       VA,             & ! [OUT]
                                       VAR_NAME,       & ! [OUT]
                                       VAR_DESC,       & ! [OUT]
                                       VAR_UNIT        ) ! [OUT]

    if ( VA > 0 ) then ! additional prognostic variables

       if ( VA > VMAX ) then
          write(*,*) 'xxx number of the prognostic variables is exceed the limit', VA, ' > ', VMAX
          call PRC_MPIstop
       endif
       allocate( PROG(KA,IA,JA,VA) )
       PROG(:,:,:,:) = UNDEF

       if( IO_L ) write(IO_FID_LOG,*)
       if( IO_L ) write(IO_FID_LOG,*) '*** [ATMOS_DYN] prognostic/diagnostic variables'
       if( IO_L ) write(IO_FID_LOG,'(1x,A,A15,A,A32,3(A))') &
                  '***       |','VARNAME        ','|', 'DESCRIPTION                     ','[', 'UNIT            ',']'
       do iv = 1, VA
          if( IO_L ) write(IO_FID_LOG,'(1x,A,i3,A,A15,A,A32,3(A))') &
                     '*** NO.',iv,'|',VAR_NAME(iv),'|',VAR_DESC(iv),'[',VAR_UNIT(iv),']'
       enddo

       if( IO_L ) write(IO_FID_LOG,*)
       if ( ATMOS_DYN_RESTART_IN_BASENAME /= '' ) then
          if( IO_L ) write(IO_FID_LOG,*) '*** Restart input?  : YES, file = ', trim(ATMOS_DYN_RESTART_IN_BASENAME)
          if( IO_L ) write(IO_FID_LOG,*) '*** Add timelabel?  : ', ATMOS_DYN_RESTART_IN_POSTFIX_TIMELABEL
       else
          if( IO_L ) write(IO_FID_LOG,*) '*** Restart input?  : NO'
       endif
       if (       ATMOS_DYN_RESTART_OUTPUT             &
            .AND. ATMOS_DYN_RESTART_OUT_BASENAME /= '' ) then
          if( IO_L ) write(IO_FID_LOG,*) '*** Restart output? : YES, file = ', trim(ATMOS_DYN_RESTART_OUT_BASENAME)
          if( IO_L ) write(IO_FID_LOG,*) '*** Add timelabel?  : ', ATMOS_DYN_RESTART_OUT_POSTFIX_TIMELABEL
       else
          if( IO_L ) write(IO_FID_LOG,*) '*** Restart output? : NO'
          ATMOS_DYN_RESTART_OUTPUT = .false.
       endif

    else ! no additional prognostic variables

       allocate( PROG(KA,IA,JA,1) ) ! for safety
       PROG(:,:,:,:) = UNDEF
       ATMOS_DYN_RESTART_OUTPUT = .false.

    endif

    return
  end subroutine ATMOS_DYN_vars_setup

  !-----------------------------------------------------------------------------
  !> HALO Communication
  subroutine ATMOS_DYN_vars_fillhalo
    use scale_comm, only: &
       COMM_vars8, &
       COMM_wait
    implicit none

    integer :: i, j, iv
    !---------------------------------------------------------------------------

    do iv = 1, VA
       do j  = JS, JE
       do i  = IS, IE
          PROG(   1:KS-1,i,j,iv) = PROG(KS,i,j,iv)
          PROG(KE+1:KA,  i,j,iv) = PROG(KE,i,j,iv)
       enddo
       enddo

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
    use scale_fileio, only: &
       FILEIO_open
    implicit none

    character(len=20)     :: timelabel
    character(len=H_LONG) :: basename
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** Input restart file (ATMOS_DYN) ***'

    if ( ATMOS_DYN_RESTART_IN_BASENAME /= '' ) then

       if ( ATMOS_DYN_RESTART_IN_POSTFIX_TIMELABEL ) then
          call TIME_gettimelabel( timelabel )
          basename = trim(ATMOS_DYN_RESTART_IN_BASENAME)//'_'//trim(timelabel)
       else
          basename = trim(ATMOS_DYN_RESTART_IN_BASENAME)
       endif

       if( IO_L ) write(IO_FID_LOG,*) '*** basename: ', trim(basename)

       call FILEIO_open( restart_fid, basename )
    else
       if( IO_L ) write(IO_FID_LOG,*) '*** restart file for ATMOS_DYN is not specified.'
    endif

    return
  end subroutine ATMOS_DYN_vars_restart_open

  !-----------------------------------------------------------------------------
  !> Read restart
  subroutine ATMOS_DYN_vars_restart_read
    use scale_fileio, only: &
       FILEIO_read, &
       FILEIO_flush
    use scale_rm_statistics, only: &
       STAT_total
    implicit none

    real(RP) :: total
    integer :: iv, i, j
    !---------------------------------------------------------------------------

    if ( restart_fid .NE. -1 ) then
       do iv = 1, VA
          call FILEIO_read( PROG(:,:,:,iv),                          & ! [OUT]
                            restart_fid, VAR_NAME(iv), 'ZXY', step=1 ) ! [IN]
       enddo

       if ( IO_PNETCDF ) then
          call FILEIO_flush( restart_fid )
          ! X/Y halos have been read from file

          ! fill K halos
          do iv = 1, VA
             do j  = 1, JA
             do i  = 1, IA
                PROG(   1:KS-1,i,j,iv) = PROG(KS,i,j,iv)
                PROG(KE+1:KA,  i,j,iv) = PROG(KE,i,j,iv)
             enddo
             enddo
          enddo
       else
          call ATMOS_DYN_vars_fillhalo
       end if

       do iv = 1, VA
          call STAT_total( total, PROG(:,:,:,iv), VAR_NAME(iv) )
       enddo
    else
       if ( IO_L ) write(IO_FID_LOG,*) '*** invalid restart file ID for ATMOS_DYN.'
    endif

    return
  end subroutine ATMOS_DYN_vars_restart_read

  !> Create restart file
  subroutine ATMOS_DYN_vars_restart_create
    use scale_time, only: &
       TIME_gettimelabel
    use scale_fileio, only: &
       FILEIO_create
    implicit none

    character(len=20)     :: timelabel
    character(len=H_LONG) :: basename

    integer :: iv
    !---------------------------------------------------------------------------

    if ( ATMOS_DYN_RESTART_OUT_BASENAME /= '' ) then

       if( IO_L ) write(IO_FID_LOG,*)
       if( IO_L ) write(IO_FID_LOG,*) '*** Output restart file (ATMOS_DYN) ***'

       if ( ATMOS_DYN_RESTART_OUT_POSTFIX_TIMELABEL ) then
          call TIME_gettimelabel( timelabel )
          basename = trim(ATMOS_DYN_RESTART_OUT_BASENAME)//'_'//trim(timelabel)
       else
          basename = trim(ATMOS_DYN_RESTART_OUT_BASENAME)
       endif

       if( IO_L ) write(IO_FID_LOG,*) '*** basename: ', trim(basename)

       call FILEIO_create( restart_fid,                                                       & ! [OUT]
                           basename, ATMOS_DYN_RESTART_OUT_TITLE, ATMOS_DYN_RESTART_OUT_DTYPE ) ! [IN]

    endif

    return
  end subroutine ATMOS_DYN_vars_restart_create

  !-----------------------------------------------------------------------------
  !> Exit netCDF define mode
  subroutine ATMOS_DYN_vars_restart_enddef
    use scale_fileio, only: &
       FILEIO_enddef
    implicit none

    if ( restart_fid .NE. -1 ) then
       call FILEIO_enddef( restart_fid ) ! [IN]
    endif

    return
  end subroutine ATMOS_DYN_vars_restart_enddef

  !-----------------------------------------------------------------------------
  !> Close restart file
  subroutine ATMOS_DYN_vars_restart_close
    use scale_fileio, only: &
       FILEIO_close
    implicit none

    if ( restart_fid .NE. -1 ) then
       call FILEIO_close( restart_fid ) ! [IN]
       restart_fid = -1
    endif

    return
  end subroutine ATMOS_DYN_vars_restart_close

  !-----------------------------------------------------------------------------
  !> Define variables in restart file
  subroutine ATMOS_DYN_vars_restart_def_var
    use scale_fileio, only: &
       FILEIO_def_var
    implicit none

    integer iv
    !---------------------------------------------------------------------------

    if ( restart_fid .NE. -1 ) then

       do iv = 1, VA
          call FILEIO_def_var( restart_fid, VAR_ID(iv), VAR_NAME(iv), VAR_DESC(iv), &
                               VAR_UNIT(iv), 'ZXY', ATMOS_DYN_RESTART_OUT_DTYPE     ) ! [IN]
       enddo

    endif

    return
  end subroutine ATMOS_DYN_vars_restart_def_var

  !-----------------------------------------------------------------------------
  !> Write variables to restart file
  subroutine ATMOS_DYN_vars_restart_write
    use scale_fileio, only: &
       FILEIO_write_var
    implicit none

    integer iv
    !---------------------------------------------------------------------------

    if ( restart_fid .NE. -1 ) then

       do iv = 1, VA
          call FILEIO_write_var( restart_fid, VAR_ID(iv), PROG(:,:,:,iv), VAR_NAME(iv), 'ZXY' ) ! [IN]
       enddo

    endif

    return
  end subroutine ATMOS_DYN_vars_restart_write

end module mod_atmos_dyn_vars
