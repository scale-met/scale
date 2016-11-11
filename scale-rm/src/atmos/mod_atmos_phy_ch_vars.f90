!-------------------------------------------------------------------------------
!> module Atmosphere / Physics Chemistry
!!
!! @par Description
!!          Container for mod_atmos_phy_ch
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2014-05-04 (H.Yashiro)    [new]
!<
!-------------------------------------------------------------------------------
#include "inc_openmp.h"
module mod_atmos_phy_ch_vars
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_stdio
  use scale_prof
  use scale_grid_index
  use scale_tracer
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: ATMOS_PHY_CH_vars_setup
  public :: ATMOS_PHY_CH_vars_fillhalo
  public :: ATMOS_PHY_CH_vars_restart_read
  public :: ATMOS_PHY_CH_vars_restart_write

  public :: ATMOS_PHY_CH_vars_restart_create
  public :: ATMOS_PHY_CH_vars_restart_open
  public :: ATMOS_PHY_CH_vars_restart_def_var
  public :: ATMOS_PHY_CH_vars_restart_enddef
  public :: ATMOS_PHY_CH_vars_restart_close

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  logical,               public :: ATMOS_PHY_CH_RESTART_OUTPUT                = .false.                !< output restart file?

  character(len=H_LONG), public :: ATMOS_PHY_CH_RESTART_IN_BASENAME           = ''                     !< Basename of the input  file
  logical,               public :: ATMOS_PHY_CH_RESTART_IN_POSTFIX_TIMELABEL  = .false.                !< Add timelabel to the basename of input  file?
  character(len=H_LONG), public :: ATMOS_PHY_CH_RESTART_OUT_BASENAME          = ''                     !< Basename of the output file
  logical,               public :: ATMOS_PHY_CH_RESTART_OUT_POSTFIX_TIMELABEL = .true.                 !< Add timelabel to the basename of output file?
  character(len=H_MID),  public :: ATMOS_PHY_CH_RESTART_OUT_TITLE             = 'ATMOS_PHY_CH restart' !< title    of the output file
  character(len=H_MID),  public :: ATMOS_PHY_CH_RESTART_OUT_DTYPE             = 'DEFAULT'              !< REAL4 or REAL8

  real(RP), public, allocatable :: ATMOS_PHY_CH_RHOQ_t(:,:,:,:) ! tendency QTRC [kg/kg/s]

  real(RP), public, allocatable :: ATMOS_PHY_CH_O3(:,:,:) ! ozone [PPM]

  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  integer,                private, parameter :: VMAX = 1       !< number of the variables
  integer,                private, parameter :: I_O3 = 1

  character(len=H_SHORT), private            :: VAR_NAME(VMAX) !< name  of the variables
  character(len=H_MID),   private            :: VAR_DESC(VMAX) !< desc. of the variables
  character(len=H_SHORT), private            :: VAR_UNIT(VMAX) !< unit  of the variables
  integer,                private            :: VAR_ID(VMAX)   !< ID    of the variables
  integer,                private            :: restart_fid = -1  ! file ID

  data VAR_NAME / 'O3' /
  data VAR_DESC / 'Ozone' /
  data VAR_UNIT / 'PPM' /

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine ATMOS_PHY_CH_vars_setup
    use scale_process, only: &
       PRC_MPIstop
    use scale_const, only: &
       UNDEF => CONST_UNDEF
    implicit none

    NAMELIST / PARAM_ATMOS_PHY_CH_VARS / &
       ATMOS_PHY_CH_RESTART_IN_BASENAME,           &
       ATMOS_PHY_CH_RESTART_IN_POSTFIX_TIMELABEL,  &
       ATMOS_PHY_CH_RESTART_OUTPUT,                &
       ATMOS_PHY_CH_RESTART_OUT_BASENAME,          &
       ATMOS_PHY_CH_RESTART_OUT_POSTFIX_TIMELABEL, &
       ATMOS_PHY_CH_RESTART_OUT_TITLE,             &
       ATMOS_PHY_CH_RESTART_OUT_DTYPE

    integer :: ierr
    integer :: iv
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '++++++ Module[VARS] / Categ[ATMOS PHY_CH] / Origin[SCALE-RM]'

    allocate( ATMOS_PHY_CH_RHOQ_t(KA,IA,JA,QA) )
    ATMOS_PHY_CH_RHOQ_t(:,:,:,:) = UNDEF

    allocate( ATMOS_PHY_CH_O3(KA,IA,JA) )
    ATMOS_PHY_CH_O3(:,:,:) = UNDEF

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_ATMOS_PHY_CH_VARS,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_ATMOS_PHY_CH_VARS. Check!'
       call PRC_MPIstop
    endif
    if( IO_LNML ) write(IO_FID_LOG,nml=PARAM_ATMOS_PHY_CH_VARS)

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** [ATMOS_PHY_CH] prognostic/diagnostic variables'
    if( IO_L ) write(IO_FID_LOG,'(1x,A,A15,A,A32,3(A))') &
               '***       |','VARNAME        ','|', 'DESCRIPTION                     ','[', 'UNIT            ',']'
    do iv = 1, VMAX
       if( IO_L ) write(IO_FID_LOG,'(1x,A,i3,A,A15,A,A32,3(A))') &
                  '*** NO.',iv,'|',VAR_NAME(iv),'|',VAR_DESC(iv),'[',VAR_UNIT(iv),']'
    enddo

    if( IO_L ) write(IO_FID_LOG,*)
    if ( ATMOS_PHY_CH_RESTART_IN_BASENAME /= '' ) then
       if( IO_L ) write(IO_FID_LOG,*) '*** Restart input?  : YES, file = ', trim(ATMOS_PHY_CH_RESTART_IN_BASENAME)
       if( IO_L ) write(IO_FID_LOG,*) '*** Add timelabel?  : ', ATMOS_PHY_CH_RESTART_IN_POSTFIX_TIMELABEL
    else
       if( IO_L ) write(IO_FID_LOG,*) '*** Restart input?  : NO'
    endif
    if (       ATMOS_PHY_CH_RESTART_OUTPUT             &
         .AND. ATMOS_PHY_CH_RESTART_OUT_BASENAME /= '' ) then
       if( IO_L ) write(IO_FID_LOG,*) '*** Restart output? : YES, file = ', trim(ATMOS_PHY_CH_RESTART_OUT_BASENAME)
       if( IO_L ) write(IO_FID_LOG,*) '*** Add timelabel?  : ', ATMOS_PHY_CH_RESTART_OUT_POSTFIX_TIMELABEL
    else
       if( IO_L ) write(IO_FID_LOG,*) '*** Restart output? : NO'
       ATMOS_PHY_CH_RESTART_OUTPUT = .false.
    endif

    return
  end subroutine ATMOS_PHY_CH_vars_setup

  !-----------------------------------------------------------------------------
  !> HALO Communication
  subroutine ATMOS_PHY_CH_vars_fillhalo
    use scale_comm, only: &
       COMM_vars8, &
       COMM_wait
    implicit none

    integer :: i, j
    !---------------------------------------------------------------------------

    do j  = JS, JE
    do i  = IS, IE
       ATMOS_PHY_CH_O3(   1:KS-1,i,j) = ATMOS_PHY_CH_O3(KS,i,j)
       ATMOS_PHY_CH_O3(KE+1:KA,  i,j) = ATMOS_PHY_CH_O3(KE,i,j)
    enddo
    enddo

    call COMM_vars8( ATMOS_PHY_CH_O3(:,:,:), 1 )
    call COMM_wait ( ATMOS_PHY_CH_O3(:,:,:), 1 )

    return
  end subroutine ATMOS_PHY_CH_vars_fillhalo

  !-----------------------------------------------------------------------------
  !> Open restart file for read
  subroutine ATMOS_PHY_CH_vars_restart_open
    use scale_time, only: &
       TIME_gettimelabel
    use scale_fileio, only: &
       FILEIO_open
    implicit none

    character(len=20)     :: timelabel
    character(len=H_LONG) :: basename
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** Input restart file (ATMOS_PHY_CH) ***'

    if ( ATMOS_PHY_CH_RESTART_IN_BASENAME /= '' ) then

       if ( ATMOS_PHY_CH_RESTART_IN_POSTFIX_TIMELABEL ) then
          call TIME_gettimelabel( timelabel )
          basename = trim(ATMOS_PHY_CH_RESTART_IN_BASENAME)//'_'//trim(timelabel)
       else
          basename = trim(ATMOS_PHY_CH_RESTART_IN_BASENAME)
       endif

       if( IO_L ) write(IO_FID_LOG,*) '*** basename: ', trim(basename)

       call FILEIO_open( restart_fid, basename )
    else
       if( IO_L ) write(IO_FID_LOG,*) '*** restart file for ATMOS_PHY_CH is not specified.'
    endif

    return
  end subroutine ATMOS_PHY_CH_vars_restart_open

  !-----------------------------------------------------------------------------
  !> Read restart
  subroutine ATMOS_PHY_CH_vars_restart_read
    use scale_fileio, only: &
       FILEIO_read, &
       FILEIO_flush
    use scale_rm_statistics, only: &
       STAT_total
    implicit none

    real(RP) :: total
    integer  :: i, j
    !---------------------------------------------------------------------------

    if ( restart_fid .NE. -1 ) then
       if( IO_L ) write(IO_FID_LOG,*) '*** basename: ', trim(ATMOS_PHY_CH_RESTART_IN_BASENAME)

       call FILEIO_read( ATMOS_PHY_CH_O3(:,:,:),                 & ! [OUT]
                         restart_fid, VAR_NAME(1), 'ZXY', step=1 ) ! [IN]

       if ( IO_AGGREGATE ) then
          call FILEIO_flush( restart_fid )
          ! X/Y halos have been read from file

          ! fill K halos
          do j  = 1, JA
          do i  = 1, IA
             ATMOS_PHY_CH_O3(   1:KS-1,i,j) = ATMOS_PHY_CH_O3(KS,i,j)
             ATMOS_PHY_CH_O3(KE+1:KA,  i,j) = ATMOS_PHY_CH_O3(KE,i,j)
          enddo
          enddo
       else
          call ATMOS_PHY_CH_vars_fillhalo
       end if

       call STAT_total( total, ATMOS_PHY_CH_O3(:,:,:), VAR_NAME(1) )
    else
       if( IO_L ) write(IO_FID_LOG,*) '*** invalid restart file for ATMOS_PHY_CH.'
    endif

    return
  end subroutine ATMOS_PHY_CH_vars_restart_read

  !-----------------------------------------------------------------------------
  !> Create restart file
  subroutine ATMOS_PHY_CH_vars_restart_create
    use scale_time, only: &
       TIME_gettimelabel
    use scale_fileio, only: &
       FILEIO_create
    implicit none

    character(len=20)     :: timelabel
    character(len=H_LONG) :: basename
    !---------------------------------------------------------------------------

    if ( ATMOS_PHY_CH_RESTART_OUT_BASENAME /= '' ) then

       if( IO_L ) write(IO_FID_LOG,*)
       if( IO_L ) write(IO_FID_LOG,*) '*** Output restart file (ATMOS_PHY_AE) ***'

       if ( ATMOS_PHY_CH_RESTART_OUT_POSTFIX_TIMELABEL ) then
          call TIME_gettimelabel( timelabel )
          basename = trim(ATMOS_PHY_CH_RESTART_OUT_BASENAME)//'_'//trim(timelabel)
       else
          basename = trim(ATMOS_PHY_CH_RESTART_OUT_BASENAME)
       endif

       if( IO_L ) write(IO_FID_LOG,*) '*** basename: ', trim(basename)

       call FILEIO_create( restart_fid,                                                             & ! [OUT]
                           basename, ATMOS_PHY_CH_RESTART_OUT_TITLE, ATMOS_PHY_CH_RESTART_OUT_DTYPE ) ! [IN]

    endif

    return
  end subroutine ATMOS_PHY_CH_vars_restart_create

  !-----------------------------------------------------------------------------
  !> Exit netCDF define mode
  subroutine ATMOS_PHY_CH_vars_restart_enddef
    use scale_fileio, only: &
       FILEIO_enddef
    implicit none

    if ( restart_fid .NE. -1 ) then
       call FILEIO_enddef( restart_fid ) ! [IN]
    endif

    return
  end subroutine ATMOS_PHY_CH_vars_restart_enddef

  !-----------------------------------------------------------------------------
  !> Close restart file
  subroutine ATMOS_PHY_CH_vars_restart_close
    use scale_fileio, only: &
       FILEIO_close
    implicit none

    if ( restart_fid .NE. -1 ) then
       call FILEIO_close( restart_fid ) ! [IN]
       restart_fid = -1
    endif

    return
  end subroutine ATMOS_PHY_CH_vars_restart_close

  !-----------------------------------------------------------------------------
  !> Write restart
  subroutine ATMOS_PHY_CH_vars_restart_def_var
    use scale_fileio, only: &
       FILEIO_def_var
    implicit none

    !---------------------------------------------------------------------------

    if ( restart_fid .NE. -1 ) then
       call FILEIO_def_var( restart_fid, VAR_ID(1), VAR_NAME(1), VAR_DESC(1), &
                            VAR_UNIT(1), 'ZXY', ATMOS_PHY_CH_RESTART_OUT_DTYPE  ) ! [IN]
    endif

    return
  end subroutine ATMOS_PHY_CH_vars_restart_def_var

  !-----------------------------------------------------------------------------
  !> Write restart
  subroutine ATMOS_PHY_CH_vars_restart_write
    use scale_fileio, only: &
       FILEIO_write_var
    implicit none

    !---------------------------------------------------------------------------

    if ( restart_fid .NE. -1 ) then
       call FILEIO_write_var( restart_fid, VAR_ID(1), ATMOS_PHY_CH_O3(:,:,:), &
                          VAR_NAME(1), 'ZXY' ) ! [IN]
    endif

    return
  end subroutine ATMOS_PHY_CH_vars_restart_write

end module mod_atmos_phy_ch_vars
