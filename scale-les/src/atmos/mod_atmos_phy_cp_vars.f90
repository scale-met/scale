!-------------------------------------------------------------------------------
!> module Atmosphere / Physics Cumulus
!!
!! @par Description
!!          Container for mod_atmos_phy_cp
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2014-05-04 (H.Yashiro)    [new]
!<
!-------------------------------------------------------------------------------
#include "inc_openmp.h"
module mod_atmos_phy_cp_vars
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
  public :: ATMOS_PHY_CP_vars_setup
  public :: ATMOS_PHY_CP_vars_fillhalo
  public :: ATMOS_PHY_CP_vars_restart_read
  public :: ATMOS_PHY_CP_vars_restart_write

  !-----------------------------------------------------------------------------
  !
  !++ included parameters
  !
  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  real(RP), public, allocatable :: ATMOS_PHY_CP_MOMZ_t(:,:,:)   ! tendency MOMZ [kg/m2/s2]
  real(RP), public, allocatable :: ATMOS_PHY_CP_MOMX_t(:,:,:)   ! tendency MOMX [kg/m2/s2]
  real(RP), public, allocatable :: ATMOS_PHY_CP_MOMY_t(:,:,:)   ! tendency MOMY [kg/m2/s2]
  real(RP), public, allocatable :: ATMOS_PHY_CP_RHOT_t(:,:,:)   ! tendency RHOT [K*kg/m3/s]
  real(RP), public, allocatable :: ATMOS_PHY_CP_QTRC_t(:,:,:,:) ! tendency QTRC [kg/kg/s]

  real(RP), public, allocatable :: ATMOS_PHY_CP_MFLX_cloudbase(:,:) ! cloud base mass flux [kg/m2/s]

  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  logical,                private, save      :: ATMOS_PHY_CP_RESTART_OUTPUT        = .false.
  character(len=H_LONG),  private, save      :: ATMOS_PHY_CP_RESTART_IN_BASENAME   = ''
  character(len=H_LONG),  private, save      :: ATMOS_PHY_CP_RESTART_OUT_BASENAME  = ''
  character(len=H_MID),   private, save      :: ATMOS_PHY_CP_RESTART_OUT_TITLE     = 'ATMOS_PHY_CP restart'
  character(len=H_MID),   private, save      :: ATMOS_PHY_CP_RESTART_OUT_DTYPE     = 'DEFAULT'              !< REAL4 or REAL8

  integer,                private, parameter :: VMAX = 1       !< number of the variables
  character(len=H_SHORT), private, save      :: VAR_NAME(VMAX) !< name  of the variables
  character(len=H_MID),   private, save      :: VAR_DESC(VMAX) !< desc. of the variables
  character(len=H_SHORT), private, save      :: VAR_UNIT(VMAX) !< unit  of the variables

  data VAR_NAME / 'MFLX_cloudbase' /
  data VAR_DESC / 'cloud base mass flux' /
  data VAR_UNIT / 'kg/m2/s' /

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine ATMOS_PHY_CP_vars_setup
    use scale_process, only: &
       PRC_MPIstop
    implicit none

    NAMELIST / PARAM_ATMOS_PHY_CP_VARS / &
       ATMOS_PHY_CP_RESTART_IN_BASENAME, &
       ATMOS_PHY_CP_RESTART_OUTPUT,      &
       ATMOS_PHY_CP_RESTART_OUT_BASENAME

    integer :: v, ierr
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[ATMOS_PHY_CP VARS]/Categ[ATMOS]'
    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_ATMOS_PHY_CP_VARS,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_ATMOS_PHY_CP_VARS. Check!'
       call PRC_MPIstop
    endif
    if( IO_L ) write(IO_FID_LOG,nml=PARAM_ATMOS_PHY_CP_VARS)

    allocate( ATMOS_PHY_CP_MOMZ_t(KA,IA,JA)    )
    allocate( ATMOS_PHY_CP_MOMX_t(KA,IA,JA)    )
    allocate( ATMOS_PHY_CP_MOMY_t(KA,IA,JA)    )
    allocate( ATMOS_PHY_CP_RHOT_t(KA,IA,JA)    )
    allocate( ATMOS_PHY_CP_QTRC_t(KA,IA,JA,QA) )

    allocate( ATMOS_PHY_CP_MFLX_cloudbase(IA,JA) )

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** [ATMOS_PHY_CP] prognostic/diagnostic variables'
    if( IO_L ) write(IO_FID_LOG,'(1x,A,A8,A,A32,3(A))') &
               '***       |',' VARNAME','|', 'DESCRIPTION                     ','[', 'UNIT            ',']'
    do v = 1, VMAX
       if( IO_L ) write(IO_FID_LOG,'(1x,A,i3,A,A8,A,A32,3(A))') &
                  '*** NO.',v,'|',trim(VAR_NAME(v)),'|',VAR_DESC(v),'[',VAR_UNIT(v),']'
    enddo

    ! restart switch
    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) 'Output...'
    if ( ATMOS_PHY_CP_RESTART_OUTPUT ) then
       if( IO_L ) write(IO_FID_LOG,*) '  Restart output (ATMOS_PHY_CP) : YES'
    else
       if( IO_L ) write(IO_FID_LOG,*) '  Restart output (ATMOS_PHY_CP) : NO'
    endif
    if( IO_L ) write(IO_FID_LOG,*)

    return
  end subroutine ATMOS_PHY_CP_vars_setup

  !-----------------------------------------------------------------------------
  !> HALO Communication
  subroutine ATMOS_PHY_CP_vars_fillhalo
    use scale_comm, only: &
       COMM_vars8, &
       COMM_wait
    implicit none
    !---------------------------------------------------------------------------

    call COMM_vars8( ATMOS_PHY_CP_MFLX_cloudbase(:,:), 1 )
    call COMM_wait ( ATMOS_PHY_CP_MFLX_cloudbase(:,:), 1 )

    return
  end subroutine ATMOS_PHY_CP_vars_fillhalo

  !-----------------------------------------------------------------------------
  !> Read restart
  subroutine ATMOS_PHY_CP_vars_restart_read
    use scale_fileio, only: &
       FILEIO_read
    use scale_const, only: &
       CONST_UNDEF
    implicit none
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** Input restart file (ATMOS_PHY_CP) ***'

    if ( ATMOS_PHY_CP_RESTART_IN_BASENAME /= '' ) then

       call FILEIO_read( ATMOS_PHY_CP_MFLX_cloudbase(:,:),                           & ! [OUT]
                         ATMOS_PHY_CP_RESTART_IN_BASENAME, VAR_NAME(1), 'XY', step=1 ) ! [IN]

       call ATMOS_PHY_CP_vars_fillhalo

    else

       if( IO_L ) write(IO_FID_LOG,*) '*** restart file for ATMOS_PHY_CP is not specified.'

       ATMOS_PHY_CP_MFLX_cloudbase(:,:) = CONST_UNDEF

    endif

    return
  end subroutine ATMOS_PHY_CP_vars_restart_read

  !-----------------------------------------------------------------------------
  !> Write restart
  subroutine ATMOS_PHY_CP_vars_restart_write
    use scale_time, only: &
       TIME_gettimelabel
    use scale_fileio, only: &
       FILEIO_write
    implicit none

    character(len=15)     :: timelabel
    character(len=H_LONG) :: basename
    !---------------------------------------------------------------------------

    if ( ATMOS_PHY_CP_RESTART_OUT_BASENAME /= '' ) then

       call TIME_gettimelabel( timelabel )
       write(basename,'(A,A,A)') trim(ATMOS_PHY_CP_RESTART_OUT_BASENAME), '_', trim(timelabel)

       if( IO_L ) write(IO_FID_LOG,*)
       if( IO_L ) write(IO_FID_LOG,*) '*** Output restart file (ATMOS_PHY_CP) ***'
       if( IO_L ) write(IO_FID_LOG,*) '*** basename: ', trim(basename)

       call FILEIO_write( ATMOS_PHY_CP_MFLX_cloudbase(:,:), basename,  ATMOS_PHY_CP_RESTART_OUT_TITLE, & ! [IN]
                          VAR_NAME(1), VAR_DESC(1), VAR_UNIT(1), 'XY', ATMOS_PHY_CP_RESTART_OUT_DTYPE  ) ! [IN]

    endif

    return
  end subroutine ATMOS_PHY_CP_vars_restart_write

end module mod_atmos_phy_cp_vars
