!-------------------------------------------------------------------------------
!> module Atmosphere / Physics Radiation
!!
!! @par Description
!!          Container for mod_atmos_phy_rd
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2014-05-04 (H.Yashiro)    [new]
!<
!-------------------------------------------------------------------------------
#include "inc_openmp.h"
module mod_atmos_phy_rd_vars
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
  public :: ATMOS_PHY_RD_vars_setup
  public :: ATMOS_PHY_RD_vars_fillhalo
  public :: ATMOS_PHY_RD_vars_restart_read
  public :: ATMOS_PHY_RD_vars_restart_write

  !-----------------------------------------------------------------------------
  !
  !++ included parameters
  !
  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  real(RP), public, allocatable :: ATMOS_PHY_RD_RHOT_t(:,:,:)   ! tendency RHOT [K*kg/m3/s]

  real(RP), public, allocatable :: ATMOS_PHY_RD_SFLX_LW_up  (:,:) ! surface upward   longwave  flux [J/m2/s]
  real(RP), public, allocatable :: ATMOS_PHY_RD_SFLX_LW_dn  (:,:) ! surface downward longwave  flux [J/m2/s]
  real(RP), public, allocatable :: ATMOS_PHY_RD_SFLX_SW_up  (:,:) ! surface upward   shortwave flux [J/m2/s]
  real(RP), public, allocatable :: ATMOS_PHY_RD_SFLX_SW_dn  (:,:) ! surface downward shortwave flux [J/m2/s]

  real(RP), public, allocatable :: ATMOS_PHY_RD_TOAFLX_LW_up(:,:) ! TOA upward   longwave  flux [J/m2/s]
  real(RP), public, allocatable :: ATMOS_PHY_RD_TOAFLX_LW_dn(:,:) ! TOA downward longwave  flux [J/m2/s]
  real(RP), public, allocatable :: ATMOS_PHY_RD_TOAFLX_SW_up(:,:) ! TOA upward   shortwave flux [J/m2/s]
  real(RP), public, allocatable :: ATMOS_PHY_RD_TOAFLX_SW_dn(:,:) ! TOA downward shortwave flux [J/m2/s]

  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  logical,                private, save      :: ATMOS_PHY_RD_RESTART_OUTPUT        = .false.
  character(len=H_LONG),  private, save      :: ATMOS_PHY_RD_RESTART_IN_BASENAME   = ''
  character(len=H_LONG),  private, save      :: ATMOS_PHY_RD_RESTART_OUT_BASENAME  = ''
  character(len=H_MID),   private, save      :: ATMOS_PHY_RD_RESTART_OUT_TITLE     = 'ATMOS_PHY_RD restart'
  character(len=H_MID),   private, save      :: ATMOS_PHY_RD_RESTART_OUT_DTYPE     = 'DEFAULT'              !< REAL4 or REAL8

  integer,                private, parameter :: VMAX = 8       !< number of the variables
  character(len=H_SHORT), private, save      :: VAR_NAME(VMAX) !< name  of the variables
  character(len=H_MID),   private, save      :: VAR_DESC(VMAX) !< desc. of the variables
  character(len=H_SHORT), private, save      :: VAR_UNIT(VMAX) !< unit  of the variables

  data VAR_NAME / 'SFLX_LW_up',   &
                  'SFLX_LW_dn',   &
                  'SFLX_SW_up',   &
                  'SFLX_SW_dn',   &
                  'TOAFLX_LW_up', &
                  'TOAFLX_LW_dn', &
                  'TOAFLX_SW_up', &
                  'TOAFLX_SW_dn'  /
  data VAR_DESC / 'surface upward   longwave  flux', &
                  'surface downward longwave  flux', &
                  'surface upward   shortwave flux', &
                  'surface downward shortwave flux', &
                  'TOA upward   longwave  flux', &
                  'TOA downward longwave  flux', &
                  'TOA upward   shortwave flux', &
                  'TOA downward shortwave flux' /
  data VAR_UNIT / 'W/m2', &
                  'W/m2', &
                  'W/m2', &
                  'W/m2', &
                  'W/m2', &
                  'W/m2', &
                  'W/m2', &
                  'W/m2'  /

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine ATMOS_PHY_RD_vars_setup
    use scale_process, only: &
       PRC_MPIstop
    implicit none

    NAMELIST / PARAM_ATMOS_PHY_RD_VARS / &
       ATMOS_PHY_RD_RESTART_IN_BASENAME, &
       ATMOS_PHY_RD_RESTART_OUTPUT,      &
       ATMOS_PHY_RD_RESTART_OUT_BASENAME

    integer :: v, ierr
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[ATMOS_PHY_RD VARS]/Categ[ATMOS]'
    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_ATMOS_PHY_RD_VARS,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_ATMOS_PHY_RD_VARS. Check!'
       call PRC_MPIstop
    endif
    if( IO_L ) write(IO_FID_LOG,nml=PARAM_ATMOS_PHY_RD_VARS)

    allocate( ATMOS_PHY_RD_RHOT_t(KA,IA,JA) )

    allocate( ATMOS_PHY_RD_SFLX_LW_up  (IA,JA) )
    allocate( ATMOS_PHY_RD_SFLX_LW_dn  (IA,JA) )
    allocate( ATMOS_PHY_RD_SFLX_SW_up  (IA,JA) )
    allocate( ATMOS_PHY_RD_SFLX_SW_dn  (IA,JA) )
    allocate( ATMOS_PHY_RD_TOAFLX_LW_up(IA,JA) )
    allocate( ATMOS_PHY_RD_TOAFLX_LW_dn(IA,JA) )
    allocate( ATMOS_PHY_RD_TOAFLX_SW_up(IA,JA) )
    allocate( ATMOS_PHY_RD_TOAFLX_SW_dn(IA,JA) )

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** [ATMOS_PHY_RD] prognostic/diagnostic variables'
    if( IO_L ) write(IO_FID_LOG,'(1x,A,A8,A,A32,3(A))') &
               '***       |',' VARNAME','|', 'DESCRIPTION                     ','[', 'UNIT            ',']'
    do v = 1, VMAX
       if( IO_L ) write(IO_FID_LOG,'(1x,A,i3,A,A8,A,A32,3(A))') &
                  '*** NO.',v,'|',trim(VAR_NAME(v)),'|',VAR_DESC(v),'[',VAR_UNIT(v),']'
    enddo

    ! restart switch
    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) 'Output...'
    if ( ATMOS_PHY_RD_RESTART_OUTPUT ) then
       if( IO_L ) write(IO_FID_LOG,*) '  Restart output (ATMOS_PHY_RD) : YES'
    else
       if( IO_L ) write(IO_FID_LOG,*) '  Restart output (ATMOS_PHY_RD) : NO'
    endif
    if( IO_L ) write(IO_FID_LOG,*)

    return
  end subroutine ATMOS_PHY_RD_vars_setup

  !-----------------------------------------------------------------------------
  !> HALO Communication
  subroutine ATMOS_PHY_RD_vars_fillhalo
    use scale_comm, only: &
       COMM_vars8, &
       COMM_wait
    implicit none
    !---------------------------------------------------------------------------

    call COMM_vars8( ATMOS_PHY_RD_SFLX_LW_up  (:,:), 1 )
    call COMM_vars8( ATMOS_PHY_RD_SFLX_LW_dn  (:,:), 2 )
    call COMM_vars8( ATMOS_PHY_RD_SFLX_SW_up  (:,:), 3 )
    call COMM_vars8( ATMOS_PHY_RD_SFLX_SW_dn  (:,:), 4 )
    call COMM_vars8( ATMOS_PHY_RD_TOAFLX_LW_up(:,:), 5 )
    call COMM_vars8( ATMOS_PHY_RD_TOAFLX_LW_dn(:,:), 6 )
    call COMM_vars8( ATMOS_PHY_RD_TOAFLX_SW_up(:,:), 7 )
    call COMM_vars8( ATMOS_PHY_RD_TOAFLX_SW_dn(:,:), 8 )
    call COMM_wait ( ATMOS_PHY_RD_SFLX_LW_up  (:,:), 1 )
    call COMM_wait ( ATMOS_PHY_RD_SFLX_LW_dn  (:,:), 2 )
    call COMM_wait ( ATMOS_PHY_RD_SFLX_SW_up  (:,:), 3 )
    call COMM_wait ( ATMOS_PHY_RD_SFLX_SW_dn  (:,:), 4 )
    call COMM_wait ( ATMOS_PHY_RD_TOAFLX_LW_up(:,:), 5 )
    call COMM_wait ( ATMOS_PHY_RD_TOAFLX_LW_dn(:,:), 6 )
    call COMM_wait ( ATMOS_PHY_RD_TOAFLX_SW_up(:,:), 7 )
    call COMM_wait ( ATMOS_PHY_RD_TOAFLX_SW_dn(:,:), 8 )

    return
  end subroutine ATMOS_PHY_RD_vars_fillhalo

  !-----------------------------------------------------------------------------
  !> Read restart
  subroutine ATMOS_PHY_RD_vars_restart_read
    use scale_fileio, only: &
       FILEIO_read
    use scale_const, only: &
       CONST_UNDEF
    implicit none
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** Input restart file (ATMOS_PHY_RD) ***'

    if ( ATMOS_PHY_RD_RESTART_IN_BASENAME /= '' ) then

       call FILEIO_read( ATMOS_PHY_RD_SFLX_LW_up(:,:),                               & ! [OUT]
                         ATMOS_PHY_RD_RESTART_IN_BASENAME, VAR_NAME(1), 'XY', step=1 ) ! [IN]
       call FILEIO_read( ATMOS_PHY_RD_SFLX_LW_dn(:,:),                               & ! [OUT]
                         ATMOS_PHY_RD_RESTART_IN_BASENAME, VAR_NAME(2), 'XY', step=1 ) ! [IN]
       call FILEIO_read( ATMOS_PHY_RD_SFLX_SW_up(:,:),                               & ! [OUT]
                         ATMOS_PHY_RD_RESTART_IN_BASENAME, VAR_NAME(3), 'XY', step=1 ) ! [IN]
       call FILEIO_read( ATMOS_PHY_RD_SFLX_SW_dn(:,:),                               & ! [OUT]
                         ATMOS_PHY_RD_RESTART_IN_BASENAME, VAR_NAME(4), 'XY', step=1 ) ! [IN]
       call FILEIO_read( ATMOS_PHY_RD_TOAFLX_LW_up(:,:),                             & ! [OUT]
                         ATMOS_PHY_RD_RESTART_IN_BASENAME, VAR_NAME(5), 'XY', step=1 ) ! [IN]
       call FILEIO_read( ATMOS_PHY_RD_TOAFLX_LW_dn(:,:),                             & ! [OUT]
                         ATMOS_PHY_RD_RESTART_IN_BASENAME, VAR_NAME(6), 'XY', step=1 ) ! [IN]
       call FILEIO_read( ATMOS_PHY_RD_TOAFLX_SW_up(:,:),                             & ! [OUT]
                         ATMOS_PHY_RD_RESTART_IN_BASENAME, VAR_NAME(7), 'XY', step=1 ) ! [IN]
       call FILEIO_read( ATMOS_PHY_RD_TOAFLX_SW_dn(:,:),                             & ! [OUT]
                         ATMOS_PHY_RD_RESTART_IN_BASENAME, VAR_NAME(8), 'XY', step=1 ) ! [IN]

       call ATMOS_PHY_RD_vars_fillhalo

    else

       if( IO_L ) write(IO_FID_LOG,*) '*** restart file for ATMOS_PHY_RD is not specified.'

       ATMOS_PHY_RD_SFLX_LW_up  (:,:) = CONST_UNDEF
       ATMOS_PHY_RD_SFLX_LW_dn  (:,:) = CONST_UNDEF
       ATMOS_PHY_RD_SFLX_SW_up  (:,:) = CONST_UNDEF
       ATMOS_PHY_RD_SFLX_SW_dn  (:,:) = CONST_UNDEF
       ATMOS_PHY_RD_TOAFLX_LW_up(:,:) = CONST_UNDEF
       ATMOS_PHY_RD_TOAFLX_LW_dn(:,:) = CONST_UNDEF
       ATMOS_PHY_RD_TOAFLX_SW_up(:,:) = CONST_UNDEF
       ATMOS_PHY_RD_TOAFLX_SW_dn(:,:) = CONST_UNDEF

    endif

    return
  end subroutine ATMOS_PHY_RD_vars_restart_read

  !-----------------------------------------------------------------------------
  !> Write restart
  subroutine ATMOS_PHY_RD_vars_restart_write
    use scale_time, only: &
       TIME_gettimelabel
    use scale_fileio, only: &
       FILEIO_write
    implicit none

    character(len=15)     :: timelabel
    character(len=H_LONG) :: basename
    !---------------------------------------------------------------------------

    if ( ATMOS_PHY_RD_RESTART_OUT_BASENAME /= '' ) then

       call TIME_gettimelabel( timelabel )
       write(basename,'(A,A,A)') trim(ATMOS_PHY_RD_RESTART_OUT_BASENAME), '_', trim(timelabel)

       if( IO_L ) write(IO_FID_LOG,*)
       if( IO_L ) write(IO_FID_LOG,*) '*** Output restart file (ATMOS_PHY_RD) ***'
       if( IO_L ) write(IO_FID_LOG,*) '*** basename: ', trim(basename)

       call FILEIO_write( ATMOS_PHY_RD_SFLX_LW_up(:,:), basename,      ATMOS_PHY_RD_RESTART_OUT_TITLE, & ! [IN]
                          VAR_NAME(1), VAR_DESC(1), VAR_UNIT(1), 'XY', ATMOS_PHY_RD_RESTART_OUT_DTYPE  ) ! [IN]
       call FILEIO_write( ATMOS_PHY_RD_SFLX_LW_dn(:,:), basename,      ATMOS_PHY_RD_RESTART_OUT_TITLE, & ! [IN]
                          VAR_NAME(2), VAR_DESC(2), VAR_UNIT(2), 'XY', ATMOS_PHY_RD_RESTART_OUT_DTYPE  ) ! [IN]
       call FILEIO_write( ATMOS_PHY_RD_SFLX_SW_up(:,:), basename,      ATMOS_PHY_RD_RESTART_OUT_TITLE, & ! [IN]
                          VAR_NAME(3), VAR_DESC(3), VAR_UNIT(3), 'XY', ATMOS_PHY_RD_RESTART_OUT_DTYPE  ) ! [IN]
       call FILEIO_write( ATMOS_PHY_RD_SFLX_SW_dn(:,:), basename,      ATMOS_PHY_RD_RESTART_OUT_TITLE, & ! [IN]
                          VAR_NAME(4), VAR_DESC(4), VAR_UNIT(4), 'XY', ATMOS_PHY_RD_RESTART_OUT_DTYPE  ) ! [IN]
       call FILEIO_write( ATMOS_PHY_RD_TOAFLX_LW_up(:,:), basename,    ATMOS_PHY_RD_RESTART_OUT_TITLE, & ! [IN]
                          VAR_NAME(5), VAR_DESC(5), VAR_UNIT(5), 'XY', ATMOS_PHY_RD_RESTART_OUT_DTYPE  ) ! [IN]
       call FILEIO_write( ATMOS_PHY_RD_TOAFLX_LW_dn(:,:), basename,    ATMOS_PHY_RD_RESTART_OUT_TITLE, & ! [IN]
                          VAR_NAME(6), VAR_DESC(6), VAR_UNIT(6), 'XY', ATMOS_PHY_RD_RESTART_OUT_DTYPE  ) ! [IN]
       call FILEIO_write( ATMOS_PHY_RD_TOAFLX_SW_up(:,:), basename,    ATMOS_PHY_RD_RESTART_OUT_TITLE, & ! [IN]
                          VAR_NAME(7), VAR_DESC(7), VAR_UNIT(7), 'XY', ATMOS_PHY_RD_RESTART_OUT_DTYPE  ) ! [IN]
       call FILEIO_write( ATMOS_PHY_RD_TOAFLX_SW_dn(:,:), basename,    ATMOS_PHY_RD_RESTART_OUT_TITLE, & ! [IN]
                          VAR_NAME(8), VAR_DESC(8), VAR_UNIT(8), 'XY', ATMOS_PHY_RD_RESTART_OUT_DTYPE  ) ! [IN]

    endif

    return
  end subroutine ATMOS_PHY_RD_vars_restart_write

end module mod_atmos_phy_rd_vars
