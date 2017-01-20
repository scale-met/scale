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
  public :: ATMOS_PHY_RD_vars_external_in

  public :: ATMOS_PHY_RD_vars_restart_create
  public :: ATMOS_PHY_RD_vars_restart_open
  public :: ATMOS_PHY_RD_vars_restart_def_var
  public :: ATMOS_PHY_RD_vars_restart_enddef
  public :: ATMOS_PHY_RD_vars_restart_close

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  logical,               public :: ATMOS_PHY_RD_RESTART_OUTPUT                = .false.                !< output restart file?

  character(len=H_LONG), public :: ATMOS_PHY_RD_RESTART_IN_BASENAME           = ''                     !< Basename of the input  file
  logical,               public :: ATMOS_PHY_RD_RESTART_IN_POSTFIX_TIMELABEL  = .false.                !< Add timelabel to the basename of input  file?
  character(len=H_LONG), public :: ATMOS_PHY_RD_RESTART_OUT_BASENAME          = ''                     !< Basename of the output file
  logical,               public :: ATMOS_PHY_RD_RESTART_OUT_POSTFIX_TIMELABEL = .true.                 !< Add timelabel to the basename of output file?
  character(len=H_MID),  public :: ATMOS_PHY_RD_RESTART_OUT_TITLE             = 'ATMOS_PHY_RD restart' !< title    of the output file
  character(len=H_SHORT), public :: ATMOS_PHY_RD_RESTART_OUT_DTYPE             = 'DEFAULT'              !< REAL4 or REAL8

  real(RP), public, allocatable :: ATMOS_PHY_RD_RHOT_t(:,:,:)   ! tendency RHOT [K*kg/m3/s]

  real(RP), public, allocatable :: ATMOS_PHY_RD_SFLX_LW_up  (:,:) ! surface upward   longwave  flux [J/m2/s]
  real(RP), public, allocatable :: ATMOS_PHY_RD_SFLX_LW_dn  (:,:) ! surface downward longwave  flux [J/m2/s]
  real(RP), public, allocatable :: ATMOS_PHY_RD_SFLX_SW_up  (:,:) ! surface upward   shortwave flux [J/m2/s]
  real(RP), public, allocatable :: ATMOS_PHY_RD_SFLX_SW_dn  (:,:) ! surface downward shortwave flux [J/m2/s]

  real(RP), public, allocatable :: ATMOS_PHY_RD_TOAFLX_LW_up(:,:) ! TOA upward   longwave  flux [J/m2/s]
  real(RP), public, allocatable :: ATMOS_PHY_RD_TOAFLX_LW_dn(:,:) ! TOA downward longwave  flux [J/m2/s]
  real(RP), public, allocatable :: ATMOS_PHY_RD_TOAFLX_SW_up(:,:) ! TOA upward   shortwave flux [J/m2/s]
  real(RP), public, allocatable :: ATMOS_PHY_RD_TOAFLX_SW_dn(:,:) ! TOA downward shortwave flux [J/m2/s]

  real(RP), public, allocatable :: ATMOS_PHY_RD_SFLX_downall(:,:,:,:) ! surface downward flux (LW/SW,direct/diffuse) [J/m2/s]

  real(RP), public, allocatable :: ATMOS_PHY_RD_solins      (:,:) ! solar insolation flux   [J/m2/s]
  real(RP), public, allocatable :: ATMOS_PHY_RD_cosSZA      (:,:) ! cos(solar zenith angle) [0-1]

  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  integer,                private, parameter :: VMAX = 12      !< number of the variables
  integer,                private, parameter :: I_SFLX_LW_up   =  1
  integer,                private, parameter :: I_SFLX_LW_dn   =  2
  integer,                private, parameter :: I_SFLX_SW_up   =  3
  integer,                private, parameter :: I_SFLX_SW_dn   =  4
  integer,                private, parameter :: I_TOAFLX_LW_up =  5
  integer,                private, parameter :: I_TOAFLX_LW_dn =  6
  integer,                private, parameter :: I_TOAFLX_SW_up =  7
  integer,                private, parameter :: I_TOAFLX_SW_dn =  8
  integer,                private, parameter :: I_SFLX_LW_dir  =  9
  integer,                private, parameter :: I_SFLX_LW_dif  = 10
  integer,                private, parameter :: I_SFLX_SW_dir  = 11
  integer,                private, parameter :: I_SFLX_SW_dif  = 12

  character(len=H_SHORT), private            :: VAR_NAME(VMAX) !< name  of the variables
  character(len=H_MID),   private            :: VAR_DESC(VMAX) !< desc. of the variables
  character(len=H_SHORT), private            :: VAR_UNIT(VMAX) !< unit  of the variables
  integer,                private            :: VAR_ID(VMAX)   !< ID    of the variables
  integer,                private            :: restart_fid = -1  ! file ID

  data VAR_NAME / 'SFLX_LW_up',   &
                  'SFLX_LW_dn',   &
                  'SFLX_SW_up',   &
                  'SFLX_SW_dn',   &
                  'TOAFLX_LW_up', &
                  'TOAFLX_LW_dn', &
                  'TOAFLX_SW_up', &
                  'TOAFLX_SW_dn', &
                  'SFLX_LW_dir',  &
                  'SFLX_LW_dif',  &
                  'SFLX_SW_dir',  &
                  'SFLX_SW_dif'   /
  data VAR_DESC / 'surface upward   longwave  flux',   &
                  'surface downward longwave  flux',   &
                  'surface upward   shortwave flux',   &
                  'surface downward shortwave flux',   &
                  'TOA upward   longwave  flux',       &
                  'TOA downward longwave  flux',       &
                  'TOA upward   shortwave flux',       &
                  'TOA downward shortwave flux',       &
                  'sfc. down. longwave  flux direct',  &
                  'sfc. down. longwave  flux diffuse', &
                  'sfc. down. shortwave flux direct',  &
                  'sfc. down. shortwave flux diffuse'  /
  data VAR_UNIT / 'W/m2', &
                  'W/m2', &
                  'W/m2', &
                  'W/m2', &
                  'W/m2', &
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
    use scale_const, only: &
       UNDEF => CONST_UNDEF
    implicit none

    NAMELIST / PARAM_ATMOS_PHY_RD_VARS / &
       ATMOS_PHY_RD_RESTART_IN_BASENAME,           &
       ATMOS_PHY_RD_RESTART_IN_POSTFIX_TIMELABEL,  &
       ATMOS_PHY_RD_RESTART_OUTPUT,                &
       ATMOS_PHY_RD_RESTART_OUT_BASENAME,          &
       ATMOS_PHY_RD_RESTART_OUT_POSTFIX_TIMELABEL, &
       ATMOS_PHY_RD_RESTART_OUT_TITLE,             &
       ATMOS_PHY_RD_RESTART_OUT_DTYPE

    integer :: ierr
    integer :: iv
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '++++++ Module[VARS] / Categ[ATMOS PHY_RD] / Origin[SCALE-RM]'

    allocate( ATMOS_PHY_RD_RHOT_t(KA,IA,JA) )
    ATMOS_PHY_RD_RHOT_t(:,:,:) = UNDEF

    allocate( ATMOS_PHY_RD_SFLX_LW_up  (IA,JA) )
    allocate( ATMOS_PHY_RD_SFLX_LW_dn  (IA,JA) )
    allocate( ATMOS_PHY_RD_SFLX_SW_up  (IA,JA) )
    allocate( ATMOS_PHY_RD_SFLX_SW_dn  (IA,JA) )
    allocate( ATMOS_PHY_RD_TOAFLX_LW_up(IA,JA) )
    allocate( ATMOS_PHY_RD_TOAFLX_LW_dn(IA,JA) )
    allocate( ATMOS_PHY_RD_TOAFLX_SW_up(IA,JA) )
    allocate( ATMOS_PHY_RD_TOAFLX_SW_dn(IA,JA) )
    ATMOS_PHY_RD_SFLX_LW_up  (:,:) = UNDEF
    ATMOS_PHY_RD_SFLX_LW_dn  (:,:) = UNDEF
    ATMOS_PHY_RD_SFLX_SW_up  (:,:) = UNDEF
    ATMOS_PHY_RD_SFLX_SW_dn  (:,:) = UNDEF
    ATMOS_PHY_RD_TOAFLX_LW_up(:,:) = UNDEF
    ATMOS_PHY_RD_TOAFLX_LW_dn(:,:) = UNDEF
    ATMOS_PHY_RD_TOAFLX_SW_up(:,:) = UNDEF
    ATMOS_PHY_RD_TOAFLX_SW_dn(:,:) = UNDEF

    allocate( ATMOS_PHY_RD_SFLX_downall(IA,JA,2,2) )
    ATMOS_PHY_RD_SFLX_downall(:,:,:,:) = UNDEF

    allocate( ATMOS_PHY_RD_solins(IA,JA) )
    allocate( ATMOS_PHY_RD_cosSZA(IA,JA) )
    ATMOS_PHY_RD_solins(:,:) = UNDEF
    ATMOS_PHY_RD_cosSZA(:,:) = UNDEF

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_ATMOS_PHY_RD_VARS,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_ATMOS_PHY_RD_VARS. Check!'
       call PRC_MPIstop
    endif
    if( IO_LNML ) write(IO_FID_LOG,nml=PARAM_ATMOS_PHY_RD_VARS)

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** [ATMOS_PHY_RD] prognostic/diagnostic variables'
    if( IO_L ) write(IO_FID_LOG,'(1x,A,A15,A,A32,3(A))') &
               '***       |','VARNAME        ','|', 'DESCRIPTION                     ','[', 'UNIT            ',']'
    do iv = 1, VMAX
       if( IO_L ) write(IO_FID_LOG,'(1x,A,i3,A,A15,A,A32,3(A))') &
                  '*** NO.',iv,'|',VAR_NAME(iv),'|',VAR_DESC(iv),'[',VAR_UNIT(iv),']'
    enddo

    if( IO_L ) write(IO_FID_LOG,*)
    if ( ATMOS_PHY_RD_RESTART_IN_BASENAME /= '' ) then
       if( IO_L ) write(IO_FID_LOG,*) '*** Restart input?  : YES, file = ', trim(ATMOS_PHY_RD_RESTART_IN_BASENAME)
       if( IO_L ) write(IO_FID_LOG,*) '*** Add timelabel?  : ', ATMOS_PHY_RD_RESTART_IN_POSTFIX_TIMELABEL
    else
       if( IO_L ) write(IO_FID_LOG,*) '*** Restart input?  : NO'
    endif
    if (       ATMOS_PHY_RD_RESTART_OUTPUT             &
         .AND. ATMOS_PHY_RD_RESTART_OUT_BASENAME /= '' ) then
       if( IO_L ) write(IO_FID_LOG,*) '*** Restart output? : YES, file = ', trim(ATMOS_PHY_RD_RESTART_OUT_BASENAME)
       if( IO_L ) write(IO_FID_LOG,*) '*** Add timelabel?  : ', ATMOS_PHY_RD_RESTART_OUT_POSTFIX_TIMELABEL
    else
       if( IO_L ) write(IO_FID_LOG,*) '*** Restart output? : NO'
       ATMOS_PHY_RD_RESTART_OUTPUT = .false.
    endif

    return
  end subroutine ATMOS_PHY_RD_vars_setup

  !-----------------------------------------------------------------------------
  !> HALO Communication
  subroutine ATMOS_PHY_RD_vars_fillhalo
    use scale_comm, only: &
       COMM_vars8, &
       COMM_wait
    implicit none

    integer :: n ,iw, id
    !---------------------------------------------------------------------------

    call COMM_vars8( ATMOS_PHY_RD_SFLX_LW_up  (:,:),  1 )
    call COMM_vars8( ATMOS_PHY_RD_SFLX_LW_dn  (:,:),  2 )
    call COMM_vars8( ATMOS_PHY_RD_SFLX_SW_up  (:,:),  3 )
    call COMM_vars8( ATMOS_PHY_RD_SFLX_SW_dn  (:,:),  4 )
    call COMM_vars8( ATMOS_PHY_RD_TOAFLX_LW_up(:,:),  5 )
    call COMM_vars8( ATMOS_PHY_RD_TOAFLX_LW_dn(:,:),  6 )
    call COMM_vars8( ATMOS_PHY_RD_TOAFLX_SW_up(:,:),  7 )
    call COMM_vars8( ATMOS_PHY_RD_TOAFLX_SW_dn(:,:),  8 )

    n = 8
    do id = 1, 2 ! direct/diffuse
    do iw = 1, 2 ! SW/LW
       n = n + 1
       call COMM_vars8( ATMOS_PHY_RD_SFLX_downall(:,:,iw,id), n )
    enddo
    enddo

    call COMM_wait ( ATMOS_PHY_RD_SFLX_LW_up  (:,:),  1 )
    call COMM_wait ( ATMOS_PHY_RD_SFLX_LW_dn  (:,:),  2 )
    call COMM_wait ( ATMOS_PHY_RD_SFLX_SW_up  (:,:),  3 )
    call COMM_wait ( ATMOS_PHY_RD_SFLX_SW_dn  (:,:),  4 )
    call COMM_wait ( ATMOS_PHY_RD_TOAFLX_LW_up(:,:),  5 )
    call COMM_wait ( ATMOS_PHY_RD_TOAFLX_LW_dn(:,:),  6 )
    call COMM_wait ( ATMOS_PHY_RD_TOAFLX_SW_up(:,:),  7 )
    call COMM_wait ( ATMOS_PHY_RD_TOAFLX_SW_dn(:,:),  8 )

    n = 8
    do id = 1, 2 ! direct/diffuse
    do iw = 1, 2 ! SW/LW
       n = n + 1
       call COMM_wait ( ATMOS_PHY_RD_SFLX_downall(:,:,iw,id), n )
    enddo
    enddo

    return
  end subroutine ATMOS_PHY_RD_vars_fillhalo

  !-----------------------------------------------------------------------------
  !> Open restart file for read
  subroutine ATMOS_PHY_RD_vars_restart_open
    use scale_time, only: &
       TIME_gettimelabel
    use scale_fileio, only: &
       FILEIO_open
    implicit none

    character(len=19)     :: timelabel
    character(len=H_LONG) :: basename
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** Input restart file (ATMOS_PHY_RD) ***'

    if ( ATMOS_PHY_RD_RESTART_IN_BASENAME /= '' ) then

       if ( ATMOS_PHY_RD_RESTART_IN_POSTFIX_TIMELABEL ) then
          call TIME_gettimelabel( timelabel )
          basename = trim(ATMOS_PHY_RD_RESTART_IN_BASENAME)//'_'//trim(timelabel)
       else
          basename = trim(ATMOS_PHY_RD_RESTART_IN_BASENAME)
       endif

       if( IO_L ) write(IO_FID_LOG,*) '*** basename: ', trim(basename)

       call FILEIO_open( restart_fid, basename )
    else
       if( IO_L ) write(IO_FID_LOG,*) '*** restart file for ATMOS_PHY_RD is not specified.'
    endif

    return
  end subroutine ATMOS_PHY_RD_vars_restart_open


  !-----------------------------------------------------------------------------
  !> Read restart
  subroutine ATMOS_PHY_RD_vars_restart_read
    use scale_fileio, only: &
       FILEIO_read, &
       FILEIO_flush
    use scale_rm_statistics, only: &
       STAT_total
    implicit none

    real(RP) :: total
    !---------------------------------------------------------------------------

    if ( restart_fid .NE. -1 ) then

       call FILEIO_read( ATMOS_PHY_RD_SFLX_LW_up(:,:),           & ! [OUT]
                         restart_fid, VAR_NAME(1) , 'XY', step=1 ) ! [IN]
       call FILEIO_read( ATMOS_PHY_RD_SFLX_LW_dn(:,:),           & ! [OUT]
                         restart_fid, VAR_NAME(2) , 'XY', step=1 ) ! [IN]
       call FILEIO_read( ATMOS_PHY_RD_SFLX_SW_up(:,:),           & ! [OUT]
                         restart_fid, VAR_NAME(3) , 'XY', step=1 ) ! [IN]
       call FILEIO_read( ATMOS_PHY_RD_SFLX_SW_dn(:,:),           & ! [OUT]
                         restart_fid, VAR_NAME(4) , 'XY', step=1 ) ! [IN]
       call FILEIO_read( ATMOS_PHY_RD_TOAFLX_LW_up(:,:),         & ! [OUT]
                         restart_fid, VAR_NAME(5) , 'XY', step=1 ) ! [IN]
       call FILEIO_read( ATMOS_PHY_RD_TOAFLX_LW_dn(:,:),         & ! [OUT]
                         restart_fid, VAR_NAME(6) , 'XY', step=1 ) ! [IN]
       call FILEIO_read( ATMOS_PHY_RD_TOAFLX_SW_up(:,:),         & ! [OUT]
                         restart_fid, VAR_NAME(7) , 'XY', step=1 ) ! [IN]
       call FILEIO_read( ATMOS_PHY_RD_TOAFLX_SW_dn(:,:),         & ! [OUT]
                         restart_fid, VAR_NAME(8) , 'XY', step=1 ) ! [IN]
       call FILEIO_read( ATMOS_PHY_RD_SFLX_downall(:,:,1,1),     & ! [OUT]
                         restart_fid, VAR_NAME(9) , 'XY', step=1 ) ! [IN]
       call FILEIO_read( ATMOS_PHY_RD_SFLX_downall(:,:,1,2),     & ! [OUT]
                         restart_fid, VAR_NAME(10), 'XY', step=1 ) ! [IN]
       call FILEIO_read( ATMOS_PHY_RD_SFLX_downall(:,:,2,1),     & ! [OUT]
                         restart_fid, VAR_NAME(11), 'XY', step=1 ) ! [IN]
       call FILEIO_read( ATMOS_PHY_RD_SFLX_downall(:,:,2,2),     & ! [OUT]
                         restart_fid, VAR_NAME(12), 'XY', step=1 ) ! [IN]

       if ( IO_AGGREGATE ) then
          call FILEIO_flush( restart_fid )
          ! halos have been read from file
       else
          call ATMOS_PHY_RD_vars_fillhalo
       end if

       call STAT_total( total, ATMOS_PHY_RD_SFLX_LW_up  (:,:),     VAR_NAME(1)  )
       call STAT_total( total, ATMOS_PHY_RD_SFLX_LW_dn  (:,:),     VAR_NAME(2)  )
       call STAT_total( total, ATMOS_PHY_RD_SFLX_SW_up  (:,:),     VAR_NAME(3)  )
       call STAT_total( total, ATMOS_PHY_RD_SFLX_SW_dn  (:,:),     VAR_NAME(4)  )
       call STAT_total( total, ATMOS_PHY_RD_TOAFLX_LW_up(:,:),     VAR_NAME(5)  )
       call STAT_total( total, ATMOS_PHY_RD_TOAFLX_LW_dn(:,:),     VAR_NAME(6)  )
       call STAT_total( total, ATMOS_PHY_RD_TOAFLX_SW_up(:,:),     VAR_NAME(7)  )
       call STAT_total( total, ATMOS_PHY_RD_TOAFLX_SW_dn(:,:),     VAR_NAME(8)  )
       call STAT_total( total, ATMOS_PHY_RD_SFLX_downall(:,:,1,1), VAR_NAME(9)  )
       call STAT_total( total, ATMOS_PHY_RD_SFLX_downall(:,:,1,2), VAR_NAME(10) )
       call STAT_total( total, ATMOS_PHY_RD_SFLX_downall(:,:,2,1), VAR_NAME(11) )
       call STAT_total( total, ATMOS_PHY_RD_SFLX_downall(:,:,2,2), VAR_NAME(12) )
    else
       if( IO_L ) write(IO_FID_LOG,*) '*** invalid restart file ID for ATMOS_PHY_RD.'
    endif

    return
  end subroutine ATMOS_PHY_RD_vars_restart_read

  !-----------------------------------------------------------------------------
  !> Input from External I/O
  subroutine ATMOS_PHY_RD_vars_external_in( &
      init_value_in  )
    implicit none

    real(RP), intent(in) :: init_value_in
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** External Input (PHY_RD) ***'

    ATMOS_PHY_RD_SFLX_LW_up  (:,:)     = init_value_in
    ATMOS_PHY_RD_SFLX_LW_dn  (:,:)     = init_value_in
    ATMOS_PHY_RD_SFLX_SW_up  (:,:)     = init_value_in
    ATMOS_PHY_RD_SFLX_SW_dn  (:,:)     = init_value_in
    ATMOS_PHY_RD_TOAFLX_LW_up(:,:)     = init_value_in
    ATMOS_PHY_RD_TOAFLX_LW_dn(:,:)     = init_value_in
    ATMOS_PHY_RD_TOAFLX_SW_up(:,:)     = init_value_in
    ATMOS_PHY_RD_TOAFLX_SW_dn(:,:)     = init_value_in
    ATMOS_PHY_RD_SFLX_downall(:,:,:,:) = init_value_in

    return
  end subroutine ATMOS_PHY_RD_vars_external_in

  !-----------------------------------------------------------------------------
  !> Create restart file
  subroutine ATMOS_PHY_RD_vars_restart_create
    use scale_time, only: &
       TIME_gettimelabel
    use scale_fileio, only: &
       FILEIO_create
    use scale_rm_statistics, only: &
       STAT_total
    implicit none

    character(len=19)     :: timelabel
    character(len=H_LONG) :: basename

    !---------------------------------------------------------------------------

    if ( ATMOS_PHY_RD_RESTART_OUT_BASENAME /= '' ) then

       if( IO_L ) write(IO_FID_LOG,*)
       if( IO_L ) write(IO_FID_LOG,*) '*** Output restart file (ATMOS_PHY_AE) ***'

       if ( ATMOS_PHY_RD_RESTART_OUT_POSTFIX_TIMELABEL ) then
          call TIME_gettimelabel( timelabel )
          basename = trim(ATMOS_PHY_RD_RESTART_OUT_BASENAME)//'_'//trim(timelabel)
       else
          basename = trim(ATMOS_PHY_RD_RESTART_OUT_BASENAME)
       endif

       if( IO_L ) write(IO_FID_LOG,*) '*** basename: ', trim(basename)

       call FILEIO_create( restart_fid,                                                             & ! [OUT]
                           basename, ATMOS_PHY_RD_RESTART_OUT_TITLE, ATMOS_PHY_RD_RESTART_OUT_DTYPE ) ! [IN]

    endif

    return
  end subroutine ATMOS_PHY_RD_vars_restart_create

  !-----------------------------------------------------------------------------
  !> Exit netCDF define mode
  subroutine ATMOS_PHY_RD_vars_restart_enddef
    use scale_fileio, only: &
       FILEIO_enddef
    implicit none

    if ( restart_fid .NE. -1 ) then
       call FILEIO_enddef( restart_fid ) ! [IN]
    endif

    return
  end subroutine ATMOS_PHY_RD_vars_restart_enddef

  !-----------------------------------------------------------------------------
  !> Close restart file
  subroutine ATMOS_PHY_RD_vars_restart_close
    use scale_fileio, only: &
       FILEIO_close
    implicit none

    if ( restart_fid .NE. -1 ) then
       call FILEIO_close( restart_fid ) ! [IN]
       restart_fid = -1
    endif

    return
  end subroutine ATMOS_PHY_RD_vars_restart_close

  !-----------------------------------------------------------------------------
  !> Define variables in restart file
  subroutine ATMOS_PHY_RD_vars_restart_def_var
    use scale_fileio, only: &
       FILEIO_def_var
    implicit none

    character(len=20)     :: timelabel
    character(len=H_LONG) :: basename

    !---------------------------------------------------------------------------

    if ( restart_fid .NE. -1 ) then

       call FILEIO_def_var( restart_fid, VAR_ID(1), VAR_NAME(1) , VAR_DESC(1) , VAR_UNIT(1) , &
                            'XY', ATMOS_PHY_RD_RESTART_OUT_DTYPE  ) ! [IN]
       call FILEIO_def_var( restart_fid, VAR_ID(2), VAR_NAME(2) , VAR_DESC(2) , VAR_UNIT(2) , &
                            'XY', ATMOS_PHY_RD_RESTART_OUT_DTYPE  ) ! [IN]
       call FILEIO_def_var( restart_fid, VAR_ID(3), VAR_NAME(3) , VAR_DESC(3) , VAR_UNIT(3) , &
                            'XY', ATMOS_PHY_RD_RESTART_OUT_DTYPE  ) ! [IN]
       call FILEIO_def_var( restart_fid, VAR_ID(4), VAR_NAME(4) , VAR_DESC(4) , VAR_UNIT(4) , &
                            'XY', ATMOS_PHY_RD_RESTART_OUT_DTYPE  ) ! [IN]
       call FILEIO_def_var( restart_fid, VAR_ID(5), VAR_NAME(5) , VAR_DESC(5) , VAR_UNIT(5) , &
                            'XY', ATMOS_PHY_RD_RESTART_OUT_DTYPE  ) ! [IN]
       call FILEIO_def_var( restart_fid, VAR_ID(6), VAR_NAME(6) , VAR_DESC(6) , VAR_UNIT(6) , &
                            'XY', ATMOS_PHY_RD_RESTART_OUT_DTYPE  ) ! [IN]
       call FILEIO_def_var( restart_fid, VAR_ID(7), VAR_NAME(7) , VAR_DESC(7) , VAR_UNIT(7) , &
                            'XY', ATMOS_PHY_RD_RESTART_OUT_DTYPE  ) ! [IN]
       call FILEIO_def_var( restart_fid, VAR_ID(8), VAR_NAME(8) , VAR_DESC(8) , VAR_UNIT(8) , &
                            'XY', ATMOS_PHY_RD_RESTART_OUT_DTYPE  ) ! [IN]
       call FILEIO_def_var( restart_fid, VAR_ID(9), VAR_NAME(9) , VAR_DESC(9) , VAR_UNIT(9) , &
                            'XY', ATMOS_PHY_RD_RESTART_OUT_DTYPE  ) ! [IN]
       call FILEIO_def_var( restart_fid, VAR_ID(10), VAR_NAME(10), VAR_DESC(10), VAR_UNIT(10), &
                            'XY', ATMOS_PHY_RD_RESTART_OUT_DTYPE  ) ! [IN]
       call FILEIO_def_var( restart_fid, VAR_ID(11), VAR_NAME(11), VAR_DESC(11), VAR_UNIT(11), &
                            'XY', ATMOS_PHY_RD_RESTART_OUT_DTYPE  ) ! [IN]
       call FILEIO_def_var( restart_fid, VAR_ID(12), VAR_NAME(12), VAR_DESC(12), VAR_UNIT(12), &
                            'XY', ATMOS_PHY_RD_RESTART_OUT_DTYPE  ) ! [IN]

    endif

    return
  end subroutine ATMOS_PHY_RD_vars_restart_def_var

  !-----------------------------------------------------------------------------
  !> Write variables to restart file
  subroutine ATMOS_PHY_RD_vars_restart_write
    use scale_fileio, only: &
       FILEIO_write_var
    use scale_rm_statistics, only: &
       STAT_total
    implicit none

    real(RP) :: total
    !---------------------------------------------------------------------------

    if ( restart_fid .NE. -1 ) then

       call ATMOS_PHY_RD_vars_fillhalo

       call STAT_total( total, ATMOS_PHY_RD_SFLX_LW_up  (:,:),     VAR_NAME(1)  )
       call STAT_total( total, ATMOS_PHY_RD_SFLX_LW_dn  (:,:),     VAR_NAME(2)  )
       call STAT_total( total, ATMOS_PHY_RD_SFLX_SW_up  (:,:),     VAR_NAME(3)  )
       call STAT_total( total, ATMOS_PHY_RD_SFLX_SW_dn  (:,:),     VAR_NAME(4)  )
       call STAT_total( total, ATMOS_PHY_RD_TOAFLX_LW_up(:,:),     VAR_NAME(5)  )
       call STAT_total( total, ATMOS_PHY_RD_TOAFLX_LW_dn(:,:),     VAR_NAME(6)  )
       call STAT_total( total, ATMOS_PHY_RD_TOAFLX_SW_up(:,:),     VAR_NAME(7)  )
       call STAT_total( total, ATMOS_PHY_RD_TOAFLX_SW_dn(:,:),     VAR_NAME(8)  )
       call STAT_total( total, ATMOS_PHY_RD_SFLX_downall(:,:,1,1), VAR_NAME(9)  )
       call STAT_total( total, ATMOS_PHY_RD_SFLX_downall(:,:,1,2), VAR_NAME(10) )
       call STAT_total( total, ATMOS_PHY_RD_SFLX_downall(:,:,2,1), VAR_NAME(11) )
       call STAT_total( total, ATMOS_PHY_RD_SFLX_downall(:,:,2,2), VAR_NAME(12) )

       call FILEIO_write_var( restart_fid, VAR_ID(1), ATMOS_PHY_RD_SFLX_LW_up(:,:), &
                              VAR_NAME(1) , 'XY'  ) ! [IN]
       call FILEIO_write_var( restart_fid, VAR_ID(2), ATMOS_PHY_RD_SFLX_LW_dn(:,:), &
                              VAR_NAME(2) , 'XY'  ) ! [IN]
       call FILEIO_write_var( restart_fid, VAR_ID(3), ATMOS_PHY_RD_SFLX_SW_up(:,:), &
                              VAR_NAME(3) , 'XY'  ) ! [IN]
       call FILEIO_write_var( restart_fid, VAR_ID(4), ATMOS_PHY_RD_SFLX_SW_dn(:,:), &
                              VAR_NAME(4) , 'XY'  ) ! [IN]
       call FILEIO_write_var( restart_fid, VAR_ID(5), ATMOS_PHY_RD_TOAFLX_LW_up(:,:), &
                              VAR_NAME(5) , 'XY'  ) ! [IN]
       call FILEIO_write_var( restart_fid, VAR_ID(6), ATMOS_PHY_RD_TOAFLX_LW_dn(:,:), &
                              VAR_NAME(6) , 'XY'  ) ! [IN]
       call FILEIO_write_var( restart_fid, VAR_ID(7), ATMOS_PHY_RD_TOAFLX_SW_up(:,:), &
                              VAR_NAME(7) , 'XY'  ) ! [IN]
       call FILEIO_write_var( restart_fid, VAR_ID(8), ATMOS_PHY_RD_TOAFLX_SW_dn(:,:), &
                              VAR_NAME(8) , 'XY'  ) ! [IN]
       call FILEIO_write_var( restart_fid, VAR_ID(9), ATMOS_PHY_RD_SFLX_downall(:,:,1,1), &
                              VAR_NAME(9) , 'XY'  ) ! [IN]
       call FILEIO_write_var( restart_fid, VAR_ID(10), ATMOS_PHY_RD_SFLX_downall(:,:,1,2), &
                              VAR_NAME(10), 'XY'  ) ! [IN]
       call FILEIO_write_var( restart_fid, VAR_ID(11), ATMOS_PHY_RD_SFLX_downall(:,:,2,1), &
                              VAR_NAME(11), 'XY'  ) ! [IN]
       call FILEIO_write_var( restart_fid, VAR_ID(12), ATMOS_PHY_RD_SFLX_downall(:,:,2,2), &
                              VAR_NAME(12), 'XY'  ) ! [IN]

    endif

    return
  end subroutine ATMOS_PHY_RD_vars_restart_write

end module mod_atmos_phy_rd_vars
