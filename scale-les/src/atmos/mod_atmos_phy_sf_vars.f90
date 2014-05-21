!-------------------------------------------------------------------------------
!> module ATMOSPHERIC Surface Variables
!!
!! @par Description
!!          Container for mod_atmos_phy_sf
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2012-03-27 (H.Yashiro)  [new]
!! @li      2013-08-31 (T.Yamaura)  [mod]
!<
!-------------------------------------------------------------------------------
module mod_atmos_phy_sf_vars
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
  public :: ATMOS_PHY_SF_vars_setup
  public :: ATMOS_PHY_SF_vars_fillhalo
  public :: ATMOS_PHY_SF_vars_restart_read
  public :: ATMOS_PHY_SF_vars_restart_write

  !-----------------------------------------------------------------------------
  !
  !++ included parameters
  !
  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  logical, public, save ::  ATMOS_PHY_SF_sw_restart = .false.

  real(RP), public, allocatable :: ATMOS_PHY_SF_DENS_t(:,:,:)   ! tendency DENS [kg/m3/s]
  real(RP), public, allocatable :: ATMOS_PHY_SF_MOMZ_t(:,:,:)   ! tendency MOMZ [kg/m2/s2]
  real(RP), public, allocatable :: ATMOS_PHY_SF_MOMX_t(:,:,:)   ! tendency MOMX [kg/m2/s2]
  real(RP), public, allocatable :: ATMOS_PHY_SF_MOMY_t(:,:,:)   ! tendency MOMY [kg/m2/s2]
  real(RP), public, allocatable :: ATMOS_PHY_SF_RHOT_t(:,:,:)   ! tendency RHOT [K*kg/m3/s]
  real(RP), public, allocatable :: ATMOS_PHY_SF_QTRC_t(:,:,:,:) ! tendency QTRC [kg/kg/s]

  real(RP), public, allocatable :: ATMOS_PHY_SF_PREC   (:,:) ! surface precipitation rate [kg/m2/s]
  real(RP), public, allocatable :: ATMOS_PHY_SF_SWD    (:,:) ! downward short-wave radiation flux (upward positive) [W/m2]
  real(RP), public, allocatable :: ATMOS_PHY_SF_LWD    (:,:) ! downward long-wave  radiation flux (upward positive) [W/m2]
  real(RP), public, allocatable :: ATMOS_PHY_SF_ZMFLX  (:,:) ! z-momentum flux [kg/m2/s]
  real(RP), public, allocatable :: ATMOS_PHY_SF_XMFLX  (:,:) ! x-momentum flux [kg/m2/s]
  real(RP), public, allocatable :: ATMOS_PHY_SF_YMFLX  (:,:) ! y-momentum flux [kg/m2/s]
  real(RP), public, allocatable :: ATMOS_PHY_SF_POTTFLX(:,:) ! potential temperature flux [K*kg/m3]
  real(RP), public, allocatable :: ATMOS_PHY_SF_SWUFLX (:,:) ! upward short-wave radiation flux (upward positive) [W/m2]
  real(RP), public, allocatable :: ATMOS_PHY_SF_LWUFLX (:,:) ! upward long-wave  radiation flux (upward positive) [W/m2]
  real(RP), public, allocatable :: ATMOS_PHY_SF_SHFLX  (:,:) ! sensible heat flux (upward positive) [W/m2]
  real(RP), public, allocatable :: ATMOS_PHY_SF_LHFLX  (:,:) ! latent   heat flux (upward positive) [W/m2]
  real(RP), public, allocatable :: ATMOS_PHY_SF_QVFLX  (:,:) ! moisture flux for atmosphere [kg/m2/s]

  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  logical,                private, save      :: ATMOS_PHY_SF_RESTART_OUTPUT        = .false.
  character(len=H_LONG),  private, save      :: ATMOS_PHY_SF_RESTART_IN_BASENAME   = ''
  character(len=H_LONG),  private, save      :: ATMOS_PHY_SF_RESTART_OUT_BASENAME  = ''
  character(len=H_MID),   private, save      :: ATMOS_PHY_SF_RESTART_OUT_TITLE     = 'ATMOS_PHY_SF restart'
  character(len=H_MID),   private, save      :: ATMOS_PHY_SF_RESTART_OUT_DTYPE     = 'DEFAULT'              !< REAL4 or REAL8

  integer,                private, parameter :: VMAX = 12      !< number of the variables
  character(len=H_SHORT), private, save      :: VAR_NAME(VMAX) !< name  of the variables
  character(len=H_MID),   private, save      :: VAR_DESC(VMAX) !< desc. of the variables
  character(len=H_SHORT), private, save      :: VAR_UNIT(VMAX) !< unit  of the variables

  data VAR_NAME / 'PREC',    &
                  'SWD',     &
                  'LWD',     &
                  'ZMFLX',   &
                  'XMFLX',   &
                  'YMFLX',   &
                  'POTTFLX', &
                  'SWUFLX',  &
                  'LWUFLX',  &
                  'SHFLX',   &
                  'LHFLX',   &
                  'QVFLX'    /
  data VAR_DESC / 'surface precipitation rate',         &
                  'downward short-wave radiation flux', &
                  'downward long-wave  radiation flux', &
                  'z-momentum flux',                    &
                  'x-momentum flux',                    &
                  'y-momentum flux',                    &
                  'potential temperature flux',         &
                  'upward short-wave radiation flux',   &
                  'upward long-wave  radiation flux',   &
                  'sensible heat flux',                 &
                  'latent   heat flux',                 &
                  'moisture flux for atmosphere'        /
  data VAR_UNIT / 'kg/m2/s', &
                  'W/m2',    &
                  'W/m2',    &
                  'kg/m2/s', &
                  'kg/m2/s', &
                  'kg/m2/s', &
                  'K*kg/m3', &
                  'W/m2',    &
                  'W/m2',    &
                  'W/m2',    &
                  'W/m2',    &
                  'kg/m2/s'  /

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine ATMOS_PHY_SF_vars_setup
    use scale_process, only: &
       PRC_MPIstop
    implicit none

    NAMELIST / PARAM_ATMOS_PHY_SF_VARS / &
       ATMOS_PHY_SF_RESTART_IN_BASENAME, &
       ATMOS_PHY_SF_RESTART_OUTPUT,      &
       ATMOS_PHY_SF_RESTART_OUT_BASENAME

    integer :: v, ierr
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[ATMOS_PHY_SF VARS]/Categ[ATMOS]'
    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_ATMOS_PHY_SF_VARS,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_ATMOS_PHY_SF_VARS. Check!'
       call PRC_MPIstop
    endif
    if( IO_L ) write(IO_FID_LOG,nml=PARAM_ATMOS_PHY_SF_VARS)

    allocate( ATMOS_PHY_SF_PREC   (IA,JA) )
    allocate( ATMOS_PHY_SF_SWD    (IA,JA) )
    allocate( ATMOS_PHY_SF_LWD    (IA,JA) )
    allocate( ATMOS_PHY_SF_ZMFLX  (IA,JA) )
    allocate( ATMOS_PHY_SF_XMFLX  (IA,JA) )
    allocate( ATMOS_PHY_SF_YMFLX  (IA,JA) )
    allocate( ATMOS_PHY_SF_POTTFLX(IA,JA) )
    allocate( ATMOS_PHY_SF_SWUFLX (IA,JA) )
    allocate( ATMOS_PHY_SF_LWUFLX (IA,JA) )
    allocate( ATMOS_PHY_SF_SHFLX  (IA,JA) )
    allocate( ATMOS_PHY_SF_LHFLX  (IA,JA) )
    allocate( ATMOS_PHY_SF_QVFLX  (IA,JA) )

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** [ATMOS_PHY_SF] prognostic/diagnostic variables'
    if( IO_L ) write(IO_FID_LOG,'(1x,A,A8,A,A32,3(A))') &
               '***       |',' VARNAME','|', 'DESCRIPTION                     ','[', 'UNIT            ',']'
    do v = 1, VMAX
       if( IO_L ) write(IO_FID_LOG,'(1x,A,i3,A,A8,A,A32,3(A))') &
                  '*** NO.',v,'|',trim(VAR_NAME(v)),'|',VAR_DESC(v),'[',VAR_UNIT(v),']'
    enddo

    ! restart switch
    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) 'Output...'
    if ( ATMOS_PHY_SF_RESTART_OUTPUT ) then
       if( IO_L ) write(IO_FID_LOG,*) '  Restart output (ATMOS_PHY_SF) : YES'
       ATMOS_PHY_SF_sw_restart = .true.
    else
       if( IO_L ) write(IO_FID_LOG,*) '  Restart output (ATMOS_PHY_SF) : NO'
       ATMOS_PHY_SF_sw_restart = .false.
    endif
    if( IO_L ) write(IO_FID_LOG,*)

    return
  end subroutine ATMOS_PHY_SF_vars_setup

  !-----------------------------------------------------------------------------
  !> Communication
  subroutine ATMOS_PHY_SF_vars_fillhalo
    use scale_comm, only: &
       COMM_vars8, &
       COMM_wait
    implicit none
    !---------------------------------------------------------------------------

    ! fill IHALO & JHALO
    call COMM_vars8( ATMOS_PHY_SF_PREC   (:,:),  1 )
    call COMM_vars8( ATMOS_PHY_SF_SWD    (:,:),  2 )
    call COMM_vars8( ATMOS_PHY_SF_LWD    (:,:),  3 )
    call COMM_vars8( ATMOS_PHY_SF_ZMFLX  (:,:),  4 )
    call COMM_vars8( ATMOS_PHY_SF_XMFLX  (:,:),  5 )
    call COMM_vars8( ATMOS_PHY_SF_YMFLX  (:,:),  6 )
    call COMM_vars8( ATMOS_PHY_SF_POTTFLX(:,:),  7 )
    call COMM_vars8( ATMOS_PHY_SF_SWUFLX (:,:),  8 )
    call COMM_vars8( ATMOS_PHY_SF_LWUFLX (:,:),  9 )
    call COMM_vars8( ATMOS_PHY_SF_SHFLX  (:,:), 10 )
    call COMM_vars8( ATMOS_PHY_SF_LHFLX  (:,:), 11 )
    call COMM_vars8( ATMOS_PHY_SF_QVFLX  (:,:), 12 )
    call COMM_wait ( ATMOS_PHY_SF_PREC   (:,:),  1 )
    call COMM_wait ( ATMOS_PHY_SF_SWD    (:,:),  2 )
    call COMM_wait ( ATMOS_PHY_SF_LWD    (:,:),  3 )
    call COMM_wait ( ATMOS_PHY_SF_ZMFLX  (:,:),  4 )
    call COMM_wait ( ATMOS_PHY_SF_XMFLX  (:,:),  5 )
    call COMM_wait ( ATMOS_PHY_SF_YMFLX  (:,:),  6 )
    call COMM_wait ( ATMOS_PHY_SF_POTTFLX(:,:),  7 )
    call COMM_wait ( ATMOS_PHY_SF_SWUFLX (:,:),  8 )
    call COMM_wait ( ATMOS_PHY_SF_LWUFLX (:,:),  9 )
    call COMM_wait ( ATMOS_PHY_SF_SHFLX  (:,:), 10 )
    call COMM_wait ( ATMOS_PHY_SF_LHFLX  (:,:), 11 )
    call COMM_wait ( ATMOS_PHY_SF_QVFLX  (:,:), 12 )

    return
  end subroutine ATMOS_PHY_SF_vars_fillhalo

  !-----------------------------------------------------------------------------
  !> Read restart
  subroutine ATMOS_PHY_SF_vars_restart_read
    use scale_fileio, only: &
       FILEIO_read
    use scale_const, only: &
       CONST_UNDEF
    implicit none
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** Input restart file (ATMOS_PHY_SF) ***'

    if ( ATMOS_PHY_SF_RESTART_IN_BASENAME /= '' ) then

       call FILEIO_read( ATMOS_PHY_SF_PREC   (:,:),                                   & ! [OUT]
                         ATMOS_PHY_SF_RESTART_IN_BASENAME, VAR_NAME( 1), 'XY', step=1 ) ! [IN]
       call FILEIO_read( ATMOS_PHY_SF_SWD    (:,:),                                   & ! [OUT]
                         ATMOS_PHY_SF_RESTART_IN_BASENAME, VAR_NAME( 2), 'XY', step=1 ) ! [IN]
       call FILEIO_read( ATMOS_PHY_SF_LWD    (:,:),                                   & ! [OUT]
                         ATMOS_PHY_SF_RESTART_IN_BASENAME, VAR_NAME( 3), 'XY', step=1 ) ! [IN]
       call FILEIO_read( ATMOS_PHY_SF_ZMFLX  (:,:),                                   & ! [OUT]
                         ATMOS_PHY_SF_RESTART_IN_BASENAME, VAR_NAME( 4), 'XY', step=1 ) ! [IN]
       call FILEIO_read( ATMOS_PHY_SF_XMFLX  (:,:),                                   & ! [OUT]
                         ATMOS_PHY_SF_RESTART_IN_BASENAME, VAR_NAME( 5), 'XY', step=1 ) ! [IN]
       call FILEIO_read( ATMOS_PHY_SF_YMFLX  (:,:),                                   & ! [OUT]
                         ATMOS_PHY_SF_RESTART_IN_BASENAME, VAR_NAME( 6), 'XY', step=1 ) ! [IN]
       call FILEIO_read( ATMOS_PHY_SF_POTTFLX(:,:),                                   & ! [OUT]
                         ATMOS_PHY_SF_RESTART_IN_BASENAME, VAR_NAME( 7), 'XY', step=1 ) ! [IN]
       call FILEIO_read( ATMOS_PHY_SF_SWUFLX (:,:),                                   & ! [OUT]
                         ATMOS_PHY_SF_RESTART_IN_BASENAME, VAR_NAME( 8), 'XY', step=1 ) ! [IN]
       call FILEIO_read( ATMOS_PHY_SF_LWUFLX (:,:),                                   & ! [OUT]
                         ATMOS_PHY_SF_RESTART_IN_BASENAME, VAR_NAME( 9), 'XY', step=1 ) ! [IN]
       call FILEIO_read( ATMOS_PHY_SF_SHFLX  (:,:),                                   & ! [OUT]
                         ATMOS_PHY_SF_RESTART_IN_BASENAME, VAR_NAME(10), 'XY', step=1 ) ! [IN]
       call FILEIO_read( ATMOS_PHY_SF_LHFLX  (:,:),                                   & ! [OUT]
                         ATMOS_PHY_SF_RESTART_IN_BASENAME, VAR_NAME(11), 'XY', step=1 ) ! [IN]
       call FILEIO_read( ATMOS_PHY_SF_QVFLX  (:,:),                                   & ! [OUT]
                         ATMOS_PHY_SF_RESTART_IN_BASENAME, VAR_NAME(12), 'XY', step=1 ) ! [IN]

       call ATMOS_PHY_SF_vars_fillhalo

    else

       if( IO_L ) write(IO_FID_LOG,*) '*** restart file for ATMOS_PHY_SF is not specified.'

       ATMOS_PHY_SF_PREC   (:,:) = CONST_UNDEF
       ATMOS_PHY_SF_SWD    (:,:) = CONST_UNDEF
       ATMOS_PHY_SF_LWD    (:,:) = CONST_UNDEF
       ATMOS_PHY_SF_ZMFLX  (:,:) = CONST_UNDEF
       ATMOS_PHY_SF_XMFLX  (:,:) = CONST_UNDEF
       ATMOS_PHY_SF_YMFLX  (:,:) = CONST_UNDEF
       ATMOS_PHY_SF_POTTFLX(:,:) = CONST_UNDEF
       ATMOS_PHY_SF_SWUFLX (:,:) = CONST_UNDEF
       ATMOS_PHY_SF_LWUFLX (:,:) = CONST_UNDEF
       ATMOS_PHY_SF_SHFLX  (:,:) = CONST_UNDEF
       ATMOS_PHY_SF_LHFLX  (:,:) = CONST_UNDEF
       ATMOS_PHY_SF_QVFLX  (:,:) = CONST_UNDEF

    endif

    return
  end subroutine ATMOS_PHY_SF_vars_restart_read


  !-----------------------------------------------------------------------------
  !> Write restart
  subroutine ATMOS_PHY_SF_vars_restart_write
    use scale_time, only: &
       TIME_gettimelabel
    use scale_fileio, only: &
       FILEIO_write
    implicit none

    character(len=15)     :: timelabel
    character(len=H_LONG) :: basename

    integer :: n
    !---------------------------------------------------------------------------

    if ( ATMOS_PHY_SF_RESTART_OUT_BASENAME /= '' ) then

       call TIME_gettimelabel( timelabel )
       write(basename,'(A,A,A)') trim(ATMOS_PHY_SF_RESTART_OUT_BASENAME), '_', trim(timelabel)

       if( IO_L ) write(IO_FID_LOG,*)
       if( IO_L ) write(IO_FID_LOG,*) '*** Output restart file (ATMOS_PHY_SF) ***'
       if( IO_L ) write(IO_FID_LOG,*) '*** basename: ', trim(basename)

       call FILEIO_write( ATMOS_PHY_SF_PREC   (:,:), basename,            ATMOS_PHY_SF_RESTART_OUT_TITLE, & ! [IN]
                          VAR_NAME( 1), VAR_DESC( 1), VAR_UNIT( 1), 'XY', ATMOS_PHY_SF_RESTART_OUT_DTYPE  ) ! [IN]
       call FILEIO_write( ATMOS_PHY_SF_SWD    (:,:), basename,            ATMOS_PHY_SF_RESTART_OUT_TITLE, & ! [IN]
                          VAR_NAME( 2), VAR_DESC( 2), VAR_UNIT( 2), 'XY', ATMOS_PHY_SF_RESTART_OUT_DTYPE  ) ! [IN]
       call FILEIO_write( ATMOS_PHY_SF_LWD    (:,:), basename,            ATMOS_PHY_SF_RESTART_OUT_TITLE, & ! [IN]
                          VAR_NAME( 3), VAR_DESC( 3), VAR_UNIT( 3), 'XY', ATMOS_PHY_SF_RESTART_OUT_DTYPE  ) ! [IN]
       call FILEIO_write( ATMOS_PHY_SF_ZMFLX  (:,:), basename,            ATMOS_PHY_SF_RESTART_OUT_TITLE, & ! [IN]
                          VAR_NAME( 4), VAR_DESC( 4), VAR_UNIT( 4), 'XY', ATMOS_PHY_SF_RESTART_OUT_DTYPE  ) ! [IN]
       call FILEIO_write( ATMOS_PHY_SF_XMFLX  (:,:), basename,            ATMOS_PHY_SF_RESTART_OUT_TITLE, & ! [IN]
                          VAR_NAME( 5), VAR_DESC( 5), VAR_UNIT( 5), 'XY', ATMOS_PHY_SF_RESTART_OUT_DTYPE  ) ! [IN]
       call FILEIO_write( ATMOS_PHY_SF_YMFLX  (:,:), basename,            ATMOS_PHY_SF_RESTART_OUT_TITLE, & ! [IN]
                          VAR_NAME( 6), VAR_DESC( 6), VAR_UNIT( 6), 'XY', ATMOS_PHY_SF_RESTART_OUT_DTYPE  ) ! [IN]
       call FILEIO_write( ATMOS_PHY_SF_POTTFLX(:,:), basename,            ATMOS_PHY_SF_RESTART_OUT_TITLE, & ! [IN]
                          VAR_NAME( 7), VAR_DESC( 7), VAR_UNIT( 7), 'XY', ATMOS_PHY_SF_RESTART_OUT_DTYPE  ) ! [IN]
       call FILEIO_write( ATMOS_PHY_SF_SWUFLX (:,:), basename,            ATMOS_PHY_SF_RESTART_OUT_TITLE, & ! [IN]
                          VAR_NAME( 8), VAR_DESC( 8), VAR_UNIT( 8), 'XY', ATMOS_PHY_SF_RESTART_OUT_DTYPE  ) ! [IN]
       call FILEIO_write( ATMOS_PHY_SF_LWUFLX (:,:), basename,            ATMOS_PHY_SF_RESTART_OUT_TITLE, & ! [IN]
                          VAR_NAME( 9), VAR_DESC( 9), VAR_UNIT( 9), 'XY', ATMOS_PHY_SF_RESTART_OUT_DTYPE  ) ! [IN]
       call FILEIO_write( ATMOS_PHY_SF_SHFLX  (:,:), basename,            ATMOS_PHY_SF_RESTART_OUT_TITLE, & ! [IN]
                          VAR_NAME(10), VAR_DESC(10), VAR_UNIT(10), 'XY', ATMOS_PHY_SF_RESTART_OUT_DTYPE  ) ! [IN]
       call FILEIO_write( ATMOS_PHY_SF_LHFLX  (:,:), basename,            ATMOS_PHY_SF_RESTART_OUT_TITLE, & ! [IN]
                          VAR_NAME(11), VAR_DESC(11), VAR_UNIT(11), 'XY', ATMOS_PHY_SF_RESTART_OUT_DTYPE  ) ! [IN]
       call FILEIO_write( ATMOS_PHY_SF_QVFLX  (:,:), basename,            ATMOS_PHY_SF_RESTART_OUT_TITLE, & ! [IN]
                          VAR_NAME(12), VAR_DESC(12), VAR_UNIT(12), 'XY', ATMOS_PHY_SF_RESTART_OUT_DTYPE  ) ! [IN]

    endif

    return
  end subroutine ATMOS_PHY_SF_vars_restart_write

end module mod_atmos_phy_sf_vars
