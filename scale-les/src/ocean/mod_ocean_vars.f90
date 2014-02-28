!-------------------------------------------------------------------------------
!> module OCEAN Variables
!!
!! @par Description
!!          Container for ocean variables
!!
!! @author Team SCALE
!!
!<
!-------------------------------------------------------------------------------
module mod_ocean_vars
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use mod_precision
  use mod_stdio
  use mod_prof
  use mod_debug
  use mod_grid_index
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: OCEAN_vars_setup
  public :: OCEAN_vars_fillhalo
  public :: OCEAN_vars_restart_read
  public :: OCEAN_vars_restart_write
  public :: OCEAN_vars_history
  public :: OCEAN_vars_total

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  real(RP), public, save, allocatable :: TW  (:,:) !< water temperature [K]
  real(RP), public, save, allocatable :: ALBW(:,:) !< surface albedo of short-wave radiation for water [0-1]
  real(RP), public, save, allocatable :: DZW (:,:) !< water depth [m]

  integer,                public, save :: I_TW   = 1
  integer,                public, save :: I_ALBW = 2
  integer,                public, save :: I_DZW  = 3
  character(len=H_SHORT), public, save :: OP_NAME(3) !< name  of the ocean variables
  character(len=H_MID),   public, save :: OP_DESC(3) !< desc. of the ocean variables
  character(len=H_SHORT), public, save :: OP_UNIT(3) !< unit  of the ocean variables

  data OP_NAME / 'TW',   &
                 'ALBW', &
                 'DZW'   /
  data OP_DESC / 'water temperature',                                &
                 'surface albedo of short-wave radiation for water', &
                 'water depth'                                       /
  data OP_UNIT / 'K',   &
                 '0-1', &
                 'm'    /

  character(len=H_SHORT),  public, save :: OCEAN_TYPE_PHY = 'OFF' !< Ocean physics type
  logical,                 public, save :: OCEAN_sw_phy           !< do ocean physics update?
  logical,                 public, save :: OCEAN_sw_restart       !< output restart?

  character(len=H_LONG),   public, save :: OCEAN_RESTART_IN_BASENAME = '' !< basename of the input file

  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  logical,               private, save :: OCEAN_RESTART_OUTPUT       = .false.                   !< output restart file?
  character(len=H_LONG), private, save :: OCEAN_RESTART_OUT_BASENAME = 'restart_out'             !< basename of the output file
  character(len=H_MID),  private, save :: OCEAN_RESTART_OUT_TITLE    = 'SCALE-LES OCEANIC VARS.' !< title    of the output file
  character(len=H_MID),  private, save :: OCEAN_RESTART_OUT_DTYPE    = 'DEFAULT'                 !< REAL4 or REAL8

  logical,               private, save :: OCEAN_VARS_CHECKRANGE      = .false.

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine OCEAN_vars_setup
    use mod_process, only: &
       PRC_MPIstop
    implicit none

    NAMELIST / PARAM_OCEAN / &
       OCEAN_TYPE_PHY

    NAMELIST / PARAM_OCEAN_VARS /  &
       OCEAN_RESTART_IN_BASENAME,  &
       OCEAN_RESTART_OUTPUT,       &
       OCEAN_RESTART_OUT_BASENAME, &
       OCEAN_RESTART_OUT_TITLE,    &
       OCEAN_VARS_CHECKRANGE

    integer :: ierr
    integer :: ip
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[OCEAN VARS]/Categ[OCEAN]'

    allocate( TW  (IA,JA) )
    allocate( ALBW(IA,JA) )
    allocate( DZW (IA,JA) )

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_OCEAN,iostat=ierr)

    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_OCEAN. Check!'
       call PRC_MPIstop
    endif
    if( IO_L ) write(IO_FID_LOG,nml=PARAM_OCEAN)

    if( IO_L ) write(IO_FID_LOG,*) '*** [OCEAN] selected components'

    if ( OCEAN_TYPE_PHY /= 'OFF' .AND. OCEAN_TYPE_PHY /= 'NONE' ) then
       if( IO_L ) write(IO_FID_LOG,*) '*** Ocean physics : ON'
       OCEAN_sw_phy = .true.
    else
       if( IO_L ) write(IO_FID_LOG,*) '*** Ocean physics : OFF'
       OCEAN_sw_phy = .false.
    endif

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[OCEAN VARS]/Categ[OCEAN]'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_OCEAN_VARS,iostat=ierr)

    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_OCEAN_VARS. Check!'
       call PRC_MPIstop
    endif
    if( IO_L ) write(IO_FID_LOG,nml=PARAM_OCEAN_VARS)

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** [OCEAN] prognostic variables'
    if( IO_L ) write(IO_FID_LOG,'(1x,A,A8,A,A32,3(A))') &
               '***       |',' VARNAME','|', 'DESCRIPTION                     ','[', 'UNIT            ',']'
    do ip = 1, 1
       if( IO_L ) write(IO_FID_LOG,'(1x,A,i3,A,A8,A,A32,3(A))') &
                  '*** NO.',ip,'|',trim(OP_NAME(ip)),'|', OP_DESC(ip),'[', OP_UNIT(ip),']'
    enddo

    if( IO_L ) write(IO_FID_LOG,*) 'Output...'
    if ( OCEAN_RESTART_OUTPUT ) then
       if( IO_L ) write(IO_FID_LOG,*) '  Ocean restart output : YES'
       OCEAN_sw_restart = .true.
    else
       if( IO_L ) write(IO_FID_LOG,*) '  Ocean restart output : NO'
       OCEAN_sw_restart = .false.
    endif
    if( IO_L ) write(IO_FID_LOG,*)

    return
  end subroutine OCEAN_vars_setup

  !-----------------------------------------------------------------------------
  !> fill HALO region of land variables
  subroutine OCEAN_vars_fillhalo
    use mod_comm, only: &
       COMM_vars8, &
       COMM_wait
    implicit none
    !---------------------------------------------------------------------------

    ! fill IHALO & JHALO
    call COMM_vars8( TW  (:,:), 1 )
    call COMM_vars8( ALBW(:,:), 2 )
    call COMM_vars8( DZW (:,:), 3 )

    call COMM_wait ( TW  (:,:), 1 )
    call COMM_wait ( ALBW(:,:), 2 )
    call COMM_wait ( DZW (:,:), 3 )

    return
  end subroutine OCEAN_vars_fillhalo

  !-----------------------------------------------------------------------------
  !> Read ocean restart
  subroutine OCEAN_vars_restart_read
    use mod_fileio, only: &
       FILEIO_read
    use mod_comm, only: &
       COMM_vars8, &
       COMM_wait
    implicit none
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** Input restart file (ocean) ***'

    call PROF_rapstart('FILE I NetCDF')

    if ( OCEAN_RESTART_IN_BASENAME /= '' ) then

       call FILEIO_read( TW(:,:),                                        & ! [OUT]
                         OCEAN_RESTART_IN_BASENAME, 'TW',   'XY', step=1 ) ! [IN]
       call FILEIO_read( ALBW(:,:),                                      & ! [OUT]
                         OCEAN_RESTART_IN_BASENAME, 'ALBW', 'XY', step=1 ) ! [IN]
       call FILEIO_read( DZW(:,:),                                       & ! [OUT]
                         OCEAN_RESTART_IN_BASENAME, 'DZW',  'XY', step=1 ) ! [IN]

       ! fill IHALO & JHALO
       call OCEAN_vars_fillhalo
    else
       if( IO_L ) write(IO_FID_LOG,*) '*** restart file for ocean is not specified.'
       TW  (:,:) = 300.0_RP
       ALBW(:,:) = 0.1_RP
       DZW (:,:) = 50.0_RP
    endif

    call PROF_rapend  ('FILE I NetCDF')

    return
  end subroutine OCEAN_vars_restart_read

  !-----------------------------------------------------------------------------
  !> Write ocean restart
  subroutine OCEAN_vars_restart_write
    use mod_time, only: &
       NOWSEC => TIME_NOWDAYSEC
    use mod_fileio, only: &
       FILEIO_write
    implicit none

    character(len=H_LONG) :: bname

    integer :: n
    !---------------------------------------------------------------------------

    call PROF_rapstart('FILE O NetCDF')

    if ( OCEAN_RESTART_OUT_BASENAME /= '' ) then

       if( IO_L ) write(IO_FID_LOG,*)
       if( IO_L ) write(IO_FID_LOG,*) '*** Output restart file (ocean) ***'

       bname = ''
       write(bname(1:15), '(F15.3)') NOWSEC
       do n = 1, 15
          if ( bname(n:n) == ' ' ) bname(n:n) = '0'
       end do
       write(bname,'(A,A,A)') trim(OCEAN_RESTART_OUT_BASENAME), '_', trim(bname)

       call FILEIO_write( TW(:,:),   bname,                         OCEAN_RESTART_OUT_TITLE, & ! [IN]
                          OP_NAME(1), OP_DESC(1), OP_UNIT(1), 'XY', OCEAN_RESTART_OUT_DTYPE  ) ! [IN]
       call FILEIO_write( ALBW(:,:), bname,                         OCEAN_RESTART_OUT_TITLE, & ! [IN]
                          OP_NAME(2), OP_DESC(2), OP_UNIT(2), 'XY', OCEAN_RESTART_OUT_DTYPE  ) ! [IN]
       call FILEIO_write( DZW(:,:),  bname,                         OCEAN_RESTART_OUT_TITLE, & ! [IN]
                          OP_NAME(3), OP_DESC(3), OP_UNIT(3), 'XY', OCEAN_RESTART_OUT_DTYPE  ) ! [IN]

    endif

    call PROF_rapend  ('FILE O NetCDF')

    return
  end subroutine OCEAN_vars_restart_write

  !-----------------------------------------------------------------------------
  !> History output set for ocean variables
  subroutine OCEAN_vars_history
    use mod_time, only: &
       TIME_DTSEC_OCEAN
    use mod_history, only: &
       HIST_in
    implicit none
    !---------------------------------------------------------------------------

    if ( OCEAN_VARS_CHECKRANGE ) then
       call VALCHECK( TW  (:,:), 0.0_RP, 1000.0_RP, OP_NAME(I_TW),   __FILE__, __LINE__ )
       call VALCHECK( ALBW(:,:), 0.0_RP,    2.0_RP, OP_NAME(I_ALBW), __FILE__, __LINE__ )
       call VALCHECK( DZW (:,:), 0.0_RP, 1.0E+6_RP, OP_NAME(I_DZW),  __FILE__, __LINE__ )
    endif

    call HIST_in( TW(:,:),   'TW',   OP_DESC(I_TW),   OP_UNIT(I_TW),   TIME_DTSEC_OCEAN )
    call HIST_in( ALBW(:,:), 'ALBW', OP_DESC(I_ALBW), OP_UNIT(I_ALBW), TIME_DTSEC_OCEAN )
    call HIST_in( DZW(:,:),  'DZW',  OP_DESC(I_DZW),  OP_UNIT(I_DZW),  TIME_DTSEC_OCEAN )

    return
  end subroutine OCEAN_vars_history

  !-----------------------------------------------------------------------------
  !> Budget monitor for ocean
  subroutine OCEAN_vars_total
    use mod_stats, only: &
       STAT_checktotal, &
       STAT_total
    implicit none

    real(RP) :: total
    !---------------------------------------------------------------------------

    if ( STAT_checktotal ) then

!       call STAT_total( total, TW(:,:),   OP_NAME(I_TW)   )
!       call STAT_total( total, ALBW(:,:), OP_NAME(I_ALBW) )
!       call STAT_total( total, DZW(:,:),  OP_NAME(I_DZW)  )

    endif

    return
  end subroutine OCEAN_vars_total

end module mod_ocean_vars
