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
  use scale_precision
  use scale_stdio
  use scale_prof
  use scale_debug
  use scale_grid_index
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
  real(RP), public, save, allocatable :: TW(:,:) !< water temperature [K]

  integer, public, parameter :: PV_NUM = 1

  integer, public, parameter :: I_TW   = 1

  character(len=H_SHORT), public, save :: PV_NAME(PV_NUM) !< name  of the ocean variables
  character(len=H_MID),   public, save :: PV_DESC(PV_NUM) !< desc. of the ocean variables
  character(len=H_SHORT), public, save :: PV_UNIT(PV_NUM) !< unit  of the ocean variables

  data PV_NAME / 'TW' /
  data PV_DESC / 'water temperature' /
  data PV_UNIT / 'K' /

  character(len=H_SHORT),  public, save :: OCEAN_TYPE_PHY = 'OFF' !< Ocean physics type
  logical,                 public, save :: OCEAN_sw_phy           !< do ocean physics update?
  logical,                 public, save :: OCEAN_sw_restart       !< output restart?

  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  logical,               private, save :: OCEAN_RESTART_OUTPUT       = .false.                   !< output restart file?
  character(len=H_LONG), private, save :: OCEAN_RESTART_IN_BASENAME  = ''                        !< basename of the input file
  character(len=H_LONG), private, save :: OCEAN_RESTART_OUT_BASENAME = ''                        !< basename of the output file
  character(len=H_MID),  private, save :: OCEAN_RESTART_OUT_TITLE    = 'SCALE-LES OCEANIC VARS.' !< title    of the output file
  character(len=H_MID),  private, save :: OCEAN_RESTART_OUT_DTYPE    = 'DEFAULT'                 !< REAL4 or REAL8

  logical,               private, save :: OCEAN_VARS_CHECKRANGE      = .false.

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine OCEAN_vars_setup
    use scale_process, only: &
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

    allocate( TW(IA,JA) )

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
    do ip = 1, PV_NUM
       if( IO_L ) write(IO_FID_LOG,'(1x,A,i3,A,A8,A,A32,3(A))') &
                  '*** NO.',ip,'|',trim(PV_NAME(ip)),'|', PV_DESC(ip),'[', PV_UNIT(ip),']'
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
    use scale_comm, only: &
       COMM_vars8, &
       COMM_wait
    implicit none
    !---------------------------------------------------------------------------

    ! fill IHALO & JHALO
    call COMM_vars8( TW(:,:), 1 )

    call COMM_wait ( TW(:,:), 1 )

    return
  end subroutine OCEAN_vars_fillhalo

  !-----------------------------------------------------------------------------
  !> Read ocean restart
  subroutine OCEAN_vars_restart_read
    use scale_fileio, only: &
       FILEIO_read
    use scale_comm, only: &
       COMM_vars8, &
       COMM_wait
    implicit none
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** Input restart file (ocean) ***'

    call PROF_rapstart('FILE I NetCDF')

    if ( OCEAN_RESTART_IN_BASENAME /= '' ) then

       call FILEIO_read( TW(:,:),                                      & ! [OUT]
                         OCEAN_RESTART_IN_BASENAME, 'TW', 'XY', step=1 ) ! [IN]

       ! fill IHALO & JHALO
       call OCEAN_vars_fillhalo
    else
       if( IO_L ) write(IO_FID_LOG,*) '*** restart file for ocean is not specified.'
       TW(:,:) = 300.0_RP
    endif

    call PROF_rapend  ('FILE I NetCDF')

    return
  end subroutine OCEAN_vars_restart_read

  !-----------------------------------------------------------------------------
  !> Write ocean restart
  subroutine OCEAN_vars_restart_write
    use scale_time, only: &
       NOWSEC => TIME_NOWDAYSEC
    use scale_fileio, only: &
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

       call FILEIO_write( TW(:,:),   bname,                         OCEAN_RESTART_OUT_TITLE,        & ! [IN]
                          PV_NAME(1), PV_DESC(1), PV_UNIT(1), 'XY', OCEAN_RESTART_OUT_DTYPE, .true. ) ! [IN]

    endif

    call PROF_rapend  ('FILE O NetCDF')

    return
  end subroutine OCEAN_vars_restart_write

  !-----------------------------------------------------------------------------
  !> History output set for ocean variables
  subroutine OCEAN_vars_history
    use scale_time, only: &
       TIME_DTSEC_OCEAN
    use scale_history, only: &
       HIST_in
    implicit none
    !---------------------------------------------------------------------------

    if ( OCEAN_VARS_CHECKRANGE ) then
       call VALCHECK( TW(:,:), 0.0_RP, 1000.0_RP, PV_NAME(I_TW), __FILE__, __LINE__ )
    endif

    call HIST_in( TW(:,:), 'TW', PV_DESC(I_TW), PV_UNIT(I_TW), TIME_DTSEC_OCEAN )

    return
  end subroutine OCEAN_vars_history

  !-----------------------------------------------------------------------------
  !> Budget monitor for ocean
  subroutine OCEAN_vars_total
    use scale_stats, only: &
       STAT_checktotal, &
       STAT_total
    implicit none

    real(RP) :: total
    !---------------------------------------------------------------------------

    if ( STAT_checktotal ) then

!       call STAT_total( total, TW(:,:), PV_NAME(I_TW) )

    endif

    return
  end subroutine OCEAN_vars_total

end module mod_ocean_vars
