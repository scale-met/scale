!-------------------------------------------------------------------------------
!> module OCEAN Variables
!!
!! @par Description
!!          Container for ocean variables
!!
!! @author Team SCALE
!! @li      2013-08-31 (T.Yamaura)  [new]
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
  public :: OCEAN_vars_external_in

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  logical, public :: OCEAN_sw_restart

  ! prognostic variables
  real(RP), public, allocatable :: TW(:,:) !< water temperature [K]

  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  logical,                private :: OCEAN_RESTART_OUTPUT       = .false.         !< output restart file?
  character(len=H_LONG),  private :: OCEAN_RESTART_IN_BASENAME  = ''              !< basename of the input file
  character(len=H_LONG),  private :: OCEAN_RESTART_OUT_BASENAME = ''              !< basename of the output file
  character(len=H_MID),   private :: OCEAN_RESTART_OUT_TITLE    = 'OCEAN restart' !< title    of the output file
  character(len=H_MID),   private :: OCEAN_RESTART_OUT_DTYPE    = 'DEFAULT'       !< REAL4 or REAL8

  logical,                private :: OCEAN_VARS_CHECKRANGE      = .false.

  integer,                private, parameter :: VMAX = 1       !< number of the variables
  integer,                private, parameter :: I_TW = 1

  character(len=H_SHORT), private            :: VAR_NAME(VMAX) !< name  of the variables
  character(len=H_MID),   private            :: VAR_DESC(VMAX) !< desc. of the variables
  character(len=H_SHORT), private            :: VAR_UNIT(VMAX) !< unit  of the variables

  data VAR_NAME / 'TW' /
  data VAR_DESC / 'water temperature' /
  data VAR_UNIT / 'K' /

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine OCEAN_vars_setup
    use scale_process, only: &
       PRC_MPIstop
    use scale_const, only: &
       UNDEF => CONST_UNDEF
    implicit none

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
    if( IO_L ) write(IO_FID_LOG,*) '++++++ Module[VARS] / Categ[OCEAN] / Origin[SCALE-LES]'

    allocate( TW(IA,JA) )
    TW(:,:) = UNDEF

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_OCEAN_VARS,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_OCEAN_VARS. Check!'
       call PRC_MPIstop
    endif
    if( IO_LNML ) write(IO_FID_LOG,nml=PARAM_OCEAN_VARS)

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** List of prognostic variables (OCEAN) ***'
    if( IO_L ) write(IO_FID_LOG,'(1x,A,A16,A,A32,3(A))') &
               '***       |','         VARNAME','|', 'DESCRIPTION                     ','[', 'UNIT            ',']'
    do ip = 1, VMAX
       if( IO_L ) write(IO_FID_LOG,'(1x,A,i3,A,A16,A,A32,3(A))') &
                  '*** NO.',ip,'|',trim(VAR_NAME(ip)),'|', VAR_DESC(ip),'[', VAR_UNIT(ip),']'
    enddo

    if( IO_L ) write(IO_FID_LOG,*)
    if ( OCEAN_RESTART_IN_BASENAME /= '' ) then
       if( IO_L ) write(IO_FID_LOG,*) '*** Restart input?  : ', trim(OCEAN_RESTART_IN_BASENAME)
    else
       if( IO_L ) write(IO_FID_LOG,*) '*** Restart input?  : NO'
    endif
    if (       OCEAN_RESTART_OUTPUT             &
         .AND. OCEAN_RESTART_OUT_BASENAME /= '' ) then
       if( IO_L ) write(IO_FID_LOG,*) '*** Restart output? : ', trim(OCEAN_RESTART_OUT_BASENAME)
       OCEAN_sw_restart = .true.
    else
       if( IO_L ) write(IO_FID_LOG,*) '*** Restart output? : NO'
       OCEAN_RESTART_OUTPUT = .false.
       OCEAN_sw_restart = .false.
    endif

    return
  end subroutine OCEAN_vars_setup

  !-----------------------------------------------------------------------------
  !> HALO Communication
  subroutine OCEAN_vars_fillhalo
    use scale_comm, only: &
       COMM_vars8, &
       COMM_wait
    implicit none
    !---------------------------------------------------------------------------

    call COMM_vars8( TW(:,:), 1 )
    call COMM_wait ( TW(:,:), 1 )

    return
  end subroutine OCEAN_vars_fillhalo

  !-----------------------------------------------------------------------------
  !> Read ocean restart
  subroutine OCEAN_vars_restart_read
    use scale_fileio, only: &
       FILEIO_read
    implicit none
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** Input restart file (OCEAN) ***'

    if ( OCEAN_RESTART_IN_BASENAME /= '' ) then
       if( IO_L ) write(IO_FID_LOG,*) '*** basename: ', trim(OCEAN_RESTART_IN_BASENAME)

       call FILEIO_read( TW(:,:),                                      & ! [OUT]
                         OCEAN_RESTART_IN_BASENAME, 'TW', 'XY', step=1 ) ! [IN]

       call OCEAN_vars_fillhalo

       call OCEAN_vars_total
    else
       if( IO_L ) write(IO_FID_LOG,*) '*** restart file for ocean is not specified.'
    endif

    return
  end subroutine OCEAN_vars_restart_read

  !-----------------------------------------------------------------------------
  !> Write ocean restart
  subroutine OCEAN_vars_restart_write
    use scale_time, only: &
       TIME_gettimelabel
    use scale_fileio, only: &
       FILEIO_write
    implicit none

    character(len=15)     :: timelabel
    character(len=H_LONG) :: basename
    !---------------------------------------------------------------------------

    if ( OCEAN_RESTART_OUT_BASENAME /= '' ) then

       call TIME_gettimelabel( timelabel )
       write(basename,'(A,A,A)') trim(OCEAN_RESTART_OUT_BASENAME), '_', trim(timelabel)

       if( IO_L ) write(IO_FID_LOG,*)
       if( IO_L ) write(IO_FID_LOG,*) '*** Output restart file (OCEAN) ***'
       if( IO_L ) write(IO_FID_LOG,*) '*** basename: ', trim(basename)

       call OCEAN_vars_total

       call FILEIO_write( TW(:,:), basename,                                    OCEAN_RESTART_OUT_TITLE, & ! [IN]
                          VAR_NAME(I_TW), VAR_DESC(I_TW), VAR_UNIT(I_TW), 'XY', OCEAN_RESTART_OUT_DTYPE  ) ! [IN]

    endif

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
       call VALCHECK( TW(:,:), 0.0_RP, 1000.0_RP, VAR_NAME(I_TW), __FILE__, __LINE__ )
    endif

    call HIST_in( TW(:,:), 'TW', VAR_DESC(I_TW), VAR_UNIT(I_TW), TIME_DTSEC_OCEAN )

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

       call STAT_total( total, TW(:,:), VAR_NAME(I_TW) )

    endif

    return
  end subroutine OCEAN_vars_total

  !-----------------------------------------------------------------------------
  !> Input from External I/O
  subroutine OCEAN_vars_external_in( &
      tw_in      ) ! (in)
    implicit none

    real(RP), intent(in) :: tw_in(:,:)
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** External Input file (ocean) ***'

    TW(:,:)   = tw_in(:,:)

    call OCEAN_vars_fillhalo

    call OCEAN_vars_total

    return
  end subroutine OCEAN_vars_external_in

end module mod_ocean_vars
