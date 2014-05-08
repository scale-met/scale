!-------------------------------------------------------------------------------
!> module URBAN Variables
!!
!! @par Description
!!          Container for urban variables
!!
!! @author Team SCALE
!!
!<
!-------------------------------------------------------------------------------
module mod_urban_vars
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
  !++ included parameters
  !
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: URBAN_vars_setup
  public :: URBAN_vars_fillhalo
  public :: URBAN_vars_restart_read
  public :: URBAN_vars_restart_write
  public :: URBAN_vars_history
  public :: URBAN_vars_total

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  character(len=H_SHORT), public, save :: URBAN_TYPE_PHY = 'OFF' !< urban physics type
  logical,                public, save :: URBAN_sw_phy           !< do urban physics update?
  logical,                public, save :: URBAN_sw_restart       !< output restart?

  ! prognostic variables

  integer,  public, parameter :: PV_NUM = 1

  character(len=H_SHORT), public, save :: PV_NAME(PV_NUM) !< name  of the urban variables
  character(len=H_MID),   public, save :: PV_DESC(PV_NUM) !< desc. of the urban variables
  character(len=H_SHORT), public, save :: PV_UNIT(PV_NUM) !< unit  of the urban variables

  data PV_NAME / 'dummy' /
  data PV_DESC / 'dummy' /
  data PV_UNIT / 'dummy' /

  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  logical,               private, save :: URBAN_RESTART_OUTPUT       = .false.                 !< output restart file?
  character(len=H_LONG), private, save :: URBAN_RESTART_IN_BASENAME  = ''                      !< basename of the restart file
  character(len=H_LONG), private, save :: URBAN_RESTART_OUT_BASENAME = ''                      !< basename of the output file
  character(len=H_MID),  private, save :: URBAN_RESTART_OUT_TITLE    = 'SCALE-LES URBAN VARS.' !< title    of the output file
  character(len=H_MID),  private, save :: URBAN_RESTART_OUT_DTYPE    = 'DEFAULT'               !< REAL4 or REAL8

  logical,               private, save :: URBAN_VARS_CHECKRANGE      = .false.

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine URBAN_vars_setup
    use scale_process, only: &
       PRC_MPIstop
    implicit none

    NAMELIST / PARAM_URBAN / &
       URBAN_TYPE_PHY

    NAMELIST / PARAM_URBAN_VARS /  &
       URBAN_RESTART_IN_BASENAME,  &
       URBAN_RESTART_OUTPUT,       &
       URBAN_RESTART_OUT_BASENAME, &
       URBAN_RESTART_OUT_TITLE,    &
       URBAN_RESTART_OUT_DTYPE,    &
       URBAN_VARS_CHECKRANGE

    integer :: ierr
    integer :: ip
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[URBAN VARS]/Categ[URBAN]'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_URBAN,iostat=ierr)

    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_URBAN. Check!'
       call PRC_MPIstop
    endif
    if( IO_L ) write(IO_FID_LOG,nml=PARAM_URBAN)

    if( IO_L ) write(IO_FID_LOG,*) '*** [URBAN] selected components'

    if ( URBAN_TYPE_PHY /= 'OFF' .AND. URBAN_TYPE_PHY /= 'NONE' ) then
       if( IO_L ) write(IO_FID_LOG,*) '*** Urban physics : ON'
       URBAN_sw_phy = .true.
    else
       if( IO_L ) write(IO_FID_LOG,*) '*** Urban physics : OFF'
       URBAN_sw_phy = .false.
    endif

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_URBAN_VARS,iostat=ierr)

    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_URBAN_VARS. Check!'
       call PRC_MPIstop
    endif
    if( IO_L ) write(IO_FID_LOG,nml=PARAM_URBAN_VARS)

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** [URBAN] prognostic variables'
    if( IO_L ) write(IO_FID_LOG,'(1x,A,A8,A,A32,3(A))') &
               '***       |',' VARNAME','|', 'DESCRIPTION                     ','[', 'UNIT            ',']'
    do ip = 1, PV_NUM
       if( IO_L ) write(IO_FID_LOG,'(1x,A,i3,A,A8,A,A32,3(A))') &
                  '*** NO.',ip,'|',trim(PV_NAME(ip)),'|', PV_DESC(ip),'[', PV_UNIT(ip),']'
    enddo

    ! restart switch
    if( IO_L ) write(IO_FID_LOG,*) 'Output...'
    if ( URBAN_RESTART_OUTPUT ) then
       if( IO_L ) write(IO_FID_LOG,*) '  Urban restart output : YES'
       URBAN_sw_restart = .true.
    else
       if( IO_L ) write(IO_FID_LOG,*) '  Urban restart output : NO'
       URBAN_sw_restart = .false.
    endif
    if( IO_L ) write(IO_FID_LOG,*)

    return
  end subroutine URBAN_vars_setup

  !-----------------------------------------------------------------------------
  !> fill HALO region of urban variables
  subroutine URBAN_vars_fillhalo
    use scale_comm, only: &
       COMM_vars8, &
       COMM_wait
    implicit none

    integer :: k
    !---------------------------------------------------------------------------

    ! fill IHALO & JHALO

    ! 2D variable

    return
  end subroutine URBAN_vars_fillhalo

  !-----------------------------------------------------------------------------
  !> Read urban restart
  subroutine URBAN_vars_restart_read
    use scale_fileio, only: &
       FILEIO_read
    use scale_const, only: &
       CONST_UNDEF
    use scale_comm, only: &
       COMM_vars8, &
       COMM_wait
    implicit none

    integer :: i, j, v
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** Input restart file (urban) ***'

    call PROF_rapstart('FILE I NetCDF')

    if ( URBAN_RESTART_IN_BASENAME /= '' ) then

       call URBAN_vars_fillhalo

       call URBAN_vars_total
    else
       if( IO_L ) write(IO_FID_LOG,*) '*** restart file for urban is not specified.'

    endif

    call PROF_rapend  ('FILE I NetCDF')

    return
  end subroutine URBAN_vars_restart_read

  !-----------------------------------------------------------------------------
  !> Write urban restart
  subroutine URBAN_vars_restart_write
    use scale_time, only: &
       NOWSEC => TIME_NOWDAYSEC
    use scale_fileio, only: &
       FILEIO_write
    implicit none

    character(len=H_LONG) :: basename

    integer :: n
    !---------------------------------------------------------------------------

    call PROF_rapstart('FILE O NetCDF')

    if ( URBAN_RESTART_OUT_BASENAME /= '' ) then

       if( IO_L ) write(IO_FID_LOG,*)
       if( IO_L ) write(IO_FID_LOG,*) '*** Output restart file (urban) ***'

       basename = ''
       write(basename(1:15), '(F15.3)') NOWSEC
       do n = 1, 15
          if ( basename(n:n) == ' ' ) basename(n:n) = '0'
       enddo
       basename = trim(URBAN_RESTART_OUT_BASENAME) // '_' // trim(basename)

    endif

    call PROF_rapend  ('FILE O NetCDF')

    call URBAN_vars_total

    return
  end subroutine URBAN_vars_restart_write

  !-----------------------------------------------------------------------------
  !> History output set for urban variables
  subroutine URBAN_vars_history
    use scale_time, only: &
       TIME_DTSEC_URBAN
    use scale_history, only: &
       HIST_in
    implicit none
    !---------------------------------------------------------------------------

    if ( URBAN_VARS_CHECKRANGE ) then
    endif

    return
  end subroutine URBAN_vars_history

  !-----------------------------------------------------------------------------
  !> Budget monitor for urban
  subroutine URBAN_vars_total
    use scale_stats, only: &
       STAT_checktotal, &
       STAT_total
    implicit none

    !real(RP) :: total
    !---------------------------------------------------------------------------

    if ( STAT_checktotal ) then

    endif

    return
  end subroutine URBAN_vars_total

end module mod_urban_vars
