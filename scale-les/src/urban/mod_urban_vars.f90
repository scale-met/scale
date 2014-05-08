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

  integer ,               public, save :: num_urb_layers = 4     !< number of urban layers

  ! prognostic variables
  real(RP), public, save, allocatable :: TR_URB (:,:)   ! Surface temperature of roof [K]
  real(RP), public, save, allocatable :: TB_URB (:,:)   ! Surface temperature of building [K]
  real(RP), public, save, allocatable :: TG_URB (:,:)   ! Surface temperature of ground [K]
  real(RP), public, save, allocatable :: TC_URB (:,:)   ! Diagnostic canopy air temperature [K]
  real(RP), public, save, allocatable :: QC_URB (:,:)   ! Diagnostic canopy humidity [-]
  real(RP), public, save, allocatable :: UC_URB (:,:)   ! Diagnostic canopy wind [m/s]
  real(RP), public, save, allocatable :: TRL_URB(:,:,:) ! temperature in layer of roof [K]
  real(RP), public, save, allocatable :: TBL_URB(:,:,:) ! temperature in layer of building [K]
  real(RP), public, save, allocatable :: TGL_URB(:,:,:) ! temperature in layer of ground [K]

  integer,  public, parameter :: PV_NUM = 9

  integer,  public, parameter :: I_TR_URB  = 1
  integer,  public, parameter :: I_TB_URB  = 2
  integer,  public, parameter :: I_TG_URB  = 3
  integer,  public, parameter :: I_TC_URB  = 4
  integer,  public, parameter :: I_QC_URB  = 5
  integer,  public, parameter :: I_UC_URB  = 6
  integer,  public, parameter :: I_TRL_URB = 7
  integer,  public, parameter :: I_TBL_URB = 8
  integer,  public, parameter :: I_TGL_URB = 9

  character(len=H_SHORT), public, save :: PV_NAME(PV_NUM) !< name  of the urban variables
  character(len=H_MID),   public, save :: PV_DESC(PV_NUM) !< desc. of the urban variables
  character(len=H_SHORT), public, save :: PV_UNIT(PV_NUM) !< unit  of the urban variables

  data PV_NAME / 'TR_URB' , &
                 'TB_URB' , &
                 'TG_URB' , &
                 'TC_URB' , &
                 'QC_URB' , &
                 'UC_URB' , &
                 'TRL_URB', &
                 'TBL_URB', &
                 'TGL_URB'  /

  data PV_DESC / 'Surface temperature of roof',       &
                 'Surface temperature of building',   &
                 'Surface temperature of ground',     &
                 'Diagnostic canopy air temperature', &
                 'Diagnostic canopy humidity',        &
                 'Diagnostic canopy wind',            &
                 'temperature in layer of roof',      &
                 'temperature in layer of building',  &
                 'temperature in layer of ground'     /

  data PV_UNIT / 'K',   &
                 'K',   &
                 'K',   &
                 'K',   &
                 '-',   &
                 'm/s', &
                 'K',   &
                 'K',   &
                 'K'    /

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
    use scale_const, only: &
       UNDEF => CONST_UNDEF
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

    ! allocate arrays
    allocate( TR_URB (IA,JA) )
    allocate( TB_URB (IA,JA) )
    allocate( TG_URB (IA,JA) )
    allocate( TC_URB (IA,JA) )
    allocate( QC_URB (IA,JA) )
    allocate( UC_URB (IA,JA) )

    allocate( TRL_URB(IA,JA,num_urb_layers) )
    allocate( TBL_URB(IA,JA,num_urb_layers) )
    allocate( TGL_URB(IA,JA,num_urb_layers) )

    ! initialize
    TR_URB (:,:) = UNDEF
    TB_URB (:,:) = UNDEF
    TG_URB (:,:) = UNDEF
    TC_URB (:,:) = UNDEF
    QC_URB (:,:) = UNDEF
    UC_URB (:,:) = UNDEF

    TRL_URB(:,:,:) = UNDEF
    TBL_URB(:,:,:) = UNDEF
    TGL_URB(:,:,:) = UNDEF

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
    do k = 1, num_urb_layers
      call COMM_vars8( TRL_URB(:,:,k), 1 )
      call COMM_vars8( TBL_URB(:,:,k), 2 )
      call COMM_vars8( TGL_URB(:,:,k), 3 )

      call COMM_wait ( TRL_URB(:,:,k), 1 )
      call COMM_wait ( TBL_URB(:,:,k), 2 )
      call COMM_wait ( TGL_URB(:,:,k), 3 )
    end do

    ! 2D variable
    call COMM_vars8( TR_URB(:,:), 4 )
    call COMM_vars8( TB_URB(:,:), 5 )
    call COMM_vars8( TG_URB(:,:), 6 )
    call COMM_vars8( TC_URB(:,:), 7 )
    call COMM_vars8( QC_URB(:,:), 8 )
    call COMM_vars8( UC_URB(:,:), 9 )

    call COMM_wait ( TR_URB(:,:), 4 )
    call COMM_wait ( TB_URB(:,:), 5 )
    call COMM_wait ( TG_URB(:,:), 6 )
    call COMM_wait ( TC_URB(:,:), 7 )
    call COMM_wait ( QC_URB(:,:), 8 )
    call COMM_wait ( UC_URB(:,:), 9 )

    return
  end subroutine URBAN_vars_fillhalo

  !-----------------------------------------------------------------------------
  !> Read urban restart
  subroutine URBAN_vars_restart_read
    use scale_fileio, only: &
       FILEIO_read
    use scale_const, only: &
       UNDEF => CONST_UNDEF
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

!       call FILEIO_read( TR_URB(:,:),                                     & ! [OUT]
!                         LAND_RESTART_IN_BASENAME, 'TR_URB', 'XY', step=1 ) ! [IN]
!       call FILEIO_read( TB_URB(:,:),                                     & ! [OUT]
!                         LAND_RESTART_IN_BASENAME, 'TB_URB', 'XY', step=1 ) ! [IN]
!       call FILEIO_read( TG_URB(:,:),                                     & ! [OUT]
!                         LAND_RESTART_IN_BASENAME, 'TG_URB', 'XY', step=1 ) ! [IN]
!       call FILEIO_read( TC_URB(:,:),                                     & ! [OUT]
!                         LAND_RESTART_IN_BASENAME, 'TC_URB', 'XY', step=1 ) ! [IN]
!       call FILEIO_read( QC_URB(:,:),                                     & ! [OUT]
!                         LAND_RESTART_IN_BASENAME, 'QC_URB', 'XY', step=1 ) ! [IN]
!       call FILEIO_read( UC_URB(:,:),                                     & ! [OUT]
!                         LAND_RESTART_IN_BASENAME, 'UC_URB', 'XY', step=1 ) ! [IN]
!
!       call FILEIO_read( TRL_URB(:,:),                                      & ! [OUT]
!                         LAND_RESTART_IN_BASENAME, 'TRL_URB', 'XYZ', step=1 ) ! [IN]
!       call FILEIO_read( TBL_URB(:,:),                                      & ! [OUT]
!                         LAND_RESTART_IN_BASENAME, 'TBL_URB', 'XYZ', step=1 ) ! [IN]
!       call FILEIO_read( TGL_URB(:,:),                                      & ! [OUT]
!                         LAND_RESTART_IN_BASENAME, 'TGL_URB', 'XYZ', step=1 ) ! [IN]

       call URBAN_vars_fillhalo

       call URBAN_vars_total
    else
       if( IO_L ) write(IO_FID_LOG,*) '*** restart file for urban is not specified.'

       TR_URB (:,:) = UNDEF
       TB_URB (:,:) = UNDEF
       TG_URB (:,:) = UNDEF
       TC_URB (:,:) = UNDEF
       QC_URB (:,:) = UNDEF
       UC_URB (:,:) = UNDEF

       TRL_URB(:,:,:) = UNDEF
       TBL_URB(:,:,:) = UNDEF
       TGL_URB(:,:,:) = UNDEF

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

!       call FILEIO_write( TR_URB(:,:),  basename,                                        URBAN_RESTART_OUT_TITLE,        & ! [IN]
!                          PV_NAME(I_TR_URB), PV_DESC(I_TR_URB), PV_UNIT(I_TR_URB), 'XY', URBAN_RESTART_OUT_DTYPE, .true. ) ! [IN]
!       call FILEIO_write( TB_URB(:,:),  basename,                                        URBAN_RESTART_OUT_TITLE,        & ! [IN]
!                          PV_NAME(I_TB_URB), PV_DESC(I_TB_URB), PV_UNIT(I_TB_URB), 'XY', URBAN_RESTART_OUT_DTYPE, .true. ) ! [IN]
!       call FILEIO_write( TG_URB(:,:),  basename,                                        URBAN_RESTART_OUT_TITLE,        & ! [IN]
!                          PV_NAME(I_TG_URB), PV_DESC(I_TG_URB), PV_UNIT(I_TG_URB), 'XY', URBAN_RESTART_OUT_DTYPE, .true. ) ! [IN]
!       call FILEIO_write( TC_URB(:,:),  basename,                                        URBAN_RESTART_OUT_TITLE,        & ! [IN]
!                          PV_NAME(I_TC_URB), PV_DESC(I_TC_URB), PV_UNIT(I_TC_URB), 'XY', URBAN_RESTART_OUT_DTYPE, .true. ) ! [IN]
!       call FILEIO_write( QC_URB(:,:),  basename,                                        URBAN_RESTART_OUT_TITLE,        & ! [IN]
!                          PV_NAME(I_QC_URB), PV_DESC(I_QC_URB), PV_UNIT(I_QC_URB), 'XY', URBAN_RESTART_OUT_DTYPE, .true. ) ! [IN]
!       call FILEIO_write( UC_URB(:,:),  basename,                                        URBAN_RESTART_OUT_TITLE,        & ! [IN]
!                          PV_NAME(I_UC_URB), PV_DESC(I_UC_URB), PV_UNIT(I_UC_URB), 'XY', URBAN_RESTART_OUT_DTYPE, .true. ) ! [IN]
!
!       call FILEIO_write( TRL_URB(:,:,:), basename,                                          URBAN_RESTART_OUT_TITLE,        & ! [IN]
!                          PV_NAME(I_TRL_URB), PV_DESC(I_TRL_URB), PV_UNIT(I_TRL_URB), 'XYZ', URBAN_RESTART_OUT_DTYPE, .true. ) ! [IN]
!       call FILEIO_write( TBL_URB(:,:,:), basename,                                          URBAN_RESTART_OUT_TITLE,        & ! [IN]
!                          PV_NAME(I_TBL_URB), PV_DESC(I_TBL_URB), PV_UNIT(I_TBL_URB), 'XYZ', URBAN_RESTART_OUT_DTYPE, .true. ) ! [IN]
!       call FILEIO_write( TGL_URB(:,:,:), basename,                                          URBAN_RESTART_OUT_TITLE,        & ! [IN]
!                          PV_NAME(I_TGL_URB), PV_DESC(I_TGL_URB), PV_UNIT(I_TGL_URB), 'XYZ', URBAN_RESTART_OUT_DTYPE, .true. ) ! [IN]

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
       call VALCHECK( TR_URB(:,:), 0.0_RP, 1000.0_RP, PV_NAME(I_TR_URB), __FILE__, __LINE__ )
       call VALCHECK( TB_URB(:,:), 0.0_RP, 1000.0_RP, PV_NAME(I_TB_URB), __FILE__, __LINE__ )
       call VALCHECK( TG_URB(:,:), 0.0_RP, 1000.0_RP, PV_NAME(I_TG_URB), __FILE__, __LINE__ )
       call VALCHECK( TC_URB(:,:), 0.0_RP, 1000.0_RP, PV_NAME(I_TC_URB), __FILE__, __LINE__ )
       call VALCHECK( QC_URB(:,:), 0.0_RP, 1000.0_RP, PV_NAME(I_QC_URB), __FILE__, __LINE__ )
       call VALCHECK( UC_URB(:,:), 0.0_RP, 1000.0_RP, PV_NAME(I_UC_URB), __FILE__, __LINE__ )

       call VALCHECK( TRL_URB(:,:,:), 0.0_RP, 1000.0_RP, PV_NAME(I_TRL_URB), __FILE__, __LINE__ )
       call VALCHECK( TBL_URB(:,:,:), 0.0_RP, 1000.0_RP, PV_NAME(I_TBL_URB), __FILE__, __LINE__ )
       call VALCHECK( TGL_URB(:,:,:), 0.0_RP, 1000.0_RP, PV_NAME(I_TGL_URB), __FILE__, __LINE__ )
    endif

    call HIST_in( TR_URB(:,:), 'TR_URB', PV_DESC(I_TR_URB), PV_UNIT(I_TR_URB), TIME_DTSEC_URBAN )
    call HIST_in( TB_URB(:,:), 'TB_URB', PV_DESC(I_TB_URB), PV_UNIT(I_TB_URB), TIME_DTSEC_URBAN )
    call HIST_in( TG_URB(:,:), 'TG_URB', PV_DESC(I_TG_URB), PV_UNIT(I_TG_URB), TIME_DTSEC_URBAN )
    call HIST_in( TC_URB(:,:), 'TC_URB', PV_DESC(I_TC_URB), PV_UNIT(I_TC_URB), TIME_DTSEC_URBAN )
    call HIST_in( QC_URB(:,:), 'QC_URB', PV_DESC(I_QC_URB), PV_UNIT(I_QC_URB), TIME_DTSEC_URBAN )
    call HIST_in( UC_URB(:,:), 'UC_URB', PV_DESC(I_UC_URB), PV_UNIT(I_UC_URB), TIME_DTSEC_URBAN )

    call HIST_in( TRL_URB(:,:,:), 'TRL_URB', PV_DESC(I_TRL_URB), PV_UNIT(I_TRL_URB), TIME_DTSEC_URBAN )
    call HIST_in( TBL_URB(:,:,:), 'TBL_URB', PV_DESC(I_TBL_URB), PV_UNIT(I_TBL_URB), TIME_DTSEC_URBAN )
    call HIST_in( TGL_URB(:,:,:), 'TGL_URB', PV_DESC(I_TGL_URB), PV_UNIT(I_TGL_URB), TIME_DTSEC_URBAN )

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
!       call STAT_total( total, TR_URB(:,:), PV_NAME(I_TR_URB) )
!       call STAT_total( total, TB_URB(:,:), PV_NAME(I_TB_URB) )
!       call STAT_total( total, TG_URB(:,:), PV_NAME(I_TG_URB) )
!       call STAT_total( total, TC_URB(:,:), PV_NAME(I_TC_URB) )
!       call STAT_total( total, QC_URB(:,:), PV_NAME(I_QC_URB) )
!       call STAT_total( total, UC_URB(:,:), PV_NAME(I_UC_URB) )

!       call STAT_total( total, TRL_URB(:,:), PV_NAME(I_TRL_URB) )
!       call STAT_total( total, TBL_URB(:,:), PV_NAME(I_TBL_URB) )
!       call STAT_total( total, TGL_URB(:,:), PV_NAME(I_TGL_URB) )
    endif

    return
  end subroutine URBAN_vars_total

end module mod_urban_vars
