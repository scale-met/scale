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

  use scale_const, only: &
     I_SW  => CONST_I_SW, &
     I_LW  => CONST_I_LW
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: OCEAN_vars_setup
  public :: OCEAN_vars_restart_read
  public :: OCEAN_vars_restart_write
  public :: OCEAN_vars_history
  public :: OCEAN_vars_total
  public :: OCEAN_vars_external_in

  public :: OCEAN_restart_create
  public :: OCEAN_restart_def_var
  public :: OCEAN_restart_enddef
  public :: OCEAN_restart_write_var
  public :: OCEAN_restart_close

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  logical,               public :: OCEAN_RESTART_OUTPUT       = .false.         !< output restart file?

  character(len=H_LONG), public :: OCEAN_RESTART_IN_BASENAME  = ''              !< basename of the input file
  character(len=H_LONG), public :: OCEAN_RESTART_OUT_BASENAME = ''              !< basename of the output file
  character(len=H_MID),  public :: OCEAN_RESTART_OUT_TITLE    = 'OCEAN restart' !< title    of the output file
  character(len=H_MID),  public :: OCEAN_RESTART_OUT_DTYPE    = 'DEFAULT'       !< REAL4 or REAL8

  ! prognostic variables
  real(RP), public, allocatable :: OCEAN_TEMP      (:,:)   !< temperature at uppermost ocean layer [K]
  real(RP), public, allocatable :: OCEAN_SFC_TEMP  (:,:)   !< ocean surface skin temperature [K]
  real(RP), public, allocatable :: OCEAN_SFC_albedo(:,:,:) !< ocean surface albedo [0-1]
  real(RP), public, allocatable :: OCEAN_SFC_Z0M   (:,:)   !< ocean surface roughness length for momentum [m]
  real(RP), public, allocatable :: OCEAN_SFC_Z0H   (:,:)   !< ocean surface roughness length for heat [m]
  real(RP), public, allocatable :: OCEAN_SFC_Z0E   (:,:)   !< ocean surface roughness length for vapor [m]

  ! tendency variables
  real(RP), public, allocatable :: OCEAN_TEMP_t      (:,:)   !< tendency of OCEAN_TEMP
  real(RP), public, allocatable :: OCEAN_SFC_TEMP_t  (:,:)   !< tendency of OCEAN_SFC_TEMP
  real(RP), public, allocatable :: OCEAN_SFC_albedo_t(:,:,:) !< tendency of OCEAN_SFC_alebdo
  real(RP), public, allocatable :: OCEAN_SFC_Z0M_t   (:,:)   !< tendency of OCEAN_SFC_Z0M
  real(RP), public, allocatable :: OCEAN_SFC_Z0H_t   (:,:)   !< tendency of OCEAN_SFC_Z0H
  real(RP), public, allocatable :: OCEAN_SFC_Z0E_t   (:,:)   !< tendency of OCEAN_SFC_Z0E

  ! surface variables for restart
  real(RP), public, allocatable :: OCEAN_SFLX_MW  (:,:) !< ocean surface w-momentum flux    [kg/m2/s]
  real(RP), public, allocatable :: OCEAN_SFLX_MU  (:,:) !< ocean surface u-momentum flux    [kg/m2/s]
  real(RP), public, allocatable :: OCEAN_SFLX_MV  (:,:) !< ocean surface v-momentum flux    [kg/m2/s]
  real(RP), public, allocatable :: OCEAN_SFLX_SH  (:,:) !< ocean surface sensible heat flux [J/m2/s]
  real(RP), public, allocatable :: OCEAN_SFLX_LH  (:,:) !< ocean surface latent heat flux   [J/m2/s]
  real(RP), public, allocatable :: OCEAN_SFLX_WH  (:,:) !< ocean surface water heat flux    [J/m2/s]
  real(RP), public, allocatable :: OCEAN_SFLX_evap(:,:) !< ocean surface water vapor flux   [kg/m2/s]

  ! diagnostic variables
  real(RP), public, allocatable :: OCEAN_U10(:,:) !< ocean surface velocity u at 10m [m/s]
  real(RP), public, allocatable :: OCEAN_V10(:,:) !< ocean surface velocity v at 10m [m/s]
  real(RP), public, allocatable :: OCEAN_T2 (:,:) !< ocean surface temperature at 2m [K]
  real(RP), public, allocatable :: OCEAN_Q2 (:,:) !< ocean surface water vapor at 2m [kg/kg]

  ! recieved atmospheric variables
  real(RP), public, allocatable :: ATMOS_TEMP     (:,:)
  real(RP), public, allocatable :: ATMOS_PRES     (:,:)
  real(RP), public, allocatable :: ATMOS_W        (:,:)
  real(RP), public, allocatable :: ATMOS_U        (:,:)
  real(RP), public, allocatable :: ATMOS_V        (:,:)
  real(RP), public, allocatable :: ATMOS_DENS     (:,:)
  real(RP), public, allocatable :: ATMOS_QV       (:,:)
  real(RP), public, allocatable :: ATMOS_PBL      (:,:)
  real(RP), public, allocatable :: ATMOS_SFC_PRES (:,:)
  real(RP), public, allocatable :: ATMOS_SFLX_LW  (:,:)
  real(RP), public, allocatable :: ATMOS_SFLX_SW  (:,:)
  real(RP), public, allocatable :: ATMOS_cosSZA   (:,:)
  real(RP), public, allocatable :: ATMOS_SFLX_prec(:,:)

  integer, public :: restart_fid = -1  ! file ID

  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  logical,                private :: OCEAN_VARS_CHECKRANGE      = .false.

  integer,                private, parameter :: VMAX        = 14 !< number of the variables
  integer,                private, parameter :: I_TEMP      =  1
  integer,                private, parameter :: I_SFC_TEMP  =  2
  integer,                private, parameter :: I_ALB_LW    =  3
  integer,                private, parameter :: I_ALB_SW    =  4
  integer,                private, parameter :: I_SFC_Z0M   =  5
  integer,                private, parameter :: I_SFC_Z0H   =  6
  integer,                private, parameter :: I_SFC_Z0E   =  7
  integer,                private, parameter :: I_SFLX_MW   =  8
  integer,                private, parameter :: I_SFLX_MU   =  9
  integer,                private, parameter :: I_SFLX_MV   = 10
  integer,                private, parameter :: I_SFLX_SH   = 11
  integer,                private, parameter :: I_SFLX_LH   = 12
  integer,                private, parameter :: I_SFLX_WH   = 13
  integer,                private, parameter :: I_SFLX_evap = 14

  character(len=H_SHORT), private            :: VAR_NAME(VMAX) !< name  of the variables
  character(len=H_MID),   private            :: VAR_DESC(VMAX) !< desc. of the variables
  character(len=H_SHORT), private            :: VAR_UNIT(VMAX) !< unit  of the variables
  integer,                private            :: VAR_ID(VMAX)   !< ID    of the variables

  data VAR_NAME / 'OCEAN_TEMP',      &
                  'OCEAN_SFC_TEMP',  &
                  'OCEAN_ALB_LW',    &
                  'OCEAN_ALB_SW',    &
                  'OCEAN_SFC_Z0M',   &
                  'OCEAN_SFC_Z0H',   &
                  'OCEAN_SFC_Z0E',   &
                  'OCEAN_SFLX_MW',   &
                  'OCEAN_SFLX_MU',   &
                  'OCEAN_SFLX_MV',   &
                  'OCEAN_SFLX_SH',   &
                  'OCEAN_SFLX_LH',   &
                  'OCEAN_SFLX_WH',   &
                  'OCEAN_SFLX_evap'  /
  data VAR_DESC / 'temperature at uppermost ocean layer', &
                  'ocean surface skin temperature',             &
                  'ocean surface albedo (longwave)',            &
                  'ocean surface albedo (shortwave)',           &
                  'ocean surface roughness length (momentum)',  &
                  'ocean surface roughness length (heat)',      &
                  'ocean surface roughness length (vapor)',     &
                  'ocean surface w-momentum flux',              &
                  'ocean surface u-momentum flux',              &
                  'ocean surface v-momentum flux',              &
                  'ocean surface sensible heat flux',           &
                  'ocean surface latent heat flux',             &
                  'ocean surface water heat flux',              &
                  'ocean surface water vapor flux'              /
  data VAR_UNIT / 'K',       &
                  'K',       &
                  '0-1',     &
                  '0-1',     &
                  'm',       &
                  'm',       &
                  'm',       &
                  'kg/m2/s', &
                  'kg/m2/s', &
                  'kg/m2/s', &
                  'J/m2/s',  &
                  'J/m2/s',  &
                  'J/m2/s',  &
                  'kg/m2/s'  /

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
    integer :: iv
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '++++++ Module[VARS] / Categ[OCEAN] / Origin[SCALE-RM]'

    allocate( OCEAN_TEMP      (IA,JA)   )
    allocate( OCEAN_SFC_TEMP  (IA,JA)   )
    allocate( OCEAN_SFC_albedo(IA,JA,2) )
    allocate( OCEAN_SFC_Z0M   (IA,JA)   )
    allocate( OCEAN_SFC_Z0H   (IA,JA)   )
    allocate( OCEAN_SFC_Z0E   (IA,JA)   )
    OCEAN_TEMP      (:,:)   = UNDEF
    OCEAN_SFC_TEMP  (:,:)   = UNDEF
    OCEAN_SFC_albedo(:,:,:) = UNDEF
    OCEAN_SFC_Z0M   (:,:)   = UNDEF
    OCEAN_SFC_Z0H   (:,:)   = UNDEF
    OCEAN_SFC_Z0E   (:,:)   = UNDEF

    allocate( OCEAN_TEMP_t      (IA,JA)   )
    allocate( OCEAN_SFC_TEMP_t  (IA,JA)   )
    allocate( OCEAN_SFC_albedo_t(IA,JA,2) )
    allocate( OCEAN_SFC_Z0M_t   (IA,JA)   )
    allocate( OCEAN_SFC_Z0H_t   (IA,JA)   )
    allocate( OCEAN_SFC_Z0E_t   (IA,JA)   )
    OCEAN_TEMP_t      (:,:)   = UNDEF
    OCEAN_SFC_TEMP_t  (:,:)   = UNDEF
    OCEAN_SFC_albedo_t(:,:,:) = UNDEF
    OCEAN_SFC_Z0M_t   (:,:)   = UNDEF
    OCEAN_SFC_Z0H_t   (:,:)   = UNDEF
    OCEAN_SFC_Z0E_t   (:,:)   = UNDEF

    allocate( OCEAN_SFLX_MW  (IA,JA) )
    allocate( OCEAN_SFLX_MU  (IA,JA) )
    allocate( OCEAN_SFLX_MV  (IA,JA) )
    allocate( OCEAN_SFLX_SH  (IA,JA) )
    allocate( OCEAN_SFLX_LH  (IA,JA) )
    allocate( OCEAN_SFLX_WH  (IA,JA) )
    allocate( OCEAN_SFLX_evap(IA,JA) )
    OCEAN_SFLX_MW  (:,:) = UNDEF
    OCEAN_SFLX_MU  (:,:) = UNDEF
    OCEAN_SFLX_MV  (:,:) = UNDEF
    OCEAN_SFLX_SH  (:,:) = UNDEF
    OCEAN_SFLX_LH  (:,:) = UNDEF
    OCEAN_SFLX_WH  (:,:) = UNDEF
    OCEAN_SFLX_evap(:,:) = UNDEF

    allocate( OCEAN_U10(IA,JA) )
    allocate( OCEAN_V10(IA,JA) )
    allocate( OCEAN_T2 (IA,JA) )
    allocate( OCEAN_Q2 (IA,JA) )
    OCEAN_U10(:,:) = UNDEF
    OCEAN_V10(:,:) = UNDEF
    OCEAN_T2 (:,:) = UNDEF
    OCEAN_Q2 (:,:) = UNDEF

    allocate( ATMOS_TEMP     (IA,JA) )
    allocate( ATMOS_PRES     (IA,JA) )
    allocate( ATMOS_W        (IA,JA) )
    allocate( ATMOS_U        (IA,JA) )
    allocate( ATMOS_V        (IA,JA) )
    allocate( ATMOS_DENS     (IA,JA) )
    allocate( ATMOS_QV       (IA,JA) )
    allocate( ATMOS_PBL      (IA,JA) )
    allocate( ATMOS_SFC_PRES (IA,JA) )
    allocate( ATMOS_SFLX_LW  (IA,JA) )
    allocate( ATMOS_SFLX_SW  (IA,JA) )
    allocate( ATMOS_cosSZA   (IA,JA) )
    allocate( ATMOS_SFLX_prec(IA,JA) )
    ATMOS_TEMP     (:,:) = UNDEF
    ATMOS_PRES     (:,:) = UNDEF
    ATMOS_W        (:,:) = UNDEF
    ATMOS_U        (:,:) = UNDEF
    ATMOS_V        (:,:) = UNDEF
    ATMOS_DENS     (:,:) = UNDEF
    ATMOS_QV       (:,:) = UNDEF
    ATMOS_PBL      (:,:) = UNDEF
    ATMOS_SFC_PRES (:,:) = UNDEF
    ATMOS_SFLX_LW  (:,:) = UNDEF
    ATMOS_SFLX_SW  (:,:) = UNDEF
    ATMOS_cosSZA   (:,:) = UNDEF
    ATMOS_SFLX_prec(:,:) = UNDEF

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
    if( IO_L ) write(IO_FID_LOG,'(1x,A,A15,A,A32,3(A))') &
               '***       |','VARNAME        ','|', 'DESCRIPTION                     ','[', 'UNIT            ',']'
    do iv = 1, VMAX
       if( IO_L ) write(IO_FID_LOG,'(1x,A,i3,A,A15,A,A32,3(A))') &
                  '*** NO.',iv,'|',VAR_NAME(iv),'|',VAR_DESC(iv),'[',VAR_UNIT(iv),']'

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
    else
       if( IO_L ) write(IO_FID_LOG,*) '*** Restart output? : NO'
       OCEAN_RESTART_OUTPUT = .false.
    endif

    return
  end subroutine OCEAN_vars_setup

  !-----------------------------------------------------------------------------
  !> Read ocean restart
  subroutine OCEAN_vars_restart_read
    use scale_fileio, only: &
       FILEIO_read
    use mod_ocean_admin, only: &
       OCEAN_sw
    implicit none
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** Input restart file (OCEAN) ***'

    if ( OCEAN_sw .and. OCEAN_RESTART_IN_BASENAME /= '' ) then
       if( IO_L ) write(IO_FID_LOG,*) '*** basename: ', trim(OCEAN_RESTART_IN_BASENAME)

       call FILEIO_read( OCEAN_TEMP(:,:),                                               & ! [OUT]
                         OCEAN_RESTART_IN_BASENAME, VAR_NAME(I_TEMP),      'XY', step=1 ) ! [IN]
       call FILEIO_read( OCEAN_SFC_TEMP(:,:),                                           & ! [OUT]
                         OCEAN_RESTART_IN_BASENAME, VAR_NAME(I_SFC_TEMP),  'XY', step=1 ) ! [IN]
       call FILEIO_read( OCEAN_SFC_albedo(:,:,I_LW),                                    & ! [OUT]
                         OCEAN_RESTART_IN_BASENAME, VAR_NAME(I_ALB_LW),    'XY', step=1 ) ! [IN]
       call FILEIO_read( OCEAN_SFC_albedo(:,:,I_SW),                                    & ! [OUT]
                         OCEAN_RESTART_IN_BASENAME, VAR_NAME(I_ALB_SW),    'XY', step=1 ) ! [IN]
       call FILEIO_read( OCEAN_SFC_Z0M(:,:),                                            & ! [OUT]
                         OCEAN_RESTART_IN_BASENAME, VAR_NAME(I_SFC_Z0M),   'XY', step=1 ) ! [IN]
       call FILEIO_read( OCEAN_SFC_Z0H(:,:),                                            & ! [OUT]
                         OCEAN_RESTART_IN_BASENAME, VAR_NAME(I_SFC_Z0H),   'XY', step=1 ) ! [IN]
       call FILEIO_read( OCEAN_SFC_Z0E(:,:),                                            & ! [OUT]
                         OCEAN_RESTART_IN_BASENAME, VAR_NAME(I_SFC_Z0E),   'XY', step=1 ) ! [IN]

       call FILEIO_read( OCEAN_SFLX_MW(:,:),                                            & ! [OUT]
                         OCEAN_RESTART_IN_BASENAME, VAR_NAME(I_SFLX_MW),   'XY', step=1 ) ! [IN]
       call FILEIO_read( OCEAN_SFLX_MU(:,:),                                            & ! [OUT]
                         OCEAN_RESTART_IN_BASENAME, VAR_NAME(I_SFLX_MU),   'XY', step=1 ) ! [IN]
       call FILEIO_read( OCEAN_SFLX_MV(:,:),                                            & ! [OUT]
                         OCEAN_RESTART_IN_BASENAME, VAR_NAME(I_SFLX_MV),   'XY', step=1 ) ! [IN]
       call FILEIO_read( OCEAN_SFLX_SH(:,:),                                            & ! [OUT]
                         OCEAN_RESTART_IN_BASENAME, VAR_NAME(I_SFLX_SH),   'XY', step=1 ) ! [IN]
       call FILEIO_read( OCEAN_SFLX_LH(:,:),                                            & ! [OUT]
                         OCEAN_RESTART_IN_BASENAME, VAR_NAME(I_SFLX_LH),   'XY', step=1 ) ! [IN]
       call FILEIO_read( OCEAN_SFLX_WH(:,:),                                            & ! [OUT]
                         OCEAN_RESTART_IN_BASENAME, VAR_NAME(I_SFLX_WH),   'XY', step=1 ) ! [IN]
       call FILEIO_read( OCEAN_SFLX_evap(:,:),                                          & ! [OUT]
                         OCEAN_RESTART_IN_BASENAME, VAR_NAME(I_SFLX_evap), 'XY', step=1 ) ! [IN]

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
    use mod_ocean_admin, only: &
       OCEAN_sw
    implicit none

    character(len=20)     :: timelabel
    character(len=H_LONG) :: basename
    !---------------------------------------------------------------------------

    if ( OCEAN_sw .and. OCEAN_RESTART_OUT_BASENAME /= '' ) then

       call TIME_gettimelabel( timelabel )
       write(basename,'(A,A,A)') trim(OCEAN_RESTART_OUT_BASENAME), '_', trim(timelabel)

       if( IO_L ) write(IO_FID_LOG,*)
       if( IO_L ) write(IO_FID_LOG,*) '*** Output restart file (OCEAN) ***'
       if( IO_L ) write(IO_FID_LOG,*) '*** basename: ', trim(basename)

       call OCEAN_vars_total

       call FILEIO_write( OCEAN_TEMP(:,:),            basename,                OCEAN_RESTART_OUT_TITLE, & ! [IN]
                          VAR_NAME(I_TEMP),           VAR_DESC(I_TEMP),        VAR_UNIT(I_TEMP),        & ! [IN]
                          'XY',                       OCEAN_RESTART_OUT_DTYPE, nohalo=.true.            ) ! [IN]
       call FILEIO_write( OCEAN_SFC_TEMP(:,:),        basename,                OCEAN_RESTART_OUT_TITLE, & ! [IN]
                          VAR_NAME(I_SFC_TEMP),       VAR_DESC(I_SFC_TEMP),    VAR_UNIT(I_SFC_TEMP),    & ! [IN]
                          'XY',                       OCEAN_RESTART_OUT_DTYPE, nohalo=.true.            ) ! [IN]
       call FILEIO_write( OCEAN_SFC_albedo(:,:,I_LW), basename,                OCEAN_RESTART_OUT_TITLE, & ! [IN]
                          VAR_NAME(I_ALB_LW),         VAR_DESC(I_ALB_LW),      VAR_UNIT(I_ALB_LW),      & ! [IN]
                          'XY',                       OCEAN_RESTART_OUT_DTYPE, nohalo=.true.            ) ! [IN]
       call FILEIO_write( OCEAN_SFC_albedo(:,:,I_SW), basename,                OCEAN_RESTART_OUT_TITLE, & ! [IN]
                          VAR_NAME(I_ALB_SW),         VAR_DESC(I_ALB_SW),      VAR_UNIT(I_ALB_SW),      & ! [IN]
                          'XY',                       OCEAN_RESTART_OUT_DTYPE, nohalo=.true.            ) ! [IN]
       call FILEIO_write( OCEAN_SFC_Z0M(:,:),         basename,                OCEAN_RESTART_OUT_TITLE, & ! [IN]
                          VAR_NAME(I_SFC_Z0M),        VAR_DESC(I_SFC_Z0M),     VAR_UNIT(I_SFC_Z0M),     & ! [IN]
                          'XY',                       OCEAN_RESTART_OUT_DTYPE, nohalo=.true.            ) ! [IN]
       call FILEIO_write( OCEAN_SFC_Z0H(:,:),         basename,                OCEAN_RESTART_OUT_TITLE, & ! [IN]
                          VAR_NAME(I_SFC_Z0H),        VAR_DESC(I_SFC_Z0H),     VAR_UNIT(I_SFC_Z0H),     & ! [IN]
                          'XY',                       OCEAN_RESTART_OUT_DTYPE, nohalo=.true.            ) ! [IN]
       call FILEIO_write( OCEAN_SFC_Z0E(:,:),         basename,                OCEAN_RESTART_OUT_TITLE, & ! [IN]
                          VAR_NAME(I_SFC_Z0E),        VAR_DESC(I_SFC_Z0E),     VAR_UNIT(I_SFC_Z0E),     & ! [IN]
                          'XY',                       OCEAN_RESTART_OUT_DTYPE, nohalo=.true.            ) ! [IN]
       call FILEIO_write( OCEAN_SFLX_MW(:,:),         basename,                OCEAN_RESTART_OUT_TITLE, & ! [IN]
                          VAR_NAME(I_SFLX_MW),        VAR_DESC(I_SFLX_MW),     VAR_UNIT(I_SFLX_MW),     & ! [IN]
                          'XY',                       OCEAN_RESTART_OUT_DTYPE, nohalo=.true.            ) ! [IN]
       call FILEIO_write( OCEAN_SFLX_MU(:,:),         basename,                OCEAN_RESTART_OUT_TITLE, & ! [IN]
                          VAR_NAME(I_SFLX_MU),        VAR_DESC(I_SFLX_MU),     VAR_UNIT(I_SFLX_MU),     & ! [IN]
                          'XY',                       OCEAN_RESTART_OUT_DTYPE, nohalo=.true.            ) ! [IN]
       call FILEIO_write( OCEAN_SFLX_MV(:,:),         basename,                OCEAN_RESTART_OUT_TITLE, & ! [IN]
                          VAR_NAME(I_SFLX_MV),        VAR_DESC(I_SFLX_MV),     VAR_UNIT(I_SFLX_MV),     & ! [IN]
                          'XY',                       OCEAN_RESTART_OUT_DTYPE, nohalo=.true.            ) ! [IN]
       call FILEIO_write( OCEAN_SFLX_SH(:,:),         basename,                OCEAN_RESTART_OUT_TITLE, & ! [IN]
                          VAR_NAME(I_SFLX_SH),        VAR_DESC(I_SFLX_SH),     VAR_UNIT(I_SFLX_SH),     & ! [IN]
                          'XY',                       OCEAN_RESTART_OUT_DTYPE, nohalo=.true.            ) ! [IN]
       call FILEIO_write( OCEAN_SFLX_LH(:,:),         basename,                OCEAN_RESTART_OUT_TITLE, & ! [IN]
                          VAR_NAME(I_SFLX_LH),        VAR_DESC(I_SFLX_LH),     VAR_UNIT(I_SFLX_LH),     & ! [IN]
                          'XY',                       OCEAN_RESTART_OUT_DTYPE, nohalo=.true.            ) ! [IN]
       call FILEIO_write( OCEAN_SFLX_WH(:,:),         basename,                OCEAN_RESTART_OUT_TITLE, & ! [IN]
                          VAR_NAME(I_SFLX_WH),        VAR_DESC(I_SFLX_WH),     VAR_UNIT(I_SFLX_WH),     & ! [IN]
                          'XY',                       OCEAN_RESTART_OUT_DTYPE, nohalo=.true.            ) ! [IN]
       call FILEIO_write( OCEAN_SFLX_evap(:,:),       basename,                OCEAN_RESTART_OUT_TITLE, & ! [IN]
                          VAR_NAME(I_SFLX_evap),      VAR_DESC(I_SFLX_evap),   VAR_UNIT(I_SFLX_evap),   & ! [IN]
                          'XY',                       OCEAN_RESTART_OUT_DTYPE, nohalo=.true.            ) ! [IN]

    endif

    return
  end subroutine OCEAN_vars_restart_write

  !-----------------------------------------------------------------------------
  !> History output set for ocean variables
  subroutine OCEAN_vars_history
    use scale_history, only: &
       HIST_in
    implicit none
    !---------------------------------------------------------------------------

    if ( OCEAN_VARS_CHECKRANGE ) then
       call VALCHECK( OCEAN_TEMP      (IS:IE,JS:JE),      0.0_RP, 1000.0_RP, VAR_NAME(I_TEMP),     &
                     __FILE__, __LINE__ )
       call VALCHECK( OCEAN_SFC_TEMP  (IS:IE,JS:JE),      0.0_RP, 1000.0_RP, VAR_NAME(I_SFC_TEMP), &
                     __FILE__, __LINE__ )
       call VALCHECK( OCEAN_SFC_albedo(IS:IE,JS:JE,I_LW), 0.0_RP,    2.0_RP, VAR_NAME(I_ALB_LW),   &
                     __FILE__, __LINE__ )
       call VALCHECK( OCEAN_SFC_albedo(IS:IE,JS:JE,I_SW), 0.0_RP,    2.0_RP, VAR_NAME(I_ALB_SW),   &
                     __FILE__, __LINE__ )
       call VALCHECK( OCEAN_SFC_Z0M   (IS:IE,JS:JE),      0.0_RP, 1000.0_RP, VAR_NAME(I_SFC_Z0M),  &
                     __FILE__, __LINE__ )
       call VALCHECK( OCEAN_SFC_Z0H   (IS:IE,JS:JE),      0.0_RP, 1000.0_RP, VAR_NAME(I_SFC_Z0H),  &
                     __FILE__, __LINE__ )
       call VALCHECK( OCEAN_SFC_Z0E   (IS:IE,JS:JE),      0.0_RP, 1000.0_RP, VAR_NAME(I_SFC_Z0E),  &
                     __FILE__, __LINE__ )
    endif

    call HIST_in( OCEAN_TEMP      (:,:),      VAR_NAME(I_TEMP),     VAR_DESC(I_TEMP),     VAR_UNIT(I_TEMP)     )
    call HIST_in( OCEAN_SFC_TEMP  (:,:),      VAR_NAME(I_SFC_TEMP), VAR_DESC(I_SFC_TEMP), VAR_UNIT(I_SFC_TEMP) )
    call HIST_in( OCEAN_SFC_albedo(:,:,I_LW), VAR_NAME(I_ALB_LW),   VAR_DESC(I_ALB_LW),   VAR_UNIT(I_ALB_LW)   )
    call HIST_in( OCEAN_SFC_albedo(:,:,I_SW), VAR_NAME(I_ALB_SW),   VAR_DESC(I_ALB_SW),   VAR_UNIT(I_ALB_SW)   )
    call HIST_in( OCEAN_SFC_Z0M   (:,:),      VAR_NAME(I_SFC_Z0M),  VAR_DESC(I_SFC_Z0M),  VAR_UNIT(I_SFC_Z0M)  )
    call HIST_in( OCEAN_SFC_Z0H   (:,:),      VAR_NAME(I_SFC_Z0H),  VAR_DESC(I_SFC_Z0H),  VAR_UNIT(I_SFC_Z0H)  )
    call HIST_in( OCEAN_SFC_Z0E   (:,:),      VAR_NAME(I_SFC_Z0E),  VAR_DESC(I_SFC_Z0E),  VAR_UNIT(I_SFC_Z0E)  )

    call HIST_in( OCEAN_SFLX_MW  (:,:), VAR_NAME(I_SFLX_MW),   VAR_DESC(I_SFLX_MW),   VAR_UNIT(I_SFLX_MW)   )
    call HIST_in( OCEAN_SFLX_MU  (:,:), VAR_NAME(I_SFLX_MU),   VAR_DESC(I_SFLX_MU),   VAR_UNIT(I_SFLX_MU)   )
    call HIST_in( OCEAN_SFLX_MV  (:,:), VAR_NAME(I_SFLX_MV),   VAR_DESC(I_SFLX_MV),   VAR_UNIT(I_SFLX_MV)   )
    call HIST_in( OCEAN_SFLX_SH  (:,:), VAR_NAME(I_SFLX_SH),   VAR_DESC(I_SFLX_SH),   VAR_UNIT(I_SFLX_SH)   )
    call HIST_in( OCEAN_SFLX_LH  (:,:), VAR_NAME(I_SFLX_LH),   VAR_DESC(I_SFLX_LH),   VAR_UNIT(I_SFLX_LH)   )
    call HIST_in( OCEAN_SFLX_WH  (:,:), VAR_NAME(I_SFLX_WH),   VAR_DESC(I_SFLX_WH),   VAR_UNIT(I_SFLX_WH)   )
    call HIST_in( OCEAN_SFLX_evap(:,:), VAR_NAME(I_SFLX_evap), VAR_DESC(I_SFLX_evap), VAR_UNIT(I_SFLX_evap) )

    return
  end subroutine OCEAN_vars_history

  !-----------------------------------------------------------------------------
  !> Budget monitor for ocean
  subroutine OCEAN_vars_total
    use scale_rm_statistics, only: &
       STATISTICS_checktotal, &
       STAT_total
    implicit none

    real(RP) :: total
    !---------------------------------------------------------------------------

    if ( STATISTICS_checktotal ) then

       call STAT_total( total, OCEAN_TEMP      (:,:),      VAR_NAME(I_TEMP)     )
       call STAT_total( total, OCEAN_SFC_TEMP  (:,:),      VAR_NAME(I_SFC_TEMP) )
       call STAT_total( total, OCEAN_SFC_albedo(:,:,I_LW), VAR_NAME(I_ALB_LW)   )
       call STAT_total( total, OCEAN_SFC_albedo(:,:,I_SW), VAR_NAME(I_ALB_SW)   )
       call STAT_total( total, OCEAN_SFC_Z0M   (:,:),      VAR_NAME(I_SFC_Z0M)  )
       call STAT_total( total, OCEAN_SFC_Z0H   (:,:),      VAR_NAME(I_SFC_Z0H)  )
       call STAT_total( total, OCEAN_SFC_Z0E   (:,:),      VAR_NAME(I_SFC_Z0E)  )

       call STAT_total( total, OCEAN_SFLX_MW  (:,:), VAR_NAME(I_SFLX_MW)   )
       call STAT_total( total, OCEAN_SFLX_MU  (:,:), VAR_NAME(I_SFLX_MU)   )
       call STAT_total( total, OCEAN_SFLX_MV  (:,:), VAR_NAME(I_SFLX_MV)   )
       call STAT_total( total, OCEAN_SFLX_SH  (:,:), VAR_NAME(I_SFLX_SH)   )
       call STAT_total( total, OCEAN_SFLX_LH  (:,:), VAR_NAME(I_SFLX_LH)   )
       call STAT_total( total, OCEAN_SFLX_WH  (:,:), VAR_NAME(I_SFLX_WH)   )
       call STAT_total( total, OCEAN_SFLX_evap(:,:), VAR_NAME(I_SFLX_evap) )

    endif

    return
  end subroutine OCEAN_vars_total

  !-----------------------------------------------------------------------------
  !> Input from External I/O
  subroutine OCEAN_vars_external_in( &
       OCEAN_TEMP_in,       &
       OCEAN_SFC_TEMP_in,   &
       OCEAN_SFC_albedo_in, &
       OCEAN_SFC_Z0M_in,    &
       OCEAN_SFC_Z0H_in,    &
       OCEAN_SFC_Z0E_in     )
    implicit none

    real(RP), intent(in) :: OCEAN_TEMP_in      (IA,JA)
    real(RP), intent(in) :: OCEAN_SFC_TEMP_in  (IA,JA)
    real(RP), intent(in) :: OCEAN_SFC_albedo_in(IA,JA,2)
    real(RP), intent(in) :: OCEAN_SFC_Z0M_in   (IA,JA)
    real(RP), intent(in) :: OCEAN_SFC_Z0H_in   (IA,JA)
    real(RP), intent(in) :: OCEAN_SFC_Z0E_in   (IA,JA)
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** External Input file (ocean) ***'

    OCEAN_TEMP      (:,:)   = OCEAN_TEMP_in      (:,:)
    OCEAN_SFC_TEMP  (:,:)   = OCEAN_SFC_TEMP_in  (:,:)
    OCEAN_SFC_albedo(:,:,:) = OCEAN_SFC_albedo_in(:,:,:)
    OCEAN_SFC_Z0M   (:,:)   = OCEAN_SFC_Z0M_in   (:,:)
    OCEAN_SFC_Z0H   (:,:)   = OCEAN_SFC_Z0H_in   (:,:)
    OCEAN_SFC_Z0E   (:,:)   = OCEAN_SFC_Z0E_in   (:,:)

    OCEAN_SFLX_MW  (:,:) = 0.0_RP
    OCEAN_SFLX_MU  (:,:) = 0.0_RP
    OCEAN_SFLX_MV  (:,:) = 0.0_RP
    OCEAN_SFLX_SH  (:,:) = 0.0_RP
    OCEAN_SFLX_LH  (:,:) = 0.0_RP
    OCEAN_SFLX_WH  (:,:) = 0.0_RP
    OCEAN_SFLX_evap(:,:) = 0.0_RP

    call OCEAN_vars_total

    return
  end subroutine OCEAN_vars_external_in

  !-----------------------------------------------------------------------------
  !> Create ocean restart file
  subroutine OCEAN_restart_create
    use scale_time, only: &
       TIME_gettimelabel
    use scale_fileio, only: &
       FILEIO_create
    use mod_ocean_admin, only: &
       OCEAN_sw
    implicit none

    character(len=20)     :: timelabel
    character(len=H_LONG) :: basename
    !---------------------------------------------------------------------------

    if ( OCEAN_sw .and. OCEAN_RESTART_OUT_BASENAME /= '' ) then

       call TIME_gettimelabel( timelabel )
       write(basename,'(A,A,A)') trim(OCEAN_RESTART_OUT_BASENAME), '_', trim(timelabel)

       if( IO_L ) write(IO_FID_LOG,*)
       if( IO_L ) write(IO_FID_LOG,*) '*** Output restart file (OCEAN) ***'
       if( IO_L ) write(IO_FID_LOG,*) '*** basename: ', trim(basename)

       call FILEIO_create(restart_fid, basename, OCEAN_RESTART_OUT_TITLE, OCEAN_RESTART_OUT_DTYPE)
    endif

    return
  end subroutine OCEAN_restart_create

  !-----------------------------------------------------------------------------
  !> Exit netCDF define mode
  subroutine OCEAN_restart_enddef
    use scale_fileio, only: &
       FILEIO_enddef
    implicit none

    if ( restart_fid .NE. -1 ) then
       call FILEIO_enddef( restart_fid ) ! [IN]
    endif

    return
  end subroutine OCEAN_restart_enddef

  !-----------------------------------------------------------------------------
  !> Close restart file
  subroutine OCEAN_restart_close
    use scale_fileio, only: &
       FILEIO_close
    implicit none

    if ( restart_fid .NE. -1 ) then
       call FILEIO_close( restart_fid ) ! [IN]
       restart_fid = -1
    endif

    return
  end subroutine OCEAN_restart_close

  !-----------------------------------------------------------------------------
  !> Define ocean variables in restart file
  subroutine OCEAN_restart_def_var
    use scale_fileio, only: &
       FILEIO_def_var
    implicit none

    !---------------------------------------------------------------------------

    if ( restart_fid .NE. -1 ) then

       call FILEIO_def_var( restart_fid, VAR_ID(I_TEMP),      VAR_NAME(I_TEMP),      VAR_DESC(I_TEMP),      &
                            VAR_UNIT(I_TEMP),      'XY', OCEAN_RESTART_OUT_DTYPE)
       call FILEIO_def_var( restart_fid, VAR_ID(I_SFC_TEMP),  VAR_NAME(I_SFC_TEMP),  VAR_DESC(I_SFC_TEMP),  &
                            VAR_UNIT(I_SFC_TEMP),  'XY', OCEAN_RESTART_OUT_DTYPE)
       call FILEIO_def_var( restart_fid, VAR_ID(I_ALB_LW),    VAR_NAME(I_ALB_LW),    VAR_DESC(I_ALB_LW),    &
                            VAR_UNIT(I_ALB_LW),    'XY', OCEAN_RESTART_OUT_DTYPE)
       call FILEIO_def_var( restart_fid, VAR_ID(I_ALB_SW),    VAR_NAME(I_ALB_SW),    VAR_DESC(I_ALB_SW),    &
                            VAR_UNIT(I_ALB_SW),    'XY', OCEAN_RESTART_OUT_DTYPE)
       call FILEIO_def_var( restart_fid, VAR_ID(I_SFC_Z0M),   VAR_NAME(I_SFC_Z0M),   VAR_DESC(I_SFC_Z0M),   &
                            VAR_UNIT(I_SFC_Z0M),   'XY', OCEAN_RESTART_OUT_DTYPE)
       call FILEIO_def_var( restart_fid, VAR_ID(I_SFC_Z0H),   VAR_NAME(I_SFC_Z0H),   VAR_DESC(I_SFC_Z0H),   &
                            VAR_UNIT(I_SFC_Z0H),   'XY', OCEAN_RESTART_OUT_DTYPE)
       call FILEIO_def_var( restart_fid, VAR_ID(I_SFC_Z0E),   VAR_NAME(I_SFC_Z0E),   VAR_DESC(I_SFC_Z0E),   &
                            VAR_UNIT(I_SFC_Z0E),   'XY', OCEAN_RESTART_OUT_DTYPE)
       call FILEIO_def_var( restart_fid, VAR_ID(I_SFLX_MW),   VAR_NAME(I_SFLX_MW),   VAR_DESC(I_SFLX_MW),   &
                            VAR_UNIT(I_SFLX_MW),   'XY', OCEAN_RESTART_OUT_DTYPE)
       call FILEIO_def_var( restart_fid, VAR_ID(I_SFLX_MU),   VAR_NAME(I_SFLX_MU),   VAR_DESC(I_SFLX_MU),   &
                             VAR_UNIT(I_SFLX_MU),  'XY', OCEAN_RESTART_OUT_DTYPE)
       call FILEIO_def_var( restart_fid, VAR_ID(I_SFLX_MV),   VAR_NAME(I_SFLX_MV),   VAR_DESC(I_SFLX_MV),   &
                            VAR_UNIT(I_SFLX_MV),   'XY', OCEAN_RESTART_OUT_DTYPE)
       call FILEIO_def_var( restart_fid, VAR_ID(I_SFLX_SH),   VAR_NAME(I_SFLX_SH),   VAR_DESC(I_SFLX_SH),   &
                            VAR_UNIT(I_SFLX_SH),   'XY', OCEAN_RESTART_OUT_DTYPE)
       call FILEIO_def_var( restart_fid, VAR_ID(I_SFLX_LH),   VAR_NAME(I_SFLX_LH),   VAR_DESC(I_SFLX_LH),   &
                            VAR_UNIT(I_SFLX_LH),   'XY', OCEAN_RESTART_OUT_DTYPE)
       call FILEIO_def_var( restart_fid, VAR_ID(I_SFLX_WH),   VAR_NAME(I_SFLX_WH),   VAR_DESC(I_SFLX_WH),   &
                            VAR_UNIT(I_SFLX_WH),   'XY', OCEAN_RESTART_OUT_DTYPE)
       call FILEIO_def_var( restart_fid, VAR_ID(I_SFLX_evap), VAR_NAME(I_SFLX_evap), VAR_DESC(I_SFLX_evap), &
                            VAR_UNIT(I_SFLX_evap), 'XY', OCEAN_RESTART_OUT_DTYPE)

    endif

    return
  end subroutine OCEAN_restart_def_var

  !-----------------------------------------------------------------------------
  !> Write ocean variables to restart file
  subroutine OCEAN_restart_write_var
    use scale_fileio, only: &
       FILEIO_write_var
    implicit none

    !---------------------------------------------------------------------------

    if ( restart_fid .NE. -1 ) then

       call OCEAN_vars_total

       call FILEIO_write_var( restart_fid, VAR_ID(I_TEMP), OCEAN_TEMP(:,:),            & ! [IN]
                              VAR_NAME(I_TEMP), 'XY', nohalo=.true.                    ) ! [IN]
       call FILEIO_write_var( restart_fid, VAR_ID(I_TEMP), OCEAN_SFC_TEMP(:,:),        & ! [IN]
                              VAR_NAME(I_SFC_TEMP), 'XY', nohalo=.true.                ) ! [IN]
       call FILEIO_write_var( restart_fid, VAR_ID(I_TEMP), OCEAN_SFC_albedo(:,:,I_LW), & ! [IN]
                              VAR_NAME(I_ALB_LW),  'XY', nohalo=.true.                 ) ! [IN]
       call FILEIO_write_var( restart_fid, VAR_ID(I_TEMP), OCEAN_SFC_albedo(:,:,I_SW), & ! [IN]
                              VAR_NAME(I_ALB_SW), 'XY',  nohalo=.true.                 ) ! [IN]
       call FILEIO_write_var( restart_fid, VAR_ID(I_TEMP), OCEAN_SFC_Z0M(:,:),         & ! [IN]
                              VAR_NAME(I_SFC_Z0M), 'XY', nohalo=.true.                 ) ! [IN]
       call FILEIO_write_var( restart_fid, VAR_ID(I_TEMP), OCEAN_SFC_Z0H(:,:),         & ! [IN]
                              VAR_NAME(I_SFC_Z0H), 'XY', nohalo=.true.                 ) ! [IN]
       call FILEIO_write_var( restart_fid, VAR_ID(I_TEMP), OCEAN_SFC_Z0E(:,:),         & ! [IN]
                              VAR_NAME(I_SFC_Z0E), 'XY', nohalo=.true.                 ) ! [IN]
       call FILEIO_write_var( restart_fid, VAR_ID(I_TEMP), OCEAN_SFLX_MW(:,:),         & ! [IN]
                              VAR_NAME(I_SFLX_MW), 'XY', nohalo=.true.                 ) ! [IN]
       call FILEIO_write_var( restart_fid, VAR_ID(I_TEMP), OCEAN_SFLX_MU(:,:),         & ! [IN]
                              VAR_NAME(I_SFLX_MU), 'XY', nohalo=.true.                 ) ! [IN]
       call FILEIO_write_var( restart_fid, VAR_ID(I_TEMP), OCEAN_SFLX_MV(:,:),         & ! [IN]
                              VAR_NAME(I_SFLX_MV), 'XY', nohalo=.true.                 ) ! [IN]
       call FILEIO_write_var( restart_fid, VAR_ID(I_TEMP), OCEAN_SFLX_SH(:,:),         & ! [IN]
                              VAR_NAME(I_SFLX_SH), 'XY', nohalo=.true.                 ) ! [IN]
       call FILEIO_write_var( restart_fid, VAR_ID(I_TEMP), OCEAN_SFLX_LH(:,:),         & ! [IN]
                              VAR_NAME(I_SFLX_LH), 'XY', nohalo=.true.                 ) ! [IN]
       call FILEIO_write_var( restart_fid, VAR_ID(I_TEMP), OCEAN_SFLX_WH(:,:),         & ! [IN]
                              VAR_NAME(I_SFLX_WH), 'XY', nohalo=.true.                 ) ! [IN]
       call FILEIO_write_var( restart_fid, VAR_ID(I_TEMP), OCEAN_SFLX_evap(:,:),       & ! [IN]
                              VAR_NAME(I_SFLX_evap), 'XY', nohalo=.true.               ) ! [IN]

    endif

    return
  end subroutine OCEAN_restart_write_var

end module mod_ocean_vars
