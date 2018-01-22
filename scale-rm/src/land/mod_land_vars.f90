!-------------------------------------------------------------------------------
!> module LAND Variables
!!
!! @par Description
!!          Container for land variables
!!
!! @author Team SCALE
!! @li      2013-08-31 (T.Yamaura)  [new]
!<
!-------------------------------------------------------------------------------
module mod_land_vars
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_stdio
  use scale_prof
  use scale_debug
  use scale_grid_index
  use scale_land_grid_index

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
  public :: LAND_vars_setup
  public :: LAND_vars_restart_read
  public :: LAND_vars_restart_write
  public :: LAND_vars_history
  public :: LAND_vars_total
  public :: LAND_vars_external_in

  public :: LAND_vars_restart_create
  public :: LAND_vars_restart_open
  public :: LAND_vars_restart_def_var
  public :: LAND_vars_restart_enddef
  public :: LAND_vars_restart_close

  public :: convert_WS2VWC

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  logical,               public :: LAND_RESTART_OUTPUT                 = .false.         !< Output restart file?

  character(len=H_LONG),  public :: LAND_RESTART_IN_BASENAME           = ''              !< Basename of the input  file
  logical,                public :: LAND_RESTART_IN_AGGREGATE                            !< Switch to use aggregate file
  logical,                public :: LAND_RESTART_IN_POSTFIX_TIMELABEL  = .false.         !< Add timelabel to the basename of input  file?
  character(len=H_LONG),  public :: LAND_RESTART_OUT_BASENAME          = ''              !< Basename of the output file
  logical,                public :: LAND_RESTART_OUT_AGGREGATE                           !< Switch to use aggregate file
  logical,                public :: LAND_RESTART_OUT_POSTFIX_TIMELABEL = .true.          !< Add timelabel to the basename of output file?
  character(len=H_MID),   public :: LAND_RESTART_OUT_TITLE             = 'LAND restart' !< Title    of the output file
  character(len=H_SHORT), public :: LAND_RESTART_OUT_DTYPE             = 'DEFAULT'       !< REAL4 or REAL8

  ! prognostic variables
  real(RP), public, allocatable :: LAND_TEMP      (:,:,:) !< temperature of each soil layer [K]
  real(RP), public, allocatable :: LAND_WATER     (:,:,:) !< moisture of each soil layer    [m3/m3]
  real(RP), public, allocatable :: LAND_SFC_TEMP  (:,:)   !< land surface skin temperature  [K]
  real(RP), public, allocatable :: LAND_SFC_albedo(:,:,:) !< land surface albedo            (0-1)

  ! for snow model
  real(RP), public, allocatable :: SNOW_SFC_TEMP  (:,:)    !< snow surface temperature     [K]
  real(RP), public, allocatable :: SNOW_SWE       (:,:)    !< snow water equivalent        [kg/m2]
  real(RP), public, allocatable :: SNOW_Depth     (:,:)    !< snow depth                   [m]
  real(RP), public, allocatable :: SNOW_Dzero     (:,:)    !< snow depth at melting point  [m]
  real(RP), public, allocatable :: SNOW_nosnowsec (:,:)    !< sec while no snow            [s]
  real(RP), public, allocatable :: LAND_type_albedo(:,:,:) !< land surface albedo (Initial) (0-1)

  ! tendency variables
  real(RP), public, allocatable :: LAND_TEMP_t      (:,:,:) !< tendency of LAND_TEMP
  real(RP), public, allocatable :: LAND_WATER_t     (:,:,:) !< tendency of LAND_WATER
  real(RP), public, allocatable :: LAND_SFC_TEMP_t  (:,:)   !< tendency of LAND_SFC_TEMP
  real(RP), public, allocatable :: LAND_SFC_albedo_t(:,:,:) !< tendency of LAND_SFC_albedo

  ! surface variables for restart
  real(RP), public, allocatable :: LAND_SFLX_MW  (:,:) !< land surface w-momentum flux    [kg/m2/s]
  real(RP), public, allocatable :: LAND_SFLX_MU  (:,:) !< land surface u-momentum flux    [kg/m2/s]
  real(RP), public, allocatable :: LAND_SFLX_MV  (:,:) !< land surface v-momentum flux    [kg/m2/s]
  real(RP), public, allocatable :: LAND_SFLX_SH  (:,:) !< land surface sensible heat flux [J/m2/s]
  real(RP), public, allocatable :: LAND_SFLX_LH  (:,:) !< land surface latent heat flux   [J/m2/s]
  real(RP), public, allocatable :: LAND_SFLX_GH  (:,:) !< land surface heat flux          [J/m2/s]
  real(RP), public, allocatable :: LAND_SFLX_evap(:,:) !< land surface water vapor flux   [kg/m2/s]

  ! diagnostic variables
  real(RP), public, allocatable :: LAND_U10(:,:) !< land surface velocity u at 10m [m/s]
  real(RP), public, allocatable :: LAND_V10(:,:) !< land surface velocity v at 10m [m/s]
  real(RP), public, allocatable :: LAND_T2 (:,:) !< land surface temperature at 2m [K]
  real(RP), public, allocatable :: LAND_Q2 (:,:) !< land surface water vapor at 2m [kg/kg]

  ! recieved atmospheric variables
  real(RP), public, allocatable :: ATMOS_TEMP     (:,:)
  real(RP), public, allocatable :: ATMOS_PRES     (:,:)
  real(RP), public, allocatable :: ATMOS_W        (:,:)
  real(RP), public, allocatable :: ATMOS_U        (:,:)
  real(RP), public, allocatable :: ATMOS_V        (:,:)
  real(RP), public, allocatable :: ATMOS_DENS     (:,:)
  real(RP), public, allocatable :: ATMOS_QV       (:,:)
  real(RP), public, allocatable :: ATMOS_PBL      (:,:)
  real(RP), public, allocatable :: ATMOS_SFC_DENS (:,:)
  real(RP), public, allocatable :: ATMOS_SFC_PRES (:,:)
  real(RP), public, allocatable :: ATMOS_SFLX_LW  (:,:)
  real(RP), public, allocatable :: ATMOS_SFLX_SW  (:,:)
  real(RP), public, allocatable :: ATMOS_cosSZA   (:,:)
  real(RP), public, allocatable :: ATMOS_SFLX_prec(:,:)

  real(RP), public, allocatable :: LAND_PROPERTY  (:,:,:) !< land surface property

  character(len=H_LONG), public :: LAND_PROPERTY_IN_FILENAME  = '' !< the file of land parameter table

  integer,  public, parameter   :: LAND_PROPERTY_nmax = 9
  integer,  public, parameter   :: I_WaterLimit       = 1 ! maximum  soil moisture           [m3/m3]
  integer,  public, parameter   :: I_WaterCritical    = 2 ! critical soil moisture           [m3/m3]
  integer,  public, parameter   :: I_StomataResist    = 3 ! stomata resistance               [1/s]
  integer,  public, parameter   :: I_ThermalCond      = 4 ! thermal conductivity for soil    [W/K/m]
  integer,  public, parameter   :: I_HeatCapacity     = 5 ! heat capacity for soil           [J/K/m3]
  integer,  public, parameter   :: I_WaterDiff        = 6 ! moisture diffusivity in the soil [m2/s]
  integer,  public, parameter   :: I_Z0M              = 7 ! roughness length for momemtum    [m]
  integer,  public, parameter   :: I_Z0H              = 8 ! roughness length for heat        [m]
  integer,  public, parameter   :: I_Z0E              = 9 ! roughness length for vapor       [m]

  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  private :: LAND_param_read

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  logical,                private :: LAND_VARS_CHECKRANGE      = .false.

  integer,                private, parameter :: VMAX        = 13 !< number of the variables
  integer,                private, parameter :: I_TEMP      =  1
  integer,                private, parameter :: I_WATER     =  2
  integer,                private, parameter :: I_WATERDS   =  3
  integer,                private, parameter :: I_SFC_TEMP  =  4
  integer,                private, parameter :: I_ALB_LW    =  5
  integer,                private, parameter :: I_ALB_SW    =  6
  integer,                private, parameter :: I_SFLX_MW   =  7
  integer,                private, parameter :: I_SFLX_MU   =  8
  integer,                private, parameter :: I_SFLX_MV   =  9
  integer,                private, parameter :: I_SFLX_SH   = 10
  integer,                private, parameter :: I_SFLX_LH   = 11
  integer,                private, parameter :: I_SFLX_GH   = 12
  integer,                private, parameter :: I_SFLX_evap = 13

  character(len=H_SHORT), private            :: VAR_NAME(VMAX) !< name  of the variables
  character(len=H_MID),   private            :: VAR_DESC(VMAX) !< desc. of the variables
  character(len=H_SHORT), private            :: VAR_UNIT(VMAX) !< unit  of the variables
  integer,                private            :: VAR_ID(VMAX)   !< ID    of the variables
  integer,                private            :: restart_fid = -1  ! file ID

  data VAR_NAME / 'LAND_TEMP',      &
                  'LAND_WATER',     &
                  'LAND_DSAT',      &
                  'LAND_SFC_TEMP',  &
                  'LAND_ALB_LW',    &
                  'LAND_ALB_SW',    &
                  'LAND_SFLX_MW',   &
                  'LAND_SFLX_MU',   &
                  'LAND_SFLX_MV',   &
                  'LAND_SFLX_SH',   &
                  'LAND_SFLX_LH',   &
                  'LAND_SFLX_GH',   &
                  'LAND_SFLX_evap'  /
  data VAR_DESC / 'temperature at each soil layer',  &
                  'moisture at each soil layer',     &
                  'degree of saturation at each soil layer', &
                  'land surface skin temperature',   &
                  'land surface albedo (longwave)',  &
                  'land surface albedo (shortwave)', &
                  'land surface w-momentum flux',    &
                  'land surface u-momentum flux',    &
                  'land surface v-momentum flux',    &
                  'land surface sensible heat flux', &
                  'land surface latent heat flux',   &
                  'land surface ground heat flux',   &
                  'land surface water vapor flux'    /
  data VAR_UNIT / 'K',       &
                  'm3/m3',   &
                  '1',       &
                  'K',       &
                  '1',       &
                  '1',       &
                  'kg/m2/s', &
                  'kg/m2/s', &
                  'kg/m2/s', &
                  'J/m2/s',  &
                  'J/m2/s',  &
                  'J/m2/s',  &
                  'kg/m2/s'  /

  real(RP), private, allocatable :: LAND_PROPERTY_table(:,:)

  integer,  private              :: LAND_QA_comm
  real(RP), private, allocatable :: work_comm(:,:,:) ! for communication

  logical,  private              :: LAND_RESTART_IN_CHECK_COORDINATES = .true.

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine LAND_vars_setup
    use scale_process, only: &
       PRC_MPIstop
    use scale_const, only: &
       UNDEF => CONST_UNDEF
    use scale_comm, only: &
       COMM_vars8, &
       COMM_wait
    use scale_landuse, only: &
       LANDUSE_index_PFT, &
       LANDUSE_PFT_nmax
    implicit none

    NAMELIST / PARAM_LAND_VARS /  &
       LAND_RESTART_IN_BASENAME,           &
       LAND_RESTART_IN_AGGREGATE,          &
       LAND_RESTART_IN_POSTFIX_TIMELABEL,  &
       LAND_RESTART_IN_CHECK_COORDINATES,  &
       LAND_RESTART_OUTPUT,                &
       LAND_RESTART_OUT_BASENAME,          &
       LAND_RESTART_OUT_AGGREGATE,         &
       LAND_RESTART_OUT_POSTFIX_TIMELABEL, &
       LAND_RESTART_OUT_TITLE,             &
       LAND_RESTART_OUT_DTYPE,             &
       LAND_VARS_CHECKRANGE

    integer :: ierr
    integer :: i, j, iv, p
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '++++++ Module[VARS] / Categ[LAND] / Origin[SCALE-RM]'

    allocate( LAND_TEMP      (LKMAX,IA,JA) )
    allocate( LAND_WATER     (LKMAX,IA,JA) )
    allocate( LAND_SFC_TEMP  (IA,JA)       )
    allocate( LAND_SFC_albedo(IA,JA,2)     )
    allocate( LAND_type_albedo(IA,JA,2)     )
    LAND_TEMP      (:,:,:) = UNDEF
    LAND_WATER     (:,:,:) = UNDEF
    LAND_SFC_TEMP  (:,:)   = UNDEF
    LAND_SFC_albedo(:,:,:) = UNDEF
    LAND_type_albedo(:,:,:) = UNDEF

    allocate( SNOW_SFC_TEMP  (IA,JA)       )
    allocate( SNOW_SWE       (IA,JA)       )
    allocate( SNOW_Depth     (IA,JA)       )
    allocate( SNOW_Dzero     (IA,JA)       )
    allocate( SNOW_nosnowsec (IA,JA)       )
    SNOW_SFC_TEMP  (:,:)   = UNDEF
    SNOW_SWE       (:,:)   = UNDEF
    SNOW_Depth     (:,:)   = UNDEF
    SNOW_Dzero     (:,:)   = UNDEF
    SNOW_nosnowsec (:,:)   = UNDEF

    allocate( LAND_TEMP_t      (LKMAX,IA,JA) )
    allocate( LAND_WATER_t     (LKMAX,IA,JA) )
    allocate( LAND_SFC_TEMP_t  (IA,JA)       )
    allocate( LAND_SFC_albedo_t(IA,JA,2)     )
    LAND_TEMP_t      (:,:,:) = UNDEF
    LAND_WATER_t     (:,:,:) = UNDEF
    LAND_SFC_TEMP_t  (:,:)   = UNDEF
    LAND_SFC_albedo_t(:,:,:) = UNDEF

    allocate( LAND_SFLX_MW  (IA,JA) )
    allocate( LAND_SFLX_MU  (IA,JA) )
    allocate( LAND_SFLX_MV  (IA,JA) )
    allocate( LAND_SFLX_SH  (IA,JA) )
    allocate( LAND_SFLX_LH  (IA,JA) )
    allocate( LAND_SFLX_GH  (IA,JA) )
    allocate( LAND_SFLX_evap(IA,JA) )
    LAND_SFLX_MW  (:,:) = UNDEF
    LAND_SFLX_MU  (:,:) = UNDEF
    LAND_SFLX_MV  (:,:) = UNDEF
    LAND_SFLX_SH  (:,:) = UNDEF
    LAND_SFLX_LH  (:,:) = UNDEF
    LAND_SFLX_GH  (:,:) = UNDEF
    LAND_SFLX_evap(:,:) = UNDEF

    allocate( LAND_U10(IA,JA) )
    allocate( LAND_V10(IA,JA) )
    allocate( LAND_T2 (IA,JA) )
    allocate( LAND_Q2 (IA,JA) )
    LAND_U10(:,:) = UNDEF
    LAND_V10(:,:) = UNDEF
    LAND_T2 (:,:) = UNDEF
    LAND_Q2 (:,:) = UNDEF

    allocate( ATMOS_TEMP     (IA,JA) )
    allocate( ATMOS_PRES     (IA,JA) )
    allocate( ATMOS_W        (IA,JA) )
    allocate( ATMOS_U        (IA,JA) )
    allocate( ATMOS_V        (IA,JA) )
    allocate( ATMOS_DENS     (IA,JA) )
    allocate( ATMOS_QV       (IA,JA) )
    allocate( ATMOS_PBL      (IA,JA) )
    allocate( ATMOS_SFC_DENS (IA,JA) )
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
    ATMOS_SFC_DENS (:,:) = UNDEF
    ATMOS_SFC_PRES (:,:) = UNDEF
    ATMOS_SFLX_LW  (:,:) = UNDEF
    ATMOS_SFLX_SW  (:,:) = UNDEF
    ATMOS_cosSZA   (:,:) = UNDEF
    ATMOS_SFLX_prec(:,:) = UNDEF

    LAND_QA_comm = LKMAX &
                 + LKMAX &
                 + 1     &
                 + 2

    allocate( work_comm(IA,JA,LAND_QA_comm) )

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_LAND_VARS,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_LAND_VARS. Check!'
       call PRC_MPIstop
    endif
    if( IO_NML ) write(IO_FID_NML,nml=PARAM_LAND_VARS)

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** List of prognostic variables (LAND) ***'
    if( IO_L ) write(IO_FID_LOG,'(1x,A,A24,A,A48,A,A12,A)') &
               '***       |', 'VARNAME                 ','|', &
               'DESCRIPTION                                     ', '[', 'UNIT        ', ']'
    do iv = 1, VMAX
       if( IO_L ) write(IO_FID_LOG,'(1x,A,I3,A,A24,A,A48,A,A12,A)') &
                  '*** NO.',iv,'|',VAR_NAME(iv),'|',VAR_DESC(iv),'[',VAR_UNIT(iv),']'
    enddo

    if( IO_L ) write(IO_FID_LOG,*)
    if ( LAND_RESTART_IN_BASENAME /= '' ) then
       if( IO_L ) write(IO_FID_LOG,*) '*** Restart input?  : YES, file = ', trim(LAND_RESTART_IN_BASENAME)
       if( IO_L ) write(IO_FID_LOG,*) '*** Add timelabel?  : ', LAND_RESTART_IN_POSTFIX_TIMELABEL
    else
       if( IO_L ) write(IO_FID_LOG,*) '*** Restart input?  : NO'
    endif
    if (       LAND_RESTART_OUTPUT             &
         .AND. LAND_RESTART_OUT_BASENAME /= '' ) then
       if( IO_L ) write(IO_FID_LOG,*) '*** Restart output? : YES, file = ', trim(LAND_RESTART_OUT_BASENAME)
       if( IO_L ) write(IO_FID_LOG,*) '*** Add timelabel?  : ', LAND_RESTART_OUT_POSTFIX_TIMELABEL
    else
       if( IO_L ) write(IO_FID_LOG,*) '*** Restart output? : NO'
       LAND_RESTART_OUTPUT = .false.
    endif

    ! Read land property table
    allocate( LAND_PROPERTY_table(LANDUSE_PFT_nmax,LAND_PROPERTY_nmax) )
    LAND_PROPERTY_table(:,:) = UNDEF

    call LAND_param_read

    ! Apply land property to 2D map
    allocate( LAND_PROPERTY(IA,JA,LAND_PROPERTY_nmax) )

    ! tentative, mosaic is off
    do p = 1, LAND_PROPERTY_nmax
    do j = JS, JE
    do i = IS, IE
       LAND_PROPERTY(i,j,p) = LAND_PROPERTY_table( LANDUSE_index_PFT(i,j,1), p )
    enddo
    enddo
    enddo

    do p = 1, LAND_PROPERTY_nmax
       call COMM_vars8( LAND_PROPERTY(:,:,p), p )
    enddo
    do p = 1, LAND_PROPERTY_nmax
       call COMM_wait ( LAND_PROPERTY(:,:,p), p )
    enddo

    return
  end subroutine LAND_vars_setup

  !-----------------------------------------------------------------------------
  !> Open land restart file for read
  subroutine LAND_vars_restart_open
    use scale_time, only: &
       TIME_gettimelabel
    use scale_file_cartesC, only: &
       FILE_CARTESC_open, &
       FILE_CARTESC_check_coordinates
    use mod_land_admin, only: &
       LAND_sw
    implicit none

    character(len=19)     :: timelabel
    character(len=H_LONG) :: basename
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** Open restart file (LAND) ***'

    if ( LAND_sw .and. LAND_RESTART_IN_BASENAME /= '' ) then

       if ( LAND_RESTART_IN_POSTFIX_TIMELABEL ) then
          call TIME_gettimelabel( timelabel )
          basename = trim(LAND_RESTART_IN_BASENAME)//'_'//trim(timelabel)
       else
          basename = trim(LAND_RESTART_IN_BASENAME)
       endif

       if( IO_L ) write(IO_FID_LOG,*) '*** basename: ', trim(basename)

       call FILE_CARTESC_open( basename, restart_fid, aggregate=LAND_RESTART_IN_AGGREGATE )

       if ( LAND_RESTART_IN_CHECK_COORDINATES ) then
          call FILE_CARTESC_check_coordinates( restart_fid, land=.true. )
       end if

    else
       if( IO_L ) write(IO_FID_LOG,*) '*** restart file for land is not specified.'
    endif

    return
  end subroutine LAND_vars_restart_open

  !-----------------------------------------------------------------------------
  !> Read land restart
  subroutine LAND_vars_restart_read
    use scale_file, only: &
       FILE_get_aggregate
    use scale_file_cartesC, only: &
       FILE_CARTESC_read, &
       FILE_CARTESC_flush
    implicit none
    !---------------------------------------------------------------------------

    if ( restart_fid /= -1 ) then
       if( IO_L ) write(IO_FID_LOG,*)
       if( IO_L ) write(IO_FID_LOG,*) '*** Read from restart file (LAND) ***'

       call FILE_CARTESC_read( restart_fid, VAR_NAME(I_TEMP),      'LXY', & ! [IN]
                               LAND_TEMP (:,:,:)                          ) ! [OUT]
       call FILE_CARTESC_read( restart_fid, VAR_NAME(I_WATER),     'LXY', & ! [IN]
                               LAND_WATER(:,:,:)                          ) ! [OUT]
       call FILE_CARTESC_read( restart_fid, VAR_NAME(I_SFC_TEMP),  'XY',  & ! [IN]
                               LAND_SFC_TEMP(:,:)                         ) ! [OUT]
       call FILE_CARTESC_read( restart_fid, VAR_NAME(I_ALB_LW),    'XY',  & ! [IN]
                               LAND_SFC_albedo(:,:,I_LW)                  ) ! [OUT]
       call FILE_CARTESC_read( restart_fid, VAR_NAME(I_ALB_SW),    'XY',  & ! [IN]
                               LAND_SFC_albedo(:,:,I_SW)                  ) ! [OUT]
       call FILE_CARTESC_read( restart_fid, VAR_NAME(I_SFLX_MW),   'XY',  & ! [IN]
                               LAND_SFLX_MW(:,:)                          ) ! [OUT]
       call FILE_CARTESC_read( restart_fid, VAR_NAME(I_SFLX_MU),   'XY',  & ! [IN]
                               LAND_SFLX_MU(:,:)                          ) ! [OUT]
       call FILE_CARTESC_read( restart_fid, VAR_NAME(I_SFLX_MV),   'XY',  & ! [IN]
                               LAND_SFLX_MV(:,:)                          ) ! [OUT]
       call FILE_CARTESC_read( restart_fid, VAR_NAME(I_SFLX_SH),   'XY',  & ! [IN]
                               LAND_SFLX_SH(:,:)                          ) ! [OUT]
       call FILE_CARTESC_read( restart_fid, VAR_NAME(I_SFLX_LH),   'XY',  & ! [IN]
                               LAND_SFLX_LH(:,:)                          ) ! [OUT]
       call FILE_CARTESC_read( restart_fid, VAR_NAME(I_SFLX_GH),   'XY',  & ! [IN]
                               LAND_SFLX_GH(:,:)                          ) ! [OUT]
       call FILE_CARTESC_read( restart_fid, VAR_NAME(I_SFLX_evap), 'XY',  & ! [IN]
                               LAND_SFLX_evap(:,:)                        ) ! [OUT]

       if( FILE_get_AGGREGATE(restart_fid) ) call FILE_CARTESC_flush( restart_fid ) ! commit all pending read requests

       !call FILE_CARTESC_read( restart_fid, 'SNOW_SFC_TEMP',      'XY',  & ! [OUT]
       !                        SNOW_SFC_TEMP(:,:)                        ) ! [IN]
       !call FILE_CARTESC_read( restart_fid, 'SNOW_SWE',           'XY',  & ! [OUT]
       !                        SNOW_SWE(:,:)                             ) ! [IN]
       !call FILE_CARTESC_read( restart_fid, 'SNOW_Depth',         'XY',  & ! [OUT]
       !                        SNOW_Depth(:,:)                           ) ! [IN]
       !call FILE_CARTESC_read( restart_fid, 'SNOW_Dzero',         'XY',  & ! [OUT]
       !                        SNOW_Dzero(:,:)                           ) ! [IN]
       !call FILE_CARTESC_read( restart_fid, 'SNOW_nosnowsec',     'XY'   & ! [OUT]
       !                        SNOW_nosnowsec(:,:)                       ) ! [IN]

       !!!!! Tentative for snow model !!!!!
       LAND_type_albedo = LAND_SFC_albedo
       SNOW_SFC_TEMP    = 273.15_RP
       SNOW_SWE         = 0.0_RP
       SNOW_Depth       = 0.0_RP
       SNOW_Dzero       = 0.0_RP
       SNOW_nosnowsec   = 0.0_RP

       call LAND_vars_total
    else
       if( IO_L ) write(IO_FID_LOG,*) '*** invalid restart file ID for land.'
    endif

    return
  end subroutine LAND_vars_restart_read

  !-----------------------------------------------------------------------------
  !> History output set for land variables
  subroutine LAND_vars_history
    use scale_file_history, only: &
       FILE_HISTORY_in
    implicit none

    real(RP) :: LAND_WATERDS(LKMAX,IA,JA)
    integer :: k, i, j
    !---------------------------------------------------------------------------

    if ( LAND_VARS_CHECKRANGE ) then
       call VALCHECK( LAND_TEMP      (:,IS:IE,JS:JE),    0.0_RP, 1000.0_RP, VAR_NAME(I_TEMP),       &
                     __FILE__, __LINE__ )
       call VALCHECK( LAND_WATER     (:,IS:IE,JS:JE),    0.0_RP, 1000.0_RP, VAR_NAME(I_WATER),      &
                     __FILE__, __LINE__ )
       call VALCHECK( LAND_SFC_TEMP  (IS:IE,JS:JE),      0.0_RP, 1000.0_RP, VAR_NAME(I_SFC_TEMP),   &
                     __FILE__, __LINE__ )
       call VALCHECK( LAND_SFC_albedo(IS:IE,JS:JE,I_LW), 0.0_RP,    2.0_RP, VAR_NAME(I_ALB_LW),     &
                     __FILE__, __LINE__ )
       call VALCHECK( LAND_SFC_albedo(IS:IE,JS:JE,I_SW), 0.0_RP,    2.0_RP, VAR_NAME(I_ALB_SW),     &
                     __FILE__, __LINE__ )

       !call VALCHECK( SNOW_SFC_TEMP (IS:IE,JS:JE),      0.0_RP, 1000.0_RP, 'SNOW_SFC_TEMP',        &
       !              __FILE__, __LINE__ )
       !call VALCHECK( SNOW_SWE      (IS:IE,JS:JE),      0.0_RP, 1000.0_RP, 'SNOW_SWE',             &
       !              __FILE__, __LINE__ )
       !call VALCHECK( SNOW_Depth    (IS:IE,JS:JE),      0.0_RP, 1000.0_RP, 'SNOW_Depth',           &
       !              __FILE__, __LINE__ )
       !call VALCHECK( SNOW_Dzero    (IS:IE,JS:JE),      0.0_RP, 1000.0_RP, 'SNOW_Dzero',           &
       !              __FILE__, __LINE__ )
    endif

    call FILE_HISTORY_in( LAND_TEMP (:,:,:), VAR_NAME(I_TEMP),  VAR_DESC(I_TEMP),  VAR_UNIT(I_TEMP),  dim_type='LXY' )
    call FILE_HISTORY_in( LAND_WATER(:,:,:), VAR_NAME(I_WATER), VAR_DESC(I_WATER), VAR_UNIT(I_WATER), dim_type='LXY' )
    do j = JS, JE
    do i = IS, IE
    do k = 1, LKMAX
       LAND_WATERDS(k,i,j) = LAND_WATER(k,i,j) / LAND_PROPERTY(i,j,I_WaterLimit)
    end do
    end do
    end do
    call FILE_HISTORY_in( LAND_WATERDS(:,:,:), VAR_NAME(I_WATERDS), VAR_DESC(I_WATERDS), VAR_UNIT(I_WATERDS), dim_type='LXY', fill_halo=.true. )


    call FILE_HISTORY_in( LAND_SFC_TEMP  (:,:),      VAR_NAME(I_SFC_TEMP),   VAR_DESC(I_SFC_TEMP),   VAR_UNIT(I_SFC_TEMP) )
    call FILE_HISTORY_in( LAND_SFC_albedo(:,:,I_LW), VAR_NAME(I_ALB_LW),     VAR_DESC(I_ALB_LW),     VAR_UNIT(I_ALB_LW)   )
    call FILE_HISTORY_in( LAND_SFC_albedo(:,:,I_SW), VAR_NAME(I_ALB_SW),     VAR_DESC(I_ALB_SW),     VAR_UNIT(I_ALB_SW)   )

    call FILE_HISTORY_in( LAND_SFLX_MW  (:,:), VAR_NAME(I_SFLX_MW),   VAR_DESC(I_SFLX_MW),   VAR_UNIT(I_SFLX_MW)   )
    call FILE_HISTORY_in( LAND_SFLX_MU  (:,:), VAR_NAME(I_SFLX_MU),   VAR_DESC(I_SFLX_MU),   VAR_UNIT(I_SFLX_MU)   )
    call FILE_HISTORY_in( LAND_SFLX_MV  (:,:), VAR_NAME(I_SFLX_MV),   VAR_DESC(I_SFLX_MV),   VAR_UNIT(I_SFLX_MV)   )
    call FILE_HISTORY_in( LAND_SFLX_SH  (:,:), VAR_NAME(I_SFLX_SH),   VAR_DESC(I_SFLX_SH),   VAR_UNIT(I_SFLX_SH)   )
    call FILE_HISTORY_in( LAND_SFLX_LH  (:,:), VAR_NAME(I_SFLX_LH),   VAR_DESC(I_SFLX_LH),   VAR_UNIT(I_SFLX_LH)   )
    call FILE_HISTORY_in( LAND_SFLX_GH  (:,:), VAR_NAME(I_SFLX_GH),   VAR_DESC(I_SFLX_GH),   VAR_UNIT(I_SFLX_GH)   )
    call FILE_HISTORY_in( LAND_SFLX_evap(:,:), VAR_NAME(I_SFLX_evap), VAR_DESC(I_SFLX_evap), VAR_UNIT(I_SFLX_evap) )

    ! snow model
    call FILE_HISTORY_in( SNOW_SFC_TEMP (:,:), 'SNOW_SFC_TEMP',  'Snow surface temperature',    'K'      )
    call FILE_HISTORY_in( SNOW_SWE      (:,:), 'SNOW_SWE',       'Snow water equivalent',       'kg/m2'  )
    call FILE_HISTORY_in( SNOW_Depth    (:,:), 'SNOW_Depth',     'Snow depth',                  'm'      )
    call FILE_HISTORY_in( SNOW_Dzero    (:,:), 'SNOW_Dzero',     'Snow depth at melting point', 'm'      )
    call FILE_HISTORY_in( SNOW_nosnowsec(:,:), 'SNOW_nosnowsec', 'Time duration without snow',  's'      )

    return
  end subroutine LAND_vars_history

  !-----------------------------------------------------------------------------
  !> Budget monitor for land
  subroutine LAND_vars_total
    use scale_rm_statistics, only: &
       STATISTICS_checktotal, &
       STAT_total
    implicit none

    real(RP) :: total

    character(len=2) :: sk
    integer          :: k
    !---------------------------------------------------------------------------

    if ( STATISTICS_checktotal ) then

       do k = LKS, LKE
          write(sk,'(I2.2)') k

          call STAT_total( total, LAND_TEMP (k,:,:), trim(VAR_NAME(I_TEMP) )//sk )
          call STAT_total( total, LAND_WATER(k,:,:), trim(VAR_NAME(I_WATER))//sk )
       enddo

       call STAT_total( total, LAND_SFC_TEMP  (:,:),      VAR_NAME(I_SFC_TEMP) )
       call STAT_total( total, LAND_SFC_albedo(:,:,I_LW), VAR_NAME(I_ALB_LW)   )
       call STAT_total( total, LAND_SFC_albedo(:,:,I_SW), VAR_NAME(I_ALB_SW)   )

       call STAT_total( total, SNOW_SFC_TEMP  (:,:),     'SNOW_SFC_TEMP'  )
       call STAT_total( total, SNOW_SWE       (:,:),     'SNOW_SWE'       )
       call STAT_total( total, SNOW_Depth     (:,:),     'SNOW_Depth'     )
       call STAT_total( total, SNOW_Dzero     (:,:),     'SNOW_Dzero'     )

    endif

    return
  end subroutine LAND_vars_total

  !-----------------------------------------------------------------------------
  !> Input from External I/O
  subroutine LAND_vars_external_in( &
      LAND_TEMP_in,      &
      LAND_WATER_in,     &
      LAND_SFC_TEMP_in,  &
      LAND_SFC_albedo_in )
    implicit none

    real(RP), intent(in) :: LAND_TEMP_in (:,:,:)
    real(RP), intent(in) :: LAND_WATER_in(:,:,:)
    real(RP), intent(in) :: LAND_SFC_TEMP_in  (IA,JA)
    real(RP), intent(in) :: LAND_SFC_albedo_in(IA,JA,2)
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** External Input (land) ***'

    LAND_TEMP      (:,:,:) = LAND_TEMP_in      (:,:,:)
    LAND_WATER     (:,:,:) = LAND_WATER_in     (:,:,:)
    LAND_SFC_TEMP  (:,:)   = LAND_SFC_TEMP_in  (:,:)
    LAND_SFC_albedo(:,:,:) = LAND_SFC_albedo_in(:,:,:)

    LAND_SFLX_MW  (:,:) = 0.0_RP
    LAND_SFLX_MU  (:,:) = 0.0_RP
    LAND_SFLX_MV  (:,:) = 0.0_RP
    LAND_SFLX_SH  (:,:) = 0.0_RP
    LAND_SFLX_LH  (:,:) = 0.0_RP
    LAND_SFLX_GH  (:,:) = 0.0_RP
    LAND_SFLX_evap(:,:) = 0.0_RP

    call LAND_vars_total

    return
  end subroutine LAND_vars_external_in

  !-----------------------------------------------------------------------------
  !> Budget monitor for land
  subroutine LAND_param_read
    use scale_process, only: &
       PRC_MPIstop
    use scale_landuse, only: &
       LANDUSE_PFT_nmax
    implicit none

    integer              :: index
    character(len=H_MID) :: description
    real(RP)             :: STRGMAX
    real(RP)             :: STRGCRT
    real(RP)             :: RSTOMA
    real(RP)             :: TCS
    real(RP)             :: HCS
    real(RP)             :: DFW
    real(RP)             :: Z0M
    real(RP)             :: Z0H
    real(RP)             :: Z0E

    NAMELIST / PARAM_LAND_PROPERTY /  &
       LAND_PROPERTY_IN_FILENAME

    NAMELIST / PARAM_LAND_DATA / &
       index,       &
       description, &
       STRGMAX,     &
       STRGCRT,     &
       RSTOMA,      &
       TCS,         &
       HCS,         &
       DFW,         &
       Z0M,         &
       Z0H,         &
       Z0E

    integer :: n
    integer :: ierr

    integer :: IO_FID_LAND_PROPERTY
    !---------------------------------------------------------------------------

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_LAND_PROPERTY,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_LAND_PROPERTY. Check!'
       call PRC_MPIstop
    endif
    if( IO_NML ) write(IO_FID_NML,nml=PARAM_LAND_PROPERTY)

    if( LAND_PROPERTY_IN_FILENAME /= '' ) then
      !--- Open land parameter file
      IO_FID_LAND_PROPERTY = IO_get_available_fid()
      open( IO_FID_LAND_PROPERTY,                     &
            file   = trim(LAND_PROPERTY_IN_FILENAME), &
            form   = 'formatted',                     &
            status = 'old',                           &
            iostat = ierr                             )

      if ( ierr /= 0 ) then
         write(*,*) 'xxx [LAND_param_read] Failed to open land parameter file! :', &
                    trim(LAND_PROPERTY_IN_FILENAME)
         call PRC_MPIstop
      else
        if( IO_L ) write(IO_FID_LOG,*)
        if( IO_L ) write(IO_FID_LOG,*) '*** Properties for each plant functional type (PFT)'
        if( IO_L ) write(IO_FID_LOG,'(12(1x,A))') '***          ', &
                                                  '                     description', &
                                                  'Max Stg.', &
                                                  'CRT Stg.', &
                                                  'Stm.Res.', &
                                                  'T condu.', &
                                                  'H capac.', &
                                                  'DFC Wat.', &
                                                  '   Z0(m)', &
                                                  '   Z0(h)', &
                                                  '   Z0(e)'

        !--- read namelist
        rewind(IO_FID_LAND_PROPERTY)

        do n = 1, LANDUSE_PFT_nmax
           ! default value
           STRGMAX =  0.2_RP
           STRGCRT =  0.1_RP
           RSTOMA  = 50.0_RP
           TCS     =  1.0_RP
           HCS     =  2.E+6_RP
           DFW     =  1.E-6_RP
           Z0M     =  0.1_RP
           Z0H     = -1.0_RP
           Z0E     = -1.0_RP

           read(IO_FID_LAND_PROPERTY,nml=PARAM_LAND_DATA,iostat=ierr)
           if ( ierr < 0 ) then !--- no more data
              exit
           elseif( ierr > 0 ) then !--- fatal error
              write(*,*) 'xxx Not appropriate names in namelist PARAM_LAND_DATA. Check!'
              call PRC_MPIstop
           endif

           if( Z0H < 0.0_RP ) then
             Z0H = Z0M / 7.4_RP ! defined by Garratt and Francey (1978)
           endif
           if( Z0E < 0.0_RP ) then
             Z0E = Z0M / 7.4_RP ! defined by Garratt and Francey (1978)
           endif

           LAND_PROPERTY_table(index,I_WaterLimit   ) = STRGMAX
           LAND_PROPERTY_table(index,I_WaterCritical) = STRGCRT
           LAND_PROPERTY_table(index,I_StomataResist) = RSTOMA
           LAND_PROPERTY_table(index,I_ThermalCond  ) = TCS
           LAND_PROPERTY_table(index,I_HeatCapacity ) = HCS
           LAND_PROPERTY_table(index,I_WaterDiff    ) = DFW
           LAND_PROPERTY_table(index,I_Z0M          ) = Z0M
           LAND_PROPERTY_table(index,I_Z0H          ) = Z0H
           LAND_PROPERTY_table(index,I_Z0E          ) = Z0E

           if( IO_L ) write(IO_FID_LOG,'(1x,A10,I3.3,1x,A32,4(1x,F6.3),2(1x,ES8.1),3(1x,F6.3))') &
                                         '*** IDX = ', index, &
                                         trim(description), &
                                         STRGMAX, &
                                         STRGCRT, &
                                         RSTOMA,  &
                                         TCS,     &
                                         HCS,     &
                                         DFW,     &
                                         Z0M,     &
                                         Z0H,     &
                                         Z0E
        enddo

      end if

      close( IO_FID_LAND_PROPERTY )

    endif

    return
  end subroutine LAND_param_read

  !-----------------------------------------------------------------------------
  !> conversion from water saturation [fraction] to volumetric water content [m3/m3]
  function convert_WS2VWC( WS, critical ) result( VWC )
    implicit none

    real(RP), intent(in) :: WS(IA,JA) ! water saturation [fraction]
    logical,  intent(in) :: critical  ! is I_WaterCritical used?

    real(RP) :: VWC(IA,JA) ! volumetric water content [m3/m3]

    ! work
    integer :: i, j, num
    !---------------------------------------------------------------------------

    if( critical ) then
      num = I_WaterCritical
    else
      num = I_WaterLimit
    end if

    do j = JS, JE
    do i = IS, IE
      VWC(i,j) = max( min( WS(i,j)*LAND_PROPERTY(i,j,num), LAND_PROPERTY(i,j,num) ), 0.0_RP )
    end do
    end do

    return
  end function convert_WS2VWC

  !-----------------------------------------------------------------------------
  !> Create land restart file
  subroutine LAND_vars_restart_create
    use scale_time, only: &
       TIME_gettimelabel
    use scale_file_cartesC, only: &
       FILE_CARTESC_create
    use mod_land_admin, only: &
       LAND_sw
    implicit none

    character(len=19)     :: timelabel
    character(len=H_LONG) :: basename
    !---------------------------------------------------------------------------

    if ( LAND_sw .and. LAND_RESTART_OUT_BASENAME /= '' ) then

       if( IO_L ) write(IO_FID_LOG,*)
       if( IO_L ) write(IO_FID_LOG,*) '*** Create restart file (LAND) ***'

       if ( LAND_RESTART_OUT_POSTFIX_TIMELABEL ) then
          call TIME_gettimelabel( timelabel )
          basename = trim(LAND_RESTART_OUT_BASENAME)//'_'//trim(timelabel)
       else
          basename = trim(LAND_RESTART_OUT_BASENAME)
       endif

       if( IO_L ) write(IO_FID_LOG,*) '*** basename: ', trim(basename)

       call FILE_CARTESC_create( &
            basename, LAND_RESTART_OUT_TITLE, LAND_RESTART_OUT_DTYPE, & ! [IN]
            restart_fid,                                              & ! [OUT]
            aggregate=LAND_RESTART_OUT_AGGREGATE                      ) ! [IN]

    endif

    return
  end subroutine LAND_vars_restart_create

  !-----------------------------------------------------------------------------
  !> Exit netCDF define mode
  subroutine LAND_vars_restart_enddef
    use scale_file_cartesC, only: &
       FILE_CARTESC_enddef
    implicit none

    if ( restart_fid /= -1 ) then
       call FILE_CARTESC_enddef( restart_fid ) ! [IN]
    endif

    return
  end subroutine LAND_vars_restart_enddef

  !-----------------------------------------------------------------------------
  !> Close restart file
  subroutine LAND_vars_restart_close
    use scale_file_cartesC, only: &
       FILE_CARTESC_close
    implicit none
    !---------------------------------------------------------------------------

    if ( restart_fid /= -1 ) then
       if( IO_L ) write(IO_FID_LOG,*)
       if( IO_L ) write(IO_FID_LOG,*) '*** Close restart file (LAND) ***'

       call FILE_CARTESC_close( restart_fid ) ! [IN]

       restart_fid = -1
    endif

    return
  end subroutine LAND_vars_restart_close

  !-----------------------------------------------------------------------------
  !> Define land variables in restart file
  subroutine LAND_vars_restart_def_var
    use scale_file_cartesC, only: &
       FILE_CARTESC_def_var
    implicit none
    integer :: i
    !---------------------------------------------------------------------------

    if ( restart_fid /= -1 ) then

       do i = I_TEMP, I_WATER
          if ( i == I_WATERDS ) cycle
          call FILE_CARTESC_def_var( restart_fid,     & ! [IN]
               VAR_NAME(i), VAR_DESC(i), VAR_UNIT(i), & ! [IN]
               'LXY', LAND_RESTART_OUT_DTYPE,         & ! [IN]
               VAR_ID(i)                              ) ! [OUT]
       end do
       do i = I_SFC_TEMP, VMAX
          call FILE_CARTESC_def_var( restart_fid,     & ! [IN]
               VAR_NAME(i), VAR_DESC(i), VAR_UNIT(i), & ! [IN]
               'XY', LAND_RESTART_OUT_DTYPE,          & ! [IN]
               VAR_ID(i)                              ) ! [OUT]
       end do

       !call FILE_CARTESC_def_var( restart_fid, 'SNOW_SFC_TEMP',  'Snow surface temperature',    &
       !                           'K',                   'XY',   LAND_RESTART_OUT_DTYPE, ?????  )
       !call FILE_CARTESC_def_var( restart_fid, 'SNOW_SWE',       'Snow water equivalent',       &
       !                           'kg/m2',               'XY',   LAND_RESTART_OUT_DTYPE, ?????  )
       !call FILE_CARTESC_def_var( restart_fid, 'SNOW_Depth',     'Snow depth',                  &
       !                           'm',                   'XY',   LAND_RESTART_OUT_DTYPE, ?????  )
       !call FILE_CARTESC_def_var( restart_fid, 'SNOW_Dzero',     'Snow depth at melting point', &
       !                           'm',                   'XY',   LAND_RESTART_OUT_DTYPE, ?????  )
       !call FILE_CARTESC_def_var( restart_fid, 'SNOW_nosnowsec', 'Time duration without snow', &
       !                           's',                 'XY',   LAND_RESTART_OUT_DTYPE,   ?????  )

    endif

    return
  end subroutine LAND_vars_restart_def_var

  !-----------------------------------------------------------------------------
  !> Write land variables to restart file
  subroutine LAND_vars_restart_write
    use scale_file_cartesC, only: &
       FILE_CARTESC_write_var
    implicit none
    !---------------------------------------------------------------------------

    if ( restart_fid /= -1 ) then

       call LAND_vars_total

       call FILE_CARTESC_write_var( restart_fid, VAR_ID(I_TEMP),      LAND_TEMP(:,:,:),          & ! [IN]
                                    VAR_NAME(I_TEMP),      'LXY',  fill_halo=.true.              ) ! [IN]
       call FILE_CARTESC_write_var( restart_fid, VAR_ID(I_WATER),     LAND_WATER(:,:,:),         & ! [IN]
                                    VAR_NAME(I_WATER),     'LXY',  fill_halo=.true.              ) ! [IN]
       call FILE_CARTESC_write_var( restart_fid, VAR_ID(I_SFC_TEMP),  LAND_SFC_TEMP(:,:),        & ! [IN]
                                    VAR_NAME(I_SFC_TEMP),  'XY',   fill_halo=.true.              ) ! [IN]
       call FILE_CARTESC_write_var( restart_fid, VAR_ID(I_ALB_LW),    LAND_SFC_albedo(:,:,I_LW), & ! [IN]
                                    VAR_NAME(I_ALB_LW),    'XY',   fill_halo=.true.              ) ! [IN]
       call FILE_CARTESC_write_var( restart_fid, VAR_ID(I_ALB_SW),    LAND_SFC_albedo(:,:,I_SW), & ! [IN]
                                    VAR_NAME(I_ALB_SW),    'XY',   fill_halo=.true.              ) ! [IN]
       call FILE_CARTESC_write_var( restart_fid, VAR_ID(I_SFLX_MW),   LAND_SFLX_MW(:,:),         & ! [IN]
                                    VAR_NAME(I_SFLX_MW),   'XY',   fill_halo=.true.              ) ! [IN]
       call FILE_CARTESC_write_var( restart_fid, VAR_ID(I_SFLX_MU),   LAND_SFLX_MU(:,:),         & ! [IN]
                                    VAR_NAME(I_SFLX_MU),   'XY',   fill_halo=.true.              ) ! [IN]
       call FILE_CARTESC_write_var( restart_fid, VAR_ID(I_SFLX_MV),   LAND_SFLX_MV(:,:),         & ! [IN]
                                    VAR_NAME(I_SFLX_MV),   'XY',   fill_halo=.true.              ) ! [IN]
       call FILE_CARTESC_write_var( restart_fid, VAR_ID(I_SFLX_SH),   LAND_SFLX_SH(:,:),         & ! [IN]
                                    VAR_NAME(I_SFLX_SH),   'XY',   fill_halo=.true.              ) ! [IN]
       call FILE_CARTESC_write_var( restart_fid, VAR_ID(I_SFLX_LH),   LAND_SFLX_LH(:,:),         & ! [IN]
                                    VAR_NAME(I_SFLX_LH),   'XY',   fill_halo=.true.              ) ! [IN]
       call FILE_CARTESC_write_var( restart_fid, VAR_ID(I_SFLX_GH),   LAND_SFLX_GH(:,:),         & ! [IN]
                                    VAR_NAME(I_SFLX_GH),   'XY',   fill_halo=.true.              ) ! [IN]
       call FILE_CARTESC_write_var( restart_fid, VAR_ID(I_SFLX_evap), LAND_SFLX_evap(:,:),       & ! [IN]
                                    VAR_NAME(I_SFLX_evap), 'XY',   fill_halo=.true.              ) ! [IN]

       !call FILE_CARTESC_write_var( restart_fid, VMAX+1, SNOW_SFC_TEMP(:,:),                    &
       !                           'SNOW_SFC_TEMP',       'XY',   fill_halo=.true.               )
       !call FILE_CARTESC_write_var( restart_fid, VMAX+2, SNOW_SWE(:,:),                         &
       !                           'SNOW_SWE',            'XY',   fill_halo=.true.               )
       !call FILE_CARTESC_write_var( restart_fid, VMAX+3, SNOW_Depth(:,:),                       &
       !                           'SNOW_Depth',          'XY',   fill_halo=.true.               )
       !call FILE_CARTESC_write_var( restart_fid, VMAX+4, SNOW_Dzero(:,:),                       &
       !                           'SNOW_Dzero',          'XY',   fill_halo=.true.               )
       !call FILE_CARTESC_write_var( restart_fid, VMAX+5, SNOW_nosnowsec(:,:),                   &
       !                           'SNOW_nosnowsec',      'XY',   fill_halo=.true.               )

    endif

    return
  end subroutine LAND_vars_restart_write

end module mod_land_vars
