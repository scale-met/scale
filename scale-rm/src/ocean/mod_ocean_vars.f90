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
  use scale_ocean_grid_cartesC_index

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

  public :: OCEAN_vars_restart_create
  public :: OCEAN_vars_restart_open
  public :: OCEAN_vars_restart_def_var
  public :: OCEAN_vars_restart_enddef
  public :: OCEAN_vars_restart_close

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  logical,               public :: OCEAN_RESTART_OUTPUT                 = .false.         !< Output restart file?

  character(len=H_LONG),  public :: OCEAN_RESTART_IN_BASENAME           = ''              !< Basename of the input  file
  logical,                public :: OCEAN_RESTART_IN_AGGREGATE                            !< Switch to use aggregate file
  logical,                public :: OCEAN_RESTART_IN_POSTFIX_TIMELABEL  = .false.         !< Add timelabel to the basename of input  file?
  character(len=H_LONG),  public :: OCEAN_RESTART_OUT_BASENAME          = ''              !< Basename of the output file
  logical,                public :: OCEAN_RESTART_OUT_AGGREGATE                           !< Switch to use aggregate file
  logical,                public :: OCEAN_RESTART_OUT_POSTFIX_TIMELABEL = .true.          !< Add timelabel to the basename of output file?
  character(len=H_MID),   public :: OCEAN_RESTART_OUT_TITLE             = 'OCEAN restart' !< Title    of the output file
  character(len=H_SHORT), public :: OCEAN_RESTART_OUT_DTYPE             = 'DEFAULT'       !< REAL4 or REAL8

  ! prognostic variables
  real(RP), public, allocatable :: OCEAN_TEMP      (:,:,:) !< ocean temperature [K]
  real(RP), public, allocatable :: OCEAN_SALT      (:,:,:) !< ocean salinity [PSU]
  real(RP), public, allocatable :: OCEAN_UVEL      (:,:,:) !< ocean zonal velocity [m/s]
  real(RP), public, allocatable :: OCEAN_VVEL      (:,:,:) !< ocean meridional velocity [m/s]

  real(RP), public, allocatable :: OCEAN_SFC_TEMP  (:,:)   !< ocean surface skin temperature [K]
  real(RP), public, allocatable :: OCEAN_SFC_albedo(:,:,:) !< ocean surface albedo (0-1)
  real(RP), public, allocatable :: OCEAN_SFC_Z0M   (:,:)   !< ocean surface roughness length for momentum [m]
  real(RP), public, allocatable :: OCEAN_SFC_Z0H   (:,:)   !< ocean surface roughness length for heat [m]
  real(RP), public, allocatable :: OCEAN_SFC_Z0E   (:,:)   !< ocean surface roughness length for vapor [m]

  ! tendency variables
  real(RP), public, allocatable :: OCEAN_TEMP_t      (:,:,:) !< tendency of ocean temperature [K/s]
  real(RP), public, allocatable :: OCEAN_SALT_t      (:,:,:) !< tendency of ocean salinity [PSU/s]
  real(RP), public, allocatable :: OCEAN_UVEL_t      (:,:,:) !< tendency of ocean zonal velocity [m/s2]
  real(RP), public, allocatable :: OCEAN_VVEL_t      (:,:,:) !< tendency of ocean meridional velocity [m/s2]

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
  real(RP), public, allocatable :: ATMOS_SFC_DENS (:,:)
  real(RP), public, allocatable :: ATMOS_SFC_PRES (:,:)
  real(RP), public, allocatable :: ATMOS_SFLX_LW  (:,:)
  real(RP), public, allocatable :: ATMOS_SFLX_SW  (:,:)
  real(RP), public, allocatable :: ATMOS_cosSZA   (:,:)
  real(RP), public, allocatable :: ATMOS_SFLX_prec(:,:)

  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  logical,                private :: OCEAN_VARS_CHECKRANGE      = .false.

  integer,                private, parameter :: VMAX        = 17 !< number of the variables 14-?27
  integer,                private, parameter :: I_TEMP      =  1
  integer,                private, parameter :: I_SALT      =  2
  integer,                private, parameter :: I_UVEL      =  3
  integer,                private, parameter :: I_VVEL      =  4
  integer,                private, parameter :: I_SFC_TEMP  =  5
  integer,                private, parameter :: I_ALB_LW    =  6
  integer,                private, parameter :: I_ALB_SW    =  7
  integer,                private, parameter :: I_SFC_Z0M   =  8
  integer,                private, parameter :: I_SFC_Z0H   =  9
  integer,                private, parameter :: I_SFC_Z0E   = 10
  integer,                private, parameter :: I_SFLX_MW   = 11
  integer,                private, parameter :: I_SFLX_MU   = 12
  integer,                private, parameter :: I_SFLX_MV   = 13
  integer,                private, parameter :: I_SFLX_SH   = 14
  integer,                private, parameter :: I_SFLX_LH   = 15
  integer,                private, parameter :: I_SFLX_WH   = 16
  integer,                private, parameter :: I_SFLX_evap = 17

  character(len=H_SHORT), private            :: VAR_NAME(VMAX) !< name  of the variables
  character(len=H_MID),   private            :: VAR_DESC(VMAX) !< desc. of the variables
  character(len=H_MID),   private            :: VAR_STDN(VMAX) !< standard name of the variables
  character(len=H_SHORT), private            :: VAR_UNIT(VMAX) !< unit  of the variables
  integer,                private            :: VAR_ID(VMAX)   !< ID    of the variables
  integer,                private            :: restart_fid = -1  ! file ID

  logical,                private            :: OCEAN_RESTART_IN_CHECK_COORDINATES = .true.

  data VAR_NAME / 'OCEAN_TEMP',      &
                  'OCEAN_SALT',      &
                  'OCEAN_UVEL',      &
                  'OCEAN_VVEL',      &
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
  data VAR_DESC / 'ocean temperature',                          &
                  'ocean salinity',                             &
                  'ocean u-velocity',                           &
                  'ocean v-velocity',                           &
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
  data VAR_STDN / 'sea_water_temperature', &
                  'sea_water_salinity', &
                  'eastward_sea_water_velocity', &
                  'northward_sea_water_velocity', &
                  'sea_surface_skin_temperature', &
                  '', &
                  '', &
                  '', &
                  '', &
                  '', &
                  '', &
                  '', &
                  '', &
                  '', &
                  '', &
                  '', &
                  '' /
  data VAR_UNIT / 'K',       &
                  'PSU',     &
                  'm/s',     &
                  'm/s',     &
                  'K',       &
                  '1',       &
                  '1',       &
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
       OCEAN_RESTART_IN_BASENAME,           &
       OCEAN_RESTART_IN_AGGREGATE,          &
       OCEAN_RESTART_IN_POSTFIX_TIMELABEL,  &
       OCEAN_RESTART_IN_CHECK_COORDINATES,  &
       OCEAN_RESTART_OUTPUT,                &
       OCEAN_RESTART_OUT_BASENAME,          &
       OCEAN_RESTART_OUT_AGGREGATE,         &
       OCEAN_RESTART_OUT_POSTFIX_TIMELABEL, &
       OCEAN_RESTART_OUT_TITLE,             &
       OCEAN_RESTART_OUT_DTYPE,             &
       OCEAN_VARS_CHECKRANGE

    integer :: ierr
    integer :: iv
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '++++++ Module[VARS] / Categ[OCEAN] / Origin[SCALE-RM]'

    allocate( OCEAN_TEMP(OKMAX,OIA,OJA) )
    allocate( OCEAN_SALT(OKMAX,OIA,OJA) )
    allocate( OCEAN_UVEL(OKMAX,OIA,OJA) )
    allocate( OCEAN_VVEL(OKMAX,OIA,OJA) )
    OCEAN_TEMP      (:,:,:) = UNDEF
    OCEAN_SALT      (:,:,:) = UNDEF
    OCEAN_UVEL      (:,:,:) = UNDEF
    OCEAN_VVEL      (:,:,:) = UNDEF

    allocate( OCEAN_SFC_TEMP  (OIA,OJA)   )
    allocate( OCEAN_SFC_albedo(OIA,OJA,2) )
    allocate( OCEAN_SFC_Z0M   (OIA,OJA)   )
    allocate( OCEAN_SFC_Z0H   (OIA,OJA)   )
    allocate( OCEAN_SFC_Z0E   (OIA,OJA)   )
    OCEAN_SFC_TEMP  (:,:)   = UNDEF
    OCEAN_SFC_albedo(:,:,:) = UNDEF
    OCEAN_SFC_Z0M   (:,:)   = UNDEF
    OCEAN_SFC_Z0H   (:,:)   = UNDEF
    OCEAN_SFC_Z0E   (:,:)   = UNDEF

    allocate( OCEAN_TEMP_t(OKMAX,OIA,OJA) )
    allocate( OCEAN_SALT_t(OKMAX,OIA,OJA) )
    allocate( OCEAN_UVEL_t(OKMAX,OIA,OJA) )
    allocate( OCEAN_VVEL_t(OKMAX,OIA,OJA) )
    OCEAN_TEMP_t(:,:,:) = UNDEF
    OCEAN_SALT_t(:,:,:) = UNDEF
    OCEAN_UVEL_t(:,:,:) = UNDEF
    OCEAN_VVEL_t(:,:,:) = UNDEF

    allocate( OCEAN_SFC_TEMP_t  (OIA,OJA)   )
    allocate( OCEAN_SFC_albedo_t(OIA,OJA,2) )
    allocate( OCEAN_SFC_Z0M_t   (OIA,OJA)   )
    allocate( OCEAN_SFC_Z0H_t   (OIA,OJA)   )
    allocate( OCEAN_SFC_Z0E_t   (OIA,OJA)   )
    OCEAN_SFC_TEMP_t  (:,:)   = UNDEF
    OCEAN_SFC_albedo_t(:,:,:) = UNDEF
    OCEAN_SFC_Z0M_t   (:,:)   = UNDEF
    OCEAN_SFC_Z0H_t   (:,:)   = UNDEF
    OCEAN_SFC_Z0E_t   (:,:)   = UNDEF

    allocate( OCEAN_SFLX_MW  (OIA,OJA) )
    allocate( OCEAN_SFLX_MU  (OIA,OJA) )
    allocate( OCEAN_SFLX_MV  (OIA,OJA) )
    allocate( OCEAN_SFLX_SH  (OIA,OJA) )
    allocate( OCEAN_SFLX_LH  (OIA,OJA) )
    allocate( OCEAN_SFLX_WH  (OIA,OJA) )
    allocate( OCEAN_SFLX_evap(OIA,OJA) )
    OCEAN_SFLX_MW  (:,:) = UNDEF
    OCEAN_SFLX_MU  (:,:) = UNDEF
    OCEAN_SFLX_MV  (:,:) = UNDEF
    OCEAN_SFLX_SH  (:,:) = UNDEF
    OCEAN_SFLX_LH  (:,:) = UNDEF
    OCEAN_SFLX_WH  (:,:) = UNDEF
    OCEAN_SFLX_evap(:,:) = UNDEF

    allocate( OCEAN_U10(OIA,OJA) )
    allocate( OCEAN_V10(OIA,OJA) )
    allocate( OCEAN_T2 (OIA,OJA) )
    allocate( OCEAN_Q2 (OIA,OJA) )
    OCEAN_U10(:,:) = UNDEF
    OCEAN_V10(:,:) = UNDEF
    OCEAN_T2 (:,:) = UNDEF
    OCEAN_Q2 (:,:) = UNDEF

    allocate( ATMOS_TEMP     (OIA,OJA) )
    allocate( ATMOS_PRES     (OIA,OJA) )
    allocate( ATMOS_W        (OIA,OJA) )
    allocate( ATMOS_U        (OIA,OJA) )
    allocate( ATMOS_V        (OIA,OJA) )
    allocate( ATMOS_DENS     (OIA,OJA) )
    allocate( ATMOS_QV       (OIA,OJA) )
    allocate( ATMOS_PBL      (OIA,OJA) )
    allocate( ATMOS_SFC_DENS (OIA,OJA) )
    allocate( ATMOS_SFC_PRES (OIA,OJA) )
    allocate( ATMOS_SFLX_LW  (OIA,OJA) )
    allocate( ATMOS_SFLX_SW  (OIA,OJA) )
    allocate( ATMOS_cosSZA   (OIA,OJA) )
    allocate( ATMOS_SFLX_prec(OIA,OJA) )
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

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_OCEAN_VARS,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_OCEAN_VARS. Check!'
       call PRC_MPIstop
    endif
    if( IO_NML ) write(IO_FID_NML,nml=PARAM_OCEAN_VARS)

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** List of prognostic variables (OCEAN) ***'
    if( IO_L ) write(IO_FID_LOG,'(1x,A,A24,A,A48,A,A12,A)') &
               '***       |', 'VARNAME                 ','|', &
               'DESCRIPTION                                     ', '[', 'UNIT        ', ']'
    do iv = 1, VMAX
       if( IO_L ) write(IO_FID_LOG,'(1x,A,I3,A,A24,A,A48,A,A12,A)') &
                  '*** NO.',iv,'|',VAR_NAME(iv),'|',VAR_DESC(iv),'[',VAR_UNIT(iv),']'

    enddo

    if( IO_L ) write(IO_FID_LOG,*)
    if ( OCEAN_RESTART_IN_BASENAME /= '' ) then
       if( IO_L ) write(IO_FID_LOG,*) '*** Restart input?  : YES, file = ', trim(OCEAN_RESTART_IN_BASENAME)
       if( IO_L ) write(IO_FID_LOG,*) '*** Add timelabel?  : ', OCEAN_RESTART_IN_POSTFIX_TIMELABEL
    else
       if( IO_L ) write(IO_FID_LOG,*) '*** Restart input?  : NO'
    endif
    if (       OCEAN_RESTART_OUTPUT             &
         .AND. OCEAN_RESTART_OUT_BASENAME /= '' ) then
       if( IO_L ) write(IO_FID_LOG,*) '*** Restart output? : YES, file = ', trim(OCEAN_RESTART_OUT_BASENAME)
       if( IO_L ) write(IO_FID_LOG,*) '*** Add timelabel?  : ', OCEAN_RESTART_OUT_POSTFIX_TIMELABEL
    else
       if( IO_L ) write(IO_FID_LOG,*) '*** Restart output? : NO'
       OCEAN_RESTART_OUTPUT = .false.
    endif

    return
  end subroutine OCEAN_vars_setup

  !-----------------------------------------------------------------------------
  !> Open ocean restart file for read
  subroutine OCEAN_vars_restart_open
    use scale_time, only: &
       TIME_gettimelabel
    use scale_file_cartesC, only: &
       FILE_CARTESC_open, &
       FILE_CARTESC_check_coordinates
    use mod_ocean_admin, only: &
       OCEAN_sw
    implicit none

    character(len=19)     :: timelabel
    character(len=H_LONG) :: basename
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** Open restart file (OCEAN) ***'

    if ( OCEAN_sw .and. OCEAN_RESTART_IN_BASENAME /= '' ) then

       if ( OCEAN_RESTART_IN_POSTFIX_TIMELABEL ) then
          call TIME_gettimelabel( timelabel )
          basename = trim(OCEAN_RESTART_IN_BASENAME)//'_'//trim(timelabel)
       else
          basename = trim(OCEAN_RESTART_IN_BASENAME)
       endif

       if( IO_L ) write(IO_FID_LOG,*) '*** basename: ', trim(basename)

       call FILE_CARTESC_open( basename, restart_fid, aggregate=OCEAN_RESTART_IN_AGGREGATE )

       if ( OCEAN_RESTART_IN_CHECK_COORDINATES ) then
          call FILE_CARTESC_check_coordinates( restart_fid )
       end if

    else
       if( IO_L ) write(IO_FID_LOG,*) '*** restart file for ocean is not specified.'
    endif

    return
  end subroutine OCEAN_vars_restart_open

  !-----------------------------------------------------------------------------
  !> Read ocean restart
  subroutine OCEAN_vars_restart_read
    use scale_file, only: &
       FILE_get_aggregate
    use scale_file_cartesC, only: &
       FILE_CARTESC_read, &
       FILE_CARTESC_flush
    implicit none
    !---------------------------------------------------------------------------

    if ( restart_fid /= -1 ) then
       if( IO_L ) write(IO_FID_LOG,*)
       if( IO_L ) write(IO_FID_LOG,*) '*** Read from restart file (OCEAN) ***'

       call FILE_CARTESC_read( restart_fid, VAR_NAME(I_TEMP),      'OXY', & ! [IN]
                               OCEAN_TEMP(:,:,:)                          ) ! [OUT]
!       call FILE_CARTESC_read( restart_fid, VAR_NAME(I_SALT),      'OXY', & ! [IN]
!                               OCEAN_SALT(:,:,:)                          ) ! [OUT]
!       call FILE_CARTESC_read( restart_fid, VAR_NAME(I_UVEL),      'OXY', & ! [IN]
!                               OCEAN_UVEL(:,:,:)                          ) ! [OUT]
!       call FILE_CARTESC_read( restart_fid, VAR_NAME(I_VVEL),      'OXY', & ! [IN]
!                               OCEAN_VVEL(:,:,:)                          ) ! [OUT]
!                         
       call FILE_CARTESC_read( restart_fid, VAR_NAME(I_SFC_TEMP),  'XY', & ! [IN]
                               OCEAN_SFC_TEMP(:,:)                       ) ! [OUT]
            
       call FILE_CARTESC_read( restart_fid, VAR_NAME(I_ALB_LW),    'XY', & ! [IN]
                               OCEAN_SFC_albedo(:,:,I_LW)                ) ! [OUT]
       call FILE_CARTESC_read( restart_fid, VAR_NAME(I_ALB_SW),    'XY', & ! [IN]
                               OCEAN_SFC_albedo(:,:,I_SW)                ) ! [OUT]
       call FILE_CARTESC_read( restart_fid, VAR_NAME(I_SFC_Z0M),   'XY', & ! [IN]
                               OCEAN_SFC_Z0M(:,:)                        ) ! [OUT]
       call FILE_CARTESC_read( restart_fid, VAR_NAME(I_SFC_Z0H),   'XY', & ! [IN]
                               OCEAN_SFC_Z0H(:,:)                        ) ! [OUT]
       call FILE_CARTESC_read( restart_fid, VAR_NAME(I_SFC_Z0E),   'XY', & ! [IN]
                               OCEAN_SFC_Z0E(:,:)                        ) ! [OUT]

       call FILE_CARTESC_read( restart_fid, VAR_NAME(I_SFLX_MW),   'XY', & ! [IN]
                               OCEAN_SFLX_MW(:,:)                        ) ! [OUT]
       call FILE_CARTESC_read( restart_fid, VAR_NAME(I_SFLX_MU),   'XY', & ! [IN]
                               OCEAN_SFLX_MU(:,:)                        ) ! [OUT]
       call FILE_CARTESC_read( restart_fid, VAR_NAME(I_SFLX_MV),   'XY', & ! [IN]
                               OCEAN_SFLX_MV(:,:)                        ) ! [OUT]
       call FILE_CARTESC_read( restart_fid, VAR_NAME(I_SFLX_SH),   'XY', & ! [IN]
                               OCEAN_SFLX_SH(:,:)                        ) ! [OUT]
       call FILE_CARTESC_read( restart_fid, VAR_NAME(I_SFLX_LH),   'XY', & ! [IN]
                               OCEAN_SFLX_LH(:,:)                        ) ! [OUT]
       call FILE_CARTESC_read( restart_fid, VAR_NAME(I_SFLX_WH),   'XY', & ! [IN]
                               OCEAN_SFLX_WH(:,:)                        ) ! [OUT]
       call FILE_CARTESC_read( restart_fid, VAR_NAME(I_SFLX_evap), 'XY', & ! [IN]
                               OCEAN_SFLX_evap(:,:)                      ) ! [OUT]
            
       if( FILE_get_AGGREGATE(restart_fid) ) call FILE_CARTESC_flush( restart_fid ) ! commit all pending read requests

       call OCEAN_vars_total
    else
       if( IO_L ) write(IO_FID_LOG,*) '*** invalid restart file ID for ocean.'
    endif

    return
  end subroutine OCEAN_vars_restart_read

  !-----------------------------------------------------------------------------
  !> History output set for ocean variables
  subroutine OCEAN_vars_history
    use scale_file_history, only: &
       FILE_HISTORY_in
    implicit none
    !---------------------------------------------------------------------------

    if ( OCEAN_VARS_CHECKRANGE ) then
       call VALCHECK( OCEAN_TEMP      (OKS:OKE,OIS:OIE,OJS:OJE), 0.0_RP, 1000.0_RP, VAR_NAME(I_TEMP), &
                     __FILE__, __LINE__ )
!        call VALCHECK( OCEAN_SALT      (OKS:OKE,OIS:OIE,OJS:OJE), 0.0_RP, 1000.0_RP, VAR_NAME(I_SALT), &
!                      __FILE__, __LINE__ )
!        call VALCHECK( OCEAN_UVEL      (OKS:OKE,OIS:OIE,OJS:OJE), 0.0_RP, 1000.0_RP, VAR_NAME(I_UVEL), &
!                      __FILE__, __LINE__ )
!        call VALCHECK( OCEAN_VVEL      (OKS:OKE,OIS:OIE,OJS:OJE), 0.0_RP, 1000.0_RP, VAR_NAME(I_VVEL), &
!                      __FILE__, __LINE__ )

       call VALCHECK( OCEAN_SFC_TEMP  (OIS:OIE,OJS:OJE),      0.0_RP, 1000.0_RP, VAR_NAME(I_SFC_TEMP), &
                     __FILE__, __LINE__ )
       call VALCHECK( OCEAN_SFC_albedo(OIS:OIE,OJS:OJE,I_LW), 0.0_RP,    2.0_RP, VAR_NAME(I_ALB_LW),   &
                     __FILE__, __LINE__ )
       call VALCHECK( OCEAN_SFC_albedo(OIS:OIE,OJS:OJE,I_SW), 0.0_RP,    2.0_RP, VAR_NAME(I_ALB_SW),   &
                     __FILE__, __LINE__ )
       call VALCHECK( OCEAN_SFC_Z0M   (OIS:OIE,OJS:OJE),      0.0_RP, 1000.0_RP, VAR_NAME(I_SFC_Z0M),  &
                     __FILE__, __LINE__ )
       call VALCHECK( OCEAN_SFC_Z0H   (OIS:OIE,OJS:OJE),      0.0_RP, 1000.0_RP, VAR_NAME(I_SFC_Z0H),  &
                     __FILE__, __LINE__ )
       call VALCHECK( OCEAN_SFC_Z0E   (OIS:OIE,OJS:OJE),      0.0_RP, 1000.0_RP, VAR_NAME(I_SFC_Z0E),  &
                     __FILE__, __LINE__ )
    endif

    call FILE_HISTORY_in( OCEAN_TEMP      (:,:,:),    VAR_NAME(I_TEMP),     VAR_DESC(I_TEMP),     VAR_UNIT(I_TEMP)    , dim_type="OXY", standard_name=VAR_STDN(I_TEMP) )
    call FILE_HISTORY_in( OCEAN_SALT      (:,:,:),    VAR_NAME(I_SALT),     VAR_DESC(I_SALT),     VAR_UNIT(I_SALT)    , dim_type="OXY", standard_name=VAR_STDN(I_SALT) )
    call FILE_HISTORY_in( OCEAN_UVEL      (:,:,:),    VAR_NAME(I_UVEL),     VAR_DESC(I_UVEL),     VAR_UNIT(I_UVEL)    , dim_type="OXY", standard_name=VAR_STDN(I_UVEL) )
    call FILE_HISTORY_in( OCEAN_VVEL      (:,:,:),    VAR_NAME(I_VVEL),     VAR_DESC(I_VVEL),     VAR_UNIT(I_VVEL)    , dim_type="OXY", standard_name=VAR_STDN(I_VVEL) )

    call FILE_HISTORY_in( OCEAN_SFC_TEMP  (:,:),      VAR_NAME(I_SFC_TEMP), VAR_DESC(I_SFC_TEMP), VAR_UNIT(I_SFC_TEMP), standard_name=VAR_STDN(I_SFC_TEMP) )
    call FILE_HISTORY_in( OCEAN_SFC_albedo(:,:,I_LW), VAR_NAME(I_ALB_LW),   VAR_DESC(I_ALB_LW),   VAR_UNIT(I_ALB_LW),  standard_name=VAR_STDN(I_ALB_LW) )
    call FILE_HISTORY_in( OCEAN_SFC_albedo(:,:,I_SW), VAR_NAME(I_ALB_SW),   VAR_DESC(I_ALB_SW),   VAR_UNIT(I_ALB_SW),  standard_name=VAR_STDN(I_ALB_LW) )
    call FILE_HISTORY_in( OCEAN_SFC_Z0M   (:,:),      VAR_NAME(I_SFC_Z0M),  VAR_DESC(I_SFC_Z0M),  VAR_UNIT(I_SFC_Z0M), standard_name=VAR_STDN(I_SFC_Z0M) )
    call FILE_HISTORY_in( OCEAN_SFC_Z0H   (:,:),      VAR_NAME(I_SFC_Z0H),  VAR_DESC(I_SFC_Z0H),  VAR_UNIT(I_SFC_Z0H), standard_name=VAR_STDN(I_SFC_Z0H) )
    call FILE_HISTORY_in( OCEAN_SFC_Z0E   (:,:),      VAR_NAME(I_SFC_Z0E),  VAR_DESC(I_SFC_Z0E),  VAR_UNIT(I_SFC_Z0E), standard_name=VAR_STDN(I_SFC_Z0H) )

    call FILE_HISTORY_in( OCEAN_SFLX_MW  (:,:), VAR_NAME(I_SFLX_MW),   VAR_DESC(I_SFLX_MW),   VAR_UNIT(I_SFLX_MW),     standard_name=VAR_STDN(I_SFLX_MW) )
    call FILE_HISTORY_in( OCEAN_SFLX_MU  (:,:), VAR_NAME(I_SFLX_MU),   VAR_DESC(I_SFLX_MU),   VAR_UNIT(I_SFLX_MU),     standard_name=VAR_STDN(I_SFLX_MU) )
    call FILE_HISTORY_in( OCEAN_SFLX_MV  (:,:), VAR_NAME(I_SFLX_MV),   VAR_DESC(I_SFLX_MV),   VAR_UNIT(I_SFLX_MV),     standard_name=VAR_STDN(I_SFLX_MV) )
    call FILE_HISTORY_in( OCEAN_SFLX_SH  (:,:), VAR_NAME(I_SFLX_SH),   VAR_DESC(I_SFLX_SH),   VAR_UNIT(I_SFLX_SH),     standard_name=VAR_STDN(I_SFLX_SH) )
    call FILE_HISTORY_in( OCEAN_SFLX_LH  (:,:), VAR_NAME(I_SFLX_LH),   VAR_DESC(I_SFLX_LH),   VAR_UNIT(I_SFLX_LH),     standard_name=VAR_STDN(I_SFLX_LH) )
    call FILE_HISTORY_in( OCEAN_SFLX_WH  (:,:), VAR_NAME(I_SFLX_WH),   VAR_DESC(I_SFLX_WH),   VAR_UNIT(I_SFLX_WH),     standard_name=VAR_STDN(I_SFLX_WH) )
    call FILE_HISTORY_in( OCEAN_SFLX_evap(:,:), VAR_NAME(I_SFLX_evap), VAR_DESC(I_SFLX_evap), VAR_UNIT(I_SFLX_evap),   standard_name=VAR_STDN(I_SFLX_evap) )

    return
  end subroutine OCEAN_vars_history

  !-----------------------------------------------------------------------------
  !> Budget monitor for ocean
  subroutine OCEAN_vars_total
    use scale_statistics, only: &
       STATISTICS_checktotal, &
       STATISTICS_total
    use scale_ocean_grid_cartesC_real, only: &
       OCEAN_GRID_CARTESC_REAL_AREA,    &
       OCEAN_GRID_CARTESC_REAL_TOTAREA, &
       OCEAN_GRID_CARTESC_REAL_VOL,     &
       OCEAN_GRID_CARTESC_REAL_TOTVOL
    implicit none

    !---------------------------------------------------------------------------

    if ( STATISTICS_checktotal ) then

       ! 3D
       call STATISTICS_total( OKA, OKS, OKE, OIA, OIS, OIE, OJA, OJS, OJE, &
                              OCEAN_TEMP(:,:,:), VAR_NAME(I_TEMP), & ! (in)
                              OCEAN_GRID_CARTESC_REAL_VOL(:,:,:),  & ! (in)
                              OCEAN_GRID_CARTESC_REAL_TOTVOL       ) ! (in)
       call STATISTICS_total( OKA, OKS, OKE, OIA, OIS, OIE, OJA, OJS, OJE, &
                              OCEAN_SALT(:,:,:), VAR_NAME(I_SALT), & ! (in)
                              OCEAN_GRID_CARTESC_REAL_VOL(:,:,:),  & ! (in)
                              OCEAN_GRID_CARTESC_REAL_TOTVOL       ) ! (in)
       call STATISTICS_total( OKA, OKS, OKE, OIA, OIS, OIE, OJA, OJS, OJE, &
                              OCEAN_UVEL(:,:,:), VAR_NAME(I_UVEL), & ! (in)
                              OCEAN_GRID_CARTESC_REAL_VOL(:,:,:),  & ! (in)
                              OCEAN_GRID_CARTESC_REAL_TOTVOL       ) ! (in)
       call STATISTICS_total( OKA, OKS, OKE, OIA, OIS, OIE, OJA, OJS, OJE, &
                              OCEAN_VVEL(:,:,:), VAR_NAME(I_VVEL), & ! (in)
                              OCEAN_GRID_CARTESC_REAL_VOL(:,:,:),  & ! (in)
                              OCEAN_GRID_CARTESC_REAL_TOTVOL       ) ! (in)

       ! 2D
       call STATISTICS_total( OIA, OIS, OIE, OJA, OJS, OJE, &
                              OCEAN_SFC_TEMP  (:,:),      VAR_NAME(I_SFC_TEMP), & ! (in)
                              OCEAN_GRID_CARTESC_REAL_AREA(:,:),                & ! (in)
                              OCEAN_GRID_CARTESC_REAL_TOTAREA                   ) ! (in)
       call STATISTICS_total( OIA, OIS, OIE, OJA, OJS, OJE, &
                              OCEAN_SFC_albedo(:,:,I_LW), VAR_NAME(I_ALB_LW),   & ! (in)
                              OCEAN_GRID_CARTESC_REAL_AREA(:,:),                & ! (in)
                              OCEAN_GRID_CARTESC_REAL_TOTAREA                   ) ! (in)
       call STATISTICS_total( OIA, OIS, OIE, OJA, OJS, OJE, &
                              OCEAN_SFC_albedo(:,:,I_SW), VAR_NAME(I_ALB_SW),   & ! (in)
                              OCEAN_GRID_CARTESC_REAL_AREA(:,:),                & ! (in)
                              OCEAN_GRID_CARTESC_REAL_TOTAREA                   ) ! (in)
       call STATISTICS_total( OIA, OIS, OIE, OJA, OJS, OJE, &
                              OCEAN_SFC_Z0M   (:,:),      VAR_NAME(I_SFC_Z0M),  & ! (in)
                              OCEAN_GRID_CARTESC_REAL_AREA(:,:),                & ! (in)
                              OCEAN_GRID_CARTESC_REAL_TOTAREA                   ) ! (in)
       call STATISTICS_total( OIA, OIS, OIE, OJA, OJS, OJE, &
                              OCEAN_SFC_Z0H   (:,:),      VAR_NAME(I_SFC_Z0H),  & ! (in)
                              OCEAN_GRID_CARTESC_REAL_AREA(:,:),                & ! (in)
                              OCEAN_GRID_CARTESC_REAL_TOTAREA                   ) ! (in)
       call STATISTICS_total( OIA, OIS, OIE, OJA, OJS, OJE, &
                              OCEAN_SFC_Z0E   (:,:),      VAR_NAME(I_SFC_Z0E),  & ! (in)
                              OCEAN_GRID_CARTESC_REAL_AREA(:,:),                & ! (in)
                              OCEAN_GRID_CARTESC_REAL_TOTAREA                   ) ! (in)

       call STATISTICS_total( OIA, OIS, OIE, OJA, OJS, OJE, &
                              OCEAN_SFLX_MW  (:,:), VAR_NAME(I_SFLX_MW),   & ! (in) 
                              OCEAN_GRID_CARTESC_REAL_AREA(:,:),           & ! (in)
                              OCEAN_GRID_CARTESC_REAL_TOTAREA              ) ! (in)
       call STATISTICS_total( OIA, OIS, OIE, OJA, OJS, OJE, &
                              OCEAN_SFLX_MU  (:,:), VAR_NAME(I_SFLX_MU),   & ! (in)
                              OCEAN_GRID_CARTESC_REAL_AREA(:,:),           & ! (in)
                              OCEAN_GRID_CARTESC_REAL_TOTAREA              ) ! (in)
       call STATISTICS_total( OIA, OIS, OIE, OJA, OJS, OJE, &
                              OCEAN_SFLX_MV  (:,:), VAR_NAME(I_SFLX_MV),   & ! (in)
                              OCEAN_GRID_CARTESC_REAL_AREA(:,:),           & ! (in)
                              OCEAN_GRID_CARTESC_REAL_TOTAREA              ) ! (in)
       call STATISTICS_total( OIA, OIS, OIE, OJA, OJS, OJE, &
                              OCEAN_SFLX_SH  (:,:), VAR_NAME(I_SFLX_SH),   & ! (in)
                              OCEAN_GRID_CARTESC_REAL_AREA(:,:),           & ! (in)
                              OCEAN_GRID_CARTESC_REAL_TOTAREA              ) ! (in)
       call STATISTICS_total( OIA, OIS, OIE, OJA, OJS, OJE, &
                              OCEAN_SFLX_LH  (:,:), VAR_NAME(I_SFLX_LH),   & ! (in)
                              OCEAN_GRID_CARTESC_REAL_AREA(:,:),           & ! (in)
                              OCEAN_GRID_CARTESC_REAL_TOTAREA              ) ! (in)
       call STATISTICS_total( OIA, OIS, OIE, OJA, OJS, OJE, &
                              OCEAN_SFLX_WH  (:,:), VAR_NAME(I_SFLX_WH),   & ! (in)
                              OCEAN_GRID_CARTESC_REAL_AREA(:,:),           & ! (in)
                              OCEAN_GRID_CARTESC_REAL_TOTAREA              ) ! (in)
       call STATISTICS_total( OIA, OIS, OIE, OJA, OJS, OJE, &
                              OCEAN_SFLX_evap(:,:), VAR_NAME(I_SFLX_evap), & ! (in)
                              OCEAN_GRID_CARTESC_REAL_AREA(:,:),           & ! (in)
                              OCEAN_GRID_CARTESC_REAL_TOTAREA              ) ! (in)

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

    real(RP), intent(in) :: OCEAN_TEMP_in      (OKMAX,OIA,OJA)
    real(RP), intent(in) :: OCEAN_SFC_TEMP_in  (OIA,OJA)
    real(RP), intent(in) :: OCEAN_SFC_albedo_in(OIA,OJA,2)
    real(RP), intent(in) :: OCEAN_SFC_Z0M_in   (OIA,OJA)
    real(RP), intent(in) :: OCEAN_SFC_Z0H_in   (OIA,OJA)
    real(RP), intent(in) :: OCEAN_SFC_Z0E_in   (OIA,OJA)
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** External Input file (ocean) ***'

    OCEAN_TEMP      (:,:,:) = OCEAN_TEMP_in      (:,:,:)

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
  subroutine OCEAN_vars_restart_create
    use scale_time, only: &
       TIME_gettimelabel
    use scale_file_cartesC, only: &
       FILE_CARTESC_create
    use mod_ocean_admin, only: &
       OCEAN_sw
    implicit none

    character(len=19)     :: timelabel
    character(len=H_LONG) :: basename
    !---------------------------------------------------------------------------

    if ( OCEAN_sw .and. OCEAN_RESTART_OUT_BASENAME /= '' ) then

       if( IO_L ) write(IO_FID_LOG,*)
       if( IO_L ) write(IO_FID_LOG,*) '*** Create restart file (OCEAN) ***'

       if ( OCEAN_RESTART_OUT_POSTFIX_TIMELABEL ) then
          call TIME_gettimelabel( timelabel )
          basename = trim(OCEAN_RESTART_OUT_BASENAME)//'_'//trim(timelabel)
       else
          basename = trim(OCEAN_RESTART_OUT_BASENAME)
       endif

       if( IO_L ) write(IO_FID_LOG,*) '*** basename: ', trim(basename)

       call FILE_CARTESC_create( &
            basename, OCEAN_RESTART_OUT_TITLE, OCEAN_RESTART_OUT_DTYPE, & ! [IN]
            restart_fid,                                                & ! [OUT]
            aggregate=OCEAN_RESTART_OUT_AGGREGATE                       ) ! [IN]

    endif

    return
  end subroutine OCEAN_vars_restart_create

  !-----------------------------------------------------------------------------
  !> Exit netCDF define mode
  subroutine OCEAN_vars_restart_enddef
    use scale_file_cartesC, only: &
       FILE_CARTESC_enddef
    implicit none

    if ( restart_fid /= -1 ) then
       call FILE_CARTESC_enddef( restart_fid ) ! [IN]
    endif

    return
  end subroutine OCEAN_vars_restart_enddef

  !-----------------------------------------------------------------------------
  !> Close restart file
  subroutine OCEAN_vars_restart_close
    use scale_file_cartesC, only: &
       FILE_CARTESC_close
    implicit none
    !---------------------------------------------------------------------------

    if ( restart_fid /= -1 ) then
       if( IO_L ) write(IO_FID_LOG,*)
       if( IO_L ) write(IO_FID_LOG,*) '*** Close restart file (OCEAN) ***'

       call FILE_CARTESC_close( restart_fid ) ! [IN]

       restart_fid = -1
    endif

    return
  end subroutine OCEAN_vars_restart_close

  !-----------------------------------------------------------------------------
  !> Define ocean variables in restart file
  subroutine OCEAN_vars_restart_def_var
    use scale_file_cartesC, only: &
       FILE_CARTESC_def_var
    implicit none
    !---------------------------------------------------------------------------

    if ( restart_fid /= -1 ) then

       call FILE_CARTESC_def_var( restart_fid, VAR_NAME(I_TEMP),      VAR_DESC(I_TEMP),      VAR_UNIT(I_TEMP),      'OXY', OCEAN_RESTART_OUT_DTYPE, &
                                  VAR_ID(I_TEMP), &
                                  standard_name=VAR_STDN(I_TEMP) )
!       call FILE_CARTESC_def_var( restart_fid, VAR_NAME(I_SALT),      VAR_DESC(I_SALT),      VAR_UNIT(I_SALT),      'OXY', OCEAN_RESTART_OUT_DTYPE, &
!                                  VAR_ID(I_SALT), &
!                                  standard_name=VAR_STDN(I_SALT) )
!       call FILE_CARTESC_def_var( restart_fid, VAR_NAME(I_UVEL),      VAR_DESC(I_UVEL),      VAR_UNIT(I_UVEL),      'OXY', OCEAN_RESTART_OUT_DTYPE, &
!                                  VAR_ID(I_UVEL), &
!                                  standard_name=VAR_STDN(I_UVEL) )
!       call FILE_CARTESC_def_var( restart_fid, VAR_NAME(I_VVEL),      VAR_DESC(I_VVEL),      VAR_UNIT(I_VVEL),      'OXY', OCEAN_RESTART_OUT_DTYPE, &
!                                  VAR_ID(I_VVEL), &
!                                  standard_name=VAR_STDN(I_VVEL) )
       call FILE_CARTESC_def_var( restart_fid, VAR_NAME(I_SFC_TEMP),  VAR_DESC(I_SFC_TEMP),  VAR_UNIT(I_SFC_TEMP),  'XY', OCEAN_RESTART_OUT_DTYPE, &
                                  VAR_ID(I_SFC_TEMP), &
                                  standard_name=VAR_STDN(I_SFC_TEMP) )
       call FILE_CARTESC_def_var( restart_fid, VAR_NAME(I_ALB_LW),    VAR_DESC(I_ALB_LW),    VAR_UNIT(I_ALB_LW),    'XY', OCEAN_RESTART_OUT_DTYPE, &
                                  VAR_ID(I_ALB_LW), &
                                  standard_name=VAR_STDN(I_ALB_LW) )
       call FILE_CARTESC_def_var( restart_fid, VAR_NAME(I_ALB_SW),    VAR_DESC(I_ALB_SW),    VAR_UNIT(I_ALB_SW),    'XY', OCEAN_RESTART_OUT_DTYPE, &
                                  VAR_ID(I_ALB_SW), &
                                  standard_name=VAR_STDN(I_ALB_SW) )
       call FILE_CARTESC_def_var( restart_fid, VAR_NAME(I_SFC_Z0M),   VAR_DESC(I_SFC_Z0M),   VAR_UNIT(I_SFC_Z0M),   'XY', OCEAN_RESTART_OUT_DTYPE, &
                                  VAR_ID(I_SFC_Z0M), &
                                  standard_name=VAR_STDN(I_SFC_Z0M) )
       call FILE_CARTESC_def_var( restart_fid, VAR_NAME(I_SFC_Z0H),   VAR_DESC(I_SFC_Z0H),   VAR_UNIT(I_SFC_Z0H),   'XY', OCEAN_RESTART_OUT_DTYPE, &
                                  VAR_ID(I_SFC_Z0H), &
                                  standard_name=VAR_STDN(I_SFC_Z0H) )
       call FILE_CARTESC_def_var( restart_fid, VAR_NAME(I_SFC_Z0E),   VAR_DESC(I_SFC_Z0E),   VAR_UNIT(I_SFC_Z0E),   'XY', OCEAN_RESTART_OUT_DTYPE, &
                                  VAR_ID(I_SFC_Z0E), &
                                  standard_name=VAR_STDN(I_SFC_Z0E) )
       call FILE_CARTESC_def_var( restart_fid, VAR_NAME(I_SFLX_MW),   VAR_DESC(I_SFLX_MW),   VAR_UNIT(I_SFLX_MW),   'XY', OCEAN_RESTART_OUT_DTYPE, &
                                  VAR_ID(I_SFLX_MW), &
                                  standard_name=VAR_STDN(I_SFLX_MW) )
       call FILE_CARTESC_def_var( restart_fid, VAR_NAME(I_SFLX_MU),   VAR_DESC(I_SFLX_MU),   VAR_UNIT(I_SFLX_MU),  'XY', OCEAN_RESTART_OUT_DTYPE, &
                                  VAR_ID(I_SFLX_MU), &
                                  standard_name=VAR_STDN(I_SFLX_MU) )
       call FILE_CARTESC_def_var( restart_fid, VAR_NAME(I_SFLX_MV),   VAR_DESC(I_SFLX_MV),   VAR_UNIT(I_SFLX_MV),   'XY', OCEAN_RESTART_OUT_DTYPE, &
                                  VAR_ID(I_SFLX_MV), &
                                  standard_name=VAR_STDN(I_SFLX_MV) )
       call FILE_CARTESC_def_var( restart_fid, VAR_NAME(I_SFLX_SH),   VAR_DESC(I_SFLX_SH),   VAR_UNIT(I_SFLX_SH),   'XY', OCEAN_RESTART_OUT_DTYPE, &
                                  VAR_ID(I_SFLX_SH), &
                                  standard_name=VAR_STDN(I_SFLX_SH) )
       call FILE_CARTESC_def_var( restart_fid, VAR_NAME(I_SFLX_LH),   VAR_DESC(I_SFLX_LH),   VAR_UNIT(I_SFLX_LH),   'XY', OCEAN_RESTART_OUT_DTYPE, &
                                  VAR_ID(I_SFLX_LH), &
                                  standard_name=VAR_STDN(I_SFLX_LH) )
       call FILE_CARTESC_def_var( restart_fid, VAR_NAME(I_SFLX_WH),   VAR_DESC(I_SFLX_WH),   VAR_UNIT(I_SFLX_WH),   'XY', OCEAN_RESTART_OUT_DTYPE, &
                                  VAR_ID(I_SFLX_WH), &
                                  standard_name=VAR_STDN(I_SFLX_WH) )
       call FILE_CARTESC_def_var( restart_fid, VAR_NAME(I_SFLX_evap), VAR_DESC(I_SFLX_evap), VAR_UNIT(I_SFLX_evap), 'XY', OCEAN_RESTART_OUT_DTYPE, &
                                  VAR_ID(I_SFLX_evap), &
                                  standard_name=VAR_STDN(I_SFLX_evap) )

    endif

    return
  end subroutine OCEAN_vars_restart_def_var

  !-----------------------------------------------------------------------------
  !> Write ocean variables to restart file
  subroutine OCEAN_vars_restart_write
    use scale_file_cartesC, only: &
       FILE_CARTESC_write_var
    implicit none
    !---------------------------------------------------------------------------

    if ( restart_fid /= -1 ) then

       call OCEAN_vars_total

       call FILE_CARTESC_write_var( restart_fid, VAR_ID(I_TEMP), OCEAN_TEMP(:,:,:),              & ! [IN]
                              VAR_NAME(I_TEMP),      'OXY', fill_halo=.true.                 ) ! [IN]
!       call FILE_CARTESC_write_var( restart_fid, VAR_ID(I_SALT), OCEAN_SALT(:,:,:),              & ! [IN]
!                              VAR_NAME(I_SALT),      'OXY', fill_halo=.true.                 ) ! [IN]
!       call FILE_CARTESC_write_var( restart_fid, VAR_ID(I_UVEL), OCEAN_UVEL(:,:,:),              & ! [IN]
!                              VAR_NAME(I_UVEL),      'OXY', fill_halo=.true.                 ) ! [IN]
!       call FILE_CARTESC_write_var( restart_fid, VAR_ID(I_VVEL), OCEAN_VVEL(:,:,:),              & ! [IN]
!                              VAR_NAME(I_VVEL),      'OXY', fill_halo=.true.                 ) ! [IN]

       call FILE_CARTESC_write_var( restart_fid, VAR_ID(I_SFC_TEMP), OCEAN_SFC_TEMP(:,:),      & ! [IN]
                              VAR_NAME(I_SFC_TEMP),  'XY', fill_halo=.true.                 ) ! [IN]
       call FILE_CARTESC_write_var( restart_fid, VAR_ID(I_ALB_LW), OCEAN_SFC_albedo(:,:,I_LW), & ! [IN]
                              VAR_NAME(I_ALB_LW),    'XY', fill_halo=.true.                 ) ! [IN]
       call FILE_CARTESC_write_var( restart_fid, VAR_ID(I_ALB_SW), OCEAN_SFC_albedo(:,:,I_SW), & ! [IN]
                              VAR_NAME(I_ALB_SW),    'XY',  fill_halo=.true.                ) ! [IN]
       call FILE_CARTESC_write_var( restart_fid, VAR_ID(I_SFC_Z0M), OCEAN_SFC_Z0M(:,:),        & ! [IN]
                              VAR_NAME(I_SFC_Z0M),   'XY', fill_halo=.true.                 ) ! [IN]
       call FILE_CARTESC_write_var( restart_fid, VAR_ID(I_SFC_Z0H), OCEAN_SFC_Z0H(:,:),        & ! [IN]
                              VAR_NAME(I_SFC_Z0H),   'XY', fill_halo=.true.                 ) ! [IN]
       call FILE_CARTESC_write_var( restart_fid, VAR_ID(I_SFC_Z0E), OCEAN_SFC_Z0E(:,:),        & ! [IN]
                              VAR_NAME(I_SFC_Z0E),   'XY', fill_halo=.true.                 ) ! [IN]
       call FILE_CARTESC_write_var( restart_fid, VAR_ID(I_SFLX_MW), OCEAN_SFLX_MW(:,:),        & ! [IN]
                              VAR_NAME(I_SFLX_MW),   'XY', fill_halo=.true.                 ) ! [IN]
       call FILE_CARTESC_write_var( restart_fid, VAR_ID(I_SFLX_MU), OCEAN_SFLX_MU(:,:),        & ! [IN]
                              VAR_NAME(I_SFLX_MU),   'XY', fill_halo=.true.                 ) ! [IN]
       call FILE_CARTESC_write_var( restart_fid, VAR_ID(I_SFLX_MV), OCEAN_SFLX_MV(:,:),        & ! [IN]
                              VAR_NAME(I_SFLX_MV),   'XY', fill_halo=.true.                 ) ! [IN]
       call FILE_CARTESC_write_var( restart_fid, VAR_ID(I_SFLX_SH), OCEAN_SFLX_SH(:,:),        & ! [IN]
                              VAR_NAME(I_SFLX_SH),   'XY', fill_halo=.true.                 ) ! [IN]
       call FILE_CARTESC_write_var( restart_fid, VAR_ID(I_SFLX_LH), OCEAN_SFLX_LH(:,:),        & ! [IN]
                              VAR_NAME(I_SFLX_LH),   'XY', fill_halo=.true.                 ) ! [IN]
       call FILE_CARTESC_write_var( restart_fid, VAR_ID(I_SFLX_WH), OCEAN_SFLX_WH(:,:),        & ! [IN]
                              VAR_NAME(I_SFLX_WH),   'XY', fill_halo=.true.                 ) ! [IN]
       call FILE_CARTESC_write_var( restart_fid, VAR_ID(I_SFLX_evap), OCEAN_SFLX_evap(:,:),    & ! [IN]
                              VAR_NAME(I_SFLX_evap), 'XY', fill_halo=.true.                 ) ! [IN]

    endif

    return
  end subroutine OCEAN_vars_restart_write

end module mod_ocean_vars
