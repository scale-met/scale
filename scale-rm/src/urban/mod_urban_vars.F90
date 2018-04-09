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
#include "scalelib.h"
module mod_urban_vars
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_stdio
  use scale_prof
  use scale_debug
  use scale_urban_grid_cartesC_index

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
  public :: URBAN_vars_setup
  public :: URBAN_vars_restart_read
  public :: URBAN_vars_restart_write
  public :: URBAN_vars_history
  public :: URBAN_vars_total
  public :: URBAN_vars_external_in

  public :: URBAN_vars_restart_create
  public :: URBAN_vars_restart_open
  public :: URBAN_vars_restart_def_var
  public :: URBAN_vars_restart_enddef
  public :: URBAN_vars_restart_close

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  logical,               public :: URBAN_RESTART_OUTPUT                 = .false.         !< Output restart file?

  character(len=H_LONG),  public :: URBAN_RESTART_IN_BASENAME           = ''              !< Basename of the input  file
  logical,                public :: URBAN_RESTART_IN_AGGREGATE                            !< Switch to use aggregate file
  logical,                public :: URBAN_RESTART_IN_POSTFIX_TIMELABEL  = .false.         !< Add timelabel to the basename of input  file?
  character(len=H_LONG),  public :: URBAN_RESTART_OUT_BASENAME          = ''              !< Basename of the output file
  logical,                public :: URBAN_RESTART_OUT_AGGREGATE                           !< Switch to use aggregate file
  logical,                public :: URBAN_RESTART_OUT_POSTFIX_TIMELABEL = .true.          !< Add timelabel to the basename of output file?
  character(len=H_MID),   public :: URBAN_RESTART_OUT_TITLE             = 'URBAN restart' !< Title    of the output file
  character(len=H_SHORT), public :: URBAN_RESTART_OUT_DTYPE             = 'DEFAULT'       !< REAL4 or REAL8

  ! prognostic variables
  real(RP), public, allocatable :: URBAN_TR   (:,:)   ! urban surface temperature of roof [K]
  real(RP), public, allocatable :: URBAN_TB   (:,:)   ! urban surface temperature of wall [K]
  real(RP), public, allocatable :: URBAN_TG   (:,:)   ! urban surface temperature of road [K]
  real(RP), public, allocatable :: URBAN_TC   (:,:)   ! urban canopy air temperature [K]
  real(RP), public, allocatable :: URBAN_QC   (:,:)   ! urban canopy humidity [kg/kg]
  real(RP), public, allocatable :: URBAN_UC   (:,:)   ! urban canopy wind [m/s]
  real(RP), public, allocatable :: URBAN_TRL  (:,:,:) ! urban temperature in layer of roof [K]
  real(RP), public, allocatable :: URBAN_TBL  (:,:,:) ! urban temperature in layer of wall [K]
  real(RP), public, allocatable :: URBAN_TGL  (:,:,:) ! urban temperature in layer of road [K]
  real(RP), public, allocatable :: URBAN_RAINR(:,:)   ! urban rain storage on roof [mm=kg/m2]
  real(RP), public, allocatable :: URBAN_RAINB(:,:)   ! urban rain storage on wall [mm=kg/m2]
  real(RP), public, allocatable :: URBAN_RAING(:,:)   ! urban rain storage on road [mm=kg/m2]
  real(RP), public, allocatable :: URBAN_ROFF (:,:)   ! urban runoff [mm=kg/m2]

  ! tendency variables
  real(RP), public, allocatable :: URBAN_TRL_t  (:,:,:) ! tendency of URBAN_TRL
  real(RP), public, allocatable :: URBAN_TBL_t  (:,:,:) ! tendency of URBAN_TBL
  real(RP), public, allocatable :: URBAN_TGL_t  (:,:,:) ! tendency of URBAN_TGL
  real(RP), public, allocatable :: URBAN_TC_t   (:,:)   ! tendency of URBAN_TC
  real(RP), public, allocatable :: URBAN_UC_t   (:,:)   ! tendency of URBAN_UC
  real(RP), public, allocatable :: URBAN_QC_t   (:,:)   ! tendency of URBAN_QC
  real(RP), public, allocatable :: URBAN_TR_t   (:,:)   ! tendency of URBAN_TR
  real(RP), public, allocatable :: URBAN_TB_t   (:,:)   ! tendency of URBAN_TB
  real(RP), public, allocatable :: URBAN_TG_t   (:,:)   ! tendency of URBAN_TG
  real(RP), public, allocatable :: URBAN_RAINR_t(:,:)   ! tendency of URBAN_RAINR
  real(RP), public, allocatable :: URBAN_RAINB_t(:,:)   ! tendency of URBAN_RAINB
  real(RP), public, allocatable :: URBAN_RAING_t(:,:)   ! tendency of URBAN_RAING
  real(RP), public, allocatable :: URBAN_ROFF_t (:,:)   ! tendency of URBAN_ROFF

  ! for restart
  real(RP), public, allocatable :: URBAN_SFC_TEMP  (:,:)   ! urban grid average of surface temperature [K]
  real(RP), public, allocatable :: URBAN_SFC_albedo(:,:,:) ! urban grid average of albedo (0-1)
  real(RP), public, allocatable :: URBAN_SFLX_MW   (:,:)   ! urban grid average of w-momentum flux [kg/m2/s]
  real(RP), public, allocatable :: URBAN_SFLX_MU   (:,:)   ! urban grid average of u-momentum flux [kg/m2/s]
  real(RP), public, allocatable :: URBAN_SFLX_MV   (:,:)   ! urban grid average of v-momentum flux [kg/m2/s]
  real(RP), public, allocatable :: URBAN_SFLX_SH   (:,:)   ! urban grid average of sensible heat flux [W/m2]
  real(RP), public, allocatable :: URBAN_SFLX_LH   (:,:)   ! urban grid average of latent heat flux [W/m2]
  real(RP), public, allocatable :: URBAN_SFLX_GH   (:,:)   ! urban grid average of ground heat flux [W/m2]
  real(RP), public, allocatable :: URBAN_SFLX_evap (:,:)   ! urban grid average of water vapor flux [kg/m2/s]

  ! diagnostic variables
  real(RP), public, allocatable :: URBAN_Z0M(:,:) ! urban grid average of rougness length (momentum) [m]
  real(RP), public, allocatable :: URBAN_Z0H(:,:) ! urban grid average of rougness length (heat) [m]
  real(RP), public, allocatable :: URBAN_Z0E(:,:) ! urban grid average of rougness length (vapor) [m]
  real(RP), public, allocatable :: URBAN_U10(:,:) ! urban grid average of velocity u at 10m [m/s]
  real(RP), public, allocatable :: URBAN_V10(:,:) ! urban grid average of velocity v at 10m [m/s]
  real(RP), public, allocatable :: URBAN_T2 (:,:) ! urban grid average of temperature at 2m [K]
  real(RP), public, allocatable :: URBAN_Q2 (:,:) ! urban grid average of water vapor at 2m [kg/kg]

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
  real(RP), public, allocatable :: ATMOS_SFLX_LW  (:,:,:)
  real(RP), public, allocatable :: ATMOS_SFLX_SW  (:,:,:)
  real(RP), public, allocatable :: ATMOS_cosSZA   (:,:)
  real(RP), public, allocatable :: ATMOS_SFLX_rain(:,:)
  real(RP), public, allocatable :: ATMOS_SFLX_snow(:,:)

  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  logical,                private :: URBAN_VARS_CHECKRANGE      = .false.

  integer,                private, parameter :: VMAX        = 23
  integer,                private, parameter :: I_TR        = 1
  integer,                private, parameter :: I_TB        = 2
  integer,                private, parameter :: I_TG        = 3
  integer,                private, parameter :: I_TC        = 4
  integer,                private, parameter :: I_QC        = 5
  integer,                private, parameter :: I_UC        = 6
  integer,                private, parameter :: I_TRL       = 7
  integer,                private, parameter :: I_TBL       = 8
  integer,                private, parameter :: I_TGL       = 9
  integer,                private, parameter :: I_RAINR     = 10
  integer,                private, parameter :: I_RAINB     = 11
  integer,                private, parameter :: I_RAING     = 12
  integer,                private, parameter :: I_ROFF      = 13
  integer,                private, parameter :: I_SFC_TEMP  = 14
  integer,                private, parameter :: I_ALB_LW    = 15
  integer,                private, parameter :: I_ALB_SW    = 16
  integer,                private, parameter :: I_SFLX_MW   = 17
  integer,                private, parameter :: I_SFLX_MU   = 18
  integer,                private, parameter :: I_SFLX_MV   = 19
  integer,                private, parameter :: I_SFLX_SH   = 20
  integer,                private, parameter :: I_SFLX_LH   = 21
  integer,                private, parameter :: I_SFLX_GH   = 22
  integer,                private, parameter :: I_SFLX_evap = 23

  character(len=H_SHORT), private            :: VAR_NAME(VMAX) !< name  of the urban variables
  character(len=H_MID),   private            :: VAR_DESC(VMAX) !< desc. of the urban variables
  character(len=H_SHORT), private            :: VAR_UNIT(VMAX) !< unit  of the urban variables
  integer,                private            :: VAR_ID(VMAX)   !< ID    of the urban variables
  integer,                private            :: restart_fid = -1  ! file ID

  logical,                private            :: URBAN_RESTART_IN_CHECK_COORDINATES = .true.

  data VAR_NAME / 'URBAN_TR' ,       &
                  'URBAN_TB' ,       &
                  'URBAN_TG' ,       &
                  'URBAN_TC' ,       &
                  'URBAN_QC' ,       &
                  'URBAN_UC' ,       &
                  'URBAN_TRL' ,      &
                  'URBAN_TBL' ,      &
                  'URBAN_TGL' ,      &
                  'URBAN_RAINR' ,    &
                  'URBAN_RAINB' ,    &
                  'URBAN_RAING' ,    &
                  'URBAN_ROFF',      &
                  'URBAN_SFC_TEMP',  &
                  'URBAN_ALB_LW',    &
                  'URBAN_ALB_SW',    &
                  'URBAN_SFLX_MW',   &
                  'URBAN_SFLX_MU',   &
                  'URBAN_SFLX_MV',   &
                  'URBAN_SFLX_SH',   &
                  'URBAN_SFLX_LH',   &
                  'URBAN_SFLX_GH',   &
                  'URBAN_SFLX_evap'  /

  data VAR_DESC / 'urban surface temperature of roof',                &
                  'urban surface temperature of wall',                &
                  'urban surface temperature of road',                &
                  'urban canopy air temperature',                     &
                  'urban canopy humidity',                            &
                  'urban canopy wind',                                &
                  'urban temperature in layer of roof',               &
                  'urban temperature in layer of wall',               &
                  'urban temperature in layer of road',               &
                  'urban rain strage on roof',                        &
                  'urban rain strage on wall',                        &
                  'urban rain strage on road',                        &
                  'urban runoff ',                                    &
                  'urban grid average of temperature',                &
                  'urban grid average of albedo LW',                  &
                  'urban grid average of albedo SW',                  &
                  'urban grid average of w-momentum flux',            &
                  'urban grid average of u-momentum flux',            &
                  'urban grid average of v-momentum flux',            &
                  'urban grid average of sensible heat flux',         &
                  'urban grid average of latent heat flux',           &
                  'urban grid average of ground heat flux',           &
                  'urban grid average of water vapor flux'            /

  data VAR_UNIT / 'K',       &
                  'K',       &
                  'K',       &
                  'K',       &
                  'kg/kg',   &
                  'm/s',     &
                  'K',       &
                  'K',       &
                  'K',       &
                  'kg/m2',   &
                  'kg/m2',   &
                  'kg/m2',   &
                  'kg/m2',   &
                  'K',       &
                  '1',       &
                  '1',       &
                  'kg/m2/s', &
                  'kg/m2/s', &
                  'kg/m2/s', &
                  'W/m2',    &
                  'W/m2',    &
                  'W/m2',    &
                  'kg/m2/s'  /

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine URBAN_vars_setup
    use scale_prc, only: &
       PRC_abort
    use scale_const, only: &
       UNDEF => CONST_UNDEF
    implicit none

    NAMELIST / PARAM_URBAN_VARS /  &
       URBAN_RESTART_IN_BASENAME,           &
       URBAN_RESTART_IN_AGGREGATE,          &
       URBAN_RESTART_IN_POSTFIX_TIMELABEL,  &
       URBAN_RESTART_IN_CHECK_COORDINATES,  &
       URBAN_RESTART_OUTPUT,                &
       URBAN_RESTART_OUT_BASENAME,          &
       URBAN_RESTART_OUT_AGGREGATE,         &
       URBAN_RESTART_OUT_POSTFIX_TIMELABEL, &
       URBAN_RESTART_OUT_TITLE,             &
       URBAN_RESTART_OUT_DTYPE,             &
       URBAN_VARS_CHECKRANGE

    integer :: ierr
    integer :: iv
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_PROGRESS(*) 'Module[VARS] / Categ[URBAN] / Origin[SCALE-RM]'

    allocate( URBAN_TR   (UIA,UJA)         )
    allocate( URBAN_TB   (UIA,UJA)         )
    allocate( URBAN_TG   (UIA,UJA)         )
    allocate( URBAN_TC   (UIA,UJA)         )
    allocate( URBAN_QC   (UIA,UJA)         )
    allocate( URBAN_UC   (UIA,UJA)         )
    allocate( URBAN_TRL  (UKS:UKE,UIA,UJA) )
    allocate( URBAN_TBL  (UKS:UKE,UIA,UJA) )
    allocate( URBAN_TGL  (UKS:UKE,UIA,UJA) )
    allocate( URBAN_RAINR(UIA,UJA)         )
    allocate( URBAN_RAINB(UIA,UJA)         )
    allocate( URBAN_RAING(UIA,UJA)         )
    allocate( URBAN_ROFF (UIA,UJA)         )
    URBAN_TR   (:,:)   = UNDEF
    URBAN_TB   (:,:)   = UNDEF
    URBAN_TG   (:,:)   = UNDEF
    URBAN_TC   (:,:)   = UNDEF
    URBAN_QC   (:,:)   = UNDEF
    URBAN_UC   (:,:)   = UNDEF
    URBAN_TRL  (:,:,:) = UNDEF
    URBAN_TBL  (:,:,:) = UNDEF
    URBAN_TGL  (:,:,:) = UNDEF
    URBAN_RAINR(:,:)   = UNDEF
    URBAN_RAINB(:,:)   = UNDEF
    URBAN_RAING(:,:)   = UNDEF
    URBAN_ROFF (:,:)   = UNDEF

    allocate( URBAN_TR_t   (UIA,UJA)         )
    allocate( URBAN_TB_t   (UIA,UJA)         )
    allocate( URBAN_TG_t   (UIA,UJA)         )
    allocate( URBAN_TC_t   (UIA,UJA)         )
    allocate( URBAN_QC_t   (UIA,UJA)         )
    allocate( URBAN_UC_t   (UIA,UJA)         )
    allocate( URBAN_TRL_t  (UKS:UKE,UIA,UJA) )
    allocate( URBAN_TBL_t  (UKS:UKE,UIA,UJA) )
    allocate( URBAN_TGL_t  (UKS:UKE,UIA,UJA) )
    allocate( URBAN_RAINR_t(UIA,UJA)         )
    allocate( URBAN_RAINB_t(UIA,UJA)         )
    allocate( URBAN_RAING_t(UIA,UJA)         )
    allocate( URBAN_ROFF_t (UIA,UJA)         )
    URBAN_TR_t   (:,:)   = UNDEF
    URBAN_TB_t   (:,:)   = UNDEF
    URBAN_TG_t   (:,:)   = UNDEF
    URBAN_TC_t   (:,:)   = UNDEF
    URBAN_QC_t   (:,:)   = UNDEF
    URBAN_UC_t   (:,:)   = UNDEF
    URBAN_TRL_t  (:,:,:) = UNDEF
    URBAN_TBL_t  (:,:,:) = UNDEF
    URBAN_TGL_t  (:,:,:) = UNDEF
    URBAN_RAINR_t(:,:)   = UNDEF
    URBAN_RAINB_t(:,:)   = UNDEF
    URBAN_RAING_t(:,:)   = UNDEF
    URBAN_ROFF_t (:,:)   = UNDEF

    allocate( URBAN_SFC_TEMP  (UIA,UJA)   )
    allocate( URBAN_SFC_albedo(UIA,UJA,2) )
    allocate( URBAN_SFLX_MW   (UIA,UJA)   )
    allocate( URBAN_SFLX_MU   (UIA,UJA)   )
    allocate( URBAN_SFLX_MV   (UIA,UJA)   )
    allocate( URBAN_SFLX_SH   (UIA,UJA)   )
    allocate( URBAN_SFLX_LH   (UIA,UJA)   )
    allocate( URBAN_SFLX_GH   (UIA,UJA)   )
    allocate( URBAN_SFLX_evap (UIA,UJA)   )
    URBAN_SFC_TEMP  (:,:)   = UNDEF
    URBAN_SFC_albedo(:,:,:) = UNDEF
    URBAN_SFLX_MW   (:,:)   = UNDEF
    URBAN_SFLX_MU   (:,:)   = UNDEF
    URBAN_SFLX_MV   (:,:)   = UNDEF
    URBAN_SFLX_SH   (:,:)   = UNDEF
    URBAN_SFLX_LH   (:,:)   = UNDEF
    URBAN_SFLX_GH   (:,:)   = UNDEF
    URBAN_SFLX_evap (:,:)   = UNDEF

    allocate( URBAN_Z0M(UIA,UJA) )
    allocate( URBAN_Z0H(UIA,UJA) )
    allocate( URBAN_Z0E(UIA,UJA) )
    allocate( URBAN_U10(UIA,UJA) )
    allocate( URBAN_V10(UIA,UJA) )
    allocate( URBAN_T2 (UIA,UJA) )
    allocate( URBAN_Q2 (UIA,UJA) )
    URBAN_Z0M(:,:) = UNDEF
    URBAN_Z0H(:,:) = UNDEF
    URBAN_Z0E(:,:) = UNDEF
    URBAN_U10(:,:) = UNDEF
    URBAN_V10(:,:) = UNDEF
    URBAN_T2 (:,:) = UNDEF
    URBAN_Q2 (:,:) = UNDEF

    allocate( ATMOS_TEMP     (UIA,UJA)   )
    allocate( ATMOS_PRES     (UIA,UJA)   )
    allocate( ATMOS_W        (UIA,UJA)   )
    allocate( ATMOS_U        (UIA,UJA)   )
    allocate( ATMOS_V        (UIA,UJA)   )
    allocate( ATMOS_DENS     (UIA,UJA)   )
    allocate( ATMOS_QV       (UIA,UJA)   )
    allocate( ATMOS_PBL      (UIA,UJA)   )
    allocate( ATMOS_SFC_DENS (UIA,UJA)   )
    allocate( ATMOS_SFC_PRES (UIA,UJA)   )
    allocate( ATMOS_SFLX_LW  (UIA,UJA,2) )
    allocate( ATMOS_SFLX_SW  (UIA,UJA,2) )
    allocate( ATMOS_cosSZA   (UIA,UJA)   )
    allocate( ATMOS_SFLX_rain(UIA,UJA)   )
    allocate( ATMOS_SFLX_snow(UIA,UJA)   )
    ATMOS_TEMP     (:,:)   = UNDEF
    ATMOS_PRES     (:,:)   = UNDEF
    ATMOS_W        (:,:)   = UNDEF
    ATMOS_U        (:,:)   = UNDEF
    ATMOS_V        (:,:)   = UNDEF
    ATMOS_DENS     (:,:)   = UNDEF
    ATMOS_QV       (:,:)   = UNDEF
    ATMOS_PBL      (:,:)   = UNDEF
    ATMOS_SFC_DENS (:,:)   = UNDEF
    ATMOS_SFC_PRES (:,:)   = UNDEF
    ATMOS_SFLX_LW  (:,:,:) = UNDEF
    ATMOS_SFLX_SW  (:,:,:) = UNDEF
    ATMOS_cosSZA   (:,:)   = UNDEF
    ATMOS_SFLX_rain(:,:)   = UNDEF
    ATMOS_SFLX_snow(:,:)   = UNDEF

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_URBAN_VARS,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("URBAN_vars_setup",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("URBAN_vars_setup",*) 'Not appropriate names in namelist PARAM_URBAN_VARS. Check!'
       call PRC_abort
    endif
    if( IO_NML ) write(IO_FID_NML,nml=PARAM_URBAN_VARS)

    LOG_NEWLINE
    LOG_INFO("URBAN_vars_setup",*) 'List of prognostic variables (URBAN) '
    LOG_INFO_CONT('(1x,A,A24,A,A48,A,A12,A)') &
               '      |', 'VARNAME                 ','|', &
               'DESCRIPTION                                     ', '[', 'UNIT        ', ']'
    do iv = 1, VMAX
       LOG_INFO_CONT('(1x,A,I3,A,A24,A,A48,A,A12,A)') &
                  'NO.',iv,'|',VAR_NAME(iv),'|',VAR_DESC(iv),'[',VAR_UNIT(iv),']'
    enddo

    LOG_NEWLINE
    if ( URBAN_RESTART_IN_BASENAME /= '' ) then
       LOG_INFO("URBAN_vars_setup",*) 'Restart input?  : YES, file = ', trim(URBAN_RESTART_IN_BASENAME)
       LOG_INFO("URBAN_vars_setup",*) 'Add timelabel?  : ', URBAN_RESTART_IN_POSTFIX_TIMELABEL
    else
       LOG_INFO("URBAN_vars_setup",*) 'Restart input?  : NO'
    endif
    if (       URBAN_RESTART_OUTPUT             &
         .AND. URBAN_RESTART_OUT_BASENAME /= '' ) then
       LOG_INFO("URBAN_vars_setup",*) 'Restart output? : YES, file = ', trim(URBAN_RESTART_OUT_BASENAME)
       LOG_INFO("URBAN_vars_setup",*) 'Add timelabel?  : ', URBAN_RESTART_OUT_POSTFIX_TIMELABEL
    else
       LOG_INFO("URBAN_vars_setup",*) 'Restart output? : NO'
       URBAN_RESTART_OUTPUT = .false.
    endif

    return
  end subroutine URBAN_vars_setup

  !-----------------------------------------------------------------------------
  !> Open urban restart file for read
  subroutine URBAN_vars_restart_open
    use scale_time, only: &
       TIME_gettimelabel
    use scale_file_cartesC, only: &
       FILE_CARTESC_open, &
       FILE_CARTESC_check_coordinates
    use mod_urban_admin, only: &
       URBAN_do
    implicit none

    character(len=19)     :: timelabel
    character(len=H_LONG) :: basename
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("URBAN_vars_restart_open",*) 'Open restart file (URBAN) '

    if ( URBAN_do .and. URBAN_RESTART_IN_BASENAME /= '' ) then

       if ( URBAN_RESTART_IN_POSTFIX_TIMELABEL ) then
          call TIME_gettimelabel( timelabel )
          basename = trim(URBAN_RESTART_IN_BASENAME)//'_'//trim(timelabel)
       else
          basename = trim(URBAN_RESTART_IN_BASENAME)
       endif

       LOG_INFO("URBAN_vars_restart_open",*) 'basename: ', trim(basename)

       call FILE_CARTESC_open( basename, restart_fid, aggregate=URBAN_RESTART_IN_AGGREGATE )

       if ( URBAN_RESTART_IN_CHECK_COORDINATES ) then
          call FILE_CARTESC_check_coordinates( restart_fid, urban=.true. )
       end if

    else
       LOG_INFO("URBAN_vars_restart_open",*) 'restart file for urban is not specified.'
    endif

    return
  end subroutine URBAN_vars_restart_open

  !-----------------------------------------------------------------------------
  !> Read urban restart
  subroutine URBAN_vars_restart_read
    use scale_prc, only: &
       PRC_abort
    use scale_file, only: &
       FILE_get_aggregate
    use scale_file_cartesC, only: &
       FILE_CARTESC_read, &
       FILE_CARTESC_flush
    implicit none
    !---------------------------------------------------------------------------

    if ( restart_fid /= -1 ) then
       LOG_NEWLINE
       LOG_INFO("URBAN_vars_restart_read",*) 'Read from restart file (URBAN) '

       call FILE_CARTESC_read( restart_fid, VAR_NAME(I_TR), 'XY', & ! [IN]
                               URBAN_TR(:,:)                      ) ! [OUT]
       call FILE_CARTESC_read( restart_fid, VAR_NAME(I_TB), 'XY', & ! [IN]
                               URBAN_TB(:,:)                      ) ! [OUT]
       call FILE_CARTESC_read( restart_fid, VAR_NAME(I_TG), 'XY', & ! [IN]
                               URBAN_TG(:,:)                      ) ! [OUT]
       call FILE_CARTESC_read( restart_fid, VAR_NAME(I_TC), 'XY', & ! [IN]
                               URBAN_TC(:,:)                      ) ! [OUT]
       call FILE_CARTESC_read( restart_fid, VAR_NAME(I_QC), 'XY', & ! [IN]
                               URBAN_QC(:,:)                      ) ! [OUT]
       call FILE_CARTESC_read( restart_fid, VAR_NAME(I_UC), 'XY', & ! [IN]
                               URBAN_UC(:,:)                      ) ! [OUT]

       call FILE_CARTESC_read( restart_fid, VAR_NAME(I_TRL), 'UXY', & ! [IN]
                               URBAN_TRL(:,:,:)                     ) ! [OUT]
       call FILE_CARTESC_read( restart_fid, VAR_NAME(I_TBL), 'UXY', & ! [IN]
                               URBAN_TBL(:,:,:)                     ) ! [OUT]
       call FILE_CARTESC_read( restart_fid, VAR_NAME(I_TGL), 'UXY', & ! [IN]
                               URBAN_TGL(:,:,:)                     ) ! [OUT]

       call FILE_CARTESC_read( restart_fid, VAR_NAME(I_RAINR), 'XY', & ! [IN]
                               URBAN_RAINR(:,:)                      ) ! [OUT]
       call FILE_CARTESC_read( restart_fid, VAR_NAME(I_RAINB), 'XY', & ! [IN]
                               URBAN_RAINB(:,:)                      ) ! [OUT]
       call FILE_CARTESC_read( restart_fid, VAR_NAME(I_RAING), 'XY', & ! [IN]
                               URBAN_RAING(:,:)                      ) ! [OUT]
       call FILE_CARTESC_read( restart_fid, VAR_NAME(I_ROFF),  'XY', & ! [IN]
                               URBAN_ROFF(:,:)                       ) ! [OUT]

       call FILE_CARTESC_read( restart_fid, VAR_NAME(I_SFC_TEMP), 'XY', & ! [IN]
                               URBAN_SFC_TEMP(:,:)                      ) ! [OUT]
       call FILE_CARTESC_read( restart_fid, VAR_NAME(I_ALB_LW),   'XY', & ! [IN]
                               URBAN_SFC_albedo(:,:,I_LW)               ) ! [OUT]
       call FILE_CARTESC_read( restart_fid, VAR_NAME(I_ALB_SW),   'XY', & ! [IN]
                               URBAN_SFC_albedo(:,:,I_SW)               ) ! [OUT]

       if( FILE_get_AGGREGATE(restart_fid) ) call FILE_CARTESC_flush( restart_fid ) ! commit all pending read requests

       call URBAN_vars_total
    else
       LOG_ERROR("URBAN_vars_restart_read",*) 'invalid restart file ID for urban.'
       call PRC_abort
    endif

    return
  end subroutine URBAN_vars_restart_read

  !-----------------------------------------------------------------------------
  !> History output set for urban variables
  subroutine URBAN_vars_history
    use scale_file_history, only: &
       FILE_HISTORY_in
    implicit none
    !---------------------------------------------------------------------------

    call PROF_rapstart('URB_History', 1)

    if ( URBAN_VARS_CHECKRANGE ) then
       call VALCHECK( URBAN_TR        (UIS:UIE,UJS:UJE),          0.0_RP, 1000.0_RP, VAR_NAME(I_TR),        &
                     __FILE__, __LINE__ )
       call VALCHECK( URBAN_TB        (UIS:UIE,UJS:UJE),          0.0_RP, 1000.0_RP, VAR_NAME(I_TB),        &
                     __FILE__, __LINE__ )
       call VALCHECK( URBAN_TG        (UIS:UIE,UJS:UJE),          0.0_RP, 1000.0_RP, VAR_NAME(I_TG),        &
                     __FILE__, __LINE__ )
       call VALCHECK( URBAN_TC        (UIS:UIE,UJS:UJE),          0.0_RP, 1000.0_RP, VAR_NAME(I_TC),        &
                     __FILE__, __LINE__ )
       call VALCHECK( URBAN_QC        (UIS:UIE,UJS:UJE),          0.0_RP, 1000.0_RP, VAR_NAME(I_QC),        &
                     __FILE__, __LINE__ )
       call VALCHECK( URBAN_UC        (UIS:UIE,UJS:UJE),          0.0_RP, 1000.0_RP, VAR_NAME(I_UC),        &
                     __FILE__, __LINE__ )
       call VALCHECK( URBAN_TRL       (:,UIS:UIE,UJS:UJE),        0.0_RP, 1000.0_RP, VAR_NAME(I_TRL),       &
                     __FILE__, __LINE__ )
       call VALCHECK( URBAN_TBL       (:,UIS:UIE,UJS:UJE),        0.0_RP, 1000.0_RP, VAR_NAME(I_TBL),       &
                     __FILE__, __LINE__ )
       call VALCHECK( URBAN_TGL       (:,UIS:UIE,UJS:UJE),        0.0_RP, 1000.0_RP, VAR_NAME(I_TGL),       &
                     __FILE__, __LINE__ )
       call VALCHECK( URBAN_RAINR     (UIS:UIE,UJS:UJE),      -1000.0_RP, 1000.0_RP, VAR_NAME(I_RAINR),     &
                     __FILE__, __LINE__ )
       call VALCHECK( URBAN_RAINB     (UIS:UIE,UJS:UJE),      -1000.0_RP, 1000.0_RP, VAR_NAME(I_RAINB),     &
                     __FILE__, __LINE__ )
       call VALCHECK( URBAN_RAING     (UIS:UIE,UJS:UJE),      -1000.0_RP, 1000.0_RP, VAR_NAME(I_RAING),     &
                     __FILE__, __LINE__ )
       call VALCHECK( URBAN_ROFF      (UIS:UIE,UJS:UJE),      -1000.0_RP, 1000.0_RP, VAR_NAME(I_ROFF),      &
                     __FILE__, __LINE__ )
       call VALCHECK( URBAN_SFC_TEMP  (UIS:UIE,UJS:UJE),          0.0_RP, 1000.0_RP, VAR_NAME(I_SFC_TEMP),  &
                     __FILE__, __LINE__ )
       call VALCHECK( URBAN_SFC_albedo(UIS:UIE,UJS:UJE,I_LW),     0.0_RP,    2.0_RP, VAR_NAME(I_ALB_LW),    &
                     __FILE__, __LINE__ )
       call VALCHECK( URBAN_SFC_albedo(UIS:UIE,UJS:UJE,I_SW),     0.0_RP,    2.0_RP, VAR_NAME(I_ALB_SW),    &
                     __FILE__, __LINE__ )
    endif

    call FILE_HISTORY_in( URBAN_TR(:,:), VAR_NAME(I_TR), VAR_DESC(I_TR), VAR_UNIT(I_TR) )
    call FILE_HISTORY_in( URBAN_TB(:,:), VAR_NAME(I_TB), VAR_DESC(I_TB), VAR_UNIT(I_TB) )
    call FILE_HISTORY_in( URBAN_TG(:,:), VAR_NAME(I_TG), VAR_DESC(I_TG), VAR_UNIT(I_TG) )
    call FILE_HISTORY_in( URBAN_TC(:,:), VAR_NAME(I_TC), VAR_DESC(I_TC), VAR_UNIT(I_TC) )
    call FILE_HISTORY_in( URBAN_QC(:,:), VAR_NAME(I_QC), VAR_DESC(I_QC), VAR_UNIT(I_QC) )
    call FILE_HISTORY_in( URBAN_UC(:,:), VAR_NAME(I_UC), VAR_DESC(I_UC), VAR_UNIT(I_UC) )

    call FILE_HISTORY_in( URBAN_TRL(:,:,:), VAR_NAME(I_TRL), VAR_DESC(I_TRL), VAR_UNIT(I_TRL), dim_type='UXY' )
    call FILE_HISTORY_in( URBAN_TBL(:,:,:), VAR_NAME(I_TBL), VAR_DESC(I_TBL), VAR_UNIT(I_TBL), dim_type='UXY' )
    call FILE_HISTORY_in( URBAN_TGL(:,:,:), VAR_NAME(I_TGL), VAR_DESC(I_TGL), VAR_UNIT(I_TGL), dim_type='UXY' )

    call FILE_HISTORY_in( URBAN_RAINR(:,:), VAR_NAME(I_RAINR), VAR_DESC(I_RAINR), VAR_UNIT(I_RAINR) )
    call FILE_HISTORY_in( URBAN_RAINB(:,:), VAR_NAME(I_RAINB), VAR_DESC(I_RAINB), VAR_UNIT(I_RAINB) )
    call FILE_HISTORY_in( URBAN_RAING(:,:), VAR_NAME(I_RAING), VAR_DESC(I_RAING), VAR_UNIT(I_RAING) )
    call FILE_HISTORY_in( URBAN_ROFF (:,:), VAR_NAME(I_ROFF),  VAR_DESC(I_ROFF),  VAR_UNIT(I_ROFF)  )

    call FILE_HISTORY_in( URBAN_SFC_TEMP  (:,:),      VAR_NAME(I_SFC_TEMP), VAR_DESC(I_SFC_TEMP), VAR_UNIT(I_SFC_TEMP) )
    call FILE_HISTORY_in( URBAN_SFC_albedo(:,:,I_LW), VAR_NAME(I_ALB_LW),   VAR_DESC(I_ALB_LW),   VAR_UNIT(I_ALB_LW)   )
    call FILE_HISTORY_in( URBAN_SFC_albedo(:,:,I_SW), VAR_NAME(I_ALB_SW),   VAR_DESC(I_ALB_SW),   VAR_UNIT(I_ALB_SW)   )

    call FILE_HISTORY_in( URBAN_SFLX_MW  (:,:), VAR_NAME(I_SFLX_MW),   VAR_DESC(I_SFLX_MW),   VAR_UNIT(I_SFLX_MW)   )
    call FILE_HISTORY_in( URBAN_SFLX_MU  (:,:), VAR_NAME(I_SFLX_MU),   VAR_DESC(I_SFLX_MU),   VAR_UNIT(I_SFLX_MU)   )
    call FILE_HISTORY_in( URBAN_SFLX_MV  (:,:), VAR_NAME(I_SFLX_MV),   VAR_DESC(I_SFLX_MV),   VAR_UNIT(I_SFLX_MV)   )
    call FILE_HISTORY_in( URBAN_SFLX_SH  (:,:), VAR_NAME(I_SFLX_SH),   VAR_DESC(I_SFLX_SH),   VAR_UNIT(I_SFLX_SH)   )
    call FILE_HISTORY_in( URBAN_SFLX_LH  (:,:), VAR_NAME(I_SFLX_LH),   VAR_DESC(I_SFLX_LH),   VAR_UNIT(I_SFLX_LH)   )
    call FILE_HISTORY_in( URBAN_SFLX_GH  (:,:), VAR_NAME(I_SFLX_GH),   VAR_DESC(I_SFLX_GH),   VAR_UNIT(I_SFLX_GH)   )
    call FILE_HISTORY_in( URBAN_SFLX_evap(:,:), VAR_NAME(I_SFLX_evap), VAR_DESC(I_SFLX_evap), VAR_UNIT(I_SFLX_evap) )

    call PROF_rapend  ('URB_History', 1)

    return
  end subroutine URBAN_vars_history

  !-----------------------------------------------------------------------------
  !> Budget monitor for urban
  subroutine URBAN_vars_total
    use scale_statistics, only: &
       STATISTICS_checktotal, &
       STATISTICS_total
    use scale_urban_grid_cartesC_real, only: &
       URBAN_GRID_CARTESC_REAL_AREA,    &
       URBAN_GRID_CARTESC_REAL_TOTAREA, &
       URBAN_GRID_CARTESC_REAL_VOL,     &
       URBAN_GRID_CARTESC_REAL_TOTVOL
    implicit none

    !---------------------------------------------------------------------------

    if ( STATISTICS_checktotal ) then

       ! 3D
       call STATISTICS_total( UKA, UKS, UKE, UIA, UIS, UIE, UJA, UJS, UJE, &
                              URBAN_TRL(:,:,:), VAR_NAME(I_TRL),  & ! (in)
                              URBAN_GRID_CARTESC_REAL_VOL(:,:,:), & ! (in)
                              URBAN_GRID_CARTESC_REAL_TOTVOL      ) ! (in)
       call STATISTICS_total( UKA, UKS, UKE, UIA, UIS, UIE, UJA, UJS, UJE, &
                              URBAN_TBL(:,:,:), VAR_NAME(I_TBL),  & ! (in)
                              URBAN_GRID_CARTESC_REAL_VOL(:,:,:), & ! (in)
                              URBAN_GRID_CARTESC_REAL_TOTVOL      ) ! (in)
       call STATISTICS_total( UKA, UKS, UKE, UIA, UIS, UIE, UJA, UJS, UJE, &
                              URBAN_TGL(:,:,:), VAR_NAME(I_TGL),  & ! (in)
                              URBAN_GRID_CARTESC_REAL_VOL(:,:,:), & ! (in)
                              URBAN_GRID_CARTESC_REAL_TOTVOL      ) ! (in)

       ! 2D
       call STATISTICS_total( UIA, UIS, UIE, UJA, UJS, UJE, &
                              URBAN_TR(:,:), VAR_NAME(I_TR),      & ! (in)
                              URBAN_GRID_CARTESC_REAL_AREA(:,:),  & ! (in)
                              URBAN_GRID_CARTESC_REAL_TOTAREA     ) ! (in)
       call STATISTICS_total( UIA, UIS, UIE, UJA, UJS, UJE, &
                              URBAN_TB(:,:), VAR_NAME(I_TB),      & ! (in)
                              URBAN_GRID_CARTESC_REAL_AREA(:,:),  & ! (in)
                              URBAN_GRID_CARTESC_REAL_TOTAREA     ) ! (in)
       call STATISTICS_total( UIA, UIS, UIE, UJA, UJS, UJE, &
                              URBAN_TG(:,:), VAR_NAME(I_TG),      & ! (in)
                              URBAN_GRID_CARTESC_REAL_AREA(:,:),  & ! (in)
                              URBAN_GRID_CARTESC_REAL_TOTAREA     ) ! (in)
       call STATISTICS_total( UIA, UIS, UIE, UJA, UJS, UJE, &
                              URBAN_TC(:,:), VAR_NAME(I_TC),      & ! (in)
                              URBAN_GRID_CARTESC_REAL_AREA(:,:),  & ! (in)
                              URBAN_GRID_CARTESC_REAL_TOTAREA     ) ! (in)
       call STATISTICS_total( UIA, UIS, UIE, UJA, UJS, UJE, &
                              URBAN_QC(:,:), VAR_NAME(I_QC),      & ! (in)
                              URBAN_GRID_CARTESC_REAL_AREA(:,:),  & ! (in)
                              URBAN_GRID_CARTESC_REAL_TOTAREA     ) ! (in)
       call STATISTICS_total( UIA, UIS, UIE, UJA, UJS, UJE, &
                              URBAN_UC(:,:), VAR_NAME(I_UC),      & ! (in)
                              URBAN_GRID_CARTESC_REAL_AREA(:,:),  & ! (in)
                              URBAN_GRID_CARTESC_REAL_TOTAREA     ) ! (in)

       call STATISTICS_total( UIA, UIS, UIE, UJA, UJS, UJE, &
                              URBAN_RAINR(:,:), VAR_NAME(I_RAINR), & ! (in)
                              URBAN_GRID_CARTESC_REAL_AREA(:,:),   & ! (in)
                              URBAN_GRID_CARTESC_REAL_TOTAREA      ) ! (in)
       call STATISTICS_total( UIA, UIS, UIE, UJA, UJS, UJE, &
                              URBAN_RAINB(:,:), VAR_NAME(I_RAINB), & ! (in)
                              URBAN_GRID_CARTESC_REAL_AREA(:,:),   & ! (in)
                              URBAN_GRID_CARTESC_REAL_TOTAREA      ) ! (in)
       call STATISTICS_total( UIA, UIS, UIE, UJA, UJS, UJE, &
                              URBAN_RAING(:,:), VAR_NAME(I_RAING), & ! (in)
                              URBAN_GRID_CARTESC_REAL_AREA(:,:),   & ! (in)
                              URBAN_GRID_CARTESC_REAL_TOTAREA      ) ! (in)
       call STATISTICS_total( UIA, UIS, UIE, UJA, UJS, UJE, &
                              URBAN_ROFF (:,:), VAR_NAME(I_ROFF),  & ! (in)
                              URBAN_GRID_CARTESC_REAL_AREA(:,:),   & ! (in)
                              URBAN_GRID_CARTESC_REAL_TOTAREA      ) ! (in)

       call STATISTICS_total( UIA, UIS, UIE, UJA, UJS, UJE, &
                              URBAN_SFC_TEMP  (:,:),      VAR_NAME(I_SFC_TEMP), & ! (in)
                              URBAN_GRID_CARTESC_REAL_AREA(:,:),                & ! (in)
                              URBAN_GRID_CARTESC_REAL_TOTAREA                   ) ! (in)
       call STATISTICS_total( UIA, UIS, UIE, UJA, UJS, UJE, &
                              URBAN_SFC_albedo(:,:,I_LW), VAR_NAME(I_ALB_LW),   & ! (in)
                              URBAN_GRID_CARTESC_REAL_AREA(:,:),                & ! (in)
                              URBAN_GRID_CARTESC_REAL_TOTAREA                   ) ! (in)
       call STATISTICS_total( UIA, UIS, UIE, UJA, UJS, UJE, &
                              URBAN_SFC_albedo(:,:,I_SW), VAR_NAME(I_ALB_SW),   & ! (in)
                              URBAN_GRID_CARTESC_REAL_AREA(:,:),                & ! (in)
                              URBAN_GRID_CARTESC_REAL_TOTAREA                   ) ! (in)

       call STATISTICS_total( UIA, UIS, UIE, UJA, UJS, UJE, &
                              URBAN_SFLX_MW  (:,:), VAR_NAME(I_SFLX_MW),   & ! (in)
                              URBAN_GRID_CARTESC_REAL_AREA(:,:),           & ! (in)
                              URBAN_GRID_CARTESC_REAL_TOTAREA              ) ! (in)
       call STATISTICS_total( UIA, UIS, UIE, UJA, UJS, UJE, &
                              URBAN_SFLX_MU  (:,:), VAR_NAME(I_SFLX_MU),   & ! (in)
                              URBAN_GRID_CARTESC_REAL_AREA(:,:),           & ! (in)
                              URBAN_GRID_CARTESC_REAL_TOTAREA              ) ! (in)
       call STATISTICS_total( UIA, UIS, UIE, UJA, UJS, UJE, &
                              URBAN_SFLX_MV  (:,:), VAR_NAME(I_SFLX_MV),   & ! (in)
                              URBAN_GRID_CARTESC_REAL_AREA(:,:),           & ! (in)
                              URBAN_GRID_CARTESC_REAL_TOTAREA              ) ! (in)
       call STATISTICS_total( UIA, UIS, UIE, UJA, UJS, UJE, &
                              URBAN_SFLX_SH  (:,:), VAR_NAME(I_SFLX_SH),   & ! (in)
                              URBAN_GRID_CARTESC_REAL_AREA(:,:),           & ! (in)
                              URBAN_GRID_CARTESC_REAL_TOTAREA              ) ! (in)
       call STATISTICS_total( UIA, UIS, UIE, UJA, UJS, UJE, &
                              URBAN_SFLX_LH  (:,:), VAR_NAME(I_SFLX_LH),   & ! (in)
                              URBAN_GRID_CARTESC_REAL_AREA(:,:),           & ! (in)
                              URBAN_GRID_CARTESC_REAL_TOTAREA              ) ! (in)
       call STATISTICS_total( UIA, UIS, UIE, UJA, UJS, UJE, &
                              URBAN_SFLX_GH  (:,:), VAR_NAME(I_SFLX_GH),   & ! (in)
                              URBAN_GRID_CARTESC_REAL_AREA(:,:),           & ! (in)
                              URBAN_GRID_CARTESC_REAL_TOTAREA              ) ! (in)
       call STATISTICS_total( UIA, UIS, UIE, UJA, UJS, UJE, &
                              URBAN_SFLX_evap(:,:), VAR_NAME(I_SFLX_evap), & ! (in)
                              URBAN_GRID_CARTESC_REAL_AREA(:,:),           & ! (in)
                              URBAN_GRID_CARTESC_REAL_TOTAREA              ) ! (in)
    endif

    return
  end subroutine URBAN_vars_total

  !-----------------------------------------------------------------------------
  !> Input from External I/O
  subroutine URBAN_vars_external_in( &
       URBAN_TC_in,        &
       URBAN_QC_in,        &
       URBAN_UC_in,        &
       URBAN_SFC_TEMP_in,  &
       URBAN_SFC_albedo_in )
    implicit none

    real(RP), intent(in) :: URBAN_TC_in        (UIA,UJA)
    real(RP), intent(in) :: URBAN_QC_in        (UIA,UJA)
    real(RP), intent(in) :: URBAN_UC_in        (UIA,UJA)
    real(RP), intent(in) :: URBAN_SFC_TEMP_in  (UIA,UJA)
    real(RP), intent(in) :: URBAN_SFC_albedo_in(UIA,UJA,2)

    integer :: k
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("URBAN_vars_external_in",*) 'External Input (coupler) '

    URBAN_TR(:,:) = URBAN_SFC_TEMP_in(:,:)
    URBAN_TB(:,:) = URBAN_SFC_TEMP_in(:,:)
    URBAN_TG(:,:) = URBAN_SFC_TEMP_in(:,:)

    URBAN_TC(:,:) = URBAN_TC_in(:,:)
    URBAN_QC(:,:) = URBAN_QC_in(:,:)
    URBAN_UC(:,:) = URBAN_UC_in(:,:)

    do k = UKS, UKE
       URBAN_TRL(k,:,:) = URBAN_SFC_TEMP_in(:,:)
       URBAN_TBL(k,:,:) = URBAN_SFC_TEMP_in(:,:)
       URBAN_TGL(k,:,:) = URBAN_SFC_TEMP_in(:,:)
    end do

    URBAN_RAINR(:,:) = 0.0_RP
    URBAN_RAINB(:,:) = 0.0_RP
    URBAN_RAING(:,:) = 0.0_RP
    URBAN_ROFF (:,:) = 0.0_RP

    URBAN_SFC_TEMP  (:,:)   = URBAN_SFC_TEMP_in  (:,:)
    URBAN_SFC_albedo(:,:,:) = URBAN_SFC_albedo_in(:,:,:)

    URBAN_Z0M      (:,:) = 2.0_RP ! tentative, will be replace in urban scheme
    URBAN_Z0H      (:,:) = 0.2_RP ! tentative, will be replace in urban scheme
    URBAN_Z0E      (:,:) = 0.2_RP ! tentative, will be replace in urban scheme
    URBAN_SFLX_MW  (:,:) = 0.0_RP
    URBAN_SFLX_MU  (:,:) = 0.0_RP
    URBAN_SFLX_MV  (:,:) = 0.0_RP
    URBAN_SFLX_SH  (:,:) = 0.0_RP
    URBAN_SFLX_LH  (:,:) = 0.0_RP
    URBAN_SFLX_GH  (:,:) = 0.0_RP
    URBAN_SFLX_evap(:,:) = 0.0_RP

    call URBAN_vars_total

    return
  end subroutine URBAN_vars_external_in

  !-----------------------------------------------------------------------------
  !> Create urban restart file
  subroutine URBAN_vars_restart_create
    use scale_time, only: &
       TIME_gettimelabel
    use scale_file_cartesC, only: &
       FILE_CARTESC_create
    use mod_urban_admin, only: &
       URBAN_do
    implicit none

    character(len=19)     :: timelabel
    character(len=H_LONG) :: basename
    !---------------------------------------------------------------------------

    if ( URBAN_do .and. URBAN_RESTART_OUT_BASENAME /= '' ) then

       LOG_NEWLINE
       LOG_INFO("URBAN_vars_restart_create",*) 'Create restart file (URBAN) '

       if ( URBAN_RESTART_OUT_POSTFIX_TIMELABEL ) then
          call TIME_gettimelabel( timelabel )
          basename = trim(URBAN_RESTART_OUT_BASENAME)//'_'//trim(timelabel)
       else
          basename = trim(URBAN_RESTART_OUT_BASENAME)
       endif

       LOG_INFO("URBAN_vars_restart_create",*) 'basename: ', trim(basename)

       call FILE_CARTESC_create( &
            basename, URBAN_RESTART_OUT_TITLE, URBAN_RESTART_OUT_DTYPE, & ! [IN]
            restart_fid,                                                & ! [OUT]
            aggregate=URBAN_RESTART_OUT_AGGREGATE                       ) ! [IN]

    endif

    return
  end subroutine URBAN_vars_restart_create

  !-----------------------------------------------------------------------------
  !> Exit netCDF define mode
  subroutine URBAN_vars_restart_enddef
    use scale_file_cartesC, only: &
       FILE_CARTESC_enddef
    implicit none

    if ( restart_fid /= -1 ) then
       call FILE_CARTESC_enddef( restart_fid ) ! [IN]
    endif

    return
  end subroutine URBAN_vars_restart_enddef

  !-----------------------------------------------------------------------------
  !> Close restart file
  subroutine URBAN_vars_restart_close
    use scale_file_cartesC, only: &
       FILE_CARTESC_close
    implicit none
    !---------------------------------------------------------------------------

    if ( restart_fid /= -1 ) then
       LOG_NEWLINE
       LOG_INFO("URBAN_vars_restart_close",*) 'Close restart file (URBAN) '

       call FILE_CARTESC_close( restart_fid ) ! [IN]

       restart_fid = -1
    endif

    return
  end subroutine URBAN_vars_restart_close

  !-----------------------------------------------------------------------------
  !> Define urban variables in restart file
  subroutine URBAN_vars_restart_def_var
    use scale_file_cartesC, only: &
       FILE_CARTESC_def_var
    implicit none
    integer :: i
    !---------------------------------------------------------------------------

    if ( restart_fid /= -1 ) then

       do i = 1, I_ALB_SW
          if ( i==I_TRL .or. i==I_TBL .or. i==I_TGL ) then ! 3D
             call FILE_CARTESC_def_var( restart_fid,     & ! [IN]
                  VAR_NAME(i), VAR_DESC(i), VAR_UNIT(i), & ! [IN]
                  'UXY', URBAN_RESTART_OUT_DTYPE,        & ! [IN]
                  VAR_ID(i)                              ) ! [OUT]
          else
             call FILE_CARTESC_def_var( restart_fid,     & ! [IN]
                  VAR_NAME(i), VAR_DESC(i), VAR_UNIT(i), & ! [IN]
                  'XY', URBAN_RESTART_OUT_DTYPE,         & ! [IN]
                  VAR_ID(i)                              ) ! [OUT]
          end if
       end do

    endif

    return
  end subroutine URBAN_vars_restart_def_var

  !-----------------------------------------------------------------------------
  !> Write urban restart
  subroutine URBAN_vars_restart_write
    use scale_file_cartesC, only: &
       FILE_CARTESC_write_var
    implicit none
    !---------------------------------------------------------------------------

    if ( restart_fid /= -1 ) then

       call URBAN_vars_total

       call FILE_CARTESC_write_var( restart_fid, VAR_ID(I_TR), URBAN_TR(:,:),                  & ! [IN]
                              VAR_NAME(I_TR), 'XY', fill_halo=.true.                        ) ! [IN]
       call FILE_CARTESC_write_var( restart_fid, VAR_ID(I_TB), URBAN_TB(:,:),                  & ! [IN]
                              VAR_NAME(I_TB), 'XY', fill_halo=.true.                        ) ! [IN]
       call FILE_CARTESC_write_var( restart_fid, VAR_ID(I_TG), URBAN_TG(:,:),                  & ! [IN]
                              VAR_NAME(I_TG), 'XY', fill_halo=.true.                        ) ! [IN]
       call FILE_CARTESC_write_var( restart_fid, VAR_ID(I_TC), URBAN_TC(:,:),                  & ! [IN]
                              VAR_NAME(I_TC), 'XY', fill_halo=.true.                        ) ! [IN]
       call FILE_CARTESC_write_var( restart_fid, VAR_ID(I_QC), URBAN_QC(:,:),                  & ! [IN]
                              VAR_NAME(I_QC), 'XY', fill_halo=.true.                        ) ! [IN]
       call FILE_CARTESC_write_var( restart_fid, VAR_ID(I_UC), URBAN_UC(:,:),                  & ! [IN]
                              VAR_NAME(I_UC), 'XY', fill_halo=.true.                        ) ! [IN]

       call FILE_CARTESC_write_var( restart_fid, VAR_ID(I_TRL), URBAN_TRL(:,:,:),              & ! [IN]
                              VAR_NAME(I_TRL), 'UXY', fill_halo=.true.                    ) ! [IN]
       call FILE_CARTESC_write_var( restart_fid, VAR_ID(I_TBL), URBAN_TBL(:,:,:),              & ! [IN]
                              VAR_NAME(I_TBL), 'UXY', fill_halo=.true.                    ) ! [IN]
       call FILE_CARTESC_write_var( restart_fid, VAR_ID(I_TGL), URBAN_TGL(:,:,:),              & ! [IN]
                              VAR_NAME(I_TGL), 'UXY', fill_halo=.true.                    ) ! [IN]

       call FILE_CARTESC_write_var( restart_fid, VAR_ID(I_RAINR), URBAN_RAINR(:,:),            & ! [IN]
                              VAR_NAME(I_RAINR), 'XY', fill_halo=.true.                     ) ! [IN]
       call FILE_CARTESC_write_var( restart_fid, VAR_ID(I_RAINB), URBAN_RAINB(:,:),            & ! [IN]
                              VAR_NAME(I_RAINB), 'XY', fill_halo=.true.                     ) ! [IN]
       call FILE_CARTESC_write_var( restart_fid, VAR_ID(I_RAING), URBAN_RAING(:,:),            & ! [IN]
                              VAR_NAME(I_RAING), 'XY', fill_halo=.true.                     ) ! [IN]
       call FILE_CARTESC_write_var( restart_fid, VAR_ID(I_ROFF), URBAN_ROFF(:,:),              & ! [IN]
                              VAR_NAME(I_ROFF),  'XY', fill_halo=.true.                     ) ! [IN]

       call FILE_CARTESC_write_var( restart_fid, VAR_ID(I_SFC_TEMP), URBAN_SFC_TEMP(:,:),      & ! [IN]
                              VAR_NAME(I_SFC_TEMP), 'XY', fill_halo=.true.                  ) ! [IN]
       call FILE_CARTESC_write_var( restart_fid, VAR_ID(I_ALB_LW), URBAN_SFC_albedo(:,:,I_LW), & ! [IN]
                              VAR_NAME(I_ALB_LW), 'XY', fill_halo=.true.                    ) ! [IN]
       call FILE_CARTESC_write_var( restart_fid, VAR_ID(I_ALB_SW), URBAN_SFC_albedo(:,:,I_SW), & ! [IN]
                              VAR_NAME(I_ALB_SW), 'XY', fill_halo=.true.                    ) ! [IN]

    endif

    return
  end subroutine URBAN_vars_restart_write

end module mod_urban_vars
