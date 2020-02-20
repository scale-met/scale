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
  use scale_io
  use scale_prof
  use scale_debug
  use scale_urban_grid_cartesC_index
  use scale_tracer
  use scale_cpl_sfc_index
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
  public :: URBAN_vars_monitor
  public :: URBAN_vars_total

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
  real(RP), public, allocatable :: URBAN_TRL  (:,:,:) ! urban temperature in layer of roof [K]
  real(RP), public, allocatable :: URBAN_TBL  (:,:,:) ! urban temperature in layer of wall [K]
  real(RP), public, allocatable :: URBAN_TGL  (:,:,:) ! urban temperature in layer of road [K]
  real(RP), public, allocatable :: URBAN_TR   (:,:)   ! urban surface temperature of roof [K]
  real(RP), public, allocatable :: URBAN_TB   (:,:)   ! urban surface temperature of wall [K]
  real(RP), public, allocatable :: URBAN_TG   (:,:)   ! urban surface temperature of road [K]
  real(RP), public, allocatable :: URBAN_TC   (:,:)   ! urban canopy air temperature [K]
  real(RP), public, allocatable :: URBAN_QC   (:,:)   ! urban canopy humidity [kg/kg]
  real(RP), public, allocatable :: URBAN_UC   (:,:)   ! urban canopy wind [m/s]
  real(RP), public, allocatable :: URBAN_RAINR(:,:)   ! urban rain storage on roof [mm=kg/m2]
  real(RP), public, allocatable :: URBAN_RAINB(:,:)   ! urban rain storage on wall [mm=kg/m2]
  real(RP), public, allocatable :: URBAN_RAING(:,:)   ! urban rain storage on road [mm=kg/m2]

  ! for restart
  real(RP), public, allocatable :: URBAN_SFC_TEMP  (:,:)     ! urban grid average of surface temperature [K]
  real(RP), public, allocatable :: URBAN_SFC_albedo(:,:,:,:) ! urban grid average of albedo (direct/diffuse,IR/near-IR/VIS) (0-1)

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

  real(RP), public, allocatable :: URBAN_ROFF (:,:)     ! urban runoff [mm/s=kg/m2/s]

  real(RP), public, allocatable :: URBAN_SFLX_MW   (:,:)     ! urban grid average of w-momentum flux [kg/m2/s]
  real(RP), public, allocatable :: URBAN_SFLX_MU   (:,:)     ! urban grid average of u-momentum flux [kg/m2/s]
  real(RP), public, allocatable :: URBAN_SFLX_MV   (:,:)     ! urban grid average of v-momentum flux [kg/m2/s]
  real(RP), public, allocatable :: URBAN_SFLX_SH   (:,:)     ! urban grid average of sensible heat flux [W/m2]
  real(RP), public, allocatable :: URBAN_SFLX_LH   (:,:)     ! urban grid average of latent heat flux [W/m2]
  real(RP), public, allocatable :: URBAN_SFLX_QTRC (:,:,:)   ! urban grid average of water vapor flux [kg/m2/s]
  real(RP), public, allocatable :: URBAN_SFLX_GH   (:,:)     ! urban grid average of ground heat flux [W/m2]

  ! given 2D variables expressing urban morphology
  real(RP), public, allocatable :: URBAN_Z0M  (:,:) ! urban grid average of rougness length (momentum) [m]
  real(RP), public, allocatable :: URBAN_Z0H  (:,:) ! urban grid average of rougness length (heat) [m]
  real(RP), public, allocatable :: URBAN_Z0E  (:,:) ! urban grid average of rougness length (vapor) [m]
  real(RP), public, allocatable :: URBAN_ZD   (:,:) ! urban grid average of displacement height [m]
  real(RP), public, allocatable :: URBAN_AH   (:,:) ! urban grid average of anthropogenic sensible heat [W/m2]
  real(RP), public, allocatable :: URBAN_AHL  (:,:) ! urban grid average of anthropogenic latent heat [W/m2]

  ! diagnostic variables
  real(RP), public, allocatable :: URBAN_Ustar(:,:) ! urban grid average of friction velocity         [m/s]
  real(RP), public, allocatable :: URBAN_Tstar(:,:) ! urban grid average of temperature scale         [K]
  real(RP), public, allocatable :: URBAN_Qstar(:,:) ! urban grid average of moisture scale            [kg/kg]
  real(RP), public, allocatable :: URBAN_Wstar(:,:) ! urban grid average of convective velocity scale [m/s]
  real(RP), public, allocatable :: URBAN_RLmo (:,:) ! urban grid average of inversed Obukhov length   [1/m]
  real(RP), public, allocatable :: URBAN_U10  (:,:) ! urban grid average of velocity u at 10m [m/s]
  real(RP), public, allocatable :: URBAN_V10  (:,:) ! urban grid average of velocity v at 10m [m/s]
  real(RP), public, allocatable :: URBAN_T2   (:,:) ! urban grid average of temperature at 2m [K]
  real(RP), public, allocatable :: URBAN_Q2   (:,:) ! urban grid average of water vapor at 2m [kg/kg]

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
  real(RP), public, allocatable :: ATMOS_SFLX_water(:,:)
  real(RP), public, allocatable :: ATMOS_SFLX_ENGI(:,:)

  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  logical,                private :: URBAN_VARS_CHECKRANGE      = .false.

  integer,                private, parameter :: VMAX              = 19
  integer,                private, parameter :: I_TRL             =  1
  integer,                private, parameter :: I_TBL             =  2
  integer,                private, parameter :: I_TGL             =  3
  integer,                private, parameter :: I_TR              =  4
  integer,                private, parameter :: I_TB              =  5
  integer,                private, parameter :: I_TG              =  6
  integer,                private, parameter :: I_TC              =  7
  integer,                private, parameter :: I_QC              =  8
  integer,                private, parameter :: I_UC              =  9
  integer,                private, parameter :: I_RAINR           = 10
  integer,                private, parameter :: I_RAINB           = 11
  integer,                private, parameter :: I_RAING           = 12
  integer,                private, parameter :: I_SFC_TEMP        = 13
  integer,                private, parameter :: I_SFC_ALB_IR_dir  = 14
  integer,                private, parameter :: I_SFC_ALB_IR_dif  = 15
  integer,                private, parameter :: I_SFC_ALB_NIR_dir = 16
  integer,                private, parameter :: I_SFC_ALB_NIR_dif = 17
  integer,                private, parameter :: I_SFC_ALB_VIS_dir = 18
  integer,                private, parameter :: I_SFC_ALB_VIS_dif = 19

  character(len=H_SHORT), private            :: VAR_NAME(VMAX) !< name  of the urban variables
  character(len=H_MID),   private            :: VAR_DESC(VMAX) !< desc. of the urban variables
  character(len=H_MID),   private            :: VAR_STDN(VMAX) !< standard name of the variables
  character(len=H_SHORT), private            :: VAR_UNIT(VMAX) !< unit  of the urban variables
  integer,                private            :: VAR_ID(VMAX)   !< ID    of the urban variables
  integer,                private            :: restart_fid = -1  ! file ID

  data VAR_NAME / 'URBAN_TRL',             &
                  'URBAN_TBL',             &
                  'URBAN_TGL',             &
                  'URBAN_TR',              &
                  'URBAN_TB',              &
                  'URBAN_TG',              &
                  'URBAN_TC',              &
                  'URBAN_QC',              &
                  'URBAN_UC',              &
                  'URBAN_RAINR',           &
                  'URBAN_RAINB',           &
                  'URBAN_RAING',           &
                  'URBAN_SFC_TEMP',        &
                  'URBAN_SFC_ALB_IR_dir',  &
                  'URBAN_SFC_ALB_IR_dif',  &
                  'URBAN_SFC_ALB_NIR_dir', &
                  'URBAN_SFC_ALB_NIR_dif', &
                  'URBAN_SFC_ALB_VIS_dir', &
                  'URBAN_SFC_ALB_VIS_dif'  /

  data VAR_DESC / 'urban temperature in layer of roof',                    &
                  'urban temperature in layer of wall',                    &
                  'urban temperature in layer of road',                    &
                  'urban surface temperature of roof',                     &
                  'urban surface temperature of wall',                     &
                  'urban surface temperature of road',                     &
                  'urban canopy air temperature',                          &
                  'urban canopy humidity',                                 &
                  'urban canopy wind',                                     &
                  'urban rain strage on roof',                             &
                  'urban rain strage on wall',                             &
                  'urban rain strage on road',                             &
                  'urban grid average of temperature',                     &
                  'urban grid average of albedo for IR (direct)',          &
                  'urban grid average of albedo for IR (diffuse)',         &
                  'urban grid average of albedo for NIR (direct)',         &
                  'urban grid average of albedo for NIR (diffuse)',        &
                  'urban grid average of albedo for VIS (direct)',         &
                  'urban grid average of albedo for VIS (diffuse)'         /

  data VAR_STDN / '', &
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
                  '', &
                  '', &
                  '', &
                  '', &
                  '', &
                  '', &
                  ''  /

  data VAR_UNIT / 'K',       &
                  'K',       &
                  'K',       &
                  'K',       &
                  'K',       &
                  'K',       &
                  'K',       &
                  'kg/kg',   &
                  'm/s',     &
                  'kg/m2',   &
                  'kg/m2',   &
                  'kg/m2',   &
                  'K',       &
                  '1',       &
                  '1',       &
                  '1',       &
                  '1',       &
                  '1',       &
                  '1'        /

  logical, private :: URBAN_RESTART_IN_CHECK_COORDINATES = .true.

  ! for monitor
  integer, parameter :: IM_TRL        = 1
  integer, parameter :: IM_TBL        = 2
  integer, parameter :: IM_TGL        = 3
  integer, parameter :: IM_TR         = 4
  integer, parameter :: IM_TB         = 5
  integer, parameter :: IM_TG         = 6
  integer, parameter :: IM_TC         = 7
  integer, parameter :: IM_UC         = 8
  integer, parameter :: IM_QC         = 9
  integer, parameter :: IM_RAINR      = 10
  integer, parameter :: IM_RAINB      = 11
  integer, parameter :: IM_RAING      = 12
  integer, parameter :: IM_ROFF       = 13
  integer, parameter :: IM_SFCWR      = 14
  integer, parameter :: IM_SFCWB      = 15
  integer, parameter :: IM_SFCWG      = 16
  integer, parameter :: IM_SFCIR      = 17
  integer, parameter :: IM_SFCIB      = 18
  integer, parameter :: IM_SFCIG      = 19
  integer, parameter :: IM_MASFLX     = 20
  integer, parameter :: IM_ENGI_S     = 21
  integer, parameter :: IM_ENGI_SR    = 22
  integer, parameter :: IM_ENGI_SB    = 23
  integer, parameter :: IM_ENGI_SG    = 24
  integer, parameter :: IM_ENGI_W     = 25
  integer, parameter :: IM_ENGI_WR    = 26
  integer, parameter :: IM_ENGI_WB    = 27
  integer, parameter :: IM_ENGI_WG    = 28
  integer, parameter :: IM_ENGSFC_GHR = 29
  integer, parameter :: IM_ENGSFC_GHB = 30
  integer, parameter :: IM_ENGSFC_GHG = 31
  integer, parameter :: IM_ENGSFC_EIR = 32
  integer, parameter :: IM_ENGSFC_EIB = 33
  integer, parameter :: IM_ENGSFC_EIG = 34
  integer, parameter :: IM_ROFF_EI    = 35
  integer, parameter :: IM_ENGFLX     = 36
  integer, parameter :: IM_max = 35
  integer, private   :: MONIT_id(IM_max)

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine URBAN_vars_setup
    use scale_prc, only: &
       PRC_abort
    use scale_const, only: &
       UNDEF => CONST_UNDEF
    use scale_monitor, only: &
       MONITOR_reg
    implicit none

    namelist / PARAM_URBAN_VARS / &
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
    LOG_INFO("URBAN_vars_setup",*) 'Setup'

    allocate( URBAN_TRL  (UKS:UKE,UIA,UJA) )
    allocate( URBAN_TBL  (UKS:UKE,UIA,UJA) )
    allocate( URBAN_TGL  (UKS:UKE,UIA,UJA) )
    allocate( URBAN_TR   (UIA,UJA)         )
    allocate( URBAN_TB   (UIA,UJA)         )
    allocate( URBAN_TG   (UIA,UJA)         )
    allocate( URBAN_TC   (UIA,UJA)         )
    allocate( URBAN_QC   (UIA,UJA)         )
    allocate( URBAN_UC   (UIA,UJA)         )
    allocate( URBAN_RAINR(UIA,UJA)         )
    allocate( URBAN_RAINB(UIA,UJA)         )
    allocate( URBAN_RAING(UIA,UJA)         )
    allocate( URBAN_ROFF (UIA,UJA)         )
    URBAN_TRL  (:,:,:) = UNDEF
    URBAN_TBL  (:,:,:) = UNDEF
    URBAN_TGL  (:,:,:) = UNDEF
    URBAN_TR   (:,:)   = UNDEF
    URBAN_TB   (:,:)   = UNDEF
    URBAN_TG   (:,:)   = UNDEF
    URBAN_TC   (:,:)   = UNDEF
    URBAN_QC   (:,:)   = UNDEF
    URBAN_UC   (:,:)   = UNDEF
    URBAN_RAINR(:,:)   = UNDEF
    URBAN_RAINB(:,:)   = UNDEF
    URBAN_RAING(:,:)   = UNDEF
    URBAN_ROFF (:,:)   = UNDEF

    allocate( URBAN_SFC_TEMP  (UIA,UJA)                     )
    allocate( URBAN_SFC_albedo(UIA,UJA,N_RAD_DIR,N_RAD_RGN) )
    URBAN_SFC_TEMP  (:,:)     = UNDEF
    URBAN_SFC_albedo(:,:,:,:) = UNDEF

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

    allocate( URBAN_SFLX_MW   (UIA,UJA)                     )
    allocate( URBAN_SFLX_MU   (UIA,UJA)                     )
    allocate( URBAN_SFLX_MV   (UIA,UJA)                     )
    allocate( URBAN_SFLX_SH   (UIA,UJA)                     )
    allocate( URBAN_SFLX_LH   (UIA,UJA)                     )
    allocate( URBAN_SFLX_GH   (UIA,UJA)                     )
    allocate( URBAN_SFLX_QTRC (UIA,UJA,QA)                  )
    URBAN_SFLX_MW   (:,:)     = UNDEF
    URBAN_SFLX_MU   (:,:)     = UNDEF
    URBAN_SFLX_MV   (:,:)     = UNDEF
    URBAN_SFLX_SH   (:,:)     = UNDEF
    URBAN_SFLX_LH   (:,:)     = UNDEF
    URBAN_SFLX_GH   (:,:)     = UNDEF
    URBAN_SFLX_QTRC (:,:,:)   = UNDEF

    allocate( URBAN_Z0M  (UIA,UJA) )
    allocate( URBAN_Z0H  (UIA,UJA) )
    allocate( URBAN_Z0E  (UIA,UJA) )
    allocate( URBAN_ZD   (UIA,UJA) )
    allocate( URBAN_AH   (UIA,UJA) )
    allocate( URBAN_AHL  (UIA,UJA) )
    allocate( URBAN_Ustar(UIA,UJA) )
    allocate( URBAN_Tstar(UIA,UJA) )
    allocate( URBAN_Qstar(UIA,UJA) )
    allocate( URBAN_Wstar(UIA,UJA) )
    allocate( URBAN_RLmo (UIA,UJA) )
    allocate( URBAN_U10  (UIA,UJA) )
    allocate( URBAN_V10  (UIA,UJA) )
    allocate( URBAN_T2   (UIA,UJA) )
    allocate( URBAN_Q2   (UIA,UJA) )
    URBAN_Z0M  (:,:) = UNDEF
    URBAN_Z0H  (:,:) = UNDEF
    URBAN_Z0E  (:,:) = UNDEF
    URBAN_ZD   (:,:) = UNDEF
    URBAN_AH   (:,:) = UNDEF
    URBAN_AHL  (:,:) = UNDEF
    URBAN_Ustar(:,:) = UNDEF
    URBAN_Tstar(:,:) = UNDEF
    URBAN_Qstar(:,:) = UNDEF
    URBAN_Wstar(:,:) = UNDEF
    URBAN_RLmo (:,:) = UNDEF
    URBAN_U10  (:,:) = UNDEF
    URBAN_V10  (:,:) = UNDEF
    URBAN_T2   (:,:) = UNDEF
    URBAN_Q2   (:,:) = UNDEF

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
    allocate( ATMOS_SFLX_water(UIA,UJA)   )
    allocate( ATMOS_SFLX_ENGI(UIA,UJA)   )
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
    ATMOS_SFLX_water(:,:)  = UNDEF
    ATMOS_SFLX_ENGI(:,:)   = UNDEF

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_URBAN_VARS,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("URBAN_vars_setup",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("URBAN_vars_setup",*) 'Not appropriate names in namelist PARAM_URBAN_VARS. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_URBAN_VARS)

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

    ! monitor
    call MONITOR_reg( 'URB_TRL',      'roof temperature',         'K m3', & ! (in)
                      MONIT_id(IM_TRL),                                   & ! (out)
                      dim_type='UXY', is_tendency=.false.                 ) ! (in)
    call MONITOR_reg( 'URB_TBL',      'wall temperature',         'K m3', & ! (in)
                      MONIT_id(IM_TBL),                                   & ! (out)
                      dim_type='UXY', is_tendency=.false.                 ) ! (in)
    call MONITOR_reg( 'URB_TGL',      'road temperature',         'K m3', & ! (in)
                      MONIT_id(IM_TGL),                                   & ! (out)
                      dim_type='UXY', is_tendency=.false.                 ) ! (in)
    call MONITOR_reg( 'URB_TR',       'roof surface temperature', 'K m2', & ! (in)
                      MONIT_id(IM_TR),                                    & ! (out)
                      dim_type='XY', is_tendency=.false.                  ) ! (in)
    call MONITOR_reg( 'URB_TB',       'wall surface temperature', 'K m2', & ! (in)
                      MONIT_id(IM_TB),                                    & ! (out)
                      dim_type='XY', is_tendency=.false.                  ) ! (in)
    call MONITOR_reg( 'URB_TG',       'road surface temperature', 'K m2', & ! (in)
                      MONIT_id(IM_TG),                                    & ! (out)
                      dim_type='XY', is_tendency=.false.                  ) ! (in)
    call MONITOR_reg( 'URB_TC',       'canopy temperature',       'K m2', & ! (in)
                      MONIT_id(IM_TC),                                    & ! (out)
                      dim_type='XY', is_tendency=.false.                  ) ! (in)
    call MONITOR_reg( 'URB_UC',       'canopy wind speed',        'm3/s', & ! (in)
                      MONIT_id(IM_UC),                                    & ! (out)
                      dim_type='XY', is_tendency=.false.                  ) ! (in)
    call MONITOR_reg( 'URB_QC',       'canopy humidity',          'kg/m', & ! (in)
                      MONIT_id(IM_QC),                                    & ! (out)
                      dim_type='XY', is_tendency=.false.                  ) ! (in)
    call MONITOR_reg( 'URB_RAINR',    'roof water',               'kg',   & ! (in)
                      MONIT_id(IM_RAINR),                                 & ! (out)
                      dim_type='XY', is_tendency=.false.                  ) ! (in)
    call MONITOR_reg( 'URB_RAINB',    'wall water',               'kg',   & ! (in)
                      MONIT_id(IM_RAINB),                                 & ! (out)
                      dim_type='XY', is_tendency=.false.                  ) ! (in)
    call MONITOR_reg( 'URB_RAING',    'road water',               'kg',   & ! (in)
                      MONIT_id(IM_RAING),                                 & ! (out)
                      dim_type='XY', is_tendency=.false.                  ) ! (in)
    call MONITOR_reg( 'URB_ROFF',     'runoff water',             'kg',   & ! (in)
                      MONIT_id(IM_ROFF),                                  & ! (out)
                      dim_type='XY', is_tendency=.false.                  ) ! (in)

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

       call FILE_CARTESC_read( restart_fid, VAR_NAME(I_TRL), 'UXY', & ! [IN]
                               URBAN_TRL(:,:,:)                     ) ! [OUT]
       call FILE_CARTESC_read( restart_fid, VAR_NAME(I_TBL), 'UXY', & ! [IN]
                               URBAN_TBL(:,:,:)                     ) ! [OUT]
       call FILE_CARTESC_read( restart_fid, VAR_NAME(I_TGL), 'UXY', & ! [IN]
                               URBAN_TGL(:,:,:)                     ) ! [OUT]

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

       call FILE_CARTESC_read( restart_fid, VAR_NAME(I_RAINR), 'XY', & ! [IN]
                               URBAN_RAINR(:,:)                      ) ! [OUT]
       call FILE_CARTESC_read( restart_fid, VAR_NAME(I_RAINB), 'XY', & ! [IN]
                               URBAN_RAINB(:,:)                      ) ! [OUT]
       call FILE_CARTESC_read( restart_fid, VAR_NAME(I_RAING), 'XY', & ! [IN]
                               URBAN_RAING(:,:)                      ) ! [OUT]

       call FILE_CARTESC_read( restart_fid, VAR_NAME(I_SFC_TEMP),        'XY',  & ! [IN]
                               URBAN_SFC_TEMP(:,:)                              ) ! [OUT]
       call FILE_CARTESC_read( restart_fid, VAR_NAME(I_SFC_ALB_IR_dir),  'XY',  & ! [IN]
                               URBAN_SFC_albedo(:,:,I_R_direct ,I_R_IR )        ) ! [OUT]
       call FILE_CARTESC_read( restart_fid, VAR_NAME(I_SFC_ALB_IR_dif),  'XY',  & ! [IN]
                               URBAN_SFC_albedo(:,:,I_R_diffuse,I_R_IR )        ) ! [OUT]
       call FILE_CARTESC_read( restart_fid, VAR_NAME(I_SFC_ALB_NIR_dir), 'XY',  & ! [IN]
                               URBAN_SFC_albedo(:,:,I_R_direct ,I_R_NIR)        ) ! [OUT]
       call FILE_CARTESC_read( restart_fid, VAR_NAME(I_SFC_ALB_NIR_dif), 'XY',  & ! [IN]
                               URBAN_SFC_albedo(:,:,I_R_diffuse,I_R_NIR)        ) ! [OUT]
       call FILE_CARTESC_read( restart_fid, VAR_NAME(I_SFC_ALB_VIS_dir), 'XY',  & ! [IN]
                               URBAN_SFC_albedo(:,:,I_R_direct ,I_R_VIS)        ) ! [OUT]
       call FILE_CARTESC_read( restart_fid, VAR_NAME(I_SFC_ALB_VIS_dif), 'XY',  & ! [IN]
                               URBAN_SFC_albedo(:,:,I_R_diffuse,I_R_VIS)        ) ! [OUT]

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
    use scale_atmos_hydrometeor, only: &
       I_QV
    use scale_landuse, only: &
       LANDUSE_exists_urban
    implicit none
    !---------------------------------------------------------------------------

    call PROF_rapstart('URB_History', 1)

    if ( URBAN_VARS_CHECKRANGE ) then
       call VALCHECK( UKA, UKS, UKE, UIA, UIS, UIE, UJA, UJS, UJE,                &
                      URBAN_TRL  (:,:,:),   0.0_RP, 1000.0_RP, VAR_NAME(I_TRL),   &
                     __FILE__, __LINE__,  mask = LANDUSE_exists_urban(:,:)        )
       call VALCHECK( UKA, UKS, UKE, UIA, UIS, UIE, UJA, UJS, UJE,                &
                      URBAN_TBL  (:,:,:),   0.0_RP, 1000.0_RP, VAR_NAME(I_TBL),   &
                     __FILE__, __LINE__,  mask = LANDUSE_exists_urban(:,:)        )
       call VALCHECK( UKA, UKS, UKE, UIA, UIS, UIE, UJA, UJS, UJE,                &
                      URBAN_TGL  (:,:,:),   0.0_RP, 1000.0_RP, VAR_NAME(I_TGL),   &
                     __FILE__, __LINE__,  mask = LANDUSE_exists_urban(:,:)        )
       call VALCHECK( UIA, UIS, UIE, UJA, UJS, UJE,                               &
                      URBAN_TR   (:,:),     0.0_RP, 1000.0_RP, VAR_NAME(I_TR),    &
                     __FILE__, __LINE__,  mask = LANDUSE_exists_urban(:,:)        )
       call VALCHECK( UIA, UIS, UIE, UJA, UJS, UJE,                               &
                      URBAN_TB   (:,:),     0.0_RP, 1000.0_RP, VAR_NAME(I_TB),    &
                     __FILE__, __LINE__,  mask = LANDUSE_exists_urban(:,:)        )
       call VALCHECK( UIA, UIS, UIE, UJA, UJS, UJE,                               &
                      URBAN_TG   (:,:),     0.0_RP, 1000.0_RP, VAR_NAME(I_TG),    &
                     __FILE__, __LINE__,  mask = LANDUSE_exists_urban(:,:)        )
       call VALCHECK( UIA, UIS, UIE, UJA, UJS, UJE,                               &
                      URBAN_TC   (:,:),     0.0_RP, 1000.0_RP, VAR_NAME(I_TC),    &
                     __FILE__, __LINE__,  mask = LANDUSE_exists_urban(:,:)        )
       call VALCHECK( UIA, UIS, UIE, UJA, UJS, UJE,                               &
                      URBAN_QC   (:,:),     0.0_RP,    1.0_RP, VAR_NAME(I_QC),    &
                     __FILE__, __LINE__,  mask = LANDUSE_exists_urban(:,:)        )
       call VALCHECK( UIA, UIS, UIE, UJA, UJS, UJE,                               &
                      URBAN_UC   (:,:),     0.0_RP,  100.0_RP, VAR_NAME(I_UC),    &
                     __FILE__, __LINE__,  mask = LANDUSE_exists_urban(:,:)        )
       call VALCHECK( UIA, UIS, UIE, UJA, UJS, UJE,                               &
                      URBAN_RAINR(:,:),     0.0_RP, 1000.0_RP, VAR_NAME(I_RAINR), &
                     __FILE__, __LINE__,  mask = LANDUSE_exists_urban(:,:)        )
       call VALCHECK( UIA, UIS, UIE, UJA, UJS, UJE,                               &
                      URBAN_RAINB(:,:),     0.0_RP, 1000.0_RP, VAR_NAME(I_RAINB), &
                     __FILE__, __LINE__,  mask = LANDUSE_exists_urban(:,:)        )
       call VALCHECK( UIA, UIS, UIE, UJA, UJS, UJE,                               &
                      URBAN_RAING(:,:),     0.0_RP, 1000.0_RP, VAR_NAME(I_RAING), &
                     __FILE__, __LINE__,  mask = LANDUSE_exists_urban(:,:)        )

       call VALCHECK( UIA, UIS, UIE, UJA, UJS, UJE,                                 &
                      URBAN_SFC_TEMP(:,:),                       0.0_RP, 1000.0_RP, &
                      VAR_NAME(I_SFC_TEMP),                     __FILE__, __LINE__, &
                      mask = LANDUSE_exists_urban(:,:)                              )

       call VALCHECK( UIA, UIS, UIE, UJA, UJS, UJE,                                 &
                      URBAN_SFC_albedo(:,:,I_R_direct ,I_R_IR ), 0.0_RP,    2.0_RP, &
                      VAR_NAME(I_SFC_ALB_IR_dir ),              __FILE__, __LINE__, &
                      mask = LANDUSE_exists_urban(:,:)                              )
       call VALCHECK( UIA, UIS, UIE, UJA, UJS, UJE,                                 &
                      URBAN_SFC_albedo(:,:,I_R_diffuse,I_R_IR ), 0.0_RP,    2.0_RP, &
                      VAR_NAME(I_SFC_ALB_IR_dif ),              __FILE__, __LINE__, &
                      mask = LANDUSE_exists_urban(:,:)                              )
       call VALCHECK( UIA, UIS, UIE, UJA, UJS, UJE,                                 &
                      URBAN_SFC_albedo(:,:,I_R_direct ,I_R_NIR), 0.0_RP,    2.0_RP, &
                      VAR_NAME(I_SFC_ALB_NIR_dir),              __FILE__, __LINE__, &
                      mask = LANDUSE_exists_urban(:,:)                              )
       call VALCHECK( UIA, UIS, UIE, UJA, UJS, UJE,                                 &
                      URBAN_SFC_albedo(:,:,I_R_diffuse,I_R_NIR), 0.0_RP,    2.0_RP, &
                      VAR_NAME(I_SFC_ALB_NIR_dif),              __FILE__, __LINE__, &
                      mask = LANDUSE_exists_urban(:,:)                              )
       call VALCHECK( UIA, UIS, UIE, UJA, UJS, UJE,                                 &
                      URBAN_SFC_albedo(:,:,I_R_direct ,I_R_VIS), 0.0_RP,    2.0_RP, &
                      VAR_NAME(I_SFC_ALB_VIS_dir),              __FILE__, __LINE__, &
                      mask = LANDUSE_exists_urban(:,:)                              )
       call VALCHECK( UIA, UIS, UIE, UJA, UJS, UJE,                                 &
                      URBAN_SFC_albedo(:,:,I_R_diffuse,I_R_VIS), 0.0_RP,    2.0_RP, &
                      VAR_NAME(I_SFC_ALB_VIS_dif),              __FILE__, __LINE__, &
                      mask = LANDUSE_exists_urban(:,:)                              )
    endif

    call FILE_HISTORY_in( URBAN_TRL(:,:,:), VAR_NAME(I_TRL), VAR_DESC(I_TRL), VAR_UNIT(I_TRL), dim_type='UXY' )
    call FILE_HISTORY_in( URBAN_TBL(:,:,:), VAR_NAME(I_TBL), VAR_DESC(I_TBL), VAR_UNIT(I_TBL), dim_type='UXY' )
    call FILE_HISTORY_in( URBAN_TGL(:,:,:), VAR_NAME(I_TGL), VAR_DESC(I_TGL), VAR_UNIT(I_TGL), dim_type='UXY' )

    call FILE_HISTORY_in( URBAN_TR(:,:), VAR_NAME(I_TR), VAR_DESC(I_TR), VAR_UNIT(I_TR) )
    call FILE_HISTORY_in( URBAN_TB(:,:), VAR_NAME(I_TB), VAR_DESC(I_TB), VAR_UNIT(I_TB) )
    call FILE_HISTORY_in( URBAN_TG(:,:), VAR_NAME(I_TG), VAR_DESC(I_TG), VAR_UNIT(I_TG) )
    call FILE_HISTORY_in( URBAN_TC(:,:), VAR_NAME(I_TC), VAR_DESC(I_TC), VAR_UNIT(I_TC) )
    call FILE_HISTORY_in( URBAN_QC(:,:), VAR_NAME(I_QC), VAR_DESC(I_QC), VAR_UNIT(I_QC) )
    call FILE_HISTORY_in( URBAN_UC(:,:), VAR_NAME(I_UC), VAR_DESC(I_UC), VAR_UNIT(I_UC) )

    call FILE_HISTORY_in( URBAN_RAINR(:,:), VAR_NAME(I_RAINR), VAR_DESC(I_RAINR), VAR_UNIT(I_RAINR) )
    call FILE_HISTORY_in( URBAN_RAINB(:,:), VAR_NAME(I_RAINB), VAR_DESC(I_RAINB), VAR_UNIT(I_RAINB) )
    call FILE_HISTORY_in( URBAN_RAING(:,:), VAR_NAME(I_RAING), VAR_DESC(I_RAING), VAR_UNIT(I_RAING) )

    call FILE_HISTORY_in( URBAN_SFC_TEMP  (:,:),                     VAR_NAME(I_SFC_TEMP),                                    &
                          VAR_DESC(I_SFC_TEMP),        VAR_UNIT(I_SFC_TEMP),        standard_name=VAR_STDN(I_SFC_TEMP)        )

    call FILE_HISTORY_in( URBAN_SFC_albedo(:,:,I_R_direct ,I_R_IR ), VAR_NAME(I_SFC_ALB_IR_dir),                              &
                          VAR_DESC(I_SFC_ALB_IR_dir),  VAR_UNIT(I_SFC_ALB_IR_dir),  standard_name=VAR_STDN(I_SFC_ALB_IR_dir)  )
    call FILE_HISTORY_in( URBAN_SFC_albedo(:,:,I_R_diffuse,I_R_IR ), VAR_NAME(I_SFC_ALB_IR_dif),                              &
                          VAR_DESC(I_SFC_ALB_IR_dif),  VAR_UNIT(I_SFC_ALB_IR_dif),  standard_name=VAR_STDN(I_SFC_ALB_IR_dif)  )
    call FILE_HISTORY_in( URBAN_SFC_albedo(:,:,I_R_direct ,I_R_NIR), VAR_NAME(I_SFC_ALB_NIR_dir),                             &
                          VAR_DESC(I_SFC_ALB_NIR_dir), VAR_UNIT(I_SFC_ALB_NIR_dir), standard_name=VAR_STDN(I_SFC_ALB_NIR_dir) )
    call FILE_HISTORY_in( URBAN_SFC_albedo(:,:,I_R_diffuse,I_R_NIR), VAR_NAME(I_SFC_ALB_NIR_dif),                             &
                          VAR_DESC(I_SFC_ALB_NIR_dif), VAR_UNIT(I_SFC_ALB_NIR_dif), standard_name=VAR_STDN(I_SFC_ALB_NIR_dif) )
    call FILE_HISTORY_in( URBAN_SFC_albedo(:,:,I_R_direct ,I_R_VIS), VAR_NAME(I_SFC_ALB_VIS_dir),                             &
                          VAR_DESC(I_SFC_ALB_VIS_dir), VAR_UNIT(I_SFC_ALB_VIS_dir), standard_name=VAR_STDN(I_SFC_ALB_VIS_dir) )
    call FILE_HISTORY_in( URBAN_SFC_albedo(:,:,I_R_diffuse,I_R_VIS), VAR_NAME(I_SFC_ALB_VIS_dif),                             &
                          VAR_DESC(I_SFC_ALB_VIS_dif), VAR_UNIT(I_SFC_ALB_VIS_dif), standard_name=VAR_STDN(I_SFC_ALB_VIS_dif) )

    call FILE_HISTORY_in( URBAN_ROFF(:,:),          'URBAN_RUNOFF',                          &
                          'urban runoff water',                                    'kg/m2/s' )

    call FILE_HISTORY_in( URBAN_SFLX_MW(:,:),       'URBAN_SFLX_MW',                         &
                          'urban grid average of w-momentum flux',                 'kg/m2/s' )
    call FILE_HISTORY_in( URBAN_SFLX_MU(:,:),       'URBAN_SFLX_MU',                         &
                          'urban grid average of u-momentum flux',                 'kg/m2/s' )
    call FILE_HISTORY_in( URBAN_SFLX_MV(:,:),       'URBAN_SFLX_MV',                         &
                          'urban grid average of v-momentum flux',                 'kg/m2/s' )
    call FILE_HISTORY_in( URBAN_SFLX_SH(:,:),       'URBAN_SFLX_SH',                         &
                          'urban grid average of sensible heat flux (upward)',     'W/m2'    )
    call FILE_HISTORY_in( URBAN_SFLX_LH(:,:),       'URBAN_SFLX_LH',                         &
                          'urban grid average of latent heat flux (upward)',       'W/m2'    )
    call FILE_HISTORY_in( URBAN_SFLX_GH(:,:),       'URBAN_SFLX_GH',                         &
                          'urban grid average of subsurface heat flux (downward)', 'W/m2'    )
    if ( I_QV > 0 ) &
    call FILE_HISTORY_in( URBAN_SFLX_QTRC(:,:,I_QV), 'URBAN_SFLX_evap',                      &
                          'urban grid average of water vapor flux (upward)',       'kg/m2/s' )

    call FILE_HISTORY_in( URBAN_Ustar(:,:), 'URBAN_Ustar',                                   &
                          'urban friction velocity',                               'm/s'     )
    call FILE_HISTORY_in( URBAN_Tstar(:,:), 'URBAN_Tstar',                                   &
                          'urban temperature scale',                               'K'       )
    call FILE_HISTORY_in( URBAN_Qstar(:,:), 'URBAN_Qstar',                                   &
                          'urban moisture scale',                                  'kg/kg'   )
    call FILE_HISTORY_in( URBAN_Wstar(:,:), 'URBAN_Wstar',                                   &
                          'urban convective velocity scale',                       'm/s'     )
    call FILE_HISTORY_in( URBAN_RLmo (:,:), 'URBAN_RLmo',                                    &
                          'urban inversed Obukhov length',                         '1/m'     )
    call FILE_HISTORY_in( URBAN_U10  (:,:), 'URBAN_U10',                                     &
                          'urban 10m x-wind',                                      'm/s'     )
    call FILE_HISTORY_in( URBAN_V10  (:,:), 'URBAN_V10',                                     &
                          'urban 10m y-wind',                                      'm/s'     )
    call FILE_HISTORY_in( URBAN_T2   (:,:), 'URBAN_T2',                                      &
                          'urban 2m temprerature',                                 'K'       )
    call FILE_HISTORY_in( URBAN_Q2   (:,:), 'URBAN_Q2',                                      &
                          'urban 2m specific humidity',                            'kg/kg'   )
    call FILE_HISTORY_in( URBAN_Z0M(:,:), 'URBAN_Z0M',                                       &
                          'urban parameter of rougness length for momentum',       'm'       )
    call FILE_HISTORY_in( URBAN_Z0H(:,:), 'URBAN_Z0H',                                       &
                          'urban parameter of rougness length for heat',           'm'       )
    call FILE_HISTORY_in( URBAN_Z0E(:,:), 'URBAN_Z0E',                                       &
                          'urban parameter of rougness length for vapor',          'm'       )
    call FILE_HISTORY_in( URBAN_ZD (:,:), 'URBAN_ZD',                                        &
                          'urban parameter of displacement height',                'm'       )
    call FILE_HISTORY_in( URBAN_AH (:,:), 'URBAN_AH',                                        &
                          'urban parameter of anthropogenic sensible heat',        'W/m2'    )
    call FILE_HISTORY_in( URBAN_AHL(:,:), 'URBAN_AHL',                                       &
                          'urban parameter of anthropogenic latent   heat',        'W/m2'    )


    call PROF_rapend  ('URB_History', 1)

    return
  end subroutine URBAN_vars_history

  !-----------------------------------------------------------------------------
  !> monitor output
  subroutine URBAN_vars_monitor
    use scale_const, only: &
       DWATR => CONST_DWATR
    use scale_atmos_hydrometeor, only: &
       CV_WATER
    use scale_monitor, only: &
       MONITOR_put
    implicit none

    real(RP)               :: WORK3D(UKA,UIA,UJA)
    real(RP)               :: WORK2D(UIA,UJA)

    integer  :: k, i, j
    !---------------------------------------------------------------------------

    call MONITOR_put( MONIT_id(IM_TRL),  URBAN_TRL(:,:,:) )
    call MONITOR_put( MONIT_id(IM_TBL),  URBAN_TBL(:,:,:) )
    call MONITOR_put( MONIT_id(IM_TGL),  URBAN_TGL(:,:,:) )

    call MONITOR_put( MONIT_id(IM_TR),  URBAN_TR(:,:) )
    call MONITOR_put( MONIT_id(IM_TB),  URBAN_TB(:,:) )
    call MONITOR_put( MONIT_id(IM_TG),  URBAN_TG(:,:) )
    call MONITOR_put( MONIT_id(IM_TC),  URBAN_TC(:,:) )
    call MONITOR_put( MONIT_id(IM_UC),  URBAN_UC(:,:) )
    if ( MONIT_id(IM_QC) > 0 ) then
       !$omp parallel do
       do j = UJS, UJE
       do i = UIS, UIE
          WORK2D(i,j) = URBAN_QC(i,j) * ATMOS_DENS(i,j)
       end do
       end do
       call MONITOR_put( MONIT_id(IM_QC), WORK2D(:,:) )
    end if

    ! mass budget
    call MONITOR_put( MONIT_id(IM_RAINR),  URBAN_RAINR(:,:) )
    call MONITOR_put( MONIT_id(IM_RAINB),  URBAN_RAINB(:,:) )
    call MONITOR_put( MONIT_id(IM_RAING),  URBAN_RAING(:,:) )
    call MONITOR_put( MONIT_id(IM_ROFF) ,  URBAN_ROFF (:,:) )
!!$    call MONITOR_put( MONIT_id(IM_SFCWR),  URBAN_SFLX_waterR(:,:) )
!!$    call MONITOR_put( MONIT_id(IM_SFCWB),  URBAN_SFLX_waterB(:,:) )
!!$    call MONITOR_put( MONIT_id(IM_SFCWG),  URBAN_SFLX_waterG(:,:) )
!!$    call MONITOR_put( MONIT_id(IM_SFCIR),  URBAN_SFLX_iceR  (:,:) )
!!$    call MONITOR_put( MONIT_id(IM_SFCIB),  URBAN_SFLX_iceB  (:,:) )
!!$    call MONITOR_put( MONIT_id(IM_SFCIG),  URBAN_SFLX_iceG  (:,:) )
!!$    if ( MONIT_id(IM_MASFLX) > 0 ) then
!!$       !$omp parallel do
!!$       do j = UJS, UJE
!!$       do i = UIS, UIE
!!$          WORK2D(i,j) = URBAN_SFLX_waterR(i,j) + URBAN_SFLX_iceR(i,j) &
!!$                      + URBAN_SFLX_waterB(i,j) + URBAN_SFLX_iceB(i,j) &
!!$                      + URBAN_SFLX_waterG(i,j) + URBAN_SFLX_iceG(i,j) &
!!$                      - URBAN_RUNOFF(i,j)
!!$       end do
!!$       end do
!!$       call MONITOR_put( MONIT_id(IM_MASFLX), WORK2D(:,:) )
!!$    end if

    ! energy budget
!!$    if ( MONIT_id(IM_ENGI_S) > 0 ) then
!!$       !$omp parallel do
!!$       do j = UJS, UJE
!!$       do i = UIS, UIE
!!$       do k = UKS, UKE
!!$          WORK3D(k,i,j) = CAPR(i,j) * URBAN_TRL(k,i,j) &
!!$                        + CAPB(i,j) * URBAN_TRB(k,i,j) &
!!$                        + CAPG(i,j) * URBAN_TRG(k,i,j)
!!$       end do
!!$       end do
!!$       end do
!!$       call MONITOR_put( MONIT_id(IM_ENGI_S), WORK3D(:,:,:) )
!!$    end if
!!$    if ( MONIT_id(IM_ENGI_SR) > 0 ) then
!!$       !$omp parallel do
!!$       do j = UJS, UJE
!!$       do i = UIS, UIE
!!$       do k = UKS, UKE
!!$          WORK3D(k,i,j) = CAPR(i,j) * URBAN_TRL(k,i,j)
!!$       end do
!!$       end do
!!$       end do
!!$       call MONITOR_put( MONIT_id(IM_ENGI_SR), WORK3D(:,:,:) )
!!$    end if
!!$    if ( MONIT_id(IM_ENGI_SB) > 0 ) then
!!$       !$omp parallel do
!!$       do j = UJS, UJE
!!$       do i = UIS, UIE
!!$       do k = UKS, UKE
!!$          WORK3D(k,i,j) = CAPB(i,j) * URBAN_TBL(k,i,j)
!!$       end do
!!$       end do
!!$       end do
!!$       call MONITOR_put( MONIT_id(IM_ENGI_SB), WORK3D(:,:,:) )
!!$    end if
!!$    if ( MONIT_id(IM_ENGI_SG) > 0 ) then
!!$       !$omp parallel do
!!$       do j = UJS, UJE
!!$       do i = UIS, UIE
!!$       do k = UKS, UKE
!!$          WORK3D(k,i,j) = CAPG(i,j) * URBAN_TGL(k,i,j)
!!$       end do
!!$       end do
!!$       end do
!!$       call MONITOR_put( MONIT_id(IM_ENGI_SG), WORK3D(:,:,:) )
!!$    end if
!!$    if ( MONIT_id(IM_ENGI_W) > 0 ) then
!!$       !$omp parallel do
!!$       do j = UJS, UJE
!!$       do i = UIS, UIE
!!$          WORK2D(i,j) = CV_WATER * ( URBAN_RAINR(k,i,j) * URBAN_TRL(UKS,i,j) &
!!$                                   + URBAN_RAINB(k,i,j) * URBAN_TBL(UKS,i,j) &
!!$                                   + URBAN_RAING(k,i,j) * URBAN_TGL(UKS,i,j) )
!!$       end do
!!$       end do
!!$       call MONITOR_put( MONIT_id(IM_ENGI_W), WORK2D(:,:) )
!!$    end if
!!$    if ( MONIT_id(IM_ENGI_WR) > 0 ) then
!!$       !$omp parallel do
!!$       do j = UJS, UJE
!!$       do i = UIS, UIE
!!$          WORK2D(i,j) = CV_WATER * URBAN_RAINR(k,i,j) * URBAN_TRL(UKS,i,j)
!!$       end do
!!$       end do
!!$       call MONITOR_put( MONIT_id(IM_ENGI_WR), WORK2D(:,:) )
!!$    end if
!!$    if ( MONIT_id(IM_ENGI_WB) > 0 ) then
!!$       !$omp parallel do
!!$       do j = UJS, UJE
!!$       do i = UIS, UIE
!!$          WORK2D(i,j) = CV_WATER * URBAN_RAINB(k,i,j) * URBAN_TBL(UKS,i,j)
!!$       end do
!!$       end do
!!$       call MONITOR_put( MONIT_id(IM_ENGI_WB), WORK2D(:,:) )
!!$    end if
!!$    if ( MONIT_id(IM_ENGI_WG) > 0 ) then
!!$       !$omp parallel do
!!$       do j = UJS, UJE
!!$       do i = UIS, UIE
!!$          WORK2D(i,j) = CV_WATER * URBAN_RAING(k,i,j) * URBAN_TGL(UKS,i,j)
!!$       end do
!!$       end do
!!$       call MONITOR_put( MONIT_id(IM_ENGI_WG), WORK2D(:,:) )
!!$    end if

!!$    call MONITOR_put( MONIT_id(IM_ENGSFC_GHR), URBAN_SFLX_GHR   (:,:) )
!!$    call MONITOR_put( MONIT_id(IM_ENGSFC_GHB), URBAN_SFLX_GHB   (:,:) )
!!$    call MONITOR_put( MONIT_id(IM_ENGSFC_GHG), URBAN_SFLX_GHG   (:,:) )
!!$    call MONITOR_put( MONIT_id(IM_ENGSFC_EIR), URBAN_SFLX_ENGIR  (:,:) )
!!$    call MONITOR_put( MONIT_id(IM_ENGSFC_EIB), URBAN_SFLX_ENGIB  (:,:) )
!!$    call MONITOR_put( MONIT_id(IM_ENGSFC_EIG), URBAN_SFLX_ENGIG  (:,:) )
!!$    call MONITOR_put( MONIT_id(IM_ROFF_EIR),   URBAN_RUNOFF_ENGIR(:,:) )
!!$    call MONITOR_put( MONIT_id(IM_ROFF_EIB),   URBAN_RUNOFF_ENGIB(:,:) )
!!$    call MONITOR_put( MONIT_id(IM_ROFF_EIG),   URBAN_RUNOFF_ENGIG(:,:) )
!!$    if ( MONIT_id(IM_ENGFLX) > 0 ) then
!!$       !$omp parallel do
!!$       do j = UJS, UJE
!!$       do i = UIS, UIE
!!$          WORK2D(i,j) = URBAN_SFLX_GHR(i,j) + URBAN_SFLX_ENGIR(i,j) &
!!$                      + URBAN_SFLX_GHB(i,j) + URBAN_SFLX_ENGIB(i,j) &
!!$                      + URBAN_SFLX_GHG(i,j) + URBAN_SFLX_ENGIG(i,j) &
!!$                      - URBAN_RUNOFF_ENGIR(i,j) &
!!$                      - URBAN_RUNOFF_ENGIB(i,j) &
!!$                      - URBAN_RUNOFF_ENGIG(i,j)
!!$       end do
!!$       end do
!!$       call MONITOR_put( MONIT_id(IM_ENGFLX), WORK2D(:,:) )
!!$    end if


    return
  end subroutine URBAN_vars_monitor

  !-----------------------------------------------------------------------------
  !> Budget monitor for urban
  subroutine URBAN_vars_total
    use scale_statistics, only: &
       STATISTICS_checktotal, &
       STATISTICS_total
    use scale_atmos_hydrometeor, only: &
       I_QV
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

       call STATISTICS_total( UIA, UIS, UIE, UJA, UJS, UJE,                                           & ! [IN]
                              URBAN_SFC_TEMP  (:,:),                     VAR_NAME(I_SFC_TEMP),        & ! [IN]
                              URBAN_GRID_CARTESC_REAL_AREA(:,:),                                      & ! [IN]
                              URBAN_GRID_CARTESC_REAL_TOTAREA                                         ) ! [IN]
       call STATISTICS_total( UIA, UIS, UIE, UJA, UJS, UJE,                                           & ! [IN]
                              URBAN_SFC_albedo(:,:,I_R_direct ,I_R_IR ), VAR_NAME(I_SFC_ALB_IR_dir),  & ! [IN]
                              URBAN_GRID_CARTESC_REAL_AREA(:,:),                                      & ! [IN]
                              URBAN_GRID_CARTESC_REAL_TOTAREA                                         ) ! [IN]
       call STATISTICS_total( UIA, UIS, UIE, UJA, UJS, UJE,                                           & ! [IN]
                              URBAN_SFC_albedo(:,:,I_R_diffuse,I_R_IR ), VAR_NAME(I_SFC_ALB_IR_dif),  & ! [IN]
                              URBAN_GRID_CARTESC_REAL_AREA(:,:),                                      & ! [IN]
                              URBAN_GRID_CARTESC_REAL_TOTAREA                                         ) ! [IN]
       call STATISTICS_total( UIA, UIS, UIE, UJA, UJS, UJE,                                           & ! [IN]
                              URBAN_SFC_albedo(:,:,I_R_direct ,I_R_NIR), VAR_NAME(I_SFC_ALB_NIR_dir), & ! [IN]
                              URBAN_GRID_CARTESC_REAL_AREA(:,:),                                      & ! [IN]
                              URBAN_GRID_CARTESC_REAL_TOTAREA                                         ) ! [IN]
       call STATISTICS_total( UIA, UIS, UIE, UJA, UJS, UJE,                                           & ! [IN]
                              URBAN_SFC_albedo(:,:,I_R_diffuse,I_R_NIR), VAR_NAME(I_SFC_ALB_NIR_dif), & ! [IN]
                              URBAN_GRID_CARTESC_REAL_AREA(:,:),                                      & ! [IN]
                              URBAN_GRID_CARTESC_REAL_TOTAREA                                         ) ! [IN]
       call STATISTICS_total( UIA, UIS, UIE, UJA, UJS, UJE,                                           & ! [IN]
                              URBAN_SFC_albedo(:,:,I_R_direct ,I_R_VIS), VAR_NAME(I_SFC_ALB_VIS_dir), & ! [IN]
                              URBAN_GRID_CARTESC_REAL_AREA(:,:),                                      & ! [IN]
                              URBAN_GRID_CARTESC_REAL_TOTAREA                                         ) ! [IN]
       call STATISTICS_total( UIA, UIS, UIE, UJA, UJS, UJE,                                           & ! [IN]
                              URBAN_SFC_albedo(:,:,I_R_diffuse,I_R_VIS), VAR_NAME(I_SFC_ALB_VIS_dif), & ! [IN]
                              URBAN_GRID_CARTESC_REAL_AREA(:,:),                                      & ! [IN]
                              URBAN_GRID_CARTESC_REAL_TOTAREA                                         ) ! [IN]
    endif

    return
  end subroutine URBAN_vars_total

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

       do i = I_TRL, I_TGL
          call FILE_CARTESC_def_var( restart_fid,     & ! [IN]
               VAR_NAME(i), VAR_DESC(i), VAR_UNIT(i), & ! [IN]
               'UXY', URBAN_RESTART_OUT_DTYPE,        & ! [IN]
               VAR_ID(i)                              ) ! [OUT]
       end do
       do i = I_TR, I_SFC_ALB_VIS_dif
          call FILE_CARTESC_def_var( restart_fid,     & ! [IN]
               VAR_NAME(i), VAR_DESC(i), VAR_UNIT(i), & ! [IN]
               'XY', URBAN_RESTART_OUT_DTYPE,         & ! [IN]
               VAR_ID(i)                              ) ! [OUT]
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

       call FILE_CARTESC_write_var( restart_fid, VAR_ID(I_TRL), URBAN_TRL(:,:,:),              & ! [IN]
                              VAR_NAME(I_TRL), 'UXY', fill_halo=.true.                    ) ! [IN]
       call FILE_CARTESC_write_var( restart_fid, VAR_ID(I_TBL), URBAN_TBL(:,:,:),              & ! [IN]
                              VAR_NAME(I_TBL), 'UXY', fill_halo=.true.                    ) ! [IN]
       call FILE_CARTESC_write_var( restart_fid, VAR_ID(I_TGL), URBAN_TGL(:,:,:),              & ! [IN]
                              VAR_NAME(I_TGL), 'UXY', fill_halo=.true.                    ) ! [IN]

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


       call FILE_CARTESC_write_var( restart_fid, VAR_ID(I_RAINR), URBAN_RAINR(:,:),            & ! [IN]
                              VAR_NAME(I_RAINR), 'XY', fill_halo=.true.                     ) ! [IN]
       call FILE_CARTESC_write_var( restart_fid, VAR_ID(I_RAINB), URBAN_RAINB(:,:),            & ! [IN]
                              VAR_NAME(I_RAINB), 'XY', fill_halo=.true.                     ) ! [IN]
       call FILE_CARTESC_write_var( restart_fid, VAR_ID(I_RAING), URBAN_RAING(:,:),            & ! [IN]
                              VAR_NAME(I_RAING), 'XY', fill_halo=.true.                     ) ! [IN]

       call FILE_CARTESC_write_var( restart_fid, VAR_ID(I_SFC_TEMP),                     & ! [IN]
                                    URBAN_SFC_TEMP(:,:),                                 & ! [IN]
                                    VAR_NAME(I_SFC_TEMP),        'XY', fill_halo=.true.  ) ! [IN]

       call FILE_CARTESC_write_var( restart_fid, VAR_ID(I_SFC_ALB_IR_dir),               & ! [IN]
                                    URBAN_SFC_albedo(:,:,I_R_direct ,I_R_IR ),           & ! [IN]
                                    VAR_NAME(I_SFC_ALB_IR_dir),  'XY',  fill_halo=.true. ) ! [IN]
       call FILE_CARTESC_write_var( restart_fid, VAR_ID(I_SFC_ALB_IR_dif),               & ! [IN]
                                    URBAN_SFC_albedo(:,:,I_R_diffuse,I_R_IR ),           & ! [IN]
                                    VAR_NAME(I_SFC_ALB_IR_dif),  'XY',  fill_halo=.true. ) ! [IN]
       call FILE_CARTESC_write_var( restart_fid, VAR_ID(I_SFC_ALB_NIR_dir),              & ! [IN]
                                    URBAN_SFC_albedo(:,:,I_R_direct ,I_R_NIR),           & ! [IN]
                                    VAR_NAME(I_SFC_ALB_NIR_dir), 'XY',  fill_halo=.true. ) ! [IN]
       call FILE_CARTESC_write_var( restart_fid, VAR_ID(I_SFC_ALB_NIR_dif),              & ! [IN]
                                    URBAN_SFC_albedo(:,:,I_R_diffuse,I_R_NIR),           & ! [IN]
                                    VAR_NAME(I_SFC_ALB_NIR_dif), 'XY',  fill_halo=.true. ) ! [IN]
       call FILE_CARTESC_write_var( restart_fid, VAR_ID(I_SFC_ALB_VIS_dir),              & ! [IN]
                                    URBAN_SFC_albedo(:,:,I_R_direct ,I_R_VIS),           & ! [IN]
                                    VAR_NAME(I_SFC_ALB_VIS_dir), 'XY',  fill_halo=.true. ) ! [IN]
       call FILE_CARTESC_write_var( restart_fid, VAR_ID(I_SFC_ALB_VIS_dif),              & ! [IN]
                                    URBAN_SFC_albedo(:,:,I_R_diffuse,I_R_VIS),           & ! [IN]
                                    VAR_NAME(I_SFC_ALB_VIS_dif), 'XY',  fill_halo=.true. ) ! [IN]

    endif

    return
  end subroutine URBAN_vars_restart_write

end module mod_urban_vars
