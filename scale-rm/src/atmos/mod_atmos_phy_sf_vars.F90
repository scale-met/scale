!-------------------------------------------------------------------------------
!> module ATMOSPHERIC Surface Variables
!!
!! @par Description
!!          Container for mod_atmos_phy_sf
!!
!! @author Team SCALE
!!
!<
!-------------------------------------------------------------------------------
#include "scalelib.h"
module mod_atmos_phy_sf_vars
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_io
  use scale_prof
  use scale_debug
  use scale_atmos_grid_cartesC_index
  use scale_tracer
  use scale_cpl_sfc_index
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: ATMOS_PHY_SF_vars_setup
  public :: ATMOS_PHY_SF_vars_finalize
  public :: ATMOS_PHY_SF_vars_fillhalo
  public :: ATMOS_PHY_SF_vars_restart_read
  public :: ATMOS_PHY_SF_vars_restart_write

  public :: ATMOS_PHY_SF_vars_restart_create
  public :: ATMOS_PHY_SF_vars_restart_open
  public :: ATMOS_PHY_SF_vars_restart_def_var
  public :: ATMOS_PHY_SF_vars_restart_enddef
  public :: ATMOS_PHY_SF_vars_restart_close

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  logical,               public :: ATMOS_PHY_SF_RESTART_OUTPUT                 = .false.                !< output restart file?

  character(len=H_LONG),  public :: ATMOS_PHY_SF_RESTART_IN_BASENAME           = ''                     !< Basename of the input  file
  logical,                public :: ATMOS_PHY_SF_RESTART_IN_AGGREGATE                                   !< Switch to use aggregate file
  logical,                public :: ATMOS_PHY_SF_RESTART_IN_POSTFIX_TIMELABEL  = .false.                !< Add timelabel to the basename of input  file?
  character(len=H_LONG),  public :: ATMOS_PHY_SF_RESTART_OUT_BASENAME          = ''                     !< Basename of the output file
  logical,                public :: ATMOS_PHY_SF_RESTART_OUT_AGGREGATE                                  !< Switch to use aggregate file
  logical,                public :: ATMOS_PHY_SF_RESTART_OUT_POSTFIX_TIMELABEL = .true.                 !< Add timelabel to the basename of output file?
  character(len=H_MID),   public :: ATMOS_PHY_SF_RESTART_OUT_TITLE             = 'ATMOS_PHY_SF restart' !< title    of the output file
  character(len=H_SHORT), public :: ATMOS_PHY_SF_RESTART_OUT_DTYPE             = 'DEFAULT'              !< REAL4 or REAL8

  real(RP), public, allocatable :: ATMOS_PHY_SF_DENS_t    (:,:)     ! tendency DENS     [    kg/m3/s]
  real(RP), public, allocatable :: ATMOS_PHY_SF_MOMZ_t    (:,:)     ! tendency MOMZ     [m/s*kg/m2/s]
  real(RP), public, allocatable :: ATMOS_PHY_SF_RHOU_t    (:,:)     ! tendency rho*U    [m/s*kg/m2/s]
  real(RP), public, allocatable :: ATMOS_PHY_SF_RHOV_t    (:,:)     ! tendency rho*V    [m/s*kg/m2/s]
  real(RP), public, allocatable :: ATMOS_PHY_SF_RHOH      (:,:)     ! diabatic heating  [J/m3/s]
  real(RP), public, allocatable :: ATMOS_PHY_SF_RHOQ_t    (:,:,:)   ! tendency rho*QTRC [    kg/kg/s]

  real(RP), public, allocatable :: ATMOS_PHY_SF_SFC_TEMP  (:,:)     ! surface skin temperature             [K]
  real(RP), public, allocatable :: ATMOS_PHY_SF_SFC_albedo(:,:,:,:) ! surface albedo (direct/diffuse,IR/near-IR/VIS) (0-1)
  real(RP), public, allocatable :: ATMOS_PHY_SF_SFC_Z0M   (:,:)     ! surface roughness length, ocean only [m]
  real(RP), public, allocatable :: ATMOS_PHY_SF_SFC_Z0H   (:,:)     ! surface roughness length, ocean only [m]
  real(RP), public, allocatable :: ATMOS_PHY_SF_SFC_Z0E   (:,:)     ! surface roughness length, ocean only [m]

  ! surface diagnostic variables
  real(RP), public, allocatable :: ATMOS_PHY_SF_SFC_DENS  (:,:)     ! surface atmosphere density  [kg/m3]
  real(RP), public, allocatable :: ATMOS_PHY_SF_SFC_PRES  (:,:)     ! surface atmosphere pressure [Pa]

  real(RP), public, allocatable :: ATMOS_PHY_SF_PREC_MASS (:,:)     ! mass flux of the precipitation [kg/m2/s]
  real(RP), public, allocatable :: ATMOS_PHY_SF_PREC_ENGI (:,:)     ! internal energy flux of the precipitation [J/m2/s]
  real(RP), public, allocatable :: ATMOS_PHY_SF_SFLX_MW   (:,:)     ! z-momentum flux (area center) [m/s*kg/m2/s]
  real(RP), public, allocatable :: ATMOS_PHY_SF_SFLX_MU   (:,:)     ! x-momentum flux (area center) [m/s*kg/m2/s]
  real(RP), public, allocatable :: ATMOS_PHY_SF_SFLX_MV   (:,:)     ! y-momentum flux (area center) [m/s*kg/m2/s]
  real(RP), public, allocatable :: ATMOS_PHY_SF_SFLX_SH   (:,:)     ! sensible heat flux [J/m2/s]
  real(RP), public, allocatable :: ATMOS_PHY_SF_SFLX_LH   (:,:)     ! latent heat flux [J/m2/s]
  real(RP), public, allocatable :: ATMOS_PHY_SF_SFLX_SHEX (:,:)     ! extra sensible heat flux [J/m2/s]
  real(RP), public, allocatable :: ATMOS_PHY_SF_SFLX_QVEX (:,:)     ! extra latent heat flux [kg/kg/m2/s]
  real(RP), public, allocatable :: ATMOS_PHY_SF_SFLX_GH   (:,:)     ! ground heat flux [J/m2/s] (downward)
  real(RP), public, allocatable, target :: ATMOS_PHY_SF_SFLX_QTRC (:,:,:) ! tracer mass flux [kg/m2/s]
  real(RP), public, allocatable :: ATMOS_PHY_SF_SFLX_ENGI (:,:)     ! internal energy flux [J/m2/s]
  real(RP), public, pointer     :: ATMOS_PHY_SF_SFLX_QV   (:,:)

  real(RP), public, allocatable :: ATMOS_PHY_SF_Ustar     (:,:)     ! friction velocity         [m/2]
  real(RP), public, allocatable :: ATMOS_PHY_SF_Tstar     (:,:)     ! temperature scale         [K]
  real(RP), public, allocatable :: ATMOS_PHY_SF_Qstar     (:,:)     ! moisture scale            [kg/kg]
  real(RP), public, allocatable :: ATMOS_PHY_SF_Wstar     (:,:)     ! convective velocity scale [m/s]

  real(RP), public, allocatable :: ATMOS_PHY_SF_U10       (:,:)     ! 10m x-wind [m/s]
  real(RP), public, allocatable :: ATMOS_PHY_SF_V10       (:,:)     ! 10m y-wind [m/s]
  real(RP), public, allocatable :: ATMOS_PHY_SF_T2        (:,:)     ! 2m temperature [K]
  real(RP), public, allocatable :: ATMOS_PHY_SF_Q2        (:,:)     ! 2m specific humidity [kg/kg]

  real(RP), public, allocatable :: ATMOS_PHY_SF_RLmo      (:,:)      ! inverse of monin-obukhov length

!  real(RP), public, allocatable :: ATMOS_PHY_SF_SFLX_QEMIS(:,:,:) ! tracer emission   flux [kg/m2/s]
!  real(RP), public, allocatable :: ATMOS_PHY_SF_SFLX_QDEP (:,:,:) ! tracer deposition flux [kg/m2/s]
!  real(RP), public, allocatable :: ATMOS_PHY_SF_SFLX_VDEP (:,:,:) ! tracer deposition velocity [m/s]

  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  integer,                private, parameter :: VMAX              = 10 !< number of the variables
  integer,                private, parameter :: I_SFC_TEMP        =  1
  integer,                private, parameter :: I_SFC_ALB_IR_dir  =  2
  integer,                private, parameter :: I_SFC_ALB_IR_dif  =  3
  integer,                private, parameter :: I_SFC_ALB_NIR_dir =  4
  integer,                private, parameter :: I_SFC_ALB_NIR_dif =  5
  integer,                private, parameter :: I_SFC_ALB_VIS_dir =  6
  integer,                private, parameter :: I_SFC_ALB_VIS_dif =  7
  integer,                private, parameter :: I_SFC_Z0M         =  8
  integer,                private, parameter :: I_SFC_Z0H         =  9
  integer,                private, parameter :: I_SFC_Z0E         = 10

  character(len=H_SHORT), private            :: VAR_NAME(VMAX) !< name  of the variables
  character(len=H_MID),   private            :: VAR_DESC(VMAX) !< desc. of the variables
  character(len=H_MID),   private            :: VAR_STDN(VMAX) !< standard name of the variables
  character(len=H_SHORT), private            :: VAR_UNIT(VMAX) !< unit  of the variables
  integer,                private            :: VAR_ID(VMAX)   !< ID    of the variables
  integer,                private            :: restart_fid = -1  ! file ID

  data VAR_NAME / 'SFC_TEMP',        &
                  'SFC_ALB_IR_dir',  &
                  'SFC_ALB_IR_dif',  &
                  'SFC_ALB_NIR_dir', &
                  'SFC_ALB_NIR_dif', &
                  'SFC_ALB_VIS_dir', &
                  'SFC_ALB_VIS_dif', &
                  'SFC_Z0M',         &
                  'SFC_Z0H',         &
                  'SFC_Z0E'          /

  data VAR_DESC / 'surface skin temperature',            &
                  'surface albedo for IR,  direct ',     &
                  'surface albedo for IR,  diffuse',     &
                  'surface albedo for NIR, direct ',     &
                  'surface albedo for NIR, diffuse',     &
                  'surface albedo for VIS, direct ',     &
                  'surface albedo for VIS, diffuse',     &
                  'surface roughness length (momentum)', &
                  'surface roughness length (heat)',     &
                  'surface roughness length (vapor)'     /

  data VAR_STDN / 'surface_temp', &
                  '', &
                  '', &
                  '', &
                  '', &
                  '', &
                  '', &
                  'surface_roughness_length_for_momentum_in_air', &
                  'surface_roughness_length_for_heat_in_air', &
                  ''  /

  data VAR_UNIT / 'K', &
                  '1', &
                  '1', &
                  '1', &
                  '1', &
                  '1', &
                  '1', &
                  'm', &
                  'm', &
                  'm'  /

  real(RP), private :: ATMOS_PHY_SF_DEFAULT_SFC_TEMP      = 300.0_RP
  real(RP), private :: ATMOS_PHY_SF_DEFAULT_SFC_albedo_LW = 0.04_RP
  real(RP), private :: ATMOS_PHY_SF_DEFAULT_SFC_albedo_SW = 0.10_RP
  real(RP), private :: ATMOS_PHY_SF_DEFAULT_SFC_Z0M       = 1E-4_RP
  real(RP), private :: ATMOS_PHY_SF_DEFAULT_SFC_Z0H       = 1E-5_RP
  real(RP), private :: ATMOS_PHY_SF_DEFAULT_SFC_Z0E       = 1E-5_RP

  real(RP), allocatable, target :: ZERO(:,:)

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine ATMOS_PHY_SF_vars_setup
    use scale_prc, only: &
       PRC_abort
    use scale_const, only: &
       UNDEF => CONST_UNDEF
    use scale_atmos_hydrometeor, only: &
       I_QV
    implicit none

    namelist / PARAM_ATMOS_PHY_SF_VARS / &
       ATMOS_PHY_SF_RESTART_IN_BASENAME,           &
       ATMOS_PHY_SF_RESTART_IN_AGGREGATE,          &
       ATMOS_PHY_SF_RESTART_IN_POSTFIX_TIMELABEL,  &
       ATMOS_PHY_SF_RESTART_OUTPUT,                &
       ATMOS_PHY_SF_RESTART_OUT_BASENAME,          &
       ATMOS_PHY_SF_RESTART_OUT_AGGREGATE,         &
       ATMOS_PHY_SF_RESTART_OUT_POSTFIX_TIMELABEL, &
       ATMOS_PHY_SF_RESTART_OUT_TITLE,             &
       ATMOS_PHY_SF_RESTART_OUT_DTYPE,             &
       ATMOS_PHY_SF_DEFAULT_SFC_TEMP,              &
       ATMOS_PHY_SF_DEFAULT_SFC_albedo_LW,         &
       ATMOS_PHY_SF_DEFAULT_SFC_albedo_SW,         &
       ATMOS_PHY_SF_DEFAULT_SFC_Z0M,               &
       ATMOS_PHY_SF_DEFAULT_SFC_Z0H,               &
       ATMOS_PHY_SF_DEFAULT_SFC_Z0E

    integer :: ierr
    integer :: iv
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("ATMOS_PHY_SF_vars_setup",*) 'Setup'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_ATMOS_PHY_SF_VARS,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("ATMOS_PHY_SF_vars_setup",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("ATMOS_PHY_SF_vars_setup",*) 'Not appropriate names in namelist PARAM_ATMOS_PHY_SF_VARS. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_ATMOS_PHY_SF_VARS)


    allocate( ATMOS_PHY_SF_DENS_t    (IA,JA)    )
    allocate( ATMOS_PHY_SF_MOMZ_t    (IA,JA)    )
    allocate( ATMOS_PHY_SF_RHOU_t    (IA,JA)    )
    allocate( ATMOS_PHY_SF_RHOV_t    (IA,JA)    )
    allocate( ATMOS_PHY_SF_RHOH      (IA,JA)    )
    allocate( ATMOS_PHY_SF_RHOQ_t    (IA,JA,QA) )
    ATMOS_PHY_SF_DENS_t    (:,:)     = UNDEF
    ATMOS_PHY_SF_MOMZ_t    (:,:)     = UNDEF
    ATMOS_PHY_SF_RHOU_t    (:,:)     = UNDEF
    ATMOS_PHY_SF_RHOV_t    (:,:)     = UNDEF
    ATMOS_PHY_SF_RHOH      (:,:)     = UNDEF
    ATMOS_PHY_SF_RHOQ_t    (:,:,:)   = UNDEF

    allocate( ATMOS_PHY_SF_SFC_TEMP  (IA,JA)                     )
    allocate( ATMOS_PHY_SF_SFC_albedo(IA,JA,N_RAD_DIR,N_RAD_RGN) )
    allocate( ATMOS_PHY_SF_SFC_Z0M   (IA,JA)                     )
    allocate( ATMOS_PHY_SF_SFC_Z0H   (IA,JA)                     )
    allocate( ATMOS_PHY_SF_SFC_Z0E   (IA,JA)                     )
    ATMOS_PHY_SF_SFC_TEMP(:,:) = ATMOS_PHY_SF_DEFAULT_SFC_TEMP
    ATMOS_PHY_SF_SFC_albedo(:,:,I_R_direct ,I_R_IR ) = ATMOS_PHY_SF_DEFAULT_SFC_albedo_LW
    ATMOS_PHY_SF_SFC_albedo(:,:,I_R_diffuse,I_R_IR ) = ATMOS_PHY_SF_DEFAULT_SFC_albedo_LW
    ATMOS_PHY_SF_SFC_albedo(:,:,I_R_direct ,I_R_NIR) = ATMOS_PHY_SF_DEFAULT_SFC_albedo_SW
    ATMOS_PHY_SF_SFC_albedo(:,:,I_R_diffuse,I_R_NIR) = ATMOS_PHY_SF_DEFAULT_SFC_albedo_SW
    ATMOS_PHY_SF_SFC_albedo(:,:,I_R_direct ,I_R_VIS) = ATMOS_PHY_SF_DEFAULT_SFC_albedo_SW
    ATMOS_PHY_SF_SFC_albedo(:,:,I_R_diffuse,I_R_VIS) = ATMOS_PHY_SF_DEFAULT_SFC_albedo_SW
    ATMOS_PHY_SF_SFC_Z0M (:,:) = ATMOS_PHY_SF_DEFAULT_SFC_Z0M
    ATMOS_PHY_SF_SFC_Z0H (:,:) = ATMOS_PHY_SF_DEFAULT_SFC_Z0H
    ATMOS_PHY_SF_SFC_Z0E (:,:) = ATMOS_PHY_SF_DEFAULT_SFC_Z0E

    allocate( ATMOS_PHY_SF_SFC_DENS  (IA,JA) )
    allocate( ATMOS_PHY_SF_SFC_PRES  (IA,JA) )
    ATMOS_PHY_SF_SFC_DENS  (:,:)     = UNDEF
    ATMOS_PHY_SF_SFC_PRES  (:,:)     = UNDEF

    allocate( ATMOS_PHY_SF_PREC_MASS (IA,JA) )
    allocate( ATMOS_PHY_SF_PREC_ENGI (IA,JA) )
    ATMOS_PHY_SF_PREC_MASS(:,:)      = UNDEF
    ATMOS_PHY_SF_PREC_ENGI(:,:)      = UNDEF

    allocate( ATMOS_PHY_SF_SFLX_MW   (IA,JA) )
    allocate( ATMOS_PHY_SF_SFLX_MU   (IA,JA) )
    allocate( ATMOS_PHY_SF_SFLX_MV   (IA,JA) )
    allocate( ATMOS_PHY_SF_SFLX_SH   (IA,JA) )
    allocate( ATMOS_PHY_SF_SFLX_LH   (IA,JA) )
    allocate( ATMOS_PHY_SF_SFLX_SHEX (IA,JA) )
    allocate( ATMOS_PHY_SF_SFLX_QVEX (IA,JA) )
    allocate( ATMOS_PHY_SF_SFLX_GH   (IA,JA) )
    allocate( ATMOS_PHY_SF_SFLX_QTRC (IA,JA,max(QA,1)) )
    allocate( ATMOS_PHY_SF_SFLX_ENGI (IA,JA) )
    ATMOS_PHY_SF_SFLX_MW   (:,:)     = UNDEF
    ATMOS_PHY_SF_SFLX_MU   (:,:)     = UNDEF
    ATMOS_PHY_SF_SFLX_MV   (:,:)     = UNDEF
    ATMOS_PHY_SF_SFLX_SH   (:,:)     = UNDEF
    ATMOS_PHY_SF_SFLX_LH   (:,:)     = UNDEF
    ATMOS_PHY_SF_SFLX_SHEX (:,:)     = UNDEF
    ATMOS_PHY_SF_SFLX_QVEX (:,:)     = UNDEF
    ATMOS_PHY_SF_SFLX_GH   (:,:)     = UNDEF
    ATMOS_PHY_SF_SFLX_QTRC (:,:,:)   = UNDEF
    ATMOS_PHY_SF_SFLX_ENGI (:,:)     = UNDEF

    allocate( ATMOS_PHY_SF_Ustar     (IA,JA) )
    allocate( ATMOS_PHY_SF_Tstar     (IA,JA) )
    allocate( ATMOS_PHY_SF_Qstar     (IA,JA) )
    allocate( ATMOS_PHY_SF_Wstar     (IA,JA) )
    ATMOS_PHY_SF_Ustar     (:,:)     = UNDEF
    ATMOS_PHY_SF_Tstar     (:,:)     = UNDEF
    ATMOS_PHY_SF_Wstar     (:,:)     = UNDEF

    allocate( ATMOS_PHY_SF_U10       (IA,JA) )
    allocate( ATMOS_PHY_SF_V10       (IA,JA) )
    allocate( ATMOS_PHY_SF_T2        (IA,JA) )
    allocate( ATMOS_PHY_SF_Q2        (IA,JA) )
    allocate( ATMOS_PHY_SF_RLmo      (IA,JA) )
    ATMOS_PHY_SF_U10       (:,:)     = UNDEF
    ATMOS_PHY_SF_V10       (:,:)     = UNDEF
    ATMOS_PHY_SF_T2        (:,:)     = UNDEF
    ATMOS_PHY_SF_Q2        (:,:)     = UNDEF
    ATMOS_PHY_SF_RLmo      (:,:)     = UNDEF


    LOG_NEWLINE
    LOG_INFO("ATMOS_PHY_SF_vars_setup",*) '[ATMOS_PHY_SF] prognostic/diagnostic variables'
    LOG_INFO_CONT('(1x,A,A24,A,A48,A,A12,A)') &
               '      |', 'VARNAME                 ','|', &
               'DESCRIPTION                                     ', '[', 'UNIT        ', ']'
    do iv = 1, VMAX
       LOG_INFO_CONT('(1x,A,I3,A,A24,A,A48,A,A12,A)') &
                  'NO.',iv,'|',VAR_NAME(iv),'|',VAR_DESC(iv),'[',VAR_UNIT(iv),']'
    enddo

    LOG_NEWLINE
    if ( ATMOS_PHY_SF_RESTART_IN_BASENAME /= '' ) then
       LOG_INFO("ATMOS_PHY_SF_vars_setup",*) 'Restart input?  : YES, file = ', trim(ATMOS_PHY_SF_RESTART_IN_BASENAME)
       LOG_INFO("ATMOS_PHY_SF_vars_setup",*) 'Add timelabel?  : ', ATMOS_PHY_SF_RESTART_IN_POSTFIX_TIMELABEL
    else
       LOG_INFO("ATMOS_PHY_SF_vars_setup",*) 'Restart input?  : NO'
    endif
    if (       ATMOS_PHY_SF_RESTART_OUTPUT             &
         .AND. ATMOS_PHY_SF_RESTART_OUT_BASENAME /= '' ) then
       LOG_INFO("ATMOS_PHY_SF_vars_setup",*) 'Restart output? : YES, file = ', trim(ATMOS_PHY_SF_RESTART_OUT_BASENAME)
       LOG_INFO("ATMOS_PHY_SF_vars_setup",*) 'Add timelabel?  : ', ATMOS_PHY_SF_RESTART_OUT_POSTFIX_TIMELABEL
    else
       LOG_INFO("ATMOS_PHY_SF_vars_setup",*) 'Restart output? : NO'
       ATMOS_PHY_SF_RESTART_OUTPUT = .false.
    endif

    if ( I_QV > 0 ) then
       ATMOS_PHY_SF_SFLX_QV => ATMOS_PHY_SF_SFLX_QTRC(:,:,I_QV)
    else
       allocate( ZERO(IA,JA) )
       ZERO(:,:) = 0.0_RP
       ATMOS_PHY_SF_SFLX_QV => ZERO
    end if

    return
  end subroutine ATMOS_PHY_SF_vars_setup

  !-----------------------------------------------------------------------------
  !> Finalize
  subroutine ATMOS_PHY_SF_vars_finalize
    implicit none
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("ATMOS_PHY_SF_vars_finalize",*) 'Finalize'

    deallocate( ATMOS_PHY_SF_DENS_t        )
    deallocate( ATMOS_PHY_SF_MOMZ_t        )
    deallocate( ATMOS_PHY_SF_RHOU_t        )
    deallocate( ATMOS_PHY_SF_RHOV_t        )
    deallocate( ATMOS_PHY_SF_RHOH          )
    deallocate( ATMOS_PHY_SF_RHOQ_t    )

    deallocate( ATMOS_PHY_SF_SFC_TEMP   )
    deallocate( ATMOS_PHY_SF_SFC_albedo )
    deallocate( ATMOS_PHY_SF_SFC_Z0M    )
    deallocate( ATMOS_PHY_SF_SFC_Z0H    )
    deallocate( ATMOS_PHY_SF_SFC_Z0E    )

    deallocate( ATMOS_PHY_SF_SFC_DENS   )
    deallocate( ATMOS_PHY_SF_SFC_PRES   )

    deallocate( ATMOS_PHY_SF_PREC_MASS  )
    deallocate( ATMOS_PHY_SF_PREC_ENGI  )


    deallocate( ATMOS_PHY_SF_SFLX_MW    )
    deallocate( ATMOS_PHY_SF_SFLX_MU    )
    deallocate( ATMOS_PHY_SF_SFLX_MV    )
    deallocate( ATMOS_PHY_SF_SFLX_SH    )
    deallocate( ATMOS_PHY_SF_SFLX_LH    )
    deallocate( ATMOS_PHY_SF_SFLX_SHEX  )
    deallocate( ATMOS_PHY_SF_SFLX_QVEX  )
    deallocate( ATMOS_PHY_SF_SFLX_GH    )
    deallocate( ATMOS_PHY_SF_SFLX_QTRC  )
    deallocate( ATMOS_PHY_SF_SFLX_ENGI  )

    deallocate( ATMOS_PHY_SF_Ustar      )
    deallocate( ATMOS_PHY_SF_Tstar      )
    deallocate( ATMOS_PHY_SF_Qstar      )
    deallocate( ATMOS_PHY_SF_Wstar      )

    deallocate( ATMOS_PHY_SF_U10        )
    deallocate( ATMOS_PHY_SF_V10        )
    deallocate( ATMOS_PHY_SF_T2         )
    deallocate( ATMOS_PHY_SF_Q2         )
    deallocate( ATMOS_PHY_SF_RLmo       )

    if ( allocated( ZERO ) ) deallocate( ZERO )
    return
  end subroutine ATMOS_PHY_SF_vars_finalize

  !-----------------------------------------------------------------------------
  !> HALO Communication
  subroutine ATMOS_PHY_SF_vars_fillhalo
    use scale_comm_cartesC, only: &
       COMM_vars8, &
       COMM_wait
    implicit none

    integer :: n ,idir, irgn
    !---------------------------------------------------------------------------

    call COMM_vars8( ATMOS_PHY_SF_SFC_TEMP (:,:), 1 )
    call COMM_vars8( ATMOS_PHY_SF_SFC_Z0M  (:,:), 2 )
    call COMM_vars8( ATMOS_PHY_SF_SFC_Z0H  (:,:), 3 )
    call COMM_vars8( ATMOS_PHY_SF_SFC_Z0E  (:,:), 4 )

    n = 4
    do irgn = I_R_IR, I_R_VIS
    do idir = I_R_direct, I_R_diffuse
       n = n + 1
       call COMM_vars8( ATMOS_PHY_SF_SFC_albedo(:,:,idir,irgn), n )
    enddo
    enddo

    call COMM_wait ( ATMOS_PHY_SF_SFC_TEMP (:,:), 1 )
    call COMM_wait ( ATMOS_PHY_SF_SFC_Z0M  (:,:), 2 )
    call COMM_wait ( ATMOS_PHY_SF_SFC_Z0H  (:,:), 3 )
    call COMM_wait ( ATMOS_PHY_SF_SFC_Z0E  (:,:), 4 )

    n = 4
    do irgn = I_R_IR, I_R_VIS
    do idir = I_R_direct, I_R_diffuse
       n = n + 1
       call COMM_wait ( ATMOS_PHY_SF_SFC_albedo(:,:,idir,irgn), n )
    enddo
    enddo

    return
  end subroutine ATMOS_PHY_SF_vars_fillhalo

  !-----------------------------------------------------------------------------
  !> Open restart file for read
  subroutine ATMOS_PHY_SF_vars_restart_open
    use scale_time, only: &
       TIME_gettimelabel
    use scale_file_cartesC, only: &
       FILE_CARTESC_open
    implicit none

    character(len=19)     :: timelabel
    character(len=H_LONG) :: basename
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("ATMOS_PHY_SF_vars_restart_open",*) 'Open restart file (ATMOS_PHY_SF) '

    if ( ATMOS_PHY_SF_RESTART_IN_BASENAME /= '' ) then

       if ( ATMOS_PHY_SF_RESTART_IN_POSTFIX_TIMELABEL ) then
          call TIME_gettimelabel( timelabel )
          basename = trim(ATMOS_PHY_SF_RESTART_IN_BASENAME)//'_'//trim(timelabel)
       else
          basename = trim(ATMOS_PHY_SF_RESTART_IN_BASENAME)
       endif

       LOG_INFO("ATMOS_PHY_SF_vars_restart_open",*) 'basename: ', trim(basename)

       call FILE_CARTESC_open( basename, restart_fid, aggregate=ATMOS_PHY_SF_RESTART_IN_AGGREGATE )

    else
       LOG_INFO("ATMOS_PHY_SF_vars_restart_open",*) 'restart file for ATMOS_PHY_SF is not specified.'
    endif

    return
  end subroutine ATMOS_PHY_SF_vars_restart_open

  !-----------------------------------------------------------------------------
  !> Read restart
  subroutine ATMOS_PHY_SF_vars_restart_read
    use scale_file, only: &
       FILE_get_aggregate
    use scale_file_cartesC, only: &
       FILE_CARTESC_read, &
       FILE_CARTESC_flush
    implicit none
    !---------------------------------------------------------------------------

    if ( restart_fid /= -1 ) then
       LOG_NEWLINE
       LOG_INFO("ATMOS_PHY_SF_vars_restart_read",*) 'Read from restart file (ATMOS_PHY_SF) '

       call FILE_CARTESC_read( restart_fid, VAR_NAME(I_SFC_TEMP),        'XY',  & ! [IN]
                               ATMOS_PHY_SF_SFC_TEMP  (:,:)                     ) ! [OUT]
       call FILE_CARTESC_read( restart_fid, VAR_NAME(I_SFC_ALB_IR_dir),  'XY',  & ! [IN]
                               ATMOS_PHY_SF_SFC_albedo(:,:,I_R_direct ,I_R_IR ) ) ! [OUT]
       call FILE_CARTESC_read( restart_fid, VAR_NAME(I_SFC_ALB_IR_dif),  'XY',  & ! [IN]
                               ATMOS_PHY_SF_SFC_albedo(:,:,I_R_diffuse,I_R_IR ) ) ! [OUT]
       call FILE_CARTESC_read( restart_fid, VAR_NAME(I_SFC_ALB_NIR_dir), 'XY',  & ! [IN]
                               ATMOS_PHY_SF_SFC_albedo(:,:,I_R_direct ,I_R_NIR) ) ! [OUT]
       call FILE_CARTESC_read( restart_fid, VAR_NAME(I_SFC_ALB_NIR_dif), 'XY',  & ! [IN]
                               ATMOS_PHY_SF_SFC_albedo(:,:,I_R_diffuse,I_R_NIR) ) ! [OUT]
       call FILE_CARTESC_read( restart_fid, VAR_NAME(I_SFC_ALB_VIS_dir), 'XY',  & ! [IN]
                               ATMOS_PHY_SF_SFC_albedo(:,:,I_R_direct ,I_R_VIS) ) ! [OUT]
       call FILE_CARTESC_read( restart_fid, VAR_NAME(I_SFC_ALB_VIS_dif), 'XY',  & ! [IN]
                               ATMOS_PHY_SF_SFC_albedo(:,:,I_R_diffuse,I_R_VIS) ) ! [OUT]
       call FILE_CARTESC_read( restart_fid, VAR_NAME(I_SFC_Z0M),         'XY',  & ! [IN]
                               ATMOS_PHY_SF_SFC_Z0M   (:,:)                     ) ! [OUT]
       call FILE_CARTESC_read( restart_fid, VAR_NAME(I_SFC_Z0H),         'XY',  & ! [IN]
                               ATMOS_PHY_SF_SFC_Z0H   (:,:)                     ) ! [OUT]
       call FILE_CARTESC_read( restart_fid, VAR_NAME(I_SFC_Z0E),         'XY',  & ! [IN]
                               ATMOS_PHY_SF_SFC_Z0E   (:,:)                     ) ! [OUT]

       if ( FILE_get_AGGREGATE(restart_fid) ) then
          call FILE_CARTESC_flush( restart_fid ) ! X/Y halos have been read from file
       else
          call ATMOS_PHY_SF_vars_fillhalo
       end if

       call ATMOS_PHY_SF_vars_check

    else
       LOG_INFO("ATMOS_PHY_SF_vars_restart_read",*) 'invalid restart file ID for ATMOS_PHY_SF.'
    endif

    return
  end subroutine ATMOS_PHY_SF_vars_restart_read

  !-----------------------------------------------------------------------------
  !> Create restart file
  subroutine ATMOS_PHY_SF_vars_restart_create
    use scale_time, only: &
       TIME_gettimelabel
    use scale_file_cartesC, only: &
       FILE_CARTESC_create
    implicit none

    character(len=19)     :: timelabel
    character(len=H_LONG) :: basename
    !---------------------------------------------------------------------------

    if ( ATMOS_PHY_SF_RESTART_OUT_BASENAME /= '' ) then

       LOG_NEWLINE
       LOG_INFO("ATMOS_PHY_SF_vars_restart_create",*) 'Create restart file (ATMOS_PHY_AE) '

       if ( ATMOS_PHY_SF_RESTART_OUT_POSTFIX_TIMELABEL ) then
          call TIME_gettimelabel( timelabel )
          basename = trim(ATMOS_PHY_SF_RESTART_OUT_BASENAME)//'_'//trim(timelabel)
       else
          basename = trim(ATMOS_PHY_SF_RESTART_OUT_BASENAME)
       endif

       LOG_INFO("ATMOS_PHY_SF_vars_restart_create",*) 'basename: ', trim(basename)

       call FILE_CARTESC_create( &
            basename, ATMOS_PHY_SF_RESTART_OUT_TITLE, ATMOS_PHY_SF_RESTART_OUT_DTYPE, & ! [IN]
            restart_fid,                                                              & ! [OUT]
            aggregate=ATMOS_PHY_SF_RESTART_OUT_AGGREGATE                              ) ! [IN]
    endif

    return
  end subroutine ATMOS_PHY_SF_vars_restart_create

  !-----------------------------------------------------------------------------
  !> Exit netCDF define mode
  subroutine ATMOS_PHY_SF_vars_restart_enddef
    use scale_file_cartesC, only: &
       FILE_CARTESC_enddef
    implicit none

    if ( restart_fid /= -1 ) then
       call FILE_CARTESC_enddef( restart_fid ) ! [IN]
    endif

    return
  end subroutine ATMOS_PHY_SF_vars_restart_enddef

  !-----------------------------------------------------------------------------
  !> Close restart file
  subroutine ATMOS_PHY_SF_vars_restart_close
    use scale_file_cartesC, only: &
       FILE_CARTESC_close
    implicit none
    !---------------------------------------------------------------------------

    if ( restart_fid /= -1 ) then
       LOG_NEWLINE
       LOG_INFO("ATMOS_PHY_SF_vars_restart_close",*) 'Close restart file (ATMOS_PHY_SF) '

       call FILE_CARTESC_close( restart_fid ) ! [IN]

       restart_fid = -1
    endif

    return
  end subroutine ATMOS_PHY_SF_vars_restart_close

  !-----------------------------------------------------------------------------
  !> Write restart
  subroutine ATMOS_PHY_SF_vars_restart_def_var
    use scale_file_cartesC, only: &
       FILE_CARTESC_def_var
    implicit none

    integer :: i
    !---------------------------------------------------------------------------

    if ( restart_fid /= -1 ) then

       do i = 1, VMAX
          call FILE_CARTESC_def_var( restart_fid,                           & ! [IN]
                                     VAR_NAME(i), VAR_DESC(i), VAR_UNIT(i), & ! [IN]
                                     'XY', ATMOS_PHY_SF_RESTART_OUT_DTYPE,  & ! [IN]
                                     VAR_ID(i),                             & ! [OUT]
                                     standard_name=VAR_STDN(i)              ) ! [IN]
       end do

    endif

    return
  end subroutine ATMOS_PHY_SF_vars_restart_def_var

  !-----------------------------------------------------------------------------
  !> Write variables to restart file
  subroutine ATMOS_PHY_SF_vars_restart_write
    use scale_file_cartesC, only: &
       FILE_CARTESC_write_var
    implicit none
    !---------------------------------------------------------------------------

    if ( restart_fid /= -1 ) then

       call ATMOS_PHY_SF_vars_fillhalo

       call ATMOS_PHY_SF_vars_check

       call FILE_CARTESC_write_var( restart_fid, VAR_ID(I_SFC_TEMP),                  & ! [IN]
                                    ATMOS_PHY_SF_SFC_TEMP  (:,:),                     & ! [IN]
                                    VAR_NAME(I_SFC_TEMP), 'XY'                        ) ! [IN]
       call FILE_CARTESC_write_var( restart_fid, VAR_ID(I_SFC_ALB_IR_dir),            & ! [IN]
                                    ATMOS_PHY_SF_SFC_albedo(:,:,I_R_direct ,I_R_IR ), & ! [IN]
                                    VAR_NAME(I_SFC_ALB_IR_dir), 'XY'                  ) ! [IN]
       call FILE_CARTESC_write_var( restart_fid, VAR_ID(I_SFC_ALB_IR_dif),            & ! [IN]
                                    ATMOS_PHY_SF_SFC_albedo(:,:,I_R_diffuse,I_R_IR ), & ! [IN]
                                    VAR_NAME(I_SFC_ALB_IR_dif), 'XY'                  ) ! [IN]
       call FILE_CARTESC_write_var( restart_fid, VAR_ID(I_SFC_ALB_NIR_dir),           & ! [IN]
                                    ATMOS_PHY_SF_SFC_albedo(:,:,I_R_direct ,I_R_NIR), & ! [IN]
                                    VAR_NAME(I_SFC_ALB_NIR_dir), 'XY'                 ) ! [IN]
       call FILE_CARTESC_write_var( restart_fid, VAR_ID(I_SFC_ALB_NIR_dif),           & ! [IN]
                                    ATMOS_PHY_SF_SFC_albedo(:,:,I_R_diffuse,I_R_NIR), & ! [IN]
                                    VAR_NAME(I_SFC_ALB_NIR_dif), 'XY'                 ) ! [IN]
       call FILE_CARTESC_write_var( restart_fid, VAR_ID(I_SFC_ALB_VIS_dir),           & ! [IN]
                                    ATMOS_PHY_SF_SFC_albedo(:,:,I_R_direct ,I_R_VIS), & ! [IN]
                                    VAR_NAME(I_SFC_ALB_VIS_dir), 'XY'                 ) ! [IN]
       call FILE_CARTESC_write_var( restart_fid, VAR_ID(I_SFC_ALB_VIS_dif),           & ! [IN]
                                    ATMOS_PHY_SF_SFC_albedo(:,:,I_R_diffuse,I_R_VIS), & ! [IN]
                                    VAR_NAME(I_SFC_ALB_VIS_dif), 'XY'                 ) ! [IN]
       call FILE_CARTESC_write_var( restart_fid, VAR_ID(I_SFC_Z0M),                   & ! [IN]
                                    ATMOS_PHY_SF_SFC_Z0M   (:,:),                     & ! [IN]
                                    VAR_NAME(I_SFC_Z0M), 'XY'                         ) ! [IN]
       call FILE_CARTESC_write_var( restart_fid, VAR_ID(I_SFC_Z0H),                   & ! [IN]
                                    ATMOS_PHY_SF_SFC_Z0H   (:,:),                     & ! [IN]
                                    VAR_NAME(I_SFC_Z0H), 'XY'                         ) ! [IN]
       call FILE_CARTESC_write_var( restart_fid, VAR_ID(I_SFC_Z0E),                   & ! [IN]
                                    ATMOS_PHY_SF_SFC_Z0E   (:,:),                     & ! [IN]
                                    VAR_NAME(I_SFC_Z0E), 'XY'                         ) ! [IN]
    endif

    return
  end subroutine ATMOS_PHY_SF_vars_restart_write

  !-----------------------------------------------------------------------------
  subroutine ATMOS_PHY_SF_vars_check
    use scale_statistics, only: &
       STATISTICS_total
    use scale_atmos_grid_cartesC_real, only: &
       ATMOS_GRID_CARTESC_REAL_AREA, &
       ATMOS_GRID_CARTESC_REAL_TOTAREA
    implicit none
    !---------------------------------------------------------------------------

    call VALCHECK( IA, IS, IE, JA, JS, JE, &
                   ATMOS_PHY_SF_SFC_TEMP(:,:),                       & ! (in)
                   0.0_RP, 1.0E3_RP, VAR_NAME(I_SFC_TEMP),           & ! (in)
                   __FILE__, __LINE__                                ) ! (in)
    call VALCHECK( IA, IS, IE, JA, JS, JE, &
                   ATMOS_PHY_SF_SFC_albedo(:,:,I_R_direct ,I_R_IR ), & ! (in)
                   0.0_RP, 1.0E0_RP, VAR_NAME(I_SFC_ALB_IR_dir),     & ! (in)
                   __FILE__, __LINE__                                ) ! (in)
    call VALCHECK( IA, IS, IE, JA, JS, JE, &
                   ATMOS_PHY_SF_SFC_albedo(:,:,I_R_diffuse,I_R_IR ), & ! (in)
                   0.0_RP, 1.0E0_RP, VAR_NAME(I_SFC_ALB_IR_dif),     & ! (in)
                   __FILE__, __LINE__                                ) ! (in)
    call VALCHECK( IA, IS, IE, JA, JS, JE, &
                   ATMOS_PHY_SF_SFC_albedo(:,:,I_R_direct ,I_R_NIR), & ! (in)
                   0.0_RP, 1.0E0_RP, VAR_NAME(I_SFC_ALB_NIR_dir),    & ! (in)
                   __FILE__, __LINE__                                ) ! (in)
    call VALCHECK( IA, IS, IE, JA, JS, JE, &
                   ATMOS_PHY_SF_SFC_albedo(:,:,I_R_diffuse,I_R_NIR), & ! (in)
                   0.0_RP, 1.0E0_RP, VAR_NAME(I_SFC_ALB_NIR_dif),    & ! (in)
                   __FILE__, __LINE__                                ) ! (in)
    call VALCHECK( IA, IS, IE, JA, JS, JE, &
                   ATMOS_PHY_SF_SFC_albedo(:,:,I_R_direct ,I_R_VIS), & ! (in)
                   0.0_RP, 1.0E0_RP, VAR_NAME(I_SFC_ALB_VIS_dir),    & ! (in)
                   __FILE__, __LINE__                                ) ! (in)
    call VALCHECK( IA, IS, IE, JA, JS, JE, &
                   ATMOS_PHY_SF_SFC_albedo(:,:,I_R_diffuse,I_R_VIS), & ! (in)
                   0.0_RP, 1.0E0_RP, VAR_NAME(I_SFC_ALB_VIS_dif),    & ! (in)
                   __FILE__, __LINE__                                ) ! (in)
    call VALCHECK( IA, IS, IE, JA, JS, JE, &
                   ATMOS_PHY_SF_SFC_Z0M(:,:),                        & ! (in)
                   0.0_RP, 1.0E2_RP, VAR_NAME(I_SFC_Z0M),            & ! (in)
                   __FILE__, __LINE__                                ) ! (in)
    call VALCHECK( IA, IS, IE, JA, JS, JE, &
                   ATMOS_PHY_SF_SFC_Z0H(:,:),                        & ! (in)
                   0.0_RP, 1.0E2_RP, VAR_NAME(I_SFC_Z0H),            & ! (in)
                   __FILE__, __LINE__                                ) ! (in)
    call VALCHECK( IA, IS, IE, JA, JS, JE, &
                   ATMOS_PHY_SF_SFC_Z0E(:,:),                        & ! (in)
                   0.0_RP, 1.0E2_RP, VAR_NAME(I_SFC_Z0E),            & ! (in)
                   __FILE__, __LINE__                                ) ! (in)

    call STATISTICS_total( IA, IS, IE, JA, JS, JE,                           & ! [IN]
                           ATMOS_PHY_SF_SFC_TEMP  (:,:),                     & ! [IN]
                           VAR_NAME(I_SFC_TEMP),                             & ! [IN]
                           ATMOS_GRID_CARTESC_REAL_AREA(:,:),                & ! [IN]
                           ATMOS_GRID_CARTESC_REAL_TOTAREA                   ) ! [IN]
    call STATISTICS_total( IA, IS, IE, JA, JS, JE,                           & ! [IN]
                           ATMOS_PHY_SF_SFC_albedo(:,:,I_R_direct ,I_R_IR ), & ! [IN]
                           VAR_NAME(I_SFC_ALB_IR_dir),                       & ! [IN]
                           ATMOS_GRID_CARTESC_REAL_AREA(:,:),                & ! [IN]
                           ATMOS_GRID_CARTESC_REAL_TOTAREA                   ) ! [IN]
    call STATISTICS_total( IA, IS, IE, JA, JS, JE,                           & ! [IN]
                           ATMOS_PHY_SF_SFC_albedo(:,:,I_R_diffuse,I_R_IR ), & ! [IN]
                           VAR_NAME(I_SFC_ALB_IR_dif),                       & ! [IN]
                           ATMOS_GRID_CARTESC_REAL_AREA(:,:),                & ! [IN]
                           ATMOS_GRID_CARTESC_REAL_TOTAREA                   ) ! [IN]
    call STATISTICS_total( IA, IS, IE, JA, JS, JE,                           & ! [IN]
                           ATMOS_PHY_SF_SFC_albedo(:,:,I_R_direct ,I_R_NIR), & ! [IN]
                           VAR_NAME(I_SFC_ALB_NIR_dir),                      & ! [IN]
                           ATMOS_GRID_CARTESC_REAL_AREA(:,:),                & ! [IN]
                           ATMOS_GRID_CARTESC_REAL_TOTAREA                   ) ! [IN]
    call STATISTICS_total( IA, IS, IE, JA, JS, JE,                           & ! [IN]
                           ATMOS_PHY_SF_SFC_albedo(:,:,I_R_diffuse,I_R_NIR), & ! [IN]
                           VAR_NAME(I_SFC_ALB_NIR_dif),                      & ! [IN]
                           ATMOS_GRID_CARTESC_REAL_AREA(:,:),                & ! [IN]
                           ATMOS_GRID_CARTESC_REAL_TOTAREA                   ) ! [IN]
    call STATISTICS_total( IA, IS, IE, JA, JS, JE,                           & ! [IN]
                           ATMOS_PHY_SF_SFC_albedo(:,:,I_R_direct ,I_R_VIS), & ! [IN]
                           VAR_NAME(I_SFC_ALB_VIS_dir),                      & ! [IN]
                           ATMOS_GRID_CARTESC_REAL_AREA(:,:),                & ! [IN]
                           ATMOS_GRID_CARTESC_REAL_TOTAREA                   ) ! [IN]
    call STATISTICS_total( IA, IS, IE, JA, JS, JE,                           & ! [IN]
                           ATMOS_PHY_SF_SFC_albedo(:,:,I_R_diffuse,I_R_VIS), & ! [IN]
                           VAR_NAME(I_SFC_ALB_VIS_dif),                      & ! [IN]
                           ATMOS_GRID_CARTESC_REAL_AREA(:,:),                & ! [IN]
                           ATMOS_GRID_CARTESC_REAL_TOTAREA                   ) ! [IN]
    call STATISTICS_total( IA, IS, IE, JA, JS, JE,                           & ! [IN]
                           ATMOS_PHY_SF_SFC_Z0M   (:,:),                     & ! [IN]
                           VAR_NAME(I_SFC_Z0M),                              & ! [IN]
                           ATMOS_GRID_CARTESC_REAL_AREA(:,:),                & ! [IN]
                           ATMOS_GRID_CARTESC_REAL_TOTAREA                   ) ! [IN]
    call STATISTICS_total( IA, IS, IE, JA, JS, JE,                           & ! [IN]
                           ATMOS_PHY_SF_SFC_Z0H   (:,:),                     & ! [IN]
                           VAR_NAME(I_SFC_Z0H),                              & ! [IN]
                           ATMOS_GRID_CARTESC_REAL_AREA(:,:),                & ! [IN]
                           ATMOS_GRID_CARTESC_REAL_TOTAREA                   ) ! [IN]
    call STATISTICS_total( IA, IS, IE, JA, JS, JE,                           & ! [IN]
                           ATMOS_PHY_SF_SFC_Z0E   (:,:),                     & ! [IN]
                           VAR_NAME(I_SFC_Z0E),                              & ! [IN]
                           ATMOS_GRID_CARTESC_REAL_AREA(:,:),                & ! [IN]
                           ATMOS_GRID_CARTESC_REAL_TOTAREA                   ) ! [IN]

    return
  end subroutine ATMOS_PHY_SF_vars_check

end module mod_atmos_phy_sf_vars
