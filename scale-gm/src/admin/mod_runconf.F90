!-------------------------------------------------------------------------------
!> Module run configuration
!!
!! @par Description
!!          This module is for managing run configuration
!!
!! @author NICAM developers, Team SCALE
!<
!-------------------------------------------------------------------------------
module mod_runconf
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use scale_precision
  use scale_io
  use scale_tracer

  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: runconf_setup

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  character(len=H_SHORT), public :: RUNNAME = ''

  !---< Component Selector >---

  !--- Dynamics
  integer,                public :: NON_HYDRO_ALPHA    = 1 ! Nonhydrostatic/hydrostatic flag
  integer,                public :: DYN_DIV_NUM        = 1
  character(len=H_SHORT), public :: TRC_ADV_TYPE       = 'MIURA2004'
  character(len=H_SHORT), public :: NDIFF_LOCATION     = 'IN_LARGE_STEP2'
  logical,                public :: FLAG_NUDGING       = .false.
  logical,                public :: THUBURN_LIM        = .true.  ! [add] 20130613 R.Yoshida

  !--- Physics
  character(len=H_SHORT), public :: ATMOS_PHY_TYPE               = 'NONE'
  character(len=H_SHORT), public :: ROUGHNESS_SEA_TYPE = 'DEFAULT'
  character(len=H_SHORT), public :: AF_TYPE            = 'NONE'

  character(len=H_SHORT), public :: OUT_FILE_TYPE      = 'DEFAULT'

  !---< tracer ID setting >---

  integer, public            :: PRG_vmax        ! total number of prognostic variables
  integer, public, parameter :: PRG_vmax0  = 6

  integer, public, parameter :: I_RHOG     =  1 ! Density x G^1/2
  integer, public, parameter :: I_RHOGVX   =  2 ! Density x G^1/2 x Horizontal velocity (X-direction)
  integer, public, parameter :: I_RHOGVY   =  3 ! Density x G^1/2 x Horizontal velocity (Y-direction)
  integer, public, parameter :: I_RHOGVZ   =  4 ! Density x G^1/2 x Horizontal velocity (Z-direction)
  integer, public, parameter :: I_RHOGW    =  5 ! Density x G^1/2 x Vertical   velocity
  integer, public, parameter :: I_RHOGE    =  6 ! Density x G^1/2 x Energy
  integer, public, parameter :: I_RHOGQstr =  7 ! tracers
  integer, public            :: I_RHOGQend = -1 !

  character(len=H_SHORT), public  :: PRG_name(PRG_vmax0)
  data PRG_name / 'rhog', 'rhogvx', 'rhogvy', 'rhogvz', 'rhogw', 'rhoge' /

  integer, public            :: DIAG_vmax       ! total number of diagnostic variables
  integer, public, parameter :: DIAG_vmax0 = 6

  integer, public, parameter :: I_pre      =  1 ! Pressure
  integer, public, parameter :: I_tem      =  2 ! Temperature
  integer, public, parameter :: I_vx       =  3 ! Horizontal velocity (X-direction)
  integer, public, parameter :: I_vy       =  4 ! Horizontal velocity (Y-direction)
  integer, public, parameter :: I_vz       =  5 ! Horizontal velocity (Z-direction)
  integer, public, parameter :: I_w        =  6 ! Vertical   velocity
  integer, public, parameter :: I_qstr     =  7 ! tracers
  integer, public            :: I_qend     = -1 !

  character(len=H_SHORT), public  :: DIAG_name(DIAG_vmax0)
  data DIAG_name / 'pre', 'tem', 'vx', 'vy', 'vz', 'w' /

  !--- No. of band for rad.
  integer, public, parameter :: NRBND     = 3
  integer, public, parameter :: NRBND_VIS = 1
  integer, public, parameter :: NRBND_NIR = 2
  integer, public, parameter :: NRBND_IR  = 3

  !--- direct/diffuse
  integer, public, parameter :: NRDIR         = 2
  integer, public, parameter :: NRDIR_DIRECT  = 1
  integer, public, parameter :: NRDIR_DIFFUSE = 2

  !--- roughness  parameter
  integer, public, parameter :: NTYPE_Z0 = 3
  integer, public, parameter :: N_Z0M    = 1
  integer, public, parameter :: N_Z0H    = 2
  integer, public, parameter :: N_Z0E    = 3

  !-----------------------------------------------------------------------------
  !
  !++ Private procedures
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine RUNCONF_setup
    use scale_prc, only: &
       PRC_abort
    use mod_ocean_admin, only: &
       OCEAN_TYPE
    use mod_land_admin, only: &
       LAND_TYPE
    use mod_cpl_admin, only: &
       CPL_sw
    implicit none

    namelist /RUNCONFPARAM/ &
       RUNNAME,            &
       NON_HYDRO_ALPHA,    &
       DYN_DIV_NUM,        &
       TRC_ADV_TYPE,       &
       NDIFF_LOCATION,     &
       FLAG_NUDGING,       &
       THUBURN_LIM,        & ! R.Yoshida 13/06/13 [add]
       ATMOS_PHY_TYPE,     & 
       ROUGHNESS_SEA_TYPE, &
       LAND_TYPE,          &
       OCEAN_TYPE,         &
       AF_TYPE,            &
       OUT_FILE_TYPE

    integer :: ierr
    !---------------------------------------------------------------------------

    !--- read parameters
    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[runconf]/Category[nhm share]'
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=RUNCONFPARAM,iostat=ierr)
    if ( ierr < 0 ) then
       if( IO_L ) write(IO_FID_LOG,*) '*** RUNCONFPARAM is not specified. use default.'
    elseif( ierr > 0 ) then
       write(*,*) 'xxx Not appropriate names in namelist RUNCONFPARAM. STOP.'
       call PRC_abort
    endif
    if( IO_NML ) write(IO_FID_NML,nml=RUNCONFPARAM)


    call RUNCONF_component_setup

    call RUNCONF_tracer_setup

    CPL_sw = .false.

    return
  end subroutine RUNCONF_setup

  !-----------------------------------------------------------------------------
  !> component check
  subroutine RUNCONF_component_setup
    implicit none
    !---------------------------------------------------------------------------

    if( THUBURN_LIM ) then ![add] 20130613 R.Yoshida
       if( IO_L ) write(IO_FID_LOG,*) 'Run with \"Thuburn Limiter\" in MIURA2004 Advection'
    else
       if( IO_L ) write(IO_FID_LOG,*) '### Without \"Thuburn Limiter\" in MIURA2004 Advection'
    endif

    return
  end subroutine RUNCONF_component_setup

  !-----------------------------------------------------------------------------
  !> tracer setup
  subroutine RUNCONF_tracer_setup
    use scale_prc, only: &
       PRC_abort
    use mod_chemvar, only: &
       CHEMVAR_setup
    use mod_atmos_admin, only: &
       ATMOS_PHY_CH_TYPE
    use scale_atmos_hydrometeor, only: &
       HYDROMETEOR_regist => ATMOS_HYDROMETEOR_regist
    implicit none

    integer :: I_QV
    integer :: v, i
    !---------------------------------------------------------------------------


    ! tentative
    if ( QA == 0 ) then
       call HYDROMETEOR_regist( 0, 0,                           & ! [IN]
                                (/ "QV" /),                     & ! [IN]
                                (/ "water vapor mass ratio" /), & ! [IN]
                                (/ "kg/kg" /),                  & ! [IN]
                                I_QV                            ) ! [OUT]
    end if

    !--- Tracer for chemistry
    if ( ATMOS_PHY_CH_TYPE == 'PASSIVE' )then
       call CHEMVAR_setup
    endif

    PRG_vmax   = PRG_vmax0  + QA
    I_RHOGQend = PRG_vmax

    DIAG_vmax  = DIAG_vmax0 + QA
    I_qend     = DIAG_vmax

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** Prognostic Tracers'
    if( IO_L ) write(IO_FID_LOG,*) '|=========================================================|'
    if( IO_L ) write(IO_FID_LOG,*) '|       :varname         :description                     |'
    do v = 1, QA
       if( IO_L ) write(IO_FID_LOG,'(1x,A,I4,A,A16,A,A,A)') '|ID=', v, ':', TRACER_name(v), ':', TRACER_desc(v),'|'
    enddo

    return
  end subroutine RUNCONF_tracer_setup

end module mod_runconf
