!-------------------------------------------------------------------------------
!> module ATMOS admin
!!
!! @par Description
!!          Atmosphere submodel administrator
!!
!! @author Team SCALE
!<
!-------------------------------------------------------------------------------
module mod_atmos_admin
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_stdio
  use scale_prof
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: ATMOS_ADMIN_setup
  public :: ATMOS_ADMIN_getscheme

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  logical,                public :: ATMOS_do          = .true. ! main switch for the model

  character(len=H_SHORT), public :: ATMOS_DYN_TYPE    = 'NONE'
  character(len=H_SHORT), public :: ATMOS_DYN_TSTEP_TRACER_TYPE  = 'FVM-HEVE'
  character(len=H_SHORT), public :: ATMOS_DYN_TSTEP_LARGE_TYPE   = 'FVM-HEVE'
  character(len=H_SHORT), public :: ATMOS_DYN_TINTEG_SHORT_TYPE  = 'RK4'
  character(len=H_SHORT), public :: ATMOS_DYN_TINTEG_TRACER_TYPE = 'RK3WS2002'
  character(len=H_SHORT), public :: ATMOS_DYN_TINTEG_LARGE_TYPE  = 'EULER'
  character(len=H_SHORT), public :: ATMOS_PHY_MP_TYPE = 'NONE'
  character(len=H_SHORT), public :: ATMOS_PHY_AE_TYPE = 'NONE'
  character(len=H_SHORT), public :: ATMOS_PHY_CH_TYPE = 'NONE'
  character(len=H_SHORT), public :: ATMOS_PHY_RD_TYPE = 'NONE'
  character(len=H_SHORT), public :: ATMOS_PHY_SF_TYPE = 'NONE'
  character(len=H_SHORT), public :: ATMOS_PHY_TB_TYPE = 'NONE'
  character(len=H_SHORT), public :: ATMOS_PHY_CP_TYPE = 'NONE'

  logical,                public :: ATMOS_USE_AVERAGE = .false.

  logical,                public :: ATMOS_sw_dyn
  logical,                public :: ATMOS_sw_phy_mp
  logical,                public :: ATMOS_sw_phy_ae
  logical,                public :: ATMOS_sw_phy_ch
  logical,                public :: ATMOS_sw_phy_rd
  logical,                public :: ATMOS_sw_phy_sf
  logical,                public :: ATMOS_sw_phy_tb
  logical,                public :: ATMOS_sw_phy_cp

  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine ATMOS_ADMIN_setup
    use scale_process, only: &
       PRC_MPIstop
    implicit none

    NAMELIST / PARAM_ATMOS / &
       ATMOS_do,          &
       ATMOS_DYN_TYPE,    &
       ATMOS_DYN_TINTEG_SHORT_TYPE, &
       ATMOS_DYN_TINTEG_TRACER_TYPE, &
       ATMOS_DYN_TINTEG_LARGE_TYPE, &
       ATMOS_PHY_MP_TYPE, &
       ATMOS_PHY_AE_TYPE, &
       ATMOS_PHY_CH_TYPE, &
       ATMOS_PHY_RD_TYPE, &
       ATMOS_PHY_SF_TYPE, &
       ATMOS_PHY_TB_TYPE, &
       ATMOS_PHY_CP_TYPE, &
       ATMOS_USE_AVERAGE

    integer :: ierr
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '++++++ Module[ADMIN] / Categ[ATMOS] / Origin[SCALE-LES]'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_ATMOS,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_ATMOS. Check!'
       call PRC_MPIstop
    endif
    if( IO_LNML ) write(IO_FID_LOG,nml=PARAM_ATMOS)

    !-----< module component check >-----

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** Atmosphere model components ***'

    if ( ATMOS_do ) then
       if( IO_L ) write(IO_FID_LOG,*) '*** Atmosphere model : ON'
    else
       if( IO_L ) write(IO_FID_LOG,*) '*** Atmosphere model : OFF'
    endif

    if( IO_L ) write(IO_FID_LOG,*) '*** Dynamics...'

    if ( ATMOS_DYN_TYPE == 'OFF' ) then
       if( IO_L ) write(IO_FID_LOG,*) '*** +Dynamical core   : OFF'
       if( IO_L ) write(IO_FID_LOG,*) '*** +Advection        : OFF'
       ATMOS_sw_dyn = .false.
    else if ( ATMOS_DYN_TYPE == 'NONE' ) then
       ! The advection is disbled
       ! The tendencies calculated by physical processed are added
       if( IO_L ) write(IO_FID_LOG,*) '*** +Dynamical core   : ON, ', trim(ATMOS_DYN_TYPE)
       if( IO_L ) write(IO_FID_LOG,*) '*** +Advection        : OFF'
       ATMOS_sw_dyn = .true.
    else
       if( IO_L ) write(IO_FID_LOG,*) '*** +Dynamical core   : ON, ', trim(ATMOS_DYN_TYPE)
       if( IO_L ) write(IO_FID_LOG,*) '*** +Advection        : ON'
       ATMOS_sw_dyn = .true.
    endif

    if ( ATMOS_sw_dyn ) then
       if( IO_L ) write(IO_FID_LOG,*) '*** +Short time step  : ', trim(ATMOS_DYN_TINTEG_SHORT_TYPE)
       if( IO_L ) write(IO_FID_LOG,*) '*** +Tracer advection : ', trim(ATMOS_DYN_TINTEG_TRACER_TYPE)
       if( IO_L ) write(IO_FID_LOG,*) '*** +Large time step  : ', trim(ATMOS_DYN_TINTEG_LARGE_TYPE)
    end if

    if( IO_L ) write(IO_FID_LOG,*) '*** Physics...'

    if ( ATMOS_PHY_MP_TYPE /= 'OFF' .AND. ATMOS_PHY_MP_TYPE /= 'NONE' ) then
       if( IO_L ) write(IO_FID_LOG,*) '*** +Cloud Microphysics   : ON, ', trim(ATMOS_PHY_MP_TYPE)
       ATMOS_sw_phy_mp = .true.
    else
       if( IO_L ) write(IO_FID_LOG,*) '*** +Cloud Microphysics   : OFF'
       ATMOS_sw_phy_mp = .false.
    endif

    if ( ATMOS_PHY_AE_TYPE /= 'OFF' .AND. ATMOS_PHY_AE_TYPE /= 'NONE' ) then
       if( IO_L ) write(IO_FID_LOG,*) '*** +Aerosol Microphysics : ON, ', trim(ATMOS_PHY_AE_TYPE)
       ATMOS_sw_phy_ae = .true.
    else
       if( IO_L ) write(IO_FID_LOG,*) '*** +Aerosol Microphysics : OFF'
       ATMOS_sw_phy_ae = .false.
    endif

    if ( ATMOS_PHY_CH_TYPE /= 'OFF' .AND. ATMOS_PHY_CH_TYPE /= 'NONE' ) then
       if( IO_L ) write(IO_FID_LOG,*) '*** +Chemistry            : ON, ', trim(ATMOS_PHY_CH_TYPE)
       ATMOS_sw_phy_ch = .true.
    else
       if( IO_L ) write(IO_FID_LOG,*) '*** +Chemistry            : OFF'
       ATMOS_sw_phy_ch = .false.
    endif

    if ( ATMOS_PHY_RD_TYPE /= 'OFF' .AND. ATMOS_PHY_RD_TYPE /= 'NONE' ) then
       if( IO_L ) write(IO_FID_LOG,*) '*** +Radiative transfer   : ON, ', trim(ATMOS_PHY_RD_TYPE)
       ATMOS_sw_phy_rd = .true.
    else
       if( IO_L ) write(IO_FID_LOG,*) '*** +Radiative transfer   : OFF'
       ATMOS_sw_phy_rd = .false.
    endif

    if ( ATMOS_PHY_SF_TYPE /= 'OFF' .AND. ATMOS_PHY_SF_TYPE /= 'NONE' ) then
       if( IO_L ) write(IO_FID_LOG,*) '*** +Surface Flux         : ON, ', trim(ATMOS_PHY_SF_TYPE)
       ATMOS_sw_phy_sf = .true.
    else
       if( IO_L ) write(IO_FID_LOG,*) '*** +Surface Flux         : OFF'
       ATMOS_sw_phy_sf = .false.
    endif

    if ( ATMOS_PHY_TB_TYPE /= 'OFF' .AND. ATMOS_PHY_TB_TYPE /= 'NONE' ) then
       if( IO_L ) write(IO_FID_LOG,*) '*** +Sub-grid Turbulence  : ON, ', trim(ATMOS_PHY_TB_TYPE)
       ATMOS_sw_phy_tb = .true.
    else
       if( IO_L ) write(IO_FID_LOG,*) '*** +Sub-grid Turbulence  : OFF'
       ATMOS_sw_phy_tb = .false.
    endif

    if ( ATMOS_PHY_CP_TYPE /= 'OFF' .AND. ATMOS_PHY_CP_TYPE /= 'NONE' ) then
       if( IO_L ) write(IO_FID_LOG,*) '*** +Convection Param.    : ON, ', trim(ATMOS_PHY_CP_TYPE)
       ATMOS_sw_phy_cp = .true.
    else
       if( IO_L ) write(IO_FID_LOG,*) '*** +Convection Param.    : OFF'
       ATMOS_sw_phy_cp = .false.
    endif

    if( IO_L ) write(IO_FID_LOG,*)
    if ( ATMOS_USE_AVERAGE ) then
       if( IO_L ) write(IO_FID_LOG,*) '*** Atmos use average?    : YES'
    else
       if( IO_L ) write(IO_FID_LOG,*) '*** Atmos use average?    : NO'
    endif

    return
  end subroutine ATMOS_ADMIN_setup

  !-----------------------------------------------------------------------------
  !> Get name of scheme for each component
  subroutine ATMOS_ADMIN_getscheme( &
       component_name, &
       scheme_name     )
    use scale_process, only: &
       PRC_MPIstop
    implicit none

    character(len=*),       intent(in)  :: component_name
    character(len=H_SHORT), intent(out) :: scheme_name
    !---------------------------------------------------------------------------

    select case (component_name)
    case ("DYN")
       scheme_name = ATMOS_DYN_TYPE
    case ("PHY_MP")
       scheme_name = ATMOS_PHY_MP_TYPE
    case ("PHY_AE")
       scheme_name = ATMOS_PHY_AE_TYPE
    case ("PHY_CH")
       scheme_name = ATMOS_PHY_CH_TYPE
    case ("PHY_RD")
       scheme_name = ATMOS_PHY_RD_TYPE
    case ("PHY_SF")
       scheme_name = ATMOS_PHY_SF_TYPE
    case ("PHY_TB")
       scheme_name = ATMOS_PHY_TB_TYPE
    case ("PHY_CP")
       scheme_name = ATMOS_PHY_CP_TYPE
    case default
       write(*,*) 'xxx Unsupported component_name. Check!', trim(component_name)
       call PRC_MPIstop
    end select

    return
  end subroutine ATMOS_ADMIN_getscheme

end module mod_atmos_admin
