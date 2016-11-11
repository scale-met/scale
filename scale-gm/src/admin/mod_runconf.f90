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
  use scale_stdio
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
  character(len=H_SHORT), public :: RAIN_TYPE          = 'DRY'
  logical,                public :: opt_2moment_water  = .false.

  character(len=H_SHORT), public :: CP_TYPE            = 'NONE'
  character(len=H_SHORT), public :: MP_TYPE            = 'NONE'
  character(len=H_SHORT), public :: RD_TYPE            = 'NONE'
  character(len=H_SHORT), public :: SF_TYPE            = 'DEFAULT'
  character(len=H_SHORT), public :: ROUGHNESS_SEA_TYPE = 'DEFAULT'
  character(len=H_SHORT), public :: OCEAN_TYPE         = 'NONE'
  character(len=H_SHORT), public :: RIV_TYPE           = 'NONE'
  character(len=H_SHORT), public :: LAND_TYPE          = 'NONE'
  character(len=H_SHORT), public :: TB_TYPE            = 'NONE'
  character(len=H_SHORT), public :: AE_TYPE            = 'NONE'
  character(len=H_SHORT), public :: CHEM_TYPE          = 'NONE'
  character(len=H_SHORT), public :: GWD_TYPE           = 'NONE'
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

  integer, public            :: TRC_vmax   =  0 ! total number of tracers

  character(len=H_SHORT), public, allocatable :: TRC_name(:) ! short name  of tracer
  character(len=H_MID),   public, allocatable :: WLABEL  (:) ! description of tracer

  integer, public            :: NQW_MAX    =  0 ! subtotal number of water mass tracers
  integer, public            :: NQW_STR    = -1 ! start index of water mass tracers
  integer, public            :: NQW_END    = -1 ! end   index of water mass tracers
  integer, public            :: I_QV       = -1 ! Water vapor
  integer, public            :: I_QC       = -1 ! Cloud water
  integer, public            :: I_QR       = -1 ! Rain
  integer, public            :: I_QI       = -1 ! Ice
  integer, public            :: I_QS       = -1 ! Snow
  integer, public            :: I_QG       = -1 ! Graupel

  integer, public            :: NNW_MAX    =  0 ! subtotal number of water number tracers
  integer, public            :: NNW_STR    = -1 ! start index of water number tracers
  integer, public            :: NNW_END    = -1 ! end   index of water number tracers
  integer, public            :: I_NC       = -1 ! Cloud water (number)
  integer, public            :: I_NR       = -1 ! Rain        (number)
  integer, public            :: I_NI       = -1 ! Ice         (number)
  integer, public            :: I_NS       = -1 ! Snow        (number)
  integer, public            :: I_NG       = -1 ! Graupel     (number)

  integer, public            :: NTB_MAX    =  0 ! subtotal number of turbulent tracers
  integer, public            :: I_TKE      = -1 ! turbulence kinetic energy
  integer, public            :: I_QKEp     = -1
  integer, public            :: I_TSQp     = -1
  integer, public            :: I_QSQp     = -1
  integer, public            :: I_COVp     = -1

  integer, public            :: NCHEM_MAX  =  0 ! subtotal number of chemical (or general purpose) tracers
  integer, public            :: NCHEM_STR  = -1 ! start index of chemical (or general purpose) tracers
  integer, public            :: NCHEM_END  = -1 ! end   index of chemical (or general purpose) tracers

  !--- specific heat of water on const pressure
  real(RP), public, allocatable :: CVW(:)
  real(RP), public, allocatable :: CPW(:)

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
    use scale_process, only: &
       PRC_MPIstop
    implicit none

    namelist /RUNCONFPARAM/ &
       RUNNAME,            &
       NON_HYDRO_ALPHA,    &
       DYN_DIV_NUM,        &
       TRC_ADV_TYPE,       &
       NDIFF_LOCATION,     &
       FLAG_NUDGING,       &
       THUBURN_LIM,        & ! R.Yoshida 13/06/13 [add]
       RAIN_TYPE,          &
       CP_TYPE,            &
       MP_TYPE,            &
       RD_TYPE,            &
       SF_TYPE,            &
       ROUGHNESS_SEA_TYPE, &
       LAND_TYPE,          &
       OCEAN_TYPE,         &
       TB_TYPE,            &
       AE_TYPE,            &
       CHEM_TYPE,          &
       GWD_TYPE,           &
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
       call PRC_MPIstop
    endif
    if( IO_LNML ) write(IO_FID_LOG,nml=RUNCONFPARAM)

    call RUNCONF_component_setup

    call RUNCONF_tracer_setup

    call RUNCONF_thermodyn_setup

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
    use scale_process, only: &
       PRC_MPIstop
    use mod_chemvar, only: &
       CHEMVAR_setup, &
       CHEM_TRC_vmax, &
       CHEM_TRC_name, &
       CHEM_TRC_desc
    implicit none

    integer :: v, i
    !---------------------------------------------------------------------------

    !--- counting tracer
    TRC_vmax = 0

    !--- Mass tracer for Water
    if    ( RAIN_TYPE == 'DRY' ) then
       NQW_MAX = 1
       I_QV    = TRC_vmax + 1
    elseif( RAIN_TYPE == 'CLOUD_PARAM' ) then
       NQW_MAX = 2
       I_QV    = TRC_vmax + 1
       I_QC    = TRC_vmax + 2
    elseif( RAIN_TYPE == 'WARM' ) then
       NQW_MAX = 3
       I_QV    = TRC_vmax + 1
       I_QC    = TRC_vmax + 2
       I_QR    = TRC_vmax + 3
    elseif( RAIN_TYPE == 'COLD' ) then
       NQW_MAX = 6
       I_QV    = TRC_vmax + 1
       I_QC    = TRC_vmax + 2
       I_QR    = TRC_vmax + 3
       I_QI    = TRC_vmax + 4
       I_QS    = TRC_vmax + 5
       I_QG    = TRC_vmax + 6
    else
       write(*,*) 'xxx RAIN_TYPE must be set to DRY,CLOUD_PARAM,WARM or COLD. STOP.'
       call PRC_MPIstop
    endif
    NQW_STR  = TRC_vmax + 1
    NQW_END  = TRC_vmax + NQW_MAX
    TRC_vmax = TRC_vmax + NQW_MAX

    !--- Number tracer for Water
    if ( opt_2moment_water ) then
       if    ( RAIN_TYPE == 'DRY' ) then
          NNW_MAX = 0
       elseif( RAIN_TYPE == 'CLOUD_PARAM' ) then
          NNW_MAX = 1
          I_NC    = TRC_vmax + 1
       elseif( RAIN_TYPE == 'WARM' ) then
          NNW_MAX = 2
          I_NC    = TRC_vmax + 1
          I_NR    = TRC_vmax + 2
       elseif( RAIN_TYPE == 'COLD' ) then
          NNW_MAX = 5
          I_NC    = TRC_vmax + 1
          I_NR    = TRC_vmax + 2
          I_NI    = TRC_vmax + 3
          I_NS    = TRC_vmax + 4
          I_NG    = TRC_vmax + 5
       endif
       NNW_STR  = TRC_vmax + min(1,NNW_MAX)
       NNW_END  = TRC_vmax + NNW_MAX
       TRC_vmax = TRC_vmax + NNW_MAX
    endif

    !--- Tracer for turbulence
    if    ( TB_TYPE == 'MY2.5' ) then
       NTB_MAX = 1
       I_TKE   = TRC_vmax + 1
    elseif( TB_TYPE == 'MYNN2.5' ) then
       NTB_MAX = 1
       I_QKEp  = TRC_vmax + 1
    elseif( TB_TYPE == 'MYNN3' ) then
       NTB_MAX = 4
       I_QKEp  = TRC_vmax + 1
       I_TSQp  = TRC_vmax + 2
       I_QSQp  = TRC_vmax + 3
       I_COVp  = TRC_vmax + 4
    endif
    TRC_vmax = TRC_vmax + NTB_MAX

    !--- Tracer for chemistry
    call CHEMVAR_setup

    if ( CHEM_TYPE == 'PASSIVE' )then
       NCHEM_MAX = CHEM_TRC_vmax
       NCHEM_STR = TRC_vmax + min(1,NCHEM_MAX)
       NCHEM_END = TRC_vmax + NCHEM_MAX
    endif
    TRC_vmax = TRC_vmax + NCHEM_MAX

    allocate( TRC_name(TRC_vmax) ) ! [add] H.Yashiro 20110819
    allocate( WLABEL  (TRC_vmax) ) ! 08/04/12 [Add] T.Mitsui
    TRC_name(:) = ""
    WLABEL  (:) = ""

    !--- Labeling
    do v = 1, TRC_vmax
       if    ( v == I_QV ) then
          TRC_name(v) = 'qv'
          WLABEL  (v) = 'VAPOR'
       elseif( v == I_QC ) then
          TRC_name(v) = 'qc'
          WLABEL  (v) = 'CLOUD'
       elseif( v == I_QR ) then
          TRC_name(v) = 'qr'
          WLABEL  (v) = 'RAIN'
       elseif( v == I_QI ) then
          TRC_name(v) = 'qi'
          WLABEL  (v) = 'ICE'
       elseif( v == I_QS ) then
          TRC_name(v) = 'qs'
          WLABEL  (v) = 'SNOW'
       elseif( v == I_QG ) then
          TRC_name(v) = 'qg'
          WLABEL  (v) = 'GRAUPEL'

       elseif( v == I_NC )then
          TRC_name(v) = 'nc'
          WLABEL  (v) = 'CLOUD_NUM'
       elseif( v == I_NR )then
          TRC_name(v) = 'nr'
          WLABEL  (v) = 'RAIN_NUM'
       elseif( v == I_NI )then
          TRC_name(v) = 'ni'
          WLABEL  (v) = 'ICE_NUM'
       elseif( v == I_NS )then
          TRC_name(v) = 'ns'
          WLABEL  (v) = 'SNOW_NUM'
       elseif( v == I_NG )then
          TRC_name(v) = 'ng'
          WLABEL  (v) = 'GRAUPEL_NUM'
       elseif( v == NCHEM_STR ) then
          do i = 1, NCHEM_MAX
             TRC_name(v+i-1) = CHEM_TRC_name(i)
             WLABEL  (v+i-1) = CHEM_TRC_desc(i)
          enddo
       endif
    enddo

    PRG_vmax   = PRG_vmax0  + TRC_vmax
    I_RHOGQend = PRG_vmax

    DIAG_vmax  = DIAG_vmax0 + TRC_vmax
    I_qend     = DIAG_vmax

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** Prognostic Tracers'
    if( IO_L ) write(IO_FID_LOG,*) '|=========================================================|'
    if( IO_L ) write(IO_FID_LOG,*) '|       :varname         :description                     |'
    do v = 1, TRC_vmax
       if( IO_L ) write(IO_FID_LOG,'(1x,A,I4,A,A16,A,A,A)') '|ID=', v, ':', TRC_name(v), ':', WLABEL(v),'|'
    enddo
    if( IO_L ) write(IO_FID_LOG,*) '|=========================================================|'
    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** thermodynamic(water) tracers'
    if( IO_L ) write(IO_FID_LOG,*) '-->', NQW_MAX, ' tracers(',NQW_STR,'-',NQW_END,')'

    return
  end subroutine RUNCONF_tracer_setup

  !-----------------------------------------------------------------------------
  !> thermodynamic setup
  subroutine RUNCONF_thermodyn_setup
    use scale_const, only: &
       CONST_THERMODYN_TYPE, &
       CPdry => CONST_CPdry, &
       CPvap => CONST_CPvap, &
       CVdry => CONST_CVdry, &
       CVvap => CONST_CVvap, &
       CL    => CONST_CL,    &
       CI    => CONST_CI
    implicit none

    integer :: v
    !---------------------------------------------------------------------------
    ! 'SIMPLE': standard approximation CVD * T
    ! 'EXACT': exact formulation
    !         -> if warm rain
    !            qd*CVD*T + qv*CVV*T + (qc+qr)*CPL*T
    !         -> if cold rain
    !            qd*CVD*T + qv*CVV*T + (qc+qr)*CPL*T
    !            + (qi+qs)*CPI*T

    !--- Heat capacity for thermodynamics
    allocate( CVW(NQW_STR:NQW_END) )
    allocate( CPW(NQW_STR:NQW_END) )

    if ( CONST_THERMODYN_TYPE == 'SIMPLE' ) then
       do v = NQW_STR, NQW_END
          if    ( v == I_QV ) then ! vapor
             CVW(v) = CVdry
             CPW(v) = CPdry
          elseif( v == I_QC ) then ! cloud
             CVW(v) = CVdry
             CPW(v) = CVdry
          elseif( v == I_QR ) then ! rain
             CVW(v) = CVdry
             CPW(v) = CVdry
          elseif( v == I_QI ) then ! ice
             CVW(v) = CVdry
             CPW(v) = CVdry
          elseif( v == I_QS ) then ! snow
             CVW(v) = CVdry
             CPW(v) = CVdry
          elseif( v == I_QG ) then ! graupel
             CVW(v) = CVdry
             CPW(v) = CVdry
          endif
       enddo
    elseif( CONST_THERMODYN_TYPE == 'SIMPLE2' ) then
       do v = NQW_STR, NQW_END
          if    ( v == I_QV ) then ! vapor
             CVW(v) = CVvap
             CPW(v) = CPvap
          elseif( v == I_QC ) then ! cloud
             CVW(v) = CPvap
             CPW(v) = CPvap
          elseif( v == I_QR ) then ! rain
             CVW(v) = CPvap
             CPW(v) = CPvap
          elseif( v == I_QI ) then ! ice
             CVW(v) = CPvap
             CPW(v) = CPvap
          elseif( v == I_QS ) then ! snow
             CVW(v) = CPvap
             CPW(v) = CPvap
          elseif( v == I_QG ) then ! graupel
             CVW(v) = CPvap
             CPW(v) = CPvap
          endif
       enddo
    elseif( CONST_THERMODYN_TYPE == 'EXACT' ) then
       do v = NQW_STR, NQW_END
          if    ( v == I_QV ) then ! vapor
             CVW(v) = CVvap
             CPW(v) = CPvap
          elseif( v == I_QC ) then ! cloud
             CVW(v) = CL
             CPW(v) = CL
          elseif( v == I_QR ) then ! rain
             CVW(v) = CL
             CPW(v) = CL
          elseif( v == I_QI ) then ! ice
             CVW(v) = CI
             CPW(v) = CI
          elseif( v == I_QS ) then ! snow
             CVW(v) = CI
             CPW(v) = CI
          elseif( v == I_QG ) then ! graupel
             CVW(v) = CI
             CPW(v) = CI
          endif
       enddo
    endif

    return
  end subroutine RUNCONF_thermodyn_setup

end module mod_runconf
