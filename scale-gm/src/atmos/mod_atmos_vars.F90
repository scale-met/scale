!-------------------------------------------------------------------------------
!> module ATMOSPHERIC Variables
!!
!! @par Description
!!          Container for atmospheric variables
!!
!! @author Team SCALE
!!
!<
!-------------------------------------------------------------------------------
#include "scalelib.h"
module mod_atmos_vars
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_io
  use scale_prof
  use scale_debug
  use scale_atmos_grid_icoA_index
  use scale_tracer

  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: ATMOS_vars_setup
  public :: ATMOS_vars_fillhalo
  public :: ATMOS_vars_restart_read
  public :: ATMOS_vars_restart_write
  public :: ATMOS_vars_restart_check
  public :: ATMOS_vars_history_setpres
  public :: ATMOS_vars_history
  public :: ATMOS_vars_total
  public :: ATMOS_vars_calc_prognostics
  public :: ATMOS_vars_calc_diagnostics_fromIcoGrid
  public :: ATMOS_vars_calc_diagnostics
  public :: ATMOS_vars_get_diagnostic
  public :: ATMOS_vars_monitor

  public :: ATMOS_vars_restart_create
  public :: ATMOS_vars_restart_open
  public :: ATMOS_vars_restart_def_var
  public :: ATMOS_vars_restart_enddef
  public :: ATMOS_vars_restart_close


  interface ATMOS_vars_get_diagnostic
     procedure ATMOS_vars_get_diagnostic_3D
     procedure ATMOS_vars_get_diagnostic_2D
     procedure ATMOS_vars_get_diagnostic_1D
  end interface ATMOS_vars_get_diagnostic

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  logical,                public :: ATMOS_RESTART_OUTPUT                = .false.         !< Output restart file?

  character(len=H_LONG),  public :: ATMOS_RESTART_IN_BASENAME           = ''              !< Basename of the input  file
  logical,                public :: ATMOS_RESTART_IN_AGGREGATE                            !< Switch to use aggregate file
  logical,                public :: ATMOS_RESTART_IN_POSTFIX_TIMELABEL  = .false.         !< Add timelabel to the basename of input  file?
  character(len=H_LONG),  public :: ATMOS_RESTART_OUT_BASENAME          = ''              !< Basename of the output file
  logical,                public :: ATMOS_RESTART_OUT_AGGREGATE                           !< Switch to use aggregate file
  logical,                public :: ATMOS_RESTART_OUT_POSTFIX_TIMELABEL = .true.          !< Add timelabel to the basename of output file?
  character(len=H_MID),   public :: ATMOS_RESTART_OUT_TITLE             = 'ATMOS restart' !< Title    of the output file
  character(len=H_SHORT), public :: ATMOS_RESTART_OUT_DTYPE             = 'DEFAULT'       !< REAL4 or REAL8

  logical,                public :: ATMOS_RESTART_CHECK                 = .false.         !< Check value consistency?
  character(len=H_LONG),  public :: ATMOS_RESTART_CHECK_BASENAME        = 'restart_check'
  real(RP),               public :: ATMOS_RESTART_CHECK_CRITERION       = 1.E-6_RP

  ! prognostic variables
  real(RP), public, allocatable         :: DENS(:,:,:,:)
  real(RP), public, allocatable         :: RHOU(:,:,:,:)
  real(RP), public, allocatable         :: RHOV(:,:,:,:)
  real(RP), public, allocatable         :: MOMZ(:,:,:,:)
  real(RP), public, allocatable         :: RHOE(:,:,:,:)
  real(RP), public, allocatable         :: RHOQ(:,:,:,:,:)
  real(RP), public, allocatable, target :: QTRC(:,:,:,:,:)

  integer,  public                      :: I_QV
  real(RP), public, pointer             :: QV(:,:,:,:)
  real(RP), public, pointer             :: QC(:,:,:,:)
  real(RP), public, pointer             :: QR(:,:,:,:)
  real(RP), public, pointer             :: QI(:,:,:,:)
  real(RP), public, pointer             :: QS(:,:,:,:)
  real(RP), public, pointer             :: QG(:,:,:,:)
  real(RP), public, pointer             :: QH(:,:,:,:)

  ! public diagnostic variables
  real(RP), public, allocatable, target :: W    (:,:,:,:) !> velocity w [m/s]
  real(RP), public, allocatable, target :: U    (:,:,:,:) !> velocity u [m/s]
  real(RP), public, allocatable, target :: V    (:,:,:,:) !> velocity v [m/s]

  real(RP), public, allocatable, target :: POTT (:,:,:,:) !> potential temperature [K]
  real(RP), public, allocatable, target :: TEMP (:,:,:,:) !> temperature           [K]
  real(RP), public, allocatable, target :: PRES (:,:,:,:) !> pressure              [Pa=J/m3]
  real(RP), public, allocatable, target :: EXNER(:,:,:,:) !> Exner function (t/pt) [1]
  real(RP), public, allocatable, target :: PHYD (:,:,:,:) !> hydrostatic pressure  [Pa=J/m3]
  real(RP), public, allocatable, target :: PHYDH(:,:,:,:) !> hydrostatic pressure  [Pa=J/m3], layer interface

  real(RP), public, allocatable, target :: Qdry (:,:,:,:) !> dry air                [1]
  real(RP), public, allocatable, target :: Rtot (:,:,:,:) !> specific gass constant [J/kg/K]
  real(RP), public, allocatable, target :: CVtot(:,:,:,:) !> specific heat          [J/kg/K]
  real(RP), public, allocatable, target :: CPtot(:,:,:,:) !> specific heat          [J/kg/K]


  ! tentative
  real(RP), public, allocatable :: CZ(:,:,:,:)
  real(RP), public, allocatable :: FZ(:,:,:,:)

  real(RP), public, allocatable :: LON(:,:,:)
  real(RP), public, allocatable :: LAT(:,:,:)

  real(RP), public, allocatable :: Z1       (:,:,:)
  real(RP), public, allocatable :: TOPO_Zsfc(:,:,:)

  real(RP), public, allocatable :: fact_ocean(:,:,:)
  real(RP), public, allocatable :: fact_land (:,:,:)
  real(RP), public, allocatable :: fact_urban(:,:,:)

  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  logical,  private :: ATMOS_VARS_CHECKRANGE    = .false.
  real(RP), private :: ATMOS_VARS_CHECKCFL_SOFT = 1.0_RP  !> if Courant number exceeds this value, put the warning message
  real(RP), private :: ATMOS_VARS_CHECKCFL_HARD = 2.0_RP  !> if Courant number exceeds this value, put the error message and stop

  type Vinfo
     character(len=H_SHORT) :: NAME
     character(len=H_MID)   :: DESC
     character(len=H_SHORT) :: UNIT
     integer                :: ndims
     character(len=H_SHORT) :: dim_type
     character(len=H_MID)   :: STDNAME
  end type Vinfo

  ! private diagnostic variables
  real(RP), allocatable, target :: LHV  (:,:,:,:) !> latent heat for vaporization [J/kg]
  real(RP), allocatable, target :: LHS  (:,:,:,:) !> latent heat for sublimation  [J/kg]
  real(RP), allocatable, target :: LHF  (:,:,:,:) !> latent heat for fusion       [J/kg]

  real(RP), allocatable, target :: POTV (:,:,:,:) !> virtual pot. temp.      [K]
  real(RP), allocatable, target :: TEML (:,:,:,:) !> liquid water temp.      [K]
  real(RP), allocatable, target :: POTL (:,:,:,:) !> liquid water pot. temp. [K]
  real(RP), allocatable, target :: POTE (:,:,:,:) !> equivalent pot. temp.   [K]

  real(RP), allocatable, target :: QTOT (:,:,:,:) !> total water content  [1]
  real(RP), allocatable, target :: QHYD (:,:,:,:) !> total hydrometeornt  [1]
  real(RP), allocatable, target :: QLIQ (:,:,:,:) !> liquid water content [1]
  real(RP), allocatable, target :: QICE (:,:,:,:) !> ice water content    [1]

  real(RP), allocatable, target :: LWP  (:,:,:)   !> liquid water path  [g/m2]
  real(RP), allocatable, target :: IWP  (:,:,:)   !> ice    water path  [g/m2]
  real(RP), allocatable, target :: PW   (:,:,:)   !> precipitable water [g/m2]

  real(RP), allocatable, target :: PREC (:,:,:)   !> surface precipitation rate CP+MP(rain+snow) [kg/m2/s]
  real(RP), allocatable, target :: RAIN (:,:,:)   !> surface rain rate CP+MP [kg/m2/s]
  real(RP), allocatable, target :: SNOW (:,:,:)   !> surface snow rate CP+MP [kg/m2/s]

  real(RP), allocatable, target :: QSAT (:,:,:,:) !> saturation specific humidity        [1]
  real(RP), allocatable, target :: RHA  (:,:,:,:) !> relative humidity (liquid+ice)      [%]
  real(RP), allocatable, target :: RHL  (:,:,:,:) !> relative humidity against to liquid [%]
  real(RP), allocatable, target :: RHI  (:,:,:,:) !> relative humidity against to ice    [%]

  real(RP), allocatable, target :: VOR  (:,:,:,:) !> vertical vorticity    [1/s]
  real(RP), allocatable, target :: DIV  (:,:,:,:) !> divergence            [1/s]
  real(RP), allocatable, target :: HDIV (:,:,:,:) !> horizontal divergence [1/s]
  real(RP), allocatable, target :: Uabs (:,:,:,:) !> absolute velocity     [m/s]

  real(RP), allocatable, target :: N2   (:,:,:,:) !> squared Brunt-Vaisala frequency [/s2]
  real(RP), allocatable, target :: PBLH (:,:,:)   !> PBL height [m]

  real(RP), allocatable, target :: MSE  (:,:,:,:) !> moist static energy [m2/s2]
  real(RP), allocatable, target :: TDEW (:,:,:,:) !> dew point [K]

  real(RP), allocatable, target :: CAPE (:,:,:)   !> CAPE       [m2/s2]
  real(RP), allocatable, target :: CIN  (:,:,:)   !> CIN        [m2/s2]
  real(RP), allocatable, target :: LCL  (:,:,:)   !> LCL height [m]
  real(RP), allocatable, target :: LFC  (:,:,:)   !> LFC height [m]
  real(RP), allocatable, target :: LNB  (:,:,:)   !> LNB height [m]

  real(RP), allocatable, target :: ENGT (:,:,:,:) !> total     energy [J/m3]
  real(RP), allocatable, target :: ENGP (:,:,:,:) !> potential energy [J/m3]
  real(RP), allocatable, target :: ENGK (:,:,:,:) !> kinetic   energy [J/m3]
  real(RP), allocatable, target :: ENGI (:,:,:,:) !> internal  energy [J/m3]

  real(RP), allocatable, target :: DENS_MEAN(:)     !> horiz. mean of density         [kg/m3]
  real(RP), allocatable, target :: W_MEAN   (:)     !> horiz. mean of w               [m/s]
  real(RP), allocatable, target :: U_MEAN   (:)     !> horiz. mean of u               [m/s]
  real(RP), allocatable, target :: V_MEAN   (:)     !> horiz. mean of v               [m/s]
  real(RP), allocatable, target :: PT_MEAN  (:)     !> horiz. mean of pot.            [K]
  real(RP), allocatable, target :: T_MEAN   (:)     !> horiz. mean of t               [K]
  real(RP), allocatable, target :: QV_MEAN  (:)     !> horiz. mean of QV
  real(RP), allocatable, target :: QHYD_MEAN(:)     !> horiz. mean of QHYD
  real(RP), allocatable, target :: QLIQ_MEAN(:)     !> horiz. mean of QLIQ
  real(RP), allocatable, target :: QICE_MEAN(:)     !> horiz. mean of QICE

  real(RP), allocatable, target :: DENS_PRIM(:,:,:,:) !> horiz. deviation of density    [kg/m3]
  real(RP), allocatable, target :: W_PRIM   (:,:,:,:) !> horiz. deviation of w          [m/s]
  real(RP), allocatable, target :: U_PRIM   (:,:,:,:) !> horiz. deviation of u          [m/s]
  real(RP), allocatable, target :: V_PRIM   (:,:,:,:) !> horiz. deviation of v          [m/s]
  real(RP), allocatable, target :: PT_PRIM  (:,:,:,:) !> horiz. deviation of pot. temp. [K]
  real(RP), allocatable, target :: W_PRIM2  (:,:,:,:) !> variance of w                  [m2/s2]
  real(RP), allocatable, target :: PT_W_PRIM(:,:,:,:) !> resolved scale heat flux       [W/s]
  real(RP), allocatable, target :: W_PRIM3  (:,:,:,:) !> skewness of w                  [m3/s3]
  real(RP), allocatable, target :: TKE_RS   (:,:,:,:) !> resolved scale TKE             [m2/s2]

  ! id of diagnostic variables
  !! public
  integer,     private, parameter :: I_W         =  1
  integer,     private, parameter :: I_U         =  2
  integer,     private, parameter :: I_V         =  3
  integer,     private, parameter :: I_POTT      =  4
  integer,     private, parameter :: I_TEMP      =  5
  integer,     private, parameter :: I_PRES      =  6
  integer,     private, parameter :: I_EXNER     =  7
  integer,     private, parameter :: I_PHYD      =  8
  integer,     private, parameter :: I_QDRY      =  9
  integer,     private, parameter :: I_RTOT      = 10
  integer,     private, parameter :: I_CVTOT     = 11
  integer,     private, parameter :: I_CPTOT     = 12
  !! private
  integer,     private, parameter :: I_LHV       = 13
  integer,     private, parameter :: I_LHS       = 14
  integer,     private, parameter :: I_LHF       = 15
  integer,     private, parameter :: I_POTV      = 16
  integer,     private, parameter :: I_TEML      = 17
  integer,     private, parameter :: I_POTL      = 18
  integer,     private, parameter :: I_POTE      = 19
  integer,     private, parameter :: I_QTOT      = 20
  integer,     private, parameter :: I_QHYD      = 21
  integer,     private, parameter :: I_QLIQ      = 22
  integer,     private, parameter :: I_QICE      = 23
  integer,     private, parameter :: I_LWP       = 24
  integer,     private, parameter :: I_IWP       = 25
  integer,     private, parameter :: I_PW        = 26
  integer,     private, parameter :: I_PREC      = 27
  integer,     private, parameter :: I_RAIN      = 28
  integer,     private, parameter :: I_SNOW      = 29
  integer,     private, parameter :: I_QSAT      = 30
  integer,     private, parameter :: I_RHA       = 31
  integer,     private, parameter :: I_RHL       = 32
  integer,     private, parameter :: I_RHI       = 33
  integer,     private, parameter :: I_VOR       = 34
  integer,     private, parameter :: I_DIV       = 35
  integer,     private, parameter :: I_HDIV      = 36
  integer,     private, parameter :: I_Uabs      = 37
  integer,     private, parameter :: I_N2        = 38
  integer,     private, parameter :: I_PBLH      = 39
  integer,     private, parameter :: I_MSE       = 40
  integer,     private, parameter :: I_TDEW      = 41
  integer,     private, parameter :: I_CAPE      = 42
  integer,     private, parameter :: I_CIN       = 43
  integer,     private, parameter :: I_LCL       = 44
  integer,     private, parameter :: I_LFC       = 45
  integer,     private, parameter :: I_LNB       = 46
  integer,     private, parameter :: I_ENGT      = 47
  integer,     private, parameter :: I_ENGP      = 48
  integer,     private, parameter :: I_ENGK      = 49
  integer,     private, parameter :: I_ENGI      = 50
  integer,     private, parameter :: I_DENS_MEAN = 51
  integer,     private, parameter :: I_W_MEAN    = 52
  integer,     private, parameter :: I_U_MEAN    = 53
  integer,     private, parameter :: I_V_MEAN    = 54
  integer,     private, parameter :: I_PT_MEAN   = 55
  integer,     private, parameter :: I_T_MEAN    = 56
  integer,     private, parameter :: I_QV_MEAN   = 57
  integer,     private, parameter :: I_QHYD_MEAN = 58
  integer,     private, parameter :: I_QLIQ_MEAN = 59
  integer,     private, parameter :: I_QICE_MEAN = 60
  integer,     private, parameter :: I_DENS_PRIM = 61
  integer,     private, parameter :: I_W_PRIM    = 62
  integer,     private, parameter :: I_U_PRIM    = 63
  integer,     private, parameter :: I_V_PRIM    = 64
  integer,     private, parameter :: I_PT_PRIM   = 65
  integer,     private, parameter :: I_W_PRIM2   = 66
  integer,     private, parameter :: I_PT_W_PRIM = 67
  integer,     private, parameter :: I_W_PRIM3   = 68
  integer,     private, parameter :: I_TKE_RS    = 69

  integer,     private, parameter :: DV_nmax     = 69
  type(Vinfo), private            :: DV_info(DV_nmax)
  logical,     private            :: DV_calclated(DV_nmax)

  data DV_info / &
       Vinfo( 'W',         'velocity w',                      'm/s',     3, 'ZXY', 'upward_air_velocity' ), &
       Vinfo( 'U',         'velocity u',                      'm/s',     3, 'ZXY', 'eastward_air_velocity' ), &
       Vinfo( 'V',         'velocity v',                      'm/s',     3, 'ZXY', 'northward_air_velocity' ), &
       Vinfo( 'PT',        'potential temp.',                 'K',       3, 'ZXY', 'air_potential_temperature' ), &
       Vinfo( 'T',         'temperature',                     'K',       3, 'ZXY', 'air_temperature' ), &
       Vinfo( 'PRES',      'pressure',                        'Pa',      3, 'ZXY', 'air_pressure' ), &
       Vinfo( 'EXNER',     'Exner function',                  '1',       3, 'ZXY', 'dimensionless_exner_function' ), &
       Vinfo( 'PHYD',      'hydrostatic pressure',            'Pa',      3, 'ZXY', '' ), &
       Vinfo( 'QDRY',      'dry air',                         'kg/kg',   3, 'ZXY', '' ), &
       Vinfo( 'RTOT',      'Total gas constant',              'J/kg/K',  3, 'ZXY', '' ), &
       Vinfo( 'CVTOT',     'Total heat capacity',             'J/kg/K',  3, 'ZXY', '' ), &
       Vinfo( 'CPTOT',     'Total heat capacity',             'J/kg/K',  3, 'ZXY', '' ), &
       Vinfo( 'LHV',       'latent heat for vaporization',    'J/kg',    3, 'ZXY', '' ), &
       Vinfo( 'LHS',       'latent heat for sublimation',     'J/kg',    3, 'ZXY', '' ), &
       Vinfo( 'LHF',       'latent heat for fusion',          'J/kg',    3, 'ZXY', '' ), &
       Vinfo( 'POTV',      'virtual potential temp.',         'K',       3, 'ZXY', '' ), &
       Vinfo( 'TEML',      'liquid water temperature',        'K',       3, 'ZXY', '' ), &
       Vinfo( 'POTL',      'liquid water potential temp.',    'K',       3, 'ZXY', '' ), &
       Vinfo( 'POTE',      'equivalent potential temp.',      'K',       3, 'ZXY', 'pseudo_equivalent_potential_temperature'  ), &
       Vinfo( 'QTOT',      'total water',                     'kg/kg',   3, 'ZXY', 'mass_fraction_of_water_in_air' ), &
       Vinfo( 'QHYD',      'total hydrometeors',              'kg/kg',   3, 'ZXY', 'mass_fraction_of_cloud_condensed_water_in_air' ), &
       Vinfo( 'QLIQ',      'total liquid water',              'kg/kg',   3, 'ZXY', '' ), &
       Vinfo( 'QICE',      'total ice water',                 'kg/kg',   3, 'ZXY', '' ), &
       Vinfo( 'LWP',       'liquid water path',               'g/m2',    2, 'XY',  'atmosphere_mass_content_of_cloud_liquid_water' ), &
       Vinfo( 'IWP',       'ice water path',                  'g/m2',    2, 'XY',  '' ), &
       Vinfo( 'PW',        'precipitable water',              'g/m2',    2, 'XY',  'atmosphere_mass_content_of_vapor' ), &
       Vinfo( 'PREC',      'surface precipitation flux',      'kg/m2/s', 2, 'XY',  'precipitation_flux' ), &
       Vinfo( 'RAIN',      'surface rain flux',               'kg/m2/s', 2, 'XY',  'rainfall_flux' ), &
       Vinfo( 'SNOW',      'surface snow flux',               'kg/m2/s', 2, 'XY',  'snowfall_flux' ), &
       Vinfo( 'QSAT',      'saturation specific humidity',    'kg/kg',   3, 'ZXY', '' ), &
       Vinfo( 'RHA',       'relative humidity(liq+ice)',      '%',       3, 'ZXY', '' ), &
       Vinfo( 'RH',        'relative humidity(liq)',          '%',       3, 'ZXY', 'relative_humidity' ), &
       Vinfo( 'RHI',       'relative humidity(ice)',          '%',       3, 'ZXY', '' ), &
       Vinfo( 'VOR',       'vertical vorticity',              '1/s',     3, 'ZXY', 'atmosphere_relative_vorticity' ), &
       Vinfo( 'DIV',       'divergence',                      '1/s',     3, 'ZXY', 'divergence_of_wind' ), &
       Vinfo( 'HDIV',      'horizontal divergence',           '1/s',     3, 'ZXY', '' ), &
       Vinfo( 'Uabs',      'absolute velocity',               'm/s',     3, 'ZXY', '' ), &
       Vinfo( 'N2',        'squared Brunt-Vaisala frequency', '1/s2',    3, 'ZXY', 'square_of_brunt_vaisala_frequency_in_air' ), &
       Vinfo( 'PBLH',      'PBL height',                      'm',       2, 'XY', 'atmosphere_boundary_layer_thickness'  ), &
       Vinfo( 'MSE',       'moist static energy',             'm2/s2',   3, 'ZXY', '' ), &
       Vinfo( 'TDEW',      'dew point',                       'K',       3, 'ZXY', 'dew_point_temperature' ), &
       Vinfo( 'CAPE',      'convective avail. pot. energy',   'm2/s2',   2, 'XY',  'atmosphere_specific_convective_available_potential_energy'  ), &
       Vinfo( 'CIN',       'convection inhibition',           'm2/s2',   2, 'XY',  '' ), &
       Vinfo( 'LCL',       'lifted condensation level',       'm',       2, 'XY',  'atmosphere_lifting_condensation_level' ), &
       Vinfo( 'LFC',       'level of free convection',        'm',       2, 'XY',  'atmosphere_level_of_free_convection' ), &
       Vinfo( 'LNB',       'level of neutral buoyancy',       'm',       2, 'XY',  '' ), &
       Vinfo( 'ENGT',      'total energy',                    'J/m3',    3, 'ZXY', '' ), &
       Vinfo( 'ENGP',      'potential energy',                'J/m3',    3, 'ZXY', '' ), &
       Vinfo( 'ENGK',      'kinetic energy',                  'J/m3',    3, 'ZXY', '' ), &
       Vinfo( 'ENGI',      'internal energy',                 'J/m3',    3, 'ZXY', '' ), &
       Vinfo( 'DENS_MEAN', 'horiz. mean of density',          'kg/m3',   1, 'Z',   '' ), &
       Vinfo( 'W_MEAN',    'horiz. mean of w',                'm/s',     1, 'Z',   '' ), &
       Vinfo( 'U_MEAN',    'horiz. mean of u',                'm/s',     1, 'Z',   '' ), &
       Vinfo( 'V_MEAN',    'horiz. mean of v',                'm/s',     1, 'Z',   '' ), &
       Vinfo( 'PT_MEAN',   'horiz. mean of pot.',             'K',       1, 'Z',   '' ), &
       Vinfo( 'T_MEAN',    'horiz. mean of t',                'K',       1, 'Z',   '' ), &
       Vinfo( 'QV_MEAN',   'horiz. mean of QV',               '1',       1, 'Z',   '' ), &
       Vinfo( 'QHYD_MEAN', 'horiz. mean of QHYD',             '1',       1, 'Z',   '' ), &
       Vinfo( 'QLIQ_MEAN', 'horiz. mean of QLIQ',             '1',       1, 'Z',   '' ), &
       Vinfo( 'QICE_MEAN', 'horiz. mean of QICE',             '1',       1, 'Z',   '' ), &
       Vinfo( 'DENS_PRIM', 'horiz. deviation of density',     'kg/m3',   3, 'ZXY', '' ), &
       Vinfo( 'W_PRIM',    'horiz. deviation of w',           'm/s',     3, 'ZXY', '' ), &
       Vinfo( 'U_PRIM',    'horiz. deviation of u',           'm/s',     3, 'ZXY', '' ), &
       Vinfo( 'V_PRIM',    'horiz. deviation of v',           'm/s',     3, 'ZXY', '' ), &
       Vinfo( 'PT_PRIM',   'horiz. deviation of pot. temp.',  'K',       3, 'ZXY', '' ), &
       Vinfo( 'W_PRIM2',   'variance of w',                   'm2/s2',   3, 'ZXY', '' ), &
       Vinfo( 'PT_W_PRIM', 'resolved scale heat flux',        'W/s',     3, 'ZXY', '' ), &
       Vinfo( 'W_PRIM3',   'skewness of w',                   'm3/s3',   3, 'ZXY', '' ), &
       Vinfo( 'TKE_RS',    'resolved scale TKE',              'm2/s2',   3, 'ZXY', '' ) /

  ! for history output and monitor
  integer, private, allocatable :: QP_HIST_id (:)       !> tracer variables
  integer, private, allocatable :: QP_MONIT_id(:)
  integer, private              :: DV_HIST_id (DV_nmax) !> diagnostic variables

  integer, private, parameter   :: IM_QDRY         =  1
  integer, private, parameter   :: IM_QTOT         =  2
  integer, private, parameter   :: IM_EVAP         =  3
  integer, private, parameter   :: IM_PREC         =  4
  integer, private, parameter   :: IM_ENGT         =  5
  integer, private, parameter   :: IM_ENGP         =  6
  integer, private, parameter   :: IM_ENGK         =  7
  integer, private, parameter   :: IM_ENGI         =  8
  integer, private, parameter   :: IM_ENGFLXT      =  9
  integer, private, parameter   :: IM_ENGSFC_SH    = 10
  integer, private, parameter   :: IM_ENGSFC_LH    = 11
  integer, private, parameter   :: IM_ENGSFC_RD    = 12
  integer, private, parameter   :: IM_ENGTOA_RD    = 13
  integer, private, parameter   :: IM_ENGSFC_LW_up = 14
  integer, private, parameter   :: IM_ENGSFC_LW_dn = 15
  integer, private, parameter   :: IM_ENGSFC_SW_up = 16
  integer, private, parameter   :: IM_ENGSFC_SW_dn = 17
  integer, private, parameter   :: IM_ENGTOA_LW_up = 18
  integer, private, parameter   :: IM_ENGTOA_LW_dn = 19
  integer, private, parameter   :: IM_ENGTOA_SW_up = 20
  integer, private, parameter   :: IM_ENGTOA_SW_dn = 21
  integer, private, parameter   :: DVM_nmax        = 21
  integer, private              :: DV_MONIT_id(DVM_nmax)


  logical,  private                      :: moist
  real(RP), private, target, allocatable :: Qe(:,:,:,:,:)  !> mass ratio of hydrometors [kg/kg]
  real(RP), private, target, allocatable :: ZERO(:,:,:,:)


  ! for restart
  integer, private :: restart_fid = -1  ! file ID
  logical, private :: ATMOS_RESTART_IN_CHECK_COORDINATES = .true.


  real(RP), private, allocatable :: WORK3D(:,:,:,:)
  real(RP), private, allocatable :: WORK2D(:,:,:)
  real(RP), private, allocatable :: WORK1D(:)

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine ATMOS_vars_setup
    use scale_const, only: &
       UNDEF => CONST_UNDEF
    use scale_prc, only: &
       PRC_abort
    use scale_file_history, only: &
       FILE_HISTORY_reg
    use scale_monitor, only: &
       MONITOR_reg
    use scale_atmos_hydrometeor, only: &
       ATMOS_HYDROMETEOR_dry, &
       N_HYD, &
       I_QV, &
       I_HC, &
       I_HR, &
       I_HI, &
       I_HS, &
       I_HG, &
       I_HH
    use mod_atmos_phy_mp_vars, only: &
       ATMOS_PHY_MP_vars_setup
    use mod_atmos_phy_ae_vars, only: &
       ATMOS_PHY_AE_vars_setup
    use mod_atmos_phy_rd_vars, only: &
       ATMOS_PHY_RD_vars_setup
    use mod_atmos_phy_sf_vars, only: &
       ATMOS_PHY_SF_vars_setup
    use mod_atmos_phy_bl_vars, only: &
       ATMOS_PHY_BL_vars_setup
    use mod_grd, only: &
       I_LON,  &
       I_LAT,  &
       GRD_s,    &
       GRD_Z,    &
       GRD_ZH,   &
       GRD_vz,   &
       GRD_ZSFC, &
       GRD_zs
    implicit none

    namelist / PARAM_ATMOS_VARS / &
       ATMOS_RESTART_IN_BASENAME,           &
       ATMOS_RESTART_IN_AGGREGATE,          &
       ATMOS_RESTART_IN_POSTFIX_TIMELABEL,  &
       ATMOS_RESTART_IN_CHECK_COORDINATES,  &
       ATMOS_RESTART_OUTPUT,                &
       ATMOS_RESTART_OUT_BASENAME,          &
       ATMOS_RESTART_OUT_AGGREGATE,         &
       ATMOS_RESTART_OUT_POSTFIX_TIMELABEL, &
       ATMOS_RESTART_OUT_TITLE,             &
       ATMOS_RESTART_OUT_DTYPE,             &
       ATMOS_RESTART_CHECK,                 &
       ATMOS_RESTART_CHECK_BASENAME,        &
       ATMOS_RESTART_CHECK_CRITERION,       &
       ATMOS_VARS_CHECKRANGE,               &
       ATMOS_VARS_CHECKCFL_SOFT,            &
       ATMOS_VARS_CHECKCFL_HARD

    integer :: ierr
    integer :: k, i, j, l, ij
    integer :: iv, iq
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("ATMOS_vars_setup",*) 'Setup'

    allocate( DENS(KA,IA,JA,ADM_lall) )
    allocate( RHOU(KA,IA,JA,ADM_lall) )
    allocate( RHOV(KA,IA,JA,ADM_lall) )
    allocate( MOMZ(KA,IA,JA,ADM_lall) )
    allocate( RHOE(KA,IA,JA,ADM_lall) )
    allocate( RHOQ(KA,IA,JA,QA,ADM_lall) )
    allocate( QTRC(KA,IA,JA,QA,ADM_lall) )
    DENS(:,:,:,:) = UNDEF
    RHOU(:,:,:,:) = UNDEF
    RHOV(:,:,:,:) = UNDEF
    MOMZ(:,:,:,:) = UNDEF
    RHOE(:,:,:,:) = UNDEF
    RHOQ(:,:,:,:,:) = UNDEF
    QTRC(:,:,:,:,:) = UNDEF


    allocate( W(KA,IA,JA,ADM_lall) )
    allocate( U(KA,IA,JA,ADM_lall) )
    allocate( V(KA,IA,JA,ADM_lall) )
    W(:,:,:,:) = UNDEF
    U(:,:,:,:) = UNDEF
    V(:,:,:,:) = UNDEF

    allocate( POTT (KA,IA,JA,ADM_lall) )
    allocate( TEMP (KA,IA,JA,ADM_lall) )
    allocate( PRES (KA,IA,JA,ADM_lall) )
    allocate( EXNER(KA,IA,JA,ADM_lall) )
    allocate( PHYD (KA,IA,JA,ADM_lall) )
    allocate( PHYDH(0:KA,IA,JA,ADM_lall) )
    POTT (:,:,:,:) = UNDEF
    TEMP (:,:,:,:) = UNDEF
    PRES (:,:,:,:) = UNDEF
    EXNER(:,:,:,:) = UNDEF
    PHYD (:,:,:,:) = UNDEF
    PHYDH(:,:,:,:) = UNDEF

    allocate( Qdry (KA,IA,JA,ADM_lall) )
    allocate( Rtot (KA,IA,JA,ADM_lall) )
    allocate( CVtot(KA,IA,JA,ADM_lall) )
    allocate( CPtot(KA,IA,JA,ADM_lall) )
    Qdry (:,:,:,:) = UNDEF
    Rtot (:,:,:,:) = UNDEF
    CVtot(:,:,:,:) = UNDEF
    CPtot(:,:,:,:) = UNDEF

    allocate( WORK3D(KA,IA,JA,ADM_lall) )
    allocate( WORK2D(   IA,JA,ADM_lall) )
    allocate( WORK1D(KA      ) )


    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_ATMOS_VARS,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("ATMOS_vars_setup",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("ATMOS_vars_setup",*) 'Not appropriate names in namelist PARAM_ATMOS_VARS. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_ATMOS_VARS)

    LOG_NEWLINE
    LOG_INFO("ATMOS_vars_setup",*) 'List of prognostic variables (ATMOS) '
    LOG_INFO_CONT('(1x,A,A24,A,A48,A,A12,A)') &
               '      |', 'VARNAME                 ','|', &
               'DESCRIPTION                                     ', '[', 'UNIT        ', ']'
    do iq = 1, QA
       LOG_INFO_CONT('(1x,A,I3,A,A24,A,A48,A,A12,A)') &
                  'NO.',5+iq,'|',TRACER_NAME(iq),'|', TRACER_DESC(iq),'[', TRACER_UNIT(iq),']'
    enddo

    LOG_NEWLINE
    if ( ATMOS_RESTART_IN_BASENAME /= '' ) then
       LOG_INFO("ATMOS_vars_setup",*) 'Restart input?  : YES, file = ', trim(ATMOS_RESTART_IN_BASENAME)
       LOG_INFO("ATMOS_vars_setup",*) 'Add timelabel?  : ', ATMOS_RESTART_IN_POSTFIX_TIMELABEL
    else
       LOG_INFO("ATMOS_vars_setup",*) 'Restart input?  : NO'
    endif
    if (       ATMOS_RESTART_OUTPUT             &
         .AND. ATMOS_RESTART_OUT_BASENAME /= '' ) then
       LOG_INFO("ATMOS_vars_setup",*) 'Restart output? : YES, file = ', trim(ATMOS_RESTART_OUT_BASENAME)
       LOG_INFO("ATMOS_vars_setup",*) 'Add timelabel?  : ', ATMOS_RESTART_OUT_POSTFIX_TIMELABEL
    else
       LOG_INFO("ATMOS_vars_setup",*) 'Restart output? : NO'
       ATMOS_RESTART_OUTPUT = .false.
    endif

    if ( ATMOS_RESTART_CHECK_BASENAME == '' ) then
       ATMOS_RESTART_CHECK = .false.
    endif

    if ( ATMOS_VARS_CHECKCFL_HARD > 0.0_RP ) then
       ATMOS_VARS_CHECKCFL_SOFT = min( ATMOS_VARS_CHECKCFL_SOFT, ATMOS_VARS_CHECKCFL_HARD )
    endif

    LOG_NEWLINE
    LOG_INFO("ATMOS_vars_setup",*) 'Check restart consistency?          : ', ATMOS_RESTART_CHECK
    LOG_INFO("ATMOS_vars_setup",*) 'Check value range of variables?     : ', ATMOS_VARS_CHECKRANGE
    if ( ATMOS_VARS_CHECKCFL_SOFT > 0.0_RP ) then
       LOG_INFO("ATMOS_vars_setup",*) 'Threshold of Courant number to warn : ', ATMOS_VARS_CHECKCFL_SOFT
    else
       LOG_INFO("ATMOS_vars_setup",*) 'Threshold of Courant number to warn : disabled'
    endif
    if ( ATMOS_VARS_CHECKCFL_HARD > 0.0_RP ) then
       LOG_INFO("ATMOS_vars_setup",*) 'Threshold of Courant number to stop : ', ATMOS_VARS_CHECKCFL_HARD
    else
       LOG_INFO("ATMOS_vars_setup",*) 'Threshold of Courant number to stop : disabled'
    endif

    call ATMOS_PHY_MP_vars_setup
    call ATMOS_PHY_AE_vars_setup
    call ATMOS_PHY_RD_vars_setup
    call ATMOS_PHY_SF_vars_setup
    call ATMOS_PHY_BL_vars_setup


    ! water content
    if ( ATMOS_HYDROMETEOR_dry ) then
       allocate( ZERO(KA,IA,JA,ADM_lall) )
!OCL XFILL
       ZERO(:,:,:,:) = 0.0_RP

       QV => ZERO
       QC => ZERO
       QR => ZERO
       QI => ZERO
       QS => ZERO
       QG => ZERO
       QH => ZERO

       moist = .false.
    else
       allocate( Qe(KA,IA,JA,ADM_lall,N_HYD) )
!OCL XFILL
       Qe(:,:,:,:,:) = UNDEF

       QV => QTRC(:,:,:,I_QV,:)
       QC => Qe(:,:,:,:,I_HC)
       QR => Qe(:,:,:,:,I_HR)
       QI => Qe(:,:,:,:,I_HI)
       QS => Qe(:,:,:,:,I_HS)
       QG => Qe(:,:,:,:,I_HG)
       QH => Qe(:,:,:,:,I_HH)

       moist = .true.
    end if


    ! tentative

    allocate( CZ(  KA,IA,JA,ADM_lall) )
    allocate( FZ(0:KA,IA,JA,ADM_lall) )
    do l = 1, ADM_lall
    do j = 1, JA
    do i = 1, IA
       ij = i + ADM_imin - 1 + ( j - 1 ) * ADM_imax
       do k = 1, KA
          CZ(k,i,j,l) = GRD_vz(ij,k,l,GRD_Z)
       end do
       do k = 0, KA-1
          FZ(k,i,j,l) = GRD_vz(ij,k+1,l,GRD_ZH)
       end do
       FZ(KA,i,j,l) = FZ(KA-1,i,j,l) + ( FZ(KA-1,i,j,l) - FZ(KA-2,i,j,l) )
    end do
    end do
    end do

    allocate( LON(IA,JA,ADM_lall) )
    allocate( LAT(IA,JA,ADM_lall) )
    do l = 1, ADM_lall
    do j = 1, JA
    do i = 1, IA
       ij = i + ADM_imin - 1 + ( j - 1 ) * ADM_iall
       LON(i,j,l) = GRD_s(ij,ADM_KNONE,l,I_LON)
       LAT(i,j,l) = GRD_s(ij,ADM_KNONE,l,I_LAT)
    end do
    end do
    end do


    allocate( Z1       (IA,JA,ADM_lall) )
    allocate( TOPO_Zsfc(IA,JA,ADM_lall) )
    do l = 1, ADM_lall
    do j = 1, JA
    do i = 1, IA
       ij = i + ADM_imin - 1 + ( j - 1 ) * ADM_iall
       TOPO_Zsfc(i,j,l) = GRD_zs(ij,ADM_KNONE,l,GRD_ZSFC)
       Z1       (i,j,l) = CZ(KS,i,j,l) - TOPO_Zsfc(i,j,l)
    end do
    end do
    end do


    allocate( fact_ocean(IA,JA,ADM_lall) )
    allocate( fact_land (IA,JA,ADM_lall) )
    allocate( fact_urban(IA,JA,ADM_lall) )
    fact_ocean(:,:,:) = 1.0_RP
    fact_land (:,:,:) = 0.0_RP
    fact_urban(:,:,:) = 0.0_RP


    DV_calclated(DV_nmax) = .false.

    !-----< history output setup >-----
    allocate( QP_HIST_id( max(QA,1) ) )
    allocate( QP_MONIT_id( max(QA,1) ) )
    QP_HIST_id (:) = -1
    QP_MONIT_id(:) = -1
    DV_HIST_id (:) = -1
    DV_MONIT_id(:) = -1


    do iv = 1, DV_nmax
       call FILE_HISTORY_reg( DV_info(iv)%NAME, DV_info(iv)%DESC, DV_info(iv)%UNIT, DV_HIST_id(iv), dim_type=DV_info(iv)%dim_type, standard_name=DV_info(iv)%STDNAME )
    end do


    !-----< monitor output setup >-----
    do iq = 1, QA
       call MONITOR_reg( TRACER_NAME(iq), TRACER_DESC(iq), TRACER_UNIT(iq), & ! (in)
                         QP_MONIT_id(iq),                                   & ! (out)
                         dim_type='ZXY', isflux=.false.                     ) ! (in)
    enddo

    call MONITOR_reg( 'QDRY',         'dry air mass',           'kg', & ! (in)
                      DV_MONIT_id(IM_QDRY),                           & ! (out)
                      dim_type='ZXY', isflux=.false.                  ) ! (in)
    call MONITOR_reg( 'QTOT',         'water mass',             'kg', & ! (in)
                      DV_MONIT_id(IM_QTOT),                           & ! (out)
                      dim_type='ZXY', isflux=.false.                  ) ! (in)
    call MONITOR_reg( 'EVAP',         'evaporation',            'kg', & ! (in)
                      DV_MONIT_id(IM_EVAP),                           & ! (out)
                      dim_type='XY', isflux=.true.                    ) ! (in)
    call MONITOR_reg( 'PRCP',         'precipitation',          'kg', & ! (in)
                      DV_MONIT_id(IM_PREC),                           & ! (out)
                      dim_type='XY', isflux=.true.                    ) ! (in)

    call MONITOR_reg( 'ENGT',         'total     energy',       'J', & ! (in)
                      DV_MONIT_id(IM_ENGT),                          & ! (out)
                      dim_type='ZXY', isflux=.false.                 ) ! (in)
    call MONITOR_reg( 'ENGP',         'potential energy',       'J', & ! (in)
                      DV_MONIT_id(IM_ENGP),                          & ! (out)
                      dim_type='ZXY', isflux=.false.                 ) ! (in)
    call MONITOR_reg( 'ENGK',         'kinetic   energy',       'J', & ! (in)
                      DV_MONIT_id(IM_ENGK),                          & ! (out)
                      dim_type='ZXY', isflux=.false.                 ) ! (in)
    call MONITOR_reg( 'ENGI',         'internal  energy',       'J', & ! (in)
                      DV_MONIT_id(IM_ENGI),                          & ! (out)
                      dim_type='ZXY', isflux=.false.                 ) ! (in)

    call MONITOR_reg( 'ENGFLXT',      'total energy flux',      'J', & ! (in)
                      DV_MONIT_id(IM_ENGFLXT),                       & ! (out)
                      dim_type='XY', isflux=.true.                   ) ! (in)
    call MONITOR_reg( 'ENGSFC_SH',    'SFC specific heat flux', 'J', & ! (in)
                      DV_MONIT_id(IM_ENGSFC_SH),                     & ! (out)
                      dim_type='XY', isflux=.true.                   ) ! (in)
    call MONITOR_reg( 'ENGSFC_LH',    'SFC latent   heat flux', 'J', & ! (in)
                      DV_MONIT_id(IM_ENGSFC_LH),                     & ! (out)
                      dim_type='XY', isflux=.true.                   ) ! (in)
    call MONITOR_reg( 'ENGSFC_RD',    'SFC net radiation flux', 'J', & ! (in)
                      DV_MONIT_id(IM_ENGSFC_RD),                     & ! (out)
                      dim_type='XY', isflux=.true.                   ) ! (in)
    call MONITOR_reg( 'ENGTOA_RD',    'TOA net radiation flux', 'J', & ! (in)
                      DV_MONIT_id(IM_ENGTOA_RD),                     & ! (out)
                      dim_type='XY', isflux=.true.                   ) ! (in)

    call MONITOR_reg( 'ENGSFC_LW_up', 'SFC LW upward   flux',   'J', & ! (in)
                      DV_MONIT_id(IM_ENGSFC_LW_up),                  & ! (out)
                      dim_type='XY', isflux=.true.                   ) ! (in)
    call MONITOR_reg( 'ENGSFC_LW_dn', 'SFC LW downward flux',   'J', & ! (in)
                      DV_MONIT_id(IM_ENGSFC_LW_dn),                  & ! (out)
                      dim_type='XY', isflux=.true.                   ) ! (in)
    call MONITOR_reg( 'ENGSFC_SW_up', 'SFC SW upward   flux',   'J', & ! (in)
                      DV_MONIT_id(IM_ENGSFC_SW_up),                  & ! (out)
                      dim_type='XY', isflux=.true.                   ) ! (in)
    call MONITOR_reg( 'ENGSFC_SW_dn', 'SFC SW downward flux',   'J', & ! (in)
                      DV_MONIT_id(IM_ENGSFC_SW_dn),                  & ! (out)
                      dim_type='XY', isflux=.true.                   ) ! (in)

    call MONITOR_reg( 'ENGTOA_LW_up', 'TOA LW upward   flux',   'J', & ! (in)
                      DV_MONIT_id(IM_ENGTOA_LW_up),                  & ! (out)
                      dim_type='XY', isflux=.true.                   ) ! (in)
    call MONITOR_reg( 'ENGTOA_LW_dn', 'TOA LW downward flux',   'J', & ! (in)
                      DV_MONIT_id(IM_ENGTOA_LW_dn),                  & ! (out)
                      dim_type='XY', isflux=.true.                   ) ! (in)
    call MONITOR_reg( 'ENGTOA_SW_up', 'TOA SW upward   flux',   'J', & ! (in)
                      DV_MONIT_id(IM_ENGTOA_SW_up),                  & ! (out)
                      dim_type='XY', isflux=.true.                   ) ! (in)
    call MONITOR_reg( 'ENGTOA_SW_dn', 'TOA SW downward flux',   'J', & ! (in)
                      DV_MONIT_id(IM_ENGTOA_SW_dn),                  & ! (out)
                      dim_type='XY', isflux=.true.                   ) ! (in)

    return
  end subroutine ATMOS_vars_setup

  !-----------------------------------------------------------------------------
  !> HALO Communication
  subroutine ATMOS_vars_fillhalo

    return
  end subroutine ATMOS_vars_fillhalo

  !-----------------------------------------------------------------------------
  !> Open restart file for reading atmospheric variables
  subroutine ATMOS_vars_restart_open

    return
  end subroutine ATMOS_vars_restart_open

  !-----------------------------------------------------------------------------
  !> Read restart of atmospheric variables
  subroutine ATMOS_vars_restart_read

    return
  end subroutine ATMOS_vars_restart_read

  !-----------------------------------------------------------------------------
  !> Set pressure for history output
  subroutine ATMOS_vars_history_setpres

!!$    call FILE_HISTORY_CARTESC_set_pres( PHYD    (:,:,:,:), & ! [IN]
!!$                                        PHYDH   (:,:,:,:), & ! [IN]
!!$                                        SFC_PRES(:,:,:)    ) ! [IN]

    return
  end subroutine ATMOS_vars_history_setpres

  !-----------------------------------------------------------------------------
  !> Check and compare between last data and sample data
  subroutine ATMOS_vars_restart_check

    return
  end subroutine ATMOS_vars_restart_check

  !-----------------------------------------------------------------------------
  !> History output set for atmospheric variables
  subroutine ATMOS_vars_history
    use scale_file_history, only: &
       FILE_HISTORY_query, &
       FILE_HISTORY_put
    use mod_atmos_phy_mp_vars, only: &
       ATMOS_PHY_MP_vars_history
    use mod_atmos_phy_ae_vars, only: &
       ATMOS_PHY_AE_vars_history
    implicit none

    logical :: do_put
    integer :: iq, iv
    !---------------------------------------------------------------------------

    call PROF_rapstart('ATM_History', 1)


    ! history output of diagnostic variables
    call FILE_HISTORY_put  ( DV_HIST_id(I_W    ), W(:,:,:,:)     )
    call FILE_HISTORY_put  ( DV_HIST_id(I_U    ), U(:,:,:,:)     )
    call FILE_HISTORY_put  ( DV_HIST_id(I_V    ), V(:,:,:,:)     )
    call FILE_HISTORY_put  ( DV_HIST_id(I_POTT ), POTT(:,:,:,:)  )
    call FILE_HISTORY_put  ( DV_HIST_id(I_TEMP ), TEMP(:,:,:,:)  )
    call FILE_HISTORY_put  ( DV_HIST_id(I_PRES ), PRES(:,:,:,:)  )

    call FILE_HISTORY_put  ( DV_HIST_id(I_EXNER), EXNER(:,:,:,:) )
    call FILE_HISTORY_put  ( DV_HIST_id(I_PHYD ), PHYD(:,:,:,:)  )

    call FILE_HISTORY_put  ( DV_HIST_id(I_QDRY ), QDRY(:,:,:,:)  )
    call FILE_HISTORY_put  ( DV_HIST_id(I_RTOT ), RTOT(:,:,:,:)  )
    call FILE_HISTORY_put  ( DV_HIST_id(I_CVTOT), CVTOT(:,:,:,:) )
    call FILE_HISTORY_put  ( DV_HIST_id(I_CPTOT), CPTOT(:,:,:,:) )

    do iv = I_CPTOT+1, DV_nmax
       if ( DV_HIST_id(iv) > 0 ) then
          call FILE_HISTORY_query( DV_HIST_id(iv), do_put )

          if ( do_put ) then
             select case( DV_info(iv)%ndims )
             case( 3 )
                call ATMOS_vars_get_diagnostic( DV_info(iv)%NAME, WORK3D(:,:,:,:) )
                call FILE_HISTORY_put( DV_HIST_id(iv), WORK3D(:,:,:,:) )
             case( 2 )
                call ATMOS_vars_get_diagnostic( DV_info(iv)%NAME, WORK2D(:,:,:) )
                call FILE_HISTORY_put( DV_HIST_id(iv), WORK2D(:,:,:) )
             case( 1 )
                call ATMOS_vars_get_diagnostic( DV_info(iv)%NAME, WORK1D(:) )
                call FILE_HISTORY_put( DV_HIST_id(iv), WORK1D(:) )
             end select
          endif
       endif
    enddo

    if ( moist ) &
         call ATMOS_PHY_MP_vars_history( DENS(:,:,:,:), TEMP(:,:,:,:), QTRC(:,:,:,:,:) )
!    if ( .false. ) then
!       call ATMOS_vars_get_diagnostic( "RH", WORK3D(:,:,:,:) )
!       call ATMOS_PHY_AE_vars_history( QTRC(:,:,:,:,:), WORK3D(:,:,:,:) )
!    end if

    call PROF_rapend  ('ATM_History', 1)

    return
  end subroutine ATMOS_vars_history

  !-----------------------------------------------------------------------------
  !> Budget monitor for atmosphere
  subroutine ATMOS_vars_total

    return
  end subroutine ATMOS_vars_total

  !-----------------------------------------------------------------------------
  !> Calc diagnostic variables from IcoA grid data
  !< outpus:  DENS, RHOU, RHOV, MOMZ, RHOE, RHOQ
  subroutine ATMOS_vars_calc_diagnostics_fromIcoGrid( &
       rhog, rhogvx, rhogvy, rhogvz, rhogw, rhoge, rhogq )
!    use scale_atmos_diagnostic_iocA, only: &
!       ATMOS_DIAGNOSTIC_ICOA_get_vel
    use mod_vmtr, only: &
       VMTR_getin_GSGAM2, &
       VMTR_getin_GSGAM2h
    implicit none

    real(RP), intent(in) :: rhog  (ADM_gall_in,ADM_kall,ADM_lall)
    real(RP), intent(in) :: rhogvx(ADM_gall_in,ADM_kall,ADM_lall)
    real(RP), intent(in) :: rhogvy(ADM_gall_in,ADM_kall,ADM_lall)
    real(RP), intent(in) :: rhogvz(ADM_gall_in,ADM_kall,ADM_lall)
    real(RP), intent(in) :: rhogw (ADM_gall_in,ADM_kall,ADM_lall)
    real(RP), intent(in) :: rhoge (ADM_gall_in,ADM_kall,ADM_lall)
    real(RP), intent(in) :: rhogq (ADM_gall_in,ADM_kall,ADM_lall,QA)

    real(RP) :: GSGAM2 (ADM_gall_in,ADM_kall,ADM_lall)
    real(RP) :: GSGAM2h(ADM_gall_in,ADM_kall,ADM_lall)

    integer :: k, i, j, ij, iq, l

    call VMTR_getin_GSGAM2 ( GSGAM2 (:,:,:) )
    call VMTR_getin_GSGAM2H( GSGAM2h(:,:,:) )

    do l = 1, ADM_lall

       !$omp parallel do private(ij)
       do j = 1, JA
       do i = 1, IA
          ij = i + ( j - 1 ) * ADM_imax
          do k = 1, KA
             DENS(k,i,j,l) = rhog (ij,k,l) / GSGAM2(ij,k,l)
             RHOE(k,i,j,l) = rhoge(ij,k,l) / GSGAM2(ij,k,l)
             do iq = 1, QA
                RHOQ(k,i,j,iq,l) = rhogq(ij,k,l,iq) / GSGAM2(ij,k,l)
             end do
          end do
       end do
       end do

    end do

    call ATMOS_DIAGNOSTIC_ICOA_get_vel( &
         KA, 1, KA, IA, 1, IA, JA, 1, JA, ADM_gall_in, ADM_imax, ADM_lall, &
         rhogvx(:,:,:), rhogvy(:,:,:), rhogvz(:,:,:), rhogw(:,:,:), & ! (in)
         GSGAM2(:,:,:), GSGAM2h(:,:,:),                             & ! (in)
         RHOU(:,:,:,:), RHOV(:,:,:,:), MOMZ(:,:,:,:)                ) ! (out)

    return
  end subroutine ATMOS_vars_calc_diagnostics_fromIcoGrid

  !-----------------------------------------------------------------------------
  !> Calc prognostic variables on IcoA grid
  !< inpus:  DENS, RHOU, RHOV, MOMZ, RHOE, RHOQ
  subroutine atmos_vars_calc_prognostics( &
       rhog, rhogvx, rhogvy, rhogvz, rhogw, rhoge, rhogq )
    use mod_vmtr, only: &
       VMTR_getin_GSGAM2, &
       VMTR_getin_GSGAM2h
    implicit none
    real(RP), intent(out) :: rhog  (ADM_gall_in,ADM_kall,ADM_lall)
    real(RP), intent(out) :: rhogvx(ADM_gall_in,ADM_kall,ADM_lall)
    real(RP), intent(out) :: rhogvy(ADM_gall_in,ADM_kall,ADM_lall)
    real(RP), intent(out) :: rhogvz(ADM_gall_in,ADM_kall,ADM_lall)
    real(RP), intent(out) :: rhogw (ADM_gall_in,ADM_kall,ADM_lall)
    real(RP), intent(out) :: rhoge (ADM_gall_in,ADM_kall,ADM_lall)
    real(RP), intent(out) :: rhogq (ADM_gall_in,ADM_kall,ADM_lall,QA)

    real(RP) :: GSGAM2 (ADM_gall_in,ADM_kall,ADM_lall)
    real(RP) :: GSGAM2h(ADM_gall_in,ADM_kall,ADM_lall)

    integer :: k, i, j, iq, ij, l

    call VMTR_getin_GSGAM2 ( GSGAM2 (:,:,:) )
    call VMTR_getin_GSGAM2H( GSGAM2h(:,:,:) )


    call ATMOS_DIAGNOSTIC_ICOA_get_vvel( &
         KA, 1, KA, IA, 1, IA, JA, 1, JA, ADM_gall_in, ADM_imax, ADM_lall, &
         RHOU(:,:,:,:), RHOV(:,:,:,:), MOMZ(:,:,:,:),              & ! [IN]
         GSGAM2(:,:,:), GSGAM2h(:,:,:),                            & ! [IN]
         rhogvx(:,:,:), rhogvy(:,:,:), rhogvz(:,:,:), rhogw(:,:,:) ) ! [OUT]

    !$omp parallel
    do l = 1, ADM_lall

       !$omp do private(ij)
       do k = 1, ADM_kall
       do j = 1, JA
       do i = 1, IA
          ij = i + ( j - 1 ) * ADM_imax
          rhog (ij,k,l) = DENS(k,i,j,l) * GSGAM2(ij,k,l)
          rhoge(ij,k,l) = RHOE(k,i,j,l) * GSGAM2(ij,k,l)
          do iq = 1, QA
             rhogq(ij,k,l,iq) = RHOQ(k,i,j,iq,l) * GSGAM2(ij,k,l)
          end do
       end do
       end do
       end do

    end do
    !$omp end parallel

    return
  end subroutine atmos_vars_calc_prognostics

  !-----------------------------------------------------------------------------
  !> Calc diagnostic variables
  !! inputs:  DENS, RHOU, RHOV, MOMZ, RHOE, RHOQ
  !< outputs: U, V, W, TEMP, POTT, PRES, EXNER, QTRC, Qdry, Rtot, CVtot, CPtot, PHYD, PHYDH
  subroutine ATMOS_vars_calc_diagnostics
    use scale_atmos_thermodyn, only: &
       ATMOS_THERMODYN_specific_heat
    use scale_atmos_diagnostic, only: &
       ATMOS_DIAGNOSTIC_get_therm_rhoe, &
       ATMOS_DIAGNOSTIC_get_phyd
    use mod_atmos_phy_mp_vars, only: &
       ATMOS_PHY_MP_vars_get_diagnostic, &
       ATMOS_PHY_MP_vars_reset_diagnostics
    use mod_atmos_phy_ae_vars, only: &
       ATMOS_PHY_AE_vars_reset_diagnostics
    implicit none

    integer :: k, i, j, iq, l

    do l = 1, ADM_lall

       !$omp parallel do
       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          U(k,i,j,l) = RHOU(k,i,j,l) / DENS(k,i,j,l)
          V(k,i,j,l) = RHOV(k,i,j,l) / DENS(k,i,j,l)
          W(k,i,j,l) = ( MOMZ(k,i,j,l) + MOMZ(k-1,i,j,l) ) * 0.5_RP / DENS(k,i,j,l)
       end do
       end do
       end do


       do iq = 1, QA
          !$omp parallel do
          do j = JS, JE
          do i = IS, IE
          do k = KS, KE
             QTRC(k,i,j,iq,l) = RHOQ(k,i,j,iq,l) / DENS(k,i,j,l)
          end do
          end do
          end do
       end do

       call ATMOS_THERMODYN_specific_heat( &
            KA, KS, KE, IA, 1, IA, JA, 1, JA, QA, &
            QTRC(:,:,:,:,l),                                             & ! (in)
            TRACER_MASS(:), TRACER_R(:), TRACER_CV(:), TRACER_CP(:),     & ! (in)
            Qdry(:,:,:,l), Rtot(:,:,:,l), CVtot(:,:,:,l), CPtot(:,:,:,l) ) ! (out)

       call ATMOS_DIAGNOSTIC_get_therm_rhoe( &
            KA, KS, KE, IA, 1, IA, JA, 1, JA, &
            DENS(:,:,:,l), RHOE(:,:,:,l),                               & ! (in)
            Rtot(:,:,:,l), CPtot(:,:,:,l), CVtot(:,:,:,l),              & ! (in)
            TEMP(:,:,:,l), POTT(:,:,:,l), PRES(:,:,:,l), EXNER(:,:,:,l) ) ! (out)

       call ATMOS_DIAGNOSTIC_get_phyd( &
            KA, KS, KE, IA, 1, IA, JA, 1, JA, &
            DENS(:,:,:,l), PRES(:,:,:,l), & ! (in)
            CZ(:,:,:,l), FZ(:,:,:,l),     & ! (in)
            PHYD(:,:,:,l), PHYDH(:,:,:,l) ) ! (out)
    end do


    call ATMOS_PHY_MP_vars_reset_diagnostics
    call ATMOS_PHY_AE_vars_reset_diagnostics

    if ( moist ) then
       call ATMOS_PHY_MP_vars_get_diagnostic( &
            DENS(:,:,:,:), TEMP(:,:,:,:), QTRC(:,:,:,:,:), & ! (in)
            Qe=Qe(:,:,:,:,:)                               ) ! (out)
    end if

    ! reset diagnostic variables
    DV_calclated(:) = .false.

    return
  end subroutine ATMOS_vars_calc_diagnostics

  ! private

  subroutine ATMOS_DIAGNOSTIC_ICOA_get_vel( &
         KA, KS, KE, IA, IS, IE, JA, JS, JE, GALL, IALL, LALL, &
         rhogvx, rhogvy, rhogvz, rhogw, &
         GSGAM2, GSGAM2h,               &
         RHOU, RHOV, MOMZ               )
    use mod_grd, only: &
       GRD_afac, &
       GRD_bfac
    use mod_gmtr, only: &
       GMTR_p_IX, &
       GMTR_p_IY, &
       GMTR_p_IZ, &
       GMTR_p_JX, &
       GMTR_p_JY, &
       GMTR_p_JZ, &
       GMTR_p
    implicit none
    integer, intent(in) :: KA, KS, KE
    integer, intent(in) :: IA, IS, IE
    integer, intent(in) :: JA, JS, JE
    integer, intent(in) :: GALL, IALL, LALL

    real(RP), intent(in) :: rhogvx(GALL,KA,LALL)
    real(RP), intent(in) :: rhogvy(GALL,KA,LALL)
    real(RP), intent(in) :: rhogvz(GALL,KA,LALL)
    real(RP), intent(in) :: rhogw (GALL,KA,LALL)
    real(RP), intent(in) :: GSGAM2 (GALL,KA,LALL)
    real(RP), intent(in) :: GSGAM2h(GALL,KA,LALL)

    real(RP), intent(out) :: RHOU(KA,IA,JA,LALL)
    real(RP), intent(out) :: RHOV(KA,IA,JA,LALL)
    real(RP), intent(out) :: MOMZ(KA,IA,JA,LALL)

    integer, parameter :: k0 = ADM_KNONE

    real(RP) :: rhovx, rhovy, rhovz
    integer :: k, i, j, l, ij, ij2

    !$omp parallel
    do l = 1, LALL
    !$omp do private(ij,ij2,rhovx,rhovy,rhovz)
    do j = JS, JE
    do i = IS, IE
       ij  = i + (j - 1) * IALL
       ij2 = i + 1 + j * (IALL+1)
       do k = KS, KE
          rhovx = rhogvx(ij,k,l) / GSGAM2(ij,k,l)
          rhovy = rhogvy(ij,k,l) / GSGAM2(ij,k,l)
          rhovz = rhogvz(ij,k,l) / GSGAM2(ij,k,l)
          RHOU(k,i,j,l) = rhovx * GMTR_p(ij2,k0,l,GMTR_p_IX) &
                        + rhovy * GMTR_p(ij2,k0,l,GMTR_p_IY) &
                        + rhovz * GMTR_p(ij2,k0,l,GMTR_p_IZ)
          RHOV(k,i,j,l) = rhovx * GMTR_p(ij2,k0,l,GMTR_p_JX) &
                        + rhovy * GMTR_p(ij2,k0,l,GMTR_p_JY) &
                        + rhovz * GMTR_p(ij2,k0,l,GMTR_p_JZ)
       end do
       do k = KS, KE-1
          MOMZ(k,i,j,l) = rhogw(ij,k+1,l) / GSGAM2h(ij,k+1,l)
       end do
       MOMZ(KE,i,j,l) = 0.0_RP
    end do
    end do
    end do
    !$omp end parallel

    return
  end subroutine ATMOS_DIAGNOSTIC_ICOA_get_vel

  subroutine ATMOS_DIAGNOSTIC_ICOA_get_vvel( &
         KA, KS, KE, IA, IS, IE, JA, JS, JE, GALL, IALL, LALL, &
         RHOU, RHOV, MOMZ, &
         GSGAM2, GSGAM2h, &
         rhogvx, rhogvy, rhogvz, rhogw )
    use mod_gmtr, only: &
       GMTR_p_IX, &
       GMTR_p_IY, &
       GMTR_p_IZ, &
       GMTR_p_JX, &
       GMTR_p_JY, &
       GMTR_p_JZ, &
       GMTR_p
    implicit none
    integer, intent(in) :: KA, KS, KE
    integer, intent(in) :: IA, IS, IE
    integer, intent(in) :: JA, JS, JE
    integer, intent(in) :: GALL, IALL, LALL

    real(RP), intent(in) :: RHOU(KA,IA,JA,LALL)
    real(RP), intent(in) :: RHOV(KA,IA,JA,LALL)
    real(RP), intent(in) :: MOMZ(KA,IA,JA,LALL)
    real(RP), intent(in) :: GSGAM2 (GALL,KA,LALL)
    real(RP), intent(in) :: GSGAM2h(GALL,KA,LALL)

    real(RP), intent(out) :: rhogvx(GALL,KA,LALL)
    real(RP), intent(out) :: rhogvy(GALL,KA,LALL)
    real(RP), intent(out) :: rhogvz(GALL,KA,LALL)
    real(RP), intent(out) :: rhogw (GALL,KA,LALL)

    integer, parameter :: k0 = ADM_KNONE
    integer :: k, i, j, l, ij, ij2

    !$omp parallel
    do l = 1, LALL
       !$omp do private(ij,ij2)
       do k = KS, KE
       do j = JS, JE
       do i = IS, IE
          ij  = i + (j - 1) * IALL
          ij2 = i + 1 + j * (IALL+1)
          rhogvx(ij,k,l) = ( RHOU(k,i,j,l) * GMTR_p(ij2,k0,l,GMTR_p_IX) &
                           + RHOV(k,i,j,l) * GMTR_p(ij2,k0,l,GMTR_p_JX) ) * GSGAM2(ij,k,l)
          rhogvy(ij,k,l) = ( RHOU(k,i,j,l) * GMTR_p(ij2,k0,l,GMTR_p_IY) &
                           + RHOV(k,i,j,l) * GMTR_p(ij2,k0,l,GMTR_p_JY) ) * GSGAM2(ij,k,l)
          rhogvz(ij,k,l) = ( RHOU(k,i,j,l) * GMTR_p(ij2,k0,l,GMTR_p_IZ) &
                           + RHOV(k,i,j,l) * GMTR_p(ij2,k0,l,GMTR_p_JZ) ) * GSGAM2(ij,k,l)
       end do
       end do
       end do
       !$omp do private(ij)
       do k = KS+1, KE
       do j = JS, JE
       do i = IS, IE
          ij = i + (j - 1) * IALL
          rhogw(ij,k,l) = MOMZ(k-1,i,j,l) * GSGAM2h(ij,k,l)
       end do
       end do
       end do
       !$omp do private(ij)
       do j = JS, JE
       do i = IS, IE
          ij = i + (j - 1) * IALL
          rhogw(ij,KS,l) = 0.0_RP
       end do
       end do
    end do
    !$omp end parallel

    return
  end subroutine ATMOS_DIAGNOSTIC_ICOA_get_vvel

  !-----------------------------------------------------------------------------
  !> get diagnostic variable 3D
  recursive subroutine ATMOS_vars_get_diagnostic_3D( &
       vname, &
       var )
    use scale_const, only: &
       GRAV  => CONST_GRAV,  &
       Rvap  => CONST_Rvap,  &
       CPdry => CONST_CPdry, &
       CVdry => CONST_CVdry
    use scale_prc, only: &
       PRC_abort
    use scale_atmos_hydrometeor, only: &
       LHVc => LHV, &
       LHFc => LHF, &
       ATMOS_HYDROMETEOR_LHV, &
       ATMOS_HYDROMETEOR_LHF, &
       ATMOS_HYDROMETEOR_LHS
    use scale_atmos_saturation, only: &
       ATMOS_SATURATION_dens2qsat_all, &
       ATMOS_SATURATION_psat_all, &
       ATMOS_SATURATION_psat_liq, &
       ATMOS_SATURATION_psat_ice, &
       ATMOS_SATURATION_tdew_liq, &
       ATMOS_SATURATION_pote
    use scale_atmos_diagnostic, only: &
       ATMOS_DIAGNOSTIC_get_potv, &
       ATMOS_DIAGNOSTIC_get_teml, &
       ATMOS_DIAGNOSTIC_get_n2
    implicit none
    character(len=*), intent(in)  :: vname
    real(RP),         intent(out) :: var(:,:,:,:)

    real(RP) :: UH  (KA,IA,JA,ADM_lall)
    real(RP) :: VH  (KA,IA,JA,ADM_lall)

    real(RP) :: WORK(KA,IA,JA,ADM_lall)

    integer :: k, i, j, l, iq

    select case ( vname )
    case ( 'W' )
       var(:,:,:,:) = W(:,:,:,:)

    case ( 'U' )
       var(:,:,:,:) = U(:,:,:,:)

    case ( 'V' )
       var(:,:,:,:) = V(:,:,:,:)

    case ( 'PT' )
       var(:,:,:,:) = POTT(:,:,:,:)

    case ( 'T' )
       var(:,:,:,:) = TEMP(:,:,:,:)

    case ( 'EXNER' )
       var(:,:,:,:) = EXNER(:,:,:,:)

    case ( 'PHYD' )
       var(:,:,:,:) = PHYD(:,:,:,:)

    case ( 'QDRY' )
       var(:,:,:,:) = QDRY(:,:,:,:)

    case ( 'RTOT' )
       var(:,:,:,:) = RTOT(:,:,:,:)

    case ( 'CVTOT' )
       var(:,:,:,:) = CVTOT(:,:,:,:)

    case ( 'CPTOT' )
       var(:,:,:,:) = CPTOT(:,:,:,:)

    case ( 'LHV' )
       if ( .not. DV_calclated(I_LHV) ) then
          call allocate_3D( LHV )
          do l = 1, ADM_lall
             call ATMOS_HYDROMETEOR_LHV( &
                  KA, KS, KE, IA, 1, IA, JA, 1, JA, &
                  TEMP(:,:,:,l), & ! (in)
                  LHV(:,:,:,l)   ) ! (out)
          end do
          DV_calclated(I_LHV) = .true.
       end if
       var(:,:,:,:) = LHV(:,:,:,:)

    case ( 'LHS' )
       if ( .not. DV_calclated(I_LHS) ) then
          call allocate_3D( LHS )
          do l = 1, ADM_lall
             call ATMOS_HYDROMETEOR_LHS( &
                  KA, KS, KE, IA, 1, IA, JA, 1, JA, &
                  TEMP(:,:,:,l), & ! (in)
                  LHS(:,:,:,l)   ) ! (out)
          end do
          DV_calclated(I_LHS) = .true.
       end if
       var(:,:,:,:) = LHS(:,:,:,:)

    case ( 'LHF' )
       if ( .not. DV_calclated(I_LHF) ) then
          call allocate_3D( LHF )
          do l = 1, ADM_lall
             call ATMOS_HYDROMETEOR_LHF( &
                  KA, KS, KE, IA, 1, IA, JA, 1, JA, &
                  TEMP(:,:,:,l), & ! (in)
                  LHF(:,:,:,l)   ) ! (out)
          end do
          DV_calclated(I_LHF) = .true.
       end if
       var(:,:,:,:) = LHF(:,:,:,:)

    case ( 'POTV' )
       if ( .not. DV_calclated(I_POTV) ) then
          call allocate_3D( POTV )
          do l = 1, ADM_lall
             call ATMOS_DIAGNOSTIC_get_potv( &
                  KA, KS, KE, IA, 1, IA, JA, 1, JA, &
                  POTT(:,:,:,l), Rtot(:,:,:,l), & ! (in)
                  POTV(:,:,:,l)                 ) ! (out)
          end do
          DV_calclated(I_POTV) = .true.
       end if
       var(:,:,:,:) = POTV(:,:,:,:)

    case ( 'TEML' )
       if ( .not. DV_calclated(I_TEML) ) then
          call allocate_3D( TEML )
          call ATMOS_vars_get_diagnostic( 'LHV', WORK3D(:,:,:,:) )
          call ATMOS_vars_get_diagnostic( 'LHS', WORK3D(:,:,:,:) )
          call ATMOS_vars_get_diagnostic( 'QLIQ', WORK3D(:,:,:,:) )
          call ATMOS_vars_get_diagnostic( 'QICE', WORK3D(:,:,:,:) )
          do l = 1, ADM_lall
             call ATMOS_DIAGNOSTIC_get_teml( &
                  KA, KS, KE, IA, 1, IA, JA, 1, JA, &
                  TEMP(:,:,:,l), LHV(:,:,:,l), LHS(:,:,:,l), & ! (in)
                  QC(:,:,:,l), QI(:,:,:,l), CPtot(:,:,:,l),  & ! (in)
                  TEML(:,:,:,l)                              ) ! (out)
          end do
          DV_calclated(I_TEML) = .true.
       end if
       var(:,:,:,:) = TEML(:,:,:,:)

    case ( 'POTL' )
       if ( .not. DV_calclated(I_POTL) ) then
          call allocate_3D( POTL )
          call ATMOS_vars_get_diagnostic( 'TEML', WORK3D(:,:,:,:) )
!OCL XFILL
          !$omp parallel do default(none) OMP_SCHEDULE_ collapse(2) &
          !$omp private(i,j,k,l) &
          !$omp shared(POTL,TEML,EXNER) &
          !$omp shared(KS,KE,IA,JA,ADM_lall)
          do l = 1, ADM_lall
          do j = 1, JA
          do i = 1, IA
          do k = KS, KE
             POTL(k,i,j,l) = TEML(k,i,j,l) / EXNER(k,i,j,l)
          enddo
          enddo
          enddo
          enddo
          DV_calclated(I_POTL) = .true.
       end if
       var(:,:,:,:) = POTL(:,:,:,:)

    case ( 'POTE' )
       if ( .not. DV_calclated(I_POTE) ) then
          call allocate_3D( POTE )
          do l = 1, ADM_lall
             call ATMOS_SATURATION_pote( &
                  KA, KS, KE, IA, 1, IA, JA, 1, JA, &
                  DENS(:,:,:,l), POTT(:,:,:,l), TEMP(:,:,:,l), QV(:,:,:,l), & ! [IN]
                  POTE(:,:,:,l)                                             ) ! [OUT]
          end do
       end if
       var(:,:,:,:) = POTE(:,:,:,:)
    case ( 'QTOT' )
       if ( .not. DV_calclated(I_QTOT) ) then
          call allocate_3D( QTOT )
          if ( moist ) then
             call ATMOS_vars_get_diagnostic( 'QHYD', WORK3D(:,:,:,:) )
!OCL XFILL
             !$omp parallel do default(none) OMP_SCHEDULE_ collapse(2) &
             !$omp private(i,j,k,l) &
             !$omp shared(QTOT,QV,QHYD) &
             !$omp shared(KS,KE,IA,JA,ADM_lall)
             do l = 1, ADM_lall
             do j = 1, JA
             do i = 1, IA
             do k = KS, KE
                QTOT(k,i,j,l) = QV(k,i,j,l) + QHYD(k,i,j,l)
             enddo
             enddo
             enddo
             enddo
          else
!OCL XFILL
             !$omp parallel do default(none) OMP_SCHEDULE_ collapse(2) &
             !$omp private(i,j,k,l) &
             !$omp shared(QTOT) &
             !$omp shared(KS,KE,IA,JA,ADM_lall)
             do l = 1, ADM_lall
             do j = 1, JA
             do i = 1, IA
             do k = KS, KE
                QTOT(k,i,j,l) = 0.0_RP
             enddo
             enddo
             enddo
             enddo
          end if
          DV_calclated(I_QTOT) = .true.
       end if
       var(:,:,:,:) = QTOT(:,:,:,:)

    case ( 'QHYD' )
       if ( .not. DV_calclated(I_QHYD) ) then
          call allocate_3D( QHYD )
          if ( moist ) then
             call ATMOS_vars_get_diagnostic( 'QLIQ', WORK3D(:,:,:,:) )
             call ATMOS_vars_get_diagnostic( 'QICE', WORK3D(:,:,:,:) )
!OCL XFILL
             !$omp parallel do default(none) OMP_SCHEDULE_ collapse(2) &
             !$omp private(i,j,k,l) &
             !$omp shared(QHYD,QLIQ,QICE) &
             !$omp shared(KS,KE,IA,JA,ADM_lall)
             do l = 1, ADM_lall
             do j = 1, JA
             do i = 1, IA
             do k = KS, KE
                QHYD(k,i,j,l) = QLIQ(k,i,j,l) + QICE(k,i,j,l)
             enddo
             enddo
             enddo
             enddo
          else
!OCL XFILL
             !$omp parallel do default(none) OMP_SCHEDULE_ collapse(2) &
             !$omp private(i,j,k,l) &
             !$omp shared(QHYD) &
             !$omp shared(KS,KE,IA,JA,ADM_lall)
             do l = 1, ADM_lall
             do j = 1, JA
             do i = 1, IA
             do k = KS, KE
                QHYD(k,i,j,l) = 0.0_RP
             enddo
             enddo
             enddo
             enddo
          end if
          DV_calclated(I_QHYD) = .true.
       end if
       var(:,:,:,:) = QHYD(:,:,:,:)

    case ( 'QLIQ' )
       if ( .not. DV_calclated(I_QLIQ) ) then
          call allocate_3D( QLIQ )
          !$omp parallel do default(none) OMP_SCHEDULE_ collapse(2) &
          !$omp private(i,j,k,l) &
          !$omp shared(QLIQ,QC,QR) &
          !$omp shared(KS,KE,IA,JA,ADM_lall)
          do l = 1, ADM_lall
          do j = 1, JA
          do i = 1, IA
          do k = KS, KE
!OCL XFILL
             QLIQ(k,i,j,l) = QC(k,i,j,l) + QR(k,i,j,l)
          enddo
          enddo
          enddo
          enddo
          DV_calclated(I_QLIQ) = .true.
       end if
       var(:,:,:,:) = QLIQ(:,:,:,:)

    case ( 'QICE' )
       if ( .not. DV_calclated(I_QICE) ) then
          call allocate_3D( QICE )
!OCL XFILL
          !$omp parallel do default(none) OMP_SCHEDULE_ collapse(2) &
          !$omp private(i,j,k,l) &
          !$omp shared(QICE,QI,QS,QG,QH) &
          !$omp shared(KS,KE,IA,JA,ADM_lall)
          do l = 1, ADM_lall
          do j = 1, JA
          do i = 1, IA
          do k = KS, KE
             QICE(k,i,j,l) = QI(k,i,j,l) + QS(k,i,j,l) + QG(k,i,j,l) + QH(k,i,j,l)
          enddo
          enddo
          enddo
          enddo
          DV_calclated(I_QICE) = .true.
       end if
       var(:,:,:,:) = QICE(:,:,:,:)

    case ( 'QSAT' )
       if ( .not. DV_calclated(I_QSAT) ) then
          call allocate_3D( QSAT )
          do l = 1, ADM_lall
             call ATMOS_SATURATION_dens2qsat_all( &
                  KA, KS, KE, IA, 1, IA, JA, 1, JA, &
                  TEMP(:,:,:,l), DENS(:,:,:,l), & ! (in)
                  QSAT(:,:,:,l)                 ) ! (out)
          end do
          DV_calclated(I_QSAT) = .true.
       end if
       var(:,:,:,:) = QSAT(:,:,:,:)

    case ( 'RHA' )
       if ( .not. DV_calclated(I_RHA) ) then
          call allocate_3D( RHA )
          if ( moist ) then
             do l = 1, ADM_lall
                call ATMOS_SATURATION_psat_all( &
                     KA, KS, KE, IA, 1, IA, JA, 1, JA, &
                     TEMP(:,:,:,l), & ! (in)
                     WORK(:,:,:,l)  ) ! (out)
             end do
!OCL XFILL
             !$omp parallel do default(none) OMP_SCHEDULE_ collapse(2) &
             !$omp private(i,j,k,l) &
             !$omp shared(RHA,DENS,QV,WORK,TEMP) &
             !$omp shared(KS,KE,IA,JA,ADM_lall)
             do l = 1, ADM_lall
             do j = 1, JA
             do i = 1, IA
             do k = KS, KE
                RHA(k,i,j,l) = DENS(k,i,j,l) * QV(k,i,j,l) &
                             / WORK(k,i,j,l) * Rvap * TEMP(k,i,j,l) &
                             * 100.0_RP
             enddo
             enddo
             enddo
             enddo
          else
!OCL XFILL
             !$omp parallel do default(none) OMP_SCHEDULE_ collapse(2) &
             !$omp private(i,j,k,l) &
             !$omp shared(RHA) &
             !$omp shared(KS,KE,IA,JA,ADM_lall)
             do l = 1, ADM_lall
             do j = 1, JA
             do i = 1, IA
             do k = KS, KE
                RHA(k,i,j,l) = 0.0_RP
             enddo
             enddo
             enddo
             enddo
          end if
          DV_calclated(I_RHA) = .true.
       end if
       var(:,:,:,:) = RHA(:,:,:,:)

    case ( 'RHL', 'RH' )
       if ( .not. DV_calclated(I_RHL) ) then
          call allocate_3D( RHL )
          if ( moist ) then
             do l = 1, ADM_lall
                call ATMOS_SATURATION_psat_liq( &
                     KA, KS, KE, IA, 1, IA, JA, 1, JA, &
                     TEMP(:,:,:,l), & ! (in)
                     WORK(:,:,:,l)  ) ! (out)
             end do
!OCL XFILL
             !$omp parallel do default(none) OMP_SCHEDULE_ collapse(2) &
             !$omp private(i,j,k,l) &
             !$omp shared(RHL,DENS,QV,WORK,TEMP) &
             !$omp shared(KS,KE,IA,JA,ADM_lall)
             do l = 1, ADM_lall
             do j = 1, JA
             do i = 1, IA
             do k = KS, KE
                RHL(k,i,j,l) = DENS(k,i,j,l) * QV(k,i,j,l) &
                             / WORK(k,i,j,l) * Rvap * TEMP(k,i,j,l) &
                             * 100.0_RP
             enddo
             enddo
             enddo
             enddo
          else
!OCL XFILL
             !$omp parallel do default(none) OMP_SCHEDULE_ collapse(2) &
             !$omp private(i,j,k,l) &
             !$omp shared(RHL) &
             !$omp shared(KS,KE,IA,JA,ADM_lall)
             do l = 1, ADM_lall
             do j = 1, JA
             do i = 1, IA
             do k = KS, KE
                RHL(k,i,j,l) = 0.0_RP
             enddo
             enddo
             enddo
             enddo
          end if
          DV_calclated(I_RHL) = .true.
       end if
       var(:,:,:,:) = RHL(:,:,:,:)

    case ( 'RHI' )
       if ( .not. DV_calclated(I_RHI) ) then
          call allocate_3D( RHI )
          if ( moist ) then
             do l = 1, ADM_lall
                call ATMOS_SATURATION_psat_ice( &
                     KA, KS, KE, IA, 1, IA, JA, 1, JA, &
                     TEMP(:,:,:,l), & ! (int)
                     WORK(:,:,:,l)  ) ! (out)
             end do
!OCL XFILL
             !$omp parallel do default(none) OMP_SCHEDULE_ collapse(2) &
             !$omp private(i,j,k,l) &
             !$omp shared(RHI,DENS,QV,WORK,TEMP) &
             !$omp shared(KS,KE,IA,JA,ADM_lall)
             do l = 1, ADM_lall
             do j = 1, JA
             do i = 1, IA
             do k = KS, KE
                RHI(k,i,j,l) = DENS(k,i,j,l) * QV(k,i,j,l) &
                             / WORK(k,i,j,l) * Rvap * TEMP(k,i,j,l) &
                             * 100.0_RP
             enddo
             enddo
             enddo
             enddo
          else
!OCL XFILL
             !$omp parallel do default(none) OMP_SCHEDULE_ collapse(2) &
             !$omp private(i,j,k,l) &
             !$omp shared(RHI) &
             !$omp shared(KS,KE,IA,JA,ADM_lall)
             do l = 1, ADM_lall
             do j = 1, JA
             do i = 1, IA
             do k = KS, KE
                RHI(k,i,j,l) = 0.0_RP
             enddo
             enddo
             enddo
             enddo
          end if
          DV_calclated(I_RHI) = .true.
       end if
       var(:,:,:,:) = RHI(:,:,:,:)

    case ( 'VOR' )
       if ( .not. DV_calclated(I_VOR) ) then
          call allocate_3D( VOR )
          write(*,*) 'xxx not supporte at this moment'
          call PRC_abort
          DV_calclated(I_VOR) = .true.
       end if
       var(:,:,:,:) = VOR(:,:,:,:)

    case ( 'DIV' )
       if ( .not. DV_calclated(I_DIV) ) then
          call allocate_3D( DIV )
          call ATMOS_vars_get_diagnostic( 'HDIV', WORK3D(:,:,:,:) )
          write(*,*) 'xxx not supporte at this moment'
          call PRC_abort
          DV_calclated(I_DIV) = .true.
       end if
       var(:,:,:,:) = DIV(:,:,:,:)

    case ( 'HDIV' )
       if ( .not. DV_calclated(I_HDIV) ) then
          call allocate_3D( HDIV )
          write(*,*) 'xxx not supporte at this moment'
          call PRC_abort
       end if
       var(:,:,:,:) = HDIV(:,:,:,:)

    case ( 'Uabs' )
       if ( .not. DV_calclated(I_Uabs) ) then
          call allocate_3D( Uabs )
!OCL XFILL
          !$omp parallel do private(k,i,j,l) OMP_SCHEDULE_ collapse(2)
          do l = 1, ADM_lall
          do j = 1, JA
          do i = 1, IA
          do k = KS, KE
             Uabs(k,i,j,l) = sqrt( W(k,i,j,l)**2 + U(k,i,j,l)**2 + V(k,i,j,l)**2 )
          enddo
          enddo
          enddo
          enddo
          DV_calclated(I_Uabs) = .true.
       end if
       var(:,:,:,:) = Uabs(:,:,:,:)

    case ( 'N2' )
       if ( .not. DV_calclated(I_N2) ) then
          call allocate_3D( N2 )
          do l = 1, ADM_lall
          call ATMOS_DIAGNOSTIC_get_n2( &
               KA, KS, KE, IA, 1, IA, JA, 1, JA, &
               POTT(:,:,:,l), Rtot(:,:,:,l), & !(in)
               CZ(:,:,:,l),                  & !(in)
               N2(:,:,:,l)                   ) ! (out)
       enddo
          DV_calclated(I_N2) = .true.
       end if
       var(:,:,:,:) = N2(:,:,:,:)

    case ( 'MSE' )
       if ( .not. DV_calclated(I_MSE) ) then
          call allocate_3D( MSE )
          call ATMOS_vars_get_diagnostic( 'LHV', WORK3D(:,:,:,:) )
!OCL XFILL
          !$omp parallel do private(k,i,j,l) OMP_SCHEDULE_ collapse(2)
          do l = 1, ADM_lall
          do j = 1, JA
          do i = 1, IA
          do k = KS, KE
             MSE(k,i,j,l) = CPTOT(k,i,j,l) * TEMP(k,i,j,l)          &
                          + GRAV * ( CZ(k,i,j,l) - FZ(KS-1,i,j,l) ) &
                          + LHV(k,i,j,l) * QV(k,i,j,l)
          enddo
          enddo
          enddo
          enddo
          DV_calclated(I_MSE) = .true.
       end if
       var(:,:,:,:) = MSE(:,:,:,:)

    case ( 'TDEW' )
       if ( .not. DV_calclated(I_TDEW) ) then
          call allocate_3D( TDEW )
          do l = 1, ADM_lall
             call ATMOS_SATURATION_tdew_liq( KA, KS, KE, IA, 1, IA, JA, 1, JA, &
                                             DENS(:,:,:,l), TEMP(:,:,:,l), QV(:,:,:,l), & ! [IN]
                                             TDEW(:,:,:,l)                              ) ! [OUT]
          enddo
          DV_calclated(I_TDEW) = .true.
       end if
       var(:,:,:,:) = TDEW(:,:,:,:)

    case ( 'ENGP' )
       if ( .not. DV_calclated(I_ENGP) ) then
          call allocate_3D( ENGP )
          !$omp parallel do private(k,i,j,l) OMP_SCHEDULE_ collapse(2)
          do l = 1, ADM_lall
          do j = 1, JA
          do i = 1, IA
          do k = KS, KE
             ENGP(k,i,j,l) = DENS(k,i,j,l) * GRAV * CZ(k,i,j,l)
          end do
          end do
          end do
          end do
          DV_calclated(I_ENGP) = .true.
       end if
       var(:,:,:,:) = ENGP(:,:,:,:)

    case ( 'ENGK' )
       if ( .not. DV_calclated(I_ENGK) ) then
          call allocate_3D( ENGK )
          !$omp parallel do private(i,j,k,l) OMP_SCHEDULE_ collapse(2)
          do l = 1, ADM_lall
          do j = 1, JA
          do i = 1, IA
          do k = KS, KE
             ENGK(k,i,j,l) = 0.5_RP * DENS(k,i,j,l) &
                         * ( W(k,i,j,l)**2 + U(k,i,j,l)**2 + V(k,i,j,l)**2 )
          end do
          end do
          end do
          end do
             DV_calclated(I_ENGK) = .true.
       end if
       var(:,:,:,:) = ENGK(:,:,:,:)

    case ( 'ENGI' )
       if ( .not. DV_calclated(I_ENGI) ) then
          call allocate_3D( ENGI )
          if ( moist ) then
             call ATMOS_vars_get_diagnostic( 'QICE', WORK3D(:,:,:,:) )
          end if
          !$omp parallel do private(i,j,k,l) OMP_SCHEDULE_ collapse(2)
          do l = 1, ADM_lall
          do j = 1, JA
          do i = 1, IA
          do k = KS, KE
             ENGI(k,i,j,l) = DENS(k,i,j,l) * QDRY(k,i,j,l) * TEMP(k,i,j,l) * CVdry
             do iq = 1, QA
                ENGI(k,i,j,l) = ENGI(k,i,j,l) &
                              + DENS(k,i,j,l) * QTRC(k,i,j,iq,l) * TEMP(k,i,j,l) * TRACER_CV(iq)
             enddo
             if ( moist ) then
                ENGI(k,i,j,l) = ENGI(k,i,j,l) &
                     + DENS(k,i,j,l) * ( QV  (k,i,j,l) * LHVc & ! Latent Heat [vapor->liquid]
                                       - QICE(k,i,j,l) * LHFc ) ! Latent Heat [ice->liquid]
             end if
          end do
          end do
          end do
          end do
          DV_calclated(I_ENGI) = .true.
       end if
       var(:,:,:,:) = ENGI(:,:,:,:)

    case ( 'ENGT' )
       if ( .not. DV_calclated(I_ENGT) ) then
          call allocate_3D( ENGT )
          call ATMOS_vars_get_diagnostic( 'ENGP', WORK3D(:,:,:,:) )
          call ATMOS_vars_get_diagnostic( 'ENGK', WORK3D(:,:,:,:) )
          call ATMOS_vars_get_diagnostic( 'ENGI', WORK3D(:,:,:,:) )
          !$omp parallel do private(i,j,k,l) OMP_SCHEDULE_ collapse(2)
          do l = 1, ADM_lall
          do j = 1, JA
          do i = 1, IA
          do k = KS, KE
             ENGT(k,i,j,l) = ENGP(k,i,j,l) + ENGK(k,i,j,l) + ENGI(k,i,j,l)
          enddo
          enddo
          enddo
          enddo
          DV_calclated(I_ENGT) = .true.
       end if
       var(:,:,:,:) = ENGT(:,:,:,:)

    case ( 'DENS_PRIM' )
       if ( .not. DV_calclated(I_DENS_PRIM) ) then
          call allocate_3D( DENS_PRIM )
          call ATMOS_vars_get_diagnostic( 'DENS_MEAN', WORK1D(:) )
!OCL XFILL
          !$omp parallel do private(i,j,k,l) OMP_SCHEDULE_ collapse(2)
          do l = 1, ADM_lall
          do j = 1, JA
          do i = 1, IA
          do k = KS, KE
             DENS_PRIM(k,i,j,l) = DENS(k,i,j,l) - DENS_MEAN(k)
          enddo
          enddo
          enddo
          enddo
          DV_calclated(I_DENS_PRIM) = .true.
       end if
       var(:,:,:,:) = DENS_PRIM(:,:,:,:)

    case ( 'W_PRIM' )
       if ( .not. DV_calclated(I_W_PRIM) ) then
          call allocate_3D( W_PRIM )
          call ATMOS_vars_get_diagnostic( 'W_MEAN', WORK1D(:) )
!OCL XFILL
          !$omp parallel do private(i,j,k,l) OMP_SCHEDULE_ collapse(2)
          do l = 1, ADM_lall
          do j = 1, JA
          do i = 1, IA
          do k = KS, KE
             W_PRIM(k,i,j,l) = W(k,i,j,l) - W_MEAN(k)
          enddo
          enddo
          enddo
          enddo
          DV_calclated(I_W_PRIM) = .true.
       end if
       var(:,:,:,:) = W_PRIM(:,:,:,:)

    case ( 'U_PRIM' )
       if ( .not. DV_calclated(I_U_PRIM) ) then
          call allocate_3D( U_PRIM )
          call ATMOS_vars_get_diagnostic( 'U_MEAN', WORK1D(:) )
!OCL XFILL
          !$omp parallel do private(i,j,k,l) OMP_SCHEDULE_ collapse(2)
          do l = 1, ADM_lall
          do j = 1, JA
          do i = 1, IA
          do k = KS, KE
             U_PRIM(k,i,j,l) = U(k,i,j,l) - U_MEAN(k)
          enddo
          enddo
          enddo
          enddo
          DV_calclated(I_U_PRIM) = .true.
       end if
       var(:,:,:,:) = U_PRIM(:,:,:,:)

    case ( 'V_PRIM' )
       if ( .not. DV_calclated(I_V_PRIM) ) then
          call allocate_3D( V_PRIM )
          call ATMOS_vars_get_diagnostic( 'V_MEAN', WORK1D(:) )
!OCL XFILL
          !$omp parallel do private(i,j,k,l) OMP_SCHEDULE_ collapse(2)
          do l = 1, ADM_lall
          do j = 1, JA
          do i = 1, IA
          do k = KS, KE
             V_PRIM(k,i,j,l) = V(k,i,j,l) - V_MEAN(k)
          enddo
          enddo
          enddo
          enddo
          DV_calclated(I_V_PRIM) = .true.
       end if
       var(:,:,:,:) = V_PRIM(:,:,:,:)

    case ( 'PT_PRIM' )
       if ( .not. DV_calclated(I_PT_PRIM) ) then
          call allocate_3D( PT_PRIM )
          call ATMOS_vars_get_diagnostic( 'PT_MEAN', WORK1D(:) )
!OCL XFILL
          !$omp parallel do private(i,j,k,l) OMP_SCHEDULE_ collapse(2)
          do l = 1, ADM_lall
          do j = 1, JA
          do i = 1, IA
          do k = KS, KE
             PT_PRIM(k,i,j,l) = POTT(k,i,j,l) - PT_MEAN(k)
          enddo
          enddo
          enddo
          enddo
          DV_calclated(I_PT_PRIM) = .true.
       end if
       var(:,:,:,:) = PT_PRIM(:,:,:,:)

    case ( 'W_PRIM2' )
       if ( .not. DV_calclated(I_W_PRIM2) ) then
          call allocate_3D( W_PRIM2 )
          call ATMOS_vars_get_diagnostic( 'W_PRIM', WORK3D(:,:,:,:) )
!OCL XFILL
          !$omp parallel do private(i,j,k,l) OMP_SCHEDULE_ collapse(2)
          do l = 1, ADM_lall
          do j = 1, JA
          do i = 1, IA
          do k = KS, KE
             W_PRIM2(k,i,j,l) = W_PRIM(k,i,j,l)**2
          enddo
          enddo
          enddo
          enddo
          DV_calclated(I_W_PRIM2) = .true.
       end if
       var(:,:,:,:) = W_PRIM2(:,:,:,:)

    case ( 'PT_W_PRIM' )
       if ( .not. DV_calclated(I_PT_W_PRIM) ) then
          call allocate_3D( PT_W_PRIM )
          call ATMOS_vars_get_diagnostic( 'W_PRIM',  WORK3D(:,:,:,:) )
          call ATMOS_vars_get_diagnostic( 'PT_PRIM', WORK3D(:,:,:,:) )
          !$omp parallel do private(i,j,k,l) OMP_SCHEDULE_ collapse(2)
          do l = 1, ADM_lall
          do j = 1, JA
          do i = 1, IA
          do k = KS, KE
             PT_W_PRIM(k,i,j,l) = W_PRIM(k,i,j,l) * PT_PRIM(k,i,j,l) * DENS(k,i,j,l) * CPdry
          enddo
          enddo
          enddo
          enddo
          DV_calclated(I_PT_W_PRIM) = .true.
       end if
       var(:,:,:,:) = PT_W_PRIM(:,:,:,:)

    case ( 'W_PRIM3' )
       if ( .not. DV_calclated(I_W_PRIM3) ) then
          call allocate_3D( W_PRIM3 )
          call ATMOS_vars_get_diagnostic( 'W_PRIM', WORK3D(:,:,:,:) )
!OCL XFILL
          !$omp parallel do private(i,j,k,l) OMP_SCHEDULE_ collapse(2)
          do l = 1, ADM_lall
          do j = 1, JA
          do i = 1, IA
          do k = KS, KE
             W_PRIM3(k,i,j,l) = W_PRIM(k,i,j,l)**3
          enddo
          enddo
          enddo
          enddo
          DV_calclated(I_W_PRIM3) = .true.
       end if
       var(:,:,:,:) = W_PRIM3(:,:,:,:)

    case ( 'TKE_RS' )
       if ( .not. DV_calclated(I_TKE_RS) ) then
          call allocate_3D( TKE_RS )
          call ATMOS_vars_get_diagnostic( 'W_PRIM', WORK3D(:,:,:,:) )
          call ATMOS_vars_get_diagnostic( 'U_PRIM', WORK3D(:,:,:,:) )
          call ATMOS_vars_get_diagnostic( 'V_PRIM', WORK3D(:,:,:,:) )
!OCL XFILL
          !$omp parallel do private(i,j,k,l) OMP_SCHEDULE_ collapse(2)
          do l = 1, ADM_lall
          do j = 1, JA
          do i = 1, IA
          do k = KS, KE
             TKE_RS(k,i,j,l) = 0.5_RP * ( W_PRIM(k,i,j,l)**2 + U_PRIM(k,i,j,l)**2 + V_PRIM(k,i,j,l)**2 )
          enddo
          enddo
          enddo
          enddo
          DV_calclated(I_TKE_RS) = .true.
       end if
       var(:,:,:,:) = TKE_RS(:,:,:,:)

    case default
       LOG_ERROR("ATMOS_vars_calc_diagnostics",*) 'name is invalid for ATMOS_vars_get_diagnostic_3D: ', trim(vname)
       call PRC_abort
    end select


    return
  end subroutine ATMOS_vars_get_diagnostic_3D

  !-----------------------------------------------------------------------------
  !> get diagnostic variable 2D
  recursive subroutine ATMOS_vars_get_diagnostic_2D( &
       vname, &
       var )
    use scale_prc, only: &
       PRC_abort
    use scale_atmos_adiabat, only: &
       ATMOS_ADIABAT_cape
    use mod_atmos_phy_mp_vars, only: &
       SFLX_rain_MP => ATMOS_PHY_MP_SFLX_rain, &
       SFLX_snow_MP => ATMOS_PHY_MP_SFLX_snow
    implicit none
    character(len=*), intent(in)  :: vname
    real(RP),         intent(out) :: var(:,:,:)

    real(RP) :: fact

    integer :: k, i, j, l, iq

    select case ( vname )
    case ( 'LWP' )
       if ( .not. DV_calclated(I_LWP) ) then
          call allocate_2D( LWP )
          call ATMOS_vars_get_diagnostic( 'QLIQ', WORK3D(:,:,:,:) )
          !$omp parallel do default(none) OMP_SCHEDULE_ collapse(2) &
          !$omp private(i,j,k,l) &
          !$omp shared(LWP,QLIQ,DENS,FZ) &
          !$omp shared(KS,KE,IA,JA,ADM_lall)
          do l = 1, ADM_lall
          do j = 1, JA
          do i = 1, IA
             LWP(i,j,l) = 0.0_RP
             do k  = KS, KE
                LWP(i,j,l) = LWP(i,j,l) &
                           + QLIQ(k,i,j,l) * DENS(k,i,j,l) * ( FZ(k,i,j,l)-FZ(k-1,i,j,l) ) * 1.E3_RP ! [kg/m2->g/m2]
             enddo
          enddo
          enddo
          enddo
          DV_calclated(I_LWP) = .true.
       end if
       var(:,:,:) = LWP(:,:,:)

    case ( 'IWP' )
       if ( .not. DV_calclated(I_IWP) ) then
          call allocate_2D( IWP )
          call ATMOS_vars_get_diagnostic( 'QICE', WORK3D(:,:,:,:) )
          !$omp parallel do default(none) OMP_SCHEDULE_ collapse(2) &
          !$omp private(i,j,k,l) &
          !$omp shared(IWP,QICE,DENS,FZ) &
          !$omp shared(KS,KE,IA,JA,ADM_lall)
          do l = 1, ADM_lall
          do j = 1, JA
          do i = 1, IA
             IWP(i,j,l) = 0.0_RP
             do k  = KS, KE
                IWP(i,j,l) = IWP(i,j,l) &
                           + QICE(k,i,j,l) * DENS(k,i,j,l) * ( FZ(k,i,j,l)-FZ(k-1,i,j,l) ) * 1.E3_RP ! [kg/m2->g/m2]
             enddo
          enddo
          enddo
          enddo
          DV_calclated(I_IWP) = .true.
       end if
       var(:,:,:) = IWP(:,:,:)

    case ( 'PW' )
       if ( .not. DV_calclated(I_PW) ) then
          call allocate_2D( PW )
          !$omp parallel do default(none) OMP_SCHEDULE_ collapse(2) &
          !$omp private(i,j,k,l) &
          !$omp shared(PW,QV,DENS,FZ) &
          !$omp shared(KS,KE,IA,JA,ADM_lall)
          do l = 1, ADM_lall
          do j = 1, JA
          do i = 1, IA
             PW(i,j,l) = 0.0_RP
             do k  = KS, KE
                PW(i,j,l) = PW(i,j,l) &
                          + QV(k,i,j,l) * DENS(k,i,j,l) * ( FZ(k,i,j,l)-FZ(k-1,i,j,l) ) * 1.E3_RP ! [kg/m2->g/m2]
             enddo
          enddo
          enddo
          enddo
          DV_calclated(I_PW) = .true.
       end if
       var(:,:,:) = PW(:,:,:)

    case ( 'PBLH' )
       if ( .not. DV_calclated(I_PBLH) ) then
          call allocate_2D( PBLH )
          call ATMOS_vars_get_diagnostic( 'POTV', WORK3D(:,:,:,:) )
          !$omp parallel do default(none) OMP_SCHEDULE_ collapse(2) &
          !$omp private(k,i,j,l) &
          !$omp private(fact) &
          !$omp shared(PBLH,POTV,CZ,FZ) &
          !$omp shared(KS,KE,IA,JA,ADM_lall)
          do l = 1, ADM_lall
          do j = 1, JA
          do i = 1, IA
             PBLH(i,j,l) = CZ(KS,i,j,l) - FZ(KS-1,i,j,l)
             do k = KS+1, KE
                if ( POTV(k,i,j,l) > POTV(KS,i,j,l) ) then
                   fact = ( POTV(KS,i,j,l) - POTV(k-1,i,j,l) ) &
                        / ( POTV(k,i,j,l)  - POTV(k-1,i,j,l) )
                   PBLH(i,j,l) = CZ(k-1,i,j,l) - FZ(KS-1,i,j,l) &
                               + fact * ( CZ(k,i,j,l) - CZ(k-1,i,j,l) )

                   exit
                endif
             enddo
          enddo
          enddo
          enddo
          DV_calclated(I_PBLH) = .true.
       end if
       var(:,:,:) = PBLH(:,:,:)

    case ( 'CAPE', 'CIN', 'LCL', 'LFC', 'LNB' )
       if ( .not. DV_calclated(I_CAPE) ) then
          call allocate_2D( CAPE )
          call allocate_2D( CIN )
          call allocate_2D( LCL )
          call allocate_2D( LFC )
          call allocate_2D( LNB )
          do l = 1, ADM_lall
             call ATMOS_ADIABAT_cape( &
                  KA, KS, KE, IA, IS, IE, JA, JS, JE, &
                  KS,                                                         & ! (in)
                  DENS(:,:,:,l), TEMP(:,:,:,l), PRES(:,:,:,l),                & ! (in)
                  QV(:,:,:,l), QC(:,:,:,l), Qdry(:,:,:,l),                    & ! (in)
                  Rtot(:,:,:,l), CPtot(:,:,:,l),                              & ! (in)
                  CZ(:,:,:,l), FZ(:,:,:,l),                                   & ! (in)
                  CAPE(:,:,l), CIN(:,:,l), LCL(:,:,l), LFC(:,:,l), LNB(:,:,l) ) ! (out)
          end do
          DV_calclated(I_CAPE) = .true.
       end if
       select case ( vname )
       case ( 'CAPE' )
          var(:,:,:) = CAPE(:,:,:)
       case ( 'CIN' )
          var(:,:,:) = CIN(:,:,:)
       case ( 'LCL' )
          var(:,:,:) = LCL(:,:,:)
       case ( 'LFC' )
          var(:,:,:) = LFC(:,:,:)
       case ( 'LNB' )
          var(:,:,:) = LNB(:,:,:)
       end select

    case ( 'PREC', 'RAIN', 'SNOW' )
       if ( .not. DV_calclated(I_PREC) ) then
          call allocate_2D( PREC )
          call allocate_2D( RAIN )
          call allocate_2D( SNOW )
          !$omp parallel do private(i,j,l) OMP_SCHEDULE_
          do l = 1, ADM_lall
          do j = 1, JA
          do i = 1, IA
             RAIN(i,j,l) = SFLX_rain_MP(i,j,l)
             SNOW(i,j,l) = SFLX_snow_MP(i,j,l)
             PREC(i,j,l) = RAIN(i,j,l) + SNOW(i,j,l)
          enddo
          enddo
          enddo
          DV_calclated(I_PREC) = .true.
       end if
       select case (vname)
       case ( 'RAIN' )
          var(:,:,:) = RAIN(:,:,:)
       case ( 'SNOW' )
          var(:,:,:) = SNOW(:,:,:)
       case ( 'PREC' )
          var(:,:,:) = PREC(:,:,:)
       end select

    case default
       LOG_ERROR("ATMOS_vars_calc_diagnostics",*) 'name is invalid for ATMOS_vars_get_diagnostic_2D: ', trim(vname)
       call PRC_abort
    end select


    return
  end subroutine ATMOS_vars_get_diagnostic_2D

  !-----------------------------------------------------------------------------
  !> get diagnostic variable 1D
  recursive subroutine ATMOS_vars_get_diagnostic_1D( &
       vname, &
       var )
    use scale_const, only: &
       CPdry => CONST_CPdry
    use scale_prc, only: &
       PRC_abort
!    use scale_comm_icoA, only: &
!       COMM_horizontal_mean
    implicit none
    character(len=*), intent(in)  :: vname
    real(RP),         intent(out) :: var(:)

    real(RP) :: WORK(KA,IA,JA,ADM_lall)

    integer :: k, i, j, l, iq

    select case ( vname )
    case ( 'DENS_MEAN' )
       if ( .not. DV_calclated(I_DENS_MEAN) ) then
          call allocate_1D( DENS_MEAN )
          write(*,*) 'xxx not supported yet'
          call PRC_abort
!          call COMM_horizontal_mean( DENS_MEAN(:), DENS(:,:,:,:) )
          DV_calclated(I_DENS_MEAN) = .true.
       end if
       var(:) = DENS_MEAN(:)

    case ( 'W_MEAN' )
       if ( .not. DV_calclated(I_W_MEAN) ) then
          write(*,*) 'xxx not supported yet'
          call PRC_abort
          call allocate_1D( W_MEAN )
          call ATMOS_vars_get_diagnostic( 'DENS_MEAN', WORK1D(:) )
!OCL XFILL
          !$omp parallel do private(i,j,k,l) OMP_SCHEDULE_ collapse(2)
          do l = 1, ADM_lall
          do j = JS, JE
          do i = IS, IE
          do k = KS, KE
             WORK(k,i,j,l) = W(k,i,j,l) * DENS(k,i,j,l)
          enddo
          enddo
          enddo
          enddo
!          call COMM_horizontal_mean( W_MEAN(:), WORK(:,:,:,:) )
          do k = KS, KE
             W_MEAN(k) = W_MEAN(k) / DENS_MEAN(k)
          enddo
          DV_calclated(I_W_MEAN) = .true.
       end if
       var(:) = W_MEAN(:)

    case ( 'U_MEAN' )
       if ( .not. DV_calclated(I_U_MEAN) ) then
          write(*,*) 'xxx not supported yet'
          call PRC_abort
          call allocate_1D( U_MEAN )
          call ATMOS_vars_get_diagnostic( 'DENS_MEAN', WORK1D(:) )
!OCL XFILL
          !$omp parallel do private(i,j,k,l) OMP_SCHEDULE_ collapse(2)
          do l = 1, ADM_lall
          do j = JS, JE
          do i = IS, IE
          do k = KS, KE
             WORK(k,i,j,l) = U(k,i,j,l) * DENS(k,i,j,l)
          enddo
          enddo
          enddo
          enddo
!          call COMM_horizontal_mean( U_MEAN(:), WORK(:,:,:,:) )
          do k = KS, KE
             U_MEAN(k) = U_MEAN(k) / DENS_MEAN(k)
          enddo
          DV_calclated(I_U_MEAN) = .true.
       end if
       var(:) = U_MEAN(:)

    case ( 'V_MEAN' )
       if ( .not. DV_calclated(I_V_MEAN) ) then
          write(*,*) 'xxx not supported yet'
          call PRC_abort
          call allocate_1D( V_MEAN )
          call ATMOS_vars_get_diagnostic( 'DENS_MEAN', WORK1D(:) )
!OCL XFILL
          !$omp parallel do private(i,j,k,l) OMP_SCHEDULE_ collapse(2)
          do l = 1, ADM_lall
          do j = JS, JE
          do i = IS, IE
          do k = KS, KE
             WORK(k,i,j,l) = V(k,i,j,l) * DENS(k,i,j,l)
          enddo
          enddo
          enddo
          enddo
!          call COMM_horizontal_mean( V_MEAN(:), WORK(:,:,:,:) )
          do k = KS, KE
             V_MEAN(k) = V_MEAN(k) / DENS_MEAN(k)
          enddo
          DV_calclated(I_V_MEAN) = .true.
       end if
       var(:) = V_MEAN(:)

    case ( 'PT_MEAN' )
       if ( .not. DV_calclated(I_PT_MEAN) ) then
          write(*,*) 'xxx not supported yet'
          call PRC_abort
          call allocate_1D( PT_MEAN )
          call ATMOS_vars_get_diagnostic( 'DENS_MEAN', WORK1D(:) )
!          call COMM_horizontal_mean( PT_MEAN(:), RHOT(:,:,:,:) )
          do k = KS, KE
             PT_MEAN(k) = PT_MEAN(k) / DENS_MEAN(k)
          enddo
          DV_calclated(I_PT_MEAN) = .true.
       end if
       var(:) = PT_MEAN(:)

    case ( 'T_MEAN' )
       if ( .not. DV_calclated(I_T_MEAN) ) then
          write(*,*) 'xxx not supported yet'
          call PRC_abort
          call allocate_1D( T_MEAN )
          call ATMOS_vars_get_diagnostic( 'DENS_MEAN', WORK1D(:) )
!OCL XFILL
          !$omp parallel do private(i,j,k,l) OMP_SCHEDULE_ collapse(2)
          do l = 1, ADM_lall
          do j = JS, JE
          do i = IS, IE
          do k = KS, KE
             WORK(k,i,j,l) = TEMP(k,i,j,l) * DENS(k,i,j,l)
          enddo
          enddo
          enddo
          enddo
!          call COMM_horizontal_mean( T_MEAN(:), WORK(:,:,:,:) )
          do k = KS, KE
             T_MEAN(k) = T_MEAN(k) / DENS_MEAN(k)
          enddo
          DV_calclated(I_T_MEAN) = .true.
       end if
       var(:) = T_MEAN(:)

    case ( 'QV_MEAN' )
       if ( .not. DV_calclated(I_QV_MEAN) ) then
          write(*,*) 'xxx not supported yet'
          call PRC_abort
          call allocate_1D( QV_MEAN )
          if ( moist ) then
             call ATMOS_vars_get_diagnostic( 'DENS_MEAN', WORK1D(:) )
!OCL XFILL
             !$omp parallel do private(i,j,k,l) OMP_SCHEDULE_ collapse(2)
             do l = 1, ADM_lall
             do j = JS, JE
             do i = IS, IE
             do k = KS, KE
                WORK(k,i,j,l) = QV(k,i,j,l) * DENS(k,i,j,l)
             enddo
             enddo
             enddo
             enddo
!             call COMM_horizontal_mean( QV_MEAN(:), WORK(:,:,:,:) )
             do k = KS, KE
                QV_MEAN(k) = QV_MEAN(k) / DENS_MEAN(k)
             enddo
          else
             !$omp parallel do private(k) OMP_SCHEDULE_
             do k = KS, KE
                QV_MEAN(k) = 0.0_RP
             enddo
          end if
          DV_calclated(I_QV_MEAN) = .true.
       end if
       var(:) = QV_MEAN(:)

    case ( 'QHYD_MEAN' )
       if ( .not. DV_calclated(I_QHYD_MEAN) ) then
          write(*,*) 'xxx not supported yet'
          call PRC_abort
          call allocate_1D( QHYD_MEAN )
          call ATMOS_vars_get_diagnostic( 'DENS_MEAN', WORK1D(:) )
          call ATMOS_vars_get_diagnostic( 'QHYD', WORK3D(:,:,:,:) )
!OCL XFILL
          !$omp parallel do private(i,j,k,l) OMP_SCHEDULE_ collapse(2)
          do l = 1, ADM_lall
          do j = JS, JE
          do i = IS, IE
          do k = KS, KE
             WORK(k,i,j,l) = QHYD(k,i,j,l) * DENS(k,i,j,l)
          enddo
          enddo
          enddo
          enddo
!          call COMM_horizontal_mean( QHYD_MEAN(:), WORK(:,:,:,:) )
          do k = KS, KE
             QHYD_MEAN(k) = QHYD_MEAN(k) / DENS_MEAN(k)
          enddo
          DV_calclated(I_QHYD_MEAN) = .true.
       end if
       var(:) = QHYD_MEAN(:)

    case ( 'QLIQ_MEAN' )
       if ( .not. DV_calclated(I_QLIQ_MEAN) ) then
          write(*,*) 'xxx not supported yet'
          call PRC_abort
          call allocate_1D( QLIQ_MEAN )
          call ATMOS_vars_get_diagnostic( 'DENS_MEAN', WORK1D(:) )
          call ATMOS_vars_get_diagnostic( 'QLIQ', WORK3D(:,:,:,:) )
!OCL XFILL
          !$omp parallel do private(i,j,k,l) OMP_SCHEDULE_ collapse(2)
          do l = 1, ADM_lall
          do j = JS, JE
          do i = IS, IE
          do k = KS, KE
             WORK(k,i,j,l) = QLIQ(k,i,j,l) * DENS(k,i,j,l)
          enddo
          enddo
          enddo
          enddo
!          call COMM_horizontal_mean( QLIQ_MEAN(:), WORK(:,:,:,:) )
          do k = KS, KE
             QLIQ_MEAN(k) = QLIQ_MEAN(k) / DENS_MEAN(k)
          enddo
          DV_calclated(I_QLIQ_MEAN) = .true.
       end if
       var(:) = QLIQ_MEAN(:)

    case ( 'QICE_MEAN' )
       if ( .not. DV_calclated(I_QICE_MEAN) ) then
          write(*,*) 'xxx not supported yet'
          call PRC_abort
          call allocate_1D( QICE_MEAN )
          call ATMOS_vars_get_diagnostic( 'DENS_MEAN', WORK1D(:) )
          call ATMOS_vars_get_diagnostic( 'QICE', WORK3D(:,:,:,:) )
!OCL XFILL
          !$omp parallel do private(i,j,k,l) OMP_SCHEDULE_ collapse(2)
          do l = 1, ADM_lall
          do j = JS, JE
          do i = IS, IE
          do k = KS, KE
             WORK(k,i,j,l) = QICE(k,i,j,l) * DENS(k,i,j,l)
          enddo
          enddo
          enddo
          enddo
!          call COMM_horizontal_mean( QICE_MEAN(:), WORK(:,:,:,:) )
          do k = KS, KE
             QICE_MEAN(k) = QICE_MEAN(k) / DENS_MEAN(k)
          enddo
          DV_calclated(I_QICE_MEAN) = .true.
       end if
       var(:) = QICE_MEAN(:)

    case default
       LOG_ERROR("ATMOS_vars_calc_diagnostics",*) 'name is invalid for ATMOS_vars_get_diagnostic_1D: ', trim(vname)
       call PRC_abort
    end select


    return
  end subroutine ATMOS_vars_get_diagnostic_1D

  !-----------------------------------------------------------------------------
  !> monitor output
  subroutine ATMOS_vars_monitor
    use scale_prc, only: &
       PRC_myrank, &
       PRC_abort
    use scale_const, only: &
       GRAV  => CONST_GRAV,  &
       CVdry => CONST_CVdry
    use scale_atmos_grid_cartesC, only: &
       RFDX => ATMOS_GRID_CARTESC_RFDX, &
       RFDY => ATMOS_GRID_CARTESC_RFDY
    use scale_atmos_grid_cartesC_metric, only: &
       MAPF => ATMOS_GRID_CARTESC_METRIC_MAPF
    use scale_statistics, only: &
       STATISTICS_checktotal, &
       STATISTICS_total,            &
       STATISTICS_detail
    use scale_monitor, only: &
       MONITOR_put
    use scale_time, only: &
       TIME_DTSEC_ATMOS_DYN
    use mod_atmos_admin, only: &
       ATMOS_DYN_TYPE
    use scale_atmos_hydrometeor, only: &
       I_QV
    use mod_atmos_phy_mp_vars, only: &
       SFLX_rain_MP => ATMOS_PHY_MP_SFLX_rain, &
       SFLX_snow_MP => ATMOS_PHY_MP_SFLX_snow
    use mod_atmos_phy_rd_vars, only: &
       SFLX_LW_up   => ATMOS_PHY_RD_SFLX_LW_up,   &
       SFLX_LW_dn   => ATMOS_PHY_RD_SFLX_LW_dn,   &
       SFLX_SW_up   => ATMOS_PHY_RD_SFLX_SW_up,   &
       SFLX_SW_dn   => ATMOS_PHY_RD_SFLX_SW_dn,   &
       TOAFLX_LW_up => ATMOS_PHY_RD_TOAFLX_LW_up, &
       TOAFLX_LW_dn => ATMOS_PHY_RD_TOAFLX_LW_dn, &
       TOAFLX_SW_up => ATMOS_PHY_RD_TOAFLX_SW_up, &
       TOAFLX_SW_dn => ATMOS_PHY_RD_TOAFLX_SW_dn
    use mod_atmos_phy_sf_vars, only: &
       SFLX_SH   => ATMOS_PHY_SF_SFLX_SH, &
       SFLX_LH   => ATMOS_PHY_SF_SFLX_LH, &
       SFLX_QTRC => ATMOS_PHY_SF_SFLX_QTRC
    implicit none

    real(RP) :: RHOQ(KA,IA,JA,ADM_lall)

    real(RP) :: ENGFLXT    (IA,JA,ADM_lall) ! total flux             [J/m2/s]
    real(RP) :: SFLX_RD_net(IA,JA,ADM_lall) ! net SFC radiation flux [J/m2/s]
    real(RP) :: TFLX_RD_net(IA,JA,ADM_lall) ! net TOA radiation flux [J/m2/s]

    real(RP)               :: WORK (KA,IA,JA,ADM_lall,3)
    character(len=H_SHORT) :: WNAME(3)
    real(RP)               :: CFLMAX

    integer  :: k, i, j, l, iq
    !---------------------------------------------------------------------------

    write(*,*) 'not supported yet'
    call PRC_abort

    !##### Mass Budget #####

    do iq = 1, QA
       !$omp parallel do private(i,j,k,l) OMP_SCHEDULE_ collapse(2)
!OCL XFILL
       do l = 1, ADM_lall
       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          RHOQ(k,i,j,l) = DENS(k,i,j,l) * QTRC(k,i,j,iq,l)
       enddo
       enddo
       enddo
       enddo

!       call MONITOR_put( QP_MONIT_id(iq), RHOQ(:,:,:,:) )
    enddo

    ! total dry airmass
    if ( DV_MONIT_id(IM_QDRY) > 0 ) then
       !$omp parallel do private(i,j,k,l) OMP_SCHEDULE_ collapse(2)
!OCL XFILL
       do l = 1, ADM_lall
       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          RHOQ(k,i,j,l) = DENS(k,i,j,l) * QDRY (k,i,j,l)
       enddo
       enddo
       enddo
       enddo
!       call MONITOR_put( DV_MONIT_id(IM_QDRY), RHOQ(:,:,:,:) )
    end if

    ! total vapor,liquid,solid tracers
    if ( DV_MONIT_id(IM_QTOT) > 0 ) then
       call ATMOS_vars_get_diagnostic( 'QTOT', WORK3D(:,:,:,:) )
       !$omp parallel do private(i,j,k,l) OMP_SCHEDULE_ collapse(2)
!OCL XFILL
       do l = 1, ADM_lall
       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          RHOQ(k,i,j,l) = DENS(k,i,j,l) * QTOT(k,i,j,l)
       enddo
       enddo
       enddo
       enddo
!       call MONITOR_put( DV_MONIT_id(IM_QTOT), RHOQ(:,:,:,:) )
    end if

    ! total evapolation
!    if ( moist ) call MONITOR_put( DV_MONIT_id(IM_EVAP), SFLX_QTRC(:,:,:,I_QV) )

    ! total precipitation
    if ( DV_MONIT_id(IM_PREC) > 0 ) then
       call ATMOS_vars_get_diagnostic( 'PREC', WORK2D(:,:,:) )
!       call MONITOR_put( DV_MONIT_id(IM_PREC), WORK2D(:,:,:) )
    end if


    !##### Energy Budget #####

    if ( DV_MONIT_id(IM_ENGT) > 0 ) then
       call ATMOS_vars_get_diagnostic( 'ENGT', WORK3D(:,:,:,:) )
!       call MONITOR_put( DV_MONIT_id(IM_ENGT), WORK3D(:,:,:,:) )
    end if
    if ( DV_MONIT_id(IM_ENGP) > 0 ) then
       call ATMOS_vars_get_diagnostic( 'ENGP', WORK3D(:,:,:,:) )
!       call MONITOR_put( DV_MONIT_id(IM_ENGP), WORK3D(:,:,:,:) )
    end if
    if ( DV_MONIT_id(IM_ENGK) > 0 ) then
       call ATMOS_vars_get_diagnostic( 'ENGK', WORK3D(:,:,:,:) )
!       call MONITOR_put( DV_MONIT_id(IM_ENGK), WORK3D(:,:,:,:) )
    end if
    if ( DV_MONIT_id(IM_ENGI) > 0 ) then
       call ATMOS_vars_get_diagnostic( 'ENGI', WORK3D(:,:,:,:) )
!       call MONITOR_put( DV_MONIT_id(IM_ENGI), WORK3D(:,:,:,:) )
    end if


    ! radiation flux
!OCL XFILL
    !$omp parallel do private(i,j,l) OMP_SCHEDULE_ collapse(2)
    do l = 1, ADM_lall
    do j = JS, JE
    do i = IS, IE
       SFLX_RD_net(i,j,l) = ( SFLX_LW_up(i,j,l) - SFLX_LW_dn(i,j,l) ) &
                          + ( SFLX_SW_up(i,j,l) - SFLX_SW_dn(i,j,l) )

       TFLX_RD_net(i,j,l) = ( TOAFLX_LW_up(i,j,l) - TOAFLX_LW_dn(i,j,l) ) &
                          + ( TOAFLX_SW_up(i,j,l) - TOAFLX_SW_dn(i,j,l) )

       ENGFLXT    (i,j,l) = SFLX_SH(i,j,l) + SFLX_LH(i,j,l) &
                          + SFLX_RD_net(i,j,l) - TFLX_RD_net(i,j,l)
    enddo
    enddo
    enddo

!!$    call MONITOR_put( DV_MONIT_id(IM_ENGFLXT),      ENGFLXT     (:,:,:) )
!!$
!!$    call MONITOR_put( DV_MONIT_id(IM_ENGSFC_SH),    SFLX_SH     (:,:,:) )
!!$    call MONITOR_put( DV_MONIT_id(IM_ENGSFC_LH),    SFLX_LH     (:,:,:) )
!!$    call MONITOR_put( DV_MONIT_id(IM_ENGSFC_RD),    SFLX_RD_net (:,:,:) )
!!$    call MONITOR_put( DV_MONIT_id(IM_ENGTOA_RD),    TFLX_RD_net (:,:,:) )
!!$
!!$    call MONITOR_put( DV_MONIT_id(IM_ENGSFC_LW_up), SFLX_LW_up  (:,:,:) )
!!$    call MONITOR_put( DV_MONIT_id(IM_ENGSFC_LW_dn), SFLX_LW_dn  (:,:,:) )
!!$    call MONITOR_put( DV_MONIT_id(IM_ENGSFC_SW_up), SFLX_SW_up  (:,:,:) )
!!$    call MONITOR_put( DV_MONIT_id(IM_ENGSFC_SW_dn), SFLX_SW_dn  (:,:,:) )
!!$
!!$    call MONITOR_put( DV_MONIT_id(IM_ENGTOA_LW_up), TOAFLX_LW_up(:,:,:) )
!!$    call MONITOR_put( DV_MONIT_id(IM_ENGTOA_LW_dn), TOAFLX_LW_dn(:,:,:) )
!!$    call MONITOR_put( DV_MONIT_id(IM_ENGTOA_SW_up), TOAFLX_SW_up(:,:,:) )
!!$    call MONITOR_put( DV_MONIT_id(IM_ENGTOA_SW_dn), TOAFLX_SW_dn(:,:,:) )



    if ( ATMOS_VARS_CHECKRANGE ) then
!OCL XFILL
       WORK(:,:,:,:,1) = W(:,:,:,:)
!OCL XFILL
       WORK(:,:,:,:,2) = U(:,:,:,:)
!OCL XFILL
       WORK(:,:,:,:,3) = V(:,:,:,:)

       WNAME(1) = "W"
       WNAME(2) = "U"
       WNAME(3) = "V"

!!$       call STATISTICS_detail( KA, KS, KE, IA, IS, IE, JA, JS, JE, 3, &
!!$                               WNAME(:), WORK(:,:,:,:,:)                )
    endif

    if (       ( ATMOS_DYN_TYPE /= 'OFF' .AND. ATMOS_DYN_TYPE /= 'NONE' )                   &
         .AND. ( ATMOS_VARS_CHECKCFL_SOFT > 0.0_RP .OR. ATMOS_VARS_CHECKCFL_HARD > 0.0_RP ) ) then
!OCL XFILL
       WORK(:,:,:,:,:) = 0.0_RP

       do l = 1, ADM_lall
       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          WORK(k,i,j,l,1) = abs(W(k,i,j,l)) &
                        * TIME_DTSEC_ATMOS_DYN / ( CZ(k+1,i,j,l) - CZ(k,i,j,l) )
          WORK(k,i,j,l,2) = abs(U(k,i,j,l)) &
                        * TIME_DTSEC_ATMOS_DYN * RFDX(i) !* MAPF(i,j,1,I_UY)
          WORK(k,i,j,l,3) = abs(V(k,i,j,l)) &
                        * TIME_DTSEC_ATMOS_DYN * RFDY(j) !* MAPF(i,j,2,I_XV)
       end do
       enddo
       enddo
       enddo

       CFLMAX = maxval( WORK(:,:,:,:,:) )

       if ( ATMOS_VARS_CHECKCFL_HARD > 0.0_RP .AND. CFLMAX > ATMOS_VARS_CHECKCFL_HARD ) then
          LOG_INFO("ATMOS_vars_monitor",*) "Courant number =", CFLMAX, " exceeded the hard limit =", ATMOS_VARS_CHECKCFL_HARD
          LOG_ERROR("ATMOS_vars_monitor",*)                     "Courant number =", CFLMAX, " exceeded the hard limit =", ATMOS_VARS_CHECKCFL_HARD
          LOG_ERROR_CONT(*)                     "Rank =", PRC_myrank
          LOG_ERROR_CONT(*)                     "Please set ATMOS_VARS_CHECKCFL_HARD in the namelist PARAM_ATMOS_VARS when you want to change the limit."

!!$          WNAME(1) = "Courant num. Z"
!!$          WNAME(2) = "Courant num. X"
!!$          WNAME(3) = "Courant num. Y"
!!$          call STATISTICS_detail( KA, KS, KE, IA, IS, IE, JA, JS, JE, 3, &
!!$                                  WNAME(:), WORK(:,:,:,:,:),               &
!!$                                  local=.true.                           )

          call PRC_abort
       endif

       if ( ATMOS_VARS_CHECKCFL_SOFT > 0.0_RP .AND. CFLMAX > ATMOS_VARS_CHECKCFL_SOFT ) then
          LOG_INFO("ATMOS_vars_monitor",*) "Courant number =", CFLMAX, " exceeded the soft limit =", ATMOS_VARS_CHECKCFL_SOFT
          LOG_ERROR("ATMOS_vars_monitor",*)                     "Courant number =", CFLMAX, " exceeded the soft limit =", ATMOS_VARS_CHECKCFL_SOFT
          LOG_ERROR_CONT(*)                     "Rank =", PRC_myrank

          WNAME(1) = "Courant num. Z"
          WNAME(2) = "Courant num. X"
          WNAME(3) = "Courant num. Y"
!!$          call STATISTICS_detail( KA, KS, KE, IA, IS, IE, JA, JS, JE, 3, &
!!$                                  WNAME(:), WORK(:,:,:,:,:),               &
!!$                                  local=.true.                           )
       endif
    endif

    return
  end subroutine ATMOS_vars_monitor

  !-----------------------------------------------------------------------------
  !> Create atmospheric restart file
  subroutine ATMOS_vars_restart_create

    return
  end subroutine ATMOS_vars_restart_create

  !-----------------------------------------------------------------------------
  !> Exit netCDF define mode
  subroutine ATMOS_vars_restart_enddef

    return
  end subroutine ATMOS_vars_restart_enddef

  !-----------------------------------------------------------------------------
  !> Close restart file
  subroutine ATMOS_vars_restart_close

    return
  end subroutine ATMOS_vars_restart_close

  !-----------------------------------------------------------------------------
  !> Define atmospheric variables in restart file
  subroutine ATMOS_vars_restart_def_var

    return
  end subroutine ATMOS_vars_restart_def_var

  !-----------------------------------------------------------------------------
  !> Write restart of atmospheric variables
  subroutine ATMOS_vars_restart_write

    return
  end subroutine ATMOS_vars_restart_write


  ! private
  subroutine allocate_3D( ary )
    use scale_const, only: &
       UNDEF => CONST_UNDEF
    real(RP), intent(inout), allocatable :: ary(:,:,:,:)

    if ( .not. allocated(ary) ) then
       allocate( ary(KA,IA,JA,ADM_lall) )
       ary(:,:,:,:) = UNDEF
    end if

    return
  end subroutine allocate_3D

  subroutine allocate_2D( ary )
    use scale_const, only: &
       UNDEF => CONST_UNDEF
    real(RP), intent(inout), allocatable :: ary(:,:,:)

    if ( .not. allocated(ary) ) then
       allocate( ary(IA,JA,ADM_lall) )
       ary(:,:,:) = UNDEF
    end if

    return
  end subroutine allocate_2D

  subroutine allocate_1D( ary )
    use scale_const, only: &
       UNDEF => CONST_UNDEF
    real(RP), intent(inout), allocatable :: ary(:)

    if ( .not. allocated(ary) ) then
       allocate( ary(KA) )
       ary(:) = UNDEF
    end if

    return
  end subroutine allocate_1D

end module mod_atmos_vars
