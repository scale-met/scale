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
  use scale_atmos_grid_cartesC_index
  use scale_index
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
  real(RP), public, target, allocatable :: DENS(:,:,:)   ! Density     [kg/m3]
  real(RP), public, target, allocatable :: MOMZ(:,:,:)   ! momentum z  [kg/m2/s]
  real(RP), public, target, allocatable :: MOMX(:,:,:)   ! momentum x  [kg/m2/s]
  real(RP), public, target, allocatable :: MOMY(:,:,:)   ! momentum y  [kg/m2/s]
  real(RP), public, target, allocatable :: RHOT(:,:,:)   ! DENS * POTT [K*kg/m3]
  real(RP), public, target, allocatable :: QTRC(:,:,:,:) ! ratio of mass of tracer to total mass[kg/kg]

  real(RP), public, target, allocatable :: DENS_avw(:,:,:)
  real(RP), public, target, allocatable :: MOMZ_avw(:,:,:)
  real(RP), public, target, allocatable :: MOMX_avw(:,:,:)
  real(RP), public, target, allocatable :: MOMY_avw(:,:,:)
  real(RP), public, target, allocatable :: RHOT_avw(:,:,:)
  real(RP), public, target, allocatable :: QTRC_avw(:,:,:,:)

  real(RP), public, pointer             :: DENS_av(:,:,:)
  real(RP), public, pointer             :: MOMZ_av(:,:,:)
  real(RP), public, pointer             :: MOMX_av(:,:,:)
  real(RP), public, pointer             :: MOMY_av(:,:,:)
  real(RP), public, pointer             :: RHOT_av(:,:,:)
  real(RP), public, pointer             :: QTRC_av(:,:,:,:)

  integer,  public                      :: I_QV
  real(RP), public, pointer             :: QV(:,:,:)
  real(RP), public, pointer             :: QC(:,:,:)
  real(RP), public, pointer             :: QR(:,:,:)
  real(RP), public, pointer             :: QI(:,:,:)
  real(RP), public, pointer             :: QS(:,:,:)
  real(RP), public, pointer             :: QG(:,:,:)
  real(RP), public, pointer             :: QH(:,:,:)

  real(RP), public, target, allocatable :: Qe(:,:,:,:)  !> mass ratio of hydrometors [kg/kg]

  ! reference state
  real(RP), public, allocatable :: DENS_ref(:,:,:)
  real(RP), public, allocatable :: POTT_ref(:,:,:)
  real(RP), public, allocatable :: TEMP_ref(:,:,:)
  real(RP), public, allocatable :: PRES_ref(:,:,:)
  real(RP), public, allocatable :: QV_ref(:,:,:)

 ! tendency by physical processes
  real(RP), public, allocatable :: DENS_tp(:,:,:)
  real(RP), public, allocatable :: MOMZ_tp(:,:,:)
  real(RP), public, allocatable :: RHOU_tp(:,:,:)
  real(RP), public, allocatable :: RHOV_tp(:,:,:)
  real(RP), public, allocatable :: RHOT_tp(:,:,:)
  real(RP), public, allocatable :: RHOH_p (:,:,:)
  real(RP), public, allocatable :: RHOQ_tp(:,:,:,:)

  ! (obsolute)
  real(RP), public, allocatable :: MOMX_tp(:,:,:)
  real(RP), public, allocatable :: MOMY_tp(:,:,:)


  ! public diagnostic variables
  real(RP), public, allocatable, target :: W    (:,:,:) !> velocity w [m/s]
  real(RP), public, allocatable, target :: U    (:,:,:) !> velocity u [m/s]
  real(RP), public, allocatable, target :: V    (:,:,:) !> velocity v [m/s]

  real(RP), public, allocatable, target :: POTT (:,:,:) !> potential temperature [K]
  real(RP), public, allocatable, target :: TEMP (:,:,:) !> temperature           [K]
  real(RP), public, allocatable, target :: PRES (:,:,:) !> pressure              [Pa=J/m3]
  real(RP), public, allocatable, target :: EXNER(:,:,:) !> Exner function (t/pt) [1]
  real(RP), public, allocatable, target :: PHYD (:,:,:) !> hydrostatic pressure  [Pa=J/m3]
  real(RP), public, allocatable, target :: PHYDH(:,:,:) !> hydrostatic pressure  [Pa=J/m3], layer interface

  real(RP), public, allocatable, target :: Qdry (:,:,:) !> dry air                [1]
  real(RP), public, allocatable, target :: Rtot (:,:,:) !> specific gass constant [J/kg/K]
  real(RP), public, allocatable, target :: CVtot(:,:,:) !> specific heat          [J/kg/K]
  real(RP), public, allocatable, target :: CPtot(:,:,:) !> specific heat          [J/kg/K]

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

  ! prognostic variables
  integer,     private, parameter   :: PV_nmax = 5
  type(Vinfo), private              :: PV_info(PV_nmax)
  integer,     private, allocatable :: PV_ID(:)

  data PV_info / &
       Vinfo( 'DENS', 'density',     'kg/m3',   3, 'ZXY',  'air_density' ), &
       Vinfo( 'MOMZ', 'momentum z',  'kg/m2/s', 3, 'ZHXY', 'upward_mass_flux_of_air' ), &
       Vinfo( 'MOMX', 'momentum x',  'kg/m2/s', 3, 'ZXHY', 'eastward_mass_flux_of_air' ), &
       Vinfo( 'MOMY', 'momentum y',  'kg/m2/s', 3, 'ZXYH', 'northward_mass_flux_of_air' ), &
       Vinfo( 'RHOT', 'rho * theta', 'kg/m3*K', 3, 'ZXY',  '' ) /


  ! private diagnostic variables
  real(RP), allocatable, target :: LHV  (:,:,:) !> latent heat for vaporization [J/kg]
  real(RP), allocatable, target :: LHS  (:,:,:) !> latent heat for sublimation  [J/kg]
  real(RP), allocatable, target :: LHF  (:,:,:) !> latent heat for fusion       [J/kg]

  real(RP), allocatable, target :: POTV (:,:,:) !> virtual pot. temp.      [K]
  real(RP), allocatable, target :: TEML (:,:,:) !> liquid water temp.      [K]
  real(RP), allocatable, target :: POTL (:,:,:) !> liquid water pot. temp. [K]
  real(RP), allocatable, target :: POTE (:,:,:) !> equivalent pot. temp.   [K]

  real(RP), allocatable, target :: QTOT (:,:,:) !> total water content  [1]
  real(RP), allocatable, target :: QHYD (:,:,:) !> total hydrometeornt  [1]
  real(RP), allocatable, target :: QLIQ (:,:,:) !> liquid water content [1]
  real(RP), allocatable, target :: QICE (:,:,:) !> ice water content    [1]

  real(RP), allocatable, target :: LWP  (:,:)   !> liquid water path  [g/m2]
  real(RP), allocatable, target :: IWP  (:,:)   !> ice    water path  [g/m2]
  real(RP), allocatable, target :: PW   (:,:)   !> precipitable water [g/m2]

  real(RP), allocatable, target :: PREC (:,:)   !> surface precipitation rate CP+MP(rain+snow) [kg/m2/s]
  real(RP), allocatable, target :: RAIN (:,:)   !> surface rain rate CP+MP [kg/m2/s]
  real(RP), allocatable, target :: SNOW (:,:)   !> surface snow rate CP+MP [kg/m2/s]

  real(RP), allocatable, target :: QSAT (:,:,:) !> saturation specific humidity        [1]
  real(RP), allocatable, target :: RHA  (:,:,:) !> relative humidity (liquid+ice)      [%]
  real(RP), allocatable, target :: RHL  (:,:,:) !> relative humidity against to liquid [%]
  real(RP), allocatable, target :: RHI  (:,:,:) !> relative humidity against to ice    [%]

  real(RP), allocatable, target :: VOR  (:,:,:) !> vertical vorticity    [1/s]
  real(RP), allocatable, target :: DIV  (:,:,:) !> divergence            [1/s]
  real(RP), allocatable, target :: HDIV (:,:,:) !> horizontal divergence [1/s]
  real(RP), allocatable, target :: Uabs (:,:,:) !> absolute velocity     [m/s]

  real(RP), allocatable, target :: N2   (:,:,:) !> squared Brunt-Vaisala frequency [/s2]
  real(RP), allocatable, target :: PBLH (:,:)   !> PBL height [m]

  real(RP), allocatable, target :: MSE  (:,:,:) !> moist static energy [m2/s2]
  real(RP), allocatable, target :: TDEW (:,:,:) !> dew point [K]

  real(RP), allocatable, target :: CAPE (:,:)   !> CAPE       [m2/s2]
  real(RP), allocatable, target :: CIN  (:,:)   !> CIN        [m2/s2]
  real(RP), allocatable, target :: LCL  (:,:)   !> LCL height [m]
  real(RP), allocatable, target :: LFC  (:,:)   !> LFC height [m]
  real(RP), allocatable, target :: LNB  (:,:)   !> LNB height [m]

  real(RP), allocatable, target :: ENGT (:,:,:) !> total     energy [J/m3]
  real(RP), allocatable, target :: ENGP (:,:,:) !> potential energy [J/m3]
  real(RP), allocatable, target :: ENGK (:,:,:) !> kinetic   energy [J/m3]
  real(RP), allocatable, target :: ENGI (:,:,:) !> internal  energy [J/m3]

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

  real(RP), allocatable, target :: DENS_PRIM(:,:,:) !> horiz. deviation of density    [kg/m3]
  real(RP), allocatable, target :: W_PRIM   (:,:,:) !> horiz. deviation of w          [m/s]
  real(RP), allocatable, target :: U_PRIM   (:,:,:) !> horiz. deviation of u          [m/s]
  real(RP), allocatable, target :: V_PRIM   (:,:,:) !> horiz. deviation of v          [m/s]
  real(RP), allocatable, target :: PT_PRIM  (:,:,:) !> horiz. deviation of pot. temp. [K]
  real(RP), allocatable, target :: W_PRIM2  (:,:,:) !> variance of w                  [m2/s2]
  real(RP), allocatable, target :: PT_W_PRIM(:,:,:) !> resolved scale heat flux       [W/s]
  real(RP), allocatable, target :: W_PRIM3  (:,:,:) !> skewness of w                  [m3/s3]
  real(RP), allocatable, target :: TKE_RS   (:,:,:) !> resolved scale TKE             [m2/s2]

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
  integer, private              :: PV_HIST_id (PV_nmax) !> prognostic variables
  integer, private              :: PV_MONIT_id(PV_nmax)
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
  real(RP), private, target, allocatable :: ZERO(:,:,:)


  ! for restart
  integer, private :: restart_fid = -1  ! file ID
  logical, private :: ATMOS_RESTART_IN_CHECK_COORDINATES = .true.


  real(RP), private, allocatable :: WORK3D(:,:,:)
  real(RP), private, allocatable :: WORK2D(:,:)
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
    use mod_atmos_admin, only: &
       ATMOS_USE_AVERAGE
    use mod_atmos_dyn_vars, only: &
       ATMOS_DYN_vars_setup
    use mod_atmos_phy_mp_vars, only: &
       ATMOS_PHY_MP_vars_setup
    use mod_atmos_phy_ae_vars, only: &
       ATMOS_PHY_AE_vars_setup
    use mod_atmos_phy_ch_vars, only: &
       ATMOS_PHY_CH_vars_setup
    use mod_atmos_phy_rd_vars, only: &
       ATMOS_PHY_RD_vars_setup
    use mod_atmos_phy_sf_vars, only: &
       ATMOS_PHY_SF_vars_setup
    use mod_atmos_phy_tb_vars, only: &
       ATMOS_PHY_TB_vars_setup
    use mod_atmos_phy_bl_vars, only: &
       ATMOS_PHY_BL_vars_setup
    use mod_atmos_phy_cp_vars, only: &
       ATMOS_PHY_CP_vars_setup
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
    integer :: iv, iq
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("ATMOS_vars_setup",*) 'Setup'

    allocate( DENS(KA,IA,JA)    )
    allocate( MOMZ(KA,IA,JA)    )
    allocate( MOMX(KA,IA,JA)    )
    allocate( MOMY(KA,IA,JA)    )
    allocate( RHOT(KA,IA,JA)    )
    allocate( QTRC(KA,IA,JA,max(QA,1)) )

    if ( ATMOS_USE_AVERAGE ) then
       allocate( DENS_avw(KA,IA,JA)    )
       allocate( MOMZ_avw(KA,IA,JA)    )
       allocate( MOMX_avw(KA,IA,JA)    )
       allocate( MOMY_avw(KA,IA,JA)    )
       allocate( RHOT_avw(KA,IA,JA)    )
       allocate( QTRC_avw(KA,IA,JA,max(QA,1)) )

       DENS_av => DENS_avw
       MOMZ_av => MOMZ_avw
       MOMX_av => MOMX_avw
       MOMY_av => MOMY_avw
       RHOT_av => RHOT_avw
       QTRC_av => QTRC_avw
    else
       DENS_av => DENS
       MOMZ_av => MOMZ
       MOMX_av => MOMX
       MOMY_av => MOMY
       RHOT_av => RHOT
       QTRC_av => QTRC
    endif

    allocate( DENS_tp(KA,IA,JA)    )
    allocate( MOMZ_tp(KA,IA,JA)    )
    allocate( RHOU_tp(KA,IA,JA)    )
    allocate( RHOV_tp(KA,IA,JA)    )
    allocate( RHOT_tp(KA,IA,JA)    )
    allocate( RHOH_p (KA,IA,JA)    )
    allocate( RHOQ_tp(KA,IA,JA,max(QA,1)) )

    allocate( W(KA,IA,JA) )
    allocate( U(KA,IA,JA) )
    allocate( V(KA,IA,JA) )
    W(:,:,:) = UNDEF
    U(:,:,:) = UNDEF
    V(:,:,:) = UNDEF

    allocate( POTT (KA,IA,JA) )
    allocate( TEMP (KA,IA,JA) )
    allocate( PRES (KA,IA,JA) )
    allocate( EXNER(KA,IA,JA) )
    allocate( PHYD (KA,IA,JA) )
    allocate( PHYDH(0:KA,IA,JA) )
    POTT (:,:,:) = UNDEF
    TEMP (:,:,:) = UNDEF
    PRES (:,:,:) = UNDEF
    EXNER(:,:,:) = UNDEF
    PHYD (:,:,:) = UNDEF
    PHYDH(:,:,:) = UNDEF

    allocate( Qdry (KA,IA,JA) )
    allocate( Rtot (KA,IA,JA) )
    allocate( CVtot(KA,IA,JA) )
    allocate( CPtot(KA,IA,JA) )
    Qdry (:,:,:) = UNDEF
    Rtot (:,:,:) = UNDEF
    CVtot(:,:,:) = UNDEF
    CPtot(:,:,:) = UNDEF

    ! obsolute
    allocate( MOMX_tp(KA,IA,JA)    )
    allocate( MOMY_tp(KA,IA,JA)    )


    MOMZ(1:KS-1,:,:) = 0.0_RP
    MOMZ(KE:KA,:,:) = 0.0_RP

    allocate( WORK3D(KA,IA,JA) )
    allocate( WORK2D(   IA,JA) )
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
    do iv = 1, PV_nmax
       LOG_INFO_CONT('(1x,A,I3,A,A24,A,A48,A,A12,A)') &
                  'NO.',iv,'|',PV_info(iv)%NAME,'|', PV_info(iv)%DESC,'[', PV_info(iv)%UNIT,']'
    enddo
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

    call ATMOS_DYN_vars_setup
    call ATMOS_PHY_MP_vars_setup
    call ATMOS_PHY_AE_vars_setup
    call ATMOS_PHY_CH_vars_setup
    call ATMOS_PHY_RD_vars_setup
    call ATMOS_PHY_SF_vars_setup
    call ATMOS_PHY_TB_vars_setup
    call ATMOS_PHY_BL_vars_setup
    call ATMOS_PHY_CP_vars_setup


    ! water content
    if ( ATMOS_HYDROMETEOR_dry ) then
       allocate( ZERO(KA,IA,JA) )
!OCL XFILL
       ZERO(:,:,:) = 0.0_RP

       QV => ZERO
       QC => ZERO
       QR => ZERO
       QI => ZERO
       QS => ZERO
       QG => ZERO
       QH => ZERO

       moist = .false.
    else
       allocate( Qe(KA,IA,JA,N_HYD) )
!OCL XFILL
       Qe(:,:,:,:) = UNDEF

       QV => QTRC_av(:,:,:,I_QV)
       QC => Qe(:,:,:,I_HC)
       QR => Qe(:,:,:,I_HR)
       QI => Qe(:,:,:,I_HI)
       QS => Qe(:,:,:,I_HS)
       QG => Qe(:,:,:,I_HG)
       QH => Qe(:,:,:,I_HH)

       moist = .true.
    end if


    DV_calclated(DV_nmax) = .false.

    !-----< history output setup >-----
    allocate( QP_HIST_id( max(QA,1) ) )
    allocate( QP_MONIT_id( max(QA,1) ) )
    PV_HIST_id (:) = -1
    PV_MONIT_id(:) = -1
    QP_HIST_id (:) = -1
    QP_MONIT_id(:) = -1
    DV_HIST_id (:) = -1
    DV_MONIT_id(:) = -1


    do iv = 1, PV_nmax
       call FILE_HISTORY_reg( PV_info(iv)%NAME, PV_info(iv)%DESC, PV_info(iv)%UNIT, PV_HIST_id(iv), dim_type=PV_info(iv)%dim_type, standard_name=PV_info(iv)%STDNAME )
    end do

    do iq = 1, QA
       call FILE_HISTORY_reg( TRACER_NAME(iq), TRACER_DESC(iq), TRACER_UNIT(iq), QP_HIST_id(iq), dim_type='ZXY' )
    enddo

    do iv = 1, DV_nmax
       call FILE_HISTORY_reg( DV_info(iv)%NAME, DV_info(iv)%DESC, DV_info(iv)%UNIT, DV_HIST_id(iv), dim_type=DV_info(iv)%dim_type, standard_name=DV_info(iv)%STDNAME )
    end do


    !-----< monitor output setup >-----
    do iv = 1, PV_nmax
       call MONITOR_reg( PV_info(iv)%NAME, PV_info(iv)%DESC, PV_info(iv)%UNIT, & ! (in)
                         PV_MONIT_id(iv),                                      & ! (out)
                         dim_type=PV_info(iv)%dim_type, isflux=.false.         ) ! (in)
    end do
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
  subroutine ATMOS_vars_fillhalo( &
       FILL_BND )
    use scale_comm_cartesC, only: &
       COMM_vars8, &
       COMM_wait
    implicit none

    logical, intent(in), optional :: FILL_BND

    logical :: FILL_BND_
    integer :: i, j, iq
    !---------------------------------------------------------------------------

    FILL_BND_ = .false.
    if ( present(FILL_BND) ) FILL_BND_ = FILL_BND

    !$omp parallel do private(i,j) OMP_SCHEDULE_ collapse(2)
    do j  = JSB, JEB
    do i  = ISB, IEB
       DENS(   1:KS-1,i,j) = DENS(KS,i,j)
       MOMZ(   1:KS-1,i,j) = MOMZ(KS,i,j)
       MOMX(   1:KS-1,i,j) = MOMX(KS,i,j)
       MOMY(   1:KS-1,i,j) = MOMY(KS,i,j)
       RHOT(   1:KS-1,i,j) = RHOT(KS,i,j)
       DENS(KE+1:KA,  i,j) = DENS(KE,i,j)
       MOMZ(KE+1:KA,  i,j) = MOMZ(KE,i,j)
       MOMX(KE+1:KA,  i,j) = MOMX(KE,i,j)
       MOMY(KE+1:KA,  i,j) = MOMY(KE,i,j)
       RHOT(KE+1:KA,  i,j) = RHOT(KE,i,j)
    enddo
    enddo

    !$omp parallel do private(i,j,iq) OMP_SCHEDULE_ collapse(3)
    do iq = 1, QA
    do j  = JSB, JEB
    do i  = ISB, IEB
       QTRC(   1:KS-1,i,j,iq) = QTRC(KS,i,j,iq)
       QTRC(KE+1:KA,  i,j,iq) = QTRC(KE,i,j,iq)
    enddo
    enddo
    enddo

    call COMM_vars8( DENS(:,:,:), 1 )
    call COMM_vars8( MOMZ(:,:,:), 2 )
    call COMM_vars8( MOMX(:,:,:), 3 )
    call COMM_vars8( MOMY(:,:,:), 4 )
    call COMM_vars8( RHOT(:,:,:), 5 )
    call COMM_wait ( DENS(:,:,:), 1, FILL_BND_ )
    call COMM_wait ( MOMZ(:,:,:), 2, FILL_BND_ )
    call COMM_wait ( MOMX(:,:,:), 3, FILL_BND_ )
    call COMM_wait ( MOMY(:,:,:), 4, FILL_BND_ )
    call COMM_wait ( RHOT(:,:,:), 5, FILL_BND_ )

    do iq = 1, QA
       call COMM_vars8( QTRC(:,:,:,iq), iq )
    enddo
    do iq = 1, QA
       call COMM_wait ( QTRC(:,:,:,iq), iq, FILL_BND_ )
    enddo

    return
  end subroutine ATMOS_vars_fillhalo

  !-----------------------------------------------------------------------------
  !> Open restart file for reading atmospheric variables
  subroutine ATMOS_vars_restart_open
    use scale_prc, only: &
       PRC_abort
    use scale_const, only: &
       GRAV  => CONST_GRAV
    use scale_time, only: &
       TIME_gettimelabel
    use scale_file_cartesC, only: &
       FILE_CARTESC_open, &
       FILE_CARTESC_check_coordinates
    use mod_atmos_admin, only: &
       ATMOS_USE_AVERAGE, &
       ATMOS_sw_dyn,      &
       ATMOS_sw_phy_mp,   &
       ATMOS_sw_phy_ae,   &
       ATMOS_sw_phy_ch,   &
       ATMOS_sw_phy_rd,   &
       ATMOS_sw_phy_sf,   &
       ATMOS_sw_phy_tb,   &
       ATMOS_sw_phy_cp
    use mod_atmos_dyn_vars, only: &
       ATMOS_DYN_vars_restart_open
    use mod_atmos_phy_mp_vars, only: &
       ATMOS_PHY_MP_vars_restart_open
    use mod_atmos_phy_ae_vars, only: &
       ATMOS_PHY_AE_vars_restart_open
    use mod_atmos_phy_ch_vars, only: &
       ATMOS_PHY_CH_vars_restart_open
    use mod_atmos_phy_rd_vars, only: &
       ATMOS_PHY_RD_vars_restart_open
    use mod_atmos_phy_sf_vars, only: &
       ATMOS_PHY_SF_vars_restart_open
    use mod_atmos_phy_tb_vars, only: &
       ATMOS_PHY_TB_vars_restart_open
    use mod_atmos_phy_cp_vars, only: &
       ATMOS_PHY_CP_vars_restart_open
    use mod_cpl_admin, only: &
       CPL_sw
    implicit none

    character(len=19)     :: timelabel
    character(len=H_LONG) :: basename
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("ATMOS_vars_restart_open",*) 'Open restart file (ATMOS) '

    if ( ATMOS_RESTART_IN_BASENAME /= '' ) then

       if ( ATMOS_RESTART_IN_POSTFIX_TIMELABEL ) then
          call TIME_gettimelabel( timelabel )
          basename = trim(ATMOS_RESTART_IN_BASENAME)//'_'//trim(timelabel)
       else
          basename = trim(ATMOS_RESTART_IN_BASENAME)
       endif

       LOG_INFO("ATMOS_vars_restart_open",*) 'basename: ', trim(basename)

       call FILE_CARTESC_open( basename, restart_fid, aggregate=ATMOS_RESTART_IN_AGGREGATE )

       if ( ATMOS_RESTART_IN_CHECK_COORDINATES ) then
          call FILE_CARTESC_check_coordinates( restart_fid, atmos=.true. )
       end if

    else
       LOG_ERROR("ATMOS_vars_restart_open",*) 'restart file for atmosphere is not specified. STOP!'
       call PRC_abort
    endif

    if ( ATMOS_USE_AVERAGE ) then
       DENS_av(:,:,:)   = DENS(:,:,:)
       MOMZ_av(:,:,:)   = MOMZ(:,:,:)
       MOMX_av(:,:,:)   = MOMX(:,:,:)
       MOMY_av(:,:,:)   = MOMY(:,:,:)
       RHOT_av(:,:,:)   = RHOT(:,:,:)
       QTRC_av(:,:,:,:) = QTRC(:,:,:,:)
    endif

    if( ATMOS_sw_dyn )    call ATMOS_DYN_vars_restart_open
    if( ATMOS_sw_phy_mp ) call ATMOS_PHY_MP_vars_restart_open
    if( ATMOS_sw_phy_ae ) call ATMOS_PHY_AE_vars_restart_open
    if( ATMOS_sw_phy_ch ) call ATMOS_PHY_CH_vars_restart_open
    if( ATMOS_sw_phy_rd ) call ATMOS_PHY_RD_vars_restart_open
    if( ATMOS_sw_phy_sf .and. (.not. CPL_sw) ) call ATMOS_PHY_SF_vars_restart_open
    if( ATMOS_sw_phy_tb ) call ATMOS_PHY_TB_vars_restart_open
    if( ATMOS_sw_phy_cp ) call ATMOS_PHY_CP_vars_restart_open

    return
  end subroutine ATMOS_vars_restart_open

  !-----------------------------------------------------------------------------
  !> Read restart of atmospheric variables
  subroutine ATMOS_vars_restart_read
    use scale_prc, only: &
       PRC_abort
    use scale_file, only: &
       FILE_get_AGGREGATE
    use scale_file_cartesC, only: &
       FILE_CARTESC_read, &
       FILE_CARTESC_flush
    use mod_atmos_admin, only: &
       ATMOS_USE_AVERAGE, &
       ATMOS_sw_dyn,      &
       ATMOS_sw_phy_mp,   &
       ATMOS_sw_phy_ae,   &
       ATMOS_sw_phy_ch,   &
       ATMOS_sw_phy_rd,   &
       ATMOS_sw_phy_sf,   &
       ATMOS_sw_phy_tb,   &
       ATMOS_sw_phy_cp
    use mod_atmos_dyn_vars, only: &
       ATMOS_DYN_vars_restart_read
    use mod_atmos_phy_mp_vars, only: &
       ATMOS_PHY_MP_vars_restart_read
    use mod_atmos_phy_ae_vars, only: &
       ATMOS_PHY_AE_vars_restart_read
    use mod_atmos_phy_ch_vars, only: &
       ATMOS_PHY_CH_vars_restart_read
    use mod_atmos_phy_rd_vars, only: &
       ATMOS_PHY_RD_vars_restart_read
    use mod_atmos_phy_sf_vars, only: &
       ATMOS_PHY_SF_vars_restart_read
    use mod_atmos_phy_tb_vars, only: &
       ATMOS_PHY_TB_vars_restart_read
    use mod_atmos_phy_cp_vars, only: &
       ATMOS_PHY_CP_vars_restart_read
    use mod_cpl_admin, only: &
       CPL_sw
    implicit none

    integer  :: i, j, iq
    !---------------------------------------------------------------------------

    if ( restart_fid /= -1 ) then
       LOG_NEWLINE
       LOG_INFO("ATMOS_vars_restart_read",*) 'Read from restart file (ATMOS) '

       call FILE_CARTESC_read( restart_fid, PV_info(I_DENS)%NAME, 'ZXY', & ! [IN]
                               DENS(:,:,:)                               ) ! [OUT]
       call FILE_CARTESC_read( restart_fid, PV_info(I_MOMZ)%NAME, 'ZHXY', & ! [IN]
                               MOMZ(:,:,:)                                ) ! [OUT]
       call FILE_CARTESC_read( restart_fid, PV_info(I_MOMX)%NAME, 'ZXHY', & ! [IN]
                               MOMX(:,:,:)                                ) ! [OUT]
       call FILE_CARTESC_read( restart_fid, PV_info(I_MOMY)%NAME, 'ZXYH', & ! [IN]
                               MOMY(:,:,:)                                ) ! [OUT]
       call FILE_CARTESC_read( restart_fid, PV_info(I_RHOT)%NAME, 'ZXY', & ! [IN]
                               RHOT(:,:,:)                               ) ! [OUT]

       do iq = 1, QA
          call FILE_CARTESC_read( restart_fid, TRACER_NAME(iq), 'ZXY', & ! [IN]
                                  QTRC(:,:,:,iq)                       ) ! [OUT]
       enddo

       if ( FILE_get_AGGREGATE(restart_fid) ) then
          call FILE_CARTESC_flush( restart_fid ) ! X/Y halos have been read from file

          ! fill k halos
          do j  = 1, JA
          do i  = 1, IA
             DENS(   1:KS-1,i,j) = DENS(KS,i,j)
             MOMZ(   1:KS-2,i,j) = MOMZ(KS-1,i,j)
             MOMX(   1:KS-1,i,j) = MOMX(KS,i,j)
             MOMY(   1:KS-1,i,j) = MOMY(KS,i,j)
             RHOT(   1:KS-1,i,j) = RHOT(KS,i,j)
             DENS(KE+1:KA,  i,j) = DENS(KE,i,j)
             MOMZ(KE+1:KA,  i,j) = MOMZ(KE,i,j)
             MOMX(KE+1:KA,  i,j) = MOMX(KE,i,j)
             MOMY(KE+1:KA,  i,j) = MOMY(KE,i,j)
             RHOT(KE+1:KA,  i,j) = RHOT(KE,i,j)
          enddo
          enddo
       else
          call ATMOS_vars_fillhalo
       end if

       call ATMOS_vars_total
    else
       LOG_ERROR("ATMOS_vars_restart_read",*) 'invalid restart file ID for atmosphere. STOP!'
       call PRC_abort
    endif

    if ( ATMOS_USE_AVERAGE ) then
       DENS_av(:,:,:)   = DENS(:,:,:)
       MOMZ_av(:,:,:)   = MOMZ(:,:,:)
       MOMX_av(:,:,:)   = MOMX(:,:,:)
       MOMY_av(:,:,:)   = MOMY(:,:,:)
       RHOT_av(:,:,:)   = RHOT(:,:,:)
       QTRC_av(:,:,:,:) = QTRC(:,:,:,:)
    endif

    if ( ATMOS_sw_dyn )    call ATMOS_DYN_vars_restart_read
    if ( ATMOS_sw_phy_mp ) call ATMOS_PHY_MP_vars_restart_read
    if ( ATMOS_sw_phy_ae ) call ATMOS_PHY_AE_vars_restart_read
    if ( ATMOS_sw_phy_ch ) call ATMOS_PHY_CH_vars_restart_read
    if ( ATMOS_sw_phy_rd ) call ATMOS_PHY_RD_vars_restart_read
    if ( ATMOS_sw_phy_sf .and. (.not. CPL_sw) ) call ATMOS_PHY_SF_vars_restart_read
    if ( ATMOS_sw_phy_tb ) call ATMOS_PHY_TB_vars_restart_read
    if ( ATMOS_sw_phy_cp ) call ATMOS_PHY_CP_vars_restart_read

    return
  end subroutine ATMOS_vars_restart_read

  !-----------------------------------------------------------------------------
  !> Set pressure for history output
  subroutine ATMOS_vars_history_setpres
    use scale_atmos_grid_cartesC_real, only: &
       REAL_CZ => ATMOS_GRID_CARTESC_REAL_CZ, &
       REAL_Z1 => ATMOS_GRID_CARTESC_REAL_Z1
    use scale_topography, only: &
       TOPO_Zsfc
    use scale_atmos_bottom, only: &
       BOTTOM_estimate => ATMOS_BOTTOM_estimate
    use scale_file_history_cartesC, only: &
       FILE_HISTORY_CARTESC_set_pres
    implicit none

    real(RP) :: SFC_DENS(IA,JA)
    real(RP) :: SFC_PRES(IA,JA)
    !---------------------------------------------------------------------------

    call BOTTOM_estimate( KA, KS, KE, IA, ISB, IEB, JA, JSB, JEB, &
                          DENS_av(:,:,:), PRES(:,:,:),                  & ! [IN]
                          REAL_CZ(:,:,:), TOPO_Zsfc(:,:), REAL_Z1(:,:), & ! [IN]
                          SFC_DENS(:,:), SFC_PRES(:,:)                  ) ! [OUT]

    call FILE_HISTORY_CARTESC_set_pres( PHYD    (:,:,:), & ! [IN]
                                        PHYDH   (:,:,:), & ! [IN]
                                        SFC_PRES(:,:)    ) ! [IN]

    return
  end subroutine ATMOS_vars_history_setpres

  !-----------------------------------------------------------------------------
  !> Check and compare between last data and sample data
  subroutine ATMOS_vars_restart_check
    use scale_prc, only: &
       PRC_myrank
    use scale_file, only: &
       FILE_get_AGGREGATE
    use scale_file_cartesC, only: &
       FILE_CARTESC_open, &
       FILE_CARTESC_read, &
       FILE_CARTESC_flush, &
       FILE_CARTESC_close
    implicit none

    real(RP) :: DENS_check(KA,IA,JA)    ! Density    [kg/m3]
    real(RP) :: MOMZ_check(KA,IA,JA)    ! momentum z [kg/s/m2]
    real(RP) :: MOMX_check(KA,IA,JA)    ! momentum x [kg/s/m2]
    real(RP) :: MOMY_check(KA,IA,JA)    ! momentum y [kg/s/m2]
    real(RP) :: RHOT_check(KA,IA,JA)    ! DENS * POTT [K*kg/m3]
    real(RP) :: QTRC_check(KA,IA,JA,QA) ! tracer mixing ratio [kg/kg]

    character(len=H_LONG) :: basename

    logical :: datacheck
    integer :: k, i, j, iq
    integer :: fid
    !---------------------------------------------------------------------------

    call PROF_rapstart('Debug')

    LOG_INFO("ATMOS_vars_restart_check",*) 'Compare last Data with ', trim(ATMOS_RESTART_CHECK_BASENAME), 'on PE=', PRC_myrank
    LOG_INFO("ATMOS_vars_restart_check",*) 'criterion = ', ATMOS_RESTART_CHECK_CRITERION
    datacheck = .true.

    basename = ATMOS_RESTART_CHECK_BASENAME

    call FILE_CARTESC_open( basename, fid )

    call FILE_CARTESC_read( fid, 'DENS', 'ZXY' , DENS_check(:,:,:) )
    call FILE_CARTESC_read( fid, 'MOMZ', 'ZHXY', MOMZ_check(:,:,:) )
    call FILE_CARTESC_read( fid, 'MOMX', 'ZXHY', MOMX_check(:,:,:) )
    call FILE_CARTESC_read( fid, 'MOMY', 'ZXYH', MOMY_check(:,:,:) )
    call FILE_CARTESC_read( fid, 'RHOT', 'ZXY' , RHOT_check(:,:,:) )
    do iq = 1, QA
       call FILE_CARTESC_read( fid, TRACER_NAME(iq), 'ZXY', QTRC_check(:,:,:,iq) )
    end do
    if ( FILE_get_AGGREGATE(fid) ) call FILE_CARTESC_flush( fid )

    call FILE_CARTESC_close( fid ) ! [IN]

    do k = KS, KE
    do j = JS, JE
    do i = IS, IE
       if ( abs( DENS(k,i,j)-DENS_check(k,i,j) ) > ATMOS_RESTART_CHECK_CRITERION ) then
          LOG_ERROR("ATMOS_vars_restart_check",*) 'there is the difference  : ', DENS(k,i,j)-DENS_check(k,i,j)
          LOG_ERROR_CONT(*) 'at (PE-id,k,i,j,varname) : ', PRC_myrank, k, i, j, 'DENS'
          datacheck = .false.
       endif
    enddo
    enddo
    enddo

    do k = KS-1, KE
    do j = JS, JE
    do i = IS, IE
       if ( abs( MOMZ(k,i,j)-MOMZ_check(k,i,j) ) > ATMOS_RESTART_CHECK_CRITERION ) then
          LOG_ERROR("ATMOS_vars_restart_check",*) 'there is the difference  : ', MOMZ(k,i,j)-MOMZ_check(k,i,j)
          LOG_ERROR_CONT(*) 'at (PE-id,k,i,j,varname) : ', PRC_myrank, k, i, j, 'MOMZ'
          datacheck = .false.
       endif
    enddo
    enddo
    enddo

    do k = KS, KE
    do j = JS, JE
    do i = IS, IE
       if ( abs( MOMX(k,i,j)-MOMX_check(k,i,j) ) > ATMOS_RESTART_CHECK_CRITERION ) then
          LOG_ERROR("ATMOS_vars_restart_check",*) 'there is the difference  : ', MOMX(k,i,j)-MOMX_check(k,i,j)
          LOG_ERROR_CONT(*) 'at (PE-id,k,i,j,varname) : ', PRC_myrank, k, i, j, 'MOMX'
          datacheck = .false.
       endif
    enddo
    enddo
    enddo

    do k = KS, KE
    do j = JS, JE
    do i = IS, IE
       if ( abs( MOMY(k,i,j)-MOMY_check(k,i,j) ) > ATMOS_RESTART_CHECK_CRITERION ) then
          LOG_ERROR("ATMOS_vars_restart_check",*) 'there is the difference  : ', MOMY(k,i,j)-MOMY_check(k,i,j)
          LOG_ERROR_CONT(*) 'at (PE-id,k,i,j,varname) : ', PRC_myrank, k, i, j, 'MOMY'
          datacheck = .false.
       endif
    enddo
    enddo
    enddo

    do k = KS, KE
    do j = JS, JE
    do i = IS, IE
       if ( abs( RHOT(k,i,j)-RHOT_check(k,i,j) ) > ATMOS_RESTART_CHECK_CRITERION ) then
          LOG_ERROR("ATMOS_vars_restart_check",*) 'there is the difference  : ', RHOT(k,i,j)-RHOT_check(k,i,j)
          LOG_ERROR_CONT(*) 'at (PE-id,k,i,j,varname) : ', PRC_myrank, k, i, j, 'RHOT'
          datacheck = .false.
       endif
    enddo
    enddo
    enddo

    do iq = 1, QA
       do k = KS, KE
       do j = JS, JE
       do i = IS, IE
          if ( abs( QTRC(k,i,j,iq)-QTRC_check(k,i,j,iq) ) > ATMOS_RESTART_CHECK_CRITERION ) then
             LOG_ERROR("ATMOS_vars_restart_check",*) 'there is the difference  : ', QTRC(k,i,j,iq)-QTRC_check(k,i,j,iq)
             LOG_ERROR_CONT(*) 'at (PE-id,k,i,j,varname) : ', PRC_myrank, k, i, j, TRACER_NAME(iq)
             datacheck = .false.
          endif
       enddo
       enddo
       enddo
    enddo

    if (datacheck) then
       LOG_INFO("ATMOS_vars_restart_check",*) 'Data Check Clear.'
    else
       LOG_INFO("ATMOS_vars_restart_check",*) 'Data Check Failed. See std. output.'
       LOG_ERROR("ATMOS_vars_restart_check",*) 'Data Check Failed.'
    endif

    call PROF_rapend('Debug')

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

    ! value check for prognostic variables
    if ( ATMOS_VARS_CHECKRANGE ) then
       call VALCHECK( DENS(:,:,:),    0.0_RP,    2.0_RP, PV_info(I_DENS)%NAME, __FILE__, __LINE__ )
       call VALCHECK( MOMZ(:,:,:), -200.0_RP,  200.0_RP, PV_info(I_MOMZ)%NAME, __FILE__, __LINE__ )
       call VALCHECK( MOMX(:,:,:), -200.0_RP,  200.0_RP, PV_info(I_MOMX)%NAME, __FILE__, __LINE__ )
       call VALCHECK( MOMY(:,:,:), -200.0_RP,  200.0_RP, PV_info(I_MOMY)%NAME, __FILE__, __LINE__ )
       call VALCHECK( RHOT(:,:,:),    0.0_RP, 1000.0_RP, PV_info(I_RHOT)%NAME, __FILE__, __LINE__ )
    endif

    ! history output of prognostic variables
    call FILE_HISTORY_put  ( PV_HIST_id(I_DENS), DENS(:,:,:) )
    call FILE_HISTORY_put  ( PV_HIST_id(I_MOMZ), MOMZ(:,:,:) )
    call FILE_HISTORY_put  ( PV_HIST_id(I_MOMX), MOMX(:,:,:) )
    call FILE_HISTORY_put  ( PV_HIST_id(I_MOMY), MOMY(:,:,:) )
    call FILE_HISTORY_put  ( PV_HIST_id(I_RHOT), RHOT(:,:,:) )
    do iq = 1, QA
       call FILE_HISTORY_put  ( QP_HIST_id(iq), QTRC(:,:,:,iq) )
    enddo


    ! history output of diagnostic variables
    call FILE_HISTORY_put  ( DV_HIST_id(I_W    ), W(:,:,:)     )
    call FILE_HISTORY_put  ( DV_HIST_id(I_U    ), U(:,:,:)     )
    call FILE_HISTORY_put  ( DV_HIST_id(I_V    ), V(:,:,:)     )
    call FILE_HISTORY_put  ( DV_HIST_id(I_POTT ), POTT(:,:,:)  )
    call FILE_HISTORY_put  ( DV_HIST_id(I_TEMP ), TEMP(:,:,:)  )
    call FILE_HISTORY_put  ( DV_HIST_id(I_PRES ), PRES(:,:,:)  )

    call FILE_HISTORY_put  ( DV_HIST_id(I_EXNER), EXNER(:,:,:) )
    call FILE_HISTORY_put  ( DV_HIST_id(I_PHYD ), PHYD(:,:,:)  )

    call FILE_HISTORY_put  ( DV_HIST_id(I_QDRY ), QDRY(:,:,:)  )
    call FILE_HISTORY_put  ( DV_HIST_id(I_RTOT ), RTOT(:,:,:)  )
    call FILE_HISTORY_put  ( DV_HIST_id(I_CVTOT), CVTOT(:,:,:) )
    call FILE_HISTORY_put  ( DV_HIST_id(I_CPTOT), CPTOT(:,:,:) )

    do iv = I_CPTOT+1, DV_nmax
       if ( DV_HIST_id(iv) > 0 ) then
          call FILE_HISTORY_query( DV_HIST_id(iv), do_put )

          if ( do_put ) then
             select case( DV_info(iv)%ndims )
             case( 3 )
                call ATMOS_vars_get_diagnostic( DV_info(iv)%NAME, WORK3D(:,:,:) )
                call FILE_HISTORY_put( DV_HIST_id(iv), WORK3D(:,:,:) )
             case( 2 )
                call ATMOS_vars_get_diagnostic( DV_info(iv)%NAME, WORK2D(:,:) )
                call FILE_HISTORY_put( DV_HIST_id(iv), WORK2D(:,:) )
             case( 1 )
                call ATMOS_vars_get_diagnostic( DV_info(iv)%NAME, WORK1D(:) )
                call FILE_HISTORY_put( DV_HIST_id(iv), WORK1D(:) )
             end select
          endif
       endif
    enddo

    if ( moist ) &
         call ATMOS_PHY_MP_vars_history( DENS_av(:,:,:), TEMP(:,:,:), QTRC_av(:,:,:,:) )
!    if ( .false. ) then
!       call ATMOS_vars_get_diagnostic( "RH", WORK3D(:,:,:) )
!       call ATMOS_PHY_AE_vars_history( QTRC_av(:,:,:,:), WORK3D(:,:,:) )
!    end if

    call PROF_rapend  ('ATM_History', 1)

    return
  end subroutine ATMOS_vars_history

  !-----------------------------------------------------------------------------
  !> Budget monitor for atmosphere
  subroutine ATMOS_vars_total
    use scale_const, only: &
       GRAV  => CONST_GRAV,  &
       CVdry => CONST_CVdry
    use scale_statistics, only: &
       STATISTICS_checktotal, &
       STATISTICS_total
    use scale_atmos_grid_cartesC_real, only: &
       ATMOS_GRID_CARTESC_REAL_VOL,       &
       ATMOS_GRID_CARTESC_REAL_TOTVOL,    &
       ATMOS_GRID_CARTESC_REAL_VOLWXY,    &
       ATMOS_GRID_CARTESC_REAL_TOTVOLWXY, &
       ATMOS_GRID_CARTESC_REAL_VOLZUY,    &
       ATMOS_GRID_CARTESC_REAL_TOTVOLZUY, &
       ATMOS_GRID_CARTESC_REAL_VOLZXV,    &
       ATMOS_GRID_CARTESC_REAL_TOTVOLZXV
    implicit none

    real(RP) :: RHOQ(KA,IA,JA)

    integer  :: i, j, k, iq
    !---------------------------------------------------------------------------

    if ( STATISTICS_checktotal ) then

       call STATISTICS_total( KA, KS, KE, IA, IS, IE, JA, JS, JE, &
                              DENS(:,:,:), PV_info(I_DENS)%NAME,    & ! (in)
                              ATMOS_GRID_CARTESC_REAL_VOL  (:,:,:), & ! (in)
                              ATMOS_GRID_CARTESC_REAL_TOTVOL        ) ! (in)
       call STATISTICS_total( KA, KS, KE, IA, IS, IE, JA, JS, JE, &
                             MOMZ(:,:,:), PV_info(I_MOMZ)%NAME,     & ! (in)
                             ATMOS_GRID_CARTESC_REAL_VOLWXY(:,:,:), & ! (in)
                             ATMOS_GRID_CARTESC_REAL_TOTVOLWXY      ) ! (in)
       call STATISTICS_total( KA, KS, KE, IA, IS, IE, JA, JS, JE, &
                              MOMX(:,:,:), PV_info(I_MOMX)%NAME,     & ! (in)
                              ATMOS_GRID_CARTESC_REAL_VOLZUY(:,:,:), & ! (in)
                              ATMOS_GRID_CARTESC_REAL_TOTVOLZUY      ) ! (in)
       call STATISTICS_total( KA, KS, KE, IA, IS, IE, JA, JS, JE, &
                              MOMY(:,:,:), PV_info(I_MOMY)%NAME,     & ! (in)
                              ATMOS_GRID_CARTESC_REAL_VOLZXV(:,:,:), & ! (in)
                              ATMOS_GRID_CARTESC_REAL_TOTVOLZXV      ) ! (in)
       call STATISTICS_total( KA, KS, KE, IA, IS, IE, JA, JS, JE, &
                              RHOT(:,:,:), PV_info(I_RHOT)%NAME,    & ! (in)
                              ATMOS_GRID_CARTESC_REAL_VOL  (:,:,:), & ! (in)
                              ATMOS_GRID_CARTESC_REAL_TOTVOL        ) ! (in)

       do iq = 1, QA
          RHOQ(:,:,:) = DENS(:,:,:) * QTRC(:,:,:,iq)

          call STATISTICS_total( KA, KS, KE, IA, IS, IE, JA, JS, JE, &
                                 RHOQ(:,:,:), TRACER_NAME(iq),       & ! (in)
                                 ATMOS_GRID_CARTESC_REAL_VOL(:,:,:), & ! (in)
                                 ATMOS_GRID_CARTESC_REAL_TOTVOL      ) ! (in)
       enddo

       call ATMOS_vars_calc_diagnostics


       RHOQ(KS:KE,IS:IE,JS:JE) = DENS(KS:KE,IS:IE,JS:JE) * QDRY (KS:KE,IS:IE,JS:JE)
       call STATISTICS_total( KA, KS, KE, IA, IS, IE, JA, JS, JE, &
                              RHOQ(:,:,:), 'QDRY',                & ! (in)
                              ATMOS_GRID_CARTESC_REAL_VOL(:,:,:), & ! (in)
                              ATMOS_GRID_CARTESC_REAL_TOTVOL      ) ! (in)

       RHOQ(KS:KE,IS:IE,JS:JE) = DENS(KS:KE,IS:IE,JS:JE) * ( 1.0_RP - QDRY (KS:KE,IS:IE,JS:JE) ) ! Qtotal
       call STATISTICS_total( KA, KS, KE, IA, IS, IE, JA, JS, JE, &
                              RHOQ(:,:,:), 'QTOT',                & ! (in)
                              ATMOS_GRID_CARTESC_REAL_VOL(:,:,:), & ! (in)
                              ATMOS_GRID_CARTESC_REAL_TOTVOL      ) ! (in)


       call ATMOS_vars_get_diagnostic( 'ENGT', WORK3D(:,:,:) )
       call STATISTICS_total( KA, KS, KE, IA, IS, IE, JA, JS, JE, &
                              WORK3D(:,:,:), 'ENGT',              & ! (in)
                              ATMOS_GRID_CARTESC_REAL_VOL(:,:,:), & ! (in)
                              ATMOS_GRID_CARTESC_REAL_TOTVOL      ) ! (in)
       call ATMOS_vars_get_diagnostic( 'ENGP', WORK3D(:,:,:) )
       call STATISTICS_total( KA, KS, KE, IA, IS, IE, JA, JS, JE, &
                              WORK3D(:,:,:), 'ENGP',              & ! (in)
                              ATMOS_GRID_CARTESC_REAL_VOL(:,:,:), & ! (in)
                              ATMOS_GRID_CARTESC_REAL_TOTVOL      ) ! (in)
       call ATMOS_vars_get_diagnostic( 'ENGK', WORK3D(:,:,:) )
       call STATISTICS_total( KA, KS, KE, IA, IS, IE, JA, JS, JE, &
                              WORK3D(:,:,:), 'ENGK',              & ! (in)
                              ATMOS_GRID_CARTESC_REAL_VOL(:,:,:), & ! (in)
                              ATMOS_GRID_CARTESC_REAL_TOTVOL      ) ! (in)
       call ATMOS_vars_get_diagnostic( 'ENGI', WORK3D(:,:,:) )
       call STATISTICS_total( KA, KS, KE, IA, IS, IE, JA, JS, JE, &
                              WORK3D(:,:,:), 'ENGI',              & ! (in)
                              ATMOS_GRID_CARTESC_REAL_VOL(:,:,:), & ! (in)
                              ATMOS_GRID_CARTESC_REAL_TOTVOL      ) ! (in)

    endif

    return
  end subroutine ATMOS_vars_total

  !-----------------------------------------------------------------------------
  !> Calc diagnostic variables
  subroutine ATMOS_vars_calc_diagnostics
    use scale_atmos_grid_cartesC_real, only: &
       REAL_CZ => ATMOS_GRID_CARTESC_REAL_CZ, &
       REAL_FZ => ATMOS_GRID_CARTESC_REAL_FZ
    use scale_atmos_thermodyn, only: &
       ATMOS_THERMODYN_specific_heat
    use scale_atmos_diagnostic_cartesC, only: &
       ATMOS_DIAGNOSTIC_CARTESC_get_vel, &
       ATMOS_DIAGNOSTIC_CARTESC_get_therm, &
       ATMOS_DIAGNOSTIC_CARTESC_get_phyd
    use scale_comm_cartesC, only: &
       COMM_vars8, &
       COMM_wait
    use scale_atmos_hydrometeor, only: &
       N_HYD
    use mod_atmos_phy_mp_vars, only: &
       ATMOS_PHY_MP_vars_get_diagnostic, &
       ATMOS_PHY_MP_vars_reset_diagnostics
    use mod_atmos_phy_ae_vars, only: &
       ATMOS_PHY_AE_vars_reset_diagnostics
    implicit none

    integer :: iq

    call ATMOS_THERMODYN_specific_heat( &
         KA, KS, KE, IA, 1, IA, JA, 1, JA, QA, &
         QTRC_av(:,:,:,:),                                        & ! (in)
         TRACER_MASS(:), TRACER_R(:), TRACER_CV(:), TRACER_CP(:), & ! (in)
         Qdry(:,:,:), Rtot(:,:,:), CVtot(:,:,:), CPtot(:,:,:)     ) ! (out)

    call ATMOS_DIAGNOSTIC_CARTESC_get_vel( &
         KA, KS, KE, IA, 1, IA, JA, 1, JA, &
         DENS_av(:,:,:), MOMZ_av(:,:,:), MOMX_av(:,:,:), MOMY_av(:,:,:), & ! (in)
         W(:,:,:), U(:,:,:), V(:,:,:)                                    ) ! (out)

    call ATMOS_DIAGNOSTIC_CARTESC_get_therm( &
         KA, KS, KE, IA, 1, IA, JA, 1, JA, &
         DENS_av(:,:,:), RHOT_av(:,:,:),                     & ! (in)
         Rtot(:,:,:), CVtot(:,:,:), CPtot(:,:,:),            & ! (in)
         POTT(:,:,:), TEMP(:,:,:), PRES(:,:,:), EXNER(:,:,:) ) ! (out)

    call ATMOS_DIAGNOSTIC_CARTESC_get_phyd( &
         KA, KS, KE, IA, 1, IA, JA, 1, JA, &
         DENS_av(:,:,:), PRES(:,:,:),    & ! (in)
         REAL_CZ(:,:,:), REAL_FZ(:,:,:), & ! (in)
         PHYD(:,:,:), PHYDH(:,:,:)       ) ! (out)

    call ATMOS_PHY_MP_vars_reset_diagnostics
    call ATMOS_PHY_AE_vars_reset_diagnostics

    if ( moist ) then
       call ATMOS_PHY_MP_vars_get_diagnostic( &
            DENS_av(:,:,:), TEMP(:,:,:), QTRC_av(:,:,:,:), & ! (in)
            Qe=Qe(:,:,:,:)                                 ) ! (out)
       do iq = 1, N_HYD
          call COMM_vars8(Qe(:,:,:,iq), iq)
       end do
       do iq = 1, N_HYD
          call COMM_wait (Qe(:,:,:,iq), iq)
       end do
    end if

    ! reset diagnostic variables
    DV_calclated(:) = .false.

    return
  end subroutine ATMOS_vars_calc_diagnostics

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
    use scale_atmos_grid_cartesC, only: &
       RCDX => ATMOS_GRID_CARTESC_RCDX, &
       RCDY => ATMOS_GRID_CARTESC_RCDY
    use scale_atmos_grid_cartesC_real, only: &
       REAL_CZ => ATMOS_GRID_CARTESC_REAL_CZ, &
       REAL_FZ => ATMOS_GRID_CARTESC_REAL_FZ
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
    use scale_atmos_diagnostic_cartesC, only: &
       ATMOS_DIAGNOSTIC_CARTESC_get_potv, &
       ATMOS_DIAGNOSTIC_CARTESC_get_teml, &
       ATMOS_DIAGNOSTIC_CARTESC_get_n2
    implicit none
    character(len=*), intent(in)  :: vname
    real(RP),         intent(out) :: var(:,:,:)

    real(RP) :: UH  (KA,IA,JA)
    real(RP) :: VH  (KA,IA,JA)

    real(RP) :: WORK(KA,IA,JA)

    integer :: k, i, j, iq

    select case ( vname )
    case ( 'W' )
       var(:,:,:) = W(:,:,:)

    case ( 'U' )
       var(:,:,:) = U(:,:,:)

    case ( 'V' )
       var(:,:,:) = V(:,:,:)

    case ( 'PT' )
       var(:,:,:) = POTT(:,:,:)

    case ( 'T' )
       var(:,:,:) = TEMP(:,:,:)

    case ( 'EXNER' )
       var(:,:,:) = EXNER(:,:,:)

    case ( 'PHYD' )
       var(:,:,:) = PHYD(:,:,:)

    case ( 'QDRY' )
       var(:,:,:) = QDRY(:,:,:)

    case ( 'RTOT' )
       var(:,:,:) = RTOT(:,:,:)

    case ( 'CVTOT' )
       var(:,:,:) = CVTOT(:,:,:)

    case ( 'CPTOT' )
       var(:,:,:) = CPTOT(:,:,:)

    case ( 'LHV' )
       if ( .not. DV_calclated(I_LHV) ) then
          call allocate_3D( LHV )
          call ATMOS_HYDROMETEOR_LHV( &
               KA, KS, KE, IA, 1, IA, JA, 1, JA, &
               TEMP(:,:,:), & ! (in)
               LHV(:,:,:)   ) ! (out)
          DV_calclated(I_LHV) = .true.
       end if
       var(:,:,:) = LHV(:,:,:)

    case ( 'LHS' )
       if ( .not. DV_calclated(I_LHS) ) then
          call allocate_3D( LHS )
          call ATMOS_HYDROMETEOR_LHS( &
               KA, KS, KE, IA, 1, IA, JA, 1, JA, &
               TEMP(:,:,:), & ! (in)
               LHS(:,:,:)   ) ! (out)
          DV_calclated(I_LHS) = .true.
       end if
       var(:,:,:) = LHS(:,:,:)

    case ( 'LHF' )
       if ( .not. DV_calclated(I_LHF) ) then
          call allocate_3D( LHF )
          call ATMOS_HYDROMETEOR_LHF( &
               KA, KS, KE, IA, 1, IA, JA, 1, JA, &
               TEMP(:,:,:), & ! (in)
               LHF(:,:,:)   ) ! (out)
          DV_calclated(I_LHF) = .true.
       end if
       var(:,:,:) = LHF(:,:,:)

    case ( 'POTV' )
       if ( .not. DV_calclated(I_POTV) ) then
          call allocate_3D( POTV )
          call ATMOS_DIAGNOSTIC_CARTESC_get_potv( &
               KA, KS, KE, IA, 1, IA, JA, 1, JA, &
               POTT(:,:,:), Rtot(:,:,:), & ! (in)
               POTV(:,:,:)               ) ! (out)
          DV_calclated(I_POTV) = .true.
       end if
       var(:,:,:) = POTV(:,:,:)

    case ( 'TEML' )
       if ( .not. DV_calclated(I_TEML) ) then
          call allocate_3D( TEML )
          call ATMOS_vars_get_diagnostic( 'LHV', WORK3D(:,:,:) )
          call ATMOS_vars_get_diagnostic( 'LHS', WORK3D(:,:,:) )
          call ATMOS_vars_get_diagnostic( 'QLIQ', WORK3D(:,:,:) )
          call ATMOS_vars_get_diagnostic( 'QICE', WORK3D(:,:,:) )
          call ATMOS_DIAGNOSTIC_CARTESC_get_teml( &
               KA, KS, KE, IA, 1, IA, JA, 1, JA, &
               TEMP(:,:,:), LHV(:,:,:), LHS(:,:,:), & ! (in)
               QC(:,:,:), QI(:,:,:), CPtot(:,:,:),  & ! (in)
               TEML(:,:,:)                          ) ! (out)
          DV_calclated(I_TEML) = .true.
       end if
       var(:,:,:) = TEML(:,:,:)

    case ( 'POTL' )
       if ( .not. DV_calclated(I_POTL) ) then
          call allocate_3D( POTL )
          call ATMOS_vars_get_diagnostic( 'TEML', WORK3D(:,:,:) )
!OCL XFILL
          !$omp parallel do default(none) OMP_SCHEDULE_ collapse(2) &
          !$omp private(i,j,k) &
          !$omp shared(POTL,TEML,EXNER) &
          !$omp shared(KS,KE,IA,JA)
          do j = 1, JA
          do i = 1, IA
          do k = KS, KE
             POTL(k,i,j) = TEML(k,i,j) / EXNER(k,i,j)
          enddo
          enddo
          enddo
          DV_calclated(I_POTL) = .true.
       end if
       var(:,:,:) = POTL(:,:,:)

    case ( 'POTE' )
       if ( .not. DV_calclated(I_POTE) ) then
          call allocate_3D( POTE )
          call ATMOS_SATURATION_pote( &
               KA, KS, KE, IA, 1, IA, JA, 1, JA, &
               DENS(:,:,:), POTT(:,:,:), TEMP(:,:,:), QV(:,:,:), & ! [IN]
               POTE(:,:,:)                                       ) ! [OUT]
       end if
       var(:,:,:) = POTE(:,:,:)
    case ( 'QTOT' )
       if ( .not. DV_calclated(I_QTOT) ) then
          call allocate_3D( QTOT )
          if ( moist ) then
             call ATMOS_vars_get_diagnostic( 'QHYD', WORK3D(:,:,:) )
!OCL XFILL
             !$omp parallel do default(none) OMP_SCHEDULE_ collapse(2) &
             !$omp private(i,j,k) &
             !$omp shared(QTOT,QV,QHYD) &
             !$omp shared(KS,KE,IA,JA)
             do j = 1, JA
             do i = 1, IA
             do k = KS, KE
                QTOT(k,i,j) = QV(k,i,j) + QHYD(k,i,j)
             enddo
             enddo
             enddo
          else
!OCL XFILL
             !$omp parallel do default(none) OMP_SCHEDULE_ collapse(2) &
             !$omp private(i,j,k) &
             !$omp shared(QTOT) &
             !$omp shared(KS,KE,IA,JA)
             do j = 1, JA
             do i = 1, IA
             do k = KS, KE
                QTOT(k,i,j) = 0.0_RP
             enddo
             enddo
             enddo
          end if
          DV_calclated(I_QTOT) = .true.
       end if
       var(:,:,:) = QTOT(:,:,:)

    case ( 'QHYD' )
       if ( .not. DV_calclated(I_QHYD) ) then
          call allocate_3D( QHYD )
          if ( moist ) then
             call ATMOS_vars_get_diagnostic( 'QLIQ', WORK3D(:,:,:) )
             call ATMOS_vars_get_diagnostic( 'QICE', WORK3D(:,:,:) )
!OCL XFILL
             !$omp parallel do default(none) OMP_SCHEDULE_ collapse(2) &
             !$omp private(i,j,k) &
             !$omp shared(QHYD,QLIQ,QICE) &
             !$omp shared(KS,KE,IA,JA)
             do j = 1, JA
             do i = 1, IA
             do k = KS, KE
                QHYD(k,i,j) = QLIQ(k,i,j) + QICE(k,i,j)
             enddo
             enddo
             enddo
          else
!OCL XFILL
             !$omp parallel do default(none) OMP_SCHEDULE_ collapse(2) &
             !$omp private(i,j,k) &
             !$omp shared(QHYD) &
             !$omp shared(KS,KE,IA,JA)
             do j = 1, JA
             do i = 1, IA
             do k = KS, KE
                QHYD(k,i,j) = 0.0_RP
             enddo
             enddo
             enddo
          end if
          DV_calclated(I_QHYD) = .true.
       end if
       var(:,:,:) = QHYD(:,:,:)

    case ( 'QLIQ' )
       if ( .not. DV_calclated(I_QLIQ) ) then
          call allocate_3D( QLIQ )
          !$omp parallel do default(none) OMP_SCHEDULE_ collapse(2) &
          !$omp private(i,j,k,iq) &
          !$omp shared(QLIQ,QC,QR) &
          !$omp shared(KS,KE,IA,JA)
          do j = 1, JA
          do i = 1, IA
          do k = KS, KE
!OCL XFILL
             QLIQ(k,i,j) = QC(k,i,j) + QR(k,i,j)
          enddo
          enddo
          enddo
          DV_calclated(I_QLIQ) = .true.
       end if
       var(:,:,:) = QLIQ(:,:,:)

    case ( 'QICE' )
       if ( .not. DV_calclated(I_QICE) ) then
          call allocate_3D( QICE )
!OCL XFILL
          !$omp parallel do default(none) OMP_SCHEDULE_ collapse(2) &
          !$omp private(i,j,k,iq) &
          !$omp shared(QICE,QI,QS,QG,QH) &
          !$omp shared(KS,KE,IA,JA)
          do j = 1, JA
          do i = 1, IA
          do k = KS, KE
             QICE(k,i,j) = QI(k,i,j) + QS(k,i,j) + QG(k,i,j) + QH(k,i,j)
          enddo
          enddo
          enddo
          DV_calclated(I_QICE) = .true.
       end if
       var(:,:,:) = QICE(:,:,:)

    case ( 'QSAT' )
       if ( .not. DV_calclated(I_QSAT) ) then
          call allocate_3D( QSAT )
          call ATMOS_SATURATION_dens2qsat_all( &
               KA, KS, KE, IA, 1, IA, JA, 1, JA, &
               TEMP(:,:,:), DENS_av(:,:,:), & ! (in)
               QSAT(:,:,:)                  ) ! (out)
          DV_calclated(I_QSAT) = .true.
       end if
       var(:,:,:) = QSAT(:,:,:)

    case ( 'RHA' )
       if ( .not. DV_calclated(I_RHA) ) then
          call allocate_3D( RHA )
          if ( moist ) then
             call ATMOS_SATURATION_psat_all( &
               KA, KS, KE, IA, 1, IA, JA, 1, JA, &
               TEMP(:,:,:), & ! (in)
               WORK(:,:,:)  ) ! (out)
!OCL XFILL
             !$omp parallel do default(none) OMP_SCHEDULE_ collapse(2) &
             !$omp private(i,j,k) &
             !$omp shared(RHA,DENS_av,QV,WORK,TEMP) &
             !$omp shared(KS,KE,IA,JA)
             do j = 1, JA
             do i = 1, IA
             do k = KS, KE
                RHA(k,i,j) = DENS_av(k,i,j) * QV(k,i,j) &
                           / WORK(k,i,j) * Rvap * TEMP(k,i,j) &
                           * 100.0_RP
             enddo
             enddo
             enddo
          else
!OCL XFILL
             !$omp parallel do default(none) OMP_SCHEDULE_ collapse(2) &
             !$omp private(i,j,k) &
             !$omp shared(RHA) &
             !$omp shared(KS,KE,IA,JA)
             do j = 1, JA
             do i = 1, IA
             do k = KS, KE
                RHA(k,i,j) = 0.0_RP
             enddo
             enddo
             enddo
          end if
          DV_calclated(I_RHA) = .true.
       end if
       var(:,:,:) = RHA(:,:,:)

    case ( 'RHL', 'RH' )
       if ( .not. DV_calclated(I_RHL) ) then
          call allocate_3D( RHL )
          if ( moist ) then
             call ATMOS_SATURATION_psat_liq( &
                  KA, KS, KE, IA, 1, IA, JA, 1, JA, &
                  TEMP(:,:,:), & ! (in)
                  WORK(:,:,:)  ) ! (out)
!OCL XFILL
             !$omp parallel do default(none) OMP_SCHEDULE_ collapse(2) &
             !$omp private(i,j,k) &
             !$omp shared(RHL,DENS_av,QV,WORK,TEMP) &
             !$omp shared(KS,KE,IA,JA)
             do j = 1, JA
             do i = 1, IA
             do k = KS, KE
                RHL(k,i,j) = DENS_av(k,i,j) * QV(k,i,j) &
                           / WORK(k,i,j) * Rvap * TEMP(k,i,j) &
                           * 100.0_RP
             enddo
             enddo
             enddo
          else
!OCL XFILL
             !$omp parallel do default(none) OMP_SCHEDULE_ collapse(2) &
             !$omp private(i,j,k) &
             !$omp shared(RHL) &
             !$omp shared(KS,KE,IA,JA)
             do j = 1, JA
             do i = 1, IA
             do k = KS, KE
                RHL(k,i,j) = 0.0_RP
             enddo
             enddo
             enddo
          end if
          DV_calclated(I_RHL) = .true.
       end if
       var(:,:,:) = RHL(:,:,:)

    case ( 'RHI' )
       if ( .not. DV_calclated(I_RHI) ) then
          call allocate_3D( RHI )
          if ( moist ) then
             call ATMOS_SATURATION_psat_ice( &
                  KA, KS, KE, IA, 1, IA, JA, 1, JA, &
                  TEMP(:,:,:), & ! (int)
                  WORK(:,:,:)  ) ! (out)
!OCL XFILL
             !$omp parallel do default(none) OMP_SCHEDULE_ collapse(2) &
             !$omp private(i,j,k) &
             !$omp shared(RHI,DENS_av,QV,WORK,TEMP) &
             !$omp shared(KS,KE,IA,JA)
             do j = 1, JA
             do i = 1, IA
             do k = KS, KE
                RHI(k,i,j) = DENS_av(k,i,j) * QV(k,i,j) &
                           / WORK(k,i,j) * Rvap * TEMP(k,i,j) &
                           * 100.0_RP
             enddo
             enddo
             enddo
          else
!OCL XFILL
             !$omp parallel do default(none) OMP_SCHEDULE_ collapse(2) &
             !$omp private(i,j,k) &
             !$omp shared(RHI) &
             !$omp shared(KS,KE,IA,JA)
             do j = 1, JA
             do i = 1, IA
             do k = KS, KE
                RHI(k,i,j) = 0.0_RP
             enddo
             enddo
             enddo
          end if
          DV_calclated(I_RHI) = .true.
       end if
       var(:,:,:) = RHI(:,:,:)

    case ( 'VOR' )
       if ( .not. DV_calclated(I_VOR) ) then
          call allocate_3D( VOR )
          !!!  to move to grid !!!
          ! at x, v, layer
          !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
!OCL XFILL
          do j = 1, JA-1
          do i = 2, IA
          do k = KS, KE
             UH(k,i,j) = 0.5_RP * ( MOMX_av(k,i  ,j) + MOMX_av(k,i  ,j+1) &
                                  + MOMX_av(k,i-1,j) + MOMX_av(k,i-1,j+1) ) &
                       / ( DENS_av(k,i,j) + DENS_av(k,i,j+1) )
          enddo
          enddo
          enddo
          ! at u, y, layer
          !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
!OCL XFILL
          do j = 2, JA
          do i = 1, IA-1
          do k = KS, KE
             VH(k,i,j) = 0.5_RP * ( MOMY_av(k,i,j  ) + MOMY_av(k,i+1,j  ) &
                                  + MOMY_av(k,i,j-1) + MOMY_av(k,i+1,j-1) ) &
                       / ( DENS_av(k,i,j) + DENS_av(k,i+1,j) )
          enddo
          enddo
          enddo
          !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
!OCL XFILL
          do j = 2, JA-1
          do i = 2, IA-1
          do k = KS, KE
             VOR(k,i,j) = ( VH(k,i,j  ) - VH(k,i-1,j  ) ) * RCDX(i) &
                        - ( UH(k,i  ,j) - UH(k,i  ,j-1) ) * RCDY(j)
          enddo
          enddo
          enddo
          !$omp parallel do private(j,k) OMP_SCHEDULE_
          do j = 1, JA
          do k = KS, KE
             VOR(k,1 ,j) = VOR(k,2   ,j)
             VOR(k,IA,j) = VOR(k,IA-1,j)
          enddo
          enddo
          !$omp parallel do private(i,k) OMP_SCHEDULE_
          do i = 1, IA
          do k = KS, KE
             VOR(k,i,1 ) = VOR(k,i,2   )
             VOR(k,i,JA) = VOR(k,i,JA-1)
          enddo
          enddo
          DV_calclated(I_VOR) = .true.
       end if
       var(:,:,:) = VOR(:,:,:)

    case ( 'DIV' )
       if ( .not. DV_calclated(I_DIV) ) then
          call allocate_3D( DIV )
          call ATMOS_vars_get_diagnostic( 'HDIV', WORK3D(:,:,:) )
          !!!! to move to grid !!!!
          !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
!OCL XFILL
          do j = 1, JA
          do i = 1, IA
          do k = KS, KE
             DIV(k,i,j) = ( MOMZ_av(k,i,j) - MOMZ_av(k-1,i  ,j  ) ) * ( REAL_FZ(k,i,j)-REAL_FZ(k-1,i,j) ) &
                        + HDIV(k,i,j)
          enddo
          enddo
          enddo
          DV_calclated(I_DIV) = .true.
       end if
       var(:,:,:) = DIV(:,:,:)

    case ( 'HDIV' )
       if ( .not. DV_calclated(I_HDIV) ) then
          call allocate_3D( HDIV )
          !!!! to move to grid !!!!
          !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
!OCL XFILL
          do j = 2, JA
          do i = 2, IA
          do k = KS, KE
             HDIV(k,i,j) = ( MOMX_av(k,i,j) - MOMX_av(k  ,i-1,j  ) ) * RCDX(i) &
                         + ( MOMY_av(k,i,j) - MOMY_av(k  ,i  ,j-1) ) * RCDY(j)
          enddo
          enddo
          enddo
          !$omp parallel do private(i,k) OMP_SCHEDULE_
          do i = 1, IA
          do k = KS, KE
             HDIV(k,i,1) = HDIV(k,i,2)
          enddo
          enddo
          !$omp parallel do private(j,k) OMP_SCHEDULE_
          do j = 1, JA
          do k = KS, KE
             HDIV(k,1,j) = HDIV(k,2,j)
          enddo
          enddo
          DV_calclated(I_HDIV) = .true.
       end if
       var(:,:,:) = HDIV(:,:,:)

    case ( 'Uabs' )
       if ( .not. DV_calclated(I_Uabs) ) then
          call allocate_3D( Uabs )
!OCL XFILL
          !$omp parallel do private(k,i,j) OMP_SCHEDULE_ collapse(2)
          do j = 1, JA
          do i = 1, IA
          do k = KS, KE
             Uabs(k,i,j) = sqrt( W(k,i,j)**2 + U(k,i,j)**2 + V(k,i,j)**2 )
          enddo
          enddo
          enddo
          DV_calclated(I_Uabs) = .true.
       end if
       var(:,:,:) = Uabs(:,:,:)

    case ( 'N2' )
       if ( .not. DV_calclated(I_N2) ) then
          call allocate_3D( N2 )
          call ATMOS_DIAGNOSTIC_CARTESC_get_n2( &
               KA, KS, KE, IA, 1, IA, JA, 1, JA, &
               POTT(:,:,:), Rtot(:,:,:), & !(in)
               REAL_CZ(:,:,:),           & !(in)
               N2(:,:,:)                 ) ! (out)
          DV_calclated(I_N2) = .true.
       end if
       var(:,:,:) = N2(:,:,:)

    case ( 'MSE' )
       if ( .not. DV_calclated(I_MSE) ) then
          call allocate_3D( MSE )
          call ATMOS_vars_get_diagnostic( 'LHV', WORK3D(:,:,:) )
!OCL XFILL
          !$omp parallel do private(k,i,j) OMP_SCHEDULE_ collapse(2)
          do j = 1, JA
          do i = 1, IA
          do k = KS, KE
             MSE(k,i,j) = CPTOT(k,i,j) * TEMP(k,i,j)                    &
                        + GRAV * ( REAL_CZ(k,i,j) - REAL_FZ(KS-1,i,j) ) &
                        + LHV(k,i,j) * QV(k,i,j)
          enddo
          enddo
          enddo
          DV_calclated(I_MSE) = .true.
       end if
       var(:,:,:) = MSE(:,:,:)

    case ( 'TDEW' )
       if ( .not. DV_calclated(I_TDEW) ) then
          call allocate_3D( TDEW )
          call ATMOS_SATURATION_tdew_liq( KA, KS, KE, IA, 1, IA, JA, 1, JA, &
                                          DENS(:,:,:), TEMP(:,:,:), QV(:,:,:), & ! [IN]
                                          TDEW(:,:,:)                          ) ! [OUT]
          DV_calclated(I_TDEW) = .true.
       end if
       var(:,:,:) = TDEW(:,:,:)

    case ( 'ENGP' )
       if ( .not. DV_calclated(I_ENGP) ) then
          call allocate_3D( ENGP )
          !$omp parallel do private(k,i,j) OMP_SCHEDULE_ collapse(2)
          do j = 1, JA
          do i = 1, IA
          do k = KS, KE
             ENGP(k,i,j) = DENS_av(k,i,j) * GRAV * REAL_CZ(k,i,j)
          end do
          end do
          end do
          DV_calclated(I_ENGP) = .true.
       end if
       var(:,:,:) = ENGP(:,:,:)

    case ( 'ENGK' )
       if ( .not. DV_calclated(I_ENGK) ) then
          call allocate_3D( ENGK )
          !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
          do j = 1, JA
          do i = 1, IA
          do k = KS, KE
             ENGK(k,i,j) = 0.5_RP * DENS_av(k,i,j) &
                         * ( W(k,i,j)**2 + U(k,i,j)**2 + V(k,i,j)**2 )
          end do
          end do
          end do
             DV_calclated(I_ENGK) = .true.
       end if
       var(:,:,:) = ENGK(:,:,:)

    case ( 'ENGI' )
       if ( .not. DV_calclated(I_ENGI) ) then
          call allocate_3D( ENGI )
          if ( moist ) then
             call ATMOS_vars_get_diagnostic( 'QICE', WORK3D(:,:,:) )
          end if
          !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
          do j = 1, JA
          do i = 1, IA
          do k = KS, KE
             ENGI(k,i,j) = DENS_av(k,i,j) * QDRY(k,i,j) * TEMP(k,i,j) * CVdry
             do iq = 1, QA
                ENGI(k,i,j) = ENGI(k,i,j) &
                            + DENS_av(k,i,j) * QTRC_av(k,i,j,iq) * TEMP(k,i,j) * TRACER_CV(iq)
             enddo
             if ( moist ) then
                ENGI(k,i,j) = ENGI(k,i,j) &
                     + DENS_av(k,i,j) * ( QV  (k,i,j) * LHVc & ! Latent Heat [vapor->liquid]
                                        - QICE(k,i,j) * LHFc ) ! Latent Heat [ice->liquid]
             end if
          end do
          end do
          end do
          DV_calclated(I_ENGI) = .true.
       end if
       var(:,:,:) = ENGI(:,:,:)

    case ( 'ENGT' )
       if ( .not. DV_calclated(I_ENGT) ) then
          call allocate_3D( ENGT )
          call ATMOS_vars_get_diagnostic( 'ENGP', WORK3D(:,:,:) )
          call ATMOS_vars_get_diagnostic( 'ENGK', WORK3D(:,:,:) )
          call ATMOS_vars_get_diagnostic( 'ENGI', WORK3D(:,:,:) )
          !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
          do j = 1, JA
          do i = 1, IA
          do k = KS, KE
             ENGT(k,i,j) = ENGP(k,i,j) + ENGK(k,i,j) + ENGI(k,i,j)
          enddo
          enddo
          enddo
          DV_calclated(I_ENGT) = .true.
       end if
       var(:,:,:) = ENGT(:,:,:)

    case ( 'DENS_PRIM' )
       if ( .not. DV_calclated(I_DENS_PRIM) ) then
          call allocate_3D( DENS_PRIM )
          call ATMOS_vars_get_diagnostic( 'DENS_MEAN', WORK1D(:) )
!OCL XFILL
          !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
          do j = 1, JA
          do i = 1, IA
          do k = KS, KE
             DENS_PRIM(k,i,j) = DENS_av(k,i,j) - DENS_MEAN(k)
          enddo
          enddo
          enddo
          DV_calclated(I_DENS_PRIM) = .true.
       end if
       var(:,:,:) = DENS_PRIM(:,:,:)

    case ( 'W_PRIM' )
       if ( .not. DV_calclated(I_W_PRIM) ) then
          call allocate_3D( W_PRIM )
          call ATMOS_vars_get_diagnostic( 'W_MEAN', WORK1D(:) )
!OCL XFILL
          !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
          do j = 1, JA
          do i = 1, IA
          do k = KS, KE
             W_PRIM(k,i,j) = W(k,i,j) - W_MEAN(k)
          enddo
          enddo
          enddo
          DV_calclated(I_W_PRIM) = .true.
       end if
       var(:,:,:) = W_PRIM(:,:,:)

    case ( 'U_PRIM' )
       if ( .not. DV_calclated(I_U_PRIM) ) then
          call allocate_3D( U_PRIM )
          call ATMOS_vars_get_diagnostic( 'U_MEAN', WORK1D(:) )
!OCL XFILL
          !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
          do j = 1, JA
          do i = 1, IA
          do k = KS, KE
             U_PRIM(k,i,j) = U(k,i,j) - U_MEAN(k)
          enddo
          enddo
          enddo
          DV_calclated(I_U_PRIM) = .true.
       end if
       var(:,:,:) = U_PRIM(:,:,:)

    case ( 'V_PRIM' )
       if ( .not. DV_calclated(I_V_PRIM) ) then
          call allocate_3D( V_PRIM )
          call ATMOS_vars_get_diagnostic( 'V_MEAN', WORK1D(:) )
!OCL XFILL
          !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
          do j = 1, JA
          do i = 1, IA
          do k = KS, KE
             V_PRIM(k,i,j) = V(k,i,j) - V_MEAN(k)
          enddo
          enddo
          enddo
          DV_calclated(I_V_PRIM) = .true.
       end if
       var(:,:,:) = V_PRIM(:,:,:)

    case ( 'PT_PRIM' )
       if ( .not. DV_calclated(I_PT_PRIM) ) then
          call allocate_3D( PT_PRIM )
          call ATMOS_vars_get_diagnostic( 'PT_MEAN', WORK1D(:) )
!OCL XFILL
          !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
          do j = 1, JA
          do i = 1, IA
          do k = KS, KE
             PT_PRIM(k,i,j) = POTT(k,i,j) - PT_MEAN(k)
          enddo
          enddo
          enddo
          DV_calclated(I_PT_PRIM) = .true.
       end if
       var(:,:,:) = PT_PRIM(:,:,:)

    case ( 'W_PRIM2' )
       if ( .not. DV_calclated(I_W_PRIM2) ) then
          call allocate_3D( W_PRIM2 )
          call ATMOS_vars_get_diagnostic( 'W_PRIM', WORK3D(:,:,:) )
!OCL XFILL
          !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
          do j = 1, JA
          do i = 1, IA
          do k = KS, KE
             W_PRIM2(k,i,j) = W_PRIM(k,i,j)**2
          enddo
          enddo
          enddo
          DV_calclated(I_W_PRIM2) = .true.
       end if
       var(:,:,:) = W_PRIM2(:,:,:)

    case ( 'PT_W_PRIM' )
       if ( .not. DV_calclated(I_PT_W_PRIM) ) then
          call allocate_3D( PT_W_PRIM )
          call ATMOS_vars_get_diagnostic( 'W_PRIM',  WORK3D(:,:,:) )
          call ATMOS_vars_get_diagnostic( 'PT_PRIM', WORK3D(:,:,:) )
          !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
          do j = 1, JA
          do i = 1, IA
          do k = KS, KE
             PT_W_PRIM(k,i,j) = W_PRIM(k,i,j) * PT_PRIM(k,i,j) * DENS_av(k,i,j) * CPdry
          enddo
          enddo
          enddo
          DV_calclated(I_PT_W_PRIM) = .true.
       end if
       var(:,:,:) = PT_W_PRIM(:,:,:)

    case ( 'W_PRIM3' )
       if ( .not. DV_calclated(I_W_PRIM3) ) then
          call allocate_3D( W_PRIM3 )
          call ATMOS_vars_get_diagnostic( 'W_PRIM', WORK3D(:,:,:) )
!OCL XFILL
          !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
          do j = 1, JA
          do i = 1, IA
          do k = KS, KE
             W_PRIM3(k,i,j) = W_PRIM(k,i,j)**3
          enddo
          enddo
          enddo
          DV_calclated(I_W_PRIM3) = .true.
       end if
       var(:,:,:) = W_PRIM3(:,:,:)

    case ( 'TKE_RS' )
       if ( .not. DV_calclated(I_TKE_RS) ) then
          call allocate_3D( TKE_RS )
          call ATMOS_vars_get_diagnostic( 'W_PRIM', WORK3D(:,:,:) )
          call ATMOS_vars_get_diagnostic( 'U_PRIM', WORK3D(:,:,:) )
          call ATMOS_vars_get_diagnostic( 'V_PRIM', WORK3D(:,:,:) )
!OCL XFILL
          !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
          do j = 1, JA
          do i = 1, IA
          do k = KS, KE
             TKE_RS(k,i,j) = 0.5_RP * ( W_PRIM(k,i,j)**2 + U_PRIM(k,i,j)**2 + V_PRIM(k,i,j)**2 )
          enddo
          enddo
          enddo
          DV_calclated(I_TKE_RS) = .true.
       end if
       var(:,:,:) = TKE_RS(:,:,:)

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
    use scale_atmos_grid_cartesC_real, only: &
       REAL_CZ => ATMOS_GRID_CARTESC_REAL_CZ, &
       REAL_FZ => ATMOS_GRID_CARTESC_REAL_FZ
    use scale_atmos_adiabat, only: &
       ATMOS_ADIABAT_cape
    use mod_atmos_phy_mp_vars, only: &
       SFLX_rain_MP => ATMOS_PHY_MP_SFLX_rain, &
       SFLX_snow_MP => ATMOS_PHY_MP_SFLX_snow
    use mod_atmos_phy_cp_vars, only: &
       SFLX_rain_CP => ATMOS_PHY_CP_SFLX_rain
    implicit none
    character(len=*), intent(in)  :: vname
    real(RP),         intent(out) :: var(:,:)

    real(RP) :: fact

    integer :: k, i, j, iq

    select case ( vname )
    case ( 'LWP' )
       if ( .not. DV_calclated(I_LWP) ) then
          call allocate_2D( LWP )
          call ATMOS_vars_get_diagnostic( 'QLIQ', WORK3D(:,:,:) )
          !$omp parallel do default(none) OMP_SCHEDULE_ collapse(2) &
          !$omp private(i,j,k) &
          !$omp shared(LWP,QLIQ,DENS_av,REAL_FZ) &
          !$omp shared(KS,KE,IA,JA)
          do j = 1, JA
          do i = 1, IA
             LWP(i,j) = 0.0_RP
             do k  = KS, KE
                LWP(i,j) = LWP(i,j) &
                         + QLIQ(k,i,j) * DENS_av(k,i,j) * ( REAL_FZ(k,i,j)-REAL_FZ(k-1,i,j) ) * 1.E3_RP ! [kg/m2->g/m2]
             enddo
          enddo
          enddo
          DV_calclated(I_LWP) = .true.
       end if
       var(:,:) = LWP(:,:)

    case ( 'IWP' )
       if ( .not. DV_calclated(I_IWP) ) then
          call allocate_2D( IWP )
          call ATMOS_vars_get_diagnostic( 'QICE', WORK3D(:,:,:) )
          !$omp parallel do default(none) OMP_SCHEDULE_ collapse(2) &
          !$omp private(i,j,k) &
          !$omp shared(IWP,QICE,DENS_av,REAL_FZ) &
          !$omp shared(KS,KE,IA,JA)
          do j = 1, JA
          do i = 1, IA
             IWP(i,j) = 0.0_RP
             do k  = KS, KE
                IWP(i,j) = IWP(i,j) &
                         + QICE(k,i,j) * DENS_av(k,i,j) * ( REAL_FZ(k,i,j)-REAL_FZ(k-1,i,j) ) * 1.E3_RP ! [kg/m2->g/m2]
             enddo
          enddo
          enddo
          DV_calclated(I_IWP) = .true.
       end if
       var(:,:) = IWP(:,:)

    case ( 'PW' )
       if ( .not. DV_calclated(I_PW) ) then
          call allocate_2D( PW )
          !$omp parallel do default(none) OMP_SCHEDULE_ collapse(2) &
          !$omp private(i,j,k) &
          !$omp shared(PW,QV,DENS_av,REAL_FZ) &
          !$omp shared(KS,KE,IA,JA)
          do j = 1, JA
          do i = 1, IA
             PW(i,j) = 0.0_RP
             do k  = KS, KE
                PW(i,j) = PW(i,j) &
                        + QV(k,i,j) * DENS_av(k,i,j) * ( REAL_FZ(k,i,j)-REAL_FZ(k-1,i,j) ) * 1.E3_RP ! [kg/m2->g/m2]
             enddo
          enddo
          enddo
          DV_calclated(I_PW) = .true.
       end if
       var(:,:) = PW(:,:)

    case ( 'PBLH' )
       if ( .not. DV_calclated(I_PBLH) ) then
          call allocate_2D( PBLH )
          call ATMOS_vars_get_diagnostic( 'POTV', WORK3D(:,:,:) )
          !$omp parallel do default(none) OMP_SCHEDULE_ collapse(2) &
          !$omp private(k,i,j) &
          !$omp private(fact) &
          !$omp shared(PBLH,POTV,REAL_CZ,REAL_FZ) &
          !$omp shared(KS,KE,IA,JA)
          do j = 1, JA
          do i = 1, IA
             PBLH(i,j) = REAL_CZ(KS,i,j) - REAL_FZ(KS-1,i,j)
             do k = KS+1, KE
                if ( POTV(k,i,j) > POTV(KS,i,j) ) then
                   fact = ( POTV(KS,i,j) - POTV(k-1,i,j) ) &
                        / ( POTV(k,i,j)  - POTV(k-1,i,j) )
                   PBLH(i,j) = REAL_CZ(k-1,i,j) - REAL_FZ(KS-1,i,j) &
                             + fact * ( REAL_CZ(k,i,j) - REAL_CZ(k-1,i,j) )

                   exit
                endif
             enddo
          enddo
          enddo
          DV_calclated(I_PBLH) = .true.
       end if
       var(:,:) = PBLH(:,:)

    case ( 'CAPE', 'CIN', 'LCL', 'LFC', 'LNB' )
       if ( .not. DV_calclated(I_CAPE) ) then
          call allocate_2D( CAPE )
          call allocate_2D( CIN )
          call allocate_2D( LCL )
          call allocate_2D( LFC )
          call allocate_2D( LNB )
          call ATMOS_ADIABAT_cape( &
               KA, KS, KE, IA, IS, IE, JA, JS, JE, &
               KS,                                               & ! (in)
               DENS_av(:,:,:), TEMP(:,:,:), PRES(:,:,:),         & ! (in)
               QV(:,:,:), QC(:,:,:), Qdry(:,:,:),                & ! (in)
               Rtot(:,:,:), CPtot(:,:,:),                        & ! (in)
               REAL_CZ(:,:,:), REAL_FZ(:,:,:),                   & ! (in)
               CAPE(:,:), CIN(:,:), LCL(:,:), LFC(:,:), LNB(:,:) ) ! (out)
          DV_calclated(I_CAPE) = .true.
       end if
       select case ( vname )
       case ( 'CAPE' )
          var(:,:) = CAPE(:,:)
       case ( 'CIN' )
          var(:,:) = CIN(:,:)
       case ( 'LCL' )
          var(:,:) = LCL(:,:)
       case ( 'LFC' )
          var(:,:) = LFC(:,:)
       case ( 'LNB' )
          var(:,:) = LNB(:,:)
       end select

    case ( 'PREC', 'RAIN', 'SNOW' )
       if ( .not. DV_calclated(I_PREC) ) then
          call allocate_2D( PREC )
          call allocate_2D( RAIN )
          call allocate_2D( SNOW )
          !$omp parallel do private(i,j) OMP_SCHEDULE_
          do j = 1, JA
          do i = 1, IA
             RAIN(i,j) = SFLX_rain_MP(i,j) + SFLX_rain_CP(i,j)
             SNOW(i,j) = SFLX_snow_MP(i,j)
             PREC(i,j) = RAIN(i,j) + SNOW(i,j)
          enddo
          enddo
          DV_calclated(I_PREC) = .true.
       end if
       select case (vname)
       case ( 'RAIN' )
          var(:,:) = RAIN(:,:)
       case ( 'SNOW' )
          var(:,:) = SNOW(:,:)
       case ( 'PREC' )
          var(:,:) = PREC(:,:)
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
    use scale_statistics, only: &
       STATISTICS_horizontal_mean
    use scale_atmos_grid_cartesC_real, only: &
       AREA => ATMOS_GRID_CARTESC_REAL_AREA
    implicit none
    character(len=*), intent(in)  :: vname
    real(RP),         intent(out) :: var(:)

    real(RP) :: WORK(KA,IA,JA)

    integer :: k, i, j, iq

    select case ( vname )
    case ( 'DENS_MEAN' )
       if ( .not. DV_calclated(I_DENS_MEAN) ) then
          call allocate_1D( DENS_MEAN )
          call STATISTICS_horizontal_mean( KA, KS, KE, IA, IS, IE, JA, JS, JE, &
                                           DENS(:,:,:), AREA(:,:), DENS_MEAN(:) )
          DV_calclated(I_DENS_MEAN) = .true.
       end if
       var(:) = DENS_MEAN(:)

    case ( 'W_MEAN' )
       if ( .not. DV_calclated(I_W_MEAN) ) then
          call allocate_1D( W_MEAN )
          call ATMOS_vars_get_diagnostic( 'DENS_MEAN', WORK1D(:) )
!OCL XFILL
          !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
          do j = JSB, JEB
          do i = ISB, IEB
          do k = KS, KE
             WORK(k,i,j) = W(k,i,j) * DENS_av(k,i,j)
          enddo
          enddo
          enddo
          call STATISTICS_horizontal_mean( KA, KS, KE, IA, IS, IE, JA, JS, JE, &
                                           WORK(:,:,:), AREA(:,:), W_MEAN(:) )
          do k = KS, KE
             W_MEAN(k) = W_MEAN(k) / DENS_MEAN(k)
          enddo
          DV_calclated(I_W_MEAN) = .true.
       end if
       var(:) = W_MEAN(:)

    case ( 'U_MEAN' )
       if ( .not. DV_calclated(I_U_MEAN) ) then
          call allocate_1D( U_MEAN )
          call ATMOS_vars_get_diagnostic( 'DENS_MEAN', WORK1D(:) )
!OCL XFILL
          !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
          do j = JSB, JEB
          do i = ISB, IEB
          do k = KS, KE
             WORK(k,i,j) = U(k,i,j) * DENS_av(k,i,j)
          enddo
          enddo
          enddo
          call STATISTICS_horizontal_mean( KA, KS, KE, IA, IS, IE, JA, JS, JE, &
                                           WORK(:,:,:), AREA(:,:), U_MEAN(:) )
          do k = KS, KE
             U_MEAN(k) = U_MEAN(k) / DENS_MEAN(k)
          enddo
          DV_calclated(I_U_MEAN) = .true.
       end if
       var(:) = U_MEAN(:)

    case ( 'V_MEAN' )
       if ( .not. DV_calclated(I_V_MEAN) ) then
          call allocate_1D( V_MEAN )
          call ATMOS_vars_get_diagnostic( 'DENS_MEAN', WORK1D(:) )
!OCL XFILL
          !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
          do j = JSB, JEB
          do i = ISB, IEB
          do k = KS, KE
             WORK(k,i,j) = V(k,i,j) * DENS_av(k,i,j)
          enddo
          enddo
          enddo
          call STATISTICS_horizontal_mean( KA, KS, KE, IA, IS, IE, JA, JS, JE, &
                                           WORK(:,:,:), AREA(:,:), V_MEAN(:) )
          do k = KS, KE
             V_MEAN(k) = V_MEAN(k) / DENS_MEAN(k)
          enddo
          DV_calclated(I_V_MEAN) = .true.
       end if
       var(:) = V_MEAN(:)

    case ( 'PT_MEAN' )
       if ( .not. DV_calclated(I_PT_MEAN) ) then
          call allocate_1D( PT_MEAN )
          call ATMOS_vars_get_diagnostic( 'DENS_MEAN', WORK1D(:) )
          call STATISTICS_horizontal_mean( KA, KS, KE, IA, IS, IE, JA, JS, JE, &
                                           RHOT(:,:,:), AREA(:,:), PT_MEAN(:) )
          do k = KS, KE
             PT_MEAN(k) = PT_MEAN(k) / DENS_MEAN(k)
          enddo
          DV_calclated(I_PT_MEAN) = .true.
       end if
       var(:) = PT_MEAN(:)

    case ( 'T_MEAN' )
       if ( .not. DV_calclated(I_T_MEAN) ) then
          call allocate_1D( T_MEAN )
          call ATMOS_vars_get_diagnostic( 'DENS_MEAN', WORK1D(:) )
!OCL XFILL
          !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
          do j = JSB, JEB
          do i = ISB, IEB
          do k = KS, KE
             WORK(k,i,j) = TEMP(k,i,j) * DENS_av(k,i,j)
          enddo
          enddo
          enddo
          call STATISTICS_horizontal_mean( KA, KS, KE, IA, IS, IE, JA, JS, JE, &
                                           WORK(:,:,:), AREA(:,:), T_MEAN(:) )
          do k = KS, KE
             T_MEAN(k) = T_MEAN(k) / DENS_MEAN(k)
          enddo
          DV_calclated(I_T_MEAN) = .true.
       end if
       var(:) = T_MEAN(:)

    case ( 'QV_MEAN' )
       if ( .not. DV_calclated(I_QV_MEAN) ) then
          call allocate_1D( QV_MEAN )
          if ( moist ) then
             call ATMOS_vars_get_diagnostic( 'DENS_MEAN', WORK1D(:) )
!OCL XFILL
             !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
             do j = JSB, JEB
             do i = ISB, IEB
             do k = KS, KE
                WORK(k,i,j) = QV(k,i,j) * DENS_av(k,i,j)
             enddo
             enddo
             enddo
             call STATISTICS_horizontal_mean( KA, KS, KE, IA, IS, IE, JA, JS, JE, &
                                              WORK(:,:,:), AREA(:,:), QV_MEAN(:) )
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
          call allocate_1D( QHYD_MEAN )
          call ATMOS_vars_get_diagnostic( 'DENS_MEAN', WORK1D(:) )
          call ATMOS_vars_get_diagnostic( 'QHYD', WORK3D(:,:,:) )
!OCL XFILL
          !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
          do j = JSB, JEB
          do i = ISB, IEB
          do k = KS, KE
             WORK(k,i,j) = QHYD(k,i,j) * DENS_av(k,i,j)
          enddo
          enddo
          enddo
          call STATISTICS_horizontal_mean( KA, KS, KE, IA, IS, IE, JA, JS, JE, &
                                           WORK(:,:,:), AREA(:,:), QHYD_MEAN(:) )
          do k = KS, KE
             QHYD_MEAN(k) = QHYD_MEAN(k) / DENS_MEAN(k)
          enddo
          DV_calclated(I_QHYD_MEAN) = .true.
       end if
       var(:) = QHYD_MEAN(:)

    case ( 'QLIQ_MEAN' )
       if ( .not. DV_calclated(I_QLIQ_MEAN) ) then
          call allocate_1D( QLIQ_MEAN )
          call ATMOS_vars_get_diagnostic( 'DENS_MEAN', WORK1D(:) )
          call ATMOS_vars_get_diagnostic( 'QLIQ', WORK3D(:,:,:) )
!OCL XFILL
          !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
          do j = JSB, JEB
          do i = ISB, IEB
          do k = KS, KE
             WORK(k,i,j) = QLIQ(k,i,j) * DENS_av(k,i,j)
          enddo
          enddo
          enddo
          call STATISTICS_horizontal_mean( KA, KS, KE, IA, IS, IE, JA, JS, JE, &
                                           WORK(:,:,:), AREA(:,:), QLIQ_MEAN(:) )
          do k = KS, KE
             QLIQ_MEAN(k) = QLIQ_MEAN(k) / DENS_MEAN(k)
          enddo
          DV_calclated(I_QLIQ_MEAN) = .true.
       end if
       var(:) = QLIQ_MEAN(:)

    case ( 'QICE_MEAN' )
       if ( .not. DV_calclated(I_QICE_MEAN) ) then
          call allocate_1D( QICE_MEAN )
          call ATMOS_vars_get_diagnostic( 'DENS_MEAN', WORK1D(:) )
          call ATMOS_vars_get_diagnostic( 'QICE', WORK3D(:,:,:) )
!OCL XFILL
          !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
          do j = JSB, JEB
          do i = ISB, IEB
          do k = KS, KE
             WORK(k,i,j) = QICE(k,i,j) * DENS_av(k,i,j)
          enddo
          enddo
          enddo
          call STATISTICS_horizontal_mean( KA, KS, KE, IA, IS, IE, JA, JS, JE, &
                                           WORK(:,:,:), AREA(:,:), QICE_MEAN(:) )
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
    use scale_atmos_grid_cartesC_real, only: &
       REAL_CZ => ATMOS_GRID_CARTESC_REAL_CZ
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
    use mod_atmos_phy_cp_vars, only: &
       SFLX_rain_CP => ATMOS_PHY_CP_SFLX_rain
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

    real(RP) :: RHOQ(KA,IA,JA)

    real(RP) :: ENGFLXT    (IA,JA) ! total flux             [J/m2/s]
    real(RP) :: SFLX_RD_net(IA,JA) ! net SFC radiation flux [J/m2/s]
    real(RP) :: TFLX_RD_net(IA,JA) ! net TOA radiation flux [J/m2/s]

    real(RP)               :: WORK (KA,IA,JA,3)
    character(len=H_SHORT) :: WNAME(3)
    real(RP)               :: CFLMAX

    integer  :: k, i, j, iq
    !---------------------------------------------------------------------------

    call MONITOR_put( PV_MONIT_id(I_DENS), DENS(:,:,:) )
    call MONITOR_put( PV_MONIT_id(I_MOMZ), MOMZ(:,:,:) )
    call MONITOR_put( PV_MONIT_id(I_MOMX), MOMX(:,:,:) )
    call MONITOR_put( PV_MONIT_id(I_MOMY), MOMY(:,:,:) )
    call MONITOR_put( PV_MONIT_id(I_RHOT), RHOT(:,:,:) )

    !##### Mass Budget #####

    do iq = 1, QA
       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
!OCL XFILL
       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          RHOQ(k,i,j) = DENS_av(k,i,j) * QTRC_av(k,i,j,iq)
       enddo
       enddo
       enddo

       call MONITOR_put( QP_MONIT_id(iq), RHOQ(:,:,:) )
    enddo

    ! total dry airmass
    if ( DV_MONIT_id(IM_QDRY) > 0 ) then
       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
!OCL XFILL
       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          RHOQ(k,i,j) = DENS(k,i,j) * QDRY (k,i,j)
       enddo
       enddo
       enddo
       call MONITOR_put( DV_MONIT_id(IM_QDRY), RHOQ(:,:,:) )
    end if

    ! total vapor,liquid,solid tracers
    if ( DV_MONIT_id(IM_QTOT) > 0 ) then
       call ATMOS_vars_get_diagnostic( 'QTOT', WORK3D(:,:,:) )
       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
!OCL XFILL
       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          RHOQ(k,i,j) = DENS(k,i,j) * QTOT(k,i,j)
       enddo
       enddo
       enddo
       call MONITOR_put( DV_MONIT_id(IM_QTOT), RHOQ(:,:,:) )
    end if

    ! total evapolation
    if ( moist ) call MONITOR_put( DV_MONIT_id(IM_EVAP), SFLX_QTRC(:,:,I_QV) )

    ! total precipitation
    if ( DV_MONIT_id(IM_PREC) > 0 ) then
       call ATMOS_vars_get_diagnostic( 'PREC', WORK2D(:,:) )
       call MONITOR_put( DV_MONIT_id(IM_PREC), WORK2D(:,:) )
    end if


    !##### Energy Budget #####

    if ( DV_MONIT_id(IM_ENGT) > 0 ) then
       call ATMOS_vars_get_diagnostic( 'ENGT', WORK3D(:,:,:) )
       call MONITOR_put( DV_MONIT_id(IM_ENGT), WORK3D(:,:,:) )
    end if
    if ( DV_MONIT_id(IM_ENGP) > 0 ) then
       call ATMOS_vars_get_diagnostic( 'ENGP', WORK3D(:,:,:) )
       call MONITOR_put( DV_MONIT_id(IM_ENGP), WORK3D(:,:,:) )
    end if
    if ( DV_MONIT_id(IM_ENGK) > 0 ) then
       call ATMOS_vars_get_diagnostic( 'ENGK', WORK3D(:,:,:) )
       call MONITOR_put( DV_MONIT_id(IM_ENGK), WORK3D(:,:,:) )
    end if
    if ( DV_MONIT_id(IM_ENGI) > 0 ) then
       call ATMOS_vars_get_diagnostic( 'ENGI', WORK3D(:,:,:) )
       call MONITOR_put( DV_MONIT_id(IM_ENGI), WORK3D(:,:,:) )
    end if


    ! radiation flux
!OCL XFILL
    !$omp parallel do private(i,j) OMP_SCHEDULE_ collapse(2)
    do j = JS, JE
    do i = IS, IE
       SFLX_RD_net(i,j) = ( SFLX_LW_up(i,j) - SFLX_LW_dn(i,j) ) &
                        + ( SFLX_SW_up(i,j) - SFLX_SW_dn(i,j) )

       TFLX_RD_net(i,j) = ( TOAFLX_LW_up(i,j) - TOAFLX_LW_dn(i,j) ) &
                        + ( TOAFLX_SW_up(i,j) - TOAFLX_SW_dn(i,j) )

       ENGFLXT    (i,j) = SFLX_SH(i,j) + SFLX_LH(i,j) &
                        + SFLX_RD_net(i,j) - TFLX_RD_net(i,j)
    enddo
    enddo

    call MONITOR_put( DV_MONIT_id(IM_ENGFLXT),      ENGFLXT     (:,:) )

    call MONITOR_put( DV_MONIT_id(IM_ENGSFC_SH),    SFLX_SH     (:,:) )
    call MONITOR_put( DV_MONIT_id(IM_ENGSFC_LH),    SFLX_LH     (:,:) )
    call MONITOR_put( DV_MONIT_id(IM_ENGSFC_RD),    SFLX_RD_net (:,:) )
    call MONITOR_put( DV_MONIT_id(IM_ENGTOA_RD),    TFLX_RD_net (:,:) )

    call MONITOR_put( DV_MONIT_id(IM_ENGSFC_LW_up), SFLX_LW_up  (:,:) )
    call MONITOR_put( DV_MONIT_id(IM_ENGSFC_LW_dn), SFLX_LW_dn  (:,:) )
    call MONITOR_put( DV_MONIT_id(IM_ENGSFC_SW_up), SFLX_SW_up  (:,:) )
    call MONITOR_put( DV_MONIT_id(IM_ENGSFC_SW_dn), SFLX_SW_dn  (:,:) )

    call MONITOR_put( DV_MONIT_id(IM_ENGTOA_LW_up), TOAFLX_LW_up(:,:) )
    call MONITOR_put( DV_MONIT_id(IM_ENGTOA_LW_dn), TOAFLX_LW_dn(:,:) )
    call MONITOR_put( DV_MONIT_id(IM_ENGTOA_SW_up), TOAFLX_SW_up(:,:) )
    call MONITOR_put( DV_MONIT_id(IM_ENGTOA_SW_dn), TOAFLX_SW_dn(:,:) )



    if ( ATMOS_VARS_CHECKRANGE ) then
!OCL XFILL
       WORK(:,:,:,1) = W(:,:,:)
!OCL XFILL
       WORK(:,:,:,2) = U(:,:,:)
!OCL XFILL
       WORK(:,:,:,3) = V(:,:,:)

       WNAME(1) = "W"
       WNAME(2) = "U"
       WNAME(3) = "V"

       call STATISTICS_detail( KA, KS, KE, IA, IS, IE, JA, JS, JE, 3, &
                               WNAME(:), WORK(:,:,:,:)                )
    endif

    if (       ( ATMOS_DYN_TYPE /= 'OFF' .AND. ATMOS_DYN_TYPE /= 'NONE' )                   &
         .AND. ( ATMOS_VARS_CHECKCFL_SOFT > 0.0_RP .OR. ATMOS_VARS_CHECKCFL_HARD > 0.0_RP ) ) then
!OCL XFILL
       WORK(:,:,:,:) = 0.0_RP

       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          WORK(k,i,j,1) = 0.5_RP * abs(MOMZ_av(k,i,j)) / ( DENS_av(k+1,i,j) + DENS_av(k,i,j) ) &
                        * TIME_DTSEC_ATMOS_DYN / ( REAL_CZ(k+1,i,j) - REAL_CZ(k,i,j) )
          WORK(k,i,j,2) = 0.5_RP * abs(MOMX_av(k,i,j)) / ( DENS_av(k,i+1,j) + DENS_av(k,i,j) ) &
                        * TIME_DTSEC_ATMOS_DYN * RFDX(i) * MAPF(i,j,1,I_UY)
          WORK(k,i,j,3) = 0.5_RP * abs(MOMY_av(k,i,j)) / ( DENS_av(k,i,j+1) + DENS_av(k,i,j) ) &
                        * TIME_DTSEC_ATMOS_DYN * RFDY(j) * MAPF(i,j,2,I_XV)
       enddo
       enddo
       enddo

       CFLMAX = maxval( WORK(:,:,:,:) )

       if ( ATMOS_VARS_CHECKCFL_HARD > 0.0_RP .AND. CFLMAX > ATMOS_VARS_CHECKCFL_HARD ) then
          LOG_INFO("ATMOS_vars_monitor",*) "Courant number =", CFLMAX, " exceeded the hard limit =", ATMOS_VARS_CHECKCFL_HARD
          LOG_ERROR("ATMOS_vars_monitor",*)                     "Courant number =", CFLMAX, " exceeded the hard limit =", ATMOS_VARS_CHECKCFL_HARD
          LOG_ERROR_CONT(*)                     "Rank =", PRC_myrank
          LOG_ERROR_CONT(*)                     "Please set ATMOS_VARS_CHECKCFL_HARD in the namelist PARAM_ATMOS_VARS when you want to change the limit."

          WNAME(1) = "Courant num. Z"
          WNAME(2) = "Courant num. X"
          WNAME(3) = "Courant num. Y"
          call STATISTICS_detail( KA, KS, KE, IA, IS, IE, JA, JS, JE, 3, &
                                  WNAME(:), WORK(:,:,:,:),               &
                                  local=.true.                           )

          call PRC_abort
       endif

       if ( ATMOS_VARS_CHECKCFL_SOFT > 0.0_RP .AND. CFLMAX > ATMOS_VARS_CHECKCFL_SOFT ) then
          LOG_INFO("ATMOS_vars_monitor",*) "Courant number =", CFLMAX, " exceeded the soft limit =", ATMOS_VARS_CHECKCFL_SOFT
          LOG_ERROR("ATMOS_vars_monitor",*)                     "Courant number =", CFLMAX, " exceeded the soft limit =", ATMOS_VARS_CHECKCFL_SOFT
          LOG_ERROR_CONT(*)                     "Rank =", PRC_myrank

          WNAME(1) = "Courant num. Z"
          WNAME(2) = "Courant num. X"
          WNAME(3) = "Courant num. Y"
          call STATISTICS_detail( KA, KS, KE, IA, IS, IE, JA, JS, JE, 3, &
                                  WNAME(:), WORK(:,:,:,:),               &
                                  local=.true.                           )
       endif
    endif

    return
  end subroutine ATMOS_vars_monitor

  !-----------------------------------------------------------------------------
  !> Create atmospheric restart file
  subroutine ATMOS_vars_restart_create
    use scale_time, only: &
       TIME_gettimelabel
    use scale_file_cartesC, only: &
       FILE_CARTESC_create
    use mod_atmos_admin, only: &
       ATMOS_sw_dyn,      &
       ATMOS_sw_phy_mp,   &
       ATMOS_sw_phy_ae,   &
       ATMOS_sw_phy_ch,   &
       ATMOS_sw_phy_rd,   &
       ATMOS_sw_phy_sf,   &
       ATMOS_sw_phy_tb,   &
       ATMOS_sw_phy_cp
    use mod_cpl_admin, only: &
       CPL_sw
    use mod_atmos_dyn_vars, only: &
       ATMOS_DYN_vars_restart_create
    use mod_atmos_phy_mp_vars, only: &
       ATMOS_PHY_MP_vars_restart_create
    use mod_atmos_phy_ae_vars, only: &
       ATMOS_PHY_AE_vars_restart_create
    use mod_atmos_phy_ch_vars, only: &
       ATMOS_PHY_CH_vars_restart_create
    use mod_atmos_phy_rd_vars, only: &
       ATMOS_PHY_RD_vars_restart_create
    use mod_atmos_phy_sf_vars, only: &
       ATMOS_PHY_SF_vars_restart_create
    use mod_atmos_phy_tb_vars, only: &
       ATMOS_PHY_TB_vars_restart_create
    use mod_atmos_phy_cp_vars, only: &
       ATMOS_PHY_CP_vars_restart_create
#ifdef SDM
    use scale_atmos_phy_mp_sdm, only: &
       sd_rest_flg_out, &
       ATMOS_PHY_MP_sdm_restart_create
    use scale_time, only: &
       NOWSEC => TIME_NOWSEC
#endif
    implicit none

    character(len=19)     :: timelabel
    character(len=H_LONG) :: basename
    !---------------------------------------------------------------------------

#ifdef SDM
    if( sd_rest_flg_out ) then
       LOG_INFO("ATMOS_vars_restart_create",*) 'Output random number for SDM '
       call ATMOS_PHY_MP_sdm_restart_create(NOWSEC)
    endif
#endif

    if ( ATMOS_RESTART_OUT_BASENAME /= '' ) then

       LOG_NEWLINE
       LOG_INFO("ATMOS_vars_restart_create",*) 'Create restart file (ATMOS) '

       if ( ATMOS_RESTART_OUT_POSTFIX_TIMELABEL ) then
          call TIME_gettimelabel( timelabel )
          basename = trim(ATMOS_RESTART_OUT_BASENAME)//'_'//trim(timelabel)
       else
          basename = trim(ATMOS_RESTART_OUT_BASENAME)
       endif

       LOG_INFO("ATMOS_vars_restart_create",*) 'basename: ', trim(basename)

       call FILE_CARTESC_create( &
            basename, ATMOS_RESTART_OUT_TITLE, ATMOS_RESTART_OUT_DTYPE, & ! [IN]
            restart_fid,                                                & ! [OUT]
            aggregate=ATMOS_RESTART_OUT_AGGREGATE                       ) ! [IN]

       allocate( PV_ID(PV_nmax+QA) )
    endif

    if( ATMOS_sw_dyn )    call ATMOS_DYN_vars_restart_create
    if( ATMOS_sw_phy_mp ) call ATMOS_PHY_MP_vars_restart_create
    if( ATMOS_sw_phy_ae ) call ATMOS_PHY_AE_vars_restart_create
    if( ATMOS_sw_phy_ch ) call ATMOS_PHY_CH_vars_restart_create
    if( ATMOS_sw_phy_rd ) call ATMOS_PHY_RD_vars_restart_create
    if( ATMOS_sw_phy_sf .and. (.not. CPL_sw) ) call ATMOS_PHY_SF_vars_restart_create
    if( ATMOS_sw_phy_tb ) call ATMOS_PHY_TB_vars_restart_create
    if( ATMOS_sw_phy_cp ) call ATMOS_PHY_CP_vars_restart_create

    return
  end subroutine ATMOS_vars_restart_create

  !-----------------------------------------------------------------------------
  !> Exit netCDF define mode
  subroutine ATMOS_vars_restart_enddef
    use scale_file_cartesC, only: &
       FILE_CARTESC_enddef
    use mod_atmos_admin, only: &
       ATMOS_sw_dyn,      &
       ATMOS_sw_phy_mp,   &
       ATMOS_sw_phy_ae,   &
       ATMOS_sw_phy_ch,   &
       ATMOS_sw_phy_rd,   &
       ATMOS_sw_phy_sf,   &
       ATMOS_sw_phy_tb,   &
       ATMOS_sw_phy_cp
    use mod_cpl_admin, only: &
       CPL_sw
    use mod_atmos_dyn_vars, only: &
       ATMOS_DYN_vars_restart_enddef
    use mod_atmos_phy_mp_vars, only: &
       ATMOS_PHY_MP_vars_restart_enddef
    use mod_atmos_phy_ae_vars, only: &
       ATMOS_PHY_AE_vars_restart_enddef
    use mod_atmos_phy_ch_vars, only: &
       ATMOS_PHY_CH_vars_restart_enddef
    use mod_atmos_phy_rd_vars, only: &
       ATMOS_PHY_RD_vars_restart_enddef
    use mod_atmos_phy_sf_vars, only: &
       ATMOS_PHY_SF_vars_restart_enddef
    use mod_atmos_phy_tb_vars, only: &
       ATMOS_PHY_TB_vars_restart_enddef
    use mod_atmos_phy_cp_vars, only: &
       ATMOS_PHY_CP_vars_restart_enddef
#ifdef SDM
    use scale_atmos_phy_mp_sdm, only: &
       sd_rest_flg_out, &
       ATMOS_PHY_MP_sdm_restart_enddef
#endif
    implicit none

    !---------------------------------------------------------------------------

#ifdef SDM
    if( sd_rest_flg_out ) then
       call ATMOS_PHY_MP_sdm_restart_enddef
    endif
#endif

    if ( restart_fid /= -1 ) then
       call FILE_CARTESC_enddef( restart_fid ) ! [IN]
    endif

    if( ATMOS_sw_dyn )    call ATMOS_DYN_vars_restart_enddef
    if( ATMOS_sw_phy_mp ) call ATMOS_PHY_MP_vars_restart_enddef
    if( ATMOS_sw_phy_ae ) call ATMOS_PHY_AE_vars_restart_enddef
    if( ATMOS_sw_phy_ch ) call ATMOS_PHY_CH_vars_restart_enddef
    if( ATMOS_sw_phy_rd ) call ATMOS_PHY_RD_vars_restart_enddef
    if( ATMOS_sw_phy_sf .and. (.not. CPL_sw) ) call ATMOS_PHY_SF_vars_restart_enddef
    if( ATMOS_sw_phy_tb ) call ATMOS_PHY_TB_vars_restart_enddef
    if( ATMOS_sw_phy_cp ) call ATMOS_PHY_CP_vars_restart_enddef

    return
  end subroutine ATMOS_vars_restart_enddef

  !-----------------------------------------------------------------------------
  !> Close restart file
  subroutine ATMOS_vars_restart_close
    use scale_file_cartesC, only: &
       FILE_CARTESC_close
    use mod_atmos_admin, only: &
       ATMOS_sw_dyn,      &
       ATMOS_sw_phy_mp,   &
       ATMOS_sw_phy_ae,   &
       ATMOS_sw_phy_ch,   &
       ATMOS_sw_phy_rd,   &
       ATMOS_sw_phy_sf,   &
       ATMOS_sw_phy_tb,   &
       ATMOS_sw_phy_cp
    use mod_cpl_admin, only: &
       CPL_sw
    use mod_atmos_dyn_vars, only: &
       ATMOS_DYN_vars_restart_close
    use mod_atmos_phy_mp_vars, only: &
       ATMOS_PHY_MP_vars_restart_close
    use mod_atmos_phy_ae_vars, only: &
       ATMOS_PHY_AE_vars_restart_close
    use mod_atmos_phy_ch_vars, only: &
       ATMOS_PHY_CH_vars_restart_close
    use mod_atmos_phy_rd_vars, only: &
       ATMOS_PHY_RD_vars_restart_close
    use mod_atmos_phy_sf_vars, only: &
       ATMOS_PHY_SF_vars_restart_close
    use mod_atmos_phy_tb_vars, only: &
       ATMOS_PHY_TB_vars_restart_close
    use mod_atmos_phy_cp_vars, only: &
       ATMOS_PHY_CP_vars_restart_close
#ifdef SDM
    use scale_atmos_phy_mp_sdm, only: &
       sd_rest_flg_out, &
       ATMOS_PHY_MP_sdm_restart_close
#endif
    implicit none
    !---------------------------------------------------------------------------

#ifdef SDM
    if( sd_rest_flg_out ) then
       call ATMOS_PHY_MP_sdm_restart_close
    endif
#endif

    if ( restart_fid /= -1 ) then
       LOG_NEWLINE
       LOG_INFO("ATMOS_vars_restart_close",*) 'Close restart file (ATMOS) '

       call FILE_CARTESC_close( restart_fid ) ! [IN]

       restart_fid = -1

       if ( allocated(PV_ID) ) deallocate( PV_ID )
    endif

    if( ATMOS_sw_dyn )    call ATMOS_DYN_vars_restart_close
    if( ATMOS_sw_phy_mp ) call ATMOS_PHY_MP_vars_restart_close
    if( ATMOS_sw_phy_ae ) call ATMOS_PHY_AE_vars_restart_close
    if( ATMOS_sw_phy_ch ) call ATMOS_PHY_CH_vars_restart_close
    if( ATMOS_sw_phy_rd ) call ATMOS_PHY_RD_vars_restart_close
    if( ATMOS_sw_phy_sf .and. (.not. CPL_sw) ) call ATMOS_PHY_SF_vars_restart_close
    if( ATMOS_sw_phy_tb ) call ATMOS_PHY_TB_vars_restart_close
    if( ATMOS_sw_phy_cp ) call ATMOS_PHY_CP_vars_restart_close

    return
  end subroutine ATMOS_vars_restart_close

  !-----------------------------------------------------------------------------
  !> Define atmospheric variables in restart file
  subroutine ATMOS_vars_restart_def_var
    use scale_file_cartesC, only: &
       FILE_CARTESC_def_var
    use mod_atmos_admin, only: &
       ATMOS_sw_dyn,      &
       ATMOS_sw_phy_mp,   &
       ATMOS_sw_phy_ae,   &
       ATMOS_sw_phy_ch,   &
       ATMOS_sw_phy_rd,   &
       ATMOS_sw_phy_sf,   &
       ATMOS_sw_phy_tb,   &
       ATMOS_sw_phy_cp
    use mod_cpl_admin, only: &
       CPL_sw
    use mod_atmos_dyn_vars, only: &
       ATMOS_DYN_vars_restart_def_var
    use mod_atmos_phy_mp_vars, only: &
       ATMOS_PHY_MP_vars_restart_def_var
    use mod_atmos_phy_ae_vars, only: &
       ATMOS_PHY_AE_vars_restart_def_var
    use mod_atmos_phy_ch_vars, only: &
       ATMOS_PHY_CH_vars_restart_def_var
    use mod_atmos_phy_rd_vars, only: &
       ATMOS_PHY_RD_vars_restart_def_var
    use mod_atmos_phy_sf_vars, only: &
       ATMOS_PHY_SF_vars_restart_def_var
    use mod_atmos_phy_tb_vars, only: &
       ATMOS_PHY_TB_vars_restart_def_var
    use mod_atmos_phy_cp_vars, only: &
       ATMOS_PHY_CP_vars_restart_def_var
#ifdef SDM
    use scale_atmos_phy_mp_sdm, only: &
       sd_rest_flg_out, &
       ATMOS_PHY_MP_sdm_restart_def_var
#endif
    implicit none

    integer iq
    !---------------------------------------------------------------------------

#ifdef SDM
    if( sd_rest_flg_out ) then
       call ATMOS_PHY_MP_sdm_restart_def_var
    endif
#endif

    if ( restart_fid /= -1 ) then

       call FILE_CARTESC_def_var( restart_fid, PV_info(I_DENS)%NAME, PV_info(I_DENS)%DESC, PV_info(I_DENS)%UNIT, 'ZXY',  ATMOS_RESTART_OUT_DTYPE, &
                                  PV_ID(I_DENS), &
                                  standard_name=PV_info(I_DENS)%STDNAME )
       call FILE_CARTESC_def_var( restart_fid, PV_info(I_MOMZ)%NAME, PV_info(I_MOMZ)%DESC, PV_info(I_MOMZ)%UNIT, 'ZHXY', ATMOS_RESTART_OUT_DTYPE, &
                                  PV_ID(I_MOMZ), &
                                  standard_name=PV_info(I_MOMZ)%STDNAME )
       call FILE_CARTESC_def_var( restart_fid, PV_info(I_MOMX)%NAME, PV_info(I_MOMX)%DESC, PV_info(I_MOMX)%UNIT, 'ZXHY', ATMOS_RESTART_OUT_DTYPE, &
                                  PV_ID(I_MOMX), &
                                  standard_name=PV_info(I_MOMX)%STDNAME )
       call FILE_CARTESC_def_var( restart_fid, PV_info(I_MOMY)%NAME, PV_info(I_MOMY)%DESC, PV_info(I_MOMY)%UNIT, 'ZXYH', ATMOS_RESTART_OUT_DTYPE, &
                                  PV_ID(I_MOMY), &
                                  standard_name=PV_info(I_MOMY)%STDNAME )
       call FILE_CARTESC_def_var( restart_fid, PV_info(I_RHOT)%NAME, PV_info(I_RHOT)%DESC, PV_info(I_RHOT)%UNIT, 'ZXY',  ATMOS_RESTART_OUT_DTYPE, &
                                  PV_ID(I_RHOT), &
                                  standard_name=PV_info(I_RHOT)%STDNAME )
       do iq = 1, QA
          call FILE_CARTESC_def_var( restart_fid, TRACER_NAME(iq), TRACER_DESC(iq), TRACER_UNIT(iq), 'ZXY',  ATMOS_RESTART_OUT_DTYPE, &
                                     PV_ID(PV_nmax+iq) )
       enddo

    endif

    if( ATMOS_sw_dyn )    call ATMOS_DYN_vars_restart_def_var
    if( ATMOS_sw_phy_mp ) call ATMOS_PHY_MP_vars_restart_def_var
    if( ATMOS_sw_phy_ae ) call ATMOS_PHY_AE_vars_restart_def_var
    if( ATMOS_sw_phy_ch ) call ATMOS_PHY_CH_vars_restart_def_var
    if( ATMOS_sw_phy_rd ) call ATMOS_PHY_RD_vars_restart_def_var
    if( ATMOS_sw_phy_sf .and. (.not. CPL_sw) ) call ATMOS_PHY_SF_vars_restart_def_var
    if( ATMOS_sw_phy_tb ) call ATMOS_PHY_TB_vars_restart_def_var
    if( ATMOS_sw_phy_cp ) call ATMOS_PHY_CP_vars_restart_def_var

    return
  end subroutine ATMOS_vars_restart_def_var

  !-----------------------------------------------------------------------------
  !> Write restart of atmospheric variables
  subroutine ATMOS_vars_restart_write
    use scale_file_cartesC, only: &
       FILE_CARTESC_write_var
    use mod_atmos_admin, only: &
       ATMOS_sw_dyn,      &
       ATMOS_sw_phy_mp,   &
       ATMOS_sw_phy_ae,   &
       ATMOS_sw_phy_ch,   &
       ATMOS_sw_phy_rd,   &
       ATMOS_sw_phy_sf,   &
       ATMOS_sw_phy_tb,   &
       ATMOS_sw_phy_cp
    use mod_cpl_admin, only: &
       CPL_sw
    use mod_atmos_dyn_vars, only: &
       ATMOS_DYN_vars_restart_write
    use mod_atmos_phy_mp_vars, only: &
       ATMOS_PHY_MP_vars_restart_write
    use mod_atmos_phy_ae_vars, only: &
       ATMOS_PHY_AE_vars_restart_write
    use mod_atmos_phy_ch_vars, only: &
       ATMOS_PHY_CH_vars_restart_write
    use mod_atmos_phy_rd_vars, only: &
       ATMOS_PHY_RD_vars_restart_write
    use mod_atmos_phy_sf_vars, only: &
       ATMOS_PHY_SF_vars_restart_write
    use mod_atmos_phy_tb_vars, only: &
       ATMOS_PHY_TB_vars_restart_write
    use mod_atmos_phy_cp_vars, only: &
       ATMOS_PHY_CP_vars_restart_write
#ifdef SDM
    use scale_atmos_phy_mp_sdm, only: &
       sd_rest_flg_out, &
       ATMOS_PHY_MP_sdm_restart_write
#endif
    implicit none

    integer iq
    !---------------------------------------------------------------------------

#ifdef SDM
    if( sd_rest_flg_out ) then
       call ATMOS_PHY_MP_sdm_restart_write
    endif
#endif

    if ( restart_fid /= -1 ) then

       call ATMOS_vars_fillhalo

       call ATMOS_vars_total

       call FILE_CARTESC_write_var( restart_fid, PV_ID(I_DENS), DENS(:,:,:), PV_info(I_DENS)%NAME, 'ZXY'  ) ! [IN]
       call FILE_CARTESC_write_var( restart_fid, PV_ID(I_MOMZ), MOMZ(:,:,:), PV_info(I_MOMZ)%NAME, 'ZHXY' ) ! [IN]
       call FILE_CARTESC_write_var( restart_fid, PV_ID(I_MOMX), MOMX(:,:,:), PV_info(I_MOMX)%NAME, 'ZXHY' ) ! [IN]
       call FILE_CARTESC_write_var( restart_fid, PV_ID(I_MOMY), MOMY(:,:,:), PV_info(I_MOMY)%NAME, 'ZXYH' ) ! [IN]
       call FILE_CARTESC_write_var( restart_fid, PV_ID(I_RHOT), RHOT(:,:,:), PV_info(I_RHOT)%NAME, 'ZXY'  ) ! [IN]

       do iq = 1, QA
          call FILE_CARTESC_write_var( restart_fid, PV_ID(PV_nmax+iq), QTRC(:,:,:,iq), TRACER_NAME(iq), 'ZXY' ) ! [IN]
       enddo

    endif

    if( ATMOS_sw_dyn )    call ATMOS_DYN_vars_restart_write
    if( ATMOS_sw_phy_mp ) call ATMOS_PHY_MP_vars_restart_write
    if( ATMOS_sw_phy_ae ) call ATMOS_PHY_AE_vars_restart_write
    if( ATMOS_sw_phy_ch ) call ATMOS_PHY_CH_vars_restart_write
    if( ATMOS_sw_phy_rd ) call ATMOS_PHY_RD_vars_restart_write
    if( ATMOS_sw_phy_sf .and. (.not. CPL_sw) ) call ATMOS_PHY_SF_vars_restart_write
    if( ATMOS_sw_phy_tb ) call ATMOS_PHY_TB_vars_restart_write
    if( ATMOS_sw_phy_cp ) call ATMOS_PHY_CP_vars_restart_write

    return
  end subroutine ATMOS_vars_restart_write


  ! private
  subroutine allocate_3D( ary )
    use scale_const, only: &
       UNDEF => CONST_UNDEF
    real(RP), intent(inout), allocatable :: ary(:,:,:)

    if ( .not. allocated(ary) ) then
       allocate( ary(KA,IA,JA) )
       ary(:,:,:) = UNDEF
    end if

    return
  end subroutine allocate_3D

  subroutine allocate_2D( ary )
    use scale_const, only: &
       UNDEF => CONST_UNDEF
    real(RP), intent(inout), allocatable :: ary(:,:)

    if ( .not. allocated(ary) ) then
       allocate( ary(IA,JA) )
       ary(:,:) = UNDEF
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
