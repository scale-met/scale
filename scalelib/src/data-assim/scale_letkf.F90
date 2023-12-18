!-------------------------------------------------------------------------------
!> module LETKF for Data-Assimilation
!!
!! @par Description
!!
!! @author Team SCALE
!!         imported from the DA system compiled by Data Assimilation Team
!!
!<
!-------------------------------------------------------------------------------
#include "scalelib.h"
module scale_letkf
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_io
  use scale_prof

  use mpi

  use scale_prc, only: &
    PRC_UNIVERSAL_IsMaster

  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: LETKF_setup
  public :: LETKF_finalize
  public :: LETKF_obs_readfile
  public :: LETKF_obs_clear
  public :: LETKF_obs_operator
  public :: LETKF_obs_initialize
  public :: LETKF_system
  public :: LETKF_param_estimation_system
  public :: LETKF_add_inflation_setup

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  integer, public, parameter :: nv3d                  = 11              ! number of 3D prognostic variables
  integer, public, parameter :: nv2d                  = 0               ! number of 2D prognostic variables
  integer, public, parameter :: nid_obs               = 16              ! number of variable types
  integer, public, parameter :: nobtype               = 24              ! number of observation report types

  integer, public, parameter :: max_obs_info_meta     = 3               ! maximum array size for type(obs_info)%meta
  integer, public, parameter :: n_qc_steps            = 2
  integer, public, parameter :: i_before_qc           = 1
  integer, public, parameter :: i_after_qc            = 2

  type, public :: obs_info
    integer               :: nobs = 0
    integer,  allocatable :: elm(:)
    real(RP), allocatable :: lon(:)
    real(RP), allocatable :: lat(:)
    real(RP), allocatable :: lev(:)
    real(RP), allocatable :: dat(:)
    real(RP), allocatable :: err(:)
    integer,  allocatable :: typ(:)
    real(RP), allocatable :: dif(:)
    real(RP)              :: meta(max_obs_info_meta) = -9.99e+33_RP
    real(RP), allocatable :: ri(:)
    real(RP), allocatable :: rj(:)
    integer,  allocatable :: rank(:)
  end type obs_info

  type, public :: obs_da_value
    integer               :: nobs = 0
    integer               :: nobs_in_key = 0
    integer,  allocatable :: set(:)
    integer,  allocatable :: idx(:)
    integer,  allocatable :: key(:)
    real(RP), allocatable :: val(:)
    real(RP), allocatable :: ensval(:,:)
    real(RP), allocatable :: eqv(:,:)          ! qv (ensemble)
    real(RP), allocatable :: qv(:)             ! qv (mean)
    real(RP), allocatable :: tm(:)             ! temp (mean)
    real(RP), allocatable :: pm(:)             ! pressure (mean)
    integer,  allocatable :: qc(:)
  end type obs_da_value

  type, public :: obs_grid_type
    integer              :: ngrd_i
    integer              :: ngrd_j
    real(RP)             :: grdspc_i
    real(RP)             :: grdspc_j
    integer              :: ngrdsch_i
    integer              :: ngrdsch_j
    integer              :: ngrdext_i
    integer              :: ngrdext_j
    integer, allocatable :: n(:,:,:)
    integer, allocatable :: ac(:,:,:)
    integer, allocatable :: tot(:)
    integer, allocatable :: n_ext(:,:)
    integer, allocatable :: ac_ext(:,:)
    integer              :: tot_ext
    integer              :: tot_sub(n_qc_steps)    ! only for diagnostic print; 1: before QC; 2: after QC
    integer              :: tot_g(n_qc_steps)      ! only for diagnostic print
    integer, allocatable :: next(:,:)              ! temporary array
  end type obs_grid_type

  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  private :: letkf_core
  private :: get_nobs
  private :: get_nobs_radar
  private :: read_obs
  private :: read_obs_radar
  private :: read_obs_radar_toshiba_pawr
  private :: read_obs_radar_toshiba_mp_pawr
  private :: str_replace
  private :: jst2utc
  private :: radar_georeference
  private :: define_grid
  private :: calc_ref_vr
  private :: obs_local
  private :: obs_local_range
  private :: obs_local_cal
  private :: relax_beta
  private :: weight_RTPP
  private :: weight_RTPS
  private :: com_distll_1
  private :: ensmean_grd
  private :: Trans_XtoY
  private :: Trans_XtoY_radar
  private :: prsadj
  private :: phys2ijk
  private :: phys2ijkz
  private :: phys2ij
  private :: itpl_2d
  private :: itpl_2d_column
  private :: itpl_3d
  private :: merge_sort_parallel
  private :: merge_sort_2threads
  private :: merge_2threads
  private :: merge
  private :: merge_sort_mpi
  private :: merge_mpi
  private :: merge_mpi_no_nest
  private :: scatter_grd_mpi
  private :: scatter_grd_mpi_all2all
  private :: gather_grd_mpi_all2all
  private :: grd_to_buf
  private :: buf_to_grd
  private :: calc_z_grd
  private :: qc_indexing_and_packing
  private :: uid_obs
  private :: uid_obs_varlocal
  private :: binary_search_i8
  private :: rank_1d_2d
  private :: rank_2d_1d
  private :: rij_rank
  private :: rij_g2l
  private :: ij_obsgrd
  private :: ij_obsgrd_ext
  private :: obs_choose
  private :: obs_choose_ext
  private :: obs_info_allocate
  private :: obs_info_deallocate
  private :: obs_da_value_allocate
  private :: obs_da_value_deallocate
  private :: obs_da_value_allreduce
  private :: obs_da_value_partial_reduce_iter
  private :: read_ens_mpi_addiinfl
  private :: monit_obs
  private :: monit_obs_mpi
  private :: monit_dep
  private :: monit_print
  private :: state_to_history

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  integer, parameter :: n_search_incr                      = 8
  integer, parameter :: nobsfilemax                        = 10
  integer, parameter :: membermax                          = 10000

  integer, parameter :: iv3d_rho                           = 1               !-- State in restart files
  integer, parameter :: iv3d_rhou                          = 2               !
  integer, parameter :: iv3d_rhov                          = 3               !
  integer, parameter :: iv3d_rhow                          = 4               !
  integer, parameter :: iv3d_rhot                          = 5               !
  integer, parameter :: iv3d_u                             = 1               !-- State for LETKF
  integer, parameter :: iv3d_v                             = 2               !
  integer, parameter :: iv3d_w                             = 3               !
  integer, parameter :: iv3d_t                             = 4               !
  integer, parameter :: iv3d_p                             = 5               !
  integer, parameter :: iv3d_q                             = 6               !
  integer, parameter :: iv3d_qc                            = 7               !
  integer, parameter :: iv3d_qr                            = 8               !
  integer, parameter :: iv3d_qi                            = 9               !
  integer, parameter :: iv3d_qs                            = 10              !
  integer, parameter :: iv3d_qg                            = 11              !

  !--- 3D, 2D diagnostic variables (in SCALE history files)
  integer, parameter :: nv3dd                              = 13
  integer, parameter :: nv2dd                              = 7
  integer, parameter :: iv3dd_u                            = 1
  integer, parameter :: iv3dd_v                            = 2
  integer, parameter :: iv3dd_w                            = 3
  integer, parameter :: iv3dd_t                            = 4
  integer, parameter :: iv3dd_p                            = 5
  integer, parameter :: iv3dd_q                            = 6
  integer, parameter :: iv3dd_qc                           = 7
  integer, parameter :: iv3dd_qr                           = 8
  integer, parameter :: iv3dd_qi                           = 9
  integer, parameter :: iv3dd_qs                           = 10
  integer, parameter :: iv3dd_qg                           = 11
  integer, parameter :: iv3dd_rh                           = 12
  integer, parameter :: iv3dd_hgt                          = 13
  integer, parameter :: iv2dd_topo                         = 1
  integer, parameter :: iv2dd_ps                           = 2
  integer, parameter :: iv2dd_rain                         = 3
  integer, parameter :: iv2dd_u10m                         = 4
  integer, parameter :: iv2dd_v10m                         = 5
  integer, parameter :: iv2dd_t2m                          = 6
  integer, parameter :: iv2dd_q2m                          = 7

  integer, parameter :: iqc_good                           = 0
  integer, parameter :: iqc_gross_err                      = 5
  integer, parameter :: iqc_ps_ter                         = 10
  integer, parameter :: iqc_ref_low                        = 11
  integer, parameter :: iqc_ref_mem                        = 12
  integer, parameter :: iqc_radar_vhi                      = 19
  integer, parameter :: iqc_out_vhi                        = 20
  integer, parameter :: iqc_out_vlo                        = 21
  integer, parameter :: iqc_obs_bad                        = 50
  integer, parameter :: iqc_otype                          = 90
  integer, parameter :: iqc_time                           = 97
  integer, parameter :: iqc_out_h                          = 98
  integer, parameter :: iqc_undef                          = 99

  integer, parameter :: nid_obs_varlocal                   = 8
  !
  ! conventional observations
  !
  integer, parameter :: id_u_obs                           = 2819
  integer, parameter :: id_v_obs                           = 2820
  integer, parameter :: id_t_obs                           = 3073
  integer, parameter :: id_tv_obs                          = 3074
  integer, parameter :: id_q_obs                           = 3330
  integer, parameter :: id_rh_obs                          = 3331
  !
  ! surface observations codes > 9999
  !
  integer, parameter :: id_ps_obs                          = 14593
  integer, parameter :: id_rain_obs                        = 19999
  integer, parameter :: id_tclon_obs                       = 99991           ! TC vital
  integer, parameter :: id_tclat_obs                       = 99992           ! TC vital
  integer, parameter :: id_tcmip_obs                       = 99993           ! TC vital
  !
  ! radar observations
  !
  integer, parameter :: id_radar_ref_obs                   = 4001
  integer, parameter :: id_radar_ref_zero_obs              = 4004
  integer, parameter :: id_radar_vr_obs                    = 4002
  integer, parameter :: id_radar_prh_obs                   = 4003
  !
  ! Himawari-8 (H08) observations
  !
  integer, parameter :: id_H08IR_obs                       = 8800

  real(RP), parameter :: dist_zero_fac                     = 3.651483717     ! SQRT(10.0d0/3.0d0) * 2.0d0
  real(RP), parameter :: dist_zero_fac_square              = 13.33333333     ! dist_zero_fac * dist_zero_fac
  real(RP), parameter :: minz                              = 0.01_RP         ! Minimum radar power.
  real(RP), parameter :: vr_min_dist                       = 8000.0_RP

  integer, parameter :: elem_uid(nid_obs) = &
     (/ id_u_obs,        id_v_obs,         id_t_obs,     id_tv_obs,        id_q_obs,              &
        id_rh_obs,       id_ps_obs,        id_rain_obs,  id_radar_ref_obs, id_radar_ref_zero_obs, &
        id_radar_vr_obs, id_radar_prh_obs, id_H08IR_obs, id_tclon_obs,     id_tclat_obs,          &
        id_tcmip_obs /)

  character(3), parameter :: obelmlist(nid_obs) = &
     (/ '  U', '  V', '  T', ' Tv', '  Q', &
        ' RH', ' PS', 'PRC', 'REF', 'RE0', &
        ' Vr', 'PRH', 'H08', 'TCX', 'TCY', &
        'TCP' /)

  character(3), parameter :: obelmlist_varlocal(nid_obs_varlocal) = &
     (/ 'WND', '  T', 'MOI', ' PS', 'PRC', &
        'TCV', 'REF', ' Vr' /)

  character(6), parameter :: obtypelist(nobtype) =        &
     (/ 'ADPUPA', 'AIRCAR', 'AIRCFT', 'SATWND', 'PROFLR', &
        'VADWND', 'SATEMP', 'ADPSFC', 'SFCSHP', 'SFCBOG', &
        'SPSSMI', 'SYNDAT', 'ERS1DA', 'GOESND', 'QKSWND', &
        'MSONET', 'GPSIPW', 'RASSDA', 'WDSATR', 'ASCATW', &
        'TMPAPR', 'PHARAD', 'H08IRB', 'TCVITL' /) ! H08

  logical :: ENS_WITH_MDET                                 = .false.         ! Run additional member of 'mdet'?

  logical :: USE_OBS(nobtype)                              = .true.
  logical :: OBS_POSTFIX_TIMELABEL                         = .false.
  logical :: OBSDA_RUN(nobsfilemax)                        = .true.
  logical :: OBSDA_OUT                                     = .false.
  character(len=H_LONG) :: OBS_IN_NAME(nobsfilemax)        = ''
  character(len=H_LONG) :: OBSDA_OUT_BASENAME              = ''
  character(len=H_LONG) :: OBSDA_MEAN_OUT_BASENAME         = ''
  character(len=H_LONG) :: OBSDA_MDET_OUT_BASENAME         = ''

  integer  :: SLOT_START                                   = 1
  integer  :: SLOT_END                                     = 1
  integer  :: SLOT_BASE                                    = 1
  real(RP) :: SLOT_TINTERVAL                               = 3600.0_RP       ! unit: [sec]

  logical  :: DEPARTURE_STAT_RADAR                         = .false.
  real(RP) :: DEPARTURE_STAT_T_RANGE                       = 0.0             ! time range within which observations are considered in the departure statistics.
                                                                             ! 0: no limit

  !--- PARAM_LETKF_VAR_LOCAL
  real(RP) :: VAR_LOCAL_UV(nv3d+nv2d)                      = 1.0_RP
  real(RP) :: VAR_LOCAL_T(nv3d+nv2d)                       = 1.0_RP
  real(RP) :: VAR_LOCAL_Q(nv3d+nv2d)                       = 1.0_RP
  real(RP) :: VAR_LOCAL_PS(nv3d+nv2d)                      = 1.0_RP
  real(RP) :: VAR_LOCAL_RAIN(nv3d+nv2d)                    = 1.0_RP
  real(RP) :: VAR_LOCAL_TC(nv3d+nv2d)                      = 1.0_RP
  real(RP) :: VAR_LOCAL_RADAR_REF(nv3d+nv2d)               = 1.0_RP
  real(RP) :: VAR_LOCAL_RADAR_VR(nv3d+nv2d)                = 1.0_RP

  character(len=H_LONG) :: INFL_MUL_IN_BASENAME            = ''
  character(len=H_LONG) :: INFL_MUL_OUT_BASENAME           = ''
  real(RP) :: INFL_MUL                                     = 1.0_RP          ! >  0: globally constant covariance inflation
                                                                             ! <= 0: use 3D inflation field from 'INFL_MUL_IN_BASENAME' file
  real(RP) :: INFL_MUL_MIN                                 = -1.0_RP         ! minimum inlfation factor (<= 0: not used)
  logical  :: INFL_MUL_ADAPTIVE                            = .false.         ! if true, outout adaptively estimated 3D inlfation field to 'INFL_MUL_OUT_BASENAME' file

  character(len=H_LONG) :: INFL_ADD_IN_BASENAME            = ''
  real(RP) :: INFL_ADD                                     = 0.0_RP          ! additive inflation
  logical  :: INFL_ADD_SHUFFLE                             = .false.         ! shuffle the additive inflation members?
  logical  :: INFL_ADD_Q_RATIO                             = .false.
  logical  :: INFL_ADD_REF_ONLY                            = .false.

  !--- PARAM_OBS_ERROR
  real(RP) :: OBSERR_U                                     = 1.0_RP
  real(RP) :: OBSERR_V                                     = 1.0_RP
  real(RP) :: OBSERR_T                                     = 1.0_RP
  real(RP) :: OBSERR_Q                                     = 0.001_RP
  real(RP) :: OBSERR_RH                                    = 10.0_RP
  real(RP) :: OBSERR_PS                                    = 100.0_RP        ! (Pa)
  real(RP) :: OBSERR_RADAR_REF                             = 5.0_RP
  real(RP) :: OBSERR_RADAR_VR                              = 3.0_RP
  real(RP) :: OBSERR_TCX                                   = 50.e+3_RP       ! (m)
  real(RP) :: OBSERR_TCY                                   = 50.e+3_RP       ! (m)
  real(RP) :: OBSERR_TCP                                   = 5.e+2_RP        ! (Pa)
  real(RP) :: OBSERR_PQ                                    = 0.001_RP        ! (kg/m3)

  ! PARAMETERS FOR RADAR DATA ASSIMILATION
  integer  :: NRADARTYPE                                   = 1               ! Currently PAWR (1) and LIDAR (2) ... not used?
  real(RP) :: RADAR_SO_SIZE_HORI                           = 1000.0_RP
  real(RP) :: RADAR_SO_SIZE_VERT                           = 1000.0_RP
  real(RP) :: RADAR_MAX_ABS_VR                             = 100.0_RP
  integer  :: RADAR_THIN_LETKF_METHOD                      = 0               ! Thinning method
                                                                             ! 0: No thinning
                                                                             ! 1: Nearest (2z)*(2x)^2 grids & their columns and
                                                                             ! rows + staggered grids
                                                                             ! x is given by RADAR_THIN_LETKF_HGRID_HNEAR
                                                                             ! z is given by RADAR_THIN_LETKF_HGRID_VNEAR
  integer  :: RADAR_THIN_LETKF_HGRID                       = 1               ! Horizontal thinning level in obs_local
  integer  :: RADAR_THIN_LETKF_VGRID                       = 1               ! Vertical thinning level in obs_local
  integer  :: RADAR_THIN_LETKF_HNEAR                       = 1               
  integer  :: RADAR_THIN_LETKF_VNEAR                       = 1               
  integer  :: RADAR_THIN_HORI                              = 1               ! Thinning horizontal interval (# of grids)
  integer  :: RADAR_THIN_VERT                              = 1               ! Thinning vertical interval (# of grids)
  logical  :: RADAR_USE_VR_STD                             = .true.          ! If we are going to use the wind std threshold within each box.
  logical  :: RADAR_BIAS_COR_RAIN                          = .false.         ! Simple bias correction for radar obs (rain)
  logical  :: RADAR_BIAS_COR_CLR                           = .false.         ! Simple bias correction for radar obs (clear sky)
  real(RP) :: RADAR_BIAS_RAIN_CONST_DBZ                    = 0.0_RP          ! Simply bias correction for radar obs (rain)
  real(RP) :: RADAR_BIAS_CLR_CONST_DBZ                     = 0.0_RP          ! Simply bias correction for radar obs (clear sky)
  real(RP) :: VR_STD_THRESHOLD                             = 2.5_RP          ! If wind variability within each superob is greather than this threshold the box is rejected.
  real(RP) :: ATTENUATION_THRESHOLD                        = 0.25_RP         ! 0.1 is 10dbz, 0.5 is aprox 5 dbz.

  logical  :: USE_RADAR_REF                                = .true.
  logical  :: USE_RADAR_VR                                 = .true.
  logical  :: USE_RADAR_PSEUDO_RH                          = .false.

  logical  :: USE_OBSERR_RADAR_REF                         = .false.
  logical  :: USE_OBSERR_RADAR_VR                          = .false.
  logical  :: RADAR_OBS_4D                                 = .false.
  logical  :: RADAR_PQV                                    = .false.         ! Pseudo qv DA for radar

  real(RP) :: RADAR_PQV_OMB                                = 25.0_RP         ! Threshold Obs-B for pseudo qv DA for radar
  real(RP) :: RADAR_REF_THRES_DBZ                          = 15.0_RP         ! Threshold of rain/no rain
  integer  :: MIN_RADAR_REF_MEMBER_OBSRAIN                 = 1               ! Minimum rainy ensemble members for assimilating rainy radar obs
  integer  :: MIN_RADAR_REF_MEMBER_OBSNORAIN               = 1               ! Minimum rainy ensemble members for assimilating clear-sky radar obs

  real(RP) :: MIN_RADAR_REF_DBZ                            = 0.0_RP          ! Minimum reflectivity
  real(RP) :: MIN_RADAR_REF_DBZ_VR                         = 5.0_RP          ! Minimum reflectivity (dBZ) for Doppler velocity observation
  real(RP) :: MIN_RADAR_REF_VR                             = 0.0_RP          ! Minimum reflectivity (Z) for Doppler velocity observation
  real(RP) :: LOW_REF_SHIFT                                = 0.0_RP

  real(RP) :: RADAR_ZMAX                                   =  99.e+3_RP      ! Height limit of radar data to be used
  real(RP) :: RADAR_ZMIN                                   = -99.e+3_RP      ! Height limit of radar data to be used
  real(RP) :: RADAR_PRH_ERROR                              =  0.1_RP         ! Obserational error for pseudo RH observations.

  integer  :: INTERPOLATION_TECHNIQUE                      = 1               ! These 2 flags affects the computation of model reflectivity and radial velocity. 
  integer  :: METHOD_REF_CALC                              = 3
  logical  :: USE_METHOD3_REF_MELT                         = .false.         ! Use radar operator considering melting (Xue et al. 2009QJRMS)
  logical  :: USE_T08_RS2014                               = .false.         ! Use RS2014 in snow obsope (must be consistent with SCALE) 
  logical  :: USE_TERMINAL_VELOCITY                        = .false.
  logical  :: USE_ATTENUATION                              = .true.          ! Consider attenuation in superobbing
  logical  :: USE_QCFLAG                                   = .true.          ! Consider or not qc flag.

  real(RP) :: RELAX_ALPHA                                  = 0.0_RP          ! RTPP relaxation parameter
  real(RP) :: RELAX_ALPHA_SPREAD                           = 0.0_RP          ! RTPS relaxation parameter
  logical  :: RELAX_TO_INFLATED_PRIOR                      = .false.         ! .true. : relaxation to multiplicatively inflated prior
                                                                             ! .false.: relaxation to original prior
  logical  :: RELAX_SPREAD_OUT                             = .false.
  character(len=H_LONG) :: RELAX_SPREAD_OUT_BASENAME       = ''

  real(RP) :: GROSS_ERROR                                  =  5.0_RP
  real(RP) :: GROSS_ERROR_RAIN                             = -1.0_RP         ! < 0: same as GROSS_ERROR
  real(RP) :: GROSS_ERROR_RADAR_REF                        = -1.0_RP         ! < 0: same as GROSS_ERROR
  real(RP) :: GROSS_ERROR_RADAR_VR                         = -1.0_RP         ! < 0: same as GROSS_ERROR
  real(RP) :: GROSS_ERROR_RADAR_PRH                        = -1.0_RP         ! < 0: same as GROSS_ERROR
  real(RP) :: GROSS_ERROR_TCX                              = -1.0_RP         ! debug ! < 0: same as GROSS_ERROR 
  real(RP) :: GROSS_ERROR_TCY                              = -1.0_RP         ! debug ! < 0: same as GROSS_ERROR
  real(RP) :: GROSS_ERROR_TCP                              = -1.0_RP         ! debug ! < 0: same as GROSS_ERROR

  real(RP) :: Q_UPDATE_TOP                                 =  0.0_RP         ! water vapor and hydrometeors are updated only below this pressure level (Pa)
  real(RP) :: Q_SPRD_MAX                                   = -1.0_RP         ! maximum q (ensemble spread)/(ensemble mean) (only effective when > 0)

  real(RP) :: BOUNDARY_BUFFER_WIDTH                        =   0.0_RP
  real(RP) :: PS_ADJUST_THRES                              = 100.0_RP

  logical :: NOBS_OUT                                      = .false.
  character(len=H_LONG) :: NOBS_OUT_BASENAME               = ''

  real(RP) :: HORI_LOCAL(nobtype) =                 &                        ! >0: localization length scale (m)
    (/  0.0_RP, -1.0_RP, -1.0_RP, -1.0_RP, -1.0_RP, &                        !  0: no localization XXX not implemented yet XXX
       -1.0_RP, -1.0_RP, -1.0_RP, -1.0_RP, -1.0_RP, &                        ! <0: same as HORI_LOCAL(1)
       -1.0_RP, -1.0_RP, -1.0_RP, -1.0_RP, -1.0_RP, &
       -1.0_RP, -1.0_RP, -1.0_RP, -1.0_RP, -1.0_RP, &
       -1.0_RP, -1.0_RP, -1.0_RP, -1.0_RP /)
  real(RP) :: VERT_LOCAL(nobtype) =                 &                        ! >0: localization length scale [ln(p) or m depends on obstype]
    (/  0.0_RP, -1.0_RP, -1.0_RP, -1.0_RP, -1.0_RP, &                        !  0: no localization
       -1.0_RP, -1.0_RP, -1.0_RP, -1.0_RP, -1.0_RP, &                        ! <0: same as VERT_LOCAL(1)
       -1.0_RP, -1.0_RP, -1.0_RP, -1.0_RP, -1.0_RP, &
       -1.0_RP, -1.0_RP, -1.0_RP, -1.0_RP, -1.0_RP, &
       -1.0_RP, -1.0_RP, -1.0_RP, -1.0_RP /)                                                       
  real(RP) :: TIME_LOCAL(nobtype) =                 &                        ! >0: localization length scale (sec) XXX not implemented yet XXX
    (/  0.0_RP, -1.0_RP, -1.0_RP, -1.0_RP, -1.0_RP, &                        !  0: no localization
       -1.0_RP, -1.0_RP, -1.0_RP, -1.0_RP, -1.0_RP, &                        ! <0: same as TIME_LOCAL(1)
       -1.0_RP, -1.0_RP, -1.0_RP, -1.0_RP, -1.0_RP, &
       -1.0_RP, -1.0_RP, -1.0_RP, -1.0_RP, -1.0_RP, &
       -1.0_RP, -1.0_RP, -1.0_RP, -1.0_RP /)

  real(RP) :: HORI_LOCAL_RADAR_OBSNOREF                    = -1.0_RP         ! <0: same as HORI_LOCAL(22=PHARAD)
  real(RP) :: HORI_LOCAL_RADAR_VR                          = -1.0_RP         ! <0: same as HORI_LOCAL(22=PHARAD)
  real(RP) :: VERT_LOCAL_RADAR_VR                          = -1.0_RP         ! <0: same as VERT_LOCAL(22=PHARAD)
  real(RP) :: VERT_LOCAL_RAIN_BASE                         = 85000.0_RP

  integer :: MAX_NOBS_PER_GRID(nobtype) =     &                              ! >0: observation number limit
    (/ 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, &                              !  0: do not limit observation numbers
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, &                              ! <0: same as MAX_NOBS_PER_GRID(1)
      -1, -1, -1, -1/)

  integer :: MAX_NOBS_PER_GRID_CRITERION                   = 1               ! 1: normalized 3D distance (from closest)
                                                                             ! 2: localization weight (from largest)
                                                                             ! 3: weighted observation error variance (from smallest)

  real(RP) :: OBS_MIN_SPACING(nobtype) =            &                        ! >0: typical minimum spacing of the obsetvation types in the densest observed area (not tuned carefully yet)
    (/  0.0_RP, -1.0_RP, -1.0_RP, -1.0_RP, -1.0_RP, &                        !     *this is only used for automatically determine OBS_SORT_GRID_SPACING. if using pre-set OBS_SORT_GRID_SPACING, this has no effect.
       -1.0_RP, -1.0_RP, -1.0_RP, -1.0_RP, -1.0_RP, &                        ! <=0: same as OBS_MIN_SPACING(1)
       -1.0_RP, -1.0_RP, -1.0_RP, -1.0_RP, -1.0_RP, &
       -1.0_RP, -1.0_RP, -1.0_RP, -1.0_RP, -1.0_RP, &
       -1.0_RP, -1.0_RP, -1.0_RP, -1.0_RP /)

  real(RP) :: OBS_SORT_GRID_SPACING(nobtype) =      &                        ! >0: optimal grid spacing for bucket sorting of observations
    (/  0.0_RP, -1.0_RP, -1.0_RP, -1.0_RP, -1.0_RP, &                        !  0: automatically determined based on HORI_LOCAL, MAX_NOBS_PER_GRID, and OBS_MIN_SPACING
       -1.0_RP, -1.0_RP, -1.0_RP, -1.0_RP, -1.0_RP, &                        ! <0: same as OBS_SORT_GRID_SPACING(1)
       -1.0_RP, -1.0_RP, -1.0_RP, -1.0_RP, -1.0_RP, &
       -1.0_RP, -1.0_RP, -1.0_RP, -1.0_RP, -1.0_RP, &
       -1.0_RP, -1.0_RP, -1.0_RP, -1.0_RP /)

  integer :: utime_obs(6) = &
    (/ -1, -1, -1, -1, -1, -1 /)

  integer :: mmean                                         = -99             ! use a value different from -1 to avoid equivalence to (my)rank_to_mem
  integer :: mmdet                                         = -99             ! when no member is corresponded to a rank/iteration
  integer :: mmdetin                                       = -99             ! 
  integer :: mmdetobs                                      = -99             ! 

  logical :: LETKF_ENTIRE_GRID_SEARCH_X                    = .false.         ! Gather all obs to analyze global constant parameters via ETKF (as in Kotsuki et al.)
  logical :: LETKF_ENTIRE_GRID_SEARCH_Y                    = .false.         !

  logical :: LETKF_DEBUG_LOG                               = .false.

  type(obs_info), allocatable :: obs(:)

  type(obs_da_value) :: obsda
  type(obs_da_value) :: obsda_sort                                           ! sorted obsda

  type(obs_grid_type), allocatable :: obsgrd(:)

  real(RP) :: MIN_RADAR_REF
  real(RP) :: RADAR_REF_THRES
  real(RP) :: RADAR_BIAS_RAIN_CONST
  real(RP) :: RADAR_BIAS_CLR_CONST

  real(RP), allocatable :: var_local(:,:)

  integer, allocatable :: n_merge(:)
  integer, allocatable :: ic_merge(:,:)

  integer :: nobstotalg
  integer :: nobstotal
  integer :: maxnobs_per_ctype

  ! combined obs type: {variable type (elm_u), report type (typ)}, allocated only when observations exist
  integer :: nctype                                                          ! number of combined obs type
  integer :: ctype_elmtyp(nid_obs,nobtype)                                   ! array of ctype for each combination of (elm_u, typ)
  integer,  allocatable :: elm_ctype(:)                                      ! array of elm  for each combined obs type
  integer,  allocatable :: elm_u_ctype(:)                                    ! array of elm_u for each combined obs type
  integer,  allocatable :: typ_ctype(:)                                      ! array of typ  for each combined obs type
  real(RP), allocatable :: hori_loc_ctype(:)                                 ! array of horizontal localization length for each combined obs type
  real(RP), allocatable :: vert_loc_ctype(:)                                 ! array of vertical localization length for each combined obs type

  ! observation monitor
  integer :: obsdep_nobs                                                     ! obsdep information
  integer,  allocatable :: obsdep_set(:)                                     !
  integer,  allocatable :: obsdep_idx(:)                                     !
  integer,  allocatable :: obsdep_qc (:)                                     !
  real(RP), allocatable :: obsdep_omb(:)                                     !
  real(RP), allocatable :: obsdep_oma(:)                                     !

  integer :: datatype
  integer :: n_merge_max
  logical :: radar_only

  integer :: nitmax
  integer :: nij1
  integer :: nij1max
  integer :: nproc_x
  integer :: nproc_y
  integer :: nens
  integer :: nensobs
  integer :: nmem
  integer :: nlon
  integer :: nlat
  integer :: nlev
  integer :: nlonh
  integer :: nlath
  integer :: nlevh
  integer :: nlong
  integer :: nlatg
  integer :: nlevall
  integer :: xhalo
  integer :: yhalo
  integer :: zhalo
  integer :: start_x
  integer :: end_x
  integer :: start_y
  integer :: end_y
  integer :: start_z
  integer :: end_z

  integer :: COMM_LCL
  integer :: NPRC_LCL
  integer :: RANK_LCL
  integer :: COMM_ENS
  integer :: NPRC_ENS
  integer :: RANK_ENS

  real(RP) :: DX
  real(RP) :: DY

  integer, allocatable :: nij1node(:)

  real(RP), allocatable :: rig1(:)
  real(RP), allocatable :: rjg1(:)
  real(RP), allocatable :: topo1(:)
  real(RP), allocatable :: hgt1(:,:)
  real(RP), allocatable :: v3dg(:,:,:,:)
  real(RP), allocatable :: v2dg(:,:,:)
  real(RP), allocatable :: v3d(:,:,:)
  real(RP), allocatable :: v2d(:,:)
  real(RP), allocatable :: topo2d(:,:)

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine LETKF_setup( &
      OBS_IN_NUM,       &
      ensemble_comm,    &
      ensemble_nprocs,  &
      ensemble_myrank,  &
      local_comm,       &
      local_nprocs,     &
      local_myrank,     &
      PRC_NUM_X,        &
      PRC_NUM_Y,        &
      KA, KS, KE,       &
      IA, IS, IE,       &
      JA, JS, JE,       &
      KMAX,             &
      IMAX,             &
      JMAX,             &
      KHALO,            &
      IHALO,            &
      JHALO,            &
      delta_x,          &
      delta_y,          &
      Zsfc              )
    use mpi
    use scale_prc, only: &
       PRC_abort
    use scale_prc_cartesC, only: &
       PRC_2Drank
    implicit none

    integer, intent(in) :: OBS_IN_NUM
    integer, intent(in) :: ensemble_comm
    integer, intent(in) :: ensemble_nprocs
    integer, intent(in) :: ensemble_myrank
    integer, intent(in) :: local_comm
    integer, intent(in) :: local_nprocs
    integer, intent(in) :: local_myrank
    integer, intent(in) :: PRC_NUM_X
    integer, intent(in) :: PRC_NUM_Y
    integer, intent(in) :: KA, KS, KE
    integer, intent(in) :: IA, IS, IE
    integer, intent(in) :: JA, JS, JE
    integer, intent(in) :: KMAX
    integer, intent(in) :: IMAX
    integer, intent(in) :: JMAX
    integer, intent(in) :: KHALO
    integer, intent(in) :: IHALO
    integer, intent(in) :: JHALO

    real(RP), intent(in) :: delta_x
    real(RP), intent(in) :: delta_y
    real(RP), intent(in) :: Zsfc(:,:)

    logical :: LETKF_DETERMINISTIC_RUN

    integer :: LETKF_MEM_NODES = 0   ! Number of nodes used for one member (0: automatically determined)

    namelist / PARAM_LETKF / &
      LETKF_DEBUG_LOG,                 &
      LETKF_DETERMINISTIC_RUN,         &
      LETKF_MEM_NODES,                 &
      SLOT_START,                      &
      SLOT_END,                        &
      SLOT_BASE,                       &
      SLOT_TINTERVAL,                  &
      DEPARTURE_STAT_RADAR,            &
      INFL_MUL,                        & 
      INFL_MUL_MIN,                    & 
      INFL_MUL_ADAPTIVE,               & 
      INFL_ADD,                        & 
      INFL_ADD_SHUFFLE,                & 
      INFL_ADD_Q_RATIO,                & 
      INFL_ADD_REF_ONLY,               & 
      RELAX_ALPHA,                     & 
      RELAX_ALPHA_SPREAD,              & 
      RELAX_TO_INFLATED_PRIOR,         & 
      GROSS_ERROR,                     & 
      GROSS_ERROR_RAIN,                & 
      GROSS_ERROR_RADAR_REF,           & 
      GROSS_ERROR_RADAR_VR,            & 
      GROSS_ERROR_RADAR_PRH,           & 
      GROSS_ERROR_TCX,                 & 
      GROSS_ERROR_TCY,                 & 
      GROSS_ERROR_TCP,                 & 
      Q_UPDATE_TOP,                    & 
      Q_SPRD_MAX,                      & 
      BOUNDARY_BUFFER_WIDTH,           & 
      HORI_LOCAL,                      &
      VERT_LOCAL,                      &
      TIME_LOCAL,                      &
      HORI_LOCAL_RADAR_OBSNOREF,       &
      HORI_LOCAL_RADAR_VR,             &
      VERT_LOCAL_RADAR_VR,             &
      MAX_NOBS_PER_GRID,               &
      MAX_NOBS_PER_GRID_CRITERION,     &
      OBS_MIN_SPACING,                 & 
      OBS_SORT_GRID_SPACING,           & 
      USE_RADAR_REF,                   & 
      USE_RADAR_VR,                    & 
      METHOD_REF_CALC,                 & 
      USE_TERMINAL_VELOCITY,           & 
      USE_OBSERR_RADAR_REF,            & 
      USE_OBSERR_RADAR_VR,             & 
      RADAR_REF_THRES_DBZ,             & 
      MIN_RADAR_REF_MEMBER_OBSRAIN,    & 
      MIN_RADAR_REF_MEMBER_OBSNORAIN,  & 
      MIN_RADAR_REF_DBZ,               & 
      MIN_RADAR_REF_DBZ_VR,            & 
      LOW_REF_SHIFT,                   & 
      RADAR_ZMAX,                      & 
      RADAR_ZMIN,                      & 
      RADAR_SO_SIZE_HORI,              & 
      RADAR_SO_SIZE_VERT,              & 
      RADAR_MAX_ABS_VR,                & 
      USE_METHOD3_REF_MELT,            & 
      RADAR_BIAS_COR_RAIN,             & 
      RADAR_BIAS_COR_CLR,              & 
      RADAR_BIAS_RAIN_CONST_DBZ,       & 
      RADAR_BIAS_CLR_CONST_DBZ,        & 
      RADAR_THIN_LETKF_METHOD,         & 
      RADAR_THIN_LETKF_HGRID,          & 
      RADAR_THIN_LETKF_VGRID,          & 
      RADAR_THIN_LETKF_HNEAR,          & 
      RADAR_THIN_LETKF_VNEAR,          & 
      RADAR_THIN_HORI,                 & 
      RADAR_THIN_VERT,                 & 
      OBSERR_U,                        &
      OBSERR_V,                        &
      OBSERR_T,                        &
      OBSERR_Q,                        &
      OBSERR_RH,                       &
      OBSERR_PS,                       &
      OBSERR_RADAR_REF,                & 
      OBSERR_RADAR_VR,                 &
      OBSERR_PQ,                       &
      LETKF_ENTIRE_GRID_SEARCH_X,      &
      LETKF_ENTIRE_GRID_SEARCH_Y

    integer :: n_mem
    integer :: n_mempn

    integer :: i, j, k, n
    integer :: ierr
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("LETKF_setup",*) 'Setup'

    if     ( RP == SP ) then
       datatype = MPI_REAL
    else if( RP == DP ) then
       datatype = MPI_DOUBLE_PRECISION
    else
       LOG_ERROR("obsope_tool",*) 'The precision has not been implemented yet:', RP
       call PRC_abort
    endif

    LETKF_DETERMINISTIC_RUN = .false.

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_LETKF,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("LETKF_setup",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("LETKF_setup",*) 'Not appropriate names in namelist PARAM_LETKF. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_LETKF)

    nproc_x = PRC_NUM_X
    nproc_y = PRC_NUM_Y
    nens    = ensemble_nprocs
    nensobs = ensemble_nprocs
    nmem    = ensemble_nprocs
    nlon    = IMAX
    nlat    = JMAX
    nlev    = KMAX
    nlonh   = IA
    nlath   = JA
    nlevh   = KA
    nlong   = IMAX * PRC_NUM_X
    nlatg   = JMAX * PRC_NUM_Y
    nlevall = nlev * nv3d + nv2d
    xhalo   = IHALO
    yhalo   = JHALO
    zhalo   = KHALO
    start_x = IS
    end_x   = IE
    start_y = JS
    end_y   = JE
    start_z = KS
    end_z   = KE
    mmean   = nens + 1

    COMM_ENS = ensemble_comm
    NPRC_ENS = ensemble_nprocs
    RANK_ENS = ensemble_myrank
    COMM_LCL = local_comm
    NPRC_LCL = local_nprocs
    RANK_LCL = local_myrank

    DX = delta_x
    DY = delta_y

    if( LETKF_DETERMINISTIC_RUN ) then ! set deterministic run
      ENS_WITH_MDET = .true.
      ! deterministic run is set to the last of ensemble member
      mmdet    = nens
      mmdetin  = nens
      mmdetobs = nens
      nmem     = nens - 1 ! member size except for deterministic run
    end if

    if( LETKF_MEM_NODES == 0 ) then
      LETKF_MEM_NODES = 1 ! (nprocs_m-1) / PPN + 1
    end if
    if( LETKF_MEM_NODES > 1 ) then
      n_mem   = NPRC_ENS / LETKF_MEM_NODES
      n_mempn = 1
    else
      n_mem   = NPRC_ENS
      n_mempn = 1 ! PPN / nprocs_m
    end if
    nitmax = ( NPRC_ENS - 1 ) / ( n_mem * n_mempn ) + 1

    i = mod( nlon*nlat, NPRC_ENS )
    nij1max = ( nlon*nlat - i ) / NPRC_ENS + 1
    if( RANK_ENS < i ) then
      nij1 = nij1max
    else
      nij1 = nij1max - 1
    end if

    allocate( obs( OBS_IN_NUM ) )
    allocate( nij1node( NPRC_ENS ) )

    do n = 1, NPRC_ENS
      if( n-1 < i ) then
        nij1node(n) = nij1max
      else
        nij1node(n) = nij1max - 1
      end if
    end do

    allocate( rig1 ( nij1       ) )
    allocate( rjg1 ( nij1       ) )
    allocate( topo1( nij1       ) )
    allocate( hgt1 ( nij1, nlev ) )

    allocate( v3dg( nlev, nlon, nlat, nv3d ) )
    allocate( v2dg(       nlon, nlat, nv2d ) )

    allocate( v3d ( nij1, nlev, nv3d ) )
    allocate( v2d ( nij1,       nv2d ) )

    do j = 1, JMAX
    do i = 1, IMAX
       v3dg(1,i,j,1) = real( i + PRC_2Drank(local_myrank,1) * IMAX + IHALO, kind=RP )
       v3dg(1,i,j,2) = real( j + PRC_2Drank(local_myrank,2) * JMAX + JHALO, kind=RP )
       v3dg(1,i,j,3) = Zsfc(i+IHALO,j+JHALO)
    end do
    end do

    call scatter_grd_mpi( mod( nens, n_mem*n_mempn ), v3dg, v2dg, v3d, v2d )

    rig1 (:) = v3d(:,1,1)
    rjg1 (:) = v3d(:,1,2)
    topo1(:) = v3d(:,1,3)

    call calc_z_grd( nij1, topo1, hgt1 )

    allocate( var_local( nv3d+nv2d, nid_obs_varlocal ) )

    if( RADAR_REF_THRES_DBZ < MIN_RADAR_REF_DBZ ) then
       RADAR_REF_THRES_DBZ = MIN_RADAR_REF_DBZ
    end if

    MIN_RADAR_REF         = 10.0_RP ** ( MIN_RADAR_REF_DBZ         / 10.0_RP )
    MIN_RADAR_REF_VR      = 10.0_RP ** ( MIN_RADAR_REF_DBZ_VR      / 10.0_RP )
    RADAR_REF_THRES       = 10.0_RP ** ( RADAR_REF_THRES_DBZ       / 10.0_RP )
    RADAR_BIAS_RAIN_CONST = 10.0_RP ** ( RADAR_BIAS_RAIN_CONST_DBZ / 10.0_RP )
    RADAR_BIAS_CLR_CONST  = 10.0_RP ** ( RADAR_BIAS_CLR_CONST_DBZ  / 10.0_RP )

    return
  end subroutine LETKF_setup

  !-----------------------------------------------------------------------------
  subroutine LETKF_finalize()
    implicit none

    deallocate( rig1  )
    deallocate( rjg1  )
    deallocate( topo1 )
    deallocate( hgt1  )
    deallocate( v3dg  )
    deallocate( v2dg  )
    deallocate( v3d   )
    deallocate( v2d   )

    deallocate( var_local )

    deallocate( obs )
    deallocate( nij1node )

    return
  end subroutine LETKF_finalize

  !-----------------------------------------------------------------------------
  subroutine LETKF_obs_readfile( &
      OBS_IN_NUM,      &
      OBS_IN_FORMAT,   &
      OBS_IN_BASENAME, &
      OBS_IN_MASKFILE  )
    use scale_prc, only: &
       PRC_abort
    use scale_time, only: &
       TIME_gettimelabel
    implicit none

    integer, intent(in) :: OBS_IN_NUM

    character(len=H_LONG), intent(in) :: OBS_IN_FORMAT(:)
    character(len=H_LONG), intent(in) :: OBS_IN_BASENAME(:)
    character(len=H_LONG), intent(in) :: OBS_IN_MASKFILE

    character(len=H_LONG) :: obsfile
    character(len=H_LONG) :: timelabel_obsfile

    logical :: err

    integer :: n
    integer :: ierr
    !---------------------------------------------------------------------------

    LOG_PROGRESS(*) 'data-assimilation / LETKF / obs / readfile'

    if( RANK_ENS == 0 .and. RANK_LCL == 0 ) then
       timelabel_obsfile = '_???????????????????.dat'
       call TIME_gettimelabel( timelabel_obsfile(2:20) )

       do n = 1, OBS_IN_NUM
          obsfile = trim(OBS_IN_BASENAME(n))//trim(timelabel_obsfile)

          if (OBS_IN_FORMAT(n) /= 'PAWR_TOSHIBA'    .and. &
              OBS_IN_FORMAT(n) /= 'MP_PAWR_TOSHIBA'       ) then
             inquire( file=obsfile, exist=err )
             if( .not. err ) then
                LOG_INFO("LETKF_obs_readfile",*) 'Warning: File (',trim(obsfile),') is not found. Skip.'
                ! skip process
                obs(n)%nobs = 0
                call obs_info_allocate(obs(n), extended=.true.)
                cycle
             end if
          end if

          select case( OBS_IN_FORMAT(n) )
          case ( 'PREPBUFR' )
            call get_nobs( obsfile, 8, obs(n)%nobs )
            call obs_info_allocate( obs(n), extended=.true. )
            call read_obs( obsfile, obs(n) )
          case ( 'RADAR' )
            call get_nobs_radar( obsfile, obs(n)%nobs, obs(n)%meta(1), obs(n)%meta(2), obs(n)%meta(3) )
            call obs_info_allocate( obs(n), extended=.true. )
            call read_obs_radar( obsfile, obs(n) )
          case ( 'PAWR_JRC' )
            LOG_ERROR("LETKF_obs_readfile",*) 'Error: This system has not been implemented yet. (OBS_IN_FORMAT(:) = PAWR_JRC)'
            call PRC_abort
          !  call read_obs_radar_jrc( obsfile, obs(n) )
          case ( 'PAWR_TOSHIBA' )
            call read_obs_radar_toshiba_pawr( obs(n), obsfile )
          case ( 'MP_PAWR_TOSHIBA' )
            call read_obs_radar_toshiba_mp_pawr( obs(n), obsfile, OBS_IN_MASKFILE )
          case default
            LOG_ERROR("LETKF_obs_readfile",*) 'Error: Unsupported observation file format.'
            call PRC_abort
          end select
       end do
    end if

    do n = 1, OBS_IN_NUM
       ! communicate obs. data to ensemble world in each local domain master at first
       if( RANK_LCL == 0 ) then
         call MPI_BCAST( obs(n)%nobs, 1, MPI_INTEGER, 0, COMM_ENS, ierr )

         if( RANK_ENS /= 0 ) then
           call obs_info_allocate( obs(n), extended=.true. )
         end if

         call MPI_BCAST( obs(n)%elm,  obs(n)%nobs,       MPI_INTEGER, 0, COMM_ENS, ierr )
         call MPI_BCAST( obs(n)%lon,  obs(n)%nobs,       datatype,    0, COMM_ENS, ierr )
         call MPI_BCAST( obs(n)%lat,  obs(n)%nobs,       datatype,    0, COMM_ENS, ierr )
         call MPI_BCAST( obs(n)%lev,  obs(n)%nobs,       datatype,    0, COMM_ENS, ierr )
         call MPI_BCAST( obs(n)%dat,  obs(n)%nobs,       datatype,    0, COMM_ENS, ierr )
         call MPI_BCAST( obs(n)%err,  obs(n)%nobs,       datatype,    0, COMM_ENS, ierr )
         call MPI_BCAST( obs(n)%typ,  obs(n)%nobs,       MPI_INTEGER, 0, COMM_ENS, ierr )
         call MPI_BCAST( obs(n)%dif,  obs(n)%nobs,       datatype,    0, COMM_ENS, ierr )
         call MPI_BCAST( obs(n)%meta, max_obs_info_meta, datatype,    0, COMM_ENS, ierr )
       end if

       ! broadcast obs. data to local domain
       call MPI_BCAST( obs(n)%nobs, 1, MPI_INTEGER, 0, COMM_LCL, ierr )

       if( RANK_LCL /= 0 ) then
         call obs_info_allocate( obs(n), extended=.true. )
       end if

       call MPI_BCAST( obs(n)%elm,  obs(n)%nobs,       MPI_INTEGER, 0, COMM_LCL, ierr )
       call MPI_BCAST( obs(n)%lon,  obs(n)%nobs,       datatype,    0, COMM_LCL, ierr )
       call MPI_BCAST( obs(n)%lat,  obs(n)%nobs,       datatype,    0, COMM_LCL, ierr )
       call MPI_BCAST( obs(n)%lev,  obs(n)%nobs,       datatype,    0, COMM_LCL, ierr )
       call MPI_BCAST( obs(n)%dat,  obs(n)%nobs,       datatype,    0, COMM_LCL, ierr )
       call MPI_BCAST( obs(n)%err,  obs(n)%nobs,       datatype,    0, COMM_LCL, ierr )
       call MPI_BCAST( obs(n)%typ,  obs(n)%nobs,       MPI_INTEGER, 0, COMM_LCL, ierr )
       call MPI_BCAST( obs(n)%dif,  obs(n)%nobs,       datatype,    0, COMM_LCL, ierr )
       call MPI_BCAST( obs(n)%meta, max_obs_info_meta, datatype,    0, COMM_LCL, ierr )

       LOG_INFO('LETKF_obs_readfile',*) 'observations input:', obs(n)%nobs
    end do

    return
  end subroutine LETKF_obs_readfile

  subroutine LETKF_obs_clear( OBS_IN_NUM )
    implicit none

    integer, intent(in) :: OBS_IN_NUM

    integer :: n
    !---------------------------------------------------------------------------

    LOG_PROGRESS(*) 'data-assimilation / LETKF / cleaning'

    do n = 1, OBS_IN_NUM
       call obs_info_deallocate( obs(n) )
    end do

    call obs_da_value_deallocate( obsda_sort )

    return
  end subroutine LETKF_obs_clear

  !-----------------------------------------------------------------------
  ! Observation operator calculation
  !-----------------------------------------------------------------------
  subroutine LETKF_obs_operator( &
      OBS_IN_NUM,    & ! [IN]
      OBS_IN_FORMAT, & ! [IN]
      U,             & ! [IN]
      V,             & ! [IN]
      W,             & ! [IN]
      TEMP,          & ! [IN]
      PRES,          & ! [IN]
      QV,            & ! [IN]
      QC,            & ! [IN]
      QR,            & ! [IN]
      QI,            & ! [IN]
      QS,            & ! [IN]
      QG,            & ! [IN]
      RH,            & ! [IN]
      HGT,           & ! [IN]
      TOPO,          & ! [IN]
      PS,            & ! [IN]
      RAIN,          & ! [IN]
      U10M,          & ! [IN]
      V10M,          & ! [IN]
      T2M,           & ! [IN]
      Q2M,           & ! [IN]
      nobs_extern    ) ! [IN]
    implicit none

    integer, intent(in) :: OBS_IN_NUM
    character(len=H_LONG), intent(in) :: OBS_IN_FORMAT(:)

    real(RP), intent(in) :: U   (nlevh,nlonh,nlath)
    real(RP), intent(in) :: V   (nlevh,nlonh,nlath)
    real(RP), intent(in) :: W   (nlevh,nlonh,nlath)
    real(RP), intent(in) :: TEMP(nlevh,nlonh,nlath)
    real(RP), intent(in) :: PRES(nlevh,nlonh,nlath)
    real(RP), intent(in) :: QV  (nlevh,nlonh,nlath)
    real(RP), intent(in) :: QC  (nlevh,nlonh,nlath)
    real(RP), intent(in) :: QR  (nlevh,nlonh,nlath)
    real(RP), intent(in) :: QI  (nlevh,nlonh,nlath)
    real(RP), intent(in) :: QS  (nlevh,nlonh,nlath)
    real(RP), intent(in) :: QG  (nlevh,nlonh,nlath)
    real(RP), intent(in) :: RH  (nlevh,nlonh,nlath)
    real(RP), intent(in) :: HGT (nlevh,nlonh,nlath)

    real(RP), intent(in) :: TOPO(nlonh,nlath)
    real(RP), intent(in) :: PS  (nlonh,nlath)
    real(RP), intent(in) :: RAIN(nlonh,nlath)
    real(RP), intent(in) :: U10M(nlonh,nlath)
    real(RP), intent(in) :: V10M(nlonh,nlath)
    real(RP), intent(in) :: T2M (nlonh,nlath)
    real(RP), intent(in) :: Q2M (nlonh,nlath)

    integer, optional, intent(in) :: nobs_extern

    type(obs_da_value) :: obsda_tmp

    integer :: i, j, k
    integer :: it, im, iof, islot, ierr
    integer :: n, nn, nsub, nmod, n1, n2

    integer :: nobs     ! observation number processed in this subroutine
    integer :: nobs_all
    integer :: nobs_max_per_file
    integer :: nobs_max_per_file_sub
    integer :: slot_nobsg

    integer :: ip, ibufs
    integer, allocatable :: cntr(:), dspr(:)
    integer, allocatable :: cnts(:), dsps(:)
    integer, allocatable :: bsn(:,:), bsna(:,:), bsnext(:,:)
    integer :: islot_time_out, islot_domain_out

    integer, allocatable :: obrank_bufs(:)
    real(RP), allocatable :: ri_bufs(:)
    real(RP), allocatable :: rj_bufs(:)

    integer, allocatable :: obset_bufs(:)
    integer, allocatable :: obidx_bufs(:)

    integer :: slot_id(SLOT_START:SLOT_END)
    real(RP) :: slot_lb(SLOT_START:SLOT_END)
    real(RP) :: slot_ub(SLOT_START:SLOT_END)

    real(RP), allocatable :: v3dg(:,:,:,:)
    real(RP), allocatable :: v2dg(:,:,:)

    real(RP) :: ril, rjl, rk, rkz

    character(len=4) :: nstr
    !-------------------------------------------------------------------------------

    LOG_PROGRESS(*) 'data-assimilation / LETKF / obs / operator'

    !-------------------------------------------------------------------------------
    ! First scan of all observation data: Compute their horizontal location and time
    !-------------------------------------------------------------------------------

    nobs_all = 0
    nobs_max_per_file = 0
    do iof = 1, OBS_IN_NUM
      if (obs(iof)%nobs > nobs_max_per_file) then
        nobs_max_per_file = obs(iof)%nobs
      end if
      if (OBSDA_RUN(iof)) then
        nobs_all = nobs_all + obs(iof)%nobs
      end if
    end do

    nobs_max_per_file_sub = (nobs_max_per_file - 1) / NPRC_LCL + 1
    allocate( obrank_bufs(nobs_max_per_file_sub) )
    allocate( ri_bufs(nobs_max_per_file_sub) )
    allocate( rj_bufs(nobs_max_per_file_sub) )

    allocate( cntr(NPRC_LCL) )
    allocate( dspr(NPRC_LCL) )

    ! Use all processes to compute the basic obsevration information
    ! (locations in model grids and the subdomains they belong to)
    !-----------------------------------------------------------------------------

    do iof = 1, OBS_IN_NUM
      if( obs(iof)%nobs > 0 ) then ! Process basic obsevration information for all observations since this information is not saved in obsda files
                                   ! when using separate observation operators; ignore the 'OBSDA_RUN' setting for this section
        nsub = obs(iof)%nobs / NPRC_LCL
        nmod = mod( obs(iof)%nobs, NPRC_LCL )
        do ip = 1, nmod
          cntr(ip) = nsub + 1
        end do
        do ip = nmod+1, NPRC_LCL
          cntr(ip) = nsub
        end do
        dspr(1) = 0
        do ip = 2, NPRC_LCL
          dspr(ip) = dspr(ip-1) + cntr(ip-1)
        end do

        obrank_bufs(:) = -1
        do ibufs = 1, cntr(RANK_LCL+1)
          n = dspr(RANK_LCL+1) + ibufs

          call phys2ij(obs(iof)%lon(n), obs(iof)%lat(n), ri_bufs(ibufs), rj_bufs(ibufs))
          call rij_rank(ri_bufs(ibufs), rj_bufs(ibufs), obrank_bufs(ibufs))
        end do

        call MPI_ALLGATHERV( obrank_bufs, cntr(RANK_LCL+1), MPI_INTEGER, obs(iof)%rank, cntr, dspr, MPI_INTEGER, COMM_LCL, ierr )
        call MPI_ALLGATHERV( ri_bufs,     cntr(RANK_LCL+1), datatype,    obs(iof)%ri,   cntr, dspr, datatype,    COMM_LCL, ierr )
        call MPI_ALLGATHERV( rj_bufs,     cntr(RANK_LCL+1), datatype,    obs(iof)%rj,   cntr, dspr, datatype,    COMM_LCL, ierr )

      end if
    end do

    deallocate(cntr, dspr)
    deallocate(obrank_bufs, ri_bufs, rj_bufs)

    ! Bucket sort of observation wrt. time slots and subdomains using the process rank 0
    !-----------------------------------------------------------------------------

    islot_time_out   = SLOT_END + 1 ! slot = SLOT_END+1 for observation not in the assimilation time window
    islot_domain_out = SLOT_END + 2 ! slot = SLOT_END+2 for observation outside of the model domain

    allocate( bsn ( SLOT_START  :SLOT_END+2, 0:NPRC_LCL-1 ) )
    allocate( bsna( SLOT_START-1:SLOT_END+2, 0:NPRC_LCL-1 ) )

    if (RANK_ENS == 0) then
      allocate ( obset_bufs(nobs_all) )
      allocate ( obidx_bufs(nobs_all) )
    end if

    if (RANK_ENS == 0 .and. RANK_LCL == 0) then
      allocate( bsnext( SLOT_START:SLOT_END+2, 0:NPRC_LCL-1 ) )
      bsn   (:,:) = 0
      bsna  (:,:) = 0
      bsnext(:,:) = 0

      do iof = 1, OBS_IN_NUM
        if (OBSDA_RUN(iof) .and. obs(iof)%nobs > 0) then
          do n = 1, obs(iof)%nobs
            if (obs(iof)%rank(n) == -1) then
              ! process the observations outside of the model domain in process rank 0
              bsn(islot_domain_out, 0) = bsn(islot_domain_out, 0) + 1
            else
              islot = ceiling(obs(iof)%dif(n) / SLOT_TINTERVAL - 0.5d0) + SLOT_BASE
              if (islot < SLOT_START .or. islot > SLOT_END) then
                islot = islot_time_out
              end if
              bsn(islot, obs(iof)%rank(n)) = bsn(islot, obs(iof)%rank(n)) + 1
            end if
          end do
        end if
      end do

      do ip = 0, NPRC_LCL-1
        if (ip > 0) then
          bsna(SLOT_START-1, ip) = bsna(SLOT_END+2, ip-1)
        end if
        do islot = SLOT_START, SLOT_END+2
          bsna(islot, ip) = bsna(islot-1, ip) + bsn(islot, ip)
        end do
        bsnext(SLOT_START:SLOT_END+2, ip) = bsna(SLOT_START-1:SLOT_END+1, ip)
      end do

      do iof = 1, OBS_IN_NUM
        if (OBSDA_RUN(iof) .and. obs(iof)%nobs > 0) then
          do n = 1, obs(iof)%nobs
            if (obs(iof)%rank(n) == -1) then
              ! process the observations outside of the model domain in process rank 0
              bsnext(islot_domain_out, 0) = bsnext(islot_domain_out, 0) + 1
              obset_bufs(bsnext(islot_domain_out, 0)) = iof
              obidx_bufs(bsnext(islot_domain_out, 0)) = n
            else
              islot = ceiling(obs(iof)%dif(n) / SLOT_TINTERVAL - 0.5d0) + SLOT_BASE
              if (islot < SLOT_START .or. islot > SLOT_END) then
                islot = islot_time_out
              end if
              bsnext(islot, obs(iof)%rank(n)) = bsnext(islot, obs(iof)%rank(n)) + 1
              obset_bufs(bsnext(islot, obs(iof)%rank(n))) = iof
              obidx_bufs(bsnext(islot, obs(iof)%rank(n))) = n
            end if
          end do
        end if
      end do

      deallocate( bsnext )

    end if

    ! Broadcast the bucket-sort observation numbers to all processes and print
    !-----------------------------------------------------------------------------

    if( RANK_LCL == 0 ) then
      call MPI_BCAST( bsn,  (SLOT_END-SLOT_START+3)*NPRC_LCL, MPI_INTEGER, 0, COMM_ENS, ierr )
      call MPI_BCAST( bsna, (SLOT_END-SLOT_START+4)*NPRC_LCL, MPI_INTEGER, 0, COMM_ENS, ierr )
    end if
    call MPI_BCAST( bsn,  (SLOT_END-SLOT_START+3)*NPRC_LCL, MPI_INTEGER, 0, COMM_LCL, ierr )
    call MPI_BCAST( bsna, (SLOT_END-SLOT_START+4)*NPRC_LCL, MPI_INTEGER, 0, COMM_LCL, ierr )

    do islot = SLOT_START, SLOT_END
      slot_id(islot) = islot - SLOT_START + 1
      slot_lb(islot) = (real(islot - SLOT_BASE, RP) - 0.5d0) * SLOT_TINTERVAL
      slot_ub(islot) = (real(islot - SLOT_BASE, RP) + 0.5d0) * SLOT_TINTERVAL
    end do

    ! Scatter the basic obsevration information to processes group {myrank_e = 0},
    ! each of which only gets the data in its own subdomain
    !-----------------------------------------------------------------------------

    nobs = bsna( SLOT_END+2, RANK_LCL ) - bsna( SLOT_START-1, RANK_LCL )

    obsda_tmp%nobs = nobs
    call obs_da_value_allocate(obsda_tmp, 0)

    if (present(nobs_extern)) then
      obsda%nobs = nobs + nobs_extern
    else
      obsda%nobs = nobs
    end if
    call obs_da_value_allocate(obsda, nitmax)

    if (RANK_ENS == 0) then
      allocate (cnts(NPRC_LCL))
      allocate (dsps(NPRC_LCL))
      do ip = 0, NPRC_LCL-1
        dsps(ip+1) = bsna(SLOT_START-1, ip)
        cnts(ip+1) = bsna(SLOT_END+2, ip) - dsps(ip+1)
      end do

      call MPI_SCATTERV(obset_bufs, cnts, dsps, MPI_INTEGER, obsda_tmp%set, cnts(RANK_LCL+1), MPI_INTEGER, 0, COMM_LCL, ierr)
      call MPI_SCATTERV(obidx_bufs, cnts, dsps, MPI_INTEGER, obsda_tmp%idx, cnts(RANK_LCL+1), MPI_INTEGER, 0, COMM_LCL, ierr)

      deallocate (cnts, dsps)
      deallocate (obset_bufs, obidx_bufs)
    end if

    ! Broadcast the basic obsevration information
    ! from processes group {myrank_e = 0} to all processes
    !-----------------------------------------------------------------------------

    call MPI_BCAST(obsda_tmp%set, nobs, MPI_INTEGER, 0, COMM_ENS, ierr)
    call MPI_BCAST(obsda_tmp%idx, nobs, MPI_INTEGER, 0, COMM_ENS, ierr)

    obsda%set(1:nobs) = obsda_tmp%set
    obsda%idx(1:nobs) = obsda_tmp%idx

    !-------------------------------------------------------------------------------
    ! Second scan of observation data in own subdomain: Compute H(x), QC, ... etc.
    !-------------------------------------------------------------------------------

    allocate ( v3dg(nlevh,nlonh,nlath,nv3dd) )
    allocate ( v2dg(      nlonh,nlath,nv2dd) )

    do j = 1, nlath
    do i = 1, nlonh
    do k = 1, nlevh
       v3dg(k,i,j,iv3dd_u  ) = U   (k,i,j)
       v3dg(k,i,j,iv3dd_v  ) = V   (k,i,j)
       v3dg(k,i,j,iv3dd_w  ) = W   (k,i,j)
       v3dg(k,i,j,iv3dd_t  ) = TEMP(k,i,j)
       v3dg(k,i,j,iv3dd_p  ) = PRES(k,i,j)
       v3dg(k,i,j,iv3dd_q  ) = QV  (k,i,j)
       v3dg(k,i,j,iv3dd_qc ) = QC  (k,i,j)
       v3dg(k,i,j,iv3dd_qr ) = QR  (k,i,j)
       v3dg(k,i,j,iv3dd_qi ) = QI  (k,i,j)
       v3dg(k,i,j,iv3dd_qs ) = QS  (k,i,j)
       v3dg(k,i,j,iv3dd_qg ) = QG  (k,i,j)
       v3dg(k,i,j,iv3dd_rh ) = RH  (k,i,j)
       v3dg(k,i,j,iv3dd_hgt) = HGT (k,i,j)
    end do
    end do
    end do
    do j = 1, nlath
    do i = 1, nlonh
       v2dg(i,j,iv2dd_topo) = TOPO(i,j)
       v2dg(i,j,iv2dd_ps  ) = PS  (i,j)
       v2dg(i,j,iv2dd_rain) = RAIN(i,j)
       v2dg(i,j,iv2dd_u10m) = U10M(i,j)
       v2dg(i,j,iv2dd_v10m) = V10M(i,j)
       v2dg(i,j,iv2dd_t2m ) = T2M (i,j)
       v2dg(i,j,iv2dd_q2m ) = Q2M (i,j)
    end do
    end do

    do it = 1, nitmax
      if (nobs > 0) then
        obsda_tmp%qc(1:nobs) = iqc_undef
      end if

      ! Observations not in the assimilation time window
      ! 
      n1 = bsna(islot_time_out-1, RANK_LCL) - bsna(SLOT_START-1, RANK_LCL) + 1
      n2 = bsna(islot_time_out,   RANK_LCL) - bsna(SLOT_START-1, RANK_LCL)
      if (n1 <= n2) then
        obsda_tmp%qc(n1:n2) = iqc_time
      end if

      ! Observations outside of the model domain
      ! 
      n1 = bsna(islot_domain_out-1, RANK_LCL) - bsna(SLOT_START-1, RANK_LCL) + 1
      n2 = bsna(islot_domain_out,   RANK_LCL) - bsna(SLOT_START-1, RANK_LCL)
      if (n1 <= n2) then
        obsda_tmp%qc(n1:n2) = iqc_out_h
      end if

      ! Valid observations: loop over time slots
      ! 
      do islot = SLOT_START, SLOT_END
        n1 = bsna(islot-1, RANK_LCL) - bsna(SLOT_START-1, RANK_LCL) + 1
        n2 = bsna(islot,   RANK_LCL) - bsna(SLOT_START-1, RANK_LCL)
        slot_nobsg = sum(bsn(islot, :))

        if (slot_nobsg <= 0) then
          cycle
        end if

        !call read_history( filename, it, islot, v3dg, v2dg )

        do nn = n1, n2
          iof = obsda_tmp%set(nn)
          n = obsda_tmp%idx(nn)

          call rij_g2l(RANK_LCL, obs(iof)%ri(n), obs(iof)%rj(n), ril, rjl)

          if (.not. USE_OBS(obs(iof)%typ(n))) then
            obsda_tmp%qc(nn) = iqc_otype
            cycle
          end if

          select case (OBS_IN_FORMAT(iof))
          !=====================================================================
          case ( 'PREPBUFR' )
          !---------------------------------------------------------------------
            call phys2ijk(v3dg(:,:,:,iv3dd_p), obs(iof)%elm(n), ril, rjl, obs(iof)%lev(n), rk, obsda_tmp%qc(nn))
            if (obsda_tmp%qc(nn) == iqc_good) then
              call Trans_XtoY(obs(iof)%elm(n), ril, rjl, rk, &
                              obs(iof)%lon(n), obs(iof)%lat(n), v3dg, v2dg, obsda_tmp%val(nn), obsda_tmp%qc(nn))
            end if
          !=====================================================================
          case ( 'RADAR', 'PAWR_TOSHIBA', 'MP_PAWR_TOSHIBA', 'PAWR_JRC', 'HIMAWARI8' )
          !---------------------------------------------------------------------
            if ( obs(iof)%lev(n) > RADAR_ZMAX .or. obs(iof)%lev(n) < RADAR_ZMIN ) then
              obsda_tmp%qc(nn) = iqc_radar_vhi
            else
              call phys2ijkz(v3dg(:,:,:,iv3dd_hgt), ril, rjl, obs(iof)%lev(n), rkz, obsda_tmp%qc(nn))
            end if
            if (obsda_tmp%qc(nn) == iqc_good) then
              call Trans_XtoY_radar(obs(iof)%elm(n), obs(iof)%meta(1), obs(iof)%meta(2), obs(iof)%meta(3), ril, rjl, rkz, &
                                    obs(iof)%lon(n), obs(iof)%lat(n), obs(iof)%lev(n), v3dg, v2dg, obsda_tmp%val(nn), obsda_tmp%qc(nn))
              !if (obsda_tmp%qc(nn) == iqc_ref_low) obsda_tmp%qc(nn) = iqc_good ! when process the observation operator, we don't care if reflectivity is too small

              call itpl_3d( v3dg(:,:,:,iv3dd_p), rkz, ril, rjl, obsda_tmp%pm(nn) )
              call itpl_3d( v3dg(:,:,:,iv3dd_t), rkz, ril, rjl, obsda_tmp%tm(nn) )
              call itpl_3d( v3dg(:,:,:,iv3dd_q), rkz, ril, rjl, obsda_tmp%qv(nn) )

            end if
          end select
        end do
      end do

      !!! (tentative) avoid a bug with GCC 5-8 compiler in L1421 (process to treat iqc_ref_low as iqc_good) !!! 
      where( obsda_tmp%qc(:) == iqc_ref_low )
        obsda_tmp%qc(:) = iqc_good
      end where

      ! Prepare variables that will need to be communicated if obsda_return is given
      ! 
      call obs_da_value_partial_reduce_iter(obsda, it, 1, nobs, obsda_tmp%val, obsda_tmp%qc, &
                                            obsda_tmp%qv, obsda_tmp%tm, obsda_tmp%pm )
    end do

    deallocate ( v3dg, v2dg )
    deallocate ( bsn, bsna )

    call obs_da_value_deallocate( obsda_tmp )

    return
  end subroutine LETKF_obs_operator

  subroutine LETKF_obs_initialize( OBS_IN_NUM, nobs_extern )
    use scale_prc, only: &
       PRC_abort
    use scale_const, only: &
      UNDEF => CONST_UNDEF, &
      TEM00 => CONST_TEM00
    IMPLICIT NONE

    integer, intent(in) :: OBS_IN_NUM
    integer, optional, intent(in) :: nobs_extern

    INTEGER :: n,i,j,ierr,im,iof,iidx

    integer :: n1, n2

    integer :: mem_ref
    real(RP) :: qvs, qdry

    integer :: it,ip
    integer :: ityp,ielm,ielm_u,ictype
    real(RP) :: target_grdspc

    integer :: myp_i,myp_j
    integer :: ip_i,ip_j

    integer :: nobs_sub(n_qc_steps),nobs_g(n_qc_steps)

    integer :: nobs_elms(nid_obs)
    integer :: nobs_elms_sum(nid_obs)

    integer :: nobs_intern
    integer :: nobs_extern_


    character(len=3) :: use_obs_print
    character(4) :: nstr

    integer :: cnts
    integer :: cntr(NPRC_LCL)
    integer :: dspr(NPRC_LCL)
    integer :: nensobs_div, nensobs_mod
    integer :: im_obs_1, im_obs_2, nensobs_part

    integer :: ns_ext, ne_ext, ns_bufr, ne_bufr
    integer :: ishift, jshift

    type(obs_da_value) :: obsbufs, obsbufr
    integer :: imin1,imax1,jmin1,jmax1,imin2,imax2,jmin2,jmax2

    real(RP),allocatable :: tmpelm(:)
    integer :: monit_nobs(nid_obs)
    real(RP) :: bias(nid_obs)
    real(RP) :: rmse(nid_obs)

    type(obs_da_value) :: obsda_ext
    logical :: ctype_use(nid_obs,nobtype)
    !-------------------------------------------------------------------------------

    LOG_PROGRESS(*) 'data-assimilation / LETKF / obs / initialize'

    if( present(nobs_extern) ) then
       nobs_extern_ = nobs_extern
    else
       nobs_extern_ = 0
    endif

    nobs_intern = obsda%nobs - nobs_extern_
    if( LETKF_DEBUG_LOG ) then
       LOG_INFO("LETKF_debug",'(1x,A,I8)') 'Internally processed observations: ', nobs_intern
       LOG_INFO("LETKF_debug",'(1x,A,I8)') 'Externally processed observations: ', nobs_extern_
       LOG_INFO("LETKF_debug",'(1x,A,I8)') 'Total                observations: ', obsda%nobs
    endif

    !-------------------------------------------------------------------------------
    ! Read externally processed observations
    !-------------------------------------------------------------------------------

    !if( nobs_extern > 0 ) then
    !  n1 = nobs_intern + 1
    !  n2 = obsda%nobs

    !  obsda_ext%nobs = nobs_extern
    !  call obs_da_value_allocate(obsda_ext,0)

    !  do it = 1, nitmax
    !    im = myrank_to_mem(it)
    !    if ((im >= 1 .and. im <= MEMBER) .or. im == mmdetin) then
    !      if (im <= MEMBER) then
    !        obsdafile = OBSDA_IN_BASENAME
    !        call filename_replace_mem(obsdafile, im)
    !      else if (im == mmean) then
    !        obsdafile = OBSDA_MEAN_IN_BASENAME
    !      else if (im == mmdet) then
    !        obsdafile = OBSDA_MDET_IN_BASENAME
    !      end if
    !      call read_obs_da(trim(obsdafile)//obsda_suffix,obsda_ext,0)

    !      if (OBSDA_OUT) then
    !        if (im <= MEMBER) then
    !          obsdafile = OBSDA_OUT_BASENAME
    !          call filename_replace_mem(obsdafile, im)
    !        else if (im == mmean) then
    !          obsdafile = OBSDA_MEAN_OUT_BASENAME
    !        else if (im == mmdet) then
    !          obsdafile = OBSDA_MDET_OUT_BASENAME
    !        end if
    !        call write_obs_da(trim(obsdafile)//obsda_suffix,obsda_ext,0,append=.true.)
    !      end if

    !      ! variables without an ensemble dimension
    !      if (it == 1) then
    !        obsda%set(n1:n2) = obsda_ext%set
    !        obsda%idx(n1:n2) = obsda_ext%idx
    !      end if

    !      call obs_da_value_partial_reduce_iter(obsda, it, n1, n2, obsda_ext%val, obsda_ext%qc, &
    !                                            obsda_ext%qv, obsda_ext%tm, obsda_ext%pm )

    !    end if ! [ (im >= 1 .and. im <= MEMBER) .or. im == mmdetin ]
    !  end do ! [ it = 1, nitmax ]

    !  call obs_da_value_deallocate(obsda_ext)

    !  ! Broadcast the observation information shared by members (e.g., grid numbers)
    !  !---------------------------------------------------------------------------

    !  if (nprocs_e > MEMBER) then
    !    call MPI_BCAST(obsda%set(n1:n2), nobs_extern, MPI_INTEGER, 0, MPI_COMM_e, ierr)
    !    call MPI_BCAST(obsda%idx(n1:n2), nobs_extern, MPI_INTEGER, 0, MPI_COMM_e, ierr)
    !  end if
    !end if

    !-------------------------------------------------------------------------------
    ! Allreduce externally processed observations
    !---------------------------------------------------------------------------

    call obs_da_value_allreduce( obsda )

    !-------------------------------------------------------------------------------
    ! Process observations and quality control (QC)
    !-------------------------------------------------------------------------------

    ! Pre-process data
    !-----------------------------------------------------------------------------

    ctype_use(:,:) = .false.
    do iof = 1, OBS_IN_NUM
    do n = 1, obs(iof)%nobs
      select case (obs(iof)%elm(n))
      case (id_radar_ref_obs)
        if (obs(iof)%dat(n) >= 0.0d0 .and. obs(iof)%dat(n) < 1.0d10) then
          if (obs(iof)%dat(n) < MIN_RADAR_REF) then
            obs(iof)%elm(n) = id_radar_ref_zero_obs
            obs(iof)%dat(n) = MIN_RADAR_REF_DBZ + LOW_REF_SHIFT
          else
            obs(iof)%dat(n) = 10.0d0 * log10(obs(iof)%dat(n))
          end if
        else
          obs(iof)%dat(n) = undef
        end if
        if (USE_OBSERR_RADAR_REF) then
          obs(iof)%err(n) = OBSERR_RADAR_REF
        end if
      case (id_radar_ref_zero_obs)
        obs(iof)%dat(n) = MIN_RADAR_REF_DBZ + LOW_REF_SHIFT
        if (USE_OBSERR_RADAR_REF) then
          obs(iof)%err(n) = OBSERR_RADAR_REF
        end if
      case (id_radar_vr_obs)
        if (USE_OBSERR_RADAR_VR) then
          obs(iof)%err(n) = OBSERR_RADAR_VR
        end if
      end select

      ! mark (elm, typ) combinations for which observations exist
      ctype_use(uid_obs(obs(iof)%elm(n)), obs(iof)%typ(n)) = .true.
    end do
    end do

    ! do this outside of the above obs loop, so these (ctype) arrays can be in ascending order
    nctype = count(ctype_use)
    if (allocated(elm_ctype     )) deallocate(elm_ctype     )
    if (allocated(elm_u_ctype   )) deallocate(elm_u_ctype   )
    if (allocated(typ_ctype     )) deallocate(typ_ctype     )
    if (allocated(hori_loc_ctype)) deallocate(hori_loc_ctype)
    if (allocated(vert_loc_ctype)) deallocate(vert_loc_ctype)
    allocate (elm_ctype     (nctype))
    allocate (elm_u_ctype   (nctype))
    allocate (typ_ctype     (nctype))
    allocate (hori_loc_ctype(nctype))
    allocate (vert_loc_ctype(nctype))
    ictype = 0
    ctype_elmtyp(:,:) = 0
    do ityp = 1, nobtype
    do ielm_u = 1, nid_obs
      if (ctype_use(ielm_u, ityp)) then
        ictype = ictype + 1
        ctype_elmtyp(ielm_u, ityp) = ictype

        elm_ctype(ictype) = elem_uid(ielm_u)
        elm_u_ctype(ictype) = ielm_u
        typ_ctype(ictype) = ityp

        ! horizontal localization
        if (elm_ctype(ictype) == id_radar_ref_zero_obs) then
          hori_loc_ctype(ictype) = HORI_LOCAL_RADAR_OBSNOREF
        else if (elm_ctype(ictype) == id_radar_vr_obs) then
          hori_loc_ctype(ictype) = HORI_LOCAL_RADAR_VR
        else
          hori_loc_ctype(ictype) = HORI_LOCAL(ityp)
        end if
        ! vertical localization
        if (elm_ctype(ictype) == id_radar_vr_obs) then
          vert_loc_ctype(ictype) = VERT_LOCAL_RADAR_VR
        else
          vert_loc_ctype(ictype) = VERT_LOCAL(ityp)
        end if
      end if ! [ ctype_use(ielm_u, ityp) ]
    end do ! [ ielm_u = 1, nid_obs ]
    end do ! [ ityp = 1, nobtype ]

    ! Compute perturbation and departure
    !  -- gross error check
    !  -- QC based on background (radar reflectivity)
    !  -- process Himawari-8 data
    !-----------------------------------------------------------------------------

    allocate(tmpelm(obsda%nobs))

    do n = 1, obsda%nobs
      IF(obsda%qc(n) > 0) CYCLE 

      iof = obsda%set(n)
      iidx = obsda%idx(n)

      tmpelm(n) = obs(iof)%elm(iidx)

      !!! ###### RADAR assimilation ######
      if (obs(iof)%elm(iidx) == id_radar_ref_obs .or. obs(iof)%elm(iidx) == id_radar_ref_zero_obs) then
        if (.not. USE_RADAR_REF) then
          obsda%qc(n) = iqc_otype
          cycle
        end if

        if (obs(iof)%dat(iidx) == undef) then
          obsda%qc(n) = iqc_obs_bad
          cycle
        end if

        !!! obsda%ensval: already converted to dBZ
        mem_ref = 0
        do i = 1, nmem
          if (obsda%ensval(i,n) > RADAR_REF_THRES_DBZ+1.0d-6 ) then
            mem_ref = mem_ref + 1
          end if
        end do

        ! Obs: Rain
        if (obs(iof)%dat(iidx) > RADAR_REF_THRES_DBZ+1.0d-6) then
          if (mem_ref < MIN_RADAR_REF_MEMBER_OBSRAIN) then
            obsda%qc(n) = iqc_ref_mem

            if ( .not. RADAR_PQV ) cycle
            ! When RADAR_PQV=True, pseudo qv obs is assimilated even if mem_ref is
            ! too small
          end if

        else
        ! Obs: No rain
          if (mem_ref < MIN_RADAR_REF_MEMBER_OBSNORAIN) then
            obsda%qc(n) = iqc_ref_mem
            cycle
          end if
        end if

      end if

      if (obs(iof)%elm(iidx) == id_radar_vr_obs) then
        if (.not. USE_RADAR_VR) then
          obsda%qc(n) = iqc_otype
          cycle
        end if
      end if
      !!! ###### end RADAR assimilation ######

      obsda%val(n) = 0.0_RP
      do i = 1, nmem
        obsda%val(n) = obsda%val(n) + obsda%ensval(i,n)
      end do
      obsda%val(n) = obsda%val(n) / real(nmem,kind=RP)

      do i = 1, nmem
        obsda%ensval(i,n) = obsda%ensval(i,n) - obsda%val(n) ! Hdx
      end do
      obsda%val(n) = obs(iof)%dat(iidx) - obsda%val(n) ! y-Hx
      if( ENS_WITH_MDET ) then
         obsda%ensval(mmdetobs,n) = obs(iof)%dat(iidx) - obsda%ensval(mmdetobs,n) ! y-Hx for deterministic run
      end if

      select case (obs(iof)%elm(iidx)) !gross error
      case (id_rain_obs)
        IF(ABS(obsda%val(n)) > GROSS_ERROR_RAIN * obs(iof)%err(iidx)) THEN
          obsda%qc(n) = iqc_gross_err
        END IF
      case (id_radar_ref_obs,id_radar_ref_zero_obs)

        if( RADAR_PQV .and. obsda%val(n) > RADAR_PQV_OMB ) then

          ! pseudo qv
          obsda%val(n) = 0.0_RP
          do i = 1, nmem
            obsda%val(n) = obsda%val(n) + obsda%eqv(i,n)
          enddo
          obsda%val(n) = obsda%val(n) / real(nmem, kind=RP)

          do i = 1, nmem
            obsda%ensval(i,n) = obsda%eqv(i,n) - obsda%val(n) ! Hdx
          enddO

          ! Tetens equation es(Pa)
          qvs = 611.2d0*exp(17.67d0*(obsda%tm(n)-TEM00)/(obsda%tm(n) - TEM00 + 243.5d0))

          ! Saturtion mixing ratio
          qvs = 0.622d0*qvs / ( obsda%pm(n) - qvs )

          obsda%val(n) = qvs - obsda%val(n) ! y-Hx

          if (ENS_WITH_MDET) then
            obsda%ensval(mmdetobs,n) = qvs - obsda%eqv(mmdetobs,n) ! y-Hx for deterministic run
          end if

          obsda%tm(n) = -1.0d0

        else
          IF(ABS(obsda%val(n)) > GROSS_ERROR_RADAR_REF * obs(iof)%err(iidx)) THEN
            obsda%qc(n) = iqc_gross_err
          END IF
        endif
      case (id_radar_vr_obs)
        IF(ABS(obsda%val(n)) > GROSS_ERROR_RADAR_VR * obs(iof)%err(iidx)) THEN
          obsda%qc(n) = iqc_gross_err
        END IF
      case (id_radar_prh_obs)
        IF(ABS(obsda%val(n)) > GROSS_ERROR_RADAR_PRH * obs(iof)%err(iidx)) THEN
          obsda%qc(n) = iqc_gross_err
        END IF
      case default
        IF(ABS(obsda%val(n)) > GROSS_ERROR * obs(iof)%err(iidx)) THEN
          obsda%qc(n) = iqc_gross_err
        END IF
      end select

    END DO

    if( LETKF_DEBUG_LOG ) then
      LOG_INFO("LETKF_debug",'(1x,A,I6,A)') 'OBSERVATIONAL DEPARTURE STATISTICS (IN THIS SUBDOMAIN #', RANK_LCL, '):'

      call monit_dep( obsda%nobs, tmpelm, obsda%val, obsda%qc, monit_nobs, bias, rmse )
      call monit_print( monit_nobs, bias, rmse )
    end if

    deallocate(tmpelm)

    !-------------------------------------------------------------------------------
    ! "Bucket sort" of observations of each combined type (with different sorting meshes)
    !-------------------------------------------------------------------------------

    if (allocated(obsgrd)) deallocate(obsgrd)
    allocate (obsgrd(nctype))

    ! Determine mesh size for bucket sort
    !-----------------------------------------------------------------------------

    do ictype = 1, nctype
      ityp = typ_ctype(ictype)

      if( OBS_SORT_GRID_SPACING(ityp) > 0.0_RP ) then
        target_grdspc = OBS_SORT_GRID_SPACING(ityp)
      else if( MAX_NOBS_PER_GRID(ityp) > 0 ) then
        target_grdspc = 0.1_RP * sqrt( real(MAX_NOBS_PER_GRID(ityp), kind=RP) ) * OBS_MIN_SPACING(ityp) ! need to be tuned
      else
        target_grdspc = hori_loc_ctype(ictype) * dist_zero_fac / 6.0_RP                ! need to be tuned
      end if
      obsgrd(ictype)%ngrd_i = min( ceiling( DX * real( nlon, kind=RP ) / target_grdspc), nlon )
      obsgrd(ictype)%ngrd_j = min( ceiling( DY * real( nlat, kind=RP ) / target_grdspc), nlat )
      obsgrd(ictype)%grdspc_i = DX * real( nlon, kind=RP ) / real( obsgrd( ictype )%ngrd_i, kind=RP )
      obsgrd(ictype)%grdspc_j = DY * real( nlat, kind=RP ) / real( obsgrd( ictype )%ngrd_j, kind=RP )
      if( LETKF_ENTIRE_GRID_SEARCH_X ) then
        obsgrd(ictype)%ngrdsch_i = ceiling( DX * nlong / obsgrd( ictype )%grdspc_i )
      else
        obsgrd(ictype)%ngrdsch_i = ceiling( hori_loc_ctype( ictype ) * dist_zero_fac / obsgrd( ictype )%grdspc_i )
      endif
      if( LETKF_ENTIRE_GRID_SEARCH_Y ) then
        obsgrd(ictype)%ngrdsch_j = ceiling( DY * nlatg / obsgrd( ictype )%grdspc_j )
      else
        obsgrd(ictype)%ngrdsch_j = ceiling( hori_loc_ctype( ictype ) * dist_zero_fac / obsgrd( ictype )%grdspc_j )
      endif
      obsgrd(ictype)%ngrdext_i = obsgrd( ictype )%ngrd_i + obsgrd( ictype )%ngrdsch_i * 2
      obsgrd(ictype)%ngrdext_j = obsgrd( ictype )%ngrd_j + obsgrd( ictype )%ngrdsch_j * 2

      allocate (obsgrd(ictype)%n (  obsgrd(ictype)%ngrd_i, obsgrd(ictype)%ngrd_j, 0:NPRC_LCL-1))
      allocate (obsgrd(ictype)%ac(0:obsgrd(ictype)%ngrd_i, obsgrd(ictype)%ngrd_j, 0:NPRC_LCL-1))
      allocate (obsgrd(ictype)%tot(0:NPRC_LCL-1))
      allocate (obsgrd(ictype)%n_ext (  obsgrd(ictype)%ngrdext_i, obsgrd(ictype)%ngrdext_j))
      allocate (obsgrd(ictype)%ac_ext(0:obsgrd(ictype)%ngrdext_i, obsgrd(ictype)%ngrdext_j))

      obsgrd(ictype)%n (:,:,:) = 0
      obsgrd(ictype)%ac(:,:,:) = 0
      obsgrd(ictype)%tot(:) = 0
      obsgrd(ictype)%n_ext (:,:) = 0
      obsgrd(ictype)%ac_ext(:,:) = 0
      obsgrd(ictype)%tot_ext = 0
      obsgrd(ictype)%tot_sub(:) = 0
      obsgrd(ictype)%tot_g(:) = 0

      allocate (obsgrd(ictype)%next(obsgrd(ictype)%ngrd_i, obsgrd(ictype)%ngrd_j))
    end do

    ! First scan: count the observation numbers in each mesh (in each subdomian)
    !-----------------------------------------------------------------------------

    do n = 1, obsda%nobs
      iof = obsda%set(n)
      iidx = obsda%idx(n)
      ictype = ctype_elmtyp(uid_obs(obs(iof)%elm(iidx)), obs(iof)%typ(iidx))

      if (obsda%qc(n) == iqc_good) then
        call ij_obsgrd( ictype, obs(iof)%ri(iidx), obs(iof)%rj(iidx), i, j )
        if (i < 1) i = 1                                          ! Assume the process assignment was correct,
        if (i > obsgrd(ictype)%ngrd_i) i = obsgrd(ictype)%ngrd_i  ! so this correction is only to remedy the round-off problem.
        if (j < 1) j = 1                                          !
        if (j > obsgrd(ictype)%ngrd_j) j = obsgrd(ictype)%ngrd_j  !

        obsgrd(ictype)%n(i,j,RANK_LCL) = obsgrd(ictype)%n(i,j,RANK_LCL) + 1
        obsgrd(ictype)%tot_sub(i_after_qc) = obsgrd(ictype)%tot_sub(i_after_qc) + 1 ! only used for diagnostic print (obs number after qc)
      end if

      obsgrd(ictype)%tot_sub(i_before_qc) = obsgrd(ictype)%tot_sub(i_before_qc) + 1 ! only used for diagnostic print (obs number before qc)
    end do

    ! Compute the accumulated numbers in each mesh
    !-----------------------------------------------------------------------------

    do ictype = 1, nctype
      if (ictype > 1) then
        obsgrd(ictype)%ac(0,1,RANK_LCL) = obsgrd(ictype-1)%ac(obsgrd(ictype-1)%ngrd_i,obsgrd(ictype-1)%ngrd_j,RANK_LCL)
      end if
      do j = 1, obsgrd(ictype)%ngrd_j
        if (j > 1) then
          obsgrd(ictype)%ac(0,j,RANK_LCL) = obsgrd(ictype)%ac(obsgrd(ictype)%ngrd_i,j-1,RANK_LCL)
        end if
        do i = 1, obsgrd(ictype)%ngrd_i
          obsgrd(ictype)%ac(i,j,RANK_LCL) = obsgrd(ictype)%ac(i-1,j,RANK_LCL) + obsgrd(ictype)%n(i,j,RANK_LCL)
        end do
      end do
      obsgrd(ictype)%next(1:obsgrd(ictype)%ngrd_i,:) = obsgrd(ictype)%ac(0:obsgrd(ictype)%ngrd_i-1,:,RANK_LCL)
    end do

    ! Second scan: save the indices of bucket-sorted observations in obsda%keys(:)
    !-----------------------------------------------------------------------------

    do n = 1, obsda%nobs
      if (obsda%qc(n) == iqc_good) then
        iof = obsda%set(n)
        iidx = obsda%idx(n)
        ictype = ctype_elmtyp(uid_obs(obs(iof)%elm(iidx)), obs(iof)%typ(iidx))

        call ij_obsgrd( ictype, obs(iof)%ri(iidx), obs(iof)%rj(iidx), i, j )
        if (i < 1) i = 1                                         ! Assume the process assignment was correct,
        if (i > obsgrd(ictype)%ngrd_i) i = obsgrd(ictype)%ngrd_i ! so this correction is only to remedy the round-off problem.
        if (j < 1) j = 1                                         !
        if (j > obsgrd(ictype)%ngrd_j) j = obsgrd(ictype)%ngrd_j !

        obsgrd(ictype)%next(i,j) = obsgrd(ictype)%next(i,j) + 1
        obsda%key(obsgrd(ictype)%next(i,j)) = n
      end if
    end do

    ! ALLREDUCE observation number information from subdomains, and compute total numbers
    !-----------------------------------------------------------------------------

    nobs_sub(:) = 0
    nobs_g(:) = 0
    do ictype = 1, nctype
      if (NPRC_LCL > 1) then
        call MPI_ALLREDUCE(MPI_IN_PLACE, obsgrd(ictype)%n, obsgrd(ictype)%ngrd_i*obsgrd(ictype)%ngrd_j*NPRC_LCL, &
                           MPI_INTEGER, MPI_SUM, COMM_LCL, ierr)
        call MPI_ALLREDUCE(MPI_IN_PLACE, obsgrd(ictype)%ac(0:obsgrd(ictype)%ngrd_i,:,:), (obsgrd(ictype)%ngrd_i+1)*obsgrd(ictype)%ngrd_j*NPRC_LCL, &
                           MPI_INTEGER, MPI_SUM, COMM_LCL, ierr)
      end if
      call MPI_ALLREDUCE(obsgrd(ictype)%tot_sub, obsgrd(ictype)%tot_g, n_qc_steps, MPI_INTEGER, MPI_SUM, COMM_LCL, ierr)

      if (ictype == 1) then
        obsgrd(ictype)%tot(:) = obsgrd(ictype)%ac(obsgrd(ictype)%ngrd_i,obsgrd(ictype)%ngrd_j,:)
      else
        obsgrd(ictype)%tot(:) = obsgrd(ictype)%ac(obsgrd(ictype)%ngrd_i,obsgrd(ictype)%ngrd_j,:) &
                              - obsgrd(ictype-1)%ac(obsgrd(ictype-1)%ngrd_i,obsgrd(ictype-1)%ngrd_j,:)
      end if

      nobs_sub(:) = nobs_sub(:) + obsgrd(ictype)%tot_sub(:)
      nobs_g(:) = nobs_g(:) + obsgrd(ictype)%tot_g(:)

      deallocate (obsgrd(ictype)%next)
    end do

    nobstotalg = nobs_g(i_after_qc) ! total obs number in the global domain (all types)

    ! Calculate observation numbers in the extended (localization) subdomain,
    ! in preparation for communicating obsetvations in the extended subdomain
    !-----------------------------------------------------------------------------

    call rank_1d_2d(RANK_LCL, myp_i, myp_j)

    do ictype = 1, nctype
      imin1 = myp_i*obsgrd(ictype)%ngrd_i+1 - obsgrd(ictype)%ngrdsch_i
      imax1 = (myp_i+1)*obsgrd(ictype)%ngrd_i + obsgrd(ictype)%ngrdsch_i
      jmin1 = myp_j*obsgrd(ictype)%ngrd_j+1 - obsgrd(ictype)%ngrdsch_j
      jmax1 = (myp_j+1)*obsgrd(ictype)%ngrd_j + obsgrd(ictype)%ngrdsch_j

      do ip = 0, NPRC_LCL-1
        call rank_1d_2d(ip, ip_i, ip_j)
        imin2 = max(1, imin1 - ip_i*obsgrd(ictype)%ngrd_i)
        imax2 = min(obsgrd(ictype)%ngrd_i, imax1 - ip_i*obsgrd(ictype)%ngrd_i)
        jmin2 = max(1, jmin1 - ip_j*obsgrd(ictype)%ngrd_j)
        jmax2 = min(obsgrd(ictype)%ngrd_j, jmax1 - ip_j*obsgrd(ictype)%ngrd_j)
        if (imin2 > imax2 .or. jmin2 > jmax2) cycle

        ishift = (ip_i - myp_i) * obsgrd(ictype)%ngrd_i + obsgrd(ictype)%ngrdsch_i
        jshift = (ip_j - myp_j) * obsgrd(ictype)%ngrd_j + obsgrd(ictype)%ngrdsch_j
        obsgrd(ictype)%n_ext(imin2+ishift:imax2+ishift, jmin2+jshift:jmax2+jshift) = obsgrd(ictype)%n(imin2:imax2, jmin2:jmax2, ip)
      end do

      if (ictype > 1) then
        obsgrd(ictype)%ac_ext(0,1) = obsgrd(ictype-1)%ac_ext(obsgrd(ictype-1)%ngrdext_i,obsgrd(ictype-1)%ngrdext_j)
      end if
      do j = 1, obsgrd(ictype)%ngrdext_j
        if (j > 1) then
          obsgrd(ictype)%ac_ext(0,j) = obsgrd(ictype)%ac_ext(obsgrd(ictype)%ngrdext_i,j-1)
        end if
        do i = 1, obsgrd(ictype)%ngrdext_i
          obsgrd(ictype)%ac_ext(i,j) = obsgrd(ictype)%ac_ext(i-1,j) + obsgrd(ictype)%n_ext(i,j)
        end do
      end do

      if (ictype == 1) then
        obsgrd(ictype)%tot_ext = obsgrd(ictype)%ac_ext(obsgrd(ictype)%ngrdext_i,obsgrd(ictype)%ngrdext_j)
      else
        obsgrd(ictype)%tot_ext = obsgrd(ictype)%ac_ext(obsgrd(ictype)%ngrdext_i,obsgrd(ictype)%ngrdext_j) &
                               - obsgrd(ictype-1)%ac_ext(obsgrd(ictype-1)%ngrdext_i,obsgrd(ictype-1)%ngrdext_j)
      end if
    end do

    if (nctype > 0) then
      nobstotal = obsgrd(nctype)%ac_ext(obsgrd(nctype)%ngrdext_i,obsgrd(nctype)%ngrdext_j) ! total obs number in the extended subdomain (all types)

      maxnobs_per_ctype = obsgrd(1)%tot_ext
      do ictype = 2, nctype
        maxnobs_per_ctype = max(maxnobs_per_ctype, obsgrd(ictype)%tot_ext)
      end do
    else
      nobstotal = 0
      maxnobs_per_ctype = 0
    end if

    ! Construct sorted obsda_sort: 
    !-----------------------------------------------------------------------------

    ! 1) Copy the observation data in own subdomains to send buffer (obsbufs) with
    ! sorted order
    !-----------------------------------------------------------------------------
    if (nctype > 0) then
      cntr(:) = obsgrd(nctype)%ac(obsgrd(nctype)%ngrd_i,obsgrd(nctype)%ngrd_j,:)
      cnts = cntr(RANK_LCL+1)
      dspr(1) = 0
      do ip = 2, NPRC_LCL 
        dspr(ip) = dspr(ip-1) + cntr(ip-1)
      end do

      nensobs_mod = mod(nensobs, NPRC_ENS)
      nensobs_div = (nensobs - nensobs_mod) / NPRC_ENS
      if (RANK_ENS < nensobs_mod) then
        im_obs_1 = (nensobs_div+1) * RANK_ENS + 1
        im_obs_2 = (nensobs_div+1) * (RANK_ENS+1)
        nensobs_part = nensobs_div + 1
      else
        im_obs_1 = (nensobs_div+1) * nensobs_mod + nensobs_div * (RANK_ENS-nensobs_mod) + 1
        im_obs_2 = (nensobs_div+1) * nensobs_mod + nensobs_div * (RANK_ENS-nensobs_mod+1)
        nensobs_part = nensobs_div
      end if

      obsbufs%nobs = nobs_sub(i_after_qc)
      call obs_da_value_allocate(obsbufs, nensobs_part)

      do n = 1, nobs_sub(i_after_qc)
        obsbufs%set(n) = obsda%set(obsda%key(n))
        obsbufs%idx(n) = obsda%idx(obsda%key(n))
        obsbufs%val(n) = obsda%val(obsda%key(n))
        obsbufs%tm(n) = obsda%tm(obsda%key(n))
        if (nensobs_part > 0) then
          obsbufs%ensval(1:nensobs_part,n) = obsda%ensval(im_obs_1:im_obs_2,obsda%key(n))
        end if
        obsbufs%qc(n) = obsda%qc(obsda%key(n))
      end do
    end if

    call obs_da_value_deallocate(obsda)

    ! 2) Communicate to get global observations;
    !    for variables with an ensemble dimension (ensval),
    !    only obtain data from partial members (nensobs_part) to save memory usage
    !-----------------------------------------------------------------------------
    if (nctype > 0) then
      obsbufr%nobs = nobs_g(i_after_qc)
      call obs_da_value_allocate(obsbufr, nensobs_part)

      call MPI_ALLGATHERV( obsbufs%set, cnts, MPI_INTEGER, obsbufr%set, cntr, dspr, MPI_INTEGER, COMM_LCL, ierr )
      call MPI_ALLGATHERV( obsbufs%idx, cnts, MPI_INTEGER, obsbufr%idx, cntr, dspr, MPI_INTEGER, COMM_LCL, ierr )
      call MPI_ALLGATHERV( obsbufs%val, cnts, datatype,    obsbufr%val, cntr, dspr, datatype,    COMM_LCL, ierr )
      call MPI_ALLGATHERV( obsbufs%tm,  cnts, datatype,    obsbufr%tm,  cntr, dspr, datatype,    COMM_LCL, ierr )
      if (nensobs_part > 0) then
        call MPI_ALLGATHERV(obsbufs%ensval, cnts*nensobs_part, datatype, obsbufr%ensval, cntr*nensobs_part, dspr*nensobs_part, datatype, COMM_LCL, ierr)
      end if
      call MPI_ALLGATHERV(obsbufs%qc, cnts, MPI_INTEGER, obsbufr%qc, cntr, dspr, MPI_INTEGER, COMM_LCL, ierr)

      call obs_da_value_deallocate(obsbufs)
    end if

    ! 3) Copy observation data within the extended (localization) subdomains
    !    from receive buffer (obsbufr) to obsda_sort; rearrange with sorted order
    !-----------------------------------------------------------------------------
    obsda_sort%nobs = nobstotal
    call obs_da_value_allocate(obsda_sort, nensobs)

    do ip = 0, NPRC_LCL-1
      call rank_1d_2d(ip, ip_i, ip_j)

      do ictype = 1, nctype
        imin1 = myp_i*obsgrd(ictype)%ngrd_i+1 - obsgrd(ictype)%ngrdsch_i
        imax1 = (myp_i+1)*obsgrd(ictype)%ngrd_i + obsgrd(ictype)%ngrdsch_i
        jmin1 = myp_j*obsgrd(ictype)%ngrd_j+1 - obsgrd(ictype)%ngrdsch_j
        jmax1 = (myp_j+1)*obsgrd(ictype)%ngrd_j + obsgrd(ictype)%ngrdsch_j

        imin2 = max(1, imin1 - ip_i*obsgrd(ictype)%ngrd_i)
        imax2 = min(obsgrd(ictype)%ngrd_i, imax1 - ip_i*obsgrd(ictype)%ngrd_i)
        jmin2 = max(1, jmin1 - ip_j*obsgrd(ictype)%ngrd_j)
        jmax2 = min(obsgrd(ictype)%ngrd_j, jmax1 - ip_j*obsgrd(ictype)%ngrd_j)
        if (imin2 > imax2 .or. jmin2 > jmax2) cycle

        ishift = (ip_i - myp_i) * obsgrd(ictype)%ngrd_i + obsgrd(ictype)%ngrdsch_i
        jshift = (ip_j - myp_j) * obsgrd(ictype)%ngrd_j + obsgrd(ictype)%ngrdsch_j

        do j = jmin2, jmax2
          ns_ext = obsgrd(ictype)%ac_ext(imin2+ishift-1,j+jshift) + 1
          ne_ext = obsgrd(ictype)%ac_ext(imax2+ishift  ,j+jshift)
          if (ns_ext > ne_ext) cycle

          ns_bufr = dspr(ip+1) + obsgrd(ictype)%ac(imin2-1,j,ip) + 1
          ne_bufr = dspr(ip+1) + obsgrd(ictype)%ac(imax2  ,j,ip)

          obsda_sort%set(ns_ext:ne_ext) = obsbufr%set(ns_bufr:ne_bufr)
          obsda_sort%idx(ns_ext:ne_ext) = obsbufr%idx(ns_bufr:ne_bufr)
          obsda_sort%val(ns_ext:ne_ext) = obsbufr%val(ns_bufr:ne_bufr)
          obsda_sort%tm(ns_ext:ne_ext) = obsbufr%tm(ns_bufr:ne_bufr)
          if (nensobs_part > 0) then
            obsda_sort%ensval(im_obs_1:im_obs_2,ns_ext:ne_ext) = obsbufr%ensval(1:nensobs_part,ns_bufr:ne_bufr)
          end if
          obsda_sort%qc(ns_ext:ne_ext) = obsbufr%qc(ns_bufr:ne_bufr)
        end do
      end do
    end do

    ! Save the keys of observations within the subdomain (excluding the localization buffer area)
    obsda_sort%nobs_in_key = 0
    do ictype = 1, nctype
      imin1 = obsgrd(ictype)%ngrdsch_i + 1
      imax1 = obsgrd(ictype)%ngrdsch_i + obsgrd(ictype)%ngrd_i
      jmin1 = obsgrd(ictype)%ngrdsch_j + 1
      jmax1 = obsgrd(ictype)%ngrdsch_j + obsgrd(ictype)%ngrd_j
      call obs_choose_ext(ictype, imin1, imax1, jmin1, jmax1, obsda_sort%nobs_in_key, obsda_sort%key)

      deallocate (obsgrd(ictype)%n)
      deallocate (obsgrd(ictype)%ac)
      deallocate (obsgrd(ictype)%n_ext)
    end do

    if (nctype > 0) then
      call obs_da_value_deallocate(obsbufr)
    end if

    ! 4) For variables with an ensemble dimension (ensval),
    !    ALLREDUCE among the ensemble dimension to obtain data of all members
    !-----------------------------------------------------------------------------
    if (NPRC_ENS > 1) then
      call MPI_ALLREDUCE(MPI_IN_PLACE, obsda_sort%ensval, nensobs*nobstotal, datatype, MPI_SUM, COMM_ENS, ierr)
    end if

    if( LETKF_DEBUG_LOG ) then
      LOG_INFO("LETKF_debug",'(1x,A,I6,A)') 'OBSERVATION COUNTS (GLOABL AND IN THIS SUBDOMAIN #', RANK_LCL, '):'
      LOG_INFO("LETKF_debug",'(1x,A)')      '====================================================================='
      LOG_INFO("LETKF_debug",'(1x,A)')      'TYPE   VAR      GLOBAL     GLOBAL  SUBDOMAIN  SUBDOMAIN EXT_SUBDOMAIN'
      LOG_INFO("LETKF_debug",'(1x,A)')      '             before QC   after QC  before QC   after QC      after QC'
      LOG_INFO("LETKF_debug",'(1x,A)')      '---------------------------------------------------------------------'
      do ictype = 1, nctype
        ityp = typ_ctype(ictype)
        ielm_u = elm_u_ctype(ictype)
        LOG_INFO("LETKF_debug",'(1x,A6,1x,A3,1x,4I11,I14)') obtypelist(ityp), obelmlist(ielm_u), &
                                                            obsgrd(ictype)%tot_g(i_before_qc),   &
                                                            obsgrd(ictype)%tot_g(i_after_qc),    &
                                                            obsgrd(ictype)%tot_sub(i_before_qc), &
                                                            obsgrd(ictype)%tot_sub(i_after_qc),  &
                                                            obsgrd(ictype)%tot_ext
      end do
      LOG_INFO("LETKF_debug",'(1x,A)')             '---------------------------------------------------------------------'
      LOG_INFO("LETKF_debug",'(1x,A,5x,4I11,I14)') 'TOTAL ', nobs_g(i_before_qc), nobs_g(i_after_qc), nobs_sub(i_before_qc), nobs_sub(i_after_qc)
      LOG_INFO("LETKF_debug",'(1x,A)')             '====================================================================='
    end if

    return
  end subroutine LETKF_obs_initialize

  !-----------------------------------------------------------------------
  ! Data Assimilation
  !-----------------------------------------------------------------------
  subroutine LETKF_system( &
      OBS_IN_NUM,    & ! [IN]
      OBS_IN_FORMAT, & ! [IN]
      U,             & ! [INOUT]
      V,             & ! [INOUT]
      W,             & ! [INOUT]
      TEMP,          & ! [INOUT]
      PRES,          & ! [INOUT]
      QV,            & ! [INOUT]
      QC,            & ! [INOUT]
      QR,            & ! [INOUT]
      QI,            & ! [INOUT]
      QS,            & ! [INOUT]
      QG             ) ! [INOUT]
    use scale_prc, only: &
      PRC_abort
    use scale_const, only: &
      UNDEF => CONST_UNDEF
    use scale_statistics, only: &
      AVERAGE => STATISTICS_AVERAGE
    use scale_matrix, only: &
      MATRIX_SOLVER_eigenvalue_decomposition
    use scale_random, only: &
      Knuth_Shuffle => RANDOM_Knuth_shuffle
    implicit none

    integer, intent(in) :: OBS_IN_NUM
    character(len=H_LONG), intent(in) :: OBS_IN_FORMAT(:)

    real(RP), intent(inout) :: U   (nlev,nlon,nlat)
    real(RP), intent(inout) :: V   (nlev,nlon,nlat)
    real(RP), intent(inout) :: W   (nlev,nlon,nlat)
    real(RP), intent(inout) :: TEMP(nlev,nlon,nlat)
    real(RP), intent(inout) :: PRES(nlev,nlon,nlat)
    real(RP), intent(inout) :: QV  (nlev,nlon,nlat)
    real(RP), intent(inout) :: QC  (nlev,nlon,nlat)
    real(RP), intent(inout) :: QR  (nlev,nlon,nlat)
    real(RP), intent(inout) :: QI  (nlev,nlon,nlat)
    real(RP), intent(inout) :: QS  (nlev,nlon,nlat)
    real(RP), intent(inout) :: QG  (nlev,nlon,nlat)

    real(RP) :: gues3d(nij1,nlev,nens+1,nv3d) ! background ensemble
    real(RP) :: gues2d(nij1,     nens+1,nv2d) ! background ensemble

    real(RP) :: anal3d(nij1,nlev,nens+1,nv3d)   ! analysis ensemble
    real(RP) :: anal2d(nij1,     nens+1,nv2d)   ! analysis ensemble

    real(RP) :: addi3d(nij1,nlev,nens+1,nv3d)
    real(RP) :: addi2d(nij1,     nens+1,nv2d)

    real(RP) :: v3dg(nlev,nlon,nlat,nv3d)
    real(RP) :: v2dg(     nlon,nlat,nv2d)

    real(RP) :: work3d(nij1,nlev,nv3d)
    real(RP) :: work2d(nij1,     nv2d)
    real(RP), allocatable :: work3da(:,:,:)     !GYL
    real(RP), allocatable :: work2da(:,:)       !GYL
    real(RP), allocatable :: work3dn(:,:,:,:)   !GYL
    real(RP), allocatable :: work2dn(:,:,:)     !GYL
    real(RP), allocatable :: work3dg(:,:,:,:)
    real(RP), allocatable :: work2dg(:,:,:)

    real(RP), allocatable :: hdxf(:,:)
    real(RP), allocatable :: rdiag(:)
    real(RP), allocatable :: rloc(:)
    real(RP), allocatable :: dep(:)
    real(RP), allocatable :: depd(:)            !GYL

    integer :: ctype_merge(nid_obs,nobtype)

    integer :: var_local_n2nc_max
    integer :: var_local_n2nc(nv3d+nv2d)
    integer :: var_local_n2n(nv3d+nv2d)
    logical :: var_local_n2n_found
    integer :: n2n, n2nc

    real(RP) :: parm
    real(RP), allocatable :: trans(:,:,:)
    real(RP), allocatable :: transm(:,:)
    real(RP), allocatable :: transmd(:,:)
    real(RP), allocatable :: pa(:,:,:)
    real(RP) :: transrlx(nmem,nmem)
    logical :: trans_done(nv3d+nv2d)

    integer :: ij,ilev,n,m,i,j,k,nobsl
    integer :: nobsl_t(nid_obs,nobtype)            !GYL
    real(RP) :: cutd_t(nid_obs,nobtype)        !GYL
    real(RP) :: beta                           !GYL
    real(RP) :: tmpinfl                        !GYL
    real(RP) :: q_mean,q_sprd                  !GYL
    real(RP) :: q_anal(nmem)                 !GYL

    integer :: mshuf,ierr                          !GYL
    integer :: ishuf(nmem)                       !GYL
    real(RP), allocatable :: addinfl_weight(:) !GYL
    real(RP) :: rdx,rdy,rdxy,ref_min_dist      !GYL
    integer :: ic,ic2,iob                          !GYL

    integer, allocatable :: search_q0(:,:,:)

    LOG_PROGRESS(*) 'data-assimilation / LETKF / system'

    ! -- prepare the first-guess data
    do j = 1, nlat
    do i = 1, nlon
    do k = 1, nlev
       v3dg(k,i,j,iv3d_u ) = U   (k,i,j)
       v3dg(k,i,j,iv3d_v ) = V   (k,i,j)
       v3dg(k,i,j,iv3d_w ) = W   (k,i,j)
       v3dg(k,i,j,iv3d_t ) = TEMP(k,i,j)
       v3dg(k,i,j,iv3d_p ) = PRES(k,i,j)
       v3dg(k,i,j,iv3d_q ) = QV  (k,i,j)
       v3dg(k,i,j,iv3d_qc) = QC  (k,i,j)
       v3dg(k,i,j,iv3d_qr) = QR  (k,i,j)
       v3dg(k,i,j,iv3d_qi) = QI  (k,i,j)
       v3dg(k,i,j,iv3d_qs) = QS  (k,i,j)
       v3dg(k,i,j,iv3d_qg) = QG  (k,i,j)
    end do
    end do
    end do
    do j = 1, nlat
    do i = 1, nlon
       v2dg(i,j,:) = UNDEF ! tentative
    end do
    end do

    if( LETKF_DEBUG_LOG ) then
       call monit_obs_mpi( OBS_IN_NUM, OBS_IN_FORMAT, v3dg, v2dg, monit_step=1 )
    endif

    call scatter_grd_mpi_all2all( 1, nens, v3dg, v2dg, gues3d(:,:,1:nens,:), gues2d(:,1:nens,:) )

    ! -- obtain the ensemble mean
    do n  = 1, nv3d
    do k  = 1, nlev
    do ij = 1, nij1
      gues3d(ij,k,mmean,n) = AVERAGE( gues3d(ij,k,1:nmem,n), UNDEF )
    end do
    end do
    end do
    do n  = 1, nv2d
    do ij = 1, nij1
      gues2d(ij,mmean,n) = AVERAGE( gues2d(ij,1:nmem,n), UNDEF )
    end do
    end do

    !
    ! Variable localization
    !
    var_local(:,1) = VAR_LOCAL_UV(:)
    var_local(:,2) = VAR_LOCAL_T(:)
    var_local(:,3) = VAR_LOCAL_Q(:)
    var_local(:,4) = VAR_LOCAL_PS(:)
    var_local(:,5) = VAR_LOCAL_RAIN(:)
    var_local(:,6) = VAR_LOCAL_TC(:)
    var_local(:,7) = VAR_LOCAL_RADAR_REF(:)
    var_local(:,8) = VAR_LOCAL_RADAR_VR(:)

    var_local_n2nc_max = 1
    var_local_n2nc(1) = 1
    var_local_n2n(1) = 1

    do n = 2, nv3d+nv2d
      var_local_n2n_found = .false.
      do i = 1, var_local_n2nc_max
        !if (maxval(abs(var_local(var_local_n2nc(i),:) - var_local(n,:))) < tiny(var_local(1,1))) then
        if (all(var_local(var_local_n2nc(i),:) == var_local(n,:))) then
          var_local_n2nc(n)   = var_local_n2nc(i)
          var_local_n2n(n)    = var_local_n2n(var_local_n2nc(n))
          var_local_n2n_found = .true.
          exit
        end if
      end do
      if( .NOT. var_local_n2n_found ) then
        var_local_n2nc_max = var_local_n2nc_max + 1
        var_local_n2nc(n)  = var_local_n2nc_max
        var_local_n2n(n)   = n
      end if
    end do

    !
    ! Observation number limit (*to be moved to namelist*)
    !
    ctype_merge(:,:) = 0
    ctype_merge(uid_obs(id_radar_ref_obs),22) = 1
    ctype_merge(uid_obs(id_radar_ref_zero_obs),22) = 1

    allocate (n_merge(nctype))
    allocate (ic_merge(nid_obs*nobtype,nctype))
    n_merge(:) = 1
    do ic = 1, nctype
      if (n_merge(ic) > 0) then
        ic_merge(1,ic) = ic
        if (ctype_merge(elm_u_ctype(ic),typ_ctype(ic)) > 0) then
          do ic2 = ic+1, nctype
            if (ctype_merge(elm_u_ctype(ic2),typ_ctype(ic2)) == ctype_merge(elm_u_ctype(ic),typ_ctype(ic))) then
              n_merge(ic) = n_merge(ic) + 1
              ic_merge(n_merge(ic),ic) = ic2
              n_merge(ic2) = 0
            end if
          end do
        end if
      end if
    end do
    n_merge_max = maxval(n_merge)

    allocate (search_q0(nctype,nv3d+1,nij1))
    search_q0(:,:,:) = 1

    radar_only = .true.
    do ic = 1, nctype
      if (obtypelist(typ_ctype(ic)) /= 'PHARAD') then
        radar_only = .false.
        exit
      end if
    end do
    !
    ! FCST PERTURBATIONS
    !
    !  .... this has been done by write_ensmean in letkf.f90
    !  CALL ensmean_grd(nens,nij1,gues3d,gues2d,mean3d,mean2d)
    DO n=1,nv3d
    DO m=1,nmem
    DO k=1,nlev
    DO i=1,nij1
       gues3d(i,k,m,n) = gues3d(i,k,m,n) - gues3d(i,k,mmean,n)
    END DO
    END DO
    END DO
    END DO
    DO n=1,nv2d
    DO m=1,nmem
    DO i=1,nij1
       gues2d(i,m,n) = gues2d(i,m,n) - gues2d(i,mmean,n)
    END DO
    END DO
    END DO

    !
    ! multiplicative inflation
    !
    IF(INFL_MUL > 0.0d0) THEN  ! fixed multiplicative inflation parameter
      work3d = INFL_MUL
      work2d = INFL_MUL
    ELSE  ! 3D parameter values are read-in
      LOG_ERROR("LETKF_system",*) 'This system has not been implemented yet. INFL_MUL must be greather than 0.0.'
      call PRC_abort
  ! not used now (TODO)
  !    allocate (work3dg(nlon,nlat,nlev,nv3d))
  !    allocate (work2dg(nlon,nlat,nv2d))
  !    IF(myrank_e == mmean_rank_e) THEN
  !      call read_restart(INFL_MUL_IN_BASENAME,work3dg,work2dg)
  !    END IF
  !
  !    CALL scatter_grd_mpi(mmean_rank_e,work3dg,work2dg,work3d,work2d)
  !
    END IF
    IF(INFL_MUL_MIN > 0.0d0) THEN
      work3d = max(work3d, INFL_MUL_MIN)
      work2d = max(work2d, INFL_MUL_MIN)
    END IF

    ! This loop cannot use OpenMP on FUGAKU (T. Honda, as of 10/16/2020)
    allocate( hdxf( nobstotal, nmem ) )
    allocate (rdiag(nobstotal))
    allocate (rloc (nobstotal))
    allocate (dep  (nobstotal))
    if( ENS_WITH_MDET ) then
       allocate (depd (nobstotal))
    end if
    allocate (trans  (nmem,nmem,var_local_n2nc_max))
    allocate (transm (nmem,     var_local_n2nc_max))
    allocate (transmd(nmem,     var_local_n2nc_max))
    allocate (pa     (nmem,nmem,var_local_n2nc_max))

    !
    ! MAIN ASSIMILATION LOOP
    !
    DO ilev=1,nlev
      DO ij=1,nij1

        trans_done(:) = .false.                                                          !GYL

        ! weight parameter based on grid locations (not for covariance inflation purpose)
        ! if the weight is zero, no analysis update is needed
        call relax_beta(rig1(ij),rjg1(ij),hgt1(ij,ilev),beta)

        if (beta == 0.0d0) then
          do n = 1, nv3d
            do m = 1, nmem
              anal3d(ij,ilev,m,n) = gues3d(ij,ilev,mmean,n) + gues3d(ij,ilev,m,n)
            end do
            if (ENS_WITH_MDET) then
              anal3d(ij,ilev,mmdet,n) = gues3d(ij,ilev,mmdet,n)
            end if
          end do
          if (ilev == 1) then
            do n = 1, nv2d
              do m = 1, nmem
                anal2d(ij,m,n) = gues2d(ij,mmean,n) + gues2d(ij,m,n)
              end do
              if (ENS_WITH_MDET) then
                anal2d(ij,mmdet,n) = gues2d(ij,mmdet,n)
              end if
            end do
          end if

          cycle
        end if

        ! update 3D variables
        DO n=1,nv3d

          n2nc = var_local_n2nc(n)
          n2n = var_local_n2n(n)

          if (gues3d(ij,ilev,mmean,iv3d_p) < Q_UPDATE_TOP .and. n >= iv3d_q .and. n <= iv3d_qg) then !GYL - Upper bound of Q update levels
            do m = 1, nmem                                                             !GYL
              anal3d(ij,ilev,m,n) = gues3d(ij,ilev,mmean,n) + gues3d(ij,ilev,m,n)        !GYL
            end do                                                                       !GYL
            if (ENS_WITH_MDET) then                                                      !GYL
              anal3d(ij,ilev,mmdet,n) = gues3d(ij,ilev,mmdet,n)                          !GYL
            end if                                                                       !GYL

            cycle                                                                        !GYL
          end if                                                                         !GYL

          if (RELAX_TO_INFLATED_PRIOR) then
            parm = work3d(ij,ilev,n)
          else
            parm = 1.0d0
          end if

          ! calculate mean and perturbation weights
          if (trans_done(n2nc)) then
            ! if weights already computed for other variables can be re-used(no variable localization), do not need to compute again
            if (INFL_MUL_ADAPTIVE) then
              work3d(ij,ilev,n) = work3d(ij,ilev,n2n)
            end if
            if (NOBS_OUT) then
              work3dn(:,ij,ilev,n) = work3dn(:,ij,ilev,n2n)
            end if

          ELSE
            ! compute weights with localized observations
            if (ENS_WITH_MDET) then                                                      !GYL
              CALL obs_local(obs(:),rig1(ij),rjg1(ij),gues3d(ij,ilev,mmean,iv3d_p),hgt1(ij,ilev),n,& !GYL
                             hdxf,rdiag,rloc,dep,nobsl,depd=depd,nobsl_t=nobsl_t,cutd_t=cutd_t,srch_q0=search_q0(:,n,ij)) !GYL
            else                                                                         !GYL
              CALL obs_local(obs(:),rig1(ij),rjg1(ij),gues3d(ij,ilev,mmean,iv3d_p),hgt1(ij,ilev),n,& !GYL
                             hdxf,rdiag,rloc,dep,nobsl,nobsl_t=nobsl_t,cutd_t=cutd_t,srch_q0=search_q0(:,n,ij)) !GYL
            end if                                                                       !GYL
            IF(RELAX_ALPHA_SPREAD /= 0.0d0) THEN                                         !GYL
              if (ENS_WITH_MDET) then                                                    !GYL
                CALL letkf_core(nmem,nobstotal,nobsl,hdxf,rdiag,rloc,dep,work3d(ij,ilev,n), & !GYL
                                trans(:,:,n2nc),transm=transm(:,n2nc),pao=pa(:,:,n2nc), & !GYL
                                rdiag_wloc=.true.,infl_update=INFL_MUL_ADAPTIVE, &       !GYL
                                depd=depd,transmd=transmd(:,n2nc))                       !GYL
              else                                                                       !GYL
                CALL letkf_core(nmem,nobstotal,nobsl,hdxf,rdiag,rloc,dep,work3d(ij,ilev,n), & !GYL
                                trans(:,:,n2nc),transm=transm(:,n2nc),pao=pa(:,:,n2nc), & !GYL
                                rdiag_wloc=.true.,infl_update=INFL_MUL_ADAPTIVE)         !GYL
              end if                                                                     !GYL
            ELSE                                                                         !GYL
              if (ENS_WITH_MDET) then                                                    !GYL
                CALL letkf_core(nmem,nobstotal,nobsl,hdxf,rdiag,rloc,dep,work3d(ij,ilev,n), & !GYL
                                trans(:,:,n2nc),transm=transm(:,n2nc),           &       !GYL
                                rdiag_wloc=.true.,infl_update=INFL_MUL_ADAPTIVE, &       !GYL
                                depd=depd,transmd=transmd(:,n2nc))                       !GYL
              else                                                                       !GYL
                CALL letkf_core(nmem,nobstotal,nobsl,hdxf,rdiag,rloc,dep,work3d(ij,ilev,n), & !GYL
                                trans(:,:,n2nc),transm=transm(:,n2nc),           &       !GYL
                                rdiag_wloc=.true.,infl_update=INFL_MUL_ADAPTIVE)         !GYL
              end if                                                                     !GYL
            END IF                                                                       !GYL
            trans_done(n2nc) = .true.                                                    !GYL
            IF(NOBS_OUT) THEN                                                            !GYL
              work3dn(:,ij,ilev,n) = real(sum(nobsl_t, dim=1),RP)                    !GYL !!! NOBS: sum over all variables for each report type
              work3dn(nobtype+1,ij,ilev,n) = real(nobsl_t(9,22),RP)                  !GYL !!! NOBS: ref
              work3dn(nobtype+2,ij,ilev,n) = real(nobsl_t(10,22),RP)                 !GYL !!! NOBS: re0
              work3dn(nobtype+3,ij,ilev,n) = real(nobsl_t(11,22),RP)                 !GYL !!! NOBS: vr
              work3dn(nobtype+4,ij,ilev,n) = real(cutd_t(9,22),RP)                   !GYL !!! CUTOFF_DIST: ref
              work3dn(nobtype+5,ij,ilev,n) = real(cutd_t(10,22),RP)                  !GYL !!! CUTOFF_DIST: re0
              work3dn(nobtype+6,ij,ilev,n) = real(cutd_t(11,22),RP)                  !GYL !!! CUTOFF_DIST: vr
            END IF                                                                       !GYL

          END IF

          ! relaxation via LETKF weight
          IF(RELAX_ALPHA /= 0.0d0) THEN                                                  !GYL - RTPP method (Zhang et al. 2004)
            CALL weight_RTPP(trans(:,:,n2nc),parm,transrlx)                              !GYL
          ELSE IF(RELAX_ALPHA_SPREAD /= 0.0d0) THEN                                      !GYL - RTPS method (Whitaker and Hamill 2012)
            IF(RELAX_SPREAD_OUT) THEN                                                    !GYL
              CALL weight_RTPS(trans(:,:,n2nc),pa(:,:,n2nc),gues3d(ij,ilev,:,n), &       !GYL
                               parm,transrlx,work3da(ij,ilev,n))                         !GYL
            ELSE                                                                         !GYL
              CALL weight_RTPS(trans(:,:,n2nc),pa(:,:,n2nc),gues3d(ij,ilev,:,n), &       !GYL
                               parm,transrlx,tmpinfl)                                    !GYL
            END IF                                                                       !GYL
          ELSE                                                                           !GYL
            transrlx = trans(:,:,n2nc)                                                   !GYL - No relaxation
          END IF                                                                         !GYL

          ! total weight matrix
          do m = 1, nmem
            do k = 1, nmem
              transrlx(k,m) = (transrlx(k,m) + transm(k,n2nc)) * beta
            end do
            transrlx(m,m) = transrlx(m,m) + (1.0_RP-beta)
          end do
          ! analysis update of members
          do m = 1, nmem
            anal3d(ij,ilev,m,n) = gues3d(ij,ilev,mmean,n)
            do k = 1, nmem
              anal3d(ij,ilev,m,n) = anal3d(ij,ilev,m,n) + gues3d(ij,ilev,k,n) * transrlx(k,m)
            end do
          end do

          if( ENS_WITH_MDET ) then
            ! analysis update of deterministic run
            anal3d(ij,ilev,mmdet,n) = 0.0_RP
            do k = 1, nmem
              anal3d(ij,ilev,mmdet,n) = anal3d(ij,ilev,mmdet,n) + gues3d(ij,ilev,k,n) * transmd(k,n2nc)
            end do
            anal3d(ij,ilev,mmdet,n) = gues3d(ij,ilev,mmdet,n) + anal3d(ij,ilev,mmdet,n) * beta
          end if

          ! limit q spread
          IF(Q_SPRD_MAX > 0.0d0 .and. n == iv3d_q) THEN                                  !GYL
            q_mean = SUM(anal3d(ij,ilev,1:nmem,n)) / REAL(nmem,RP)               !GYL
            q_sprd = 0.0d0                                                               !GYL
            DO m=1,nmem                                                                !GYL
              q_anal(m) = anal3d(ij,ilev,m,n) - q_mean                                   !GYL
              q_sprd = q_sprd + q_anal(m)**2                                             !GYL
            END DO                                                                       !GYL

            if ( q_mean > 0.0_RP ) then
              q_sprd = SQRT(q_sprd / REAL(nmem-1,RP)) / q_mean                       !GYL
              IF(q_sprd > Q_SPRD_MAX) THEN                                                 !GYL
                DO m=1,nmem                                                              !GYL
                  anal3d(ij,ilev,m,n) = q_mean + q_anal(m) * Q_SPRD_MAX / q_sprd           !GYL
                END DO                                                                     !GYL
              END IF                                                                       !GYL
            endif
          END IF                                                                         !GYL

        END DO ! [ n=1,nv3d ]

        ! update 2D variables at ilev = 1
        IF(ilev == 1) THEN 

          DO n=1,nv2d

            n2nc = var_local_n2nc(nv3d+n)
            n2n = var_local_n2n(nv3d+n)

            if (RELAX_TO_INFLATED_PRIOR) then
              parm = work2d(ij,n)
            else
              parm = 1.0d0
            end if

            ! calculate mean and perturbation weights
            if (trans_done(n2nc)) then
              ! if weights already computed for other variables can be re-used(no variable localization), do not need to compute again
              IF(n2n <= nv3d) then
                if (INFL_MUL_ADAPTIVE) then
                  work2d(ij,n) = work3d(ij,ilev,n2n)
                end if
                if (NOBS_OUT) then
                  work2dn(:,ij,n) = work3dn(:,ij,ilev,n2n)
                end if
              else
                if (INFL_MUL_ADAPTIVE) then
                  work2d(ij,n) = work2d(ij,n2n-nv3d)
                end if
                if (NOBS_OUT) then
                  work2dn(:,ij,n) = work2dn(:,ij,n2n-nv3d)
                end if
              end if

            ELSE
              ! compute weights with localized observations
              if (ENS_WITH_MDET) then                                                    !GYL
                CALL obs_local(obs(:),rig1(ij),rjg1(ij),gues3d(ij,ilev,mmean,iv3d_p),hgt1(ij,ilev),nv3d+n,hdxf,rdiag,rloc,dep,nobsl,depd=depd,nobsl_t=nobsl_t,cutd_t=cutd_t,srch_q0=search_q0(:,nv3d+1,ij))
              else                                                                       !GYL
                CALL obs_local(obs(:),rig1(ij),rjg1(ij),gues3d(ij,ilev,mmean,iv3d_p),hgt1(ij,ilev),nv3d+n,hdxf,rdiag,rloc,dep,nobsl,nobsl_t=nobsl_t,cutd_t=cutd_t,srch_q0=search_q0(:,nv3d+1,ij))
              end if                                                                     !GYL
              IF(RELAX_ALPHA_SPREAD /= 0.0d0) THEN                                       !GYL
                if (ENS_WITH_MDET) then                                                  !GYL
                  CALL letkf_core(nmem,nobstotal,nobsl,hdxf,rdiag,rloc,dep,work2d(ij,n), & !GYL
                                  trans(:,:,n2nc),transm=transm(:,n2nc),pao=pa(:,:,n2nc), & !GYL
                                  rdiag_wloc=.true.,infl_update=INFL_MUL_ADAPTIVE, &     !GYL
                                  depd=depd,transmd=transmd(:,n2nc))                     !GYL
                else                                                                     !GYL
                  CALL letkf_core(nmem,nobstotal,nobsl,hdxf,rdiag,rloc,dep,work2d(ij,n), & !GYL
                                  trans(:,:,n2nc),transm=transm(:,n2nc),pao=pa(:,:,n2nc), & !GYL
                                  rdiag_wloc=.true.,infl_update=INFL_MUL_ADAPTIVE)       !GYL
                end if                                                                   !GYL
              ELSE                                                                       !GYL
                if (ENS_WITH_MDET) then                                                  !GYL
                  CALL letkf_core(nmem,nobstotal,nobsl,hdxf,rdiag,rloc,dep,work2d(ij,n), & !GYL
                                  trans(:,:,n2nc),transm=transm(:,n2nc),           &     !GYL
                                  rdiag_wloc=.true.,infl_update=INFL_MUL_ADAPTIVE, &     !GYL
                                  depd=depd,transmd=transmd(:,n2nc))                     !GYL
                else                                                                     !GYL
                  CALL letkf_core(nmem,nobstotal,nobsl,hdxf,rdiag,rloc,dep,work2d(ij,n), & !GYL
                                  trans(:,:,n2nc),transm=transm(:,n2nc),           &     !GYL
                                  rdiag_wloc=.true.,infl_update=INFL_MUL_ADAPTIVE)       !GYL
                end if                                                                   !GYL
              END IF                                                                     !GYL
              trans_done(n2nc) = .true.                                                  !GYL
              IF(NOBS_OUT) THEN                                                          !GYL
                work2dn(:,ij,n) = real(sum(nobsl_t,dim=1),RP)                        !GYL !!! NOBS: sum over all variables for each report type
              END IF                                                                     !GYL

            END IF

            ! relaxation via LETKF weight
            IF(RELAX_ALPHA /= 0.0d0) THEN                                              !GYL - RTPP method (Zhang et al. 2004)
              CALL weight_RTPP(trans(:,:,n2nc),parm,transrlx)                          !GYL
            ELSE IF(RELAX_ALPHA_SPREAD /= 0.0d0) THEN                                  !GYL - RTPS method (Whitaker and Hamill 2012)
              IF(RELAX_SPREAD_OUT) THEN                                                !GYL
                CALL weight_RTPS(trans(:,:,n2nc),pa(:,:,n2nc),gues2d(ij,:,n), &        !GYL
                                 parm,transrlx,work2da(ij,n))                          !GYL
              ELSE                                                                     !GYL
                CALL weight_RTPS(trans(:,:,n2nc),pa(:,:,n2nc),gues2d(ij,:,n), &        !GYL
                                 parm,transrlx,tmpinfl)                                !GYL
              END IF                                                                   !GYL
            ELSE                                                                       !GYL
              transrlx = trans(:,:,n2nc)                                               !GYL - No relaxation
            END IF                                                                     !GYL

            ! total weight matrix
            do m = 1, nmem
              do k = 1, nmem
                transrlx(k,m) = (transrlx(k,m) + transm(k,n2nc)) * beta
              end do
              transrlx(m,m) = transrlx(m,m) + (1.0_RP-beta)
            end do

            ! analysis update of members
            do m = 1, nmem
              anal2d(ij,m,n) = gues2d(ij,mmean,n)
              do k = 1, nmem
                anal2d(ij,m,n) = anal2d(ij,m,n) + gues2d(ij,k,n) * transrlx(k,m)
              end do
            end do

            ! analysis update of deterministic run
            if (ENS_WITH_MDET) then
              anal2d(ij,mmdet,n) = 0.0_RP
              do k = 1, nmem
                anal2d(ij,mmdet,n) = anal2d(ij,mmdet,n) + gues2d(ij,k,n) * transmd(k,n2nc)
              end do
              anal2d(ij,mmdet,n) = gues2d(ij,mmdet,n) + anal2d(ij,mmdet,n) * beta
            end if

          END DO
        END IF
      END DO
    END DO ! [ ilev=1,nlev ]

    deallocate (hdxf,rdiag,rloc,dep)
    if (ENS_WITH_MDET) then
      deallocate (depd)
    end if
    deallocate (trans,transm,transmd,pa)

    deallocate (n_merge,ic_merge)
    deallocate (search_q0)

    IF (allocated(work3dg)) deallocate (work3dg)
    IF (allocated(work2dg)) deallocate (work2dg)

    !
    ! Additive inflation
    !
    IF(INFL_ADD > 0.0d0) THEN

      if (INFL_ADD_Q_RATIO) then
        work3d(:,:,:) = gues3d(:,:,mmean,:)
      else
        work3d(:,:,:) = 1.0d0
      end if

      allocate (addinfl_weight(nij1))
      if (INFL_ADD_REF_ONLY) then
        addinfl_weight(:) = 0.0d0
        ic = ctype_elmtyp(uid_obs(id_radar_ref_obs), 22)
        if (ic > 0) then
          do ij = 1, nij1
            ref_min_dist = 1.0d33
            !!!!!! save this (ref_min_dist) information when doing DA
            do iob = obsgrd(ic)%ac_ext(0, 1) + 1, obsgrd(ic)%ac_ext(obsgrd(ic)%ngrdext_i, obsgrd(ic)%ngrdext_j)
              rdx = (rig1(ij) - obs(obsda_sort%set(iob))%ri(obsda_sort%idx(iob))) * DX
              rdy = (rjg1(ij) - obs(obsda_sort%set(iob))%rj(obsda_sort%idx(iob))) * DY
              rdxy = rdx*rdx + rdy*rdy
              if (rdxy < ref_min_dist) then
                ref_min_dist = rdxy
              end if
            end do

            ref_min_dist = ref_min_dist / (hori_loc_ctype(ic) * hori_loc_ctype(ic))
            if (ref_min_dist <= dist_zero_fac_square) then
              addinfl_weight(ij) = EXP(-0.5d0 * ref_min_dist)
            end if
          end do
        end if
      else
        addinfl_weight(:) = 1.0d0
      end if

      LOG_INFO("LETKF_system",*) 'Additive covariance inflation, parameter:', INFL_ADD
      if (INFL_ADD_SHUFFLE) then
        if (RANK_ENS== 0) then
          call Knuth_Shuffle(nens, ishuf)
        end if
        call MPI_BCAST(ishuf, nens, MPI_INTEGER, 0, COMM_ENS, ierr)
        LOG_INFO("LETKF_system",*) 'shuffle members: on'
        LOG_INFO("LETKF_system",*) 'shuffle sequence:', ishuf
      end if

      DO n=1,nv3d
      DO m=1,nens
        if (INFL_ADD_SHUFFLE) then
          mshuf = ishuf(m)
        else
          mshuf = m
        end if
        if (n == iv3d_q .or. n == iv3d_qc .or. n == iv3d_qr .or. n == iv3d_qi .or. n == iv3d_qs .or. n == iv3d_qg) then
          DO k=1,nlev
          DO i=1,nij1
            anal3d(i,k,m,n) = anal3d(i,k,m,n) &
              & + addi3d(i,k,mshuf,n) * INFL_ADD * addinfl_weight(i) * work3d(i,k,n)
          END DO
          END DO
        else
          DO k=1,nlev
          DO i=1,nij1
            anal3d(i,k,m,n) = anal3d(i,k,m,n) &
              & + addi3d(i,k,mshuf,n) * INFL_ADD * addinfl_weight(i)
          END DO
          END DO
        end if
      END DO
      END DO
      DO n=1,nv2d
      DO m=1,nens
        if (INFL_ADD_SHUFFLE) then
          mshuf = ishuf(m)
        else
          mshuf = m
        end if
        DO i=1,nij1
          anal2d(i,m,n) = anal2d(i,m,n) + addi2d(i,mshuf,n) * INFL_ADD * addinfl_weight(i)
        END DO
      END DO
      END DO

      deallocate (addinfl_weight)

    END IF

    ! -- prepare the output values
    call gather_grd_mpi_all2all( 1, nens, anal3d(:,:,1:nens,:), anal2d(:,1:nens,:), v3dg, v2dg )

    do j = 1, nlat
    do i = 1, nlon
    do k = 1, nlev
       U   (k,i,j) = v3dg(k,i,j,iv3d_u )
       V   (k,i,j) = v3dg(k,i,j,iv3d_v )
       W   (k,i,j) = v3dg(k,i,j,iv3d_w )
       TEMP(k,i,j) = v3dg(k,i,j,iv3d_t )
       PRES(k,i,j) = v3dg(k,i,j,iv3d_p )
       QV  (k,i,j) = v3dg(k,i,j,iv3d_q )
       QC  (k,i,j) = v3dg(k,i,j,iv3d_qc)
       QR  (k,i,j) = v3dg(k,i,j,iv3d_qr)
       QI  (k,i,j) = v3dg(k,i,j,iv3d_qi)
       QS  (k,i,j) = v3dg(k,i,j,iv3d_qs)
       QG  (k,i,j) = v3dg(k,i,j,iv3d_qg)
    end do
    end do
    end do

    if( LETKF_DEBUG_LOG ) then
       call monit_obs_mpi( OBS_IN_NUM, OBS_IN_FORMAT, v3dg, v2dg, monit_step=2 )
    endif

    return
  end subroutine LETKF_system

  !-----------------------------------------------------------------------
  ! Data Assimilation (Parameter estmation for global constant parameters)
  !-----------------------------------------------------------------------
  subroutine LETKF_param_estimation_system( &
      PEST_PMAX, &
      PEST_VAR0  )
    use scale_const, only: &
      UNDEF => CONST_UNDEF
    use scale_statistics, only: &
      AVERAGE => STATISTICS_AVERAGE
    implicit none

    integer, intent(in) :: PEST_PMAX

    real(RP), intent(inout) :: PEST_VAR0(nens,PEST_PMAX)

    real(RP) :: gues0d(nens+1,PEST_PMAX)
    real(RP) :: anal0d(nens+1,PEST_PMAX)

    real(RP) :: work0d(PEST_PMAX)

    real(RP),allocatable :: hdxf(:,:)
    real(RP),allocatable :: rdiag(:)
    real(RP),allocatable :: rloc(:)
    real(RP),allocatable :: dep(:)
    real(RP),allocatable :: depd(:)            !gyl

    real(RP) :: parm
    real(RP),allocatable :: trans(:,:)
    real(RP),allocatable :: transm(:)
    real(RP),allocatable :: transmd(:)
    real(RP),allocatable :: pa(:,:)
    real(RP) :: transrlx(nmem,nmem)

    integer :: n,m,k,nobsl
    real(RP) :: beta                           !gyl
    real(RP) :: tmpinfl                        !gyl
    integer :: ierr

    LOG_PROGRESS(*) 'data-assimilation / LETKF / parameter estimation'

    ! -- prepare the first-guess data
    do n = 1, PEST_PMAX
    do m = 1, nens
      gues0d(m,n) = PEST_VAR0(m,n)
    end do
    end do

    ! -- obtain the ensemble mean
    do n = 1, PEST_PMAX
      gues0d(mmean,n) = AVERAGE( gues0d(1:nmem,n), UNDEF )
    end do

    !
    ! No variable localization
    !

    !
    ! FCST PERTURBATIONS
    !

    do n = 1, PEST_PMAX
    do m = 1, nmem
      gues0d(m,n) = gues0d(m,n) - gues0d(mmean,n)
    end do
    end do

    allocate (hdxf (nobstotal,nmem))
    allocate (rdiag(nobstotal))
    allocate (rloc (nobstotal))
    allocate (dep  (nobstotal))
    if (ENS_WITH_MDET) then
      allocate (depd (nobstotal))
      allocate (transmd(nmem))
    end if
    allocate (trans (nmem,nmem))
    allocate (transm(nmem))
    allocate (pa    (nmem,nmem))

    !
    ! MAIN ASSIMILATION LOOP
    !
    parm = 1.0_RP
    beta = 1.0_RP

    work0d = 1.0_RP

    DO n = 1, PEST_PMAX
      ! calculate mean and perturbation weights
      ! compute weights with localized observations
      if (ENS_WITH_MDET) then
        CALL obs_pest_etkf(hdxf, rdiag, rloc, dep, nobsl, depd=depd)
      else
        CALL obs_pest_etkf(hdxf, rdiag, rloc, dep, nobsl)
      end if

      if (ENS_WITH_MDET) then
        CALL letkf_core(nmem,nobstotal,nobsl,hdxf,rdiag,rloc,dep,work0d(n), &
                        trans(:,:),transm=transm(:),pao=pa(:,:), &
                        rdiag_wloc=.true.,&
                        depd=depd,transmd=transmd(:))
      else
        CALL letkf_core(nmem,nobstotal,nobsl,hdxf,rdiag,rloc,dep,work0d(n), &
                        trans(:,:),transm=transm(:),pao=pa(:,:), &
                        rdiag_wloc=.true.)
      end if

      ! relaxation via LETKF weight
      CALL weight_RTPS_const(trans(:,:),pa(:,:),gues0d(:,n),transrlx,tmpinfl)

      ! total weight matrix
      do m = 1, nmem
        do k = 1, nmem
          transrlx(k,m) = (transrlx(k,m) + transm(k)) * beta
        end do
        transrlx(m,m) = transrlx(m,m) + (1.0_RP-beta)
      end do

      ! analysis update of members
      do m = 1, nmem
        anal0d(m,n) = gues0d(mmean,n)
        do k = 1, nmem
          anal0d(m,n) = anal0d(m,n) + gues0d(k,n) * transrlx(k,m)
        end do
      end do

      ! analysis update of deterministic run
      if (ENS_WITH_MDET) then
        anal0d(mmdet,n) = 0.0_RP
        do k = 1, nmem
          anal0d(mmdet,n) = anal0d(mmdet,n) + gues0d(k,n) * transmd(k)
        end do
        anal0d(mmdet,n) = gues0d(mmdet,n) + anal0d(mmdet,n) * beta
      end if

    enddo

    do n = 1, PEST_PMAX
    do m = 1, nens
      PEST_VAR0(m,n) = anal0d(m,n)
    end do
    end do

    deallocate (hdxf,rdiag,rloc,dep)
    if (ENS_WITH_MDET) then
      deallocate (depd,transmd)
    end if
    deallocate (trans,transm,pa)

    return
  end subroutine LETKF_param_estimation_system

  !-----------------------------------------------------------------------
  ! Setup of additive inflation
  !-----------------------------------------------------------------------
  subroutine LETKF_add_inflation_setup( &
      addi3d, &
      addi2d  )
    implicit none

    real(RP), intent(out) :: addi3d(:,:,:,:)
    real(RP), intent(out) :: addi2d(:,:,:)
    integer :: i, k, m, n

    call read_ens_mpi_addiinfl(addi3d, addi2d)

    call ensmean_grd(nens, nens, nij1, addi3d, addi2d)

    do n = 1, nv3d
    do m = 1, nens
    do k = 1, nlev
    do i = 1, nij1
      addi3d(i,k,m,n) = addi3d(i,k,m,n) - addi3d(i,k,mmean,n)
    end do
    end do
    end do
    end do

    do n = 1, nv2d
    do m = 1, nens
    do i = 1, nij1
      addi2d(i,m,n) = addi2d(i,m,n) - addi2d(i,mmean,n)
    end do
    end do
    end do

    return
  end subroutine LETKF_add_inflation_setup

  ! =======================================================================
  !  Main Subroutine of LETKF Core
  !   INPUT
  !     ne               : ensemble size                                           !GYL
  !     nobs             : array size, but only first nobsl elements are used
  !     nobsl            : total number of observation assimilated at the point
  !     hdxb(nobs,ne)    : obs operator times fcst ens perturbations
  !     rdiag(nobs)      : observation error variance
  !     rloc(nobs)       : localization weighting function
  !     dep(nobs)        : observation departure (yo-Hxb)
  !     parm_infl        : covariance inflation parameter
  !     rdiag_wloc       : (optional) flag indicating that rdiag = rdiag / rloc    !GYL
  !     infl_update      : (optional) flag to return updated inflation parameter   !GYL
  !     depd(nobs)       : observation departure (yo-Hxb) for deterministic run    !GYL
  !   OUTPUT
  !     parm_infl        : updated covariance inflation parameter
  !     trans(ne,ne)     : transformation matrix
  !     transm(ne)       : (optional) transformation matrix mean                   !GYL
  !     pao(ne,ne)       : (optional) analysis covariance matrix in ensemble space !GYL
  !     transmd(ne)      : (optional) transformation matrix mean for deterministic run !GYL
  ! =======================================================================
  ! 
  !  [PURPOSE:] Local Ensemble Transform Kalman Filtering (LETKF)
  !             Model Independent Core Module
  ! 
  !  [REFERENCES:]
  !   [1] Ott et al., 2004: A local ensemble Kalman filter for atmospheric
  !     data assimilation. Tellus, 56A, 415-428.
  !   [2] Hunt et al., 2007: Efficient Data Assimilation for Spatiotemporal
  !     Chaos: A Local Ensemble Transform Kalman Filter. Physica D, 230,
  !     112-126.
  ! 
  !  [HISTORY:]
  !   01/21/2009 Takemasa Miyoshi  Created at U. of Maryland, College Park
  ! 
  ! =======================================================================
  subroutine letkf_core(ne,nobs,nobsl,hdxb,rdiag,rloc,dep,parm_infl,trans,transm,pao,rdiag_wloc,infl_update,depd,transmd)
    use scale_sort, only: &
      SORT_quickselect_desc_arg, &
      SORT_quickselect_arg
    use scale_matrix, only: &
      MATRIX_SOLVER_eigenvalue_decomposition
    implicit none
    INTEGER,INTENT(IN) :: ne                      !GYL
    INTEGER,INTENT(IN) :: nobs
    INTEGER,INTENT(IN) :: nobsl
    REAL(RP),INTENT(IN) :: hdxb(1:nobs,1:ne)
    REAL(RP),INTENT(IN) :: rdiag(1:nobs)
    REAL(RP),INTENT(IN) :: rloc(1:nobs)
    REAL(RP),INTENT(IN) :: dep(1:nobs)
    REAL(RP),INTENT(INOUT) :: parm_infl
    REAL(RP),INTENT(OUT) :: trans(ne,ne)
    REAL(RP),INTENT(OUT),OPTIONAL :: transm(ne)
    REAL(RP),INTENT(OUT),OPTIONAL :: pao(ne,ne)
    LOGICAL,INTENT(IN),OPTIONAL :: rdiag_wloc     !GYL
    LOGICAL,INTENT(IN),OPTIONAL :: infl_update    !GYL
    REAL(RP),INTENT(IN),OPTIONAL :: depd(1:nobs) !GYL
    REAL(RP),INTENT(OUT),OPTIONAL :: transmd(ne)  !GYL

    REAL(RP) :: hdxb_rinv(nobsl,ne)
    REAL(RP) :: eivec(ne,ne)
    REAL(RP) :: eival(ne)
    REAL(RP) :: pa(ne,ne)
    REAL(RP) :: work1(ne,ne)
    REAL(RP) :: work2(ne)
    REAL(RP) :: work2d(ne)
    REAL(RP) :: work3(ne)
    REAL(RP) :: rho
    REAL(RP) :: parm(4),sigma_o,gain
    REAL(RP),PARAMETER :: sigma_b = 0.04d0 !error stdev of parm_infl
    LOGICAL :: rdiag_wloc_
    LOGICAL :: infl_update_
    INTEGER :: i,j,k

    real(RP) :: work1_evp(ne,ne)
    real(RP) :: eivec_evp(ne,ne)
    real(RP) :: eival_evp(ne)

#ifdef DA
    rdiag_wloc_ = .FALSE.                               !GYL
    IF(present(rdiag_wloc)) rdiag_wloc_ = rdiag_wloc    !GYL
    infl_update_ = .FALSE.                              !GYL
    IF(present(infl_update)) infl_update_ = infl_update !GYL

    IF(nobsl == 0) THEN
      trans = 0.0
      DO i=1,ne
        trans(i,i) = SQRT(parm_infl)
      END DO
      IF (PRESENT(transm)) THEN
        transm = 0.0
      END IF
      IF (PRESENT(transmd)) THEN
        transmd = 0.0
      END IF
      IF (PRESENT(pao)) THEN
        pao = 0.0
        DO i=1,ne
          pao(i,i) = parm_infl / REAL(ne-1,kind=RP)
        END DO
      END IF
      RETURN
    END IF
  !-----------------------------------------------------------------------
  !  Rinv hdxb
  !-----------------------------------------------------------------------
    IF(rdiag_wloc_) THEN                            !GYL
      DO j=1,ne                                     !GYL
        DO i=1,nobsl                                !GYL
          hdxb_rinv(i,j) = hdxb(i,j) / rdiag(i)     !GYL
        END DO                                      !GYL
      END DO                                        !GYL
    ELSE                                            !GYL
      DO j=1,ne
        DO i=1,nobsl
          hdxb_rinv(i,j) = hdxb(i,j) / rdiag(i) * rloc(i)
        END DO
      END DO
    END IF                                          !GYL

  !-----------------------------------------------------------------------
  !  hdxb^T Rinv hdxb
  !-----------------------------------------------------------------------
    if( RP == SP ) then
      CALL sgemm('t','n',ne,ne,nobsl,1.0e0,hdxb_rinv,nobsl,hdxb(1:nobsl,:),nobsl,0.0e0,work1,ne)
    else
      CALL dgemm('t','n',ne,ne,nobsl,1.0d0,hdxb_rinv,nobsl,hdxb(1:nobsl,:),nobsl,0.0d0,work1,ne)
    end if
  !  DO j=1,ne
  !    DO i=1,ne
  !      work1(i,j) = hdxb_rinv(1,i) * hdxb(1,j)
  !      DO k=2,nobsl
  !        work1(i,j) = work1(i,j) + hdxb_rinv(k,i) * hdxb(k,j)
  !      END DO
  !    END DO
  !  END DO
  !-----------------------------------------------------------------------
  !  hdxb^T Rinv hdxb + (m-1) I / rho (covariance inflation)
  !-----------------------------------------------------------------------
    rho = 1.0d0 / parm_infl
    DO i=1,ne
      work1(i,i) = work1(i,i) + REAL(ne-1,kind=RP) * rho
    END DO
  !-----------------------------------------------------------------------
  !  eigenvalues and eigenvectors of [ hdxb^T Rinv hdxb + (m-1) I ]
  !-----------------------------------------------------------------------
    work1_evp = real( work1, kind=RP )
    call MATRIX_SOLVER_eigenvalue_decomposition( ne, work1_evp, eival_evp, eivec_evp )
    eival = real( eival_evp, kind=RP )
    eivec = real( eivec_evp, kind=RP )
  !-----------------------------------------------------------------------
  !  Pa = [ hdxb^T Rinv hdxb + (m-1) I ]inv
  !-----------------------------------------------------------------------
    DO j=1,ne
      DO i=1,ne
        work1(i,j) = eivec(i,j) / eival(j)
      END DO
    END DO
    if( RP == SP ) then
      CALL sgemm('n','t',ne,ne,ne,1.0e0,work1,ne,eivec,ne,0.0e0,pa,ne)
    else
      CALL dgemm('n','t',ne,ne,ne,1.0d0,work1,ne,eivec,ne,0.0d0,pa,ne)
    end if
  !  DO j=1,ne
  !    DO i=1,ne
  !      pa(i,j) = work1(i,1) * eivec(j,1)
  !      DO k=2,ne
  !        pa(i,j) = pa(i,j) + work1(i,k) * eivec(j,k)
  !      END DO
  !    END DO
  !  END DO
  !-----------------------------------------------------------------------
  !  hdxb_rinv^T dep
  !-----------------------------------------------------------------------
    if( RP == SP ) then
      CALL sgemv('t',nobsl,ne,1.0e0,hdxb_rinv,nobsl,dep,1,0.0e0,work2,1)
    else
      CALL dgemv('t',nobsl,ne,1.0d0,hdxb_rinv,nobsl,dep,1,0.0d0,work2,1)
    end if
  !  DO i=1,ne
  !    work2(i) = hdxb_rinv(1,i) * dep(1)
  !    DO j=2,nobsl
  !      work2(i) = work2(i) + hdxb_rinv(j,i) * dep(j)
  !    END DO
  !  END DO
    IF (PRESENT(depd) .AND. PRESENT(transmd)) THEN
      if( RP == SP ) then
        CALL sgemv('t',nobsl,ne,1.0e0,hdxb_rinv,nobsl,depd,1,0.0e0,work2d,1)
      else
        CALL dgemv('t',nobsl,ne,1.0d0,hdxb_rinv,nobsl,depd,1,0.0d0,work2d,1)
      end if
  !    DO i=1,ne
  !       work2d(i) = hdxb_rinv(1,i) * depd(1)
  !      DO j=2,nobsl
  !        work2d(i) = work2d(i) + hdxb_rinv(j,i) * depd(j)
  !      END DO
  !    END DO
    END IF
  !-----------------------------------------------------------------------
  !  Pa hdxb_rinv^T dep
  !-----------------------------------------------------------------------
    if( RP == SP ) then
      CALL sgemv('n',ne,ne,1.0e0,pa,ne,work2,1,0.0e0,work3,1)
    else
      CALL dgemv('n',ne,ne,1.0d0,pa,ne,work2,1,0.0d0,work3,1)
    end if
  !  DO i=1,ne
  !    work3(i) = 0.0_RP
  !  end DO
  !  !$omp parallel do reduction(+:work3)
  !  DO j=1,ne
  !     DO i=1,ne
  !        work3(i) = work3(i) + pa(i,j) * work2(j)
  !    END DO
  !  END DO
    IF (PRESENT(depd) .AND. PRESENT(transmd)) THEN
      if( RP == SP ) then
        CALL sgemv('n',ne,ne,1.0e0,pa,ne,work2d,1,0.0e0,transmd,1)
      else
        CALL dgemv('n',ne,ne,1.0d0,pa,ne,work2d,1,0.0d0,transmd,1)
      end if
  !    DO i=1,ne
  !      transmd(i) = 0.0_RP
  !    end DO
  !    !$omp parallel do reduction(+:transmd)
  !    DO j=1,ne
  !       DO i=1,ne
  !          transmd(i) = transmd(i) + pa(i,j) * work2d(j)
  !       END DO
  !    END DO
    END IF

  !!$!-----------------------------------------------------------------------
  !!$!  Pa hdxb_rinv^T
  !!$!-----------------------------------------------------------------------
  !!$#ifdef SINGLELETKF
  !!$  CALL sgemm('n','t',ne,nobsl,ne,1.0e0,pa,ne,hdxb_rinv,&
  !!$    & nobsl,0.0e0,work2,ne)
  !!$#else
  !!$  CALL dgemm('n','t',ne,nobsl,ne,1.0d0,pa,ne,hdxb_rinv,&
  !!$    & nobsl,0.0d0,work2,ne)
  !!$#endif
  !!$!  DO j=1,nobsl
  !!$!    DO i=1,ne
  !!$!      work2(i,j) = pa(i,1) * hdxb_rinv(j,1)
  !!$!      DO k=2,ne
  !!$!        work2(i,j) = work2(i,j) + pa(i,k) * hdxb_rinv(j,k)
  !!$!      END DO
  !!$!    END DO
  !!$!  END DO
  !!$!-----------------------------------------------------------------------
  !!$!  Pa hdxb_rinv^T dep
  !!$!-----------------------------------------------------------------------
  !!$  DO i=1,ne
  !!$    work3(i) = work2(i,1) * dep(1)
  !!$    DO j=2,nobsl
  !!$      work3(i) = work3(i) + work2(i,j) * dep(j)
  !!$    END DO
  !!$  END DO
  !!$  IF (PRESENT(depd) .AND. PRESENT(transmd)) THEN       !GYL
  !!$    DO i=1,ne                                          !GYL
  !!$      transmd(i) = work2(i,1) * depd(1)                !GYL
  !!$      DO j=2,nobsl                                     !GYL
  !!$        transmd(i) = transmd(i) + work2(i,j) * depd(j) !GYL
  !!$      END DO                                           !GYL
  !!$    END DO                                             !GYL
  !!$  END IF                                               !GYL
  !-----------------------------------------------------------------------
  !  T = sqrt[(m-1)Pa]
  !-----------------------------------------------------------------------
    DO j=1,ne
      rho = SQRT( REAL(ne-1,kind=RP) / eival(j) )
      DO i=1,ne
        work1(i,j) = eivec(i,j) * rho
      END DO
    END DO
    if( RP == SP ) then
      CALL sgemm('n','t',ne,ne,ne,1.0e0,work1,ne,eivec,ne,0.0e0,trans,ne)
    else
      CALL dgemm('n','t',ne,ne,ne,1.0d0,work1,ne,eivec,ne,0.0d0,trans,ne)
    end if
  !  DO j=1,ne
  !    DO i=1,ne
  !      trans(i,j) = work1(i,1) * eivec(j,1)
  !      DO k=2,ne
  !        trans(i,j) = trans(i,j) + work1(i,k) * eivec(j,k)
  !      END DO
  !    END DO
  !  END DO
  !-----------------------------------------------------------------------
  !  T + Pa hdxb_rinv^T dep
  !-----------------------------------------------------------------------
    IF (PRESENT(transm)) THEN                !GYL - if transm is present,
      transm = work3                         !GYL - return both trans and transm without adding them
    ELSE                                     !GYL
      DO j=1,ne
        DO i=1,ne
          trans(i,j) = trans(i,j) + work3(i)
        END DO
      END DO
    END IF                                   !GYL
    IF (PRESENT(pao)) pao = pa               !GYL

    IF (.NOT. infl_update_) RETURN           !GYL - skip the following if no inflation update is required
  !-----------------------------------------------------------------------
  !  Inflation estimation
  !-----------------------------------------------------------------------
    parm = 0.0d0
    IF(rdiag_wloc_) THEN                            !GYL
      DO i=1,nobsl                                  !GYL
        parm(1) = parm(1) + dep(i)*dep(i)/rdiag(i)  !GYL
      END DO                                        !GYL
    ELSE                                            !GYL
      DO i=1,nobsl
        parm(1) = parm(1) + dep(i)*dep(i)/rdiag(i) * rloc(i)
      END DO
    END IF                                          !GYL
    DO j=1,ne
      DO i=1,nobsl
        parm(2) = parm(2) + hdxb_rinv(i,j) * hdxb(i,j)
      END DO
    END DO
    parm(2) = parm(2) / REAL(ne-1,kind=RP)
    parm(3) = SUM(rloc(1:nobsl))
    parm(4) = (parm(1)-parm(3))/parm(2) - parm_infl
  !  sigma_o = 1.0d0/REAL(nobsl,kind=RP)/MAXVAL(rloc(1:nobsl))
    sigma_o = 2.0d0/parm(3)*((parm_infl*parm(2)+parm(3))/parm(2))**2
    gain = sigma_b**2 / (sigma_o + sigma_b**2)
    parm_infl = parm_infl + gain * parm(4)
#endif

    return
  end subroutine letkf_core

  !-----------------------------------------------------------------------
  ! Basic modules for observation input
  !-----------------------------------------------------------------------
  subroutine get_nobs(cfile,nrec,nn)
    implicit none

    character(*),intent(in) :: cfile
    integer,intent(in) :: nrec
    integer,intent(out) :: nn
    integer :: iunit
    logical :: ex
    integer :: sz

    nn = 0
    iunit = IO_get_available_fid()

    inquire(file=cfile,exist=ex)
    if(ex) then
      open(unit=iunit,file=cfile,form='unformatted',access='stream')
      inquire(unit=iunit, size=sz)
      if (mod(sz, SP * (nrec+2) ) /= 0) then
        LOG_ERROR('get_nobs',*) ': reading error -- skipped: ', cfile
        close(unit=iunit)
        return
      end if
      nn = sz / (SP * (nrec+2) )
      close(unit=iunit)
    else
      LOG_INFO('get_nobs',*) 'file not found -- skipped: ', cfile
    end if

    return
  end subroutine get_nobs

  subroutine get_nobs_radar(cfile,nn,radarlon,radarlat,radarz)
    implicit none
    character(*),intent(in) :: cfile
    integer,intent(out) :: nn
    real(RP),intent(out) :: radarlon,radarlat,radarz
    integer :: nrec
    real(SP) :: tmp
    integer :: ios
    integer :: iunit
    logical :: ex
    integer :: sz

    if(RADAR_OBS_4D) then
      nrec = 8
    else
      nrec = 7
    end if
    nn = 0
    iunit = IO_get_available_fid()

    radarlon = 0.0
    radarlat = 0.0
    radarz   = 0.0

    inquire(file=cfile,exist=ex)
    if(ex) then
      open(unit=iunit,file=cfile,form='unformatted',access='stream')
      read(unit=iunit,iostat=ios) tmp ! dummy
      read(unit=iunit,iostat=ios) tmp
      if(ios /= 0) then
        LOG_ERROR('get_nobs_radar',*) ': reading error -- skipped: ', cfile
        close(unit=iunit)
        return
      end if
      radarlon = real(tmp,kind=RP)
      read(unit=iunit,iostat=ios) tmp ! dummy
      read(unit=iunit,iostat=ios) tmp ! dummy
      read(unit=iunit,iostat=ios) tmp
      if(ios /= 0) then
        LOG_ERROR('get_nobs_radar',*) ': reading error -- skipped: ', cfile
        close(unit=iunit)
        return
      end if
      radarlat = real(tmp,kind=RP)
      read(unit=iunit,iostat=ios) tmp ! dummy
      read(unit=iunit,iostat=ios) tmp ! dummy
      read(unit=iunit,iostat=ios) tmp
      if(ios /= 0) then
        LOG_ERROR('get_nobs_radar',*) ': reading error -- skipped: ', cfile
        close(unit=iunit)
        return
      end if
      radarz = real(tmp,kind=RP)
      read(unit=iunit,iostat=ios) tmp ! dummy

      ! get file size by INQUIRE statement... may not work for some older fortran compilers
      inquire(unit=iunit, size=sz)
      sz = sz - SP * (1+2) * 3 ! substract the radar data header
      if (mod(sz, SP * (nrec+2)) /= 0) then
        LOG_ERROR('get_nobs_radar',*) ': reading error -- skipped: ', cfile
        close(unit=iunit)
        return
      end if
      nn = sz / (SP * (nrec+2))

      close(unit=iunit)
    else
      LOG_INFO('get_nobs_radar',*) 'file not found -- skipped: ', cfile
    end if

    return
  end subroutine get_nobs_radar

  subroutine read_obs(cfile,obs)
    use scale_const, only: &
        PI => CONST_PI
    use scale_mapprojection, only: &
        mapprojection_lonlat2xy
    implicit none

    character(*),intent(in) :: cfile
    type(obs_info),intent(inout) :: obs
    real(SP) :: wk(8)
    real(RP) :: x_RP, y_RP
    real(SP) :: tmp
    integer :: n,iunit

    iunit = IO_get_available_fid()
    open(unit=iunit,file=cfile,form='unformatted',access='stream')
    do n = 1, obs%nobs
      read(unit=iunit) tmp
      read(unit=iunit) wk(:)
      read(unit=iunit) tmp
      select case(nint(wk(1)))
      case(id_u_obs)
        wk(4) = wk(4) * 100.0 ! hpa -> pa
      case(id_v_obs)
        wk(4) = wk(4) * 100.0 ! hpa -> pa
      case(id_t_obs)
        wk(4) = wk(4) * 100.0 ! hpa -> pa
      case(id_tv_obs)
        wk(4) = wk(4) * 100.0 ! hpa -> pa
      case(id_q_obs)
        wk(4) = wk(4) * 100.0 ! hpa -> pa
      case(id_ps_obs)
        wk(5) = wk(5) * 100.0 ! hpa -> pa
        wk(6) = wk(6) * 100.0 ! hpa -> pa
      case(id_rh_obs)
        wk(4) = wk(4) * 100.0 ! hpa -> pa
        wk(5) = wk(5) * 0.01 ! percent input
        wk(6) = wk(6) * 0.01 ! percent input
      case(id_tcmip_obs)
        wk(4) = wk(4) * 100.0 ! hpa -> pa
        wk(5) = wk(5) * 100.0 ! hpa -> pa
        wk(6) = real( obserr_tcp,kind=SP )
      case(id_tclon_obs)
        call mapprojection_lonlat2xy( real(wk(2),kind=RP)*PI/180.0_RP,&
                                      real(wk(3),kind=RP)*PI/180.0_RP,&
                                      x_RP, y_RP )
        wk(4) = wk(4) * 100.0 ! hpa -> pa
        wk(5) = real( x_RP, kind=SP )
        wk(6) = real( obserr_tcx, kind=SP )
      case(id_tclat_obs)
        call mapprojection_lonlat2xy( real(wk(2),kind=RP)*PI/180.0_RP,&
                                      real(wk(3),kind=RP)*PI/180.0_RP,&
                                      x_RP, y_RP )
        wk(4) = wk(4) * 100.0 ! hpa -> pa
        wk(5) = real( y_RP, kind=SP )
        wk(6) = real( obserr_tcy, kind=SP )
      end select
      obs%elm(n) = nint( wk(1) )
      obs%lon(n) = real( wk(2), kind=RP )
      obs%lat(n) = real( wk(3), kind=RP )
      obs%lev(n) = real( wk(4), kind=RP )
      obs%dat(n) = real( wk(5), kind=RP )
      obs%err(n) = real( wk(6), kind=RP )
      obs%typ(n) = nint( wk(7) )
      obs%dif(n) = real( wk(8), kind=RP )
    end do
    close(unit=iunit)

    return
  end subroutine read_obs

  subroutine read_obs_radar(cfile,obs)
    implicit none
    character(*),intent(in) :: cfile
    type(obs_info),intent(inout) :: obs
    real(SP) :: wk(8)
    integer :: nrec
    real(SP) :: tmp
    integer :: n,iunit,ios

    if(RADAR_OBS_4D) then
      nrec = 8
    else
      nrec = 7
    end if
    iunit = IO_get_available_fid()
    open(unit=iunit,file=cfile,form='unformatted',access='stream')
    read(unit=iunit,iostat=ios) tmp ! dummy
    read(unit=iunit,iostat=ios) tmp ! radarlon
    read(unit=iunit,iostat=ios) tmp ! dummy
    read(unit=iunit,iostat=ios) tmp ! dummy
    read(unit=iunit,iostat=ios) tmp ! radarlat
    read(unit=iunit,iostat=ios) tmp ! dummy
    read(unit=iunit,iostat=ios) tmp ! dummy
    read(unit=iunit,iostat=ios) tmp ! radarz
    read(unit=iunit,iostat=ios) tmp ! dummy
    do n = 1, obs%nobs
      read(unit=iunit) tmp ! dummy
      read(unit=iunit) wk(1:nrec)
      read(unit=iunit) tmp ! dummy
      obs%elm(n) = nint(wk(1))
      obs%lon(n) = real(wk(2),kind=RP)
      obs%lat(n) = real(wk(3),kind=RP)
      obs%lev(n) = real(wk(4),kind=RP)
      obs%dat(n) = real(wk(5),kind=RP)
      obs%err(n) = real(wk(6),kind=RP)
      obs%typ(n) = 22
      if(RADAR_OBS_4D) then
        obs%dif(n) = real(wk(8),kind=RP)
      else
        obs%dif(n) = 0.0
      end if
    end do
    close(unit=iunit)

    return
  end subroutine read_obs_radar

  subroutine read_obs_radar_toshiba_pawr( obs, cfile )
    use iso_c_binding
    use scale_da_read_pawr_toshiba, only: &
      read_toshiba => DA_read_pawr_toshiba, &
      c_pawr_header, &
      RDIM,          &
      AZDIM,         &
      ELDIM
    implicit none

    character(len=*), intent(in) :: cfile
    type(obs_info), intent(out) :: obs

    type(obs_info) :: obs_ref

    integer, parameter :: n_type = 3
    character(len=4), parameter :: file_type_sfx(n_type) = &
      (/'.ze ', '.vr ', '.qcf'/)
    logical, parameter :: input_is_dbz = .true.

    type(c_pawr_header) :: hd(n_type)
    real(kind=c_float), allocatable, save :: rtdat(:, :, :, :)
    real(kind=c_float), allocatable, save :: az(:, :, :)
    real(kind=c_float), allocatable, save :: el(:, :, :)
    integer :: j, ierr, ierr2
    character(len=3) :: fname
    integer, save::i=0

    real(RP), allocatable :: ze(:, :, :), vr(:, :, :), qcflag(:, :, :), attenuation(:, :, :), rrange(:)
    real(RP), allocatable :: radlon(:, :, :), radlat(:, :, :), radz(:, :, :)
    real(RP), allocatable :: lon(:), lat(:), z(:)
    integer(8), allocatable :: grid_index(:), grid_count_ze(:), grid_count_vr(:)
    real(RP), allocatable :: grid_ze(:), grid_vr(:)
    real(RP), allocatable :: grid_lon_ze(:), grid_lat_ze(:), grid_z_ze(:)
    real(RP), allocatable :: grid_lon_vr(:),  grid_lat_vr(:),  grid_z_vr(:)

    character(len=1024) :: input_fname(n_type)
    integer ia, ir, ie
    real(RP) :: dlon, dlat
    integer(8) nobs_sp

    integer,save :: na, nr, ne
    real(RP),save :: lon0, lat0, z0
    real(RP),save :: missing
    integer,save :: range_res

    real(RP) :: max_obs_ze , min_obs_ze , max_obs_vr , min_obs_vr 
    integer :: nobs_ze, nobs_vr
    integer(8) :: idx, n, n_ref
    integer :: pos
    integer, parameter :: int1 = selected_int_kind(1) !1-BYTE INT
    integer(kind = int1) :: tmp_qcf, valid_qcf

    integer,parameter :: qcf_mask(8)=(/ 0, 1, 1, 1, 1, 0, 0, 0 /) !! valid, shadow, clutter possible, clutter certain, interference, range sidelobe /

    integer::qcf_count(0:255)

    character(len=90) :: plotname
    character(len=19) :: timelabel

    integer :: ii, jj, kk

    real(SP), allocatable :: ref3d(:,:,:)
    real(SP), allocatable :: vr3d(:,:,:)
    character(len=255) :: filename
    integer :: irec, iunit, iolen
    integer :: k

    integer ::  nlon, nlat, nlev

    RADAR_SO_SIZE_HORI = max( real( DX, kind=RP ), RADAR_SO_SIZE_HORI )
    RADAR_SO_SIZE_HORI = max( real( DY, kind=RP ), RADAR_SO_SIZE_HORI )

    if( .not. allocated(rtdat) ) allocate(rtdat(RDIM, AZDIM, ELDIM, n_type))
    if( .not. allocated(az)    ) allocate(az(AZDIM, ELDIM, n_type))
    if( .not. allocated(el)    ) allocate(el(AZDIM, ELDIM, n_type))

    do j = 1, n_type
      input_fname(j) = trim(cfile)
      call str_replace(input_fname(j), '.<type>', trim(file_type_sfx(j)), pos)
      ierr = read_toshiba( input_fname(j), hd(j), az(:,:,j), el(:,:,j), rtdat(:,:,:,j) )
      if( ierr /= 0 ) then
        LOG_INFO("LETKF_obs_readfile",*) 'Warning: File (',trim(input_fname(j)),') cannot be read. Skip.'
        obs%nobs = 0
        return
      endif
    end do

    ! Set obs information
    lon0 = hd(1)%longitude
    lat0 = hd(1)%latitude
    z0   = hd(1)%altitude

    missing   = real( hd(1)%mesh_offset, kind=RP )
    range_res = hd(1)%range_res

    nr = hd(1)%range_num
    na = hd(1)%sector_num
    ne = hd(1)%el_num

    call jst2utc( hd(1)%s_yr, hd(1)%s_mn, hd(1)%s_dy, hd(1)%s_hr, hd(1)%s_mi, hd(1)%s_sc, 0.0_DP, utime_obs )

    allocate( ze         (na, nr, ne) )
    allocate( vr         (na, nr, ne) )
    allocate( qcflag     (na, nr, ne) )
    allocate( attenuation(na, nr, ne) )

    valid_qcf = 0
    do j = 1, 8  
      if(qcf_mask(j) > 0) valid_qcf = ibset(valid_qcf, j - 1) 
    end do

    do ie = 1, ne
       do ir = 1, nr
       do ia = 1, na
         ze(ia,ir,ie) = rtdat(ir,ia,ie,1)
         vr(ia,ir,ie) = rtdat(ir,ia,ie,2)                                                                 

         if( vr(ia,ir,ie) > RADAR_MAX_ABS_VR .or. vr(ia,ir,ie) < -RADAR_MAX_ABS_VR ) vr(ia,ir,ie) = missing

         tmp_qcf = int(rtdat(ir, ia, ie, 3), int1)
         if(iand(valid_qcf, tmp_qcf) == 0) then      
           qcflag(ia, ir, ie) = 0.0_RP !valid
         else
           qcflag(ia, ir, ie) = 1000.0_RP !invalid
         end if
         attenuation(ia, ir, ie) = 1.0_RP !not implemented yet
       end do
       end do
    end do
    deallocate(rtdat)

    allocate(rrange(nr))
    do ir = 1, nr
       rrange(ir) = (dble(ir) - 0.5_RP) * range_res
    end do

    allocate( radlon(na, nr, ne) )
    allocate( radlat(na, nr, ne) )
    allocate( radz  (na, nr, ne) )

    call radar_georeference( lon0, lat0, z0, na, nr, ne,   & ! input
                             real( az(:, 1, 1), kind=RP ), & ! input (assume ordinary scan strategy)
                             rrange,                       & ! input (assume ordinary scan strategy)
                             real( el(1, :, 1), kind=RP ), & ! input (assume ordinary scan strategy)
                             radlon, radlat, radz          ) ! output  

    call define_grid( lon0, lat0, nr, rrange, rrange(nr), RADAR_ZMAX,             & ! input
                      RADAR_SO_SIZE_HORI, RADAR_SO_SIZE_HORI, RADAR_SO_SIZE_VERT, & ! input
                      dlon, dlat, nlon, nlat, nlev, lon, lat, z                   ) ! output

    call radar_superobing( na, nr, ne, radlon, radlat, radz, ze, vr,                      & ! input spherical
                           qcflag, attenuation,                                           & ! input spherical
                           nlon, nlat, nlev, lon, lat, z, dlon, dlat, RADAR_SO_SIZE_VERT, & ! input cartesian
                           missing, input_is_dbz,                                         & ! input param
                           lon0, lat0,                                                    &
                           nobs_sp, grid_index,                                           & ! output array info
                           grid_ze, grid_lon_ze, grid_lat_ze, grid_z_ze, grid_count_ze,   & ! output ze
                           grid_vr, grid_lon_vr, grid_lat_vr, grid_z_vr, grid_count_vr    ) ! output vr

    if( allocated( ze          ) ) deallocate(ze)
    if( allocated( vr          ) ) deallocate(vr)
    if( allocated( qcflag      ) ) deallocate(qcflag)
    if( allocated( attenuation ) ) deallocate(attenuation)
    if( allocated( rrange      ) ) deallocate(rrange)

    obs%meta(1) = lon0
    obs%meta(2) = lat0
    obs%meta(3) = z0

    obs%nobs = 0
    obs_ref%nobs = 0

    do idx = 1, nobs_sp
      ! Thinning
      ii = nint( ( grid_lon_ze(idx) - lon(1) ) / dlon ) + 1
      jj = nint( ( grid_lat_ze(idx) - lat(1) ) / dlat ) + 1
      kk = nint( ( grid_z_ze(idx) - z(1) ) / RADAR_SO_SIZE_VERT ) + 1

      if ( mod(ii, RADAR_THIN_HORI) /= 0 .or. mod(jj, RADAR_THIN_HORI) /= 0 .or. &
           mod(kk, RADAR_THIN_VERT) /= 0 ) cycle

      if (grid_count_ze(idx) > 0) then
        obs%nobs = obs%nobs + 1

        ! Count refrectivity obs ( > MIN_RADAR_REF ) below RADAR_ZMAX
        if ( grid_ze(idx) > MIN_RADAR_REF .and. grid_z_ze(idx) < RADAR_ZMAX ) then
          obs_ref%nobs = obs_ref%nobs + 1
        end if
      end if
      if (grid_count_vr(idx) > 0) then
        obs%nobs = obs%nobs + 1
      end if
    end do

    call obs_info_allocate( obs,     extended=.true. )
    call obs_info_allocate( obs_ref, extended=.true. )

    n = 0
    n_ref = 0
    nobs_ze = 0
    nobs_vr = 0
    min_obs_ze = huge(1.0d0)
    max_obs_ze = -huge(1.0d0)
    min_obs_vr = huge(1.0d0)
    max_obs_vr = -huge(1.0d0)

    do idx = 1, nobs_sp
      ii = nint( ( grid_lon_ze(idx) - lon(1) ) / dlon ) + 1
      jj = nint( ( grid_lat_ze(idx) - lat(1) ) / dlat ) + 1
      kk = nint( ( grid_z_ze(idx) - z(1) ) / RADAR_SO_SIZE_VERT ) + 1

      ! Thinning
      if ( mod(ii, RADAR_THIN_HORI) /= 0 .or. mod(jj, RADAR_THIN_HORI) /= 0 .or. &
           mod(kk, RADAR_THIN_VERT) /= 0 ) cycle

      if (grid_count_ze(idx) > 0) then
        n = n + 1
        obs%elm(n) = id_radar_ref_obs
        obs%lon(n) = grid_lon_ze(idx)
        obs%lat(n) = grid_lat_ze(idx)
        obs%lev(n) = grid_z_ze(idx)
        obs%dat(n) = grid_ze(idx)
        ! Add RADAR_BIAS_CONST_DBZ in dBZ
        if ( RADAR_BIAS_COR_RAIN .and. grid_ze(idx) > MIN_RADAR_REF ) then 
          obs%dat(n) = grid_ze(idx) * RADAR_BIAS_RAIN_CONST
        elseif ( RADAR_BIAS_COR_CLR .and. grid_ze(idx) < MIN_RADAR_REF )  then
          obs%dat(n) = grid_ze(idx) * RADAR_BIAS_CLR_CONST
        endif
        obs%err(n) = OBSERR_RADAR_REF
        obs%typ(n) = 22
        obs%dif(n) = 0.0d0
        nobs_ze = nobs_ze + 1
        if (grid_ze(idx) > max_obs_ze) max_obs_ze = grid_ze(idx)
        if (grid_ze(idx) < min_obs_ze) min_obs_ze = grid_ze(idx)

        if ( grid_ze(idx) > MIN_RADAR_REF .and. grid_z_ze(idx) < RADAR_ZMAX ) then
          n_ref = n_ref + 1
          obs_ref%elm(n_ref) = id_radar_ref_obs
          obs_ref%lon(n_ref) = grid_lon_ze(idx)
          obs_ref%lat(n_ref) = grid_lat_ze(idx)
          obs_ref%lev(n_ref) = grid_z_ze(idx)
          if ( RADAR_BIAS_COR_RAIN ) then
            ! Add RADAR_BIAS_CONST_DBZ in dBZ
            obs_ref%dat(n_ref) = grid_ze(idx) * RADAR_BIAS_RAIN_CONST
          else
            obs_ref%dat(n_ref) = grid_ze(idx)
          end if
        end if
      end if

      if (grid_count_vr(idx) > 0) then
        n = n + 1
        obs%elm(n) = id_radar_vr_obs
        obs%lon(n) = grid_lon_ze(idx)
        obs%lat(n) = grid_lat_ze(idx)
        obs%lev(n) = grid_z_ze(idx)
        obs%dat(n) = grid_vr(idx)
        obs%err(n) = OBSERR_RADAR_VR
        obs%typ(n) = 22
        obs%dif(n) = 0.0d0
        nobs_vr = nobs_vr + 1
        if (grid_vr(idx) > max_obs_vr) max_obs_vr = grid_vr(idx)
        if (grid_vr(idx) < min_obs_vr) min_obs_vr = grid_vr(idx)
      end if
    end do

    if( allocated(radlon) ) deallocate(radlon)
    if( allocated(radlat) ) deallocate(radlat)
    if( allocated(radz  ) ) deallocate(radz)
    deallocate(az, el)

    call obs_info_deallocate( obs_ref )

    if( allocated( grid_index    ) ) deallocate(grid_index)
    if( allocated( grid_ze       ) ) deallocate(grid_ze)
    if( allocated( grid_lon_ze   ) ) deallocate(grid_lon_ze)
    if( allocated( grid_lat_ze   ) ) deallocate(grid_lat_ze)
    if( allocated( grid_z_ze     ) ) deallocate(grid_z_ze)
    if( allocated( grid_count_ze ) ) deallocate(grid_count_ze)
    if( allocated( grid_vr       ) ) deallocate(grid_vr)
    if( allocated( grid_lon_vr   ) ) deallocate(grid_lon_vr)
    if( allocated( grid_lat_vr   ) ) deallocate(grid_lat_vr)
    if( allocated( grid_z_vr     ) ) deallocate(grid_z_vr)
    if( allocated( grid_count_vr ) ) deallocate(grid_count_vr)

    return
  end subroutine read_obs_radar_toshiba_pawr

  subroutine read_obs_radar_toshiba_mp_pawr( obs, cfile, maskfile )
    use iso_c_binding
    use scale_da_read_mp_pawr_toshiba, only: &
      read_toshiba_mpr => DA_read_mp_pawr_toshiba, &
      c_mppawr_header, &
      RDIM,            &
      AZDIM,           &
      ELDIM
    implicit none

    character(len=*), intent(in) :: cfile
    character(len=*), intent(in) :: maskfile
    type(obs_info), intent(out) :: obs

    type(obs_info) :: obs_ref

    integer, parameter :: n_type = 2
    character(len=4), parameter :: file_type_sfx(n_type) = (/'.ze', '.vr'/)
    logical, parameter :: input_is_dbz = .true.
    integer, parameter :: opt_verbose = 0 !!! for MP-PAWR toshiba format    

    integer :: access !FILE INQUIRY
    integer :: ios
    integer(4), save :: shadow_na, shadow_ne
    integer(4), allocatable, save :: tmpshadow(:)
    integer(2), allocatable, save :: shadow(:,:)
    real(8),save :: shadow_del_az

    type(c_mppawr_header) :: hd(n_type)
    real(kind=c_float), allocatable, save :: rtdat(:, :, :, :)
    real(kind=c_float), allocatable, save :: az(:, :, :)
    real(kind=c_float), allocatable, save :: el(:, :, :)
    integer :: j, ierr, ierr2
    character(len=3) :: fname
    integer, save::i=0

    real(RP), allocatable :: ze(:, :, :), vr(:, :, :), qcflag(:, :, :), attenuation(:, :, :), rrange(:)
    real(RP), allocatable :: radlon(:, :, :), radlat(:, :, :), radz(:, :, :)
    real(RP), allocatable :: lon(:), lat(:), z(:)
    integer(8), allocatable :: grid_index(:), grid_count_ze(:), grid_count_vr(:)
    real(RP), allocatable :: grid_ze(:), grid_vr(:)
    real(RP), allocatable :: grid_lon_ze(:), grid_lat_ze(:), grid_z_ze(:)
    real(RP), allocatable :: grid_lon_vr(:),  grid_lat_vr(:),  grid_z_vr(:)

    character(len=1024) :: input_fname(n_type)
    integer ia, ir, ie
    real(RP) :: dlon, dlat
    integer(8) nobs_sp

    integer,save :: na, nr, ne
    real(RP),save :: lon0, lat0, z0
    real(RP),save :: missing
    integer,save :: range_res

    real(RP) :: max_obs_ze , min_obs_ze , max_obs_vr , min_obs_vr 
    integer :: nobs_ze, nobs_vr
    integer(8) :: idx, n, n_ref
    integer :: pos
    integer, parameter :: int1 = selected_int_kind(1) !1-BYTE INT
    integer(kind = int1) :: tmp_qcf, valid_qcf

    integer,parameter :: qcf_mask(8)=(/ 0, 1, 1, 1, 1, 0, 0, 0 /) !! valid, shadow, clutter possible, clutter certain, interference, range sidelobe /

    integer::qcf_count(0:255)

    character(len=90) :: plotname
    character(len=19) :: timelabel

    integer :: ii, jj, kk

    real(SP), allocatable :: ref3d(:,:,:)
    real(SP), allocatable :: vr3d(:,:,:)
    character(len=255) :: filename
    integer :: irec, iunit, iolen
    integer :: k

    integer ::  nlon, nlat, nlev

    RADAR_SO_SIZE_HORI = max( real( DX, kind=RP ), RADAR_SO_SIZE_HORI )
    RADAR_SO_SIZE_HORI = max( real( DY, kind=RP ), RADAR_SO_SIZE_HORI )

    if( .not. allocated(rtdat) ) allocate(rtdat(RDIM, AZDIM, ELDIM, n_type))
    if( .not. allocated(az)    ) allocate(az(AZDIM, ELDIM, n_type))
    if( .not. allocated(el)    ) allocate(el(AZDIM, ELDIM, n_type))

    if( trim(maskfile) /= '' .and. .NOT. allocated( shadow ) ) then
      iunit = IO_get_available_fid()
      open(iunit, file=trim(maskfile), status="old", access="stream", form="unformatted")
      read(iunit,iostat=ios) shadow_na, shadow_ne
      allocate( shadow( shadow_na, shadow_ne ) )
      if ( ios == 0 )then
        read(iunit,iostat=ios) shadow
        close(iunit)
      else
        write(6,'(3A)') 'file ',trim(maskfile) ,' not found or unsupported format.'
        stop 1
      end if 
    end if

    do j = 1, n_type
      input_fname(j) = trim(cfile)
      call str_replace(input_fname(j), '.<type>', trim(file_type_sfx(j)), pos)
      ierr = read_toshiba_mpr( input_fname(j), opt_verbose, hd(j), az(:,:,j), el(:,:,j), rtdat(:,:,:,j) )
      if( ierr /= 0 ) then
        LOG_INFO("LETKF_obs_readfile",*) 'Warning: File (',trim(input_fname(j)),') cannot be read. Skip.'
        obs%nobs = 0
        return
      endif
    end do

    ! Set obs information
    lon0 = hd(1)%longitude
    lat0 = hd(1)%latitude
    z0   = hd(1)%altitude

    missing   = real( hd(1)%mesh_offset, kind=RP )
    range_res = hd(1)%range_res

    nr = hd(1)%range_num
    na = hd(1)%ray_num
    ne = hd(1)%el_num

    call jst2utc( hd(1)%s_yr, hd(1)%s_mn, hd(1)%s_dy, hd(1)%s_hr, hd(1)%s_mi, hd(1)%s_sc, 0.0_DP, utime_obs )

    allocate( ze         (na, nr, ne) )
    allocate( vr         (na, nr, ne) )
    allocate( qcflag     (na, nr, ne) )
    allocate( attenuation(na, nr, ne) )

    shadow_del_az = 360.0_RP / shadow_na
    if ( trim(maskfile) /= '' .and. .not. allocated(tmpshadow) ) then
       allocate(tmpshadow(na))
    end if

    do ie = 1, ne
       do ir = 1, nr
       do ia = 1, na
          ze(ia,ir,ie) = rtdat(ir,ia,ie,1)
          vr(ia,ir,ie) = rtdat(ir,ia,ie,2)                                                                 

          if( vr(ia,ir,ie) > RADAR_MAX_ABS_VR .or. vr(ia,ir,ie) < -RADAR_MAX_ABS_VR ) vr(ia,ir,ie) = missing
       end do
       end do

       do ir = 1, nr
       do ia = 1, na
          qcflag(ia,ir,ie) = 0.0_RP !valid
       end do
       end do
       if( trim(maskfile) /= '' .and. allocated(shadow) ) then
          if( ie <= shadow_ne ) then
             do ia = 1, na
                tmpshadow(ia) = shadow( min(shadow_na,nint(az(ia, ie, 1) / shadow_del_az) + 1), ie )
             end do
             do ir = 1, nr
             do ia = 1, na
                if( tmpshadow(ia) /= 0 .and. ir >= tmpshadow(ia) ) then
                   qcflag(ia, ir, ie) = 1000.0_RP  !invalid
                end if
             end do
             end do
          end if
       end if

       do ir = 1, nr
       do ia = 1, na
          attenuation(ia, ir, ie) = 1.0_RP !not implemented yet
       end do
       end do
    end do
    deallocate(rtdat)

    allocate(rrange(nr))
    do ir = 1, nr
       rrange(ir) = (dble(ir) - 0.5_RP) * range_res
    end do

    allocate( radlon(na, nr, ne) )
    allocate( radlat(na, nr, ne) )
    allocate( radz  (na, nr, ne) )

    call radar_georeference( lon0, lat0, z0, na, nr, ne,   & ! input
                             real( az(:, 1, 1), kind=RP ), & ! input (assume ordinary scan strategy)
                             rrange,                       & ! input (assume ordinary scan strategy)
                             real( el(1, :, 1), kind=RP ), & ! input (assume ordinary scan strategy)
                             radlon, radlat, radz          ) ! output  

    call define_grid( lon0, lat0, nr, rrange, rrange(nr), RADAR_ZMAX,             & ! input
                      RADAR_SO_SIZE_HORI, RADAR_SO_SIZE_HORI, RADAR_SO_SIZE_VERT, & ! input
                      dlon, dlat, nlon, nlat, nlev, lon, lat, z                   ) ! output

    call radar_superobing( na, nr, ne, radlon, radlat, radz, ze, vr,                      & ! input spherical
                           qcflag, attenuation,                                           & ! input spherical
                           nlon, nlat, nlev, lon, lat, z, dlon, dlat, RADAR_SO_SIZE_VERT, & ! input cartesian
                           missing, input_is_dbz,                                         & ! input param
                           lon0, lat0,                                                    &
                           nobs_sp, grid_index,                                           & ! output array info
                           grid_ze, grid_lon_ze, grid_lat_ze, grid_z_ze, grid_count_ze,   & ! output ze
                           grid_vr, grid_lon_vr, grid_lat_vr, grid_z_vr, grid_count_vr    ) ! output vr

    if( allocated( ze          ) ) deallocate(ze)
    if( allocated( vr          ) ) deallocate(vr)
    if( allocated( qcflag      ) ) deallocate(qcflag)
    if( allocated( attenuation ) ) deallocate(attenuation)
    if( allocated( rrange      ) ) deallocate(rrange)

    obs%meta(1) = lon0
    obs%meta(2) = lat0
    obs%meta(3) = z0

    obs%nobs = 0
    obs_ref%nobs = 0

    do idx = 1, nobs_sp
      ! Thinning
      ii = nint( ( grid_lon_ze(idx) - lon(1) ) / dlon ) + 1
      jj = nint( ( grid_lat_ze(idx) - lat(1) ) / dlat ) + 1
      kk = nint( ( grid_z_ze(idx) - z(1) ) / RADAR_SO_SIZE_VERT ) + 1

      if ( mod(ii, RADAR_THIN_HORI) /= 0 .or. mod(jj, RADAR_THIN_HORI) /= 0 .or. &
           mod(kk, RADAR_THIN_VERT) /= 0 ) cycle

      if (grid_count_ze(idx) > 0) then
        obs%nobs = obs%nobs + 1

        ! Count refrectivity obs ( > MIN_RADAR_REF ) below RADAR_ZMAX
        if ( grid_ze(idx) > MIN_RADAR_REF .and. grid_z_ze(idx) < RADAR_ZMAX ) then
          obs_ref%nobs = obs_ref%nobs + 1
        end if
      end if
      if (grid_count_vr(idx) > 0) then
        obs%nobs = obs%nobs + 1
      end if
    end do

    call obs_info_allocate( obs,     extended=.true. )
    call obs_info_allocate( obs_ref, extended=.true. )

    n = 0
    n_ref = 0
    nobs_ze = 0
    nobs_vr = 0
    min_obs_ze = huge(1.0d0)
    max_obs_ze = -huge(1.0d0)
    min_obs_vr = huge(1.0d0)
    max_obs_vr = -huge(1.0d0)

    do idx = 1, nobs_sp
      ii = nint( ( grid_lon_ze(idx) - lon(1) ) / dlon ) + 1
      jj = nint( ( grid_lat_ze(idx) - lat(1) ) / dlat ) + 1
      kk = nint( ( grid_z_ze(idx) - z(1) ) / RADAR_SO_SIZE_VERT ) + 1

      ! Thinning
      if ( mod(ii, RADAR_THIN_HORI) /= 0 .or. mod(jj, RADAR_THIN_HORI) /= 0 .or. &
           mod(kk, RADAR_THIN_VERT) /= 0 ) cycle

      if (grid_count_ze(idx) > 0) then
        n = n + 1
        obs%elm(n) = id_radar_ref_obs
        obs%lon(n) = grid_lon_ze(idx)
        obs%lat(n) = grid_lat_ze(idx)
        obs%lev(n) = grid_z_ze(idx)
        obs%dat(n) = grid_ze(idx)
        ! Add RADAR_BIAS_CONST_DBZ in dBZ
        if ( RADAR_BIAS_COR_RAIN .and. grid_ze(idx) > MIN_RADAR_REF ) then 
          obs%dat(n) = grid_ze(idx) * RADAR_BIAS_RAIN_CONST
        elseif ( RADAR_BIAS_COR_CLR .and. grid_ze(idx) < MIN_RADAR_REF )  then
          obs%dat(n) = grid_ze(idx) * RADAR_BIAS_CLR_CONST
        endif
        obs%err(n) = OBSERR_RADAR_REF
        obs%typ(n) = 22
        obs%dif(n) = 0.0d0
        nobs_ze = nobs_ze + 1
        if (grid_ze(idx) > max_obs_ze) max_obs_ze = grid_ze(idx)
        if (grid_ze(idx) < min_obs_ze) min_obs_ze = grid_ze(idx)

        if ( grid_ze(idx) > MIN_RADAR_REF .and. grid_z_ze(idx) < RADAR_ZMAX ) then
          n_ref = n_ref + 1
          obs_ref%elm(n_ref) = id_radar_ref_obs
          obs_ref%lon(n_ref) = grid_lon_ze(idx)
          obs_ref%lat(n_ref) = grid_lat_ze(idx)
          obs_ref%lev(n_ref) = grid_z_ze(idx)
          if ( RADAR_BIAS_COR_RAIN ) then
            ! Add RADAR_BIAS_CONST_DBZ in dBZ
            obs_ref%dat(n_ref) = grid_ze(idx) * RADAR_BIAS_RAIN_CONST
          else
            obs_ref%dat(n_ref) = grid_ze(idx)
          end if
        end if
      end if

      if (grid_count_vr(idx) > 0) then
        n = n + 1
        obs%elm(n) = id_radar_vr_obs
        obs%lon(n) = grid_lon_ze(idx)
        obs%lat(n) = grid_lat_ze(idx)
        obs%lev(n) = grid_z_ze(idx)
        obs%dat(n) = grid_vr(idx)
        obs%err(n) = OBSERR_RADAR_VR
        obs%typ(n) = 22
        obs%dif(n) = 0.0d0
        nobs_vr = nobs_vr + 1
        if (grid_vr(idx) > max_obs_vr) max_obs_vr = grid_vr(idx)
        if (grid_vr(idx) < min_obs_vr) min_obs_vr = grid_vr(idx)
      end if
    end do

    if( allocated(radlon) ) deallocate(radlon)
    if( allocated(radlat) ) deallocate(radlat)
    if( allocated(radz  ) ) deallocate(radz)
    deallocate(az, el)

    call obs_info_deallocate( obs_ref )

    if( allocated( grid_index    ) ) deallocate(grid_index)
    if( allocated( grid_ze       ) ) deallocate(grid_ze)
    if( allocated( grid_lon_ze   ) ) deallocate(grid_lon_ze)
    if( allocated( grid_lat_ze   ) ) deallocate(grid_lat_ze)
    if( allocated( grid_z_ze     ) ) deallocate(grid_z_ze)
    if( allocated( grid_count_ze ) ) deallocate(grid_count_ze)
    if( allocated( grid_vr       ) ) deallocate(grid_vr)
    if( allocated( grid_lon_vr   ) ) deallocate(grid_lon_vr)
    if( allocated( grid_lat_vr   ) ) deallocate(grid_lat_vr)
    if( allocated( grid_z_vr     ) ) deallocate(grid_z_vr)
    if( allocated( grid_count_vr ) ) deallocate(grid_count_vr)

    return
  end subroutine read_obs_radar_toshiba_mp_pawr

  !-------------------------------------------------------------------------------
  ! Replace the first occurrence of 'oldsub' in 'str' with 'newsub';
  ! note that 'str' will be left-adjusted no matter whether 'oldsub' is found
  !-------------------------------------------------------------------------------
  ! [INPUT]
  !   str    : input string
  !   oldsub : old substring to be replaced
  !   newsub : new substring
  ! [OUTPUT]
  !   str    : output string with substring replaced
  !   pos    : the start position of the replaced substring; if not found, return 0
  !-------------------------------------------------------------------------------
  subroutine str_replace(str, oldsub, newsub, pos)
    use scale_prc, only: &
      PRC_abort
    implicit none
    character(len=*), intent(inout) :: str
    character(len=*), intent(in) :: oldsub
    character(len=*), intent(in) :: newsub
    integer, intent(out) :: pos
    integer :: str_lent, oldsub_len, newsub_len, shift

    str = adjustl(str)
    str_lent = len_trim(str)
    oldsub_len = len(oldsub)
    newsub_len = len(newsub)

    pos = index(str, oldsub)
    if (pos >= 1) then
      shift = newsub_len - oldsub_len
      if (shift > 0) then
        if (str_lent+shift > len(str)) then
          LOG_ERROR('str_replace', *) "The length of 'str' string is not enough for substitution."
          call PRC_abort
        end if
        str(pos+oldsub_len:str_lent+shift) = adjustr(str(pos+oldsub_len:str_lent+shift))
      else if (shift < 0) then
        str(pos+newsub_len:pos+oldsub_len-1) = repeat(' ', 0-shift)
        str(pos+newsub_len:str_lent) = adjustl(str(pos+newsub_len:str_lent))
      end if
      str(pos:pos+newsub_len-1) = newsub
    end if

    return
  end subroutine str_replace

  subroutine jst2utc(jyear, jmonth, jday, jhour, jminute, jsecond, jtime_ms, utime)
    use scale_calendar, only: &
        CALENDAR_date2daysec, &
        CALENDAR_daysec2date, &
        CALENDAR_adjust_daysec
    implicit none

    integer, intent(in) :: jyear, jmonth, jday
    integer, intent(in) :: jhour, jminute, jsecond
    real(DP) :: jtime_ms
    integer, intent(out) :: utime(6)
    integer :: jtime(6)
    integer :: absday
    real(DP) :: abssec, utime_ms

    jtime(1) = jyear
    jtime(2) = jmonth
    jtime(3) = jday
    jtime(4) = jhour
    jtime(5) = jminute
    jtime(6) = jsecond

    call CALENDAR_date2daysec( absday,       & ! [OUT]
                               abssec,       & ! [OUT]
                               jtime,        & ! [IN]
                               jtime_ms,     & ! [IN]
                               0             ) ! [IN]

    abssec = abssec - real(3600*9, kind=DP)

    call CALENDAR_adjust_daysec( absday,   & ! [INOUT]
                                 abssec )    ! [INOUT]

    call CALENDAR_daysec2date( utime,   & ! [OUT]
                               utime_ms, & ! [OUT]
                               absday,      & ! [IN]
                               abssec,      & ! [IN]
                               0            ) ! [IN]

    return
  end subroutine jst2utc

  !=======================================================================
  !
  ! [PURPOSE:] Thinning of radar data
  !
  ! [HISTORY:] This version produce a superobing of the observations but
  ! the data is stored in azimuth , range , elevation. Conversion to the 
  ! lat , lon , z is performed by the observation operator.
  !
  ! Modified to produce 1D list of superobservations with lon, lat, z
  !
  !=======================================================================
  subroutine radar_georeference(lon0, lat0, z0, na, nr, ne, azimuth, rrange, elevation, radlon, radlat, radz)
    use scale_const, only: &
      RADIUS => CONST_RADIUS, &
      D2R    => CONST_D2R,    &
      R2D    => CONST_R2D
    implicit none

    real(RP), parameter :: ke = 4.0_RP / 3.0_RP ! factor for the effective radius of the earth

    real(RP) :: ke_Re ! effective radius of the earth [m]

    real(RP), intent(in) :: lon0, lat0, z0
    integer,  intent(in) :: na, nr, ne
    real(RP), intent(in) :: azimuth(na), rrange(nr), elevation(ne)
    real(RP), intent(out) :: radlon(na, nr, ne), radlat(na, nr, ne)
    real(RP), intent(out) :: radz(na, nr, ne)

    real(RP) sin_elev_ke_Re_2, cos_elev_div_ke_Re, tmpdist
    real(RP) cdist, sdist, sinll1, cosll1, sinll1_cdist, cosll1_sdist, cosll1_cdist, sinll1_sdist
    real(RP) :: azimuth_rad, sin_azim(na), cos_azim(na)
    integer :: ia, ir, ie

    ke_Re = ke * RADIUS 

    ! This code is copied from juan ruiz's qc code and modified
    sinll1 = sin(lat0 * D2R)
    cosll1 = cos(lat0 * D2R)
    do ia = 1, na
       azimuth_rad = azimuth(ia) * D2R
       sin_azim(ia) = sin(azimuth_rad)
       cos_azim(ia) = cos(azimuth_rad)
    end do

    do ie = 1, ne
       sin_elev_ke_Re_2 = sin(elevation(ie) * D2R) * ke_Re * 2
       cos_elev_div_ke_Re = cos(elevation(ie) * D2R) / ke_Re
       do ir = 1, nr
          ! Perform standard height beam heigth computation.
          radz(:, ir, ie) = z0 + sqrt(rrange(ir) * (rrange(ir) + sin_elev_ke_Re_2) + ke_Re * ke_Re) - ke_Re
          tmpdist = ke * asin(rrange(ir) * cos_elev_div_ke_Re)
          if (tmpdist .eq. 0d0) then
             radlon(1:na, ir, ie) = lon0
             radlat(1:na, ir, ie) = lat0
          else
             cdist = cos(tmpdist)
             sdist = sin(tmpdist)
             sinll1_cdist = sinll1 * cdist
             cosll1_sdist = cosll1 * sdist
             cosll1_cdist = cosll1 * cdist
             sinll1_sdist = sinll1 * sdist
             do ia = 1, na
                radlat(ia, ir, ie) = asin(sinll1_cdist + cosll1_sdist * cos_azim(ia)) * R2D
                radlon(ia, ir, ie) = lon0 + atan2(sdist * sin_azim(ia), cosll1_cdist - sinll1_sdist * cos_azim(ia)) * R2D
             end do
          end if
       end do
    end do

    return
  end subroutine radar_georeference

  !-----------------------------------------------------------------------
  ! Main superobing routine
  !-----------------------------------------------------------------------
  subroutine radar_superobing(na, nr, ne, radlon, radlat, radz, ze, vr, & ! input spherical
       &                       qcflag, attenuation, & ! input spherical 2
       &                       nlon, nlat, nlev, lon, lat, z, dlon, dlat, dz, & ! input cartesian
       &                       missing, input_is_dbz, & ! input param
       &                       lon0, lat0, & ! input param
       &                       nobs_sp, grid_index, & ! output array info
       &                       grid_ref, grid_lon_ref, grid_lat_ref, grid_z_ref, grid_count_ref, &  ! output ze
       &                       grid_vr, grid_lon_vr, grid_lat_vr, grid_z_vr, grid_count_vr )       ! output vr
    use scale_sort, only: &
       SORT_uniq_int_sorted
    implicit none

    integer, intent(in) :: na, nr, ne ! array size of the spherical grid
    real(RP), intent(in), dimension(na, nr, ne) :: radlon, radlat, radz ! georeference
    real(RP), intent(in) :: ze(na, nr, ne), vr(na, nr, ne) ! main data
    real(RP), intent(in) :: qcflag(na, nr, ne), attenuation(na, nr, ne) ! additional info
    integer, intent(in) :: nlon, nlat, nlev ! array size of the cartesian grid
    real(RP), intent(in) :: lon(nlon), lat(nlat), z(nlev)
    real(RP), intent(in) :: dlon, dlat, dz, missing
    logical, intent(in) :: input_is_dbz
    real(RP), intent(in) :: lon0, lat0

    integer(8), intent(out) :: nobs_sp
    integer(8), allocatable, intent(out) :: grid_index(:), grid_count_ref(:), grid_count_vr(:)
    REAL(RP), allocatable, intent(out) :: grid_ref(:), grid_vr(:)
    REAL(RP), allocatable, intent(out) :: grid_lon_ref(:), grid_lat_ref(:), grid_z_ref(:) !Lat, Lon and Z weighted by the observations.
    REAL(RP), allocatable, intent(out) :: grid_lon_vr(:),  grid_lat_vr(:), grid_z_vr(:)  !Lat, Lon and Z weighted by the observations.

    REAL(RP), ALLOCATABLE :: grid_w_vr(:)
    REAL(RP), ALLOCATABLE :: grid_meanvr(:), grid_stdvr(:) !Non weighted radial velocity and its dispersion within each box.

    integer(8), allocatable :: packed_grid(:), pack2uniq(:), nobs_each_elev(:)
    real(RP), allocatable :: packed_data(:, :)
    logical, allocatable :: packed_attn(:)

    REAL(RP) :: count_inv, tmpstdvr, tmpmeanvr
    integer(8) :: idx, jdx, nobs, sidx(ne), eidx(ne)

    integer :: i
    !integer i, e0, e1, ne_mpi
    !integer, allocatable :: j_mpi(:, :), sendcount(:), recvcount(:), recvoffset(:)
    !integer(8), allocatable :: sendbuf(:), nobs_each_mpi(:)
    !integer(8), allocatable :: packed_grid_mpi(:)
    !real(RP), allocatable :: packed_data_mpi(:, :)
    !logical, allocatable :: packed_attn_mpi(:)

    integer :: ierr

    !allocate(sendcount(0:(mpiprocs - 1)), recvcount(0:(mpiprocs - 1)), recvoffset(0:(mpiprocs - 1)))
    !allocate(j_mpi(2, 0:(mpiprocs - 1)))
    !call set_mpi_div(j_mpi, mpiprocs, int(ne, 8))
    !e0 = j_mpi(1, myrank)
    !e1 = j_mpi(2, myrank)
    !ne_mpi = e1 - e0 + 1

    ! AVERAGE DATA AND INCLUDE OBSERVATIONA ERROR.
    ! We will compute the i,j,k for each radar grid point and box average the
    ! data.

    ! QC, Indexing, and packing simultaneously
    allocate(nobs_each_elev(ne))
    !if(mpiprocs > 1) then
    !   allocate(packed_grid_mpi(na * nr * ne_mpi), &
    !        &   packed_data_mpi(5, na * nr * ne_mpi), &
    !        &   packed_attn_mpi(na * nr * ne_mpi), &
    !        &   nobs_each_mpi(e0:e1))
    !   call qc_indexing_and_packing( &
    !        & na, nr, ne_mpi, ze(:, :, e0:e1), vr(:, :, e0:e1), & ! input spherical
    !        & radlon(:, :, e0:e1), radlat(:, :, e0:e1), radz(:, :, e0:e1), & ! input spherical
    !        & qcflag(:, :, e0:e1), input_is_dbz, attenuation(:, :, e0:e1), & ! input spherical
    !        & nlon, nlat, nlev, lon, lat, z, dlon, dlat, dz, & ! input cartesian
    !        & missing, & ! input param
    !        & lon0, lat0, &
    !        & nobs_each_mpi, packed_grid_mpi, packed_data_mpi, packed_attn_mpi) ! output

    !   sendcount(:) = ne_mpi
    !   recvcount(0:(mpiprocs - 1)) = j_mpi(2, 0:(mpiprocs - 1)) - j_mpi(1, 0:(mpiprocs - 1)) + 1
    !   recvoffset = j_mpi(1, :) - 1 !START FROM 0
    !   call MPI_Allgatherv(nobs_each_mpi, sendcount(0), MPI_INTEGER8, &
    !        &              nobs_each_elev, recvcount, recvoffset, MPI_INTEGER8, &
    !        &              comm, ierr)
    !   deallocate(nobs_each_mpi)
    !else
       allocate(packed_grid(na * nr * ne), &
            &   packed_data(5, na * nr * ne), &
            &   packed_attn(na * nr * ne))
       call qc_indexing_and_packing( &
            & na, nr, ne, ze, vr, radlon, radlat, radz, & ! input spherical
            & qcflag, input_is_dbz, attenuation, & ! input spherical
            & nlon, nlat, nlev, lon, lat, z, dlon, dlat, dz, & ! input cartesian
            & missing, & ! input param
            & lon0, lat0, &
            & nobs_each_elev, packed_grid, packed_data, packed_attn) ! output
    !end if
    nobs = sum(nobs_each_elev)
    sidx(1) = 1
    do i = 1, ne - 1
       sidx(i + 1) = sidx(i) + nobs_each_elev(i)
    end do !i
    eidx = sidx + nobs_each_elev - 1
    deallocate(nobs_each_elev)
    !! MPI packed data
    !if(mpiprocs > 1) then
    !   sendcount = eidx(e1) - sidx(e0) + 1
    !   do i = 0, mpiprocs - 1
    !      recvoffset(i) = sidx(j_mpi(1, i)) - 1 !START FROM 0
    !      recvcount(i) = eidx(j_mpi(2, i)) - sidx(j_mpi(1, i)) + 1
    !   end do
    !   allocate(packed_grid(nobs), packed_data(5, nobs), packed_attn(nobs))
    !   ! NEEDED BY ALL
    !   call MPI_Allgatherv(packed_grid_mpi, sendcount(0), MPI_INTEGER8, &
    !        &              packed_grid, recvcount, recvoffset, MPI_INTEGER8, &
    !        &              comm, ierr)
    !   ! ONLY NEEDED BY ROOT
    !   call MPI_Gatherv(packed_attn_mpi, sendcount(0), MPI_LOGICAL, &
    !        &           packed_attn, recvcount, recvoffset, MPI_LOGICAL, &
    !        &           0, comm, ierr)
    !   call MPI_Gatherv(packed_data_mpi, sendcount(0) * 5, datatype, &
    !        &           packed_data, recvcount * 5, recvoffset * 5, datatype, &
    !        &           0, comm, ierr)
    !   deallocate(packed_grid_mpi, packed_data_mpi, packed_attn_mpi)

    !end if

    ! Sort index array
    allocate(grid_index(nobs))
    do idx = 1, nobs
       grid_index(idx) = packed_grid(idx)
    end do
    !if(mpiprocs > 1) then
    !   call merge_sort_mpi(nobs, grid_index, comm) !ONLY RANK 0 RETURNS DATA
    !   call MPI_Bcast(grid_index, int(nobs), MPI_INTEGER8, 0, comm, ierr)
    !else
       call merge_sort_parallel(nobs, grid_index)
    !end if

    ! Unique superobs (nobs_sp)
    call SORT_uniq_int_sorted(nobs, grid_index, nobs_sp)

    ! Inverted indexing
    allocate(pack2uniq(nobs))
    !call set_mpi_div(j_mpi, mpiprocs, nobs)
    !allocate(sendbuf(j_mpi(2, myrank) - j_mpi(1, myrank) + 1))
    !do idx = j_mpi(1, myrank), j_mpi(2, myrank)
    !   sendbuf(idx - j_mpi(1, myrank) + 1) = binary_search_i8(nobs_sp, grid_index, packed_grid(idx))
    do idx = 1, nobs
       pack2uniq(idx) = binary_search_i8(nobs_sp, grid_index, packed_grid(idx))
    end do
    !sendcount(:) = j_mpi(2, myrank) - j_mpi(1, myrank) + 1
    !recvcount(0:(mpiprocs - 1)) = j_mpi(2, 0:(mpiprocs - 1)) - j_mpi(1, 0:(mpiprocs - 1)) + 1
    !recvoffset = j_mpi(1, :) - 1
    !! ONLY NEEDED BY ROOT
    !call MPI_Gatherv(sendbuf, sendcount(0), MPI_INTEGER8, &
    !     &           pack2uniq, recvcount, recvoffset, MPI_INTEGER8, &
    !     &           0, comm, ierr)
    !deallocate(j_mpi, sendbuf, sendcount, recvcount, recvoffset) !END OF MPI

    ! Allocate output arrays
    allocate(grid_ref(nobs_sp), grid_vr(nobs_sp))
    allocate(grid_count_ref(nobs_sp), grid_count_vr(nobs_sp))
    allocate(grid_lon_ref(nobs_sp), grid_lat_ref(nobs_sp), grid_z_ref(nobs_sp))
    allocate(grid_lon_vr(nobs_sp) , grid_lat_vr(nobs_sp) , grid_z_vr(nobs_sp))
    allocate(grid_w_vr(nobs_sp))
    allocate(grid_meanvr(nobs_sp), grid_stdvr(nobs_sp))

    !if(myrank > 0) return !ONLY RANK 0 DOES THE REMAINING WORK

    !Initialize arrays
    do idx = 1, nobs_sp
       grid_count_ref(idx) = 0
       grid_count_vr(idx)  = 0
       grid_ref(idx)       = 0.0d0
       grid_vr(idx)        = 0.0d0
       grid_w_vr(idx)      = 0.0d0
       grid_lon_ref(idx)   = 0.0d0
       grid_lat_ref(idx)   = 0.0d0
       grid_z_ref(idx)     = 0.0d0
       grid_lon_vr(idx)    = 0.0d0
       grid_lat_vr(idx)    = 0.0d0
       grid_z_vr(idx)      = 0.0d0
       grid_meanvr(idx)    = 0.0d0
       grid_stdvr(idx)     = 0.0d0
    end do

    do jdx = 1, nobs
       idx = pack2uniq(jdx)
       ! PROCESS REFLECITIVITY
       ! use attenuation estimates / ignore estimates
       if(packed_attn(jdx)) then
          grid_ref(idx) = grid_ref(idx) + packed_data(1, jdx)
          grid_count_ref(idx) = grid_count_ref(idx) + 1
          grid_lon_ref(idx) = grid_lon_ref(idx) + packed_data(3, jdx)
          grid_lat_ref(idx) = grid_lat_ref(idx) + packed_data(4, jdx)
          grid_z_ref(idx) = grid_z_ref(idx) + packed_data(5, jdx)
       ENDIF

       ! CONSIDER THE RADIAL VELOCITY
       ! Wind will be averaged using an average weighted by returned power.
       ! (this should reduce noise). 
       IF(packed_data(2, jdx) .GT. missing) THEN !PROCESS WIND
          grid_w_vr(idx) = grid_w_vr(idx) + packed_data(1, jdx)
          grid_count_vr(idx) = grid_count_vr(idx) + 1
          grid_meanvr(idx) = grid_meanvr(idx) + packed_data(2, jdx)
          grid_stdvr(idx) = grid_stdvr(idx) + packed_data(2, jdx) ** 2
          grid_vr(idx) = grid_vr(idx) + packed_data(2, jdx) * packed_data(1, jdx)
          grid_lon_vr(idx) = grid_lon_vr(idx) + packed_data(3, jdx) * packed_data(1, jdx)
          grid_lat_vr(idx) = grid_lat_vr(idx) + packed_data(4, jdx) * packed_data(1, jdx)
          grid_z_vr(idx) = grid_z_vr(idx) + packed_data(5, jdx) * packed_data(1, jdx)
       ENDIF
    ENDDO !jdx

    ! Average data and write observation file (FOR LETKF)
    DO idx = 1, nobs_sp
       IF(grid_count_ref(idx) > 0) THEN  !Process reflectivity
          count_inv = 1.0d0 / dble(grid_count_ref(idx))
          grid_ref(idx)     = grid_ref(idx)     * count_inv
          grid_lon_ref(idx) = grid_lon_ref(idx) * count_inv
          grid_lat_ref(idx) = grid_lat_ref(idx) * count_inv
          grid_z_ref(idx)   = grid_z_ref(idx)   * count_inv
       ENDIF

       IF(grid_count_vr(idx) > 0) THEN
          count_inv = 1.0d0 / grid_w_vr(idx)
          grid_vr(idx)     = grid_vr(idx)     * count_inv
          grid_lon_vr(idx) = grid_lon_vr(idx) * count_inv
          grid_lat_vr(idx) = grid_lat_vr(idx) * count_inv
          grid_z_vr(idx)   = grid_z_vr(idx)   * count_inv

          ! If variability within a box is big then we may have:
          ! -small scale strong phenomena (tornado!)
          ! -noise in the data.
          ! In any case averaging the data is not a god idea so this data
          ! can be rejected a priori.
          IF( RADAR_USE_VR_STD ) THEN
             count_inv = 1.0d0 / dble(grid_count_vr(idx))
             tmpmeanvr = grid_meanvr(idx) * count_inv
             tmpstdvr = grid_stdvr(idx)  * count_inv
             tmpstdvr = SQRT(tmpstdvr - (tmpmeanvr ** 2))
             IF(tmpstdvr > vr_std_threshold) grid_count_vr(idx) = 0 !Reject the observation.
          ENDIF
       end IF
    end do

    return
  end SUBROUTINE radar_superobing

  !-----------------------------------------------------------------------
  ! Define grid
  !-----------------------------------------------------------------------
  subroutine define_grid(lon0, lat0, nr, rrange, maxrange, maxz, dx, dy, dz, & ! input
       &                 dlon, dlat, nlon, nlat, nlev, lon, lat, z)            ! output
    use scale_const, only: &
      RADIUS => CONST_RADIUS, &
      D2R    => CONST_D2R,    &
      R2D    => CONST_R2D
    implicit none
    real(RP), intent(in) :: lon0, lat0
    integer, intent(in) :: nr
    real(RP), intent(in) :: rrange(nr)
    real(RP), intent(in) :: maxrange, maxz, dx, dy, dz
    real(RP), intent(out) :: dlon, dlat
    integer, intent(out) :: nlon, nlat, nlev
    real(RP), allocatable, intent(out) :: lon(:), lat(:), z(:)
    real(RP) :: tmpmaxrange

    integer :: i, j, k

    ! Translate DX into an appropiate DLON and DLAT.
    ! Hopefully nobody will put a radar at the pole.
    dlon = R2D * dx / (cos(lat0 * D2R) * RADIUS)
    dlat = R2D * dx / RADIUS

    ! Compute possible value for NLON in order to cover the maximum radar range.
    tmpmaxrange = maxrange
    tmpmaxrange = min(tmpmaxrange, rrange(nr)) 
    tmpmaxrange = 2.0 * tmpmaxrange
    nlon = CEILING(tmpmaxrange / dx)
    nlat = CEILING(tmpmaxrange / dy)
    nlev = CEILING(maxz / dz)

    allocate(lon(nlon), lat(nlat), z(nlev))

    ! Force grid dimensions to be odd
    IF( MOD(nlon,2) == 0 ) nlon = nlon + 1
    IF( MOD(nlat,2) == 0 ) nlat = nlat + 1

    do i = 1, nlon
       lon(i) = lon0 + dlon * (-1.0 - (nlon - 1.0) / 2.0 + i)
    end do

    do j = 1, nlat
       lat(j) = lat0 + dlat * (-1.0 - (nlat - 1.0) / 2.0 + j)
    end DO

    do k = 1, nlev
       z(k) = dz * (k - 1)
    end do

  end subroutine define_grid

  !-----------------------------------------------------------------------
  ! Compute radar reflectivity and radial wind.
  ! Radial wind computations for certain methods depend on model reflectivity
  ! so both functions has been merged into a single one.
  ! First reflectivity is computed, and the the radial velocity is computed.
  !-----------------------------------------------------------------------
  subroutine calc_ref_vr(qv,qc,qr,qci,qs,qg,u,v,w,t,p,az,elev,ref,vr)
    use scale_const, only: &
      PI    => CONST_PI,   &
      GRAV  => CONST_GRAV, &
      D2R   => CONST_D2R,  &
      R2D   => CONST_R2D,  &
      Rdry  => CONST_Rdry
    implicit none
    real(RP), intent(in) :: qv        !Water vapor
    real(RP), intent(in) :: qc,qr     !Cloud and rain water
    real(RP), intent(in) :: qci,qs,qg !Cloud ice, snow and graupel
    real(RP), intent(in) :: t,p       !Temperature and pressure.
    real(RP), intent(in) :: u,v,w     !velocities with respecto to earth.
    real(RP), intent(inout) :: ref      !Reflectivity.
    real(RP)              :: ro
    real(RP), intent(in) :: az     !Azimuth respect to the radar.
    real(RP), intent(in) :: elev   !Elevation angle respect to the surface.
    real(RP), intent(inout) :: vr    !Radial velocity.
    real(RP)  :: qms , qmg !Melting species concentration (METHOD_REF_CALC 3)
    real(RP)  :: qt        !Total condensate mixing ratio (METHOD_REF_CALC 1)
    real(RP)  :: zr , zs , zg !Rain, snow and graupel's reflectivities.
    real(RP)  :: wr , ws , wg !Rain, snow and graupel's mean terminal velocities.
    real(RP)  :: wt           !Total mean terminal velocity.
    real(RP)  :: nor, nos, nog !Rain, snow and graupel's intercepting parameters.
    real(RP)  :: ror, ros, rog , roi !Rain, snow and graupel, ice densities.
    real(RP)  :: a,b,c,d,Cd    !Constant for fall speed computations.
    real(RP)  :: cr_t08,cs_t08,cg_t08,dr_t08,ds_t08,dg_t08    !Constant for fall speed computations (Tomita2008)
    real(RP)  :: cf, pip , roo
    real(RP)  :: ki2 , kr2
    real(RP)  :: lr , ls , lg
    real(RP)  :: tmp_factor , rofactor
    real(RP)  :: p0
    real(RP)  :: Fs, Fg , zms , zmg , fws , fwg !METHOD_REF_CALC 3
    real(RP)  :: qrp , qsp , qgp
    real(RP)  :: maxf                     !Maximum mixture relative concentration. (METHOD_REF_CALC 3)

    real(RP) , parameter :: qeps = 1.0d-20 !Avoid overflow

    real(RP) , parameter :: as_RS14 = 6.9d-2 
    real(RP)  :: tc , MOMs_0bs

    !Note: While equivalent reflectivity is assumed to be independent of the radar, in 
    !practice short wavelengths as those associated with K band radars migh frequently
    !experience Mie scattering. In that case, the equivalent reflectivity is not longer
    !radar independent and an appropiate relationship between the forecasted concentrations
    !and the reflectivity should be used.

    !This model reflectivity won't be lower than this value.

    !Initialize reflectivities
    zr=0.0d0
    zs=0.0d0
    zg=0.0d0
    zms=0.0d0
    zmg=0.0d0
    ref=0.0d0

    !Compute air density (all methods use this)

    ro =  p / ( Rdry * t)


    !Begin computation of reflectivity and vr

    if (METHOD_REF_CALC == 1) then

      !WRF method: See for example Sugimoto et al. Evaluation of the Performance of Ra
      !dial Velocity Assimilation with WRF-3DVAR System and Simulated Multiple-Doppler
      !Radar Data
      !Only rain is used to estimate the terminal velocity of hidrometeors.
      !Only rain is used to compute equivalent reflectivity.
      !Marshall-Palmer distributions are assumed to find the relationship bestween
      !concentration and equivalent reflectivity.
      ! Sun and Crook 1997 , 1998.
      !Derived for C-band radars.
      !Attenuation is not computed.

      !Reflectivity
      nor=8.0d6      ![m^-4]
      ror=1000.0d0   ![Kg/m3]
      pip= PI ** 1.75 !factor
      cf =10.0d18 * 72 !factor
      p0=1.0d5            !Reference pressure.

      qt=qr + qs + qg  !Assume that the total condensate is rain water
                       !But ignore cloud ice and cloud water

      IF( qt .GT. qeps )THEN
      ref = cf * ( ( ro * qt )**1.75 )
      ref = ref / ( pip * ( nor ** 0.75 ) * ( ror ** 1.75 ) )
     !ref= 2.04d4 *( ( ro * qt * 1.0d3 ) ** 1.75 ) !Original Sun and Crook expresion.
      ELSE
      ref=0.0d0
      ENDIF

      !Radial wind

      IF ( qt .GT. qeps )THEN
      a=(p0/p)**0.4
      wt = 5.40d0 * a * ( qt ** 0.125 )
      ELSE
      wt=0d0
      ENDIF


    else if (METHOD_REF_CALC == 2) then
      !Observation operator from Tong and Xue 2006, 2008 a and b.
      !Based on Smith et al 1975.
      !It includes reflectivity contribution by all the microphisical species.
      !is assumes Marshall and Palmer distributions.
      !Based on C band radars.
      !Attenuation is not computed.
      nor=8.0d6      ![m^-4]
      nos=3.0d6      ![m^-4]
      nog=4.0d4      ![m^-4]
      ror=1000.0d0   ![Kg/m3]
      ros=100.0d0    ![Kg/m3]
      rog=913.0d0    ![Kg/m3] 
      roi=917.0d0    ![Kg/m3]
      roo=1.0d0      ![Kg/m3] Surface air density.
      ki2=0.176d0    !Dielectric factor for ice.
      kr2=0.930d0    !Dielectric factor for water.
      pip= PI ** 1.75 !factor
      cf =1.0d18 * 720 !factor 

      IF( qr .GT. qeps )THEN
      zr= cf * ( ( ro * qr )**1.75 )
      zr= zr / ( pip * ( nor ** 0.75 ) * ( ror ** 1.75 ) )
      ENDIF
      !The contribution of snow depends on temperature (bright band effect)
      IF( qs .GT. qeps )THEN
      IF ( t <= 273.16 )THEN
       zs = cf * ki2 * ( ros ** 0.25 ) * ( ( ro * qs ) ** 1.75 )
       zs = zs / ( pip * kr2 * ( nos ** 0.75  ) * ( roi ** 2 ) )
      ELSE
      !WARNING: This formulation has to be checked the paper says that 
       !ros instead of roi should be used in thes expresion, but that 
       !leads to unrealistic jumps in the model derived reflectivity.
       zs = cf * ( ( ro * qs ) ** 1.75 )
       zs = zs / ( pip * ( nos ** 0.75 ) * ( roi ** 1.75 ) )
      ENDIF
      ENDIF

      !Only dry graupel contribution is ussed.
      IF( qg .GT. qeps )THEN
      zg= ( cf / ( pip * ( nog ** 0.75) * ( rog ** 1.75 ) ) ) ** 0.95
      zg= zg * ( ( ro * qg ) ** 1.6625 )
      ENDIF

      ref = zr + zs + zg

      !Compute reflectivity weigthed terminal velocity.
      !Lin et al 1983.
      IF( ref > 0.0d0 )THEN
      !There are hidrometeors, compute their terminal velocity.
      !Change units to be consistent with Lin et al 1983 and
      !to obtain wt in m/s
      nor=nor*1e-3      ![cm^-4]
      nos=nos*1e-3      ![cm^-4]
      nog=nog*1e-3      ![cm^-4]
      ror=ror*1e-3        ![g/cm3]
      ros=ros*1e-3        ![g/cm3]
      rog=rog*1e-3      ![g/cm3] 
      roo=roo*1e-3      ![g/cm3] Surface air density.
      ro= ro*1e-3

      a=2115d0   ![cm**1-b / s]
      b=0.8d0
      c=152.93d0 ![cm**1-b / s]
      d=0.25d0
      Cd=0.6d0

      rofactor= ( roo / ro  ) ** 0.25
      if(qr > qeps )then
        CALL com_gamma( 4.0_RP + b , tmp_factor )
        lr= ( PI * ror * nor / ( ro * qr ) ) ** 0.25
        wr= a * tmp_factor / ( 6.0d0 * ( lr ** b ) )
        wr= 1.0d-2*wr * rofactor
      else
        wr = 0.0d0
      endif

      if(qs > qeps )then
        CALL com_gamma( 4.0_RP + d , tmp_factor )
        ls= ( PI * ros * nos / ( ro * qs ) ) ** 0.25
        ws= c * tmp_factor / ( 6.0d0 * ( ls ** d ) )
        ws= 1.0d-2*ws * rofactor
      else
        ws = 0.0d0
      endif

      if(qg > qeps )then
        CALL com_gamma( 4.5_RP , tmp_factor )
        lg= ( PI * rog * nog / ( ro * qg ) ) ** 0.25
        wg= tmp_factor * ( ( ( 4.0d0 * GRAV * 100.0d0 * rog )/( 3.0d0 * Cd * ro ) ) ** 0.5 )
        wg= 1.0d-2*wg / ( 6.0d0 * ( lg ** 0.5 ) )
      else
        wg = 0.0d0
      endif

      !Reflectivity weighted terminal velocity. 
      wt = ( wr * zr + ws * zs + wg * zg )/ ( zr + zs + zg )

      ELSE

      wt=0.0d0

      ENDIF

    else if (METHOD_REF_CALC == 3) then
      !Observation operator from Xue et al 2007
      !Asumes power law between Z and mass concentration of different 
      !hydrometeor categories. 
      !Derived for X-Band radars tacking into account Mie scattering.
      !Includes a computation of the attenuation. (k)

      MAXF=0.5d0

      !First we need to compute the mixtures between rain, snow and hail
      !Following Jung et al 2007 eq (2) and (3)
      Fg=0.0d0
      Fs=0.0d0
      fwg=0.0d0
      fws=0.0d0
      IF( qr .GT. qeps .AND. qg .GT. qeps)THEN
        Fg=MAXF * ( min( qr/qg , qg/qr ) )**(1.0d0/3.0d0)
        fwg= qr / ( qr + qg )
      ENDIF
      IF( qr .GT. qeps .AND. qs .GT. qeps)THEN
        Fs=MAXF * ( min( qr/qs , qs/qr ) )**(1.0d0/3.0d0)
        fws= qr / ( qr + qs )
      ENDIF

      if ( .not. USE_METHOD3_REF_MELT ) then
        Fs = 0.0_RP
        Fg = 0.0_RP
      endif 

      !Correct the rain, snow and hail mixing ratios assuming
      !that we have a mixture due to melting.

      qrp=(1.0d0-Fs-Fg)*qr

      qsp=(1.0d0-Fs)*qs

      qgp=(1.0d0-Fg)*qg

      !Compute the new species resulting from melting.

      qms=Fs * (qr + qs) !Melting snow concentration.

      qmg=Fg * (qr + qg) !Melting hail concentration.

      !Compute reflectivities for each species including the melting species.

      IF( qrp .GT. qeps)THEN
      zr= 2.53d4 * ( ro * qrp * 1.0d3 )**1.84
      ENDIF
      IF( qsp .GT. qeps)THEN
      zs= 3.48d3 * ( ro * qsp * 1.0d3 )**1.66
      ENDIF
      IF( qgp .GT. qeps)THEN
      zg= 5.54d3 * ( ro * qgp * 1.0d3 )**1.70  !!! graupel (A.Amemiya 2019)
      ENDIF
      IF( qms .GT. qeps )THEN
      zms=( 0.00491 + 5.75*fws - 5.588*(fws**2) )*1.0d5
      zms= zms * ( ro * qms * 1.0d3 )**( 1.67 - 0.202*fws + 0.398*(fws**2) )

      ENDIF
      IF( qmg .GT. qeps )THEN
      zmg=( 0.0358 + 5.27*fwg -9.51*(fwg**2) + 4.68 *(fwg**3) )*1.0d5
      zmg= zmg * ( ro * qmg * 1.0d3 )**( 1.70 + 0.020*fwg + 0.287 * (fwg**2) - 0.186*(fwg**3) ) !!! graupel (A. Amemiya 2020)
       ENDIF

      ref = zr +  zg  + zs + zms + zmg

      !Compute reflectivity weigthed terminal velocity.
      !Lin et al 1983. (The distribution parameters are 
      !consistent with the work of Jung et al 2007)

      !!! graupel paramters and terminal velocity equations are modified to be
      !!! consistent with Tomita2008 default settings

      IF( ref > 0.0d0 )THEN
        !There are hidrometeors, compute their terminal velocity.
        !Units according to Lin et al 1983.
        nor=8.0d-2      ![cm^-4]
        nos=3.0d-2      ![cm^-4]
        nog=4.0d-2      ![cm^-4]          !!!!! Tomita 2008
        ror=1.0d0        ![g/cm3]
        ros=0.1d0        ![g/cm3]
        rog=0.400d0      ![g/cm3]         !!!!! Tomita 2008
        roo=1.28 * 0.001d0      ![g/cm3] Surface air density.  !!! Tomita 2008 
        ro=1.0d-3 * ro
        Cd=0.6d0         !!!!! drag_g in SCALE

        cr_t08=130.0d0 ![m**1-b / s]     !!! SCALE default
        dr_t08=0.5d0
        cs_t08=4.84d0 ![m**1-b / s]
        ds_t08=0.25d0
        dg_t08=0.5d0  !!! SCALE default


        rofactor= ( roo / ro  ) ** 0.5

        IF ( qr .GT. qeps )THEN
        lr= ( PI * ror * nor / ( ro * qr ) ) ** 0.25 !!! [cm ^-1]
        CALL com_gamma( 4.0_RP + dr_t08 , tmp_factor )
        wr= cr_t08 * tmp_factor / ( 6.0d0 * ( ( lr * 1.0e2 ) ** dr_t08 ) ) !!! [m/s]
        wr= wr * rofactor
        ELSE
        wr=0.0d0
        ENDIF

        IF( qs .GT. qeps )THEN
          IF  ( USE_T08_RS2014 ) then
            tc=min(-0.1_RP, t-273.15_RP)
            MOMs_0bs = exp(log(10.0_RP) * loga_(tc,2.0_RP) + log(ro * qs / as_RS14) * b_(tc,2.0_RP) ) 
            ws = cs_t08 * rofactor * exp(log(10.0_RP) * loga_(tc,2.25_RP) + log(ro * qs / as_RS14) * b_(tc,2.25_RP) ) / MOMs_0bs
          ELSE
            ls= ( PI * ros * nos / ( ro * qs ) ) ** 0.25 !!! [cm ^-1]
            CALL com_gamma( 4.0_RP + ds_t08 , tmp_factor )
            ws= cs_t08 * tmp_factor / ( 6.0d0 * ( ( ls * 1.0e2 ) ** ds_t08 ) ) !!! [m/s]
            ws= ws * rofactor
          ENDIF
        ELSE
          ws=0.0d0
        ENDIF


        IF ( qg .GT. qeps )THEN
        lg= ( PI * rog * nog / ( ro * qg ) ) ** 0.25
        CALL com_gamma( 4.0_RP + dg_t08 , tmp_factor )
        wg = tmp_factor * ( ( ( 4.0d0 * GRAV * rog )/( 3.0d0 * Cd * roo ) ) ** 0.5 )  !!! fixed 2021.5.14
        wg = wg / ( 6.0d0 * ( ( lg * 1.0e2 ) ** dg_t08 ) )   !!! [m/s]
        wg= wg * rofactor
        ELSE
        wg=0.0d0
        ENDIF

        !Reflectivity weighted terminal velocity. 
        !The melting species are assumed to fail as their non-melting counterparts.
        !however this might not be the case for melting snow.
        wt = ( wr * zr + ws *  zs +  ws * zms + wg *  zg + wg * zmg ) / ( zr + zs + zg + zms + zmg )

      ELSE

        wt=0.0d0

      ENDIF

    else  !IF OVER DIFFERENT OPTIONS

      WRITE(6,*)'[Error] Not recognized method for radar reflectivity and wind computation'
      STOP

    end if ! [METHOD_REF_CALC == ?]


    !Compute radial velocity
    !WRITE(6,*)'ICRV',u,v,w,wt,az,elev
    vr = u * cos(elev*D2R) * sin(az*D2R)
    vr = vr + v * cos(elev*D2R) * cos(az*D2R)
    IF( USE_TERMINAL_VELOCITY )THEN
      vr = vr + (w - wt)*sin(elev*D2R)
    ELSE
      vr = vr + (w)*sin(elev*D2R)
    ENDIF

    !!! NaN check 
    if (.not.(abs(vr) .lt. 1.0e20)) then
      write(*,*) '***ERROR*** vr is NaN' 
      write(*,*) vr
      write(*,*) wr,wg,ws
      write(*,*) lr,lg,ls
      write(*,*) zr,zg,zs,zms,zmg
      write(*,*) u,v,w,wt
      WRITE(6,*) elev , az , D2R
      stop
    elseif (.not. (abs(ref) .lt. 1.0e20))  then
      write(*,*) '***ERROR*** ref is NaN' 
      write(*,*) ref
      write(*,*) zr,  zg, zs, zms, zmg
      stop 
    elseif (ref < 0.0)  then
      write(*,*) '***ERROR*** ref is negative' 
      write(*,*) ref
      write(*,*) zr,  zg, zs, zms, zmg
      stop
    end if  

    return

  contains 
    !-----------------------------------------------------------------------
    ! GAMMA FUNCTION
    !-----------------------------------------------------------------------
    !==================================================
    !       Purpose: Compute the gamma function 但(x)
    !       Input :  x  --- Argument of 但(x)
    !                       ( x is not equal to 0,-1,-2,炭炭炭 )
    !       Output:  GA --- 但(x)
    !==================================================
    ! Proff. Jianming Jin
    ! Department of Electrical and Computer Engineering
    ! University of Illinois
    subroutine com_gamma(x,ga)
      implicit none
      real(RP) :: x , ga
      real(RP) :: g(26)
      real(RP) :: z , r , gr
      integer  :: m1, k , m

      !dimension g(26)
      if (x.eq.int(x)) then
         if (x.gt.0.0d0) then
            ga=1.0d0
            m1=x-1
            do k=2,m1
               ga=ga*k
            end do
         else
            ga=huge(1.0_RP)
         endif
      else

         if (abs(x).gt.1.0_RP) then
            z=abs(x)
            m=int(z)
            r=1.0_RP
            do k=1,m
               r=r*(z-k)
            end do
            z=z-m
         else
            z=x
         endif
         g(:) = (/               &
           1.0,                  &
           0.5772156649015329e0, &
          -0.6558780715202538e0, &
          -0.420026350340952e-1, &
           0.1665386113822915e0, &
          -0.421977345555443e-1, &
          -0.96219715278770e-2,  &
           0.72189432466630e-2,  &
          -0.11651675918591e-2,  &
          -0.2152416741149e-3,   &
           0.1280502823882e-3,   &
          -0.201348547807e-4,    &
          -0.12504934821e-5,     &
           0.11330272320e-5,     &
          -0.2056338417e-6,      &
           0.61160950e-8,        &
           0.50020075e-8,        &
          -0.11812746e-8,        &
           0.1043427e-9,         &
           0.77823e-11,          &
          -0.36968e-11,          &
           0.51e-12,             &
          -0.206e-13,            &   
          -0.54e-14,             &
           0.14e-14,             &
           0.1e-15               /)
         gr=g(26)
         do k=25,1,-1
            gr=gr*z+g(k)
         end do
         ga=1.0d0/(gr*z)
         if (abs(x).gt.1.0_RP) then
            ga=ga*r
            if (x.lt.0.0_RP) ga=-pi/(x*ga*sin(pi*x))
         endif
      endif

      return
    end subroutine com_gamma

    function loga_(tems, nm)
      real(RP) :: loga_
      real(RP) :: tems
      real(RP) :: nm
      real(RP),  parameter   :: coef_a01 =  5.065339_RP
      real(RP),  parameter   :: coef_a02 = -0.062659_RP
      real(RP),  parameter   :: coef_a03 = -3.032362_RP
      real(RP),  parameter   :: coef_a04 =  0.029469_RP
      real(RP),  parameter   :: coef_a05 = -0.000285_RP
      real(RP),  parameter   :: coef_a06 =  0.31255_RP
      real(RP),  parameter   :: coef_a07 =  0.000204_RP
      real(RP),  parameter   :: coef_a08 =  0.003199_RP
      real(RP),  parameter   :: coef_a09 =  0.0_RP
      real(RP),  parameter   :: coef_a10 = -0.015952_RP
      real(RP) :: coef_at(4)     

      coef_at(1) = coef_a01 + tems * ( coef_a02 + tems * ( coef_a05 + tems * coef_a09 ) )
      coef_at(2) = coef_a03 + tems * ( coef_a04 + tems *   coef_a07 )
      coef_at(3) = coef_a06 + tems *   coef_a08
      coef_at(4) = coef_a10
      loga_  = coef_at(1) + nm * ( coef_at(2) + nm * ( coef_at(3) + nm * coef_at(4) ) )

    end function loga_

    function b_(tems, nm)
      real(RP) :: b_
      real(RP) :: tems
      real(RP) :: nm
      real(RP), parameter   :: coef_b01 =  0.476221_RP
      real(RP), parameter   :: coef_b02 = -0.015896_RP
      real(RP),  parameter   :: coef_b03 =  0.165977_RP
      real(RP),  parameter   :: coef_b04 =  0.007468_RP
      real(RP),  parameter   :: coef_b05 = -0.000141_RP
      real(RP),  parameter   :: coef_b06 =  0.060366_RP
      real(RP),  parameter   :: coef_b07 =  0.000079_RP
      real(RP),  parameter   :: coef_b08 =  0.000594_RP
      real(RP),  parameter   :: coef_b09 =  0.0_RP
      real(RP),  parameter   :: coef_b10 = -0.003577_RP
      real(RP) :: coef_bt(4)     

      coef_bt(1) = coef_b01 + tems * ( coef_b02 + tems * ( coef_b05 + tems * coef_b09 ) )
      coef_bt(2) = coef_b03 + tems * ( coef_b04 + tems *   coef_b07 )
      coef_bt(3) = coef_b06 + tems *   coef_b08
      coef_bt(4) = coef_b10
      b_  = coef_bt(1) + nm * ( coef_bt(2) + nm * ( coef_bt(3) + nm * coef_bt(4) ) )

     end function b_

  end subroutine calc_ref_vr

  !-------------------------------------------------------------------------------
  ! Find local observations to be used for a targeted grid
  !-------------------------------------------------------------------------------
  ! [INPUT]
  !   ri      : horizontal i-grid cooridnate of the targeted grid
  !   rj      : horizontal j-grid cooridnate of the targeted grid
  !   rlev    : vertical pressure of the targeted grid
  !   rz      : vertical height   of the targeted grid
  !   nvar    : variable index of the targeted grid
  !   srch_q0 : (optional) first guess of the multiplier of incremental search distances
  ! [OUT]
  !   hdxf    : fcstast ensemble perturbations in the observation space
  !   rdiag   : localization-weighted observation error variances
  !   rloc    : localization weights
  !   dep     : observation departure
  !   nobsl   : number of valid observations (in hdxf, rdiag, rloc, dep)
  !   depd    : (optional) observation departure for the deterministic run
  !   nobsl_t : (optional) number of assimilated observations wrt. observation variables/types
  !   cutd_t  : (optional) cutoff distance of assimilated observations wrt. observation variables/types
  !   srch_q0 : (optional) revised first guess of the multiplier of incremental search distances for the next call
  !-------------------------------------------------------------------------------
  subroutine obs_local(obs, ri, rj, rlev, rz, nvar, hdxf, rdiag, rloc, dep, nobsl, depd, nobsl_t, cutd_t, srch_q0)
    use scale_sort, only: &
      SORT_quickselect_desc_arg, &
      SORT_quickselect_arg
    implicit none

    type(obs_info), intent(in) :: obs(:)
    real(RP), intent(in) :: ri, rj, rlev, rz
    integer, intent(in) :: nvar
    real(RP), intent(out) :: hdxf(nobstotal,nmem)
    real(RP), intent(out) :: rdiag(nobstotal)
    real(RP), intent(out) :: rloc(nobstotal)
    real(RP), intent(out) :: dep(nobstotal)
    integer, intent(out) :: nobsl
    real(RP), intent(out), optional :: depd(nobstotal)
    integer, intent(out), optional :: nobsl_t(nid_obs,nobtype)
    real(RP), intent(out), optional :: cutd_t(nid_obs,nobtype)
    integer, intent(inout), optional :: srch_q0(nctype)

    integer, allocatable :: nobs_use(:)
    integer, allocatable :: nobs_use2(:)
    real(RP), allocatable :: dist_tmp(:)
    real(RP), allocatable :: rloc_tmp(:)
    real(RP), allocatable :: rdiag_tmp(:)

    real(RP) :: nrloc, nrdiag
    real(RP) :: ndist_dummy
    integer :: iob, ityp, ielm
    integer :: imin, imax, jmin, jmax
    integer :: ic, ic2, icm
    integer :: n, nn, nn_prev
    integer :: nobsl_prev, nobsl_incr
    integer :: nobsl_max_master
    integer :: ielm_u_master, ityp_master

    integer :: q
    logical :: loop
    integer :: nn_steps(n_merge_max+1)
    logical :: reach_cutoff
    integer :: imin_cutoff(n_merge_max), imax_cutoff(n_merge_max)
    integer :: jmin_cutoff(n_merge_max), jmax_cutoff(n_merge_max)
    real(RP) :: search_incr(n_merge_max)
    real(RP) :: search_incr_i(n_merge_max), search_incr_j(n_merge_max)
    real(RP) :: dist_cutoff_fac, dist_cutoff_fac_square

    real(RP) :: cutd
    integer :: nobsl_cm(n_merge_max)

    !-----------------------------------------------------------------------------
    ! Initialize
    !-----------------------------------------------------------------------------

    nobsl = 0
    if (present(nobsl_t)) then
      nobsl_t(:,:) = 0
    end if
    if (present(cutd_t)) then
      cutd_t(:,:) = 0.0d0
      if (MAX_NOBS_PER_GRID_CRITERION == 1) then
        do ic = 1, nctype
          cutd_t(elm_u_ctype(ic),typ_ctype(ic)) = hori_loc_ctype(ic) * dist_zero_fac
        end do
      end if
    end if

    if (nobstotal == 0) then
      return
    end if

    if (maxval(MAX_NOBS_PER_GRID(:)) > 0) then
      allocate (nobs_use (nobstotal))
      allocate (nobs_use2(nobstotal))
      allocate (rloc_tmp (nobstotal))
      allocate (rdiag_tmp(nobstotal))
      if (MAX_NOBS_PER_GRID_CRITERION == 1) then
        allocate (dist_tmp (nobstotal))
      end if
  !    dist_tmp(:) = -1.0d6
      rloc_tmp(:) = -1.0d6
  !    rdiag_tmp(:) = -1.0d6
    else
      allocate (nobs_use(maxnobs_per_ctype))
    end if

    !-----------------------------------------------------------------------------
    ! For each observation type,
    ! do rough data search by a rectangle using the sorting mesh, and then
    ! do precise data search by normalized 3D distance and variable localization.
    !-----------------------------------------------------------------------------

    do ic = 1, nctype
      if (n_merge(ic) == 0) then
        if (present(cutd_t)) then
          cutd_t(elm_u_ctype(ic),typ_ctype(ic)) = 0.0d0
        end if
        cycle
      end if

      nobsl_max_master = MAX_NOBS_PER_GRID(typ_ctype(ic)) ! Use the number limit setting of the "master" obs type for all group of obs types
      ielm_u_master = elm_u_ctype(ic)                     ! Count observation numbers    at the "master" obs type for all group of obs types
      ityp_master = typ_ctype(ic)                         !

      if (nobsl_max_master <= 0) then
      !---------------------------------------------------------------------------
      ! When obs number limit is not enabled,
      ! directly prepare (hdxf, dep, depd, rdiag, rloc) output.
      !---------------------------------------------------------------------------

        nobsl_prev = nobsl

        do icm = 1, n_merge(ic)
          ic2 = ic_merge(icm,ic)
          ielm = elm_ctype(ic2)
          ityp = typ_ctype(ic2)

          if (obsgrd(ic2)%tot_ext > 0) then
            nn = 0
            call obs_local_range(ic2, ri, rj, imin, imax, jmin, jmax)
            call obs_choose_ext(ic2, imin, imax, jmin, jmax, nn, nobs_use)

            do n = 1, nn
              iob = nobs_use(n)

              call obs_local_cal(obs(:),ri, rj, rlev, rz, nvar, iob, ic2, ndist_dummy, nrloc, nrdiag)
              if (nrloc == 0.0d0) cycle

              nobsl = nobsl + 1
              hdxf(nobsl,:) = obsda_sort%ensval(1:nmem,iob)
              rdiag(nobsl) = nrdiag
              rloc(nobsl) = nrloc
              dep(nobsl) = obsda_sort%val(iob)
              if (present(depd)) then
                depd(nobsl) = obsda_sort%ensval(mmdetobs,iob)
              end if
            end do ! [ n = 1, nn ]
          end if ! [ obsgrd(ic2)%tot_ext > 0 ]

          if (present(nobsl_t)) then
            nobsl_t(elm_u_ctype(ic2),ityp) = nobsl - nobsl_prev
          end if
        end do ! [ do icm = 1, n_merge(ic) ]

      !---------------------------------------------------------------------------
      else if (MAX_NOBS_PER_GRID_CRITERION == 1) then
      !---------------------------------------------------------------------------
      ! When obs number limit is enabled and the sorting criterion is simply distance,
      ! try the incremental observation location search
      ! (within the localization cutoff area) before conduct selection.
      !---------------------------------------------------------------------------

        nn = 0
        do icm = 1, n_merge(ic)
          ic2 = ic_merge(icm,ic)

          if (obsgrd(ic2)%tot_ext > 0) then
            call obs_local_range(ic2, ri, rj, imin, imax, jmin, jmax)
            call obs_choose_ext(ic2, imin, imax, jmin, jmax, nn, nobs_use)
          end if ! [ obsgrd(ic2)%tot_ext > 0 ]
        end do ! [ do icm = 1, n_merge(ic) ]

        if (nn == 0) cycle

        ! Determine the search_incr based on the master obs ctype:
        ! (zero-weight distance) / (# of search increment), but not smaller than the sorting mesh size
        search_incr(1) = hori_loc_ctype(ic) * dist_zero_fac / real(n_search_incr, RP)  ! (unit: m)
        search_incr(1) = max(search_incr(1), obsgrd(ic)%grdspc_i, obsgrd(ic)%grdspc_j)     ! (unit: m)

        do icm = 1, n_merge(ic)
          ic2 = ic_merge(icm,ic)
          ! Determine the search_incr of the other obs ctypes based on its horizontal localization scale
          ! relative to that of the master obs ctype
          if (icm > 1) then
            search_incr(icm) = search_incr(1) / hori_loc_ctype(ic) * hori_loc_ctype(ic2)
          end if
          search_incr_i(icm) = search_incr(icm) / DX ! (unit: x-grid)
          search_incr_j(icm) = search_incr(icm) / DY ! (unit: y-grid)
          call obs_local_range(ic2, ri, rj, imin_cutoff(icm), imax_cutoff(icm), jmin_cutoff(icm), jmax_cutoff(icm))
        end do

        nobsl_incr = 0
        if (present(srch_q0)) then
          q = srch_q0(ic) - 1
        else
          q = 0
        end if
        loop = .true.

        do while (loop)
          q = q + 1
          nn = 0
          reach_cutoff = .true.

          do icm = 1, n_merge(ic)
            ic2 = ic_merge(icm,ic)
            nn_steps(icm) = nn

            if (obsgrd(ic2)%tot_ext > 0) then
              call ij_obsgrd_ext(ic2, ri-search_incr_i(icm)*q, rj-search_incr_j(icm)*q, imin, jmin)
              call ij_obsgrd_ext(ic2, ri+search_incr_i(icm)*q, rj+search_incr_j(icm)*q, imax, jmax)

              if (imin <= imin_cutoff(icm) .and. imax >= imax_cutoff(icm) .and. &
                  jmin <= jmin_cutoff(icm) .and. jmax >= jmax_cutoff(icm)) then
                imin = imin_cutoff(icm)
                imax = imax_cutoff(icm)
                jmin = jmin_cutoff(icm)
                jmax = jmax_cutoff(icm)
              else
                reach_cutoff = .false.
              end if

              call obs_choose_ext(ic2, imin, imax, jmin, jmax, nn, nobs_use)
            end if ! [ obsgrd(ic2)%tot_ext > 0 ]
          end do ! [ do icm = 1, n_merge(ic) ]

          nn_steps(n_merge(ic)+1) = nn

          if ((.not. reach_cutoff) .and. nn < nobsl_max_master) cycle
          if (reach_cutoff) then
            loop = .false.
            if (nn == 0) exit
          end if

          nobsl_incr = 0

          ! Determine the cutoff fraction in this incremental search,
          ! which should be the same for all obs ctypes based on its definition
          dist_cutoff_fac = search_incr(1) * q / hori_loc_ctype(ic)
          dist_cutoff_fac_square = dist_cutoff_fac * dist_cutoff_fac

          do icm = 1, n_merge(ic)
            ic2 = ic_merge(icm,ic)

            if (nn_steps(icm+1) > nn_steps(icm)) then
              do n = nn_steps(icm)+1, nn_steps(icm+1)
                iob = nobs_use(n)

                if (rloc_tmp(iob) == 0.0d0) cycle

                if (rloc_tmp(iob) < 0.0d0) then
                  call obs_local_cal(obs(:),ri, rj, rlev, rz, nvar, iob, ic2, dist_tmp(iob), rloc_tmp(iob), rdiag_tmp(iob))
                  if (rloc_tmp(iob) == 0.0d0) cycle
                end if

                if (.not. reach_cutoff) then
                  if (dist_tmp(iob) > dist_cutoff_fac_square) cycle
                end if

                nobsl_incr = nobsl_incr + 1
                nobs_use2(nobsl_incr) = iob
              end do
            end if ! [ nn_steps(icm+1) > nn_steps(icm) ]
          end do ! [ do icm = 1, n_merge(ic) ]

          if (nobsl_incr >= nobsl_max_master) loop = .false.
        end do ! [ loop ]

        if (present(srch_q0)) then
          if (q == srch_q0(ic) .and. nobsl_incr > nobsl_max_master * 3) then ! when (nobsl_incr >= nobsl_max_master) too soon, decrease srch_q0
            srch_q0(ic) = q - 1
          else if (q > srch_q0(ic)) then ! when (nobsl_incr >= nobsl_max_master) too late, increase srch_q0
            srch_q0(ic) = q
          end if
        end if

        if (nobsl_incr == 0) cycle

        if (nobsl_incr > nobsl_max_master) then
          call SORT_quickselect_arg(dist_tmp, nobs_use2, 1, nobsl_incr, nobsl_max_master)
          nobsl_incr = nobsl_max_master
        end if

        do n = 1, nobsl_incr
          nobsl = nobsl + 1
          iob = nobs_use2(n)
          hdxf(nobsl,:) = obsda_sort%ensval(1:nmem,iob)
          rdiag(nobsl) = rdiag_tmp(iob)
          rloc(nobsl) = rloc_tmp(iob)
          dep(nobsl) = obsda_sort%val(iob)
          if (present(depd)) then
            depd(nobsl) = obsda_sort%ensval(mmdetobs,iob)
          end if
        end do

        if (present(nobsl_t)) then
          nobsl_t(ielm_u_master,ityp_master) = nobsl_incr
        end if
        if (present(cutd_t)) then
          if (nobsl_incr == nobsl_max_master) then
            cutd_t(ielm_u_master,ityp_master) = hori_loc_ctype(ic) * sqrt(dist_tmp(nobs_use2(nobsl_incr)))
          end if
        end if

      !---------------------------------------------------------------------------
      else
      !---------------------------------------------------------------------------
      ! When obs number limit is enabled and the sorting criterion is NOT distance,
      ! conduct selection to all observations within the localization cutoff area.
      !---------------------------------------------------------------------------

        nn = 0
        nobsl_incr = 0

        do icm = 1, n_merge(ic)
          ic2 = ic_merge(icm,ic)

          if (obsgrd(ic2)%tot_ext > 0) then
            nn_prev = nn
            call obs_local_range(ic2, ri, rj, imin, imax, jmin, jmax)
            call obs_choose_ext(ic2, imin, imax, jmin, jmax, nn, nobs_use)

            do n = nn_prev+1, nn
              iob = nobs_use(n)

              call obs_local_cal(obs(:),ri, rj, rlev, rz, nvar, iob, ic2, ndist_dummy, rloc_tmp(iob), rdiag_tmp(iob))
              if (rloc_tmp(iob) == 0.0d0) cycle

              nobsl_incr = nobsl_incr + 1
              nobs_use2(nobsl_incr) = iob
            end do
          end if ! [ obsgrd(ic2)%tot_ext > 0 ]
        end do ! [ do icm = 1, n_merge(ic) ]

        if (nobsl_incr == 0) cycle

        if (nobsl_incr > nobsl_max_master) then
          if (MAX_NOBS_PER_GRID_CRITERION == 2) then
            call SORT_quickselect_desc_arg(rloc_tmp, nobs_use2, 1, nobsl_incr, nobsl_max_master)
          else if (MAX_NOBS_PER_GRID_CRITERION == 3) then
            call SORT_quickselect_arg(rdiag_tmp, nobs_use2, 1, nobsl_incr, nobsl_max_master)
          else
            write (6, '(A,I6)') "[Error] Unsupported 'MAX_NOBS_PER_GRID_CRITERION':", MAX_NOBS_PER_GRID_CRITERION
            stop 99
          end if
          nobsl_incr = nobsl_max_master
        end if

        do n = 1, nobsl_incr
          nobsl = nobsl + 1
          iob = nobs_use2(n)
          hdxf(nobsl,:) = obsda_sort%ensval(1:nmem,iob)
          rdiag(nobsl) = rdiag_tmp(iob)
          rloc(nobsl) = rloc_tmp(iob)
          dep(nobsl) = obsda_sort%val(iob)
          if (present(depd)) then
            depd(nobsl) = obsda_sort%ensval(mmdetobs,iob)
          end if
        end do

        if (present(nobsl_t)) then
          nobsl_t(ielm_u_master,ityp_master) = nobsl_incr
        end if
        if (present(cutd_t)) then
          if (nobsl_incr == nobsl_max_master) then
            if (MAX_NOBS_PER_GRID_CRITERION == 2) then
              cutd_t(ielm_u_master,ityp_master) = rloc_tmp(nobs_use2(nobsl_incr))
            else if (MAX_NOBS_PER_GRID_CRITERION == 3) then
              cutd_t(ielm_u_master,ityp_master) = rdiag_tmp(nobs_use2(nobsl_incr))
            end if
          end if
        end if

      !---------------------------------------------------------------------------
      end if

    end do ! [ ic = 1, nctype ]

    !-----------------------------------------------------------------------------
    ! Finalize
    !-----------------------------------------------------------------------------

    deallocate (nobs_use)
    if (maxval(MAX_NOBS_PER_GRID(:)) > 0) then
      deallocate (nobs_use2)
      deallocate (rloc_tmp)
      deallocate (rdiag_tmp)
      if (MAX_NOBS_PER_GRID_CRITERION == 1) then
        deallocate (dist_tmp)
      end if
    end if

    return
  end subroutine obs_local

  !-------------------------------------------------------------------------------
  ! Calculate the range of the rectangle that covers the (horizontal) localization
  ! cut-off length in the extended subdomain, given the observation type
  !-------------------------------------------------------------------------------
  subroutine obs_local_range(ctype, ri, rj, imin, imax, jmin, jmax)
    implicit none
    integer, intent(in) :: ctype
    real(RP), intent(in) :: ri, rj
    integer, intent(out) :: imin, imax, jmin, jmax

    real(RP) :: dist_zero_i, dist_zero_j


    if (nlon==1) then !!! 2-D
      dist_zero_i = 0.0_RP
    else
      dist_zero_i = hori_loc_ctype(ctype) * dist_zero_fac / DX
    end if
    dist_zero_j = hori_loc_ctype(ctype) * dist_zero_fac / DY
    call ij_obsgrd_ext(ctype, ri - dist_zero_i, rj - dist_zero_j, imin, jmin)
    call ij_obsgrd_ext(ctype, ri + dist_zero_i, rj + dist_zero_j, imax, jmax)

    return
  end subroutine obs_local_range

  !-------------------------------------------------------------------------------
  ! Subroutine for main calculation of obs_local
  !-------------------------------------------------------------------------------
  subroutine obs_local_cal(obs, ri, rj, rlev, rz, nvar, iob, ic, ndist, nrloc, nrdiag)
    implicit none
    type(obs_info), intent(in) :: obs(:)
    real(RP), intent(in) :: ri, rj, rlev, rz ! coordinate of the targeted model grid
    integer, intent(in) :: nvar         ! index of targeted model variable
    integer, intent(in) :: iob          ! index of observation in obsda_sort
    integer, intent(in) :: ic           ! observation combined type
    real(RP), intent(out) :: ndist  ! normalized 3D distance SQUARE       (in case of rejected obs: -1.)
    real(RP), intent(out) :: nrloc  ! localization weight                 (in case of rejected obs:  0.)
    real(RP), intent(out) :: nrdiag ! weighted observation error variance (in case of rejected obs: -1.)

    integer :: obelm           ! observation variable type
    integer :: obtyp           ! observation report type
    integer :: obset
    integer :: obidx
    real(RP) :: rdx, rdy, rdz
    real(RP) :: nd_h, nd_v ! normalized horizontal/vertical distances

    integer :: di, dj, dk

    nrloc = 0.0d0
    nrdiag = -1.0d0
    ndist = -1.0d0

    obelm = elm_ctype(ic)
    obtyp = typ_ctype(ic)
    obset = obsda_sort%set(iob)
    obidx = obsda_sort%idx(iob)
    !
    ! Calculate variable localization
    !
    if (nvar > 0) then  ! use variable localization only when nvar > 0
      nrloc = var_local(nvar,uid_obs_varlocal(obelm))

      !--- reject obs by variable localization
      if (nrloc < tiny(var_local)) then
        nrloc = 0.0d0
        return
      end if
    end if
    !
    ! Calculate normalized vertical distances
    !
    if (vert_loc_ctype(ic) == 0.0d0) then
      nd_v = 0.0d0                                                            ! no vertical localization
    else if (obelm == id_ps_obs) then
      nd_v = ABS(LOG(obs(obset)%dat(obidx)) - LOG(rlev)) / vert_loc_ctype(ic) ! for ps, use observed ps value for the base of vertical localization
    else if (obelm == id_rain_obs) then
      nd_v = ABS(LOG(VERT_LOCAL_RAIN_BASE) - LOG(rlev)) / vert_loc_ctype(ic)  ! for rain, use VERT_LOCAL_RAIN_BASE for the base of vertical localization
    else if (obtyp == 22) then ! obtypelist(obtyp) == 'PHARAD'
      nd_v = ABS(obs(obset)%lev(obidx) - rz) / vert_loc_ctype(ic)             ! for PHARAD, use z-coordinate for vertical localization
    else
      nd_v = ABS(LOG(obs(obset)%lev(obidx)) - LOG(rlev)) / vert_loc_ctype(ic)
    end if

    !--- reject obs by normalized vertical distance
    !    (do this first because there is large possibility to reject obs by the vertical distrance)
    if (nd_v > dist_zero_fac) then
      nrloc = 0.0d0
      return
    end if
    !
    ! Calculate normalized horizontal distances
    !
    rdx = (ri - obs(obset)%ri(obidx)) * DX
    rdy = (rj - obs(obset)%rj(obidx)) * DY
    nd_h = sqrt(rdx*rdx + rdy*rdy) / hori_loc_ctype(ic)

    !--- reject obs by normalized horizontal distance
    if (nd_h > dist_zero_fac) then
      nrloc = 0.0d0
      return
    end if
    !
    ! Calculate (normalized 3D distances)^2
    !
    ndist = nd_h * nd_h + nd_v * nd_v

    !--- reject obs by normalized 3D distance
    if (ndist > dist_zero_fac_square) then
      nrloc = 0.0d0
      ndist = -1.0d0
      return
    end if

    if ( obtyp == 22 .and. ( RADAR_THIN_LETKF_METHOD > 0 ) ) then ! obtypelist(obtyp) == 'PHARAD'
      rdz = obs(obset)%lev(obidx) - rz 


      select case( RADAR_THIN_LETKF_METHOD )
      case( 1 )
        ! Pick up nearest 8 obs (2x2x2)
        ! and then choose every HGRID/VGRID
        di = int( abs( rdx / RADAR_SO_SIZE_HORI ) )
        dj = int( abs( rdy / RADAR_SO_SIZE_HORI ) )
        dk = int( abs( obs(obset)%lev(obidx) - rz ) / RADAR_SO_SIZE_VERT ) 

        if ( ( mod( di, RADAR_THIN_LETKF_HGRID ) /= 0 .or. &
               mod( dj, RADAR_THIN_LETKF_HGRID ) /= 0 .or. &
               mod( dk, RADAR_THIN_LETKF_VGRID ) /= 0 ) .and. &
              ( ( di >= RADAR_THIN_LETKF_HNEAR ) .or. &
                ( dj >= RADAR_THIN_LETKF_HNEAR ) .or. &
                ( dk >= RADAR_THIN_LETKF_VNEAR ) ) ) then
          nrloc = 0.0d0
          ndist = -1.0d0
          return
        endif
      case( 2 )
        ! Pick up nearest 1 obs 
        ! and then choose every HGRID/VGRID
        di = nint( rdx / RADAR_SO_SIZE_HORI )
        dj = nint( rdy / RADAR_SO_SIZE_HORI )
        dk = nint( obs(obset)%lev(obidx) - rz ) / RADAR_SO_SIZE_VERT 

        if ( mod( di, RADAR_THIN_LETKF_HGRID ) /= 0 .or. &
             mod( dj, RADAR_THIN_LETKF_HGRID ) /= 0 .or. &
             mod( dk, RADAR_THIN_LETKF_VGRID ) /= 0 ) then
          nrloc = 0.0d0
          ndist = -1.0d0
          return
        endif

      case default
        ! No thinning
      end select

    endif 


    !
    ! Calculate observational localization
    !
    nrloc = nrloc * EXP(-0.5d0 * ndist)
    !
    ! Calculate (observation variance / localization)
    !
    nrdiag = obs(obset)%err(obidx) * obs(obset)%err(obidx) / nrloc
    if ( RADAR_PQV .and. obelm == id_radar_ref_obs .and. obsda_sort%tm(iob) < 0.0d0 ) then
      nrdiag = OBSERR_PQ**2 / nrloc
    endif

    return
  end subroutine obs_local_cal

  !-------------------------------------------------------------------------------
  ! Find local observations to be used for a targeted grid
  !-------------------------------------------------------------------------------
  ! [INPUT]
  !   srch_q0 : (optional) first guess of the multiplier of incremental search distances
  ! [OUT]
  !   hdxf    : fcstast ensemble perturbations in the observation space
  !   rdiag   : localization-weighted observation error variances
  !   rloc    : localization weights
  !   dep     : observation departure
  !   nobsl   : number of valid observations (in hdxf, rdiag, rloc, dep)
  !   depd    : (optional) observation departure for the deterministic run
  !-------------------------------------------------------------------------------
  subroutine obs_pest_etkf(hdxf, rdiag, rloc, dep, nobsl, depd)
    implicit none

    real(RP), intent(out) :: hdxf(nobstotal,nmem)
    real(RP), intent(out) :: rdiag(nobstotal)
    real(RP), intent(out) :: rloc(nobstotal)
    real(RP), intent(out) :: dep(nobstotal)
    integer, intent(out) :: nobsl
    real(RP), intent(out), optional :: depd(nobstotal)

    integer :: obset
    integer :: obidx

    real(RP) :: nrloc, nrdiag
    integer :: iob
    integer :: n, nn

    !-----------------------------------------------------------------------------
    ! For each observation type,
    ! do rough data search by a rectangle using the sorting mesh, and then
    ! do precise data search by normalized 3D distance and variable localization.
    !-----------------------------------------------------------------------------

    nobsl = 0
    nrloc = 1.0 ! ETKF ! No localization
    nn = obsda_sort%nobs ! obsda_sort%nobs = nobstotal

    do n = 1, nn
      iob = n
      obset = obsda_sort%set(n)
      obidx = obsda_sort%idx(n)
      nrdiag = obs(obset)%err(obidx) * obs(obset)%err(obidx) / nrloc
      if (nrloc == 0.0) cycle

      nobsl = nobsl + 1
      hdxf(nobsl,:) = obsda_sort%ensval(1:nmem,iob)
      rdiag(nobsl) = nrdiag
      rloc(nobsl) = nrloc
      dep(nobsl) = obsda_sort%val(iob)
      if (present(depd)) then
        depd(nobsl) = obsda_sort%ensval(mmdetobs,iob)
      end if
    end do

    return
  end subroutine obs_pest_etkf

  !-------------------------------------------------------------------------------
  ! Relaxation parameter based on grid locations (not for covariance inflation purpose)
  !-------------------------------------------------------------------------------
  subroutine relax_beta(ri, rj, rz, beta)
    implicit none
    real(RP), intent(in) :: ri, rj, rz
    real(RP), intent(out) :: beta
    real(RP) :: dist_bdy

    beta = 1.0d0
    !
    ! Upper bound of updates when RADAR_ZMAX is set and only radar observations are assimilated
    !
    if (radar_only .and. rz > RADAR_ZMAX + max(VERT_LOCAL(22), VERT_LOCAL_RADAR_VR) * dist_zero_fac) then
      beta = 0.0d0
      return
    end if
    !
    ! Boundary buffer
    !
    if (BOUNDARY_BUFFER_WIDTH > 0.0d0) then
      dist_bdy = min(min(ri-xhalo, nlong+xhalo+1-ri) * DX, &
                     min(rj-yhalo, nlatg+yhalo+1-rj) * DY) / BOUNDARY_BUFFER_WIDTH
      if (dist_bdy < 1.0d0) then
        beta = max(dist_bdy, 0.0d0)
      end if
    end if

    return
  end subroutine relax_beta

  !-------------------------------------------------------------------------------
  ! Relaxation via LETKF weight - RTPP method
  !-------------------------------------------------------------------------------
  subroutine weight_RTPP(w, infl, wrlx)
    implicit none
    real(RP), intent(in) :: w(nmem,nmem)
    real(RP), intent(in) :: infl
    real(RP), intent(out) :: wrlx(nmem,nmem)
    integer :: m

    wrlx = (1.0d0 - RELAX_ALPHA) * w
    do m = 1, nmem
      wrlx(m,m) = wrlx(m,m) + RELAX_ALPHA * sqrt(infl)
    end do

    return
  end subroutine weight_RTPP

  !-------------------------------------------------------------------------------
  ! Relaxation via LETKF weight - RTPS method
  !-------------------------------------------------------------------------------
  subroutine weight_RTPS(w, pa, xb, infl, wrlx, infl_out)
    implicit none
    real(RP), intent(in) :: w(nmem,nmem)
    real(RP), intent(in) :: pa(nmem,nmem)
    real(RP), intent(in) :: xb(nmem)
    real(RP), intent(in) :: infl
    real(RP), intent(out) :: wrlx(nmem,nmem)
    real(RP), intent(out) :: infl_out
    real(RP) :: var_g, var_a
    integer :: m, k

    var_g = 0.0d0
    var_a = 0.0d0
    do m = 1, nmem
      var_g = var_g + xb(m) * xb(m)
      do k = 1, nmem
        var_a = var_a + xb(k) * pa(k,m) * xb(m)
      end do
    end do
    if (var_g > 0.0d0 .and. var_a > 0.0d0) then
      infl_out = RELAX_ALPHA_SPREAD * sqrt(var_g * infl / (var_a * real(nmem-1,RP))) &  ! Whitaker and Hamill 2012
               - RELAX_ALPHA_SPREAD + 1.0d0                                                   !
      wrlx = w * infl_out
    else
      wrlx = w
      infl_out = 1.0d0
    end if

    return
  end subroutine weight_RTPS

  !-------------------------------------------------------------------------------
  ! Relaxation via LETKF weight - RTPS method for parameter estimation (alpha = 1.0)
  !-------------------------------------------------------------------------------
  subroutine weight_RTPS_const(w, pa, xb, wrlx, infl_out)
    implicit none
    real(RP), intent(in) :: w(nmem,nmem)
    real(RP), intent(in) :: pa(nmem,nmem)
    real(RP), intent(in) :: xb(nmem)
    real(RP), intent(out) :: wrlx(nmem,nmem)
    real(RP), intent(out) :: infl_out
    real(RP) :: var_g, var_a
    integer :: m, k

    real(RP), parameter :: RTPS_const = 1.0d0

    var_g = 0.0d0
    var_a = 0.0d0
    do m = 1, nmem
      var_g = var_g + xb(m) * xb(m)
      do k = 1, nmem
        var_a = var_a + xb(k) * pa(k,m) * xb(m)
      end do
    end do
    if (var_g > 0.0d0 .and. var_a > 0.0d0) then
      infl_out = RTPS_const * sqrt(var_g * 1.0d0 / (var_a * real(nmem-1,kind=RP))) &  ! Whitaker and Hamill 2012
               - RTPS_const + 1.0d0                                                   !
      wrlx = w * infl_out
    else
      wrlx = w
      infl_out = 1.0d0
    end if

    return
  end subroutine weight_RTPS_const

  !-----------------------------------------------------------------------
  ! DISTANCE BETWEEN TWO POINTS (LONa,LATa)-(LONb,LATb)
  !-----------------------------------------------------------------------
  subroutine com_distll_1(alon,alat,blon,blat,dist)
    use scale_const, only: &
      PI => CONST_PI, &
      RE => CONST_RADIUS
    implicit none
    real(RP), intent(in)  :: alon
    real(RP), intent(in)  :: alat
    real(RP), intent(in)  :: blon
    real(RP), intent(in)  :: blat
    real(RP), intent(out) :: dist
    real(RP), parameter   :: r180 = 1.0_RP/180.0_RP
    real(RP) :: lon1,lon2,lat1,lat2
    real(RP) :: cosd

    lon1 = alon * pi * r180
    lon2 = blon * pi * r180
    lat1 = alat * pi * r180
    lat2 = blat * pi * r180

    cosd = sin(lat1)*sin(lat2) + cos(lat1)*cos(lat2)*cos(lon2-lon1)
    cosd = min( 1._RP,cosd)
    cosd = max(-1._RP,cosd)

    dist = acos( cosd ) * re

    return
  end subroutine com_distll_1

  !-------------------------------------------------------------------------------
  ! Calculate ensemble mean (on scattered grids)
  !-------------------------------------------------------------------------------
  ! [INPUT]
  !   mem                     : ensemble size
  !   nens                    : ensemble demension of state variables
  !   nij                     : scattered grid numbers
  !   v3d(nij,nlev,nens,nv3d) : 3D ensemble state variables (on scattered grids)
  !                             inputted by (:,:,1..mem,:)
  !   v2d(nij,     nens,nv3d) : 2D ensemble state variables (on scattered grids)
  !                             inputted by (:,  1..mem,:)
  ! [OUTPUT]
  !   v3d(nij,nlev,nens,nv3d) : ensemble mean of 3D state variables (on scattered grids)
  !                             outputted by (:,:,mem+1,:)
  !   v2d(nij,     nens,nv3d) : ensemble mean of 2D state variables (on scattered grids)
  !                             outputted by (:  ,mem+1,:)
  !-------------------------------------------------------------------------------
  subroutine ensmean_grd(mem, nens, nij, v3d, v2d)
    implicit none
    integer,  intent(in)    :: mem
    integer,  intent(in)    :: nens
    integer,  intent(in)    :: nij
    real(RP), intent(inout) :: v3d(nij,nlev,nens,nv3d)
    real(RP), intent(inout) :: v2d(nij,nens,nv2d)

    integer :: i, k, m, n, mmean
    !---------------------------------------------------------------------

    mmean = mem + 1

    do n = 1, nv3d
    do k = 1, nlev
    do i = 1, nij
       v3d(i,k,mmean,n) = v3d(i,k,1,n)
       do m = 2, mem
          v3d(i,k,mmean,n) = v3d(i,k,mmean,n) + v3d(i,k,m,n)
       end do
       v3d(i,k,mmean,n) = v3d(i,k,mmean,n) / real(mem, kind=RP)
    end do
    end do
    end do
    do n = 1, nv2d
    do i = 1, nij
       v2d(i,mmean,n) = v2d(i,1,n)
       do m = 2, mem
          v2d(i,mmean,n) = v2d(i,mmean,n) + v2d(i,m,n)
       end do
       v2d(i,mmean,n) = v2d(i,mmean,n) / real(mem, kind=RP)
    end do
    end do

    return
  end subroutine ensmean_grd

  !-----------------------------------------------------------------------
  ! Transformation from model variables to an observation
  !
  ! stggrd: grid type of u and v
  !  0: non-staggered grid
  !  1: staggered grid
  !-----------------------------------------------------------------------
  subroutine Trans_XtoY(elm,ri,rj,rk,lon,lat,v3d,v2d,yobs,qc,stggrd)
    use scale_const, only: &
      UNDEF => CONST_UNDEF, &
      FVIRT => CONST_EPSTVap
    implicit none
    integer,  intent(in)  :: elm
    real(RP), intent(in)  :: ri,rj,rk
    real(RP), intent(in)  :: lon,lat
    real(RP), intent(in)  :: v3d(:,:,:,:)
    real(RP), intent(in)  :: v2d(  :,:,:)
    real(RP), intent(out) :: yobs
    integer,  intent(out) :: qc
    integer,  intent(in), optional :: stggrd
    real(RP) :: u,v,t,q,topo

    integer :: stggrd_ = 0
    if (present(stggrd)) stggrd_ = stggrd

    yobs = undef
    qc = iqc_good

    select case (elm)
    case(id_u_obs,id_v_obs)  ! U,V
      if (stggrd_ == 1) then
        call itpl_3d(v3d(:,:,:,iv3dd_u),rk,ri-0.5_RP,rj,u)  !###### should modity itpl_3d to prevent '1.0' problem....??
        call itpl_3d(v3d(:,:,:,iv3dd_v),rk,ri,rj-0.5_RP,v)  !######
      else
        call itpl_3d(v3d(:,:,:,iv3dd_u),rk,ri,rj,u)
        call itpl_3d(v3d(:,:,:,iv3dd_v),rk,ri,rj,v)
      end if
      if (elm == id_u_obs) then
        yobs = u
      else
        yobs = v
      end if
    case(id_t_obs)  ! T
      call itpl_3d(v3d(:,:,:,iv3dd_t),rk,ri,rj,yobs)
    case(id_tv_obs)  ! Tv
      call itpl_3d(v3d(:,:,:,iv3dd_t),rk,ri,rj,yobs)
      call itpl_3d(v3d(:,:,:,iv3dd_q),rk,ri,rj,q)
      yobs = yobs * (1.0d0 + fvirt * q)
    case(id_q_obs)  ! Q
      call itpl_3d(v3d(:,:,:,iv3dd_q),rk,ri,rj,yobs)
    case(id_ps_obs) ! PS
      call itpl_2d(v2d(:,:,iv2dd_t2m),ri,rj,t)
      call itpl_2d(v2d(:,:,iv2dd_q2m),ri,rj,q)
      call itpl_2d(v2d(:,:,iv2dd_topo),ri,rj,topo)
      call itpl_2d(v2d(:,:,iv2dd_ps),ri,rj,yobs)
      call prsadj(yobs,rk-topo,t,q)
      if( abs(rk-topo) > PS_ADJUST_THRES ) then
        qc = iqc_ps_ter
      end if
    case(id_rh_obs) ! RH
      call itpl_3d(v3d(:,:,:,iv3dd_rh),rk,ri,rj,yobs)
    case default
      qc = iqc_otype
    end select

    return
  end subroutine Trans_XtoY

  subroutine Trans_XtoY_radar(elm,radar_lon,radar_lat,radar_z,ri,rj,rk,lon,lat,lev,v3d,v2d,yobs,qc,stggrd)
    use scale_const, only: &
      UNDEF => CONST_UNDEF, &
      D2R   => CONST_D2R,   &
      R2D   => CONST_R2D
    use scale_mapprojection, only: &
      MAPPROJECTION_rotcoef
    implicit none

    integer,  intent(in)  :: elm
    real(RP), intent(in)  :: ri,rj,rk,radar_lon,radar_lat,radar_z !!!!! Use only, ri, rj, rk eventually... (radar_lon,lat,z in ri,rj,rk)
    real(RP), intent(in)  :: lon,lat,lev
    real(RP), intent(in)  :: v3d(:,:,:,:)
    real(RP), intent(in)  :: v2d(:,:,:)
    real(RP), intent(out) :: yobs
    integer,  intent(out) :: qc
    integer,  intent(in), optional :: stggrd

    integer :: stggrd_ = 0

    real(RP) :: qvr,qcr,qrr,qir,qsr,qgr,ur,vr,wr,tr,pr !,rhr
    real(RP) :: dist , dlon , dlat , az , elev , radar_ref,radar_rv

    real(RP) :: utmp, vtmp
    real(RP) :: rotc(1,1,2)
    real(RP) :: lon_tmp(1,1),lat_tmp(1,1)

    if (present(stggrd)) stggrd_ = stggrd

    yobs = undef
    qc = iqc_good

    if (stggrd_ == 1) then
      CALL itpl_3d(v3d(:,:,:,iv3dd_u),rk,ri-0.5_RP,rj,ur)  !###### should modity itpl_3d to prevent '1.0' problem....??
      CALL itpl_3d(v3d(:,:,:,iv3dd_v),rk,ri,rj-0.5_RP,vr)  !######
      CALL itpl_3d(v3d(:,:,:,iv3dd_w),rk-0.5_RP,ri,rj,wr)  !######
    else
      CALL itpl_3d(v3d(:,:,:,iv3dd_u),rk,ri,rj,ur)
      CALL itpl_3d(v3d(:,:,:,iv3dd_v),rk,ri,rj,vr)
      CALL itpl_3d(v3d(:,:,:,iv3dd_w),rk,ri,rj,wr)
    end if
    CALL itpl_3d(v3d(:,:,:,iv3dd_t),rk,ri,rj,tr)
    CALL itpl_3d(v3d(:,:,:,iv3dd_p),rk,ri,rj,pr)
    CALL itpl_3d(v3d(:,:,:,iv3dd_q),rk,ri,rj,qvr)
    CALL itpl_3d(v3d(:,:,:,iv3dd_qc),rk,ri,rj,qcr)
    CALL itpl_3d(v3d(:,:,:,iv3dd_qr),rk,ri,rj,qrr)
    CALL itpl_3d(v3d(:,:,:,iv3dd_qi),rk,ri,rj,qir)
    CALL itpl_3d(v3d(:,:,:,iv3dd_qs),rk,ri,rj,qsr)
    CALL itpl_3d(v3d(:,:,:,iv3dd_qg),rk,ri,rj,qgr)

    ! Compute az and elevation for the current observation.
    ! Simple approach (TODO: implement a more robust computation)

    ! Azimuth
    dlon = lon-radar_lon
    dlat = lat-radar_lat
    if ( dlon == 0.0_RP .and. dlat == 0.0_RP  )then
      qc = iqc_out_h
      return
    else
      az = R2D*atan2(dlon*cos(radar_lat*D2R),dlat)
    endif
    if( az < 0 ) az = 360.0_RP + az
    ! elevation
    call com_distll_1(lon,lat,radar_lon,radar_lat,dist)
    elev = R2D*atan2(lev-radar_z,dist)

    lon_tmp(1,1) = lon*D2R
    lat_tmp(1,1) = lat*D2R
    call MAPPROJECTION_rotcoef(1, 1, 1, 1, 1, 1, &
                               lon_tmp(:,:),lat_tmp(:,:),rotc(:,:,1),rotc(:,:,2))

    utmp = ur
    vtmp = vr
    ur = utmp * rotc(1,1,1) - vtmp * rotc(1,1,2)
    vr = utmp * rotc(1,1,2) + vtmp * rotc(1,1,1)

    ! Check that the azimuth and elevation angles are within the expected range.
    ! Some grid points may be at the radar location.

    CALL calc_ref_vr(qvr,qcr,qrr,qir,qsr,qgr,ur,vr,wr,tr,pr,az,elev,radar_ref,radar_rv)

    select case (elm)
    case(id_radar_ref_obs,id_radar_ref_zero_obs)
      if (radar_ref < MIN_RADAR_REF) then
        qc = iqc_ref_low
        yobs = MIN_RADAR_REF_DBZ + LOW_REF_SHIFT  !!! even if the above qc is bad, still return the value
      else
        yobs = 10.0_RP * log10(radar_ref)
      end if
    case(id_radar_vr_obs)
      if (radar_ref < MIN_RADAR_REF) then
        qc = iqc_ref_low
      end if
      yobs = radar_rv  !!! even if the above qc is bad, still return the value
    case default
      qc = iqc_otype
    end select

    return
  end subroutine Trans_XtoY_radar

  !-----------------------------------------------------------------------
  ! Pressure adjustment for a different height level
  !-----------------------------------------------------------------------
  subroutine prsadj(p,dz,t,q)
    use scale_const, only: &
      GG => CONST_GRAV, &
      RD => CONST_Rdry
    implicit none
    real(RP), intent(inout) :: p
    real(RP), intent(in) :: dz ! height difference (target - original) [m]
    real(RP), intent(in) :: t  ! temperature [K] at original level
    real(RP), intent(in) :: q  ! humidity [kg/kg] at original level
    real(RP), parameter  :: gamma=5.0d-3 ! lapse rate [K/m]
    real(RP) :: tv

    if(dz /= 0) then
      tv = t * (1.0d0 + 0.608d0 * q)
      p = p * ((-gamma*dz+tv)/tv)**(gg/(gamma*rd)) !tv is at original level
    end if

    return
  end subroutine prsadj

  !-----------------------------------------------------------------------
  ! Coordinate conversion (find rk in pressure levels)
  !
  ! rk = 0.0d0  : surface observation
  !-----------------------------------------------------------------------
  subroutine phys2ijk(p_full,elem,ri,rj,rlev,rk,qc)
    use scale_const, only: &
      UNDEF => CONST_UNDEF
    implicit none

    real(RP), intent(in)  :: p_full(:,:,:)
    integer,  intent(in)  :: elem
    real(RP), intent(in)  :: ri
    real(RP), intent(in)  :: rj
    real(RP), intent(in)  :: rlev ! pressure levels (for 3D variable only)
    real(RP), intent(out) :: rk
    integer,  intent(out) :: qc

    real(RP) :: ak
    real(RP) :: lnps(nlevh,nlonh,nlath)
    real(RP) :: plev(nlevh)
    real(RP) :: ptmp
    integer  :: i,j,k, ii, jj, ks

    qc = iqc_good
    !
    ! rlev -> rk
    !
    if( ri < 1.0_RP .or. ri > nlonh .or. rj < 1.0_RP .or. rj > nlath ) then
      rk = undef
      qc = iqc_out_h
      return ! [Warning] observation is outside of the horizontal domain
    end if
    !
    if(elem > 9999) then ! surface observation
      rk = rlev
    else
      !
      ! horizontal interpolation
      !
      i = ceiling(ri)
      j = ceiling(rj)
      !
      ! Find the lowest valid level
      !
      ks = 1+zhalo
      do jj = j-1, j
      do ii = max(i-1,1), min(i,nlong+xhalo)
        do k=1+zhalo,nlev+zhalo
              if( p_full(k,ii,jj) >= 0.0_RP ) exit
        end do
        if (k > ks) ks = k
      end do
      end do

      lnps(:,i-1:i,j-1:j) = log(p_full(:,i-1:i,j-1:j))
      call itpl_2d_column(lnps,ri,rj,plev)
      !
      ! Log pressure
      !
      rk = log(rlev)
      !
      ! determine if rk is within bound.
      !
      if(rk < plev(nlev+zhalo)) then
        call itpl_2d(p_full(nlev+zhalo,:,:),ri,rj,ptmp)
        rk = undef
        qc = iqc_out_vhi
        return ! [Warning] observation is too high: ptop=', ptmp, ', lev=', rlev, ', elem=', elem
      end if
      if(rk > plev(ks)) then
        call itpl_2d(p_full(ks,:,:),ri,rj,ptmp)
        rk = undef
        qc = iqc_out_vlo
        return ! '[Warning] observation is too low: pbottom=', ptmp, ', lev=', rlev, ', elem=', elem
      end if
      !
      ! find rk
      !
      do k=ks+1,nlev+zhalo
        if(plev(k) < rk) exit ! assuming descending order of plev
      end do

      if (k == nlev+zhalo+1) then
        ak = 0.99
        rk = real(k-2,kind=RP) + ak
      else
        ak = (rk - plev(k-1)) / (plev(k) - plev(k-1))
        rk = real(k-1,kind=RP) + ak
      end if
    end if

    return
  end subroutine phys2ijk

  !-----------------------------------------------------------------------
  ! Coordinate conversion (find rk in height levels)
  !
  ! rk = 0.0d0  : surface observation
  !-----------------------------------------------------------------------
  subroutine phys2ijkz(z_full,ri,rj,rlev,rk,qc)
    use scale_const, only: &
      UNDEF => CONST_UNDEF
    implicit none

    real(RP), intent(in)  :: z_full(:,:,:)
    real(RP), intent(in)  :: ri
    real(RP), intent(in)  :: rj
    real(RP), intent(in)  :: rlev ! height levels
    real(RP), intent(out) :: rk
    integer,  intent(out) :: qc

    real(RP) :: ak
    real(RP) :: zlev(nlevh)
    real(RP) :: ztmp
    integer :: i,j,k, ii, jj, ks

    qc = iqc_good
    !
    ! rlev -> rk
    !
    if( ri < 1.0_RP .or. ri > nlonh .or. rj < 1.0_RP .or. rj > nlath ) then
      rk = undef
      qc = iqc_out_h
      return ! '[Warning] observation is outside of the horizontal domain'
    end if
    !
    ! horizontal interpolation
    !
    i = ceiling(ri)
    j = ceiling(rj)
    !
    ! Find the lowest valid level
    !
    ks = 1+zhalo
    do jj = j-1, j
    do ii = max(i-1,1), min(i,nlong+xhalo)
      do k=1+zhalo,nlev+zhalo
        if( z_full(k,ii,jj) > -300.0_RP .and. z_full(k,ii,jj) < 10000.0_RP ) exit
      end do
      if (k > ks) ks = k
    end do
    end do

    call itpl_2d_column(z_full,ri,rj,zlev)

    !
    ! determine if rlev is within bound.
    !
    if( rlev > zlev(nlev+zhalo) ) then
      call itpl_2d(z_full(nlev+zhalo,:,:),ri,rj,ztmp)
      rk = undef
      qc = iqc_out_vhi
      return ! '[Warning] observation is too high: ztop=', ztmp, ', lev=', rlev
    end if
    if(rlev < zlev(ks)) then
      call itpl_2d(z_full(ks,:,:),ri,rj,ztmp)
      rk = undef
      qc = iqc_out_vlo
      return ! '[Warning] observation is too low: zbottom=', ztmp, ', lev=', rlev
    end if
    !
    ! find rk
    !
    do k=ks+1,nlev+zhalo
      if(zlev(k) > rlev) exit ! assuming ascending order of zlev
    end do

    if (k == nlev+zhalo+1) then
      ak = 0.99
      rk = real(k-2,kind=RP) + ak
    else
      ak = (rlev - zlev(k-1)) / (zlev(k) - zlev(k-1))
      rk = real(k-1,kind=RP) + ak
    end if

    return
  end subroutine phys2ijkz

  !-----------------------------------------------------------------------
  ! Coordinate conversion
  !-----------------------------------------------------------------------
  subroutine phys2ij(rlon,rlat,rig,rjg)
    use scale_const, only: &
        PI => CONST_PI
    use scale_atmos_grid_cartesC, only: &
        CXG => ATMOS_GRID_CARTESC_CXG, &
        CYG => ATMOS_GRID_CARTESC_CYG
    use scale_mapprojection, only: &
        MAPPROJECTION_lonlat2xy
    implicit none
    real(RP), intent(in)  :: rlon
    real(RP), intent(in)  :: rlat
    real(RP), intent(out) :: rig
    real(RP), intent(out) :: rjg
    real(RP) :: rig_RP
    real(RP) :: rjg_RP
    !
    ! rlon,rlat -> ri,rj
    !
    call MAPPROJECTION_lonlat2xy( real(rlon*PI/180.0_RP, kind=RP), &
                                  real(rlat*PI/180.0_RP, kind=RP), rig_RP, rjg_RP )
    rig = real((rig_RP - CXG(1)) / DX, kind=RP) + 1.0_RP
    rjg = real((rjg_RP - CYG(1)) / DY, kind=RP) + 1.0_RP

    if (nlonh==1) then !!! adjustment : ideal 2-D case
      rig = 1.0_RP
    end if

    return
  end subroutine phys2ij

  !-----------------------------------------------------------------------
  ! Interpolation
  !-----------------------------------------------------------------------
  subroutine itpl_2d(var,ri,rj,var5)
    implicit none
    real(RP), intent(in)  :: var(:,:)
    real(RP), intent(in)  :: ri
    real(RP), intent(in)  :: rj
    real(RP), intent(out) :: var5
    real(RP) :: ai,aj
    integer :: i,j

    i = ceiling(ri)
    ai = ri - real(i-1,kind=RP)
    j = ceiling(rj)
    aj = rj - real(j-1,kind=RP)

    if (nlonh==1) then
    var5 = var(i  ,j-1) * (1-aj) &
       & + var(i  ,j  ) *    aj
    else
    var5 = var(i-1,j-1) * (1-ai) * (1-aj) &
       & + var(i  ,j-1) *    ai  * (1-aj) &
       & + var(i-1,j  ) * (1-ai) *    aj  &
       & + var(i  ,j  ) *    ai  *    aj
    end if

    return
  end subroutine itpl_2d

  subroutine itpl_2d_column(var,ri,rj,var5)
    implicit none
    real(RP), intent(in)  :: var(:,:,:)
    real(RP), intent(in)  :: ri
    real(RP), intent(in)  :: rj
    real(RP), intent(out) :: var5(:)

    real(RP) :: ai,aj
    integer :: i,j

    i = ceiling(ri)
    ai = ri - real(i-1,kind=RP)
    j = ceiling(rj)
    aj = rj - real(j-1,kind=RP)

    if (nlonh==1) then
    var5(:) = var(:,i  ,j-1) * (1-aj) &
          & + var(:,i  ,j  ) *    aj
    else
    var5(:) = var(:,i-1,j-1) * (1-ai) * (1-aj) &
          & + var(:,i  ,j-1) *    ai  * (1-aj) &
          & + var(:,i-1,j  ) * (1-ai) *    aj  &
          & + var(:,i  ,j  ) *    ai  *    aj
    end if

    return
  end subroutine itpl_2d_column

  subroutine itpl_3d(var,rk,ri,rj,var5)
    implicit none
    real(RP), intent(in)  :: var(:,:,:)
    real(RP), intent(in)  :: ri
    real(RP), intent(in)  :: rj
    real(RP), intent(in)  :: rk
    real(RP), intent(out) :: var5
    real(RP) :: ai,aj,ak
    integer :: i,j,k

    i = ceiling(ri)
    ai = ri - real(i-1,kind=RP)
    j = ceiling(rj)
    aj = rj - real(j-1,kind=RP)
    k = ceiling(rk)
    ak = rk - real(k-1,kind=RP)

    if (nlonh==1) then
    var5 = var(k-1,i  ,j-1) * (1-aj) * (1-ak) &
       & + var(k-1,i  ,j  ) *    aj  * (1-ak) &
       & + var(k,  i  ,j-1) * (1-aj) *    ak  &
       & + var(k,  i  ,j  ) *    aj  *    ak
    else
    var5 = var(k-1,i-1,j-1) * (1-ai) * (1-aj) * (1-ak) &
       & + var(k-1,i  ,j-1) *    ai  * (1-aj) * (1-ak) &
       & + var(k-1,i-1,j  ) * (1-ai) *    aj  * (1-ak) &
       & + var(k-1,i  ,j  ) *    ai  *    aj  * (1-ak) &
       & + var(k,  i-1,j-1) * (1-ai) * (1-aj) *    ak  &
       & + var(k,  i  ,j-1) *    ai  * (1-aj) *    ak  &
       & + var(k,  i-1,j  ) * (1-ai) *    aj  *    ak  &
       & + var(k,  i  ,j  ) *    ai  *    aj  *    ak
    end if

    return
  end subroutine itpl_3d

  subroutine merge_sort_parallel(n, array)
    use omp_lib
    integer(8), intent(in) :: n
    integer(8), intent(inout) :: array(n)
    logical omp_nested
    integer maxnest

#ifdef DA
    omp_nested = omp_get_nested()
    maxnest = floor(log(dble(omp_get_max_threads())) / log(2.0d0))
    call omp_set_nested(.true.)
    call merge_sort_2threads(n, array, 0, maxnest)
    call omp_set_nested(omp_nested)
#endif

  end subroutine merge_sort_parallel

  recursive subroutine merge_sort_2threads(n, array, nest, maxnest)
    use scale_sort, only: &
      SORT_quicksort
    implicit none
    integer(8), intent(in) :: n
    integer, intent(in) :: nest, maxnest
    integer(8), intent(inout) :: array(n)
    integer(8) :: asize(2)
    integer(8), allocatable :: tmpary(:)
    integer, parameter:: nmin = 4

    asize = n / 2
    if(mod(n,2_8) .ne. 0) asize(1) = asize(1) + 1

    if(nest < maxnest) then
       allocate(tmpary(n))
       tmpary(1:asize(1)) = array(1:asize(1))
       if(asize(1) > nmin) then
          call merge_sort_2threads(asize(1), tmpary(1:asize(1)), nest + 1, maxnest)
       else
          call SORT_quicksort(asize(1), tmpary(1:asize(1)))
       end if
       tmpary((asize(1) + 1):n) = array((asize(1) + 1):n)
       if(asize(2) > nmin) then
          call merge_sort_2threads(asize(2), tmpary((asize(1) + 1):n), nest + 1, maxnest)
       else
          call SORT_quicksort(asize(2), tmpary((asize(1) + 1):n))
       end if
       call merge_2threads(asize(1), tmpary(1:asize(1)), asize(2), tmpary((asize(1) + 1):n), n, array, nest, maxnest)
    else
       allocate(tmpary(n))
       tmpary = array
       call SORT_quicksort(asize(1), tmpary(1:asize(1)))
       call SORT_quicksort(asize(2), tmpary((asize(1) + 1):n))
       call merge(asize(1), tmpary(1:asize(1)), asize(2), tmpary((asize(1) + 1):n), n, array)
    end if
  end subroutine merge_sort_2threads

  recursive subroutine merge_2threads(n1, ary1, n2, ary2, n3, ary3, nest, maxnest)
    integer(8), intent(in) :: n1, ary1(n1), n2, ary2(n2), n3
    integer, intent(in) :: nest, maxnest
    integer(8), intent(inout) :: ary3(n3)
    integer(8) k, m
    integer, parameter :: threshold = 10

    if(nest >= maxnest .or. n1 < threshold .or. n2 < threshold) then
       call merge(n1, ary1, n2, ary2, n3, ary3)
       return
    end if

    k = n1 / 2
    m = binary_search_i8(n2, ary2, ary1(k))
    call merge_2threads(k, ary1(1:k), m, ary2(1:m), k + m, ary3(1:(k + m)), nest + 1, maxnest)
    call merge_2threads(n1 - k, ary1((k + 1):n1), n2 - m, ary2((m + 1):n2), n3 - (k + m), ary3((k + m + 1):n3), nest + 1, maxnest)

  end subroutine merge_2threads

  subroutine merge(n1, ary1, n2, ary2, n3, ary3)
    integer(8), intent(in) :: n1, ary1(n1), n2, ary2(n2), n3
    integer(8), intent(inout) :: ary3(n3)
    integer(8) i, j, k
    i = 1
    j = 1
    do k = 1, n3
       if(ary1(i) < ary2(j)) then
          ary3(k) = ary1(i)
          if(i == n1) then
             ary3((k + 1):n3) = ary2(j:n2)
             exit
          end if
          i = i + 1
       else
          ary3(k) = ary2(j)
          if(j == n2) then
             ary3((k + 1):n3) = ary1(i:n1)
             exit
          end if
          j = j + 1
       end if
    end do
  end subroutine merge

  recursive subroutine merge_sort_mpi(n, array, comm)
    use scale_sort, only: &
      SORT_quicksort
    implicit none
    !ONLY RANK 0 RETURNS RESULT
    integer(8), intent(in) :: n
    integer, intent(in) :: comm
    integer(8), intent(inout) :: array(n)
    integer(8) :: asize(2)
    integer(8), allocatable :: tmpary(:)
    integer :: parent, child, procs(2), newcomm, tag, nprocs, rank
    integer :: ierr

    call mpi_comm_size(comm, nprocs, ierr)
    call mpi_comm_rank(comm, rank, ierr)

    procs(1) = nprocs / 2
    procs(2) = nprocs - procs(1)
    parent = 0
    child = parent + procs(1)

    asize = n / 2
    if(mod(n,2_8) .ne. 0) asize(1) = asize(1) + 1

    allocate(tmpary(n))
    if(rank < child) then
       call MPI_comm_split(comm, 0, rank, newcomm, ierr)
       tmpary(1:asize(1)) = array(1:asize(1))
       if(procs(1) > 1) then
          call merge_sort_mpi(asize(1), tmpary(1:asize(1)), newcomm)
       else
          call SORT_quicksort(asize(1), tmpary(1:asize(1)))
       end if
    else
       call MPI_comm_split(comm, 1, rank, newcomm, ierr)
       tmpary((asize(1) + 1):n) = array((asize(1) + 1):n)
       if(procs(2) > 1) then
          call merge_sort_mpi(asize(2), tmpary((asize(1) + 1):n), newcomm)
       else
          call SORT_quicksort(asize(2), tmpary((asize(1) + 1):n))
       end if
    end if

    call merge_mpi_no_nest(asize(1), tmpary(1:asize(1)), asize(2), tmpary((asize(1) + 1):n), n, array, parent, child, comm)

  end subroutine merge_sort_mpi

  recursive subroutine merge_mpi(n1, ary1, n2, ary2, n3, ary3, parent, child, nprocs, comm)
    implicit none

    integer(8), intent(in) :: n1, ary1(n1), n2, ary2(n2), n3
    integer, intent(in) :: parent, child, nprocs, comm
    integer(8), intent(inout) :: ary3(n3)
    integer(8) k, m, pivot
    integer, parameter :: threshold = 10
    integer procs(2), grandchild1, grandchild2, rank, newcomm
    integer, parameter :: tag_p2c = 12345
    integer, parameter :: tag_p2g = 23451
    integer, parameter :: tag_c2g = 34512
    integer, parameter :: tag_p2h = 45123
    integer, parameter :: tag_c2h = 51234

    integer :: ierr

    call mpi_comm_rank(comm, rank, ierr)
    procs(1) = child - parent
    procs(2) = nprocs - procs(1)

    k = n1 / 2
    if(rank == parent) pivot = ary1(k)
    call MPI_bcast(pivot, 1, MPI_INTEGER8, parent, comm, ierr)
    if(rank == child) m = binary_search_i8(n2, ary2, pivot)
    call MPI_bcast(m, 1, MPI_INTEGER8, child, comm, ierr)

    if(procs(1) > 1 .and. n1 >= threshold) then
       grandchild1 = parent + (procs(1) / 2)
    else
       grandchild1 = parent
    end if
    if(procs(2) > 1 .and. n2 >= threshold) then
       grandchild2 = child + (procs(2) / 2)
    else
       grandchild2 = child
    end if

    if(rank >= parent .and. rank < child) then
       call MPI_comm_split(comm, 0, rank, newcomm, ierr)

       if(rank == parent) then
          if(rank == grandchild1) then
             call MPI_sendrecv(ary1((k + 1):n1), int(n1 - k), MPI_INTEGER8, grandchild2, tag_p2h, &
                  &            ary2(1:m), int(m), MPI_INTEGER8, child, tag_c2g, &
                  &            comm, MPI_STATUS_IGNORE, ierr)
          else
             call MPI_send(ary1((k + 1):n1), int(n1 - k,4), MPI_INTEGER8, grandchild2, &
                  &        tag_p2h, comm, ierr)
          end if
       else if(rank == grandchild1) then
          call MPI_recv(ary2(1:m), int(m), MPI_INTEGER8, child, &
               &        tag_c2g, comm, MPI_STATUS_IGNORE, ierr)
       end if

       if(procs(1) > 1 .and. n1 >= threshold) then
          call merge_mpi(k, ary1(1:k), m, ary2(1:m), k + m, ary3(1:(k + m)), 0, grandchild1 - parent, procs(1), newcomm)
       else
          call merge(k, ary1(1:k), m, ary2(1:m), k + m, ary3(1:(k + m)))
       end if
       if(rank == parent) then
          call MPI_recv(ary3((k + m + 1):n3), int(n3 - (k + m)), MPI_INTEGER8, child, &
               &        tag_p2c, comm, MPI_STATUS_IGNORE, ierr)
       end if
    else if(rank >= child .and. rank < parent + nprocs) then
       call MPI_comm_split(comm, 1, rank, newcomm, ierr)

       if(rank == child) then
          if(rank == grandchild2) then
             call MPI_sendrecv(ary2(1:m), int(m), MPI_INTEGER8, grandchild1, tag_c2g, &
                  &            ary1((k + 1):n1), int(n1 - k), MPI_INTEGER8, parent, tag_p2h, &
                  &            comm, MPI_STATUS_IGNORE, ierr)
          else
             call MPI_send(ary2(1:m), int(m,4), MPI_INTEGER8, grandchild1, &
                  &        tag_c2g, comm, ierr)
          end if
       else if(rank == grandchild2) then
          call MPI_recv(ary1((k + 1):n1), int(n1 - k), MPI_INTEGER8, parent, &
               &        tag_p2h, comm, MPI_STATUS_IGNORE, ierr)
       end if

       if(procs(2) > 1 .and. n2 >= threshold) then
          call merge_mpi(n2 - m, ary2((m + 1):n2), n1 - k, ary1((k + 1):n1), & !NOTE: FLIP SIDE
               &         n3 - (k + m), ary3((k + m + 1):n3), 0, grandchild2 - child, procs(2), newcomm)
       else
          call merge(n1 - k, ary1((k + 1):n1), n2 - m, ary2((m + 1):n2), &
               &     n3 - (k + m), ary3((k + m + 1):n3))
       end if
       if(rank == child) then
          call MPI_send(ary3((k + m + 1):n3), int(n3 - (k + m)), MPI_INTEGER8, parent, &
               &        tag_p2c, comm, ierr)
       end if
    else
       call MPI_comm_split(comm, 3, rank, newcomm, ierr) !SHOULD NOT HAPPEN (BUG)
       write(*, *) "something wrong in merge_mpi"
    end if
  end subroutine merge_mpi

  recursive subroutine merge_mpi_no_nest(n1, ary1, n2, ary2, n3, ary3, parent, child, comm)
    implicit none

    integer(8), intent(in) :: n1, ary1(n1), n2, ary2(n2), n3
    integer, intent(in) :: parent, child, comm
    integer(8), intent(inout) :: ary3(n3)
    integer(8) k, m, pivot
    integer, parameter :: threshold = 10
    integer rank
    integer, parameter :: tag_p2c = 12345

    integer(8) pivot_dummy

    integer :: ierr

    pivot_dummy=10

    call mpi_comm_rank(comm, rank, ierr)

    k = n1 / 2
    if(rank == parent) then
       pivot = ary1(k)

       call MPI_send(pivot, 1, MPI_INTEGER8, child, tag_p2c, comm, ierr)

       call MPI_recv(m, 1, MPI_INTEGER8, child, tag_p2c, comm, MPI_STATUS_IGNORE, ierr)
       call MPI_sendrecv(ary1((k + 1):n1), int(n1 - k), MPI_INTEGER8, child, tag_p2c, &
            &            ary2(1:m), int(m), MPI_INTEGER8, child, tag_p2c, &
            &            comm, MPI_STATUS_IGNORE, ierr)
       call merge(k, ary1(1:k), m, ary2(1:m), k + m, ary3(1:(k + m)))
       call MPI_recv(ary3((k + m + 1):n3), int(n3 - (k + m)), MPI_INTEGER8, child, &
               &        tag_p2c, comm, MPI_STATUS_IGNORE, ierr)
    else if(rank == child) then
       call MPI_recv(pivot, 1, MPI_INTEGER8, parent, tag_p2c, comm, MPI_STATUS_IGNORE, ierr)
       m = binary_search_i8(n2, ary2, pivot)
       call MPI_send(m, 1, MPI_INTEGER8, parent, tag_p2c, comm, ierr)
       call MPI_sendrecv(ary2(1:m), int(m), MPI_INTEGER8, parent, tag_p2c, &
            &            ary1((k + 1):n1), int(n1 - k), MPI_INTEGER8, parent, tag_p2c, &
            &            comm, MPI_STATUS_IGNORE, ierr)
       call merge(n1 - k, ary1((k + 1):n1), n2 - m, ary2((m + 1):n2), &
            &     n3 - (k + m), ary3((k + m + 1):n3))
       call MPI_send(ary3((k + m + 1):n3), int(n3 - (k + m)), MPI_INTEGER8, parent, &
            &        tag_p2c, comm, ierr)
    end if
  end subroutine merge_mpi_no_nest

  !-------------------------------------------------------------------------------
  ! Scatter gridded data to processes (nrank -> all)
  !-------------------------------------------------------------------------------
  subroutine scatter_grd_mpi(nrank,v3dg,v2dg,v3d,v2d)
    integer,intent(in) :: nrank
    real(RP),intent(in) :: v3dg(nlev,nlon,nlat,nv3d)
    real(RP),intent(in) :: v2dg(nlon,nlat,nv2d)
    real(RP),intent(out) :: v3d(nij1,nlev,nv3d)
    real(RP),intent(out) :: v2d(nij1,nv2d)
    real(RP) :: bufs(nij1max,nlevall,NPRC_ENS)
    real(RP) :: bufr(nij1max,nlevall)
    integer :: j,k,n,ierr,ns,nr

    ns = nij1max * nlevall
    nr = ns
    if( RANK_ENS == nrank ) then
      j=0
      do n=1,nv3d
      do k=1,nlev
        j = j+1
        call grd_to_buf( NPRC_ENS, v3dg(k,:,:,n), bufs(:,j,:) )
      end do
      end do

      do n=1,nv2d
        j = j+1
        call grd_to_buf( NPRC_ENS, v2dg(:,:,n),bufs(:,j,:) )
      end do
    end if

    CALL MPI_SCATTER( bufs, ns, datatype, bufr, nr, datatype, nrank, COMM_ENS, ierr )

    j=0
    do n=1,nv3d
    do k=1,nlev
      j = j+1
      v3d(:,k,n) = real( bufr(1:nij1,j), kind=RP )
    end do
    end do

    do n=1,nv2d
      j = j+1
      v2d(:,n) = real( bufr(1:nij1,j), kind=RP )
    end do

    return
  end subroutine scatter_grd_mpi

  !-------------------------------------------------------------------------------
  ! Scatter gridded data using MPI_ALLTOALL(V) (mstart~mend -> all)
  !-------------------------------------------------------------------------------
  subroutine scatter_grd_mpi_all2all(mstart,mend,v3dg,v2dg,v3d,v2d)
    integer,intent(in) :: mstart,mend
    real(RP),intent(in) :: v3dg(nlev,nlon,nlat,nv3d)
    real(RP),intent(in) :: v2dg(nlon,nlat,nv2d)
    real(RP),intent(inout) :: v3d(nij1,nlev,nens,nv3d)
    real(RP),intent(inout) :: v2d(nij1,nens,nv2d)
    real(RP) :: bufs(nij1max,nlevall,NPRC_ENS)
    real(RP) :: bufr(nij1max,nlevall,NPRC_ENS)
    integer :: k,n,j,m,mcount,ierr
    integer :: ns(NPRC_ENS),nst(NPRC_ENS),nr(NPRC_ENS),nrt(NPRC_ENS)

    mcount = mend - mstart + 1

    if(RANK_ENS < mcount) then
      j = 0
      do n = 1, nv3d
      do k = 1, nlev
        j = j+1
        call grd_to_buf(NPRC_ENS,v3dg(k,:,:,n),bufs(:,j,:))
      end do
      end do
      do n = 1, nv2d
        j = j+1
        call grd_to_buf(NPRC_ENS,v2dg(:,:,n),bufs(:,j,:))
      end do
    end if

    if(mcount == NPRC_ENS) then
      call MPI_ALLTOALL(bufs, nij1max*nlevall, datatype, &
                        bufr, nij1max*nlevall, datatype, COMM_ENS, ierr)
    else
      call set_all2allv_counts(mcount,nij1max*nlevall,NPRC_ENS,nr,nrt,ns,nst)
      call MPI_ALLTOALLV(bufs, ns, nst, datatype, &
                         bufr, nr, nrt, datatype, COMM_ENS, ierr)
    end if

    do m = mstart, mend
      j = 0
      do n = 1, nv3d
      do k = 1, nlev
        j = j+1
        v3d(:,k,m,n) = real(bufr(1:nij1,j,m-mstart+1),kind=RP)
      end do
      end do
      do n = 1, nv2d
        j = j+1
        v2d(:,m,n) = real(bufr(1:nij1,j,m-mstart+1),kind=RP)
      end do
    end do

    return
  end subroutine scatter_grd_mpi_all2all

  !-------------------------------------------------------------------------------
  ! Gather gridded data using MPI_ALLTOALL(V) (all -> mstart~mend)
  !-------------------------------------------------------------------------------
  subroutine gather_grd_mpi_all2all(mstart,mend,v3d,v2d,v3dg,v2dg)
    integer,intent(in) :: mstart,mend
    real(RP),intent(in) :: v3d(nij1,nlev,nens,nv3d)
    real(RP),intent(in) :: v2d(nij1,nens,nv2d)
    real(RP),intent(out) :: v3dg(nlev,nlon,nlat,nv3d)
    real(RP),intent(out) :: v2dg(nlon,nlat,nv2d)
    real(RP) :: bufs(nij1max,nlevall,NPRC_ENS)
    real(RP) :: bufr(nij1max,nlevall,NPRC_ENS)
    integer :: k,n,j,m,mcount,ierr
    integer :: ns(NPRC_ENS),nst(NPRC_ENS),nr(NPRC_ENS),nrt(NPRC_ENS)

    mcount = mend - mstart + 1

    do m = mstart, mend
      j = 0
      do n = 1, nv3d
      do k = 1, nlev
        j = j+1
        bufs(1:nij1,j,m-mstart+1) = real(v3d(:,k,m,n),kind=RP)
      end do
      end do
      do n=1,nv2d
        j = j+1
        bufs(1:nij1,j,m-mstart+1) = real(v2d(:,m,n),kind=RP)
      end do
    end do

    if(mcount == NPRC_ENS) then
      call MPI_ALLTOALL(bufs, nij1max*nlevall, datatype, &
                        bufr, nij1max*nlevall, datatype, COMM_ENS, ierr)
    else
      call set_all2allv_counts(mcount,nij1max*nlevall,NPRC_ENS,ns,nst,nr,nrt)
      call MPI_ALLTOALLV(bufs, ns, nst, datatype, &
                         bufr, nr, nrt, datatype, COMM_ENS, ierr)
    end if

    if(RANK_ENS < mcount) then
      j = 0
      do n = 1, nv3d
      do k = 1, nlev
        j = j+1
        call buf_to_grd(NPRC_ENS,bufr(:,j,:),v3dg(k,:,:,n))
      end do
      end do
      do n = 1, nv2d
        j = j+1
        call buf_to_grd(NPRC_ENS,bufr(:,j,:),v2dg(:,:,n))
      end do
    end if

    return
  end subroutine gather_grd_mpi_all2all

  !-------------------------------------------------------------------------------
  ! Set the send/recieve counts of MPI_ALLTOALLV
  !-------------------------------------------------------------------------------
  subroutine set_all2allv_counts(mcount,ngpblock,np,n_ens,nt_ens,n_mem,nt_mem)
    integer,intent(in) :: mcount,ngpblock
    integer,intent(in) :: np
    integer,intent(out) :: n_ens(np),nt_ens(np),n_mem(np),nt_mem(np)
    integer :: p

    n_ens = 0
    nt_ens = 0
    n_mem = 0
    nt_mem = 0
    do p = 1, mcount
      n_ens(p) = ngpblock
      if(RANK_ENS+1 == p) then
        n_mem(:) = ngpblock
      end if
    end do
    do p = 2, np
      nt_ens(p) = nt_ens(p-1) + n_ens(p-1)
      nt_mem(p) = nt_mem(p-1) + n_mem(p-1)
    end do

    return
  end subroutine set_all2allv_counts

  !-------------------------------------------------------------------------------
  ! gridded data -> buffer
  !-------------------------------------------------------------------------------
  subroutine grd_to_buf(np,grd,buf)
    use scale_const, only: &
      UNDEF => CONST_UNDEF
    integer,intent(in) :: np
    real(RP),intent(in) :: grd(nlon,nlat)
    real(RP),intent(out) :: buf(nij1max,np)
    integer :: i,j,m,ilon,ilat

    do m = 1, np
    do i = 1, nij1node(m)
      j = m-1 + np * (i-1)
      ilon = mod(j,nlon) + 1
      ilat = (j-ilon+1) / nlon + 1

      buf(i,m) = grd(ilon,ilat)
    end do
    end do

    do m = 1,np
      if( nij1node(m) < nij1max ) buf(nij1max,m) = UNDEF
    end do

    return
  end subroutine grd_to_buf

  !-------------------------------------------------------------------------------
  ! buffer -> gridded data
  !-------------------------------------------------------------------------------
  subroutine buf_to_grd(np,buf,grd)
    integer,intent(in) :: np
    real(RP),intent(in) :: buf(nij1max,np)
    real(RP),intent(out) :: grd(nlon,nlat)
    integer :: i,j,m,ilon,ilat

    do m = 1,np
    do i = 1,nij1node(m)
      j = m-1 + np * (i-1)
      ilon = mod(j,nlon) + 1
      ilat = (j-ilon+1) / nlon + 1
      grd(ilon,ilat) = buf(i,m)
    end do
    end do

    return
  end subroutine buf_to_grd

  !-------------------------------------------------------------------------------
  ! Calculate 3D height coordinate given the topography height (on scattered grids)
  !-------------------------------------------------------------------------------
  ! [INPUT]
  !   nij         : scattered grid numbers
  !   topo(nij)   : topography height (on scattered grids)
  ! [OUTPUT]
  !   z(nij,nlev) : 3D height coordinate (on scattered grids)
  !-------------------------------------------------------------------------------
  subroutine calc_z_grd(nij, topo, z)
    use scale_atmos_grid_cartesC, only: &
       ATMOS_GRID_CARTESC_CZ, &
       ATMOS_GRID_CARTESC_FZ
    use scale_atmos_grid_cartesC_index, only: &
       KHALO, KS, KE
    implicit none

    integer, intent(in) :: nij
    real(RP), intent(in) :: topo(nij)
    real(RP), intent(out) :: z(nij,nlev)
    real(RP) :: ztop
    integer :: k, i

    ztop = ATMOS_GRID_CARTESC_FZ(KE) - ATMOS_GRID_CARTESC_FZ(KS-1)
    do k = 1, nlev
    do i = 1, nij
      z(i,k) = (ztop - topo(i)) / ztop * ATMOS_GRID_CARTESC_CZ(k+KHALO) + topo(i)
    end do
    end do

    return
  end subroutine calc_z_grd

  subroutine qc_indexing_and_packing(&
       & na, nr, ne, ze, vr, radlon, radlat, radz, & ! input spherical
       & qcflag, input_is_dbz, attenuation, & ! input spherical
       & nlon, nlat, nlev, lon, lat, z, dlon, dlat, dz, & ! input cartesian
       & missing, & ! input param
       & lon0, lat0, &  ! input param
       & nobs_each_elev, packed_grid, packed_data, packed_attn) ! output
    use scale_const, only: &
      RADIUS => CONST_RADIUS, &
      D2R    => CONST_D2R
    implicit none

    integer, intent(in) :: na, nr, ne ! array size of the spherical grid
    real(RP), intent(in) :: ze(na * nr, ne), vr(na * nr, ne) ! main data
    real(RP), intent(in), dimension(na * nr, ne) :: radlon, radlat, radz ! georeference
    real(RP), intent(in) :: qcflag(na * nr, ne) ! additional info
    logical, intent(in) :: input_is_dbz
    real(RP), intent(in) :: lon0, lat0
    real(RP), intent(in), dimension(na * nr, ne) :: attenuation
    integer, intent(in) :: nlon, nlat, nlev ! array size of the cartesian grid
    real(RP), intent(in) :: lon(nlon), lat(nlat), z(nlev)
    real(RP), intent(in) :: dlon, dlat, dz, missing
    integer(8), intent(out) :: nobs_each_elev(ne)
    integer(8), intent(out) :: packed_grid(na * nr * ne) !MAX DATA SIZE
    real(RP), intent(out) :: packed_data(5, na * nr * ne) !MAX DATA SIZE
    logical, intent(out) :: packed_attn(na * nr * ne) !MAX DATA SIZE
    REAL(RP) :: ri, rj, rk, dlon_inv, dlat_inv, dz_inv, qced_ze, qced_vr
    INTEGER :: i, ie
    integer(8) :: idx

    real(RP) :: lon_coef, lat_coef, dist_i, dist_j
    real(RP) :: vr_min_dist_square

    lon_coef = D2R * (cos(lat0 * D2R) * RADIUS)
    lat_coef = D2R * RADIUS
    vr_min_dist_square = vr_min_dist * vr_min_dist

    dlon_inv = 1.0d0 / dlon
    dlat_inv = 1.0d0 / dlat
    dz_inv = 1.0d0 / dz

    nobs_each_elev = 0
    idx = 1
    do ie = 1, ne
    do i = 1, na * nr
       !QC
       qced_ze = ze(i, ie)
       qced_vr = vr(i, ie)

       ! We will work with reflectivity (not in dbz) if we have dbz as input
       ! then transform it
       IF(input_is_dbz .and. (qced_ze .ne. missing)) qced_ze = 10.0d0 ** (qced_ze * 0.1d0)

       ! Missing values not associated with clutter will be asigned a minimum
       ! reflectivity value. 
       if(qced_ze .LE. missing) then
          qced_ze = minz
          qced_vr = missing !added by Otsuka
       end if
       if(qced_ze .LT. minz) then
          qced_ze = minz
       end if
       if(qced_ze .LT. MIN_RADAR_REF_VR) then         
         qced_vr = missing !added by Otsuka  
       end if                                

       ! If dealing with qced real data it will be a good idea to remove all
       ! values detected as non weather echoes in the qc algorithm.
       ! This can also be useful to deal with simulated tophographyc shading
       ! in OSSES.
       ! We need reflectivity to average vr observations. Remove vr
       ! observations where the reflectivity is missing.
       if(USE_QCFLAG .and. (qcflag(i, ie) .GE. 900.0d0)) then
         qced_ze = missing
         qced_vr = missing
       endif

       if(qced_ze == missing) cycle

       ! Get i,j,k very simple approach since we are assuming a regular
       ! lat/lon/z grid.
       ri = (radlon(i, ie) - lon(1)) * dlon_inv
       rj = (radlat(i, ie) - lat(1)) * dlat_inv
       rk = (radz(i, ie)   - z(1)  ) * dz_inv

       if(ri < 0 .or. ri >= nlon .or. rj < 0 .or. rj >= nlat .or. rk < 0 .or. rk >= nlev) cycle

       if (qced_vr > missing) then
         if (radz(i, ie) < RADAR_ZMIN) then
            qced_vr = missing
         else
            dist_i = (radlon(i, ie) - lon0) * lon_coef
            dist_j = (radlat(i, ie) - lat0) * lat_coef
            if (dist_i * dist_i + dist_j * dist_j < vr_min_dist_square) then
              qced_vr = missing
            end if
         end if
       end if

       packed_grid(idx)    = aint(ri) + (aint(rj) + aint(rk) * nlat) * nlon + 1
       packed_data(1, idx) = qced_ze
       packed_data(2, idx) = qced_vr
       packed_data(3, idx) = radlon(i, ie)
       packed_data(4, idx) = radlat(i, ie)
       packed_data(5, idx) = radz(i, ie)
       packed_attn(idx)    = ((.not. USE_ATTENUATION) .or. &
            &                 (attenuation(i, ie) > attenuation_threshold)) !seems wrong but not yet confirmed
       idx = idx + 1
       nobs_each_elev(ie) = nobs_each_elev(ie) + 1

    end do ! i
    end do ! ie
  end subroutine qc_indexing_and_packing

  !-----------------------------------------------------------------------
  ! Convert a raw obsID to a sequential obsID (1 - nid_obs)
  !-----------------------------------------------------------------------
  integer function uid_obs(id_obs)
    implicit none
    integer, intent(in)  :: id_obs
    !---------------------------------------------------------------------

    select case(id_obs)
    case(id_u_obs)
      uid_obs = 1
    case(id_v_obs)
      uid_obs = 2
    case(id_t_obs)
      uid_obs = 3
    case(id_tv_obs)
      uid_obs = 4
    case(id_q_obs)
      uid_obs = 5
    case(id_rh_obs)
      uid_obs = 6
    case(id_ps_obs)
      uid_obs = 7
    case(id_rain_obs)
      uid_obs = 8
    case(id_radar_ref_obs)
      uid_obs = 9
    case(id_radar_ref_zero_obs)
      uid_obs = 10
    case(id_radar_vr_obs)
      uid_obs = 11
    case(id_radar_prh_obs)
      uid_obs = 12
    case(id_h08ir_obs) ! H08
      uid_obs = 13     ! H08
    case(id_tclon_obs)
      uid_obs = 14
    case(id_tclat_obs)
      uid_obs = 15
    case(id_tcmip_obs)
      uid_obs = 16
    case default
      uid_obs = -1     ! error
    end select

    return
  end function uid_obs

  !-----------------------------------------------------------------------
  ! Convert a raw obsID to a sequential obsID for variable localization (1 - nid_obs_verlocal)
  !-----------------------------------------------------------------------
  integer function uid_obs_varlocal(id_obs)
    implicit none
    integer, intent(in)  :: id_obs
    !---------------------------------------------------------------------

    select case(id_obs)
    case(id_u_obs, id_v_obs)
      uid_obs_varlocal = 1
    case(id_t_obs, id_tv_obs)
      uid_obs_varlocal = 2
    case(id_q_obs, id_rh_obs)
      uid_obs_varlocal = 3
    case(id_ps_obs)
      uid_obs_varlocal = 4
    case(id_rain_obs)
      uid_obs_varlocal = 5
    case(id_tclon_obs, id_tclat_obs, id_tcmip_obs)
      uid_obs_varlocal = 6
    case(id_radar_ref_obs, id_radar_ref_zero_obs, id_radar_prh_obs)
      uid_obs_varlocal = 7
    case(id_radar_vr_obs)
      uid_obs_varlocal = 8
    case(id_h08ir_obs)      ! H08
      uid_obs_varlocal = 9  ! H08
    case default
      uid_obs_varlocal = -1 ! error
    end select
  end function uid_obs_varlocal

  function binary_search_i8(n, ary, val)
    integer(8) binary_search_i8
    integer(8), intent(in) :: n
    integer(8), intent(in) :: ary(n), val
    integer(8) pivot, nmax, nmin
    nmin = 1
    nmax = n

    if(ary(1) < ary(n)) then
       do while(nmax > nmin + 1)
          pivot = nmin + (nmax - nmin) / 2
          if(val == ary(pivot)) then
             nmin = pivot
             exit
          else if(val < ary(pivot)) then
             nmax = pivot
          else
             nmin = pivot
          end if
       end do
    else
       do while(nmax > nmin + 1)
          pivot = nmin + (nmax - nmin) / 2
          if(val == ary(pivot)) then
             nmin = pivot
             exit
          else if(val > ary(pivot)) then
             nmax = pivot
          else
             nmin = pivot
          end if
       end do
    end if
    binary_search_i8 = nmin
  end function binary_search_i8

  !-------------------------------------------------------------------------------
  ! Convert 1D rank of process to 2D rank
  !-------------------------------------------------------------------------------
  subroutine rank_1d_2d(rank, rank_i, rank_j)
    use scale_prc_cartesC, only: PRC_2Drank
    implicit none
    integer, intent(in)  :: rank
    integer, intent(out) :: rank_i, rank_j
    !---------------------------------------------------------------------

    rank_i = PRC_2Drank(rank,1)
    rank_j = PRC_2Drank(rank,2)

    return  
  end subroutine rank_1d_2d

  !-------------------------------------------------------------------------------
  ! Convert 2D rank of process to 1D rank
  !-------------------------------------------------------------------------------
  subroutine rank_2d_1d(rank_i, rank_j, rank)
    implicit none
    integer, intent(in)  :: rank_i, rank_j
    integer, intent(out) :: rank

    rank = rank_j * nproc_x + rank_i

    return  
  end subroutine rank_2d_1d

  !-------------------------------------------------------------------------------
  ! Given <real> global grid coordinates (i,j), return the 1D rank of process 
  ! * HALO grids are used
  !-------------------------------------------------------------------------------
  ! [INPUT]
  !   ig, jg : global grid coordinates
  ! [OUTPUT]
  !   rank   : the 1D rank of process where the grid resides;
  !            * return -1 if the grid is outside of the global domain
  !-------------------------------------------------------------------------------
  subroutine rij_rank(ig, jg, rank)
    implicit none
    real(RP), intent(in)  :: ig
    real(RP), intent(in)  :: jg
    integer,  intent(out) :: rank
    integer :: rank_i, rank_j

    if (ig < real(start_x,kind=RP) .or. ig > real(nlon*nproc_x+xhalo,kind=RP) .or. &
        jg < real(start_y,kind=RP) .or. jg > real(nlat*nproc_y+yhalo,kind=RP)) then
      !!! exception : 2-D ideal case
      if (.not.(start_x == end_x .and. (ig < real(start_x,kind=RP) .or. ig > real(nlong+xhalo,kind=RP) ))) then
        rank = -1
        return
      end if
      rank = -1
      return
    end if

    rank_i = ceiling((ig-real(xhalo,kind=RP)-0.5_RP) / real(nlon,kind=RP)) - 1
    rank_j = ceiling((jg-real(yhalo,kind=RP)-0.5_RP) / real(nlat,kind=RP)) - 1
    call rank_2d_1d(rank_i, rank_j, rank)

    return
  end subroutine rij_rank

  !-------------------------------------------------------------------------------
  ! Convert <real> global grid coordinates (i,j) to local given the 1D rank of process
  !-------------------------------------------------------------------------------
  subroutine rij_g2l(rank, ig, jg, il, jl)
    implicit none
    integer,  intent(in)  :: rank
    real(RP), intent(in)  :: ig
    real(RP), intent(in)  :: jg
    real(RP), intent(out) :: il
    real(RP), intent(out) :: jl

    integer :: rank_i, rank_j
    !---------------------------------------------------------------------

    call rank_1d_2d( rank, rank_i, rank_j )
    il = ig - real(rank_i * nlon, kind=RP)
    jl = jg - real(rank_j * nlat, kind=RP)

    return  
  end subroutine rij_g2l

  !-----------------------------------------------------------------------
  ! Convert grid (i,j) values to obsgrid (ogi, ogj) sorting mesh
  !-----------------------------------------------------------------------
  subroutine ij_obsgrd( ctype, ri, rj, ogi, ogj )
    implicit none
    integer, intent(in) :: ctype
    real(RP), intent(in) :: ri, rj
    integer, intent(out) :: ogi, ogj
    real(RP) :: ril, rjl

    call rij_g2l( RANK_LCL, ri, rj, ril, rjl )
    ogi = ceiling( ( ril - real( xhalo, kind=RP ) - 0.5 ) * real( obsgrd(ctype)%ngrd_i, kind=RP ) / real( nlon, kind=RP ) )
    ogj = ceiling( ( rjl - real( yhalo, kind=RP ) - 0.5 ) * real( obsgrd(ctype)%ngrd_j, kind=RP ) / real( nlat, kind=RP ) )

    return
  end subroutine ij_obsgrd

  !-----------------------------------------------------------------------
  ! Convert grid (i,j) values to obsgrid (ogi, ogj) sorting mesh in the extended subdomain
  !-----------------------------------------------------------------------
  subroutine ij_obsgrd_ext( ctype, ri, rj, ogi, ogj )
    implicit none
    integer,  intent(in)  :: ctype
    real(RP), intent(in)  :: ri, rj
    integer,  intent(out) :: ogi, ogj

    real(RP) :: ril, rjl
    !---------------------------------------------------------------------

    call rij_g2l( RANK_LCL, ri, rj, ril, rjl )
    ogi = ceiling( ( ril - real( xhalo, kind=RP ) - 0.5 ) * real( obsgrd(ctype)%ngrd_i, kind=RP ) / real( nlon, kind=RP ) ) &
        + obsgrd(ctype)%ngrdsch_i
    ogj = ceiling( ( rjl - real( yhalo, kind=RP ) - 0.5 ) * real( obsgrd(ctype)%ngrd_j, kind=RP ) / real( nlat, kind=RP ) ) &
        + obsgrd(ctype)%ngrdsch_j

    return
  end subroutine ij_obsgrd_ext

  !-----------------------------------------------------------------------
  ! Choose observations in a rectangle using the bucket sort results
  !-----------------------------------------------------------------------
  subroutine obs_choose(ctype, proc, imin, imax, jmin, jmax, nn, nobs_use)
    implicit none
    integer, intent(in) :: ctype
    integer, intent(in) :: proc
    integer, intent(in) :: imin, imax, jmin, jmax
    integer, intent(inout) :: nn
    integer, intent(inout), optional :: nobs_use(:)
    integer :: n, j

    if (imin > imax .or. jmin > jmax) return
    if (obsgrd(ctype)%tot(proc) == 0) return

    do j = jmin, jmax
      if (present(nobs_use)) then
        do n = obsgrd(ctype)%ac(imin-1,j,proc)+1, obsgrd(ctype)%ac(imax,j,proc)
          nn = nn + 1
          nobs_use(nn) = n
        end do
      else
        nn = nn + obsgrd(ctype)%ac(imax,j,proc) - obsgrd(ctype)%ac(imin-1,j,proc)
      end if
    end do

    return
  end subroutine obs_choose

  !-----------------------------------------------------------------------
  ! Choose observations in a rectangle using the bucket sort results in the extended subdomain
  !-----------------------------------------------------------------------
  subroutine obs_choose_ext(ctype, imin, imax, jmin, jmax, nn, nobs_use)
    implicit none
    integer, intent(in)    :: ctype
    integer, intent(in)    :: imin, imax, jmin, jmax
    integer, intent(inout) :: nn
    integer, intent(out), optional :: nobs_use(:)

    integer :: n, j
    !---------------------------------------------------------------------

    if (imin > imax .or. jmin > jmax) return
    if (obsgrd(ctype)%tot_ext == 0) return

    do j = jmin, jmax
      if (present(nobs_use)) then
        do n = obsgrd(ctype)%ac_ext(imin-1,j)+1, obsgrd(ctype)%ac_ext(imax,j)
          nn = nn + 1
          nobs_use(nn) = n
        end do
      else
        nn = nn + obsgrd(ctype)%ac_ext(imax,j) - obsgrd(ctype)%ac_ext(imin-1,j)
      end if
    end do

    return
  end subroutine obs_choose_ext

  subroutine obs_info_allocate(obs, extended)
    implicit none
    type(obs_info),intent(inout) :: obs
    logical, optional, intent(in) :: extended

    call obs_info_deallocate(obs)

    allocate( obs%elm (obs%nobs) )
    allocate( obs%lon (obs%nobs) )
    allocate( obs%lat (obs%nobs) )
    allocate( obs%lev (obs%nobs) )
    allocate( obs%dat (obs%nobs) )
    allocate( obs%err (obs%nobs) )
    allocate( obs%typ (obs%nobs) )
    allocate( obs%dif (obs%nobs) )

    obs%elm = 0
    obs%lon = 0.0d0
    obs%lat = 0.0d0
    obs%lev = 0.0d0
    obs%dat = 0.0d0
    obs%err = 0.0d0
    obs%typ = 0
    obs%dif = 0.0d0

    if (present(extended)) then
      if (extended) then
        allocate( obs%ri (obs%nobs) )
        allocate( obs%rj (obs%nobs) )
        allocate( obs%rank (obs%nobs) )

        obs%ri = 0.0d0
        obs%rj = 0.0d0
        obs%rank = -1
      end if
    end if

    return
  end subroutine obs_info_allocate

  subroutine obs_info_deallocate(obs)
    implicit none
    type(obs_info),intent(inout) :: obs

    if( allocated(obs%elm) ) deallocate(obs%elm)
    if( allocated(obs%lon) ) deallocate(obs%lon)
    if( allocated(obs%lat) ) deallocate(obs%lat)
    if( allocated(obs%lev) ) deallocate(obs%lev)
    if( allocated(obs%dat) ) deallocate(obs%dat)
    if( allocated(obs%err) ) deallocate(obs%err)
    if( allocated(obs%typ) ) deallocate(obs%typ)
    if( allocated(obs%dif) ) deallocate(obs%dif)
    if( allocated(obs%ri ) ) deallocate(obs%ri)
    if( allocated(obs%rj ) ) deallocate(obs%rj)
    if( allocated(obs%rank)) deallocate(obs%rank)

    return
  end subroutine obs_info_deallocate

  subroutine obs_da_value_allocate( obsda, member )
    implicit none
    type(obs_da_value), intent(inout) :: obsda
    integer, intent(in) :: member

    call obs_da_value_deallocate( obsda )

    allocate( obsda%set (obsda%nobs) )
    allocate( obsda%idx (obsda%nobs) )
    allocate( obsda%key (obsda%nobs) )
    allocate( obsda%val (obsda%nobs) )
    allocate( obsda%qc  (obsda%nobs) )

    allocate( obsda%tm (obsda%nobs) )
    allocate( obsda%pm (obsda%nobs) )
    allocate( obsda%qv (obsda%nobs) )

    obsda%nobs_in_key = 0
    obsda%idx = 0
    obsda%key = 0
    obsda%val = 0.0d0
    obsda%qc = 0

    obsda%tm = 0.0d0
    obsda%pm = 0.0d0
    obsda%qv = 0.0d0

    if (member > 0) then
      allocate( obsda%ensval(member,obsda%nobs) )
      obsda%ensval = 0.0d0

      allocate( obsda%eqv(member,obsda%nobs) )
      obsda%eqv = 0.0d0
    end if

    return
  end subroutine obs_da_value_allocate

  subroutine obs_da_value_deallocate( obsda )
    implicit none
    type(obs_da_value), intent(inout) :: obsda

    obsda%nobs_in_key = 0

    if( allocated(obsda%set   ) ) deallocate(obsda%set   )
    if( allocated(obsda%idx   ) ) deallocate(obsda%idx   )
    if( allocated(obsda%key   ) ) deallocate(obsda%key   )
    if( allocated(obsda%val   ) ) deallocate(obsda%val   )
    if( allocated(obsda%ensval) ) deallocate(obsda%ensval)
    if( allocated(obsda%qc    ) ) deallocate(obsda%qc    )

    if( allocated(obsda%tm    ) ) deallocate(obsda%tm)
    if( allocated(obsda%pm    ) ) deallocate(obsda%pm)
    if( allocated(obsda%eqv   ) ) deallocate(obsda%eqv)
    if( allocated(obsda%qv    ) ) deallocate(obsda%qv)

    return
  end subroutine obs_da_value_deallocate

  subroutine obs_da_value_allreduce( obsda )
    implicit none
    type(obs_da_value), intent(inout) :: obsda

    real(RP), allocatable :: ensval_bufs(:,:)
    real(RP), allocatable :: ensval_bufr(:,:)
    real(RP), allocatable :: ensval_bufs2(:,:)
    real(RP), allocatable :: ensval_bufr2(:,:)
    integer :: cnts
    integer :: cntr(NPRC_ENS)
    integer :: dspr(NPRC_ENS)
    integer :: current_shape(2)
    integer :: ie, it, im, imb, ierr

    if( obsda%nobs <= 0 ) then
      return
    end if

    ! variables with an ensemble dimension
    cntr(:) = 0
    do ie = 1, NPRC_ENS
      cntr(ie) = cntr(ie) + 1
    end do
    allocate( ensval_bufs(obsda%nobs, cntr(RANK_ENS+1)) )
    allocate( ensval_bufr(obsda%nobs, NPRC_ENS) )
    allocate( ensval_bufs2(obsda%nobs, cntr(RANK_ENS+1)) )
    allocate( ensval_bufr2(obsda%nobs, NPRC_ENS) )

    do im = 1, cntr(RANK_ENS+1)
      ensval_bufs(:,im) = obsda%ensval(im,:)
      ensval_bufs2(:,im) = obsda%eqv(im,:)
    end do

    cntr(:) = cntr(:) * obsda%nobs
    cnts = cntr(RANK_ENS+1)
    dspr(1) = 0
    do ie = 2, NPRC_ENS
      dspr(ie) = dspr(ie-1) + cntr(ie-1)
    end do

    call MPI_ALLGATHERV( ensval_bufs,  cnts, datatype, ensval_bufr,  cntr, dspr, datatype, COMM_ENS, ierr )
    call MPI_ALLGATHERV( ensval_bufs2, cnts, datatype, ensval_bufr2, cntr, dspr, datatype, COMM_ENS, ierr )

    current_shape = shape(obsda%ensval)
    if (current_shape(1) < NPRC_ENS) then
      deallocate (obsda%ensval)
      allocate (obsda%ensval(NPRC_ENS, obsda%nobs))
      deallocate (obsda%eqv)
      allocate (obsda%eqv(NPRC_ENS, obsda%nobs))
    end if

    do ie = 1, NPRC_ENS
      obsda%ensval(ie,:) = ensval_bufr(:,ie)
      obsda%eqv(ie,:) = ensval_bufr2(:,ie)
    end do
    deallocate(ensval_bufs, ensval_bufr)
    deallocate(ensval_bufs2, ensval_bufr2)

    ! variables without an ensemble dimension
    if( NPRC_ENS > 1 ) then
      call MPI_ALLREDUCE( MPI_IN_PLACE, obsda%qc(:), obsda%nobs, MPI_INTEGER, MPI_MAX, COMM_ENS, ierr )
      call MPI_ALLREDUCE( MPI_IN_PLACE, obsda%tm(:), obsda%nobs, datatype,    MPI_SUM, COMM_ENS, ierr )
      call MPI_ALLREDUCE( MPI_IN_PLACE, obsda%pm(:), obsda%nobs, datatype,    MPI_SUM, COMM_ENS, ierr )
    end if
    obsda%tm = obsda%tm / real( NPRC_ENS, kind=RP )
    obsda%pm = obsda%pm / real( NPRC_ENS, kind=RP )

    return
  end subroutine obs_da_value_allreduce

  subroutine obs_da_value_partial_reduce_iter(obsda, iter, nstart, nobs, ensval, qc, eqv, tm, pm )
    implicit none
    type(obs_da_value), intent(inout) :: obsda
    integer,  intent(in) :: iter
    integer,  intent(in) :: nstart
    integer,  intent(in) :: nobs
    real(RP), intent(in) :: ensval(nobs)
    integer,  intent(in) :: qc(nobs)
    real(RP), intent(in) :: eqv(nobs)
    real(RP), intent(in) :: tm(nobs)
    real(RP), intent(in) :: pm(nobs)

    integer :: nend

    if (nobs <= 0) then
      return
    end if
    nend = nstart + nobs - 1

    ! variables with an ensemble dimension
    obsda%ensval(iter,nstart:nend) = ensval
    obsda%eqv(iter,nstart:nend) = eqv

    ! variables without an ensemble dimension
    obsda%qc(nstart:nend) = max(obsda%qc(nstart:nend), qc)
    ! only consider tm & pm from members, not from the mean
    obsda%tm(nstart:nend) = obsda%tm(nstart:nend) + tm
    obsda%pm(nstart:nend) = obsda%pm(nstart:nend) + pm

    return
  end subroutine obs_da_value_partial_reduce_iter

  !-------------------------------------------------------------------------------
  ! Read ensemble additive inflation parameter and distribute to processes
  !-------------------------------------------------------------------------------
  subroutine read_ens_mpi_addiinfl(v3d, v2d)
    implicit none
    real(RP), intent(out) :: v3d(:,:,:,:)
    real(RP), intent(out) :: v2d(:,:,:)

    character(len=H_LONG) :: filename
    real(RP) :: v3dg(nlev,nlon,nlat,nv3d)
    real(RP) :: v2dg(nlon,nlat,nv2d)
    integer :: it, im, mstart, mend
    !---------------------------------------------------------------------

    ! [TODO]: restructure this subroutine
    v3d = 0.0_RP
    v2d = 0.0_RP
    !
    !do it = 1, nitmax
    !  im = myrank_to_mem(it)

    !  ! Note: read all members
    !  if (im >= 1 .and. im <= MEMBER) then
    !    filename = INFL_ADD_IN_BASENAME
    !    call filename_replace_mem(filename, im)
    !    call read_restart(filename, v3dg, v2dg)
    !  end if

    !  mstart = 1 + (it-1)*nprocs_e
    !  mend = min(it*nprocs_e, MEMBER)
    !  if (mstart <= mend) then
    !    call scatter_grd_mpi_alltoall(mstart, mend, v3dg, v2dg, v3d, v2d)
    !  end if
    !end do

    return
  end subroutine read_ens_mpi_addiinfl

  !-----------------------------------------------------------------------
  ! Monitor observation departure by giving the v3dg,v2dg data
  !-----------------------------------------------------------------------
  subroutine monit_obs( OBS_IN_NUM, OBS_IN_FORMAT, v3dg, v2dg, nobs, bias, rmse, monit_type, use_key, step )
    use scale_const, only: &
      UNDEF => CONST_UNDEF
    use scale_prc, only: &
      PRC_myrank
    implicit none

    integer, intent(in) :: OBS_IN_NUM
    character(len=H_LONG), intent(in) :: OBS_IN_FORMAT(:)

    real(RP), intent(in)  :: v3dg(nlev,nlon,nlat,nv3d)
    real(RP), intent(in)  :: v2dg(nlon,nlat,nv2d)
    integer,  intent(out) :: nobs(nid_obs)
    real(RP), intent(out) :: bias(nid_obs)
    real(RP), intent(out) :: rmse(nid_obs)
    logical,  intent(out) :: monit_type(nid_obs)
    logical,  intent(in)  :: use_key
    integer,  intent(in)  :: step

    real(RP) :: v3dgh(nlevh,nlonh,nlath,nv3dd)
    real(RP) :: v2dgh(nlonh,nlath,nv2dd)
    integer  :: nnobs
    integer  :: n,nn
    integer  :: iset,iidx
    real(RP) :: ril,rjl,rk,rkz

    real(RP), allocatable :: oelm(:)
    real(RP), allocatable :: ohx(:)
    integer,  allocatable :: oqc(:)

    call state_to_history(v3dg, v2dg, v3dgh, v2dgh)

    if (use_key) then
      nnobs = obsda_sort%nobs_in_key
    else
      nnobs = obsda_sort%nobs
    end if

    allocate (oelm(nnobs))
    allocate (ohx(nnobs))
    allocate (oqc(nnobs))

    if (step == 1) then
      obsdep_nobs = nnobs
      allocate (obsdep_set(obsdep_nobs))
      allocate (obsdep_idx(obsdep_nobs))
      allocate (obsdep_qc (obsdep_nobs))
      allocate (obsdep_omb(obsdep_nobs))
      allocate (obsdep_oma(obsdep_nobs))
    end if

    oqc = -1

    do n = 1, nnobs

      if (use_key) then
        nn = obsda_sort%key(n)
      else
        nn = n
      end if

      iset = obsda_sort%set(nn)
      iidx = obsda_sort%idx(nn)

      if (step == 1) then
        obsdep_set(n) = iset
        obsdep_idx(n) = iidx
      end if

      oelm(n) = obs(iset)%elm(iidx)
      call rij_g2l(PRC_myrank, obs(iset)%ri(iidx), obs(iset)%rj(iidx), ril, rjl)

      if (DEPARTURE_STAT_T_RANGE <= 0.0d0 .or. &
          abs(obs(iset)%dif(iidx)) <= DEPARTURE_STAT_T_RANGE) then

        oqc(n) = iqc_otype

        select case (OBS_IN_FORMAT(iset))
        !=========================================================================
        case ( 'PREPBUFR' )
        !-------------------------------------------------------------------------
          call phys2ijk(v3dgh(:,:,:,iv3dd_p),obs(iset)%elm(iidx), &
                        ril,rjl,obs(iset)%lev(iidx),rk,oqc(n))
          if (oqc(n) == iqc_good) then
            call Trans_XtoY(obs(iset)%elm(iidx),ril,rjl,rk, &
                            obs(iset)%lon(iidx),obs(iset)%lat(iidx), &
                            v3dgh,v2dgh,ohx(n),oqc(n),stggrd=1)
          end if
        !=========================================================================
        case ( 'RADAR', 'PAWR_TOSHIBA', 'MP_PAWR_TOSHIBA', 'PAWR_JRC', 'HIMAWARI8' )
        !-------------------------------------------------------------------------
          if (DEPARTURE_STAT_RADAR) then
            call phys2ijkz(v3dgh(:,:,:,iv3dd_hgt),ril,rjl,obs(iset)%lev(iidx),rkz,oqc(n))
            if (oqc(n) == iqc_good) then
              call Trans_XtoY_radar(obs(iset)%elm(iidx),obs(iset)%meta(1), &
                                    obs(iset)%meta(2),obs(iset)%meta(3),ril,rjl,rkz, &
                                    obs(iset)%lon(iidx),obs(iset)%lat(iidx), &
                                    obs(iset)%lev(iidx),v3dgh,v2dgh,ohx(n),oqc(n),stggrd=1)
              if (oqc(n) == iqc_ref_low) oqc(n) = iqc_good ! when process the observation operator, we don't care if reflectivity is 350
            end if
          end if
        !=========================================================================
        end select

        if (oqc(n) == iqc_good) then
          ohx(n) = obs(iset)%dat(iidx) - ohx(n)
        else
          ohx(n) = undef
        end if

        if (step == 1) then
          obsdep_qc(n) = oqc(n)
          obsdep_omb(n) = ohx(n)
        else if (step == 2) then
          if (obsdep_qc(n) == iqc_good) then ! Use the QC value of y_a only if the QC of y_b is good
            obsdep_qc(n) = oqc(n)            !
          end if                             !
          obsdep_oma(n) = ohx(n)
        end if

      end if ! [ DEPARTURE_STAT_T_RANGE <= 0.0d0 .or. &
             !   abs(obs(iset)%dif(iidx)) <= DEPARTURE_STAT_T_RANGE ]

    end do ! [ n = 1, nnobs ]


    call monit_dep(nnobs,oelm,ohx,oqc,nobs,bias,rmse)

    monit_type = .false.
    monit_type(uid_obs(id_u_obs)) = .true.
    monit_type(uid_obs(id_v_obs)) = .true.
    monit_type(uid_obs(id_t_obs)) = .true.
    monit_type(uid_obs(id_tv_obs)) = .true.
    monit_type(uid_obs(id_q_obs)) = .true.
    monit_type(uid_obs(id_rh_obs)) = .true.
    monit_type(uid_obs(id_ps_obs)) = .true.
    if (DEPARTURE_STAT_RADAR) then
      monit_type(uid_obs(id_radar_ref_obs)) = .true.
      monit_type(uid_obs(id_radar_ref_zero_obs)) = .true.
      monit_type(uid_obs(id_radar_vr_obs)) = .true.
    end if

    deallocate (oelm)
    deallocate (ohx)
    deallocate (oqc)

    return
  end subroutine monit_obs

  !-------------------------------------------------------------------------------
  ! MPI driver for monitoring observation departure statistics
  !-------------------------------------------------------------------------------
  subroutine monit_obs_mpi( OBS_IN_NUM, OBS_IN_FORMAT, v3dg, v2dg, monit_step, timelabel )
    use scale_const, only: &
      UNDEF => CONST_UNDEF
    implicit none

    integer, intent(in) :: OBS_IN_NUM
    character(len=H_LONG), intent(in) :: OBS_IN_FORMAT(:)

    real(RP), intent(in) :: v3dg(nlev,nlon,nlat,nv3d)
    real(RP), intent(in) :: v2dg(nlon,nlat,nv2d)
    integer,  intent(in) :: monit_step
    character(15), intent(in), optional :: timelabel

    integer  :: nobs(nid_obs)
    integer  :: nobs_g(nid_obs)
    real(RP) :: bias(nid_obs)
    real(RP) :: bias_g(nid_obs)
    real(RP) :: rmse(nid_obs)
    real(RP) :: rmse_g(nid_obs)
    logical  :: monit_type(nid_obs)
    integer  :: obsdep_g_nobs
    integer,  allocatable :: obsdep_g_set(:)
    integer,  allocatable :: obsdep_g_idx(:)
    integer,  allocatable :: obsdep_g_qc(:)
    real(RP), allocatable :: obsdep_g_omb(:)
    real(RP), allocatable :: obsdep_g_oma(:)
    integer  :: cnts
    integer  :: cntr(NPRC_LCL)
    integer  :: dspr(NPRC_LCL)
    integer  :: i, ip, ierr

    if( RANK_ENS == 0 ) then
      call monit_obs( OBS_IN_NUM, OBS_IN_FORMAT, v3dg, v2dg, nobs, bias, rmse, monit_type, .true., monit_step )

      do i = 1, nid_obs
        if (monit_type(i)) then
          nobs_g(i) = nobs(i)
          if (nobs(i) == 0) then
            bias_g(i) = 0.0d0
            rmse_g(i) = 0.0d0
          else
            bias_g(i) = bias(i) * real( nobs(i), kind=RP )
            rmse_g(i) = rmse(i) * rmse(i) * real( nobs(i), kind=RP )
          end if
        end if
      end do

      if( NPRC_LCL > 1 ) then
        call MPI_ALLREDUCE(MPI_IN_PLACE, nobs_g, nid_obs, MPI_INTEGER, MPI_SUM, COMM_LCL, ierr)
        call MPI_ALLREDUCE(MPI_IN_PLACE, bias_g, nid_obs, datatype,    MPI_SUM, COMM_LCL, ierr)
        call MPI_ALLREDUCE(MPI_IN_PLACE, rmse_g, nid_obs, datatype,    MPI_SUM, COMM_LCL, ierr)
      end if

      do i = 1, nid_obs
        if (monit_type(i)) then
          if (nobs_g(i) == 0) then
            bias_g(i) = undef
            rmse_g(i) = undef
          else
            bias_g(i) = bias_g(i) / REAL(nobs_g(i),kind=RP)
            rmse_g(i) = sqrt(rmse_g(i) / REAL(nobs_g(i),kind=RP))
          end if
        else
          nobs_g(i) = -1
          bias_g(i) = undef
          rmse_g(i) = undef
        end if
      end do

      if (monit_step == 2) then
        cnts = obsdep_nobs
        cntr = 0
        cntr(RANK_LCL+1) = cnts
        call MPI_ALLREDUCE( MPI_IN_PLACE, cntr, NPRC_LCL, MPI_INTEGER, MPI_SUM, COMM_LCL, ierr )
        dspr = 0
        do ip = 1, NPRC_LCL-1
          dspr(ip+1) = dspr(ip) + cntr(ip)
        end do

        obsdep_g_nobs = dspr(NPRC_LCL) + cntr(NPRC_LCL)
        allocate (obsdep_g_set(obsdep_g_nobs))
        allocate (obsdep_g_idx(obsdep_g_nobs))
        allocate (obsdep_g_qc (obsdep_g_nobs))
        allocate (obsdep_g_omb(obsdep_g_nobs))
        allocate (obsdep_g_oma(obsdep_g_nobs))

        if (obsdep_g_nobs > 0) then
          call MPI_GATHERV(obsdep_set, cnts, MPI_INTEGER, obsdep_g_set, cntr, dspr, MPI_INTEGER, 0, COMM_LCL, ierr)
          call MPI_GATHERV(obsdep_idx, cnts, MPI_INTEGER, obsdep_g_idx, cntr, dspr, MPI_INTEGER, 0, COMM_LCL, ierr)
          call MPI_GATHERV(obsdep_qc,  cnts, MPI_INTEGER, obsdep_g_qc,  cntr, dspr, MPI_INTEGER, 0, COMM_LCL, ierr)
          call MPI_GATHERV(obsdep_omb, cnts, datatype,    obsdep_g_omb, cntr, dspr, datatype,    0, COMM_LCL, ierr)
          call MPI_GATHERV(obsdep_oma, cnts, datatype,    obsdep_g_oma, cntr, dspr, datatype,    0, COMM_LCL, ierr)
        end if

        deallocate (obsdep_g_set)
        deallocate (obsdep_g_idx)
        deallocate (obsdep_g_qc )
        deallocate (obsdep_g_omb)
        deallocate (obsdep_g_oma)

      end if ! [ OBSDEP_OUT .and. monit_step == 2 ]

      if (monit_step == 2) then
        deallocate (obsdep_set)
        deallocate (obsdep_idx)
        deallocate (obsdep_qc )
        deallocate (obsdep_omb)
        deallocate (obsdep_oma)
      end if
    end if ! [ myrank_e == mmean_rank_e ]

    call MPI_BCAST(nobs,       nid_obs, MPI_INTEGER, 0, COMM_ENS, ierr)
    call MPI_BCAST(bias,       nid_obs, datatype,    0, COMM_ENS, ierr)
    call MPI_BCAST(rmse,       nid_obs, datatype,    0, COMM_ENS, ierr)
    call MPI_BCAST(nobs_g,     nid_obs, MPI_INTEGER, 0, COMM_ENS, ierr)
    call MPI_BCAST(bias_g,     nid_obs, datatype,    0, COMM_ENS, ierr)
    call MPI_BCAST(rmse_g,     nid_obs, datatype,    0, COMM_ENS, ierr)
    call MPI_BCAST(monit_type, nid_obs, MPI_LOGICAL, 0, COMM_ENS, ierr)

    if( monit_step == 1 ) then
      LOG_INFO("LETKF_debug",'(1x,A)') 'OBSERVATIONAL DEPARTURE STATISTICS [GUESS] (IN THIS SUBDOMAIN):'
    else if (monit_step == 2) then
      LOG_INFO("LETKF_debug",'(1x,A)') 'OBSERVATIONAL DEPARTURE STATISTICS [ANALYSIS] (IN THIS SUBDOMAIN):'
    end if
    call monit_print(nobs, bias, rmse, monit_type)

    if (monit_step == 1) then
      LOG_INFO("LETKF_debug",'(1x,A)') 'OBSERVATIONAL DEPARTURE STATISTICS [GUESS] (GLOBAL):'
    else if (monit_step == 2) then
      LOG_INFO("LETKF_debug",'(1x,A)') 'OBSERVATIONAL DEPARTURE STATISTICS [ANALYSIS] (GLOBAL):'
    end if
    call monit_print(nobs_g, bias_g, rmse_g, monit_type)

    return
  end subroutine monit_obs_mpi

  !
  ! monit_obs is ported into obs/obs_tools.f90
  !
  SUBROUTINE monit_dep(nn,elm,dep,qc,nobs,bias,rmse)
    use scale_const, only: &
      UNDEF => CONST_UNDEF
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: nn
    REAL(RP),INTENT(IN) :: elm(nn)
    REAL(RP),INTENT(IN) :: dep(nn)
    INTEGER,INTENT(IN) :: qc(nn)
    INTEGER,INTENT(OUT) :: nobs(nid_obs)
    REAL(RP),INTENT(OUT) :: bias(nid_obs)
    REAL(RP),INTENT(OUT) :: rmse(nid_obs)
    INTEGER :: n,i,ielm

    nobs = 0
    bias = 0.0d0
    rmse = 0.0d0

    DO n=1,nn
      IF(qc(n) /= iqc_good) CYCLE

      ielm = NINT(elm(n))
      if (ielm == id_tv_obs) then ! compute Tv as T
        ielm = id_t_obs
      end if
      if (ielm == id_radar_ref_zero_obs) then ! compute RE0 as REF
        ielm = id_radar_ref_obs
      end if

      i = uid_obs(ielm)
      nobs(i) = nobs(i) + 1
      bias(i) = bias(i) + dep(n)
      rmse(i) = rmse(i) + dep(n)**2
    END DO

    DO i = 1, nid_obs
      IF(nobs(i) == 0) THEN
        bias(i) = undef
        rmse(i) = undef
      ELSE
        bias(i) = bias(i) / REAL(nobs(i),kind=RP)
        rmse(i) = SQRT(rmse(i) / REAL(nobs(i),kind=RP))
      END IF
    END DO

    RETURN
  END SUBROUTINE monit_dep

  !-----------------------------------------------------------------------
  ! Monitor departure
  !-----------------------------------------------------------------------
  SUBROUTINE monit_print(nobs,bias,rmse,monit_type)
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: nobs(nid_obs)
    REAL(RP),INTENT(IN) :: bias(nid_obs)
    REAL(RP),INTENT(IN) :: rmse(nid_obs)
    LOGICAL,INTENT(IN),OPTIONAL :: monit_type(nid_obs)

    character(12) :: var_show(nid_obs)
    character(12) :: nobs_show(nid_obs)
    character(12) :: bias_show(nid_obs)
    character(12) :: rmse_show(nid_obs)

    integer :: i, n
    character(4) :: nstr

    logical :: monit_type_(nid_obs)

    monit_type_ = .true.
    if (present(monit_type)) monit_type_ = monit_type

    n = 0
    do i = 1, nid_obs
      if (monit_type_(i) .and. i /= uid_obs(id_tv_obs) .and. i /= uid_obs(id_radar_ref_zero_obs)) then
        n = n + 1
        write(var_show(n),'(A12)') obelmlist(i)
        write(nobs_show(n),'(I12)') nobs(i)
        if (nobs(i) > 0) then
          write(bias_show(n),'(ES12.3)') bias(i)
          write(rmse_show(n),'(ES12.3)') rmse(i)
        else
          write(bias_show(n),'(A12)') 'N/A'
          write(rmse_show(n),'(A12)') 'N/A'
        end if
      end if
    end do
    write(nstr, '(I4)') n

    LOG_INFO("LETKF_debug",'(1x,A,' // trim(nstr) // "('============'))") '======'
    LOG_INFO("LETKF_debug",'(1x,A,' // trim(nstr) // 'A)'               ) '      ', var_show(1:n)
    LOG_INFO("LETKF_debug",'(1x,A,' // trim(nstr) // "('------------'))") '------'
    LOG_INFO("LETKF_debug",'(1x,A,' // trim(nstr) // 'A)'               ) 'BIAS  ', bias_show(1:n)
    LOG_INFO("LETKF_debug",'(1x,A,' // trim(nstr) // 'A)'               ) 'RMSE  ', rmse_show(1:n)
    LOG_INFO("LETKF_debug",'(1x,A,' // trim(nstr) // 'A)'               ) 'NUMBER', nobs_show(1:n)
    LOG_INFO("LETKF_debug",'(1x,A,' // trim(nstr) // "('============'))") '======'

    RETURN
  END SUBROUTINE monit_print

  !-------------------------------------------------------------------------------
  ! Transform the LETKF state variables to the variables in SCALE history files
  ! (with HALO), so that they can be used for observation operator calculation
  !-------------------------------------------------------------------------------
  ! [INPUT]
  !   v3dg, v2dg   : 3D, 2D state variables
  !   topo         : topography
  ! [OUTPUT]
  !   v3dgh, v2dgh : 3D, 2D SCALE history variables
  !-------------------------------------------------------------------------------
  subroutine state_to_history(v3dg, v2dg, v3dgh, v2dgh)
    use scale_comm_cartesC, only: &
        COMM_vars8, &
        COMM_wait
    use scale_atmos_grid_cartesC_metric, only: &
        ROTC => ATMOS_GRID_CARTESC_METRIC_ROTC
    use scale_atmos_grid_cartesC_real, only: &
        REAL_CZ => ATMOS_GRID_CARTESC_REAL_CZ
    use scale_const, only: &
        UNDEF => CONST_UNDEF, &
        Rdry  => CONST_Rdry,  &
        Rvap  => CONST_Rvap
    use scale_atmos_saturation, only: &
        ATMOS_SATURATION_psat_all
    implicit none

    real(RP), intent(in) :: v3dg(nlev,nlon,nlat,nv3d)
    real(RP), intent(in) :: v2dg(nlon,nlat,nv2d)
    real(RP), intent(out) :: v3dgh(nlevh,nlonh,nlath,nv3dd)
    real(RP), intent(out) :: v2dgh(nlonh,nlath,nv2dd)
    real(RP) :: v3dgh_RP(nlevh,nlonh,nlath,nv3dd)
    real(RP) :: v2dgh_RP(nlonh,nlath,nv2dd)

    integer :: i, j, k, iv3d, iv2d

    real(RP) :: utmp, vtmp
    real(RP) :: qdry, Rtot
    real(RP) :: psat(nlevh,nlonh,nlath)

    ! Variables that can be directly copied
    !---------------------------------------------------------
    v3dgh_RP(:,:,:,:) = UNDEF
    v2dgh_RP(:,:,:) = UNDEF

    v3dgh_RP(start_z:end_z,start_x:end_x,start_y:end_y,iv3dd_u) = v3dg(:,:,:,iv3d_u)
    v3dgh_RP(start_z:end_z,start_x:end_x,start_y:end_y,iv3dd_v) = v3dg(:,:,:,iv3d_v)
    v3dgh_RP(start_z:end_z,start_x:end_x,start_y:end_y,iv3dd_w) = v3dg(:,:,:,iv3d_w)
    v3dgh_RP(start_z:end_z,start_x:end_x,start_y:end_y,iv3dd_t) = v3dg(:,:,:,iv3d_t)
    v3dgh_RP(start_z:end_z,start_x:end_x,start_y:end_y,iv3dd_p) = v3dg(:,:,:,iv3d_p)
    v3dgh_RP(start_z:end_z,start_x:end_x,start_y:end_y,iv3dd_q) = v3dg(:,:,:,iv3d_q)
    v3dgh_RP(start_z:end_z,start_x:end_x,start_y:end_y,iv3dd_qc) = v3dg(:,:,:,iv3d_qc)
    v3dgh_RP(start_z:end_z,start_x:end_x,start_y:end_y,iv3dd_qr) = v3dg(:,:,:,iv3d_qr)
    v3dgh_RP(start_z:end_z,start_x:end_x,start_y:end_y,iv3dd_qi) = v3dg(:,:,:,iv3d_qi)
    v3dgh_RP(start_z:end_z,start_x:end_x,start_y:end_y,iv3dd_qs) = v3dg(:,:,:,iv3d_qs)
    v3dgh_RP(start_z:end_z,start_x:end_x,start_y:end_y,iv3dd_qg) = v3dg(:,:,:,iv3d_qg)

    ! Rotate U/V (model coord. wind) and obtain Umet/Vmet (true zonal/meridional wind)
    !-------------
    do j = start_y, end_y
    do i = start_x, end_x
    do k = start_z, end_z
      utmp = v3dgh_RP(k,i,j,iv3d_u)
      vtmp = v3dgh_RP(k,i,j,iv3d_v)

      v3dgh_RP(k,i,j,iv3d_u) = utmp * ROTC(i,j,1) - vtmp * ROTC(i,j,2)
      v3dgh_RP(k,i,j,iv3d_v) = utmp * ROTC(i,j,2) + vtmp * ROTC(i,j,1)
    enddo
    enddo
    enddo

    ! RH
    !---------------------------------------------------------

    call ATMOS_SATURATION_psat_all( nlevh, start_z, end_z,   & ! (in)
                                    nlonh, start_x, end_x,   & ! (in)
                                    nlath, start_y, end_y,   & ! (in)
                                    v3dgh_RP(:,:,:,iv3dd_t), & ! (in)
                                    psat(:,:,:)              ) ! (out)

    !$omp parallel do private(k,i,j,qdry,Rtot) schedule(static) collapse(2)
    do j = start_y, end_y
    do i = start_x, end_x
    do k = start_z, end_z
      qdry  = 1.0
      do iv3d = iv3dd_q, iv3dd_qg ! loop over all moisture variables
        qdry  = qdry - v3dgh_RP(k,i,j,iv3d)
      enddo
      Rtot  = Rdry  * qdry + Rvap * v3dgh_RP(k,i,j,iv3dd_q)

      v3dgh_RP(k,i,j,iv3dd_rh) =  v3dgh_RP(k,i,j,iv3dd_q) * v3dgh_RP(k,i,j,iv3dd_p) / psat(k,i,j) * Rvap / Rtot 
    end do
    end do
    end do

    ! Calculate height based the the topography and vertical coordinate
    !---------------------------------------------------------

    v3dgh_RP(start_z:end_z,start_x:end_x,start_y:end_y,iv3dd_hgt) = REAL_CZ(start_z:end_z,start_x:end_x,start_y:end_y)

    ! Surface variables: use the 1st level as the surface (although it is not)
    !---------------------------------------------------------

    v2dgh_RP(start_x:end_x,start_y:end_y,iv2dd_topo) = REAL_CZ(start_z,start_x:end_x,start_y:end_y)

    v2dgh_RP(start_x:end_x,start_y:end_y,iv2dd_ps)   = v3dg(1,1:nlon,1:nlat,iv3d_p)
    v2dgh_RP(start_x:end_x,start_y:end_y,iv2dd_u10m) = v3dg(1,1:nlon,1:nlat,iv3d_u)
    v2dgh_RP(start_x:end_x,start_y:end_y,iv2dd_v10m) = v3dg(1,1:nlon,1:nlat,iv3d_v)
    v2dgh_RP(start_x:end_x,start_y:end_y,iv2dd_t2m)  = v3dg(1,1:nlon,1:nlat,iv3d_t)
    v2dgh_RP(start_x:end_x,start_y:end_y,iv2dd_q2m)  = v3dg(1,1:nlon,1:nlat,iv3d_q)

    !v2dgh_RP(start_x:end_x,start_y:end_y,iv2dd_rain) = [[No way]]

    v2dgh_RP(start_x:end_x,start_y:end_y,iv2dd_ps) = v3dg(1,:,:,iv3d_p)
    v2dgh_RP(start_x:end_x,start_y:end_y,iv2dd_u10m) = v3dg(1,:,:,iv3d_u)
    v2dgh_RP(start_x:end_x,start_y:end_y,iv2dd_v10m) = v3dg(1,:,:,iv3d_v)
    v2dgh_RP(start_x:end_x,start_y:end_y,iv2dd_t2m) = v3dg(1,:,:,iv3d_t)
    v2dgh_RP(start_x:end_x,start_y:end_y,iv2dd_q2m) = v3dg(1,:,:,iv3d_q)

    !v2dgh_RP(start_x:end_x,start_y:end_y,iv2dd_rain) = [[No way]]

    do iv3d = 1, nv3dd
    do j = start_y, end_y
    do i = start_x, end_x
      v3dgh_RP(      1:start_z-1,i,j,iv3d) = v3dgh_RP(start_z,i,j,iv3d)
      v3dgh_RP(end_z+1:nlevh    ,i,j,iv3d) = v3dgh_RP(end_z  ,i,j,iv3d)
    end do
    end do
    end do

    ! Communicate the lateral halo areas
    !---------------------------------------------------------

    do iv3d = 1, nv3dd
      call COMM_vars8( v3dgh_RP(:,:,:,iv3d), iv3d )
    end do
    do iv3d = 1, nv3dd
      call COMM_wait ( v3dgh_RP(:,:,:,iv3d), iv3d )
    end do

    do iv2d = 1, nv2dd
      call COMM_vars8( v2dgh_RP(:,:,iv2d), iv2d )
    end do
    do iv2d = 1, nv2dd
      call COMM_wait ( v2dgh_RP(:,:,iv2d), iv2d )
    end do

    v3dgh = real(v3dgh_RP, kind=RP)
    v2dgh = real(v2dgh_RP, kind=RP)

    return
  end subroutine state_to_history

end module scale_letkf
