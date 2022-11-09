!-------------------------------------------------------------------------------
!> module Communication CartesianC nesting
!!
!! @par Description
!!          Grid module for nesting system
!!
!! @author Team SCALE
!!
!<
!-------------------------------------------------------------------------------
#include "scalelib.h"
module scale_comm_cartesC_nest
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use mpi
  use scale_precision
  use scale_io
  use scale_prof
  use scale_debug
  use scale_atmos_grid_cartesC_index
  use scale_index
  use scale_tracer
  use scale_file_h
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: COMM_CARTESC_NEST_setup
  public :: COMM_CARTESC_NEST_domain_regist_file
  public :: COMM_CARTESC_NEST_parent_info
  public :: COMM_CARTESC_NEST_domain_shape
  public :: COMM_CARTESC_NEST_nestdown_send
  public :: COMM_CARTESC_NEST_nestdown_recv
  public :: COMM_CARTESC_NEST_recvwait_issue_send
  public :: COMM_CARTESC_NEST_recvwait_issue_recv
  public :: COMM_CARTESC_NEST_recv_cancel_send
  public :: COMM_CARTESC_NEST_recv_cancel_recv
  public :: COMM_CARTESC_NEST_test_send
  public :: COMM_CARTESC_NEST_test_recv
  public :: COMM_CARTESC_NEST_finalize

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  type, public :: domain_info
     integer  :: prc_num_x
     integer  :: prc_num_y
     integer  :: KMAX
     integer  :: KHALO
     integer  :: IMAX
     integer  :: IHALO
     integer  :: JMAX
     integer  :: JHALO
     integer  :: OKMAX
     integer  :: LKMAX
     integer  :: UKMAX
     logical  :: periodic_x
     logical  :: periodic_y
     real(RP), allocatable :: latlon_catalogue(:,:,:)
     integer, allocatable :: tile_id(:)
     integer  :: tile_num_x
     integer  :: tile_num_y
     character(len=FILE_HLONG) :: basename
  end type domain_info


  integer,  public :: COMM_CARTESC_NEST_Filiation(10)    !< index of parent-daughter relation (p>0, d<0)
  integer,  public :: HANDLING_NUM                       !< handing number of nesting relation

  integer,  public :: COMM_CARTESC_NEST_INTERP_LEVEL        = 5       !< horizontal interpolation level
  integer,  public :: COMM_CARTESC_NEST_INTERP_WEIGHT_ORDER = 2       !< horizontal interpolation weight order

  logical,  public :: USE_NESTING              = .false.
  logical,  public :: ONLINE_IAM_PARENT        = .false. !< a flag to say "I am a parent"
  logical,  public :: ONLINE_IAM_DAUGHTER      = .false. !< a flag to say "I am a daughter"
  integer,  public :: ONLINE_DOMAIN_NUM        = 1
  logical,  public :: ONLINE_USE_VELZ          = .false.
  logical,  public :: ONLINE_NO_ROTATE         = .false.
  logical,  public :: ONLINE_BOUNDARY_USE_QHYD = .false.

  logical,  public :: ONLINE_RECV_DIAGQHYD = .false.
  logical,  public :: ONLINE_SEND_DIAGQHYD = .false.
  integer,  public :: ONLINE_RECV_QA = 0   !< number of tracer received from the parent domain
  integer,  public :: ONLINE_SEND_QA = 0   !< number of tracer sent to the daughter domain

  real(DP), public :: ONLINE_PARENT_DTSEC   !< parent DT [sec]
  integer,  public :: ONLINE_PARENT_NSTEP   !< parent nsteps

  integer,  public :: ONLINE_DAUGHTER_nprocs = -1

  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  private :: COMM_CARTESC_NEST_ping
  private :: COMM_CARTESC_NEST_parentsize
  private :: COMM_CARTESC_NEST_catalogue
  private :: COMM_CARTESC_NEST_setup_nestdown
  private :: COMM_CARTESC_NEST_importgrid_nestdown
  private :: COMM_CARTESC_NEST_intercomm_nestdown
  private :: COMM_CARTESC_NEST_issuer_of_receive
  private :: COMM_CARTESC_NEST_issuer_of_wait

  interface COMM_CARTESC_NEST_intercomm_nestdown
     module procedure COMM_CARTESC_NEST_intercomm_nestdown_3D
  end interface COMM_CARTESC_NEST_intercomm_nestdown

  interface COMM_CARTESC_NEST_issuer_of_receive
     module procedure COMM_CARTESC_NEST_issuer_of_receive_3D
  end interface COMM_CARTESC_NEST_issuer_of_receive

  interface COMM_CARTESC_NEST_issuer_of_wait
     module procedure COMM_CARTESC_NEST_issuer_of_wait_3D
  end interface COMM_CARTESC_NEST_issuer_of_wait

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  integer, private, parameter        :: MAX_DINFO = 3
  integer, private                   :: num_dom = 0
  type(domain_info), private, target :: dom_info(MAX_DINFO)
  integer,  private   :: I_PARENT = -1

  real(RP), private              :: latlon_local  (4,2)      !< local latlon info [rad]

  integer,  private              :: COMM_CARTESC_NEST_TILE_ALL            !< NUM of TILEs in the local node
  integer,  private              :: COMM_CARTESC_NEST_TILE_ALLMAX_p       !< MAXNUM of TILEs among whole processes for parent
  integer,  private              :: COMM_CARTESC_NEST_TILE_ALLMAX_d       !< MAXNUM of TILEs among whole processes for daughter
  integer,  private, allocatable :: COMM_CARTESC_NEST_TILE_LIST_p(:,:)    !< relationship list in whole system for parent
  integer,  private, allocatable :: COMM_CARTESC_NEST_TILE_LIST_d(:,:)    !< relationship list in whole system for daughter
  integer,  private, allocatable :: COMM_CARTESC_NEST_TILE_LIST_YP(:)     !< yellow-page of daughter targets for parent
  integer,  private              :: NUM_YP                   !< page number of yellow-page

  integer(8),            private :: ONLINE_WAIT_LIMIT         !< limit times of waiting loop in "COMM_CARTESC_NEST_waitall"
  logical,               private :: ONLINE_DAUGHTER_USE_VELZ
  logical,               private :: ONLINE_DAUGHTER_NO_ROTATE
  logical,               private :: ONLINE_AGGRESSIVE_COMM

  integer, private :: TILEAL_KA
  integer, private :: TILEAL_IA
  integer, private :: TILEAL_JA

  integer,  parameter :: I_LON    = 1
  integer,  parameter :: I_LAT    = 2

  integer,  parameter :: I_MIN = 1
  integer,  parameter :: I_MAX = 2

  integer,  parameter :: I_SCLR   = 1                       !< interpolation kinds of grid point (scalar)
  integer,  parameter :: I_ZSTG   = 2                       !< interpolation kinds of grid point (z-axis staggered)
  integer,  parameter :: I_XSTG   = 3                       !< interpolation kinds of grid point (x-axis staggered)
  integer,  parameter :: I_YSTG   = 4                       !< interpolation kinds of grid point (y-axis staggered)

  integer,  parameter :: itp_ng   = 4                       !< # of interpolation kinds of grid point
  integer,  private   :: itp_nh   = 4                       !< # of interpolation kinds of horizontal direction

  integer,  parameter :: tag_lon   = 1
  integer,  parameter :: tag_lat   = 2
  integer,  parameter :: tag_lonuy = 3
  integer,  parameter :: tag_latuy = 4
  integer,  parameter :: tag_lonxv = 5
  integer,  parameter :: tag_latxv = 6
  integer,  parameter :: tag_cz    = 7
  integer,  parameter :: tag_fz    = 8

  integer,  parameter :: tag_dens = 1
  integer,  parameter :: tag_momz = 2
  integer,  parameter :: tag_momx = 3
  integer,  parameter :: tag_momy = 4
  integer,  parameter :: tag_rhot = 5
  integer,  parameter :: tag_qx   = 6

  integer,  parameter :: order_tag_comm = 100000
  integer,  parameter :: order_tag_var  = 1000
  ! intercomm tag id:  IC | VAR |  YP
  ! (total: 6columns)   X   X X   X X X

  integer,  private, parameter :: interp_search_divnum = 10

  integer,  private   :: INTERCOMM_ID(2)

  integer,  private              :: max_isu                ! maximum number of receive/wait issue
  integer,  private              :: max_rq    = 1000       ! maximum number of req: tentative approach
  integer,  private              :: rq_ctl_p               ! for control request id (counting)
  integer,  private              :: rq_ctl_d               ! for control request id (counting)
  integer,  private              :: rq_tot_p               ! for control request id (total number)
  integer,  private              :: rq_tot_d               ! for control request id (total number)
  integer,  private, allocatable :: ireq_p(:)              ! buffer of request-id for parent
  integer,  private, allocatable :: ireq_d(:)              ! buffer of request-id for daughter
  integer,  private, allocatable :: call_order(:)          ! calling order from parent

  real(RP), private, allocatable :: recvbuf_3D(:,:,:,:)    ! buffer of receiver: 3D (with HALO)

  real(RP), private, allocatable :: buffer_ref_LON  (:,:)  ! buffer of communicator: LON
  real(RP), private, allocatable :: buffer_ref_LONUY(:,:)  ! buffer of communicator: LONUY
  real(RP), private, allocatable :: buffer_ref_LONXV(:,:)  ! buffer of communicator: LONXV
  real(RP), private, allocatable :: buffer_ref_LAT  (:,:)  ! buffer of communicator: LAT
  real(RP), private, allocatable :: buffer_ref_LATUY(:,:)  ! buffer of communicator: LATUY
  real(RP), private, allocatable :: buffer_ref_LATXV(:,:)  ! buffer of communicator: LATXV
  real(RP), private, allocatable :: buffer_ref_CZ  (:,:,:) ! buffer of communicator: CZ
  real(RP), private, allocatable :: buffer_ref_FZ  (:,:,:) ! buffer of communicator: FZ


  real(RP), private, allocatable :: buffer_ref_3D (:,:,:)  ! buffer of communicator: 3D data      (with HALO)

  real(RP), private, allocatable :: org_DENS(:,:,:)        ! buffer of communicator: DENS
  real(RP), private, allocatable :: org_MOMZ(:,:,:)        ! buffer of communicator: MOMZ
  real(RP), private, allocatable :: org_MOMX(:,:,:)        ! buffer of communicator: MOMX
  real(RP), private, allocatable :: org_MOMY(:,:,:)        ! buffer of communicator: MOMY
  real(RP), private, allocatable :: org_U_ll(:,:,:)        ! buffer of communicator: U_ll
  real(RP), private, allocatable :: org_V_ll(:,:,:)        ! buffer of communicator: V_ll
  real(RP), private, allocatable :: org_RHOT(:,:,:)        ! buffer of communicator: RHOT
  real(RP), private, allocatable :: org_QTRC(:,:,:,:)      ! buffer of communicator: QTRC

  integer,  private, allocatable :: igrd (:,:,:,:)         ! interpolation target grids in x-axis
  integer,  private, allocatable :: jgrd (:,:,:,:)         ! interpolation target grids in y-axis
  real(RP), private, allocatable :: hfact(:,:,:,:)         ! interpolation factor for horizontal direction
  integer,  private, allocatable :: kgrd (:,:,:,:,:,:)     ! interpolation target grids in z-axis
  real(RP), private, allocatable :: vfact(:,  :,:,:,:)     ! interpolation factor for vertical direction

  integer(8), private :: nwait_p, nwait_d, nrecv, nsend

  character(len=H_SHORT) :: MP_TYPE

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine COMM_CARTESC_NEST_setup ( &
       QA_MP,     &
       MP_TYPE_in )
    use scale_file, only: &
       FILE_open,          &
       FILE_read,          &
       FILE_get_attribute, &
       FILE_get_shape
    use scale_const, only: &
       D2R => CONST_D2R
    use scale_time, only: &
       TIME_NSTEP, &
       TIME_DTSEC
    use scale_prc, only: &
       PRC_abort,            &
       PRC_GLOBAL_domainID,  &
       PRC_IsMaster,         &
       PRC_INTERCOMM_PARENT, &
       PRC_INTERCOMM_CHILD
    use scale_interp, only: &
       INTERP_setup,   &
       INTERP_factor3d
    use scale_comm_cartesC, only: &
       COMM_Bcast
    use scale_atmos_grid_cartesC, only: &
       ATMOS_GRID_CARTESC_CX, &
       ATMOS_GRID_CARTESC_FX, &
       ATMOS_GRID_CARTESC_CY, &
       ATMOS_GRID_CARTESC_FY
    use scale_atmos_grid_cartesC_real, only: &
       ATMOS_GRID_CARTESC_REAL_LON,   &
       ATMOS_GRID_CARTESC_REAL_LAT,   &
       ATMOS_GRID_CARTESC_REAL_LONUY, &
       ATMOS_GRID_CARTESC_REAL_LONXV, &
       ATMOS_GRID_CARTESC_REAL_LONUV, &
       ATMOS_GRID_CARTESC_REAL_LATUY, &
       ATMOS_GRID_CARTESC_REAL_LATXV, &
       ATMOS_GRID_CARTESC_REAL_LATUV, &
       ATMOS_GRID_CARTESC_REAL_CZ,    &
       ATMOS_GRID_CARTESC_REAL_FZ
    use scale_atmos_hydrometeor, only: &
       ATMOS_HYDROMETEOR_dry
    use scale_mapprojection, only: &
       MAPPROJECTION_lonlat2xy
    implicit none

    integer,          intent(in) :: QA_MP
    character(len=*), intent(in) :: MP_TYPE_in

    character(len=H_SHORT) :: COMM_CARTESC_NEST_INTERP_TYPE = 'LINEAR' ! "LINEAR" or "DIST-WEIGHT"
                                                                       !   LINEAR     : bi-linear interpolation
                                                                       !   DIST-WEIGHT: distance-weighted mean of the nearest N-neighbors

    real(RP), allocatable :: X_ref(:,:)
    real(RP), allocatable :: Y_ref(:,:)

    integer :: ONLINE_SPECIFIED_MAXRQ = 0
    integer :: n, i, j
    integer :: fid, ierr
    integer :: parent_id

    logical :: flag_parent
    logical :: flag_child

    integer :: nprocs
    logical :: parent_periodic_x
    logical :: parent_periodic_y

    logical :: error

    namelist / PARAM_COMM_CARTESC_NEST / &
       ONLINE_DOMAIN_NUM,        &
       ONLINE_IAM_PARENT,        &
       ONLINE_IAM_DAUGHTER,      &
       ONLINE_USE_VELZ,          &
       ONLINE_NO_ROTATE,         &
       ONLINE_BOUNDARY_USE_QHYD, &
       ONLINE_AGGRESSIVE_COMM,   &
       ONLINE_WAIT_LIMIT,        &
       ONLINE_SPECIFIED_MAXRQ,   &
       COMM_CARTESC_NEST_INTERP_TYPE,  &
       COMM_CARTESC_NEST_INTERP_LEVEL, &
       COMM_CARTESC_NEST_INTERP_WEIGHT_ORDER

    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("COMM_CARTESC_NEST_setup",*) 'Setup'

    flag_child  = PRC_INTERCOMM_PARENT /= MPI_COMM_NULL ! exist parent, so work as a child
    flag_parent = PRC_INTERCOMM_CHILD  /= MPI_COMM_NULL ! exist child, so work as a parent

    nwait_p = 0
    nwait_d = 0
    nrecv = 0
    nsend = 0

    HANDLING_NUM           = 0
    COMM_CARTESC_NEST_Filiation(:)      = 0
    ONLINE_WAIT_LIMIT      = 999999999
    ONLINE_AGGRESSIVE_COMM = .true.

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_COMM_CARTESC_NEST,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("COMM_CARTESC_NEST_setup",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("COMM_CARTESC_NEST_setup",*) 'Not appropriate names in namelist PARAM_COMM_CARTESC_NEST. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_COMM_CARTESC_NEST)

    PRC_GLOBAL_domainID = ONLINE_DOMAIN_NUM


    if ( ONLINE_IAM_DAUGHTER .or. ONLINE_IAM_PARENT ) then
       USE_NESTING = .true.
    endif

    call INTERP_setup( COMM_CARTESC_NEST_INTERP_WEIGHT_ORDER ) ! [IN]

    select case ( COMM_CARTESC_NEST_INTERP_TYPE )
    case ( 'LINEAR' )
       itp_nh = 4
    case ( 'DIST-WEIGHT' )
       itp_nh = COMM_CARTESC_NEST_INTERP_LEVEL
    case default
       LOG_ERROR("COMM_CARTESC_NEST_setup",*) 'Unsupported type of COMM_CARTESC_NEST_INTERP_TYPE : ', trim(COMM_CARTESC_NEST_INTERP_TYPE)
       LOG_ERROR_CONT(*) '       It must be "LINEAR" or "DIST-WEIGHT"'
       call PRC_abort
    end select


    latlon_local(I_MIN,I_LON) = minval(ATMOS_GRID_CARTESC_REAL_LONUV(:,:)) / D2R
    latlon_local(I_MAX,I_LON) = maxval(ATMOS_GRID_CARTESC_REAL_LONUV(:,:)) / D2R
    latlon_local(I_MIN,I_LAT) = minval(ATMOS_GRID_CARTESC_REAL_LATUV(:,:)) / D2R
    latlon_local(I_MAX,I_LAT) = maxval(ATMOS_GRID_CARTESC_REAL_LATUV(:,:)) / D2R

    if ( .not. USE_NESTING ) return


    DEBUG_DOMAIN_NUM = ONLINE_DOMAIN_NUM
    if( ONLINE_SPECIFIED_MAXRQ > max_rq ) max_rq = ONLINE_SPECIFIED_MAXRQ

    allocate( ireq_p(max_rq)     )
    allocate( ireq_d(max_rq)     )
    allocate( call_order(max_rq) )
    ireq_p(:) = MPI_REQUEST_NULL
    ireq_d(:) = MPI_REQUEST_NULL


    ! ONLINE_(RECV|SEND)_QA can be modified according to the configuration in the other side
    ! See COMM_CARTESC_NEST_parentsize

    if( ONLINE_BOUNDARY_USE_QHYD ) then
       MP_TYPE = MP_TYPE_in
       ONLINE_RECV_QA = QA_MP
    elseif ( ATMOS_HYDROMETEOR_dry ) then
       MP_TYPE = "DRY"
       ONLINE_RECV_QA = 0
    else
       MP_TYPE = "QV"
       ONLINE_RECV_QA = 1
    endif

    if ( ATMOS_HYDROMETEOR_dry ) then
       MP_TYPE = "DRY"
       ONLINE_SEND_QA = 0
    else if ( MP_TYPE_in == "NONE" ) then
       MP_TYPE = "QV"
       ONLINE_SEND_QA = 1
    else
       MP_TYPE = MP_TYPE_in
       ONLINE_SEND_QA = QA_MP
    end if

    LOG_INFO("COMM_CARTESC_NEST_setup",*) "flag_parent", flag_parent, "flag_child", flag_child
    LOG_INFO("COMM_CARTESC_NEST_setup",*) "ONLINE_IAM_PARENT", ONLINE_IAM_PARENT, "ONLINE_IAM_DAUGHTER", ONLINE_IAM_DAUGHTER

    if( flag_parent ) then ! must do first before daughter processes
       !-------------------------------------------------
       if ( .NOT. ONLINE_IAM_PARENT ) then
          LOG_ERROR("COMM_CARTESC_NEST_setup",*) '[NEST_setup] Parent Flag from launcher is not consistent with namelist!'
          LOG_ERROR_CONT(*) 'PARENT - domain : ', ONLINE_DOMAIN_NUM
          call PRC_abort
       endif

       HANDLING_NUM = 1 !HANDLING_NUM + 1
       INTERCOMM_ID(HANDLING_NUM) = ONLINE_DOMAIN_NUM
       COMM_CARTESC_NEST_Filiation(INTERCOMM_ID(HANDLING_NUM)) = 1

       LOG_INFO("COMM_CARTESC_NEST_setup",'(1x,A,I2,A)') 'Online Nesting - PARENT [INTERCOMM_ID:', &
                                                        INTERCOMM_ID(HANDLING_NUM), ' ]'
       LOG_INFO("COMM_CARTESC_NEST_setup",*) 'Online Nesting - INTERCOMM :', PRC_INTERCOMM_CHILD

       call COMM_CARTESC_NEST_ping( HANDLING_NUM )
       call COMM_CARTESC_NEST_parentsize( HANDLING_NUM )
       call COMM_CARTESC_NEST_catalogue( HANDLING_NUM )
       call MPI_BARRIER(PRC_INTERCOMM_CHILD, ierr)

       LOG_INFO("COMM_CARTESC_NEST_setup",'(1x,A)'     ) 'Informations of Daughter Domain'
       LOG_INFO_CONT('(1x,A,I6)'  ) '--- DAUGHTER_nprocs :', ONLINE_DAUGHTER_nprocs
       LOG_INFO_CONT('(1x,A,I6)  ') 'Limit Num. NCOMM req.   :', max_rq

       allocate( org_DENS(KA,IA,JA) )
       allocate( org_MOMZ(KA,IA,JA) )
       allocate( org_MOMX(KA,IA,JA) )
       allocate( org_MOMY(KA,IA,JA) )
       allocate( org_U_ll(KA,IA,JA) )
       allocate( org_V_ll(KA,IA,JA) )
       allocate( org_RHOT(KA,IA,JA) )
       allocate( org_QTRC(KA,IA,JA,max(ONLINE_RECV_QA,1)) )

       call COMM_CARTESC_NEST_setup_nestdown( HANDLING_NUM )

       !---------------------------------- end of parent routines
    endif


    if( flag_child ) then
       !-------------------------------------------------
       if ( .NOT. ONLINE_IAM_DAUGHTER ) then
          LOG_ERROR("COMM_CARTESC_NEST_setup",*) '[NEST_setup] Child Flag from launcher is not consistent with namelist!'
          LOG_ERROR_CONT(*) 'DAUGHTER - domain : ', ONLINE_DOMAIN_NUM
          call PRC_abort
       endif

       HANDLING_NUM = 2 !HANDLING_NUM + 1
       INTERCOMM_ID(HANDLING_NUM) = ONLINE_DOMAIN_NUM - 1
       COMM_CARTESC_NEST_Filiation(INTERCOMM_ID(HANDLING_NUM)) = -1

       LOG_INFO("COMM_CARTESC_NEST_setup",'(1x,A,I2,A)') 'Online Nesting - DAUGHTER [INTERCOMM_ID:', &
                                                        INTERCOMM_ID(HANDLING_NUM), ' ]'
       LOG_INFO("COMM_CARTESC_NEST_setup",*) 'Online Nesting - INTERCOMM :', PRC_INTERCOMM_PARENT

       num_dom = num_dom + 1
       I_PARENT = num_dom
       if ( I_PARENT > MAX_DINFO ) then
          LOG_ERROR("COMM_CARTESC_NEST_setup",*) 'number of domain exeeds the limit'
          call PRC_abort
       end if

       call COMM_CARTESC_NEST_ping( HANDLING_NUM )

       call COMM_CARTESC_NEST_parentsize( HANDLING_NUM )

       nprocs = dom_info(I_PARENT)%prc_num_x * dom_info(I_PARENT)%prc_num_y
       allocate( dom_info(I_PARENT)%latlon_catalogue(nprocs,2,2) )
       call COMM_CARTESC_NEST_catalogue( HANDLING_NUM )
       call MPI_BARRIER(PRC_INTERCOMM_PARENT, ierr)

       call COMM_CARTESC_NEST_domain_relate( I_PARENT )

       TILEAL_KA = dom_info(I_PARENT)%KMAX + dom_info(I_PARENT)%KHALO * 2
       TILEAL_IA = dom_info(I_PARENT)%IMAX * dom_info(I_PARENT)%tile_num_x
       TILEAL_JA = dom_info(I_PARENT)%JMAX * dom_info(I_PARENT)%tile_num_y

       LOG_INFO("COMM_CARTESC_NEST_setup",'(1x,A)'     ) 'Informations of Parent Domain'
       LOG_INFO_CONT('(1x,A,I6)'  ) '--- PARENT_PRC_nprocs   :', nprocs
       LOG_INFO_CONT('(1x,A,I6)'  ) '--- PARENT_PRC_NUM_X    :', dom_info(I_PARENT)%prc_num_x
       LOG_INFO_CONT('(1x,A,I6)'  ) '--- PARENT_PRC_NUM_Y    :', dom_info(I_PARENT)%prc_num_y
       LOG_INFO_CONT('(1x,A,I6)'  ) '--- PARENT_KMAX         :', dom_info(I_PARENT)%KMAX
       LOG_INFO_CONT('(1x,A,I6)'  ) '--- PARENT_IMAX         :', dom_info(I_PARENT)%IMAX
       LOG_INFO_CONT('(1x,A,I6)'  ) '--- PARENT_JMAX         :', dom_info(I_PARENT)%JMAX
       LOG_INFO_CONT('(1x,A,F9.3)') '--- PARENT_DTSEC        :', ONLINE_PARENT_DTSEC
       LOG_INFO_CONT('(1x,A,I6)'  ) '--- PARENT_NSTEP        :', ONLINE_PARENT_NSTEP
       LOG_INFO_CONT('(1x,A)'     ) 'Informations of Daughter Domain [me]'
       LOG_INFO_CONT('(1x,A,F9.3)') '--- DAUGHTER_DTSEC      :', TIME_DTSEC
       LOG_INFO_CONT('(1x,A,I6)'  ) '--- DAUGHTER_NSTEP      :', TIME_NSTEP
       LOG_INFO_CONT('(1x,A)'     ) 'Informations of Target Tiles'
       LOG_INFO_CONT('(1x,A,I6)'  ) '--- TILEALL_KA      :', TILEAL_KA
       LOG_INFO_CONT('(1x,A,I6)'  ) '--- TILEALL_IA      :', TILEAL_IA
       LOG_INFO_CONT('(1x,A,I6)'  ) '--- TILEALL_JA      :', TILEAL_JA
       LOG_INFO_CONT('(1x,A,I6)  ') 'Limit Num. NCOMM req. :', max_rq

       allocate( buffer_ref_LON  (TILEAL_IA,TILEAL_JA) )
       allocate( buffer_ref_LONUY(TILEAL_IA,TILEAL_JA) )
       allocate( buffer_ref_LONXV(TILEAL_IA,TILEAL_JA) )
       allocate( buffer_ref_LAT  (TILEAL_IA,TILEAL_JA) )
       allocate( buffer_ref_LATUY(TILEAL_IA,TILEAL_JA) )
       allocate( buffer_ref_LATXV(TILEAL_IA,TILEAL_JA) )

       allocate( buffer_ref_CZ(TILEAL_KA,TILEAL_IA,TILEAL_JA) )
       allocate( buffer_ref_FZ(TILEAL_KA,TILEAL_IA,TILEAL_JA) )

       allocate( buffer_ref_3D(TILEAL_KA,TILEAL_IA,TILEAL_JA) )

       allocate( igrd (IA,JA,itp_nh,itp_ng) )
       allocate( jgrd (IA,JA,itp_nh,itp_ng) )
       allocate( hfact(IA,JA,itp_nh,itp_ng) )
       allocate( kgrd (KA,2,IA,JA,itp_nh,itp_ng) )
       allocate( vfact(KA,  IA,JA,itp_nh,itp_ng) )

       call COMM_CARTESC_NEST_setup_nestdown( HANDLING_NUM )


       select case ( COMM_CARTESC_NEST_INTERP_TYPE )
       case ( 'LINEAR' )

          allocate( X_ref(TILEAL_IA,TILEAL_JA) )
          allocate( Y_ref(TILEAL_IA,TILEAL_JA) )

          ! for scalar points
          call MAPPROJECTION_lonlat2xy( TILEAL_IA, 1, TILEAL_IA, &
                                        TILEAL_JA, 1, TILEAL_JA, &
                                        buffer_ref_LON(:,:),   & ! [IN]
                                        buffer_ref_LAT(:,:),   & ! [IN]
                                        X_ref(:,:), Y_ref(:,:) ) ! [OUT]
          call INTERP_factor3d( TILEAL_KA, KHALO+1, TILEAL_KA-KHALO, &
                                TILEAL_IA, TILEAL_JA,                &
                                KA, KS, KE, IA, JA,                  &
                                X_ref(:,:), Y_ref(:,:),            & ! [IN]
                                buffer_ref_CZ (:,:,:),             & ! [IN]
                                ATMOS_GRID_CARTESC_CX(:),          & ! [IN]
                                ATMOS_GRID_CARTESC_CY(:),          & ! [IN]
                                ATMOS_GRID_CARTESC_REAL_CZ(:,:,:), & ! [IN]
                                igrd (    :,:,:,I_SCLR),           & ! [OUT]
                                jgrd (    :,:,:,I_SCLR),           & ! [OUT]
                                hfact(    :,:,:,I_SCLR),           & ! [OUT]
                                kgrd (:,:,:,:,:,I_SCLR),           & ! [OUT]
                                vfact(:,  :,:,:,I_SCLR)            ) ! [OUT]

          ! for z staggered points
          call INTERP_factor3d( TILEAL_KA+1, KHALO+1, TILEAL_KA+1-KHALO, &
                                TILEAL_IA, TILEAL_JA,                    &
                                KA, KS, KE, IA, JA,                      &
                                X_ref(:,:), Y_ref(:,:),               & ! [IN]
                                buffer_ref_FZ (:,:,:),                & ! [IN]
                                ATMOS_GRID_CARTESC_CX(:),             & ! [IN]
                                ATMOS_GRID_CARTESC_CY(:),             & ! [IN]
                                ATMOS_GRID_CARTESC_REAL_FZ(1:KA,:,:), & ! [IN]
                                igrd (    :,:,:,I_ZSTG),              & ! [OUT]
                                jgrd (    :,:,:,I_ZSTG),              & ! [OUT]
                                hfact(    :,:,:,I_ZSTG),              & ! [OUT]
                                kgrd (:,:,:,:,:,I_ZSTG),              & ! [OUT]
                                vfact(:,  :,:,:,I_ZSTG)               ) ! [OUT]

          ! for x staggered points
          call MAPPROJECTION_lonlat2xy( TILEAL_IA, 1, TILEAL_IA, &
                                        TILEAL_JA, 1, TILEAL_JA, &
                                        buffer_ref_LONUY(:,:),   & ! [IN]
                                        buffer_ref_LATUY(:,:),   & ! [IN]
                                        X_ref(:,:), Y_ref(:,:)   ) ! [OUT]
          call INTERP_factor3d( TILEAL_KA, KHALO+1, TILEAL_KA-KHALO, &
                                TILEAL_IA, TILEAL_JA,                &
                                KA, KS, KE, IA, JA,                  &
                                X_ref(:,:), Y_ref(:,:),            & ! [IN]
                                buffer_ref_CZ  (:,:,:),            & ! [IN]
                                ATMOS_GRID_CARTESC_FX(1:IA),       & ! [IN]
                                ATMOS_GRID_CARTESC_CY(:),          & ! [IN]
                                ATMOS_GRID_CARTESC_REAL_CZ(:,:,:), & ! [IN]
                                igrd (    :,:,:,I_XSTG),           & ! [OUT]
                                jgrd (    :,:,:,I_XSTG),           & ! [OUT]
                                hfact(    :,:,:,I_XSTG),           & ! [OUT]
                                kgrd (:,:,:,:,:,I_XSTG),           & ! [OUT]
                                vfact(:,  :,:,:,I_XSTG)            ) ! [OUT]

          ! for y staggered points
          call MAPPROJECTION_lonlat2xy( TILEAL_IA, 1, TILEAL_IA, &
                                        TILEAL_JA, 1, TILEAL_JA, &
                                        buffer_ref_LONXV(:,:),   & ! [IN]
                                        buffer_ref_LATXV(:,:),   & ! [IN]
                                        X_ref(:,:), Y_ref(:,:)   ) ! [OUT]
          call INTERP_factor3d( TILEAL_KA, KHALO+1, TILEAL_KA-KHALO, &
                                TILEAL_IA, TILEAL_JA,                &
                                KA, KS, KE, IA, JA,                  &
                                X_ref(:,:), Y_ref(:,:),            & ! [IN]
                                buffer_ref_CZ  (:,:,:),            & ! [IN]
                                ATMOS_GRID_CARTESC_CX(:),          & ! [IN]
                                ATMOS_GRID_CARTESC_FY(1:JA),       & ! [IN]
                                ATMOS_GRID_CARTESC_REAL_CZ(:,:,:), & ! [IN]
                                igrd (    :,:,:,I_YSTG),           & ! [OUT]
                                jgrd (    :,:,:,I_YSTG),           & ! [OUT]
                                hfact(    :,:,:,I_YSTG),           & ! [OUT]
                                kgrd (:,:,:,:,:,I_YSTG),           & ! [OUT]
                                vfact(:,  :,:,:,I_YSTG)            ) ! [OUT]

          deallocate( X_ref, Y_ref )

       case ( 'DIST-WEIGHT' )

          ! for scalar points
          call INTERP_factor3d( itp_nh,                              &
                                TILEAL_KA, KHALO+1, TILEAL_KA-KHALO, &
                                TILEAL_IA, TILEAL_JA,                &
                                KA, KS, KE, IA, JS,                  &
                                buffer_ref_LON(:,:),                & ! [IN]
                                buffer_ref_LAT(:,:),                & ! [IN]
                                buffer_ref_CZ (:,:,:),              & ! [IN]
                                ATMOS_GRID_CARTESC_REAL_LON(:,:),   & ! [IN]
                                ATMOS_GRID_CARTESC_REAL_LAT(:,:),   & ! [IN]
                                ATMOS_GRID_CARTESC_REAL_CZ (:,:,:), & ! [IN]
                                igrd (    :,:,:,I_SCLR),            & ! [OUT]
                                jgrd (    :,:,:,I_SCLR),            & ! [OUT]
                                hfact(    :,:,:,I_SCLR),            & ! [OUT]
                                kgrd (:,:,:,:,:,I_SCLR),            & ! [OUT]
                                vfact(:,  :,:,:,I_SCLR)             ) ! [OUT]

          ! for z staggered points
          call INTERP_factor3d( itp_nh,                            &
                                TILEAL_KA, KHALO, TILEAL_KA-KHALO, &
                                TILEAL_IA, TILEAL_JA,              &
                                KA, KS, KE, IA, JA,                &
                                buffer_ref_LON(:,:),                   & ! [IN]
                                buffer_ref_LAT(:,:),                   & ! [IN]
                                buffer_ref_FZ (:,:,:),                 & ! [IN]
                                ATMOS_GRID_CARTESC_REAL_LON(:,:),      & ! [IN]
                                ATMOS_GRID_CARTESC_REAL_LAT(:,:),      & ! [IN]
                                ATMOS_GRID_CARTESC_REAL_FZ (1:KA,:,:), & ! [IN]
                                igrd (    :,:,:,I_ZSTG),               & ! [OUT]
                                jgrd (    :,:,:,I_ZSTG),               & ! [OUT]
                                hfact(    :,:,:,I_ZSTG),               & ! [OUT]
                                kgrd (:,:,:,:,:,I_ZSTG),               & ! [OUT]
                                vfact(:,  :,:,:,I_ZSTG)                ) ! [OUT]

          ! for x staggered points
          call INTERP_factor3d( itp_nh,                              &
                                TILEAL_KA, KHALO+1, TILEAL_KA-KHALO, &
                                TILEAL_IA, TILEAL_JA,                &
                                KA, KS, KE, IA, JA,                  &
                                buffer_ref_LONUY(:,:),                    & ! [IN]
                                buffer_ref_LATUY(:,:),                    & ! [IN]
                                buffer_ref_CZ  (:,:,:),                   & ! [IN]
                                ATMOS_GRID_CARTESC_REAL_LONUY(1:IA,1:JA), & ! [IN]
                                ATMOS_GRID_CARTESC_REAL_LATUY(1:IA,1:JA), & ! [IN]
                                ATMOS_GRID_CARTESC_REAL_CZ(:,:,:),        & ! [IN]
                                igrd (    :,:,:,I_XSTG),                  & ! [OUT]
                                jgrd (    :,:,:,I_XSTG),                  & ! [OUT]
                                hfact(    :,:,:,I_XSTG),                  & ! [OUT]
                                kgrd (:,:,:,:,:,I_XSTG),                  & ! [OUT]
                                vfact(:,  :,:,:,I_XSTG)                   ) ! [OUT]

          ! for y staggered points
          call INTERP_factor3d( itp_nh,                              &
                                TILEAL_KA, KHALO+1, TILEAL_KA-KHALO, &
                                TILEAL_IA, TILEAL_JA,                &
                                KA, KS, KE, IA, JA,                  &
                                buffer_ref_LONXV(:,:),                    & ! [IN]
                                buffer_ref_LATXV(:,:),                    & ! [IN]
                                buffer_ref_CZ  (:,:,:),                   & ! [IN]
                                ATMOS_GRID_CARTESC_REAL_LONXV(1:IA,1:JA), & ! [IN]
                                ATMOS_GRID_CARTESC_REAL_LATXV(1:IA,1:JA), & ! [IN]
                                ATMOS_GRID_CARTESC_REAL_CZ(:,:,:),        & ! [IN]
                                igrd (    :,:,:,I_YSTG),                  & ! [OUT]
                                jgrd (    :,:,:,I_YSTG),                  & ! [OUT]
                                hfact(    :,:,:,I_YSTG),                  & ! [OUT]
                                kgrd (:,:,:,:,:,I_YSTG),                  & ! [OUT]
                                vfact(:,  :,:,:,I_YSTG)                   ) ! [OUT]

       end select

       !---------------------------------- end of child routines
    end if

    !LOG_INFO("COMM_CARTESC_NEST_setup",'(1x,A,I2)') 'Number of Related Domains :', HANDLING_NUM
    !if ( HANDLING_NUM > 2 ) then
    !   f( IO_L ) LOG_ERROR("COMM_CARTESC_NEST_setup",*) 'Too much handing domains (up to 2)'
    !   call PRC_abort
    !endif

    return
  end subroutine COMM_CARTESC_NEST_setup

  !-----------------------------------------------------------------------------
  !> offline setup
  subroutine COMM_CARTESC_NEST_domain_regist_file( &
       dom_id, &
       PARENT_BASENAME,  &
       PARENT_PRC_NUM_X, &
       PARENT_PRC_NUM_Y, &
       LATLON_CATALOGUE_FNAME    )
    use scale_prc, only: &
       PRC_IsMaster, &
       PRC_abort
    use scale_file, only: &
       FILE_open, &
       FILE_get_attribute, &
       FILE_get_shape, &
       FILE_read
    use scale_comm_cartesC, only: &
       COMM_Bcast
    integer, intent(out) :: dom_id

    character(len=*), intent(in) :: PARENT_BASENAME
    integer,          intent(in), optional :: PARENT_PRC_NUM_X
    integer,          intent(in), optional :: PARENT_PRC_NUM_Y
    character(len=*), intent(in), optional :: LATLON_CATALOGUE_FNAME !< metadata files for lat-lon domain for all processes

    type(domain_info), pointer :: dinfo

    integer :: nprocs
    integer :: pnum_x(1), pnum_y(1)
    integer :: imaxg(1), jmaxg(1)
    integer :: dims(1), dims2(1)
    integer :: halos(2)
    integer :: parent_x, parent_xh
    integer :: parent_y, parent_yh

    real(RP), allocatable :: work(:,:), work_uv(:,:), minmax(:,:,:)

    character(len=H_LONG) :: fname
    integer :: fid
    integer :: parent_id
    logical :: existed, error
    integer :: ierr

    integer :: i, j, n

    do n = 1, num_dom
       if ( dom_info(n)%basename == parent_basename ) then
          dom_id = n
          return
       end if
    end do

    num_dom = num_dom + 1
    dom_id = num_dom
    if ( dom_id > MAX_DINFO ) then
       LOG_ERROR("COMM_CARTESC_NEST_domain_regist_file",*) 'number of domains exceed the limit'
       call PRC_abort
    end if
    dinfo => dom_info(dom_id)
    dinfo%basename = parent_basename

    if ( PRC_IsMaster ) then
       call FILE_open( PARENT_BASENAME, & ! (in)
                       fid,                     & ! (out)
                       aggregate = .false.      ) ! (in)

       call FILE_get_attribute( fid, "global", "scale_atmos_grid_cartesC_index_imaxg", &
                                imaxg(:), existed=existed                        )
       if ( existed ) then
          call FILE_get_attribute( fid, "global", "scale_cartesC_prc_num_x", &
                                   pnum_x(:)                                 )
          dinfo%prc_num_x = pnum_x(1)
          dinfo%IMAX = imaxg(1) / pnum_x(1)

          call FILE_get_attribute( fid, "x", "halo_global", & ! (in)
                                   halos(:)                 ) ! (out)
          dinfo%IHALO = halos(1)
          call FILE_get_attribute( fid, "global", "scale_cartesC_prc_periodic_x", &
                                   dinfo%periodic_x                               )

          call FILE_get_attribute( fid, "global", "scale_atmos_grid_cartesC_index_jmaxg", &
                                   jmaxg(:)                                               )
          call FILE_get_attribute( fid, "global", "scale_cartesC_prc_num_y", &
                                   pnum_y(:)                                 )
          dinfo%prc_num_y = pnum_y(1)
          dinfo%JMAX = jmaxg(1) / pnum_y(1)
          call FILE_get_attribute( fid, "y", "halo_global", & ! (in)
                                   halos(:)                 ) ! (out)
          dinfo%JHALO = halos(1)
          call FILE_get_attribute( fid, "global", "scale_cartesC_prc_periodic_y", &
                                   dinfo%periodic_y                               )

       else
          ! for old file (for backward compatibility)

          if ( present(PARENT_PRC_NUM_X) .and. present(PARENT_PRC_NUM_Y) ) then
             dinfo%prc_num_x = PARENT_PRC_NUM_X
             dinfo%prc_num_y = PARENT_PRC_NUM_Y
          else
             LOG_ERROR("COMM_CARTESC_NEST_domain_regist_file",*) 'PARENT_PRC_NUM_(X|Y) is needed for files generated by the older version'
             call PRC_abort
          end if

          call FILE_get_shape( fid, "CX", dims(:) )
          call FILE_get_shape( fid, "x", dims2(:) )
          dinfo%IHALO = dims2(1) + IHALO * 2 - dims(1) ! assume IHALO is the same
          dinfo%IMAX = dims(1) - dinfo%IHALO*2

          call FILE_get_shape( fid, "CY", dims(:) )
          call FILE_get_shape( fid, "y", dims2(:) )
          dinfo%JHALO = dims2(1) + JHALO * 2 - dims(1) ! assume JHALO is the same
          dinfo%JMAX = dims(1) - dinfo%JHALO*2

          dinfo%periodic_x = .false.
          dinfo%periodic_y = .false.

       endif

       call FILE_get_attribute( fid, "global", "scale_atmos_grid_cartesC_index_kmax", &
                                dims(:), existed=existed                        )
       if ( existed ) then
          dinfo%KMAX = dims(1)
       else
          call FILE_get_shape( fid, "z", dims(:), error=error )
          if ( error ) then
             dinfo%KMAX = 0
          else
             dinfo%KMAX = dims(1)
          endif
       end if
       dinfo%KHALO = 0

       call FILE_get_attribute( fid, "global", "scale_ocean_grid_cartesC_index_kmax", &
                                dims(:), existed=existed                        )
       if ( existed ) then
          dinfo%OKMAX = dims(1)
       else
          call FILE_get_shape( fid, "oz", dims(:), error=error )
          if ( error ) then
             dinfo%OKMAX = 0
          else
             dinfo%OKMAX = dims(1)
          endif
       end if

       call FILE_get_attribute( fid, "global", "scale_land_grid_cartesC_index_kmax", &
                                dims(:), existed=existed                       )
       if ( existed ) then
          dinfo%LKMAX = dims(1)
       else
          call FILE_get_shape( fid, "lz", dims(:), error=error )
          if ( error ) then
             dinfo%LKMAX = 0
          else
             dinfo%LKMAX = dims(1)
          endif
       end if

       call FILE_get_attribute( fid, "global", "scale_urban_grid_cartesC_index_kmax", &
                                dims(:), existed=existed                       )
       if ( existed ) then
          dinfo%UKMAX = dims(1)
       else
          call FILE_get_shape( fid, "uz", dims(:), error=error )
          if ( error ) then
             dinfo%UKMAX = 0
          else
             dinfo%UKMAX = dims(1)
          endif
       end if

    end if ! master node


    call COMM_Bcast( dinfo%prc_num_x )
    call COMM_Bcast( dinfo%prc_num_y )
    call COMM_Bcast( dinfo%KMAX )
    call COMM_Bcast( dinfo%OKMAX )
    call COMM_Bcast( dinfo%LKMAX )
    call COMM_Bcast( dinfo%UKMAX )
    call COMM_Bcast( dinfo%IMAX )
    call COMM_Bcast( dinfo%JMAX )
    call COMM_Bcast( dinfo%IHALO )
    call COMM_Bcast( dinfo%JHALO )
    call COMM_Bcast( dinfo%periodic_x )
    call COMM_Bcast( dinfo%periodic_y )

    !--- latlon catalogue
    nprocs = dinfo%prc_num_x * dinfo%prc_num_y
    allocate( dinfo%latlon_catalogue(nprocs,2,2) )

    if ( PRC_IsMaster ) then

       existed =  present(LATLON_CATALOGUE_FNAME)
       if ( existed ) then
          existed = LATLON_CATALOGUE_FNAME /= ""
       end if
       if ( existed ) then
          ! read from catalogue file

          fid = IO_get_available_fid()
          call IO_get_fname(fname, LATLON_CATALOGUE_FNAME)
          open( fid,                  &
                file   = fname,       &
                form   = 'formatted', &
                status = 'old',       &
                iostat = ierr         )

          if ( ierr /= 0 ) then
             LOG_ERROR("COMM_CARTESC_NEST_domain_regist_file",*) 'cannot open latlon-catalogue file!: ', trim(fname)
             call PRC_abort
          endif

          do i = 1, nprocs
             read(fid,'(i8,4f32.24)',iostat=ierr) parent_id, &
                                                  dinfo%latlon_catalogue(i,I_MIN,I_LON), dinfo%latlon_catalogue(i,I_MAX,I_LON), & ! LON: MIN, MAX
                                                  dinfo%latlon_catalogue(i,I_MIN,I_LAT), dinfo%latlon_catalogue(i,I_MAX,I_LAT)    ! LAT: MIN, MAX
             if ( ierr /= 0 .or. i /= parent_id ) then
                LOG_ERROR("COMM_CARTESC_NEST_domain_regist_file",*) 'catalogue file is invalid, ', trim(fname)
                call PRC_abort
             end if
             if ( ierr /= 0 ) exit
          enddo
          close(fid)

       else
          ! read from netcdf file

          allocate( minmax(nprocs,2,2) )

          n = 1
          do j = 1, dinfo%prc_num_y
          do i = 1, dinfo%prc_num_x
             call FILE_open( PARENT_BASENAME,     & ! (in)
                             fid,                 & ! (out)
                             aggregate = .false., & ! (in)
                             rankid    = n-1      ) ! (in)

             call FILE_get_shape( fid, "xh", dims(:) )
             parent_xh = dims(1)
             call FILE_get_shape( fid, "yh", dims(:) )
             parent_yh = dims(1)
             allocate( work_uv( parent_xh, parent_yh ) )

             if ( dinfo%periodic_x .or. dinfo%periodic_y ) then
                call FILE_get_shape( fid, "x", dims(:) )
                parent_x = dims(1)
                call FILE_get_shape( fid, "y", dims(:) )
                parent_y = dims(1)
             end if

             call FILE_read( fid, "lon_uv", work_uv(:,:) )
             dinfo%latlon_catalogue(n,I_MIN,I_LON) = minval( work_uv(:,:) )
             dinfo%latlon_catalogue(n,I_MAX,I_LON) = maxval( work_uv(:,:) )

             if ( i > 1 ) then
                dinfo%latlon_catalogue(n,I_MIN,I_LON) = min( dinfo%latlon_catalogue(n,I_MIN,I_LON), minmax(n-1,I_MIN,I_LON) )
                dinfo%latlon_catalogue(n,I_MAX,I_LON) = max( dinfo%latlon_catalogue(n,I_MAX,I_LON), minmax(n-1,I_MAX,I_LON) )
             else
                if ( dinfo%periodic_x ) then
                   allocate( work( parent_x, parent_yh ) )
                   call FILE_read( fid, "lon_xv", work(:,:) )
                   ! This assumes an equally spaced grid
                   work(1,:) = work(1,:) * 2.0_RP - work_uv(1,:)
                   dinfo%latlon_catalogue(n,I_MIN,I_LON) = min( dinfo%latlon_catalogue(n,I_MIN,I_LON), minval( work ) )
                   dinfo%latlon_catalogue(n,I_MAX,I_LON) = max( dinfo%latlon_catalogue(n,I_MAX,I_LON), maxval( work ) )
                   deallocate( work )
                end if
             end if
             minmax(n,I_MIN,I_LON) = minval( work_uv(parent_xh,:) )
             minmax(n,I_MAX,I_LON) = maxval( work_uv(parent_xh,:) )

             call FILE_read( fid, "lat_uv", work_uv(:,:) )
             dinfo%latlon_catalogue(n,I_MIN,I_LAT) = minval( work_uv(:,:) )
             dinfo%latlon_catalogue(n,I_MAX,I_LAT) = maxval( work_uv(:,:) )

             if ( j > 1 ) then
                dinfo%latlon_catalogue(n,I_MIN,I_LAT) = min( dinfo%latlon_catalogue(n,I_MIN,I_LAT), minmax(n-dinfo%prc_num_x,I_MIN,I_LAT) )
                dinfo%latlon_catalogue(n,I_MAX,I_LAT) = max( dinfo%latlon_catalogue(n,I_MAX,I_LAT), minmax(n-dinfo%prc_num_x,I_MAX,I_LAT) )
             else
                if ( dinfo%periodic_y ) then
                   allocate( work( parent_xh, parent_y ) )
                   call FILE_read( fid, "lat_uy", work(:,:) )
                   ! This assumes an equally spaced grid
                   work(:,1) = work(:,1) * 2.0_RP - work_uv(:,1)
                   dinfo%latlon_catalogue(n,I_MIN,I_LAT) = min( dinfo%latlon_catalogue(n,I_MIN,I_LAT), minval( work(:,1) ) )
                   dinfo%latlon_catalogue(n,I_MAX,I_LAT) = max( dinfo%latlon_catalogue(n,I_MAX,I_LAT), maxval( work(:,1) ) )
                   deallocate( work )
                end if
             end if
             minmax(n,I_MIN,I_LAT) = minval( work_uv(:,parent_yh) )
             minmax(n,I_MAX,I_LAT) = maxval( work_uv(:,parent_yh) )

             deallocate( work_uv )

             n = n + 1
          enddo
          enddo

          deallocate( minmax )

       endif

    end if ! master node

    call COMM_Bcast( nprocs, 2, 2, dinfo%latlon_catalogue(:,:,:) )

    call COMM_CARTESC_NEST_domain_relate(dom_id)


  end subroutine COMM_CARTESC_NEST_domain_regist_file

  !-----------------------------------------------------------------------------
  !> Solve relationship between ParentDomain & Daughter Domain
  subroutine COMM_CARTESC_NEST_domain_relate( &
       dom_id )
    use scale_prc, only: &
       PRC_myrank,   &
       PRC_abort
    implicit none

    integer, intent(in) :: dom_id !< id number of domain information

    type(domain_info), pointer :: dinfo

    integer :: nprocs
    integer :: x_min, x_max
    integer :: y_min, y_max
    logical :: hit(2,2)
    real(RP) :: dx, dy
    integer :: p, i, j
    !---------------------------------------------------------------------------

    if ( dom_id < 1 .or. dom_id > num_dom ) then
       LOG_ERROR("COMM_CARTESC_NEST_domain_relate",*) "domain id is invalid: ", dom_id
       call PRC_abort
    end if

    dinfo => dom_info(dom_id)
    nprocs = dinfo%prc_num_x * dinfo%prc_num_y

    x_min = dinfo%prc_num_x
    x_max = -1
    y_min = dinfo%prc_num_y
    y_max = -1
    hit(:,:) = .false.

    do p = 1, nprocs
       dx = ( dinfo%latlon_catalogue(p,I_MAX,I_LON) - dinfo%latlon_catalogue(p,I_MIN,I_LON) ) / dinfo%IMAX
       dy = ( dinfo%latlon_catalogue(p,I_MAX,I_LAT) - dinfo%latlon_catalogue(p,I_MIN,I_LAT) ) / dinfo%JMAX
       if ( ( (     latlon_local(I_MIN,I_LON) >= dinfo%latlon_catalogue(p,I_MIN,I_LON) - dx &
              .AND. latlon_local(I_MIN,I_LON) <= dinfo%latlon_catalogue(p,I_MAX,I_LON) + dx ) .OR. &
              (     latlon_local(I_MAX,I_LON) >= dinfo%latlon_catalogue(p,I_MIN,I_LON) - dx &
              .AND. latlon_local(I_MAX,I_LON) <= dinfo%latlon_catalogue(p,I_MAX,I_LON) + dx ) .OR. &
              (     dinfo%latlon_catalogue(p,I_MIN,I_LON) >= latlon_local(I_MIN,I_LON) - dx &
              .AND. dinfo%latlon_catalogue(p,I_MIN,I_LON) <= latlon_local(I_MAX,I_LON) + dx ) .OR. &
              (     dinfo%latlon_catalogue(p,I_MAX,I_LON) >= latlon_local(I_MIN,I_LON) - dx &
              .AND. dinfo%latlon_catalogue(p,I_MAX,I_LON) <= latlon_local(I_MAX,I_LON) + dx ) ) .AND. &
            ( (     latlon_local(I_MIN,I_LAT) >= dinfo%latlon_catalogue(p,I_MIN,I_LAT) - dy &
              .AND. latlon_local(I_MIN,I_LAT) <= dinfo%latlon_catalogue(p,I_MAX,I_LAT) + dy ) .OR. &
              (     latlon_local(I_MAX,I_LAT) >= dinfo%latlon_catalogue(p,I_MIN,I_LAT) - dy &
              .AND. latlon_local(I_MAX,I_LAT) <= dinfo%latlon_catalogue(p,I_MAX,I_LAT) + dy ) .OR. &
              (     dinfo%latlon_catalogue(p,I_MIN,I_LAT) >= latlon_local(I_MIN,I_LAT) - dy &
              .AND. dinfo%latlon_catalogue(p,I_MIN,I_LAT) <= latlon_local(I_MAX,I_LAT) + dy ) .OR. &
              (     dinfo%latlon_catalogue(p,I_MAX,I_LAT) >= latlon_local(I_MIN,I_LAT) - dy &
              .AND. dinfo%latlon_catalogue(p,I_MAX,I_LAT) <= latlon_local(I_MAX,I_LAT) + dy ) ) ) then
          if ( dinfo%latlon_catalogue(p,I_MIN,I_LON) <= latlon_local(I_MIN,I_LON) ) hit(I_MIN,I_LON) = .true.
          if ( dinfo%latlon_catalogue(p,I_MAX,I_LON) >= latlon_local(I_MAX,I_LON) ) hit(I_MAX,I_LON) = .true.
          if ( dinfo%latlon_catalogue(p,I_MIN,I_LAT) <= latlon_local(I_MIN,I_LAT) ) hit(I_MIN,I_LAT) = .true.
          if ( dinfo%latlon_catalogue(p,I_MAX,I_LAT) >= latlon_local(I_MAX,I_LAT) ) hit(I_MAX,I_LAT) = .true.
          i = mod(p-1, dinfo%prc_num_x)
          j = (p-1) / dinfo%prc_num_x
          if ( i < x_min ) x_min = i
          if ( i > x_max ) x_max = i
          if ( j < y_min ) y_min = j
          if ( j > y_max ) y_max = j
       end if
    end do

    if ( .not. ( hit(I_MIN,I_LON) .and. hit(I_MAX,I_LON) .and. hit(I_MIN,I_LAT) .and. hit(I_MAX,I_LAT) ) ) then
       LOG_ERROR("COMM_CARTESC_NEST_domain_relate",*) 'region of daughter domain is larger than that of parent'
       LOG_ERROR_CONT(*)              '                                  at rank:', PRC_myrank, ' of domain:', ONLINE_DOMAIN_NUM
       LOG_ERROR_CONT(*) 'LON MIN: ',hit(I_MIN,I_LON), ', LON MAX: ',hit(I_MAX,I_LON), ', LAT MIN: ',hit(I_MIN,I_LAT), ', LAT MAX: ',hit(I_MAX,I_LAT)
       LOG_ERROR_CONT('(A,F12.6,1x,F12.6)') 'daughter local (me) MIN-MAX: LON=', &
            latlon_local(I_MIN,I_LON), latlon_local(I_MAX,I_LON)
       do p = 1, nprocs
          LOG_ERROR_CONT('(A,I5,A,F12.6,1x,F12.6)') '     parent (', p,') MIN-MAX: LON=', &
                     dinfo%latlon_catalogue(p,I_MIN,I_LON) ,dinfo%latlon_catalogue(p,I_MAX,I_LON)
       enddo
       LOG_ERROR_CONT('(A,F12.6,1x,F12.6)') 'daughter local (me): MIN-MAX LAT=', &
            latlon_local(I_MIN,I_LAT), latlon_local(I_MAX,I_LAT)
       do p = 1, nprocs
          LOG_ERROR_CONT('(A,I5,A,F12.6,1x,F12.6)') '     parent (', p,') MIN-MAX: LAT=', &
                     dinfo%latlon_catalogue(p,I_MIN,I_LAT) ,dinfo%latlon_catalogue(p,I_MAX,I_LAT)
       enddo
       call PRC_abort
    end if



    dinfo%tile_num_x = x_max - x_min + 1
    dinfo%tile_num_y = y_max - y_min + 1

    allocate( dinfo%tile_id(dinfo%tile_num_x * dinfo%tile_num_y) )

    LOG_INFO("COMM_CARTESC_NEST_domain_relate",'(1x,A)') 'NEST: target process tile in parent domain'
    p = 1
    do j = 1, dinfo%tile_num_y
    do i = 1, dinfo%tile_num_x
       dinfo%tile_id(p) = x_min + i - 1 + (y_min + j - 1) * dinfo%prc_num_x
       LOG_INFO_CONT('(1x,A,I4,A,I6)') '(', p, ') target mpi-process:', dinfo%tile_id(p)
       p = p + 1
    enddo
    enddo

    return
  end subroutine COMM_CARTESC_NEST_domain_relate

  !-----------------------------------------------------------------------------
  !> Return infomation of parent domain (for offline)
  subroutine COMM_CARTESC_NEST_parent_info( &
       dom_id, &
       KMAX,  &
       LKMAX, &
       IMAXG, &
       JMAXG, &
       num_tile, &
       tile_id)
    use scale_prc, only: &
       PRC_abort
    integer, intent(in) :: dom_id

    integer, intent(out), optional :: KMAX
    integer, intent(out), optional :: LKMAX
    integer, intent(out), optional :: IMAXG
    integer, intent(out), optional :: JMAXG
    integer, intent(out), optional :: num_tile
    integer, intent(out), optional :: tile_id(:)

    integer :: i

    !---------------------------------------------------------------------------

    if ( dom_id < 1 .or. dom_id > num_dom ) then
       LOG_ERROR("COMM_CARTESC_NEST_domina_shape",*) 'domain id is invalid: ', dom_id
       call PRC_abort
    end if

    if ( present(KMAX) ) &
         KMAX = dom_info(dom_id)%KMAX

    if ( present(LKMAX) ) &
         LKMAX = dom_info(dom_id)%LKMAX

    if ( present(IMAXG) ) &
         IMAXG = dom_info(dom_id)%IMAX * dom_info(dom_id)%tile_num_x

    if ( present(JMAXG) ) &
         JMAXG = dom_info(dom_id)%JMAX * dom_info(dom_id)%tile_num_y

    if ( present(num_tile) ) &
         NUM_TILE = dom_info(dom_id)%tile_num_x * dom_info(dom_id)%tile_num_y

    if ( present(tile_id) ) then
       do i = 1, min( size(tile_id), size(dom_info(dom_id)%tile_id) )
          tile_id(i) = dom_info(dom_id)%tile_id(i)
       end do
    end if

    return
  end subroutine COMM_CARTESC_NEST_parent_info

  !-----------------------------------------------------------------------------
  !> Return shape of ParentDomain at the specified rank (for offline)
  !  including definition array size with BND or not in Parent domain
  subroutine COMM_CARTESC_NEST_domain_shape ( &
       tilei,     &
       tilej,     &
       cxs, cxe,  &
       cys, cye,  &
       pxs, pxe,  &
       pys, pye,  &
       dom_id,    &
       iloc,      &
       xstg, ystg )
    use scale_prc, only: &
       PRC_abort
    implicit none

    integer, intent(out) :: tilei, tilej
    integer, intent(out) :: cxs, cxe, cys, cye
    integer, intent(out) :: pxs, pxe, pys, pye

    integer, intent(in) :: dom_id
    integer, intent(in) :: iloc   ! rank number; start from 1

    logical, intent(in), optional :: xstg
    logical, intent(in), optional :: ystg

    type(domain_info), pointer :: dinfo

    integer :: hdl = 1      ! handler number
    integer :: rank
    integer :: xloc,  yloc
    integer :: xlocg, ylocg ! location over whole domain
    logical :: xstg_, ystg_
    !---------------------------------------------------------------------------

    if ( dom_id < 1 .or. dom_id > num_dom ) then
       LOG_ERROR("COMM_CARTESC_NEST_domina_shape",*) 'domain id is invalid: ', dom_id
       call PRC_abort
    end if

    if ( present(xstg) ) then
       xstg_ = xstg
    else
       xstg_ = .false.
    end if
    if ( present(ystg) ) then
       ystg_ = ystg
    else
       ystg_ = .false.
    end if

    dinfo => dom_info(dom_id)

    rank  = dinfo%tile_id(iloc)
    xloc  = mod( iloc-1, dinfo%tile_num_x ) + 1
    yloc  = int( real(iloc-1) / real(dinfo%tile_num_x) ) + 1
    xlocg = mod( rank, dinfo%prc_num_x ) + 1
    ylocg = int( real(rank) / real(dinfo%prc_num_x) ) + 1
    tilei = dinfo%IMAX
    tilej = dinfo%JMAX

    cxs   = tilei * (xloc-1) + 1
    cxe   = tilei * xloc
    cys   = tilej * (yloc-1) + 1
    cye   = tilej * yloc
    pxs   = 1
    pxe   = tilei
    pys   = 1
    pye   = tilej

    if ( .not. dinfo%periodic_x )  then
       if ( xlocg == 1 ) then ! BND_W
          tilei = tilei + dinfo%IHALO
          pxs = pxs + dinfo%IHALO
          pxe = pxe + dinfo%IHALO
       endif
       if ( xlocg == dinfo%prc_num_x ) then ! BND_E
          tilei = tilei + dinfo%IHALO
       endif

       if ( xstg_ ) then ! staggarded grid
          if ( xlocg == 1 ) then ! BND_W
             tilei = tilei + 1
             if ( dinfo%IHALO > 0 ) then
                pxs = pxs - 1
             else
                pxe = pxe + 1
             end if
          else
             cxs = cxs + 1
          end if
          cxe = cxe + 1
       end if
    end if

    if ( .not. dinfo%periodic_y ) then
       if ( ylocg == 1 ) then ! BND_S
          tilej = tilej + dinfo%JHALO
          pys = pys + dinfo%JHALO
          pye = pye + dinfo%JHALO
       endif
       if ( ylocg == dinfo%prc_num_y ) then ! BND_N
          tilej = tilej + dinfo%JHALO
       endif

       if ( ystg_ ) then ! staggarded grid
          if ( ylocg == 1 ) then ! BND_W
             tilej = tilej + 1
             if ( dinfo%JHALO > 0 ) then
                pys = pys - 1
             else
                pye = pye + 1
             end if
          else
             cys = cys + 1
          end if
          cye = cye + 1
       end if
    end if

    return
  end subroutine COMM_CARTESC_NEST_domain_shape

  !-----------------------------------------------------------------------------
  !> Get parent domain size
  subroutine COMM_CARTESC_NEST_parentsize( &
       HANDLE )
    use scale_prc, only: &
       PRC_abort,            &
       PRC_nprocs,           &
       PRC_myrank,           &
       PRC_IsMaster,         &
       PRC_INTERCOMM_PARENT, &
       PRC_INTERCOMM_CHILD
    use scale_prc_cartesC, only: &
       PRC_NUM_X,   &
       PRC_NUM_Y
    use scale_time, only: &
       TIME_NSTEP, &
       TIME_DTSEC
    use scale_comm_cartesC, only: &
       COMM_bcast
    use scale_atmos_hydrometeor, only: &
       N_HYD
    implicit none

    integer, intent(in) :: HANDLE !< id number of nesting relation in this process target

    real(RP) :: buffer
    integer, parameter :: ileng = 10
    integer  :: datapack(ileng)

    integer                :: QA_OTHERSIDE
    character(len=H_SHORT) :: MP_TYPE_OTHERSIDE

    integer :: ireq1, ireq2, ireq3, ireq4, ireq5, ireq6
    integer :: ierr1, ierr2, ierr3, ierr4, ierr5, ierr6
    integer :: istatus(MPI_STATUS_SIZE)
    integer :: tag
    !---------------------------------------------------------------------------

    if( .NOT. USE_NESTING ) return

    tag   = INTERCOMM_ID(HANDLE) * 100

    if ( COMM_CARTESC_NEST_Filiation( INTERCOMM_ID(HANDLE) ) > 0 ) then

       !##### parent ####

       datapack( 1) = PRC_NUM_X
       datapack( 2) = PRC_NUM_Y
       datapack( 3) = KMAX
       datapack( 4) = KHALO
       datapack( 5) = IMAX
       datapack( 6) = IHALO
       datapack( 7) = JMAX
       datapack( 8) = JHALO
       datapack( 9) = TIME_NSTEP
       datapack(10) = ONLINE_SEND_QA
       buffer       = TIME_DTSEC

       if ( PRC_IsMaster ) then
          ! from daughter to parent
          call MPI_IRECV(ONLINE_DAUGHTER_NPROCS, 1, MPI_INTEGER, PRC_myrank, tag+3, PRC_INTERCOMM_CHILD, ireq4, ierr4)
          call MPI_IRECV(QA_OTHERSIDE, 1, MPI_INTEGER, PRC_myrank, tag+4, PRC_INTERCOMM_CHILD, ireq5, ierr5)
          call MPI_IRECV(MP_TYPE_OTHERSIDE, H_SHORT, MPI_CHARACTER, PRC_myrank, tag+5, PRC_INTERCOMM_CHILD, ireq6, ierr6)

          ! from parent to daughter
          call MPI_ISEND(datapack, ileng, MPI_INTEGER, PRC_myrank, tag, PRC_INTERCOMM_CHILD, ireq1, ierr1)
          call MPI_ISEND(buffer, 1, MPI_DOUBLE_PRECISION, PRC_myrank, tag+1, PRC_INTERCOMM_CHILD, ireq2, ierr2)
          call MPI_ISEND(MP_TYPE, H_SHORT, MPI_CHARACTER, PRC_myrank, tag+2, PRC_INTERCOMM_CHILD, ireq3, ierr3)
          call MPI_WAIT(ireq1, istatus, ierr1)
          call MPI_WAIT(ireq2, istatus, ierr2)
          call MPI_WAIT(ireq3, istatus, ierr3)

          call MPI_WAIT(ireq4, istatus, ierr4)
          call MPI_WAIT(ireq5, istatus, ierr5)
          call MPI_WAIT(ireq6, istatus, ierr6)
       end if

       call COMM_bcast(ONLINE_DAUGHTER_NPROCS)
       call COMM_bcast(QA_OTHERSIDE)
       call COMM_bcast(MP_TYPE_OTHERSIDE)

       if ( MP_TYPE == MP_TYPE_OTHERSIDE .and. ONLINE_SEND_QA == QA_OTHERSIDE ) then
          ONLINE_SEND_DIAGQHYD = .false.
       else if ( MP_TYPE_OTHERSIDE == "DRY" .or. MP_TYPE == "DRY" ) then
          ONLINE_SEND_QA = 0
          ONLINE_SEND_DIAGQHYD = .false.
       else if ( MP_TYPE_OTHERSIDE == "QV" .or. MP_TYPE == "QV" ) then
          ONLINE_SEND_QA = 1
          ONLINE_SEND_DIAGQHYD = .false.
       else
          LOG_INFO("COMM_CARTESC_NEST_parentsize",*) 'Hydrometeor will be diagnosed on children side'
          LOG_INFO("COMM_CARTESC_NEST_parentsize",*) 'MP type      (remote,local) = ', trim(MP_TYPE_OTHERSIDE), ", ", trim(MP_TYPE)
          LOG_INFO("COMM_CARTESC_NEST_parentsize",*) 'Number of QA (remote,local) = ', QA_OTHERSIDE, ONLINE_SEND_QA
          ONLINE_SEND_QA = N_HYD + 1 ! QV + hydrometeors
          ONLINE_SEND_DIAGQHYD = .true.
       endif


    elseif( COMM_CARTESC_NEST_Filiation( INTERCOMM_ID(HANDLE) ) < 0 ) then

       !##### child ####

       if ( PRC_IsMaster ) then
          ! from parent to daughter
          call MPI_IRECV(datapack, ileng, MPI_INTEGER, PRC_myrank, tag, PRC_INTERCOMM_PARENT, ireq1, ierr1)
          call MPI_IRECV(buffer, 1, MPI_DOUBLE_PRECISION, PRC_myrank, tag+1, PRC_INTERCOMM_PARENT, ireq2, ierr2)
          call MPI_IRECV(MP_TYPE_OTHERSIDE, H_SHORT, MPI_CHARACTER, PRC_myrank, tag+2, PRC_INTERCOMM_PARENT, ireq3, ierr3)

          ! from daughter to parent
          call MPI_ISEND(PRC_nprocs, 1, MPI_INTEGER, PRC_myrank, tag+3, PRC_INTERCOMM_PARENT, ireq4, ierr4)
          call MPI_ISEND(ONLINE_RECV_QA, 1, MPI_INTEGER, PRC_myrank, tag+4, PRC_INTERCOMM_PARENT, ireq5, ierr5)
          call MPI_ISEND(MP_TYPE, H_SHORT, MPI_CHARACTER, PRC_myrank, tag+5, PRC_INTERCOMM_PARENT, ireq6, ierr6)

          call MPI_WAIT(ireq4, istatus, ierr4)
          call MPI_WAIT(ireq5, istatus, ierr5)
          call MPI_WAIT(ireq6, istatus, ierr6)

          call MPI_WAIT(ireq1, istatus, ierr1)
          call MPI_WAIT(ireq2, istatus, ierr2)
          call MPI_WAIT(ireq3, istatus, ierr3)

       endif
       call COMM_bcast(ileng, datapack)
       call COMM_bcast(buffer)
       call COMM_bcast(MP_TYPE_OTHERSIDE)

       dom_info(I_PARENT)%prc_num_x = datapack( 1)
       dom_info(I_PARENT)%prc_num_y = datapack( 2)
       dom_info(I_PARENT)%KMAX      = datapack( 3)
       dom_info(I_PARENT)%KHALO     = datapack( 4)
       dom_info(I_PARENT)%IMAX      = datapack( 5)
       dom_info(I_PARENT)%IHALO     = datapack( 6)
       dom_info(I_PARENT)%JMAX      = datapack( 7)
       dom_info(I_PARENT)%JHALO     = datapack( 8)
       ONLINE_PARENT_NSTEP          = datapack( 9)
       QA_OTHERSIDE                 = datapack(10)
       ONLINE_PARENT_DTSEC          = buffer

       if ( MP_TYPE == MP_TYPE_OTHERSIDE .and. ONLINE_RECV_QA == QA_OTHERSIDE ) then
          ONLINE_RECV_DIAGQHYD = .false.
       else if ( MP_TYPE == "DRY" ) then
          ONLINE_RECV_QA = 0
          ONLINE_RECV_DIAGQHYD = .false.
       else if ( MP_TYPE == "QV" ) then
          ONLINE_RECV_QA = 1
          ONLINE_RECV_DIAGQHYD = .false.
       else
          LOG_INFO("COMM_CARTESC_NEST_parentsize",*) 'Hydrometeor will be diagnosed on this side'
          LOG_INFO("COMM_CARTESC_NEST_parentsize",*) 'MP type      (remote,local) = ', trim(MP_TYPE_OTHERSIDE), ", ", trim(MP_TYPE)
          LOG_INFO("COMM_CARTESC_NEST_parentsize",*) 'Number of QA (remote,local) = ', QA_OTHERSIDE, ONLINE_RECV_QA
          ONLINE_RECV_QA = N_HYD + 1 ! QV + hydrometeors
          ONLINE_RECV_DIAGQHYD = .true.
       endif

    else
       LOG_ERROR("COMM_CARTESC_NEST_parentsize",*) '[COMM_CARTESC_NEST_parentsize] internal error'
       call PRC_abort
    endif

    return
  end subroutine COMM_CARTESC_NEST_parentsize

  !-----------------------------------------------------------------------------
  !> Get parent latlon catalogue
  subroutine COMM_CARTESC_NEST_catalogue( &
       HANDLE )
    use scale_prc, only: &
       PRC_abort,            &
       PRC_nprocs,           &
       PRC_myrank,           &
       PRC_IsMaster,         &
       PRC_INTERCOMM_PARENT, &
       PRC_INTERCOMM_CHILD
    use scale_comm_cartesC, only: &
       COMM_datatype,  &
       COMM_bcast
    use scale_atmos_grid_cartesC_real, only: &
       ATMOS_GRID_CARTESC_REAL_DOMAIN_CATALOGUE
    implicit none

    integer, intent(in) :: HANDLE !< id number of nesting relation in this process target

    integer :: nprocs
    integer :: ireq, ierr, ileng
    integer :: istatus(MPI_STATUS_SIZE)
    integer :: tag
    !---------------------------------------------------------------------------

    if( .NOT. USE_NESTING ) return

    tag = INTERCOMM_ID(HANDLE) * 100

    if ( COMM_CARTESC_NEST_Filiation( INTERCOMM_ID(HANDLE) ) > 0 ) then

       !##### parent ####

       ileng = PRC_nprocs * 2 * 2

       if ( PRC_IsMaster ) then
          call MPI_ISEND(ATMOS_GRID_CARTESC_REAL_DOMAIN_CATALOGUE, ileng, COMM_datatype, PRC_myrank, tag, PRC_INTERCOMM_CHILD, ireq, ierr)
          call MPI_WAIT(ireq, istatus, ierr)
       endif

    elseif( COMM_CARTESC_NEST_Filiation( INTERCOMM_ID(HANDLE) ) < 0 ) then

       !##### child ####

       nprocs = dom_info(I_PARENT)%prc_num_x * dom_info(I_PARENT)%prc_num_y
       ileng = nprocs * 2 * 2

       if ( PRC_IsMaster ) then
          call MPI_IRECV(dom_info(I_PARENT)%latlon_catalogue, ileng, COMM_datatype, PRC_myrank, tag, PRC_INTERCOMM_PARENT, ireq, ierr)
          call MPI_WAIT(ireq, istatus, ierr)
       endif
       call COMM_bcast( nprocs, 2, 2, dom_info(I_PARENT)%latlon_catalogue )

    else
       LOG_ERROR("COMM_CARTESC_NEST_catalogue",*) 'internal error'
       call PRC_abort
    endif

    return
  end subroutine COMM_CARTESC_NEST_catalogue

  !-----------------------------------------------------------------------------
  !> Check Communication Inter-domains
  subroutine COMM_CARTESC_NEST_ping( &
       HANDLE )
    use scale_prc, only: &
       PRC_abort,            &
       PRC_myrank,           &
       PRC_IsMaster,         &
       PRC_INTERCOMM_PARENT, &
       PRC_INTERCOMM_CHILD
    use scale_comm_cartesC, only: &
       COMM_bcast
    implicit none

    integer, intent(in) :: HANDLE !< id number of nesting relation in this process target

    integer :: ping, pong
    integer :: ireq1, ireq2, ierr1, ierr2
    integer :: istatus(MPI_STATUS_SIZE)
    integer :: tag
    logical :: ping_error
    !---------------------------------------------------------------------------

    if( .NOT. USE_NESTING ) return

    tag        = INTERCOMM_ID(HANDLE) * 100
    ping_error = .false.

    if ( COMM_CARTESC_NEST_Filiation( INTERCOMM_ID(HANDLE) ) > 0 ) then

       !##### parent ####

       ping = ONLINE_DOMAIN_NUM
       pong = 0

       if ( PRC_IsMaster ) then
          call MPI_IRECV(pong, 1, MPI_INTEGER, PRC_myrank, tag+2, PRC_INTERCOMM_CHILD, ireq2, ierr2)
          call MPI_ISEND(ping, 1, MPI_INTEGER, PRC_myrank, tag+1, PRC_INTERCOMM_CHILD, ireq1, ierr1)
          call MPI_WAIT(ireq1, istatus, ierr1)
          call MPI_WAIT(ireq2, istatus, ierr2)
       endif

       call COMM_bcast(pong)

       if ( pong /= INTERCOMM_ID(HANDLE)+1 ) ping_error = .true.

    elseif( COMM_CARTESC_NEST_Filiation( INTERCOMM_ID(HANDLE) ) < 0 ) then

       !##### child ####

       ping = ONLINE_DOMAIN_NUM
       pong = 0

       if ( PRC_IsMaster ) then
          call MPI_IRECV(pong, 1, MPI_INTEGER, PRC_myrank, tag+1, PRC_INTERCOMM_PARENT, ireq2, ierr2)
          call MPI_ISEND(ping, 1, MPI_INTEGER, PRC_myrank, tag+2, PRC_INTERCOMM_PARENT, ireq1, ierr1)
          call MPI_WAIT(ireq1, istatus, ierr1)
          call MPI_WAIT(ireq2, istatus, ierr2)
       endif

       call COMM_bcast(pong)

       if ( pong /= INTERCOMM_ID(HANDLE) ) ping_error = .true.

    else
       LOG_ERROR("COMM_CARTESC_NEST_ping",*) 'internal error'
       call PRC_abort
    endif

    if ( ping_error ) then
       LOG_ERROR("COMM_CARTESC_NEST_ping",*) 'ping destination error'
       call PRC_abort
    endif

    return
  end subroutine COMM_CARTESC_NEST_ping

  !-----------------------------------------------------------------------------
  !> Inter-domain communication setup for nestdown
  subroutine COMM_CARTESC_NEST_setup_nestdown( &
       HANDLE )
    use scale_prc, only: &
       PRC_abort,            &
       PRC_myrank,           &
       PRC_IsMaster,         &
       PRC_INTERCOMM_PARENT, &
       PRC_INTERCOMM_CHILD
    use scale_prc, only: &
       PRC_nprocs
    use scale_comm_cartesC, only: &
       COMM_world, &
       COMM_bcast
    implicit none

    integer, intent(in) :: HANDLE !< id number of nesting relation in this process target

    integer, allocatable :: buffer_LIST   (:)
    integer, allocatable :: buffer_ALLLIST(:)

    integer :: parent_ka
    integer :: parent_ia
    integer :: parent_ja

    integer :: ireq, ierr, ileng
    integer :: istatus(MPI_STATUS_SIZE)
    integer :: tag, target_rank

    integer :: i, j, k
    !---------------------------------------------------------------------------

    if( .NOT. USE_NESTING ) return

    tag = INTERCOMM_ID(HANDLE) * 100

    if ( COMM_CARTESC_NEST_Filiation( INTERCOMM_ID(HANDLE) ) > 0 ) then

       !##### parent ####

       if ( PRC_IsMaster ) then
          call MPI_IRECV(COMM_CARTESC_NEST_TILE_ALLMAX_p, 1, MPI_INTEGER, PRC_myrank, tag+1, PRC_INTERCOMM_CHILD, ireq, ierr)
          call MPI_WAIT(ireq, istatus, ierr)
       endif
       call COMM_bcast(COMM_CARTESC_NEST_TILE_ALLMAX_p)

       allocate( COMM_CARTESC_NEST_TILE_LIST_p (COMM_CARTESC_NEST_TILE_ALLMAX_p,ONLINE_DAUGHTER_nprocs) )
       allocate( COMM_CARTESC_NEST_TILE_LIST_YP(COMM_CARTESC_NEST_TILE_ALLMAX_p*ONLINE_DAUGHTER_nprocs) )

       ileng = COMM_CARTESC_NEST_TILE_ALLMAX_p * ONLINE_DAUGHTER_nprocs
       if ( PRC_IsMaster ) then
          call MPI_IRECV(COMM_CARTESC_NEST_TILE_LIST_p, ileng, MPI_INTEGER, PRC_myrank, tag+2, PRC_INTERCOMM_CHILD, ireq, ierr)
          call MPI_WAIT(ireq, istatus, ierr)
       endif
       call COMM_bcast(COMM_CARTESC_NEST_TILE_ALLMAX_p, ONLINE_DAUGHTER_nprocs, COMM_CARTESC_NEST_TILE_LIST_p)

       COMM_CARTESC_NEST_TILE_LIST_YP(:) = -1

       k = 0
       do j = 1, ONLINE_DAUGHTER_nprocs
       do i = 1, COMM_CARTESC_NEST_TILE_ALLMAX_p
          if ( COMM_CARTESC_NEST_TILE_LIST_p(i,j) == PRC_myrank ) then
             k = k + 1
             COMM_CARTESC_NEST_TILE_LIST_YP(k) = j - 1  !rank number is started from 1
          endif
       enddo
       enddo
       NUM_YP = k

       LOG_INFO("COMM_CARTESC_NEST_setup_nestdown",'(A,I5,A,I5)') "[P]   Num YP =",NUM_YP,"  Num TILE(MAX) =",COMM_CARTESC_NEST_TILE_ALLMAX_p

       if ( PRC_IsMaster ) then
          call MPI_IRECV(ONLINE_DAUGHTER_USE_VELZ, 1, MPI_LOGICAL, PRC_myrank, tag+3, PRC_INTERCOMM_CHILD, ireq, ierr)
          call MPI_WAIT(ireq, istatus, ierr)
       endif
       call COMM_bcast(ONLINE_DAUGHTER_USE_VELZ)

       LOG_INFO("COMM_CARTESC_NEST_setup_nestdown",'(1x,A,L2)') 'NEST: ONLINE_DAUGHTER_USE_VELZ =', ONLINE_DAUGHTER_USE_VELZ

       if ( PRC_IsMaster ) then
          call MPI_IRECV(ONLINE_DAUGHTER_NO_ROTATE, 1, MPI_LOGICAL, PRC_myrank, tag+4, PRC_INTERCOMM_CHILD, ireq, ierr)
          call MPI_WAIT(ireq, istatus, ierr)
       endif
       call COMM_bcast(ONLINE_DAUGHTER_NO_ROTATE)

       if( ONLINE_NO_ROTATE .neqv. ONLINE_DAUGHTER_NO_ROTATE ) then
          LOG_ERROR("COMM_CARTESC_NEST_setup_nestdown",*) 'Flag of NO_ROTATE is not consistent with the child domain'
          LOG_ERROR_CONT(*) 'ONLINE_NO_ROTATE = ', ONLINE_NO_ROTATE
          LOG_ERROR_CONT(*) 'ONLINE_DAUGHTER_NO_ROTATE =', ONLINE_DAUGHTER_NO_ROTATE
          call PRC_abort
       endif
       LOG_INFO("COMM_CARTESC_NEST_setup_nestdown",'(1x,A,L2)') 'NEST: ONLINE_DAUGHTER_NO_ROTATE =', ONLINE_DAUGHTER_NO_ROTATE

       call COMM_CARTESC_NEST_importgrid_nestdown( HANDLE )

       do i = 1, NUM_YP
          target_rank = COMM_CARTESC_NEST_TILE_LIST_YP(i)
          call MPI_ISEND(i, 1, MPI_INTEGER, target_rank, tag+5, PRC_INTERCOMM_CHILD, ireq, ierr)
          call MPI_WAIT(ireq, istatus, ierr)
       enddo

       call MPI_BARRIER(PRC_INTERCOMM_CHILD, ierr)

    elseif( COMM_CARTESC_NEST_Filiation( INTERCOMM_ID(HANDLE) ) < 0 ) then

       !##### child ####

       COMM_CARTESC_NEST_TILE_ALL = size( dom_info(I_PARENT)%TILE_ID(:) ) ! should be equal to "NEST_TILE_NUM_X*NEST_TILE_NUM_Y"
       call MPI_Allreduce( COMM_CARTESC_NEST_TILE_ALL,      &
                           COMM_CARTESC_NEST_TILE_ALLMAX_d, &
                           1,                  &
                           MPI_INTEGER,        &
                           MPI_MAX,            &
                           COMM_world,         &
                           ierr                )
       LOG_INFO("COMM_CARTESC_NEST_setup_nestdown",'(A,I5,A,I5)') "[D]   Num YP =",COMM_CARTESC_NEST_TILE_ALL,"  Num TILE(MAX) =",COMM_CARTESC_NEST_TILE_ALLMAX_d

       if ( PRC_IsMaster ) then
          call MPI_ISEND(COMM_CARTESC_NEST_TILE_ALLMAX_d, 1, MPI_INTEGER, PRC_myrank, tag+1, PRC_INTERCOMM_PARENT, ireq, ierr)
          call MPI_WAIT(ireq, istatus, ierr)
       endif

       parent_ka = dom_info(I_PARENT)%KMAX + dom_info(I_PARENT)%KHALO * 2
       parent_ia = dom_info(I_PARENT)%IMAX + dom_info(I_PARENT)%IHALO * 2
       parent_ja = dom_info(I_PARENT)%JMAX + dom_info(I_PARENT)%JHALO * 2

       max_isu = 4 + ONLINE_RECV_QA
       if ( ONLINE_USE_VELZ ) max_isu = max_isu + 1
       max_isu = COMM_CARTESC_NEST_TILE_ALL * max_isu
       allocate( recvbuf_3D( parent_ka, parent_ia, parent_ja, max_isu ) )

       allocate( buffer_LIST   (COMM_CARTESC_NEST_TILE_ALLMAX_d)            )
       allocate( buffer_ALLLIST(COMM_CARTESC_NEST_TILE_ALLMAX_d*PRC_nprocs)   )
       allocate( COMM_CARTESC_NEST_TILE_LIST_d(COMM_CARTESC_NEST_TILE_ALLMAX_d,PRC_nprocs) )

       do i = 1, COMM_CARTESC_NEST_TILE_ALLMAX_d
          if ( i <= COMM_CARTESC_NEST_TILE_ALL ) then
             buffer_LIST(i) = dom_info(I_PARENT)%TILE_ID(i)
          else
             buffer_LIST(i) = -1
          endif
       enddo

       ileng = COMM_CARTESC_NEST_TILE_ALLMAX_d
       call MPI_Allgather( buffer_LIST(:),     &
                           ileng,              &
                           MPI_INTEGER,        &
                           buffer_ALLLIST(:),  &
                           ileng,              &
                           MPI_INTEGER,        &
                           COMM_world,         &
                           ierr                )
       k = 1
       do j = 1, PRC_nprocs
       do i = 1, COMM_CARTESC_NEST_TILE_ALLMAX_d
          COMM_CARTESC_NEST_TILE_LIST_d(i,j) = buffer_ALLLIST(k)
          k = k + 1
       enddo
       enddo

       deallocate( buffer_LIST    )
       deallocate( buffer_ALLLIST )

       ileng = COMM_CARTESC_NEST_TILE_ALLMAX_d * PRC_nprocs
       if ( PRC_IsMaster ) then
          call MPI_ISEND(COMM_CARTESC_NEST_TILE_LIST_d, ileng, MPI_INTEGER, PRC_myrank, tag+2, PRC_INTERCOMM_PARENT, ireq, ierr)
          call MPI_WAIT(ireq, istatus, ierr)
       endif

       if ( PRC_IsMaster ) then
          call MPI_ISEND(ONLINE_USE_VELZ, 1, MPI_LOGICAL, PRC_myrank, tag+3, PRC_INTERCOMM_PARENT, ireq, ierr)
          call MPI_WAIT(ireq, istatus, ierr)
       endif

       if ( PRC_IsMaster ) then
          call MPI_ISEND(ONLINE_NO_ROTATE, 1, MPI_LOGICAL, PRC_myrank, tag+4, PRC_INTERCOMM_PARENT, ireq, ierr)
          call MPI_WAIT(ireq, istatus, ierr)
       endif
       call COMM_bcast(ONLINE_DAUGHTER_NO_ROTATE)

       call COMM_CARTESC_NEST_importgrid_nestdown( HANDLE )

       do i = 1, COMM_CARTESC_NEST_TILE_ALL
          target_rank = COMM_CARTESC_NEST_TILE_LIST_d(i,PRC_myrank+1)
          call MPI_IRECV( call_order(i), 1, MPI_INTEGER, target_rank, tag+5, PRC_INTERCOMM_PARENT, ireq, ierr )
          call MPI_WAIT(ireq, istatus, ierr)
       enddo

       call MPI_BARRIER(PRC_INTERCOMM_PARENT, ierr)
    else
       LOG_ERROR("COMM_CARTESC_NEST_setup_nestdown",*) 'internal error'
       call PRC_abort
    endif

    if( NUM_YP * 16 > max_rq .OR. COMM_CARTESC_NEST_TILE_ALL * 16 > max_rq ) then ! 16 = dyn:5 + qtrc:11
       LOG_ERROR("COMM_CARTESC_NEST_setup_nestdown",*) 'internal error (overflow number of ireq)'
       LOG_ERROR_CONT(*) 'NUM_YP x 16        = ', NUM_YP * 16
       LOG_ERROR_CONT(*) 'COMM_CARTESC_NEST_TILE_ALL x 16 = ', COMM_CARTESC_NEST_TILE_ALL * 16
       LOG_ERROR_CONT(*) 'max_rq             = ', max_rq
       call PRC_abort
    endif

    return
  end subroutine COMM_CARTESC_NEST_setup_nestdown

  !-----------------------------------------------------------------------------
  !> Grid Data transfer from parent to daughter: nestdown
  subroutine COMM_CARTESC_NEST_importgrid_nestdown( &
       HANDLE )
    use scale_prc, only: &
       PRC_myrank,           &
       PRC_abort,            &
       PRC_INTERCOMM_PARENT, &
       PRC_INTERCOMM_CHILD
    use scale_atmos_grid_cartesC_real, only: &
       ATMOS_GRID_CARTESC_REAL_LON,   &
       ATMOS_GRID_CARTESC_REAL_LAT,   &
       ATMOS_GRID_CARTESC_REAL_LONUY, &
       ATMOS_GRID_CARTESC_REAL_LONXV, &
       ATMOS_GRID_CARTESC_REAL_LATUY, &
       ATMOS_GRID_CARTESC_REAL_LATXV, &
       ATMOS_GRID_CARTESC_REAL_CZ,    &
       ATMOS_GRID_CARTESC_REAL_FZ
    use scale_comm_cartesC, only: &
       COMM_datatype
    implicit none

    integer, intent(in) :: HANDLE !< id number of nesting relation in this process target

    integer :: parent_KA
    integer :: parent_IA, parent_IS, parent_IE, parent_IMAX
    integer :: parent_JA, parent_JS, parent_JE, parent_JMAX

    integer  :: ierr, ileng
    integer  :: istatus(MPI_STATUS_SIZE)
    integer  :: tag, tagbase, target_rank
    integer  :: rq_str, rq_end, rq_tot

    integer  :: xloc, yloc
    integer  :: xs, xe
    integer  :: ys, ye

    real(RP) :: max_ref, max_loc

    real(RP), allocatable :: sendbuf_2D(:,:,:)
    real(RP), allocatable :: sendbuf_3D(:,:,:,:)
    real(RP), allocatable :: recvbuf_2D(:,:,:)

    integer  :: i, k, rq
    !---------------------------------------------------------------------------

    if( .NOT. USE_NESTING ) return

    tagbase = INTERCOMM_ID(HANDLE) * 100
    rq      = 0

    if ( COMM_CARTESC_NEST_Filiation( INTERCOMM_ID(HANDLE) ) > 0 ) then

       !##### parent [send issue] #####

       allocate( sendbuf_2D(     IA, JA, 4 ) )
       allocate( sendbuf_3D( KA, IA, JA, 1 ) )

       do i = 1, NUM_YP
          ! send data to multiple daughter processes
          target_rank = COMM_CARTESC_NEST_TILE_LIST_YP(i)

          rq_str = rq + 1

          rq = rq + 1
          ileng = IA * JA
          tag   = tagbase + tag_lon
          call MPI_ISEND(ATMOS_GRID_CARTESC_REAL_LON, ileng, COMM_datatype, target_rank, tag, PRC_INTERCOMM_CHILD, ireq_p(rq), ierr)

          rq = rq + 1
          ileng = IA * JA
          tag   = tagbase + tag_lat
          call MPI_ISEND(ATMOS_GRID_CARTESC_REAL_LAT, ileng, COMM_datatype, target_rank, tag, PRC_INTERCOMM_CHILD, ireq_p(rq), ierr)

          sendbuf_2D(:,:,1) = ATMOS_GRID_CARTESC_REAL_LONUY(1:IA,1:JA)
          rq = rq + 1
          ileng = IA * JA
          tag   = tagbase + tag_lonuy
          call MPI_ISEND(sendbuf_2D(:,:,1), ileng, COMM_datatype, target_rank, tag, PRC_INTERCOMM_CHILD, ireq_p(rq), ierr)

          sendbuf_2D(:,:,2) = ATMOS_GRID_CARTESC_REAL_LATUY(1:IA,1:JA)
          rq = rq + 1
          ileng = IA * JA
          tag   = tagbase + tag_latuy
          call MPI_ISEND(sendbuf_2D(:,:,2), ileng, COMM_datatype, target_rank, tag, PRC_INTERCOMM_CHILD, ireq_p(rq), ierr)

          sendbuf_2D(:,:,3) = ATMOS_GRID_CARTESC_REAL_LONXV(1:IA,1:JA)
          rq = rq + 1
          ileng = IA * JA
          tag   = tagbase + tag_lonxv
          call MPI_ISEND(sendbuf_2D(:,:,3), ileng, COMM_datatype, target_rank, tag, PRC_INTERCOMM_CHILD, ireq_p(rq), ierr)

          sendbuf_2D(:,:,4) = ATMOS_GRID_CARTESC_REAL_LATXV(1:IA,1:JA)
          rq = rq + 1
          ileng = IA * JA
          tag   = tagbase + tag_latxv
          call MPI_ISEND(sendbuf_2D(:,:,4), ileng, COMM_datatype, target_rank, tag, PRC_INTERCOMM_CHILD, ireq_p(rq), ierr)

          rq = rq + 1
          ileng = KA * IA * JA
          tag   = tagbase + tag_cz
          call MPI_ISEND(ATMOS_GRID_CARTESC_REAL_CZ, ileng, COMM_datatype, target_rank, tag, PRC_INTERCOMM_CHILD, ireq_p(rq), ierr)

          sendbuf_3D(:,:,:,1) = ATMOS_GRID_CARTESC_REAL_FZ(1:,:,:)
          rq = rq + 1
          ileng = KA * IA * JA
          tag   = tagbase + tag_fz
          call MPI_ISEND(sendbuf_3D(:,:,:,1), ileng, COMM_datatype, target_rank, tag, PRC_INTERCOMM_CHILD, ireq_p(rq), ierr)

          rq_end = rq
          rq_tot = rq_end - rq_str + 1

          call COMM_CARTESC_NEST_waitall( rq_tot, ireq_p(rq_str:rq_end) )
       enddo

       deallocate( sendbuf_2D )
       deallocate( sendbuf_3D )

    elseif( COMM_CARTESC_NEST_Filiation( INTERCOMM_ID(HANDLE) ) < 0 ) then

       !##### child [recv & wait issue] #####

       parent_KA = dom_info(I_PARENT)%KMAX + dom_info(I_PARENT)%KHALO * 2

       parent_IMAX = dom_info(I_PARENT)%IMAX
       parent_IA = parent_IMAX + dom_info(I_PARENT)%IHALO * 2
       parent_IS = dom_info(I_PARENT)%IHALO + 1
       parent_IE = parent_IMAX + dom_info(I_PARENT)%IHALO

       parent_JMAX = dom_info(I_PARENT)%JMAX
       parent_JA = parent_JMAX + dom_info(I_PARENT)%JHALO * 2
       parent_JS = dom_info(I_PARENT)%JHALO + 1
       parent_JE = parent_JMAX + dom_info(I_PARENT)%JHALO

       allocate( recvbuf_2D( PARENT_IA, PARENT_JA, 6 ) )

       do i = 1, COMM_CARTESC_NEST_TILE_ALL
          ! receive data from multiple parent tiles
          target_rank = COMM_CARTESC_NEST_TILE_LIST_d(i,PRC_myrank+1)

          xloc = mod( i-1, dom_info(I_PARENT)%tile_num_x ) + 1
          yloc = int( real(i-1) / real(dom_info(I_PARENT)%tile_num_x) ) + 1

          xs = PARENT_IMAX * (xloc-1) + 1
          xe = PARENT_IMAX *  xloc
          ys = PARENT_JMAX * (yloc-1) + 1
          ye = PARENT_JMAX *  yloc

          rq_str = rq + 1

          rq = rq + 1
          ileng = PARENT_IA * PARENT_JA
          tag   = tagbase + tag_lon
          call MPI_IRECV(recvbuf_2D(:,:,tag_lon), ileng, COMM_datatype, target_rank, tag, PRC_INTERCOMM_PARENT, ireq_d(rq), ierr)

          rq = rq + 1
          ileng = PARENT_IA * PARENT_JA
          tag   = tagbase + tag_lat
          call MPI_IRECV(recvbuf_2D(:,:,tag_lat), ileng, COMM_datatype, target_rank, tag, PRC_INTERCOMM_PARENT, ireq_d(rq), ierr)

          rq = rq + 1
          ileng = PARENT_IA * PARENT_JA
          tag   = tagbase + tag_lonuy
          call MPI_IRECV(recvbuf_2D(:,:,tag_lonuy), ileng, COMM_datatype, target_rank, tag, PRC_INTERCOMM_PARENT, ireq_d(rq), ierr)

          rq = rq + 1
          ileng = PARENT_IA * PARENT_JA
          tag   = tagbase + tag_latuy
          call MPI_IRECV(recvbuf_2D(:,:,tag_latuy), ileng, COMM_datatype, target_rank, tag, PRC_INTERCOMM_PARENT, ireq_d(rq), ierr)

          rq = rq + 1
          ileng = PARENT_IA * PARENT_JA
          tag   = tagbase + tag_lonxv
          call MPI_IRECV(recvbuf_2D(:,:,tag_lonxv), ileng, COMM_datatype, target_rank, tag, PRC_INTERCOMM_PARENT, ireq_d(rq), ierr)

          rq = rq + 1
          ileng = PARENT_IA * PARENT_JA
          tag   = tagbase + tag_latxv
          call MPI_IRECV(recvbuf_2D(:,:,tag_latxv), ileng, COMM_datatype, target_rank, tag, PRC_INTERCOMM_PARENT, ireq_d(rq), ierr)

          rq = rq + 1
          ileng = PARENT_KA * PARENT_IA * PARENT_JA
          tag   = tagbase + tag_cz
          call MPI_IRECV(recvbuf_3D(:,:,:,tag_cz), ileng, COMM_datatype, target_rank, tag, PRC_INTERCOMM_PARENT, ireq_d(rq), ierr)

          rq = rq + 1
          ileng = PARENT_KA * PARENT_IA * PARENT_JA
          tag   = tagbase + tag_fz
          call MPI_IRECV(recvbuf_3D(:,:,:,tag_fz), ileng, COMM_datatype, target_rank, tag, PRC_INTERCOMM_PARENT, ireq_d(rq), ierr)

          rq_end = rq
          rq_tot = rq_end - rq_str + 1

          call COMM_CARTESC_NEST_waitall( rq_tot, ireq_d(rq_str:rq_end) )

          buffer_ref_LON  (xs:xe,ys:ye)  = recvbuf_2D(PARENT_IS:PARENT_IE,PARENT_JS:PARENT_JE,tag_lon  )
          buffer_ref_LAT  (xs:xe,ys:ye)  = recvbuf_2D(PARENT_IS:PARENT_IE,PARENT_JS:PARENT_JE,tag_lat  )
          buffer_ref_LONUY(xs:xe,ys:ye)  = recvbuf_2D(PARENT_IS:PARENT_IE,PARENT_JS:PARENT_JE,tag_lonuy)
          buffer_ref_LATUY(xs:xe,ys:ye)  = recvbuf_2D(PARENT_IS:PARENT_IE,PARENT_JS:PARENT_JE,tag_latuy)
          buffer_ref_LONXV(xs:xe,ys:ye)  = recvbuf_2D(PARENT_IS:PARENT_IE,PARENT_JS:PARENT_JE,tag_lonxv)
          buffer_ref_LATXV(xs:xe,ys:ye)  = recvbuf_2D(PARENT_IS:PARENT_IE,PARENT_JS:PARENT_JE,tag_latxv)

          do k = 1, PARENT_KA
             buffer_ref_CZ(k,xs:xe,ys:ye)  = recvbuf_3D(k,PARENT_IS:PARENT_IE,PARENT_JS:PARENT_JE,tag_cz)
             buffer_ref_FZ(k,xs:xe,ys:ye)  = recvbuf_3D(k,PARENT_IS:PARENT_IE,PARENT_JS:PARENT_JE,tag_fz)
          enddo
       enddo

       ! check domain compatibility
       max_ref = maxval( buffer_ref_FZ(:,:,:) )
       max_loc = maxval( ATMOS_GRID_CARTESC_REAL_FZ(KS-1:KE,:,:) ) ! HALO + 1
       if ( max_ref < max_loc ) then
          LOG_ERROR("COMM_CARTESC_NEST_importgrid_nestdown",*) 'REQUESTED DOMAIN IS TOO MUCH BROAD'
          LOG_ERROR_CONT(*) '-- VERTICAL direction over the limit'
          LOG_ERROR_CONT(*) '-- reference max: ', max_ref
          LOG_ERROR_CONT(*) '--     local max: ', max_loc
          call PRC_abort
       endif

       deallocate( recvbuf_2D )

    else
       LOG_ERROR("COMM_CARTESC_NEST_importgrid_nestdown",*) 'internal error'
       call PRC_abort
    endif

    return
  end subroutine COMM_CARTESC_NEST_importgrid_nestdown

  !-----------------------------------------------------------------------------
  !> Boundary data transfer from parent to daughter: nestdown (parent side)
  subroutine COMM_CARTESC_NEST_nestdown_send( &
       DENS_send, &
       MOMZ_send, &
       MOMX_send, &
       MOMY_send, &
       RHOT_send, &
       QTRC_send  )
    use scale_prc, only: &
       PRC_abort, &
       PRC_INTERCOMM_PARENT, &
       PRC_INTERCOMM_CHILD
    use scale_comm_cartesC, only: &
       COMM_vars8, &
       COMM_wait
    use scale_atmos_grid_cartesC_metric, only: &
       ROTC => ATMOS_GRID_CARTESC_METRIC_ROTC
    implicit none

    real(RP), intent(in)  :: DENS_send(KA,IA,JA)
    real(RP), intent(in)  :: MOMZ_send(KA,IA,JA)
    real(RP), intent(in)  :: MOMX_send(KA,IA,JA)
    real(RP), intent(in)  :: MOMY_send(KA,IA,JA)
    real(RP), intent(in)  :: RHOT_send(KA,IA,JA)
    real(RP), intent(in)  :: QTRC_send(KA,IA,JA,ONLINE_SEND_QA)

    integer, parameter :: HANDLE = 1

    real(RP) :: WORK1_send(KA,IA,JA)
    real(RP) :: WORK2_send(KA,IA,JA)
    real(RP) :: u_on_map, v_on_map

    real(RP) :: dummy(1,1,1)
    integer  :: tagbase, tagcomm
    integer  :: isu_tag

    integer  :: ierr
    integer  :: i, j, k, iq
    !---------------------------------------------------------------------------

    if( .NOT. USE_NESTING ) return

    tagcomm = INTERCOMM_ID(HANDLE) * order_tag_comm

    if ( COMM_CARTESC_NEST_Filiation( INTERCOMM_ID(HANDLE) ) > 0 ) then

       !##### parent [send issue] #####

       call PROF_rapstart('NEST_total_P', 2)
       call PROF_rapstart('NEST_pack_P', 2)

       nsend = nsend + 1
       LOG_INFO("COMM_CARTESC_NEST_nestdown",'(1X,A,I5,A)') "CONeP[P] send( ", nsend, " )"

       ! to keep values at that time by finish of sending process
!OCL XFILL
       org_DENS(:,:,:) = DENS_send(:,:,:)
!OCL XFILL
       org_MOMZ(:,:,:) = MOMZ_send(:,:,:)
!OCL XFILL
       org_MOMX(:,:,:) = MOMX_send(:,:,:)
!OCL XFILL
       org_MOMY(:,:,:) = MOMY_send(:,:,:)
!OCL XFILL
       org_RHOT(:,:,:) = RHOT_send(:,:,:)
       do iq = 1, ONLINE_SEND_QA
!OCL XFILL
          org_QTRC(:,:,:,iq) = QTRC_send(:,:,:,iq)
       enddo

       !*** request control
       !--- do not change the calling order below;
       !--- it should be consistent with the order in "COMM_CARTESC_NEST_recvwait_issue"
       rq_ctl_p = 0

       if ( .NOT. ONLINE_DAUGHTER_NO_ROTATE ) then
          ! from staggered point to scalar point
          do j = 1, JA
          do i = 2, IA
          do k = 1, KA
             WORK1_send(k,i,j) = ( org_MOMX(k,i-1,j) + org_MOMX(k,i,j) ) * 0.5_RP
          enddo
          enddo
          enddo

          do j = 1, JA
          do k = 1, KA
             WORK1_send(k,1,j) = org_MOMX(k,1,j)
          enddo
          enddo

          call COMM_vars8( WORK1_send(:,:,:), 1 )

          do j = 2, JA
          do i = 1, IA
          do k = 1, KA
             WORK2_send(k,i,j) = ( org_MOMY(k,i,j-1) + org_MOMY(k,i,j) ) * 0.5_RP
          enddo
          enddo
          enddo

          do i = 1, IA
          do k = 1, KA
             WORK2_send(k,i,1) = org_MOMY(k,i,1)
          enddo
          enddo

          call COMM_vars8( WORK2_send(:,:,:), 2 )

          call COMM_wait ( WORK1_send(:,:,:), 1, .false. )
          call COMM_wait ( WORK2_send(:,:,:), 2, .false. )

          ! rotation from map-projected field to latlon field
          do j = 1, JA
          do i = 1, IA
          do k = 1, KA
             u_on_map = WORK1_send(k,i,j) / org_DENS(k,i,j)
             v_on_map = WORK2_send(k,i,j) / org_DENS(k,i,j)

             org_U_ll(k,i,j) = u_on_map * ROTC(i,j,1) - v_on_map * ROTC(i,j,2)
             org_V_ll(k,i,j) = u_on_map * ROTC(i,j,2) + v_on_map * ROTC(i,j,1)
          enddo
          enddo
          enddo
       endif

       tagbase = tagcomm + tag_dens*order_tag_var
       call COMM_CARTESC_NEST_intercomm_nestdown( org_DENS(:,:,:),         & ! [IN]
                                          dummy   (:,:,:),         & ! [OUT]
                                          tagbase, I_SCLR, HANDLE, & ! [IN]
                                          isu_tag,                 & ! [INOUT]
                                          flag_dens = .true.       ) ! [IN]

       tagbase = tagcomm + tag_momz*order_tag_var
       if ( ONLINE_DAUGHTER_USE_VELZ ) then
          call COMM_CARTESC_NEST_intercomm_nestdown( org_MOMZ(:,:,:),         & ! [IN]
                                             dummy   (:,:,:),         & ! [OUT]
                                             tagbase, I_ZSTG, HANDLE, & ! [IN]
                                             isu_tag                  ) ! [INOUT]
       endif

       tagbase = tagcomm + tag_momx*order_tag_var
       if ( ONLINE_DAUGHTER_NO_ROTATE ) then
          call COMM_CARTESC_NEST_intercomm_nestdown( org_MOMX(:,:,:),         & ! [IN]
                                             dummy   (:,:,:),         & ! [OUT]
                                             tagbase, I_XSTG, HANDLE, & ! [IN]
                                             isu_tag                  ) ! [INOUT]
       else
          call COMM_CARTESC_NEST_intercomm_nestdown( org_U_ll(:,:,:),         & ! [IN]
                                             dummy   (:,:,:),         & ! [OUT]
                                             tagbase, I_SCLR, HANDLE, & ! [IN]
                                             isu_tag                  ) ! [INOUT]
       endif

       tagbase = tagcomm + tag_momy*order_tag_var
       if ( ONLINE_DAUGHTER_NO_ROTATE ) then
          call COMM_CARTESC_NEST_intercomm_nestdown( org_MOMY(:,:,:),         & ! [IN]
                                             dummy   (:,:,:),         & ! [OUT]
                                             tagbase, I_YSTG, HANDLE, & ! [IN]
                                             isu_tag                  ) ! [INOUT]
       else
          call COMM_CARTESC_NEST_intercomm_nestdown( org_V_ll(:,:,:),         & ! [IN]
                                             dummy   (:,:,:),         & ! [OUT]
                                             tagbase, I_SCLR, HANDLE, & ! [IN]
                                             isu_tag                  ) ! [INOUT]
       endif

       tagbase = tagcomm + tag_rhot*order_tag_var
       call COMM_CARTESC_NEST_intercomm_nestdown( org_RHOT(:,:,:),         & ! [IN]
                                          dummy   (:,:,:),         & ! [OUT]
                                          tagbase, I_SCLR, HANDLE, & ! [IN]
                                          isu_tag                  ) ! [INOUT]

       do iq = 1, ONLINE_SEND_QA
          tagbase = tagcomm + (tag_qx*10+iq)*order_tag_var
          call COMM_CARTESC_NEST_intercomm_nestdown( org_QTRC(:,:,:,iq),      & ! [IN]
                                             dummy   (:,:,:),         & ! [OUT]
                                             tagbase, I_SCLR, HANDLE, & ! [IN]
                                             isu_tag                  ) ! [INOUT]
       enddo

       rq_tot_p = rq_ctl_p

       call PROF_rapend  ('NEST_pack_P', 2)
       call PROF_rapend  ('NEST_total_P', 2)

    else
       LOG_ERROR("COMM_CARTESC_NEST_nestdown_send",*) 'internal error'
       call PRC_abort
    endif

    return
  end subroutine COMM_CARTESC_NEST_nestdown_send
  !-----------------------------------------------------------------------------
  !> Boundary data transfer from parent to daughter: nestdown (daughter side)
  subroutine COMM_CARTESC_NEST_nestdown_recv( &
       DENS_recv, &
       VELZ_recv, &
       VELX_recv, &
       VELY_recv, &
       POTT_recv, &
       QTRC_recv  )
    use scale_prc, only: &
       PRC_abort, &
       PRC_INTERCOMM_PARENT, &
       PRC_INTERCOMM_CHILD
    use scale_comm_cartesC, only: &
       COMM_vars8, &
       COMM_wait
    use scale_atmos_grid_cartesC_metric, only: &
       ROTC => ATMOS_GRID_CARTESC_METRIC_ROTC
    implicit none

    real(RP), intent(out) :: DENS_recv(KA,IA,JA)
    real(RP), intent(out) :: VELZ_recv(KA,IA,JA)
    real(RP), intent(out) :: VELX_recv(KA,IA,JA)
    real(RP), intent(out) :: VELY_recv(KA,IA,JA)
    real(RP), intent(out) :: POTT_recv(KA,IA,JA)
    real(RP), intent(out) :: QTRC_recv(KA,IA,JA,ONLINE_RECV_QA)

    integer, parameter :: HANDLE = 2

    real(RP) :: WORK1_recv(KA,IA,JA)
    real(RP) :: WORK2_recv(KA,IA,JA)
    real(RP) :: U_ll_recv (KA,IA,JA)
    real(RP) :: V_ll_recv (KA,IA,JA)
    real(RP) :: u_on_map, v_on_map

    real(RP) :: dummy(1,1,1)
    integer  :: tagbase, tagcomm
    integer  :: isu_tag

    integer  :: ierr
    integer  :: i, j, k, iq
    !---------------------------------------------------------------------------

    if( .NOT. USE_NESTING ) return

    tagcomm = INTERCOMM_ID(HANDLE) * order_tag_comm

    if( COMM_CARTESC_NEST_Filiation( INTERCOMM_ID(HANDLE) ) < 0 ) then

       !##### child [wait issue] #####

       call PROF_rapstart('NEST_total_C', 2)
       call PROF_rapstart('NEST_wait_C', 2)

       nwait_d = nwait_d + 1
       !LOG_INFO("COMM_CARTESC_NEST_nestdown",'(1X,A,I5,A)') "NestIDC [C]: que wait ( ", nwait_d, " )"

       !*** reset issue tag and request control
       !--- do not change the calling order below;
       !--- it should be consistent with the order in "COMM_CARTESC_NEST_recvwait_issue"
       isu_tag  = 0

       call COMM_CARTESC_NEST_waitall( rq_tot_d, ireq_d )

       if ( ONLINE_AGGRESSIVE_COMM ) then
          ! nothing to do
       else
          call MPI_BARRIER(PRC_INTERCOMM_PARENT, ierr)
       endif

       call PROF_rapend  ('NEST_wait_C', 2)
       call PROF_rapstart('NEST_unpack_C', 2)

       tagbase = tagcomm + tag_dens*order_tag_var
       call COMM_CARTESC_NEST_intercomm_nestdown( dummy     (:,:,:),       & ! [IN]
                                          WORK1_recv(:,:,:),       & ! [OUT]
                                          tagbase, I_SCLR, HANDLE, & ! [IN]
                                          isu_tag,                 & ! [INOUT]
                                          flag_dens = .true.       ) ! [IN]
!OCL XFILL
       do j = 1, JA
       do i = 1, IA
       do k = KS, KE
          DENS_recv(k,i,j) = WORK1_recv(k,i,j)
       enddo
       enddo
       enddo

       call COMM_vars8( DENS_recv, 1 )

       tagbase = tagcomm + tag_momz*order_tag_var
       if ( ONLINE_USE_VELZ ) then
          call COMM_CARTESC_NEST_intercomm_nestdown( dummy     (:,:,:),       & ! [IN]
                                             WORK2_recv(:,:,:),       & ! [OUT]
                                             tagbase, I_ZSTG, HANDLE, & ! [IN]
                                             isu_tag                  ) ! [INOUT]
!OCL XFILL
          do j = 1, JA
          do i = 1, IA
          do k = KS, KE-1
             VELZ_recv(k,i,j) = WORK2_recv(k,i,j) / ( WORK1_recv(k,i,j) + WORK1_recv(k+1,i,j) ) * 2.0_RP
          enddo
          enddo
          enddo

          do j = 1, JA
          do i = 1, IA
             VELZ_recv(KS-1,i,j) = 0.0_RP
             VELZ_recv(KE  ,i,j) = 0.0_RP
          enddo
          enddo
       endif

       call COMM_wait ( DENS_recv, 1, .false. )

       tagbase = tagcomm + tag_momx*order_tag_var
       if ( ONLINE_NO_ROTATE ) then
          ! U_ll_recv receives MOMX
          call COMM_CARTESC_NEST_intercomm_nestdown( dummy     (:,:,:),       & ! [IN]
                                             WORK1_recv(:,:,:),       & ! [OUT]
                                             tagbase, I_XSTG, HANDLE, & ! [IN]
                                             isu_tag                  ) ! [INOUT]
       else
          ! U_ll_recv receives MOMX/DENS
          call COMM_CARTESC_NEST_intercomm_nestdown( dummy    (:,:,:),        & ! [IN]
                                             U_ll_recv(:,:,:),        & ! [OUT]
                                             tagbase, I_SCLR, HANDLE, & ! [IN]
                                             isu_tag                  ) ! [INOUT]
       endif

       tagbase = tagcomm + tag_momy*order_tag_var
       if ( ONLINE_NO_ROTATE ) then
          ! V_ll_recv receives MOMY
          call COMM_CARTESC_NEST_intercomm_nestdown( dummy     (:,:,:),       & ! [IN]
                                             WORK2_recv(:,:,:),       & ! [OUT]
                                             tagbase, I_YSTG, HANDLE, & ! [IN]
                                             isu_tag                  ) ! [INOUT]
       else
          ! V_ll_recv receives MOMY/DENS
          call COMM_CARTESC_NEST_intercomm_nestdown( dummy    (:,:,:),        & ! [IN]
                                             V_ll_recv(:,:,:),        & ! [OUT]
                                             tagbase, I_SCLR, HANDLE, & ! [IN]
                                             isu_tag                  ) ! [INOUT]
       endif

       if ( ONLINE_NO_ROTATE ) then

!OCL XFILL
          do j = 1, JA
          do i = 1, IA-1
          do k = KS, KE
             VELX_recv(k,i,j) = WORK1_recv(k,i,j) / ( DENS_recv(k,i+1,j) + DENS_recv(k,i,j) ) * 2.0_RP
          enddo
          enddo
          enddo

!OCL XFILL
          do j = 1, JA
          do k = KS, KE
             VELX_recv(k,IA,j) = WORK1_recv(k,IA,j) / DENS_recv(k,IA,j)
          enddo
          enddo

          call COMM_vars8( VELX_recv, 2 )

!OCL XFILL
          do j = 1, JA-1
          do i = 1, IA
          do k = KS, KE
             VELY_recv(k,i,j) = WORK2_recv(k,i,j) / ( DENS_recv(k,i,j+1) + DENS_recv(k,i,j) ) * 2.0_RP
          enddo
          enddo
          enddo

!OCL XFILL
          do i = 1, IA
          do k = KS, KE
             VELY_recv(k,i,JA) = WORK2_recv(k,i,JA) / DENS_recv(k,i,JA)
          enddo
          enddo

          call COMM_vars8( VELY_recv, 3 )

          call COMM_wait ( VELX_recv, 2, .false. )
          call COMM_wait ( VELY_recv, 3, .false. )

       else ! rotate

          ! rotation from latlon field to map-projected field
!OCL XFILL
          do j = 1, JA
          do i = 1, IA
          do k = KS, KE
             WORK1_recv(k,i,j) =  U_ll_recv(k,i,j) * ROTC(i,j,1) + V_ll_recv(k,i,j) * ROTC(i,j,2)
             WORK2_recv(k,i,j) = -U_ll_recv(k,i,j) * ROTC(i,j,2) + V_ll_recv(k,i,j) * ROTC(i,j,1)
          enddo
          enddo
          enddo

          ! from scalar point to staggered point
!OCL XFILL
          do j = 1, JA
          do i = 1, IA-1
          do k = KS, KE
             VELX_recv(k,i,j) = ( WORK1_recv(k,i+1,j) + WORK1_recv(k,i,j) ) * 0.5_RP
          enddo
          enddo
          enddo

!OCL XFILL
          do j = 1, JA
          do k = KS, KE
             VELX_recv(k,IA,j) = WORK1_recv(k,IA,j)
          enddo
          enddo

          call COMM_vars8( VELX_recv, 2 )

!OCL XFILL
          do j = 1, JA-1
          do i = 1, IA
          do k = KS, KE
             VELY_recv(k,i,j) = ( WORK2_recv(k,i,j+1) + WORK2_recv(k,i,j) ) * 0.5_RP
          enddo
          enddo
          enddo

!OCL XFILL
          do i = 1, IA
          do k = KS, KE
             VELY_recv(k,i,JA) = WORK2_recv(k,i,JA)
          enddo
          enddo

          call COMM_vars8( VELY_recv, 3 )

          call COMM_wait ( VELX_recv, 2, .false. )
          call COMM_wait ( VELY_recv, 3, .false. )

       endif

       tagbase = tagcomm + tag_rhot*order_tag_var
       call COMM_CARTESC_NEST_intercomm_nestdown( dummy     (:,:,:),       & ! [IN]
                                          WORK1_recv(:,:,:),       & ! [OUT]
                                          tagbase, I_SCLR, HANDLE, & ! [IN]
                                          isu_tag                  ) ! [INOUT]
!OCL XFILL
       do j = 1, JA
       do i = 1, IA
       do k = KS, KE
          POTT_recv(k,i,j) = WORK1_recv(k,i,j) / DENS_recv(k,i,j)
       enddo
       enddo
       enddo

       do iq = 1, ONLINE_RECV_QA
          tagbase = tagcomm + (tag_qx*10+iq)*order_tag_var
          call COMM_CARTESC_NEST_intercomm_nestdown( dummy     (:,:,:),       & ! [IN]
                                             WORK1_recv(:,:,:),       & ! [OUT]
                                             tagbase, I_SCLR, HANDLE, & ! [IN]
                                             isu_tag                  ) ! [INOUT]
!OCL XFILL
          do j = 1, JA
          do i = 1, IA
          do k = KS, KE
             QTRC_recv(k,i,j,iq) = WORK1_recv(k,i,j)
          enddo
          enddo
          enddo
       enddo

       call PROF_rapend  ('NEST_unpack_C', 2)
       call PROF_rapend  ('NEST_total_C', 2)

    else
       LOG_ERROR("COMM_CARTESC_NEST_nestdown_recv",*) 'internal error'
       call PRC_abort
    endif

    return
  end subroutine COMM_CARTESC_NEST_nestdown_recv

  !-----------------------------------------------------------------------------
  !> Sub-command for data transfer from parent to daughter: nestdown (parent side)
  subroutine COMM_CARTESC_NEST_recvwait_issue_send
    use scale_prc, only: &
       PRC_abort, &
       PRC_INTERCOMM_CHILD
    implicit none

    integer, parameter :: HANDLE = 1

    integer :: isu_tag
    integer :: tagbase, tagcomm
    integer :: ierr
    integer :: iq
    !---------------------------------------------------------------------------

    if( .NOT. USE_NESTING ) return

    tagcomm = INTERCOMM_ID(HANDLE) * order_tag_comm

    if ( COMM_CARTESC_NEST_Filiation( INTERCOMM_ID(HANDLE) ) > 0 ) then

       !##### parent [wait issue] #####

       call PROF_rapstart('NEST_total_P', 2)
       call PROF_rapstart('NEST_wait_P', 2)

       nwait_p = nwait_p + 1
       !LOG_INFO("COMM_CARTESC_NEST_recvwait_issue",'(1X,A,I5,A)') "NestIDC [P]: que wait ( ", nwait_p, " )"

       call COMM_CARTESC_NEST_issuer_of_wait( HANDLE )

       if ( ONLINE_AGGRESSIVE_COMM ) then
          ! nothing to do
       else
          call MPI_BARRIER(PRC_INTERCOMM_CHILD, ierr)
       endif

       call PROF_rapend  ('NEST_wait_P', 2)
       call PROF_rapend  ('NEST_total_P', 2)

    else
       LOG_ERROR("COMM_CARTESC_NEST_recvwait_issue_send",*) 'internal error'
       call PRC_abort
    endif

    return
  end subroutine COMM_CARTESC_NEST_recvwait_issue_send

  !-----------------------------------------------------------------------------
  !> Sub-command for data transfer from parent to daughter: nestdown (daughter side)
  subroutine COMM_CARTESC_NEST_recvwait_issue_recv
    use scale_prc, only: &
       PRC_abort, &
       PRC_INTERCOMM_CHILD
    implicit none

    integer, parameter :: HANDLE = 2

    integer :: isu_tag
    integer :: tagbase, tagcomm
    integer :: ierr
    integer :: iq
    !---------------------------------------------------------------------------

    if( .NOT. USE_NESTING ) return

    tagcomm = INTERCOMM_ID(HANDLE) * order_tag_comm

    if( COMM_CARTESC_NEST_Filiation( INTERCOMM_ID(HANDLE) ) < 0 ) then

       !##### child [receive issue] #####

       call PROF_rapstart('NEST_total_C', 2)

       nrecv = nrecv + 1
       LOG_INFO("COMM_CARTESC_NEST_recvwait_issue_recv",'(1X,A,I5,A)') "NestIDC [C]: que recv ( ", nrecv, " )"

       !*** reset issue tag and request control
       !--- do not change the calling order below;
       !--- it should be consistent with the order in "COMM_CARTESC_NEST_nestdown"
       isu_tag  = 0
       rq_ctl_d = 0

       tagbase = tagcomm + tag_dens*order_tag_var
       call COMM_CARTESC_NEST_issuer_of_receive( tagbase, I_SCLR, HANDLE, isu_tag )

       tagbase = tagcomm + tag_momz*order_tag_var
       if ( ONLINE_USE_VELZ ) then
          call COMM_CARTESC_NEST_issuer_of_receive( tagbase, I_ZSTG, HANDLE, isu_tag )
       endif

       tagbase = tagcomm + tag_momx*order_tag_var
       if ( ONLINE_NO_ROTATE ) then
          call COMM_CARTESC_NEST_issuer_of_receive( tagbase, I_XSTG, HANDLE, isu_tag )
       else
          call COMM_CARTESC_NEST_issuer_of_receive( tagbase, I_SCLR, HANDLE, isu_tag )
       endif

       tagbase = tagcomm + tag_momy*order_tag_var
       if ( ONLINE_NO_ROTATE ) then
          call COMM_CARTESC_NEST_issuer_of_receive( tagbase, I_YSTG, HANDLE, isu_tag )
       else
          call COMM_CARTESC_NEST_issuer_of_receive( tagbase, I_SCLR, HANDLE, isu_tag )
       endif

       tagbase = tagcomm + tag_rhot*order_tag_var
       call COMM_CARTESC_NEST_issuer_of_receive( tagbase, I_SCLR, HANDLE, isu_tag )

       do iq = 1, ONLINE_RECV_QA
          tagbase = tagcomm + (tag_qx*10+iq)*order_tag_var
          call COMM_CARTESC_NEST_issuer_of_receive( tagbase, I_SCLR, HANDLE, isu_tag )
       enddo

       rq_tot_d = rq_ctl_d

       call PROF_rapend('NEST_total_C', 2)

    else
       LOG_ERROR("COMM_CARTESC_NEST_recvwait_issue_recv",*) 'internal error'
       call PRC_abort
    endif

    return
  end subroutine COMM_CARTESC_NEST_recvwait_issue_recv

  !-----------------------------------------------------------------------------
  !> Sub-command for data transfer from parent to daughter: nestdown (parent side)
  subroutine COMM_CARTESC_NEST_recv_cancel_send
    use scale_prc, only: &
       PRC_abort
    implicit none

    integer, parameter :: HANDLE = 1

    !logical :: flag
    !integer :: istatus(MPI_STATUS_SIZE)

    integer :: rq
    integer :: ierr
    !---------------------------------------------------------------------------

    if( .NOT. USE_NESTING ) return

    if ( COMM_CARTESC_NEST_Filiation( INTERCOMM_ID(HANDLE) ) > 0 ) then

       !##### parent #####
       ! Nothing to do

    else
       LOG_ERROR("COMM_CARTESC_NEST_recv_cancel_send",*) 'internal error'
       call PRC_abort
    endif

    return
  end subroutine COMM_CARTESC_NEST_recv_cancel_send

  !-----------------------------------------------------------------------------
  !> Sub-command for data transfer from parent to daughter: nestdown (daughter side)
  subroutine COMM_CARTESC_NEST_recv_cancel_recv
    use scale_prc, only: &
       PRC_abort
    implicit none

    integer, parameter :: HANDLE = 2

    !logical :: flag
    !integer :: istatus(MPI_STATUS_SIZE)

    integer :: rq
    integer :: ierr
    !---------------------------------------------------------------------------

    if( .NOT. USE_NESTING ) return

    if( COMM_CARTESC_NEST_Filiation( INTERCOMM_ID(HANDLE) ) < 0 ) then

       !##### child #####

       LOG_INFO("COMM_CARTESC_NEST_recv_cancel_recv",'(1X,A,I5,A)') "NestIDC [C]: CANCEL recv ( ", nrecv, " )"

       do rq = 1, rq_tot_d
          if ( ireq_d(rq) /= MPI_REQUEST_NULL ) then

             call MPI_CANCEL(ireq_d(rq), ierr)

!              call MPI_TEST_CANCELLED(istatus, flag, ierr)
!              if ( .NOT. flag ) then
!                 LOG_ERROR("COMM_CARTESC_NEST_recv_cancel_recv",*) 'receive actions do not cancelled, req = ', rq
!              endif
          endif
       enddo

    else
       LOG_ERROR("COMM_CARTESC_NEST_recv_cancel_recv",*) 'internal error'
       call PRC_abort
    endif

    return
  end subroutine COMM_CARTESC_NEST_recv_cancel_recv

  !-----------------------------------------------------------------------------
  !> Inter-communication from parent to daughter: nestdown
  subroutine COMM_CARTESC_NEST_intercomm_nestdown_3D( &
       pvar,     &
       dvar,     &
       tagbase,  &
       id_stag,  &
       HANDLE,   &
       isu_tag,  &
       flag_dens )
    use scale_prc, only: &
       PRC_abort, &
       PRC_INTERCOMM_PARENT, &
       PRC_INTERCOMM_CHILD
    use scale_comm_cartesC, only: &
       COMM_datatype
    use scale_interp, only: &
       INTERP_interp3d
    use scale_atmos_grid_cartesC_real, &
       REAL_CZ => ATMOS_GRID_CARTESC_REAL_CZ, &
       REAL_FZ => ATMOS_GRID_CARTESC_REAL_FZ
    implicit none

    real(RP), intent(in)    :: pvar(:,:,:) !< variable from parent domain (PARENT_KA,PARENT_IA,PARENT_JA / 1,1,1)
    real(RP), intent(out)   :: dvar(:,:,:) !< variable to daughter domain (1,1,1 / MY_KA,MY_IA,MY_JA)
    integer,  intent(in)    :: tagbase     !< communication tag of the variable
    integer,  intent(in)    :: id_stag     !< id of staggered grid option
    integer,  intent(in)    :: HANDLE      !< id number of nesting relation in this process target
    integer,  intent(inout) :: isu_tag     !< tag for receive buffer

    logical , intent(in), optional :: flag_dens !< flag of logarithmic interpolation for density

    integer :: tile_num_x

    integer :: ileng, tag, target_rank

    integer :: parent_KA

    integer :: xloc, yloc
    integer :: gxs, gxe, gys, gye ! for large  domain
    integer :: pxs, pxe, pys, pye ! for parent domain
    integer :: zs, ze

    integer :: ig, rq, yp
    logical :: no_zstag
    logical :: spline
    logical :: logarithmic

    integer :: ierr
    integer :: i, j
    !---------------------------------------------------------------------------

    if( .NOT. USE_NESTING ) return

    logarithmic = .false.
    spline = .false.
    if ( present(flag_dens) ) then
       if( flag_dens ) then
          logarithmic = .true.
          spline = .true.
       end if
    endif

    if    ( id_stag == I_SCLR ) then
       no_zstag = .true.
       ig       = I_SCLR
    elseif( id_stag == I_ZSTG ) then
       no_zstag = .false.
       ig       = I_ZSTG
    elseif( id_stag == I_XSTG ) then
       no_zstag = .true.
       ig       = I_XSTG
    elseif( id_stag == I_YSTG ) then
       no_zstag = .true.
       ig       = I_YSTG
    endif

    if ( COMM_CARTESC_NEST_Filiation( INTERCOMM_ID(HANDLE) ) > 0 ) then
       !##### parent #####

       ileng = KA * IA * JA

       rq = rq_ctl_p

       do yp = 1, NUM_YP
          rq = rq + 1

          ! send data to multiple daughter processes
          target_rank = COMM_CARTESC_NEST_TILE_LIST_YP(yp)
          tag         = tagbase + yp

          call MPI_ISEND( pvar,                &
                          ileng,               &
                          COMM_datatype,       &
                          target_rank,         &
                          tag,                 &
                          PRC_INTERCOMM_CHILD, &
                          ireq_p(rq),          &
                          ierr                 )

          dvar(:,:,:) = -1.0_RP  ! input as a dummy value
       enddo

       rq_ctl_p = rq

    elseif( COMM_CARTESC_NEST_Filiation( INTERCOMM_ID(HANDLE) ) < 0 ) then

       !##### child #####

       parent_KA = dom_info(I_PARENT)%KMAX + dom_info(I_PARENT)%KHALO * 2
       ileng = parent_KA &
             * ( dom_info(I_PARENT)%IMAX + dom_info(I_PARENT)%IHALO * 2 ) &
             * ( dom_info(I_PARENT)%JMAX + dom_info(I_PARENT)%JHALO * 2 )

       tile_num_x = dom_info(I_PARENT)%tile_num_x

       zs = 1
       ze = parent_KA

       pxs = dom_info(I_PARENT)%IHALO + 1
       pxe = dom_info(I_PARENT)%IMAX + dom_info(I_PARENT)%IHALO
       pys = dom_info(I_PARENT)%JHALO + 1
       pye = dom_info(I_PARENT)%JMAX + dom_info(I_PARENT)%JHALO

       rq = rq_ctl_d

       do yp = 1, COMM_CARTESC_NEST_TILE_ALL
          rq = rq + 1

          xloc = mod( yp-1, dom_info(I_PARENT)%TILE_NUM_X ) + 1
          yloc = int( real(yp-1) / real(dom_info(I_PARENT)%TILE_NUM_X) ) + 1

          gxs = dom_info(I_PARENT)%IMAX * (xloc-1) + 1
          gxe = dom_info(I_PARENT)%IMAX * xloc
          gys = dom_info(I_PARENT)%JMAX * (yloc-1) + 1
          gye = dom_info(I_PARENT)%JMAX * yloc

          isu_tag = isu_tag + 1

          if ( isu_tag > max_isu ) then
             LOG_ERROR("COMM_CARTESC_NEST_intercomm_nestdown_3D",*) 'Exceeded maximum issue'
             LOG_ERROR_CONT(*) 'isu_tag  = ', isu_tag
             call PRC_abort
          endif

!OCL XFILL
          buffer_ref_3D(zs:ze,gxs:gxe,gys:gye) = recvbuf_3D(zs:ze,pxs:pxe,pys:pye,isu_tag)

       enddo

       rq_ctl_d = rq

       if ( no_zstag ) then
          call INTERP_interp3d( itp_nh,                              &
                                TILEAL_KA, KHALO+1, TILEAL_KA-KHALO, &
                                TILEAL_IA, TILEAL_JA,                &
                                KA, KS, KE, IA, JA,                  &
                                igrd         (    :,:,:,ig), & ! [IN]
                                jgrd         (    :,:,:,ig), & ! [IN]
                                hfact        (    :,:,:,ig), & ! [IN]
                                kgrd         (:,:,:,:,:,ig), & ! [IN]
                                vfact        (:,  :,:,:,ig), & ! [IN]
                                buffer_ref_CZ(:,:,:),        & ! [IN]
                                REAL_CZ      (:,:,:),        & ! [IN]
                                buffer_ref_3D(:,:,:),        & ! [IN]
                                dvar         (:,:,:),        & ! [OUT]
                                spline = spline,             & ! [IN, optional]
                                logwgt = logarithmic         ) ! [IN, optional]

       else
          call INTERP_interp3d( itp_nh,                            &
                                TILEAL_KA, KHALO, TILEAL_KA-KHALO, &
                                TILEAL_IA, TILEAL_JA,              &
                                KA, KS, KE, IA, JA,                &
                                igrd         (    :,:,:,ig), & ! [IN]
                                jgrd         (    :,:,:,ig), & ! [IN]
                                hfact        (    :,:,:,ig), & ! [IN]
                                kgrd         (:,:,:,:,:,ig), & ! [IN]
                                vfact        (:,  :,:,:,ig), & ! [IN]
                                buffer_ref_FZ(:,:,:),        & ! [IN]
                                REAL_FZ      (1:,:,:),       & ! [IN]
                                buffer_ref_3D(:,:,:),        & ! [IN]
                                dvar         (:,:,:),        & ! [OUT]
                                spline = spline,             & ! [IN, optional]
                                logwgt = logarithmic         ) ! [IN, optional]
       endif

       do j = 1, JA
       do i = 1, IA
          dvar(   1:KS-1,i,j) = 0.0_RP
          dvar(KE+1:KA  ,i,j) = 0.0_RP
       enddo
       enddo

    else
       LOG_ERROR("COMM_CARTESC_NEST_intercomm_nestdown_3D",*) 'internal error'
       call PRC_abort
    endif

    return
  end subroutine COMM_CARTESC_NEST_intercomm_nestdown_3D

  !-----------------------------------------------------------------------------
  !> [substance of issuer] Inter-communication from parent to daughter: nestdown
  subroutine COMM_CARTESC_NEST_issuer_of_receive_3D( &
       tagbase, &
       id_stag, &
       HANDLE,  &
       isu_tag  )
    use scale_prc, only: &
       PRC_myrank, &
       PRC_abort,  &
       PRC_INTERCOMM_PARENT
    use scale_comm_cartesC, only: &
       COMM_datatype
    implicit none

    integer, intent(in)    :: tagbase  !< communication tag of the variable
    integer, intent(in)    :: id_stag  !< id of staggered grid option
    integer, intent(in)    :: HANDLE   !< id number of nesting relation in this process target
    integer, intent(inout) :: isu_tag  !< tag for receive buffer

    integer :: ierr, ileng
    integer :: tag, target_rank

    integer :: rq, yp
    !---------------------------------------------------------------------------

    if( .NOT. USE_NESTING ) return

    if ( COMM_CARTESC_NEST_Filiation( INTERCOMM_ID(HANDLE) ) > 0 ) then

       !##### parent #####
       ! nothing to do

    elseif( COMM_CARTESC_NEST_Filiation( INTERCOMM_ID(HANDLE) ) < 0 ) then

       !##### child #####

       ileng = ( dom_info(I_PARENT)%KMAX + dom_info(I_PARENT)%KHALO * 2 ) &
             * ( dom_info(I_PARENT)%IMAX + dom_info(I_PARENT)%IHALO * 2 ) &
             * ( dom_info(I_PARENT)%JMAX + dom_info(I_PARENT)%JHALO * 2 )

       rq = rq_ctl_d

       do yp = 1, COMM_CARTESC_NEST_TILE_ALL
          rq = rq + 1

          target_rank = COMM_CARTESC_NEST_TILE_LIST_d(yp,PRC_myrank+1)
          tag         = tagbase + call_order(yp)

          isu_tag = isu_tag + 1

          if ( isu_tag > max_isu ) then
             LOG_ERROR("COMM_CARTESC_NEST_issuer_of_receive_3D",*) 'Exceeded maximum issue'
             LOG_ERROR_CONT(*) 'isu_tag  = ', isu_tag
             call PRC_abort
          endif

          recvbuf_3D(:,:,:,isu_tag) = 0.0_RP

          call MPI_IRECV( recvbuf_3D(:,:,:,isu_tag), &
                          ileng,                     &
                          COMM_datatype,             &
                          target_rank,               &
                          tag,                       &
                          PRC_INTERCOMM_PARENT,      &
                          ireq_d(rq),                &
                          ierr                       )

       enddo

       rq_ctl_d = rq

    else
       LOG_ERROR("COMM_CARTESC_NEST_issuer_of_receive_3D",*) 'internal error'
       call PRC_abort
    endif

    return
  end subroutine COMM_CARTESC_NEST_issuer_of_receive_3D

  !-----------------------------------------------------------------------------
  !> [substance of issuer] Inter-communication from parent to daughter: nestdown
  subroutine COMM_CARTESC_NEST_issuer_of_wait_3D( &
       HANDLE )
    use scale_prc, only: &
       PRC_abort
    implicit none

    integer, intent(in) :: HANDLE  !< id number of nesting relation in this process target
    !---------------------------------------------------------------------------

    if( .NOT. USE_NESTING ) return

    if ( COMM_CARTESC_NEST_Filiation( INTERCOMM_ID(HANDLE) ) > 0 ) then

       !##### parent #####
       call COMM_CARTESC_NEST_waitall( rq_tot_p, ireq_p(:) )

    elseif( COMM_CARTESC_NEST_Filiation( INTERCOMM_ID(HANDLE) ) < 0 ) then

       !##### child #####
       ! nothing to do

    else
       LOG_ERROR("COMM_CARTESC_NEST_issuer_of_wait_3D",*) 'internal error'
       call PRC_abort
    endif

    return
  end subroutine COMM_CARTESC_NEST_issuer_of_wait_3D

  !-----------------------------------------------------------------------------
  !> [substance of comm_wait] Inter-communication
  subroutine COMM_CARTESC_NEST_waitall( &
       req_count, &
       ireq       )
    use scale_prc, only: &
       PRC_abort
    implicit none

    integer, intent(in)    :: req_count
    integer, intent(inout) :: ireq(max_rq)

    integer :: i
    integer :: ierr
    integer :: istatus(MPI_STATUS_SIZE,req_count)
    integer :: req_count2
    integer :: ireq2(max_rq)

!    logical    :: flag = .false.
!    integer(8) :: num  = 0
    !---------------------------------------------------------------------------

    if( .NOT. USE_NESTING ) return

    req_count2 = 0
    do i = 1, req_count
       if ( ireq(i) /= MPI_REQUEST_NULL ) then
          req_count2 = req_count2 + 1
          ireq2(req_count2) = ireq(i)
       endif
    enddo

    if( req_count2 /= 0 ) call MPI_WAITALL( req_count2, ireq2(1:req_count2), istatus, ierr )

!    do while ( .NOT. flag )
!       num = num + 1
!       call MPI_TESTALL( req_count, ireq, flag, istatus, ierr )
!
!       if ( num > ONLINE_WAIT_LIMIT ) then
!          LOG_ERROR("COMM_CARTESC_NEST_waitall",'(1x,A)') 'over the limit of waiting time [NESTCOM]'
!          LOG_ERROR_CONT('(1x,A)') 'over the limit of waiting time [NESTCOM]'
!          call PRC_abort
!       endif
!    enddo

    return
  end subroutine COMM_CARTESC_NEST_waitall

  !-----------------------------------------------------------------------------
  !> [check communication status] Inter-communication (parent side)
  subroutine COMM_CARTESC_NEST_test_send
    use scale_prc, only: &
       PRC_abort
    implicit none

    integer, parameter :: HANDLE = 1

    integer :: istatus(MPI_STATUS_SIZE)
    integer :: ierr
    logical :: flag
    !---------------------------------------------------------------------------

    if( .NOT. USE_NESTING ) return

    if ( COMM_CARTESC_NEST_Filiation( INTERCOMM_ID(HANDLE) ) > 0 ) then

       !##### parent #####

       call PROF_rapstart('NEST_test_P', 2)
       if ( rq_ctl_p > 0 ) call MPI_TEST(ireq_p(1), flag, istatus, ierr)
       call PROF_rapend('NEST_test_P', 2)

    else
       LOG_ERROR("COMM_CARTESC_NEST_test_send",*) 'error'
       call PRC_abort
    endif

    return
  end subroutine COMM_CARTESC_NEST_test_send

  !-----------------------------------------------------------------------------
  !> [check communication status] Inter-communication (daughter side)
  subroutine COMM_CARTESC_NEST_test_recv
    use scale_prc, only: &
       PRC_abort
    implicit none

    integer, parameter :: HANDLE = 2

    integer :: istatus(MPI_STATUS_SIZE)
    integer :: ierr
    logical :: flag
    !---------------------------------------------------------------------------

    if( .NOT. USE_NESTING ) return

    if( COMM_CARTESC_NEST_Filiation( INTERCOMM_ID(HANDLE) ) < 0 ) then

       !##### child #####

       call PROF_rapstart('NEST_test_C', 2)
       if ( rq_ctl_d > 0 ) call MPI_TEST(ireq_d(1), flag, istatus, ierr)
       call PROF_rapend('NEST_test_C', 2)

    else
       LOG_ERROR("COMM_CARTESC_NEST_test_recv",*) 'error'
       call PRC_abort
    endif

    return
  end subroutine COMM_CARTESC_NEST_test_recv

  !-----------------------------------------------------------------------------
  !> finalize
  subroutine COMM_CARTESC_NEST_finalize
    implicit none

    integer :: i

    do i = 1, num_dom
       if ( allocated( dom_info(i)%latlon_catalogue ) ) deallocate( dom_info(i)%latlon_catalogue )
       if ( allocated( dom_info(i)%tile_id ) )          deallocate( dom_info(i)%tile_id )
       num_dom = 0
    end do

    if ( allocated( COMM_CARTESC_NEST_TILE_LIST_p ) )  deallocate( COMM_CARTESC_NEST_TILE_LIST_p )
    if ( allocated( COMM_CARTESC_NEST_TILE_LIST_d ) )  deallocate( COMM_CARTESC_NEST_TILE_LIST_d )
    if ( allocated( COMM_CARTESC_NEST_TILE_LIST_YP ) ) deallocate( COMM_CARTESC_NEST_TILE_LIST_YP )

    if ( allocated( ireq_p ) ) deallocate( ireq_p )
    if ( allocated( ireq_d ) ) deallocate( ireq_d )

    if ( allocated( call_order ) ) deallocate( call_order )
    if ( allocated( recvbuf_3D ) ) deallocate( recvbuf_3D )

    if ( allocated( buffer_ref_LON ) )   deallocate( buffer_ref_LON )
    if ( allocated( buffer_ref_LONUY ) ) deallocate( buffer_ref_LONUY )
    if ( allocated( buffer_ref_LONXV ) ) deallocate( buffer_ref_LONXV )
    if ( allocated( buffer_ref_LAT ) )   deallocate( buffer_ref_LAT )
    if ( allocated( buffer_ref_LATUY ) ) deallocate( buffer_ref_LATUY )
    if ( allocated( buffer_ref_LATXV ) ) deallocate( buffer_ref_LATXV )
    if ( allocated( buffer_ref_CZ ) )    deallocate( buffer_ref_CZ )
    if ( allocated( buffer_ref_FZ ) )    deallocate( buffer_ref_FZ )
    if ( allocated( buffer_ref_3D ) )    deallocate( buffer_ref_3D )

    if ( allocated( org_DENS ) ) deallocate( org_DENS )
    if ( allocated( org_MOMZ ) ) deallocate( org_MOMZ )
    if ( allocated( org_MOMX ) ) deallocate( org_MOMX )
    if ( allocated( org_MOMY ) ) deallocate( org_MOMY )
    if ( allocated( org_U_ll ) ) deallocate( org_U_ll )
    if ( allocated( org_V_ll ) ) deallocate( org_V_ll )
    if ( allocated( org_RHOT ) ) deallocate( org_RHOT )
    if ( allocated( org_QTRC ) ) deallocate( org_QTRC )

    if ( allocated( igrd ) )  deallocate( igrd )
    if ( allocated( jgrd ) )  deallocate( jgrd )
    if ( allocated( hfact ) ) deallocate( hfact )
    if ( allocated( kgrd ) )  deallocate( kgrd )
    if ( allocated( vfact ) ) deallocate( vfact )

    return
  end subroutine COMM_CARTESC_NEST_finalize

end module scale_comm_cartesC_nest
