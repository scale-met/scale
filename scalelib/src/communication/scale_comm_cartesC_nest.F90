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
module scale_comm_cartesC_nest
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use mpi
  use scale_precision
  use scale_stdio
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
  public :: COMM_CARTESC_NEST_setup
  public :: COMM_CARTESC_NEST_domain_relate
  public :: COMM_CARTESC_NEST_domain_shape
  public :: COMM_CARTESC_NEST_nestdown
  public :: COMM_CARTESC_NEST_recvwait_issue
  public :: COMM_CARTESC_NEST_recv_cancel
  public :: COMM_CARTESC_NEST_test
  public :: COMM_CARTESC_NEST_disconnect

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  integer,  public              :: INTERCOMM_PARENT                   ! inter-communicator to parent
  integer,  public              :: INTERCOMM_DAUGHTER                 ! inter-communicator to daughter

  integer,  public              :: COMM_CARTESC_NEST_Filiation(10)    !< index of parent-daughter relation (p>0, d<0)
  integer,  public              :: HANDLING_NUM                       !< handing number of nesting relation
  integer,  public              :: COMM_CARTESC_NEST_TILE_NUM_X       !< parent tile number in x-direction
  integer,  public              :: COMM_CARTESC_NEST_TILE_NUM_Y       !< parent tile number in y-direction
  integer,  public, allocatable :: COMM_CARTESC_NEST_TILE_ID(:)       !< parent tile real id

  integer,  public              :: PARENT_KMAX(2)                     !< parent max number in z-direction
  integer,  public              :: PARENT_IMAX(2)                     !< parent max number in x-direction
  integer,  public              :: PARENT_JMAX(2)                     !< parent max number in y-direction
  integer,  public              :: PARENT_KA(2)                       !< parent max number in z-direction (with halo)
  integer,  public              :: PARENT_IA(2)                       !< parent max number in x-direction (with halo)
  integer,  public              :: PARENT_JA(2)                       !< parent max number in y-direction (with halo)
  integer,  public              :: PARENT_OKMAX(2)                    !< parent max number in oz-direction
  integer,  public              :: PARENT_LKMAX(2)                    !< parent max number in lz-direction
  real(DP), public              :: PARENT_DTSEC(2)                    !< parent DT [sec]
  integer,  public              :: PARENT_NSTEP(2)                    !< parent step [number]

  integer,  public              :: DAUGHTER_KMAX(2)                   !< daughter max number in z-direction
  integer,  public              :: DAUGHTER_IMAX(2)                   !< daughter max number in x-direction
  integer,  public              :: DAUGHTER_JMAX(2)                   !< daughter max number in y-direction
  integer,  public              :: DAUGHTER_KA(2)                     !< daughter max number in z-direction (with halo)
  integer,  public              :: DAUGHTER_IA(2)                     !< daughter max number in x-direction (with halo)
  integer,  public              :: DAUGHTER_JA(2)                     !< daughter max number in y-direction (with halo)
  integer,  public              :: DAUGHTER_OKMAX(2)                  !< daughter max number in oz-direction
  integer,  public              :: DAUGHTER_LKMAX(2)                  !< daughter max number in lz-direction
  real(DP), public              :: DAUGHTER_DTSEC(2)                  !< daughter DT [sec]
  integer,  public              :: DAUGHTER_NSTEP(2)                  !< daughter steps [number]

  integer,  public              :: PRNT_KS(2)                         !< start index in z-direction in parent
  integer,  public              :: PRNT_KE(2)                         !< end index   in z-direction in parent
  integer,  public              :: PRNT_IS(2)                         !< start index in x-direction in parent
  integer,  public              :: PRNT_IE(2)                         !< end index   in x-direction in parent
  integer,  public              :: PRNT_JS(2)                         !< start index in y-direction in parent
  integer,  public              :: PRNT_JE(2)                         !< end index   in y-direction in parent

  integer,  public              :: DATR_KS(2)                         !< start index in z-direction in daughter
  integer,  public              :: DATR_KE(2)                         !< end index   in z-direction in daughter
  integer,  public              :: DATR_IS(2)                         !< start index in x-direction in daughter
  integer,  public              :: DATR_IE(2)                         !< end index   in x-direction in daughter
  integer,  public              :: DATR_JS(2)                         !< start index in y-direction in daughter
  integer,  public              :: DATR_JE(2)                         !< end index   in y-direction in daughter

  integer,  public              :: TILEAL_KA(2)                       !< cells of all tiles in z-direction
  integer,  public              :: TILEAL_IA(2)                       !< cells of all tiles in x-direction
  integer,  public              :: TILEAL_JA(2)                       !< cells of all tiles in y-direction

  integer,  public              :: COMM_CARTESC_NEST_BND_QA              = 1       !< number of tracer treated in nesting system
  integer,  public              :: COMM_CARTESC_NEST_INTERP_LEVEL        = 5       !< horizontal interpolation level
  integer,  public              :: COMM_CARTESC_NEST_INTERP_WEIGHT_ORDER = 2       !< horizontal interpolation weight order

  logical,  public              :: USE_NESTING              = .false.
  logical,  public              :: OFFLINE                  = .false.
  logical,  public              :: ONLINE_IAM_PARENT        = .false. !< a flag to say "I am a parent"
  logical,  public              :: ONLINE_IAM_DAUGHTER      = .false. !< a flag to say "I am a daughter"
  integer,  public              :: ONLINE_DOMAIN_NUM        = 1
  logical,  public              :: ONLINE_USE_VELZ          = .false.
  logical,  public              :: ONLINE_NO_ROTATE         = .false.
  logical,  public              :: ONLINE_BOUNDARY_USE_QHYD = .false.
  logical,  public              :: ONLINE_BOUNDARY_DIAGQNUM = .false.

  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  private :: COMM_CARTESC_NEST_parentsize
  private :: COMM_CARTESC_NEST_catalogue
  private :: COMM_CARTESC_NEST_ping
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
  real(RP), private, allocatable :: latlon_catalog(:,:,:)    !< parent latlon catalog [rad]
  real(RP), private              :: corner_loc(4,2)          !< local corner location [rad]

  integer,  private              :: PARENT_PRC_NUM_X(2)      !< MPI processes in x-direction in parent
  integer,  private              :: PARENT_PRC_NUM_Y(2)      !< MPI processes in y-direction in parent
  integer,  private              :: PARENT_PRC_nprocs(2)     !< MPI total processes in parent

  integer,  private              :: DAUGHTER_PRC_NUM_X(2)    !< MPI processes in x-direction in daughter
  integer,  private              :: DAUGHTER_PRC_NUM_Y(2)    !< MPI processes in y-direction in daughter
  integer,  private              :: DAUGHTER_PRC_nprocs(2)   !< MPI total processes in daughter

  integer,  private              :: COMM_CARTESC_NEST_TILE_ALL            !< NUM of TILEs in the local node
  integer,  private              :: COMM_CARTESC_NEST_TILE_ALLMAX_p       !< MAXNUM of TILEs among whole processes for parent
  integer,  private              :: COMM_CARTESC_NEST_TILE_ALLMAX_d       !< MAXNUM of TILEs among whole processes for daughter
  integer,  private, allocatable :: COMM_CARTESC_NEST_TILE_LIST_p(:,:)    !< relationship list in whole system for parent
  integer,  private, allocatable :: COMM_CARTESC_NEST_TILE_LIST_d(:,:)    !< relationship list in whole system for daughter
  integer,  private, allocatable :: COMM_CARTESC_NEST_TILE_LIST_YP(:)     !< yellow-page of daughter targets for parent
  integer,  private              :: NUM_YP                   !< page number of yellow-page

  character(len=H_LONG), private :: OFFLINE_PARENT_BASENAME   !< parent file base name
  integer,               private :: OFFLINE_PARENT_PRC_NUM_X  !< MPI processes in x-direction in parent [for namelist]
  integer,               private :: OFFLINE_PARENT_PRC_NUM_Y  !< MPI processes in y-direction in parent [for namelist]
  integer,               private :: OFFLINE_PARENT_KMAX       !< parent max number in z-direction [for namelist]
  integer,               private :: OFFLINE_PARENT_IMAX       !< parent max number in x-direction [for namelist]
  integer,               private :: OFFLINE_PARENT_JMAX       !< parent max number in y-direction [for namelist]
  integer,               private :: OFFLINE_PARENT_LKMAX      !< parent max number in lz-direction [for namelist]
  integer,               private :: OFFLINE_PARENT_OKMAX      !< parent max number in oz-direction [for namelist]
  integer(8),            private :: ONLINE_WAIT_LIMIT         !< limit times of waiting loop in "COMM_CARTESC_NEST_waitall"
  logical,               private :: ONLINE_DAUGHTER_USE_VELZ
  logical,               private :: ONLINE_DAUGHTER_NO_ROTATE
  logical,               private :: ONLINE_AGGRESSIVE_COMM

  integer,  parameter :: I_LON    = 1
  integer,  parameter :: I_LAT    = 2

  integer,  parameter :: I_NW     = 1
  integer,  parameter :: I_NE     = 2
  integer,  parameter :: I_SW     = 3
  integer,  parameter :: I_SE     = 4
  integer,  parameter :: I_BNDQA  = 20                      !< tentative approach (prefixed allocate size)

  integer,  parameter :: I_SCLR   = 1                       !< interpolation kinds of grid point (scalar)
  integer,  parameter :: I_ZSTG   = 2                       !< interpolation kinds of grid point (z-axis staggered)
  integer,  parameter :: I_XSTG   = 3                       !< interpolation kinds of grid point (x-axis staggered)
  integer,  parameter :: I_YSTG   = 4                       !< interpolation kinds of grid point (y-axis staggered)

  integer,  parameter :: itp_ng   = 4                       !< # of interpolation kinds of grid point
  integer,  private   :: itp_nh   = 4                       !< # of interpolation kinds of horizontal direction
  integer,  private   :: itp_nv   = 2                       !< # of interpolation kinds of vertical direction

  integer,  parameter :: tag_lon  = 1
  integer,  parameter :: tag_lat  = 2
  integer,  parameter :: tag_lonu = 3
  integer,  parameter :: tag_latu = 4
  integer,  parameter :: tag_lonv = 5
  integer,  parameter :: tag_latv = 6
  integer,  parameter :: tag_cz   = 7
  integer,  parameter :: tag_fz   = 8

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

  integer,  private, parameter   :: max_isu   = 100        ! maximum number of receive/wait issue
  integer,  private, parameter   :: max_isuf  = 20         ! maximum number of receive/wait issue (z-stag)
  integer,  private, parameter   :: max_bndqa = 12         ! maximum number of QA in boundary: tentative approach
  integer,  private              :: max_rq    = 1000       ! maximum number of req: tentative approach
  integer,  private              :: rq_ctl_p               ! for control request id (counting)
  integer,  private              :: rq_ctl_d               ! for control request id (counting)
  integer,  private              :: rq_tot_p               ! for control request id (total number)
  integer,  private              :: rq_tot_d               ! for control request id (total number)
  integer,  private, allocatable :: ireq_p(:)              ! buffer of request-id for parent
  integer,  private, allocatable :: ireq_d(:)              ! buffer of request-id for daughter
  integer,  private, allocatable :: call_order(:)          ! calling order from parent

  real(RP), private, allocatable :: buffer_2D  (:,:)       ! buffer of communicator: 2D (with HALO)
  real(RP), private, allocatable :: buffer_3D  (:,:,:)     ! buffer of communicator: 3D (with HALO)
  real(RP), private, allocatable :: buffer_3DF (:,:,:)     ! buffer of communicator: 3D-Kface (with HALO)
  real(RP), private, allocatable :: recvbuf_3D (:,:,:,:)   ! buffer of receiver: 3D (with HALO)
  real(RP), private, allocatable :: recvbuf_3DF(:,:,:,:)   ! buffer of receiver: 3D-Kface (with HALO)

  real(RP), private, allocatable :: buffer_ref_LON (:,:)   ! buffer of communicator: LON
  real(RP), private, allocatable :: buffer_ref_LONU(:,:)   ! buffer of communicator: LONU
  real(RP), private, allocatable :: buffer_ref_LONV(:,:)   ! buffer of communicator: LONV
  real(RP), private, allocatable :: buffer_ref_LAT (:,:)   ! buffer of communicator: LAT
  real(RP), private, allocatable :: buffer_ref_LATU(:,:)   ! buffer of communicator: LATU
  real(RP), private, allocatable :: buffer_ref_LATV(:,:)   ! buffer of communicator: LATV
  real(RP), private, allocatable :: buffer_ref_CZ  (:,:,:) ! buffer of communicator: CZ
  real(RP), private, allocatable :: buffer_ref_FZ  (:,:,:) ! buffer of communicator: FZ


  real(RP), private, allocatable :: buffer_ref_3D (:,:,:)  ! buffer of communicator: 3D data      (with HALO)
  real(RP), private, allocatable :: buffer_ref_3DF(:,:,:)  ! buffer of communicator: 3D at z-Face (with HALO)

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
  real(RP), private, allocatable :: vfact(:,:,:,:,:,:)     ! interpolation factor for vertical direction

  integer(8), private :: nwait_p, nwait_d, nrecv, nsend

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine COMM_CARTESC_NEST_setup ( &
       inter_parent, &
       inter_child   )
    use scale_file, only: &
       FILE_open, &
       FILE_get_attribute, &
       FILE_get_shape
    use scale_const, only: &
       D2R => CONST_D2R
    use scale_process, only: &
       PRC_MPIstop,         &
       PRC_GLOBAL_domainID, &
       PRC_IsMaster
    use scale_rm_process, only: &
       PRC_HAS_W, &
       PRC_HAS_E, &
       PRC_HAS_S, &
       PRC_HAS_N
    use scale_interp, only: &
       INTRP_setup,   &
       INTRP_factor3d
    use scale_comm, only: &
       COMM_Bcast
    use scale_atmos_grid_cartesC_real, only: &
       ATMOS_GRID_CARTESC_REAL_LON,   &
       ATMOS_GRID_CARTESC_REAL_LAT,   &
       ATMOS_GRID_CARTESC_REAL_LONU,  &
       ATMOS_GRID_CARTESC_REAL_LONV,  &
       ATMOS_GRID_CARTESC_REAL_LONUV, &
       ATMOS_GRID_CARTESC_REAL_LATU,  &
       ATMOS_GRID_CARTESC_REAL_LATV,  &
       ATMOS_GRID_CARTESC_REAL_LATUV, &
       ATMOS_GRID_CARTESC_REAL_CZ,    &
       ATMOS_GRID_CARTESC_REAL_FZ
    use scale_atmos_hydrometeor, only: &
       I_QV
    use scale_atmos_phy_mp, only: &
       QA_MP
    implicit none

    integer, intent(in), optional :: inter_parent
    integer, intent(in), optional :: inter_child

    !< metadata files for lat-lon domain for all processes
    character(len=H_LONG)  :: LATLON_CATALOGUE_FNAME = 'latlon_domain_catalogue.txt'

    integer :: ONLINE_SPECIFIED_MAXRQ = 0
    integer :: i
    integer :: fid, ierr
    integer :: parent_id

    logical :: flag_parent = .false.
    logical :: flag_child  = .false.

    integer :: imaxg(1), jmaxg(1)
    integer :: pnum_x(1), pnum_y(1)
    integer :: dims(1)
    logical :: error, existed

    namelist / PARAM_COMM_CARTESC_NEST /      &
       LATLON_CATALOGUE_FNAME,   &
       OFFLINE_PARENT_BASENAME,  &
       OFFLINE_PARENT_PRC_NUM_X, &
       OFFLINE_PARENT_PRC_NUM_Y, &
       ONLINE_DOMAIN_NUM,        &
       ONLINE_IAM_PARENT,        &
       ONLINE_IAM_DAUGHTER,      &
       ONLINE_USE_VELZ,          &
       ONLINE_NO_ROTATE,         &
       ONLINE_BOUNDARY_USE_QHYD, &
       ONLINE_BOUNDARY_DIAGQNUM, &
       ONLINE_AGGRESSIVE_COMM,   &
       ONLINE_WAIT_LIMIT,        &
       ONLINE_SPECIFIED_MAXRQ,   &
       COMM_CARTESC_NEST_INTERP_LEVEL,        &
       COMM_CARTESC_NEST_INTERP_WEIGHT_ORDER

    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '++++++ Module[GRID_NEST] / Categ[ATMOS-RM GRID] / Origin[SCALElib]'

    if( inter_parent /= MPI_COMM_NULL ) flag_child  = .true. ! exist parent, so work as a child
    if( inter_child  /= MPI_COMM_NULL ) flag_parent = .true. ! exist child, so work as a parent

    OFFLINE_PARENT_BASENAME = ""

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
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_COMM_CARTESC_NEST. Check!'
       call PRC_MPIstop
    endif
    if( IO_NML ) write(IO_FID_NML,nml=PARAM_COMM_CARTESC_NEST)

    PRC_GLOBAL_domainID = ONLINE_DOMAIN_NUM

    if ( OFFLINE_PARENT_BASENAME /= "" ) then

       OFFLINE = .true.
       USE_NESTING = .true.

       if ( PRC_IsMaster ) then
          call FILE_open( OFFLINE_PARENT_BASENAME, & ! (in)
                          fid,                     & ! (out)
                          aggregate = .false.      ) ! (in)

          call FILE_get_attribute( fid, "global", "scale_atmos_grid_cartesC_imaxg", &
                                   imaxg(:), existed=existed                        )
          if ( existed ) then
             call FILE_get_attribute( fid, "global", "scale_cartesC_prc_num_x", &
                                      pnum_x(:), existed=existed                )
          end if
          if ( existed ) then
             OFFLINE_PARENT_IMAX = imaxg(1) / pnum_x(1)
          else
             ! for old file
             call FILE_get_shape( fid, "CX", dims(:) )
             OFFLINE_PARENT_IMAX = dims(1)-IHALO*2
          end if

          call FILE_get_attribute( fid, "global", "scale_atmos_grid_cartesC_jmaxg", &
                                   jmaxg(:), existed=existed                        )
          if ( existed ) then
             call FILE_get_attribute( fid, "global", "scale_cartesC_prc_num_y", &
                                      pnum_y(:), existed=existed                )
          end if
          if ( existed ) then
             OFFLINE_PARENT_JMAX = jmaxg(1) / pnum_y(1)
          else
             ! for old file
             call FILE_get_shape( fid, "CY", dims(:) )
             OFFLINE_PARENT_JMAX = dims(1)-JHALO*2
          end if

          call FILE_get_attribute( fid, "global", "scale_atmos_grid_cartesC_kmax", &
                                   dims(:), existed=existed                        )
          if ( existed ) then
             OFFLINE_PARENT_KMAX = dims(1)
          else
             call FILE_get_shape( fid, "z", dims(:), error=error )
             if ( error ) then
                OFFLINE_PARENT_KMAX = 0
             else
                OFFLINE_PARENT_KMAX = dims(1)
             endif
          end if

          call FILE_get_attribute( fid, "global", "scale_ocean_grid_cartesC_kmax", &
                                   dims(:), existed=existed                        )
          if ( existed ) then
             OFFLINE_PARENT_OKMAX = dims(1)
          else
             ! for old file
             call FILE_get_shape( fid, "oz", dims(:), error=error )
             if ( error ) then
                OFFLINE_PARENT_OKMAX = 0
             else
                OFFLINE_PARENT_OKMAX = dims(1)
             endif
          end if

          call FILE_get_attribute( fid, "global", "scale_land_grid_cartesC_kmax", &
                                   dims(:), existed=existed                       )
          if ( existed ) then
             OFFLINE_PARENT_LKMAX = dims(1)
          else
             ! for old file
             call FILE_get_shape( fid, "lz", dims(:), error=error )
             if ( error ) then
                OFFLINE_PARENT_LKMAX = 0
             else
                OFFLINE_PARENT_LKMAX = dims(1)
             endif
          end if

       endif
       call COMM_Bcast( OFFLINE_PARENT_IMAX  )
       call COMM_Bcast( OFFLINE_PARENT_JMAX  )
       call COMM_Bcast( OFFLINE_PARENT_KMAX  )
       call COMM_Bcast( OFFLINE_PARENT_OKMAX )
       call COMM_Bcast( OFFLINE_PARENT_LKMAX )
    endif

    if ( ONLINE_IAM_DAUGHTER .or. ONLINE_IAM_PARENT ) then

       if ( OFFLINE ) then
          write(*,*) 'xxx OFFLINE and ONLINE cannot be use at the same time'
          call PRC_MPIstop
       endif

       USE_NESTING = .true.
    endif

    call INTRP_setup( interp_search_divnum,    & ! [IN]
                      COMM_CARTESC_NEST_INTERP_WEIGHT_ORDER ) ! [IN]

    itp_nh = COMM_CARTESC_NEST_INTERP_LEVEL
    itp_nv = 2

    DEBUG_DOMAIN_NUM = ONLINE_DOMAIN_NUM
    if( ONLINE_SPECIFIED_MAXRQ > max_rq ) max_rq = ONLINE_SPECIFIED_MAXRQ

    allocate( ireq_p(max_rq)     )
    allocate( ireq_d(max_rq)     )
    allocate( call_order(max_rq) )
    ireq_p(:) = MPI_REQUEST_NULL
    ireq_d(:) = MPI_REQUEST_NULL

    if ( USE_NESTING ) then

       if ( OFFLINE .OR. ONLINE_IAM_DAUGHTER ) then
          corner_loc(I_NW,I_LON) = ATMOS_GRID_CARTESC_REAL_LONUV( 0,JA) / D2R
          corner_loc(I_NE,I_LON) = ATMOS_GRID_CARTESC_REAL_LONUV(IA,JA) / D2R
          corner_loc(I_SW,I_LON) = ATMOS_GRID_CARTESC_REAL_LONUV( 0, 0) / D2R
          corner_loc(I_SE,I_LON) = ATMOS_GRID_CARTESC_REAL_LONUV(IA, 0) / D2R
          corner_loc(I_NW,I_LAT) = ATMOS_GRID_CARTESC_REAL_LATUV( 0,JA) / D2R
          corner_loc(I_NE,I_LAT) = ATMOS_GRID_CARTESC_REAL_LATUV(IA,JA) / D2R
          corner_loc(I_SW,I_LAT) = ATMOS_GRID_CARTESC_REAL_LATUV( 0, 0) / D2R
          corner_loc(I_SE,I_LAT) = ATMOS_GRID_CARTESC_REAL_LATUV(IA, 0) / D2R
       endif

       if ( OFFLINE ) then

         HANDLING_NUM = 1
         PARENT_PRC_NUM_X(HANDLING_NUM) = OFFLINE_PARENT_PRC_NUM_X
         PARENT_PRC_NUM_Y(HANDLING_NUM) = OFFLINE_PARENT_PRC_NUM_Y
         PARENT_KMAX(HANDLING_NUM)      = OFFLINE_PARENT_KMAX
         PARENT_IMAX(HANDLING_NUM)      = OFFLINE_PARENT_IMAX
         PARENT_JMAX(HANDLING_NUM)      = OFFLINE_PARENT_JMAX
         PARENT_OKMAX(HANDLING_NUM)     = OFFLINE_PARENT_OKMAX
         PARENT_LKMAX(HANDLING_NUM)     = OFFLINE_PARENT_LKMAX

         PARENT_PRC_nprocs(HANDLING_NUM) = PARENT_PRC_NUM_X(HANDLING_NUM) * PARENT_PRC_NUM_Y(HANDLING_NUM)
         allocate( latlon_catalog(PARENT_PRC_nprocs(HANDLING_NUM),4,2) )

         !--- read latlon catalogue
         fid = IO_get_available_fid()
         open( fid,                                    &
               file   = trim(LATLON_CATALOGUE_FNAME),  &
               form   = 'formatted',                   &
               status = 'old',                         &
               iostat = ierr                           )

         if ( ierr /= 0 ) then
            write(*,*) 'xxx [NEST_setup] cannot open latlon-catalogue file!'
            call PRC_MPIstop
         endif

         do i = 1, PARENT_PRC_nprocs(HANDLING_NUM)
            read(fid,'(i8,8f32.24)',iostat=ierr) parent_id, &
                                                 latlon_catalog(i,I_NW,I_LON), latlon_catalog(i,I_NE,I_LON), & ! LON: NW, NE
                                                 latlon_catalog(i,I_SW,I_LON), latlon_catalog(i,I_SE,I_LON), & ! LON: SW, SE
                                                 latlon_catalog(i,I_NW,I_LAT), latlon_catalog(i,I_NE,I_LAT), & ! LAT: NW, NE
                                                 latlon_catalog(i,I_SW,I_LAT), latlon_catalog(i,I_SE,I_LAT)    ! LAT: SW, SE
            if ( i /= parent_id ) then
               write(*,*) 'xxx [NEST_setup] internal error: parent mpi id'
               call PRC_MPIstop
            endif
            if ( ierr /= 0 ) exit
         enddo
         close(fid)

         call COMM_CARTESC_NEST_domain_relate(HANDLING_NUM)

      else ! ONLINE RELATIONSHIP
!         if ( present(flag_parent) .AND. present(flag_child) ) then
!            if( IO_L ) write(IO_FID_LOG,'(1x,A)') &
!                       '*** Setup Online Nesting Inter-Domain Communicator (IDC)'
!         else
!            write(*,*) 'xxx Internal Error:'
!            write(*,*) 'xxx The flag_parent and flag_child are needed.'
!            write(*,*) '    domain: ', ONLINE_DOMAIN_NUM
!            call PRC_MPIstop
!         endif

         if( ONLINE_BOUNDARY_USE_QHYD ) then
            COMM_CARTESC_NEST_BND_QA = QA_MP
         elseif ( I_QV > 0 ) then
            COMM_CARTESC_NEST_BND_QA = 1
         else
            COMM_CARTESC_NEST_BND_QA = 0
         endif

         if( IO_L ) write(IO_FID_LOG,*) "flag_parent", flag_parent, "flag_child", flag_child
         if( IO_L ) write(IO_FID_LOG,*) "ONLINE_IAM_PARENT", ONLINE_IAM_PARENT, "ONLINE_IAM_DAUGHTER", ONLINE_IAM_DAUGHTER

         if( flag_parent ) then ! must do first before daughter processes
         !-------------------------------------------------
            if ( .NOT. ONLINE_IAM_PARENT ) then
               write(*,*) 'xxx [NEST_setup] Parent Flag from launcher is not consistent with namelist!'
               write(*,*) 'xxx PARENT - domain : ', ONLINE_DOMAIN_NUM
               call PRC_MPIstop
            endif

            HANDLING_NUM = 1 !HANDLING_NUM + 1
            INTERCOMM_ID(HANDLING_NUM) = ONLINE_DOMAIN_NUM
            COMM_CARTESC_NEST_Filiation(INTERCOMM_ID(HANDLING_NUM)) = 1

            INTERCOMM_DAUGHTER = inter_child
            if( IO_L ) write(IO_FID_LOG,'(1x,A,I2,A)') '*** Online Nesting - PARENT [INTERCOMM_ID:', &
                                                        INTERCOMM_ID(HANDLING_NUM), ' ]'
            if( IO_L ) write(IO_FID_LOG,*) '*** Online Nesting - INTERCOMM :', INTERCOMM_DAUGHTER

            call COMM_CARTESC_NEST_ping( HANDLING_NUM )

            call COMM_CARTESC_NEST_parentsize( HANDLING_NUM )

            call COMM_CARTESC_NEST_catalogue( HANDLING_NUM )
            call MPI_BARRIER(INTERCOMM_DAUGHTER, ierr)

            PARENT_KA(HANDLING_NUM)   = PARENT_KMAX(HANDLING_NUM)   + KHALO * 2
            PARENT_IA(HANDLING_NUM)   = PARENT_IMAX(HANDLING_NUM)   + IHALO * 2
            PARENT_JA(HANDLING_NUM)   = PARENT_JMAX(HANDLING_NUM)   + JHALO * 2
            DAUGHTER_KA(HANDLING_NUM) = DAUGHTER_KMAX(HANDLING_NUM) + KHALO * 2
            DAUGHTER_IA(HANDLING_NUM) = DAUGHTER_IMAX(HANDLING_NUM) + IHALO * 2
            DAUGHTER_JA(HANDLING_NUM) = DAUGHTER_JMAX(HANDLING_NUM) + JHALO * 2
            TILEAL_KA(HANDLING_NUM)   = 0
            TILEAL_IA(HANDLING_NUM)   = 0
            TILEAL_JA(HANDLING_NUM)   = 0

            if( IO_L ) write(IO_FID_LOG,'(1x,A)'     ) '***  Informations of Parent Domain [me]'
            if( IO_L ) write(IO_FID_LOG,'(1x,A,I6)'  ) '***  --- PARENT_PRC_nprocs   :', PARENT_PRC_nprocs(HANDLING_NUM)
            if( IO_L ) write(IO_FID_LOG,'(1x,A,I6)'  ) '***  --- PARENT_PRC_NUM_X    :', PARENT_PRC_NUM_X(HANDLING_NUM)
            if( IO_L ) write(IO_FID_LOG,'(1x,A,I6)'  ) '***  --- PARENT_PRC_NUM_Y    :', PARENT_PRC_NUM_Y(HANDLING_NUM)
            if( IO_L ) write(IO_FID_LOG,'(1x,A,I6)'  ) '***  --- PARENT_KMAX         :', PARENT_KMAX(HANDLING_NUM)
            if( IO_L ) write(IO_FID_LOG,'(1x,A,I6)'  ) '***  --- PARENT_IMAX         :', PARENT_IMAX(HANDLING_NUM)
            if( IO_L ) write(IO_FID_LOG,'(1x,A,I6)'  ) '***  --- PARENT_JMAX         :', PARENT_JMAX(HANDLING_NUM)
            if( IO_L ) write(IO_FID_LOG,'(1x,A,F9.3)') '***  --- PARENT_DTSEC        :', PARENT_DTSEC(HANDLING_NUM)
            if( IO_L ) write(IO_FID_LOG,'(1x,A,I6)  ') '***  --- PARENT_NSTEP        :', PARENT_NSTEP(HANDLING_NUM)
            if( IO_L ) write(IO_FID_LOG,'(1x,A)'     ) '***  Informations of Daughter Domain'
            if( IO_L ) write(IO_FID_LOG,'(1x,A,I6)'  ) '***  --- DAUGHTER_PRC_nprocs :', DAUGHTER_PRC_nprocs(HANDLING_NUM)
            if( IO_L ) write(IO_FID_LOG,'(1x,A,I6)'  ) '***  --- DAUGHTER_PRC_NUM_X  :', DAUGHTER_PRC_NUM_X(HANDLING_NUM)
            if( IO_L ) write(IO_FID_LOG,'(1x,A,I6)'  ) '***  --- DAUGHTER_PRC_NUM_Y  :', DAUGHTER_PRC_NUM_Y(HANDLING_NUM)
            if( IO_L ) write(IO_FID_LOG,'(1x,A,I6)'  ) '***  --- DAUGHTER_KMAX       :', DAUGHTER_KMAX(HANDLING_NUM)
            if( IO_L ) write(IO_FID_LOG,'(1x,A,I6)'  ) '***  --- DAUGHTER_IMAX       :', DAUGHTER_IMAX(HANDLING_NUM)
            if( IO_L ) write(IO_FID_LOG,'(1x,A,I6)'  ) '***  --- DAUGHTER_JMAX       :', DAUGHTER_JMAX(HANDLING_NUM)
            if( IO_L ) write(IO_FID_LOG,'(1x,A,F9.3)') '***  --- DAUGHTER_DTSEC      :', DAUGHTER_DTSEC(HANDLING_NUM)
            if( IO_L ) write(IO_FID_LOG,'(1x,A,I6)  ') '***  --- DAUGHTER_NSTEP      :', DAUGHTER_NSTEP(HANDLING_NUM)
            if( IO_L ) write(IO_FID_LOG,'(1x,A,I6)  ') '***  Limit Num. NCOMM req.   :', max_rq

            allocate( org_DENS(PARENT_KA(HANDLING_NUM),PARENT_IA(HANDLING_NUM),PARENT_JA(HANDLING_NUM))           )
            allocate( org_MOMZ(PARENT_KA(HANDLING_NUM),PARENT_IA(HANDLING_NUM),PARENT_JA(HANDLING_NUM))           )
            allocate( org_MOMX(PARENT_KA(HANDLING_NUM),PARENT_IA(HANDLING_NUM),PARENT_JA(HANDLING_NUM))           )
            allocate( org_MOMY(PARENT_KA(HANDLING_NUM),PARENT_IA(HANDLING_NUM),PARENT_JA(HANDLING_NUM))           )
            allocate( org_U_ll(PARENT_KA(HANDLING_NUM),PARENT_IA(HANDLING_NUM),PARENT_JA(HANDLING_NUM))           )
            allocate( org_V_ll(PARENT_KA(HANDLING_NUM),PARENT_IA(HANDLING_NUM),PARENT_JA(HANDLING_NUM))           )
            allocate( org_RHOT(PARENT_KA(HANDLING_NUM),PARENT_IA(HANDLING_NUM),PARENT_JA(HANDLING_NUM))           )
            allocate( org_QTRC(PARENT_KA(HANDLING_NUM),PARENT_IA(HANDLING_NUM),PARENT_JA(HANDLING_NUM),max_bndqa) )

            call COMM_CARTESC_NEST_setup_nestdown( HANDLING_NUM )

         !---------------------------------- end of parent routines
         endif


         if( flag_child ) then
         !-------------------------------------------------
            if ( .NOT. ONLINE_IAM_DAUGHTER ) then
               write(*,*) 'xxx [NEST_setup] Child Flag from launcher is not consistent with namelist!'
               write(*,*) 'xxx DAUGHTER - domain : ', ONLINE_DOMAIN_NUM
               call PRC_MPIstop
            endif

            HANDLING_NUM = 2 !HANDLING_NUM + 1
            INTERCOMM_ID(HANDLING_NUM) = ONLINE_DOMAIN_NUM - 1
            COMM_CARTESC_NEST_Filiation(INTERCOMM_ID(HANDLING_NUM)) = -1

            INTERCOMM_PARENT = inter_parent
            if( IO_L ) write(IO_FID_LOG,'(1x,A,I2,A)') '*** Online Nesting - DAUGHTER [INTERCOMM_ID:', &
                                                        INTERCOMM_ID(HANDLING_NUM), ' ]'
            if( IO_L ) write(IO_FID_LOG,*) '*** Online Nesting - INTERCOMM :', INTERCOMM_PARENT

            call COMM_CARTESC_NEST_ping( HANDLING_NUM )

            call COMM_CARTESC_NEST_parentsize( HANDLING_NUM )

            allocate( latlon_catalog(PARENT_PRC_nprocs(HANDLING_NUM),4,2) )
            call COMM_CARTESC_NEST_catalogue( HANDLING_NUM )
            call MPI_BARRIER(INTERCOMM_PARENT, ierr)

            call COMM_CARTESC_NEST_domain_relate( HANDLING_NUM )

            PARENT_KA  (HANDLING_NUM) = PARENT_KMAX  (HANDLING_NUM) + KHALO * 2
            PARENT_IA  (HANDLING_NUM) = PARENT_IMAX  (HANDLING_NUM) + IHALO * 2
            PARENT_JA  (HANDLING_NUM) = PARENT_JMAX  (HANDLING_NUM) + JHALO * 2
            DAUGHTER_KA(HANDLING_NUM) = DAUGHTER_KMAX(HANDLING_NUM) + KHALO * 2
            DAUGHTER_IA(HANDLING_NUM) = DAUGHTER_IMAX(HANDLING_NUM) + IHALO * 2
            DAUGHTER_JA(HANDLING_NUM) = DAUGHTER_JMAX(HANDLING_NUM) + JHALO * 2
            TILEAL_KA  (HANDLING_NUM) = PARENT_KA    (HANDLING_NUM)
            TILEAL_IA  (HANDLING_NUM) = PARENT_IMAX  (HANDLING_NUM) * COMM_CARTESC_NEST_TILE_NUM_X
            TILEAL_JA  (HANDLING_NUM) = PARENT_JMAX  (HANDLING_NUM) * COMM_CARTESC_NEST_TILE_NUM_Y

            if( IO_L ) write(IO_FID_LOG,'(1x,A)'     ) '***  Informations of Parent Domain'
            if( IO_L ) write(IO_FID_LOG,'(1x,A,I6)'  ) '***  --- PARENT_PRC_nprocs   :', PARENT_PRC_nprocs(HANDLING_NUM)
            if( IO_L ) write(IO_FID_LOG,'(1x,A,I6)'  ) '***  --- PARENT_PRC_NUM_X    :', PARENT_PRC_NUM_X(HANDLING_NUM)
            if( IO_L ) write(IO_FID_LOG,'(1x,A,I6)'  ) '***  --- PARENT_PRC_NUM_Y    :', PARENT_PRC_NUM_Y(HANDLING_NUM)
            if( IO_L ) write(IO_FID_LOG,'(1x,A,I6)'  ) '***  --- PARENT_KMAX         :', PARENT_KMAX(HANDLING_NUM)
            if( IO_L ) write(IO_FID_LOG,'(1x,A,I6)'  ) '***  --- PARENT_IMAX         :', PARENT_IMAX(HANDLING_NUM)
            if( IO_L ) write(IO_FID_LOG,'(1x,A,I6)'  ) '***  --- PARENT_JMAX         :', PARENT_JMAX(HANDLING_NUM)
            if( IO_L ) write(IO_FID_LOG,'(1x,A,F9.3)') '***  --- PARENT_DTSEC        :', PARENT_DTSEC(HANDLING_NUM)
            if( IO_L ) write(IO_FID_LOG,'(1x,A,I6)'  ) '***  --- PARENT_NSTEP        :', PARENT_NSTEP(HANDLING_NUM)
            if( IO_L ) write(IO_FID_LOG,'(1x,A)'     ) '***  Informations of Daughter Domain [me]'
            if( IO_L ) write(IO_FID_LOG,'(1x,A,I6)'  ) '***  --- DAUGHTER_PRC_nprocs :', DAUGHTER_PRC_nprocs(HANDLING_NUM)
            if( IO_L ) write(IO_FID_LOG,'(1x,A,I6)'  ) '***  --- DAUGHTER_PRC_NUM_X  :', DAUGHTER_PRC_NUM_X(HANDLING_NUM)
            if( IO_L ) write(IO_FID_LOG,'(1x,A,I6)'  ) '***  --- DAUGHTER_PRC_NUM_Y  :', DAUGHTER_PRC_NUM_Y(HANDLING_NUM)
            if( IO_L ) write(IO_FID_LOG,'(1x,A,I6)'  ) '***  --- DAUGHTER_KMAX       :', DAUGHTER_KMAX(HANDLING_NUM)
            if( IO_L ) write(IO_FID_LOG,'(1x,A,I6)'  ) '***  --- DAUGHTER_IMAX       :', DAUGHTER_IMAX(HANDLING_NUM)
            if( IO_L ) write(IO_FID_LOG,'(1x,A,I6)'  ) '***  --- DAUGHTER_JMAX       :', DAUGHTER_JMAX(HANDLING_NUM)
            if( IO_L ) write(IO_FID_LOG,'(1x,A,F9.3)') '***  --- DAUGHTER_DTSEC      :', DAUGHTER_DTSEC(HANDLING_NUM)
            if( IO_L ) write(IO_FID_LOG,'(1x,A,I6)'  ) '***  --- DAUGHTER_NSTEP      :', DAUGHTER_NSTEP(HANDLING_NUM)
            if( IO_L ) write(IO_FID_LOG,'(1x,A)'     ) '***  Informations of Target Tiles'
            if( IO_L ) write(IO_FID_LOG,'(1x,A,I6)'  ) '***  --- TILEALL_KA      :', TILEAL_KA(HANDLING_NUM)
            if( IO_L ) write(IO_FID_LOG,'(1x,A,I6)'  ) '***  --- TILEALL_IA      :', TILEAL_IA(HANDLING_NUM)
            if( IO_L ) write(IO_FID_LOG,'(1x,A,I6)'  ) '***  --- TILEALL_JA      :', TILEAL_JA(HANDLING_NUM)
            if( IO_L ) write(IO_FID_LOG,'(1x,A,I6)  ') '***  Limit Num. NCOMM req. :', max_rq

            allocate( buffer_2D  (                            PARENT_IA(HANDLING_NUM), PARENT_JA(HANDLING_NUM) ) )
            allocate( buffer_3D  (   PARENT_KA(HANDLING_NUM), PARENT_IA(HANDLING_NUM), PARENT_JA(HANDLING_NUM) ) )
            allocate( buffer_3DF ( 0:PARENT_KA(HANDLING_NUM), PARENT_IA(HANDLING_NUM), PARENT_JA(HANDLING_NUM) ) )

            allocate( recvbuf_3D (   PARENT_KA(HANDLING_NUM), PARENT_IA(HANDLING_NUM), PARENT_JA(HANDLING_NUM), max_isu  ) )
            allocate( recvbuf_3DF( 0:PARENT_KA(HANDLING_NUM), PARENT_IA(HANDLING_NUM), PARENT_JA(HANDLING_NUM), max_isuf ) )

            allocate( buffer_ref_LON (                          TILEAL_IA(HANDLING_NUM),TILEAL_JA(HANDLING_NUM)) )
            allocate( buffer_ref_LONU(                          TILEAL_IA(HANDLING_NUM),TILEAL_JA(HANDLING_NUM)) )
            allocate( buffer_ref_LONV(                          TILEAL_IA(HANDLING_NUM),TILEAL_JA(HANDLING_NUM)) )
            allocate( buffer_ref_LAT (                          TILEAL_IA(HANDLING_NUM),TILEAL_JA(HANDLING_NUM)) )
            allocate( buffer_ref_LATU(                          TILEAL_IA(HANDLING_NUM),TILEAL_JA(HANDLING_NUM)) )
            allocate( buffer_ref_LATV(                          TILEAL_IA(HANDLING_NUM),TILEAL_JA(HANDLING_NUM)) )
            allocate( buffer_ref_CZ  (  TILEAL_KA(HANDLING_NUM),TILEAL_IA(HANDLING_NUM),TILEAL_JA(HANDLING_NUM)) )
            allocate( buffer_ref_FZ  (0:TILEAL_KA(HANDLING_NUM),TILEAL_IA(HANDLING_NUM),TILEAL_JA(HANDLING_NUM)) )

            allocate( buffer_ref_3D  (  TILEAL_KA(HANDLING_NUM),TILEAL_IA(HANDLING_NUM),TILEAL_JA(HANDLING_NUM)) )
            allocate( buffer_ref_3DF (0:TILEAL_KA(HANDLING_NUM),TILEAL_IA(HANDLING_NUM),TILEAL_JA(HANDLING_NUM)) )

            allocate( igrd (                                 DAUGHTER_IA(HANDLING_NUM),DAUGHTER_JA(HANDLING_NUM),itp_nh,itp_ng) )
            allocate( jgrd (                                 DAUGHTER_IA(HANDLING_NUM),DAUGHTER_JA(HANDLING_NUM),itp_nh,itp_ng) )
            allocate( hfact(                                 DAUGHTER_IA(HANDLING_NUM),DAUGHTER_JA(HANDLING_NUM),itp_nh,itp_ng) )
            allocate( kgrd (DAUGHTER_KA(HANDLING_NUM),itp_nv,DAUGHTER_IA(HANDLING_NUM),DAUGHTER_JA(HANDLING_NUM),itp_nh,itp_ng) )
            allocate( vfact(DAUGHTER_KA(HANDLING_NUM),itp_nv,DAUGHTER_IA(HANDLING_NUM),DAUGHTER_JA(HANDLING_NUM),itp_nh,itp_ng) )

            call COMM_CARTESC_NEST_setup_nestdown( HANDLING_NUM )


            ! for scalar points
            call INTRP_factor3d( itp_nh,                            & ! [IN]
                                 TILEAL_KA  (HANDLING_NUM),         & ! [IN]
                                 1,                                 & ! [IN]
                                 TILEAL_KA  (HANDLING_NUM),         & ! [IN]
                                 TILEAL_IA  (HANDLING_NUM),         & ! [IN]
                                 TILEAL_JA  (HANDLING_NUM),         & ! [IN]
                                 buffer_ref_LON (:,:),              & ! [IN]
                                 buffer_ref_LAT (:,:),              & ! [IN]
                                 buffer_ref_CZ  (:,:,:),            & ! [IN]
                                 DAUGHTER_KA(HANDLING_NUM),         & ! [IN]
                                 DATR_KS    (HANDLING_NUM),         & ! [IN]
                                 DATR_KE    (HANDLING_NUM),         & ! [IN]
                                 DAUGHTER_IA(HANDLING_NUM),         & ! [IN]
                                 DAUGHTER_JA(HANDLING_NUM),         & ! [IN]
                                 ATMOS_GRID_CARTESC_REAL_LON       (:,:),              & ! [IN]
                                 ATMOS_GRID_CARTESC_REAL_LAT       (:,:),              & ! [IN]
                                 ATMOS_GRID_CARTESC_REAL_CZ        (:,:,:),            & ! [IN]
                                 igrd           (    :,:,:,I_SCLR), & ! [OUT]
                                 jgrd           (    :,:,:,I_SCLR), & ! [OUT]
                                 hfact          (    :,:,:,I_SCLR), & ! [OUT]
                                 kgrd           (:,:,:,:,:,I_SCLR), & ! [OUT]
                                 vfact          (:,:,:,:,:,I_SCLR)  ) ! [OUT]

            ! for z staggered points
            call INTRP_factor3d( itp_nh,                            & ! [IN]
                                 TILEAL_KA  (HANDLING_NUM)+1,       & ! [IN]
                                 1,                                 & ! [IN]
                                 TILEAL_KA  (HANDLING_NUM)+1,       & ! [IN]
                                 TILEAL_IA  (HANDLING_NUM),         & ! [IN]
                                 TILEAL_JA  (HANDLING_NUM),         & ! [IN]
                                 buffer_ref_LON (:,:),              & ! [IN]
                                 buffer_ref_LAT (:,:),              & ! [IN]
                                 buffer_ref_FZ  (:,:,:),            & ! [IN]
                                 DAUGHTER_KA(HANDLING_NUM),         & ! [IN]
                                 DATR_KS    (HANDLING_NUM),         & ! [IN]
                                 DATR_KE    (HANDLING_NUM),         & ! [IN]
                                 DAUGHTER_IA(HANDLING_NUM),         & ! [IN]
                                 DAUGHTER_JA(HANDLING_NUM),         & ! [IN]
                                 ATMOS_GRID_CARTESC_REAL_LON       (:,:),              & ! [IN]
                                 ATMOS_GRID_CARTESC_REAL_LAT       (:,:),              & ! [IN]
                                 ATMOS_GRID_CARTESC_REAL_FZ        (1:KA,:,:),         & ! [IN]
                                 igrd           (    :,:,:,I_ZSTG), & ! [OUT]
                                 jgrd           (    :,:,:,I_ZSTG), & ! [OUT]
                                 hfact          (    :,:,:,I_ZSTG), & ! [OUT]
                                 kgrd           (:,:,:,:,:,I_ZSTG), & ! [OUT]
                                 vfact          (:,:,:,:,:,I_ZSTG)  ) ! [OUT]

            ! for x staggered points
            call INTRP_factor3d( itp_nh,                            & ! [IN]
                                 TILEAL_KA  (HANDLING_NUM),         & ! [IN]
                                 1,                                 & ! [IN]
                                 TILEAL_KA  (HANDLING_NUM),         & ! [IN]
                                 TILEAL_IA  (HANDLING_NUM),         & ! [IN]
                                 TILEAL_JA  (HANDLING_NUM),         & ! [IN]
                                 buffer_ref_LONU(:,:),              & ! [IN]
                                 buffer_ref_LATV(:,:),              & ! [IN]
                                 buffer_ref_CZ  (:,:,:),            & ! [IN]
                                 DAUGHTER_KA(HANDLING_NUM),         & ! [IN]
                                 DATR_KS    (HANDLING_NUM),         & ! [IN]
                                 DATR_KE    (HANDLING_NUM),         & ! [IN]
                                 DAUGHTER_IA(HANDLING_NUM),         & ! [IN]
                                 DAUGHTER_JA(HANDLING_NUM),         & ! [IN]
                                 ATMOS_GRID_CARTESC_REAL_LONU      (1:IA,1:JA),        & ! [IN]
                                 ATMOS_GRID_CARTESC_REAL_LATU      (1:IA,1:JA),        & ! [IN]
                                 ATMOS_GRID_CARTESC_REAL_CZ        (:,:,:),            & ! [IN]
                                 igrd           (    :,:,:,I_XSTG), & ! [OUT]
                                 jgrd           (    :,:,:,I_XSTG), & ! [OUT]
                                 hfact          (    :,:,:,I_XSTG), & ! [OUT]
                                 kgrd           (:,:,:,:,:,I_XSTG), & ! [OUT]
                                 vfact          (:,:,:,:,:,I_XSTG)  ) ! [OUT]

            ! for y staggered points
            call INTRP_factor3d( itp_nh,                            & ! [IN]
                                 TILEAL_KA  (HANDLING_NUM),         & ! [IN]
                                 1,                                 & ! [IN]
                                 TILEAL_KA  (HANDLING_NUM),         & ! [IN]
                                 TILEAL_IA  (HANDLING_NUM),         & ! [IN]
                                 TILEAL_JA  (HANDLING_NUM),         & ! [IN]
                                 buffer_ref_LONV(:,:),              & ! [IN]
                                 buffer_ref_LATV(:,:),              & ! [IN]
                                 buffer_ref_CZ  (:,:,:),            & ! [IN]
                                 DAUGHTER_KA(HANDLING_NUM),         & ! [IN]
                                 DATR_KS    (HANDLING_NUM),         & ! [IN]
                                 DATR_KE    (HANDLING_NUM),         & ! [IN]
                                 DAUGHTER_IA(HANDLING_NUM),         & ! [IN]
                                 DAUGHTER_JA(HANDLING_NUM),         & ! [IN]
                                 ATMOS_GRID_CARTESC_REAL_LONV      (1:IA,1:JA),        & ! [IN]
                                 ATMOS_GRID_CARTESC_REAL_LATV      (1:IA,1:JA),        & ! [IN]
                                 ATMOS_GRID_CARTESC_REAL_CZ        (:,:,:),            & ! [IN]
                                 igrd           (    :,:,:,I_YSTG), & ! [OUT]
                                 jgrd           (    :,:,:,I_YSTG), & ! [OUT]
                                 hfact          (    :,:,:,I_YSTG), & ! [OUT]
                                 kgrd           (:,:,:,:,:,I_YSTG), & ! [OUT]
                                 vfact          (:,:,:,:,:,I_YSTG)  ) ! [OUT]

            deallocate( buffer_2D  )
            deallocate( buffer_3D  )
            deallocate( buffer_3DF )

         else
            ONLINE_USE_VELZ = .false.
         endif

         !if( IO_L ) write(IO_FID_LOG,'(1x,A,I2)') '*** Number of Related Domains :', HANDLING_NUM
         !if ( HANDLING_NUM > 2 ) then
         !   f( IO_L ) write(*,*) 'xxx Too much handing domains (up to 2)'
         !   call PRC_MPIstop
         !endif

      endif !--- OFFLINE or NOT

    endif !--- USE_NESTING

    return
  end subroutine COMM_CARTESC_NEST_setup

  !-----------------------------------------------------------------------------
  !> Solve relationship between ParentDomain & Daughter Domain
  subroutine COMM_CARTESC_NEST_domain_relate( &
       HANDLE )
    use scale_process, only: &
       PRC_myrank,   &
       PRC_MPIstop
    implicit none

    integer, intent(in) :: HANDLE !< id number of nesting relation in this process target

    logical              :: hit = .false.
    integer, allocatable :: pd_tile_num(:,:)

    real(RP) :: wid_lon, wid_lat
    integer  :: pd_sw_tile
    integer  :: pd_ne_tile
    integer  :: i, j, k
    !---------------------------------------------------------------------------

    if( .NOT. USE_NESTING ) return

    allocate( pd_tile_num(0:PARENT_PRC_nprocs(HANDLE)-1,2) )

    k = 0 ! MPI process number starts from zero
    do j = 1, PARENT_PRC_NUM_Y(HANDLE)
    do i = 1, PARENT_PRC_NUM_X(HANDLE)
       pd_tile_num(k,1) = i
       pd_tile_num(k,2) = j
       k = k + 1
    enddo
    enddo

    !--- SW search
    hit = .false.
    do i = 1, PARENT_PRC_nprocs(HANDLE)
       wid_lon = abs((latlon_catalog(i,I_SW,I_LON) - latlon_catalog(i,I_SE,I_LON)) &
                      / real( PARENT_IMAX(HANDLE)-1, kind=RP )) * 0.8_RP
       wid_lat = abs((latlon_catalog(i,I_SW,I_LAT) - latlon_catalog(i,I_NW,I_LAT)) &
                      / real( PARENT_JMAX(HANDLE)-1, kind=RP )) * 0.8_RP

       if ( corner_loc(I_SW,I_LON) >= min(latlon_catalog(i,I_SW,I_LON),latlon_catalog(i,I_NW,I_LON))-wid_lon .AND. &
            corner_loc(I_SW,I_LAT) >= min(latlon_catalog(i,I_SW,I_LAT),latlon_catalog(i,I_SE,I_LAT))-wid_lat .AND. &
            corner_loc(I_SW,I_LON) <= max(latlon_catalog(i,I_NE,I_LON),latlon_catalog(i,I_SE,I_LON))+wid_lon .AND. &
            corner_loc(I_SW,I_LAT) <= max(latlon_catalog(i,I_NE,I_LAT),latlon_catalog(i,I_NW,I_LAT))+wid_lat ) then

          pd_sw_tile = i-1 ! MPI process number starts from zero
          hit = .true.
          exit ! exit loop
       endif
    enddo
    if ( .NOT. hit ) then
       write(*,*) 'xxx [NEST_domain_relate] region of daughter domain is larger than that of parent: SW search'
       write(*,*) '                                  at rank:', PRC_myrank, ' of domain:', ONLINE_DOMAIN_NUM
       if( IO_L ) write(IO_FID_LOG,'(1x,A)') 'xxx region of daughter domain is larger than that of parent: SW search'
       if( IO_L ) write(IO_FID_LOG,*) ' grid width: half width in lat:', wid_lat, ' half width in lon:', wid_lon
       if( IO_L ) write(IO_FID_LOG,'(1x,A,F12.6)') '    daughter local (me): LON=',corner_loc(I_SW,I_LON)
       do i = 1, PARENT_PRC_nprocs(HANDLE)
          if( IO_L ) write(IO_FID_LOG,'(1x,A,F12.6,1x,F12.6)') '     parent local SW-NE: LON=', &
                     latlon_catalog(i,I_SW,I_LON) ,latlon_catalog(i,I_NE,I_LON)
       enddo
       if( IO_L ) write(IO_FID_LOG,'(1x,A,F12.6)') '    daughter local (me): LAT=',corner_loc(I_SW,I_LAT)
       do i = 1, PARENT_PRC_nprocs(HANDLE)
          if( IO_L ) write(IO_FID_LOG,'(1x,A,F12.6,1x,F12.6)') '     parent local SW-NE: LAT=', &
                     latlon_catalog(i,I_SW,I_LAT) ,latlon_catalog(i,I_NE,I_LAT)
       enddo
       call PRC_MPIstop
    endif

    !--- NE search
    hit = .false.
    do i = PARENT_PRC_nprocs(HANDLE), 1, -1
       wid_lon = abs((latlon_catalog(i,I_NW,I_LON) - latlon_catalog(i,I_NE,I_LON)) &
                      / real( PARENT_IMAX(HANDLE)-1, kind=RP )) * 0.8_RP
       wid_lat = abs((latlon_catalog(i,I_SE,I_LAT) - latlon_catalog(i,I_NE,I_LAT)) &
                      / real( PARENT_JMAX(HANDLE)-1, kind=RP )) * 0.8_RP

       if ( corner_loc(I_NE,I_LON) >= min(latlon_catalog(i,I_SW,I_LON),latlon_catalog(i,I_NW,I_LON))-wid_lon .AND. &
            corner_loc(I_NE,I_LAT) >= min(latlon_catalog(i,I_SW,I_LAT),latlon_catalog(i,I_SE,I_LAT))-wid_lat .AND. &
            corner_loc(I_NE,I_LON) <= max(latlon_catalog(i,I_NE,I_LON),latlon_catalog(i,I_SE,I_LON))+wid_lon .AND. &
            corner_loc(I_NE,I_LAT) <= max(latlon_catalog(i,I_NE,I_LAT),latlon_catalog(i,I_NW,I_LAT))+wid_lat ) then

          pd_ne_tile = i-1 ! MPI process number starts from zero
          hit = .true.
          exit ! exit loop
       endif
    enddo
    if ( .NOT. hit ) then
       write(*,*) 'xxx [NEST_domain_relate] region of daughter domain is larger than that of parent: NE search'
       write(*,*) '                                  at rank:', PRC_myrank, ' of domain:', ONLINE_DOMAIN_NUM
       if( IO_L ) write(IO_FID_LOG,'(1x,A)') 'xxx region of daughter domain is larger than that of parent: NE search'
       if( IO_L ) write(IO_FID_LOG,*) ' grid width: half width in lat:', wid_lat, ' half width in lon:', wid_lon
       if( IO_L ) write(IO_FID_LOG,'(1x,A,F12.6)') '    daughter local (me): LON=',corner_loc(I_NE,I_LON)
       do i = 1, PARENT_PRC_nprocs(HANDLE)
          if( IO_L ) write(IO_FID_LOG,'(1x,A,F12.6,1x,F12.6)') '     parent local SW-NE: LON=', &
                     latlon_catalog(i,I_SW,I_LON) ,latlon_catalog(i,I_NE,I_LON)
       enddo
       if( IO_L ) write(IO_FID_LOG,'(1x,A,F12.6)') '    daughter local (me): LAT=',corner_loc(I_NE,I_LAT)
       do i = 1, PARENT_PRC_nprocs(HANDLE)
          if( IO_L ) write(IO_FID_LOG,'(1x,A,F12.6,1x,F12.6)') '     parent local SW-NE: LAT=', &
                     latlon_catalog(i,I_SW,I_LAT) ,latlon_catalog(i,I_NE,I_LAT)
       enddo
       call PRC_MPIstop
    endif

    COMM_CARTESC_NEST_TILE_NUM_X = pd_tile_num(pd_ne_tile,1) - pd_tile_num(pd_sw_tile,1) + 1
    COMM_CARTESC_NEST_TILE_NUM_Y = pd_tile_num(pd_ne_tile,2) - pd_tile_num(pd_sw_tile,2) + 1

    allocate( COMM_CARTESC_NEST_TILE_ID( COMM_CARTESC_NEST_TILE_NUM_X*COMM_CARTESC_NEST_TILE_NUM_Y ) )

    if( IO_L ) write(IO_FID_LOG,'(1x,A)') '*** NEST: target process tile in parent domain'
    k = 1
    do j = 1, COMM_CARTESC_NEST_TILE_NUM_Y
    do i = 1, COMM_CARTESC_NEST_TILE_NUM_X
       COMM_CARTESC_NEST_TILE_ID(k) = pd_sw_tile + (i-1) + PARENT_PRC_NUM_X(HANDLE)*(j-1)
       if( IO_L ) write(IO_FID_LOG,'(1x,A,I4,A,I6)') '    (', k, ') target mpi-process:', COMM_CARTESC_NEST_TILE_ID(k)
       k = k + 1
    enddo
    enddo

    return
  end subroutine COMM_CARTESC_NEST_domain_relate

  !-----------------------------------------------------------------------------
  !> Return shape of ParentDomain at the specified rank (for offline)
  !  including definition array size with BND or not in Parent domain
  subroutine COMM_CARTESC_NEST_domain_shape ( &
       tilei,    &
       tilej,    &
       cxs, cxe, &
       cys, cye, &
       pxs, pxe, &
       pys, pye, &
       iloc      )
    implicit none

    integer, intent(out) :: tilei, tilej
    integer, intent(out) :: cxs, cxe, cys, cye
    integer, intent(out) :: pxs, pxe, pys, pye
    integer, intent(in)  :: iloc                ! rank number; start from 1

    integer :: hdl = 1      ! handler number
    integer :: rank
    integer :: xloc,  yloc
    integer :: xlocg, ylocg ! location over whole domain
    !---------------------------------------------------------------------------

    if( .NOT. USE_NESTING ) return

    rank  = COMM_CARTESC_NEST_TILE_ID(iloc)
    xloc  = mod( iloc-1, COMM_CARTESC_NEST_TILE_NUM_X ) + 1
    yloc  = int( real(iloc-1) / real(COMM_CARTESC_NEST_TILE_NUM_X) ) + 1
    xlocg = mod( rank, OFFLINE_PARENT_PRC_NUM_X ) + 1
    ylocg = int( real(rank) / real(OFFLINE_PARENT_PRC_NUM_X) ) + 1
    tilei = PARENT_IMAX(hdl)
    tilej = PARENT_JMAX(hdl)

    cxs   = tilei * (xloc-1) + 1
    cxe   = tilei * xloc
    cys   = tilej * (yloc-1) + 1
    cye   = tilej * yloc
    pxs   = 1
    pxe   = tilei
    pys   = 1
    pye   = tilej

    if ( xlocg == 1 ) then ! BND_W
       tilei = tilei + 2
       pxs = pxs + 2
       pxe = pxe + 2
    endif
    if ( xlocg == OFFLINE_PARENT_PRC_NUM_X ) then ! BND_E
       tilei = tilei + 2
    endif
    if ( ylocg == 1 ) then ! BND_S
       tilej = tilej + 2
       pys = pys + 2
       pye = pye + 2
    endif
    if ( ylocg == OFFLINE_PARENT_PRC_NUM_Y ) then ! BND_N
       tilej = tilej + 2
    endif

    return
  end subroutine COMM_CARTESC_NEST_domain_shape

  !-----------------------------------------------------------------------------
  !> Get parent domain size
  subroutine COMM_CARTESC_NEST_parentsize( &
       HANDLE )
    use scale_process, only: &
       PRC_MPIstop, &
       PRC_nprocs,  &
       PRC_myrank,  &
       PRC_IsMaster
    use scale_rm_process, only: &
       PRC_NUM_X,   &
       PRC_NUM_Y
    use scale_time, only: &
       TIME_NSTEP, &
       TIME_DTSEC
    use scale_comm, only: &
       COMM_bcast
    implicit none

    integer, intent(in) :: HANDLE !< id number of nesting relation in this process target

    real(RP) :: buffer
    integer  :: datapack(14)
    integer  :: QA_OTHERSIDE
    integer  :: ireq1, ireq2, ierr1, ierr2, ileng
    integer  :: istatus(MPI_STATUS_SIZE)
    integer  :: tag
    !---------------------------------------------------------------------------

    if( .NOT. USE_NESTING ) return

    tag   = INTERCOMM_ID(HANDLE) * 100
    ileng = 14

    if ( COMM_CARTESC_NEST_Filiation( INTERCOMM_ID(HANDLE) ) > 0 ) then

       !##### parent ####

       ! from parent to daughter
       datapack( 1) = PRC_nprocs
       datapack( 2) = PRC_NUM_X
       datapack( 3) = PRC_NUM_Y
       datapack( 4) = KMAX
       datapack( 5) = IMAX
       datapack( 6) = JMAX
       datapack( 7) = KS
       datapack( 8) = KE
       datapack( 9) = IS
       datapack(10) = IE
       datapack(11) = JS
       datapack(12) = JE
       datapack(13) = TIME_NSTEP
       datapack(14) = COMM_CARTESC_NEST_BND_QA
       buffer       = TIME_DTSEC

       if ( PRC_IsMaster ) then
          call MPI_ISEND(datapack, ileng, MPI_INTEGER, PRC_myrank, tag, INTERCOMM_DAUGHTER, ireq1, ierr1)
          call MPI_ISEND(buffer, 1, MPI_DOUBLE_PRECISION, PRC_myrank, tag+1, INTERCOMM_DAUGHTER, ireq2, ierr2)
          call MPI_WAIT(ireq1, istatus, ierr1)
          call MPI_WAIT(ireq2, istatus, ierr2)
       endif

       PARENT_PRC_nprocs(HANDLE) = datapack( 1)
       PARENT_PRC_NUM_X (HANDLE) = datapack( 2)
       PARENT_PRC_NUM_Y (HANDLE) = datapack( 3)
       PARENT_KMAX      (HANDLE) = datapack( 4)
       PARENT_IMAX      (HANDLE) = datapack( 5)
       PARENT_JMAX      (HANDLE) = datapack( 6)
       PRNT_KS          (HANDLE) = datapack( 7)
       PRNT_KE          (HANDLE) = datapack( 8)
       PRNT_IS          (HANDLE) = datapack( 9)
       PRNT_IE          (HANDLE) = datapack(10)
       PRNT_JS          (HANDLE) = datapack(11)
       PRNT_JE          (HANDLE) = datapack(12)
       PARENT_NSTEP     (HANDLE) = datapack(13)
       PARENT_DTSEC     (HANDLE) = buffer

       ! from daughter to parent
       if ( PRC_IsMaster ) then
          call MPI_IRECV(datapack, ileng, MPI_INTEGER, PRC_myrank, tag+2, INTERCOMM_DAUGHTER, ireq1, ierr1)
          call MPI_IRECV(buffer, 1, MPI_DOUBLE_PRECISION, PRC_myrank, tag+3, INTERCOMM_DAUGHTER, ireq2, ierr2)
          call MPI_WAIT(ireq1, istatus, ierr1)
          call MPI_WAIT(ireq2, istatus, ierr2)
       endif
       call COMM_bcast(datapack, ileng)
       call COMM_bcast(buffer)

       DAUGHTER_PRC_nprocs(HANDLE) = datapack( 1)
       DAUGHTER_PRC_NUM_X (HANDLE) = datapack( 2)
       DAUGHTER_PRC_NUM_Y (HANDLE) = datapack( 3)
       DAUGHTER_KMAX      (HANDLE) = datapack( 4)
       DAUGHTER_IMAX      (HANDLE) = datapack( 5)
       DAUGHTER_JMAX      (HANDLE) = datapack( 6)
       DATR_KS            (HANDLE) = datapack( 7)
       DATR_KE            (HANDLE) = datapack( 8)
       DATR_IS            (HANDLE) = datapack( 9)
       DATR_IE            (HANDLE) = datapack(10)
       DATR_JS            (HANDLE) = datapack(11)
       DATR_JE            (HANDLE) = datapack(12)
       DAUGHTER_NSTEP     (HANDLE) = datapack(13)
       QA_OTHERSIDE                = datapack(14)
       DAUGHTER_DTSEC     (HANDLE) = buffer


    elseif( COMM_CARTESC_NEST_Filiation( INTERCOMM_ID(HANDLE) ) < 0 ) then

       !##### child ####

       ! from parent to daughter
       if ( PRC_IsMaster ) then
          call MPI_IRECV(datapack, ileng, MPI_INTEGER, PRC_myrank, tag, INTERCOMM_PARENT, ireq1, ierr1)
          call MPI_IRECV(buffer, 1, MPI_DOUBLE_PRECISION, PRC_myrank, tag+1, INTERCOMM_PARENT, ireq2, ierr2)
          call MPI_WAIT(ireq1, istatus, ierr1)
          call MPI_WAIT(ireq2, istatus, ierr2)
       endif
       call COMM_bcast(datapack, ileng)
       call COMM_bcast(buffer)

       PARENT_PRC_nprocs(HANDLE) = datapack( 1)
       PARENT_PRC_NUM_X (HANDLE) = datapack( 2)
       PARENT_PRC_NUM_Y (HANDLE) = datapack( 3)
       PARENT_KMAX      (HANDLE) = datapack( 4)
       PARENT_IMAX      (HANDLE) = datapack( 5)
       PARENT_JMAX      (HANDLE) = datapack( 6)
       PRNT_KS          (HANDLE) = datapack( 7)
       PRNT_KE          (HANDLE) = datapack( 8)
       PRNT_IS          (HANDLE) = datapack( 9)
       PRNT_IE          (HANDLE) = datapack(10)
       PRNT_JS          (HANDLE) = datapack(11)
       PRNT_JE          (HANDLE) = datapack(12)
       PARENT_NSTEP     (HANDLE) = datapack(13)
       QA_OTHERSIDE              = datapack(14)
       PARENT_DTSEC     (HANDLE) = buffer

       ! from daughter to parent
       datapack( 1) = PRC_nprocs
       datapack( 2) = PRC_NUM_X
       datapack( 3) = PRC_NUM_Y
       datapack( 4) = KMAX
       datapack( 5) = IMAX
       datapack( 6) = JMAX
       datapack( 7) = KS
       datapack( 8) = KE
       datapack( 9) = IS
       datapack(10) = IE
       datapack(11) = JS
       datapack(12) = JE
       datapack(13) = TIME_NSTEP
       datapack(14) = COMM_CARTESC_NEST_BND_QA
       buffer       = TIME_DTSEC

       if ( PRC_IsMaster ) then
          call MPI_ISEND(datapack, ileng, MPI_INTEGER, PRC_myrank, tag+2, INTERCOMM_PARENT, ireq1, ierr1)
          call MPI_ISEND(buffer, 1, MPI_DOUBLE_PRECISION, PRC_myrank, tag+3, INTERCOMM_PARENT, ireq2, ierr2)
          call MPI_WAIT(ireq1, istatus, ierr1)
          call MPI_WAIT(ireq2, istatus, ierr2)
       endif

       DAUGHTER_PRC_nprocs(HANDLE) = datapack( 1)
       DAUGHTER_PRC_NUM_X (HANDLE) = datapack( 2)
       DAUGHTER_PRC_NUM_Y (HANDLE) = datapack( 3)
       DAUGHTER_KMAX      (HANDLE) = datapack( 4)
       DAUGHTER_IMAX      (HANDLE) = datapack( 5)
       DAUGHTER_JMAX      (HANDLE) = datapack( 6)
       DATR_KS            (HANDLE) = datapack( 7)
       DATR_KE            (HANDLE) = datapack( 8)
       DATR_IS            (HANDLE) = datapack( 9)
       DATR_IE            (HANDLE) = datapack(10)
       DATR_JS            (HANDLE) = datapack(11)
       DATR_JE            (HANDLE) = datapack(12)
       DAUGHTER_NSTEP     (HANDLE) = datapack(13)
       DAUGHTER_DTSEC     (HANDLE) = buffer
    else
       write(*,*) 'xxx [COMM_CARTESC_NEST_parentsize] internal error'
       call PRC_MPIstop
    endif

    if ( ONLINE_BOUNDARY_DIAGQNUM ) then
       if( IO_L ) write(IO_FID_LOG,*) '*** Number concentration of hydrometeor will be diagnosed'
       if( IO_L ) write(IO_FID_LOG,*) '*** Number of QA (remote,local) = ', QA_OTHERSIDE, COMM_CARTESC_NEST_BND_QA
       COMM_CARTESC_NEST_BND_QA = min(QA_OTHERSIDE, COMM_CARTESC_NEST_BND_QA)
    else
       if ( QA_OTHERSIDE /= COMM_CARTESC_NEST_BND_QA ) then
          write(*,*) 'xxx [COMM_CARTESC_NEST_parentsize] NUMBER of QA are not matched!'
          write(*,*) 'xxx check a flag of ONLINE_BOUNDARY_USE_QHYD.'
          write(*,*) 'xxx Number of QA (remote,local) = ', QA_OTHERSIDE, COMM_CARTESC_NEST_BND_QA
          call PRC_MPIstop
       endif
    endif

    return
  end subroutine COMM_CARTESC_NEST_parentsize

  !-----------------------------------------------------------------------------
  !> Get parent latlon catalogue
  subroutine COMM_CARTESC_NEST_catalogue( &
       HANDLE  )
    use scale_process, only: &
       PRC_MPIstop, &
       PRC_nprocs,  &
       PRC_myrank,  &
       PRC_IsMaster
    use scale_comm, only: &
       COMM_datatype,  &
       COMM_bcast
    use scale_atmos_grid_cartesC_real, only: &
       ATMOS_GRID_CARTESC_REAL_DOMAIN_CATALOGUE
    implicit none

    integer, intent(in) :: HANDLE !< id number of nesting relation in this process target

    integer :: ireq, ierr, ileng
    integer :: istatus(MPI_STATUS_SIZE)
    integer :: tag
    !---------------------------------------------------------------------------

    if( .NOT. USE_NESTING ) return

    tag = INTERCOMM_ID(HANDLE) * 100

    if ( COMM_CARTESC_NEST_Filiation( INTERCOMM_ID(HANDLE) ) > 0 ) then

       !##### parent ####

       ileng = PRC_nprocs * 4 * 2

       if ( PRC_IsMaster ) then
          call MPI_ISEND(ATMOS_GRID_CARTESC_REAL_DOMAIN_CATALOGUE, ileng, COMM_datatype, PRC_myrank, tag, INTERCOMM_DAUGHTER, ireq, ierr)
          call MPI_WAIT(ireq, istatus, ierr)
       endif

    elseif( COMM_CARTESC_NEST_Filiation( INTERCOMM_ID(HANDLE) ) < 0 ) then

       !##### child ####

       ileng = PARENT_PRC_nprocs(HANDLE) * 4 * 2

       if ( PRC_IsMaster ) then
          call MPI_IRECV(latlon_catalog, ileng, COMM_datatype, PRC_myrank, tag, INTERCOMM_PARENT, ireq, ierr)
          call MPI_WAIT(ireq, istatus, ierr)
       endif
       call COMM_bcast( latlon_catalog, PARENT_PRC_nprocs(HANDLE), 4, 2 )

    else
       write(*,*) 'xxx [COMM_CARTESC_NEST_catalogue] internal error'
       call PRC_MPIstop
    endif

    return
  end subroutine COMM_CARTESC_NEST_catalogue

  !-----------------------------------------------------------------------------
  !> Check Communication Inter-domains
  subroutine COMM_CARTESC_NEST_ping( &
       HANDLE )
    use scale_process, only: &
       PRC_MPIstop, &
       PRC_myrank,  &
       PRC_IsMaster
    use scale_comm, only: &
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
          call MPI_ISEND(ping, 1, MPI_INTEGER, PRC_myrank, tag+1, INTERCOMM_DAUGHTER, ireq1, ierr1)
          call MPI_IRECV(pong, 1, MPI_INTEGER, PRC_myrank, tag+2, INTERCOMM_DAUGHTER, ireq2, ierr2)
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
          call MPI_ISEND(ping, 1, MPI_INTEGER, PRC_myrank, tag+2, INTERCOMM_PARENT, ireq1, ierr1)
          call MPI_IRECV(pong, 1, MPI_INTEGER, PRC_myrank, tag+1, INTERCOMM_PARENT, ireq2, ierr2)
          call MPI_WAIT(ireq1, istatus, ierr1)
          call MPI_WAIT(ireq2, istatus, ierr2)
       endif

       call COMM_bcast(pong)

       if ( pong /= INTERCOMM_ID(HANDLE) ) ping_error = .true.

    else
       write(*,*) 'xxx [COMM_CARTESC_NEST_ping] internal error'
       call PRC_MPIstop
    endif

    if ( ping_error ) then
       write(*,*) 'xxx [COMM_CARTESC_NEST_ping] ping destination error'
       call PRC_MPIstop
    endif

    return
  end subroutine COMM_CARTESC_NEST_ping

  !-----------------------------------------------------------------------------
  !> Inter-domain communication setup for nestdown
  subroutine COMM_CARTESC_NEST_setup_nestdown( &
       HANDLE )
    use scale_process, only: &
       PRC_MPIstop, &
       PRC_myrank,  &
       PRC_IsMaster
    use scale_comm, only: &
       COMM_world, &
       COMM_bcast
    implicit none

    integer, intent(in) :: HANDLE !< id number of nesting relation in this process target

    integer, allocatable :: buffer_LIST   (:)
    integer, allocatable :: buffer_ALLLIST(:)

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
          call MPI_IRECV(COMM_CARTESC_NEST_TILE_ALLMAX_p, 1, MPI_INTEGER, PRC_myrank, tag+1, INTERCOMM_DAUGHTER, ireq, ierr)
          call MPI_WAIT(ireq, istatus, ierr)
       endif
       call COMM_bcast(COMM_CARTESC_NEST_TILE_ALLMAX_p)

       allocate( COMM_CARTESC_NEST_TILE_LIST_p (COMM_CARTESC_NEST_TILE_ALLMAX_p,DAUGHTER_PRC_nprocs(HANDLE)) )
       allocate( COMM_CARTESC_NEST_TILE_LIST_YP(COMM_CARTESC_NEST_TILE_ALLMAX_p*DAUGHTER_PRC_nprocs(HANDLE)) )

       ileng = COMM_CARTESC_NEST_TILE_ALLMAX_p*DAUGHTER_PRC_nprocs(HANDLE)
       if ( PRC_IsMaster ) then
          call MPI_IRECV(COMM_CARTESC_NEST_TILE_LIST_p, ileng, MPI_INTEGER, PRC_myrank, tag+2, INTERCOMM_DAUGHTER, ireq, ierr)
          call MPI_WAIT(ireq, istatus, ierr)
       endif
       call COMM_bcast(COMM_CARTESC_NEST_TILE_LIST_p, COMM_CARTESC_NEST_TILE_ALLMAX_p, DAUGHTER_PRC_nprocs(HANDLE))

       COMM_CARTESC_NEST_TILE_LIST_YP(:) = -1

       k = 0
       do j = 1, DAUGHTER_PRC_nprocs(HANDLE)
       do i = 1, COMM_CARTESC_NEST_TILE_ALLMAX_p
          if ( COMM_CARTESC_NEST_TILE_LIST_p(i,j) == PRC_myrank ) then
             k = k + 1
             COMM_CARTESC_NEST_TILE_LIST_YP(k) = j - 1  !rank number is started from 1
          endif
       enddo
       enddo
       NUM_YP = k

       if( IO_L ) write(IO_FID_LOG,'(A,I5,A,I5)') "[P]   Num YP =",NUM_YP,"  Num TILE(MAX) =",COMM_CARTESC_NEST_TILE_ALLMAX_p

       if ( PRC_IsMaster ) then
          call MPI_IRECV(ONLINE_DAUGHTER_USE_VELZ, 1, MPI_LOGICAL, PRC_myrank, tag+3, INTERCOMM_DAUGHTER, ireq, ierr)
          call MPI_WAIT(ireq, istatus, ierr)
       endif
       call COMM_bcast(ONLINE_DAUGHTER_USE_VELZ)

       if( IO_L ) write(IO_FID_LOG,'(1x,A,L2)') '*** NEST: ONLINE_DAUGHTER_USE_VELZ =', ONLINE_DAUGHTER_USE_VELZ

       if ( PRC_IsMaster ) then
          call MPI_IRECV(ONLINE_DAUGHTER_NO_ROTATE, 1, MPI_LOGICAL, PRC_myrank, tag+4, INTERCOMM_DAUGHTER, ireq, ierr)
          call MPI_WAIT(ireq, istatus, ierr)
       endif
       call COMM_bcast(ONLINE_DAUGHTER_NO_ROTATE)

       if( ONLINE_NO_ROTATE .neqv. ONLINE_DAUGHTER_NO_ROTATE ) then
          write(*,*) 'xxx [COMM_CARTESC_NEST_setup_nestdown] Flag of NO_ROTATE is not consistent with the child domain'
          if( IO_L ) write(IO_FID_LOG,*) 'xxx ONLINE_NO_ROTATE = ', ONLINE_NO_ROTATE
          if( IO_L ) write(IO_FID_LOG,*) 'xxx ONLINE_DAUGHTER_NO_ROTATE =', ONLINE_DAUGHTER_NO_ROTATE
          call PRC_MPIstop
       endif
       if( IO_L ) write(IO_FID_LOG,'(1x,A,L2)') '*** NEST: ONLINE_DAUGHTER_NO_ROTATE =', ONLINE_DAUGHTER_NO_ROTATE

       call COMM_CARTESC_NEST_importgrid_nestdown( HANDLE )

       do i = 1, NUM_YP
          target_rank = COMM_CARTESC_NEST_TILE_LIST_YP(i)
          call MPI_ISEND(i, 1, MPI_INTEGER, target_rank, tag+5, INTERCOMM_DAUGHTER, ireq, ierr)
          call MPI_WAIT(ireq, istatus, ierr)
       enddo

       call MPI_BARRIER(INTERCOMM_DAUGHTER, ierr)

    elseif( COMM_CARTESC_NEST_Filiation( INTERCOMM_ID(HANDLE) ) < 0 ) then

       !##### child ####

       COMM_CARTESC_NEST_TILE_ALL = size( COMM_CARTESC_NEST_TILE_ID(:) ) ! should be equal to "NEST_TILE_NUM_X*NEST_TILE_NUM_Y"
       call MPI_Allreduce( COMM_CARTESC_NEST_TILE_ALL,      &
                           COMM_CARTESC_NEST_TILE_ALLMAX_d, &
                           1,                  &
                           MPI_INTEGER,        &
                           MPI_MAX,            &
                           COMM_world,         &
                           ierr                )
       if( IO_L ) write(IO_FID_LOG,'(A,I5,A,I5)') "[D]   Num YP =",COMM_CARTESC_NEST_TILE_ALL,"  Num TILE(MAX) =",COMM_CARTESC_NEST_TILE_ALLMAX_d

       if ( PRC_IsMaster ) then
          call MPI_ISEND(COMM_CARTESC_NEST_TILE_ALLMAX_d, 1, MPI_INTEGER, PRC_myrank, tag+1, INTERCOMM_PARENT, ireq, ierr)
          call MPI_WAIT(ireq, istatus, ierr)
       endif

       allocate( buffer_LIST   (COMM_CARTESC_NEST_TILE_ALLMAX_d)            )
       allocate( buffer_ALLLIST(COMM_CARTESC_NEST_TILE_ALLMAX_d*DAUGHTER_PRC_nprocs(HANDLE))   )
       allocate( COMM_CARTESC_NEST_TILE_LIST_d(COMM_CARTESC_NEST_TILE_ALLMAX_d,DAUGHTER_PRC_nprocs(HANDLE)) )

       do i = 1, COMM_CARTESC_NEST_TILE_ALLMAX_d
          if ( i <= COMM_CARTESC_NEST_TILE_ALL ) then
             buffer_LIST(i) = COMM_CARTESC_NEST_TILE_ID(i)
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
       do j = 1, DAUGHTER_PRC_nprocs(HANDLE)
       do i = 1, COMM_CARTESC_NEST_TILE_ALLMAX_d
          COMM_CARTESC_NEST_TILE_LIST_d(i,j) = buffer_ALLLIST(k)
          k = k + 1
       enddo
       enddo

       deallocate( buffer_LIST    )
       deallocate( buffer_ALLLIST )

       ileng = COMM_CARTESC_NEST_TILE_ALLMAX_d*DAUGHTER_PRC_nprocs(HANDLE)
       if ( PRC_IsMaster ) then
          call MPI_ISEND(COMM_CARTESC_NEST_TILE_LIST_d, ileng, MPI_INTEGER, PRC_myrank, tag+2, INTERCOMM_PARENT, ireq, ierr)
          call MPI_WAIT(ireq, istatus, ierr)
       endif

       if ( PRC_IsMaster ) then
          call MPI_ISEND(ONLINE_USE_VELZ, 1, MPI_LOGICAL, PRC_myrank, tag+3, INTERCOMM_PARENT, ireq, ierr)
          call MPI_WAIT(ireq, istatus, ierr)
       endif

       if ( PRC_IsMaster ) then
          call MPI_ISEND(ONLINE_NO_ROTATE, 1, MPI_LOGICAL, PRC_myrank, tag+4, INTERCOMM_PARENT, ireq, ierr)
          call MPI_WAIT(ireq, istatus, ierr)
       endif
       call COMM_bcast(ONLINE_DAUGHTER_NO_ROTATE)

       call COMM_CARTESC_NEST_importgrid_nestdown( HANDLE )

       do i = 1, COMM_CARTESC_NEST_TILE_ALL
          target_rank = COMM_CARTESC_NEST_TILE_LIST_d(i,PRC_myrank+1)
          call MPI_IRECV( call_order(i), 1, MPI_INTEGER, target_rank, tag+5, INTERCOMM_PARENT, ireq, ierr )
          call MPI_WAIT(ireq, istatus, ierr)
       enddo

       call MPI_BARRIER(INTERCOMM_PARENT, ierr)
    else
       write(*,*) 'xxx [COMM_CARTESC_NEST_setup_nestdown] internal error'
       call PRC_MPIstop
    endif

    if( NUM_YP * 16 > max_rq .OR. COMM_CARTESC_NEST_TILE_ALL * 16 > max_rq ) then ! 16 = dyn:5 + qtrc:11
       write(*,*) 'xxx [COMM_CARTESC_NEST_setup_nestdown] internal error (overflow number of ireq)'
       write(*,*) 'xxx NUM_YP x 16        = ', NUM_YP * 16
       write(*,*) 'xxx COMM_CARTESC_NEST_TILE_ALL x 16 = ', COMM_CARTESC_NEST_TILE_ALL * 16
       write(*,*) 'xxx max_rq             = ', max_rq
       call PRC_MPIstop
    endif

    return
  end subroutine COMM_CARTESC_NEST_setup_nestdown

  !-----------------------------------------------------------------------------
  !> Grid Data transfer from parent to daughter: nestdown
  subroutine COMM_CARTESC_NEST_importgrid_nestdown( &
       HANDLE )
    use scale_process, only: &
       PRC_myrank,  &
       PRC_MPIstop
    use scale_atmos_grid_cartesC_real, only: &
       ATMOS_GRID_CARTESC_REAL_LON,  &
       ATMOS_GRID_CARTESC_REAL_LAT,  &
       ATMOS_GRID_CARTESC_REAL_LONU, &
       ATMOS_GRID_CARTESC_REAL_LONV, &
       ATMOS_GRID_CARTESC_REAL_LATU, &
       ATMOS_GRID_CARTESC_REAL_LATV, &
       ATMOS_GRID_CARTESC_REAL_CZ,   &
       ATMOS_GRID_CARTESC_REAL_FZ
    use scale_comm, only: &
       COMM_datatype
    implicit none

    integer, intent(in) :: HANDLE !< id number of nesting relation in this process target

    integer  :: ierr, ileng
    integer  :: istatus(MPI_STATUS_SIZE)
    integer  :: tag, tagbase, target_rank

    integer  :: xloc, yloc
    integer  :: xs, xe
    integer  :: ys, ye

    real(RP) :: max_ref, max_loc

    integer  :: i, k, rq
    !---------------------------------------------------------------------------

    if( .NOT. USE_NESTING ) return

    tagbase = INTERCOMM_ID(HANDLE) * 100
    rq      = 0

    if ( COMM_CARTESC_NEST_Filiation( INTERCOMM_ID(HANDLE) ) > 0 ) then

       !##### parent [send issue] #####

       do i = 1, NUM_YP
          ! send data to multiple daughter processes
          target_rank = COMM_CARTESC_NEST_TILE_LIST_YP(i)

          rq = rq + 1
          ileng = PARENT_IA(HANDLE) * PARENT_JA(HANDLE)
          tag   = tagbase + tag_lon
          call MPI_ISEND(ATMOS_GRID_CARTESC_REAL_LON, ileng, COMM_datatype, target_rank, tag, INTERCOMM_DAUGHTER, ireq_p(rq), ierr)
          call MPI_WAIT(ireq_p(rq), istatus, ierr)

          rq = rq + 1
          ileng = PARENT_IA(HANDLE) * PARENT_JA(HANDLE)
          tag   = tagbase + tag_lat
          call MPI_ISEND(ATMOS_GRID_CARTESC_REAL_LAT, ileng, COMM_datatype, target_rank, tag, INTERCOMM_DAUGHTER, ireq_p(rq), ierr)
          call MPI_WAIT(ireq_p(rq), istatus, ierr)

          rq = rq + 1
          ileng = PARENT_IA(HANDLE) * PARENT_JA(HANDLE)
          tag   = tagbase + tag_lonu
          call MPI_ISEND(ATMOS_GRID_CARTESC_REAL_LONU(1:IA,1:JA), ileng, COMM_datatype, target_rank, tag, INTERCOMM_DAUGHTER, ireq_p(rq), ierr)
          call MPI_WAIT(ireq_p(rq), istatus, ierr)

          rq = rq + 1
          ileng = PARENT_IA(HANDLE) * PARENT_JA(HANDLE)
          tag   = tagbase + tag_latu
          call MPI_ISEND(ATMOS_GRID_CARTESC_REAL_LATU(1:IA,1:JA), ileng, COMM_datatype, target_rank, tag, INTERCOMM_DAUGHTER, ireq_p(rq), ierr)
          call MPI_WAIT(ireq_p(rq), istatus, ierr)

          rq = rq + 1
          ileng = PARENT_IA(HANDLE) * PARENT_JA(HANDLE)
          tag   = tagbase + tag_lonv
          call MPI_ISEND(ATMOS_GRID_CARTESC_REAL_LONV(1:IA,1:JA), ileng, COMM_datatype, target_rank, tag, INTERCOMM_DAUGHTER, ireq_p(rq), ierr)
          call MPI_WAIT(ireq_p(rq), istatus, ierr)

          rq = rq + 1
          ileng = PARENT_IA(HANDLE) * PARENT_JA(HANDLE)
          tag   = tagbase + tag_latv
          call MPI_ISEND(ATMOS_GRID_CARTESC_REAL_LATV(1:IA,1:JA), ileng, COMM_datatype, target_rank, tag, INTERCOMM_DAUGHTER, ireq_p(rq), ierr)
          call MPI_WAIT(ireq_p(rq), istatus, ierr)

          rq = rq + 1
          ileng = PARENT_KA(HANDLE) * PARENT_IA(HANDLE) * PARENT_JA(HANDLE)
          tag   = tagbase + tag_cz
          call MPI_ISEND(ATMOS_GRID_CARTESC_REAL_CZ, ileng, COMM_datatype, target_rank, tag, INTERCOMM_DAUGHTER, ireq_p(rq), ierr)
          call MPI_WAIT(ireq_p(rq), istatus, ierr)

          rq = rq + 1
          ileng = (PARENT_KA(HANDLE)+1) * PARENT_IA(HANDLE) * PARENT_JA(HANDLE)
          tag   = tagbase + tag_fz
          call MPI_ISEND(ATMOS_GRID_CARTESC_REAL_FZ, ileng, COMM_datatype, target_rank, tag, INTERCOMM_DAUGHTER, ireq_p(rq), ierr)
          call MPI_WAIT(ireq_p(rq), istatus, ierr)
       enddo

    elseif( COMM_CARTESC_NEST_Filiation( INTERCOMM_ID(HANDLE) ) < 0 ) then

       !##### child [recv & wait issue] #####

       do i = 1, COMM_CARTESC_NEST_TILE_ALL
          ! receive data from multiple parent tiles
          target_rank = COMM_CARTESC_NEST_TILE_LIST_d(i,PRC_myrank+1)

          xloc = mod( i-1, COMM_CARTESC_NEST_TILE_NUM_X ) + 1
          yloc = int( real(i-1) / real(COMM_CARTESC_NEST_TILE_NUM_X) ) + 1

          xs = PARENT_IMAX(HANDLE) * (xloc-1) + 1
          xe = PARENT_IMAX(HANDLE) *  xloc
          ys = PARENT_JMAX(HANDLE) * (yloc-1) + 1
          ye = PARENT_JMAX(HANDLE) *  yloc

          rq = rq + 1
          ileng = PARENT_IA(HANDLE) * PARENT_JA(HANDLE)
          tag   = tagbase + tag_lon
          call MPI_IRECV(buffer_2D, ileng, COMM_datatype, target_rank, tag, INTERCOMM_PARENT, ireq_d(rq), ierr)
          call MPI_WAIT(ireq_d(rq), istatus, ierr)
          buffer_ref_LON(xs:xe,ys:ye)  = buffer_2D(PRNT_IS(HANDLE):PRNT_IE(HANDLE),PRNT_JS(HANDLE):PRNT_JE(HANDLE))

          rq = rq + 1
          ileng = PARENT_IA(HANDLE) * PARENT_JA(HANDLE)
          tag   = tagbase + tag_lat
          call MPI_IRECV(buffer_2D, ileng, COMM_datatype, target_rank, tag, INTERCOMM_PARENT, ireq_d(rq), ierr)
          call MPI_WAIT(ireq_d(rq), istatus, ierr)
          buffer_ref_LAT(xs:xe,ys:ye)  = buffer_2D(PRNT_IS(HANDLE):PRNT_IE(HANDLE),PRNT_JS(HANDLE):PRNT_JE(HANDLE))

          rq = rq + 1
          ileng = PARENT_IA(HANDLE) * PARENT_JA(HANDLE)
          tag   = tagbase + tag_lonu
          call MPI_IRECV(buffer_2D, ileng, COMM_datatype, target_rank, tag, INTERCOMM_PARENT, ireq_d(rq), ierr)
          call MPI_WAIT(ireq_d(rq), istatus, ierr)
          buffer_ref_LONU(xs:xe,ys:ye)  = buffer_2D(PRNT_IS(HANDLE):PRNT_IE(HANDLE),PRNT_JS(HANDLE):PRNT_JE(HANDLE))

          rq = rq + 1
          ileng = PARENT_IA(HANDLE) * PARENT_JA(HANDLE)
          tag   = tagbase + tag_latu
          call MPI_IRECV(buffer_2D, ileng, COMM_datatype, target_rank, tag, INTERCOMM_PARENT, ireq_d(rq), ierr)
          call MPI_WAIT(ireq_d(rq), istatus, ierr)
          buffer_ref_LATU(xs:xe,ys:ye)  = buffer_2D(PRNT_IS(HANDLE):PRNT_IE(HANDLE),PRNT_JS(HANDLE):PRNT_JE(HANDLE))

          rq = rq + 1
          ileng = PARENT_IA(HANDLE) * PARENT_JA(HANDLE)
          tag   = tagbase + tag_lonv
          call MPI_IRECV(buffer_2D, ileng, COMM_datatype, target_rank, tag, INTERCOMM_PARENT, ireq_d(rq), ierr)
          call MPI_WAIT(ireq_d(rq), istatus, ierr)
          buffer_ref_LONV(xs:xe,ys:ye)  = buffer_2D(PRNT_IS(HANDLE):PRNT_IE(HANDLE),PRNT_JS(HANDLE):PRNT_JE(HANDLE))

          rq = rq + 1
          ileng = PARENT_IA(HANDLE) * PARENT_JA(HANDLE)
          tag   = tagbase + tag_latv
          call MPI_IRECV(buffer_2D, ileng, COMM_datatype, target_rank, tag, INTERCOMM_PARENT, ireq_d(rq), ierr)
          call MPI_WAIT(ireq_d(rq), istatus, ierr)
          buffer_ref_LATV(xs:xe,ys:ye)  = buffer_2D(PRNT_IS(HANDLE):PRNT_IE(HANDLE),PRNT_JS(HANDLE):PRNT_JE(HANDLE))

          rq = rq + 1
          ileng = PARENT_KA(HANDLE) * PARENT_IA(HANDLE) * PARENT_JA(HANDLE)
          tag   = tagbase + tag_cz
          call MPI_IRECV(buffer_3D, ileng, COMM_datatype, target_rank, tag, INTERCOMM_PARENT, ireq_d(rq), ierr)
          call MPI_WAIT(ireq_d(rq), istatus, ierr)
          do k = 1, PARENT_KA(HANDLE)
             buffer_ref_CZ(k,xs:xe,ys:ye)  = buffer_3D(k,PRNT_IS(HANDLE):PRNT_IE(HANDLE),PRNT_JS(HANDLE):PRNT_JE(HANDLE))
          enddo

          rq = rq + 1
          ileng = (PARENT_KA(HANDLE)+1) * PARENT_IA(HANDLE) * PARENT_JA(HANDLE)
          tag   = tagbase + tag_fz
          call MPI_IRECV(buffer_3DF,ileng, COMM_datatype, target_rank, tag, INTERCOMM_PARENT, ireq_d(rq), ierr)
          call MPI_WAIT(ireq_d(rq), istatus, ierr)
          do k = 0, PARENT_KA(HANDLE)
             buffer_ref_FZ(k,xs:xe,ys:ye)  = buffer_3DF(k,PRNT_IS(HANDLE):PRNT_IE(HANDLE),PRNT_JS(HANDLE):PRNT_JE(HANDLE))
          enddo
       enddo

       ! check domain compatibility
       max_ref = maxval( buffer_ref_FZ(:,:,:) )
       max_loc = maxval( ATMOS_GRID_CARTESC_REAL_FZ(KS-1:KE,:,:) ) ! HALO + 1
       if ( max_ref < max_loc ) then
          write(*,*) 'xxx [COMM_CARTESC_NEST_importgrid_nestdown] REQUESTED DOMAIN IS TOO MUCH BROAD'
          write(*,*) 'xxx -- VERTICAL direction over the limit'
          write(*,*) 'xxx -- reference max: ', max_ref
          write(*,*) 'xxx --     local max: ', max_loc
          call PRC_MPIstop
       endif

    else
       write(*,*) 'xxx [COMM_CARTESC_NEST_importgrid_nestdown] internal error'
       call PRC_MPIstop
    endif

    return
  end subroutine COMM_CARTESC_NEST_importgrid_nestdown

  !-----------------------------------------------------------------------------
  !> Boundary data transfer from parent to daughter: nestdown
  subroutine COMM_CARTESC_NEST_nestdown( &
       HANDLE,    &
       BND_QA,    &
       DENS_send, &
       MOMZ_send, &
       MOMX_send, &
       MOMY_send, &
       RHOT_send, &
       QTRC_send, &
       DENS_recv, &
       VELZ_recv, &
       VELX_recv, &
       VELY_recv, &
       POTT_recv, &
       QTRC_recv  )
    use scale_process, only: &
       PRC_MPIstop
    use scale_comm, only: &
       COMM_vars8, &
       COMM_wait
    use scale_atmos_grid_cartesC_metric, only: &
       ROTC => ATMOS_GRID_CARTESC_METRIC_ROTC
    implicit none

    integer,  intent(in)    :: HANDLE !< id number of nesting relation in this process target
    integer,  intent(in)    :: BND_QA !< num of tracer
    real(RP), intent(in)    :: DENS_send(PARENT_KA  (HANDLE),PARENT_IA  (HANDLE),PARENT_JA  (HANDLE))
    real(RP), intent(in)    :: MOMZ_send(PARENT_KA  (HANDLE),PARENT_IA  (HANDLE),PARENT_JA  (HANDLE))
    real(RP), intent(in)    :: MOMX_send(PARENT_KA  (HANDLE),PARENT_IA  (HANDLE),PARENT_JA  (HANDLE))
    real(RP), intent(in)    :: MOMY_send(PARENT_KA  (HANDLE),PARENT_IA  (HANDLE),PARENT_JA  (HANDLE))
    real(RP), intent(in)    :: RHOT_send(PARENT_KA  (HANDLE),PARENT_IA  (HANDLE),PARENT_JA  (HANDLE))
    real(RP), intent(in)    :: QTRC_send(PARENT_KA  (HANDLE),PARENT_IA  (HANDLE),PARENT_JA  (HANDLE),BND_QA)
    real(RP), intent(inout) :: DENS_recv(DAUGHTER_KA(HANDLE),DAUGHTER_IA(HANDLE),DAUGHTER_JA(HANDLE))
    real(RP), intent(inout) :: VELZ_recv(DAUGHTER_KA(HANDLE),DAUGHTER_IA(HANDLE),DAUGHTER_JA(HANDLE))
    real(RP), intent(inout) :: VELX_recv(DAUGHTER_KA(HANDLE),DAUGHTER_IA(HANDLE),DAUGHTER_JA(HANDLE))
    real(RP), intent(inout) :: VELY_recv(DAUGHTER_KA(HANDLE),DAUGHTER_IA(HANDLE),DAUGHTER_JA(HANDLE))
    real(RP), intent(inout) :: POTT_recv(DAUGHTER_KA(HANDLE),DAUGHTER_IA(HANDLE),DAUGHTER_JA(HANDLE))
    real(RP), intent(inout) :: QTRC_recv(DAUGHTER_KA(HANDLE),DAUGHTER_IA(HANDLE),DAUGHTER_JA(HANDLE),BND_QA)

    real(RP) :: WORK1_send(PARENT_KA  (HANDLE),PARENT_IA  (HANDLE),PARENT_JA  (HANDLE))
    real(RP) :: WORK2_send(PARENT_KA  (HANDLE),PARENT_IA  (HANDLE),PARENT_JA  (HANDLE))
    real(RP) :: WORK1_recv(DAUGHTER_KA(HANDLE),DAUGHTER_IA(HANDLE),DAUGHTER_JA(HANDLE))
    real(RP) :: WORK2_recv(DAUGHTER_KA(HANDLE),DAUGHTER_IA(HANDLE),DAUGHTER_JA(HANDLE))
    real(RP) :: U_ll_recv (DAUGHTER_KA(HANDLE),DAUGHTER_IA(HANDLE),DAUGHTER_JA(HANDLE))
    real(RP) :: V_ll_recv (DAUGHTER_KA(HANDLE),DAUGHTER_IA(HANDLE),DAUGHTER_JA(HANDLE))
    real(RP) :: u_on_map, v_on_map

    real(RP) :: dummy(1,1,1)
    integer  :: tagbase, tagcomm
    integer  :: isu_tag, isu_tagf

    integer  :: ierr
    integer  :: i, j, k, iq
    !---------------------------------------------------------------------------

    if( .NOT. USE_NESTING ) return

    if ( BND_QA > I_BNDQA ) then
       write(*,*) 'xxx [COMM_CARTESC_NEST_nestdown] internal error: BND_QA is larger than I_BNDQA'
       call PRC_MPIstop
    elseif( BND_QA > max_bndqa ) then
       write(*,*) 'xxx [COMM_CARTESC_NEST_nestdown] internal error: BND_QA is larger than max_bndqa'
       call PRC_MPIstop
    endif

    tagcomm = INTERCOMM_ID(HANDLE) * order_tag_comm

    if ( COMM_CARTESC_NEST_Filiation( INTERCOMM_ID(HANDLE) ) > 0 ) then

       !##### parent [send issue] #####

       call PROF_rapstart('NEST_total_P', 2)
       call PROF_rapstart('NEST_pack_P', 2)

       nsend = nsend + 1
       if( IO_L ) write(IO_FID_LOG,'(1X,A,I5,A)') "*** CONeP[P] send( ", nsend, " )"

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
       do iq = 1, BND_QA
!OCL XFILL
          org_QTRC(:,:,:,iq) = QTRC_send(:,:,:,iq)
       enddo

       !*** request control
       !--- do not change the calling order below;
       !--- it should be consistent with the order in "COMM_CARTESC_NEST_recvwait_issue"
       rq_ctl_p = 0

       if ( .NOT. ONLINE_DAUGHTER_NO_ROTATE ) then
          ! from staggered point to scalar point
          do j = 1, PARENT_JA(HANDLE)
          do i = 2, PARENT_IA(HANDLE)
          do k = 1, PARENT_KA(HANDLE)
             WORK1_send(k,i,j) = ( org_MOMX(k,i-1,j) + org_MOMX(k,i,j) ) * 0.5_RP
          enddo
          enddo
          enddo

          do j = 1, PARENT_JA(HANDLE)
          do k = 1, PARENT_KA(HANDLE)
             WORK1_send(k,1,j) = org_MOMX(k,1,j)
          enddo
          enddo

          call COMM_vars8( WORK1_send(:,:,:), 1 )

          do j = 2, PARENT_JA(HANDLE)
          do i = 1, PARENT_IA(HANDLE)
          do k = 1, PARENT_KA(HANDLE)
             WORK2_send(k,i,j) = ( org_MOMY(k,i,j-1) + org_MOMY(k,i,j) ) * 0.5_RP
          enddo
          enddo
          enddo

          do i = 1, PARENT_IA(HANDLE)
          do k = 1, PARENT_KA(HANDLE)
             WORK2_send(k,i,1) = org_MOMY(k,i,1)
          enddo
          enddo

          call COMM_vars8( WORK2_send(:,:,:), 2 )

          call COMM_wait ( WORK1_send(:,:,:), 1, .false. )
          call COMM_wait ( WORK2_send(:,:,:), 2, .false. )

          ! rotation from map-projected field to latlon field
          do j = 1, PARENT_JA(HANDLE)
          do i = 1, PARENT_IA(HANDLE)
          do k = 1, PARENT_KA(HANDLE)
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
                                          isu_tag, isu_tagf,       & ! [INOUT]
                                          flag_dens = .true.       ) ! [IN]

       tagbase = tagcomm + tag_momz*order_tag_var
       if ( ONLINE_DAUGHTER_USE_VELZ ) then
          call COMM_CARTESC_NEST_intercomm_nestdown( org_MOMZ(:,:,:),         & ! [IN]
                                             dummy   (:,:,:),         & ! [OUT]
                                             tagbase, I_ZSTG, HANDLE, & ! [IN]
                                             isu_tag, isu_tagf        ) ! [INOUT]
       endif

       tagbase = tagcomm + tag_momx*order_tag_var
       if ( ONLINE_DAUGHTER_NO_ROTATE ) then
          call COMM_CARTESC_NEST_intercomm_nestdown( org_MOMX(:,:,:),         & ! [IN]
                                             dummy   (:,:,:),         & ! [OUT]
                                             tagbase, I_XSTG, HANDLE, & ! [IN]
                                             isu_tag, isu_tagf        ) ! [INOUT]
       else
          call COMM_CARTESC_NEST_intercomm_nestdown( org_U_ll(:,:,:),         & ! [IN]
                                             dummy   (:,:,:),         & ! [OUT]
                                             tagbase, I_SCLR, HANDLE, & ! [IN]
                                             isu_tag, isu_tagf        ) ! [INOUT]
       endif

       tagbase = tagcomm + tag_momy*order_tag_var
       if ( ONLINE_DAUGHTER_NO_ROTATE ) then
          call COMM_CARTESC_NEST_intercomm_nestdown( org_MOMY(:,:,:),         & ! [IN]
                                             dummy   (:,:,:),         & ! [OUT]
                                             tagbase, I_YSTG, HANDLE, & ! [IN]
                                             isu_tag, isu_tagf        ) ! [INOUT]
       else
          call COMM_CARTESC_NEST_intercomm_nestdown( org_V_ll(:,:,:),         & ! [IN]
                                             dummy   (:,:,:),         & ! [OUT]
                                             tagbase, I_SCLR, HANDLE, & ! [IN]
                                             isu_tag, isu_tagf        ) ! [INOUT]
       endif

       tagbase = tagcomm + tag_rhot*order_tag_var
       call COMM_CARTESC_NEST_intercomm_nestdown( org_RHOT(:,:,:),         & ! [IN]
                                          dummy   (:,:,:),         & ! [OUT]
                                          tagbase, I_SCLR, HANDLE, & ! [IN]
                                          isu_tag, isu_tagf        ) ! [INOUT]

       do iq = 1, BND_QA
          tagbase = tagcomm + (tag_qx*10+iq)*order_tag_var
          call COMM_CARTESC_NEST_intercomm_nestdown( org_QTRC(:,:,:,iq),      & ! [IN]
                                             dummy   (:,:,:),         & ! [OUT]
                                             tagbase, I_SCLR, HANDLE, & ! [IN]
                                             isu_tag, isu_tagf        ) ! [INOUT]
       enddo

       rq_tot_p = rq_ctl_p

       call PROF_rapend  ('NEST_pack_P', 2)
       call PROF_rapend  ('NEST_total_P', 2)

    elseif( COMM_CARTESC_NEST_Filiation( INTERCOMM_ID(HANDLE) ) < 0 ) then

       !##### child [wait issue] #####

       call PROF_rapstart('NEST_total_C', 2)
       call PROF_rapstart('NEST_wait_C', 2)

       nwait_d = nwait_d + 1
       !if( IO_L ) write(IO_FID_LOG,'(1X,A,I5,A)') "*** NestIDC [C]: que wait ( ", nwait_d, " )"

       !*** reset issue tag and request control
       !--- do not change the calling order below;
       !--- it should be consistent with the order in "COMM_CARTESC_NEST_recvwait_issue"
       isu_tag  = 0
       isu_tagf = 0

       call COMM_CARTESC_NEST_waitall( rq_tot_d, ireq_d )

       if ( ONLINE_AGGRESSIVE_COMM ) then
          ! nothing to do
       else
          call MPI_BARRIER(INTERCOMM_PARENT, ierr)
       endif

       call PROF_rapend  ('NEST_wait_C', 2)
       call PROF_rapstart('NEST_unpack_C', 2)

       tagbase = tagcomm + tag_dens*order_tag_var
       call COMM_CARTESC_NEST_intercomm_nestdown( dummy     (:,:,:),       & ! [IN]
                                          WORK1_recv(:,:,:),       & ! [OUT]
                                          tagbase, I_SCLR, HANDLE, & ! [IN]
                                          isu_tag, isu_tagf,       & ! [INOUT]
                                          flag_dens = .true.       ) ! [IN]
!OCL XFILL
       do j = 1, DAUGHTER_JA(HANDLE)
       do i = 1, DAUGHTER_IA(HANDLE)
       do k = DATR_KS(HANDLE), DATR_KE(HANDLE)
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
                                             isu_tag, isu_tagf        ) ! [INOUT]
!OCL XFILL
          do j = 1, DAUGHTER_JA(HANDLE)
          do i = 1, DAUGHTER_IA(HANDLE)
          do k = DATR_KS(HANDLE), DATR_KE(HANDLE)-1
             VELZ_recv(k,i,j) = WORK2_recv(k,i,j) / ( WORK1_recv(k,i,j) + WORK1_recv(k+1,i,j) ) * 2.0_RP
          enddo
          enddo
          enddo

          do j = 1, DAUGHTER_JA(HANDLE)
          do i = 1, DAUGHTER_IA(HANDLE)
             VELZ_recv(DATR_KS(HANDLE)-1,i,j) = 0.0_RP
             VELZ_recv(DATR_KE(HANDLE)  ,i,j) = 0.0_RP
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
                                             isu_tag, isu_tagf        ) ! [INOUT]
       else
          ! U_ll_recv receives MOMX/DENS
          call COMM_CARTESC_NEST_intercomm_nestdown( dummy    (:,:,:),        & ! [IN]
                                             U_ll_recv(:,:,:),        & ! [OUT]
                                             tagbase, I_SCLR, HANDLE, & ! [IN]
                                             isu_tag, isu_tagf        ) ! [INOUT]
       endif

       tagbase = tagcomm + tag_momy*order_tag_var
       if ( ONLINE_NO_ROTATE ) then
          ! V_ll_recv receives MOMY
          call COMM_CARTESC_NEST_intercomm_nestdown( dummy     (:,:,:),       & ! [IN]
                                             WORK2_recv(:,:,:),       & ! [OUT]
                                             tagbase, I_YSTG, HANDLE, & ! [IN]
                                             isu_tag, isu_tagf        ) ! [INOUT]
       else
          ! V_ll_recv receives MOMY/DENS
          call COMM_CARTESC_NEST_intercomm_nestdown( dummy    (:,:,:),        & ! [IN]
                                             V_ll_recv(:,:,:),        & ! [OUT]
                                             tagbase, I_SCLR, HANDLE, & ! [IN]
                                             isu_tag, isu_tagf        ) ! [INOUT]
       endif

       if ( ONLINE_NO_ROTATE ) then

!OCL XFILL
          do j = 1, DAUGHTER_JA(HANDLE)
          do i = 1, DAUGHTER_IA(HANDLE)-1
          do k = DATR_KS(HANDLE), DATR_KE(HANDLE)
             VELX_recv(k,i,j) = WORK1_recv(k,i,j) / ( DENS_recv(k,i+1,j) + DENS_recv(k,i,j) ) * 2.0_RP
          enddo
          enddo
          enddo

          i = DAUGHTER_IA(HANDLE)
!OCL XFILL
          do j = 1, DAUGHTER_JA(HANDLE)
          do k = DATR_KS(HANDLE), DATR_KE(HANDLE)
             VELX_recv(k,i,j) = WORK1_recv(k,i,j) / DENS_recv(k,i,j)
          enddo
          enddo

          call COMM_vars8( VELX_recv, 2 )

!OCL XFILL
          do j = 1, DAUGHTER_JA(HANDLE)-1
          do i = 1, DAUGHTER_IA(HANDLE)
          do k = DATR_KS(HANDLE), DATR_KE(HANDLE)
             VELY_recv(k,i,j) = WORK2_recv(k,i,j) / ( DENS_recv(k,i,j+1) + DENS_recv(k,i,j) ) * 2.0_RP
          enddo
          enddo
          enddo

          j = DAUGHTER_JA(HANDLE)
!OCL XFILL
          do i = 1, DAUGHTER_IA(HANDLE)
          do k = DATR_KS(HANDLE), DATR_KE(HANDLE)
             VELY_recv(k,i,j) = WORK2_recv(k,i,j) / DENS_recv(k,i,j)
          enddo
          enddo

          call COMM_vars8( VELY_recv, 3 )

          call COMM_wait ( VELX_recv, 2, .false. )
          call COMM_wait ( VELY_recv, 3, .false. )

       else ! rotate

          ! rotation from latlon field to map-projected field
!OCL XFILL
          do j = 1, DAUGHTER_JA(HANDLE)
          do i = 1, DAUGHTER_IA(HANDLE)
          do k = DATR_KS(HANDLE), DATR_KE(HANDLE)
             WORK1_recv(k,i,j) =  U_ll_recv(k,i,j) * ROTC(i,j,1) + V_ll_recv(k,i,j) * ROTC(i,j,2)
             WORK2_recv(k,i,j) = -U_ll_recv(k,i,j) * ROTC(i,j,2) + V_ll_recv(k,i,j) * ROTC(i,j,1)
          enddo
          enddo
          enddo

          ! from scalar point to staggered point
!OCL XFILL
          do j = 1, DAUGHTER_JA(HANDLE)
          do i = 1, DAUGHTER_IA(HANDLE)-1
          do k = DATR_KS(HANDLE), DATR_KE(HANDLE)
             VELX_recv(k,i,j) = ( WORK1_recv(k,i+1,j) + WORK1_recv(k,i,j) ) * 0.5_RP
          enddo
          enddo
          enddo

          i = DAUGHTER_IA(HANDLE)
!OCL XFILL
          do j = 1, DAUGHTER_JA(HANDLE)
          do k = DATR_KS(HANDLE), DATR_KE(HANDLE)
             VELX_recv(k,i,j) = WORK1_recv(k,i,j)
          enddo
          enddo

          call COMM_vars8( VELX_recv, 2 )

!OCL XFILL
          do j = 1, DAUGHTER_JA(HANDLE)-1
          do i = 1, DAUGHTER_IA(HANDLE)
          do k = DATR_KS(HANDLE), DATR_KE(HANDLE)
             VELY_recv(k,i,j) = ( WORK2_recv(k,i,j+1) + WORK2_recv(k,i,j) ) * 0.5_RP
          enddo
          enddo
          enddo

          j = DAUGHTER_JA(HANDLE)
!OCL XFILL
          do i = 1, DAUGHTER_IA(HANDLE)
          do k = DATR_KS(HANDLE), DATR_KE(HANDLE)
             VELY_recv(k,i,j) = WORK2_recv(k,i,j)
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
                                          isu_tag, isu_tagf        ) ! [INOUT]
!OCL XFILL
       do j = 1, DAUGHTER_JA(HANDLE)
       do i = 1, DAUGHTER_IA(HANDLE)
       do k = DATR_KS(HANDLE), DATR_KE(HANDLE)
          POTT_recv(k,i,j) = WORK1_recv(k,i,j) / DENS_recv(k,i,j)
       enddo
       enddo
       enddo

       do iq = 1, BND_QA
          tagbase = tagcomm + (tag_qx*10+iq)*order_tag_var
          call COMM_CARTESC_NEST_intercomm_nestdown( dummy     (:,:,:),       & ! [IN]
                                             WORK1_recv(:,:,:),       & ! [OUT]
                                             tagbase, I_SCLR, HANDLE, & ! [IN]
                                             isu_tag, isu_tagf        ) ! [INOUT]
!OCL XFILL
          do j = 1, DAUGHTER_JA(HANDLE)
          do i = 1, DAUGHTER_IA(HANDLE)
          do k = DATR_KS(HANDLE), DATR_KE(HANDLE)
             QTRC_recv(k,i,j,iq) = WORK1_recv(k,i,j)
          enddo
          enddo
          enddo
       enddo

       call PROF_rapend  ('NEST_unpack_C', 2)
       call PROF_rapend  ('NEST_total_C', 2)

    else
       write(*,*) 'xxx [COMM_CARTESC_NEST_nestdown] internal error'
       call PRC_MPIstop
    endif

    return
  end subroutine COMM_CARTESC_NEST_nestdown

  !-----------------------------------------------------------------------------
  !> Sub-command for data transfer from parent to daughter: nestdown
  subroutine COMM_CARTESC_NEST_recvwait_issue( &
       HANDLE, &
       BND_QA  )
    use scale_process, only: &
       PRC_MPIstop
    implicit none

    integer, intent(in) :: HANDLE !< id number of nesting relation in this process target
    integer, intent(in) :: BND_QA !< num of tracer in online-nesting

    integer :: isu_tag, isu_tagf
    integer :: tagbase, tagcomm
    integer :: ierr
    integer :: iq
    !---------------------------------------------------------------------------

    if( .NOT. USE_NESTING ) return

    if ( BND_QA > I_BNDQA ) then
       write(*,*) 'xxx [COMM_CARTESC_NEST_recvwait_issue] internal error: about BND_QA'
       call PRC_MPIstop
    endif

    tagcomm = INTERCOMM_ID(HANDLE) * order_tag_comm

    if ( COMM_CARTESC_NEST_Filiation( INTERCOMM_ID(HANDLE) ) > 0 ) then

       !##### parent [wait issue] #####

       call PROF_rapstart('NEST_total_P', 2)
       call PROF_rapstart('NEST_wait_P', 2)

       nwait_p = nwait_p + 1
       !if( IO_L ) write(IO_FID_LOG,'(1X,A,I5,A)') "*** NestIDC [P]: que wait ( ", nwait_p, " )"

       call COMM_CARTESC_NEST_issuer_of_wait( HANDLE )

       if ( ONLINE_AGGRESSIVE_COMM ) then
          ! nothing to do
       else
          call MPI_BARRIER(INTERCOMM_DAUGHTER, ierr)
       endif

       call PROF_rapend  ('NEST_wait_P', 2)
       call PROF_rapend  ('NEST_total_P', 2)

    elseif( COMM_CARTESC_NEST_Filiation( INTERCOMM_ID(HANDLE) ) < 0 ) then

       !##### child [receive issue] #####

       call PROF_rapstart('NEST_total_C', 2)

       nrecv = nrecv + 1
       if( IO_L ) write(IO_FID_LOG,'(1X,A,I5,A)') "*** NestIDC [C]: que recv ( ", nrecv, " )"

       !*** reset issue tag and request control
       !--- do not change the calling order below;
       !--- it should be consistent with the order in "COMM_CARTESC_NEST_nestdown"
       isu_tag  = 0
       isu_tagf = 0
       rq_ctl_d = 0

       tagbase = tagcomm + tag_dens*order_tag_var
       call COMM_CARTESC_NEST_issuer_of_receive( tagbase, I_SCLR, HANDLE, isu_tag, isu_tagf )

       tagbase = tagcomm + tag_momz*order_tag_var
       if ( ONLINE_USE_VELZ ) then
          call COMM_CARTESC_NEST_issuer_of_receive( tagbase, I_ZSTG, HANDLE, isu_tag, isu_tagf )
       endif

       tagbase = tagcomm + tag_momx*order_tag_var
       if ( ONLINE_NO_ROTATE ) then
          call COMM_CARTESC_NEST_issuer_of_receive( tagbase, I_XSTG, HANDLE, isu_tag, isu_tagf )
       else
          call COMM_CARTESC_NEST_issuer_of_receive( tagbase, I_SCLR, HANDLE, isu_tag, isu_tagf )
       endif

       tagbase = tagcomm + tag_momy*order_tag_var
       if ( ONLINE_NO_ROTATE ) then
          call COMM_CARTESC_NEST_issuer_of_receive( tagbase, I_YSTG, HANDLE, isu_tag, isu_tagf )
       else
          call COMM_CARTESC_NEST_issuer_of_receive( tagbase, I_SCLR, HANDLE, isu_tag, isu_tagf )
       endif

       tagbase = tagcomm + tag_rhot*order_tag_var
       call COMM_CARTESC_NEST_issuer_of_receive( tagbase, I_SCLR, HANDLE, isu_tag, isu_tagf )

       do iq = 1, BND_QA
          tagbase = tagcomm + (tag_qx*10+iq)*order_tag_var
          call COMM_CARTESC_NEST_issuer_of_receive( tagbase, I_SCLR, HANDLE, isu_tag, isu_tagf )
       enddo

       rq_tot_d = rq_ctl_d

       call PROF_rapend('NEST_total_C', 2)

    else
       write(*,*) 'xxx [COMM_CARTESC_NEST_recvwait_issue] internal error'
       call PRC_MPIstop
    endif

    return
  end subroutine COMM_CARTESC_NEST_recvwait_issue

  !-----------------------------------------------------------------------------
  !> Sub-command for data transfer from parent to daughter: nestdown
  subroutine COMM_CARTESC_NEST_recv_cancel( &
       HANDLE )
    use scale_process, only: &
       PRC_MPIstop
    implicit none

    integer, intent(in) :: HANDLE !< id number of nesting relation in this process target

    !logical :: flag
    !integer :: istatus(MPI_STATUS_SIZE)

    integer :: rq
    integer :: ierr
    !---------------------------------------------------------------------------

    if( .NOT. USE_NESTING ) return

    if ( COMM_CARTESC_NEST_Filiation( INTERCOMM_ID(HANDLE) ) > 0 ) then

       !##### parent #####
       ! Nothing to do

    elseif( COMM_CARTESC_NEST_Filiation( INTERCOMM_ID(HANDLE) ) < 0 ) then

       !##### child #####

       if( IO_L ) write(IO_FID_LOG,'(1X,A,I5,A)') "*** NestIDC [C]: CANCEL recv ( ", nrecv, " )"

       do rq = 1, rq_tot_d
          if ( ireq_d(rq) /= MPI_REQUEST_NULL ) then

             call MPI_CANCEL(ireq_d(rq), ierr)

!              call MPI_TEST_CANCELLED(istatus, flag, ierr)
!              if ( .NOT. flag ) then
!                 write(*,*) 'xxx receive actions do not cancelled, req = ', rq
!              endif
          endif
       enddo

    else
       write(*,*) 'xxx [COMM_CARTESC_NEST_recv_cancel] internal error'
       call PRC_MPIstop
    endif

    return
  end subroutine COMM_CARTESC_NEST_recv_cancel

  !-----------------------------------------------------------------------------
  !> Inter-communication from parent to daughter: nestdown
  subroutine COMM_CARTESC_NEST_intercomm_nestdown_3D( &
       pvar,     &
       dvar,     &
       tagbase,  &
       id_stag,  &
       HANDLE,   &
       isu_tag,  &
       isu_tagf, &
       flag_dens )
    use scale_process, only: &
       PRC_MPIstop
    use scale_comm, only: &
       COMM_datatype
    use scale_interp, only: &
       INTRP_interp3d
    implicit none

    real(RP), intent(in)    :: pvar(:,:,:) !< variable from parent domain (PARENT_KA,PARENT_IA,PARENT_JA / 1,1,1)
    real(RP), intent(out)   :: dvar(:,:,:) !< variable to daughter domain (1,1,1 / MY_KA,MY_IA,MY_JA)
    integer,  intent(in)    :: tagbase     !< communication tag of the variable
    integer,  intent(in)    :: id_stag     !< id of staggered grid option
    integer,  intent(in)    :: HANDLE      !< id number of nesting relation in this process target
    integer,  intent(inout) :: isu_tag     !< tag for receive buffer
    integer,  intent(inout) :: isu_tagf    !< tag for receive buffer

    logical , intent(in), optional :: flag_dens !< flag of logarithmic interpolation for density

    integer :: ileng, tag, target_rank

    integer :: xloc, yloc
    integer :: gxs, gxe, gys, gye ! for large  domain
    integer :: pxs, pxe, pys, pye ! for parent domain
    integer :: zs, ze

    integer :: ig, rq, yp
    logical :: no_zstag    = .true.
    logical :: logarithmic = .false.

    integer :: ierr
    integer :: i, j
    !---------------------------------------------------------------------------

    if( .NOT. USE_NESTING ) return

    logarithmic = .false.
    if ( present(flag_dens) ) then
       if( flag_dens ) logarithmic = .true.
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

    if ( no_zstag ) then
       ileng = (PARENT_KA(HANDLE)  ) * PARENT_IA(HANDLE) * PARENT_JA(HANDLE)
    else
       ileng = (PARENT_KA(HANDLE)+1) * PARENT_IA(HANDLE) * PARENT_JA(HANDLE)
    endif

    if ( COMM_CARTESC_NEST_Filiation( INTERCOMM_ID(HANDLE) ) > 0 ) then

       !##### parent #####
       rq = rq_ctl_p

       do yp = 1, NUM_YP
          rq = rq + 1

          ! send data to multiple daughter processes
          target_rank = COMM_CARTESC_NEST_TILE_LIST_YP(yp)
          tag         = tagbase + yp

          call MPI_ISEND( pvar,               &
                          ileng,              &
                          COMM_datatype,      &
                          target_rank,        &
                          tag,                &
                          INTERCOMM_DAUGHTER, &
                          ireq_p(rq),         &
                          ierr                )

          dvar(:,:,:) = -1.0_RP  ! input as a dummy value
       enddo

       rq_ctl_p = rq

    elseif( COMM_CARTESC_NEST_Filiation( INTERCOMM_ID(HANDLE) ) < 0 ) then

       !##### child #####
       rq = rq_ctl_d

       do yp = 1, COMM_CARTESC_NEST_TILE_ALL
          rq = rq + 1

          xloc = mod( yp-1, COMM_CARTESC_NEST_TILE_NUM_X ) + 1
          yloc = int( real(yp-1) / real(COMM_CARTESC_NEST_TILE_NUM_X) ) + 1

          gxs = PARENT_IMAX(HANDLE) * (xloc-1) + 1
          gxe = PARENT_IMAX(HANDLE) * xloc
          gys = PARENT_JMAX(HANDLE) * (yloc-1) + 1
          gye = PARENT_JMAX(HANDLE) * yloc

          pxs = PRNT_IS(HANDLE)
          pxe = PRNT_IE(HANDLE)
          pys = PRNT_JS(HANDLE)
          pye = PRNT_JE(HANDLE)

          if ( no_zstag ) then
             isu_tag = isu_tag + 1

             zs = 1
             ze = PARENT_KA(HANDLE)
!OCL XFILL
             buffer_ref_3D(zs:ze,gxs:gxe,gys:gye) = recvbuf_3D(zs:ze,pxs:pxe,pys:pye,isu_tag)
          else
             isu_tagf = isu_tagf + 1

             zs = 0
             ze = PARENT_KA(HANDLE)
!OCL XFILL
             buffer_ref_3DF(zs:ze,gxs:gxe,gys:gye) = recvbuf_3DF(zs:ze,pxs:pxe,pys:pye,isu_tagf)
          endif

          if ( isu_tag > max_isu .OR. isu_tagf > max_isuf ) then
             write(*,*) 'xxx [COMM_CARTESC_NEST_intercomm_nestdown_3D] Exceeded maximum issue'
             write(*,*) 'xxx isu_tag  = ', isu_tag
             write(*,*) 'xxx isu_tagf = ', isu_tagf
             call PRC_MPIstop
          endif

       enddo

       rq_ctl_d = rq

       if ( no_zstag ) then
          call INTRP_interp3d( itp_nh,                       & ! [IN]
                               TILEAL_KA  (HANDLE),          & ! [IN]
                               TILEAL_IA  (HANDLE),          & ! [IN]
                               TILEAL_JA  (HANDLE),          & ! [IN]
                               DAUGHTER_KA(HANDLE),          & ! [IN]
                               DATR_KS    (HANDLE),          & ! [IN]
                               DATR_KE    (HANDLE),          & ! [IN]
                               DAUGHTER_IA(HANDLE),          & ! [IN]
                               DAUGHTER_JA(HANDLE),          & ! [IN]
                               igrd          (    :,:,:,ig), & ! [IN]
                               jgrd          (    :,:,:,ig), & ! [IN]
                               hfact         (    :,:,:,ig), & ! [IN]
                               kgrd          (:,:,:,:,:,ig), & ! [IN]
                               vfact         (:,:,:,:,:,ig), & ! [IN]
                               buffer_ref_3D (:,:,:),        & ! [INOUT]
                               dvar          (:,:,:),        & ! [OUT]
                               logwgt = logarithmic          ) ! [IN]
       else
          call INTRP_interp3d( itp_nh,                       & ! [IN]
                               TILEAL_KA  (HANDLE)+1,        & ! [IN]
                               TILEAL_IA  (HANDLE),          & ! [IN]
                               TILEAL_JA  (HANDLE),          & ! [IN]
                               DAUGHTER_KA(HANDLE),          & ! [IN]
                               DATR_KS    (HANDLE),          & ! [IN]
                               DATR_KE    (HANDLE),          & ! [IN]
                               DAUGHTER_IA(HANDLE),          & ! [IN]
                               DAUGHTER_JA(HANDLE),          & ! [IN]
                               igrd          (    :,:,:,ig), & ! [IN]
                               jgrd          (    :,:,:,ig), & ! [IN]
                               hfact         (    :,:,:,ig), & ! [IN]
                               kgrd          (:,:,:,:,:,ig), & ! [IN]
                               vfact         (:,:,:,:,:,ig), & ! [IN]
                               buffer_ref_3DF(:,:,:),        & ! [INOUT]
                               dvar          (:,:,:),        & ! [OUT]
                               logwgt = logarithmic          ) ! [IN]
       endif

       do j = 1, DAUGHTER_JA(HANDLE)
       do i = 1, DAUGHTER_IA(HANDLE)
          dvar(                1:DATR_KS    (HANDLE)-1,i,j) = 0.0_RP
          dvar(DATR_KE(HANDLE)+1:DAUGHTER_KA(HANDLE)  ,i,j) = 0.0_RP
       enddo
       enddo

    else
       write(*,*) 'xxx [COMM_CARTESC_NEST_intercomm_nestdown_3D] internal error'
       call PRC_MPIstop
    endif

    return
  end subroutine COMM_CARTESC_NEST_intercomm_nestdown_3D

  !-----------------------------------------------------------------------------
  !> [substance of issuer] Inter-communication from parent to daughter: nestdown
  subroutine COMM_CARTESC_NEST_issuer_of_receive_3D( &
       tagbase, &
       id_stag, &
       HANDLE,  &
       isu_tag, &
       isu_tagf )
    use scale_process, only: &
       PRC_myrank,  &
       PRC_MPIstop
    use scale_comm, only: &
       COMM_datatype
    implicit none

    integer, intent(in)    :: tagbase  !< communication tag of the variable
    integer, intent(in)    :: id_stag  !< id of staggered grid option
    integer, intent(in)    :: HANDLE   !< id number of nesting relation in this process target
    integer, intent(inout) :: isu_tag  !< tag for receive buffer
    integer, intent(inout) :: isu_tagf !< tag for receive buffer

    integer :: ierr, ileng
    integer :: tag, target_rank

    integer :: ig, rq, yp
    logical :: no_zstag    = .true.
    !---------------------------------------------------------------------------

    if( .NOT. USE_NESTING ) return

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

    if ( no_zstag ) then
       ileng = (PARENT_KA(HANDLE)  ) * PARENT_IA(HANDLE) * PARENT_JA(HANDLE)
    else
       ileng = (PARENT_KA(HANDLE)+1) * PARENT_IA(HANDLE) * PARENT_JA(HANDLE)
    endif

    if ( COMM_CARTESC_NEST_Filiation( INTERCOMM_ID(HANDLE) ) > 0 ) then

       !##### parent #####
       ! nothing to do

    elseif( COMM_CARTESC_NEST_Filiation( INTERCOMM_ID(HANDLE) ) < 0 ) then

       !##### child #####
       rq = rq_ctl_d

       do yp = 1, COMM_CARTESC_NEST_TILE_ALL
          rq = rq + 1

          target_rank = COMM_CARTESC_NEST_TILE_LIST_d(yp,PRC_myrank+1)
          tag         = tagbase + call_order(yp)

          if ( no_zstag ) then
             isu_tag = isu_tag + 1

             recvbuf_3D(:,:,:,isu_tag) = 0.0_RP

             call MPI_IRECV( recvbuf_3D(:,:,:,isu_tag), &
                             ileng,                     &
                             COMM_datatype,             &
                             target_rank,               &
                             tag,                       &
                             INTERCOMM_PARENT,          &
                             ireq_d(rq),                &
                             ierr                       )
          else
             isu_tagf = isu_tagf + 1

             recvbuf_3DF(:,:,:,isu_tagf) = 0.0_RP

             call MPI_IRECV( recvbuf_3DF(:,:,:,isu_tagf), &
                             ileng,                       &
                             COMM_datatype,               &
                             target_rank,                 &
                             tag,                         &
                             INTERCOMM_PARENT,            &
                             ireq_d(rq),                  &
                             ierr                         )
          endif

       enddo

       if ( isu_tag > max_isu .OR. isu_tagf > max_isuf ) then
          write(*,*) 'xxx [COMM_CARTESC_NEST_issuer_of_receive_3D] Exceeded maximum issue'
          write(*,*) 'xxx isu_tag  = ', isu_tag
          write(*,*) 'xxx isu_tagf = ', isu_tagf
          call PRC_MPIstop
       endif

       rq_ctl_d = rq

    else
       write(*,*) 'xxx [COMM_CARTESC_NEST_issuer_of_receive_3D] internal error'
       call PRC_MPIstop
    endif

    return
  end subroutine COMM_CARTESC_NEST_issuer_of_receive_3D

  !-----------------------------------------------------------------------------
  !> [substance of issuer] Inter-communication from parent to daughter: nestdown
  subroutine COMM_CARTESC_NEST_issuer_of_wait_3D( &
       HANDLE )
    use scale_process, only: &
       PRC_MPIstop
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
       write(*,*) 'xxx [COMM_CARTESC_NEST_issuer_of_wait_3D] internal error'
       call PRC_MPIstop
    endif

    return
  end subroutine COMM_CARTESC_NEST_issuer_of_wait_3D

  !-----------------------------------------------------------------------------
  !> [substance of comm_wait] Inter-communication
  subroutine COMM_CARTESC_NEST_waitall( &
       req_count, &
       ireq       )
    use scale_process, only: &
       PRC_MPIstop
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
!          if( IO_L ) write(IO_FID_LOG,'(1x,A)') '*** ERROR: over the limit of waiting time [NESTCOM]'
!          write(*,'(1x,A)') '*** ERROR: over the limit of waiting time [NESTCOM]'
!          call PRC_MPIstop
!       endif
!    enddo

    return
  end subroutine COMM_CARTESC_NEST_waitall

  !-----------------------------------------------------------------------------
  !> [check communication status] Inter-communication
  subroutine COMM_CARTESC_NEST_test( &
       HANDLE )
    use scale_process, only: &
       PRC_MPIstop
    implicit none

    integer, intent(in) :: HANDLE !< id number of nesting relation in this process target

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

    elseif( COMM_CARTESC_NEST_Filiation( INTERCOMM_ID(HANDLE) ) < 0 ) then

       !##### child #####

       call PROF_rapstart('NEST_test_C', 2)
       if ( rq_ctl_d > 0 ) call MPI_TEST(ireq_d(1), flag, istatus, ierr)
       call PROF_rapend('NEST_test_C', 2)

    else
       write(*,*) 'xxx [COMM_CARTESC_NEST_test] error'
       call PRC_MPIstop
    endif

    return
  end subroutine COMM_CARTESC_NEST_test

  !-----------------------------------------------------------------------------
  !> [finalize: disconnect] Inter-communication
  subroutine COMM_CARTESC_NEST_disconnect
    use scale_process, only: &
       PRC_GLOBAL_COMM_WORLD
    implicit none

    integer :: ierr
    !---------------------------------------------------------------------------

    if( .NOT. USE_NESTING ) return

    if( IO_L ) write(IO_FID_LOG,'(1x,A)') '*** Waiting finish of whole processes'
    call MPI_BARRIER(PRC_GLOBAL_COMM_WORLD, ierr)

    if ( ONLINE_IAM_PARENT ) then
       !if( IO_L ) write(IO_FID_LOG,'(1x,A)') '*** Waiting finish of whole processes as a parent'
       !call MPI_BARRIER(INTERCOMM_DAUGHTER, ierr)
       call MPI_COMM_FREE(INTERCOMM_DAUGHTER, ierr)
       if( IO_L ) write(IO_FID_LOG,'(1x,A)') '*** Disconnected communication with child'
    endif

    if ( ONLINE_IAM_DAUGHTER ) then
       !if( IO_L ) write(IO_FID_LOG,'(1x,A)') '*** Waiting finish of whole processes as a child'
       !call MPI_BARRIER(INTERCOMM_PARENT, ierr)
       call MPI_COMM_FREE(INTERCOMM_PARENT, ierr)
       if( IO_L ) write(IO_FID_LOG,'(1x,A)') '*** Disconnected communication with parent'
    endif

    return
  end subroutine COMM_CARTESC_NEST_disconnect

end module scale_comm_cartesC_nest
