!-------------------------------------------------------------------------------
!> module GRID (nesting system)
!!
!! @par Description
!!          Grid module for nesting system
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2014-07-28 (R.Yoshida)  [new]
!! @li      2014-09-05 (R.Yoshida)  [add] online communication system
!!
!<
!-------------------------------------------------------------------------------
module scale_grid_nest
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use mpi              ![external]
  use scale_precision
  use scale_stdio
  use scale_prof
  use scale_debug
  use scale_grid_index
  use scale_index
  use scale_tracer
  use scale_const, only: &
     r_in_m => CONST_RADIUS
  use scale_interpolation_nest
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: NEST_setup
  public :: NEST_domain_relate
  public :: NEST_COMM_nestdown
  public :: NEST_COMM_recvwait_issue
  public :: NEST_COMM_recv_cancel
  public :: NEST_COMM_test
  public :: NEST_COMM_disconnect

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  integer,  public              :: INTERCOMM_PARENT     ! inter-communicator to parent
  integer,  public              :: INTERCOMM_DAUGHTER   ! inter-communicator to daughter

  integer,  public              :: NEST_Filiation(10)   !< index of parent-daughter relation (p>0, d<0)
  integer,  public              :: HANDLING_NUM         !< handing number of nesting relation
  integer,  public              :: NEST_TILE_NUM_X      !< parent tile number in x-direction
  integer,  public              :: NEST_TILE_NUM_Y      !< parent tile number in y-direction
  integer,  public, allocatable :: NEST_TILE_ID(:)      !< parent tile real id

  integer,  public              :: PARENT_KMAX(2)       !< parent max number in z-direction
  integer,  public              :: PARENT_IMAX(2)       !< parent max number in x-direction
  integer,  public              :: PARENT_JMAX(2)       !< parent max number in y-direction
  integer,  public              :: PARENT_KA(2)         !< parent max number in z-direction (with halo)
  integer,  public              :: PARENT_IA(2)         !< parent max number in x-direction (with halo)
  integer,  public              :: PARENT_JA(2)         !< parent max number in y-direction (with halo)
  integer,  public              :: PARENT_LKMAX(2)      !< parent max number in lz-direction
  real(DP), public              :: PARENT_DTSEC(2)      !< parent DT [sec]
  integer,  public              :: PARENT_NSTEP(2)      !< parent step [number]

  integer,  public              :: DAUGHTER_KMAX(2)     !< daughter max number in z-direction
  integer,  public              :: DAUGHTER_IMAX(2)     !< daughter max number in x-direction
  integer,  public              :: DAUGHTER_JMAX(2)     !< daughter max number in y-direction
  integer,  public              :: DAUGHTER_KA(2)       !< daughter max number in z-direction (with halo)
  integer,  public              :: DAUGHTER_IA(2)       !< daughter max number in x-direction (with halo)
  integer,  public              :: DAUGHTER_JA(2)       !< daughter max number in y-direction (with halo)
  integer,  public              :: DAUGHTER_LKMAX(2)    !< daughter max number in lz-direction
  real(DP), public              :: DAUGHTER_DTSEC(2)    !< daughter DT [sec]
  integer,  public              :: DAUGHTER_NSTEP(2)    !< daughter steps [number]

  integer,  public              :: PRNT_KS(2)           !< start index in z-direction in parent
  integer,  public              :: PRNT_KE(2)           !< end index   in z-direction in parent
  integer,  public              :: PRNT_IS(2)           !< start index in x-direction in parent
  integer,  public              :: PRNT_IE(2)           !< end index   in x-direction in parent
  integer,  public              :: PRNT_JS(2)           !< start index in y-direction in parent
  integer,  public              :: PRNT_JE(2)           !< end index   in y-direction in parent

  integer,  public              :: DATR_KS(2)           !< start index in z-direction in daughter
  integer,  public              :: DATR_KE(2)           !< end index   in z-direction in daughter
  integer,  public              :: DATR_IS(2)           !< start index in x-direction in daughter
  integer,  public              :: DATR_IE(2)           !< end index   in x-direction in daughter
  integer,  public              :: DATR_JS(2)           !< start index in y-direction in daughter
  integer,  public              :: DATR_JE(2)           !< end index   in y-direction in daughter

  integer,  public              :: TILEAL_KA(2)         !< cells of all tiles in z-direction
  integer,  public              :: TILEAL_IA(2)         !< cells of all tiles in x-direction
  integer,  public              :: TILEAL_JA(2)         !< cells of all tiles in y-direction

  integer,  public              :: NEST_BND_QA = 1      !< number of tracer treated in nesting system
  integer,  public              :: NEST_INTERP_LEVEL = 3 !< horizontal interpolation level

  logical,  public              :: USE_NESTING          = .false.
  logical,  public              :: OFFLINE              = .true.
  logical,  public              :: ONLINE_IAM_PARENT    = .false.   !< a flag to say "I am a parent"
  logical,  public              :: ONLINE_IAM_DAUGHTER  = .false.   !< a flag to say "I am a daughter"
  integer,  public              :: ONLINE_DOMAIN_NUM    = 1
  logical,  public              :: ONLINE_USE_VELZ      = .false.
  logical,  public              :: ONLINE_NO_ROTATE     = .false.
  logical,  public              :: ONLINE_BOUNDARY_USE_QHYD = .false.

  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  private :: NEST_COMM_parentsize
  private :: NEST_COMM_catalogue
  private :: NEST_COMM_ping
  private :: NEST_COMM_setup_nestdown
  private :: NEST_COMM_importgrid_nestdown
  private :: NEST_COMM_intercomm_nestdown
  private :: NEST_COMM_issuer_of_receive
  private :: NEST_COMM_issuer_of_wait
!  private :: NEST_latlonz_interporation_fact

  interface NEST_COMM_intercomm_nestdown
     module procedure NEST_COMM_intercomm_nestdown_3D
  end interface NEST_COMM_intercomm_nestdown

  interface NEST_COMM_issuer_of_receive
     module procedure NEST_COMM_issuer_of_receive_3D
  end interface NEST_COMM_issuer_of_receive

  interface NEST_COMM_issuer_of_wait
     module procedure NEST_COMM_issuer_of_wait_3D
  end interface NEST_COMM_issuer_of_wait

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  real(RP), private, allocatable :: latlon_catalog(:,:,:)    !< parent latlon catalog [rad]
  real(RP), private              :: corner_loc(4,2)          !< local corner location [rad]

  integer, private               :: PARENT_PRC_NUM_X(2)      !< MPI processes in x-direction in parent
  integer, private               :: PARENT_PRC_NUM_Y(2)      !< MPI processes in y-direction in parent
  integer, private               :: PARENT_PRC_nmax(2)       !< MPI total processes in parent

  integer, private               :: DAUGHTER_PRC_NUM_X(2)    !< MPI processes in x-direction in daughter
  integer, private               :: DAUGHTER_PRC_NUM_Y(2)    !< MPI processes in y-direction in daughter
  integer, private               :: DAUGHTER_PRC_nmax(2)     !< MPI total processes in daughter

  integer, private               :: NEST_TILE_ALL            !< NUM of TILEs in the local node
  integer, private               :: NEST_TILE_ALLMAX_p       !< MAXNUM of TILEs among whole processes for parent
  integer, private               :: NEST_TILE_ALLMAX_d       !< MAXNUM of TILEs among whole processes for daughter
  integer, private, allocatable  :: NEST_TILE_LIST_p(:,:)    !< relationship list in whole system for parent
  integer, private, allocatable  :: NEST_TILE_LIST_d(:,:)    !< relationship list in whole system for daughter
  integer, private, allocatable  :: NEST_TILE_LIST_YP(:)     !< yellow-page of daughter targets for parent
  integer, private               :: NUM_YP                   !< page number of yellow-page

  integer, private               :: OFFLINE_PARENT_PRC_NUM_X !< MPI processes in x-direction in parent [for namelist]
  integer, private               :: OFFLINE_PARENT_PRC_NUM_Y !< MPI processes in y-direction in parent [for namelist]
  integer, private               :: OFFLINE_PARENT_KMAX      !< parent max number in z-direction [for namelist]
  integer, private               :: OFFLINE_PARENT_IMAX      !< parent max number in x-direction [for namelist]
  integer, private               :: OFFLINE_PARENT_JMAX      !< parent max number in y-direction [for namelist]
  integer, private               :: OFFLINE_PARENT_LKMAX     !< parent max number in lz-direction [for namelist]
  integer(8), private            :: ONLINE_WAIT_LIMIT        !< limit times of waiting loop in "NEST_COMM_waitall"
  logical, private               :: ONLINE_DAUGHTER_USE_VELZ
  logical, private               :: ONLINE_DAUGHTER_NO_ROTATE
  logical, private               :: ONLINE_AGGRESSIVE_COMM

  integer, parameter :: I_LON    = 1
  integer, parameter :: I_LAT    = 2

  integer, parameter :: I_NW     = 1
  integer, parameter :: I_NE     = 2
  integer, parameter :: I_SW     = 3
  integer, parameter :: I_SE     = 4
  integer, parameter :: I_BNDQA  = 20                      !< tentative approach (prefixed allocate size)

  integer, parameter :: I_SCLR   = 1                       !< interpolation kinds of grid point (scalar)
  integer, parameter :: I_ZSTG   = 2                       !< interpolation kinds of grid point (z-axis staggered)
  integer, parameter :: I_XSTG   = 3                       !< interpolation kinds of grid point (x-axis staggered)
  integer, parameter :: I_YSTG   = 4                       !< interpolation kinds of grid point (y-axis staggered)

  integer, parameter :: itp_ng   = 4                       !< # of interpolation kinds of grid point
  integer, private   :: itp_nh   = 3                       !< # of interpolation kinds of horizontal direction
  integer, private   :: itp_nv   = 2                       !< # of interpolation kinds of vertical direction

  integer, parameter :: tag_lon  = 1
  integer, parameter :: tag_lat  = 2
  integer, parameter :: tag_lonx = 3
  integer, parameter :: tag_latx = 4
  integer, parameter :: tag_lony = 5
  integer, parameter :: tag_laty = 6
  integer, parameter :: tag_cz   = 7
  integer, parameter :: tag_fz   = 8

  integer, parameter :: tag_dens = 1
  integer, parameter :: tag_momz = 2
  integer, parameter :: tag_momx = 3
  integer, parameter :: tag_momy = 4
  integer, parameter :: tag_rhot = 5
  integer, parameter :: tag_qx   = 6

  integer, parameter :: order_tag_comm = 100000
  integer, parameter :: order_tag_var  = 1000
  ! intercomm tag id:  IC | VAR |  YP
  ! (total: 6columns)   X   X X   X X X

!  real(RP), private, parameter :: large_number_one   = 9.999E+15_RP
!  real(RP), private, parameter :: large_number_two   = 8.888E+15_RP
!  real(RP), private, parameter :: large_number_three = 7.777E+15_RP
  integer,  private            :: interp_search_divnum

  integer, private   :: INTERCOMM_ID(2)

  integer, private, parameter :: max_isu   = 100             ! maximum number of receive/wait issue
  integer, private, parameter :: max_isuf  = 20              ! maximum number of receive/wait issue (z-stag)
  integer, private, parameter :: max_bndqa = 12              ! maximum number of QA in boundary: tentative approach
  integer, private            :: max_rq    = 1000            ! maximum number of req: tentative approach
  integer, private            :: rq_ctl_p                    ! for control request id (counting)
  integer, private            :: rq_ctl_d                    ! for control request id (counting)
  integer, private            :: rq_tot_p                    ! for control request id (total number)
  integer, private            :: rq_tot_d                    ! for control request id (total number)
  integer, private, allocatable :: ireq_p(:)                 ! buffer of request-id for parent
  integer, private, allocatable :: ireq_d(:)                 ! buffer of request-id for daughter
  integer, private, allocatable :: call_order(:)             ! calling order from parent

  real(RP), private, allocatable :: buffer_2D (:,:)          ! buffer of communicator: 2D (with HALO)
  real(RP), private, allocatable :: buffer_3D (:,:,:)        ! buffer of communicator: 3D (with HALO)
  real(RP), private, allocatable :: buffer_3DF(:,:,:)        ! buffer of communicator: 3D-Kface (with HALO)
  real(RP), private, allocatable :: recvbuf_3D (:,:,:,:)     ! buffer of receiver: 3D (with HALO)
  real(RP), private, allocatable :: recvbuf_3DF(:,:,:,:)     ! buffer of receiver: 3D-Kface (with HALO)

  real(RP), private, allocatable :: buffer_ref_LON (:,:)     ! buffer of communicator: DENS (with HALO)
  real(RP), private, allocatable :: buffer_ref_LONX(:,:)     ! buffer of communicator: VELZ (with HALO)
  real(RP), private, allocatable :: buffer_ref_LONY(:,:)     ! buffer of communicator: VELX (with HALO)
  real(RP), private, allocatable :: buffer_ref_LAT (:,:)     ! buffer of communicator: VELY (with HALO)
  real(RP), private, allocatable :: buffer_ref_LATX(:,:)     ! buffer of communicator: POTT (with HALO)
  real(RP), private, allocatable :: buffer_ref_LATY(:,:)     ! buffer of communicator: POTT (with HALO)
  real(RP), private, allocatable :: buffer_ref_CZ  (:,:,:)   ! buffer of communicator: VELY (with HALO)
  real(RP), private, allocatable :: buffer_ref_FZ  (:,:,:)   ! buffer of communicator: VELY (with HALO)

  !real(RP), private, allocatable :: buffer_ref_2D (:,:)      ! buffer of communicator: 2D data (with HALO)
  real(RP), private, allocatable :: buffer_ref_3D (:,:,:)    ! buffer of communicator: 3D data (with HALO)
  real(RP), private, allocatable :: buffer_ref_3DF(:,:,:)    ! buffer of communicator: 3D at z-Face (with HALO)
  real(RP), private, allocatable :: u_llp(:,:,:)
  real(RP), private, allocatable :: v_llp(:,:,:)

  real(RP), private, allocatable :: org_DENS(:,:,:)          ! buffer to keep values at that time in parent
  real(RP), private, allocatable :: org_MOMZ(:,:,:)          ! buffer to keep values at that time in parent
  real(RP), private, allocatable :: org_MOMX(:,:,:)          ! buffer to keep values at that time in parent
  real(RP), private, allocatable :: org_MOMY(:,:,:)          ! buffer to keep values at that time in parent
  real(RP), private, allocatable :: org_RHOT(:,:,:)          ! buffer to keep values at that time in parent
  real(RP), private, allocatable :: org_QTRC(:,:,:,:)        ! buffer to keep values at that time in parent

  real(RP), private, allocatable :: hfact(:,:,:,:)           ! interpolation factor for horizontal direction
  real(RP), private, allocatable :: vfact(:,:,:,:,:,:)       ! interpolation factor for vertical direction
  integer,  private, allocatable :: kgrd (:,:,:,:,:,:)       ! interpolation target grids in z-axis
  integer,  private, allocatable :: igrd (:,:,:,:)           ! interpolation target grids in x-axis
  integer,  private, allocatable :: jgrd (:,:,:,:)           ! interpolation target grids in y-axis

  integer(8), private :: nwait_p, nwait_d, nrecv, nsend
  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine NEST_setup ( &
      icomm_parent,  &
      icomm_child,   &
      flag_parent,   &
      flag_child     )
    use scale_const, only: &
       D2R => CONST_D2R
    use scale_process, only: &
       PRC_nmax,          &
       PRC_master,        &
       PRC_myrank,        &
       PRC_MPIstop,       &
       PRC_HAS_W,         &
       PRC_HAS_E,         &
       PRC_HAS_S,         &
       PRC_HAS_N
    use scale_grid_real, only: &
       REAL_LONXY,           &
       REAL_LATXY,           &
       MY_LON  => REAL_LON,  &
       MY_LAT  => REAL_LAT,  &
       MY_LONX => REAL_LONX, &
       MY_LATX => REAL_LATX, &
       MY_LONY => REAL_LONY, &
       MY_LATY => REAL_LATY, &
       MY_CZ   => REAL_CZ,   &
       MY_FZ   => REAL_FZ,   &
       p_latlon_catalog => REAL_DOMAIN_CATALOGUE
    use scale_comm, only: &
       COMM_world
    implicit none

    integer, intent(in), optional :: icomm_parent
    integer, intent(in), optional :: icomm_child
    logical, intent(in), optional :: flag_parent
    logical, intent(in), optional :: flag_child

    !< metadata files for lat-lon domain for all processes
    character(len=H_LONG)  :: LATLON_CATALOGUE_FNAME = 'latlon_domain_catalogue.txt'

    integer :: ONLINE_SPECIFIED_MAXRQ = 0
    integer :: i
    integer :: fid, ierr
    integer :: parent_id
    integer, allocatable :: errcodes(:)

    integer :: ims, ime
    integer :: jms, jme

    character(2) :: dom_num

    namelist / PARAM_NEST /      &
       USE_NESTING,              &
       LATLON_CATALOGUE_FNAME,   &
       OFFLINE_PARENT_PRC_NUM_X, &
       OFFLINE_PARENT_PRC_NUM_Y, &
       OFFLINE_PARENT_KMAX,      &
       OFFLINE_PARENT_IMAX,      &
       OFFLINE_PARENT_JMAX,      &
       OFFLINE_PARENT_LKMAX,     &
       OFFLINE,                  &
       ONLINE_DOMAIN_NUM,        &
       ONLINE_IAM_PARENT,        &
       ONLINE_IAM_DAUGHTER,      &
       ONLINE_USE_VELZ,          &
       ONLINE_NO_ROTATE,         &
       ONLINE_BOUNDARY_USE_QHYD, &
       ONLINE_AGGRESSIVE_COMM,   &
       ONLINE_WAIT_LIMIT,        &
       ONLINE_SPECIFIED_MAXRQ,   &
       NEST_INTERP_LEVEL

    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*) ''
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[NEST]/Categ[GRID]'

    nwait_p = 0
    nwait_d = 0
    nrecv = 0
    nsend = 0

    HANDLING_NUM           = 0
    NEST_Filiation(:)      = 0
    ONLINE_WAIT_LIMIT      = 999999999
    ONLINE_AGGRESSIVE_COMM = .false.
    interp_search_divnum   = 10

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_NEST,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_NEST. Check!'
       call PRC_MPIstop
    endif
    if( IO_LNML ) write(IO_FID_LOG,nml=PARAM_NEST)

    call INTRPNEST_setup ( interp_search_divnum, NEST_INTERP_LEVEL, OFFLINE )
    itp_nh = int( NEST_INTERP_LEVEL )
    itp_nv = 2

    ! only for register
    if ( ONLINE_IAM_PARENT .or. ONLINE_IAM_DAUGHTER ) then
       call PROF_rapstart('NESTCOM total parent')
       call PROF_rapend  ('NESTCOM total parent')
       call PROF_rapstart('NESTCOM send parent')
       call PROF_rapend  ('NESTCOM send parent')
       call PROF_rapstart('NESTCOM test parent')
       call PROF_rapend  ('NESTCOM test parent')
       call PROF_rapstart('NESTCOM wait parent')
       call PROF_rapend  ('NESTCOM wait parent')
       call PROF_rapstart('NESTCOM total child')
       call PROF_rapend  ('NESTCOM total child')
       call PROF_rapstart('NESTCOM recv child')
       call PROF_rapend  ('NESTCOM recv child')
       call PROF_rapstart('NESTCOM test child')
       call PROF_rapend  ('NESTCOM test child')
       call PROF_rapstart('NESTCOM wait child')
       call PROF_rapend  ('NESTCOM wait child')
       call PROF_rapstart('NESTCOM interp')
       call PROF_rapend  ('NESTCOM interp')
    endif

    DEBUG_DOMAIN_NUM = ONLINE_DOMAIN_NUM
    if( ONLINE_SPECIFIED_MAXRQ > max_rq ) max_rq = ONLINE_SPECIFIED_MAXRQ
    allocate( ireq_p(max_rq)     )
    allocate( ireq_d(max_rq)     )
    allocate( call_order(max_rq) )

    if( USE_NESTING ) then

       if ( OFFLINE .or. ONLINE_IAM_DAUGHTER ) then

          ims = IS-1
          ime = IE
          jms = JS-1
          jme = JE
          if ( .not. PRC_HAS_W ) ims = 1
          if ( .not. PRC_HAS_E ) ime = IA
          if ( .not. PRC_HAS_S ) jms = 1
          if ( .not. PRC_HAS_N ) jme = JA
          corner_loc(I_NW,I_LON) = REAL_LONXY(ims,jme) / D2R
          corner_loc(I_NE,I_LON) = REAL_LONXY(ime,jme) / D2R
          corner_loc(I_SW,I_LON) = REAL_LONXY(ims,jms) / D2R
          corner_loc(I_SE,I_LON) = REAL_LONXY(ime,jms) / D2R
          corner_loc(I_NW,I_LAT) = REAL_LATXY(ims,jme) / D2R
          corner_loc(I_NE,I_LAT) = REAL_LATXY(ime,jme) / D2R
          corner_loc(I_SW,I_LAT) = REAL_LATXY(ims,jms) / D2R
          corner_loc(I_SE,I_LAT) = REAL_LATXY(ime,jms) / D2R
       end if

      if( OFFLINE ) then
         HANDLING_NUM = 1
         PARENT_PRC_NUM_X(HANDLING_NUM) = OFFLINE_PARENT_PRC_NUM_X
         PARENT_PRC_NUM_Y(HANDLING_NUM) = OFFLINE_PARENT_PRC_NUM_Y
         PARENT_KMAX(HANDLING_NUM)      = OFFLINE_PARENT_KMAX
         PARENT_IMAX(HANDLING_NUM)      = OFFLINE_PARENT_IMAX
         PARENT_JMAX(HANDLING_NUM)      = OFFLINE_PARENT_JMAX
         PARENT_LKMAX(HANDLING_NUM)     = OFFLINE_PARENT_LKMAX

         PARENT_PRC_nmax(HANDLING_NUM) = PARENT_PRC_NUM_X(HANDLING_NUM) * PARENT_PRC_NUM_Y(HANDLING_NUM)
         allocate( latlon_catalog(PARENT_PRC_nmax(HANDLING_NUM),4,2) )

         !--- read latlon catalogue
         fid = IO_get_available_fid()
         open( fid,                                    &
               file   = trim(LATLON_CATALOGUE_FNAME),  &
               form   = 'formatted',                   &
               status = 'old',                         &
               iostat = ierr                           )

         if ( ierr /= 0 ) then
            write(*,*) 'xxx cannot open latlon-catalogue file!'
            call PRC_MPIstop
         endif

         do i = 1, PARENT_PRC_nmax(HANDLING_NUM)
            read(fid,'(i8,8f32.24)',iostat=ierr) parent_id, &
                                                 latlon_catalog(i,I_NW,I_LON), latlon_catalog(i,I_NE,I_LON), & ! LON: NW, NE
                                                 latlon_catalog(i,I_SW,I_LON), latlon_catalog(i,I_SE,I_LON), & ! LON: SW, SE
                                                 latlon_catalog(i,I_NW,I_LAT), latlon_catalog(i,I_NE,I_LAT), & ! LAT: NW, NE
                                                 latlon_catalog(i,I_SW,I_LAT), latlon_catalog(i,I_SE,I_LAT)    ! LAT: SW, SE
            if ( i /= parent_id ) then
               if( IO_L ) write(*,*) 'xxx internal error: parent mpi id'
               call PRC_MPIstop
            endif
            if ( ierr /= 0 ) exit
         enddo
         close(fid)

         call NEST_domain_relate(HANDLING_NUM)

      else ! ONLINE RELATIONSHIP
         if ( present(flag_parent) .and. present(flag_child) ) then
            if( IO_L ) write(IO_FID_LOG,'(1x,A)') &
                       '*** Setup Online Nesting Inter-Domain Communicator (IDC)'
         else
            write(*,*) 'xxx Internal Error:'
            write(*,*) 'xxx The flag_parent and flag_child are needed.'
            write(*,*) '    domain: ', ONLINE_DOMAIN_NUM
            call PRC_MPIstop
         endif

         if( ONLINE_BOUNDARY_USE_QHYD ) then
            NEST_BND_QA = QA
         else
            NEST_BND_QA = I_QV
         endif

         if( flag_parent ) then ! must do first before daughter processes
         !-------------------------------------------------
            if ( .NOT. ONLINE_IAM_PARENT ) then
               write(*,*) 'xxx Flag from launcher is not consistent with namelist!'
               write(*,*) '    PARENT - domain: ', ONLINE_DOMAIN_NUM
               call PRC_MPIstop
            endif

            HANDLING_NUM = 1 !HANDLING_NUM + 1
            INTERCOMM_ID(HANDLING_NUM) = ONLINE_DOMAIN_NUM
            NEST_Filiation(INTERCOMM_ID(HANDLING_NUM)) = 1

            INTERCOMM_DAUGHTER = icomm_child
            if( IO_L ) write(IO_FID_LOG,'(1x,A,I2,A)') '*** Online Nesting - PARENT [INTERCOMM_ID:', &
                                                        INTERCOMM_ID(HANDLING_NUM), ' ]'
            if( IO_L ) write(IO_FID_LOG,*) '*** Online Nesting - INTERCOMM :', INTERCOMM_DAUGHTER

            call NEST_COMM_ping( HANDLING_NUM )

            call NEST_COMM_parentsize( HANDLING_NUM )

            call NEST_COMM_catalogue( HANDLING_NUM )
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

            if ( .NOT. ONLINE_NO_ROTATE ) then
               allocate( u_llp(PARENT_KA(HANDLING_NUM), PARENT_IA(HANDLING_NUM), PARENT_JA(HANDLING_NUM) ) )
               allocate( v_llp(PARENT_KA(HANDLING_NUM), PARENT_IA(HANDLING_NUM), PARENT_JA(HANDLING_NUM) ) )
               u_llp(:,:,:) = 0.0_RP
               v_llp(:,:,:) = 0.0_RP
            endif

            if( IO_L ) write(IO_FID_LOG,'(1x,A)'     ) '***  Informations of Parent Domain [me]'
            if( IO_L ) write(IO_FID_LOG,'(1x,A,I6)'  ) '***  --- PARENT_PRC_nmax   :', PARENT_PRC_nmax(HANDLING_NUM)
            if( IO_L ) write(IO_FID_LOG,'(1x,A,I6)'  ) '***  --- PARENT_PRC_NUM_X  :', PARENT_PRC_NUM_X(HANDLING_NUM)
            if( IO_L ) write(IO_FID_LOG,'(1x,A,I6)'  ) '***  --- PARENT_PRC_NUM_Y  :', PARENT_PRC_NUM_Y(HANDLING_NUM)
            if( IO_L ) write(IO_FID_LOG,'(1x,A,I6)'  ) '***  --- PARENT_KMAX       :', PARENT_KMAX(HANDLING_NUM)
            if( IO_L ) write(IO_FID_LOG,'(1x,A,I6)'  ) '***  --- PARENT_IMAX       :', PARENT_IMAX(HANDLING_NUM)
            if( IO_L ) write(IO_FID_LOG,'(1x,A,I6)'  ) '***  --- PARENT_JMAX       :', PARENT_JMAX(HANDLING_NUM)
            if( IO_L ) write(IO_FID_LOG,'(1x,A,F9.3)') '***  --- PARENT_DTSEC      :', PARENT_DTSEC(HANDLING_NUM)
            if( IO_L ) write(IO_FID_LOG,'(1x,A,I6)  ') '***  --- PARENT_NSTEP      :', PARENT_NSTEP(HANDLING_NUM)
            if( IO_L ) write(IO_FID_LOG,'(1x,A)'     ) '***  Informations of Daughter Domain'
            if( IO_L ) write(IO_FID_LOG,'(1x,A,I6)'  ) '***  --- DAUGHTER_PRC_nmax :', DAUGHTER_PRC_nmax(HANDLING_NUM)
            if( IO_L ) write(IO_FID_LOG,'(1x,A,I6)'  ) '***  --- DAUGHTER_PRC_NUM_X:', DAUGHTER_PRC_NUM_X(HANDLING_NUM)
            if( IO_L ) write(IO_FID_LOG,'(1x,A,I6)'  ) '***  --- DAUGHTER_PRC_NUM_Y:', DAUGHTER_PRC_NUM_Y(HANDLING_NUM)
            if( IO_L ) write(IO_FID_LOG,'(1x,A,I6)'  ) '***  --- DAUGHTER_KMAX     :', DAUGHTER_KMAX(HANDLING_NUM)
            if( IO_L ) write(IO_FID_LOG,'(1x,A,I6)'  ) '***  --- DAUGHTER_IMAX     :', DAUGHTER_IMAX(HANDLING_NUM)
            if( IO_L ) write(IO_FID_LOG,'(1x,A,I6)'  ) '***  --- DAUGHTER_JMAX     :', DAUGHTER_JMAX(HANDLING_NUM)
            if( IO_L ) write(IO_FID_LOG,'(1x,A,F9.3)') '***  --- DAUGHTER_DTSEC    :', DAUGHTER_DTSEC(HANDLING_NUM)
            if( IO_L ) write(IO_FID_LOG,'(1x,A,I6)  ') '***  --- DAUGHTER_NSTEP    :', DAUGHTER_NSTEP(HANDLING_NUM)
            if( IO_L ) write(IO_FID_LOG,'(1x,A,I6)  ') '***  Limit Num. NCOMM req. :', max_rq

            allocate( org_DENS(PARENT_KA(HANDLING_NUM),PARENT_IA(HANDLING_NUM),PARENT_JA(HANDLING_NUM))           )
            allocate( org_MOMZ(PARENT_KA(HANDLING_NUM),PARENT_IA(HANDLING_NUM),PARENT_JA(HANDLING_NUM))           )
            allocate( org_MOMX(PARENT_KA(HANDLING_NUM),PARENT_IA(HANDLING_NUM),PARENT_JA(HANDLING_NUM))           )
            allocate( org_MOMY(PARENT_KA(HANDLING_NUM),PARENT_IA(HANDLING_NUM),PARENT_JA(HANDLING_NUM))           )
            allocate( org_RHOT(PARENT_KA(HANDLING_NUM),PARENT_IA(HANDLING_NUM),PARENT_JA(HANDLING_NUM))           )
            allocate( org_QTRC(PARENT_KA(HANDLING_NUM),PARENT_IA(HANDLING_NUM),PARENT_JA(HANDLING_NUM),max_bndqa) )

            call NEST_COMM_setup_nestdown( HANDLING_NUM )

         !---------------------------------- end of parent routines
         endif


         if( flag_child ) then
         !-------------------------------------------------
            if ( .NOT. ONLINE_IAM_DAUGHTER ) then
               write(*,*) 'xxx Flag from launcher is not consistent with namelist!'
               write(*,*) '    DAUGHTER - domain: ', ONLINE_DOMAIN_NUM
               call PRC_MPIstop
            endif

            HANDLING_NUM = 2 !HANDLING_NUM + 1
            INTERCOMM_ID(HANDLING_NUM) = ONLINE_DOMAIN_NUM - 1
            NEST_Filiation(INTERCOMM_ID(HANDLING_NUM)) = -1

            INTERCOMM_PARENT = icomm_parent
            if( IO_L ) write(IO_FID_LOG,'(1x,A,I2,A)') '*** Online Nesting - DAUGHTER [INTERCOMM_ID:', &
                                                        INTERCOMM_ID(HANDLING_NUM), ' ]'
            if( IO_L ) write(IO_FID_LOG,*) '*** Online Nesting - INTERCOMM :', INTERCOMM_PARENT

            call NEST_COMM_ping( HANDLING_NUM )

            call NEST_COMM_parentsize( HANDLING_NUM )

            allocate( latlon_catalog(PARENT_PRC_nmax(HANDLING_NUM),4,2) )
            call NEST_COMM_catalogue( HANDLING_NUM )
            call MPI_BARRIER(INTERCOMM_PARENT, ierr)

            call NEST_domain_relate( HANDLING_NUM )

            PARENT_KA(HANDLING_NUM)   = PARENT_KMAX(HANDLING_NUM)   + KHALO * 2
            PARENT_IA(HANDLING_NUM)   = PARENT_IMAX(HANDLING_NUM)   + IHALO * 2
            PARENT_JA(HANDLING_NUM)   = PARENT_JMAX(HANDLING_NUM)   + JHALO * 2
            DAUGHTER_KA(HANDLING_NUM) = DAUGHTER_KMAX(HANDLING_NUM) + KHALO * 2
            DAUGHTER_IA(HANDLING_NUM) = DAUGHTER_IMAX(HANDLING_NUM) + IHALO * 2
            DAUGHTER_JA(HANDLING_NUM) = DAUGHTER_JMAX(HANDLING_NUM) + JHALO * 2
            TILEAL_KA(HANDLING_NUM)   = PARENT_KA(HANDLING_NUM)
            TILEAL_IA(HANDLING_NUM)   = PARENT_IMAX(HANDLING_NUM) * NEST_TILE_NUM_X
            TILEAL_JA(HANDLING_NUM)   = PARENT_JMAX(HANDLING_NUM) * NEST_TILE_NUM_Y

            if( IO_L ) write(IO_FID_LOG,'(1x,A)'     ) '***  Informations of Parent Domain'
            if( IO_L ) write(IO_FID_LOG,'(1x,A,I6)'  ) '***  --- PARENT_PRC_nmax   :', PARENT_PRC_nmax(HANDLING_NUM)
            if( IO_L ) write(IO_FID_LOG,'(1x,A,I6)'  ) '***  --- PARENT_PRC_NUM_X  :', PARENT_PRC_NUM_X(HANDLING_NUM)
            if( IO_L ) write(IO_FID_LOG,'(1x,A,I6)'  ) '***  --- PARENT_PRC_NUM_Y  :', PARENT_PRC_NUM_Y(HANDLING_NUM)
            if( IO_L ) write(IO_FID_LOG,'(1x,A,I6)'  ) '***  --- PARENT_KMAX       :', PARENT_KMAX(HANDLING_NUM)
            if( IO_L ) write(IO_FID_LOG,'(1x,A,I6)'  ) '***  --- PARENT_IMAX       :', PARENT_IMAX(HANDLING_NUM)
            if( IO_L ) write(IO_FID_LOG,'(1x,A,I6)'  ) '***  --- PARENT_JMAX       :', PARENT_JMAX(HANDLING_NUM)
            if( IO_L ) write(IO_FID_LOG,'(1x,A,F9.3)') '***  --- PARENT_DTSEC      :', PARENT_DTSEC(HANDLING_NUM)
            if( IO_L ) write(IO_FID_LOG,'(1x,A,I6)'  ) '***  --- PARENT_NSTEP      :', PARENT_NSTEP(HANDLING_NUM)
            if( IO_L ) write(IO_FID_LOG,'(1x,A)'     ) '***  Informations of Daughter Domain [me]'
            if( IO_L ) write(IO_FID_LOG,'(1x,A,I6)'  ) '***  --- DAUGHTER_PRC_nmax :', DAUGHTER_PRC_nmax(HANDLING_NUM)
            if( IO_L ) write(IO_FID_LOG,'(1x,A,I6)'  ) '***  --- DAUGHTER_PRC_NUM_X:', DAUGHTER_PRC_NUM_X(HANDLING_NUM)
            if( IO_L ) write(IO_FID_LOG,'(1x,A,I6)'  ) '***  --- DAUGHTER_PRC_NUM_Y:', DAUGHTER_PRC_NUM_Y(HANDLING_NUM)
            if( IO_L ) write(IO_FID_LOG,'(1x,A,I6)'  ) '***  --- DAUGHTER_KMAX     :', DAUGHTER_KMAX(HANDLING_NUM)
            if( IO_L ) write(IO_FID_LOG,'(1x,A,I6)'  ) '***  --- DAUGHTER_IMAX     :', DAUGHTER_IMAX(HANDLING_NUM)
            if( IO_L ) write(IO_FID_LOG,'(1x,A,I6)'  ) '***  --- DAUGHTER_JMAX     :', DAUGHTER_JMAX(HANDLING_NUM)
            if( IO_L ) write(IO_FID_LOG,'(1x,A,F9.3)') '***  --- DAUGHTER_DTSEC    :', DAUGHTER_DTSEC(HANDLING_NUM)
            if( IO_L ) write(IO_FID_LOG,'(1x,A,I6)'  ) '***  --- DAUGHTER_NSTEP    :', DAUGHTER_NSTEP(HANDLING_NUM)
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
            allocate( buffer_ref_LONX(                          TILEAL_IA(HANDLING_NUM),TILEAL_JA(HANDLING_NUM)) )
            allocate( buffer_ref_LONY(                          TILEAL_IA(HANDLING_NUM),TILEAL_JA(HANDLING_NUM)) )
            allocate( buffer_ref_LAT (                          TILEAL_IA(HANDLING_NUM),TILEAL_JA(HANDLING_NUM)) )
            allocate( buffer_ref_LATX(                          TILEAL_IA(HANDLING_NUM),TILEAL_JA(HANDLING_NUM)) )
            allocate( buffer_ref_LATY(                          TILEAL_IA(HANDLING_NUM),TILEAL_JA(HANDLING_NUM)) )
            allocate( buffer_ref_CZ  (  PARENT_KA(HANDLING_NUM),TILEAL_IA(HANDLING_NUM),TILEAL_JA(HANDLING_NUM)) )
            allocate( buffer_ref_FZ  (0:PARENT_KA(HANDLING_NUM),TILEAL_IA(HANDLING_NUM),TILEAL_JA(HANDLING_NUM)) )

            !allocate( buffer_ref_2D  (                          TILEAL_IA(HANDLING_NUM),TILEAL_JA(HANDLING_NUM)) )
            allocate( buffer_ref_3D  (  PARENT_KA(HANDLING_NUM),TILEAL_IA(HANDLING_NUM),TILEAL_JA(HANDLING_NUM)) )
            allocate( buffer_ref_3DF (0:PARENT_KA(HANDLING_NUM),TILEAL_IA(HANDLING_NUM),TILEAL_JA(HANDLING_NUM)) )

            allocate( hfact(            DAUGHTER_IA(HANDLING_NUM),DAUGHTER_JA(HANDLING_NUM),itp_nh,       itp_ng) )
            allocate( vfact(DAUGHTER_KA(HANDLING_NUM),DAUGHTER_IA(HANDLING_NUM),DAUGHTER_JA(HANDLING_NUM),itp_nh,itp_nv,itp_ng) )
            allocate( igrd (            DAUGHTER_IA(HANDLING_NUM),DAUGHTER_JA(HANDLING_NUM),itp_nh,       itp_ng) )
            allocate( jgrd (            DAUGHTER_IA(HANDLING_NUM),DAUGHTER_JA(HANDLING_NUM),itp_nh,       itp_ng) )
            allocate( kgrd (DAUGHTER_KA(HANDLING_NUM),DAUGHTER_IA(HANDLING_NUM),DAUGHTER_JA(HANDLING_NUM),itp_nh,itp_nv,itp_ng) )

            call NEST_COMM_setup_nestdown( HANDLING_NUM )


            ! for scalar points
!            call NEST_latlonz_interporation_fact( hfact          (:,:,:,I_SCLR),     & ! [OUT]
!                                                  vfact          (:,:,:,:,:,I_SCLR), & ! [OUT]
!                                                  kgrd           (:,:,:,:,:,I_SCLR), & ! [OUT]
!                                                  igrd           (:,:,:,I_SCLR),     & ! [OUT]
!                                                  jgrd           (:,:,:,I_SCLR),     & ! [OUT]
!                                                  MY_CZ          (:,:,:),            & ! [IN]
!                                                  MY_LAT         (:,:),              & ! [IN]
!                                                  MY_LON         (:,:),              & ! [IN]
!                                                  buffer_ref_CZ  (:,:,:),            & ! [IN]
!                                                  buffer_ref_LAT (:,:),              & ! [IN]
!                                                  buffer_ref_LON (:,:),              & ! [IN]
!                                                  TILEAL_KA(HANDLING_NUM),           & ! [IN]
!                                                  TILEAL_IA(HANDLING_NUM),           & ! [IN]
!                                                  TILEAL_JA(HANDLING_NUM),           & ! [IN]
!                                                  HANDLING_NUM                       ) ! [IN]
            call INTRPNEST_interp_fact_llz( hfact          (:,:,:,I_SCLR),     & ! [OUT]
                                            vfact          (:,:,:,:,:,I_SCLR), & ! [OUT]
                                            kgrd           (:,:,:,:,:,I_SCLR), & ! [OUT]
                                            igrd           (:,:,:,I_SCLR),     & ! [OUT]
                                            jgrd           (:,:,:,I_SCLR),     & ! [OUT]
                                            MY_CZ          (:,:,:),            & ! [IN]
                                            MY_LAT         (:,:),              & ! [IN]
                                            MY_LON         (:,:),              & ! [IN]
                                            DATR_KS(HANDLING_NUM),             & ! [IN]
                                            DATR_KE(HANDLING_NUM),             & ! [IN]
                                            DAUGHTER_IA(HANDLING_NUM),         & ! [IN]
                                            DAUGHTER_JA(HANDLING_NUM),         & ! [IN]
                                            buffer_ref_CZ  (:,:,:),            & ! [IN]
                                            buffer_ref_LAT (:,:),              & ! [IN]
                                            buffer_ref_LON (:,:),              & ! [IN]
                                            TILEAL_KA(HANDLING_NUM),           & ! [IN]
                                            TILEAL_IA(HANDLING_NUM),           & ! [IN]
                                            TILEAL_JA(HANDLING_NUM)            ) ! [IN]


            ! for z staggered points
!            call NEST_latlonz_interporation_fact( hfact          (:,:,:,I_ZSTG),     & ! [OUT]
!                                                  vfact          (:,:,:,:,:,I_ZSTG), & ! [OUT]
!                                                  kgrd           (:,:,:,:,:,I_ZSTG), & ! [OUT]
!                                                  igrd           (:,:,:,I_ZSTG),     & ! [OUT]
!                                                  jgrd           (:,:,:,I_ZSTG),     & ! [OUT]
!                                                  MY_FZ          (:,:,:),            & ! [IN]
!                                                  MY_LAT         (:,:),              & ! [IN]
!                                                  MY_LON         (:,:),              & ! [IN]
!                                                  buffer_ref_FZ  (:,:,:),            & ! [IN]
!                                                  buffer_ref_LAT (:,:),              & ! [IN]
!                                                  buffer_ref_LON (:,:),              & ! [IN]
!                                                  TILEAL_KA(HANDLING_NUM)+1,         & ! [IN]
!                                                  TILEAL_IA(HANDLING_NUM),           & ! [IN]
!                                                  TILEAL_JA(HANDLING_NUM),           & ! [IN]
!                                                  HANDLING_NUM                       ) ! [IN]
            call INTRPNEST_interp_fact_llz( hfact          (:,:,:,I_ZSTG),     & ! [OUT]
                                            vfact          (:,:,:,:,:,I_ZSTG), & ! [OUT]
                                            kgrd           (:,:,:,:,:,I_ZSTG), & ! [OUT]
                                            igrd           (:,:,:,I_ZSTG),     & ! [OUT]
                                            jgrd           (:,:,:,I_ZSTG),     & ! [OUT]
                                            MY_FZ          (:,:,:),            & ! [IN]
                                            MY_LAT         (:,:),              & ! [IN]
                                            MY_LON         (:,:),              & ! [IN]
                                            DATR_KS(HANDLING_NUM),             & ! [IN]
                                            DATR_KE(HANDLING_NUM),             & ! [IN]
                                            DAUGHTER_IA(HANDLING_NUM),         & ! [IN]
                                            DAUGHTER_JA(HANDLING_NUM),         & ! [IN]
                                            buffer_ref_FZ  (:,:,:),            & ! [IN]
                                            buffer_ref_LAT (:,:),              & ! [IN]
                                            buffer_ref_LON (:,:),              & ! [IN]
                                            TILEAL_KA(HANDLING_NUM)+1,         & ! [IN]
                                            TILEAL_IA(HANDLING_NUM),           & ! [IN]
                                            TILEAL_JA(HANDLING_NUM)            ) ! [IN]


            ! for x staggered points
!            call NEST_latlonz_interporation_fact( hfact          (:,:,:,I_XSTG),     & ! [OUT]
!                                                  vfact          (:,:,:,:,:,I_XSTG), & ! [OUT]
!                                                  kgrd           (:,:,:,:,:,I_XSTG), & ! [OUT]
!                                                  igrd           (:,:,:,I_XSTG),     & ! [OUT]
!                                                  jgrd           (:,:,:,I_XSTG),     & ! [OUT]
!                                                  MY_CZ          (:,:,:),            & ! [IN]
!                                                  MY_LATX        (:,:),              & ! [IN]
!                                                  MY_LONX        (:,:),              & ! [IN]
!                                                  buffer_ref_CZ  (:,:,:),            & ! [IN]
!                                                  buffer_ref_LATX(:,:),              & ! [IN]
!                                                  buffer_ref_LONX(:,:),              & ! [IN]
!                                                  TILEAL_KA(HANDLING_NUM),           & ! [IN]
!                                                  TILEAL_IA(HANDLING_NUM),           & ! [IN]
!                                                  TILEAL_JA(HANDLING_NUM),           & ! [IN]
!                                                  HANDLING_NUM                       ) ! [IN]
            call INTRPNEST_interp_fact_llz( hfact          (:,:,:,I_XSTG),     & ! [OUT]
                                            vfact          (:,:,:,:,:,I_XSTG), & ! [OUT]
                                            kgrd           (:,:,:,:,:,I_XSTG), & ! [OUT]
                                            igrd           (:,:,:,I_XSTG),     & ! [OUT]
                                            jgrd           (:,:,:,I_XSTG),     & ! [OUT]
                                            MY_CZ          (:,:,:),            & ! [IN]
                                            MY_LATX        (:,:),              & ! [IN]
                                            MY_LONX        (:,:),              & ! [IN]
                                            DATR_KS(HANDLING_NUM),             & ! [IN]
                                            DATR_KE(HANDLING_NUM),             & ! [IN]
                                            DAUGHTER_IA(HANDLING_NUM),         & ! [IN]
                                            DAUGHTER_JA(HANDLING_NUM),         & ! [IN]
                                            buffer_ref_CZ  (:,:,:),            & ! [IN]
                                            buffer_ref_LATX(:,:),              & ! [IN]
                                            buffer_ref_LONX(:,:),              & ! [IN]
                                            TILEAL_KA(HANDLING_NUM),           & ! [IN]
                                            TILEAL_IA(HANDLING_NUM),           & ! [IN]
                                            TILEAL_JA(HANDLING_NUM)            ) ! [IN]

            ! for y staggered points
!            call NEST_latlonz_interporation_fact( hfact          (:,:,:,I_YSTG),     & ! [OUT]
!                                                  vfact          (:,:,:,:,:,I_YSTG), & ! [OUT]
!                                                  kgrd           (:,:,:,:,:,I_YSTG), & ! [OUT]
!                                                  igrd           (:,:,:,I_YSTG),     & ! [OUT]
!                                                  jgrd           (:,:,:,I_YSTG),     & ! [OUT]
!                                                  MY_CZ          (:,:,:),            & ! [IN]
!                                                  MY_LATY        (:,:),              & ! [IN]
!                                                  MY_LONY        (:,:),              & ! [IN]
!                                                  buffer_ref_CZ  (:,:,:),            & ! [IN]
!                                                  buffer_ref_LATY(:,:),              & ! [IN]
!                                                  buffer_ref_LONY(:,:),              & ! [IN]
!                                                  TILEAL_KA(HANDLING_NUM),           & ! [IN]
!                                                  TILEAL_IA(HANDLING_NUM),           & ! [IN]
!                                                  TILEAL_JA(HANDLING_NUM),           & ! [IN]
!                                                  HANDLING_NUM                       ) ! [IN]
            call INTRPNEST_interp_fact_llz( hfact          (:,:,:,I_YSTG),     & ! [OUT]
                                            vfact          (:,:,:,:,:,I_YSTG), & ! [OUT]
                                            kgrd           (:,:,:,:,:,I_YSTG), & ! [OUT]
                                            igrd           (:,:,:,I_YSTG),     & ! [OUT]
                                            jgrd           (:,:,:,I_YSTG),     & ! [OUT]
                                            MY_CZ          (:,:,:),            & ! [IN]
                                            MY_LATY        (:,:),              & ! [IN]
                                            MY_LONY        (:,:),              & ! [IN]
                                            DATR_KS(HANDLING_NUM),             & ! [IN]
                                            DATR_KE(HANDLING_NUM),             & ! [IN]
                                            DAUGHTER_IA(HANDLING_NUM),         & ! [IN]
                                            DAUGHTER_JA(HANDLING_NUM),         & ! [IN]
                                            buffer_ref_CZ  (:,:,:),            & ! [IN]
                                            buffer_ref_LATY(:,:),              & ! [IN]
                                            buffer_ref_LONY(:,:),              & ! [IN]
                                            TILEAL_KA(HANDLING_NUM),           & ! [IN]
                                            TILEAL_IA(HANDLING_NUM),           & ! [IN]
                                            TILEAL_JA(HANDLING_NUM)            ) ! [IN]

            deallocate( buffer_2D  )
            deallocate( buffer_3D  )
            deallocate( buffer_3DF )

         !---------------------------------- end of daughter routines
         else
            ONLINE_USE_VELZ = .false.
         endif

         !if( IO_L ) write(IO_FID_LOG,'(1x,A,I2)') '*** Number of Related Domains :', HANDLING_NUM
         !if ( HANDLING_NUM > 2 ) then
         !   if( IO_L ) write(*,*) 'xxx Too much handing domains (up to 2)'
         !   call PRC_MPIstop
         !endif

      endif !--- OFFLINE or NOT

    endif !--- USE_NESTING

    return
  end subroutine NEST_setup

  !-----------------------------------------------------------------------------
  !> Solve relationship between ParentDomain & Daughter Domain
  subroutine NEST_domain_relate( &
      HANDLE  )
    use scale_const, only: &
       EPS => CONST_EPS
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

    allocate( pd_tile_num(0:PARENT_PRC_nmax(HANDLE)-1,2) )

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
    do i = 1, PARENT_PRC_nmax(HANDLE)
       wid_lon = abs((latlon_catalog(i,I_SW,I_LON) - latlon_catalog(i,I_SE,I_LON)) &
                      / real( PARENT_IMAX(HANDLE)-1, kind=RP )) * 0.8_RP
       wid_lat = abs((latlon_catalog(i,I_SW,I_LAT) - latlon_catalog(i,I_NW,I_LAT)) &
                      / real( PARENT_JMAX(HANDLE)-1, kind=RP )) * 0.8_RP

       if ( corner_loc(I_SW,I_LON) >= min(latlon_catalog(i,I_SW,I_LON),latlon_catalog(i,I_NW,I_LON))-wid_lon .and. &
            corner_loc(I_SW,I_LAT) >= min(latlon_catalog(i,I_SW,I_LAT),latlon_catalog(i,I_SE,I_LAT))-wid_lat .and. &
            corner_loc(I_SW,I_LON) <= max(latlon_catalog(i,I_NE,I_LON),latlon_catalog(i,I_SE,I_LON))+wid_lon .and. &
            corner_loc(I_SW,I_LAT) <= max(latlon_catalog(i,I_NE,I_LAT),latlon_catalog(i,I_NW,I_LAT))+wid_lat ) then

          pd_sw_tile = i-1 ! MPI process number starts from zero
          hit = .true.
          exit ! exit loop
       endif
    enddo
    if ( .NOT. hit ) then
       write(*,*) 'xxx region of daughter domain is larger than that of parent: SW search'
       write(*,*) '    at rank:', PRC_myrank, ' of domain:', ONLINE_DOMAIN_NUM
       if( IO_L ) write(IO_FID_LOG,'(1x,A)') 'xxx region of daughter domain is larger than that of parent: SW search'
       if( IO_L ) write(IO_FID_LOG,*) ' grid width: half width in lat:', wid_lat, ' half width in lon:', wid_lon
       if( IO_L ) write(IO_FID_LOG,'(1x,A,F12.6)') '    daughter local (me): LON=',corner_loc(I_SW,I_LON)
       do i = 1, PARENT_PRC_nmax(HANDLE)
          if( IO_L ) write(IO_FID_LOG,'(1x,A,F12.6,1x,F12.6)') '     parent local SW-NE: LON=', &
                     latlon_catalog(i,I_SW,I_LON) ,latlon_catalog(i,I_NE,I_LON)
       enddo
       if( IO_L ) write(IO_FID_LOG,'(1x,A,F12.6)') '    daughter local (me): LAT=',corner_loc(I_SW,I_LAT)
       do i = 1, PARENT_PRC_nmax(HANDLE)
          if( IO_L ) write(IO_FID_LOG,'(1x,A,F12.6,1x,F12.6)') '     parent local SW-NE: LAT=', &
                     latlon_catalog(i,I_SW,I_LAT) ,latlon_catalog(i,I_NE,I_LAT)
       enddo
       call PRC_MPIstop
    endif

    !--- NE search
    hit = .false.
    do i = PARENT_PRC_nmax(HANDLE), 1, -1
       wid_lon = abs((latlon_catalog(i,I_NW,I_LON) - latlon_catalog(i,I_NE,I_LON)) &
                      / real( PARENT_IMAX(HANDLE)-1, kind=RP )) * 0.8_RP
       wid_lat = abs((latlon_catalog(i,I_SE,I_LAT) - latlon_catalog(i,I_NE,I_LAT)) &
                      / real( PARENT_JMAX(HANDLE)-1, kind=RP )) * 0.8_RP

       if ( corner_loc(I_NE,I_LON) >= min(latlon_catalog(i,I_SW,I_LON),latlon_catalog(i,I_NW,I_LON))-wid_lon .and. &
            corner_loc(I_NE,I_LAT) >= min(latlon_catalog(i,I_SW,I_LAT),latlon_catalog(i,I_SE,I_LAT))-wid_lat .and. &
            corner_loc(I_NE,I_LON) <= max(latlon_catalog(i,I_NE,I_LON),latlon_catalog(i,I_SE,I_LON))+wid_lon .and. &
            corner_loc(I_NE,I_LAT) <= max(latlon_catalog(i,I_NE,I_LAT),latlon_catalog(i,I_NW,I_LAT))+wid_lat ) then

          pd_ne_tile = i-1 ! MPI process number starts from zero
          hit = .true.
          exit ! exit loop
       endif
    enddo
    if ( .NOT. hit ) then
       write(*,*) 'xxx region of daughter domain is larger than that of parent: NE search'
       write(*,*) '    at rank:', PRC_myrank, ' of domain:', ONLINE_DOMAIN_NUM
       if( IO_L ) write(IO_FID_LOG,'(1x,A)') 'xxx region of daughter domain is larger than that of parent: NE search'
       if( IO_L ) write(IO_FID_LOG,*) ' grid width: half width in lat:', wid_lat, ' half width in lon:', wid_lon
       if( IO_L ) write(IO_FID_LOG,'(1x,A,F12.6)') '    daughter local (me): LON=',corner_loc(I_NE,I_LON)
       do i = 1, PARENT_PRC_nmax(HANDLE)
          if( IO_L ) write(IO_FID_LOG,'(1x,A,F12.6,1x,F12.6)') '     parent local SW-NE: LON=', &
                     latlon_catalog(i,I_SW,I_LON) ,latlon_catalog(i,I_NE,I_LON)
       enddo
       if( IO_L ) write(IO_FID_LOG,'(1x,A,F12.6)') '    daughter local (me): LAT=',corner_loc(I_NE,I_LAT)
       do i = 1, PARENT_PRC_nmax(HANDLE)
          if( IO_L ) write(IO_FID_LOG,'(1x,A,F12.6,1x,F12.6)') '     parent local SW-NE: LAT=', &
                     latlon_catalog(i,I_SW,I_LAT) ,latlon_catalog(i,I_NE,I_LAT)
       enddo
       call PRC_MPIstop
    endif

    NEST_TILE_NUM_X = pd_tile_num(pd_ne_tile,1) - pd_tile_num(pd_sw_tile,1) + 1
    NEST_TILE_NUM_Y = pd_tile_num(pd_ne_tile,2) - pd_tile_num(pd_sw_tile,2) + 1

    allocate( NEST_TILE_ID( NEST_TILE_NUM_X*NEST_TILE_NUM_Y ) )

    if( IO_L ) write(IO_FID_LOG,'(1x,A)') '*** NEST: target process tile in parent domain'
    k = 1
    do j = 1, NEST_TILE_NUM_Y
    do i = 1, NEST_TILE_NUM_X
       NEST_TILE_ID(k) = pd_sw_tile + (i-1) + PARENT_PRC_NUM_X(HANDLE)*(j-1)
       if( IO_L ) write(IO_FID_LOG,'(1x,A,I4,A,I6)') '    (', k, ') target mpi-process:', NEST_TILE_ID(k)
       k = k + 1
    enddo
    enddo

    return
  end subroutine NEST_domain_relate

  !-----------------------------------------------------------------------------
  !> Get parent domain size
  subroutine NEST_COMM_parentsize( &
      HANDLE  )
    use scale_process, only: &
       PRC_myrank,  &
       PRC_master,  &
       PRC_nmax,    &
       PRC_NUM_X,   &
       PRC_NUM_Y,   &
       PRC_MPIstop
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

    tag   = INTERCOMM_ID(HANDLE) * 100
    ileng = 14

    if ( NEST_Filiation( INTERCOMM_ID(HANDLE) ) > 0 ) then !--- parent
       ! from parent to daughter
       datapack( 1) = PRC_nmax
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
       datapack(14) = NEST_BND_QA
       buffer       = TIME_DTSEC

       if ( PRC_myrank == PRC_master ) then
          call MPI_ISEND(datapack, ileng, MPI_INTEGER, PRC_myrank, tag, INTERCOMM_DAUGHTER, ireq1, ierr1)
          call MPI_ISEND(buffer, 1, MPI_DOUBLE_PRECISION, PRC_myrank, tag+1, INTERCOMM_DAUGHTER, ireq2, ierr2)
          call MPI_WAIT(ireq1, istatus, ierr1)
          call MPI_WAIT(ireq2, istatus, ierr2)
       endif

       PARENT_PRC_nmax(HANDLE)  = datapack( 1)
       PARENT_PRC_NUM_X(HANDLE) = datapack( 2)
       PARENT_PRC_NUM_Y(HANDLE) = datapack( 3)
       PARENT_KMAX(HANDLE)      = datapack( 4)
       PARENT_IMAX(HANDLE)      = datapack( 5)
       PARENT_JMAX(HANDLE)      = datapack( 6)
       PRNT_KS(HANDLE)          = datapack( 7)
       PRNT_KE(HANDLE)          = datapack( 8)
       PRNT_IS(HANDLE)          = datapack( 9)
       PRNT_IE(HANDLE)          = datapack(10)
       PRNT_JS(HANDLE)          = datapack(11)
       PRNT_JE(HANDLE)          = datapack(12)
       PARENT_NSTEP(HANDLE)     = datapack(13)
       PARENT_DTSEC(HANDLE)     = buffer

       ! from daughter to parent
       if ( PRC_myrank == PRC_master ) then
          call MPI_IRECV(datapack, ileng, MPI_INTEGER, PRC_myrank, tag+2, INTERCOMM_DAUGHTER, ireq1, ierr1)
          call MPI_IRECV(buffer, 1, MPI_DOUBLE_PRECISION, PRC_myrank, tag+3, INTERCOMM_DAUGHTER, ireq2, ierr2)
          call MPI_WAIT(ireq1, istatus, ierr1)
          call MPI_WAIT(ireq2, istatus, ierr2)
       endif
       call COMM_bcast(datapack, ileng)
       call COMM_bcast(buffer)

       DAUGHTER_PRC_nmax(HANDLE)  = datapack( 1)
       DAUGHTER_PRC_NUM_X(HANDLE) = datapack( 2)
       DAUGHTER_PRC_NUM_Y(HANDLE) = datapack( 3)
       DAUGHTER_KMAX(HANDLE)      = datapack( 4)
       DAUGHTER_IMAX(HANDLE)      = datapack( 5)
       DAUGHTER_JMAX(HANDLE)      = datapack( 6)
       DATR_KS(HANDLE)            = datapack( 7)
       DATR_KE(HANDLE)            = datapack( 8)
       DATR_IS(HANDLE)            = datapack( 9)
       DATR_IE(HANDLE)            = datapack(10)
       DATR_JS(HANDLE)            = datapack(11)
       DATR_JE(HANDLE)            = datapack(12)
       DAUGHTER_NSTEP(HANDLE)     = datapack(13)
       QA_OTHERSIDE               = datapack(14)
       DAUGHTER_DTSEC(HANDLE)     = buffer


    elseif ( NEST_Filiation( INTERCOMM_ID(HANDLE) ) < 0 ) then !--- daughter
       ! from parent to daughter
       if ( PRC_myrank == PRC_master ) then
          call MPI_IRECV(datapack, ileng, MPI_INTEGER, PRC_myrank, tag, INTERCOMM_PARENT, ireq1, ierr1)
          call MPI_IRECV(buffer, 1, MPI_DOUBLE_PRECISION, PRC_myrank, tag+1, INTERCOMM_PARENT, ireq2, ierr2)
          call MPI_WAIT(ireq1, istatus, ierr1)
          call MPI_WAIT(ireq2, istatus, ierr2)
       endif
       call COMM_bcast(datapack, ileng)
       call COMM_bcast(buffer)

       PARENT_PRC_nmax(HANDLE)  = datapack( 1)
       PARENT_PRC_NUM_X(HANDLE) = datapack( 2)
       PARENT_PRC_NUM_Y(HANDLE) = datapack( 3)
       PARENT_KMAX(HANDLE)      = datapack( 4)
       PARENT_IMAX(HANDLE)      = datapack( 5)
       PARENT_JMAX(HANDLE)      = datapack( 6)
       PRNT_KS(HANDLE)          = datapack( 7)
       PRNT_KE(HANDLE)          = datapack( 8)
       PRNT_IS(HANDLE)          = datapack( 9)
       PRNT_IE(HANDLE)          = datapack(10)
       PRNT_JS(HANDLE)          = datapack(11)
       PRNT_JE(HANDLE)          = datapack(12)
       PARENT_NSTEP(HANDLE)     = datapack(13)
       QA_OTHERSIDE             = datapack(14)
       PARENT_DTSEC(HANDLE)     = buffer

       ! from daughter to parent
       datapack( 1) = PRC_nmax
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
       datapack(14) = NEST_BND_QA
       buffer       = TIME_DTSEC

       if ( PRC_myrank == PRC_master ) then
          call MPI_ISEND(datapack, ileng, MPI_INTEGER, PRC_myrank, tag+2, INTERCOMM_PARENT, ireq1, ierr1)
          call MPI_ISEND(buffer, 1, MPI_DOUBLE_PRECISION, PRC_myrank, tag+3, INTERCOMM_PARENT, ireq2, ierr2)
          call MPI_WAIT(ireq1, istatus, ierr1)
          call MPI_WAIT(ireq2, istatus, ierr2)
       endif

       DAUGHTER_PRC_nmax(HANDLE)  = datapack( 1)
       DAUGHTER_PRC_NUM_X(HANDLE) = datapack( 2)
       DAUGHTER_PRC_NUM_Y(HANDLE) = datapack( 3)
       DAUGHTER_KMAX(HANDLE)      = datapack( 4)
       DAUGHTER_IMAX(HANDLE)      = datapack( 5)
       DAUGHTER_JMAX(HANDLE)      = datapack( 6)
       DATR_KS(HANDLE)            = datapack( 7)
       DATR_KE(HANDLE)            = datapack( 8)
       DATR_IS(HANDLE)            = datapack( 9)
       DATR_IE(HANDLE)            = datapack(10)
       DATR_JS(HANDLE)            = datapack(11)
       DATR_JE(HANDLE)            = datapack(12)
       DAUGHTER_NSTEP(HANDLE)     = datapack(13)
       DAUGHTER_DTSEC(HANDLE)     = buffer
    else
       write(*,*) 'xxx internal error [parentsize: nest/grid]'
       call PRC_MPIstop
    endif

    if( QA_OTHERSIDE /= NEST_BND_QA ) then
       write(*,*) 'xxx ERROR: NUMBER of QA are not matched! [parentsize: nest/grid]'
       write(*,*) 'xxx check a flag of ONLINE_BOUNDARY_USE_QHYD.', QA_OTHERSIDE, NEST_BND_QA
       call PRC_MPIstop
    endif

    return
  end subroutine NEST_COMM_parentsize

  !-----------------------------------------------------------------------------
  !> Get parent latlon catalogue
  subroutine NEST_COMM_catalogue( &
      HANDLE  )
    use scale_process, only: &
       PRC_myrank,  &
       PRC_master,  &
       PRC_nmax,    &
       PRC_MPIstop
    use scale_grid_real, only: &
       REAL_DOMAIN_CATALOGUE
    use scale_comm, only: &
       COMM_datatype,  &
       COMM_bcast
    implicit none

    integer, intent(in) :: HANDLE !< id number of nesting relation in this process target

    integer :: ireq, ierr, ileng
    integer :: istatus(MPI_STATUS_SIZE)
    integer :: tag
    !---------------------------------------------------------------------------

    tag = INTERCOMM_ID(HANDLE) * 100

    if ( NEST_Filiation( INTERCOMM_ID(HANDLE) ) > 0 ) then !--- parent
       ileng = PRC_nmax * 4 * 2
       if ( PRC_myrank == PRC_master ) then
          call MPI_ISEND(REAL_DOMAIN_CATALOGUE, ileng, COMM_datatype, PRC_myrank, tag, INTERCOMM_DAUGHTER, ireq, ierr)
          call MPI_WAIT(ireq, istatus, ierr)
       endif

    elseif ( NEST_Filiation( INTERCOMM_ID(HANDLE) ) < 0 ) then !--- daughter
       ileng = PARENT_PRC_nmax(HANDLE) * 4 * 2
       if ( PRC_myrank == PRC_master ) then
          call MPI_IRECV(latlon_catalog, ileng, COMM_datatype, PRC_myrank, tag, INTERCOMM_PARENT, ireq, ierr)
          call MPI_WAIT(ireq, istatus, ierr)
       endif
       call COMM_bcast( latlon_catalog, PARENT_PRC_nmax(HANDLE), 4, 2 )
    else
       if( IO_L ) write(*,*) 'xxx internal error [nest/grid]'
       call PRC_MPIstop
    endif

    return
  end subroutine NEST_COMM_catalogue

  !-----------------------------------------------------------------------------
  !> Check Communication Inter-domains
  subroutine NEST_COMM_ping( &
      HANDLE  )
    use scale_process, only: &
       PRC_myrank,  &
       PRC_master,  &
       PRC_MPIstop
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

    tag        = INTERCOMM_ID(HANDLE) * 100
    ping_error = .false.

    if ( NEST_Filiation( INTERCOMM_ID(HANDLE) ) > 0 ) then !--- parent
       ping = ONLINE_DOMAIN_NUM
       pong = 0

       if ( PRC_myrank == PRC_master ) then
          call MPI_ISEND(ping, 1, MPI_INTEGER, PRC_myrank, tag+1, INTERCOMM_DAUGHTER, ireq1, ierr1)
          call MPI_IRECV(pong, 1, MPI_INTEGER, PRC_myrank, tag+2, INTERCOMM_DAUGHTER, ireq2, ierr2)
          call MPI_WAIT(ireq1, istatus, ierr1)
          call MPI_WAIT(ireq2, istatus, ierr2)
       endif

       call COMM_bcast(pong)

       if ( pong /= INTERCOMM_ID(HANDLE)+1 ) ping_error = .true.

    elseif ( NEST_Filiation( INTERCOMM_ID(HANDLE) ) < 0 ) then !--- daughter
       ping = ONLINE_DOMAIN_NUM
       pong = 0

       if ( PRC_myrank == PRC_master ) then
          call MPI_ISEND(ping, 1, MPI_INTEGER, PRC_myrank, tag+2, INTERCOMM_PARENT, ireq1, ierr1)
          call MPI_IRECV(pong, 1, MPI_INTEGER, PRC_myrank, tag+1, INTERCOMM_PARENT, ireq2, ierr2)
          call MPI_WAIT(ireq1, istatus, ierr1)
          call MPI_WAIT(ireq2, istatus, ierr2)
       endif

       call COMM_bcast(pong)

       if ( pong /= INTERCOMM_ID(HANDLE) ) ping_error = .true.

    else
       if( IO_L ) write(*,*) 'xxx internal error [nest/grid]'
       call PRC_MPIstop
    endif

    if ( ping_error ) then
       if( IO_L ) write(*,*) 'xxx ping destination error [nest/grid]'
       call PRC_MPIstop
    endif

    return
  end subroutine NEST_COMM_ping

  !-----------------------------------------------------------------------------
  !> Inter-domain communication setup for nestdown
  subroutine NEST_COMM_setup_nestdown( &
      HANDLE  )
    use scale_process, only: &
       PRC_myrank,  &
       PRC_master,  &
       PRC_MPIstop
    use scale_grid_real, only: &
       REAL_DOMAIN_CATALOGUE
    use scale_comm, only: &
       COMM_datatype,  &
       COMM_world,     &
       COMM_bcast
    implicit none

    integer, intent(in) :: HANDLE !< id number of nesting relation in this process target

    integer, allocatable :: buffer_LIST(:)
    integer, allocatable :: buffer_ALLLIST(:)

    integer :: ireq, ierr, ileng
    integer :: istatus(MPI_STATUS_SIZE)
    integer :: tag, target_rank

    integer :: i, j, k
    !---------------------------------------------------------------------------

    tag = INTERCOMM_ID(HANDLE) * 100

    if ( NEST_Filiation( INTERCOMM_ID(HANDLE) ) > 0 ) then
    !--------------------------------------------------- parent

       if ( PRC_myrank == PRC_master ) then
          call MPI_IRECV(NEST_TILE_ALLMAX_p, 1, MPI_INTEGER, PRC_myrank, tag+1, INTERCOMM_DAUGHTER, ireq, ierr)
          call MPI_WAIT(ireq, istatus, ierr)
       endif
       call COMM_bcast(NEST_TILE_ALLMAX_p)

       allocate( NEST_TILE_LIST_p (NEST_TILE_ALLMAX_p,DAUGHTER_PRC_nmax(HANDLE)) )
       allocate( NEST_TILE_LIST_YP(NEST_TILE_ALLMAX_p*DAUGHTER_PRC_nmax(HANDLE)) )

       ileng = NEST_TILE_ALLMAX_p*DAUGHTER_PRC_nmax(HANDLE)
       if ( PRC_myrank == PRC_master ) then
          call MPI_IRECV(NEST_TILE_LIST_p, ileng, MPI_INTEGER, PRC_myrank, tag+2, INTERCOMM_DAUGHTER, ireq, ierr)
          call MPI_WAIT(ireq, istatus, ierr)
       endif
       call COMM_bcast(NEST_TILE_LIST_p, NEST_TILE_ALLMAX_p, DAUGHTER_PRC_nmax(HANDLE))

       NEST_TILE_LIST_YP(:) = -1

       k = 0
       do j = 1, DAUGHTER_PRC_nmax(HANDLE)
       do i = 1, NEST_TILE_ALLMAX_p
          if ( NEST_TILE_LIST_p(i,j) == PRC_myrank ) then
             k = k + 1
             NEST_TILE_LIST_YP(k) = j - 1  !rank number is started from 1
          endif
       enddo
       enddo
       NUM_YP = k

       if( IO_L ) write(IO_FID_LOG,'(A,I5,A,I5)') "[P]   Num YP =",NUM_YP,"  Num TILE(MAX) =",NEST_TILE_ALLMAX_p

       if ( PRC_myrank == PRC_master ) then
          call MPI_IRECV(ONLINE_DAUGHTER_USE_VELZ, 1, MPI_LOGICAL, PRC_myrank, tag+3, INTERCOMM_DAUGHTER, ireq, ierr)
          call MPI_WAIT(ireq, istatus, ierr)
       endif
       call COMM_bcast(ONLINE_DAUGHTER_USE_VELZ)

       if( IO_L ) write(IO_FID_LOG,'(1x,A,L2)') '*** NEST: ONLINE_DAUGHTER_USE_VELZ =', ONLINE_DAUGHTER_USE_VELZ

       if ( PRC_myrank == PRC_master ) then
          call MPI_IRECV(ONLINE_DAUGHTER_NO_ROTATE, 1, MPI_LOGICAL, PRC_myrank, tag+4, INTERCOMM_DAUGHTER, ireq, ierr)
          call MPI_WAIT(ireq, istatus, ierr)
       endif
       call COMM_bcast(ONLINE_DAUGHTER_NO_ROTATE)

       if( IO_L ) write(IO_FID_LOG,'(1x,A,L2)') '*** NEST: ONLINE_DAUGHTER_NO_ROTATE =', ONLINE_DAUGHTER_NO_ROTATE

       call NEST_COMM_importgrid_nestdown( HANDLE )

       do i = 1, NUM_YP
          target_rank = NEST_TILE_LIST_YP(i)
          call MPI_ISEND(i, 1, MPI_INTEGER, target_rank, tag+5, INTERCOMM_DAUGHTER, ireq, ierr)
          call MPI_WAIT(ireq, istatus, ierr)
       enddo

       call MPI_BARRIER(INTERCOMM_DAUGHTER, ierr)

    elseif ( NEST_Filiation( INTERCOMM_ID(HANDLE) ) < 0 ) then
    !--------------------------------------------------- daughter

       NEST_TILE_ALL = size( NEST_TILE_ID(:) ) ! should be equal to "NEST_TILE_NUM_X*NEST_TILE_NUM_Y"
       call MPI_Allreduce( NEST_TILE_ALL,      &
                           NEST_TILE_ALLMAX_d, &
                           1,                  &
                           MPI_INTEGER,        &
                           MPI_MAX,            &
                           COMM_world,         &
                           ierr                )
       if( IO_L ) write(IO_FID_LOG,'(A,I5,A,I5)') "[D]   Num YP =",NEST_TILE_ALL,"  Num TILE(MAX) =",NEST_TILE_ALLMAX_d

       if ( PRC_myrank == PRC_master ) then
          call MPI_ISEND(NEST_TILE_ALLMAX_d, 1, MPI_INTEGER, PRC_myrank, tag+1, INTERCOMM_PARENT, ireq, ierr)
          call MPI_WAIT(ireq, istatus, ierr)
       endif

       allocate( buffer_LIST   (NEST_TILE_ALLMAX_d)            )
       allocate( buffer_ALLLIST(NEST_TILE_ALLMAX_d*DAUGHTER_PRC_nmax(HANDLE))   )
       allocate( NEST_TILE_LIST_d(NEST_TILE_ALLMAX_d,DAUGHTER_PRC_nmax(HANDLE)) )

       do i = 1, NEST_TILE_ALLMAX_d
          if ( i <= NEST_TILE_ALL ) then
             buffer_LIST(i) = NEST_TILE_ID(i)
          else
             buffer_LIST(i) = -1
          endif
       enddo

       ileng = NEST_TILE_ALLMAX_d
       call MPI_Allgather( buffer_LIST(:),     &
                           ileng,              &
                           MPI_INTEGER,        &
                           buffer_ALLLIST(:),  &
                           ileng,              &
                           MPI_INTEGER,        &
                           COMM_world,         &
                           ierr                )
       k = 1
       do j = 1, DAUGHTER_PRC_nmax(HANDLE)
       do i = 1, NEST_TILE_ALLMAX_d
          NEST_TILE_LIST_d(i,j) = buffer_ALLLIST(k)
          k = k + 1
       enddo
       enddo

       deallocate( buffer_LIST    )
       deallocate( buffer_ALLLIST )

       ileng = NEST_TILE_ALLMAX_d*DAUGHTER_PRC_nmax(HANDLE)
       if ( PRC_myrank == PRC_master ) then
          call MPI_ISEND(NEST_TILE_LIST_d, ileng, MPI_INTEGER, PRC_myrank, tag+2, INTERCOMM_PARENT, ireq, ierr)
          call MPI_WAIT(ireq, istatus, ierr)
       endif

       if ( PRC_myrank == PRC_master ) then
          call MPI_ISEND(ONLINE_USE_VELZ, 1, MPI_LOGICAL, PRC_myrank, tag+3, INTERCOMM_PARENT, ireq, ierr)
          call MPI_WAIT(ireq, istatus, ierr)
       endif

       if ( PRC_myrank == PRC_master ) then
          call MPI_ISEND(ONLINE_NO_ROTATE, 1, MPI_LOGICAL, PRC_myrank, tag+4, INTERCOMM_PARENT, ireq, ierr)
          call MPI_WAIT(ireq, istatus, ierr)
       endif
       call COMM_bcast(ONLINE_DAUGHTER_NO_ROTATE)

       call NEST_COMM_importgrid_nestdown( HANDLE )

       do i = 1, NEST_TILE_ALL
          target_rank = NEST_TILE_LIST_d(i,PRC_myrank+1)
          call MPI_IRECV( call_order(i), 1, MPI_INTEGER, target_rank, tag+5, INTERCOMM_PARENT, ireq, ierr )
          call MPI_WAIT(ireq, istatus, ierr)
       enddo

       call MPI_BARRIER(INTERCOMM_PARENT, ierr)
    else
    !---------------------------------------------------
       if( IO_L ) write(*,*) 'xxx internal error [nest/grid]'
       call PRC_MPIstop
    endif

    if( NUM_YP * 16 > max_rq .or. NEST_TILE_ALL * 16 > max_rq ) then ! 16 = dyn:5 + qtrc:11
       write(*,*) 'xxx internal error (overflow number of ireq) [nest/grid]'
       write(*,*) '    NUM_YP x 16        = ', NUM_YP * 16
       write(*,*) '    NEST_TILE_ALL x 16 = ', NEST_TILE_ALL * 16
       write(*,*) '    max_rq             = ', max_rq
       call PRC_MPIstop
    endif

    return
  end subroutine NEST_COMM_setup_nestdown

  !-----------------------------------------------------------------------------
  !> Grid Data transfer from parent to daughter: nestdown
  subroutine NEST_COMM_importgrid_nestdown( &
      HANDLE  )
    use scale_process, only: &
       PRC_myrank,  &
       PRC_nmax,    &
       PRC_MPIstop
    use scale_grid_real, only: &
       REAL_LON,    &
       REAL_LAT,    &
       REAL_LONX,   &
       REAL_LONY,   &
       REAL_LATX,   &
       REAL_LATY,   &
       REAL_CZ,     &
       REAL_FZ
    use scale_comm, only: &
       COMM_datatype
    implicit none

    integer, intent(in) :: HANDLE !< id number of nesting relation in this process target

    integer :: ierr, ileng
    integer :: istatus(MPI_STATUS_SIZE)
    integer :: tag, tagbase, target_rank

    integer :: xloc, yloc
    integer :: xs, xe
    integer :: ys, ye

    real(RP) :: max_ref, max_loc

    integer :: i, k, rq
    !---------------------------------------------------------------------------

    tagbase = INTERCOMM_ID(HANDLE) * 100
    rq      = 0

    if ( NEST_Filiation( INTERCOMM_ID(HANDLE) ) > 0 ) then
    !--------------------------------------------------- parent
       do i = 1, NUM_YP
          ! send data to multiple daughter processes
          target_rank = NEST_TILE_LIST_YP(i)

          rq = rq + 1
          ileng = PARENT_IA(HANDLE) * PARENT_JA(HANDLE)
          tag = tagbase + tag_lon
          call MPI_ISEND(REAL_LON, ileng, COMM_datatype, target_rank, tag, INTERCOMM_DAUGHTER, ireq_p(rq), ierr)
          call MPI_WAIT(ireq_p(rq), istatus, ierr)

          rq = rq + 1
          ileng = PARENT_IA(HANDLE) * PARENT_JA(HANDLE)
          tag = tagbase + tag_lat
          call MPI_ISEND(REAL_LAT, ileng, COMM_datatype, target_rank, tag, INTERCOMM_DAUGHTER, ireq_p(rq), ierr)
          call MPI_WAIT(ireq_p(rq), istatus, ierr)

          rq = rq + 1
          ileng = PARENT_IA(HANDLE) * PARENT_JA(HANDLE)
          tag = tagbase + tag_lonx
          call MPI_ISEND(REAL_LONX, ileng, COMM_datatype, target_rank, tag, INTERCOMM_DAUGHTER, ireq_p(rq), ierr)
          call MPI_WAIT(ireq_p(rq), istatus, ierr)

          rq = rq + 1
          ileng = PARENT_IA(HANDLE) * PARENT_JA(HANDLE)
          tag = tagbase + tag_latx
          call MPI_ISEND(REAL_LATX, ileng, COMM_datatype, target_rank, tag, INTERCOMM_DAUGHTER, ireq_p(rq), ierr)
          call MPI_WAIT(ireq_p(rq), istatus, ierr)

          rq = rq + 1
          ileng = PARENT_IA(HANDLE) * PARENT_JA(HANDLE)
          tag = tagbase + tag_lony
          call MPI_ISEND(REAL_LONY, ileng, COMM_datatype, target_rank, tag, INTERCOMM_DAUGHTER, ireq_p(rq), ierr)
          call MPI_WAIT(ireq_p(rq), istatus, ierr)

          rq = rq + 1
          ileng = PARENT_IA(HANDLE) * PARENT_JA(HANDLE)
          tag = tagbase + tag_laty
          call MPI_ISEND(REAL_LATY, ileng, COMM_datatype, target_rank, tag, INTERCOMM_DAUGHTER, ireq_p(rq), ierr)
          call MPI_WAIT(ireq_p(rq), istatus, ierr)

          rq = rq + 1
          ileng = PARENT_KA(HANDLE) * PARENT_IA(HANDLE) * PARENT_JA(HANDLE)
          tag = tagbase + tag_cz
          call MPI_ISEND(REAL_CZ, ileng, COMM_datatype, target_rank, tag, INTERCOMM_DAUGHTER, ireq_p(rq), ierr)
          call MPI_WAIT(ireq_p(rq), istatus, ierr)

          rq = rq + 1
          ileng = (PARENT_KA(HANDLE)+1) * PARENT_IA(HANDLE) * PARENT_JA(HANDLE)
          tag = tagbase + tag_fz
          call MPI_ISEND(REAL_FZ, ileng, COMM_datatype, target_rank, tag, INTERCOMM_DAUGHTER, ireq_p(rq), ierr)
          call MPI_WAIT(ireq_p(rq), istatus, ierr)
       enddo


    elseif ( NEST_Filiation( INTERCOMM_ID(HANDLE) ) < 0 ) then
    !--------------------------------------------------- daughter
       do i = 1, NEST_TILE_ALL
          ! receive data from multiple parent tiles
          target_rank = NEST_TILE_LIST_d(i,PRC_myrank+1)

          xloc = mod( i-1, NEST_TILE_NUM_X ) + 1
          yloc = int( real(i-1) / real(NEST_TILE_NUM_X) ) + 1

          xs = PARENT_IMAX(HANDLE) * (xloc-1) + 1
          xe = PARENT_IMAX(HANDLE) * xloc
          ys = PARENT_JMAX(HANDLE) * (yloc-1) + 1
          ye = PARENT_JMAX(HANDLE) * yloc

          rq = rq + 1
          ileng = PARENT_IA(HANDLE) * PARENT_JA(HANDLE)
          tag = tagbase + tag_lon
          call MPI_IRECV(buffer_2D, ileng, COMM_datatype, target_rank, tag, INTERCOMM_PARENT, ireq_d(rq), ierr)
          call MPI_WAIT(ireq_d(rq), istatus, ierr)
          buffer_ref_LON(xs:xe,ys:ye)  = buffer_2D(PRNT_IS(HANDLE):PRNT_IE(HANDLE),PRNT_JS(HANDLE):PRNT_JE(HANDLE))

          rq = rq + 1
          ileng = PARENT_IA(HANDLE) * PARENT_JA(HANDLE)
          tag = tagbase + tag_lat
          call MPI_IRECV(buffer_2D, ileng, COMM_datatype, target_rank, tag, INTERCOMM_PARENT, ireq_d(rq), ierr)
          call MPI_WAIT(ireq_d(rq), istatus, ierr)
          buffer_ref_LAT(xs:xe,ys:ye)  = buffer_2D(PRNT_IS(HANDLE):PRNT_IE(HANDLE),PRNT_JS(HANDLE):PRNT_JE(HANDLE))

          rq = rq + 1
          ileng = PARENT_IA(HANDLE) * PARENT_JA(HANDLE)
          tag = tagbase + tag_lonx
          call MPI_IRECV(buffer_2D, ileng, COMM_datatype, target_rank, tag, INTERCOMM_PARENT, ireq_d(rq), ierr)
          call MPI_WAIT(ireq_d(rq), istatus, ierr)
          buffer_ref_LONX(xs:xe,ys:ye)  = buffer_2D(PRNT_IS(HANDLE):PRNT_IE(HANDLE),PRNT_JS(HANDLE):PRNT_JE(HANDLE))

          rq = rq + 1
          ileng = PARENT_IA(HANDLE) * PARENT_JA(HANDLE)
          tag = tagbase + tag_latx
          call MPI_IRECV(buffer_2D, ileng, COMM_datatype, target_rank, tag, INTERCOMM_PARENT, ireq_d(rq), ierr)
          call MPI_WAIT(ireq_d(rq), istatus, ierr)
          buffer_ref_LATX(xs:xe,ys:ye)  = buffer_2D(PRNT_IS(HANDLE):PRNT_IE(HANDLE),PRNT_JS(HANDLE):PRNT_JE(HANDLE))

          rq = rq + 1
          ileng = PARENT_IA(HANDLE) * PARENT_JA(HANDLE)
          tag = tagbase + tag_lony
          call MPI_IRECV(buffer_2D, ileng, COMM_datatype, target_rank, tag, INTERCOMM_PARENT, ireq_d(rq), ierr)
          call MPI_WAIT(ireq_d(rq), istatus, ierr)
          buffer_ref_LONY(xs:xe,ys:ye)  = buffer_2D(PRNT_IS(HANDLE):PRNT_IE(HANDLE),PRNT_JS(HANDLE):PRNT_JE(HANDLE))

          rq = rq + 1
          ileng = PARENT_IA(HANDLE) * PARENT_JA(HANDLE)
          tag = tagbase + tag_laty
          call MPI_IRECV(buffer_2D, ileng, COMM_datatype, target_rank, tag, INTERCOMM_PARENT, ireq_d(rq), ierr)
          call MPI_WAIT(ireq_d(rq), istatus, ierr)
          buffer_ref_LATY(xs:xe,ys:ye)  = buffer_2D(PRNT_IS(HANDLE):PRNT_IE(HANDLE),PRNT_JS(HANDLE):PRNT_JE(HANDLE))

          rq = rq + 1
          ileng = PARENT_KA(HANDLE) * PARENT_IA(HANDLE) * PARENT_JA(HANDLE)
          tag = tagbase + tag_cz
          call MPI_IRECV(buffer_3D, ileng, COMM_datatype, target_rank, tag, INTERCOMM_PARENT, ireq_d(rq), ierr)
          call MPI_WAIT(ireq_d(rq), istatus, ierr)
          do k = 1, PARENT_KA(HANDLE)
             buffer_ref_CZ(k,xs:xe,ys:ye)  = buffer_3D(k,PRNT_IS(HANDLE):PRNT_IE(HANDLE),PRNT_JS(HANDLE):PRNT_JE(HANDLE))
          enddo

          rq = rq + 1
          ileng = (PARENT_KA(HANDLE)+1) * PARENT_IA(HANDLE) * PARENT_JA(HANDLE)
          tag = tagbase + tag_fz
          call MPI_IRECV(buffer_3DF,ileng, COMM_datatype, target_rank, tag, INTERCOMM_PARENT, ireq_d(rq), ierr)
          call MPI_WAIT(ireq_d(rq), istatus, ierr)
          do k = 0, PARENT_KA(HANDLE)
             buffer_ref_FZ(k,xs:xe,ys:ye)  = buffer_3DF(k,PRNT_IS(HANDLE):PRNT_IE(HANDLE),PRNT_JS(HANDLE):PRNT_JE(HANDLE))
          enddo
       enddo

       ! check domain compatibility
       max_ref = maxval( buffer_ref_FZ(:,:,:) )
       max_loc = maxval( REAL_FZ(KS-1:KE,:,:) ) ! HALO + 1
       if ( max_ref < max_loc ) then
          write(*,*) 'xxx ERROR: REQUESTED DOMAIN IS TOO MUCH BROAD'
          write(*,*) 'xxx -- VERTICAL direction over the limit'
          write(*,*) 'xxx -- reference max: ', max_ref
          write(*,*) 'xxx --     local max: ', max_loc
          call PRC_MPIstop
       endif

    else
    !---------------------------------------------------
       write(*,*) 'xxx internal error [nest/grid]'
       call PRC_MPIstop
    endif

    return
  end subroutine NEST_COMM_importgrid_nestdown

  !-----------------------------------------------------------------------------
  !> Boundary data transfer from parent to daughter: nestdown
  subroutine NEST_COMM_nestdown( &
      HANDLE,              & ! [in   ]
      BND_QA,              & ! [in   ]
      ipt_DENS,            & ! [in   ]
      ipt_MOMZ,            & ! [in   ]
      ipt_MOMX,            & ! [in   ]
      ipt_MOMY,            & ! [in   ]
      ipt_RHOT,            & ! [in   ]
      ipt_QTRC,            & ! [in   ]
      interped_ref_DENS,   & ! [inout]
      interped_ref_VELZ,   & ! [inout]
      interped_ref_VELX,   & ! [inout]
      interped_ref_VELY,   & ! [inout]
      interped_ref_POTT,   & ! [inout]
      interped_ref_QTRC    ) ! [inout]
    use scale_process, only: &
       PRC_myrank,  &
       PRC_nmax,    &
       PRC_MPIstop
    use scale_grid_real, only: &
       REAL_DOMAIN_CATALOGUE
    use scale_comm, only: &
       COMM_vars8,  &
       COMM_wait,   &
       COMM_world,  &
       COMM_datatype
    use scale_gridtrans, only: &
       rotc => GTRANS_ROTC
    implicit none

    integer,  intent(in)    :: HANDLE        !< id number of nesting relation in this process target
    integer,  intent(in)    :: BND_QA        !< num of tracer

    real(RP), intent(in   ) :: ipt_DENS(PARENT_KA(HANDLE),PARENT_IA(HANDLE),PARENT_JA(HANDLE))
    real(RP), intent(in   ) :: ipt_MOMZ(PARENT_KA(HANDLE),PARENT_IA(HANDLE),PARENT_JA(HANDLE))
    real(RP), intent(in   ) :: ipt_MOMX(PARENT_KA(HANDLE),PARENT_IA(HANDLE),PARENT_JA(HANDLE))
    real(RP), intent(in   ) :: ipt_MOMY(PARENT_KA(HANDLE),PARENT_IA(HANDLE),PARENT_JA(HANDLE))
    real(RP), intent(in   ) :: ipt_RHOT(PARENT_KA(HANDLE),PARENT_IA(HANDLE),PARENT_JA(HANDLE))
    real(RP), intent(in   ) :: ipt_QTRC(PARENT_KA(HANDLE),PARENT_IA(HANDLE),PARENT_JA(HANDLE),BND_QA)

    real(RP), intent(inout) :: interped_ref_DENS(DAUGHTER_KA(HANDLE),DAUGHTER_IA(HANDLE),DAUGHTER_JA(HANDLE))
    real(RP), intent(inout) :: interped_ref_VELZ(DAUGHTER_KA(HANDLE),DAUGHTER_IA(HANDLE),DAUGHTER_JA(HANDLE))
    real(RP), intent(inout) :: interped_ref_VELX(DAUGHTER_KA(HANDLE),DAUGHTER_IA(HANDLE),DAUGHTER_JA(HANDLE))
    real(RP), intent(inout) :: interped_ref_VELY(DAUGHTER_KA(HANDLE),DAUGHTER_IA(HANDLE),DAUGHTER_JA(HANDLE))
    real(RP), intent(inout) :: interped_ref_POTT(DAUGHTER_KA(HANDLE),DAUGHTER_IA(HANDLE),DAUGHTER_JA(HANDLE))
    real(RP), intent(inout) :: interped_ref_QTRC(DAUGHTER_KA(HANDLE),DAUGHTER_IA(HANDLE),DAUGHTER_JA(HANDLE),BND_QA)

    real(RP) :: dummy(1,1,1)
    real(RP) :: dens (DAUGHTER_KA(HANDLE),DAUGHTER_IA(HANDLE),DAUGHTER_JA(HANDLE))
    real(RP) :: u_lld(DAUGHTER_KA(HANDLE),DAUGHTER_IA(HANDLE),DAUGHTER_JA(HANDLE))
    real(RP) :: v_lld(DAUGHTER_KA(HANDLE),DAUGHTER_IA(HANDLE),DAUGHTER_JA(HANDLE))
    real(RP) :: work1(DAUGHTER_KA(HANDLE),DAUGHTER_IA(HANDLE),DAUGHTER_JA(HANDLE))

    real(RP) :: u_on_map, v_on_map
    integer :: ierr
    integer :: tagbase, tagcomm
    integer :: isu_tag, isu_tagf
    integer :: i, j, k, iq

    integer, parameter :: cosin = 1
    integer, parameter :: sine  = 2
    !---------------------------------------------------------------------------

    if ( BND_QA > I_BNDQA ) then
       if( IO_L ) write(*,*) 'xxx internal error: BND_QA is larger than I_BNDQA [nest/grid]'
       call PRC_MPIstop
    endif
    if ( BND_QA > max_bndqa ) then
       if( IO_L ) write(*,*) 'xxx internal error: BND_QA is larger than max_bndqa [nest/grid]'
       call PRC_MPIstop
    endif

    tagcomm = INTERCOMM_ID(HANDLE) * order_tag_comm

    if ( NEST_Filiation( INTERCOMM_ID(HANDLE) ) > 0 ) then
    !-------------------------------------------------------- parent [send issue]
       call PROF_rapstart('NESTCOM total parent')

       nsend = nsend + 1
       if( IO_L ) write(IO_FID_LOG,'(1X,A,I5,A)') "*** NestIDC [P]: que send ( ", nsend, " )"

       ! to keep values at that time by finish of sending process
       org_DENS(:,:,:) = ipt_DENS(:,:,:)
       org_MOMZ(:,:,:) = ipt_MOMZ(:,:,:)
       org_MOMX(:,:,:) = ipt_MOMX(:,:,:)
       org_MOMY(:,:,:) = ipt_MOMY(:,:,:)
       org_RHOT(:,:,:) = ipt_RHOT(:,:,:)
       do iq = 1, BND_QA
          org_QTRC(:,:,:,iq) = ipt_QTRC(:,:,:,iq)
       enddo


       !*** request control
       !--- do not change the calling order below;
       !--- it should be consistent with the order in "NEST_COMM_recvwait_issue"
       rq_ctl_p = 0

       if ( .not. ONLINE_DAUGHTER_NO_ROTATE ) then
          do j = PRNT_JS(HANDLE), PRNT_JE(HANDLE)
          do i = PRNT_IS(HANDLE), PRNT_IE(HANDLE)
          do k = PRNT_KS(HANDLE), PRNT_KE(HANDLE)
             u_on_map = org_MOMX(k,i,j) / ( org_DENS(k,i+1,j) + org_DENS(k,i,j) ) * 2.0_RP
             v_on_map = org_MOMY(k,i,j) / ( org_DENS(k,i,j+1) + org_DENS(k,i,j) ) * 2.0_RP

             u_llp(k,i,j) = u_on_map * rotc(i,j,cosin) - v_on_map * rotc(i,j,sine )
             v_llp(k,i,j) = u_on_map * rotc(i,j,sine ) + v_on_map * rotc(i,j,cosin)
          enddo
          enddo
          enddo
       end if

       tagbase = tagcomm + tag_dens*order_tag_var
       call NEST_COMM_intercomm_nestdown( org_DENS, dummy, tagbase, I_SCLR, HANDLE, isu_tag, isu_tagf, .true. )

       tagbase = tagcomm + tag_momz*order_tag_var
       if ( ONLINE_DAUGHTER_USE_VELZ ) then
          call NEST_COMM_intercomm_nestdown( org_MOMZ, dummy, tagbase, I_ZSTG, HANDLE, isu_tag, isu_tagf )
       end if

       tagbase = tagcomm + tag_momx*order_tag_var
       if ( ONLINE_DAUGHTER_NO_ROTATE ) then
          call NEST_COMM_intercomm_nestdown( org_MOMX, dummy, tagbase, I_XSTG, HANDLE, isu_tag, isu_tagf )
       else
          call NEST_COMM_intercomm_nestdown( u_llp, dummy, tagbase, I_SCLR, HANDLE, isu_tag, isu_tagf )
       end if

       tagbase = tagcomm + tag_momy*order_tag_var
       if ( ONLINE_DAUGHTER_NO_ROTATE ) then
          call NEST_COMM_intercomm_nestdown( org_MOMY, dummy, tagbase, I_YSTG, HANDLE, isu_tag, isu_tagf )
       else
          call NEST_COMM_intercomm_nestdown( v_llp, dummy, tagbase, I_SCLR, HANDLE, isu_tag, isu_tagf )
       end if

       tagbase = tagcomm + tag_rhot*order_tag_var
       call NEST_COMM_intercomm_nestdown( org_RHOT, dummy, tagbase, I_SCLR, HANDLE, isu_tag, isu_tagf )

       do iq = 1, BND_QA
          tagbase = tagcomm + (tag_qx*10+iq)*order_tag_var
          call NEST_COMM_intercomm_nestdown( org_QTRC(:,:,:,iq), dummy, tagbase, I_SCLR, HANDLE, isu_tag, isu_tagf )
       enddo

       rq_tot_p = rq_ctl_p

       call PROF_rapend('NESTCOM total parent')

    elseif ( NEST_Filiation( INTERCOMM_ID(HANDLE) ) < 0 ) then
    !-------------------------------------------------------- daughter [wait issue]
       call PROF_rapstart('NESTCOM total child')

       nwait_d = nwait_d + 1
       !if( IO_L ) write(IO_FID_LOG,'(1X,A,I5,A)') "*** NestIDC [C]: que wait ( ", nwait_d, " )"

       !*** reset issue tag and request control
       !--- do not change the calling order below;
       !--- it should be consistent with the order in "NEST_COMM_recvwait_issue"
       isu_tag  = 0
       isu_tagf = 0

       call PROF_rapstart('NESTCOM wait child')
       call NEST_COMM_waitall( rq_tot_d, ireq_d )
       if ( ONLINE_AGGRESSIVE_COMM ) then
          ! nothing to do
       else
          call MPI_BARRIER(INTERCOMM_PARENT, ierr)
       endif
       call PROF_rapend  ('NESTCOM wait child')

       tagbase = tagcomm + tag_dens*order_tag_var
       call NEST_COMM_intercomm_nestdown( dummy, dens, tagbase, I_SCLR, HANDLE, isu_tag, isu_tagf, .true. )
       do j = 1, DAUGHTER_JA(HANDLE)
       do i = 1, DAUGHTER_IA(HANDLE)
       do k = DATR_KS(HANDLE), DATR_KE(HANDLE)
          interped_ref_DENS(k,i,j) = dens(k,i,j)
       enddo
       enddo
       enddo
       call COMM_vars8( interped_ref_DENS, 1 )
       call COMM_wait ( interped_ref_DENS, 1, .false. )

       tagbase = tagcomm + tag_momz*order_tag_var
       if ( ONLINE_USE_VELZ ) then
          call NEST_COMM_intercomm_nestdown( dummy, work1, tagbase, I_ZSTG, HANDLE, isu_tag, isu_tagf )
          do j = 1, DAUGHTER_JA(HANDLE)
          do i = 1, DAUGHTER_IA(HANDLE)
          do k = DATR_KS(HANDLE), DATR_KE(HANDLE)-1
             interped_ref_VELZ(k,i,j) = work1(k,i,j) / ( dens(k,i,j) + dens(k+1,i,j) ) * 2.0_RP
          enddo
          enddo
          enddo
       end if

       tagbase = tagcomm + tag_momx*order_tag_var
       if ( ONLINE_NO_ROTATE ) then
          call NEST_COMM_intercomm_nestdown( dummy, u_lld, tagbase, I_XSTG, HANDLE, isu_tag, isu_tagf )
       else
          call NEST_COMM_intercomm_nestdown( dummy, u_lld, tagbase, I_SCLR, HANDLE, isu_tag, isu_tagf )
       endif
       tagbase = tagcomm + tag_momy*order_tag_var
       if ( ONLINE_NO_ROTATE ) then
          call NEST_COMM_intercomm_nestdown( dummy, v_lld, tagbase, I_YSTG, HANDLE, isu_tag, isu_tagf )
       else
          call NEST_COMM_intercomm_nestdown( dummy, v_lld, tagbase, I_SCLR, HANDLE, isu_tag, isu_tagf )
       endif

       if ( ONLINE_NO_ROTATE ) then
          do j = 1, DAUGHTER_JA(HANDLE)
          do i = 1, DAUGHTER_IA(HANDLE)-1
          do k = DATR_KS(HANDLE), DATR_KE(HANDLE)
             interped_ref_VELX(k,i,j) = u_lld(k,i,j) &
                                      / ( interped_ref_DENS(k,i+1,j) + interped_ref_DENS(k,i,j) ) * 2.0_RP
          enddo
          enddo
          enddo
          do j = 1, DAUGHTER_JA(HANDLE)
          do k = DATR_KS(HANDLE), DATR_KE(HANDLE)
             interped_ref_VELX(k,DAUGHTER_IA(HANDLE),j) = u_lld(k,DAUGHTER_IA(HANDLE),j) &
                                                        / interped_ref_DENS(k,DAUGHTER_IA(HANDLE),j)
          enddo
          enddo
          call COMM_vars8( interped_ref_VELX, 2 )
          do j = 1, DAUGHTER_JA(HANDLE)-1
          do i = 1, DAUGHTER_IA(HANDLE)
          do k = DATR_KS(HANDLE), DATR_KE(HANDLE)
             interped_ref_VELY(k,i,j) = v_lld(k,i,j) &
                                      / ( interped_ref_DENS(k,i,j+1) + interped_ref_DENS(k,i,j) ) * 2.0_RP
          enddo
          enddo
          enddo
          do i = 1, DAUGHTER_IA(HANDLE)
          do k = DATR_KS(HANDLE), DATR_KE(HANDLE)
             interped_ref_VELY(k,i,DAUGHTER_JA(HANDLE)) = v_lld(k,i,DAUGHTER_JA(HANDLE)) &
                                                        / interped_ref_DENS(k,i,DAUGHTER_JA(HANDLE))
          enddo
          enddo
          call COMM_vars8( interped_ref_VELY, 3 )
          call COMM_wait ( interped_ref_VELX, 2, .false. )
          call COMM_wait ( interped_ref_VELY, 3, .false. )
       else ! rotate
          do j = 1, DAUGHTER_JA(HANDLE)
          do i = 1, DAUGHTER_IA(HANDLE)
          do k = DATR_KS(HANDLE), DATR_KE(HANDLE)
             interped_ref_VELX(k,i,j) =   u_lld(k,i,j) * rotc(i,j,cosin) + v_lld(k,i,j) * rotc(i,j,sine )
             interped_ref_VELY(k,i,j) = - u_lld(k,i,j) * rotc(i,j,sine ) + v_lld(k,i,j) * rotc(i,j,cosin)
          enddo
          enddo
          enddo
       end if

       tagbase = tagcomm + tag_rhot*order_tag_var
       call NEST_COMM_intercomm_nestdown( dummy, work1, tagbase, I_SCLR, HANDLE, isu_tag, isu_tagf )
       do j = 1, DAUGHTER_JA(HANDLE)
       do i = 1, DAUGHTER_IA(HANDLE)
       do k = DATR_KS(HANDLE), DATR_KE(HANDLE)
          interped_ref_POTT(k,i,j) = work1(k,i,j) / interped_ref_DENS(k,i,j)
       enddo
       enddo
       enddo

       do iq = 1, BND_QA
          tagbase = tagcomm + (tag_qx*10+iq)*order_tag_var
          call NEST_COMM_intercomm_nestdown( dummy, work1, tagbase, I_SCLR, HANDLE, isu_tag, isu_tagf )
          do j = 1, DAUGHTER_JA(HANDLE)
          do i = 1, DAUGHTER_IA(HANDLE)
          do k = DATR_KS(HANDLE), DATR_KE(HANDLE)
             interped_ref_QTRC(k,i,j,iq) = work1(k,i,j)
          enddo
          enddo
          enddo
       enddo

       call PROF_rapend('NESTCOM total child')
    else
       write(*,*) 'xxx internal error [nestdown: nest/grid]'
       call PRC_MPIstop
    endif

    return
  end subroutine NEST_COMM_nestdown

  !-----------------------------------------------------------------------------
  !> Sub-command for data transfer from parent to daughter: nestdown
  subroutine NEST_COMM_recvwait_issue( &
      HANDLE,     & ! [in]
      BND_QA      ) ! [in]
    use scale_process, only: &
       PRC_MPIstop
    use scale_comm, only: &
       COMM_world
    implicit none

    integer, intent(in) :: HANDLE    !< id number of nesting relation in this process target
    integer, intent(in) :: BND_QA    !< num of tracer in online-nesting

    integer :: isu_tag, isu_tagf
    integer :: tagbase, tagcomm
    integer :: ierr
    integer :: iq
    !---------------------------------------------------------------------------

    if ( BND_QA > I_BNDQA ) then
       write(*,*) 'xxx internal error: about BND_QA [nest/grid]'
       call PRC_MPIstop
    endif

    tagcomm = INTERCOMM_ID(HANDLE) * order_tag_comm

    if ( NEST_Filiation( INTERCOMM_ID(HANDLE) ) > 0 ) then
    !-------------------------------------------------------- parent [wait issue]
       call PROF_rapstart('NESTCOM total parent')
       nwait_p = nwait_p + 1
       !if( IO_L ) write(IO_FID_LOG,'(1X,A,I5,A)') "*** NestIDC [P]: que wait ( ", nwait_p, " )"

       call PROF_rapstart('NESTCOM wait parent')
       call NEST_COMM_issuer_of_wait( HANDLE )

       if ( ONLINE_AGGRESSIVE_COMM ) then
          ! nothing to do
       else
          call MPI_BARRIER(INTERCOMM_DAUGHTER, ierr)
       endif
       call PROF_rapend  ('NESTCOM wait parent')

       call PROF_rapend('NESTCOM total parent')
    elseif ( NEST_Filiation( INTERCOMM_ID(HANDLE) ) < 0 ) then
    !-------------------------------------------------------- daughter [receive issue]
       call PROF_rapstart('NESTCOM total child')
       nrecv = nrecv + 1
       if( IO_L ) write(IO_FID_LOG,'(1X,A,I5,A)') "*** NestIDC [C]: que recv ( ", nrecv, " )"

       !*** reset issue tag and request control
       !--- do not change the calling order below;
       !--- it should be consistent with the order in "NEST_COMM_nestdown"
       isu_tag  = 0
       isu_tagf = 0
       rq_ctl_d = 0

       tagbase = tagcomm + tag_dens*order_tag_var
       call NEST_COMM_issuer_of_receive( tagbase, I_SCLR, HANDLE, isu_tag, isu_tagf, .true. )

       tagbase = tagcomm + tag_momz*order_tag_var
       if ( ONLINE_USE_VELZ ) then
          call NEST_COMM_issuer_of_receive( tagbase, I_ZSTG, HANDLE, isu_tag, isu_tagf )
       end if

       tagbase = tagcomm + tag_momx*order_tag_var
       if ( ONLINE_NO_ROTATE ) then
          call NEST_COMM_issuer_of_receive( tagbase, I_XSTG, HANDLE, isu_tag, isu_tagf )
       else
          call NEST_COMM_issuer_of_receive( tagbase, I_SCLR, HANDLE, isu_tag, isu_tagf )
       endif

       tagbase = tagcomm + tag_momy*order_tag_var
       if ( ONLINE_NO_ROTATE ) then
          call NEST_COMM_issuer_of_receive( tagbase, I_YSTG, HANDLE, isu_tag, isu_tagf )
       else
          call NEST_COMM_issuer_of_receive( tagbase, I_SCLR, HANDLE, isu_tag, isu_tagf )
       endif

       tagbase = tagcomm + tag_rhot*order_tag_var
       call NEST_COMM_issuer_of_receive( tagbase, I_SCLR, HANDLE, isu_tag, isu_tagf )

       do iq = 1, BND_QA
          tagbase = tagcomm + (tag_qx*10+iq)*order_tag_var
          call NEST_COMM_issuer_of_receive( tagbase, I_SCLR, HANDLE, isu_tag, isu_tagf )
       enddo

       rq_tot_d = rq_ctl_d

       call PROF_rapend('NESTCOM total child')
    else
       write(*,*) 'xxx internal error [issue: nest/grid]'
       call PRC_MPIstop
    endif

    return
  end subroutine NEST_COMM_recvwait_issue

  !-----------------------------------------------------------------------------
  !> Sub-command for data transfer from parent to daughter: nestdown
  subroutine NEST_COMM_recv_cancel( &
      HANDLE     ) ! [in]
    use scale_process, only: &
       PRC_myrank,  &
       PRC_nmax,    &
       PRC_MPIstop
    implicit none

    integer, intent(in) :: HANDLE   !< id number of nesting relation in this process target

    integer :: ierr
    integer :: istatus(MPI_STATUS_SIZE)

    integer :: rq
    logical :: flag
    !---------------------------------------------------------------------------

    if ( NEST_Filiation( INTERCOMM_ID(HANDLE) ) > 0 ) then
    !-------------------------------------------------------- parent
       !--- Nothing to do

    elseif ( NEST_Filiation( INTERCOMM_ID(HANDLE) ) < 0 ) then
    !-------------------------------------------------------- daughter [receive issue]
       if( IO_L ) write(IO_FID_LOG,'(1X,A,I5,A)') "*** NestIDC [C]: CANCEL recv ( ", nrecv, " )"
       do rq = 1, rq_tot_d
          if ( ireq_d(rq) .ne. MPI_REQUEST_NULL ) then
             call MPI_CANCEL(ireq_d(rq), ierr)
             call MPI_TEST_CANCELLED(istatus, flag, ierr)
             !if ( .NOT. flag ) then
             !   write(IO_FID_LOG,*) 'XXX ERROR: receive actions do not cancelled, req = ', rq
             !endif
          endif
       enddo

    else
       write(*,*) 'xxx internal error [cancel: nest/grid]'
       call PRC_MPIstop
    endif

    return
  end subroutine NEST_COMM_recv_cancel

  !-----------------------------------------------------------------------------
  !> Inter-communication from parent to daughter: nestdown
  subroutine NEST_COMM_intercomm_nestdown_3D( &
      pvar,      & ! [in ]
      dvar,      & ! [out]
      tagbase,   & ! [in ]
      id_stag,   & ! [in ]
      HANDLE,    & ! [in ]
      isu_tag,   & ! [inout]
      isu_tagf,  & ! [inout]
      flag_dens  ) ! [in ]: optional
    use scale_process, only: &
       PRC_myrank,  &
       PRC_nmax,    &
       PRC_MPIstop
    use scale_comm, only: &
       COMM_datatype
    implicit none

    real(RP), intent(in)    :: pvar(:,:,:)        !< variable from parent domain (PARENT_KA,PARENT_IA,PARENT_JA / 1,1,1)
    real(RP), intent(out)   :: dvar(:,:,:)        !< variable to daughter domain (1,1,1 / MY_KA,MY_IA,MY_JA)
    integer , intent(in)    :: tagbase            !< communication tag of the variable
    integer , intent(in)    :: id_stag            !< id of staggered grid option
    integer , intent(in)    :: HANDLE             !< id number of nesting relation in this process target
    integer , intent(inout) :: isu_tag            !< tag for receive buffer
    integer , intent(inout) :: isu_tagf           !< tag for receive buffer

    logical , intent(in), optional  :: flag_dens  !< flag of logarithmic interpolation for density

    integer :: ierr, ileng
    integer :: tag, target_rank

    integer :: xloc, yloc
    integer :: xs, xe
    integer :: ys, ye

    integer :: i, j, k
    integer :: ig, rq, yp
    logical :: no_zstag    = .true.
    logical :: logarithmic = .false.
    !---------------------------------------------------------------------------
    logarithmic = .false.
    if ( present(flag_dens) ) then
    if ( flag_dens ) then
       logarithmic = .true.
    endif
    endif

    if ( id_stag == I_SCLR )     then
       no_zstag = .true.
       ig       = I_SCLR
    elseif ( id_stag == I_ZSTG ) then
       no_zstag = .false.
       ig       = I_ZSTG
    elseif ( id_stag == I_XSTG ) then
       no_zstag = .true.
       ig       = I_XSTG
    elseif ( id_stag == I_YSTG ) then
       no_zstag = .true.
       ig       = I_YSTG
    endif

    if ( no_zstag ) then
       ileng = PARENT_KA(HANDLE) * PARENT_IA(HANDLE) * PARENT_JA(HANDLE)
    else
       ileng = (PARENT_KA(HANDLE)+1) * PARENT_IA(HANDLE) * PARENT_JA(HANDLE)
    endif

    if ( NEST_Filiation( INTERCOMM_ID(HANDLE) ) > 0 ) then
    !--------------------------------------------------- parent
       rq = rq_ctl_p

       do yp = 1, NUM_YP
          rq = rq + 1

          call PROF_rapstart('NESTCOM send parent')
          ! send data to multiple daughter processes
          target_rank = NEST_TILE_LIST_YP(yp)
          tag = tagbase  + yp
          call MPI_ISEND(pvar, ileng, COMM_datatype, target_rank, tag, INTERCOMM_DAUGHTER, ireq_p(rq), ierr)
          call PROF_rapend('NESTCOM send parent')

          dvar(:,:,:) = -1.0_RP  ! input as a dummy value
       enddo

       rq_ctl_p = rq

    elseif ( NEST_Filiation( INTERCOMM_ID(HANDLE) ) < 0 ) then
    !--------------------------------------------------- daughter
       rq = rq_ctl_d

       do yp = 1, NEST_TILE_ALL ! YP Loop
          rq = rq + 1

          xloc = mod( yp-1, NEST_TILE_NUM_X ) + 1
          yloc = int( real(yp-1) / real(NEST_TILE_NUM_X) ) + 1

          xs = PARENT_IMAX(HANDLE) * (xloc-1) + 1
          xe = PARENT_IMAX(HANDLE) * xloc
          ys = PARENT_JMAX(HANDLE) * (yloc-1) + 1
          ye = PARENT_JMAX(HANDLE) * yloc

          if ( no_zstag ) then
             isu_tag = isu_tag + 1

             if ( .not. logarithmic ) then
                ! linear interpolation
                do k = 1, PARENT_KA(HANDLE)
                   buffer_ref_3D(k,xs:xe,ys:ye) &
                   = recvbuf_3D(k,PRNT_IS(HANDLE):PRNT_IE(HANDLE),PRNT_JS(HANDLE):PRNT_JE(HANDLE),isu_tag)
                enddo
             else
                ! logarithmic weighted interpolation
                do k = 1, PARENT_KA(HANDLE)
                   buffer_ref_3D(k,xs:xe,ys:ye) &
                   = log( recvbuf_3D(k,PRNT_IS(HANDLE):PRNT_IE(HANDLE),PRNT_JS(HANDLE):PRNT_JE(HANDLE),isu_tag) )
                enddo
             endif
          else
             isu_tagf = isu_tagf + 1

             do k = 0, PARENT_KA(HANDLE)
                buffer_ref_3DF(k,xs:xe,ys:ye) &
                = recvbuf_3DF(k,PRNT_IS(HANDLE):PRNT_IE(HANDLE),PRNT_JS(HANDLE):PRNT_JE(HANDLE),isu_tagf)
             enddo
          endif

          if ( isu_tag > max_isu .or. isu_tagf > max_isuf ) then
             write(*,*) 'xxx Exceeded maximum issue [intercomm: nest/grid]'
             write(*,*) 'xxx isu_tag  = ', isu_tag
             write(*,*) 'xxx isu_tagf = ', isu_tagf
             call PRC_MPIstop
          endif

       enddo ! YP Loop
       rq_ctl_d = rq

       call PROF_rapstart('NESTCOM interp')
       dvar(:,:,:) = 0.0_RP 

       if ( no_zstag ) then
!          if ( .not. logarithmic ) then
!             ! linear interpolation
!             do j = 1, DAUGHTER_JA(HANDLE)
!             do i = 1, DAUGHTER_IA(HANDLE)
!             do k = DATR_KS(HANDLE), DATR_KE(HANDLE)
!                dvar(k,i,j) = buffer_ref_3D(kgrd(k,i,j,1,1,ig),igrd(i,j,1,ig),jgrd(i,j,1,ig)) &
!                            * hfact(i,j,1,ig) * vfact(k,i,j,1,1,ig)                           &
!                            + buffer_ref_3D(kgrd(k,i,j,2,1,ig),igrd(i,j,2,ig),jgrd(i,j,2,ig)) &
!                            * hfact(i,j,2,ig) * vfact(k,i,j,2,1,ig)                           &
!                            + buffer_ref_3D(kgrd(k,i,j,3,1,ig),igrd(i,j,3,ig),jgrd(i,j,3,ig)) &
!                            * hfact(i,j,3,ig) * vfact(k,i,j,3,1,ig)                           &
!                            + buffer_ref_3D(kgrd(k,i,j,1,2,ig),igrd(i,j,1,ig),jgrd(i,j,1,ig)) &
!                            * hfact(i,j,1,ig) * vfact(k,i,j,1,2,ig)                           &
!                            + buffer_ref_3D(kgrd(k,i,j,2,2,ig),igrd(i,j,2,ig),jgrd(i,j,2,ig)) &
!                            * hfact(i,j,2,ig) * vfact(k,i,j,2,2,ig)                           &
!                            + buffer_ref_3D(kgrd(k,i,j,3,2,ig),igrd(i,j,3,ig),jgrd(i,j,3,ig)) &
!                            * hfact(i,j,3,ig) * vfact(k,i,j,3,2,ig)
!             end do
!             end do
!             end do
!          else
!             ! logarithmic weighted interpolation
!             do j = 1, DAUGHTER_JA(HANDLE)
!             do i = 1, DAUGHTER_IA(HANDLE)
!             do k = DATR_KS(HANDLE), DATR_KE(HANDLE)
!                dvar(k,i,j) = exp( buffer_ref_3D(kgrd(k,i,j,1,1,ig),igrd(i,j,1,ig),jgrd(i,j,1,ig)) &
!                                 * hfact(i,j,1,ig) * vfact(k,i,j,1,1,ig)                           &
!                                 + buffer_ref_3D(kgrd(k,i,j,2,1,ig),igrd(i,j,2,ig),jgrd(i,j,2,ig)) &
!                                 * hfact(i,j,2,ig) * vfact(k,i,j,2,1,ig)                           &
!                                 + buffer_ref_3D(kgrd(k,i,j,3,1,ig),igrd(i,j,3,ig),jgrd(i,j,3,ig)) &
!                                 * hfact(i,j,3,ig) * vfact(k,i,j,3,1,ig)                           &
!                                 + buffer_ref_3D(kgrd(k,i,j,1,2,ig),igrd(i,j,1,ig),jgrd(i,j,1,ig)) &
!                                 * hfact(i,j,1,ig) * vfact(k,i,j,1,2,ig)                           &
!                                 + buffer_ref_3D(kgrd(k,i,j,2,2,ig),igrd(i,j,2,ig),jgrd(i,j,2,ig)) &
!                                 * hfact(i,j,2,ig) * vfact(k,i,j,2,2,ig)                           &
!                                 + buffer_ref_3D(kgrd(k,i,j,3,2,ig),igrd(i,j,3,ig),jgrd(i,j,3,ig)) &
!                                 * hfact(i,j,3,ig) * vfact(k,i,j,3,2,ig) )
!             end do
!             end do
!             end do
!          endif
          call INTRPNEST_interp_3d( dvar,                 &
                                    buffer_ref_3D,        &
                                    hfact(:,:,:,ig),      &
                                    vfact(:,:,:,:,:,ig),  &
                                    kgrd(:,:,:,:,:,ig),   &
                                    igrd(:,:,:,ig),       &
                                    jgrd(:,:,:,ig),       &
                                    DAUGHTER_IA(HANDLE),  &
                                    DAUGHTER_JA(HANDLE),  &
                                    DATR_KS(HANDLE),      &
                                    DATR_KE(HANDLE),      &
                                    logarithmic           )
       else
          ! linear interpolation (z-staggered)
!          do j = 1, DAUGHTER_JA(HANDLE)
!          do i = 1, DAUGHTER_IA(HANDLE)
!          do k = DATR_KS(HANDLE), DATR_KE(HANDLE)
!             dvar(k,i,j) = buffer_ref_3DF(kgrd(k,i,j,1,1,ig),igrd(i,j,1,ig),jgrd(i,j,1,ig)) &
!                         * hfact(i,j,1,ig) * vfact(k,i,j,1,1,ig)                            &
!                         + buffer_ref_3DF(kgrd(k,i,j,2,1,ig),igrd(i,j,2,ig),jgrd(i,j,2,ig)) &
!                         * hfact(i,j,2,ig) * vfact(k,i,j,2,1,ig)                            &
!                         + buffer_ref_3DF(kgrd(k,i,j,3,1,ig),igrd(i,j,3,ig),jgrd(i,j,3,ig)) &
!                         * hfact(i,j,3,ig) * vfact(k,i,j,3,1,ig)                            &
!                         + buffer_ref_3DF(kgrd(k,i,j,1,2,ig),igrd(i,j,1,ig),jgrd(i,j,1,ig)) &
!                         * hfact(i,j,1,ig) * vfact(k,i,j,1,2,ig)                            &
!                         + buffer_ref_3DF(kgrd(k,i,j,2,2,ig),igrd(i,j,2,ig),jgrd(i,j,2,ig)) &
!                         * hfact(i,j,2,ig) * vfact(k,i,j,2,2,ig)                            &
!                         + buffer_ref_3DF(kgrd(k,i,j,3,2,ig),igrd(i,j,3,ig),jgrd(i,j,3,ig)) &
!                         * hfact(i,j,3,ig) * vfact(k,i,j,3,2,ig)
!          end do
!          end do
!          end do
          call INTRPNEST_interp_3d( dvar,                 &
                                    buffer_ref_3DF,       &
                                    hfact(:,:,:,ig),      &
                                    vfact(:,:,:,:,:,ig),  &
                                    kgrd(:,:,:,:,:,ig),   &
                                    igrd(:,:,:,ig),       &
                                    jgrd(:,:,:,ig),       &
                                    DAUGHTER_IA(HANDLE),  &
                                    DAUGHTER_JA(HANDLE),  &
                                    DATR_KS(HANDLE),      &
                                    DATR_KE(HANDLE),      &
                                    logarithmic           )
       endif

       call PROF_rapend('NESTCOM interp')

    else
    !---------------------------------------------------
       write(*,*) 'xxx internal error [nest/grid]'
       call PRC_MPIstop
    endif

    return
  end subroutine NEST_COMM_intercomm_nestdown_3D

  !-----------------------------------------------------------------------------
  !> [substance of issuer] Inter-communication from parent to daughter: nestdown
  subroutine NEST_COMM_issuer_of_receive_3D( &
      tagbase,   & ! [in ]
      id_stag,   & ! [in ]
      HANDLE,    & ! [in ]
      isu_tag,   & ! [inout]
      isu_tagf,  & ! [inout]
      flag_dens  ) ! [in ]: optional
    use scale_process, only: &
       PRC_myrank,  &
       PRC_nmax,    &
       PRC_MPIstop
    use scale_comm, only: &
       COMM_datatype
    implicit none

    integer , intent(in) :: tagbase              !< communication tag of the variable
    integer , intent(in) :: id_stag              !< id of staggered grid option
    integer , intent(in) :: HANDLE               !< id number of nesting relation in this process target
    integer , intent(inout) :: isu_tag           !< tag for receive buffer
    integer , intent(inout) :: isu_tagf          !< tag for receive buffer
    logical , intent(in), optional :: flag_dens  !< flag of logarithmic interpolation for density

    integer :: ierr, ileng
    integer :: tag, target_rank

    integer :: ig, rq, yp
    logical :: no_zstag    = .true.
    logical :: logarithmic = .false.
    !---------------------------------------------------------------------------

    logarithmic = .false.
    if ( present(flag_dens) ) then
    if ( flag_dens ) then
       logarithmic = .true.
    endif
    endif

    if ( id_stag == I_SCLR )     then
       no_zstag = .true.
       ig       = I_SCLR
    elseif ( id_stag == I_ZSTG ) then
       no_zstag = .false.
       ig       = I_ZSTG
    elseif ( id_stag == I_XSTG ) then
       no_zstag = .true.
       ig       = I_XSTG
    elseif ( id_stag == I_YSTG ) then
       no_zstag = .true.
       ig       = I_YSTG
    endif

    if ( no_zstag ) then
       ileng = PARENT_KA(HANDLE) * PARENT_IA(HANDLE) * PARENT_JA(HANDLE)
    else
       ileng = (PARENT_KA(HANDLE)+1) * PARENT_IA(HANDLE) * PARENT_JA(HANDLE)
    endif

    if ( NEST_Filiation( INTERCOMM_ID(HANDLE) ) > 0 ) then
    !--------------------------------------------------- parent
       ! nothing to do
       ! rq = rq_ctl_p
       ! rq_ctl_p = rq
    elseif ( NEST_Filiation( INTERCOMM_ID(HANDLE) ) < 0 ) then
    !--------------------------------------------------- daughter
       rq = rq_ctl_d

       do yp = 1, NEST_TILE_ALL ! YP Loop
          rq = rq + 1

          call PROF_rapstart('NESTCOM recv child')
          target_rank = NEST_TILE_LIST_d(yp,PRC_myrank+1)

          tag = tagbase + call_order(yp)
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
             call PROF_rapend('NESTCOM recv child')
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
             call PROF_rapend('NESTCOM recv child')
          endif

       enddo ! YP Loop

       if ( isu_tag > max_isu .or. isu_tagf > max_isuf ) then
          write(*,*) 'xxx Exceeded maximum issue [receive: nest/grid]'
          write(*,*) 'xxx isu_tag  = ', isu_tag
          write(*,*) 'xxx isu_tagf = ', isu_tagf
          call PRC_MPIstop
       endif

       rq_ctl_d = rq
    else
    !---------------------------------------------------
       write(*,*) 'xxx internal error [receive: nest/grid]'
       call PRC_MPIstop
    endif

    return
  end subroutine NEST_COMM_issuer_of_receive_3D

  !-----------------------------------------------------------------------------
  !> [substance of issuer] Inter-communication from parent to daughter: nestdown
  subroutine NEST_COMM_issuer_of_wait_3D( &
      HANDLE    ) ! [in]
    use scale_process, only: &
       PRC_MPIstop
    implicit none

    integer, intent(in) :: HANDLE  !< id number of nesting relation in this process target
    !---------------------------------------------------------------------------

    if ( NEST_Filiation( INTERCOMM_ID(HANDLE) ) > 0 ) then
    !--------------------------------------------------- parent
       call NEST_COMM_waitall( rq_tot_p, ireq_p )

    elseif ( NEST_Filiation( INTERCOMM_ID(HANDLE) ) < 0 ) then
    !--------------------------------------------------- daughter
       ! nothing to do
    else
    !---------------------------------------------------
       write(*,*) 'xxx internal error [wait: nest/grid]'
       call PRC_MPIstop
    endif

    return
  end subroutine NEST_COMM_issuer_of_wait_3D

  !-----------------------------------------------------------------------------
  !> [substance of comm_wait] Inter-communication
  subroutine NEST_COMM_waitall( &
      req_count,  & ! [in]
      ireq        ) ! [in]
    use scale_process, only: &
       PRC_MPIstop
    implicit none

    integer, intent(in)    :: req_count
    integer, intent(inout) :: ireq(max_rq)

    integer :: ierr
    integer :: istatus(MPI_STATUS_SIZE,req_count)

    logical    :: flag
    integer(8) :: num
    !---------------------------------------------------------------------------
    num  = 0
    flag = .false.

    do while ( .not. flag )
       num = num + 1
       call MPI_TESTALL( req_count, ireq, flag, istatus, ierr )

       if ( num > ONLINE_WAIT_LIMIT ) then
          if( IO_L ) write(IO_FID_LOG,'(1x,A)') '*** ERROR: over the limit of waiting time [NESTCOM]'
          write(*,'(1x,A)') '*** ERROR: over the limit of waiting time [NESTCOM]'
          call PRC_MPIstop
       endif
    enddo

    return
  end subroutine NEST_COMM_waitall

  !-----------------------------------------------------------------------------
  !> [check communication status] Inter-communication
  subroutine NEST_COMM_test( &
      HANDLE    ) ! [in]
    use scale_process, only: &
       PRC_MPIstop
    implicit none

    integer, intent(in) :: HANDLE  !< id number of nesting relation in this process target

    integer :: istatus(MPI_STATUS_SIZE)
    integer :: ierr
    logical :: flag
    !---------------------------------------------------------------------------

    if ( NEST_Filiation( INTERCOMM_ID(HANDLE) ) > 0 ) then
    !--------------------------------------------------- parent
       call PROF_rapstart('NESTCOM test parent')
       if ( rq_ctl_p > 0 ) call MPI_TEST(ireq_p(1), flag, istatus, ierr)
       call PROF_rapend('NESTCOM test parent')

    elseif ( NEST_Filiation( INTERCOMM_ID(HANDLE) ) < 0 ) then
    !--------------------------------------------------- daughter
       call PROF_rapstart('NESTCOM test child')
       if ( rq_ctl_d > 0 ) call MPI_TEST(ireq_d(1), flag, istatus, ierr)
       call PROF_rapend('NESTCOM test child')
    else
    !---------------------------------------------------
       write(*,*) 'xxx internal error [test: nest/grid]'
       call PRC_MPIstop
    endif

    return
  end subroutine NEST_COMM_test

  !-----------------------------------------------------------------------------
  !> [finalize: disconnect] Inter-communication
  subroutine NEST_COMM_disconnect ( )
    use scale_process, only: &
       GLOBAL_COMM_WORLD
    implicit none

    integer :: ierr
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,'(1x,A)') '*** Waiting finish of whole processes'
    call MPI_BARRIER(GLOBAL_COMM_WORLD, ierr)

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
  end subroutine NEST_COMM_disconnect

end module scale_grid_nest
