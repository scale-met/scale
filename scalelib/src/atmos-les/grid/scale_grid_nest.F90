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
  use scale_grid_index
  use scale_const, only: &
     r_in_m => CONST_RADIUS
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

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  integer,  public              :: NEST_Filiation(10)   !< index of parent-daughter relation
  integer,  public              :: NEST_HANDLING_NUM    !< parent tile number in y-direction
  integer,  public              :: NEST_TILE_NUM_X      !< parent tile number in x-direction
  integer,  public              :: NEST_TILE_NUM_Y      !< parent tile number in y-direction
  integer,  public, allocatable :: NEST_TILE_ID(:)      !< parent tile real id

  integer,  public              :: PARENT_KMAX          !< parent max number in z-direction
  integer,  public              :: PARENT_IMAX          !< parent max number in x-direction
  integer,  public              :: PARENT_JMAX          !< parent max number in y-direction
  integer,  public              :: PARENT_KA            !< parent max number in z-direction (with halo)
  integer,  public              :: PARENT_IA            !< parent max number in x-direction (with halo)
  integer,  public              :: PARENT_JA            !< parent max number in y-direction (with halo)
  integer,  public              :: PARENT_LKMAX         !< parent max number in lz-direction
  real(DP), public              :: PARENT_DTSEC         !< parent DT [sec]

  integer,  public              :: DAUGHTER_KMAX        !< daughter max number in z-direction
  integer,  public              :: DAUGHTER_IMAX        !< daughter max number in x-direction
  integer,  public              :: DAUGHTER_JMAX        !< daughter max number in y-direction
  integer,  public              :: DAUGHTER_KA          !< daughter max number in z-direction (with halo)
  integer,  public              :: DAUGHTER_IA          !< daughter max number in x-direction (with halo)
  integer,  public              :: DAUGHTER_JA          !< daughter max number in y-direction (with halo)
  integer,  public              :: DAUGHTER_LKMAX       !< daughter max number in lz-direction
  real(DP), public              :: DAUGHTER_DTSEC       !< daughter DT [sec]

  integer,  public              :: PRNT_KS              !< start index in z-direction in parent
  integer,  public              :: PRNT_KE              !< end index   in z-direction in parent
  integer,  public              :: PRNT_IS              !< start index in x-direction in parent
  integer,  public              :: PRNT_IE              !< end index   in x-direction in parent
  integer,  public              :: PRNT_JS              !< start index in y-direction in parent
  integer,  public              :: PRNT_JE              !< end index   in y-direction in parent

  integer,  public              :: DATR_KS              !< start index in z-direction in daughter
  integer,  public              :: DATR_KE              !< end index   in z-direction in daughter
  integer,  public              :: DATR_IS              !< start index in x-direction in daughter
  integer,  public              :: DATR_IE              !< end index   in x-direction in daughter
  integer,  public              :: DATR_JS              !< start index in y-direction in daughter
  integer,  public              :: DATR_JE              !< end index   in y-direction in daughter

  integer,  public              :: TILEAL_KA            !< cells of all tiles in z-direction
  integer,  public              :: TILEAL_IA            !< cells of all tiles in x-direction
  integer,  public              :: TILEAL_JA            !< cells of all tiles in y-direction

  logical,  public              :: USE_NESTING        = .false.
  logical,  public              :: OFFLINE            = .true.
  logical,  public              :: ONLINE_PARENT      = .false.
  logical,  public              :: ONLINE_DAUGHTER    = .false.
  integer,  public              :: ONLINE_DOMAIN_NUM  = 1

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

  interface NEST_COMM_intercomm_nestdown
     module procedure NEST_COMM_intercomm_nestdown_3D
  end interface NEST_COMM_intercomm_nestdown

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  real(RP), private, allocatable :: latlon_catalog(:,:,:)  !< parent latlon catalog [rad]
  real(RP), private              :: corner_loc(4,2)        !< local corner location [rad]

  integer, private               :: PARENT_PRC_NUM_X       !< MPI processes in x-direction in parent
  integer, private               :: PARENT_PRC_NUM_Y       !< MPI processes in y-direction in parent
  integer, private               :: PARENT_PRC_nmax        !< MPI total processes in parent

  integer, private               :: DAUGHTER_PRC_NUM_X     !< MPI processes in x-direction in daughter
  integer, private               :: DAUGHTER_PRC_NUM_Y     !< MPI processes in y-direction in daughter
  integer, private               :: DAUGHTER_PRC_nmax      !< MPI total processes in daughter

  integer, private               :: NEST_TILE_ALL          !< NUM of TILEs in the local node
  integer, private               :: NEST_TILE_ALLMAX       !< MAXNUM of TILEs among whole processes
  integer, private, allocatable  :: NEST_TILE_LIST(:,:)    !< relationship list in whole system for daughter
  integer, private, allocatable  :: NEST_TILE_LIST_YP(:)   !< yellow-page of daughter targets for parent
  integer, private               :: NUM_YP                 !< page number of yellow-page

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
  integer, parameter :: itp_nh   = 3                       !< # of interpolation kinds of horizontal direction
  integer, parameter :: itp_nv   = 2                       !< # of interpolation kinds of vertical direction

  integer, parameter :: tag_lon  = 1
  integer, parameter :: tag_lat  = 2
  integer, parameter :: tag_lonx = 3
  integer, parameter :: tag_latx = 4
  integer, parameter :: tag_lony = 5
  integer, parameter :: tag_laty = 6
  integer, parameter :: tag_cz   = 7
  integer, parameter :: tag_fz   = 8
  integer, parameter :: tag_dens = 9
  integer, parameter :: tag_momx = 10
  integer, parameter :: tag_momy = 11
  integer, parameter :: tag_momz = 12
  integer, parameter :: tag_rhot = 13
  integer, parameter :: tag_qx   = 20

  real(RP), private, parameter :: large_number_one   = 9.999E+15_RP
  real(RP), private, parameter :: large_number_two   = 8.888E+15_RP
  real(RP), private, parameter :: large_number_three = 7.777E+15_RP
  integer,  private            :: interp_search_divnum = 10

  integer, private   :: INTERCOMM_ID
  integer, private   :: INTERCOMM_PARENT
  integer, private   :: INTERCOMM_DAUGHTER

  real(RP), private, allocatable :: buffer_2D (:,:)          ! buffer of cummunicator: 2D (with HALO)
  real(RP), private, allocatable :: buffer_3D (:,:,:)        ! buffer of cummunicator: 3D (with HALO)
  real(RP), private, allocatable :: buffer_3DF(:,:,:)        ! buffer of cummunicator: 3D-Kface (with HALO)

  real(RP), private, allocatable :: buffer_ref_LON (:,:)     ! buffer of cummunicator: DENS (with HALO)
  real(RP), private, allocatable :: buffer_ref_LONX(:,:)     ! buffer of cummunicator: VELZ (with HALO)
  real(RP), private, allocatable :: buffer_ref_LONY(:,:)     ! buffer of cummunicator: VELX (with HALO)
  real(RP), private, allocatable :: buffer_ref_LAT (:,:)     ! buffer of cummunicator: VELY (with HALO)
  real(RP), private, allocatable :: buffer_ref_LATX(:,:)     ! buffer of cummunicator: POTT (with HALO)
  real(RP), private, allocatable :: buffer_ref_LATY(:,:)     ! buffer of cummunicator: POTT (with HALO)
  real(RP), private, allocatable :: buffer_ref_CZ  (:,:,:)   ! buffer of cummunicator: VELY (with HALO)
  real(RP), private, allocatable :: buffer_ref_FZ  (:,:,:)   ! buffer of cummunicator: VELY (with HALO)

  real(RP), private, allocatable :: buffer_ref_2D (:,:)      ! buffer of cummunicator: 2D data (with HALO)
  real(RP), private, allocatable :: buffer_ref_3D (:,:,:)    ! buffer of cummunicator: 3D data (with HALO)
  real(RP), private, allocatable :: buffer_ref_3DF(:,:,:)    ! buffer of cummunicator: 3D at z-Face (with HALO)

  real(RP), private, allocatable :: hfact(:,:,:,:)           ! interpolation factor for horizontal direction
  real(RP), private, allocatable :: vfact(:,:,:,:,:,:)       ! interpolation factor for vertical direction
  integer,  private, allocatable :: kgrd (:,:,:,:,:,:)       ! interpolation target grids in z-axis
  integer,  private, allocatable :: igrd (:,:,:,:)           ! interpolation target grids in x-axis
  integer,  private, allocatable :: jgrd (:,:,:,:)           ! interpolation target grids in y-axis

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine NEST_setup
    use scale_const, only: &
       D2R => CONST_D2R
    use scale_process, only: &
       PRC_nmax,   &
       PRC_master, &
       PRC_myrank, &
       PRC_MPIstop
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

    !< metadata files for lat-lon domain for all processes
    character(len=H_LONG) :: LATLON_CATALOGUE_FNAME = 'latlon_domain_catalogue.txt'
    character(len=H_MID)  :: cmd                    = './scale-les'
    character(20)         :: argv(2)

    integer :: i
    integer :: fid, ierr
    integer :: parent_id
    integer :: istatus(MPI_STATUS_SIZE)
    integer :: errcodes(1:PRC_nmax)

    character(2) :: dom_num

    namelist / PARAM_NEST /    &
       USE_NESTING,            &
       PARENT_PRC_NUM_X,       &
       PARENT_PRC_NUM_Y,       &
       PARENT_KMAX,            &
       PARENT_IMAX,            &
       PARENT_JMAX,            &
       PARENT_LKMAX,           &
       LATLON_CATALOGUE_FNAME, &
       OFFLINE,                &
       ONLINE_DOMAIN_NUM,      &
       ONLINE_PARENT,          &
       ONLINE_DAUGHTER

    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[NEST]/Categ[GRID]'

    argv(1) = 'run.conf'
    argv(2) = ''

    NEST_HANDLING_NUM = 0
    NEST_Filiation(:) = 0

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

    if( USE_NESTING ) then

      corner_loc(I_NW,I_LON) = REAL_LONXY(IS-1,JE  ) / D2R
      corner_loc(I_NE,I_LON) = REAL_LONXY(IE  ,JE  ) / D2R
      corner_loc(I_SW,I_LON) = REAL_LONXY(IS-1,JS-1) / D2R
      corner_loc(I_SE,I_LON) = REAL_LONXY(IE  ,JS-1) / D2R
      corner_loc(I_NW,I_LAT) = REAL_LATXY(IS-1,JE  ) / D2R
      corner_loc(I_NE,I_LAT) = REAL_LATXY(IE  ,JE  ) / D2R
      corner_loc(I_SW,I_LAT) = REAL_LATXY(IS-1,JS-1) / D2R
      corner_loc(I_SE,I_LAT) = REAL_LATXY(IE  ,JS-1) / D2R

      if( OFFLINE ) then

         PARENT_PRC_nmax = PARENT_PRC_NUM_X * PARENT_PRC_NUM_Y
         allocate( latlon_catalog(PARENT_PRC_nmax,4,2) )

         !--- read latlon catalogue
         fid = IO_get_available_fid()
         open( fid,                                    &
               file   = trim(LATLON_CATALOGUE_FNAME),  &
               form   = 'formatted',                   &
               status = 'old',                         &
               iostat = ierr                           )

         if ( ierr /= 0 ) then
            if( IO_L ) write(*,*) 'xxx cannot open latlon-catalogue file!'
            call PRC_MPIstop
         endif

         do i = 1, PARENT_PRC_nmax
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

         call NEST_domain_relate

      else ! ONLINE RELATIONSHIP
         if( ONLINE_PARENT ) then
         !-------------------------------------------------
            INTERCOMM_ID = ONLINE_DOMAIN_NUM
            NEST_Filiation(INTERCOMM_ID) = 1

            ! Launch Daughter Domain
            write(dom_num,'(I2.2)') ONLINE_DOMAIN_NUM+1
            argv(1) = 'run.d'//dom_num//'.conf'
            if( IO_L ) write(IO_FID_LOG,'(1x,A,I2,A)') '*** Launch Daughter Domain [INTERCOMM_ID:', INTERCOMM_ID, ' ]'
            if( IO_L ) write(IO_FID_LOG,'(1x,A,A,A,A)') '*** Launch Command: ', trim(cmd), ' ', trim(argv(1))
            call MPI_COMM_SPAWN( trim(cmd),          &
                                 argv,               &
                                 PRC_nmax,           &
                                 MPI_INFO_NULL,      &
                                 PRC_master,         &
                                 COMM_world,         &
                                 INTERCOMM_DAUGHTER, &
                                 errcodes,           &
                                 ierr                )

            call NEST_COMM_ping( INTERCOMM_ID )

            call NEST_COMM_parentsize( INTERCOMM_ID )

            allocate( latlon_catalog(PARENT_PRC_nmax,4,2) )
            call NEST_COMM_catalogue( INTERCOMM_ID )
            call MPI_BARRIER(INTERCOMM_DAUGHTER, ierr)

            PARENT_KA   = PARENT_KMAX   + KHALO * 2
            PARENT_IA   = PARENT_IMAX   + IHALO * 2
            PARENT_JA   = PARENT_JMAX   + JHALO * 2
            DAUGHTER_KA = DAUGHTER_KMAX + KHALO * 2
            DAUGHTER_IA = DAUGHTER_IMAX + IHALO * 2
            DAUGHTER_JA = DAUGHTER_JMAX + JHALO * 2
            TILEAL_KA   = 0
            TILEAL_IA   = 0
            TILEAL_JA   = 0

            if( IO_L ) write(IO_FID_LOG,'(1x,A)'     ) '***  Informations of Parent Domain'
            if( IO_L ) write(IO_FID_LOG,'(1x,A,I6)'  ) '***  --- PARENT_PRC_nmax   :', PARENT_PRC_nmax
            if( IO_L ) write(IO_FID_LOG,'(1x,A,I6)'  ) '***  --- PARENT_PRC_NUM_X  :', PARENT_PRC_NUM_X
            if( IO_L ) write(IO_FID_LOG,'(1x,A,I6)'  ) '***  --- PARENT_PRC_NUM_Y  :', PARENT_PRC_NUM_Y
            if( IO_L ) write(IO_FID_LOG,'(1x,A,I6)'  ) '***  --- PARENT_KMAX       :', PARENT_KMAX
            if( IO_L ) write(IO_FID_LOG,'(1x,A,I6)'  ) '***  --- PARENT_IMAX       :', PARENT_IMAX
            if( IO_L ) write(IO_FID_LOG,'(1x,A,I6)'  ) '***  --- PARENT_JMAX       :', PARENT_JMAX
            if( IO_L ) write(IO_FID_LOG,'(1x,A,F9.3)') '***  --- PARENT_DTSEC      :', PARENT_DTSEC
            if( IO_L ) write(IO_FID_LOG,'(1x,A)'     ) '***  Informations of Daughter Domain'
            if( IO_L ) write(IO_FID_LOG,'(1x,A,I6)'  ) '***  --- DAUGHTER_PRC_nmax :', DAUGHTER_PRC_nmax
            if( IO_L ) write(IO_FID_LOG,'(1x,A,I6)'  ) '***  --- DAUGHTER_PRC_NUM_X:', DAUGHTER_PRC_NUM_X
            if( IO_L ) write(IO_FID_LOG,'(1x,A,I6)'  ) '***  --- DAUGHTER_PRC_NUM_Y:', DAUGHTER_PRC_NUM_Y
            if( IO_L ) write(IO_FID_LOG,'(1x,A,I6)'  ) '***  --- DAUGHTER_KMAX     :', DAUGHTER_KMAX
            if( IO_L ) write(IO_FID_LOG,'(1x,A,I6)'  ) '***  --- DAUGHTER_IMAX     :', DAUGHTER_IMAX
            if( IO_L ) write(IO_FID_LOG,'(1x,A,I6)'  ) '***  --- DAUGHTER_JMAX     :', DAUGHTER_JMAX
            if( IO_L ) write(IO_FID_LOG,'(1x,A,F9.3)') '***  --- DAUGHTER_DTSEC    :', DAUGHTER_DTSEC

            call NEST_COMM_setup_nestdown

         elseif( ONLINE_DAUGHTER ) then
         !-------------------------------------------------
            INTERCOMM_ID = ONLINE_DOMAIN_NUM - 1
            NEST_Filiation(INTERCOMM_ID) = -1
            NEST_HANDLING_NUM = NEST_HANDLING_NUM + 1

            call MPI_COMM_GET_PARENT( INTERCOMM_PARENT, ierr )
            if( IO_L ) write(IO_FID_LOG,'(1x,A,I2,A)') '*** Activated Daughter Domain [INTERCOMM_ID:', INTERCOMM_ID, ' ]'

            call NEST_COMM_ping( INTERCOMM_ID )

            call NEST_COMM_parentsize( INTERCOMM_ID )

            allocate( latlon_catalog(PARENT_PRC_nmax,4,2) )
            call NEST_COMM_catalogue( INTERCOMM_ID )
            call MPI_BARRIER(INTERCOMM_PARENT, ierr)

            call NEST_domain_relate

            PARENT_KA   = PARENT_KMAX   + KHALO * 2
            PARENT_IA   = PARENT_IMAX   + IHALO * 2
            PARENT_JA   = PARENT_JMAX   + JHALO * 2
            DAUGHTER_KA = DAUGHTER_KMAX + KHALO * 2
            DAUGHTER_IA = DAUGHTER_IMAX + IHALO * 2
            DAUGHTER_JA = DAUGHTER_JMAX + JHALO * 2
            TILEAL_KA   = PARENT_KA
            TILEAL_IA   = PARENT_IA * NEST_TILE_NUM_X
            TILEAL_JA   = PARENT_JA * NEST_TILE_NUM_Y

            if( IO_L ) write(IO_FID_LOG,'(1x,A)'     ) '***  Informations of Parent Domain'
            if( IO_L ) write(IO_FID_LOG,'(1x,A,I6)'  ) '***  --- PARENT_PRC_nmax   :', PARENT_PRC_nmax
            if( IO_L ) write(IO_FID_LOG,'(1x,A,I6)'  ) '***  --- PARENT_PRC_NUM_X  :', PARENT_PRC_NUM_X
            if( IO_L ) write(IO_FID_LOG,'(1x,A,I6)'  ) '***  --- PARENT_PRC_NUM_Y  :', PARENT_PRC_NUM_Y
            if( IO_L ) write(IO_FID_LOG,'(1x,A,I6)'  ) '***  --- PARENT_KMAX       :', PARENT_KMAX
            if( IO_L ) write(IO_FID_LOG,'(1x,A,I6)'  ) '***  --- PARENT_IMAX       :', PARENT_IMAX
            if( IO_L ) write(IO_FID_LOG,'(1x,A,I6)'  ) '***  --- PARENT_JMAX       :', PARENT_JMAX
            if( IO_L ) write(IO_FID_LOG,'(1x,A,F9.3)') '***  --- PARENT_DTSEC      :', PARENT_DTSEC
            if( IO_L ) write(IO_FID_LOG,'(1x,A)'     ) '***  Informations of Daughter Domain'
            if( IO_L ) write(IO_FID_LOG,'(1x,A,I6)'  ) '***  --- DAUGHTER_PRC_nmax :', DAUGHTER_PRC_nmax
            if( IO_L ) write(IO_FID_LOG,'(1x,A,I6)'  ) '***  --- DAUGHTER_PRC_NUM_X:', DAUGHTER_PRC_NUM_X
            if( IO_L ) write(IO_FID_LOG,'(1x,A,I6)'  ) '***  --- DAUGHTER_PRC_NUM_Y:', DAUGHTER_PRC_NUM_Y
            if( IO_L ) write(IO_FID_LOG,'(1x,A,I6)'  ) '***  --- DAUGHTER_KMAX     :', DAUGHTER_KMAX
            if( IO_L ) write(IO_FID_LOG,'(1x,A,I6)'  ) '***  --- DAUGHTER_IMAX     :', DAUGHTER_IMAX
            if( IO_L ) write(IO_FID_LOG,'(1x,A,I6)'  ) '***  --- DAUGHTER_JMAX     :', DAUGHTER_JMAX
            if( IO_L ) write(IO_FID_LOG,'(1x,A,F9.3)') '***  --- DAUGHTER_DTSEC    :', DAUGHTER_DTSEC
            if( IO_L ) write(IO_FID_LOG,'(1x,A)'     ) '***  Informations of Target Tiles'
            if( IO_L ) write(IO_FID_LOG,'(1x,A,I6)'  ) '***  --- TILEALL_KA      :', TILEAL_KA
            if( IO_L ) write(IO_FID_LOG,'(1x,A,I6)'  ) '***  --- TILEALL_IA      :', TILEAL_IA
            if( IO_L ) write(IO_FID_LOG,'(1x,A,I6)'  ) '***  --- TILEALL_JA      :', TILEAL_JA

            allocate( buffer_2D  (              PARENT_IA, PARENT_JA ) )
            allocate( buffer_3D  (   PARENT_KA, PARENT_IA, PARENT_JA ) )
            allocate( buffer_3DF ( 0:PARENT_KA, PARENT_IA, PARENT_JA ) )

            allocate( buffer_ref_LON (            TILEAL_IA,TILEAL_JA) )
            allocate( buffer_ref_LONX(            TILEAL_IA,TILEAL_JA) )
            allocate( buffer_ref_LONY(            TILEAL_IA,TILEAL_JA) )
            allocate( buffer_ref_LAT (            TILEAL_IA,TILEAL_JA) )
            allocate( buffer_ref_LATX(            TILEAL_IA,TILEAL_JA) )
            allocate( buffer_ref_LATY(            TILEAL_IA,TILEAL_JA) )
            allocate( buffer_ref_CZ  (  PARENT_KA,TILEAL_IA,TILEAL_JA) )
            allocate( buffer_ref_FZ  (0:PARENT_KA,TILEAL_IA,TILEAL_JA) )

            allocate( buffer_ref_2D  (            TILEAL_IA,TILEAL_JA) )
            allocate( buffer_ref_3D  (  PARENT_KA,TILEAL_IA,TILEAL_JA) )
            allocate( buffer_ref_3DF (0:PARENT_KA,TILEAL_IA,TILEAL_JA) )

            allocate( hfact(            DAUGHTER_IA,DAUGHTER_JA,itp_nh,       itp_ng) )
            allocate( vfact(DAUGHTER_KA,DAUGHTER_IA,DAUGHTER_JA,itp_nh,itp_nv,itp_ng) )
            allocate( igrd (            DAUGHTER_IA,DAUGHTER_JA,itp_nh,       itp_ng) )
            allocate( jgrd (            DAUGHTER_IA,DAUGHTER_JA,itp_nh,       itp_ng) )
            allocate( kgrd (DAUGHTER_KA,DAUGHTER_IA,DAUGHTER_JA,itp_nh,itp_nv,itp_ng) )

            call NEST_COMM_setup_nestdown

            ! for scalar points
            call latlonz_interporation_fact( hfact          (:,:,:,I_SCLR),     & ! [OUT]
                                             vfact          (:,:,:,:,:,I_SCLR), & ! [OUT]
                                             kgrd           (:,:,:,:,:,I_SCLR), & ! [OUT]
                                             igrd           (:,:,:,I_SCLR),     & ! [OUT]
                                             jgrd           (:,:,:,I_SCLR),     & ! [OUT]
                                             MY_CZ          (:,:,:),            & ! [IN]
                                             MY_LAT         (:,:),              & ! [IN]
                                             MY_LON         (:,:),              & ! [IN]
                                             buffer_ref_CZ  (:,:,:),            & ! [IN]
                                             buffer_ref_LAT (:,:),              & ! [IN]
                                             buffer_ref_LON (:,:),              & ! [IN]
                                             TILEAL_KA,                         & ! [IN]
                                             TILEAL_IA,                         & ! [IN]
                                             TILEAL_JA                          ) ! [IN]

            ! for z staggered points
            call latlonz_interporation_fact( hfact          (:,:,:,I_ZSTG),     & ! [OUT]
                                             vfact          (:,:,:,:,:,I_ZSTG), & ! [OUT]
                                             kgrd           (:,:,:,:,:,I_ZSTG), & ! [OUT]
                                             igrd           (:,:,:,I_ZSTG),     & ! [OUT]
                                             jgrd           (:,:,:,I_ZSTG),     & ! [OUT]
                                             MY_FZ          (:,:,:),            & ! [IN]
                                             MY_LAT         (:,:),              & ! [IN]
                                             MY_LON         (:,:),              & ! [IN]
                                             buffer_ref_FZ  (:,:,:),            & ! [IN]
                                             buffer_ref_LAT (:,:),              & ! [IN]
                                             buffer_ref_LON (:,:),              & ! [IN]
                                             TILEAL_KA,                         & ! [IN]
                                             TILEAL_IA,                         & ! [IN]
                                             TILEAL_JA                          ) ! [IN]

            ! for x staggered points
            call latlonz_interporation_fact( hfact          (:,:,:,I_XSTG),     & ! [OUT]
                                             vfact          (:,:,:,:,:,I_XSTG), & ! [OUT]
                                             kgrd           (:,:,:,:,:,I_XSTG), & ! [OUT]
                                             igrd           (:,:,:,I_XSTG),     & ! [OUT]
                                             jgrd           (:,:,:,I_XSTG),     & ! [OUT]
                                             MY_CZ          (:,:,:),            & ! [IN]
                                             MY_LATX        (:,:),              & ! [IN]
                                             MY_LONX        (:,:),              & ! [IN]
                                             buffer_ref_CZ  (:,:,:),            & ! [IN]
                                             buffer_ref_LATX(:,:),              & ! [IN]
                                             buffer_ref_LONX(:,:),              & ! [IN]
                                             TILEAL_KA,                         & ! [IN]
                                             TILEAL_IA,                         & ! [IN]
                                             TILEAL_JA                          ) ! [IN]

            ! for y staggered points
            call latlonz_interporation_fact( hfact          (:,:,:,I_YSTG),     & ! [OUT]
                                             vfact          (:,:,:,:,:,I_YSTG), & ! [OUT]
                                             kgrd           (:,:,:,:,:,I_YSTG), & ! [OUT]
                                             igrd           (:,:,:,I_YSTG),     & ! [OUT]
                                             jgrd           (:,:,:,I_YSTG),     & ! [OUT]
                                             MY_CZ          (:,:,:),            & ! [IN]
                                             MY_LATY        (:,:),              & ! [IN]
                                             MY_LONY        (:,:),              & ! [IN]
                                             buffer_ref_CZ  (:,:,:),            & ! [IN]
                                             buffer_ref_LATY(:,:),              & ! [IN]
                                             buffer_ref_LONY(:,:),              & ! [IN]
                                             TILEAL_KA,                         & ! [IN]
                                             TILEAL_IA,                         & ! [IN]
                                             TILEAL_JA                          ) ! [IN]

         !-------------------------------------------------
         endif

      endif

    endif !USE_NESTING

    return
  end subroutine NEST_setup

  !-----------------------------------------------------------------------------
  !> Solve relationship between ParentDomain & Daughter Domain
  subroutine NEST_domain_relate
    use scale_const, only: &
       EPS => CONST_EPS
    use scale_process, only: &
       PRC_MPIstop
    implicit none

    logical :: hit

    integer, allocatable :: pd_tile_num(:,:)

    integer :: pd_sw_tile
    integer :: pd_ne_tile
    integer :: i, j, ii, jj, k
    !---------------------------------------------------------------------------

    allocate( pd_tile_num(0:PARENT_PRC_nmax-1,2) )

    k = 0 ! MPI process number starts from zero
    do j = 1, PARENT_PRC_NUM_Y
    do i = 1, PARENT_PRC_NUM_X
       pd_tile_num(k,1) = i
       pd_tile_num(k,2) = j
       k = k + 1
    enddo
    enddo

    !--- SW search
    hit = .false.
    do i = 1, PARENT_PRC_nmax
       if ( corner_loc(I_SW,I_LON) > latlon_catalog(i,I_SW,I_LON) .and. &
            corner_loc(I_SW,I_LAT) > latlon_catalog(i,I_SW,I_LAT) .and. &
            corner_loc(I_SW,I_LON) < latlon_catalog(i,I_NE,I_LON) .and. &
            corner_loc(I_SW,I_LAT) < latlon_catalog(i,I_NE,I_LAT) .or.  &
            abs( corner_loc(I_SW,I_LON) - latlon_catalog(i,I_SW,I_LON) ) < EPS .and. &
            abs( corner_loc(I_SW,I_LAT) - latlon_catalog(i,I_SW,I_LAT) ) < EPS       ) then

          pd_sw_tile = i-1 ! MPI process number starts from zero
          hit = .true.
          exit ! exit loop
       endif
    enddo
    if ( .NOT. hit ) then
       if( IO_L ) write(*,*) 'xxx domain mismatch between parent and daughter: SW search'
       call PRC_MPIstop
    endif

    !--- NE search
    hit = .false.
    do i = 1, PARENT_PRC_nmax
       if ( corner_loc(I_NE,I_LON) > latlon_catalog(i,I_SW,I_LON) .and. &
            corner_loc(I_NE,I_LAT) > latlon_catalog(i,I_SW,I_LAT) .and. &
            corner_loc(I_NE,I_LON) < latlon_catalog(i,I_NE,I_LON) .and. &
            corner_loc(I_NE,I_LAT) < latlon_catalog(i,I_NE,I_LAT) .or.  &
            abs( corner_loc(I_NE,I_LON) - latlon_catalog(i,I_NE,I_LON) ) < EPS .and. &
            abs( corner_loc(I_NE,I_LAT) - latlon_catalog(i,I_NE,I_LAT) ) < EPS       ) then

          pd_ne_tile = i-1 ! MPI process number starts from zero
          hit = .true.
          exit ! exit loop
       endif
    enddo
    if ( .NOT. hit ) then
       if( IO_L ) write(*,*) 'xxx domain mismatch between parent and daughter: NE search'
       call PRC_MPIstop
    endif

    NEST_TILE_NUM_X = pd_tile_num(pd_ne_tile,1) - pd_tile_num(pd_sw_tile,1) + 1
    NEST_TILE_NUM_Y = pd_tile_num(pd_ne_tile,2) - pd_tile_num(pd_sw_tile,2) + 1

    allocate( NEST_TILE_ID( NEST_TILE_NUM_X*NEST_TILE_NUM_Y ) )

    k = 1
    do j = 1, NEST_TILE_NUM_Y
    do i = 1, NEST_TILE_NUM_X
       NEST_TILE_ID(k) = pd_sw_tile + (i-1) + PARENT_PRC_NUM_X*(j-1)
       k = k + 1
    enddo
    enddo

    return
  end subroutine NEST_domain_relate

  !-----------------------------------------------------------------------------
  !> Get parent domain size
  subroutine NEST_COMM_parentsize( &
      INTERCOMM_ID  )
    use scale_process, only: &
       PRC_myrank,  &
       PRC_nmax,    &
       PRC_NUM_X,   &
       PRC_NUM_Y,   &
       PRC_MPIstop
    use scale_time, only: &
       TIME_DTSEC
    implicit none

    integer, intent(in) :: INTERCOMM_ID !< id of domain-intercommunication

    real(DP) :: buffer
    integer  :: datapack(12)
    integer  :: parent, daughter
    integer  :: ireq1, ireq2, ierr1, ierr2, ileng
    integer  :: istatus(MPI_STATUS_SIZE)
    integer  :: tag
    !---------------------------------------------------------------------------

    tag        = INTERCOMM_ID * 100
    parent     = INTERCOMM_ID
    daughter   = INTERCOMM_ID + 1
    ileng = 12

    if ( ONLINE_DOMAIN_NUM == parent ) then
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
       buffer      = TIME_DTSEC

       call MPI_ISEND(datapack, ileng, MPI_INTEGER, PRC_myrank, tag, INTERCOMM_DAUGHTER, ireq1, ierr1)
       call MPI_ISEND(buffer, 1, MPI_DOUBLE_PRECISION, PRC_myrank, tag+1, INTERCOMM_DAUGHTER, ireq2, ierr2)
       call MPI_WAIT(ireq1, istatus, ierr1)
       call MPI_WAIT(ireq2, istatus, ierr2)

       PARENT_PRC_nmax  = datapack( 1)
       PARENT_PRC_NUM_X = datapack( 2)
       PARENT_PRC_NUM_Y = datapack( 3)
       PARENT_KMAX      = datapack( 4)
       PARENT_IMAX      = datapack( 5)
       PARENT_JMAX      = datapack( 6)
       PRNT_KS          = datapack( 7)
       PRNT_KE          = datapack( 8)
       PRNT_IS          = datapack( 9)
       PRNT_IE          = datapack(10)
       PRNT_JS          = datapack(11)
       PRNT_JE          = datapack(12)
       PARENT_DTSEC     = buffer

       ! from daughter to parent
       call MPI_IRECV(datapack, ileng, MPI_INTEGER, PRC_myrank, tag+2, INTERCOMM_DAUGHTER, ireq1, ierr1)
       call MPI_IRECV(buffer, 1, MPI_DOUBLE_PRECISION, PRC_myrank, tag+3, INTERCOMM_DAUGHTER, ireq2, ierr2)
       call MPI_WAIT(ireq1, istatus, ierr1)
       call MPI_WAIT(ireq2, istatus, ierr2)

       DAUGHTER_PRC_nmax  = datapack( 1)
       DAUGHTER_PRC_NUM_X = datapack( 2)
       DAUGHTER_PRC_NUM_Y = datapack( 3)
       DAUGHTER_KMAX      = datapack( 4)
       DAUGHTER_IMAX      = datapack( 5)
       DAUGHTER_JMAX      = datapack( 6)
       DATR_KS            = datapack( 7)
       DATR_KE            = datapack( 8)
       DATR_IS            = datapack( 9)
       DATR_IE            = datapack(10)
       DATR_JS            = datapack(11)
       DATR_JE            = datapack(12)
       DAUGHTER_DTSEC     = buffer

    elseif ( ONLINE_DOMAIN_NUM == daughter ) then
       ! from parent to daughter
       call MPI_IRECV(datapack, ileng, MPI_INTEGER, PRC_myrank, tag, INTERCOMM_PARENT, ireq1, ierr1)
       call MPI_IRECV(buffer, 1, MPI_DOUBLE_PRECISION, PRC_myrank, tag+1, INTERCOMM_PARENT, ireq2, ierr2)
       call MPI_WAIT(ireq1, istatus, ierr1)
       call MPI_WAIT(ireq2, istatus, ierr2)

       PARENT_PRC_nmax  = datapack( 1)
       PARENT_PRC_NUM_X = datapack( 2)
       PARENT_PRC_NUM_Y = datapack( 3)
       PARENT_KMAX      = datapack( 4)
       PARENT_IMAX      = datapack( 5)
       PARENT_JMAX      = datapack( 6)
       PRNT_KS          = datapack( 7)
       PRNT_KE          = datapack( 8)
       PRNT_IS          = datapack( 9)
       PRNT_IE          = datapack(10)
       PRNT_JS          = datapack(11)
       PRNT_JE          = datapack(12)
       PARENT_DTSEC     = buffer

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
       buffer       = TIME_DTSEC

       call MPI_ISEND(datapack, ileng, MPI_INTEGER, PRC_myrank, tag+2, INTERCOMM_PARENT, ireq1, ierr1)
       call MPI_ISEND(buffer, 1, MPI_DOUBLE_PRECISION, PRC_myrank, tag+3, INTERCOMM_PARENT, ireq2, ierr2)
       call MPI_WAIT(ireq1, istatus, ierr1)
       call MPI_WAIT(ireq2, istatus, ierr2)

       DAUGHTER_PRC_nmax  = datapack( 1)
       DAUGHTER_PRC_NUM_X = datapack( 2)
       DAUGHTER_PRC_NUM_Y = datapack( 3)
       DAUGHTER_KMAX      = datapack( 4)
       DAUGHTER_IMAX      = datapack( 5)
       DAUGHTER_JMAX      = datapack( 6)
       DATR_KS            = datapack( 7)
       DATR_KE            = datapack( 8)
       DATR_IS            = datapack( 9)
       DATR_IE            = datapack(10)
       DATR_JS            = datapack(11)
       DATR_JE            = datapack(12)
       DAUGHTER_DTSEC     = buffer
    else
       if( IO_L ) write(*,*) 'xxx internal error [nest/grid]'
       call PRC_MPIstop
    endif

    return
  end subroutine NEST_COMM_parentsize

  !-----------------------------------------------------------------------------
  !> Get parent latlon catalogue
  subroutine NEST_COMM_catalogue( &
      INTERCOMM_ID  )
    use scale_process, only: &
       PRC_myrank,  &
       PRC_nmax,    &
       PRC_MPIstop
    use scale_grid_real, only: &
       REAL_DOMAIN_CATALOGUE
    use scale_comm, only: &
       COMM_datatype
    implicit none

    integer, intent(in) :: INTERCOMM_ID !< id of domain-intercommunication

    integer :: parent, daughter
    integer :: ireq, ierr, ileng
    integer :: istatus(MPI_STATUS_SIZE)
    integer :: tag
    !---------------------------------------------------------------------------

    tag        = INTERCOMM_ID * 100
    parent     = INTERCOMM_ID
    daughter   = INTERCOMM_ID + 1

    if ( ONLINE_DOMAIN_NUM == parent ) then
       ileng = PRC_nmax * 4 * 2
       call MPI_ISEND(REAL_DOMAIN_CATALOGUE, ileng, COMM_datatype, PRC_myrank, tag, INTERCOMM_DAUGHTER, ireq, ierr)
       call MPI_WAIT(ireq, istatus, ierr)

    elseif ( ONLINE_DOMAIN_NUM == daughter ) then
       ileng = PARENT_PRC_nmax * 4 * 2
       call MPI_IRECV(latlon_catalog, ileng, COMM_datatype, PRC_myrank, tag, INTERCOMM_PARENT, ireq, ierr)
       call MPI_WAIT(ireq, istatus, ierr)

    else
       if( IO_L ) write(*,*) 'xxx internal error [nest/grid]'
       call PRC_MPIstop
    endif

    return
  end subroutine NEST_COMM_catalogue

  !-----------------------------------------------------------------------------
  !> Check Communication Inter-domains
  subroutine NEST_COMM_ping( &
      INTERCOMM_ID  )
    use scale_process, only: &
       PRC_myrank,  &
       PRC_MPIstop
    implicit none

    integer, intent(in) :: INTERCOMM_ID !< id of domain-intercommunication

    integer :: ping, pong
    integer :: parent, daughter
    integer :: ireq1, ireq2, ierr1, ierr2
    integer :: istatus(MPI_STATUS_SIZE)
    integer :: tag
    logical :: ping_error
    !---------------------------------------------------------------------------

    tag        = INTERCOMM_ID * 100
    parent     = INTERCOMM_ID
    daughter   = INTERCOMM_ID + 1
    ping_error = .false.

    if ( ONLINE_DOMAIN_NUM == parent ) then
       ping = ONLINE_DOMAIN_NUM
       pong = 0

       call MPI_ISEND(ping, 1, MPI_INTEGER, PRC_myrank, tag+1, INTERCOMM_DAUGHTER, ireq1, ierr1)
       call MPI_IRECV(pong, 1, MPI_INTEGER, PRC_myrank, tag+2, INTERCOMM_DAUGHTER, ireq2, ierr2)
       call MPI_WAIT(ireq1, istatus, ierr1)
       call MPI_WAIT(ireq2, istatus, ierr2)

       if ( pong /= daughter ) ping_error = .true.

    elseif ( ONLINE_DOMAIN_NUM == daughter ) then
       ping = ONLINE_DOMAIN_NUM
       pong = 0

       call MPI_ISEND(ping, 1, MPI_INTEGER, PRC_myrank, tag+2, INTERCOMM_PARENT, ireq1, ierr1)
       call MPI_IRECV(pong, 1, MPI_INTEGER, PRC_myrank, tag+1, INTERCOMM_PARENT, ireq2, ierr2)
       call MPI_WAIT(ireq1, istatus, ierr1)
       call MPI_WAIT(ireq2, istatus, ierr2)

       if ( pong /= parent ) ping_error = .true.

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
  subroutine NEST_COMM_setup_nestdown
    use scale_process, only: &
       PRC_myrank,  &
       PRC_nmax,    &
       PRC_MPIstop
    use scale_grid_real, only: &
       REAL_DOMAIN_CATALOGUE
    use scale_comm, only: &
       COMM_datatype,  &
       COMM_world
    implicit none

    integer, allocatable :: buffer_LIST(:)
    integer, allocatable :: buffer_ALLLIST(:)

    integer :: parent, daughter
    integer :: ireq, ierr, ileng
    integer :: istatus(MPI_STATUS_SIZE)
    integer :: tag

    integer :: i, j, k
    !---------------------------------------------------------------------------

    tag        = INTERCOMM_ID * 100
    parent     = INTERCOMM_ID
    daughter   = INTERCOMM_ID + 1

    if ( ONLINE_DOMAIN_NUM == parent ) then
    !---------------------------------------------------

       call MPI_IRECV(NEST_TILE_ALLMAX, 1, MPI_INTEGER, PRC_myrank, tag+1, INTERCOMM_DAUGHTER, ireq, ierr)
       call MPI_WAIT(ireq, istatus, ierr)

       allocate( NEST_TILE_LIST   (NEST_TILE_ALLMAX,PRC_nmax) )
       allocate( NEST_TILE_LIST_YP(NEST_TILE_ALLMAX*PRC_nmax) )

       ileng = NEST_TILE_ALLMAX*PRC_nmax
       call MPI_IRECV(NEST_TILE_LIST, ileng, MPI_INTEGER, PRC_myrank, tag+2, INTERCOMM_DAUGHTER, ireq, ierr)
       call MPI_WAIT(ireq, istatus, ierr)

       NEST_TILE_LIST_YP(:) = -1

       k = 0
       do j = 1, PRC_nmax
       do i = 1, NEST_TILE_ALLMAX
          if ( NEST_TILE_LIST(i,j) == PRC_myrank ) then
             k = k + 1
             NEST_TILE_LIST_YP(k) = j - 1  !rank number is started from 1
          endif
       enddo
       enddo
       NUM_YP = k

       if( IO_L ) write(IO_FID_LOG,'(A,I5,A,I5)') "   Num YP =",NUM_YP,"  Num TILE(MAX) =",NEST_TILE_ALLMAX

       call NEST_COMM_importgrid_nestdown

       call MPI_BARRIER(INTERCOMM_DAUGHTER, ierr)
    elseif ( ONLINE_DOMAIN_NUM == daughter ) then
    !---------------------------------------------------

       NEST_TILE_ALL = size( NEST_TILE_ID(:) ) ! should be equal to "NEST_TILE_NUM_X*NEST_TILE_NUM_Y"
       call MPI_Allreduce( NEST_TILE_ALL,    &
                           NEST_TILE_ALLMAX, &
                           1,                &
                           MPI_INTEGER,      &
                           MPI_MAX,          &
                           COMM_world,       &
                           ierr              )
       if( IO_L ) write(IO_FID_LOG,'(A,I5,A,I5)') "   Num YP =",NEST_TILE_ALL,"  Num TILE(MAX) =",NEST_TILE_ALLMAX

       call MPI_ISEND(NEST_TILE_ALLMAX, 1, MPI_INTEGER, PRC_myrank, tag+1, INTERCOMM_PARENT, ireq, ierr)
       call MPI_WAIT(ireq, istatus, ierr)

       allocate( buffer_LIST   (NEST_TILE_ALLMAX)          )
       allocate( buffer_ALLLIST(NEST_TILE_ALLMAX*PRC_nmax) )
       allocate( NEST_TILE_LIST(NEST_TILE_ALLMAX,PRC_nmax) )

       do i = 1, NEST_TILE_ALLMAX
          if ( i <= NEST_TILE_ALL ) then
             buffer_LIST(i) = NEST_TILE_ID(i)
          else
             buffer_LIST(i) = -1
          endif
       enddo

       ileng = NEST_TILE_ALLMAX
       call MPI_Allgather( buffer_LIST(:),     &
                           ileng,              &
                           MPI_INTEGER,        &
                           buffer_ALLLIST(:),  &
                           ileng,              &
                           MPI_INTEGER,        &
                           COMM_world,         &
                           ierr                )
       k = 1
       do j = 1, PRC_nmax
       do i = 1, NEST_TILE_ALLMAX
          NEST_TILE_LIST(i,j) = buffer_ALLLIST(k)
          k = k + 1
       enddo
       enddo

       deallocate( buffer_LIST    )
       deallocate( buffer_ALLLIST )

       ileng = NEST_TILE_ALLMAX*PRC_nmax
       call MPI_ISEND(NEST_TILE_LIST, ileng, MPI_INTEGER, PRC_myrank, tag+2, INTERCOMM_PARENT, ireq, ierr)
       call MPI_WAIT(ireq, istatus, ierr)

       call NEST_COMM_importgrid_nestdown

       call MPI_BARRIER(INTERCOMM_PARENT, ierr)
    else
    !---------------------------------------------------
       if( IO_L ) write(*,*) 'xxx internal error [nest/grid]'
       call PRC_MPIstop
    endif

    if( NUM_YP > 1000 .or. NEST_TILE_ALL > 1000 ) then
       write(*,*) 'xxx internal error (overflow number of ireq) [nest/grid]'
       call PRC_MPIstop
    endif

    return
  end subroutine NEST_COMM_setup_nestdown

  !-----------------------------------------------------------------------------
  !> Grid Data transfer from parent to daughter: nestdown
  subroutine NEST_COMM_importgrid_nestdown
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

    integer :: parent, daughter
    integer :: ireq(1000)  ! tentative approach
    integer :: ierr(1000)  ! tentative approach
    integer :: ileng
    integer :: istatus(MPI_STATUS_SIZE)
    integer :: tag, tagbase, target_rank

    integer :: xloc, yloc
    integer :: xs, xe
    integer :: ys, ye

    integer :: i, j, k, l
    !---------------------------------------------------------------------------

    tagbase    = INTERCOMM_ID * 100
    parent     = INTERCOMM_ID
    daughter   = INTERCOMM_ID + 1
    l          = 0

    if( NUM_YP*8 > 1000 .or. NEST_TILE_ALL*8 > 1000 ) then
       write(*,*) 'xxx internal error (overflow number of ireq) [nest/grid]'
       call PRC_MPIstop
    endif

    if ( ONLINE_DOMAIN_NUM == parent ) then
    !---------------------------------------------------
       do i = 1, NUM_YP
          ! send data to multiple daughter processes
          target_rank = NEST_TILE_LIST_YP(i)

          l = l + 1
          ileng = PARENT_IA * PARENT_JA
          tag = tagbase + tag_lon
          call MPI_ISEND(REAL_LON, ileng, COMM_datatype, target_rank, tag, INTERCOMM_DAUGHTER, ireq(l), ierr(l))
          call MPI_WAIT(ireq(l), istatus, ierr(l))

          l = l + 1
          ileng = PARENT_IA * PARENT_JA
          tag = tagbase + tag_lat
          call MPI_ISEND(REAL_LAT, ileng, COMM_datatype, target_rank, tag, INTERCOMM_DAUGHTER, ireq(l), ierr(l))
          call MPI_WAIT(ireq(l), istatus, ierr(l))

          l = l + 1
          ileng = PARENT_IA * PARENT_JA
          tag = tagbase + tag_lonx
          call MPI_ISEND(REAL_LONX, ileng, COMM_datatype, target_rank, tag, INTERCOMM_DAUGHTER, ireq(l), ierr(l))
          call MPI_WAIT(ireq(l), istatus, ierr(l))

          l = l + 1
          ileng = PARENT_IA * PARENT_JA
          tag = tagbase + tag_latx
          call MPI_ISEND(REAL_LATX, ileng, COMM_datatype, target_rank, tag, INTERCOMM_DAUGHTER, ireq(l), ierr(l))
          call MPI_WAIT(ireq(l), istatus, ierr(l))

          l = l + 1
          ileng = PARENT_IA * PARENT_JA
          tag = tagbase + tag_lony
          call MPI_ISEND(REAL_LONY, ileng, COMM_datatype, target_rank, tag, INTERCOMM_DAUGHTER, ireq(l), ierr(l))
          call MPI_WAIT(ireq(l), istatus, ierr(l))

          l = l + 1
          ileng = PARENT_IA * PARENT_JA
          tag = tagbase + tag_laty
          call MPI_ISEND(REAL_LATY, ileng, COMM_datatype, target_rank, tag, INTERCOMM_DAUGHTER, ireq(l), ierr(l))
          call MPI_WAIT(ireq(l), istatus, ierr(l))

          l = l + 1
          ileng = PARENT_KA * PARENT_IA * PARENT_JA
          tag = tagbase + tag_cz
          call MPI_ISEND(REAL_CZ, ileng, COMM_datatype, target_rank, tag, INTERCOMM_DAUGHTER, ireq(l), ierr(l))
          call MPI_WAIT(ireq(l), istatus, ierr(l))

          l = l + 1
          ileng = (PARENT_KA+1) * PARENT_IA * PARENT_JA
          tag = tagbase + tag_fz
          call MPI_ISEND(REAL_FZ, ileng, COMM_datatype, target_rank, tag, INTERCOMM_DAUGHTER, ireq(l), ierr(l))
          call MPI_WAIT(ireq(l), istatus, ierr(l))
       enddo
    elseif ( ONLINE_DOMAIN_NUM == daughter ) then
    !---------------------------------------------------
       do i = 1, NEST_TILE_ALL
          ! receive data from multiple parent tiles
          target_rank = NEST_TILE_LIST(i,PRC_myrank+1)

          xloc = mod( i-1, NEST_TILE_NUM_X ) + 1
          yloc = int( real(i-1) / real(NEST_TILE_NUM_X) ) + 1

          xs = PARENT_IA * (xloc-1) + 1
          xe = PARENT_IA * xloc
          ys = PARENT_JA * (yloc-1) + 1
          ye = PARENT_JA * yloc

          l = l + 1
          ileng = PARENT_IA * PARENT_JA
          tag = tagbase + tag_lon
          call MPI_IRECV(buffer_2D, ileng, COMM_datatype, target_rank, tag, INTERCOMM_PARENT, ireq(l), ierr(l))
          call MPI_WAIT(ireq(l), istatus, ierr(l))
          buffer_ref_LON(xs:xe,ys:ye)  = buffer_2D(:,:)

          l = l + 1
          ileng = PARENT_IA * PARENT_JA
          tag = tagbase + tag_lat
          call MPI_IRECV(buffer_2D, ileng, COMM_datatype, target_rank, tag, INTERCOMM_PARENT, ireq(l), ierr(l))
          call MPI_WAIT(ireq(l), istatus, ierr(l))
          buffer_ref_LAT(xs:xe,ys:ye)  = buffer_2D(:,:)

          l = l + 1
          ileng = PARENT_IA * PARENT_JA
          tag = tagbase + tag_lonx
          call MPI_IRECV(buffer_2D, ileng, COMM_datatype, target_rank, tag, INTERCOMM_PARENT, ireq(l), ierr(l))
          call MPI_WAIT(ireq(l), istatus, ierr(l))
          buffer_ref_LONX(xs:xe,ys:ye)  = buffer_2D(:,:)

          l = l + 1
          ileng = PARENT_IA * PARENT_JA
          tag = tagbase + tag_latx
          call MPI_IRECV(buffer_2D, ileng, COMM_datatype, target_rank, tag, INTERCOMM_PARENT, ireq(l), ierr(l))
          call MPI_WAIT(ireq(l), istatus, ierr(l))
          buffer_ref_LATX(xs:xe,ys:ye)  = buffer_2D(:,:)

          l = l + 1
          ileng = PARENT_IA * PARENT_JA
          tag = tagbase + tag_lony
          call MPI_IRECV(buffer_2D, ileng, COMM_datatype, target_rank, tag, INTERCOMM_PARENT, ireq(l), ierr(l))
          call MPI_WAIT(ireq(l), istatus, ierr(l))
          buffer_ref_LONY(xs:xe,ys:ye)  = buffer_2D(:,:)

          l = l + 1
          ileng = PARENT_IA * PARENT_JA
          tag = tagbase + tag_laty
          call MPI_IRECV(buffer_2D, ileng, COMM_datatype, target_rank, tag, INTERCOMM_PARENT, ireq(l), ierr(l))
          call MPI_WAIT(ireq(l), istatus, ierr(l))
          buffer_ref_LATY(xs:xe,ys:ye)  = buffer_2D(:,:)

          l = l + 1
          ileng = PARENT_KA * PARENT_IA * PARENT_JA
          tag = tagbase + tag_cz
          call MPI_IRECV(buffer_3D, ileng, COMM_datatype, target_rank, tag, INTERCOMM_PARENT, ireq(l), ierr(l))
          call MPI_WAIT(ireq(l), istatus, ierr(l))
          do k = 1, KA
             buffer_ref_CZ(k,xs:xe,ys:ye)  = buffer_3D(k,:,:)
          enddo

          l = l + 1
          ileng = (PARENT_KA+1) * PARENT_IA * PARENT_JA
          tag = tagbase + tag_fz
          call MPI_IRECV(buffer_3DF,ileng, COMM_datatype, target_rank, tag, INTERCOMM_PARENT, ireq(l), ierr(l))
          call MPI_WAIT(ireq(l), istatus, ierr(l))
          do k = 0, KA
             buffer_ref_FZ(k,xs:xe,ys:ye)  = buffer_3DF(k,:,:)
          enddo
       enddo
    else
    !---------------------------------------------------
       if( IO_L ) write(*,*) 'xxx internal error [nest/grid]'
       call PRC_MPIstop
    endif

    return
  end subroutine NEST_COMM_importgrid_nestdown

  !-----------------------------------------------------------------------------
  !> Boundary data transfer from parent to daughter: nestdown
  subroutine NEST_COMM_nestdown( &
      BND_QA,              & ! [in   ]
      org_DENS,            & ! [in   ]
      org_MOMX,            & ! [in   ]
      org_MOMY,            & ! [in   ]
      org_RHOT,            & ! [in   ]
      org_QTRC,            & ! [in   ]
      interped_ref_DENS,   & ! [inout]
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
       COMM_datatype
    implicit none

    integer,  intent(in)    :: BND_QA                                                        !< num of tracer

    real(RP), intent(in   ) :: org_DENS(PARENT_KA,PARENT_IA,PARENT_JA)                       !< (KMAX,IMAX,JMAX)
    real(RP), intent(in   ) :: org_MOMX(PARENT_KA,PARENT_IA,PARENT_JA)                       !< (KMAX,IMAX,JMAX)
    real(RP), intent(in   ) :: org_MOMY(PARENT_KA,PARENT_IA,PARENT_JA)                       !< (KMAX,IMAX,JMAX)
    real(RP), intent(in   ) :: org_RHOT(PARENT_KA,PARENT_IA,PARENT_JA)                       !< (KMAX,IMAX,JMAX)
    real(RP), intent(in   ) :: org_QTRC(PARENT_KA,PARENT_IA,PARENT_JA,BND_QA)                !< (KMAX,IMAX,JMAX,QA)
    real(RP), intent(inout) :: interped_ref_DENS(DAUGHTER_KA,DAUGHTER_IA,DAUGHTER_JA)        !< (KMAX,IMAX,JMAX)
    real(RP), intent(inout) :: interped_ref_VELX(DAUGHTER_KA,DAUGHTER_IA,DAUGHTER_JA)        !< (KMAX,IMAX,JMAX)
    real(RP), intent(inout) :: interped_ref_VELY(DAUGHTER_KA,DAUGHTER_IA,DAUGHTER_JA)        !< (KMAX,IMAX,JMAX)
    real(RP), intent(inout) :: interped_ref_POTT(DAUGHTER_KA,DAUGHTER_IA,DAUGHTER_JA)        !< (KMAX,IMAX,JMAX)
    real(RP), intent(inout) :: interped_ref_QTRC(DAUGHTER_KA,DAUGHTER_IA,DAUGHTER_JA,BND_QA) !< (KMAX,IMAX,JMAX,QA)

    real(RP) :: dummy(1,1,1)
    real(RP) :: dens (DAUGHTER_KA,  DAUGHTER_IA,DAUGHTER_JA)
    real(RP) :: work1(DAUGHTER_KA,  DAUGHTER_IA,DAUGHTER_JA)
    real(RP) :: work2(DAUGHTER_KA+1,DAUGHTER_IA,DAUGHTER_JA)

    integer :: parent, daughter
    integer :: ireq, ierr, ileng
    integer :: istatus(MPI_STATUS_SIZE)
    integer :: tag, tagbase
    integer :: i, j, k, iq
    !---------------------------------------------------------------------------

    if ( BND_QA > I_BNDQA ) then
       if( IO_L ) write(*,*) 'xxx internal error: about BND_QA [nest/grid]'
       call PRC_MPIstop
    endif

    tagbase    = INTERCOMM_ID * 100
    parent     = INTERCOMM_ID
    daughter   = INTERCOMM_ID + 1

    if ( ONLINE_DOMAIN_NUM == parent ) then
       call MPI_BARRIER(INTERCOMM_DAUGHTER, ierr) ![start] inter-domain communication

       tag = tagbase + tag_dens
       call NEST_COMM_intercomm_nestdown( org_DENS, dummy, tag, I_SCLR, .true. )

       tag = tagbase + tag_momx
       call NEST_COMM_intercomm_nestdown( org_MOMX, dummy, tag, I_XSTG )

       tag = tagbase + tag_momy
       call NEST_COMM_intercomm_nestdown( org_MOMY, dummy, tag, I_YSTG )

       tag = tagbase + tag_rhot
       call NEST_COMM_intercomm_nestdown( org_RHOT, dummy, tag, I_SCLR )

       do iq = 1, BND_QA
          tag = tagbase+ tag_qx + iq
          call NEST_COMM_intercomm_nestdown( org_QTRC(:,:,:,iq), dummy, tag, I_SCLR )
       enddo

       call MPI_BARRIER(INTERCOMM_DAUGHTER, ierr) ![ end ] inter-domain communication
    elseif ( ONLINE_DOMAIN_NUM == daughter ) then
       call MPI_BARRIER(INTERCOMM_PARENT, ierr) ![start] inter-domain communication

       tag = tagbase + tag_dens
       call NEST_COMM_intercomm_nestdown( dummy, dens, tag, I_SCLR, .true. )
       do j = DATR_JS, DATR_JE
       do i = DATR_IS, DATR_IE
       do k = DATR_KS, DATR_KE
       interped_ref_DENS(k,i,j) = dens(k,i,j)
       enddo
       enddo
       enddo

       tag = tagbase + tag_momx
       call NEST_COMM_intercomm_nestdown( dummy, work1, tag, I_XSTG )
       do j = DATR_JS, DATR_JE
       do i = DATR_IS, DATR_IE
       do k = DATR_KS, DATR_KE
          interped_ref_VELX(k,i,j) = work1(k,i,j) / ( dens(k,i+1,j) + dens(k,i,j) ) * 2.0_RP
       enddo
       enddo
       enddo

       tag = tagbase + tag_momy
       call NEST_COMM_intercomm_nestdown( dummy, work1, tag, I_YSTG )
       do j = DATR_JS, DATR_JE
       do i = DATR_IS, DATR_IE
       do k = DATR_KS, DATR_KE
          interped_ref_VELY(k,i,j) = work1(k,i,j) / ( dens(k,i+1,j) + dens(k,i,j) ) * 2.0_RP
       enddo
       enddo
       enddo

       tag = tagbase + tag_rhot
       call NEST_COMM_intercomm_nestdown( dummy, work1, tag, I_SCLR )
       do j = DATR_JS, DATR_JE
       do i = DATR_IS, DATR_IE
       do k = DATR_KS, DATR_KE
          interped_ref_POTT(k,i,j) = work1(k,i,j) / dens(k,i,j)
       enddo
       enddo
       enddo

       do iq = 1, BND_QA
          tag = tagbase + tag_qx + iq
          call NEST_COMM_intercomm_nestdown( dummy, work1, tag, I_SCLR )
          do j = DATR_JS, DATR_JE
          do i = DATR_IS, DATR_IE
          do k = DATR_KS, DATR_KE
             interped_ref_QTRC(k,i,j,iq) = work1(k,i,j)
          enddo
          enddo
          enddo
       enddo

       call MPI_BARRIER(INTERCOMM_PARENT, ierr) ![ end ] inter-domain communication
    else
       if( IO_L ) write(*,*) 'xxx internal error [nest/grid]'
       call PRC_MPIstop
    endif

    return
  end subroutine NEST_COMM_nestdown

  !-----------------------------------------------------------------------------
  !> Inter-communication from parent to daughter: nestdown
  subroutine NEST_COMM_intercomm_nestdown_3D( &
      pvar,      & ! [in ]
      dvar,      & ! [out]
      tag_var,   & ! [in ]
      id_stag,   & ! [in ]
      flag_dens  ) ! [in ]: optional
    use scale_process, only: &
       PRC_myrank,  &
       PRC_nmax,    &
       PRC_MPIstop
    use scale_comm, only: &
       COMM_datatype
    implicit none

    real(RP), intent(in)  :: pvar(:,:,:)          !< variable from parent domain (PARENT_KA,PARENT_IA,PARENT_JA / 1,1,1)
    real(RP), intent(out) :: dvar(:,:,:)          !< variable to daughter domain (1,1,1 / MY_KA,MY_IA,MY_JA)
    integer , intent(in)  :: tag_var              !< communication tag of the variable
    integer , intent(in)  :: id_stag              !< id of staggered grid option

    logical , intent(in), optional  :: flag_dens  !< flag of logarithmic interpolation for density

    integer :: parent, daughter
    integer :: ireq(1000)  ! tentative approach
    integer :: ierr(1000)  ! tentative approach
    integer :: ileng
    integer :: istatus(MPI_STATUS_SIZE)
    integer :: tag, target_rank

    integer :: xloc, yloc
    integer :: xs, xe
    integer :: ys, ye

    integer :: i, j, k, l, yp
    integer :: ig

    logical :: no_zstag    = .true.
    logical :: logarithmic = .false.
    !---------------------------------------------------------------------------

    parent   = INTERCOMM_ID
    daughter = INTERCOMM_ID + 1
    l        = 0

    logarithmic = .false.
    if ( present(flag_dens) .and. flag_dens ) then
       logarithmic = .true.
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
       ileng = PARENT_KA * PARENT_IA * PARENT_JA
    else
       ileng = (PARENT_KA+1) * PARENT_IA * PARENT_JA
    endif
    tag = tag_var

    if ( ONLINE_DOMAIN_NUM == parent ) then
    !---------------------------------------------------
       do yp = 1, NUM_YP
          ! send data to multiple daughter processes
          target_rank = NEST_TILE_LIST_YP(yp)

          l = l + 1
          call MPI_ISEND(pvar, ileng, COMM_datatype, target_rank, tag, INTERCOMM_DAUGHTER, ireq(l), ierr(l))
          call MPI_WAIT(ireq(l), istatus, ierr(l))

          dvar(:,:,:) = -1.0_RP  ! input as a dummy value
       enddo
    elseif ( ONLINE_DOMAIN_NUM == daughter ) then
    !---------------------------------------------------
       do yp = 1, NEST_TILE_ALL
          ! receive data from multiple parent tiles
          target_rank = NEST_TILE_LIST(yp,PRC_myrank+1)

          xloc = mod( yp-1, NEST_TILE_NUM_X ) + 1
          yloc = int( real(yp-1) / real(NEST_TILE_NUM_X) ) + 1

          xs = PARENT_IA * (xloc-1) + 1
          xe = PARENT_IA * xloc
          ys = PARENT_JA * (yloc-1) + 1
          ye = PARENT_JA * yloc

          l = l + 1
          if ( no_zstag ) then
             call MPI_IRECV(buffer_3D, ileng, COMM_datatype, target_rank, tag, INTERCOMM_PARENT, ireq(l), ierr(l))
             call MPI_WAIT(ireq(l), istatus, ierr(l))

             if ( .not. logarithmic ) then
                ! linear interpolation
                do k = 1, PARENT_KA
                   buffer_ref_3D(k,xs:xe,ys:ye)  = buffer_3D(k,1:PARENT_IA,1:PARENT_JA)
                enddo

                do j = DATR_JS-1, DATR_JE+1
                do i = DATR_IS-1, DATR_IE+1
                do k = DATR_KS-1, DATR_KE+1
                   dvar(k,i,j) = buffer_ref_3D(kgrd(k,i,j,1,1,ig),igrd(i,j,1,ig),jgrd(i,j,1,ig)) &
                               * hfact(i,j,1,ig) * vfact(k,i,j,1,1,ig)                           &
                               + buffer_ref_3D(kgrd(k,i,j,2,1,ig),igrd(i,j,2,ig),jgrd(i,j,2,ig)) &
                               * hfact(i,j,2,ig) * vfact(k,i,j,2,1,ig)                           &
                               + buffer_ref_3D(kgrd(k,i,j,3,1,ig),igrd(i,j,3,ig),jgrd(i,j,3,ig)) &
                               * hfact(i,j,3,ig) * vfact(k,i,j,3,1,ig)                           &
                               + buffer_ref_3D(kgrd(k,i,j,1,2,ig),igrd(i,j,1,ig),jgrd(i,j,1,ig)) &
                               * hfact(i,j,1,ig) * vfact(k,i,j,1,2,ig)                           &
                               + buffer_ref_3D(kgrd(k,i,j,2,2,ig),igrd(i,j,2,ig),jgrd(i,j,2,ig)) &
                               * hfact(i,j,2,ig) * vfact(k,i,j,2,2,ig)                           &
                               + buffer_ref_3D(kgrd(k,i,j,3,2,ig),igrd(i,j,3,ig),jgrd(i,j,3,ig)) &
                               * hfact(i,j,3,ig) * vfact(k,i,j,3,2,ig)
                end do
                end do
                end do
             else
                ! logarithmic weighted interpolation
                do k = 1, PARENT_KA
                   buffer_ref_3D(k,xs:xe,ys:ye)  = log( buffer_3D(k,1:PARENT_IA,1:PARENT_JA) )
                enddo

                do j = DATR_JS-1, DATR_JE+1
                do i = DATR_IS-1, DATR_IE+1
                do k = DATR_KS-1, DATR_KE+1
                   dvar(k,i,j) = exp( buffer_ref_3D(kgrd(k,i,j,1,1,ig),igrd(i,j,1,ig),jgrd(i,j,1,ig)) &
                                    * hfact(i,j,1,ig) * vfact(k,i,j,1,1,ig)                           &
                                    + buffer_ref_3D(kgrd(k,i,j,2,1,ig),igrd(i,j,2,ig),jgrd(i,j,2,ig)) &
                                    * hfact(i,j,2,ig) * vfact(k,i,j,2,1,ig)                           &
                                    + buffer_ref_3D(kgrd(k,i,j,3,1,ig),igrd(i,j,3,ig),jgrd(i,j,3,ig)) &
                                    * hfact(i,j,3,ig) * vfact(k,i,j,3,1,ig)                           &
                                    + buffer_ref_3D(kgrd(k,i,j,1,2,ig),igrd(i,j,1,ig),jgrd(i,j,1,ig)) &
                                    * hfact(i,j,1,ig) * vfact(k,i,j,1,2,ig)                           &
                                    + buffer_ref_3D(kgrd(k,i,j,2,2,ig),igrd(i,j,2,ig),jgrd(i,j,2,ig)) &
                                    * hfact(i,j,2,ig) * vfact(k,i,j,2,2,ig)                           &
                                    + buffer_ref_3D(kgrd(k,i,j,3,2,ig),igrd(i,j,3,ig),jgrd(i,j,3,ig)) &
                                    * hfact(i,j,3,ig) * vfact(k,i,j,3,2,ig) )
                end do
                end do
                end do
             endif
          else
             call MPI_IRECV(buffer_3DF, ileng, COMM_datatype, target_rank, tag, INTERCOMM_PARENT, ireq(l), ierr(l))
             call MPI_WAIT(ireq(l), istatus, ierr(l))
             do k = 0, PARENT_KA
                buffer_ref_3DF(k,xs:xe,ys:ye)  = buffer_3DF(k,1:PARENT_IA,1:PARENT_JA)
             enddo

             do j = DATR_JS-1, DATR_JE+1
             do i = DATR_IS-1, DATR_IE+1
             do k = DATR_KS-1, DATR_KE+1
                dvar(k,i,j) = buffer_ref_3DF(kgrd(k,i,j,1,1,ig),igrd(i,j,1,ig),jgrd(i,j,1,ig)) &
                            * hfact(i,j,1,ig) * vfact(k,i,j,1,1,ig)                            &
                            + buffer_ref_3DF(kgrd(k,i,j,2,1,ig),igrd(i,j,2,ig),jgrd(i,j,2,ig)) &
                            * hfact(i,j,2,ig) * vfact(k,i,j,2,1,ig)                            &
                            + buffer_ref_3DF(kgrd(k,i,j,3,1,ig),igrd(i,j,3,ig),jgrd(i,j,3,ig)) &
                            * hfact(i,j,3,ig) * vfact(k,i,j,3,1,ig)                            &
                            + buffer_ref_3DF(kgrd(k,i,j,1,2,ig),igrd(i,j,1,ig),jgrd(i,j,1,ig)) &
                            * hfact(i,j,1,ig) * vfact(k,i,j,1,2,ig)                            &
                            + buffer_ref_3DF(kgrd(k,i,j,2,2,ig),igrd(i,j,2,ig),jgrd(i,j,2,ig)) &
                            * hfact(i,j,2,ig) * vfact(k,i,j,2,2,ig)                            &
                            + buffer_ref_3DF(kgrd(k,i,j,3,2,ig),igrd(i,j,3,ig),jgrd(i,j,3,ig)) &
                            * hfact(i,j,3,ig) * vfact(k,i,j,3,2,ig)
             end do
             end do
             end do
          endif
       enddo

    else
    !---------------------------------------------------
       if( IO_L ) write(*,*) 'xxx internal error [nest/grid]'
       call PRC_MPIstop
    endif

    return
  end subroutine NEST_COMM_intercomm_nestdown_3D


  !-----------------------------------------------------------------------------
  ! WITHOUT TIME DIMENSION
  subroutine latlonz_interporation_fact( &
      hfact,      & ! (out)
      vfact,      & ! (out)
      kgrd,       & ! (out)
      igrd,       & ! (out)
      jgrd,       & ! (out)
      myhgt,      & ! (in)
      mylat,      & ! (in)
      mylon,      & ! (in)
      inhgt,      & ! (in)
      inlat,      & ! (in)
      inlon,      & ! (in)
      nz,         & ! (in)
      nx,         & ! (in)
      ny,         & ! (in)
      landgrid    ) ! (in)
    implicit none

    real(RP), intent(out) :: hfact(:,:,:)
    real(RP), intent(out) :: vfact(:,:,:,:,:)
    integer,  intent(out) :: kgrd (:,:,:,:,:)
    integer,  intent(out) :: igrd (:,:,:)
    integer,  intent(out) :: jgrd (:,:,:)

    real(RP), intent(in)  :: myhgt(:,:,:)
    real(RP), intent(in)  :: mylat(:,:)
    real(RP), intent(in)  :: mylon(:,:)
    real(RP), intent(in)  :: inhgt(:,:,:)
    real(RP), intent(in)  :: inlat(:,:)
    real(RP), intent(in)  :: inlon(:,:)
    integer,  intent(in)  :: nz
    integer,  intent(in)  :: nx
    integer,  intent(in)  :: ny

    logical,  intent(in), optional :: landgrid

    real(RP) :: distance
    real(RP) :: denom
    real(RP) :: dist(itp_nh)
    integer :: i, j, k, ii, jj, kk
    integer :: idx
    integer :: istart, iend, iinc, blk_i
    integer :: jstart, jend, jinc, blk_j
    integer :: kstart, kend
    logical :: lndgrd
    !---------------------------------------------------------------------------

    lndgrd = .false.
    if ( present(landgrid) .and. landgrid ) then
       lndgrd = .true.
    endif

    hfact(:,:,:) = 0.0_RP
    vfact(:,:,:,:,:) = 0.0_RP

    do j = DATR_JS-1, DATR_JE+1
    do i = DATR_IS-1, DATR_IE+1
       ! nearest block search
       iinc = (nx + 1) / interp_search_divnum
       jinc = (ny + 1) / interp_search_divnum
       dist(1) = large_number_one
       jj = 1 + (jinc/2)
       do while (jj <= ny)
          ii = 1 + (iinc/2)
          do while (ii <= nx)
             distance = haversine( mylat(i,j),mylon(i,j),inlat(ii,jj),inlon(ii,jj) )

             if( distance < dist(1) )then
                dist(1) = distance
                blk_i = ii
                blk_j = jj
             endif
             ii = ii + iinc
          enddo
          jj = jj + jinc
       enddo
       istart = blk_i - (iinc/2) - 1
       if( istart < 1 ) istart = 1
       iend   = blk_i + (iinc/2) + 1
       if( iend  > nx ) iend   = nx
       jstart = blk_j - (jinc/2) - 1
       if( jstart < 1 ) jstart = 1
       jend   = blk_j + (jinc/2) + 1
       if( jend  > ny ) jend   = ny

       ! main search
       dist(1) = large_number_three
       dist(2) = large_number_two
       dist(3) = large_number_one
       do jj = jstart, jend
       do ii = istart, iend
          distance = haversine( mylat(i,j),mylon(i,j),inlat(ii,jj),inlon(ii,jj) )
          if ( distance <= dist(1) ) then
             dist(3) = dist(2);     igrd(i,j,3) = igrd(i,j,2);  jgrd(i,j,3) = jgrd(i,j,2)
             dist(2) = dist(1);     igrd(i,j,2) = igrd(i,j,1);  jgrd(i,j,2) = jgrd(i,j,1)
             dist(1) = distance;    igrd(i,j,1) = ii;           jgrd(i,j,1) = jj
          elseif ( dist(1) < distance .and. distance <= dist(2) ) then
             dist(3) = dist(2);     igrd(i,j,3) = igrd(i,j,2);  jgrd(i,j,3) = jgrd(i,j,2)
             dist(2) = distance;    igrd(i,j,2) = ii;           jgrd(i,j,2) = jj
          elseif ( dist(2) < distance .and. distance <= dist(3) ) then
             dist(3) = distance;    igrd(i,j,3) = ii;           jgrd(i,j,3) = jj
          endif
       enddo
       enddo
       if( dist(1)==0.0_RP )then
          hfact(i,j,1) = 1.0_RP
          hfact(i,j,2) = 0.0_RP
          hfact(i,j,3) = 0.0_RP
       else
          denom = 1.0_RP / ( (1.0_RP/dist(1)) + (1.0_RP/dist(2)) + (1.0_RP/dist(3)) )
          hfact(i,j,1) = ( 1.0_RP/dist(1) ) * denom
          hfact(i,j,2) = ( 1.0_RP/dist(2) ) * denom
          hfact(i,j,3) = ( 1.0_RP/dist(3) ) * denom
       endif

       !if ( lndgrd ) then
       !   kstart = 1;     kend = LKMAX
       !else
          kstart = DATR_KS-1;  kend = DATR_KE+1
       !endif
       do idx = 1, itp_nh
          ii = igrd(i,j,idx)
          jj = jgrd(i,j,idx)
          do k = kstart, kend
             dist(1) = large_number_two
             dist(2) = large_number_one
             do kk = 1, nz
                distance = abs( myhgt(k,i,j) - inhgt(kk,ii,jj) )
                if ( distance <= dist(1) ) then
                   dist(2) = dist(1);     kgrd(k,i,j,idx,2) = kgrd(k,i,j,idx,1)
                   dist(1) = distance;    kgrd(k,i,j,idx,1) = kk
                elseif ( dist(1) < distance .and. distance <= dist(2) ) then
                   dist(2) = distance;    kgrd(k,i,j,idx,2) = kk
                endif
             enddo
             if( dist(1)==0.0_RP )then
                vfact(k,i,j,idx,1) = 1.0_RP
                vfact(k,i,j,idx,2) = 0.0_RP
             else
                denom = 1.0_RP / ( (1.0_RP/dist(1)) + (1.0_RP/dist(2)) )
                vfact(k,i,j,idx,1) = ( 1.0_RP/dist(1) ) * denom
                vfact(k,i,j,idx,2) = ( 1.0_RP/dist(2) ) * denom
             endif
          enddo
       enddo
    enddo
    enddo

    return
  end subroutine latlonz_interporation_fact


  !-----------------------------------------------------------------------------
  ! Haversine Formula (from R.W. Sinnott, "Virtues of the Haversine",
  ! Sky and Telescope, vol. 68, no. 2, 1984, p. 159):
  function haversine( &
      la0,       &
      lo0,       &
      la,        &
      lo )       &
      result( d )
    implicit none
    real(RP), intent(in) :: la0, lo0, la, lo   ! la,la0: Lat, lo,lo0: Lon; [rad]
    real(RP) :: d, dlon, dlat, work1, work2
    !---------------------------------------------------------------------------

    ! output unit : [m]
    dlon = lo0 - lo
    dlat = la0 - la
    work1 = (sin(dlat/2.0_RP))**2.0_RP + &
            cos(la0) * cos(la) * (sin(dlon/2.0_RP))**2.0_RP
    work2 = 2.0_RP * asin(min( 1.0_RP, sqrt(work1) ))
    d = r_in_m * work2

  end function haversine

end module scale_grid_nest
