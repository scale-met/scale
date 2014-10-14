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

  integer,  public              :: DAUGHTER_KMAX(2)     !< daughter max number in z-direction
  integer,  public              :: DAUGHTER_IMAX(2)     !< daughter max number in x-direction
  integer,  public              :: DAUGHTER_JMAX(2)     !< daughter max number in y-direction
  integer,  public              :: DAUGHTER_KA(2)       !< daughter max number in z-direction (with halo)
  integer,  public              :: DAUGHTER_IA(2)       !< daughter max number in x-direction (with halo)
  integer,  public              :: DAUGHTER_JA(2)       !< daughter max number in y-direction (with halo)
  integer,  public              :: DAUGHTER_LKMAX(2)    !< daughter max number in lz-direction
  real(DP), public              :: DAUGHTER_DTSEC(2)    !< daughter DT [sec]

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

  logical,  public              :: USE_NESTING          = .false.
  logical,  public              :: OFFLINE              = .true.
  logical,  public              :: ONLINE_IAM_PARENT    = .false.
  logical,  public              :: ONLINE_IAM_DAUGHTER  = .false.
  integer,  public              :: ONLINE_DOMAIN_NUM    = 1
  integer,  public              :: ONLINE_DAUGHTER_PRC  = 1

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

  integer, private               :: PARENT_PRC_NUM_X(2)    !< MPI processes in x-direction in parent
  integer, private               :: PARENT_PRC_NUM_Y(2)    !< MPI processes in y-direction in parent
  integer, private               :: PARENT_PRC_nmax(2)     !< MPI total processes in parent

  integer, private               :: DAUGHTER_PRC_NUM_X(2)  !< MPI processes in x-direction in daughter
  integer, private               :: DAUGHTER_PRC_NUM_Y(2)  !< MPI processes in y-direction in daughter
  integer, private               :: DAUGHTER_PRC_nmax(2)   !< MPI total processes in daughter

  integer, private               :: NEST_TILE_ALL          !< NUM of TILEs in the local node
  integer, private               :: NEST_TILE_ALLMAX_p     !< MAXNUM of TILEs among whole processes for parent
  integer, private               :: NEST_TILE_ALLMAX_d     !< MAXNUM of TILEs among whole processes for daughter
  integer, private, allocatable  :: NEST_TILE_LIST_p(:,:)  !< relationship list in whole system for parent
  integer, private, allocatable  :: NEST_TILE_LIST_d(:,:)  !< relationship list in whole system for daughter
  integer, private, allocatable  :: NEST_TILE_LIST_YP(:)   !< yellow-page of daughter targets for parent
  integer, private               :: NUM_YP                 !< page number of yellow-page

  integer, private               :: OFFLINE_PARENT_PRC_NUM_X !< MPI processes in x-direction in parent [for namelist]
  integer, private               :: OFFLINE_PARENT_PRC_NUM_Y !< MPI processes in y-direction in parent [for namelist]
  integer, private               :: OFFLINE_PARENT_KMAX    !< parent max number in z-direction [for namelist]
  integer, private               :: OFFLINE_PARENT_IMAX    !< parent max number in x-direction [for namelist]
  integer, private               :: OFFLINE_PARENT_JMAX    !< parent max number in y-direction [for namelist]
  integer, private               :: OFFLINE_PARENT_LKMAX   !< parent max number in lz-direction [for namelist]

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

  integer, private   :: INTERCOMM_ID(2)
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
    integer, allocatable :: errcodes(:)

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
       ONLINE_DAUGHTER_PRC

    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[NEST]/Categ[GRID]'

    argv(1) = 'run.conf'
    argv(2) = ''

    HANDLING_NUM = 0
    NEST_Filiation(:) = 0
    ONLINE_DAUGHTER_PRC = PRC_nmax

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

    allocate ( errcodes(1:ONLINE_DAUGHTER_PRC) )

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
            if( IO_L ) write(*,*) 'xxx cannot open latlon-catalogue file!'
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
         if( IO_L ) write(IO_FID_LOG,'(1x,A)') '*** Setup Online Nesting'

         if( ONLINE_IAM_PARENT ) then ! must do first before daughter processes
         !-------------------------------------------------
            HANDLING_NUM = 1 !HANDLING_NUM + 1
            INTERCOMM_ID(HANDLING_NUM) = ONLINE_DOMAIN_NUM
            NEST_Filiation(INTERCOMM_ID(HANDLING_NUM)) = 1

            ! Launch Daughter Domain
            write(dom_num,'(I2.2)') ONLINE_DOMAIN_NUM+1
            argv(1) = 'run.d'//dom_num//'.conf'
            if( IO_L ) write(IO_FID_LOG,'(1x,A,I2,A)') '*** Launch Daughter Domain [INTERCOMM_ID:', INTERCOMM_ID(HANDLING_NUM), ' ]'
            if( IO_L ) write(IO_FID_LOG,'(1x,A,A,A,A)') '*** Launch Command: ', trim(cmd), ' ', trim(argv(1))
            if( IO_L ) write(IO_FID_LOG,'(1x,A,I5   )') '*** Number of Daughter Processes: ', ONLINE_DAUGHTER_PRC
            call MPI_COMM_SPAWN( trim(cmd),           &
                                 argv,                &
                                 ONLINE_DAUGHTER_PRC, &
                                 MPI_INFO_NULL,       &
                                 PRC_master,          &
                                 COMM_world,          &
                                 INTERCOMM_DAUGHTER,  &
                                 errcodes,            &
                                 ierr                 )

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

            if( IO_L ) write(IO_FID_LOG,'(1x,A)'     ) '***  Informations of Parent Domain [me]'
            if( IO_L ) write(IO_FID_LOG,'(1x,A,I6)'  ) '***  --- PARENT_PRC_nmax   :', PARENT_PRC_nmax(HANDLING_NUM)
            if( IO_L ) write(IO_FID_LOG,'(1x,A,I6)'  ) '***  --- PARENT_PRC_NUM_X  :', PARENT_PRC_NUM_X(HANDLING_NUM)
            if( IO_L ) write(IO_FID_LOG,'(1x,A,I6)'  ) '***  --- PARENT_PRC_NUM_Y  :', PARENT_PRC_NUM_Y(HANDLING_NUM)
            if( IO_L ) write(IO_FID_LOG,'(1x,A,I6)'  ) '***  --- PARENT_KMAX       :', PARENT_KMAX(HANDLING_NUM)
            if( IO_L ) write(IO_FID_LOG,'(1x,A,I6)'  ) '***  --- PARENT_IMAX       :', PARENT_IMAX(HANDLING_NUM)
            if( IO_L ) write(IO_FID_LOG,'(1x,A,I6)'  ) '***  --- PARENT_JMAX       :', PARENT_JMAX(HANDLING_NUM)
            if( IO_L ) write(IO_FID_LOG,'(1x,A,F9.3)') '***  --- PARENT_DTSEC      :', PARENT_DTSEC(HANDLING_NUM)
            if( IO_L ) write(IO_FID_LOG,'(1x,A)'     ) '***  Informations of Daughter Domain'
            if( IO_L ) write(IO_FID_LOG,'(1x,A,I6)'  ) '***  --- DAUGHTER_PRC_nmax :', DAUGHTER_PRC_nmax(HANDLING_NUM)
            if( IO_L ) write(IO_FID_LOG,'(1x,A,I6)'  ) '***  --- DAUGHTER_PRC_NUM_X:', DAUGHTER_PRC_NUM_X(HANDLING_NUM)
            if( IO_L ) write(IO_FID_LOG,'(1x,A,I6)'  ) '***  --- DAUGHTER_PRC_NUM_Y:', DAUGHTER_PRC_NUM_Y(HANDLING_NUM)
            if( IO_L ) write(IO_FID_LOG,'(1x,A,I6)'  ) '***  --- DAUGHTER_KMAX     :', DAUGHTER_KMAX(HANDLING_NUM)
            if( IO_L ) write(IO_FID_LOG,'(1x,A,I6)'  ) '***  --- DAUGHTER_IMAX     :', DAUGHTER_IMAX(HANDLING_NUM)
            if( IO_L ) write(IO_FID_LOG,'(1x,A,I6)'  ) '***  --- DAUGHTER_JMAX     :', DAUGHTER_JMAX(HANDLING_NUM)
            if( IO_L ) write(IO_FID_LOG,'(1x,A,F9.3)') '***  --- DAUGHTER_DTSEC    :', DAUGHTER_DTSEC(HANDLING_NUM)

            call NEST_COMM_setup_nestdown( HANDLING_NUM )

         !---------------------------------- parent routines
         endif


         if( ONLINE_IAM_DAUGHTER ) then
         !-------------------------------------------------
            HANDLING_NUM = 2 !HANDLING_NUM + 1
            INTERCOMM_ID(HANDLING_NUM) = ONLINE_DOMAIN_NUM - 1
            NEST_Filiation(INTERCOMM_ID(HANDLING_NUM)) = -1

            call MPI_COMM_GET_PARENT( INTERCOMM_PARENT, ierr )
            if( IO_L ) write(IO_FID_LOG,'(1x,A,I2,A)') &
            '*** Activated Daughter Domain [INTERCOMM_ID:', INTERCOMM_ID(HANDLING_NUM), ' ]'

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
            if( IO_L ) write(IO_FID_LOG,'(1x,A)'     ) '***  Informations of Daughter Domain [me]'
            if( IO_L ) write(IO_FID_LOG,'(1x,A,I6)'  ) '***  --- DAUGHTER_PRC_nmax :', DAUGHTER_PRC_nmax(HANDLING_NUM)
            if( IO_L ) write(IO_FID_LOG,'(1x,A,I6)'  ) '***  --- DAUGHTER_PRC_NUM_X:', DAUGHTER_PRC_NUM_X(HANDLING_NUM)
            if( IO_L ) write(IO_FID_LOG,'(1x,A,I6)'  ) '***  --- DAUGHTER_PRC_NUM_Y:', DAUGHTER_PRC_NUM_Y(HANDLING_NUM)
            if( IO_L ) write(IO_FID_LOG,'(1x,A,I6)'  ) '***  --- DAUGHTER_KMAX     :', DAUGHTER_KMAX(HANDLING_NUM)
            if( IO_L ) write(IO_FID_LOG,'(1x,A,I6)'  ) '***  --- DAUGHTER_IMAX     :', DAUGHTER_IMAX(HANDLING_NUM)
            if( IO_L ) write(IO_FID_LOG,'(1x,A,I6)'  ) '***  --- DAUGHTER_JMAX     :', DAUGHTER_JMAX(HANDLING_NUM)
            if( IO_L ) write(IO_FID_LOG,'(1x,A,F9.3)') '***  --- DAUGHTER_DTSEC    :', DAUGHTER_DTSEC(HANDLING_NUM)
            if( IO_L ) write(IO_FID_LOG,'(1x,A)'     ) '***  Informations of Target Tiles'
            if( IO_L ) write(IO_FID_LOG,'(1x,A,I6)'  ) '***  --- TILEALL_KA      :', TILEAL_KA(HANDLING_NUM)
            if( IO_L ) write(IO_FID_LOG,'(1x,A,I6)'  ) '***  --- TILEALL_IA      :', TILEAL_IA(HANDLING_NUM)
            if( IO_L ) write(IO_FID_LOG,'(1x,A,I6)'  ) '***  --- TILEALL_JA      :', TILEAL_JA(HANDLING_NUM)

            allocate( buffer_2D  (                            PARENT_IA(HANDLING_NUM), PARENT_JA(HANDLING_NUM) ) )
            allocate( buffer_3D  (   PARENT_KA(HANDLING_NUM), PARENT_IA(HANDLING_NUM), PARENT_JA(HANDLING_NUM) ) )
            allocate( buffer_3DF ( 0:PARENT_KA(HANDLING_NUM), PARENT_IA(HANDLING_NUM), PARENT_JA(HANDLING_NUM) ) )

            allocate( buffer_ref_LON (                          TILEAL_IA(HANDLING_NUM),TILEAL_JA(HANDLING_NUM)) )
            allocate( buffer_ref_LONX(                          TILEAL_IA(HANDLING_NUM),TILEAL_JA(HANDLING_NUM)) )
            allocate( buffer_ref_LONY(                          TILEAL_IA(HANDLING_NUM),TILEAL_JA(HANDLING_NUM)) )
            allocate( buffer_ref_LAT (                          TILEAL_IA(HANDLING_NUM),TILEAL_JA(HANDLING_NUM)) )
            allocate( buffer_ref_LATX(                          TILEAL_IA(HANDLING_NUM),TILEAL_JA(HANDLING_NUM)) )
            allocate( buffer_ref_LATY(                          TILEAL_IA(HANDLING_NUM),TILEAL_JA(HANDLING_NUM)) )
            allocate( buffer_ref_CZ  (  PARENT_KA(HANDLING_NUM),TILEAL_IA(HANDLING_NUM),TILEAL_JA(HANDLING_NUM)) )
            allocate( buffer_ref_FZ  (0:PARENT_KA(HANDLING_NUM),TILEAL_IA(HANDLING_NUM),TILEAL_JA(HANDLING_NUM)) )

            allocate( buffer_ref_2D  (                          TILEAL_IA(HANDLING_NUM),TILEAL_JA(HANDLING_NUM)) )
            allocate( buffer_ref_3D  (  PARENT_KA(HANDLING_NUM),TILEAL_IA(HANDLING_NUM),TILEAL_JA(HANDLING_NUM)) )
            allocate( buffer_ref_3DF (0:PARENT_KA(HANDLING_NUM),TILEAL_IA(HANDLING_NUM),TILEAL_JA(HANDLING_NUM)) )

            allocate( hfact(            DAUGHTER_IA(HANDLING_NUM),DAUGHTER_JA(HANDLING_NUM),itp_nh,       itp_ng) )
            allocate( vfact(DAUGHTER_KA(HANDLING_NUM),DAUGHTER_IA(HANDLING_NUM),DAUGHTER_JA(HANDLING_NUM),itp_nh,itp_nv,itp_ng) )
            allocate( igrd (            DAUGHTER_IA(HANDLING_NUM),DAUGHTER_JA(HANDLING_NUM),itp_nh,       itp_ng) )
            allocate( jgrd (            DAUGHTER_IA(HANDLING_NUM),DAUGHTER_JA(HANDLING_NUM),itp_nh,       itp_ng) )
            allocate( kgrd (DAUGHTER_KA(HANDLING_NUM),DAUGHTER_IA(HANDLING_NUM),DAUGHTER_JA(HANDLING_NUM),itp_nh,itp_nv,itp_ng) )

            call NEST_COMM_setup_nestdown( HANDLING_NUM )

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
                                             TILEAL_KA(HANDLING_NUM),           & ! [IN]
                                             TILEAL_IA(HANDLING_NUM),           & ! [IN]
                                             TILEAL_JA(HANDLING_NUM),           & ! [IN]
                                             HANDLING_NUM                       ) ! [IN]

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
                                             TILEAL_KA(HANDLING_NUM)+1,         & ! [IN]
                                             TILEAL_IA(HANDLING_NUM),           & ! [IN]
                                             TILEAL_JA(HANDLING_NUM),           & ! [IN]
                                             HANDLING_NUM                       ) ! [IN]

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
                                             TILEAL_KA(HANDLING_NUM),           & ! [IN]
                                             TILEAL_IA(HANDLING_NUM),           & ! [IN]
                                             TILEAL_JA(HANDLING_NUM),           & ! [IN]
                                             HANDLING_NUM                       ) ! [IN]

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
                                             TILEAL_KA(HANDLING_NUM),           & ! [IN]
                                             TILEAL_IA(HANDLING_NUM),           & ! [IN]
                                             TILEAL_JA(HANDLING_NUM),           & ! [IN]
                                             HANDLING_NUM                       ) ! [IN]

         !---------------------------------- daughter routines
         endif

         if( IO_L ) write(IO_FID_LOG,'(1x,A,I2)') '*** Number of Related Domains :', HANDLING_NUM
         if ( HANDLING_NUM > 2 ) then
            if( IO_L ) write(*,*) 'xxx Too much handing domains (up to 2)'
            call PRC_MPIstop
         endif

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

    logical :: hit

    integer, allocatable :: pd_tile_num(:,:)

    real(RP) :: eps_lon, eps_lat
    integer  :: pd_sw_tile
    integer  :: pd_ne_tile
    integer  :: i, j, ii, jj, k
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
       eps_lon = abs((latlon_catalog(i,I_SW,I_LON) - latlon_catalog(i,I_SE,I_LON)) &
                      / real( PARENT_IMAX(HANDLE), kind=RP )) * 0.5_RP
       eps_lat = abs((latlon_catalog(i,I_SW,I_LAT) - latlon_catalog(i,I_NW,I_LAT)) &
                      / real( PARENT_JMAX(HANDLE), kind=RP )) * 0.5_RP

       if ( corner_loc(I_SW,I_LON) > latlon_catalog(i,I_SW,I_LON) .and. &
            corner_loc(I_SW,I_LAT) > latlon_catalog(i,I_SW,I_LAT) .and. &
            corner_loc(I_SW,I_LON) < latlon_catalog(i,I_NE,I_LON) .and. &
            corner_loc(I_SW,I_LAT) < latlon_catalog(i,I_NE,I_LAT) &
            .or.  &
            abs( corner_loc(I_SW,I_LON) - latlon_catalog(i,I_SW,I_LON) ) < eps_lon .and. &
            corner_loc(I_SW,I_LAT) > latlon_catalog(i,I_SW,I_LAT) .and. &
            corner_loc(I_SW,I_LAT) < latlon_catalog(i,I_NE,I_LAT) &
            .or.  &
            corner_loc(I_SW,I_LON) > latlon_catalog(i,I_SW,I_LON) .and. &
            corner_loc(I_SW,I_LON) < latlon_catalog(i,I_NE,I_LON) .and. &
            abs( corner_loc(I_SW,I_LAT) - latlon_catalog(i,I_SW,I_LAT) ) < eps_lat &
            .or.  &
            abs( corner_loc(I_SW,I_LON) - latlon_catalog(i,I_SW,I_LON) ) < eps_lon .and. &
            abs( corner_loc(I_SW,I_LAT) - latlon_catalog(i,I_SW,I_LAT) ) < eps_lat       ) then

          pd_sw_tile = i-1 ! MPI process number starts from zero
          hit = .true.
          exit ! exit loop
       endif
    enddo
    if ( .NOT. hit ) then
       write(*,*) 'xxx domain mismatch between parent and daughter: SW search'
       write(*,*) '    at rank:', PRC_myrank, ' of domain:', ONLINE_DOMAIN_NUM
       if( IO_L ) write(IO_FID_LOG,'(1x,A)') 'xxx domain mismatch between parent and daughter: SW search'
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
    do i = 1, PARENT_PRC_nmax(HANDLE)
       eps_lon = abs((latlon_catalog(i,I_NW,I_LON) - latlon_catalog(i,I_NE,I_LON)) &
                      / real( PARENT_IMAX(HANDLE), kind=RP )) * 0.5_RP
       eps_lat = abs((latlon_catalog(i,I_SE,I_LAT) - latlon_catalog(i,I_NE,I_LAT)) &
                      / real( PARENT_JMAX(HANDLE), kind=RP )) * 0.5_RP

       if ( corner_loc(I_NE,I_LON) > latlon_catalog(i,I_SW,I_LON) .and. &
            corner_loc(I_NE,I_LAT) > latlon_catalog(i,I_SW,I_LAT) .and. &
            corner_loc(I_NE,I_LON) < latlon_catalog(i,I_NE,I_LON) .and. &
            corner_loc(I_NE,I_LAT) < latlon_catalog(i,I_NE,I_LAT) &
            .or.  &
            abs( corner_loc(I_NE,I_LON) - latlon_catalog(i,I_NE,I_LON) ) < eps_lon .and. &
            corner_loc(I_NE,I_LAT) > latlon_catalog(i,I_SW,I_LAT) .and. &
            corner_loc(I_NE,I_LAT) < latlon_catalog(i,I_NE,I_LAT) &
            .or.  &
            corner_loc(I_NE,I_LON) > latlon_catalog(i,I_SW,I_LON) .and. &
            corner_loc(I_NE,I_LON) < latlon_catalog(i,I_NE,I_LON) .and. &
            abs( corner_loc(I_NE,I_LAT) - latlon_catalog(i,I_NE,I_LAT) ) < eps_lat &
            .or.  &
            abs( corner_loc(I_NE,I_LON) - latlon_catalog(i,I_NE,I_LON) ) < eps_lon .and. &
            abs( corner_loc(I_NE,I_LAT) - latlon_catalog(i,I_NE,I_LAT) ) < eps_lat       ) then

          pd_ne_tile = i-1 ! MPI process number starts from zero
          hit = .true.
          exit ! exit loop
       endif
    enddo
    if ( .NOT. hit ) then
       write(*,*) 'xxx domain mismatch between parent and daughter: NE search'
       write(*,*) '    at rank:', PRC_myrank, ' of domain:', ONLINE_DOMAIN_NUM
       if( IO_L ) write(IO_FID_LOG,'(1x,A)') 'xxx domain mismatch between parent and daughter: NE search'
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

    k = 1
    do j = 1, NEST_TILE_NUM_Y
    do i = 1, NEST_TILE_NUM_X
       NEST_TILE_ID(k) = pd_sw_tile + (i-1) + PARENT_PRC_NUM_X(HANDLE)*(j-1)
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
       PRC_nmax,    &
       PRC_NUM_X,   &
       PRC_NUM_Y,   &
       PRC_MPIstop
    use scale_time, only: &
       TIME_DTSEC
    implicit none

    integer, intent(in) :: HANDLE !< id number of nesting relation in this process target

    real(DP) :: buffer
    integer  :: datapack(12)
    integer  :: ireq1, ireq2, ierr1, ierr2, ileng
    integer  :: istatus(MPI_STATUS_SIZE)
    integer  :: tag
    !---------------------------------------------------------------------------

    tag   = INTERCOMM_ID(HANDLE) * 100
    ileng = 12

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
       buffer      = TIME_DTSEC

       call MPI_ISEND(datapack, ileng, MPI_INTEGER, PRC_myrank, tag, INTERCOMM_DAUGHTER, ireq1, ierr1)
       call MPI_ISEND(buffer, 1, MPI_DOUBLE_PRECISION, PRC_myrank, tag+1, INTERCOMM_DAUGHTER, ireq2, ierr2)
       call MPI_WAIT(ireq1, istatus, ierr1)
       call MPI_WAIT(ireq2, istatus, ierr2)

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
       PARENT_DTSEC(HANDLE)     = buffer

       ! from daughter to parent
       call MPI_IRECV(datapack, ileng, MPI_INTEGER, PRC_myrank, tag+2, INTERCOMM_DAUGHTER, ireq1, ierr1)
       call MPI_IRECV(buffer, 1, MPI_DOUBLE_PRECISION, PRC_myrank, tag+3, INTERCOMM_DAUGHTER, ireq2, ierr2)
       call MPI_WAIT(ireq1, istatus, ierr1)
       call MPI_WAIT(ireq2, istatus, ierr2)

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
       DAUGHTER_DTSEC(HANDLE)     = buffer


    elseif ( NEST_Filiation( INTERCOMM_ID(HANDLE) ) < 0 ) then !--- daughter
       ! from parent to daughter
       call MPI_IRECV(datapack, ileng, MPI_INTEGER, PRC_myrank, tag, INTERCOMM_PARENT, ireq1, ierr1)
       call MPI_IRECV(buffer, 1, MPI_DOUBLE_PRECISION, PRC_myrank, tag+1, INTERCOMM_PARENT, ireq2, ierr2)
       call MPI_WAIT(ireq1, istatus, ierr1)
       call MPI_WAIT(ireq2, istatus, ierr2)

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
       buffer       = TIME_DTSEC

       call MPI_ISEND(datapack, ileng, MPI_INTEGER, PRC_myrank, tag+2, INTERCOMM_PARENT, ireq1, ierr1)
       call MPI_ISEND(buffer, 1, MPI_DOUBLE_PRECISION, PRC_myrank, tag+3, INTERCOMM_PARENT, ireq2, ierr2)
       call MPI_WAIT(ireq1, istatus, ierr1)
       call MPI_WAIT(ireq2, istatus, ierr2)

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
       DAUGHTER_DTSEC(HANDLE)     = buffer
    else
       if( IO_L ) write(*,*) 'xxx internal error [nest/grid]'
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
       PRC_nmax,    &
       PRC_MPIstop
    use scale_grid_real, only: &
       REAL_DOMAIN_CATALOGUE
    use scale_comm, only: &
       COMM_datatype
    implicit none

    integer, intent(in) :: HANDLE !< id number of nesting relation in this process target

    integer :: ireq, ierr, ileng
    integer :: istatus(MPI_STATUS_SIZE)
    integer :: tag
    !---------------------------------------------------------------------------

    tag = INTERCOMM_ID(HANDLE) * 100

    if ( NEST_Filiation( INTERCOMM_ID(HANDLE) ) > 0 ) then !--- parent
       ileng = PRC_nmax * 4 * 2
       call MPI_ISEND(REAL_DOMAIN_CATALOGUE, ileng, COMM_datatype, PRC_myrank, tag, INTERCOMM_DAUGHTER, ireq, ierr)
       call MPI_WAIT(ireq, istatus, ierr)

    elseif ( NEST_Filiation( INTERCOMM_ID(HANDLE) ) < 0 ) then !--- daughter
       ileng = PARENT_PRC_nmax(HANDLE) * 4 * 2
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
      HANDLE  )
    use scale_process, only: &
       PRC_myrank,  &
       PRC_master,  &
       PRC_MPIstop
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
       PRC_nmax,    &
       PRC_MPIstop
    use scale_grid_real, only: &
       REAL_DOMAIN_CATALOGUE
    use scale_comm, only: &
       COMM_datatype,  &
       COMM_world
    implicit none

    integer, intent(in) :: HANDLE !< id number of nesting relation in this process target

    integer, allocatable :: buffer_LIST(:)
    integer, allocatable :: buffer_ALLLIST(:)

    integer :: ireq, ierr, ileng
    integer :: istatus(MPI_STATUS_SIZE)
    integer :: tag

    integer :: i, j, k
    !---------------------------------------------------------------------------

    tag = INTERCOMM_ID(HANDLE) * 100

    if ( NEST_Filiation( INTERCOMM_ID(HANDLE) ) > 0 ) then
    !--------------------------------------------------- parent

       call MPI_IRECV(NEST_TILE_ALLMAX_p, 1, MPI_INTEGER, PRC_myrank, tag+1, INTERCOMM_DAUGHTER, ireq, ierr)
       call MPI_WAIT(ireq, istatus, ierr)

       allocate( NEST_TILE_LIST_p (NEST_TILE_ALLMAX_p,PRC_nmax) )
       allocate( NEST_TILE_LIST_YP(NEST_TILE_ALLMAX_p*PRC_nmax) )

       ileng = NEST_TILE_ALLMAX_p*PRC_nmax
       call MPI_IRECV(NEST_TILE_LIST_p, ileng, MPI_INTEGER, PRC_myrank, tag+2, INTERCOMM_DAUGHTER, ireq, ierr)
       call MPI_WAIT(ireq, istatus, ierr)

       NEST_TILE_LIST_YP(:) = -1

       k = 0
       do j = 1, PRC_nmax
       do i = 1, NEST_TILE_ALLMAX_p
          if ( NEST_TILE_LIST_p(i,j) == PRC_myrank ) then
             k = k + 1
             NEST_TILE_LIST_YP(k) = j - 1  !rank number is started from 1
          endif
       enddo
       enddo
       NUM_YP = k

       if( IO_L ) write(IO_FID_LOG,'(A,I5,A,I5)') "   Num YP =",NUM_YP,"  Num TILE(MAX) =",NEST_TILE_ALLMAX_p

       call NEST_COMM_importgrid_nestdown( HANDLE )

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
       if( IO_L ) write(IO_FID_LOG,'(A,I5,A,I5)') "   Num YP =",NEST_TILE_ALL,"  Num TILE(MAX) =",NEST_TILE_ALLMAX_d

       call MPI_ISEND(NEST_TILE_ALLMAX_d, 1, MPI_INTEGER, PRC_myrank, tag+1, INTERCOMM_PARENT, ireq, ierr)
       call MPI_WAIT(ireq, istatus, ierr)

       allocate( buffer_LIST   (NEST_TILE_ALLMAX_d)            )
       allocate( buffer_ALLLIST(NEST_TILE_ALLMAX_d*PRC_nmax)   )
       allocate( NEST_TILE_LIST_d(NEST_TILE_ALLMAX_d,PRC_nmax) )

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
       do j = 1, PRC_nmax
       do i = 1, NEST_TILE_ALLMAX_d
          NEST_TILE_LIST_d(i,j) = buffer_ALLLIST(k)
          k = k + 1
       enddo
       enddo

       deallocate( buffer_LIST    )
       deallocate( buffer_ALLLIST )

       ileng = NEST_TILE_ALLMAX_d*PRC_nmax
       call MPI_ISEND(NEST_TILE_LIST_d, ileng, MPI_INTEGER, PRC_myrank, tag+2, INTERCOMM_PARENT, ireq, ierr)
       call MPI_WAIT(ireq, istatus, ierr)

       call NEST_COMM_importgrid_nestdown( HANDLE )

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

    tagbase = INTERCOMM_ID(HANDLE) * 100
    l       = 0

    if( NUM_YP*8 > 1000 .or. NEST_TILE_ALL*8 > 1000 ) then
       write(*,*) 'xxx internal error (overflow number of ireq) [nest/grid]'
       call PRC_MPIstop
    endif

    if ( NEST_Filiation( INTERCOMM_ID(HANDLE) ) > 0 ) then
    !--------------------------------------------------- parent
       do i = 1, NUM_YP
          ! send data to multiple daughter processes
          target_rank = NEST_TILE_LIST_YP(i)

          l = l + 1
          ileng = PARENT_IA(HANDLE) * PARENT_JA(HANDLE)
          tag = tagbase + tag_lon
          call MPI_ISEND(REAL_LON, ileng, COMM_datatype, target_rank, tag, INTERCOMM_DAUGHTER, ireq(l), ierr(l))
          call MPI_WAIT(ireq(l), istatus, ierr(l))

          l = l + 1
          ileng = PARENT_IA(HANDLE) * PARENT_JA(HANDLE)
          tag = tagbase + tag_lat
          call MPI_ISEND(REAL_LAT, ileng, COMM_datatype, target_rank, tag, INTERCOMM_DAUGHTER, ireq(l), ierr(l))
          call MPI_WAIT(ireq(l), istatus, ierr(l))

          l = l + 1
          ileng = PARENT_IA(HANDLE) * PARENT_JA(HANDLE)
          tag = tagbase + tag_lonx
          call MPI_ISEND(REAL_LONX, ileng, COMM_datatype, target_rank, tag, INTERCOMM_DAUGHTER, ireq(l), ierr(l))
          call MPI_WAIT(ireq(l), istatus, ierr(l))

          l = l + 1
          ileng = PARENT_IA(HANDLE) * PARENT_JA(HANDLE)
          tag = tagbase + tag_latx
          call MPI_ISEND(REAL_LATX, ileng, COMM_datatype, target_rank, tag, INTERCOMM_DAUGHTER, ireq(l), ierr(l))
          call MPI_WAIT(ireq(l), istatus, ierr(l))

          l = l + 1
          ileng = PARENT_IA(HANDLE) * PARENT_JA(HANDLE)
          tag = tagbase + tag_lony
          call MPI_ISEND(REAL_LONY, ileng, COMM_datatype, target_rank, tag, INTERCOMM_DAUGHTER, ireq(l), ierr(l))
          call MPI_WAIT(ireq(l), istatus, ierr(l))

          l = l + 1
          ileng = PARENT_IA(HANDLE) * PARENT_JA(HANDLE)
          tag = tagbase + tag_laty
          call MPI_ISEND(REAL_LATY, ileng, COMM_datatype, target_rank, tag, INTERCOMM_DAUGHTER, ireq(l), ierr(l))
          call MPI_WAIT(ireq(l), istatus, ierr(l))

          l = l + 1
          ileng = PARENT_KA(HANDLE) * PARENT_IA(HANDLE) * PARENT_JA(HANDLE)
          tag = tagbase + tag_cz
          call MPI_ISEND(REAL_CZ, ileng, COMM_datatype, target_rank, tag, INTERCOMM_DAUGHTER, ireq(l), ierr(l))
          call MPI_WAIT(ireq(l), istatus, ierr(l))

          l = l + 1
          ileng = (PARENT_KA(HANDLE)+1) * PARENT_IA(HANDLE) * PARENT_JA(HANDLE)
          tag = tagbase + tag_fz
          call MPI_ISEND(REAL_FZ, ileng, COMM_datatype, target_rank, tag, INTERCOMM_DAUGHTER, ireq(l), ierr(l))
          call MPI_WAIT(ireq(l), istatus, ierr(l))
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

          l = l + 1
          ileng = PARENT_IA(HANDLE) * PARENT_JA(HANDLE)
          tag = tagbase + tag_lon
          call MPI_IRECV(buffer_2D, ileng, COMM_datatype, target_rank, tag, INTERCOMM_PARENT, ireq(l), ierr(l))
          call MPI_WAIT(ireq(l), istatus, ierr(l))
          buffer_ref_LON(xs:xe,ys:ye)  = buffer_2D(PRNT_IS(HANDLE):PRNT_IE(HANDLE),PRNT_JS(HANDLE):PRNT_JE(HANDLE))

          l = l + 1
          ileng = PARENT_IA(HANDLE) * PARENT_JA(HANDLE)
          tag = tagbase + tag_lat
          call MPI_IRECV(buffer_2D, ileng, COMM_datatype, target_rank, tag, INTERCOMM_PARENT, ireq(l), ierr(l))
          call MPI_WAIT(ireq(l), istatus, ierr(l))
          buffer_ref_LAT(xs:xe,ys:ye)  = buffer_2D(PRNT_IS(HANDLE):PRNT_IE(HANDLE),PRNT_JS(HANDLE):PRNT_JE(HANDLE))

          l = l + 1
          ileng = PARENT_IA(HANDLE) * PARENT_JA(HANDLE)
          tag = tagbase + tag_lonx
          call MPI_IRECV(buffer_2D, ileng, COMM_datatype, target_rank, tag, INTERCOMM_PARENT, ireq(l), ierr(l))
          call MPI_WAIT(ireq(l), istatus, ierr(l))
          buffer_ref_LONX(xs:xe,ys:ye)  = buffer_2D(PRNT_IS(HANDLE):PRNT_IE(HANDLE),PRNT_JS(HANDLE):PRNT_JE(HANDLE))

          l = l + 1
          ileng = PARENT_IA(HANDLE) * PARENT_JA(HANDLE)
          tag = tagbase + tag_latx
          call MPI_IRECV(buffer_2D, ileng, COMM_datatype, target_rank, tag, INTERCOMM_PARENT, ireq(l), ierr(l))
          call MPI_WAIT(ireq(l), istatus, ierr(l))
          buffer_ref_LATX(xs:xe,ys:ye)  = buffer_2D(PRNT_IS(HANDLE):PRNT_IE(HANDLE),PRNT_JS(HANDLE):PRNT_JE(HANDLE))

          l = l + 1
          ileng = PARENT_IA(HANDLE) * PARENT_JA(HANDLE)
          tag = tagbase + tag_lony
          call MPI_IRECV(buffer_2D, ileng, COMM_datatype, target_rank, tag, INTERCOMM_PARENT, ireq(l), ierr(l))
          call MPI_WAIT(ireq(l), istatus, ierr(l))
          buffer_ref_LONY(xs:xe,ys:ye)  = buffer_2D(PRNT_IS(HANDLE):PRNT_IE(HANDLE),PRNT_JS(HANDLE):PRNT_JE(HANDLE))

          l = l + 1
          ileng = PARENT_IA(HANDLE) * PARENT_JA(HANDLE)
          tag = tagbase + tag_laty
          call MPI_IRECV(buffer_2D, ileng, COMM_datatype, target_rank, tag, INTERCOMM_PARENT, ireq(l), ierr(l))
          call MPI_WAIT(ireq(l), istatus, ierr(l))
          buffer_ref_LATY(xs:xe,ys:ye)  = buffer_2D(PRNT_IS(HANDLE):PRNT_IE(HANDLE),PRNT_JS(HANDLE):PRNT_JE(HANDLE))

          l = l + 1
          ileng = PARENT_KA(HANDLE) * PARENT_IA(HANDLE) * PARENT_JA(HANDLE)
          tag = tagbase + tag_cz
          call MPI_IRECV(buffer_3D, ileng, COMM_datatype, target_rank, tag, INTERCOMM_PARENT, ireq(l), ierr(l))
          call MPI_WAIT(ireq(l), istatus, ierr(l))
          do k = 1, KA
             buffer_ref_CZ(k,xs:xe,ys:ye)  = buffer_3D(k,PRNT_IS(HANDLE):PRNT_IE(HANDLE),PRNT_JS(HANDLE):PRNT_JE(HANDLE))
          enddo

          l = l + 1
          ileng = (PARENT_KA(HANDLE)+1) * PARENT_IA(HANDLE) * PARENT_JA(HANDLE)
          tag = tagbase + tag_fz
          call MPI_IRECV(buffer_3DF,ileng, COMM_datatype, target_rank, tag, INTERCOMM_PARENT, ireq(l), ierr(l))
          call MPI_WAIT(ireq(l), istatus, ierr(l))
          do k = 0, KA
             buffer_ref_FZ(k,xs:xe,ys:ye)  = buffer_3DF(k,PRNT_IS(HANDLE):PRNT_IE(HANDLE),PRNT_JS(HANDLE):PRNT_JE(HANDLE))
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
      HANDLE,              & ! [in   ]
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
    use scale_gridtrans, only: &
       rotc => GTRANS_ROTC
    implicit none

    integer,  intent(in)    :: HANDLE        !< id number of nesting relation in this process target
    integer,  intent(in)    :: BND_QA        !< num of tracer

    real(RP), intent(in   ) :: org_DENS(PARENT_KA(HANDLE),PARENT_IA(HANDLE),PARENT_JA(HANDLE))
    real(RP), intent(in   ) :: org_MOMX(PARENT_KA(HANDLE),PARENT_IA(HANDLE),PARENT_JA(HANDLE))
    real(RP), intent(in   ) :: org_MOMY(PARENT_KA(HANDLE),PARENT_IA(HANDLE),PARENT_JA(HANDLE))
    real(RP), intent(in   ) :: org_RHOT(PARENT_KA(HANDLE),PARENT_IA(HANDLE),PARENT_JA(HANDLE))
    real(RP), intent(in   ) :: org_QTRC(PARENT_KA(HANDLE),PARENT_IA(HANDLE),PARENT_JA(HANDLE),BND_QA)

    real(RP), intent(inout) :: interped_ref_DENS(DAUGHTER_KA(HANDLE),DAUGHTER_IA(HANDLE),DAUGHTER_JA(HANDLE))
    real(RP), intent(inout) :: interped_ref_VELX(DAUGHTER_KA(HANDLE),DAUGHTER_IA(HANDLE),DAUGHTER_JA(HANDLE))
    real(RP), intent(inout) :: interped_ref_VELY(DAUGHTER_KA(HANDLE),DAUGHTER_IA(HANDLE),DAUGHTER_JA(HANDLE))
    real(RP), intent(inout) :: interped_ref_POTT(DAUGHTER_KA(HANDLE),DAUGHTER_IA(HANDLE),DAUGHTER_JA(HANDLE))
    real(RP), intent(inout) :: interped_ref_QTRC(DAUGHTER_KA(HANDLE),DAUGHTER_IA(HANDLE),DAUGHTER_JA(HANDLE),BND_QA)

    real(RP) :: dummy(1,1,1)
    real(RP) :: u_llp(PARENT_KA(HANDLE),  PARENT_IA(HANDLE),  PARENT_JA(HANDLE)  )
    real(RP) :: v_llp(PARENT_KA(HANDLE),  PARENT_IA(HANDLE),  PARENT_JA(HANDLE)  )
    real(RP) :: dens (DAUGHTER_KA(HANDLE),DAUGHTER_IA(HANDLE),DAUGHTER_JA(HANDLE))
    real(RP) :: u_lld(DAUGHTER_KA(HANDLE),DAUGHTER_IA(HANDLE),DAUGHTER_JA(HANDLE))
    real(RP) :: v_lld(DAUGHTER_KA(HANDLE),DAUGHTER_IA(HANDLE),DAUGHTER_JA(HANDLE))
    real(RP) :: work1(DAUGHTER_KA(HANDLE),DAUGHTER_IA(HANDLE),DAUGHTER_JA(HANDLE))

    real(RP) :: u_on_map, v_on_map
    integer :: ireq, ierr, ileng
    integer :: istatus(MPI_STATUS_SIZE)
    integer :: tag, tagbase
    integer :: i, j, k, iq

    integer, parameter :: cosin = 1
    integer, parameter :: sine  = 2
    !---------------------------------------------------------------------------

    if ( BND_QA > I_BNDQA ) then
       if( IO_L ) write(*,*) 'xxx internal error: about BND_QA [nest/grid]'
       call PRC_MPIstop
    endif

    tagbase = INTERCOMM_ID(HANDLE) * 100

    if ( NEST_Filiation( INTERCOMM_ID(HANDLE) ) > 0 ) then
    !-------------------------------------------------------- parent
       call MPI_BARRIER(INTERCOMM_DAUGHTER, ierr) ![start] inter-domain communication

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

       tag = tagbase + tag_dens
       call NEST_COMM_intercomm_nestdown( org_DENS, dummy, tag, I_SCLR, HANDLE, .true. )

       tag = tagbase + tag_momx
       call NEST_COMM_intercomm_nestdown( u_llp, dummy, tag, I_XSTG, HANDLE )

       tag = tagbase + tag_momy
       call NEST_COMM_intercomm_nestdown( v_llp, dummy, tag, I_YSTG, HANDLE )

       tag = tagbase + tag_rhot
       call NEST_COMM_intercomm_nestdown( org_RHOT, dummy, tag, I_SCLR, HANDLE )

       do iq = 1, BND_QA
          tag = tagbase+ tag_qx + iq
          call NEST_COMM_intercomm_nestdown( org_QTRC(:,:,:,iq), dummy, tag, I_SCLR, HANDLE )
       enddo

       call MPI_BARRIER(INTERCOMM_DAUGHTER, ierr) ![ end ] inter-domain communication


    elseif ( NEST_Filiation( INTERCOMM_ID(HANDLE) ) < 0 ) then
    !-------------------------------------------------------- daughter
       call MPI_BARRIER(INTERCOMM_PARENT, ierr) ![start] inter-domain communication

       tag = tagbase + tag_dens
       call NEST_COMM_intercomm_nestdown( dummy, dens, tag, I_SCLR, HANDLE, .true. )
       do j = DATR_JS(HANDLE), DATR_JE(HANDLE)
       do i = DATR_IS(HANDLE), DATR_IE(HANDLE)
       do k = DATR_KS(HANDLE), DATR_KE(HANDLE)
       interped_ref_DENS(k,i,j) = dens(k,i,j)
       enddo
       enddo
       enddo

       tag = tagbase + tag_momx
       call NEST_COMM_intercomm_nestdown( dummy, u_lld, tag, I_XSTG, HANDLE )
       tag = tagbase + tag_momy
       call NEST_COMM_intercomm_nestdown( dummy, v_lld, tag, I_YSTG, HANDLE )
       do j = DATR_JS(HANDLE), DATR_JE(HANDLE)
       do i = DATR_IS(HANDLE), DATR_IE(HANDLE)
       do k = DATR_KS(HANDLE), DATR_KE(HANDLE)
          interped_ref_VELX(k,i,j) =   u_lld(k,i,j) * rotc(i,j,cosin) + v_lld(k,i,j) * rotc(i,j,sine )
          interped_ref_VELY(k,i,j) = - u_lld(k,i,j) * rotc(i,j,sine ) + v_lld(k,i,j) * rotc(i,j,cosin)

       enddo
       enddo
       enddo

       tag = tagbase + tag_rhot
       call NEST_COMM_intercomm_nestdown( dummy, work1, tag, I_SCLR, HANDLE )
       do j = DATR_JS(HANDLE), DATR_JE(HANDLE)
       do i = DATR_IS(HANDLE), DATR_IE(HANDLE)
       do k = DATR_KS(HANDLE), DATR_KE(HANDLE)
          interped_ref_POTT(k,i,j) = work1(k,i,j) / dens(k,i,j)
       enddo
       enddo
       enddo

       do iq = 1, BND_QA
          tag = tagbase + tag_qx + iq
          call NEST_COMM_intercomm_nestdown( dummy, work1, tag, I_SCLR, HANDLE )
          do j = DATR_JS(HANDLE), DATR_JE(HANDLE)
          do i = DATR_IS(HANDLE), DATR_IE(HANDLE)
          do k = DATR_KS(HANDLE), DATR_KE(HANDLE)
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
      HANDLE,    & ! [in ]
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
    integer , intent(in)  :: HANDLE               !< id number of nesting relation in this process target

    logical , intent(in), optional  :: flag_dens  !< flag of logarithmic interpolation for density

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

    l = 0

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
       ileng = PARENT_KA(HANDLE) * PARENT_IA(HANDLE) * PARENT_JA(HANDLE)
    else
       ileng = (PARENT_KA(HANDLE)+1) * PARENT_IA(HANDLE) * PARENT_JA(HANDLE)
    endif
    tag = tag_var


    if ( NEST_Filiation( INTERCOMM_ID(HANDLE) ) > 0 ) then
    !--------------------------------------------------- parent
       do yp = 1, NUM_YP
          ! send data to multiple daughter processes
          target_rank = NEST_TILE_LIST_YP(yp)

          l = l + 1
          call MPI_ISEND(pvar, ileng, COMM_datatype, target_rank, tag, INTERCOMM_DAUGHTER, ireq(l), ierr(l))
          call MPI_WAIT(ireq(l), istatus, ierr(l))

          dvar(:,:,:) = -1.0_RP  ! input as a dummy value
       enddo


    elseif ( NEST_Filiation( INTERCOMM_ID(HANDLE) ) < 0 ) then
    !--------------------------------------------------- daughter
       do yp = 1, NEST_TILE_ALL ! YP Loop
          ! receive data from multiple parent tiles
          target_rank = NEST_TILE_LIST_d(yp,PRC_myrank+1)

          xloc = mod( yp-1, NEST_TILE_NUM_X ) + 1
          yloc = int( real(yp-1) / real(NEST_TILE_NUM_X) ) + 1

          xs = PARENT_IMAX(HANDLE) * (xloc-1) + 1
          xe = PARENT_IMAX(HANDLE) * xloc
          ys = PARENT_JMAX(HANDLE) * (yloc-1) + 1
          ye = PARENT_JMAX(HANDLE) * yloc

          l = l + 1
          if ( no_zstag ) then
             call MPI_IRECV(buffer_3D, ileng, COMM_datatype, target_rank, tag, INTERCOMM_PARENT, ireq(l), ierr(l))
             call MPI_WAIT(ireq(l), istatus, ierr(l))

             if ( .not. logarithmic ) then
                ! linear interpolation
                do k = 1, PARENT_KA(HANDLE)
                   buffer_ref_3D(k,xs:xe,ys:ye)  = buffer_3D(k,PRNT_IS(HANDLE):PRNT_IE(HANDLE),PRNT_JS(HANDLE):PRNT_JE(HANDLE))
                enddo
             else
                ! logarithmic weighted interpolation
                do k = 1, PARENT_KA(HANDLE)
                   buffer_ref_3D(k,xs:xe,ys:ye)  = log( buffer_3D(k,PRNT_IS(HANDLE):PRNT_IE(HANDLE),PRNT_JS(HANDLE):PRNT_JE(HANDLE)) )
                enddo
             endif
          else
             call MPI_IRECV(buffer_3DF, ileng, COMM_datatype, target_rank, tag, INTERCOMM_PARENT, ireq(l), ierr(l))
             call MPI_WAIT(ireq(l), istatus, ierr(l))
             do k = 0, PARENT_KA(HANDLE)
                buffer_ref_3DF(k,xs:xe,ys:ye)  = buffer_3DF(k,PRNT_IS(HANDLE):PRNT_IE(HANDLE),PRNT_JS(HANDLE):PRNT_JE(HANDLE))
             enddo
          endif
       enddo ! YP Loop

          if ( no_zstag ) then
             if ( .not. logarithmic ) then
                ! linear interpolation
                do j = DATR_JS(HANDLE), DATR_JE(HANDLE)
                do i = DATR_IS(HANDLE), DATR_IE(HANDLE)
                do k = DATR_KS(HANDLE), DATR_KE(HANDLE)
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
                do j = DATR_JS(HANDLE), DATR_JE(HANDLE)
                do i = DATR_IS(HANDLE), DATR_IE(HANDLE)
                do k = DATR_KS(HANDLE), DATR_KE(HANDLE)
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
             ! linear interpolation (z-staggered)
             do j = DATR_JS(HANDLE), DATR_JE(HANDLE)
             do i = DATR_IS(HANDLE), DATR_IE(HANDLE)
             do k = DATR_KS(HANDLE), DATR_KE(HANDLE)
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
      HANDLE,     & ! (in)
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

    integer,  intent(in)  :: HANDLE

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

    do j = DATR_JS(HANDLE)-1, DATR_JE(HANDLE)+1
    do i = DATR_IS(HANDLE)-1, DATR_IE(HANDLE)+1
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
          kstart = DATR_KS(HANDLE);  kend = DATR_KE(HANDLE)
       !endif
       do idx = 1, itp_nh
          ii = igrd(i,j,idx)
          jj = jgrd(i,j,idx)
          do k = kstart, kend
             dist(1) = large_number_two
             dist(2) = large_number_one
             do kk = 1+KHALO, nz-KHALO
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
