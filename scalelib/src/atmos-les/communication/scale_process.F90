!-------------------------------------------------------------------------------
!> module PROCESS
!!
!! @par Description
!!          MPI/non-MPI management module
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2011-10-11 (R.Yoshida)  [new]
!! @li      2011-11-11 (H.Yashiro)  [mod] Integrate to SCALE-LES ver.3
!!
!<
module scale_process
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use mpi
  use gtool_file, only: &
     FileCloseAll
  use scale_precision
  use scale_stdio
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: PRC_MPIstart
  public :: PRC_MPIsetup
  public :: PRC_NOMPIstart
  public :: PRC_MPIstop
  public :: PRC_MPIfinish
  public :: PRC_setup
  public :: PRC_MPIsplit
  public :: PRC_MPItime
  public :: PRC_MPItimestat
  public :: PRC_BULKsetup

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  !-----------------------------------------------------------------------------
  !                          [ communicator system ]
  !    MPI_COMM_WORLD
  !          |
  ! MASTER_COMM_WORLD --split--> BULK_COMM_WORLD
  !                                     |
  !                            GLOBAL_COMM_WORLD --split--> LOCAL_COMM_WORLD
  !-----------------------------------------------------------------------------                    

  integer, public, parameter   :: PRC_master = 0      !< master node

  integer, public, parameter   :: max_depth = 1000    !< max depth of domains
  integer, public, parameter   :: split_root = 0      !< root process in each color

  integer, public              :: MASTER_COMM_WORLD   !< master (original) communicator
  integer, public              :: GLOBAL_COMM_WORLD   !< global communicator
  integer, public              :: ABORT_COMM_WORLD    !< communicator for aborting
  integer, public              :: LOCAL_COMM_WORLD    !< local communicator (split)
  integer, public              :: MASTER_myrank       !< myrank in master communicator
  integer, public              :: MASTER_nmax         !< process num in master communicator
  integer, public              :: GLOBAL_myrank       !< myrank in global communicator
  integer, public              :: GLOBAL_nmax         !< process num in global communicator
  logical, public              :: MASTER_LOG

  integer, public              :: PRC_myrank = 0   !< my node ID (Local)
  integer, public              :: PRC_nmax   = 1   !< total number of processors (Local)
  integer, public              :: PRC_NUM_X  = 1   !< x length of 2D processor topology
  integer, public              :: PRC_NUM_Y  = 1   !< y length of 2D processor topology

  integer, public, allocatable :: PRC_2Drank(:,:)  !< node index in 2D topology

  integer, public            :: PRC_next(8) = -1 !< node ID of 8 neighbour process

  integer, public, parameter :: PRC_W  = 1       !< [node direction] west
  integer, public, parameter :: PRC_N  = 2       !< [node direction] north
  integer, public, parameter :: PRC_E  = 3       !< [node direction] east
  integer, public, parameter :: PRC_S  = 4       !< [node direction] south
  integer, public, parameter :: PRC_NW = 5       !< [node direction] northwest
  integer, public, parameter :: PRC_NE = 6       !< [node direction] northeast
  integer, public, parameter :: PRC_SW = 7       !< [node direction] southwest
  integer, public, parameter :: PRC_SE = 8       !< [node direction] southeast

  logical, public            :: PRC_HAS_W
  logical, public            :: PRC_HAS_N
  logical, public            :: PRC_HAS_E
  logical, public            :: PRC_HAS_S

  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  private :: PRC_MPIcoloring
  private :: PRC_sort_ascd

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  integer, private, parameter :: abort_code   = -1 !< mpi abort code in error handler
!  integer, private, parameter :: abort_code_p = 2 !< mpi abort code in error handler from parent
!  integer, private, parameter :: abort_code_d = 3 !< mpi abort code in error handler from daughter
  logical, private :: PRC_mpi_alive   = .false. !< whether MPI is alive or not?

  integer, private :: new_abort
  integer, private :: handler_status

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Start MPI
  subroutine PRC_MPIstart
    implicit none

    integer :: ierr
    !---------------------------------------------------------------------------

    call MPI_Init(ierr)
    PRC_mpi_alive = .true.
    MASTER_COMM_WORLD = MPI_COMM_WORLD
    ABORT_COMM_WORLD  = MPI_COMM_WORLD

    call MPI_COMM_CREATE_ERRHANDLER( PRC_MPI_errorhandler, new_abort, ierr )
    call MPI_COMM_SET_ERRHANDLER   ( MPI_COMM_WORLD,       new_abort, ierr )
    call MPI_COMM_GET_ERRHANDLER   ( MPI_COMM_WORLD,  handler_status, ierr )

    call MPI_Comm_size(MPI_COMM_WORLD,MASTER_nmax,  ierr)
    call MPI_Comm_rank(MPI_COMM_WORLD,MASTER_myrank,ierr)

    if ( MASTER_myrank == PRC_master ) then
       MASTER_LOG = .true.
    else
       MASTER_LOG = .false.
    endif

    return
  end subroutine PRC_MPIstart

  !-----------------------------------------------------------------------------
  !> Setup MPI
  subroutine PRC_MPIsetup( &
      MY_COMM_WORLD   ) ! [in]
    implicit none

    integer, intent(in) :: MY_COMM_WORLD

    character(len=H_LONG)  :: fname ! name of logfile for each process
    character(len=H_SHORT) :: info

    integer :: ierr
    !---------------------------------------------------------------------------

    LOCAL_COMM_WORLD = MY_COMM_WORLD
    call MPI_COMM_RANK(LOCAL_COMM_WORLD,PRC_myrank,ierr)
    call MPI_COMM_SIZE(LOCAL_COMM_WORLD,PRC_nmax,  ierr)

    if ( .not. IO_LOG_SUPPRESS ) then
       if ( PRC_myrank == PRC_master ) then ! master node
          IO_L = .true.
       else
          IO_L = IO_LOG_ALLNODE
       endif
    endif

    if ( IO_LOG_NML_SUPPRESS ) then
       IO_LNML = .false.
    else
       IO_LNML = IO_L
    endif

    if ( IO_L ) then

       !--- Open logfile
       if ( IO_LOG_BASENAME == IO_STDOUT ) then
          IO_FID_LOG = IO_FID_STDOUT
       else
          IO_FID_LOG = IO_get_available_fid()
          call IO_make_idstr(fname,trim(IO_LOG_BASENAME),'pe',PRC_myrank)
          open( unit   = IO_FID_LOG,  &
                file   = trim(fname), &
                form   = 'formatted', &
                iostat = ierr         )
          if ( ierr /= 0 ) then
             write(*,*) 'xxx File open error! :', trim(fname)
             call PRC_MPIstop
          endif
       endif

       write(IO_FID_LOG,*) ''
       write(IO_FID_LOG,*) '                                               -+++++++++;              '
       write(IO_FID_LOG,*) '                                             ++++++++++++++=            '
       write(IO_FID_LOG,*) '                                           ++++++++++++++++++-          '
       write(IO_FID_LOG,*) '                                          +++++++++++++++++++++         '
       write(IO_FID_LOG,*) '                                        .+++++++++++++++++++++++        '
       write(IO_FID_LOG,*) '                                        +++++++++++++++++++++++++       '
       write(IO_FID_LOG,*) '                                       +++++++++++++++++++++++++++      '
       write(IO_FID_LOG,*) '                                      =++++++=x######+++++++++++++;     '
       write(IO_FID_LOG,*) '                                     .++++++X####XX####++++++++++++     '
       write(IO_FID_LOG,*) '         =+xxx=,               ++++  +++++=##+       .###++++++++++-    '
       write(IO_FID_LOG,*) '      ,xxxxxxxxxx-            +++++.+++++=##           .##++++++++++    '
       write(IO_FID_LOG,*) '     xxxxxxxxxxxxx           -+++x#;+++++#+              ##+++++++++.   '
       write(IO_FID_LOG,*) '    xxxxxxxx##xxxx,          ++++# +++++XX                #+++++++++-   '
       write(IO_FID_LOG,*) '   xxxxxxx####+xx+x         ++++#.++++++#                  #+++++++++   '
       write(IO_FID_LOG,*) '  +xxxxxX#X   #Xx#=        =+++#x=++++=#.                  x=++++++++   '
       write(IO_FID_LOG,*) '  xxxxxx#,    x###        .++++#,+++++#=                    x++++++++   '
       write(IO_FID_LOG,*) ' xxxxxx#.                 ++++# +++++x#                     #++++++++   '
       write(IO_FID_LOG,*) ' xxxxxx+                 ++++#-+++++=#                      #++++++++   '
       write(IO_FID_LOG,*) ',xxxxxX                 -+++XX-+++++#,                      +++++++++   '
       write(IO_FID_LOG,*) '=xxxxxX                .++++#.+++++#x                       -++++++++   '
       write(IO_FID_LOG,*) '+xxxxx=                ++++#.++++++#                        ++++++++#   '
       write(IO_FID_LOG,*) 'xxxxxx;               ++++#+=++++=#                         ++++++++#   '
       write(IO_FID_LOG,*) 'xxxxxxx              ,+++x#,+++++#-                        ;++++++++-   '
       write(IO_FID_LOG,*) '#xxxxxx              +++=# +++++xX                         ++++++++#    '
       write(IO_FID_LOG,*) 'xxxxxxxx            ++++#-+++++=#                         +++++++++X    '
       write(IO_FID_LOG,*) '-+xxxxxx+          ++++X#-++++=#.            -++;        =++++++++#     '
       write(IO_FID_LOG,*) ' #xxxxxxxx.      .+++++# +++++#x            =++++-      +++++++++XX     '
       write(IO_FID_LOG,*) ' #xxxxxxxxxx=--=++++++#.++++++#             ++++++    -+++++++++x#      '
       write(IO_FID_LOG,*) '  #+xxxxxxxxxx+++++++#x=++++=#              ++++++;=+++++++++++x#       '
       write(IO_FID_LOG,*) '  =#+xxxxxxxx+++++++##,+++++#=             =++++++++++++++++++##.       '
       write(IO_FID_LOG,*) '   X#xxxxxxxx++++++## +++++x#              ;x++++++++++++++++##.        '
       write(IO_FID_LOG,*) '    x##+xxxx+++++x## +++++=#                ##++++++++++++x##X          '
       write(IO_FID_LOG,*) '     ,###Xx+++x###x -++++=#,                .####x+++++X####.           '
       write(IO_FID_LOG,*) '       -########+   -#####x                   .X#########+.             '
       write(IO_FID_LOG,*) '           .,.      ......                         .,.                  '
       write(IO_FID_LOG,*) '                                                                        '
       write(IO_FID_LOG,*) '    .X#######     +###-        =###+           ###x         x########   '
       write(IO_FID_LOG,*) '   .#########    ######X      #######         ####        .#########x   '
       write(IO_FID_LOG,*) '   ####+++++=   X#######.    -#######x       .###;        ####x+++++.   '
       write(IO_FID_LOG,*) '   ###          ###= ####    #### x###       ####        -###.          '
       write(IO_FID_LOG,*) '  .###         ####   ###+  X###   ###X     =###.        ####           '
       write(IO_FID_LOG,*) '   ###-       ;###,        .###+   -###     ####        x##########+    '
       write(IO_FID_LOG,*) '   +####x     ####         ####     ####   ####         ###########.    '
       write(IO_FID_LOG,*) '    x######. =###          ###,     .###-  ###+        x###--------     '
       write(IO_FID_LOG,*) '      =##### X###         -###       #### ,###         ####             '
       write(IO_FID_LOG,*) '        .###=x###;        .###+       ###X ###X        ####.            '
       write(IO_FID_LOG,*) ' ###########; ###########+ ###########     ########### ,##########.     '
       write(IO_FID_LOG,*) '-###########  ,##########,  #########X      ##########  +#########      '
       write(IO_FID_LOG,*) ',,,,,,,,,,.     ,,,,,,,,,    .,,,,,,,.       .,,,,,,,,    ,,,,,,,,      '
       write(IO_FID_LOG,*) '                                                                        '
       write(IO_FID_LOG,*) '     SCALE : Scalable Computing by Advanced Library and Environment     '
       write(IO_FID_LOG,*) ''
       write(IO_FID_LOG,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
       write(IO_FID_LOG,*) '+                   LES-scale Numerical Weather Model                  +'
       write(IO_FID_LOG,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
       write(IO_FID_LOG,*) trim(H_LIBNAME)
       write(IO_FID_LOG,*) trim(H_MODELNAME)
       write(IO_FID_LOG,*) ''
       write(IO_FID_LOG,*) '++++++ Start MPI'
       write(IO_FID_LOG,*) '*** LOCAL_COMM_WORLD        : ', LOCAL_COMM_WORLD
       write(IO_FID_LOG,*) '*** total process  [LOCAL]  : ', PRC_nmax
       write(IO_FID_LOG,*) '*** master rank    [LOCAL]  : ', PRC_master
       write(IO_FID_LOG,*) '*** my process ID  [LOCAL]  : ', PRC_myrank
       write(IO_FID_LOG,*) '*** GLOBAL_COMM_WORLD       : ', GLOBAL_COMM_WORLD
       write(IO_FID_LOG,*) '*** total process  [GLOBAL] : ', GLOBAL_nmax
       write(IO_FID_LOG,*) '*** my process ID  [GLOBAL] : ', GLOBAL_myrank
       write(IO_FID_LOG,*) '*** MASTER_COMM_WORLD       : ', MASTER_COMM_WORLD
       write(IO_FID_LOG,*) '*** total process  [MASTER] : ', MASTER_nmax
       write(IO_FID_LOG,*) '*** my process ID  [MASTER] : ', MASTER_myrank
       write(IO_FID_LOG,*) '*** ABORT_COMM_WORLD        : ', ABORT_COMM_WORLD
       write(IO_FID_LOG,*) ''
       write(IO_FID_LOG,*) '++++++ Module[STDIO] / Categ[IO] / Origin[SCALElib]'
       write(IO_FID_LOG,*) ''
       write(IO_FID_LOG,*) '*** Open config file, FID = ', IO_FID_CONF
       write(IO_FID_LOG,*) '*** Open log    file, FID = ', IO_FID_LOG
       write(IO_FID_LOG,*) '*** basename of log file  = ', trim(IO_LOG_BASENAME)
       write(IO_FID_LOG,*) '*** detailed log output   = ', IO_LNML
       write(IO_FID_LOG,*) '*** Created Error Handler: rt=', new_abort,' st=', handler_status

    else

       if ( PRC_myrank == PRC_master ) then ! master node
          write(*,*) '*** Log report is suppressed.'
       endif

    endif

    return
  end subroutine PRC_MPIsetup

  !-----------------------------------------------------------------------------
  !> Dummy subroutine of MPIstart
  subroutine PRC_NOMPIstart
    implicit none

    character(len=H_LONG) :: fname ! name of logfile

    integer :: ierr
    !---------------------------------------------------------------------------

    PRC_nmax   = 1
    PRC_myrank = 0
    PRC_mpi_alive = .false.

    if ( .not. IO_LOG_SUPPRESS ) then
       if ( PRC_myrank == PRC_master ) then ! master node
          IO_L = .true.
       else
          IO_L = IO_LOG_ALLNODE
       endif
    endif

    if ( IO_L ) then

       !--- Open logfile
       if ( IO_LOG_BASENAME == IO_STDOUT ) then
          IO_FID_LOG = IO_FID_STDOUT
       else
          IO_FID_LOG = IO_get_available_fid()
          call IO_make_idstr(fname,'LOG','pe',PRC_myrank)
          open( unit   = IO_FID_LOG,  &
                file   = trim(fname), &
                form   = 'formatted', &
                iostat = ierr         )
          if ( ierr /= 0 ) then
             write(*,*) 'xxx File open error! :', trim(fname)
             call PRC_MPIstop
          endif
       endif

       write(IO_FID_LOG,*)
       write(IO_FID_LOG,*) '     SCALE : Scalable Computing by Advanced Library and Environment     '
       write(IO_FID_LOG,*)
       write(IO_FID_LOG,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
       write(IO_FID_LOG,*) '+                   LES-scale Numerical Weather Model                  +'
       write(IO_FID_LOG,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
       write(IO_FID_LOG,*) trim(H_LIBNAME)
       write(IO_FID_LOG,*) trim(H_MODELNAME)
       write(IO_FID_LOG,*)
       write(IO_FID_LOG,*) '++++++ Start WITHOUT MPI'
       write(IO_FID_LOG,*)
       write(IO_FID_LOG,*) '++++++ Module[STDIO] / Categ[IO] / Origin[SCALElib]'
       write(IO_FID_LOG,*)
       write(IO_FID_LOG,*) '*** Open config file, FID = ', IO_FID_CONF
       write(IO_FID_LOG,*) '*** Open log    file, FID = ', IO_FID_LOG
       write(IO_FID_LOG,*) '*** detailed log output   = ', IO_LNML

    endif

    return
  end subroutine PRC_NOMPIstart

  !-----------------------------------------------------------------------------
  !> Setup BULK JOB
  subroutine PRC_BULKsetup( &
      BULK_COMM_WORLD,   & ! [in]
      ABORT_ALL_JOBS     ) ! [in]
    implicit none

    integer, intent(in) :: BULK_COMM_WORLD
    logical, intent(in) :: ABORT_ALL_JOBS

    integer :: ierr
    !---------------------------------------------------------------------------


    GLOBAL_COMM_WORLD = BULK_COMM_WORLD
    call MPI_Comm_size(GLOBAL_COMM_WORLD,GLOBAL_nmax,  ierr)
    call MPI_Comm_rank(GLOBAL_COMM_WORLD,GLOBAL_myrank,ierr)

    if ( ABORT_ALL_JOBS ) then
       ABORT_COMM_WORLD = MPI_COMM_WORLD
    else
       ABORT_COMM_WORLD = GLOBAL_COMM_WORLD
    endif

    return
  end subroutine PRC_BULKsetup

  !-----------------------------------------------------------------------------
  !> Abort MPI
  subroutine PRC_MPIstop
    implicit none

    integer :: ierr
    !---------------------------------------------------------------------------

    if ( IO_L ) then
       write(IO_FID_LOG,*) ''
       write(IO_FID_LOG,*) '++++++ PRC_MPIstop', PRC_myrank
       write(IO_FID_LOG,*) ''
    end if

    write(*,*) ''
    write(*,*) '++++++ PRC_MPIstop', PRC_myrank
    write(*,*) ''

    if ( PRC_mpi_alive ) then
        call MPI_COMM_CALL_ERRHANDLER(MPI_COMM_WORLD, abort_code, ierr)
    endif

    stop
  end subroutine PRC_MPIstop

  !-----------------------------------------------------------------------------
  !> Stop MPI peacefully
  subroutine PRC_MPIfinish
    implicit none

    character(len=H_SHORT) :: request = 'STOP'

    integer :: ierr
    !---------------------------------------------------------------------------

    ! Stop MPI
    if ( PRC_mpi_alive ) then
       if ( IO_L ) then
          write(IO_FID_LOG,*)
          write(IO_FID_LOG,*) '++++++ Stop MPI'
          write(IO_FID_LOG,*) '*** Broadcast STOP signal'
       endif

       ! free splitted communicator
       if ( LOCAL_COMM_WORLD .ne. GLOBAL_COMM_WORLD ) then
          call MPI_COMM_FREE( LOCAL_COMM_WORLD, ierr )
       endif

       call MPI_BCAST( request,        &
                       H_SHORT,      &
                       MPI_CHARACTER,  & !--- type
                       PRC_master,     & !--- source rank
                       MPI_COMM_WORLD, &
                       ierr            )
       call MPI_Barrier(MPI_COMM_WORLD,ierr)
       call MPI_Finalize(ierr)
       if( IO_L ) write(IO_FID_LOG,*) '*** MPI is peacefully finalized'
    endif

    ! Close logfile, configfile
    if ( IO_L ) then
       if ( IO_FID_LOG /= IO_FID_STDOUT ) close(IO_FID_LOG)
    endif
    close(IO_FID_CONF)

    ! Stop program
    stop
  end subroutine PRC_MPIfinish

  !-----------------------------------------------------------------------------
  !> Setup Processor topology
  subroutine PRC_setup
    implicit none

    logical :: PRC_PERIODIC_X   = .true.  !< periodic condition or not (X)?
    logical :: PRC_PERIODIC_Y   = .true.  !< periodic condition or not (Y)?
    logical :: PRC_CART_REORDER = .false. !< flag for rank reordering over the cartesian map

    namelist / PARAM_PRC / &
       PRC_NUM_X,      &
       PRC_NUM_Y,      &
       PRC_PERIODIC_X, &
       PRC_PERIODIC_Y, &
       PRC_CART_REORDER

    logical :: period(2)
    integer :: divide(2)
    integer :: coords_W(2)
    integer :: coords_N(2)
    integer :: coords_E(2)
    integer :: coords_S(2)
    integer :: next_coords(2)
    integer :: iptbl
    integer :: next(8)

    integer :: ierr
    integer :: p
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '++++++ Module[PROCESS] / Categ[ATMOS-LES COMM] / Origin[SCALElib]'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_PRC,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_PRC. Check!'
       call PRC_MPIstop
    endif
    if( IO_LNML ) write(IO_FID_LOG,nml=PARAM_PRC)

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** Process Allocation ***'
    if( IO_L ) write(IO_FID_LOG,*) '*** No. of Node   :', PRC_NUM_X," x ",PRC_NUM_Y

    if ( PRC_NUM_X*PRC_NUM_Y /= PRC_nmax ) then
       write(*,*) 'xxx total number of node does not match that requested. Check!'
       call PRC_MPIstop
    endif

    if ( mod(PRC_nmax,PRC_NUM_X) /= 0 ) then
       write(*,*) 'xxx number of requested node cannot devide to 2D. Check!'
       call PRC_MPIstop
    endif

    ! set communication topology
    allocate( PRC_2Drank(-1:PRC_nmax-1,2) )
    PRC_2Drank(:,:) = -1

    do p = 0, PRC_nmax-1
       PRC_2Drank(p,1) = mod(p,PRC_NUM_X)
       PRC_2Drank(p,2) = (p-PRC_2Drank(p,1)) / PRC_NUM_X
    enddo

    divide(1) = PRC_NUM_Y
    divide(2) = PRC_NUM_X
    period(1) = PRC_PERIODIC_Y
    period(2) = PRC_PERIODIC_X
    if ( PRC_mpi_alive ) then
       call MPI_CART_CREATE(LOCAL_COMM_WORLD,2,divide,period,PRC_CART_REORDER,iptbl,ierr)
       call MPI_CART_SHIFT(iptbl,0,1,PRC_next(PRC_S),PRC_next(PRC_N),ierr) ! next rank search Down/Up
       call MPI_CART_SHIFT(iptbl,1,1,PRC_next(PRC_W),PRC_next(PRC_E),ierr) ! next rank search Left/Right

       ! get neighbor_coordinates
       PRC_HAS_W = PRC_next(PRC_W) /= MPI_PROC_NULL
       if( PRC_HAS_W ) call MPI_CART_COORDS(iptbl,PRC_next(PRC_W),2,coords_W,ierr)
       PRC_HAS_N = PRC_next(PRC_N) /= MPI_PROC_NULL
       if( PRC_HAS_N ) call MPI_CART_COORDS(iptbl,PRC_next(PRC_N),2,coords_N,ierr)
       PRC_HAS_E = PRC_next(PRC_E) /= MPI_PROC_NULL
       if( PRC_HAS_E ) call MPI_CART_COORDS(iptbl,PRC_next(PRC_E),2,coords_E,ierr)
       PRC_HAS_S = PRC_next(PRC_S) /= MPI_PROC_NULL
       if( PRC_HAS_S ) call MPI_CART_COORDS(iptbl,PRC_next(PRC_S),2,coords_S,ierr)
       ! next rank search NorthWest
       if (      .not. PRC_HAS_N &
            .OR. .not. PRC_HAS_W ) then
          PRC_next(PRC_NW) = MPI_PROC_NULL
       else
          next_coords(1) = coords_N(1)
          next_coords(2) = coords_W(2)
          call MPI_CART_RANK(iptbl, next_coords, PRC_next(PRC_NW), ierr)
       endif
       ! next rank search NorthEast
       if (      .not. PRC_HAS_N &
            .OR. .not. PRC_HAS_E ) then
          PRC_next(PRC_NE) = MPI_PROC_NULL
       else
          next_coords(1) = coords_N(1)
          next_coords(2) = coords_E(2)
          call MPI_CART_RANK(iptbl, next_coords, PRC_next(PRC_NE), ierr)
       endif
       ! next rank search SouthWest
       if (      .not. PRC_HAS_S &
            .OR. .not. PRC_HAS_W ) then
          PRC_next(PRC_SW) = MPI_PROC_NULL
       else
          next_coords(1) = coords_S(1)
          next_coords(2) = coords_W(2)
          call MPI_CART_RANK(iptbl, next_coords, PRC_next(PRC_SW), ierr)
       endif
       ! next rank search SouthEast
       if (      .not. PRC_HAS_S &
            .OR. .not. PRC_HAS_E ) then
          PRC_next(PRC_SE) = MPI_PROC_NULL
       else
          next_coords(1) = coords_S(1)
          next_coords(2) = coords_E(2)
          call MPI_CART_RANK(iptbl, next_coords, PRC_next(PRC_SE), ierr)
       endif
    endif

    next(:) = max(PRC_next(:),-1) ! avoid if MPI_PROC_NULL < -1

    if( IO_L ) write(IO_FID_LOG,*) '*** Node topology :'
    if( IO_L ) write(IO_FID_LOG,'(1x,A,I5,A,I5,A,I5,A,A,I5,A,I5,A,I5,A,A,I5,A,I5,A,I5,A)') &
    '***  NW(',next(PRC_NW),',',PRC_2Drank(next(PRC_NW),1),',',PRC_2Drank(next(PRC_NW),2),')', &
      ' -  N(',next(PRC_N) ,',',PRC_2Drank(next(PRC_N) ,1),',',PRC_2Drank(next(PRC_N) ,2),')', &
      ' - NE(',next(PRC_NE),',',PRC_2Drank(next(PRC_NE),1),',',PRC_2Drank(next(PRC_NE),2),')'
    if( IO_L ) write(IO_FID_LOG,'(1x,A)') '***                                  |'
    if( IO_L ) write(IO_FID_LOG,'(1x,A,I5,A,I5,A,I5,A,A,I5,A,I5,A,I5,A,A,I5,A,I5,A,I5,A)') &
    '***   W(',next(PRC_W),',',PRC_2Drank(next(PRC_W),1),',',PRC_2Drank(next(PRC_W),2),')', &
      ' -  P(',PRC_myrank ,',',PRC_2Drank(PRC_myrank, 1),',',PRC_2Drank(PRC_myrank, 2),')', &
      ' -  E(',next(PRC_E),',',PRC_2Drank(next(PRC_E),1),',',PRC_2Drank(next(PRC_E),2),')'
    if( IO_L ) write(IO_FID_LOG,'(1x,A)') '***                                  |'
    if( IO_L ) write(IO_FID_LOG,'(1x,A,I5,A,I5,A,I5,A,A,I5,A,I5,A,I5,A,A,I5,A,I5,A,I5,A)') &
    '***  SW(',next(PRC_SW),',',PRC_2Drank(next(PRC_SW),1),',',PRC_2Drank(next(PRC_SW),2),')', &
      ' -  S(',next(PRC_S) ,',',PRC_2Drank(next(PRC_S) ,1),',',PRC_2Drank(next(PRC_S) ,2),')', &
      ' - SE(',next(PRC_SE),',',PRC_2Drank(next(PRC_SE),1),',',PRC_2Drank(next(PRC_SE),2),')'

    return
  end subroutine PRC_setup

  !-----------------------------------------------------------------------------
  !> MPI Communicator Split
  subroutine PRC_MPIsplit( &
      ORG_COMM,         & ! [in ]
      NUM_DOMAIN,       & ! [in ]
      PRC_DOMAINS,      & ! [in ]
      CONF_FILES,       & ! [in ]
      LOG_SPLIT,        & ! [in ]
      bulk_split,       & ! [in ]
      INTRA_COMM,       & ! [out]
      inter_parent,     & ! [out]
      inter_child,      & ! [out]
      fname_local       ) ! [out]
    implicit none

    integer, intent(in)  :: ORG_COMM
    integer, intent(in)  :: NUM_DOMAIN
    integer, intent(in)  :: PRC_DOMAINS(:)
    character(len=H_LONG), intent(in) :: CONF_FILES(:)
    logical, intent(in)  :: LOG_SPLIT
    logical, intent(in)  :: bulk_split

    integer, intent(out) :: intra_comm
    integer, intent(out) :: inter_parent
    integer, intent(out) :: inter_child
    character(len=H_LONG), intent(out) :: fname_local

    integer :: PARENT_COL(max_depth)        ! parent color number
    integer :: CHILD_COL(max_depth)         ! child  color number
    integer :: PRC_ROOT(0:max_depth)        ! root process in the color
    integer, allocatable :: COLOR_LIST(:)   ! member list in each color
    integer, allocatable :: KEY_LIST(:)     ! local process number in each color

    integer :: total_nmax
    integer :: ORG_myrank  ! my rank number in the original communicator
    integer :: ORG_nmax    ! total rank number in the original communicator

    logical :: do_create_p(max_depth)
    logical :: do_create_c(max_depth)
    logical :: color_reorder = .false.

    integer :: COL_NMAX(0:max_depth)
    character(len=H_LONG) :: COL_FILE(0:max_depth)
    character(4) :: col_num

    integer :: i
    integer :: itag, ierr
    !---------------------------------------------------------------------------

    INTRA_COMM       = ORG_COMM
    inter_parent     = MPI_COMM_NULL
    inter_child      = MPI_COMM_NULL
    fname_local      = CONF_FILES(1)

    if ( NUM_DOMAIN > 1 ) then ! multi domain run
       call MPI_COMM_RANK(ORG_COMM,ORG_myrank,ierr)
       call MPI_COMM_SIZE(ORG_COMM,ORG_nmax,  ierr)
       allocate ( COLOR_LIST(0:ORG_nmax-1) )
       allocate ( KEY_LIST  (0:ORG_nmax-1) )

       total_nmax = 0
       do i = 1, NUM_DOMAIN
          total_nmax = total_nmax + PRC_DOMAINS(i)
       enddo
       if ( total_nmax .ne. ORG_nmax ) then
          if ( MASTER_LOG ) write (*,*) ""
          if ( MASTER_LOG ) write (*,*) "ERROR: MPI PROCESS NUMBER is INCONSISTENT"
          if ( MASTER_LOG ) write (*,*) "REQUESTED NPROCS = ", total_nmax, "  LAUNCHED NPROCS = ", ORG_nmax
          call PRC_MPIstop
       endif

       if ( bulk_split ) then
          color_reorder = .false.
       else
          color_reorder = .true.
       endif
       call PRC_MPIcoloring( ORG_COMM,     & ! [in ]
                             NUM_DOMAIN,   & ! [in ]
                             PRC_DOMAINS,  & ! [in ]
                             CONF_FILES,   & ! [in ]
                             color_reorder, & ! [in ]
                             LOG_SPLIT,     & ! [in ] 
                             COLOR_LIST,   & ! [out]
                             PRC_ROOT,     & ! [out]
                             KEY_LIST,     & ! [out]
                             PARENT_COL,   & ! [out]
                             CHILD_COL,    & ! [out]
                             COL_FILE      ) ! [out]


       ! split comm_world
       call MPI_COMM_SPLIT(ORG_COMM,               &
                           COLOR_LIST(ORG_myrank), &
                           KEY_LIST(ORG_myrank),   &
                           INTRA_COMM, ierr)
       if ( bulk_split ) then
          write(col_num,'(I4.4)') COLOR_LIST(ORG_myrank)
          fname_local = col_num
       else
          fname_local = COL_FILE(COLOR_LIST(ORG_myrank))
       endif

       ! set parent-child relationship
       do_create_p(:) = .false.
       do_create_c(:) = .false.
       if ( .NOT. bulk_split ) then
          do i = 1, NUM_DOMAIN-1
             if ( MASTER_LOG ) write ( *, '(1X,A,I4)' ) "relationship: ", i
             if ( MASTER_LOG ) write ( *, '(1X,A,I4,A,I4)' ) &
                               "--- parent color = ", PARENT_COL(i), "  child color = ", CHILD_COL(i)
             if ( COLOR_LIST(ORG_myrank) == PARENT_COL(i) ) then
                do_create_p(i) = .true.
             elseif ( COLOR_LIST(ORG_myrank) == CHILD_COL(i) ) then
                do_create_c(i) = .true.
             endif
          enddo
       endif

       ! create inter-commnunicator
       inter_parent = MPI_COMM_NULL
       inter_child  = MPI_COMM_NULL
       if ( .NOT. bulk_split ) then
          do i = 1, NUM_DOMAIN-1
             itag = i*100
             if ( do_create_p(i) ) then ! as a parent
                call MPI_INTERCOMM_CREATE( INTRA_COMM, split_root,           &
                                           ORG_COMM,   PRC_ROOT(CHILD_COL(i)), &
                                           itag, inter_child,  ierr)
             elseif ( do_create_c(i) ) then ! as a child
                call MPI_INTERCOMM_CREATE( INTRA_COMM, split_root,           &
                                           ORG_COMM,   PRC_ROOT(PARENT_COL(i)), &
                                           itag, inter_parent, ierr)
             endif
             call MPI_BARRIER(ORG_COMM, ierr)
          enddo
       endif

       deallocate ( COLOR_LIST, KEY_LIST )

    elseif ( NUM_DOMAIN == 1 ) then ! single domain run
       if ( MASTER_LOG ) write (*,*) "*** a single comunicator"
    else
       if ( MASTER_LOG ) write (*,*) "ERROR: REQUESTED DOMAIN NUMBER IS NOT ACCEPTABLE"
       call PRC_MPIstop
    endif

    return
  end subroutine PRC_MPIsplit

  !-----------------------------------------------------------------------------
  !> Set color and keys for COMM_SPLIT
  subroutine PRC_MPIcoloring( &
      ORG_COMM,        & ! [in ]
      NUM_DOMAIN,      & ! [in ]
      PRC_DOMAINS,     & ! [in ]
      CONF_FILES,      & ! [in ]
      color_reorder,   & ! [in ]
      LOG_SPLIT,       & ! [in ]
      COLOR_LIST,      & ! [out]
      PRC_ROOT,        & ! [out]
      KEY_LIST,        & ! [out]
      PARENT_COL,      & ! [out]
      CHILD_COL,       & ! [out]
      COL_FILE         ) ! [out]
    implicit none

    integer, intent(in)  :: ORG_COMM
    integer, intent(in)  :: NUM_DOMAIN
    integer, intent(in)  :: PRC_DOMAINS(:)
    character(len=H_LONG), intent(in) :: CONF_FILES(:)
    logical, intent(in)  :: color_reorder
    logical, intent(in)  :: LOG_SPLIT
    integer, intent(out) :: COLOR_LIST(:)             ! member list in each color
    integer, intent(out) :: PRC_ROOT(0:max_depth)     ! root process in each color
    integer, intent(out) :: KEY_LIST(:)               ! local process number in each color
    integer, intent(out) :: PARENT_COL(:)             ! parent color number
    integer, intent(out) :: CHILD_COL(:)              ! child  color number
    character(len=H_LONG), intent(out) :: COL_FILE(0:max_depth) ! conf file in each color

    integer :: touch(max_depth)
    integer :: PRC_ORDER(max_depth)                   ! reordered number of process
    integer :: ORDER2DOM(max_depth)                   ! get domain number by order number
    integer :: DOM2ORDER(max_depth)                   ! get order number by domain number
    integer :: DOM2COL(max_depth)                     ! get color number by domain number
    integer :: COL2DOM(0:max_depth)                   ! get domain number by color number
    integer :: RO_PRC_DOMAINS(max_depth)              ! reordered values
    integer :: RO_DOM2COL(max_depth)                  ! reordered values
    integer :: RO_PARENT_COL(max_depth)               ! reordered values
    integer :: RO_CHILD_COL(max_depth)                ! reordered values
    character(len=H_LONG) :: RO_CONF_FILES(max_depth) ! reordered values

    integer :: ORG_nmax   ! parent domain number
    integer :: id_parent  ! parent domain number
    integer :: id_child   ! child domain number
    integer :: dnum, nprc, order, key
    integer :: i, j
    integer :: ierr
    !---------------------------------------------------------------------------

    call MPI_COMM_SIZE(ORG_COMM,ORG_nmax, ierr)

    if ( color_reorder ) then
       !--- make color order
       !    domain num is counted from 1
       !    color num  is counted from 0
       touch(:) = -1
       PRC_ORDER(:) = PRC_DOMAINS(:)
       call PRC_sort_ascd( PRC_ORDER(1:NUM_DOMAIN), 1, NUM_DOMAIN )
   
       do i = 1, NUM_DOMAIN
       do j = 1, NUM_DOMAIN
          if ( PRC_DOMAINS(i) .eq. PRC_ORDER(j) .and. touch(j) < 0 ) then
             DOM2COL(i         ) = j - 1  ! domain_num --> color_num
             COL2DOM(DOM2COL(i)) = i      ! color_num  --> domain_num
             touch(j) = 1
             exit
          endif
       enddo
       enddo
   
       PARENT_COL(:) = -1
       CHILD_COL(:)  = -1
       do i = 1, NUM_DOMAIN
          id_parent = i - 1
          id_child  = i + 1
   
          if ( 1 <= id_parent .and. id_parent <= NUM_DOMAIN ) then
             PARENT_COL(i) = DOM2COL(id_parent)
          endif
          if ( 1 <= id_child  .and. id_child  <= NUM_DOMAIN ) then
             CHILD_COL(i) = DOM2COL(id_child)
          endif
   
          if ( MASTER_LOG .and. LOG_SPLIT ) then
             write( *, '(1X,A,I2,1X,A,I2,2(2X,A,I2))' )  &
                  "DOMAIN: ", i, "MY_COL: ", DOM2COL(i), &
                  "PARENT: COL= ", PARENT_COL(i), "CHILD: COL= ", CHILD_COL(i)
          endif
       enddo
     
       !--- reorder following color order
       do i = 1, NUM_DOMAIN
          dnum = COL2DOM(i-1)
          ORDER2DOM(i)      = dnum
          DOM2ORDER(dnum)   = i
          RO_PRC_DOMAINS(i) = PRC_DOMAINS(dnum)
          RO_DOM2COL(dnum)  = DOM2COL(dnum)
          RO_CONF_FILES(i)  = CONF_FILES(dnum)
          RO_PARENT_COL(i)  = PARENT_COL(dnum)
          RO_CHILD_COL(i)   = CHILD_COL (dnum)
       enddo
     
       !--- set relationship by ordering of relationship number
       PARENT_COL(:)   = -1
       CHILD_COL(:)    = -1
       do i = 1, NUM_DOMAIN-1
          PARENT_COL(i) = RO_PARENT_COL( DOM2ORDER(i+1) ) ! from child to parent
          CHILD_COL(i)  = RO_CHILD_COL ( DOM2ORDER(i)   ) ! from parent to child
       enddo
     
       do i = 1, NUM_DOMAIN
          if ( MASTER_LOG ) write ( *, * ) ""
          if ( MASTER_LOG ) write ( *, '(1X,A,I2,A,I5)' ) "ORDER (",i,") -> DOMAIN: ", ORDER2DOM(i)
          if ( MASTER_LOG ) write ( *, '(1X,A,I1,A,I5)' ) "NUM PRC_DOMAINS(",i,")  = ", RO_PRC_DOMAINS(i)
          if ( MASTER_LOG ) write ( *, '(1X,A,I1,A,I3)' ) "MY COLOR(",i,") = ", RO_DOM2COL(ORDER2DOM(i))
          if ( MASTER_LOG ) write ( *, '(1X,A,I1,A,I3)' ) "PARENT COLOR(",i,") = ", RO_PARENT_COL(i)
          if ( MASTER_LOG ) write ( *, '(1X,A,I1,A,I3)' ) "CHILD COLOR(",i,") = ", RO_CHILD_COL(i)
          if ( MASTER_LOG ) write ( *, '(1X,A,I1,A,A)'  ) "CONF_FILES(",i,")    = ", trim(RO_CONF_FILES(i))
       enddo
       if ( MASTER_LOG ) write ( *, * ) ""
   
       do i=1, max_depth
          COL_FILE(i-1) = RO_CONF_FILES(i) ! final copy
       enddo

    else !--- without reordering of colors
       ORDER2DOM(:)      = -1
       RO_DOM2COL(:)     = -1
       RO_PRC_DOMAINS(:) = -1
       RO_PRC_DOMAINS(:) = -1
       RO_PARENT_COL(:)  = -1
       RO_CHILD_COL(:)   = -1
   
       do i = 1, NUM_DOMAIN
          ORDER2DOM(i)      = i
          RO_DOM2COL(i)     = i-1
          RO_PRC_DOMAINS(i) = PRC_DOMAINS(i)
          RO_CONF_FILES(i)  = CONF_FILES(i)
       enddo
   
       do i = 1, NUM_DOMAIN
          id_parent = i - 1
          id_child  = i + 1
   
          if ( 1 <= id_parent .and. id_parent <= NUM_DOMAIN ) then
             RO_PARENT_COL(i) = RO_DOM2COL(id_parent)
          endif
          if ( 1 <= id_child  .and. id_child  <= NUM_DOMAIN ) then
             RO_CHILD_COL(i) = RO_DOM2COL(id_child)
          endif
       enddo
   
       ! make relationship
       do i = 1, NUM_DOMAIN-1
          PARENT_COL(i) = RO_PARENT_COL(i+1) ! from child to parent
          CHILD_COL(i)  = RO_CHILD_COL (i  ) ! from parent to child
       enddo

    endif

    ! make a process table
    order = 1
    key   = 0
    nprc  = RO_PRC_DOMAINS(order)
    PRC_ROOT(:) = -999

    do i = 0, ORG_nmax-1
       COLOR_LIST(i+1) = RO_DOM2COL(ORDER2DOM(order))
       KEY_LIST(i+1)   = key
       if ( key == 0 ) then
          PRC_ROOT(COLOR_LIST(i+1)) = i
          COL_FILE(COLOR_LIST(i+1)) = RO_CONF_FILES(order)
       endif
       if ( LOG_SPLIT .and. MASTER_LOG ) then
          write ( *, '(1X,4(A,I5))' ) "PE:", i, "   COLOR:", COLOR_LIST(i+1), &
                "   KEY:", KEY_LIST(i+1), "   PRC_ROOT:", PRC_ROOT(COLOR_LIST(i+1))
       endif
       key = key + 1
       if ( key >= nprc ) then
          order = order + 1
          key   = 0
          nprc  = RO_PRC_DOMAINS(order)
       endif
    enddo

    return
  end subroutine PRC_MPIcoloring
  !-----------------------------------------------------------------------------

  !-----------------------------------------------------------------------------
  !> quicksort (ascending order)
  recursive subroutine PRC_sort_ascd(a, top, bottom)
    implicit none
    integer, intent(inout) :: a(:)
    integer, intent(in)    :: top, bottom
    integer :: i, j, cnt, trg
    !---------------------------------------------------------------------------
    cnt = a( (top+bottom) / 2 )
    i = top; j = bottom
    do
       do while ( a(i) > cnt ) !ascending evaluation
          i = i + 1
       enddo
       do while ( cnt > a(j) ) !ascending evaluation
          j = j - 1
       enddo
       if ( i >= j ) exit
       trg = a(i);  a(i) = a(j);  a(j) = trg
       i = i + 1
       j = j - 1
    enddo
    if ( top < i-1    ) call PRC_sort_ascd( a, top, i-1    )
    if ( j+1 < bottom ) call PRC_sort_ascd( a, j+1, bottom )
    return
  end subroutine PRC_sort_ascd
  !-----------------------------------------------------------------------------

  !-----------------------------------------------------------------------------
  !> Get MPI time
  !> @return time
  function PRC_MPItime() result(time)
    implicit none

    real(DP) :: time
    !---------------------------------------------------------------------------

    if ( PRC_mpi_alive ) then
       time = real(MPI_WTIME(), kind=DP)
    else
       call cpu_time(time)
    endif

  end function PRC_MPItime

  !-----------------------------------------------------------------------------
  !> Calc global statistics for timer
  subroutine PRC_MPItimestat( &
      avgvar, &
      maxvar, &
      minvar, &
      maxidx, &
      minidx, &
      var     )
    implicit none

    real(DP), intent(out) :: avgvar(:) !< average
    real(DP), intent(out) :: maxvar(:) !< maximum
    real(DP), intent(out) :: minvar(:) !< minimum
    integer,  intent(out) :: maxidx(:) !< index of maximum
    integer,  intent(out) :: minidx(:) !< index of minimum
    real(DP), intent(in)  :: var(:)    !< values for statistics

    real(DP), allocatable :: statval(:,:)
    integer               :: vsize

    real(DP) :: totalvar
    integer  :: ierr
    integer  :: v, p
    !---------------------------------------------------------------------------

    vsize = size(var(:))

    allocate( statval(vsize,0:PRC_nmax-1) )
    statval(:,:) = 0.0_DP

    do v = 1, vsize
       statval(v,PRC_myrank) = var(v)
    enddo

    ! MPI broadcast
    do p = 0, PRC_nmax-1
       call MPI_Bcast( statval(1,p),         &
                       vsize,                &
                       MPI_DOUBLE_PRECISION, &
                       p,                    &
                       LOCAL_COMM_WORLD,     &
                       ierr                  )
    enddo

    do v = 1, vsize
       totalvar = 0.0_DP
       do p = 0, PRC_nmax-1
          totalvar = totalvar + statval(v,p)
       enddo
       avgvar(v) = totalvar / PRC_nmax

       maxvar(v)   = maxval(statval(v,:))
       minvar(v)   = minval(statval(v,:))
       maxidx(v:v) = maxloc(statval(v,:))
       minidx(v:v) = minloc(statval(v,:))
    enddo

    deallocate( statval )

    return
  end subroutine PRC_MPItimestat

  !-----------------------------------------------------------------------------
  !> MPI Error Handler
  subroutine PRC_MPI_errorhandler( &
      comm,     &
      errcode   )
    implicit none

    ! attributes are needed to be the same with COMM_ERRHANDLER function
    integer :: comm    !< MPI communicator
    integer :: errcode !< error code

    character(len=MPI_MAX_ERROR_STRING) :: msg
    integer :: len
    integer :: ierr
    !---------------------------------------------------------------------------

    ! Print Error Messages
    if ( PRC_mpi_alive ) then
          ! flush 1kbyte
       if ( IO_L ) then
          write(IO_FID_LOG,'(32A32)') '                                '
          write(IO_FID_LOG,*) '++++++ Abort MPI'
          write(IO_FID_LOG,*) ''
       end if

       if ( errcode .eq. abort_code ) then ! called from PRC_MPIstop
       elseif ( errcode <= MPI_ERR_LASTCODE ) then
          call MPI_ERROR_STRING(errcode, msg, len, ierr)
          if ( IO_L ) write(IO_FID_LOG,*) '++++++ ', errcode, trim(msg)
          write(*,*) '++++++ ', errcode, trim(msg)
       else
          if ( IO_L ) write(IO_FID_LOG,*) '++++++ Unexpected error code', errcode
          write(*,*) '++++++ Unexpected error code', errcode
       endif

       if ( comm .ne. ABORT_COMM_WORLD ) then
          if ( IO_L ) write(IO_FID_LOG,*) '++++++ Unexpected communicator'
          write(*,*) '++++++ Unexpected communicator'
       endif
       if ( IO_L ) write(IO_FID_LOG,*) ''
       write(*,*) ''
    endif

    call FileCloseAll

    ! Close logfile, configfile
    if ( IO_L ) then
       if ( IO_FID_LOG /= IO_FID_STDOUT ) close(IO_FID_LOG)
    endif
    close(IO_FID_CONF)

    ! Abort MPI
    if ( PRC_mpi_alive ) then
       call sleep(5)
       call MPI_ABORT(ABORT_COMM_WORLD, abort_code, ierr)
    endif

    stop
  end subroutine PRC_MPI_errorhandler

end module scale_process
!-------------------------------------------------------------------------------
