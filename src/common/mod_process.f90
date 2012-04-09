!-------------------------------------------------------------------------------
!> module PROCESS
!!
!! @par Description
!!          MPI/non-MPI management module
!!
!! @author H.Tomita and SCALE developpers
!!
!! @par History
!! @li      2011-10-11 (R.Yoshida) [new]
!! @li      2011-11-11 (H.Yashiro) [mod] Integrate to SCALE3
!!
!<
!-------------------------------------------------------------------------------
module mod_process
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use mpi
  use mod_stdio, only : &
     IO_FID_LOG, &
     IO_L
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: PRC_MPIstart
  public :: PRC_NOMPIstart
  public :: PRC_MPIstop
  public :: PRC_setup
  public :: PRC_MPItime
  public :: PRC_MPItimestat

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  integer, public,         parameter :: PRC_master = 0   !< master node

  integer, public,              save :: PRC_myrank = 0   !< my node ID
  integer, public,              save :: PRC_nmax   = 1   !< total number of processors
  integer, public,              save :: PRC_NUM_X  = 1   !< x length of 2D processor topology
  integer, public,              save :: PRC_NUM_Y  = 1   !< y length of 2D processor topology
  integer, public, allocatable, save :: PRC_2Drank(:,:)  !< node index in 2D topology

  integer, public,              save :: PRC_next(8) = -1 !< node ID of 8 neighbour process
  integer, public,         parameter :: PRC_W  = 1       !< [node direction] west
  integer, public,         parameter :: PRC_N  = 2       !< [node direction] north
  integer, public,         parameter :: PRC_E  = 3       !< [node direction] east
  integer, public,         parameter :: PRC_S  = 4       !< [node direction] south
  integer, public,         parameter :: PRC_NW = 5       !< [node direction] northwest
  integer, public,         parameter :: PRC_NE = 6       !< [node direction] northeast
  integer, public,         parameter :: PRC_SW = 7       !< [node direction] southwest
  integer, public,         parameter :: PRC_SE = 8       !< [node direction] southeast

  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  logical, private, save :: PRC_mpi_alive = .false. !< whether MPI is alive or not?

  !-----------------------------------------------------------------------------
contains

  !-----------------------------------------------------------------------------
  !> Start MPI
  !-----------------------------------------------------------------------------
  subroutine PRC_MPIstart
    use mod_stdio, only : &
       IO_FID_CONF,          &
       IO_FILECHR,           &
       IO_get_available_fid, &
       IO_make_idstr,        &
       IO_LOG_SUPPRESS,      &
       IO_LOG_ALLNODE
    implicit none

    character(len=IO_FILECHR) :: fname !< name of logfile for each process

    integer :: ierr
    !---------------------------------------------------------------------------

    call MPI_Init(ierr)
    call MPI_Comm_size(MPI_COMM_WORLD,PRC_nmax,  ierr)
    call MPI_Comm_rank(MPI_COMM_WORLD,PRC_myrank,ierr)

    PRC_mpi_alive = .true.

    if ( .not. IO_LOG_SUPPRESS ) then
       if ( PRC_myrank == PRC_master ) then ! master node
          IO_L = .true.
       else
          IO_L = IO_LOG_ALLNODE
       endif
    endif

    if ( IO_L ) then

       !--- Open logfile
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

       write(IO_FID_LOG,*)
       write(IO_FID_LOG,*) '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
       write(IO_FID_LOG,*) '+ SCALE: Scalable Computing by Advanced Library and Environment +'
       write(IO_FID_LOG,*) '+ SCALE-LES ver.3 (SCALE3): LES-scale Numerical weather model   +'
       write(IO_FID_LOG,*) '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
       write(IO_FID_LOG,*)
       write(IO_FID_LOG,*) '++++++ Start MPI'
       write(IO_FID_LOG,*) '*** total process : ', PRC_nmax
       write(IO_FID_LOG,*) '*** master rank   : ', PRC_master
       write(IO_FID_LOG,*) '*** my process ID : ', PRC_myrank
       write(IO_FID_LOG,*)
       write(IO_FID_LOG,*) '+++ Module[STDIO]/Categ[COMMON]'
       write(IO_FID_LOG,*) '*** Open config file, FID =', IO_FID_CONF
       write(IO_FID_LOG,*) '*** Open log    file, FID =', IO_FID_LOG
    else
       if ( PRC_myrank == PRC_master ) then ! master node
          write(*,*) '*** Log report is suppressed.'
       endif
    endif

    return
  end subroutine PRC_MPIstart

  !-----------------------------------------------------------------------------
  !> Dummy subroutine of MPIstart
  !-----------------------------------------------------------------------------
  subroutine PRC_NOMPIstart
    use mod_stdio, only : &
       IO_FID_CONF,          &
       IO_FILECHR,           &
       IO_get_available_fid, &
       IO_make_idstr,        &
       IO_LOG_SUPPRESS,      &
       IO_LOG_ALLNODE
    implicit none

    character(len=IO_FILECHR) :: fname

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

       write(IO_FID_LOG,*)
       write(IO_FID_LOG,*) '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
       write(IO_FID_LOG,*) '+ SCALE: Scalable Computing by Advanced Library and Environment +'
       write(IO_FID_LOG,*) '+ SCALE-LES ver.3 (SCALE3): LES-scale Numerical weather model   +'
       write(IO_FID_LOG,*) '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
       write(IO_FID_LOG,*)
       write(IO_FID_LOG,*) '++++++ Start WITHOUT MPI'
       write(IO_FID_LOG,*)
       write(IO_FID_LOG,*) '+++ Module[STDIO]/Categ[COMMON]'
       write(IO_FID_LOG,*) '*** Open config file, FID =', IO_FID_CONF
       write(IO_FID_LOG,*) '*** Open log    file, FID =', IO_FID_LOG

    endif

    return
  end subroutine PRC_NOMPIstart

  !-----------------------------------------------------------------------------
  !> Stop MPI
  !-----------------------------------------------------------------------------
  subroutine PRC_MPIstop
    use mod_stdio, only : &
       IO_FID_CONF, &
       IO_SYSCHR
    implicit none

    character(len=IO_SYSCHR) :: request = 'STOP'

    integer :: ierr
    !---------------------------------------------------------------------------

    ! Stop MPI
    if ( PRC_mpi_alive ) then
       if( IO_L ) write(IO_FID_LOG,*)
       if( IO_L ) write(IO_FID_LOG,*) '++++++ Stop MPI'
       if( IO_L ) write(IO_FID_LOG,*) '*** Broadcast STOP signal'
       call MPI_BCAST( request,        &
                       IO_SYSCHR,      &
                       MPI_CHARACTER,  & !--- type
                       PRC_master,     & !--- source rank
                       MPI_COMM_WORLD, &
                       ierr            )
       call MPI_Barrier(MPI_COMM_WORLD,ierr)
       call MPI_Finalize(ierr)
       if( IO_L ) write(IO_FID_LOG,*) '*** MPI is normaly finalized'
    endif

    ! Close logfile, configfile
    if ( IO_L ) then
       close(IO_FID_LOG)
    endif
    close(IO_FID_CONF)

    ! Stop program
    stop
  end subroutine PRC_MPIstop

  !-----------------------------------------------------------------------------
  !> Setup Processor topology
  !-----------------------------------------------------------------------------
  subroutine PRC_setup
    use mod_stdio, only: &
       IO_FID_CONF
    implicit none

    logical :: PRC_PERIODIC_X = .true. !< periodic condition or not (X)?
    logical :: PRC_PERIODIC_Y = .true. !< periodic condition or not (Y)? 

    namelist / PARAM_PRC / &
       PRC_NUM_X,      &
       PRC_NUM_Y,      &
       PRC_PERIODIC_X, &
       PRC_PERIODIC_Y

    logical :: period(2)
    integer :: divide(2)
    integer :: coords_W(2)
    integer :: coords_N(2)
    integer :: coords_E(2)
    integer :: coords_S(2)
    integer :: next_coords(2)
    integer :: iptbl

    integer :: ierr
    integer :: p
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[PRC]/Categ[COMMON]'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_PRC,iostat=ierr)

    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_PRC. Check!'
       call PRC_MPIstop
    endif
    if( IO_L ) write(IO_FID_LOG,nml=PARAM_PRC)

    if( IO_L ) write(IO_FID_LOG,*) '*** Process Allocation ***'
    if( IO_L ) write(IO_FID_LOG,*) '*** No. of Node :', PRC_NUM_X," x ",PRC_NUM_Y

    if ( PRC_NUM_X*PRC_NUM_Y /= PRC_nmax ) then
       write(*,*) 'xxx total number of node does not match that requested. Check!'
       call PRC_MPIstop
    endif

    if ( mod(PRC_nmax,PRC_NUM_X) /= 0 ) then
       write(*,*) 'xxx number of requested node cannot devide to 2D. Check!'
       call PRC_MPIstop
    endif

    ! set communication topology
    allocate( PRC_2Drank(-1:PRC_nmax-1,2) ); PRC_2Drank(:,:) = -1

    do p = 0, PRC_nmax-1
       PRC_2Drank(p,1) = mod(p,PRC_NUM_X)
       PRC_2Drank(p,2) = (p-PRC_2Drank(p,1)) / PRC_NUM_X
    enddo

    divide(1) = PRC_NUM_Y
    divide(2) = PRC_NUM_X
    period(1) = PRC_PERIODIC_Y
    period(2) = PRC_PERIODIC_X
    if ( PRC_mpi_alive ) then
       call MPI_CART_CREATE(MPI_COMM_WORLD,2,divide,period,.false.,iptbl,ierr)
       call MPI_CART_SHIFT(iptbl,0,1,PRC_next(PRC_S),PRC_next(PRC_N),ierr) ! next rank search Down/Up
       call MPI_CART_SHIFT(iptbl,1,1,PRC_next(PRC_W),PRC_next(PRC_E),ierr) ! next rank search Left/Right

       ! get neighbor_coordinates
       call MPI_CART_COORDS(iptbl,PRC_next(PRC_W),2,coords_W,ierr)
       call MPI_CART_COORDS(iptbl,PRC_next(PRC_N),2,coords_N,ierr)
       call MPI_CART_COORDS(iptbl,PRC_next(PRC_E),2,coords_E,ierr)
       call MPI_CART_COORDS(iptbl,PRC_next(PRC_S),2,coords_S,ierr)
       ! next rank search NorthWest
       next_coords(1) = coords_N(1)
       next_coords(2) = coords_W(2)
       call MPI_CART_RANK(iptbl, next_coords, PRC_next(PRC_NW), ierr)
       ! next rank search NorthEast
       next_coords(1) = coords_N(1)
       next_coords(2) = coords_E(2)
       call MPI_CART_RANK(iptbl, next_coords, PRC_next(PRC_NE), ierr)
       ! next rank search SouthWest
       next_coords(1) = coords_S(1)
       next_coords(2) = coords_W(2)
       call MPI_CART_RANK(iptbl, next_coords, PRC_next(PRC_SW), ierr)
       ! next rank search SouthEast
       next_coords(1) = coords_S(1)
       next_coords(2) = coords_E(2)
       call MPI_CART_RANK(iptbl, next_coords, PRC_next(PRC_SE), ierr)
    endif

    if( IO_L ) write(IO_FID_LOG,*) '*** Node topology ***'
    if( IO_L ) write(IO_FID_LOG,'(1x,A,I5,A,I5,A,I5,A,A,I5,A,I5,A,I5,A,A,I5,A,I5,A,I5,A)') &
    '***  NW(',PRC_next(PRC_NW),',',PRC_2Drank(PRC_next(PRC_NW),1),',',PRC_2Drank(PRC_next(PRC_NW),2),')', &
      ' -  N(',PRC_next(PRC_N) ,',',PRC_2Drank(PRC_next(PRC_N),1) ,',',PRC_2Drank(PRC_next(PRC_N),2) ,')', &
      ' - NE(',PRC_next(PRC_NE),',',PRC_2Drank(PRC_next(PRC_NE),1),',',PRC_2Drank(PRC_next(PRC_NE),2),')'
    if( IO_L ) write(IO_FID_LOG,'(1x,A)') '***                                  |'
    if( IO_L ) write(IO_FID_LOG,'(1x,A,I5,A,I5,A,I5,A,A,I5,A,I5,A,I5,A,A,I5,A,I5,A,I5,A)') &
    '***   W(',PRC_next(PRC_W),',',PRC_2Drank(PRC_next(PRC_W),1),',',PRC_2Drank(PRC_next(PRC_W),2),')', &
      ' -  P(',PRC_myrank     ,',',PRC_2Drank(PRC_myrank,     1),',',PRC_2Drank(PRC_myrank,     2),')', &
      ' -  E(',PRC_next(PRC_E),',',PRC_2Drank(PRC_next(PRC_E),1),',',PRC_2Drank(PRC_next(PRC_E),2),')'
    if( IO_L ) write(IO_FID_LOG,'(1x,A)') '***                                  |'
    if( IO_L ) write(IO_FID_LOG,'(1x,A,I5,A,I5,A,I5,A,A,I5,A,I5,A,I5,A,A,I5,A,I5,A,I5,A)') &
    '***  SW(',PRC_next(PRC_SW),',',PRC_2Drank(PRC_next(PRC_SW),1),',',PRC_2Drank(PRC_next(PRC_SW),2),')', &
      ' -  S(',PRC_next(PRC_S) ,',',PRC_2Drank(PRC_next(PRC_S),1) ,',',PRC_2Drank(PRC_next(PRC_S),2) ,')', &
      ' - SE(',PRC_next(PRC_SE),',',PRC_2Drank(PRC_next(PRC_SE),1),',',PRC_2Drank(PRC_next(PRC_SE),2),')'

    return
  end subroutine PRC_setup

  !-----------------------------------------------------------------------------
  !> get MPI time
  !> @return time
  !-----------------------------------------------------------------------------
  function PRC_MPItime() result(time)
    implicit none

    real(8) :: time
    !---------------------------------------------------------------------------

    if ( PRC_mpi_alive ) then
       time = MPI_WTIME()
    else
       call cpu_time(time)
    endif

  end function PRC_MPItime

  !-----------------------------------------------------------------------------
  subroutine PRC_MPItimestat( &
      avgvar, &
      maxvar, &
      minvar, &
      maxidx, &
      minidx, &
      var     )
    implicit none

    real(8), intent(out) :: avgvar(:)
    real(8), intent(out) :: maxvar(:)
    real(8), intent(out) :: minvar(:)
    integer, intent(out) :: maxidx(:)
    integer, intent(out) :: minidx(:)
    real(8), intent(in)  :: var(:)

    real(8), allocatable :: statval(:,:)
    integer              :: vsize

    real(8) :: totalvar
    integer :: ierr
    integer :: v, p
    !---------------------------------------------------------------------------

    vsize = size(var(:))

    allocate( statval(vsize,0:PRC_nmax-1) )
    statval(:,:) = 0.D0

    do v = 1, vsize
       statval(v,PRC_myrank) = var(v)
    enddo

    ! MPI broadcast
    do p = 0, PRC_nmax-1
       call MPI_Bcast( statval(1,p),         &
                       vsize,                &
                       MPI_DOUBLE_PRECISION, &
                       p,                    &
                       MPI_COMM_WORLD,       &
                       ierr                  )
    enddo

    do v = 1, vsize
       totalvar = 0.D0
       do p = 0, PRC_nmax-1
          totalvar = totalvar + statval(v,p)
       enddo
       avgvar(v) = totalvar / PRC_nmax

       maxvar(v)   = maxval(statval(v,:))
       minvar(v)   = minval(statval(v,:))
       maxidx(v:v) = maxloc(statval(v,:))
       minidx(v:v) = minloc(statval(v,:))
    enddo

    return
  end subroutine PRC_MPItimestat

end module mod_process
