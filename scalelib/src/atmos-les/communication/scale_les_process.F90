!-------------------------------------------------------------------------------
!> module LES PROCESS
!!
!! @par Description
!!          MPI process management module for LES model
!!
!! @author Team SCALE
!!
!<
module scale_les_process
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
  public :: PRC_setup

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  !-----------------------------------------------------------------------------
  integer, public, parameter   :: PRC_W  = 1       !< [node direction] west
  integer, public, parameter   :: PRC_N  = 2       !< [node direction] north
  integer, public, parameter   :: PRC_E  = 3       !< [node direction] east
  integer, public, parameter   :: PRC_S  = 4       !< [node direction] south
  integer, public, parameter   :: PRC_NW = 5       !< [node direction] northwest
  integer, public, parameter   :: PRC_NE = 6       !< [node direction] northeast
  integer, public, parameter   :: PRC_SW = 7       !< [node direction] southwest
  integer, public, parameter   :: PRC_SE = 8       !< [node direction] southeast

  integer, public              :: PRC_NUM_X  = 1   !< x length of 2D processor topology
  integer, public              :: PRC_NUM_Y  = 1   !< y length of 2D processor topology

  integer, public, allocatable :: PRC_2Drank(:,:)  !< node index in 2D topology
  integer, public              :: PRC_next(8) = -1 !< node ID of 8 neighbour process

  logical, public              :: PRC_HAS_W
  logical, public              :: PRC_HAS_N
  logical, public              :: PRC_HAS_E
  logical, public              :: PRC_HAS_S

  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup Processor topology
  subroutine PRC_setup
    use scale_process, only: &
       PRC_MPIstop,              &
       PRC_masterrank,           &
       PRC_mpi_alive,            &
       PRC_ABORT_COMM_WORLD,     &
       PRC_UNIVERSAL_COMM_WORLD, &
       PRC_UNIVERSAL_nprocs,     &
       PRC_UNIVERSAL_myrank,     &
       PRC_UNIVERSAL_IsMaster,   &
       PRC_GLOBAL_COMM_WORLD,    &
       PRC_GLOBAL_nprocs,        &
       PRC_GLOBAL_myrank,        &
       PRC_GLOBAL_IsMaster,      &
       PRC_LOCAL_COMM_WORLD,     &
       PRC_nprocs,               &
       PRC_myrank,               &
       PRC_IsMaster
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

    if ( IO_L ) then
       write(IO_FID_LOG,*)           ''
       write(IO_FID_LOG,*)           '++++++ Start MPI'
       write(IO_FID_LOG,'(1x,A,I8)') '*** UNIVERSAL_COMM_WORLD        : ', PRC_UNIVERSAL_COMM_WORLD
       write(IO_FID_LOG,'(1x,A,I8)') '*** total process [UNIVERSAL]   : ', PRC_UNIVERSAL_nprocs
       write(IO_FID_LOG,'(1x,A,I8)') '*** my process ID [UNIVERSAL]   : ', PRC_UNIVERSAL_myrank
       write(IO_FID_LOG,'(1x,A,L8)') '*** master rank?  [UNIVERSAL]   : ', PRC_UNIVERSAL_IsMaster
       write(IO_FID_LOG,'(1x,A,I8)') '*** GLOBAL_COMM_WORLD           : ', PRC_GLOBAL_COMM_WORLD
       write(IO_FID_LOG,'(1x,A,I8)') '*** total process [GLOBAL]      : ', PRC_GLOBAL_nprocs
       write(IO_FID_LOG,'(1x,A,I8)') '*** my process ID [GLOBAL]      : ', PRC_GLOBAL_myrank
       write(IO_FID_LOG,'(1x,A,L8)') '*** master rank?  [GLOBAL]      : ', PRC_GLOBAL_IsMaster
       write(IO_FID_LOG,'(1x,A,I8)') '*** LOCAL_COMM_WORLD            : ', PRC_LOCAL_COMM_WORLD
       write(IO_FID_LOG,'(1x,A,I8)') '*** total process [LOCAL]       : ', PRC_nprocs
       write(IO_FID_LOG,'(1x,A,I8)') '*** my process ID [LOCAL]       : ', PRC_myrank
       write(IO_FID_LOG,'(1x,A,L8)') '*** master rank?  [LOCAL]       : ', PRC_IsMaster
       write(IO_FID_LOG,'(1x,A,I8)') '*** ABORT_COMM_WORLD            : ', PRC_ABORT_COMM_WORLD
       write(IO_FID_LOG,'(1x,A,I8)') '*** master rank ID [each world] : ', PRC_masterrank
    endif

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

    if ( PRC_NUM_X*PRC_NUM_Y /= PRC_nprocs ) then
       write(*,*) 'xxx total number of node does not match that requested. Check!'
       call PRC_MPIstop
    endif

    if ( mod(PRC_nprocs,PRC_NUM_X) /= 0 ) then
       write(*,*) 'xxx number of requested node cannot devide to 2D. Check!'
       call PRC_MPIstop
    endif

    ! set communication topology
    allocate( PRC_2Drank(-1:PRC_nprocs-1,2) )
    PRC_2Drank(:,:) = -1

    do p = 0, PRC_nprocs-1
       PRC_2Drank(p,1) = mod(p,PRC_NUM_X)
       PRC_2Drank(p,2) = (p-PRC_2Drank(p,1)) / PRC_NUM_X
    enddo

    divide(1) = PRC_NUM_Y
    divide(2) = PRC_NUM_X
    period(1) = PRC_PERIODIC_Y
    period(2) = PRC_PERIODIC_X
    if ( PRC_mpi_alive ) then
       call MPI_CART_CREATE(PRC_LOCAL_COMM_WORLD,2,divide,period,PRC_CART_REORDER,iptbl,ierr)
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

end module scale_les_process
!-------------------------------------------------------------------------------
