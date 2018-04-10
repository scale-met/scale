!-------------------------------------------------------------------------------
!> module process / cartesC
!!
!! @par Description
!!          MPI process management module for Cartesian-C grid
!!
!! @author Team SCALE
!!
!<
#include "scalelib.h"
module scale_prc_cartesC
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use mpi
  use scale_precision
  use scale_io
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: PRC_CARTESC_setup

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

  logical, public              :: PRC_PERIODIC_X   = .true.  !< periodic condition or not (X)?
  logical, public              :: PRC_PERIODIC_Y   = .true.  !< periodic condition or not (Y)?

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
  subroutine PRC_CARTESC_setup
    use scale_prc, only: &
       PRC_abort,              &
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

    logical :: PRC_CART_REORDER = .false. !< flag for rank reordering over the cartesian map

    namelist / PARAM_PRC_CARTESC / &
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

    LOG_NEWLINE
    LOG_INFO("PRC_CARTESC_setup",*) 'Setup'

    if ( IO_L ) then
       LOG_NEWLINE
       LOG_PROGRESS(*) 'start MPI'
       LOG_NEWLINE
       LOG_INFO("PRC_CARTESC_setup",*)            'Process information '
       LOG_INFO_CONT('(1x,A,I12)') 'UNIVERSAL_COMM_WORLD        : ', PRC_UNIVERSAL_COMM_WORLD
       LOG_INFO_CONT('(1x,A,I12)') 'total process [UNIVERSAL]   : ', PRC_UNIVERSAL_nprocs
       LOG_INFO_CONT('(1x,A,I12)') 'my process ID [UNIVERSAL]   : ', PRC_UNIVERSAL_myrank
       LOG_INFO_CONT('(1x,A,L12)') 'master rank?  [UNIVERSAL]   : ', PRC_UNIVERSAL_IsMaster
       LOG_INFO_CONT('(1x,A,I12)') 'GLOBAL_COMM_WORLD           : ', PRC_GLOBAL_COMM_WORLD
       LOG_INFO_CONT('(1x,A,I12)') 'total process [GLOBAL]      : ', PRC_GLOBAL_nprocs
       LOG_INFO_CONT('(1x,A,I12)') 'my process ID [GLOBAL]      : ', PRC_GLOBAL_myrank
       LOG_INFO_CONT('(1x,A,L12)') 'master rank?  [GLOBAL]      : ', PRC_GLOBAL_IsMaster
       LOG_INFO_CONT('(1x,A,I12)') 'LOCAL_COMM_WORLD            : ', PRC_LOCAL_COMM_WORLD
       LOG_INFO_CONT('(1x,A,I12)') 'total process [LOCAL]       : ', PRC_nprocs
       LOG_INFO_CONT('(1x,A,I12)') 'my process ID [LOCAL]       : ', PRC_myrank
       LOG_INFO_CONT('(1x,A,L12)') 'master rank?  [LOCAL]       : ', PRC_IsMaster
       LOG_INFO_CONT('(1x,A,I12)') 'ABORT_COMM_WORLD            : ', PRC_ABORT_COMM_WORLD
       LOG_INFO_CONT('(1x,A,I12)') 'master rank ID [each world] : ', PRC_masterrank
    endif

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_PRC_CARTESC,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("PRC_CARTESC_setup",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("PRC_CARTESC_setup",*) 'Not appropriate names in namelist PARAM_PRC_CARTESC. Check!'
       call PRC_abort
    endif
    if( IO_NML ) write(IO_FID_NML,nml=PARAM_PRC_CARTESC)

    LOG_NEWLINE
    LOG_INFO("PRC_CARTESC_setup",*) 'Process allocation '
    LOG_INFO_CONT(*) 'No. of Node   :', PRC_NUM_X," x ",PRC_NUM_Y

    if ( PRC_NUM_X*PRC_NUM_Y /= PRC_nprocs ) then
       LOG_ERROR("PRC_CARTESC_setup",*) 'total number of node does not match that requested. Check!'
       call PRC_abort
    endif

    if ( mod(PRC_nprocs,PRC_NUM_X) /= 0 ) then
       LOG_ERROR("PRC_CARTESC_setup",*) 'number of requested node cannot devide to 2D. Check!'
       call PRC_abort
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
       if (      .NOT. PRC_HAS_N &
            .OR. .NOT. PRC_HAS_W ) then
          PRC_next(PRC_NW) = MPI_PROC_NULL
       else
          next_coords(1) = coords_N(1)
          next_coords(2) = coords_W(2)
          call MPI_CART_RANK(iptbl, next_coords, PRC_next(PRC_NW), ierr)
       endif
       ! next rank search NorthEast
       if (      .NOT. PRC_HAS_N &
            .OR. .NOT. PRC_HAS_E ) then
          PRC_next(PRC_NE) = MPI_PROC_NULL
       else
          next_coords(1) = coords_N(1)
          next_coords(2) = coords_E(2)
          call MPI_CART_RANK(iptbl, next_coords, PRC_next(PRC_NE), ierr)
       endif
       ! next rank search SouthWest
       if (      .NOT. PRC_HAS_S &
            .OR. .NOT. PRC_HAS_W ) then
          PRC_next(PRC_SW) = MPI_PROC_NULL
       else
          next_coords(1) = coords_S(1)
          next_coords(2) = coords_W(2)
          call MPI_CART_RANK(iptbl, next_coords, PRC_next(PRC_SW), ierr)
       endif
       ! next rank search SouthEast
       if (      .NOT. PRC_HAS_S &
            .OR. .NOT. PRC_HAS_E ) then
          PRC_next(PRC_SE) = MPI_PROC_NULL
       else
          next_coords(1) = coords_S(1)
          next_coords(2) = coords_E(2)
          call MPI_CART_RANK(iptbl, next_coords, PRC_next(PRC_SE), ierr)
       endif
       call MPI_COMM_FREE(iptbl,ierr)
    endif

    next(:) = max(PRC_next(:),-1) ! avoid if MPI_PROC_NULL < -1

    LOG_INFO("PRC_CARTESC_setup",*) 'Node topology :'
    LOG_INFO_CONT('(1x,A,I5,A,I5,A,I5,A,A,I5,A,I5,A,I5,A,A,I5,A,I5,A,I5,A)') &
       'NW(',next(PRC_NW),',',PRC_2Drank(next(PRC_NW),1),',',PRC_2Drank(next(PRC_NW),2),')', &
    ' -  N(',next(PRC_N) ,',',PRC_2Drank(next(PRC_N) ,1),',',PRC_2Drank(next(PRC_N) ,2),')', &
    ' - NE(',next(PRC_NE),',',PRC_2Drank(next(PRC_NE),1),',',PRC_2Drank(next(PRC_NE),2),')'
    LOG_INFO_CONT('(1x,A)') '             |                       |                       |'
    LOG_INFO_CONT('(1x,A,I5,A,I5,A,I5,A,A,I5,A,I5,A,I5,A,A,I5,A,I5,A,I5,A)') &
       ' W(',next(PRC_W),',',PRC_2Drank(next(PRC_W),1),',',PRC_2Drank(next(PRC_W),2),')', &
    ' -  P(',PRC_myrank ,',',PRC_2Drank(PRC_myrank, 1),',',PRC_2Drank(PRC_myrank, 2),')', &
    ' -  E(',next(PRC_E),',',PRC_2Drank(next(PRC_E),1),',',PRC_2Drank(next(PRC_E),2),')'
    LOG_INFO_CONT('(1x,A)') '             |                       |                       |'
    LOG_INFO_CONT('(1x,A,I5,A,I5,A,I5,A,A,I5,A,I5,A,I5,A,A,I5,A,I5,A,I5,A)') &
       'SW(',next(PRC_SW),',',PRC_2Drank(next(PRC_SW),1),',',PRC_2Drank(next(PRC_SW),2),')', &
    ' -  S(',next(PRC_S) ,',',PRC_2Drank(next(PRC_S) ,1),',',PRC_2Drank(next(PRC_S) ,2),')', &
    ' - SE(',next(PRC_SE),',',PRC_2Drank(next(PRC_SE),1),',',PRC_2Drank(next(PRC_SE),2),')'

    return
  end subroutine PRC_CARTESC_setup

end module scale_prc_cartesC
