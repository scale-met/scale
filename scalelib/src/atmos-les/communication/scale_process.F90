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
  use scale_precision
  use scale_stdio
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ included parameters
  !
#include "scalelib.h"

  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: PRC_MPIstart
  public :: PRC_NOMPIstart
  public :: PRC_MPIstop
  public :: PRC_MPIfinish
  public :: PRC_setup
  public :: PRC_MPItime
  public :: PRC_MPItimestat

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  integer, public, parameter :: PRC_master = 0   !< master node

  integer, public            :: PRC_myrank = 0   !< my node ID
  integer, public            :: PRC_nmax   = 1   !< total number of processors
  integer, public            :: PRC_NUM_X  = 1   !< x length of 2D processor topology
  integer, public            :: PRC_NUM_Y  = 1   !< y length of 2D processor topology

  integer, public, allocatable, save :: PRC_2Drank(:,:)  !< node index in 2D topology

  integer, public            :: PRC_next(8) = -1 !< node ID of 8 neighbour process

  integer, public, parameter :: PRC_W  = 1       !< [node direction] west
  integer, public, parameter :: PRC_N  = 2       !< [node direction] north
  integer, public, parameter :: PRC_E  = 3       !< [node direction] east
  integer, public, parameter :: PRC_S  = 4       !< [node direction] south
  integer, public, parameter :: PRC_NW = 5       !< [node direction] northwest
  integer, public, parameter :: PRC_NE = 6       !< [node direction] northeast
  integer, public, parameter :: PRC_SW = 7       !< [node direction] southwest
  integer, public, parameter :: PRC_SE = 8       !< [node direction] southeast

  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  logical, private :: PRC_mpi_alive = .false. !< whether MPI is alive or not?

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Start MPI
  subroutine PRC_MPIstart
    implicit none

    character(len=H_LONG) :: fname ! name of logfile for each process

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
       call IO_make_idstr(fname,trim(IO_LOG_BASENAME),'pe',PRC_myrank)
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
       write(IO_FID_LOG,*) '+ SCALE-LES ver. '//VERSION//' : LES-scale Numerical weather model   +'
       write(IO_FID_LOG,*) '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
       write(IO_FID_LOG,*)
       write(IO_FID_LOG,*) '++++++ Start MPI'
       write(IO_FID_LOG,*) '*** total process : ', PRC_nmax
       write(IO_FID_LOG,*) '*** master rank   : ', PRC_master
       write(IO_FID_LOG,*) '*** my process ID : ', PRC_myrank
       write(IO_FID_LOG,*)
       write(IO_FID_LOG,*) '+++ Module[STDIO]/Categ[COMMON]'
       write(IO_FID_LOG,*) '*** Open config file, FID = ', IO_FID_CONF
       write(IO_FID_LOG,*) '*** Open log    file, FID = ', IO_FID_LOG
       write(IO_FID_LOG,*) '*** basename of log file  = ', trim(IO_LOG_BASENAME)

    else

       if ( PRC_myrank == PRC_master ) then ! master node
          write(*,*) '*** Log report is suppressed.'
       endif

    endif

    return
  end subroutine PRC_MPIstart

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
       write(IO_FID_LOG,*) '+ SCALE-LES ver. '//VERSION//' : LES-scale Numerical weather model   +'
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
  !> Abort MPI
  subroutine PRC_MPIstop
    implicit none

    integer :: ierr
    !---------------------------------------------------------------------------

    ! flush 1kbyte
    if( IO_L ) write(IO_FID_LOG,'(32A32)') '                                '

    if ( PRC_mpi_alive ) then
       if( IO_L ) write(IO_FID_LOG,*)
       if( IO_L ) write(IO_FID_LOG,*) '++++++ Abort MPI'
       if( IO_L ) close(IO_FID_LOG)
       call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
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
       if( IO_L ) write(IO_FID_LOG,*)
       if( IO_L ) write(IO_FID_LOG,*) '++++++ Stop MPI'
       if( IO_L ) write(IO_FID_LOG,*) '*** Broadcast STOP signal'
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
       close(IO_FID_LOG)
    endif
    close(IO_FID_CONF)

    ! Stop program
    stop
  end subroutine PRC_MPIfinish

  !-----------------------------------------------------------------------------
  !> Setup Processor topology
  subroutine PRC_setup
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
    integer :: next(8)

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
       call MPI_CART_CREATE(MPI_COMM_WORLD,2,divide,period,.false.,iptbl,ierr)
       call MPI_CART_SHIFT(iptbl,0,1,PRC_next(PRC_S),PRC_next(PRC_N),ierr) ! next rank search Down/Up
       call MPI_CART_SHIFT(iptbl,1,1,PRC_next(PRC_W),PRC_next(PRC_E),ierr) ! next rank search Left/Right

       ! get neighbor_coordinates
       if( PRC_next(PRC_W) /= MPI_PROC_NULL ) call MPI_CART_COORDS(iptbl,PRC_next(PRC_W),2,coords_W,ierr)
       if( PRC_next(PRC_N) /= MPI_PROC_NULL ) call MPI_CART_COORDS(iptbl,PRC_next(PRC_N),2,coords_N,ierr)
       if( PRC_next(PRC_E) /= MPI_PROC_NULL ) call MPI_CART_COORDS(iptbl,PRC_next(PRC_E),2,coords_E,ierr)
       if( PRC_next(PRC_S) /= MPI_PROC_NULL ) call MPI_CART_COORDS(iptbl,PRC_next(PRC_S),2,coords_S,ierr)
       ! next rank search NorthWest
       if (      PRC_next(PRC_N) == MPI_PROC_NULL &
            .OR. PRC_next(PRC_W) == MPI_PROC_NULL ) then
          PRC_next(PRC_NW) = MPI_PROC_NULL
       else
          next_coords(1) = coords_N(1)
          next_coords(2) = coords_W(2)
          call MPI_CART_RANK(iptbl, next_coords, PRC_next(PRC_NW), ierr)
       endif
       ! next rank search NorthEast
       if (      PRC_next(PRC_N) == MPI_PROC_NULL &
            .OR. PRC_next(PRC_E) == MPI_PROC_NULL ) then
          PRC_next(PRC_NE) = MPI_PROC_NULL
       else
          next_coords(1) = coords_N(1)
          next_coords(2) = coords_E(2)
          call MPI_CART_RANK(iptbl, next_coords, PRC_next(PRC_NE), ierr)
       endif
       ! next rank search SouthWest
       if (      PRC_next(PRC_S) == MPI_PROC_NULL &
            .OR. PRC_next(PRC_W) == MPI_PROC_NULL ) then
          PRC_next(PRC_SW) = MPI_PROC_NULL
       else
          next_coords(1) = coords_S(1)
          next_coords(2) = coords_W(2)
          call MPI_CART_RANK(iptbl, next_coords, PRC_next(PRC_SW), ierr)
       endif
       ! next rank search SouthEast
       if (      PRC_next(PRC_S) == MPI_PROC_NULL &
            .OR. PRC_next(PRC_E) == MPI_PROC_NULL ) then
          PRC_next(PRC_SE) = MPI_PROC_NULL
       else
          next_coords(1) = coords_S(1)
          next_coords(2) = coords_E(2)
          call MPI_CART_RANK(iptbl, next_coords, PRC_next(PRC_SE), ierr)
       endif
    endif

    next(:) = max(PRC_next(:),-1) ! avoid if MPI_PROC_NULL < -1

    if( IO_L ) write(IO_FID_LOG,*) '*** Node topology ***'
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
                       MPI_COMM_WORLD,       &
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

end module scale_process
!-------------------------------------------------------------------------------
