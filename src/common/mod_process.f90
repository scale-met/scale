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
  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  integer, public,              save :: PRC_master = 0  !< master node

  integer, public,              save :: PRC_myrank = 0  !< my node ID
  integer, public,              save :: PRC_nmax   = 1  !< total number of processors
  integer, public,              save :: PRC_NUM_X  = 1
  integer, public,              save :: PRC_NUM_Y  = 1
  integer, public, allocatable, save :: PRC_2Drank(:,:)

  integer, public,              save :: PRC_next(4) = -1
  integer, public,         parameter :: PRC_W = 1
  integer, public,         parameter :: PRC_N = 2
  integer, public,         parameter :: PRC_E = 3
  integer, public,         parameter :: PRC_S = 4

  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  logical, private, save :: PRC_mpi_alive = .false.
  !-----------------------------------------------------------------------------
contains

  !-----------------------------------------------------------------------------
  !> Start MPI
  !-----------------------------------------------------------------------------
  subroutine PRC_MPIstart
    use mod_stdio, only : &
       IO_get_available_fid, &
       IO_make_idstr,        &
       IO_FILECHR,           &
       IO_FID_CONF,          &
       IO_FID_LOG,           &
       IO_L
    implicit none

    character(len=IO_FILECHR) :: fname

    integer :: ierr
    !---------------------------------------------------------------------------

    call MPI_Init(ierr)
    call MPI_Comm_size(MPI_COMM_WORLD,PRC_nmax,  ierr)
    call MPI_Comm_rank(MPI_COMM_WORLD,PRC_myrank,ierr)

    PRC_mpi_alive = .true.

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
    write(IO_FID_LOG,*) '+ Numerical model for LES-scale weather                         +'
    write(IO_FID_LOG,*) '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
    write(IO_FID_LOG,*)
    write(IO_FID_LOG,*) '++++++ Start MPI'
    write(IO_FID_LOG,*) '*** total process : ', PRC_nmax
    write(IO_FID_LOG,*) '*** master rank   : ', PRC_master
    write(IO_FID_LOG,*) '*** my process ID : ', PRC_myrank
    write(IO_FID_LOG,*)
    write(IO_FID_LOG,*) '+++ Module[STDIO]/Categ[COMMON]'
    write(IO_FID_LOG,*) '+++ Open config file, FID =', IO_FID_CONF
    write(IO_FID_LOG,*) '+++ Open log    file, FID =', IO_FID_LOG

    if ( .not. IO_L ) then
       write(IO_FID_LOG,*) '+++ Following log message is suppressed.'
    endif

    return
  end subroutine PRC_MPIstart

  !-----------------------------------------------------------------------------
  !> Start MPI
  !-----------------------------------------------------------------------------
  subroutine PRC_NOMPIstart
    use mod_stdio, only : &
       IO_get_available_fid, &
       IO_make_idstr,        &
       IO_FILECHR,           &
       IO_FID_CONF,          &
       IO_FID_LOG,           &
       IO_L
    implicit none

    character(len=IO_FILECHR) :: fname

    integer :: ierr
    !---------------------------------------------------------------------------

    PRC_nmax   = 1
    PRC_myrank = 0
    PRC_mpi_alive = .false.

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
    write(IO_FID_LOG,*) '+ Numerical model for LES-scale weather                         +'
    write(IO_FID_LOG,*) '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
    write(IO_FID_LOG,*)
    write(IO_FID_LOG,*) '++++++ Start Without MPI'
    write(IO_FID_LOG,*) '*** total process : ', PRC_nmax
    write(IO_FID_LOG,*) '*** master rank   : ', PRC_master
    write(IO_FID_LOG,*) '*** my process ID : ', PRC_myrank
    write(IO_FID_LOG,*)
    write(IO_FID_LOG,*) '+++ Module[STDIO]/Categ[COMMON]'
    write(IO_FID_LOG,*) '+++ Open config file, FID =', IO_FID_CONF
    write(IO_FID_LOG,*) '+++ Open log    file, FID =', IO_FID_LOG

    if ( .not. IO_L ) then
       write(IO_FID_LOG,*) '+++ Following log message is suppressed.'
    endif

    return
  end subroutine PRC_NOMPIstart

  !-----------------------------------------------------------------------------
  !> Stop MPI
  !-----------------------------------------------------------------------------
  subroutine PRC_MPIstop
    use mod_stdio, only : &
       IO_SYSCHR,  &
       IO_FID_LOG, &
       IO_L
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

    ! Stop program
    close(IO_FID_LOG)
    stop

    return
  end subroutine PRC_MPIstop


  !-----------------------------------------------------------------------------
  !> Setup Precessor numbers
  !-----------------------------------------------------------------------------
  subroutine PRC_setup
    use mod_stdio, only: &
       IO_FID_CONF, &
       IO_FID_LOG,  &
       IO_L
    implicit none

    logical :: PRC_PERIODIC_X = .true.
    logical :: PRC_PERIODIC_Y = .true.

    NAMELIST / PARAM_PRC / &
       PRC_NUM_X,      &
       PRC_NUM_Y,      &
       PRC_PERIODIC_X, &
       PRC_PERIODIC_Y

    logical :: period(2)
    integer :: divide(2)
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
       call MPI_CART_SHIFT(iptbl,0,1,PRC_next(PRC_N),PRC_next(PRC_S),ierr)  ! next rank search Left/Right
       call MPI_CART_SHIFT(iptbl,1,1,PRC_next(PRC_W),PRC_next(PRC_E),ierr)  ! next rank search Up/Down
    endif

    if( IO_L ) write(IO_FID_LOG,*) '*** Node topology ***'
    if( IO_L ) write(IO_FID_LOG,'(1x,A,A,I5,A,I5,A,I5,A)') '***                         ', &
         'N(',PRC_next(PRC_N),',',PRC_2Drank(PRC_next(PRC_N),1) ,',',PRC_2Drank(PRC_next(PRC_N),2) ,')'
    if( IO_L ) write(IO_FID_LOG,'(1x,A)') '***                                 |'
    if( IO_L ) write(IO_FID_LOG,'(1x,A,I5,A,I5,A,I5,A,A,I5,A,I5,A,I5,A,A,I5,A,I5,A,I5,A)') &
    '***  W(',PRC_next(PRC_W),',',PRC_2Drank(PRC_next(PRC_W),1),',',PRC_2Drank(PRC_next(PRC_W),2),')', &
      ' - P(',PRC_myrank     ,',',PRC_2Drank(PRC_myrank,     1),',',PRC_2Drank(PRC_myrank,     2),')', &
      ' - E(',PRC_next(PRC_E),',',PRC_2Drank(PRC_next(PRC_E),1),',',PRC_2Drank(PRC_next(PRC_E),2),')'
    if( IO_L ) write(IO_FID_LOG,'(1x,A)') '***                                 |'
    if( IO_L ) write(IO_FID_LOG,'(1x,A,A,I5,A,I5,A,I5,A)') '***                         ', &
         'S(',PRC_next(PRC_S),',',PRC_2Drank(PRC_next(PRC_S),1) ,',',PRC_2Drank(PRC_next(PRC_S),2) ,')'

    return
  end subroutine PRC_setup

  !-----------------------------------------------------------------------------
  !> get MPI time
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

end module mod_process
