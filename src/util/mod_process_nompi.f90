!-------------------------------------------------------------------------------
!> module PROCESS NOMPI
!!
!! @par Description
!!          Dummy module for non-MPI environment (no use, no include)
!!
!! @author H.Tomita and SCALE developpers
!!
!! @par History
!! @li      2011-11-11 (H.Yashiro) [new]
!!
!<
!-------------------------------------------------------------------------------
module mod_process
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
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
  public :: PRC_MPIstop
  public :: PRC_MPItime
  public :: PRC_MPItimestat

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  include "scale-les.h"

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
  !> Dummy subroutine of MPIstart
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
  end subroutine PRC_MPIstart

  !-----------------------------------------------------------------------------
  !> Stop MPI
  !-----------------------------------------------------------------------------
  subroutine PRC_MPIstop
    use mod_stdio, only : &
       IO_FID_CONF, &
       IO_SYSCHR
    implicit none

!    character(len=IO_SYSCHR) :: request = 'STOP'
!
!    integer :: ierr
    !---------------------------------------------------------------------------

    ! Stop MPI
!    if ( PRC_mpi_alive ) then
!       if( IO_L ) write(IO_FID_LOG,*)
!       if( IO_L ) write(IO_FID_LOG,*) '++++++ Stop MPI'
!       if( IO_L ) write(IO_FID_LOG,*) '*** Broadcast STOP signal'
!       call MPI_BCAST( request,        &
!                       IO_SYSCHR,      &
!                       MPI_CHARACTER,  & !--- type
!                       PRC_master,     & !--- source rank
!                       MPI_COMM_WORLD, &
!                       ierr            )
!       call MPI_Barrier(MPI_COMM_WORLD,ierr)
!       call MPI_Finalize(ierr)
!       if( IO_L ) write(IO_FID_LOG,*) '*** MPI is normaly finalized'
!    endif

    ! Close logfile, configfile
    if ( IO_L ) then
       close(IO_FID_LOG)
    endif
    close(IO_FID_CONF)

    ! Stop program
    stop
  end subroutine PRC_MPIstop

  !-----------------------------------------------------------------------------
  !> get MPI time
  !> @return time
  !-----------------------------------------------------------------------------
  function PRC_MPItime() result(time)
    implicit none

    real(8) :: time
    !---------------------------------------------------------------------------

!    if ( PRC_mpi_alive ) then
!       time = MPI_WTIME()
!    else
       call cpu_time(time)
!    endif

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
!    integer :: ierr
    integer :: v, p
    !---------------------------------------------------------------------------

    vsize = size(var(:))

    allocate( statval(vsize,0:PRC_nmax-1) )
    statval(:,:) = 0.D0

    do v = 1, vsize
       statval(v,PRC_myrank) = var(v)
    enddo

    ! MPI broadcast
!    do p = 0, PRC_nmax-1
!       call MPI_Bcast( statval(1,p),         &
!                       vsize,                &
!                       MPI_DOUBLE_PRECISION, &
!                       p,                    &
!                       MPI_COMM_WORLD,       &
!                       ierr                  )
!    enddo

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
