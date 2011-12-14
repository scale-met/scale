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
  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !

  integer, public, save :: PRC_master = 0 !< master node

  integer, public,              save :: PRC_myrank = 0 !< my node ID
  integer, public,              save :: PRC_nmax   = 1 !< total number of processors
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

!    call MPI_Init(ierr)
!    call MPI_Comm_size(MPI_COMM_WORLD,PRC_nmax,  ierr)
!    call MPI_Comm_rank(MPI_COMM_WORLD,PRC_myrank,ierr)

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
    write(IO_FID_LOG,*) '++++++ Start WITHOUT MPI'
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
  !> Stop MPI
  !-----------------------------------------------------------------------------
  subroutine PRC_MPIstop
    use mod_stdio, only : &
       IO_SYSCHR,  &
       IO_FID_LOG, &
       IO_L
    implicit none

!    character(len=IO_SYSCHR) :: request = 'STOP'

!    integer :: ierr
    !---------------------------------------------------------------------------

    ! Stop MPI
    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '++++++ Stop Program'
!    if( IO_L ) write(IO_FID_LOG,*) '++++++ Stop MPI'
!    if( IO_L ) write(IO_FID_LOG,*) '*** Broadcast STOP signal'
!    call MPI_BCAST( request,        &
!                    IO_SYSCHR,      &
!                    MPI_CHARACTER,  & !--- type
!                    PRC_master,     & !--- source rank
!                    MPI_COMM_WORLD, &
!                    ierr            )
!    call MPI_Barrier(MPI_COMM_WORLD,ierr)
!    call MPI_Finalize(ierr)
!    if( IO_L ) write(IO_FID_LOG,*) '*** MPI is normaly finalized'

    ! Stop program
    close(IO_FID_LOG)
    stop

    return
  end subroutine PRC_MPIstop

  !-----------------------------------------------------------------------------
  !> get MPI time
  !-----------------------------------------------------------------------------
  function PRC_MPItime() result(time)
    implicit none

    real(8) :: time
    !---------------------------------------------------------------------------

    call cpu_time(time)

  end function PRC_MPItime

end module mod_process
