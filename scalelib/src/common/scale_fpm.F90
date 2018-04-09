!-------------------------------------------------------------------------------
!> module FPM
!!
!! @par Description
!!          Failure Process Management module
!!
!! @author Team SCALE
!!
!<
module scale_fpm
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use mpi
  use scale_stdio
  !-----------------------------------------------------------------------------
  implicit none
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: FPM_Init
  public :: FPM_Polling

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  integer, public :: FPM_MAX_FAILURE  = 1
  integer, public :: FPM_POLLING_FREQ = 5
  logical, public :: FPM_alive

  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  integer, private, parameter :: FPM_MANAGER_MASTER = 0  !< rank number of master

  integer, private :: FPM_UNIVERSAL_COMM           !< universal communicator
  integer, private :: FPM_unv_myproc               !< my proc in universal world
  integer, private :: FPM_unv_nprocs               !< num proc in universal world
  integer, private :: FPM_GLOBAL_COMM              !< global communicator
  integer, private :: FPM_glb_myproc               !< my proc in global world
  integer, private :: FPM_glb_nprocs               !< num proc in global world
  integer, private :: FPM_LOCAL_COMM               !< local communicator
  integer, private :: FPM_LOCAL_MASTER             !< master proc in local communicator
  integer, private :: FPM_lcl_myproc               !< my proc in local world
  integer, private :: FPM_lcl_nprocs               !< num proc in local world

  integer, private :: FPM_MANAGER_COMM             !< manager communicator
  integer, private :: FPM_num_member               !< number of manager members

  logical, private :: FPM_master                   !< flag of master manager
  logical, private :: FPM_manager                  !< flag of manager
  logical, private, allocatable :: FPM_running(:)  !< flag of running status
  logical, private, allocatable :: FPM_lcl_running(:)  !< flag of running status

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------

  !-----------------------------------------------------------------------------
  !> Initialize FPM
  subroutine FPM_Init( &
      max_failure,     &
      polling_freq,    &
      universal_comm,  &
      global_comm,     &
      local_comm,      &
      num_member,      &
      global_root,     &
      use_fpm          )
    implicit none
    integer, intent(in) :: max_failure     !< threshold of failure procs
    integer, intent(in) :: polling_freq    !< polling frequency per time step (0: no polling)
    integer, intent(in) :: universal_comm  !< communicator
    integer, intent(in) :: global_comm     !< communicator
    integer, intent(in) :: local_comm      !< communicator
    integer, intent(in) :: num_member      !< number of members (NUM_BULKJOB)
    integer, intent(in) :: global_root(:)  !< root ranks of global comms
    logical, intent(in) :: use_fpm         !< fpm switch

    integer, allocatable :: manager_list(:)
    integer, allocatable :: exclude_list(:)
    integer :: num_exclude
    integer :: group_univ               !< group ID for universal world
    integer :: group_manager            !< group ID for manager world
    integer :: i, j, k
    integer :: ierr

    !---------------------------------------------------------------------------

    FPM_alive       = use_fpm
    FPM_master      = .false.
    FPM_manager     = .false.
    FPM_num_member  = num_member
    FPM_MAX_FAILURE = max_failure
    FPM_POLLING_FREQ = polling_freq

    if ( FPM_alive ) then
       FPM_UNIVERSAL_COMM = universal_comm
       call MPI_COMM_RANK( FPM_UNIVERSAL_COMM, FPM_unv_myproc, ierr )
       call MPI_COMM_SIZE( FPM_UNIVERSAL_COMM, FPM_unv_nprocs, ierr )
       FPM_GLOBAL_COMM    = global_comm
       call MPI_COMM_RANK( FPM_GLOBAL_COMM, FPM_glb_myproc, ierr )
       call MPI_COMM_SIZE( FPM_GLOBAL_COMM, FPM_glb_nprocs, ierr )
       FPM_LOCAL_COMM     = local_comm
       FPM_LOCAL_MASTER   = 0
       call MPI_COMM_RANK( FPM_LOCAL_COMM, FPM_lcl_myproc, ierr )
       call MPI_COMM_SIZE( FPM_LOCAL_COMM, FPM_lcl_nprocs, ierr )

       if ( FPM_unv_myproc == FPM_MANAGER_MASTER ) FPM_master = .true.
       if ( FPM_master ) write(*,*) ''
       if ( FPM_master ) write(*,*) '*** Failure Procs Manager: available'
       if ( FPM_master ) write(*,*) '*** Threshold of Failure Procs = ', FPM_MAX_FAILURE
       if ( FPM_master ) then
          if ( FPM_POLLING_FREQ > 0 ) then
             write(*,*) '*** FPM Polling Frequency per DT = ', FPM_POLLING_FREQ
          else
             write(*,*) '*** FPM: NO Polling'
          endif
       endif

       ! create manager communicator
       allocate( manager_list(FPM_num_member) )
       do i=1, FPM_num_member
          manager_list(i) = global_root(i)
          if ( FPM_unv_myproc == manager_list(i) ) FPM_manager = .true.
       enddo

       num_exclude = FPM_unv_nprocs - FPM_num_member
       allocate( exclude_list(num_exclude) )
       j = 1
       k = 1
       do i=0, FPM_unv_nprocs-1
          if ( i == manager_list(j) ) then
             if ( j < FPM_num_member ) j = j + 1
          else
             exclude_list(k) = i
             if ( k < num_exclude ) k = k + 1
          endif
       enddo

       call MPI_COMM_GROUP( FPM_UNIVERSAL_COMM, &
                            group_univ,         &
                            ierr                )
       call MPI_GROUP_EXCL( group_univ,         &
                            num_exclude,        &
                            exclude_list,       &
                            group_manager,      &
                            ierr                )
       call MPI_COMM_CREATE( FPM_UNIVERSAL_COMM, &
                             group_manager,      &
                             FPM_MANAGER_COMM,   &
                             ierr                )

       allocate( FPM_running    (FPM_num_member) )
       allocate( FPM_lcl_running(FPM_lcl_nprocs) )
       FPM_running    (:) = .true.
       FPM_lcl_running(:) = .true.
    endif

  end subroutine FPM_Init

  !-----------------------------------------------------------------------------
  !> Main system of FPM
  subroutine FPM_Polling( &
      run_stat,    &
      stop_signal  )
    implicit none
    logical, intent(in ) :: run_stat     !< running status
    logical, intent(out) :: stop_signal  !< exit sign

    integer :: sendcounts, recvcounts
    integer :: failcount
    integer :: i
    integer :: ierr

    logical :: local_stat
    logical :: sendbuff
    logical :: force_stop
    !---------------------------------------------------------------------------

    sendcounts  = 1
    recvcounts  = 1
    stop_signal = .false.
    local_stat  = .true.

    ! participants level
    sendbuff = run_stat
    call MPI_GATHER( sendbuff,           &
                     sendcounts,         &
                     MPI_LOGICAL,        &
                     FPM_lcl_running(:), &
                     recvcounts,         &
                     MPI_LOGICAL,        &
                     FPM_LOCAL_MASTER,   &
                     FPM_LOCAL_COMM,     &
                     ierr                )

    ! manager level
    !-------------------------------------------<<<
    if ( FPM_manager ) then
       do i=1, FPM_lcl_nprocs
          if ( .NOT. FPM_lcl_running(i) ) then
             local_stat = .false.
             exit
          endif
       enddo

       !call MPI_BARRIER(FPM_MANAGER_COMM, ierr)
       sendbuff = local_stat
       call MPI_GATHER( sendbuff,           &
                        sendcounts,         &
                        MPI_LOGICAL,        &
                        FPM_running(:),     &
                        recvcounts,         &
                        MPI_LOGICAL,        &
                        FPM_MANAGER_MASTER, &
                        FPM_MANAGER_COMM,   &
                        ierr                )

       ! master level
       !=======================================<<<
       if ( FPM_master ) then
          failcount = 0
          do i=1, FPM_num_member
             if ( .NOT. FPM_running(i) ) then
                failcount = failcount + 1
             endif
          enddo

          if ( failcount >= FPM_MAX_FAILURE ) then
             stop_signal = .true.
          else
             stop_signal = .false.
          endif
       endif
       !========================================>>>
    endif
    !------------------------------------------->>>

    ! broadcast signal
    call MPI_BCAST( stop_signal,        &
                    sendcounts,         &
                    MPI_LOGICAL,        &
                    FPM_MANAGER_MASTER, &
                    FPM_UNIVERSAL_COMM, &
                    ierr                )

  end subroutine FPM_Polling

end module scale_fpm !failure_process_management
!-------------------------------------------------------------------------------
