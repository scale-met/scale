!-------------------------------------------------------------------------------
!> module PROCESS
!!
!! @par Description
!!          MPI/non-MPI management module
!!
!! @author Team SCALE
!!
!<
#include "scalelib.h"
module scale_prc
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use mpi
  use scale_precision
  use scale_io
  use scale_fpm, only: &
     FPM_alive,  &
     FPM_Polling
  use scale_sigvars
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: PRC_MPIstart
  public :: PRC_UNIVERSAL_setup
  public :: PRC_GLOBAL_setup
  public :: PRC_LOCAL_setup
  public :: PRC_SINGLECOM_setup
  public :: PRC_abort
  public :: PRC_MPIfinish
  public :: PRC_MPIsplit_bulk
  public :: PRC_MPIsplit_nest

  public :: PRC_MPIbarrier
  public :: PRC_MPItime
  public :: PRC_MPItimestat

  public :: PRC_set_file_closer

  public :: PRC_ERRHANDLER_setup

  abstract interface
     subroutine closer( skip_abort )
       logical, intent(in), optional :: skip_abort
     end subroutine closer
  end interface

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  !-----------------------------------------------------------------------------
  !                          [ communicator system ]
  !    MPI_COMM_WORLD
  !          |
  ! PRC_UNIVERSAL_COMM_WORLD --split--> BULK_COMM_WORLD
  !                                     |
  !                            PRC_GLOBAL_COMM_WORLD --split--> PRC_LOCAL_COMM_WORLD
  !-----------------------------------------------------------------------------
  integer, public, parameter :: PRC_masterrank  = 0             !< master process in each communicator
  integer, public, parameter :: PRC_DOMAIN_nlim = 10000         !< max depth of domains
  integer, public, parameter :: PRC_COMM_NULL   = MPI_COMM_NULL

  ! universal world
  integer, public :: PRC_UNIVERSAL_COMM_WORLD = -1      !< original communicator
  integer, public :: PRC_UNIVERSAL_myrank     = -1      !< myrank         in universal communicator
  integer, public :: PRC_UNIVERSAL_nprocs     = -1      !< process num    in universal communicator
  logical, public :: PRC_UNIVERSAL_IsMaster   = .false. !< master process in universal communicator?

  integer, public :: PRC_UNIVERSAL_jobID      = 0       !< my job ID      in universal communicator

  ! global world
  integer, public :: PRC_GLOBAL_COMM_WORLD    = -1      !< global communicator
  integer, public :: PRC_GLOBAL_myrank        = -1      !< myrank         in global communicator
  integer, public :: PRC_GLOBAL_nprocs        = -1      !< process num    in global communicator
  logical, public :: PRC_GLOBAL_IsMaster      = .false. !< master process in global communicator?

  integer, public :: PRC_GLOBAL_domainID      = 0       !< my domain ID   in global communicator
  integer, public :: PRC_GLOBAL_ROOT(PRC_DOMAIN_nlim)   !< root processes in global members

  ! local world
  integer, public :: PRC_LOCAL_COMM_WORLD     = -1      !< local communicator
  integer, public :: PRC_nprocs               = 1       !< myrank         in local communicator
  integer, public :: PRC_myrank               = 0       !< process num    in local communicator
  logical, public :: PRC_IsMaster             = .false. !< master process in local communicator?
  !$acc declare create(PRC_myrank)

  ! inter-domain world
  integer, public :: PRC_INTERCOMM_parent     = MPI_COMM_NULL !< communicator between this rank and parent domain
  integer, public :: PRC_INTERCOMM_child      = MPI_COMM_NULL !< communicator between this rank and child  domain

  ! error handling
  logical, public :: PRC_mpi_alive = .false.            !< MPI is alive?
  integer, public :: PRC_UNIVERSAL_handler              !< error handler  in universal communicator
  integer, public :: PRC_ABORT_COMM_WORLD               !< communicator for aborting
  integer, public :: PRC_ABORT_handler                  !< error handler communicator for aborting

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
  integer, private, parameter :: PRC_ABORT_code = -1     !< MPI abort code in error handler
!  integer, private, parameter :: PRC_ABORT_code_p = 2 !< mpi abort code in error handler from parent
!  integer, private, parameter :: PRC_ABORT_code_d = 3 !< mpi abort code in error handler from daughter

  procedure(closer), pointer :: PRC_FILE_Closer => NULL()

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Start MPI
  subroutine PRC_MPIstart( &
       comm )
    implicit none

    integer, intent(out) :: comm ! communicator

    integer :: ierr
    !---------------------------------------------------------------------------

    call MPI_Init(ierr)

    PRC_mpi_alive = .true.
!    PRC_UNIVERSAL_handler = MPI_ERRHANDLER_NULL
!    call MPI_COMM_CREATE_ERRHANDLER( PRC_MPI_errorhandler, PRC_UNIVERSAL_handler, ierr )

    comm = MPI_COMM_WORLD
    PRC_ABORT_COMM_WORLD = comm

    return
  end subroutine PRC_MPIstart

  !-----------------------------------------------------------------------------
  !> setup MPI in universal communicator
  subroutine PRC_UNIVERSAL_setup( &
       comm,    &
       nprocs,  &
       myrank,  &
       ismaster )
    implicit none

    integer, intent(in)  :: comm     ! communicator
    integer, intent(out) :: nprocs   ! number of procs in this communicator
    integer, intent(out) :: myrank   ! myrank          in this communicator
    logical, intent(out) :: ismaster ! master process  in this communicator?

    integer :: ierr
    !---------------------------------------------------------------------------

    PRC_UNIVERSAL_COMM_WORLD = comm

    call MPI_Comm_size(PRC_UNIVERSAL_COMM_WORLD,PRC_UNIVERSAL_nprocs,ierr)
    call MPI_Comm_rank(PRC_UNIVERSAL_COMM_WORLD,PRC_UNIVERSAL_myrank,ierr)

    if ( PRC_UNIVERSAL_myrank == PRC_masterrank ) then
       PRC_UNIVERSAL_IsMaster = .true.
    else
       PRC_UNIVERSAL_IsMaster = .false.
    endif

    nprocs   = PRC_UNIVERSAL_nprocs
    myrank   = PRC_UNIVERSAL_myrank
    ismaster = PRC_UNIVERSAL_IsMaster



!    PRC_ABORT_COMM_WORLD = PRC_UNIVERSAL_COMM_WORLD
!
!    call MPI_Comm_set_errhandler(PRC_ABORT_COMM_WORLD,PRC_UNIVERSAL_handler,ierr)
!    call MPI_Comm_get_errhandler(PRC_ABORT_COMM_WORLD,PRC_ABORT_handler    ,ierr)

    return
  end subroutine PRC_UNIVERSAL_setup

  !-----------------------------------------------------------------------------
  !> setup MPI in global communicator
  subroutine PRC_GLOBAL_setup( &
       abortall, &
       comm      )
    implicit none

    logical, intent(in)  :: abortall ! abort all jobs?
    integer, intent(in)  :: comm     ! communicator

    integer :: ierr
    !---------------------------------------------------------------------------

    PRC_GLOBAL_COMM_WORLD = comm

    call MPI_Comm_size(PRC_GLOBAL_COMM_WORLD,PRC_GLOBAL_nprocs,ierr)
    call MPI_Comm_rank(PRC_GLOBAL_COMM_WORLD,PRC_GLOBAL_myrank,ierr)

    if ( PRC_GLOBAL_myrank == PRC_masterrank ) then
       PRC_GLOBAL_IsMaster = .true.
    else
       PRC_GLOBAL_IsMaster = .false.
    endif

!    if ( .NOT. abortall ) then
!       PRC_ABORT_COMM_WORLD = PRC_GLOBAL_COMM_WORLD
!
!       call MPI_COMM_SET_ERRHANDLER(PRC_ABORT_COMM_WORLD,PRC_UNIVERSAL_handler,ierr)
!       call MPI_COMM_GET_ERRHANDLER(PRC_ABORT_COMM_WORLD,PRC_ABORT_handler    ,ierr)
!    endif

    return
  end subroutine PRC_GLOBAL_setup

  !-----------------------------------------------------------------------------
  !> Setup MPI in local communicator
  subroutine PRC_LOCAL_setup( &
       comm,    &
       myrank,  &
       ismaster )
    implicit none

    integer, intent(in)  :: comm     ! communicator
    integer, intent(out) :: myrank   ! myrank         in this communicator
    logical, intent(out) :: ismaster ! master process in this communicator?

    integer :: ierr
    !---------------------------------------------------------------------------

    PRC_LOCAL_COMM_WORLD = comm

    call MPI_COMM_RANK(PRC_LOCAL_COMM_WORLD,PRC_myrank,ierr)
    call MPI_COMM_SIZE(PRC_LOCAL_COMM_WORLD,PRC_nprocs,ierr)
    !$acc update device(PRC_myrank)

    if ( PRC_myrank == PRC_masterrank ) then
       PRC_IsMaster = .true.
    else
       PRC_IsMaster = .false.
    endif

    myrank   = PRC_myrank
    ismaster = PRC_IsMaster

    return
  end subroutine PRC_LOCAL_setup

  !-----------------------------------------------------------------------------
  !> Setup MPI single communicator (not use universal-global-local setting)
  subroutine PRC_SINGLECOM_setup( &
       comm,    &
       nprocs,  &
       myrank,  &
       ismaster )
    implicit none

    integer, intent(in)  :: comm     ! communicator
    integer, intent(out) :: nprocs   ! number of procs
    integer, intent(out) :: myrank   ! myrank
    logical, intent(out) :: ismaster ! master process?

    integer :: ierr
    !---------------------------------------------------------------------------

    call MPI_Comm_size(comm,nprocs,ierr)
    call MPI_Comm_rank(comm,myrank,ierr)

    if ( myrank == PRC_masterrank ) then
       ismaster = .true.
    else
       ismaster = .false.
    endif

    PRC_UNIVERSAL_COMM_WORLD = comm
    PRC_UNIVERSAL_nprocs     = nprocs
    PRC_UNIVERSAL_myrank     = myrank
    PRC_UNIVERSAL_IsMaster   = ismaster

    PRC_GLOBAL_COMM_WORLD    = comm
    PRC_GLOBAL_nprocs        = nprocs
    PRC_GLOBAL_myrank        = myrank
    PRC_GLOBAL_IsMaster      = ismaster

    PRC_LOCAL_COMM_WORLD     = comm
    PRC_nprocs               = nprocs
    PRC_myrank               = myrank
    PRC_IsMaster             = ismaster
    !$acc update device(PRC_myrank)



    PRC_ABORT_COMM_WORLD = comm

!    call MPI_Comm_set_errhandler(PRC_ABORT_COMM_WORLD,PRC_UNIVERSAL_handler,ierr)
!    call MPI_Comm_get_errhandler(PRC_ABORT_COMM_WORLD,PRC_ABORT_handler    ,ierr)

    return
  end subroutine PRC_SINGLECOM_setup

  !-----------------------------------------------------------------------------
  !> Setup MPI error handler
  subroutine PRC_ERRHANDLER_setup( &
       use_fpm, &
       master   )
    implicit none

    logical, intent(in) :: use_fpm ! fpm switch
    logical, intent(in) :: master  ! master flag

    integer :: ierr
    !---------------------------------------------------------------------------

    call MPI_COMM_CREATE_ERRHANDLER(PRC_MPI_errorhandler,PRC_UNIVERSAL_handler,ierr)

    call MPI_COMM_SET_errhandler(PRC_ABORT_COMM_WORLD,PRC_UNIVERSAL_handler,ierr)
    call MPI_COMM_GET_errhandler(PRC_ABORT_COMM_WORLD,PRC_ABORT_handler    ,ierr)

    if ( PRC_UNIVERSAL_handler /= PRC_ABORT_handler ) then
       if( PRC_UNIVERSAL_IsMaster ) write(*,*) ""
       if( PRC_UNIVERSAL_IsMaster ) write(*,*) "ERROR: MPI HANDLER is INCONSISTENT"
       if( PRC_UNIVERSAL_IsMaster ) write(*,*) "     PRC_UNIVERSAL_handler = ", PRC_UNIVERSAL_handler
       if( PRC_UNIVERSAL_IsMaster ) write(*,*) "     PRC_ABORT_handler     = ", PRC_ABORT_handler
       call PRC_abort
    endif

    if ( use_fpm ) then
       call SIGVARS_Get_all( master )
       call signal( SIGINT,  PRC_abort )
!       call signal( SIGQUIT, PRC_abort )
!       call signal( SIGABRT, PRC_abort )
!       call signal( SIGFPE,  PRC_abort )
!       call signal( SIGSEGV, PRC_abort )
!       call signal( SIGTERM, PRC_abort )
    endif

    return
  end subroutine PRC_ERRHANDLER_setup

  !-----------------------------------------------------------------------------
  !> Abort Process
  subroutine PRC_abort
    implicit none

    integer :: ierr
    !---------------------------------------------------------------------------

    if ( PRC_mpi_alive ) then
       ! tentative approach; input "PRC_UNIVERSAL_COMM_WORLD".
       call MPI_COMM_CALL_ERRHANDLER(PRC_UNIVERSAL_COMM_WORLD,PRC_ABORT_code,ierr)
    endif

    stop
  end subroutine PRC_abort

  !-----------------------------------------------------------------------------
  !> Stop MPI peacefully
  subroutine PRC_MPIfinish
    implicit none

    integer :: ierr
    logical :: sign_status
    logical :: sign_exit
    !---------------------------------------------------------------------------

    call MPI_BARRIER(PRC_GLOBAL_COMM_WORLD, ierr)
    if ( PRC_INTERCOMM_child /= MPI_COMM_NULL ) &
         call MPI_Comm_free(PRC_INTERCOMM_child, ierr)
    if ( PRC_INTERCOMM_parent /= MPI_COMM_NULL ) &
         call MPI_Comm_free(PRC_INTERCOMM_PARENT, ierr)

    ! FPM polling
    if ( FPM_alive ) then
       sign_status = .false.
       sign_exit   = .false.
       do while ( .NOT. sign_exit )
          call FPM_Polling( sign_status, sign_exit )
       enddo
    endif

    if (PRC_UNIVERSAL_handler .NE. MPI_ERRHANDLER_NULL) then
        call MPI_Errhandler_free(PRC_UNIVERSAL_handler, ierr)
    endif
    if (PRC_ABORT_handler .NE. MPI_ERRHANDLER_NULL) then
        call MPI_Errhandler_free(PRC_ABORT_handler, ierr)
    endif

    ! Stop MPI
    if ( PRC_mpi_alive ) then
       ! free splitted communicator
       if ( PRC_LOCAL_COMM_WORLD  /= PRC_GLOBAL_COMM_WORLD ) then
          call MPI_Comm_free(PRC_LOCAL_COMM_WORLD,ierr)
       endif

       call MPI_Barrier(PRC_UNIVERSAL_COMM_WORLD,ierr)

       call MPI_Finalize(ierr)
    endif

    return
  end subroutine PRC_MPIfinish

  !-----------------------------------------------------------------------------
  !> MPI Communicator Split (bulk job)
  subroutine PRC_MPIsplit_bulk( &
      ORG_COMM_WORLD, &
      NUM_BULKJOB,    &
      PRC_BULKJOB,    &
      debug,          &
      SUB_COMM_WORLD, &
      ID_BULKJOB      )
    implicit none

    integer, intent(in)  :: ORG_COMM_WORLD ! communicator (original group)
    integer, intent(in)  :: NUM_BULKJOB    ! number of bulk jobs
    integer, intent(in)  :: PRC_BULKJOB(:) ! number of ranks in subgroup communicator
    logical, intent(in)  :: debug          ! log-output for mpi splitting?
    integer, intent(out) :: SUB_COMM_WORLD ! communicator (new subgroup)
    integer, intent(out) :: ID_BULKJOB     ! bulk job id

    integer :: ORG_myrank ! my rank         in the original communicator
    integer :: ORG_nrank  ! number of ranks in the original communicator
    integer :: sum_nrank

    integer, allocatable  :: prc2color (:)                 ! color id      for each process
    integer, allocatable  :: prc2key   (:)                 ! local rank id for each process
    integer               :: COL_domain(0:PRC_DOMAIN_nlim) ! domain id    of this color
    integer               :: COL_master(0:PRC_DOMAIN_nlim) ! master rank  of this color
    integer               :: COL_parent(0:PRC_DOMAIN_nlim) ! parent color of this color
    integer               :: COL_child (0:PRC_DOMAIN_nlim) ! child  color of this color

    logical :: color_reorder = .false. ! dummy

    integer :: i, color
    integer :: itag, ierr
    !---------------------------------------------------------------------------


    if ( NUM_BULKJOB == 1 ) then ! single domain run

       SUB_COMM_WORLD      = ORG_COMM_WORLD
       ID_BULKJOB          = 0
       PRC_UNIVERSAL_jobID = 0

    elseif( NUM_BULKJOB > 1 ) then ! multi domain run

       call MPI_COMM_RANK(ORG_COMM_WORLD,ORG_myrank,ierr)
       call MPI_COMM_SIZE(ORG_COMM_WORLD,ORG_nrank ,ierr)

       sum_nrank = sum(PRC_BULKJOB(1:NUM_BULKJOB))

       if ( sum_nrank /= ORG_nrank ) then
          if ( PRC_UNIVERSAL_IsMaster ) then
             LOG_ERROR("PRC_MPIsplit",*) "MPI PROCESS NUMBER is INCONSISTENT"
             LOG_ERROR_CONT(*)           " REQUESTED NPROCS = ", sum_nrank, "  LAUNCHED NPROCS = ", ORG_nrank
          endif
          call PRC_abort
       endif

       allocate( prc2color(0:ORG_nrank-1) )
       allocate( prc2key  (0:ORG_nrank-1) )

       call PRC_MPIcoloring( ORG_COMM_WORLD, & ! [IN]
                             ORG_nrank,      & ! [IN]
                             NUM_BULKJOB,    & ! [IN]
                             PRC_BULKJOB(:), & ! [IN]
                             color_reorder,  & ! [IN]
                             .true.,         & ! [IN]
                             debug,          & ! [IN]
                             prc2color (:),  & ! [OUT]
                             prc2key   (:),  & ! [OUT]
                             COL_domain(:),  & ! [OUT]
                             COL_master(:),  & ! [OUT]
                             COL_parent(:),  & ! [OUT]
                             COL_child (:)   ) ! [OUT]

       ! split comm_world
       call MPI_COMM_SPLIT( ORG_COMM_WORLD,        & ! [IN]
                            prc2color(ORG_myrank), & ! [IN]
                            prc2key  (ORG_myrank), & ! [IN]
                            SUB_COMM_WORLD,        & ! [OUT]
                            ierr                   ) ! [OUT]

       do i = 1, NUM_BULKJOB
          color = i-1
          PRC_GLOBAL_ROOT(i) = COL_master(color)
       enddo

       ID_BULKJOB          = prc2color(ORG_myrank) ! color = bulk id
       PRC_UNIVERSAL_jobID = prc2color(ORG_myrank) ! color = bulk id

       deallocate( prc2color )
       deallocate( prc2key   )

    else
       if ( PRC_UNIVERSAL_IsMaster ) then
          LOG_ERROR("PRC_MPIsplit",*) "REQUESTED DOMAIN NUMBER IS NOT ACCEPTABLE"
       endif
       call PRC_abort
    endif

    return
  end subroutine PRC_MPIsplit_bulk

  !-----------------------------------------------------------------------------
  !> MPI Communicator Split (nesting)
  subroutine PRC_MPIsplit_nest( &
      ORG_COMM_WORLD,   &
      NUM_DOMAIN,       &
      PRC_DOMAIN,       &
      debug,            &
      color_reorder,    &
      SUB_COMM_WORLD,   &
      ID_DOMAIN         )
    implicit none

    integer,               intent(in)  :: ORG_COMM_WORLD   ! communicator (original group)
    integer,               intent(in)  :: NUM_DOMAIN       ! number of bulk jobs
    integer,               intent(in)  :: PRC_DOMAIN(:)    ! number of ranks in subgroup communicator
    logical,               intent(in)  :: debug            ! log-output for mpi splitting?
    logical,               intent(in)  :: color_reorder    ! reorder
    integer,               intent(out) :: SUB_COMM_WORLD   ! communicator (new subgroup)
    integer,               intent(out) :: ID_DOMAIN        ! domain id

    integer :: ORG_myrank ! my rank         in the original communicator
    integer :: ORG_nrank  ! number of ranks in the original communicator
    integer :: sum_nrank

    integer, allocatable  :: prc2color (:)                 ! color id      for each process
    integer, allocatable  :: prc2key   (:)                 ! local rank id for each process
    integer               :: COL_domain(0:PRC_DOMAIN_nlim) ! domain id    of this color
    integer               :: COL_master(0:PRC_DOMAIN_nlim) ! master rank  of this color
    integer               :: COL_parent(0:PRC_DOMAIN_nlim) ! parent color of this color
    integer               :: COL_child (0:PRC_DOMAIN_nlim) ! child  color of this color

    integer :: i, color
    integer :: itag, ierr
    !---------------------------------------------------------------------------

    if ( NUM_DOMAIN == 1 ) then ! single domain run

       SUB_COMM_WORLD      = ORG_COMM_WORLD
       ID_DOMAIN           = 1
       PRC_GLOBAL_domainID = 1

    elseif( NUM_DOMAIN > 1 ) then ! multi domain run

       call MPI_COMM_RANK(ORG_COMM_WORLD,ORG_myrank,ierr)
       call MPI_COMM_SIZE(ORG_COMM_WORLD,ORG_nrank ,ierr)

       sum_nrank = sum(PRC_DOMAIN(1:NUM_DOMAIN))

       if ( sum_nrank /= ORG_nrank ) then
          if ( PRC_UNIVERSAL_IsMaster ) then
             LOG_ERROR("PRC_MPIsplit",*) "MPI PROCESS NUMBER is INCONSISTENT"
             LOG_ERROR_CONT(*)           " REQUESTED NPROCS = ", sum_nrank, "  LAUNCHED NPROCS = ", ORG_nrank
          endif
          call PRC_abort
       endif

       allocate( prc2color(0:ORG_nrank-1) )
       allocate( prc2key  (0:ORG_nrank-1) )

       call PRC_MPIcoloring( ORG_COMM_WORLD, & ! [IN]
                             ORG_nrank,      & ! [IN]
                             NUM_DOMAIN,     & ! [IN]
                             PRC_DOMAIN(:),  & ! [IN]
                             color_reorder,  & ! [IN]
                             .false.,        & ! [IN]
                             debug,          & ! [IN]
                             prc2color (:),  & ! [OUT]
                             prc2key   (:),  & ! [OUT]
                             COL_domain(:),  & ! [OUT]
                             COL_master(:),  & ! [OUT]
                             COL_parent(:),  & ! [OUT]
                             COL_child (:)   ) ! [OUT]

       ! split comm_world
       call MPI_COMM_SPLIT( ORG_COMM_WORLD,        & ! [IN]
                            prc2color(ORG_myrank), & ! [IN]
                            prc2key  (ORG_myrank), & ! [IN]
                            SUB_COMM_WORLD,        & ! [OUT]
                            ierr                   ) ! [OUT]

       ID_DOMAIN           = COL_domain(prc2color(ORG_myrank)) ! color /= domain id
       PRC_GLOBAL_domainID = COL_domain(prc2color(ORG_myrank)) ! color /= domain id

       if ( PRC_UNIVERSAL_IsMaster ) then
          write(*,*)
          write(*,*) "INFO [PRC_MPIsplit] Inter-domain relationship information"
          do i = 1, NUM_DOMAIN
             color = i-1
             if ( COL_parent(color) >= 0 ) then
                write(*,'(5x,A,I2.2)')          "Relationship No.", i
                write(*,'(5x,A,I2.2,A,A,I2.2)') "Parent color: ", COL_parent(color), " <=> ", &
                                                "Child  color: ", color
             endif
          enddo
       endif

       ! create inter-communicator
       do i = 1, NUM_DOMAIN
          color = i-1
          if ( COL_parent(color) >= 0 ) then
             itag  = i*100

             if    ( prc2color(ORG_myrank) == COL_parent(color) ) then ! as a parent

                call MPI_INTERCOMM_CREATE( SUB_COMM_WORLD, PRC_masterrank,    &
                                           ORG_COMM_WORLD, COL_master(color), &
                                           itag, PRC_INTERCOMM_child, ierr    )

             elseif( prc2color(ORG_myrank) == color ) then ! as a child

                call MPI_INTERCOMM_CREATE( SUB_COMM_WORLD, PRC_masterrank,                &
                                           ORG_COMM_WORLD, COL_master(COL_parent(color)), &
                                           itag, PRC_INTERCOMM_parent, ierr               )

             endif

             call MPI_BARRIER(ORG_COMM_WORLD,ierr)
          endif
       enddo

       deallocate( prc2color )
       deallocate( prc2key   )

    else
       if ( PRC_UNIVERSAL_IsMaster ) then
          write(*,*)"ERROR [RPC_MPIsplit] REQUESTED DOMAIN NUMBER IS NOT ACCEPTABLE"
       endif
       call PRC_abort
    endif

    return
  end subroutine PRC_MPIsplit_nest

  !-----------------------------------------------------------------------------
  !> Set color and keys for COMM_SPLIT
  subroutine PRC_MPIcoloring( &
      ORG_COMM_WORLD, &
      ORG_nrank,      &
      NUM_SUBGROUP,   &
      PRC_SUBGROUP,   &
      color_reorder,  &
      bulkjob,        &
      debug,          &
      prc2color,      &
      prc2key,        &
      COL_domain,     &
      COL_master,     &
      COL_parent,     &
      COL_child       )
    implicit none

    integer, intent(in)  :: ORG_COMM_WORLD
    integer, intent(in)  :: ORG_nrank
    integer, intent(in)  :: NUM_SUBGROUP
    integer, intent(in)  :: PRC_SUBGROUP(:)
    logical, intent(in)  :: color_reorder
    logical, intent(in)  :: bulkjob
    logical, intent(in)  :: debug
    integer, intent(out) :: prc2color (0:ORG_nrank-1)     ! color id      for each process
    integer, intent(out) :: prc2key   (0:ORG_nrank-1)     ! local rank id for each process
    integer, intent(out) :: COL_domain(0:PRC_DOMAIN_nlim) ! domain id    of the color
    integer, intent(out) :: COL_master(0:PRC_DOMAIN_nlim) ! master rank  of the color
    integer, intent(out) :: COL_parent(0:PRC_DOMAIN_nlim) ! parent color of the color
    integer, intent(out) :: COL_child (0:PRC_DOMAIN_nlim) ! child  color of the color

    integer :: PRC_REORDER(  NUM_SUBGROUP)   ! reordered  number of process
    integer :: DOM2COL    (  NUM_SUBGROUP)   ! get color  number by domain number
    integer :: COL2DOM    (0:NUM_SUBGROUP-1) ! get domain number by color  number
    logical :: touch      (0:NUM_SUBGROUP-1)

    integer :: id_parent  ! parent domain number
    integer :: id_child   ! child  domain number
    integer :: nprc, key
    integer :: i, domain, color, p
    integer :: ierr
    !---------------------------------------------------------------------------

    if ( color_reorder .AND. .NOT. bulkjob ) then ! with reordering of colors

       PRC_REORDER(1:NUM_SUBGROUP) = PRC_SUBGROUP(1:NUM_SUBGROUP)
       call PRC_sort_ascd( PRC_REORDER(1:NUM_SUBGROUP), 1, NUM_SUBGROUP )

       touch(:) = .false.
       do domain = 1, NUM_SUBGROUP
       do i      = NUM_SUBGROUP, 1, -1
          color = i-1 ! counted from 0
          if ( PRC_SUBGROUP(domain) == PRC_REORDER(i) .AND. ( .NOT. touch(color) ) ) then
             DOM2COL(domain) = color  ! domain_num -> color_num
             COL2DOM(color)  = domain ! color_num  -> domain_num
             touch  (color)  = .true.
             exit
          endif
       enddo ! order(=color+1) loop
       enddo ! domain loop

    else ! without reordering of colors

       do domain = 1, NUM_SUBGROUP
          color = domain-1 ! counted from 0
          DOM2COL(domain) = color  ! domain_num -> color_num
          COL2DOM(color)  = domain ! color_num  -> domain_num
       enddo ! domain loop

    endif

    ! make relationship
    COL_parent(:) = -1
    COL_child (:) = -1
    if ( .NOT. bulkjob ) then
       do i = 1, NUM_SUBGROUP
          id_parent = i - 1
          id_child  = i + 1

          if ( id_parent >= 1 .AND. id_parent <= NUM_SUBGROUP ) then
             COL_parent(DOM2COL(i)) = DOM2COL(id_parent)
          endif
          if ( id_child  >= 1 .AND. id_child  <= NUM_SUBGROUP ) then
             COL_child (DOM2COL(i)) = DOM2COL(id_child)
          endif

          if ( debug .AND. PRC_UNIVERSAL_IsMaster ) then
             write(*,'(4(A,I2))') &
             "DOMAIN: ", i, ", MY COL: ", DOM2COL(i), ", PARENT COL:", COL_parent(i), ", CHILD COL:", COL_child(i)
          endif
       enddo ! domain loop
    endif

    if ( PRC_UNIVERSAL_IsMaster ) then
       write(*,*)
       write(*,*)           'INFO [PRC_MPIcoloring] Domain information'
       write(*,'(5x,A,L2)')      'Reordering? : ', color_reorder
       do i = 1, NUM_SUBGROUP
          color  = i-1
          domain = COL2DOM(color)
          write(*,*)
          write(*,'(5x,2(A,I2.2))') "Order No. ", i, " -> Domain No. ", domain
          write(*,'(5x,A,I5)')      "Number of process      = ", PRC_SUBGROUP(domain)
          write(*,'(5x,A,I5)')      "Color of this   domain = ", color
          if ( COL_parent(color) >= 0 ) then
             write(*,'(5x,A,I5)')   "Color of parent domain = ", COL_parent(color)
          else
             write(*,'(5x,A)'   )   "Color of parent domain = no parent"
          endif
          if ( COL_child(color) >= 0 ) then
             write(*,'(5x,A,I5)')   "Color of child  domain = ", COL_child(color)
          else
             write(*,'(5x,A)'   )   "Color of child  domain = no child"
          endif
       enddo ! order(=color+1) loop
    endif

    ! make a process table
    color  = 0
    key    = 0
    domain = COL2DOM(color)
    nprc   = PRC_SUBGROUP(domain)
    COL_master(:) = -999

    do p = 0, ORG_nrank-1
       prc2color(p) = color
       prc2key  (p) = key

       if ( key == 0 ) then ! master rank
          COL_domain(color) = domain
          COL_master(color) = p
       endif

       if ( debug .AND. PRC_UNIVERSAL_IsMaster ) then
          write(*,'(5x,4(A,I5))') &
          "PE:", p, " COLOR:", prc2color(p), " KEY:", prc2key(p), " COL_master:", COL_master(color)
       endif

       key = key + 1
       if ( key >= nprc ) then
          color = color + 1
          key   = 0
          if ( color < NUM_SUBGROUP ) then ! ignore last
             domain = COL2DOM(color)
             nprc   = PRC_SUBGROUP(domain)
          endif
       endif
    enddo

    return
  end subroutine PRC_MPIcoloring

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
  !> Barrier MPI
  subroutine PRC_MPIbarrier
    implicit none

    integer  :: ierr
    !---------------------------------------------------------------------------

    if ( PRC_mpi_alive ) then
       call MPI_barrier(PRC_LOCAL_COMM_WORLD,ierr)
    endif

  end subroutine PRC_MPIbarrier

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
       call CPU_TIME(time)
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

    allocate( statval(vsize,0:PRC_nprocs-1) )
    statval(:,:) = 0.0_DP

    do v = 1, vsize
       statval(v,PRC_myrank) = var(v)
    enddo

    ! MPI broadcast
    do p = 0, PRC_nprocs-1
       call MPI_Bcast( statval(1,p),         &
                       vsize,                &
                       MPI_DOUBLE_PRECISION, &
                       p,                    &
                       PRC_LOCAL_COMM_WORLD, &
                       ierr                  )
    enddo

    do v = 1, vsize
       totalvar = 0.0_DP
       do p = 0, PRC_nprocs-1
          totalvar = totalvar + statval(v,p)
       enddo
       avgvar(v) = totalvar / PRC_nprocs

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
    logical :: sign_status
    logical :: sign_exit
    !---------------------------------------------------------------------------

!print *, "into errhandler:", PRC_UNIVERSAL_myrank

    ! FPM polling
    if ( FPM_alive ) then
       sign_status = .false.
       sign_exit   = .false.
       do while ( .NOT. sign_exit )
          call FPM_Polling( sign_status, sign_exit )
       enddo
    endif

    ! Print Error Messages
    if ( PRC_mpi_alive ) then
          ! flush 1kbyte
       if ( IO_L ) then
          LOG_PROGRESS(*) 'abort MPI'
          flush(IO_FID_LOG)
       endif

       if ( PRC_IsMaster ) then
          write(*,*) '+++++ BULK   ID       : ', PRC_UNIVERSAL_jobID
          write(*,*) '+++++ DOMAIN ID       : ', PRC_GLOBAL_domainID
          write(*,*) '+++++ MASTER LOCATION : ', PRC_UNIVERSAL_myrank,'/',PRC_UNIVERSAL_nprocs
          write(*,*) '+++++ GLOBAL LOCATION : ', PRC_GLOBAL_myrank,'/',PRC_GLOBAL_nprocs
          write(*,*) '+++++ LOCAL  LOCATION : ', PRC_myrank,'/',PRC_nprocs
          write(*,*) ''
       endif

       if    ( errcode == PRC_ABORT_code ) then ! called from PRC_abort
          ! do nothing
       elseif( errcode <= MPI_ERR_LASTCODE ) then
          call MPI_ERROR_STRING(errcode, msg, len, ierr)
          if( IO_L ) write(IO_FID_LOG,*) '+++++ ', errcode, trim(msg)
          write(*,*)                     '+++++ ', errcode, trim(msg)
       else
          if( IO_L ) write(IO_FID_LOG,*) '+++++ Unexpected error code', errcode
          write(*,*)                     '+++++ Unexpected error code', errcode
       endif

       if ( comm /= PRC_ABORT_COMM_WORLD ) then
          if( IO_L ) write(IO_FID_LOG,*) '+++++ Unexpected communicator'
          write(*,*)                     '+++++ Unexpected communicator'
       endif
       if( IO_L ) write(IO_FID_LOG,*) ''
       write(*,*)                     ''
    endif

    if ( associated( PRC_FILE_CLOSER ) ) call PRC_FILE_Closer( .true. )

    ! Close logfile, configfile
    if ( IO_L ) then
       if( IO_FID_LOG /= IO_FID_STDOUT ) close(IO_FID_LOG)
    endif
    close(IO_FID_CONF)

    ! Abort MPI
    if ( PRC_mpi_alive ) then
       call sleep(5)
       call MPI_ABORT(PRC_ABORT_COMM_WORLD, PRC_ABORT_code, ierr)
    endif

    stop
  end subroutine PRC_MPI_errorhandler

  subroutine PRC_set_file_closer( routine )
    procedure(closer) :: routine

    PRC_FILE_Closer => routine

    return
  end subroutine PRC_set_file_closer

end module scale_prc
