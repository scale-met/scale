!-------------------------------------------------------------------------------
!> module SCALE
!!
!! @par Description
!!          initialization
!!
!<
module scale
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_stdio
  use scale_process
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: SCALE_init
  public :: SCALE_finalize

  ! from scale_process
  public :: PRC_abort
  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  ! from scale_precision
  public :: RP
  public :: DP
  public :: SP

  ! from scale_io
  public :: H_SHORT
  public :: H_MID
  public :: H_LONG
  public :: IO_FID_CONF
  public :: IO_FID_LOG
  public :: IO_FID_NML
  public :: IO_L
  public :: IO_NML
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
  !> Initialization
  subroutine SCALE_init( &
       app_name )
    use dc_log, only: &
       LogInit
    use scale_process, only: &
       PRC_MPIstart, &
       PRC_SINGLECOM_setup
    use scale_const, only: &
       CONST_setup
    use scale_prof, only: &
       PROF_setup, &
       PROF_rapstart
    implicit none
    character(len=*), optional :: app_name !> application name

    character(len=H_SHORT) :: name

    ! for prc setup
    integer :: comm
    integer :: nprocs
    integer :: myrank
    logical :: ismaster

    if ( present(app_name) ) then
       name = app_name
    else
       name = "SCALE APPLICATION"
    end if


    ! start MPI
    call PRC_MPIstart( comm ) ! [OUT]

    ! setup MPI communicator
    call PRC_SINGLECOM_setup( comm,    & ! [IN]
                              nprocs,  & ! [OUT]
                              myrank,  & ! [OUT]
                              ismaster ) ! [OUT]

    ! setup scale_io
    call IO_setup( name, allow_noconf = .true. )

    ! setup log
    call IO_LOG_setup( myrank, ismaster )
    call LogInit( IO_FID_CONF,       &
                  IO_FID_LOG, IO_L,  &
                  IO_FID_NML, IO_NML )

    ! setup profiler
    call PROF_setup

    ! setup constants
    call CONST_setup

    call PROF_rapstart( 'Main', 0 )

    return
  end subroutine SCALE_init

  subroutine SCALE_finalize
    use scale_file, only: &
       FILE_Close_All
    use scale_process, only: &
       PRC_MPIfinish
    use scale_prof, only: &
       PROF_rapend, &
       PROF_rapreport

    call PROF_rapend( 'Main', 0 )

    call FILE_Close_All
    call PROF_rapreport

    ! stop mpi
    call PRC_MPIfinish

    return
  end subroutine SCALE_finalize
end module scale
