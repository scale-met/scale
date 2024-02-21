!-------------------------------------------------------------------------------
!> module Communication for Ensemble system
!!
!! @par Description
!!          MPI Communication module for Ensemble system
!!
!! @author Team SCALE
!!
!<
!-------------------------------------------------------------------------------
#include "scalelib.h"
module scale_comm_ensemble
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use mpi
  use scale_precision
  use scale_io
  use scale_prof
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: COMM_ENSEMBLE_setup
  public :: COMM_ENSEMBLE_finalize

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  integer, public :: COMM_ENSEMBLE_world
  integer, public :: COMM_ENSEMBLE_myrank
  integer, public :: COMM_ENSEMBLE_nprocs

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
  !> Setup
  subroutine COMM_ENSEMBLE_setup
    use scale_prc, only: &
       PRC_UNIVERSAL_COMM_WORLD, &
       PRC_GLOBAL_myrank,        &
       PRC_GLOBAL_nprocs,        &
       PRC_LOCAL_COMM_WORLD,     &
       PRC_nprocs,               &
       PRC_myrank,               &
       PRC_abort
    implicit none

    integer, allocatable :: ensemble2color(:)
    integer, allocatable :: ensemble2key  (:)

    integer :: n
    integer :: myrank, nprocs
    integer :: ierr
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("COMM_ENSEMBLE_setup",*) 'Setup'

    !!--- read namelist
    !rewind(IO_FID_CONF)
    !read(IO_FID_CONF,nml=PARAM_COMM_ENSEMBLE,iostat=ierr)
    !if( ierr < 0 ) then !--- missing
    !   LOG_INFO("COMM_ENSEMBLE_setup",*) 'Not found namelist. Default used.'
    !elseif( ierr > 0 ) then !--- fatal error
    !   LOG_ERROR("COMM_ENSEMBLE_setup",*) 'Not appropriate names in namelist PARAM_COMM_ENSEMBLE. Check!'
    !   call PRC_abort
    !endif
    !LOG_NML(PARAM_COMM_ENSEMBLE)

    ! [global] -> [ensemble]
    if( PRC_GLOBAL_nprocs < 2 ) then
       LOG_ERROR("COMM_ENSEMBLE_setup",*) 'PRC_GLOBAL_nprocs must be greater than or equal to 2:', PRC_GLOBAL_nprocs
       call PRC_abort
    end if

    allocate( ensemble2color(0:PRC_GLOBAL_nprocs-1) )
    allocate( ensemble2key  (0:PRC_GLOBAL_nprocs-1) )

    do n = 0, PRC_GLOBAL_nprocs-1
       ensemble2color(n) = n
       ensemble2key  (n) = n
    end do

    call MPI_COMM_SPLIT( PRC_UNIVERSAL_COMM_WORLD,          & ! [IN]
                         ensemble2color(PRC_GLOBAL_myrank), & ! [IN]
                         ensemble2key  (PRC_GLOBAL_myrank), & ! [IN]
                         COMM_ENSEMBLE_world,               & ! [OUT]
                         ierr                               ) ! [OUT]

    call MPI_COMM_RANK( COMM_ENSEMBLE_world, COMM_ENSEMBLE_myrank, ierr )
    call MPI_COMM_SIZE( COMM_ENSEMBLE_world, COMM_ENSEMBLE_nprocs, ierr )

    deallocate( ensemble2color )
    deallocate( ensemble2key   )

    return
  end subroutine COMM_ENSEMBLE_setup

  !-----------------------------------------------------------------------------
  !> finalize
  subroutine COMM_ENSEMBLE_finalize
    implicit none

    integer :: ierr
    !---------------------------------------------------------------------------

    call MPI_COMM_FREE( COMM_ENSEMBLE_world, ierr )

    return
  end subroutine COMM_ENSEMBLE_finalize

end module scale_comm_ensemble
