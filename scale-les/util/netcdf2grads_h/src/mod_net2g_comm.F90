!-------------------------------------------------------------------------------------------
!> module NET2G comm
!!
!! @par Description
!!          MPI coomunication module for post-process of scale
!!
!! @author Team SCALE
!!
!! @par History
!! @li  2015-02-03 (R.Yoshida)  original
!!
!<
!-------------------------------------------------------------------------------------------
module mod_net2g_comm
  !-----------------------------------------------------------------------------------------
  !
  !++ used modules
  !
  use mpi

  use mod_net2g_error

  !-----------------------------------------------------------------------------------------
  implicit none
  private
  !++ included parameters
#include "inc_net2g.h"
  !-----------------------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: comm_initialize
  public :: comm_setup
  public :: comm_finalize
  public :: comm_gather_grid
  public :: comm_gather_vars

  !-----------------------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  integer, public :: irank
  integer, public :: isize

  !-----------------------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  !-----------------------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  real(SP),allocatable :: sendbuf(:,:)
  real(SP),allocatable :: sendbuf_gx(:)
  real(SP),allocatable :: sendbuf_gy(:)
!  real(SP),allocatable :: recvbuf(:,:)
!  real(SP),allocatable :: cx_gather(:)
!  real(SP),allocatable :: cy_gather(:)
!  real(SP),allocatable :: cdx_gather(:)
!  real(SP),allocatable :: cdy_gather(:)

  !-----------------------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------------------


  !> communication: initialize
  !-----------------------------------------------------------------------------------------
  subroutine comm_initialize()
    implicit none

    integer :: ierr
    !---------------------------------------------------------------------------

    call MPI_INIT( ierr )
    call MPI_COMM_SIZE( MPI_COMM_WORLD, isize, ierr )
    call MPI_COMM_RANK( MPI_COMM_WORLD, irank, ierr )
    if ( irank == master ) write( *, * ) "+++ MPI COMM: Corrective Initialize"

    return
  end subroutine comm_initialize

  !> communication: setup
  !-----------------------------------------------------------------------------------------
  subroutine comm_setup( &
      mnxp,      & ! [in]
      mnyp,      & ! [in]
      nxgp,      & ! [in]
      nygp,      & ! [in]
      nmnge      ) ! [in]
!      tproc      ) ! [in]
    implicit none

    integer, intent(in) :: mnxp, mnyp
    integer, intent(in) :: nxgp, nygp
    integer, intent(in) :: nmnge
!    integer, intent(in) :: tproc
    !---------------------------------------------------------------------------

    allocate( sendbuf     (mnxp,       mnyp*nmnge) )
    allocate( sendbuf_gx  (nxgp*nmnge            ) )
    allocate( sendbuf_gy  (nygp*nmnge            ) )

!    if ( irank == master ) then
!       allocate( recvbuf   (mnxp,             mnyp*nmnge*tproc) )
!       allocate( cx_gather (nxgp*nmnge*tproc                  ) )
!       allocate( cy_gather (nygp*nmnge*tproc                  ) )
!       allocate( cdx_gather(nxgp*nmnge*tproc                  ) )
!       allocate( cdy_gather(nygp*nmnge*tproc                  ) )
!    else
!       allocate( recvbuf   (1,                1               ) )
!       allocate( cx_gather (1                                 ) )
!       allocate( cy_gather (1                                 ) )
!       allocate( cdx_gather(1                                 ) )
!       allocate( cdy_gather(1                                 ) )
!    endif

    return
  end subroutine comm_setup

  !> communication: finalize
  !-----------------------------------------------------------------------------------------
  subroutine comm_finalize()
    implicit none

    integer :: ierr
    !---------------------------------------------------------------------------

    if ( irank == master ) write( *, * ) "+++ MPI COMM: Corrective Finalize"
    call MPI_FINALIZE( ierr )

    return
  end subroutine comm_finalize

  !> communication: gather [1D]
  !-----------------------------------------------------------------------------------------
  subroutine comm_gather_grid( &
      nxgp,        & ! [in ]
      nygp,        & ! [in ]
      nmnge,       & ! [in ]
      p_cx,        & ! [in ]
      p_cdx,       & ! [in ]
      p_cy,        & ! [in ]
      p_cdy,       & ! [in ]
      cx_gather,   & ! [out]
      cdx_gather,  & ! [out]
      cy_gather,   & ! [out]
      cdy_gather   ) ! [out]
    implicit none

    integer,  intent(in)  :: nxgp, nygp
    integer,  intent(in)  :: nmnge
    real(SP), intent(in)  :: p_cx(:), p_cdx(:)
    real(SP), intent(in)  :: p_cy(:), p_cdy(:)
    real(SP), intent(out) :: cx_gather(:), cdx_gather(:)
    real(SP), intent(out) :: cy_gather(:), cdy_gather(:)

    integer :: sendcounts
    integer :: recvcounts
    integer :: iix, jjy
    integer :: ierr
    !---------------------------------------------------------------------------

    ! grids for x-direction
    sendcounts = nxgp*nmnge
    recvcounts = nxgp*nmnge

    sendbuf_gx(:) = UNDEF_SP
    do iix = 1, nxgp*nmnge
       sendbuf_gx(iix) = real( p_cx(iix) )
    enddo
    call MPI_GATHER( sendbuf_gx(:),  &
                     sendcounts,     &
                     MPI_REAL,       &
                     cx_gather(:),   &
                     recvcounts,     &
                     MPI_REAL,       &
                     master,         &
                     MPI_COMM_WORLD, &
                     ierr            )

    sendbuf_gx(:) = UNDEF_SP
    do iix = 1, nxgp*nmnge
       sendbuf_gx(iix) = real( p_cdx(iix) )
    enddo
    call MPI_GATHER( sendbuf_gx(:),  &
                     sendcounts,     &
                     MPI_REAL,       &
                     cdx_gather(:),  &
                     recvcounts,     &
                     MPI_REAL,       &
                     master,         &
                     MPI_COMM_WORLD, &
                     ierr            )

    ! grids for y-direction
    sendcounts = nygp*nmnge
    recvcounts = nygp*nmnge

    sendbuf_gy(:) = UNDEF_SP
    do jjy = 1, nygp*nmnge
       sendbuf_gy(jjy) = real( p_cy(jjy) )
    enddo
    call MPI_GATHER( sendbuf_gy(:),  &
                     sendcounts,     &
                     MPI_REAL,       &
                     cy_gather(:),   &
                     recvcounts,     &
                     MPI_REAL,       &
                     master,         &
                     MPI_COMM_WORLD, &
                     ierr            )

    sendbuf_gy(:) = UNDEF_SP
    do jjy = 1, nygp*nmnge
       sendbuf_gy(jjy) = real( p_cdy(jjy) )
    enddo
    call MPI_GATHER( sendbuf_gy(:),  &
                     sendcounts,     &
                     MPI_REAL,       &
                     cdy_gather(:),  &
                     recvcounts,     &
                     MPI_REAL,       &
                     master,         &
                     MPI_COMM_WORLD, &
                     ierr            )

    return
  end subroutine comm_gather_grid

  !> communication: gather [2D]
  !-----------------------------------------------------------------------------------------
  subroutine comm_gather_vars( &
      mnxp,      & ! [in ]
      mnyp,      & ! [in ]
      nmnge,     & ! [in ]
      invar,     & ! [in ]
      recvbuf    ) ! [out]
    implicit none

    integer,  intent(in)  :: mnxp, mnyp
    integer,  intent(in)  :: nmnge
    real(SP), intent(in)  :: invar (:,:)
    real(SP), intent(out) :: recvbuf(:,:)

    integer :: sendcounts
    integer :: recvcounts
    integer :: iix, jjy
    integer :: ierr
    !---------------------------------------------------------------------------

    sendcounts = mnxp * mnyp * nmnge 
    recvcounts = mnxp * mnyp * nmnge 

    sendbuf(:,:) = UNDEF_SP
    do jjy = 1, mnyp*nmnge 
    do iix = 1, mnxp
       sendbuf(iix,jjy) = invar(iix,jjy)
    enddo
    enddo

    call MPI_GATHER( sendbuf(:,:), &
                     sendcounts,      &
                     MPI_REAL,        &
                     recvbuf(:,:), &
                     recvcounts,      &
                     MPI_REAL,        &
                     master,          &
                     MPI_COMM_WORLD,  &
                     ierr             )

    return
  end subroutine comm_gather_vars

end module mod_net2g_comm
