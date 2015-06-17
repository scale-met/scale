!-------------------------------------------------------------------------------------------
!> module NET2G error
!!
!! @par Description
!!          Error handling module for post-process of scale
!!
!! @author Team SCALE
!!
!! @par History
!! @li  2015-02-03 (R.Yoshida)  original
!!
!<
!-------------------------------------------------------------------------------------------
module mod_net2g_error
  !-----------------------------------------------------------------------------------------
  !
  !++ used modules
  !
  use mpi
  use netcdf

  use mod_net2g_vars

  !-----------------------------------------------------------------------------------------
  implicit none
  private
  !++ included parameters
#include "inc_net2g.h"
  !-----------------------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: handle_err
  public :: err_abort

  !-----------------------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  !-----------------------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  !-----------------------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  integer, private :: ierr

  !-----------------------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------------------

  !> error handler for netcdf90 system
  !---------------------------------------------------------------------------
  subroutine handle_err( &
      istat,  & ! [in]
      nline   ) ! [in]
    implicit none

    integer, intent(in) :: istat
    integer, intent(in) :: nline
    !---------------------------------------------------------------------------

    if ( LOUT ) write (*, *) nf90_strerror(istat)
    call err_abort( -1, nline )

    stop
  end subroutine handle_err


  !> process abort: error handing
  !---------------------------------------------------------------------------
  subroutine err_abort( &
      ecode,  & ! [in]
      nline   ) ! [in]
    implicit none

    integer, intent(in) :: ecode
    integer, intent(in) :: nline
    !---------------------------------------------------------------------------

    if ( ecode == err_internal ) then
       write (*, *) "##### ERROR: internal error"
    elseif ( ecode == err_netcdf ) then
       write (*, *) "##### ERROR: netcdf error"
    endif

!    write (*, *) "***** Abort: by rank =", irank
    write (*, *) "***** Abort: at Line =", nline

    if ( LOUT ) close ( FID_LOG )
    call MPI_ABORT( MPI_COMM_WORLD, ecode, ierr )

    stop
  end subroutine err_abort


end module mod_net2g_error
