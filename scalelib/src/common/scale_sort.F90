!-------------------------------------------------------------------------------
!> module SORT
!!
!! @par Description
!!          Sort data
!!
!! @author Team SCALE
!!
!<
#include "scalelib.h"
module scale_sort
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_io
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: SORT_exec

  interface SORT_exec
     module procedure SORT_exec_without_idx
     module procedure SORT_exec_with_idxs
     module procedure SORT_exec_with_idx
  end interface SORT_exec

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
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
  !> bubble sort
!OCL SERIAL
  subroutine SORT_exec_with_idxs( &
      npoints,      &
      val,          &
      idx_i, idx_j, &
      reverse       )
    !$acc routine seq
    implicit none
    integer,  intent(in)    :: npoints                ! number of interpolation points
    real(RP), intent(inout) :: val  (npoints)         ! value to sort
    integer,  intent(inout) :: idx_i(npoints)         ! i-index
    integer,  intent(inout) :: idx_j(npoints)         ! j-index

    logical,  intent(in), optional :: reverse

    real(RP) :: sig
    integer  :: itmp
    integer  :: jtmp
    real(RP) :: vtmp

    integer  :: n1, n2
    !---------------------------------------------------------------------------

    sig = 1.0_RP
    if ( present(reverse) ) then
       if ( reverse ) sig = -1.0_RP
    end if

    do n1 = 1, npoints-1
    do n2 = n1+1, npoints
       if ( val(n1) * sig > val(n2) * sig ) then
          itmp      = idx_i(n1)
          jtmp      = idx_j(n1)
          vtmp      = val  (n1)

          idx_i(n1) = idx_i(n2)
          idx_j(n1) = idx_j(n2)
          val  (n1) = val  (n2)

          idx_i(n2) = itmp
          idx_j(n2) = jtmp
          val  (n2) = vtmp
       endif
    enddo
    enddo

    return
  end subroutine SORT_exec_with_idxs

!OCL SERIAL
  subroutine SORT_exec_with_idx( &
      npoints,    &
      val, index, &
      reverse     )
    !$acc routine seq
    implicit none
    integer,  intent(in)    :: npoints                ! number of interpolation points
    real(RP), intent(inout) :: val  (npoints)         ! value to sort
    integer,  intent(inout) :: index(npoints)         ! index

    logical,  intent(in), optional :: reverse

    real(RP) :: sig
    integer  :: itmp
    real(RP) :: vtmp

    integer  :: n1, n2
    !---------------------------------------------------------------------------

    sig = 1.0_RP
    if ( present(reverse) ) then
       if ( reverse ) sig = -1.0_RP
    end if

    do n1 = 1, npoints-1
    do n2 = n1+1, npoints
       if ( val(n1) * sig > val(n2) * sig ) then
          itmp    = index(n1)
          vtmp    = val  (n1)

          index(n1) = index(n2)
          val  (n1) = val  (n2)

          index(n2) = itmp
          val  (n2) = vtmp
       endif
    enddo
    enddo

    return
  end subroutine SORT_exec_with_idx

!OCL SERIAL
  subroutine SORT_exec_without_idx( &
      npoints, &
      val,     &
      reverse  )
    !$acc routine seq
    implicit none
    integer,  intent(in)    :: npoints                ! number of interpolation points
    real(RP), intent(inout) :: val  (npoints)         ! value to sort

    logical,  intent(in), optional :: reverse

    real(RP) :: sig
    real(RP) :: vtmp

    integer  :: n1, n2
    !---------------------------------------------------------------------------

    sig = 1.0_RP
    if ( present(reverse) ) then
       if ( reverse ) sig = -1.0_RP
    end if

    do n1 = 1, npoints-1
    do n2 = n1+1, npoints
       if ( val(n1) * sig > val(n2) * sig ) then
          vtmp      = val  (n1)

          val  (n1) = val  (n2)

          val  (n2) = vtmp
       endif
    enddo
    enddo

    return
  end subroutine SORT_exec_without_idx

end module scale_sort
