!-------------------------------------------------------------------------------
!> module TRACER / none
!!
!! @par Description
!!          Tracer none module
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2015-03-27 (Y.Sato)   [new] imported from scale_tracer_suzuki10.f90
!!
!<
!-------------------------------------------------------------------------------
module scale_aetracer_none
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_stdio
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: AETRACER_none_setup

  include "inc_aetracer_none.h"

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  subroutine AETRACER_none_setup
    use scale_process, only: &
      PRC_MPIstop

    return
  end subroutine AETRACER_none_setup

end module scale_aetracer_none
