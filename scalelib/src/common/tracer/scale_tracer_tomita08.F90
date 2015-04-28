!-------------------------------------------------------------------------------
!> module TRACER / tomita08
!!
!! @par Description
!!          Tracer tomita08 module
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2013-12-05 (S.Nishizawa)   [new] imported from inc_tracer_tomita08.f90
!!
!<
!-------------------------------------------------------------------------------
module scale_tracer_tomita08
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
  public :: TRACER_tomita08_setup

  include "inc_tracer_tomita08.h"

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  subroutine TRACER_tomita08_setup
    implicit none

    return
  end subroutine TRACER_tomita08_setup

end module scale_tracer_tomita08
