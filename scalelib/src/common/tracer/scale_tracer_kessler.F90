!-------------------------------------------------------------------------------
!> module TRACER / kessler
!!
!! @par Description
!!          Tracer kessler module
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2013-12-05 (S.Nishizawa)   [new] imported from inc_tracer_kessler.f90
!!
!<
!-------------------------------------------------------------------------------
module scale_tracer_kessler
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
  public :: TRACER_kessler_setup

  include "inc_tracer_kessler.h"

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  subroutine TRACER_kessler_setup
    implicit none

    return
  end subroutine TRACER_kessler_setup

end module scale_tracer_kessler
