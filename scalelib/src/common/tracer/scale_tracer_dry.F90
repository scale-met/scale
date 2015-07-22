!-------------------------------------------------------------------------------
!> module TRACER / dry
!!
!! @par Description
!!          Tracer dry module
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2013-12-05 (S.Nishizawa)   [new] imported from inc_tracer_dry.f90
!!
!<
!-------------------------------------------------------------------------------
module scale_tracer_dry
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
  public :: TRACER_dry_setup

  include "inc_tracer_dry.h"

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  subroutine TRACER_dry_setup
    implicit none

    return
  end subroutine TRACER_dry_setup

end module scale_tracer_dry
