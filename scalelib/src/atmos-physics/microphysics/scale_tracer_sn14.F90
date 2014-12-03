!-------------------------------------------------------------------------------
!> module TRACER / sn14
!!
!! @par Description
!!          Tracer sn14 module
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2013-12-05 (S.Nishizawa)   [new] imported from inc_tracer_sn14.f90
!!
!<
!-------------------------------------------------------------------------------
module scale_tracer_sn14
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
  public :: TRACER_sn14_setup

  include "inc_tracer_sn14.h"

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  subroutine TRACER_sn14_setup
    implicit none

    return
  end subroutine TRACER_sn14_setup

end module scale_tracer_sn14
