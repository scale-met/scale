!-------------------------------------------------------------------------------
!> module PRECISION
!!
!! @par Description
!!          precision module
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2013-12-03 (S.Nishizawa) [new]
!!
!<
#include "scalelib.h"
module scale_precision
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  !-----------------------------------------------------------------------------
  integer, public, parameter :: SP = kind(0.0E0) ! Single Precision
  integer, public, parameter :: DP = kind(0.0D0) ! Double Precision

  integer, public, parameter :: SP_PREC = precision(0.E0)
  integer, public, parameter :: DP_PREC = precision(0.D0)

#ifdef SINGLE
  integer, public, parameter :: RP      = SP      ! single precision
  integer, public, parameter :: RP_PREC = SP_PREC
#else
  integer, public, parameter :: RP      = DP      ! double precision
  integer, public, parameter :: RP_PREC = DP_PREC
#endif

  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  !-----------------------------------------------------------------------------
end module scale_precision
