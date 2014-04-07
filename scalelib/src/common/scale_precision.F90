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
module scale_precision
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use dc_types, only: &
     dc_SP => SP, &
     dc_DP => DP
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ included parameters
  !
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  !-----------------------------------------------------------------------------
  integer, public, parameter :: SP = dc_SP ! Single Precision: kind(0.E0)
  integer, public, parameter :: DP = dc_DP ! Double Precision: kind(0.D0)
#ifdef SINGLE
  integer, public, parameter :: RP = SP ! single precision
#else
  integer, public, parameter :: RP = DP ! double precision
#endif
  !
  !-----------------------------------------------------------------------------
  !++ Private procedure
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  !-----------------------------------------------------------------------------
end module scale_precision
!-------------------------------------------------------------------------------
