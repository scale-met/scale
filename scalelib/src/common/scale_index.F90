!-------------------------------------------------------------------------------
!> module Index
!!
!! @par Description
!!          Index module
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2014-07-09 (S.Nishizawa)   [new]
!!
!<
!-------------------------------------------------------------------------------
module scale_index
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  !-----------------------------------------------------------------------------
  use scale_stdio
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
  integer, public, parameter :: I_DENS = 1
  integer, public, parameter :: I_MOMZ = 2
  integer, public, parameter :: I_MOMX = 3
  integer, public, parameter :: I_MOMY = 4
  integer, public, parameter :: I_RHOT = 5
  integer, public, parameter :: I_QTRC = 6

end module scale_index
