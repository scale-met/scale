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
  integer, public, parameter :: I_DENS     = 1
  integer, public, parameter :: I_MOMZ     = 2
  integer, public, parameter :: I_MOMX     = 3
  integer, public, parameter :: I_MOMY     = 4
  integer, public, parameter :: I_RHOT     = 5
  integer, public, parameter :: I_QTRC     = 6

  integer, public, parameter :: I_BND_DENS = 1 ! reference density     [kg/m3]
  integer, public, parameter :: I_BND_VELZ = 2 ! reference momentum (z) [m/s]
  integer, public, parameter :: I_BND_VELX = 3 ! reference momentum (x) [m/s]
  integer, public, parameter :: I_BND_VELY = 4 ! reference momentum (y) [m/s]
  integer, public, parameter :: I_BND_POTT = 5 ! reference mass-weighted potential temperature [K]
  integer, public, parameter :: I_BND_QV   = 6 ! reference water vapor [kg/kg]

  integer, public, parameter :: I_BND_SIZE = 6 ! size of array for reference quantities

end module scale_index
