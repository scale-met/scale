!-------------------------------------------------------------------------------
!> module Aerosol
!!
!! @par Description
!!          Aerosol module
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2016-08-09 (S.Nishizawa)   [new]
!!
!<
!-------------------------------------------------------------------------------
module scale_atmos_aerosol
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_stdio
  use scale_grid_index
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: ATMOS_AEROSOL_setup

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  integer, public, parameter :: N_AE = 9
  integer, public, parameter :: I_AD  = 1 !< dust-like
  integer, public, parameter :: I_ASO = 2 !< soot
  integer, public, parameter :: I_AVA = 3 !< volcanic ash
  integer, public, parameter :: I_AS  = 4 !< sulfate (H2SO4)
  integer, public, parameter :: I_AR  = 5 !< rural
  integer, public, parameter :: I_ASS = 6 !< sea salt
  integer, public, parameter :: I_AU  = 7 !< urban
  integer, public, parameter :: I_AT  = 8 !< tropo
  integer, public, parameter :: I_AOC = 9 !< yellow dust

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
  !> Setup
  subroutine ATMOS_AEROSOL_setup
    implicit none
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '++++++ Module[AEROSOL] / Categ[ATMOS SHARE] / Origin[SCALElib]'

    return
  end subroutine ATMOS_AEROSOL_setup

end module scale_atmos_aerosol
