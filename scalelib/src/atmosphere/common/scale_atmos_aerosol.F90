!-------------------------------------------------------------------------------
!> module atmosphere / aerosol
!!
!! @par Description
!!          Aerosol module
!!
!! @author Team SCALE
!!
!<
!-------------------------------------------------------------------------------
#include "scalelib.h"
module scale_atmos_aerosol
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

  character(len=H_SHORT), public, parameter :: AE_NAME(N_AE) = &
       (/ "AD ", "ASO", "AVA", "AS ", "AR ", "ASS", "AU ", "AT ", "AOC" /)
  character(len=H_MID),   public, parameter :: AE_DESC(N_AE) = &
       (/ "dust-like      ", &
          "soot           ", &
          "volcanic ash   ", &
          "sulfate (H2SO4)", &
          "rural          ", &
          "sea salt       ", &
          "urban          ", &
          "tropo          ", &
          "yellow dust    " /)

  real(RP), parameter, private :: rhod_ae = 1.83_RP ! particle density [g/cm3] sulfate assumed
  real(RP), parameter, public :: AE_DENS(N_AE) = (/ rhod_ae, rhod_ae, rhod_ae, rhod_ae, rhod_ae, rhod_ae, rhod_ae, rhod_ae, rhod_ae /)


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

    LOG_NEWLINE
    LOG_INFO("ATMOS_AEROSOL_setup",*) 'Setup'

    return
  end subroutine ATMOS_AEROSOL_setup

end module scale_atmos_aerosol
