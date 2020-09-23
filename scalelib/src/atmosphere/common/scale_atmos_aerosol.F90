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
  use scale_io
  use scale_prof
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
  integer, public, parameter :: N_AE  = 7
  integer, public, parameter :: I_A01 = 1 ! Soil dust
  integer, public, parameter :: I_A02 = 2 ! Carbonacerous (BC/OC=0.3)
  integer, public, parameter :: I_A03 = 3 ! Carbonacerous (BC/OC=0.15)
  integer, public, parameter :: I_A04 = 4 ! Carbonacerous (BC/OC=0.)
  integer, public, parameter :: I_A05 = 5 ! Black carbon
  integer, public, parameter :: I_A06 = 6 ! Sulfate
  integer, public, parameter :: I_A07 = 7 ! Sea salt

  character(len=H_SHORT), public, parameter :: AE_NAME(N_AE) = &
       (/ "A01", "A02", "A03", "A04", "A05", "A06", "A07" /)
  character(len=H_MID),   public, parameter :: AE_DESC(N_AE) = &
       (/ "Soil dust                 ", &  ! Soil dust
          "Carbonacerous (BC/OC=0.3) ", &  ! Carbonacerous (BC/OC=0.3)
          "Carbonacerous (BC/OC=0.15)", &  ! Carbonacerous (BC/OC=0.15)
          "Carbonacerous (BC/OC=0.)  ", &  ! Carbonacerous (BC/OC=0.)
          "Black carbon              ", &  ! Black carbon
          "Sulfate                   ", &  ! Sulfate
          "Sea salt                  "  /) ! Sea salt

  real(RP), public, parameter :: AE_DENS(N_AE) = & ! aerosol density [kg/m3]
       (/ 2.50E+3_RP, &  ! Soil dust
          1.43E+3_RP, &  ! Carbonacerous (BC/OC=0.3)
          1.46E+3_RP, &  ! Carbonacerous (BC/OC=0.15)
          1.50E+3_RP, &  ! Carbonacerous (BC/OC=0.)
          1.25E+3_RP, &  ! Black carbon
          1.77E+3_RP, &  ! Sulfate
          2.20E+3_RP  /) ! Sea salt

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
