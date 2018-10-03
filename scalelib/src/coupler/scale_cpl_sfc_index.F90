!-------------------------------------------------------------------------------
!> module coupler / surface-atmospehre
!!
!! @par Description
!!         Common index among surface submodels and atmosphere
!!
!! @author Team SCALE
!<
!-------------------------------------------------------------------------------
#include "scalelib.h"
module scale_cpl_sfc_index
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
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
  integer,  public, parameter :: N_RAD_RGN   = 3 ! number of the radiation band region
  integer,  public, parameter :: I_R_IR      = 1 ! infrared         lamda > 4.0um
  integer,  public, parameter :: I_R_NIR     = 2 ! near-IR  4.0um > lamda > 0.7um
  integer,  public, parameter :: I_R_VIS     = 3 ! visible  0.7um > lamda

  real(RP), public, parameter :: RAD_RGN_boundary_VIS = 0.7_RP ! [um]
  real(RP), public, parameter :: RAD_RGN_boundary_IR  = 4.0_RP ! [um]

  integer,  public, parameter :: N_RAD_DIR   = 2 ! number of the radiation direction
  integer,  public, parameter :: I_R_direct  = 1
  integer,  public, parameter :: I_R_diffuse = 2

end module scale_cpl_sfc_index
