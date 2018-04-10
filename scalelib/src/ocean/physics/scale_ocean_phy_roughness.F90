!-------------------------------------------------------------------------------
!> module ocean / physics / surface roughness length
!!
!! @par Description
!!          surface roughness length common module
!!
!! @author Team SCALE
!!
!<
!-------------------------------------------------------------------------------
#include "scalelib.h"
module scale_ocean_phy_roughness
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
  public :: OCEAN_PHY_ROUGHNESS_setup

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  real(RP), public :: OCEAN_PHY_ROUGHNESS_visck     = 1.5E-5_RP ! kinematic viscosity
  real(RP), public :: OCEAN_PHY_ROUGHNESS_Ustar_min = 1.0E-3_RP ! minimum fiction velocity
  real(RP), public :: OCEAN_PHY_ROUGHNESS_Z0M_min   = 1.0E-5_RP ! minimum roughness length for momentum [m]
  real(RP), public :: OCEAN_PHY_ROUGHNESS_Z0H_min   = 1.0E-5_RP ! minimum roughness length for heat     [m]
  real(RP), public :: OCEAN_PHY_ROUGHNESS_Z0E_min   = 1.0E-5_RP ! minimum roughness length for moisture [m]

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  subroutine OCEAN_PHY_ROUGHNESS_setup
    use scale_prc, only: &
       PRC_abort
    implicit none

    NAMELIST / PARAM_OCEAN_PHY_ROUGHNESS / &
       OCEAN_PHY_ROUGHNESS_visck,     &
       OCEAN_PHY_ROUGHNESS_Ustar_min, &
       OCEAN_PHY_ROUGHNESS_Z0M_min,   &
       OCEAN_PHY_ROUGHNESS_Z0H_min,   &
       OCEAN_PHY_ROUGHNESS_Z0E_min

    integer :: ierr
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("OCEAN_PHY_ROUGHNESS_setup",*) 'Setup'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_OCEAN_PHY_ROUGHNESS,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("OCEAN_PHY_ROUGHNESS_setup",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("OCEAN_PHY_ROUGHNESS_setup",*) 'Not appropriate names in namelist PARAM_OCEAN_PHY_ROUGHNESS. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_OCEAN_PHY_ROUGHNESS)

    return
  end subroutine OCEAN_PHY_ROUGHNESS_setup

end module scale_ocean_phy_roughness
