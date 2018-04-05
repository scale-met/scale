!-------------------------------------------------------------------------------
!> module Coupler admin
!!
!! @par Description
!!          Coupler submodel administrator
!!
!! @author Team SCALE
!<
!-------------------------------------------------------------------------------
#include "scalelib.h"
module mod_cpl_admin
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_stdio
  use scale_prof
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: CPL_ADMIN_setup
  public :: CPL_ADMIN_getscheme

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  logical, public :: CPL_sw ! do coupler calculation?

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
  subroutine CPL_ADMIN_setup
    use scale_prc, only: &
       PRC_abort
    use mod_atmos_admin, only: &
       ATMOS_PHY_SF_TYPE, &
       ATMOS_sw_phy_sf
    use mod_ocean_admin, only: &
       OCEAN_do
    use mod_land_admin, only: &
       LAND_do
    use mod_urban_admin, only: &
       URBAN_do
    implicit none
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_PROGRESS(*) 'Module[ADMIN] / Categ[CPL] / Origin[SCALE-RM]'

    !-----< module component check >-----

    LOG_NEWLINE
    LOG_INFO("CPL_ADMIN_setup",*) 'Coupler components '

    ! Atmos-Ocean/Land/Urban Switch
    if ( OCEAN_do .OR. LAND_do .OR. URBAN_do ) then
       CPL_sw = .true.
    else
       CPL_sw = .false.
    endif

    ! Check Atmos_Surface setting
    if ( CPL_sw ) then
       LOG_INFO("CPL_ADMIN_setup",*) 'Coupler : ON'

       if ( ATMOS_PHY_SF_TYPE == 'COUPLE' ) then
          ! do nothing
       elseif( ATMOS_PHY_SF_TYPE == 'NONE' ) then
          LOG_INFO("CPL_ADMIN_setup",*) '-> Surface Flux Type is forced to change from NONE to COUPLE.'
          ! overwrite
          ATMOS_PHY_SF_TYPE = 'COUPLE'
          ATMOS_sw_phy_sf   = .true.
       else
          LOG_INFO("CPL_ADMIN_setup",*) 'Surface Flux : ', trim(ATMOS_PHY_SF_TYPE)
          LOG_INFO("CPL_ADMIN_setup",*) 'xxx Setting conflicts between coupler and surface flux! STOP.'
          LOG_ERROR("CPL_ADMIN_setup",*) 'Setting conflicts between coupler and surface flux! STOP.'
          call PRC_abort
       endif
    else
       LOG_INFO("CPL_ADMIN_setup",*) 'Coupler : OFF'
    endif

    return
  end subroutine CPL_ADMIN_setup

  !-----------------------------------------------------------------------------
  !> Get name of scheme for each component
  subroutine CPL_ADMIN_getscheme
    implicit none
    !---------------------------------------------------------------------------

    return
  end subroutine CPL_ADMIN_getscheme

end module mod_cpl_admin
