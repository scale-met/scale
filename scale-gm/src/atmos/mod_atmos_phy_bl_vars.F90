!-------------------------------------------------------------------------------
!> module atmosphere / physics / PBL
!!
!! @par Description
!!          Container for mod_atmos_phy_bl
!!
!! @author Team SCALE
!!
!<
!-------------------------------------------------------------------------------
#include "scalelib.h"
module mod_atmos_phy_bl_vars
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_io
  use scale_prof
  use scale_atmos_grid_icoA_index
  use scale_tracer
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: ATMOS_PHY_BL_vars_setup

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  integer,  public              :: QS, QE

  real(RP), public, allocatable :: ATMOS_PHY_BL_Zi(:,:,:)        ! depth of the PBL
  real(RP), public, allocatable :: ATMOS_PHY_BL_SFLX_BUOY(:,:,:) ! surface flux of buoyancy
  real(RP), public, allocatable :: ATMOS_PHY_BL_QL(:,:,:,:)      ! cloud water
  real(RP), public, allocatable :: ATMOS_PHY_BL_cldfrac(:,:,:,:) ! cloud fraction

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
  subroutine ATMOS_PHY_BL_vars_setup
    use scale_prc, only: &
       PRC_abort
    use scale_const, only: &
       UNDEF => CONST_UNDEF
    implicit none

    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("ATMOS_PHY_BL_vars_setup",*) 'Setup'

    allocate( ATMOS_PHY_BL_Zi       (IA,JA,ADM_lall) )
    allocate( ATMOS_PHY_BL_SFLX_BUOY(IA,JA,ADM_lall) )
    allocate( ATMOS_PHY_BL_QL       (KA,IA,JA,ADM_lall) )
    allocate( ATMOS_PHY_BL_cldfrac  (KA,IA,JA,ADM_lall) )
    ATMOS_PHY_BL_Zi       (:,:,:)   = UNDEF
    ATMOS_PHY_BL_SFLX_BUOY(:,:,:)   = UNDEF
    ATMOS_PHY_BL_QL       (:,:,:,:) = UNDEF
    ATMOS_PHY_BL_cldfrac  (:,:,:,:) = UNDEF

    return
  end subroutine ATMOS_PHY_BL_vars_setup

end module mod_atmos_phy_bl_vars
