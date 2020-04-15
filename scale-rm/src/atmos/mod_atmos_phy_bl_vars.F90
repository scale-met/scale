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
  use scale_atmos_grid_cartesC_index
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

  real(RP), public, allocatable :: ATMOS_PHY_BL_RHOU_t(:,:,:)   ! tendency RHOU [kg/m2/s2]
  real(RP), public, allocatable :: ATMOS_PHY_BL_RHOV_t(:,:,:)   ! tendency RHOV [kg/m2/s2]
  real(RP), public, allocatable :: ATMOS_PHY_BL_RHOT_t(:,:,:)   ! tendency RHOT [K*kg/m3/s]

  real(RP), public, allocatable, target :: ATMOS_PHY_BL_RHOQ_t(:,:,:,:) ! tendency rho*QTRC [kg/kg/s]

  real(RP), public, allocatable :: ATMOS_PHY_BL_Zi    (:,:)     ! depth of the PBL

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

    allocate( ATMOS_PHY_BL_RHOU_t(KA,IA,JA)    )
    allocate( ATMOS_PHY_BL_RHOV_t(KA,IA,JA)    )
    allocate( ATMOS_PHY_BL_RHOT_t(KA,IA,JA)    )
    allocate( ATMOS_PHY_BL_RHOQ_t(KA,IA,JA,QA) )
    ATMOS_PHY_BL_RHOU_t(:,:,:)   = UNDEF
    ATMOS_PHY_BL_RHOV_t(:,:,:)   = UNDEF
    ATMOS_PHY_BL_RHOT_t(:,:,:)   = UNDEF
    ATMOS_PHY_BL_RHOQ_t(:,:,:,:) = UNDEF

    allocate( ATMOS_PHY_BL_Zi(IA,JA) )
    ATMOS_PHY_BL_Zi(:,:) = 100.0_RP ! tentative

    return
  end subroutine ATMOS_PHY_BL_vars_setup

end module mod_atmos_phy_bl_vars
