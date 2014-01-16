!-------------------------------------------------------------------------------
!> module LAND driver
!!
!! @par Description
!!          Land module driver
!!
!! @author Team SCALE
!! @li      2013-08-31 (T.Yamaura)  [new]
!!
!<
!-------------------------------------------------------------------------------
module mod_land
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use mod_precision
  use mod_stdio
  use mod_prof
  use mod_grid_index
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: LAND_setup
  public :: LAND_step

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
  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine LAND_setup
    use mod_land_vars, only: &
       sw_phy => LAND_sw_phy,  &
       LAND_vars_setup,        &
       LAND_vars_restart_read
    use mod_land_phy_bucket, only: &
       LAND_PHY_setup
    implicit none
    !---------------------------------------------------------------------------

    call LAND_vars_setup

    call LAND_vars_restart_read

    if ( sw_phy ) call LAND_PHY_setup

    return
  end subroutine LAND_setup

  !-----------------------------------------------------------------------------
  !> Land step
  subroutine LAND_step
    use mod_land_vars, only: &
       sw_phy => LAND_sw_phy, &
       SFLX_GH,            &
       SFLX_PREC,          &
       SFLX_QVLnd,         &
       TG,                 &
       QvEfc,              &
       I_EMIT,             &
       I_ALB,              &
       I_TCS,              &
       I_DZg,              &
       I_Z0M,              &
       I_Z0H,              &
       I_Z0E,              &
       LAND_PROPERTY,      &
       LAND_vars_fillhalo, &
       LAND_vars_history
    use mod_land_phy_bucket, only: &
       LAND_PHY
    use mod_cpl_vars, only: &
       CPL_putLnd,      &
       CPL_getCPL2Lnd,  &
       sw_AtmLnd => CPL_sw_AtmLnd
    implicit none
    !---------------------------------------------------------------------------

    !########## from Coupler ##########
    if ( sw_AtmLnd ) then
       call CPL_getCPL2Lnd( SFLX_GH   (:,:), & ! [OUT]
                            SFLX_PREC (:,:), & ! [OUT]
                            SFLX_QVLnd(:,:)  ) ! [OUT]
    endif

    !########## Physics ##########
    call PROF_rapstart('LND Physics')
    if ( sw_phy ) then
       call LAND_PHY
    endif
    call PROF_rapend  ('LND Physics')

    !########## Fill HALO ##########
    call LAND_vars_fillhalo

    !########## History & Monitor ##########
    call PROF_rapstart('LND History')
       call LAND_vars_history
    call PROF_rapend  ('LND History')

    !########## to Coupler ##########
    if ( sw_AtmLnd ) then
       call CPL_putLnd( TG           (:,:),        & ! [IN]
                        QvEfc        (:,:),        & ! [IN]
                        LAND_PROPERTY(:,:,I_EMIT), & ! [IN]
                        LAND_PROPERTY(:,:,I_ALB),  & ! [IN]
                        LAND_PROPERTY(:,:,I_TCS),  & ! [IN]
                        LAND_PROPERTY(:,:,I_DZg),  & ! [IN]
                        LAND_PROPERTY(:,:,I_Z0M),  & ! [IN]
                        LAND_PROPERTY(:,:,I_Z0H),  & ! [IN]
                        LAND_PROPERTY(:,:,I_Z0E)   ) ! [IN]
    endif

    return
  end subroutine LAND_step

end module mod_land
