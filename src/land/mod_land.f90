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
  use mod_stdio, only: &
     IO_FID_LOG,  &
     IO_L
  use mod_time, only: &
     TIME_rapstart, &
     TIME_rapend
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: LAND_setup
  public :: LAND_step

  !
  !++ included parameters
  !
include "inc_precision.h"
include "inc_index.h"
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
       SFLX_GH,       &
       SFLX_PREC,     &
       SFLX_QVLnd,    &
       TG,            &
       QvEfc,         &
       EMIT,          &
       ALB, TCS, DZg, &
       Z00, Z0R, Z0S, &
       Zt0, ZtR, ZtS, &
       Ze0, ZeR, ZeS, &
       sw_phy => LAND_sw_phy, &
       LAND_vars_fillhalo,    &
       LAND_vars_history
    use mod_land_phy_bucket, only: &
       LAND_PHY
    use mod_cpl_vars, only: &
       sw_AtmLnd => CPL_sw_AtmLnd
    use mod_cpl_atmos_land, only: &
       CPL_AtmLnd_putLnd,      &
       CPL_AtmLnd_getDat2Lnd,  &
       CPL_AtmLnd_flushDat2Lnd
    implicit none
integer::i,j
    !---------------------------------------------------------------------------

    !########## Surface Flux ##########
    if ( sw_AtmLnd ) then
       call CPL_AtmLnd_getDat2Lnd( &
          SFLX_GH, SFLX_PREC, SFLX_QVLnd )
       call CPL_AtmLnd_flushDat2Lnd
    endif

    !########## Physics ##########
    call TIME_rapstart('LND Physics')
    if ( sw_phy ) then
       call LAND_PHY
    endif
    call TIME_rapend  ('LND Physics')

    !########## Fill HALO ##########
    call LAND_vars_fillhalo

    !########## History & Monitor ##########
    call TIME_rapstart('LND History')
       call LAND_vars_history
    call TIME_rapend  ('LND History')

    !########## for Coupler ##########
    if ( sw_AtmLnd ) then
       call CPL_AtmLnd_putLnd( &
          TG, QvEfc, EMIT, ALB, TCS, DZg,             &
          Z00, Z0R, Z0S, Zt0, ZtR, ZtS, Ze0, ZeR, ZeS )
    endif

    return
  end subroutine LAND_step

end module mod_land
