!-------------------------------------------------------------------------------
!> module Atmosphere
!!
!! @par Description
!!          Atmosphere module
!!
!! @author H.Tomita and SCALE developpers
!!
!! @par History
!! @li      2011-11-11 (H.Yashiro) [new]
!! @li      2012-02-15 (H.Yashiro) [add] Microphysics
!! @li      2012-03-23 (H.Yashiro) [mod] Cleaning
!!
!<
!-------------------------------------------------------------------------------
module mod_atmos
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
  public :: ATMOS_setup
  public :: ATMOS_step

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
  !> Setup atmosphere
  !-----------------------------------------------------------------------------
  subroutine ATMOS_setup
    use mod_atmos_vars, only: &
       ATMOS_vars_setup,  &
       ATMOS_vars_restart_read
    use mod_atmos_refstate, only: &
       ATMOS_REFSTATE_setup
    use mod_atmos_boundary, only: &
       ATMOS_BOUNDARY_setup
    use mod_atmos_thermodyn, only: &
       ATMOS_THRRMODYN_setup
    use mod_atmos_dyn, only: &
       ATMOS_DYN_setup
    use mod_atmos_phy_sf, only: &
       ATMOS_PHY_SF_setup
    use mod_atmos_phy_tb, only: &
       ATMOS_PHY_TB_setup
!    use mod_atmos_phy_mp, only: &
!       ATMOS_PHY_MP_setup
!    use mod_atmos_phy_rd, only: &
!       ATMOS_PHY_RD_setup
    implicit none
    !---------------------------------------------------------------------------

    call ATMOS_THRRMODYN_setup

    call ATMOS_vars_setup

    call ATMOS_vars_restart_read

    call ATMOS_REFSTATE_setup

    call ATMOS_BOUNDARY_setup

    call ATMOS_DYN_setup

    call ATMOS_PHY_SF_setup

    call ATMOS_PHY_TB_setup

!    call ATMOS_PHY_MP_setup
!
!    call ATMOS_PHY_RD_setup

    return
  end subroutine ATMOS_setup

  !-----------------------------------------------------------------------------
  !> advance atmospheric state
  !-----------------------------------------------------------------------------
  subroutine ATMOS_step
    use mod_time, only: &
       do_dyn    => TIME_DOATMOS_DYN,    &
       do_phy_tb => TIME_DOATMOS_PHY_TB, &
       do_phy_mp => TIME_DOATMOS_PHY_MP, &
       do_phy_rd => TIME_DOATMOS_PHY_RD
    use mod_atmos_vars, only: &
       sw_dyn    => ATMOS_sw_dyn,    &
       sw_phy_tb => ATMOS_sw_phy_tb, &
       sw_phy_mp => ATMOS_sw_phy_mp, &
       sw_phy_rd => ATMOS_sw_phy_rd, &
       ATMOS_vars_history
    use mod_atmos_dyn, only: &
       ATMOS_DYN
    use mod_atmos_phy_sf, only: &
       ATMOS_PHY_SF
    use mod_atmos_phy_tb, only: &
       ATMOS_PHY_TB
!    use mod_atmos_phy_mp, only: &
!       ATMOS_PHY_MP
!    use mod_atmos_phy_rd, only: &
!       ATMOS_PHY_RD
    implicit none
    !---------------------------------------------------------------------------

    !########## Dynamics ##########
    call TIME_rapstart('Dynamics')
#ifdef _FPCOLL_
call START_COLLECTION("Dynamics")
#endif
    if ( sw_dyn .AND. do_dyn ) then
       call ATMOS_DYN
    endif
#ifdef _FPCOLL_
call STOP_COLLECTION  ("Dynamics")
#endif
    call TIME_rapend  ('Dynamics')

    !########## Turbulence ##########

    call TIME_rapstart('Turbulence')
#ifdef _FPCOLL_
call START_COLLECTION("Turbulence")
#endif
    if ( sw_phy_tb .AND. do_phy_tb ) then
       call ATMOS_PHY_SF
       call ATMOS_PHY_TB
    endif
#ifdef _FPCOLL_
call STOP_COLLECTION  ("Turbulence")
#endif
    call TIME_rapend  ('Turbulence')

    !########## Microphysics ##########
    call TIME_rapstart('Microphysics')
#ifdef _FPCOLL_
call START_COLLECTION("Microphysics")
#endif
    if ( sw_phy_mp .AND. do_phy_mp ) then
!       call ATMOS_PHY_MP
    endif
#ifdef _FPCOLL_
call STOP_COLLECTION  ("Microphysics")
#endif
    call TIME_rapend  ('Microphysics')

    !########## Radiation ##########
    call TIME_rapstart('Radiation')
#ifdef _FPCOLL_
call START_COLLECTION("Radiation")
#endif
    if ( sw_phy_rd .AND. do_phy_rd ) then
!       call ATMOS_PHY_RD
    endif
#ifdef _FPCOLL_
call STOP_COLLECTION  ("Radiation")
#endif
    call TIME_rapend  ('Radiation')

    !########## History&Monitor ##########
    call TIME_rapstart('History')
#ifdef _FPCOLL_
call START_COLLECTION("History")
#endif
       call ATMOS_vars_history
#ifdef _FPCOLL_
call STOP_COLLECTION  ("History")
#endif
    call TIME_rapend  ('History')

    return
  end subroutine ATMOS_step

end module mod_atmos
