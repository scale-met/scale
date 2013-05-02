!-------------------------------------------------------------------------------
!> module OCEAN driver
!!
!! @par Description
!!          Ocean module
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2011-12-11 (H.Yashiro)  [new]
!! @li      2012-03-23 (H.Yashiro)  [mod] FIXEDSST
!!
!<
!-------------------------------------------------------------------------------
module mod_ocean
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
  public :: OCEAN_setup
  public :: OCEAN_step

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
  subroutine OCEAN_setup
    use mod_ocean_vars, only: &
       OCEAN_vars_setup,        &
       OCEAN_vars_restart_read, &
       OCEAN_sw_phy
    use mod_ocean_phy, only: &
       OCEAN_PHY_SST_setup
    implicit none
    !---------------------------------------------------------------------------

    call OCEAN_vars_setup

    call OCEAN_vars_restart_read

    if ( OCEAN_sw_phy ) call OCEAN_PHY_SST_setup

    return
  end subroutine OCEAN_setup

  !-----------------------------------------------------------------------------
  !> Ocean step
  subroutine OCEAN_step
    use mod_ocean_phy, only: &
       OCEAN_PHY_SST
    use mod_ocean_vars, only: &
!       OCEAN_vars_history, &
       OCEAN_sw_phy
    implicit none
    !---------------------------------------------------------------------------

    call TIME_rapstart('OCN Physics')
    if ( OCEAN_sw_phy ) then
       call OCEAN_PHY_SST
    endif
    call TIME_rapend  ('OCN Physics')

    return
  end subroutine OCEAN_step

end module mod_ocean
