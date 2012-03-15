!-------------------------------------------------------------------------------
!> module OCEAN
!!
!! @par Description
!!          Ocean module
!!
!! @author H.Tomita and SCALE developpers
!!
!! @par History
!! @li      2011-12-11 (H.Yashiro) [new]
!!
!<
!-------------------------------------------------------------------------------
module mod_ocean
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
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
  !> Setup ocean
  !-----------------------------------------------------------------------------
  subroutine OCEAN_setup
    use mod_stdio, only: &
       IO_FID_LOG,  &
       IO_L
    use mod_ocean_vars, only: &
       OCEAN_vars_setup,  &
       OCEAN_vars_restart_read
    implicit none
    !---------------------------------------------------------------------------

    call OCEAN_vars_setup

!    call OCEAN_vars_restart_read

    return
  end subroutine OCEAN_setup

  !-----------------------------------------------------------------------------
  !> advance ocean state
  !-----------------------------------------------------------------------------
  subroutine OCEAN_step
    implicit none
    !---------------------------------------------------------------------------

    ! Do nothing yet

    return
  end subroutine OCEAN_step

end module mod_ocean
