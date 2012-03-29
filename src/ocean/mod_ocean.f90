!-------------------------------------------------------------------------------
!> module OCEAN
!!
!! @par Description
!!          Ocean module
!!
!! @author H.Tomita and SCALE developpers
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
  !> Setup ocean
  !-----------------------------------------------------------------------------
  subroutine OCEAN_setup
    use mod_ocean_vars, only: &
       OCEAN_vars_setup, &
       OCEAN_vars_restart_read
    use mod_ocean_sf, only: &
       OCEAN_FIXEDSST_setup
    implicit none
    !---------------------------------------------------------------------------

    call OCEAN_vars_setup

    call OCEAN_vars_restart_read

    call OCEAN_FIXEDSST_setup

    return
  end subroutine OCEAN_setup

  !-----------------------------------------------------------------------------
  !> advance ocean state
  !-----------------------------------------------------------------------------
  subroutine OCEAN_step
    use mod_ocean_vars, only: &
       sw_sf => OCEAN_sw_sf
    use mod_ocean_sf, only: &
       OCEAN_FIXEDSST
    implicit none
    !---------------------------------------------------------------------------

    call TIME_rapstart('Ocean')
    if ( sw_sf ) then
       call OCEAN_FIXEDSST
    endif
    call TIME_rapend  ('Ocean')

    return
  end subroutine OCEAN_step

end module mod_ocean
