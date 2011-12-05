!-------------------------------------------------------------------------------
!> module Ocean
!!
!! @par Description
!!          Ocean module
!!
!! @author H.Tomita and SCALE developpers
!!
!! @par History
!! @li      2011-11-11 (H.Yashiro) [new] Imported from SCALE-LES ver.2
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
    use mod_grid, only: &
       IA   => GRID_IA, &
       JA   => GRID_JA
    use mod_ocean_vars, only: &
       OCEAN_vars_get, &
       OCEAN_vars_put
    implicit none

    real(8) :: sst(IA,JA,1)    ! Sea Surface Temperatue [K]
    !---------------------------------------------------------------------------

!    call OCEAN_vars_get( sst )

    ! Do nothing yet

!    call OCEAN_vars_put( sst )

    return
  end subroutine OCEAN_step

end module mod_ocean
