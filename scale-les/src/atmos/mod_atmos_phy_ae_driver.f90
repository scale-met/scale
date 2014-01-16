!-------------------------------------------------------------------------------
!> module ATMOSPHERE / Physics Aerosol Microphysics
!!
!! @par Description
!!          Aerosol Microphysics wrapper
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2013-12-06 (S.Nishizawa)  [new]
!!
!<
!-------------------------------------------------------------------------------
module mod_atmos_phy_ae_driver
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use mod_precision
  use mod_index
  use mod_stdio
  use mod_prof
  use mod_tracer
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: ATMOS_PHY_AE_driver_setup
  public :: ATMOS_PHY_AE_driver

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
  subroutine ATMOS_PHY_AE_driver_setup( AE_TYPE )
    use mod_stdio, only: &
       IO_FID_LOG, &
       IO_L, &
       IO_SYSCHR
    use mod_atmos_phy_ae, only: &
       ATMOS_PHY_AE_setup
    implicit none
    character(len=IO_SYSCHR), intent(in) :: AE_TYPE

    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[Physics-AE]/Categ[ATMOS]'

    call ATMOS_PHY_AE_setup( AE_TYPE )

    call ATMOS_PHY_AE_driver( .true., .false. )

    return
  end subroutine ATMOS_PHY_AE_driver_setup

  !-----------------------------------------------------------------------------
  !> Aerosol Microphysics
  subroutine ATMOS_PHY_AE_driver( update_flag, history_flag )
    use mod_atmos_phy_ae, only: &
       ATMOS_PHY_AE
    implicit none
    logical, intent(in) :: update_flag
    logical, intent(in) :: history_flag

    if ( update_flag ) then
       call ATMOS_PHY_AE
    end if

    return
  end subroutine ATMOS_PHY_AE_driver

end module mod_atmos_phy_ae_driver
