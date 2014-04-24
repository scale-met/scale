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
  use scale_precision
  use scale_stdio
  use scale_prof
  use scale_grid_index
  use scale_tracer
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
    use scale_atmos_phy_ae, only: &
       ATMOS_PHY_AE_setup
    implicit none

    character(len=H_SHORT), intent(in) :: AE_TYPE

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
    use scale_atmos_phy_ae, only: &
       ATMOS_PHY_AE
    use mod_atmos_vars, only: &
       ATMOS_vars_fillhalo, &
       ATMOS_vars_total, &
       DENS, &
       MOMZ, &
       MOMX, &
       MOMY, &
       RHOT, &
       QTRC
    implicit none
    logical, intent(in) :: update_flag
    logical, intent(in) :: history_flag

    if ( update_flag ) then
       call ATMOS_PHY_AE( &
            DENS, &
            MOMZ, &
            MOMX, &
            MOMY, &
            RHOT, &
            QTRC  )

       call ATMOS_vars_fillhalo

       call ATMOS_vars_total
    endif

    return
  end subroutine ATMOS_PHY_AE_driver

end module mod_atmos_phy_ae_driver
