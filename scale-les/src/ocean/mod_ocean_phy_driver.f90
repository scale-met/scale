!-------------------------------------------------------------------------------
!> module OCEAN / Physics
!!
!! @par Description
!!          ocean physics module
!!
!! @author Team SCALE
!!
!<
!-------------------------------------------------------------------------------
module mod_ocean_phy_driver
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_stdio
  use scale_prof
  use scale_debug
  use scale_grid_index
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: OCEAN_PHY_driver_setup
  public :: OCEAN_PHY_driver

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
  subroutine OCEAN_PHY_driver_setup
    use scale_ocean_phy, only: &
       OCEAN_PHY_setup
    use mod_ocean_admin, only: &
       OCEAN_TYPE, &
       OCEAN_sw
    implicit none
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '++++++ Module[DRIVER] / Categ[OCEAN PHY] / Origin[SCALE-LES]'

    if ( OCEAN_sw ) then

       ! setup library component
       call OCEAN_PHY_setup( OCEAN_TYPE )

       ! run once (only for the diagnostic value)
       call PROF_rapstart('OCN Physics', 1)
       call OCEAN_PHY_driver( update_flag = .true. )
       call PROF_rapend  ('OCN Physics', 1)

    else
       if( IO_L ) write(IO_FID_LOG,*) '*** this component is never called.'
    endif

    return
  end subroutine OCEAN_PHY_driver_setup

  !-----------------------------------------------------------------------------
  !> Driver
  subroutine OCEAN_PHY_driver( update_flag )
    use scale_time, only: &
       dt => TIME_DTSEC_OCEAN
    use scale_statistics, only: &
       STATISTICS_checktotal, &
       STAT_total
    use scale_history, only: &
       HIST_in
    use scale_ocean_phy, only: &
       OCEAN_PHY
    use mod_ocean_vars, only: &
       OCEAN_TEMP,      &
       OCEAN_TEMP_t,    &
       OCEAN_SFLX_heat, &
       OCEAN_SFLX_prec, &
       OCEAN_SFLX_evap
    implicit none

    logical, intent(in) :: update_flag

    real(RP) :: total ! dummy
    !---------------------------------------------------------------------------

    if ( update_flag ) then

       call OCEAN_PHY( OCEAN_TEMP     (:,:), & ! [IN]
                       OCEAN_SFLX_heat(:,:), & ! [IN]
                       OCEAN_SFLX_prec(:,:), & ! [IN]
                       OCEAN_SFLX_evap(:,:), & ! [IN]
                       OCEAN_TEMP_t   (:,:)  ) ! [OUT]

       call HIST_in( OCEAN_TEMP_t(:,:), 'OCEAN_TEMP_t', 'SST tendency', 'K' )

    endif

    if ( STATISTICS_checktotal ) then
       call STAT_total( total, OCEAN_TEMP_t(:,:), 'OCEAN_TEMP_t' )
    endif

    return
  end subroutine OCEAN_PHY_driver

end module mod_ocean_phy_driver
