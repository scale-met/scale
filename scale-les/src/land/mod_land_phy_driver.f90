!-------------------------------------------------------------------------------
!> module LAND / Physics
!!
!! @par Description
!!          land physics module
!!
!! @author Team SCALE
!<
!-------------------------------------------------------------------------------
module mod_land_phy_driver
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_stdio
  use scale_prof
  use scale_debug
  use scale_grid_index
  use scale_land_grid_index
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: LAND_PHY_driver_setup
  public :: LAND_PHY_driver

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
  subroutine LAND_PHY_driver_setup
    use scale_land_phy, only: &
       LAND_PHY_setup
    use mod_land_admin, only: &
       LAND_TYPE, &
       LAND_sw
    implicit none
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '++++++ Module[DRIVER] / Categ[LAND PHY] / Origin[SCALE-LES]'

    if ( LAND_sw ) then

       ! setup library component
       call LAND_PHY_setup( LAND_TYPE )

       ! run once (only for the diagnostic value)
       call LAND_PHY_driver( update_flag=.true., history_flag=.false. )

    else
       if( IO_L ) write(IO_FID_LOG,*) '*** this component is never called.'
    endif

    return
  end subroutine LAND_PHY_driver_setup

  !-----------------------------------------------------------------------------
  !> Driver
  subroutine LAND_PHY_driver( update_flag, history_flag )
    use scale_time, only: &
       dt => TIME_DTSEC_LAND
    use scale_land_grid, only: &
       GRID_LCDZ
    use scale_statistics, only: &
       STATISTICS_checktotal, &
       STAT_total
    use scale_history, only: &
       HIST_in
    use scale_land_phy, only: &
       LAND_PHY
    use mod_land_vars, only: &
       I_WaterLimit,   &
       I_ThermalCond,  &
       I_HeatCapacity, &
       I_WaterDiff,    &
       LAND_TEMP,      &
       LAND_WATER,     &
       LAND_TEMP_t,    &
       LAND_WATER_t,   &
       LAND_PROPERTY
    use mod_cpl_admin, only: &
       CPL_sw_AtmLnd
    use mod_cpl_vars, only: &
       CPL_getLnd
    implicit none

    logical, intent(in) :: update_flag
    logical, intent(in) :: history_flag

    real(RP) :: FLX_heat  (IA,JA)
    real(RP) :: FLX_precip(IA,JA)
    real(RP) :: FLX_evap  (IA,JA)

    real(RP) :: total ! dummy

    character(len=2) :: sk
    integer          :: k
    !---------------------------------------------------------------------------

    if ( update_flag ) then

       if ( CPL_sw_AtmLnd ) then
          call CPL_getLnd( FLX_heat  (:,:), & ! [OUT]
                           FLX_precip(:,:), & ! [OUT]
                           FLX_evap  (:,:)  ) ! [OUT]
       endif

       call LAND_PHY( LAND_TEMP    (:,:,:),              & ! [IN]
                      LAND_WATER   (:,:,:),              & ! [IN]
                      LAND_PROPERTY(:,:,I_WaterLimit),   & ! [IN]
                      LAND_PROPERTY(:,:,I_ThermalCond),  & ! [IN]
                      LAND_PROPERTY(:,:,I_HeatCapacity), & ! [IN]
                      LAND_PROPERTY(:,:,I_WaterDiff),    & ! [IN]
                      FLX_heat     (:,:),                & ! [IN]
                      FLX_precip   (:,:),                & ! [IN]
                      FLX_evap     (:,:),                & ! [IN]
                      GRID_LCDZ    (:),                  & ! [IN]
                      LAND_TEMP_t  (:,:,:),              & ! [OUT]
                      LAND_WATER_t (:,:,:)               ) ! [OUT]

       if ( history_flag ) then
          call HIST_in( LAND_TEMP_t (:,:,:), 'LAND_TEMP_t',  'Soil temperature tendency', 'K',     zdim='land' )
          call HIST_in( LAND_WATER_t(:,:,:), 'LAND_WATER_t', 'Soil moisture    tendency', 'm3/m3', zdim='land' )

          call HIST_in( LAND_PROPERTY(:,:,1), 'LP_WaterLimit   ', 'LAND PROPERTY, WaterLimit   ', 'm3/m3'  )
          call HIST_in( LAND_PROPERTY(:,:,2), 'LP_WaterCritical', 'LAND PROPERTY, WaterCritical', 'm3/m3'  )
          call HIST_in( LAND_PROPERTY(:,:,3), 'LP_ThermalCond  ', 'LAND PROPERTY, ThermalCond  ', 'W/K/m'  )
          call HIST_in( LAND_PROPERTY(:,:,4), 'LP_HeatCapacity ', 'LAND PROPERTY, HeatCapacity ', 'J/K/m3' )
          call HIST_in( LAND_PROPERTY(:,:,5), 'LP_WaterDiff    ', 'LAND PROPERTY, WaterDiff    ', 'm2/s'   )
          call HIST_in( LAND_PROPERTY(:,:,6), 'LP_Z0M          ', 'LAND PROPERTY, Z0M          ', 'm'      )
          call HIST_in( LAND_PROPERTY(:,:,7), 'LP_Z0H          ', 'LAND PROPERTY, Z0H          ', 'm'      )
          call HIST_in( LAND_PROPERTY(:,:,8), 'LP_Z0E          ', 'LAND PROPERTY, Z0E          ', 'm'      )
       endif

    endif

    if ( STATISTICS_checktotal ) then
       do k = LKS, LKE
          write(sk,'(I2.2)') k

          call STAT_total( total, LAND_TEMP_t (k,:,:), 'LAND_TEMP_t'//sk  )
          call STAT_total( total, LAND_WATER_t(k,:,:), 'LAND_WATER_t'//sk )
       enddo
    endif

    return
  end subroutine LAND_PHY_driver

end module mod_land_phy_driver
