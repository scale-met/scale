!-------------------------------------------------------------------------------
!> module LAND / Physics
!!
!! @par Description
!!          land physics module
!!
!! @author Team SCALE
!<
!-------------------------------------------------------------------------------
module scale_land_phy
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_stdio
  use scale_prof
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: LAND_PHY_setup

  abstract interface
     subroutine lnd( &
           LAND_TEMP_t,       &
           LAND_WATER_t,      &
           LAND_TEMP,         &
           LAND_WATER,        &
           LAND_WaterLimit,   &
           LAND_ThermalCond,  &
           LAND_HeatCapacity, &
           LAND_WaterDiff,    &
           LAND_SFLX_GH,      &
           LAND_SFLX_prec,    &
           LAND_SFLX_evap,    &
           CDZ,               &
           dt                 )
       use scale_precision
       use scale_grid_index
       use scale_land_grid_index
       implicit none

       real(RP), intent(out) :: LAND_TEMP_t      (LKMAX,IA,JA)
       real(RP), intent(out) :: LAND_WATER_t     (LKMAX,IA,JA)

       real(RP), intent(in)  :: LAND_TEMP        (LKMAX,IA,JA)
       real(RP), intent(in)  :: LAND_WATER       (LKMAX,IA,JA)
       real(RP), intent(in)  :: LAND_WaterLimit  (IA,JA)
       real(RP), intent(in)  :: LAND_ThermalCond (IA,JA)
       real(RP), intent(in)  :: LAND_HeatCapacity(IA,JA)
       real(RP), intent(in)  :: LAND_WaterDiff   (IA,JA)
       real(RP), intent(in)  :: LAND_SFLX_GH     (IA,JA)
       real(RP), intent(in)  :: LAND_SFLX_prec   (IA,JA)
       real(RP), intent(in)  :: LAND_SFLX_evap   (IA,JA)
       real(RP), intent(in)  :: CDZ              (LKMAX)
       real(RP), intent(in)  :: dt
     end subroutine lnd
  end interface
  procedure(lnd), pointer :: LAND_PHY => NULL()
  public :: LAND_PHY

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
  subroutine LAND_PHY_setup( LAND_TYPE )
    use scale_process, only: &
       PRC_MPIstop
    use scale_land_phy_slab, only: &
       LAND_PHY_SLAB_setup, &
       LAND_PHY_SLAB
    use scale_land_phy_matsiro, only: &
       LAND_PHY_matsiro_setup, &
       LAND_PHY_matsiro
    implicit none

    character(len=*), intent(in) :: LAND_TYPE
    !---------------------------------------------------------------------------

    select case ( LAND_TYPE )
    case ( 'CONST' )
       call LAND_PHY_SLAB_setup( LAND_TYPE )
       LAND_PHY => LAND_PHY_SLAB
    case ( 'SLAB' )
       call LAND_PHY_SLAB_setup( LAND_TYPE )
       LAND_PHY => LAND_PHY_SLAB
    case ( 'MATSIRO' )
       call LAND_PHY_MATSIRO_setup( LAND_TYPE )
       LAND_PHY => LAND_PHY_MATSIRO
    case default
       write(*,*) 'xxx invalid Land type(', trim(LAND_TYPE), '). CHECK!'
       call PRC_MPIstop
    end select

    return
  end subroutine LAND_PHY_setup

end module scale_land_phy
