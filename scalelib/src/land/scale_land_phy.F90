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
     subroutine land( &
          LAND_TEMP,         &
          LAND_WATER,        &
          LAND_WaterLimit,   &
          LAND_ThermalCond,  &
          LAND_HeatCapacity, &
          LAND_WaterDiff,    &
          FLX_heat,          &
          FLX_precip,        &
          FLX_evap,          &
          CDZ,               &
          LAND_TEMP_t,       &
          LAND_WATER_t       )
       use scale_precision
       use scale_grid_index
       use scale_land_grid_index
       implicit none

       real(RP), intent(in)  :: LAND_TEMP        (LKMAX,IA,JA)
       real(RP), intent(in)  :: LAND_WATER       (LKMAX,IA,JA)
       real(RP), intent(in)  :: LAND_WaterLimit  (IA,JA)
       real(RP), intent(in)  :: LAND_ThermalCond (IA,JA)
       real(RP), intent(in)  :: LAND_HeatCapacity(IA,JA)
       real(RP), intent(in)  :: LAND_WaterDiff   (IA,JA)
       real(RP), intent(in)  :: FLX_heat         (IA,JA)
       real(RP), intent(in)  :: FLX_precip       (IA,JA)
       real(RP), intent(in)  :: FLX_evap         (IA,JA)
       real(RP), intent(in)  :: CDZ              (LKMAX)
       real(RP), intent(out) :: LAND_TEMP_t      (LKMAX,IA,JA)
       real(RP), intent(out) :: LAND_WATER_t     (LKMAX,IA,JA)
     end subroutine land
  end interface
  procedure(land), pointer :: LAND_PHY => NULL()
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
#define EXTM(pre, name, post) pre ## name ## post
#define NAME(pre, name, post) EXTM(pre, name, post)
#ifdef RD
    use NAME(scale_land_phy_, LAND,), only: &
       NAME(LAND_PHY_, LAND, _setup), &
       NAME(LAND_PHY_, LAND,)
#else
    use scale_land_phy_bucket, only: &
       LAND_PHY_bucket_setup, &
       LAND_PHY_bucket
    use scale_land_phy_matsiro, only: &
       LAND_PHY_matsiro_setup, &
       LAND_PHY_matsiro
#endif
    implicit none

    character(len=*), intent(in) :: LAND_TYPE
    !---------------------------------------------------------------------------

    select case ( LAND_TYPE )
    case ( 'BUCKET' )
       call LAND_PHY_bucket_setup( LAND_TYPE )
       LAND_PHY => LAND_PHY_bucket
    case ( 'MATSIRO' )
       call LAND_PHY_matsiro_setup( LAND_TYPE )
       LAND_PHY => LAND_PHY_matsiro
    case default
       write(*,*) 'xxx invalid Land type(', trim(LAND_TYPE), '). CHECK!'
       call PRC_MPIstop
    end select

    return
  end subroutine LAND_PHY_setup

end module scale_land_phy
