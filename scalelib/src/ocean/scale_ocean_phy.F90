!-------------------------------------------------------------------------------
!> module OCEAN / Physics
!!
!! @par Description
!!          ocean physics module
!!
!! @author Team SCALE
!<
!-------------------------------------------------------------------------------
module scale_ocean_phy
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
  public :: OCEAN_PHY_setup

  abstract interface
     subroutine ocn( &
           OCEAN_TEMP_t,    &
           OCEAN_TEMP,      &
           OCEAN_SFLX_WH,   &
           OCEAN_SFLX_prec, &
           OCEAN_SFLX_evap, &
           dt               )
       use scale_precision
       use scale_ocean_grid_cartesC_index
       implicit none

       real(RP), intent(out) :: OCEAN_TEMP_t   (OKMAX,OIA,OJA)
       real(RP), intent(in)  :: OCEAN_TEMP     (OKMAX,OIA,OJA)
       real(RP), intent(in)  :: OCEAN_SFLX_WH  (OIA,OJA)
       real(RP), intent(in)  :: OCEAN_SFLX_prec(OIA,OJA)
       real(RP), intent(in)  :: OCEAN_SFLX_evap(OIA,OJA)
       real(DP), intent(in)  :: dt
     end subroutine ocn
  end interface
  procedure(ocn), pointer :: OCEAN_PHY => NULL()
  public :: OCEAN_PHY

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
  subroutine OCEAN_PHY_setup( OCEAN_TYPE )
    use scale_prc, only: &
       PRC_abort
    use scale_ocean_phy_const, only: &
       OCEAN_PHY_CONST_setup, &
       OCEAN_PHY_CONST
    use scale_ocean_phy_slab, only: &
       OCEAN_PHY_SLAB_setup, &
       OCEAN_PHY_SLAB
    implicit none

    character(len=*), intent(in) :: OCEAN_TYPE
    !---------------------------------------------------------------------------

    select case( OCEAN_TYPE )
    case( 'CONST' )
       call OCEAN_PHY_CONST_setup( OCEAN_TYPE )
       OCEAN_PHY => OCEAN_PHY_CONST
    case( 'SLAB' )
       call OCEAN_PHY_SLAB_setup( OCEAN_TYPE )
       OCEAN_PHY => OCEAN_PHY_SLAB
    case default
       write(*,*) 'xxx invalid Ocean type(', trim(OCEAN_TYPE), '). CHECK!'
       call PRC_abort
    end select

    return
  end subroutine OCEAN_PHY_setup

end module scale_ocean_phy
