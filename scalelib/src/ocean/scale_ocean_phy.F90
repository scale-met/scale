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
       use scale_grid_index
       implicit none

       real(RP), intent(out) :: OCEAN_TEMP_t   (IA,JA)
       real(RP), intent(in)  :: OCEAN_TEMP     (IA,JA)
       real(RP), intent(in)  :: OCEAN_SFLX_WH  (IA,JA)
       real(RP), intent(in)  :: OCEAN_SFLX_prec(IA,JA)
       real(RP), intent(in)  :: OCEAN_SFLX_evap(IA,JA)
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
    use scale_process, only: &
       PRC_MPIstop
    use scale_ocean_phy_slab, only: &
       OCEAN_PHY_SLAB_setup, &
       OCEAN_PHY_SLAB
    use scale_ocean_phy_file, only: &
       OCEAN_PHY_FILE_setup, &
       OCEAN_PHY_FILE
    implicit none

    character(len=*), intent(in) :: OCEAN_TYPE
    !---------------------------------------------------------------------------

    select case( OCEAN_TYPE )
    case( 'CONST' )
       call OCEAN_PHY_SLAB_setup( OCEAN_TYPE )
       OCEAN_PHY => OCEAN_PHY_SLAB
    case( 'SLAB' )
       call OCEAN_PHY_SLAB_setup( OCEAN_TYPE )
       OCEAN_PHY => OCEAN_PHY_SLAB
    case( 'FILE' )
       call OCEAN_PHY_FILE_setup( OCEAN_TYPE )
       OCEAN_PHY => OCEAN_PHY_FILE
    case default
       write(*,*) 'xxx invalid Ocean type(', trim(OCEAN_TYPE), '). CHECK!'
       call PRC_MPIstop
    end select

    return
  end subroutine OCEAN_PHY_setup

end module scale_ocean_phy
