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
          OCEAN_TEMP,   &
          FLX_heat,     &
          FLX_precip,   &
          FLX_evap,     &
          OCEAN_TEMP_t  )
       use scale_precision
       use scale_grid_index
       implicit none

       real(RP), intent(in)  :: OCEAN_TEMP  (IA,JA)
       real(RP), intent(in)  :: FLX_heat    (IA,JA)
       real(RP), intent(in)  :: FLX_precip  (IA,JA)
       real(RP), intent(in)  :: FLX_evap    (IA,JA)
       real(RP), intent(out) :: OCEAN_TEMP_t(IA,JA)
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
#define EXTM(pre, name, post) pre ## name ## post
#define NAME(pre, name, post) EXTM(pre, name, post)
#ifdef OCN
    use NAME(scale_ocean_phy_, OCN,), only: &
       NAME(OCEAN_PHY_, OCN, _setup), &
       NAME(OCEAN_PHY_, OCN,)
#else
    use scale_ocean_phy_slab, only: &
       OCEAN_PHY_slab_setup, &
       OCEAN_PHY_slab
#endif
    implicit none

    character(len=*), intent(in) :: OCEAN_TYPE
    !---------------------------------------------------------------------------

    select case ( OCEAN_TYPE )
    case ( 'SLAB' )
       call OCEAN_PHY_slab_setup( OCEAN_TYPE )
       OCEAN_PHY => OCEAN_PHY_slab
    case default
       write(*,*) 'xxx invalid Ocean type(', trim(OCEAN_TYPE), '). CHECK!'
       call PRC_MPIstop
    end select

    return
  end subroutine OCEAN_PHY_setup

end module scale_ocean_phy
