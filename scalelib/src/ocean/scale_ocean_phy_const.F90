!-------------------------------------------------------------------------------
!> module OCEAN / Physics Constant model
!!
!! @par Description
!!          ocean physics module, constant model
!!
!! @author Team SCALE
!!
!<
!-------------------------------------------------------------------------------
module scale_ocean_phy_const
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
  public :: OCEAN_PHY_CONST_setup
  public :: OCEAN_PHY_CONST

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
  subroutine OCEAN_PHY_CONST_setup( OCEAN_TYPE )
    implicit none

    character(len=*), intent(in) :: OCEAN_TYPE
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '++++++ Module[CONST] / Categ[OCEAN PHY] / Origin[SCALElib]'

    return
  end subroutine OCEAN_PHY_CONST_setup

  !-----------------------------------------------------------------------------
  !> Slab ocean model
  subroutine OCEAN_PHY_CONST( &
       OCEAN_TEMP_t,    &
       OCEAN_TEMP,      &
       OCEAN_SFLX_WH,   &
       OCEAN_SFLX_prec, &
       OCEAN_SFLX_evap, &
       dt               )
    implicit none

    real(RP), intent(out) :: OCEAN_TEMP_t   (IA,JA)
    real(RP), intent(in)  :: OCEAN_TEMP     (IA,JA)
    real(RP), intent(in)  :: OCEAN_SFLX_WH  (IA,JA)
    real(RP), intent(in)  :: OCEAN_SFLX_prec(IA,JA)
    real(RP), intent(in)  :: OCEAN_SFLX_evap(IA,JA)
    real(DP), intent(in)  :: dt
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*) '*** Ocean physics step: Const'

    OCEAN_TEMP_t(:,:) = 0.0_RP

    return
  end subroutine OCEAN_PHY_CONST

end module scale_ocean_phy_const
