!-------------------------------------------------------------------------------
!> module ATMOSPHERE / Physics Surface fluxes
!!
!! @par Description
!!          Flux from/to bottom wall of atmosphere (surface)
!!          Bulk Method
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2013-08-31 (T.Yamaura)  [new]
!<
!-------------------------------------------------------------------------------
module mod_atmos_phy_sf
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use mod_precision
  use mod_index
  use mod_tracer
  use mod_stdio, only: &
     IO_FID_LOG,  &
     IO_L
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: ATMOS_PHY_SF_setup
  public :: ATMOS_PHY_SF

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
  !
  !-----------------------------------------------------------------------------
  subroutine ATMOS_PHY_SF_setup
    implicit none

    !---------------------------------------------------------------------------

    return
  end subroutine ATMOS_PHY_SF_setup

  !-----------------------------------------------------------------------------
  ! calculation flux
  !-----------------------------------------------------------------------------
  subroutine ATMOS_PHY_SF( &
       update_flag, &
       history_flag )
    use mod_const, only: &
       CPdry  => CONST_CPdry,  &
       LH0    => CONST_LH0
    use mod_grid, only: &
       RCDZ => GRID_RCDZ, &
       RFDZ => GRID_RFDZ
    use mod_atmos_vars, only: &
       DENS,    &
       RHOT,    &
       DENS_tp, &
       MOMZ_tp, &
       MOMX_tp, &
       MOMY_tp, &
       RHOT_tp, &
       QTRC_tp
    use mod_atmos_vars_sf, only: &
       SFLX_MOMZ, &
       SFLX_MOMX, &
       SFLX_MOMY, &
       SFLX_SH,   &
       SFLX_LH,   &
       SFLX_QVAtm
    implicit none

    logical, intent(in) :: update_flag
    logical, intent(in), optional :: history_flag

    integer :: i, j

    if( IO_L ) write(IO_FID_LOG,*) '*** Physics step: Surface'

    do j = JS, JE
    do i = IS, IE
       RHOT_tp(KS,i,j) = RHOT_tp(KS,i,j) &
            + ( SFLX_SH(i,j)/CPdry &
              + SFLX_QVAtm(i,j) * RHOT(KS,i,j) / DENS(KS,i,j) &
              ) * RCDZ(KS)
       DENS_tp(KS,i,j) = DENS_tp(KS,i,j) &
            + SFLX_QVAtm(i,j)  * RCDZ(KS)
       MOMZ_tp(KS,i,j) = MOMZ_tp(KS,i,j) &
            + SFLX_MOMZ(i,j)   * RFDZ(KS)
       MOMX_tp(KS,i,j) = MOMX_tp(KS,i,j) &
            + SFLX_MOMX(i,j)   * RCDZ(KS)
       MOMY_tp(KS,i,j) = MOMY_tp(KS,i,j) &
            + SFLX_MOMY(i,j)   * RCDZ(KS)
       QTRC_tp(KS,i,j,I_QV) = QTRC_tp(KS,i,j,I_QV) &
            + SFLX_QVAtm(i,j)  * RCDZ(KS)
    enddo
    enddo

    return
  end subroutine ATMOS_PHY_SF

end module mod_atmos_phy_sf
