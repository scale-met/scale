!-------------------------------------------------------------------------------
!> module ATMOSPHERE / Physics Radiation
!!
!! @par Description
!!          Common module for Radiation
!!          Heating Rate
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2013-02-06 (H.Yashiro) [new]
!!
!<
!-------------------------------------------------------------------------------
module mod_atmos_phy_rd_common
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use mod_stdio, only: &
     IO_FID_LOG,  &
     IO_L
  use mod_time, only: &
     TIME_rapstart, &
     TIME_rapend
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ included parameters
  !
  include "inc_precision.h"
  include 'inc_index.h'
  include 'inc_tracer.h'

  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: ATMOS_PHY_RD_heating

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

  ! for direction
  integer,  private, parameter :: I_up = 1
  integer,  private, parameter :: I_dn = 2
  ! for band region
  integer,  private, parameter :: I_LW = 1
  integer,  private, parameter :: I_SW = 2

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Calc heating rate
  subroutine ATMOS_PHY_RD_heating( &
       RHOE_t,  &
       RHOE,    &
       flux_rad )
    use mod_time, only: &
       dt => TIME_DTSEC_ATMOS_PHY_RD
    use mod_grid, only: &
       RCDZ => GRID_RCDZ
    implicit none

    real(RP), intent(out)   :: RHOE_t  (KA,IA,JA,2)
    real(RP), intent(inout) :: RHOE    (KA,IA,JA)
    real(RP), intent(in)    :: flux_rad(KA,IA,JA,2,2)

    integer :: k, i, j
    !---------------------------------------------------------------------------

    do j = JS, JE
    do i = IS, IE
    do k = KS, KE

       RHOE_t(k,i,j,I_LW) = ( ( flux_rad(k,i,j,I_LW,I_dn) - flux_rad(k-1,i,j,I_LW,I_dn) ) &
                            - ( flux_rad(k,i,j,I_LW,I_up) - flux_rad(k-1,i,j,I_LW,I_up) ) ) * RCDZ(k)

       RHOE_t(k,i,j,I_SW) = ( ( flux_rad(k,i,j,I_SW,I_dn) - flux_rad(k-1,i,j,I_SW,I_dn) ) &
                            - ( flux_rad(k,i,j,I_SW,I_up) - flux_rad(k-1,i,j,I_SW,I_up) ) ) * RCDZ(k)

       RHOE(k,i,j) = RHOE(k,i,j) + dt * ( RHOE_t(k,i,j,I_LW) + RHOE_t(k,i,j,I_SW) )

    enddo
    enddo
    enddo

    return
  end subroutine ATMOS_PHY_RD_heating

end module mod_atmos_phy_rd_common
