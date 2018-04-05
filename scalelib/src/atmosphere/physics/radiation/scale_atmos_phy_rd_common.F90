!-------------------------------------------------------------------------------
!> module atmosphere / physics / radiation / common
!!
!! @par Description
!!          Common module for Radiation
!!
!! @author Team SCALE
!!
!<
!-------------------------------------------------------------------------------
#include "scalelib.h"
module scale_atmos_phy_rd_common
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
  public :: ATMOS_PHY_RD_calc_heating

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  ! for direction
  integer,  public, parameter :: I_up = 1
  integer,  public, parameter :: I_dn = 2
  ! for band region
  integer,  public, parameter :: I_LW = 1
  integer,  public, parameter :: I_SW = 2
  ! for band region
  integer,  public, parameter :: I_direct  = 1
  integer,  public, parameter :: I_diffuse = 2

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
  !> Calc heating rate
  subroutine ATMOS_PHY_RD_calc_heating( &
       KA, KS, KE, IA, IS, IE, JA, JS, JE, &
       flux_rad,     &
       DENS, TEMP,   &
       CVtot,        &
       FZ,           &
       RHOH,         &
       TEMP_t        )
    implicit none
    integer, intent(in) :: KA, KS, KE
    integer, intent(in) :: IA, IS, IE
    integer, intent(in) :: JA, JS, JE

    real(RP), intent(in)  :: flux_rad(KA,IA,JA,2,2)
    real(RP), intent(in)  :: DENS    (KA,IA,JA)
    real(RP), intent(in)  :: TEMP    (KA,IA,JA)
    real(RP), intent(in)  :: CVtot   (KA,IA,JA)
    real(RP), intent(in)  :: FZ      (0:KA,IA,JA)

    real(RP), intent(out) :: RHOH(KA,IA,JA)

    real(RP), intent(out), optional :: TEMP_t(KA,IA,JA,3)

    real(RP) :: RHOH_LW, RHOH_SW

    integer :: k, i, j
    !---------------------------------------------------------------------------

    !$omp parallel do default(none) OMP_SCHEDULE_ collapse(2) &
    !$omp private(i,j,k, &
    !$omp         RHOH_LW,RHOH_SW) &
    !$omp shared(RHOH,TEMP_t,flux_rad,DENS,CVtot,FZ, &
    !$omp        KS,KE,IS,IE,JS,JE)
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE

       RHOH_LW = ( ( flux_rad(k,i,j,I_LW,I_dn) - flux_rad(k-1,i,j,I_LW,I_dn) ) &
                 - ( flux_rad(k,i,j,I_LW,I_up) - flux_rad(k-1,i,j,I_LW,I_up) ) &
                 ) / ( FZ(k,i,j) - FZ(k-1,i,j) )

       RHOH_SW = ( ( flux_rad(k,i,j,I_SW,I_dn) - flux_rad(k-1,i,j,I_SW,I_dn) ) &
                 - ( flux_rad(k,i,j,I_SW,I_up) - flux_rad(k-1,i,j,I_SW,I_up) ) &
                 ) / ( FZ(k,i,j) - FZ(k-1,i,j) )

       RHOH(k,i,j) = RHOH_LW + RHOH_SW

       TEMP_t(k,i,j,I_LW) = RHOH_LW / DENS(k,i,j) / CVtot(k,i,j) * 86400.0_RP ! [K/day]
       TEMP_t(k,i,j,I_SW) = RHOH_SW / DENS(k,i,j) / CVtot(k,i,j) * 86400.0_RP ! [K/day]

       TEMP_t(k,i,j,3)    = TEMP_t(k,i,j,I_LW) + TEMP_t(k,i,j,I_SW)
    enddo
    enddo
    enddo

    return
  end subroutine ATMOS_PHY_RD_calc_heating

end module scale_atmos_phy_rd_common
