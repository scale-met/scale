!-------------------------------------------------------------------------------
!> module atmosphere / bottom boundary extrapolation
!!
!! @par Description
!!          Bottom boundary treatment of model domain
!!          Extrapolation of DENS, PRES, TEMP
!!
!! @author Team SCALE
!<
!-------------------------------------------------------------------------------
#include "scalelib.h"
module scale_atmos_bottom
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_io
  use scale_prof
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: ATMOS_BOTTOM_estimate

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
  !> Calc bottom boundary of atmosphere (just above surface)
  subroutine ATMOS_BOTTOM_estimate( &
       KA, KS, KE, IA, IS, IE, JA, JS, JE, &
       DENS, PRES, QV, &
       SFC_TEMP,       &
       Z1,             &
       SFC_DENS, SFC_PRES )
    use scale_const, only: &
       GRAV    => CONST_GRAV, &
       Rdry    => CONST_Rdry, &
       EPSTvap => CONST_EPSTvap
    implicit none
    integer, intent(in) :: KA, KS, KE
    integer, intent(in) :: IA, IS, IE
    integer, intent(in) :: JA, JS, JE

    real(RP), intent(in)  :: DENS    (KA,IA,JA)
    real(RP), intent(in)  :: PRES    (KA,IA,JA)
    real(RP), intent(in)  :: QV      (KA,IA,JA)
    real(RP), intent(in)  :: SFC_TEMP(IA,JA)
    real(RP), intent(in)  :: Z1      (IA,JA)
    real(RP), intent(out) :: SFC_DENS(IA,JA)
    real(RP), intent(out) :: SFC_PRES(IA,JA)

    real(RP) :: Rtot

    integer :: i, j
    !---------------------------------------------------------------------------

    !$omp parallel do OMP_SCHEDULE_ &
    !$omp private (Rtot)
    do j = JS, JE
    do i = IS, IE
       Rtot = Rdry * ( 1.0_RP + EPSTvap * QV(KS,i,j) )
       ! ( PRES(KS) - ( SFC_DENS * Rtot * SFC_TEMP ) ) / Z1 = - GRAV * ( DENS(KS) + SFC_DENS ) * 0.5
       SFC_DENS(i,j) = ( PRES(KS,i,j) + GRAV * Z1(i,j) * DENS(KS,i,j) * 0.5_RP ) &
                     / ( Rtot * SFC_TEMP(i,j) - GRAV * Z1(i,j) * 0.5_RP )
       SFC_PRES(i,j) = SFC_DENS(i,j) * Rtot * SFC_TEMP(i,j)
    enddo
    enddo

    return
  end subroutine ATMOS_BOTTOM_estimate

  !-----------------------------------------------------------------------------
  function lagrange_interp( p, x, y ) result(q)
    implicit none

    real(RP), intent(in) :: p
    real(RP), intent(in) :: x(3), y(3)
    real(RP)             :: q
    !---------------------------------------------------------------------------

    q = ( (p-x(2)) * (p-x(3)) ) / ( (x(1)-x(2)) * (x(1)-x(3)) ) * y(1) &
      + ( (p-x(1)) * (p-x(3)) ) / ( (x(2)-x(1)) * (x(2)-x(3)) ) * y(2) &
      + ( (p-x(1)) * (p-x(2)) ) / ( (x(3)-x(1)) * (x(3)-x(2)) ) * y(3)

    return
  end function lagrange_interp

end module scale_atmos_bottom
