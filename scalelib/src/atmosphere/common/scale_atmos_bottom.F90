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
#include "inc_openmp.h"
module scale_atmos_bottom
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
       DENS, PRES,        &
       CZ, Zsfc, Z1,      &
       SFC_DENS, SFC_PRES )
    use scale_const, only: &
       GRAV  => CONST_GRAV
    implicit none
    integer, intent(in) :: KA, KS, KE
    integer, intent(in) :: IA, IS, IE
    integer, intent(in) :: JA, JS, JE

    real(RP), intent(in)  :: DENS    (KA,IA,JA)
    real(RP), intent(in)  :: PRES    (KA,IA,JA)
    real(RP), intent(in)  :: CZ      (KA,IA,JA)
    real(RP), intent(in)  :: Zsfc    (IA,JA)
    real(RP), intent(in)  :: Z1      (IA,JA)
    real(RP), intent(out) :: SFC_DENS(IA,JA)
    real(RP), intent(out) :: SFC_PRES(IA,JA)

    integer :: i, j
    !---------------------------------------------------------------------------

    ! estimate surface density (extrapolation)
    !$omp parallel do OMP_SCHEDULE_
    do j = JS, JE
    do i = IS, IE
       SFC_DENS(i,j) = lagrange_interp( Zsfc(i,j),         & ! [IN]
                                        CZ  (KS:KS+2,i,j), & ! [IN]
                                        DENS(KS:KS+2,i,j)  ) ! [IN]
    enddo
    enddo

    ! estimate surface pressure (hydrostatic balance)
    !$omp parallel do OMP_SCHEDULE_
    do j = JS, JE
    do i = IS, IE
       SFC_PRES(i,j) = PRES(KS,i,j) &
                     + 0.5_RP * ( SFC_DENS(i,j) + DENS(KS,i,j) ) * GRAV * Z1(i,j)
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
