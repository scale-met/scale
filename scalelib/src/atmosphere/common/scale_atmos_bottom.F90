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
       FZ,             &
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
    real(RP), intent(in)  :: FZ      (0:KA,IA,JA)
    real(RP), intent(out) :: SFC_DENS(IA,JA)
    real(RP), intent(out) :: SFC_PRES(IA,JA)

    real(RP) :: Rtot
    real(RP) :: dz1, dz2
    real(RP) :: F2H1, F2H2
    real(RP) :: LogP1, LogP2
    real(RP) :: PRESH

    integer :: i, j
    !---------------------------------------------------------------------------

    !$acc data copyin(DENS, PRES, QV, SFC_TEMP, FZ) copyout(SFC_DENS, SFC_PRES)

    !$omp parallel do OMP_SCHEDULE_ &
    !$omp private (Rtot,dz1,dz2,F2H1,F2H2,LogP1,LogP2,PRESH)
    !$acc kernels
    do j = JS, JE
    do i = IS, IE
       ! interpolate half-level pressure
       dz1 = FZ(KS+1,i,j) - FZ(KS,i,j)
       dz2 = FZ(KS,i,j) - FZ(KS-1,i,j)

       F2H1 = dz2 / ( dz1 + dz2 )
       F2H2 = dz1 / ( dz1 + dz2 )

       LogP1 = log( PRES(KS+1,i,j) )
       LogP2 = log( PRES(KS,i,j) )

       PRESH = exp( F2H1 * LogP1 + F2H2 * LogP2 )

       Rtot = Rdry * ( 1.0_RP + EPSTvap * QV(KS,i,j) )
       ! ( PRESH - SFC_PRES ) / dz2 = - GRAV * DENS(KS)
       SFC_PRES(i,j) = PRESH + GRAV * DENS(KS,i,j) * dz2
       SFC_DENS(i,j) = SFC_PRES(i,j) / ( Rtot * SFC_TEMP(i,j) )
    enddo
    enddo
    !$acc end kernels

    !$acc end data

    return
  end subroutine ATMOS_BOTTOM_estimate

  !-----------------------------------------------------------------------------
  function lagrange_interp( p, x, y ) result(q)
    !$acc routine seq
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
