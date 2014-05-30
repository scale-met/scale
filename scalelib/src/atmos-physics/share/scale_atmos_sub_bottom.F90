!-------------------------------------------------------------------------------
!> module ATMOSPHERE / Bottom boundary treatment
!!
!! @par Description
!!          Bottom boundary treatment of model domain
!!          Extrapolation of DENS, PRES, TEMP
!!
!! @author Team SCALE
!<
!-------------------------------------------------------------------------------
module scale_atmos_bottom
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_stdio
  use scale_prof
  use scale_grid_index
  use scale_tracer
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
       DENS,     &
       PRES,     &
       CZ,       &
       Zsfc,     &
       Z1,       &
       SFC_DENS, &
       SFC_PRES  )
    use scale_const, only: &
       GRAV  => CONST_GRAV
    implicit none

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
    do j = 1, JA
    do i = 1, IA
       SFC_DENS(i,j) = lagrange_interp( Zsfc(i,j),         & ! [IN]
                                        CZ  (KS:KS+2,i,j), & ! [IN]
                                        DENS(KS:KS+2,i,j)  ) ! [IN]
    enddo
    enddo

    ! estimate surface pressure (hydrostatic balance)
    do j = 1, JA
    do i = 1, IA
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
