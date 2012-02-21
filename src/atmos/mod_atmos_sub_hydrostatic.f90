!-------------------------------------------------------------------------------
!> module hydrostatic
!!
!! @par Description
!!          make hydrostatic profile in the model
!!
!! @author H.Tomita and SCALE developpers
!!
!! @par History
!! @li      2012-2-20 (H.Yashiro) [new] Extract from the tool of Y.Miyamoto
!!
!<
!-------------------------------------------------------------------------------
module mod_atmos_hydrostatic
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use mod_stdio, only: &
     IO_FID_LOG, &
     IO_L
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: hydro_buildrho

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
  !> Buildup density from surface
  !-----------------------------------------------------------------------------
  subroutine hydro_buildrho( &
      dens,     &
      pres,     &
      pott,     &
      pres_sfc, &
      pott_sfc  )
    use mod_const, only : &
       GRAV   => CONST_GRAV,   &
       Rdry   => CONST_Rdry,   &
       RovCP  => CONST_RovCP,  &
       CPovR  => CONST_CPovR,  &
       RovCV  => CONST_RovCV,  &
       CVovCP => CONST_CVovCP, &
       CPovCV => CONST_CPovCV, &
       Pstd   => CONST_Pstd
    use mod_grid, only : &
       KA  => GRID_KA, &
       KS  => GRID_KS, &
       KE  => GRID_KE, &
       CZ  => GRID_CZ, &
       FDZ => GRID_FDZ
    implicit none

    real(8), intent(out) :: dens(KA) !< density [kg/m3]
    real(8), intent(out) :: pres(KA) !< pressure [Pa]
    real(8), intent(in)  :: pott(KA) !< potential temperature [K]
    real(8), intent(in)  :: pres_sfc !< surface pressure [Pa]
    real(8), intent(in)  :: pott_sfc !< surface potential temperature [K]

    real(8) :: dens_sfc
    real(8) :: RovP, dens_s, dhyd, dgrd

    real(8), parameter :: criteria = 1.D-10
    integer, parameter :: itelim = 100

    integer :: k, ite
    !---------------------------------------------------------------------------

    RovP  = Rdry / Pstd**RovCP

    ! make density at surface
    dens_sfc = Pstd / Rdry / pott(KS) * ( pres_sfc/Pstd )**CVovCP

    ! make density at lowermost cell center
    k = KS

    dens_s  = 0.D0
    dens(k) = dens_sfc ! first guess

    do ite = 1, itelim
       if ( abs(dens(k)-dens_s) <= criteria ) exit

       dens_s = dens(k)

       dhyd = - ( RovP*pott(k) *dens_s   )**CPovCV / CZ(k) &
              + ( RovP*pott_sfc*dens_sfc )**CPovCV / CZ(k) &
              - 0.5D0 * GRAV * ( dens_s + dens_sfc )

       dgrd = - ( RovP*pott(k)           )**CPovCV / CZ(k) &
              * CPovCV * dens_s**RovCV                     &
              - 0.5D0 * GRAV

       dens(k) = dens_s - dhyd/dgrd
    enddo

    if ( ite > itelim ) then
       if( IO_L ) write(IO_FID_LOG,*) 'xxx iteration not converged!', k, ite, dens(k), dens_s, dhyd, dgrd
    endif

    ! make density
    do k = KS+1, KE

       dens_s  = 0.D0
       dens(k) = dens(k-1)

       do ite = 1, itelim
          if ( abs(dens(k)-dens_s) <= criteria ) exit

          dens_s = dens(k)

          dhyd = - ( RovP*pott(k  )*dens_s    )**CPovCV / FDZ(k-1) &
                 + ( RovP*pott(k-1)*dens(k-1) )**CPovCV / FDZ(k-1) &
                 - 0.5D0 * GRAV * ( dens_s + dens(k-1) )

          dgrd = - ( RovP*pott(k)             )**CPovCV / FDZ(k-1) &
                 * CPovCV * dens_s**RovCV &
                 - 0.5D0 * GRAV

          dens(k) = dens_s - dhyd/dgrd
       enddo

       if ( ite > itelim ) then
          if( IO_L ) write(IO_FID_LOG,*) 'xxx iteration not converged!', k, ite, dens(k), dens_s, dhyd, dgrd
       endif
    enddo

    do k = KS, KE
       pres(k) = ( dens(k) * Rdry * pott(k) )**CPovCV * Pstd**(-RovCV)
    enddo

    return
  end subroutine hydro_buildrho

end module mod_atmos_hydrostatic
