!-------------------------------------------------------------------------------
!> module COUPLER / Surface fluxes
!!
!! @par Description
!!          Surface flux with Bulk Method
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2014-02-21 (T.Yamaura)  [new]
!<
!-------------------------------------------------------------------------------
module mod_cpl_sfcflux
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use mod_precision
  use mod_stdio
  use mod_grid_index
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: CPL_sfcflux

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
  subroutine CPL_sfcflux( &
      RES, DRES,                                              & ! (out)
      XMFLX, YMFLX, ZMFLX,                                    & ! (out)
      SWUFLX, LWUFLX, SHFLX, LHFLX, GHFLX,                    & ! (out)
      DZ, TS,                                                 & ! (in)
      DENS, MOMX, MOMY, MOMZ, RHOS, PRES, ATMP, QV, SWD, LWD, & ! (in)
      TG, QVEF, EMIT, ALB, TCS, DZG, Z0M, Z0H, Z0E            ) ! (in)
    use mod_const, only: &
      GRAV   => CONST_GRAV,  &
      CPdry  => CONST_CPdry, &
      Rvap   => CONST_Rvap,  &
      STB    => CONST_STB,   &
      LH0    => CONST_LH0,   &
      P00    => CONST_PRE00
    use mod_atmos_saturation, only: &
      qsat => ATMOS_SATURATION_pres2qsat_all
    use mod_cpl_bulkcoef, only: &
      CPL_bulkcoef
    implicit none

    ! argument
    real(RP), intent(out) :: RES   (IA,JA) ! residual in the equation of heat balance
    real(RP), intent(out) :: DRES  (IA,JA) ! d(residual) / d(Ts)
    real(RP), intent(out) :: XMFLX (IA,JA) ! x-momentum flux at the surface [kg/m2/s]
    real(RP), intent(out) :: YMFLX (IA,JA) ! y-momentum flux at the surface [kg/m2/s]
    real(RP), intent(out) :: ZMFLX (IA,JA) ! z-momentum flux at the surface [kg/m2/s]
    real(RP), intent(out) :: SWUFLX(IA,JA) ! upward shortwave flux at the surface [W/m2]
    real(RP), intent(out) :: LWUFLX(IA,JA) ! upward longwave flux at the surface [W/m2]
    real(RP), intent(out) :: SHFLX (IA,JA) ! sensible heat flux at the surface [W/m2]
    real(RP), intent(out) :: LHFLX (IA,JA) ! latent heat flux at the surface [W/m2]
    real(RP), intent(out) :: GHFLX (IA,JA) ! ground heat flux at the surface [W/m2]

    real(RP), intent(in) :: DZ  (IA,JA) ! height from the surface to the lowest atmospheric layer [m]
    real(RP), intent(in) :: TS  (IA,JA) ! skin temperature [K]

    real(RP), intent(in) :: DENS(IA,JA) ! air density at the lowest atmospheric layer [kg/m3]
    real(RP), intent(in) :: MOMX(IA,JA) ! momentum x at the lowest atmospheric layer [kg/m2/s]
    real(RP), intent(in) :: MOMY(IA,JA) ! momentum y at the lowest atmospheric layer [kg/m2/s]
    real(RP), intent(in) :: MOMZ(IA,JA) ! momentum z at the lowest atmospheric layer [kg/m2/s]
    real(RP), intent(in) :: RHOS(IA,JA) ! air density at the sruface [kg/m3]
    real(RP), intent(in) :: PRES(IA,JA) ! pressure at the surface [Pa]
    real(RP), intent(in) :: ATMP(IA,JA) ! air temperature at the surface [K]
    real(RP), intent(in) :: QV  (IA,JA) ! ratio of water vapor mass to total mass at the lowest atmospheric layer [kg/kg]
    real(RP), intent(in) :: SWD (IA,JA) ! downward short-wave radiation flux at the surface (upward positive) [W/m2]
    real(RP), intent(in) :: LWD (IA,JA) ! downward long-wave radiation flux at the surface (upward positive) [W/m2]

    real(RP), intent(in) :: TG  (IA,JA) ! soil temperature [K]
    real(RP), intent(in) :: QVEF(IA,JA) ! efficiency of evaporation [no unit]
    real(RP), intent(in) :: EMIT(IA,JA) ! emissivity in long-wave radiation [no unit]
    real(RP), intent(in) :: ALB (IA,JA) ! surface albedo in short-wave radiation [no unit]
    real(RP), intent(in) :: TCS (IA,JA) ! thermal conductivity for soil [W/m/K]
    real(RP), intent(in) :: DZG (IA,JA) ! soil depth [m]
    real(RP), intent(in) :: Z0M (IA,JA) ! roughness length for momemtum [m]
    real(RP), intent(in) :: Z0H (IA,JA) ! roughness length for heat [m]
    real(RP), intent(in) :: Z0E (IA,JA) ! roughness length for vapor [m]

    ! constant
    real(RP), parameter :: dTS    =   1.0E-8_RP ! delta TS
    real(RP), parameter :: U_minM =   0.0_RP    ! minimum U_abs for u,v,w
    real(RP), parameter :: U_minH =   0.0_RP    !                   T
    real(RP), parameter :: U_minE =   0.0_RP    !                   q
    real(RP), parameter :: U_maxM = 100.0_RP    ! maximum U_abs for u,v,w
    real(RP), parameter :: U_maxH = 100.0_RP    !                   T
    real(RP), parameter :: U_maxE = 100.0_RP    !                   q

    ! work
    real(RP) :: Uabs ! absolute velocity at the lowest atmospheric layer [m/s]
    real(RP) :: Cm, Ch, Ce ! bulk transfer coeff. [no unit]
    real(RP) :: dCm, dCh, dCe

    real(RP) :: SQV, dSQV ! saturation water vapor mixing ratio at surface [kg/kg]
    real(RP) :: dLWUFLX, dGHFLX, dSHFLX, dLHFLX

    integer :: i, j
    !---------------------------------------------------------------------------

    ! at (u, y, layer)
    do j = JS, JE
    do i = IS, IE
      Uabs = sqrt( &
             ( 0.5_RP * ( MOMZ(i,j) + MOMZ(i+1,j)                               ) )**2 &
           + ( 2.0_RP *   MOMX(i,j)                                               )**2 &
           + ( 0.5_RP * ( MOMY(i,j-1) + MOMY(i,j) + MOMY(i+1,j-1) + MOMY(i+1,j) ) )**2 &
           ) / ( DENS(i,j) + DENS(i+1,j) )

      call CPL_bulkcoef( &
          Cm, Ch, Ce,                           & ! (out)
          ( ATMP(i,j) + ATMP(i+1,j) ) * 0.5_RP, & ! (in)
          ( TS  (i,j) + TS  (i+1,j) ) * 0.5_RP, & ! (in)
          DZ(i,j), Uabs,                        & ! (in)
          Z0M(i,j), Z0H(i,j), Z0E(i,j)          ) ! (in)

      XMFLX(i,j) = - Cm * min(max(Uabs,U_minM),U_maxM) * MOMX(i,j)
    enddo
    enddo

    ! at (x, v, layer)
    do j = JS, JE
    do i = IS, IE
      Uabs = sqrt( &
             ( 0.5_RP * ( MOMZ(i,j) + MOMZ(i,j+1)                               ) )**2 &
           + ( 0.5_RP * ( MOMX(i-1,j) + MOMX(i,j) + MOMX(i-1,j+1) + MOMX(i,j+1) ) )**2 &
           + ( 2.0_RP *   MOMY(i,j)                                               )**2 &
           ) / ( DENS(i,j) + DENS(i,j+1) )

      call CPL_bulkcoef( &
          Cm, Ch, Ce,                           & ! (out)
          ( ATMP(i,j) + ATMP(i,j+1) ) * 0.5_RP, & ! (in)
          ( TS  (i,j) + TS  (i,j+1) ) * 0.5_RP, & ! (in)
          DZ(i,j), Uabs,                        & ! (in)
          Z0M(i,j), Z0H(i,j), Z0E(i,j)          ) ! (in)

      YMFLX(i,j) = - Cm * min(max(Uabs,U_minM),U_maxM) * MOMY(i,j)
    enddo
    enddo

    ! at cell center
    do j = JS-1, JE+1
    do i = IS-1, IE+1
      Uabs = sqrt( &
             ( MOMZ(i,j)               )**2 &
           + ( MOMX(i-1,j) + MOMX(i,j) )**2 &
           + ( MOMY(i,j-1) + MOMY(i,j) )**2 &
           ) / DENS(i,j) * 0.5_RP

      call CPL_bulkcoef( &
          Cm, Ch, Ce,                  & ! (out)
          ATMP(i,j), TS(i,j),          & ! (in)
          DZ(i,j), Uabs,               & ! (in)
          Z0M(i,j), Z0H(i,j), Z0E(i,j) ) ! (in)

      ZMFLX(i,j) = - Cm * min(max(Uabs,U_minM),U_maxM) * MOMZ(i,j) * 0.5_RP

      ! saturation at the surface
      call qsat( SQV, TS(i,j), PRES(i,j) )

      SHFLX (i,j) = CPdry * min(max(Uabs,U_minH),U_maxH) * RHOS(i,j) * Ch * ( TS(i,j) - ATMP(i,j) )
      LHFLX (i,j) = LH0   * min(max(Uabs,U_minE),U_maxE) * RHOS(i,j) * QVEF(i,j) * Ce * ( SQV - QV(i,j) )
      GHFLX (i,j) = -2.0_RP * TCS(i,j) * ( TS(i,j) - TG(i,j)  ) / DZG(i,j)
      SWUFLX(i,j) = ALB(i,j) * SWD(i,j)
      LWUFLX(i,j) = EMIT(i,j) * STB * TS(i,j)**4

      ! calculation for residual
      RES(i,j) = SWD(i,j) - SWUFLX(i,j) + LWD(i,j) - LWUFLX(i,j) - SHFLX(i,j) - LHFLX(i,j) + GHFLX(i,j)

      call CPL_bulkcoef( &
          dCm, dCh, dCe,               & ! (out)
          ATMP(i,j), TS(i,j)+dTS,      & ! (in)
          DZ(i,j), Uabs,               & ! (in)
          Z0M(i,j), Z0H(i,j), Z0E(i,j) ) ! (in)

      call qsat( dSQV, TS(i,j)+dTS, PRES(i,j) )

      dSHFLX  = CPdry * min(max(Uabs,U_minH),U_maxH) * RHOS(i,j) * ( (dCh-Ch)/dTS * ( TS(i,j) - ATMP(i,j) ) + Ch )
      dLHFLX  = LH0   * min(max(Uabs,U_minE),U_maxE) * RHOS(i,j) * QVEF(i,j) * ( (dCe-Ce)/dTS * ( SQV - QV(i,j) ) + Ce * (dSQV-SQV)/dTS )
      dGHFLX  = -2.0_RP * TCS(i,j) / DZG(i,j)
      dLWUFLX = 4.0_RP * EMIT(i,j) * STB * TS(i,j)**3

      ! calculation for d(residual)/dTS
      DRES(i,j) = - dLWUFLX - dSHFLX - dLHFLX + dGHFLX
    enddo
    enddo

    return
  end subroutine CPL_sfcflux

end module mod_cpl_sfcflux
