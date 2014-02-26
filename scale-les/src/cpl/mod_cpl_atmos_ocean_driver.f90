!-------------------------------------------------------------------------------
!> module COUPLER / Atmosphere-Ocean Driver
!!
!! @par Description
!!          Coupler driver: atmosphere-ocean
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2014-02-26 (T.Yamaura)  [new]
!<
!-------------------------------------------------------------------------------
module mod_cpl_atmos_ocean_driver
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
  public :: CPL_AtmOcn_driver_setup
  public :: CPL_AtmOcn_driver

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

  subroutine CPL_AtmOcn_driver_setup
    use mod_atmos_phy_sf_driver, only: &
       ATMOS_PHY_SF_driver_final
!    use mod_ocean_phy_fixed, only: &
!       OCEAN_PHY_driver_final
    use mod_cpl_vars, only: &
       CPL_TYPE_AtmOcn,    &
       CPL_flushAtm,       &
       CPL_flushOcn,       &
       CPL_flushCPL
    use mod_cpl_atmos_ocean, only: &
       CPL_AtmOcn_setup
    implicit none
    !---------------------------------------------------------------------------

    call CPL_flushAtm
    call CPL_flushOcn
    call CPL_flushCPL

    call ATMOS_PHY_SF_driver_final
!    call OCEAN_PHY_driver_final

    call CPL_AtmOcn_setup( CPL_TYPE_AtmOcn )
    call CPL_AtmOcn_driver( .false. )

    return
  end subroutine CPL_AtmOcn_driver_setup

  subroutine CPL_AtmOcn_driver( update_flag )
    use mod_grid_real, only: &
       CZ => REAL_CZ, &
       FZ => REAL_FZ
    use mod_cpl_vars, only: &
       CPL_AtmOcn_putCPL,     &
       CPL_AtmOcn_getAtm2CPL, &
       CPL_AtmOcn_getOcn2CPL, &
       SST
    use mod_cpl_atmos_ocean, only: &
       CPL_AtmOcn
    implicit none

    ! argument
    logical, intent(in) :: update_flag

    ! work
    real(RP) :: XMFLX (IA,JA) ! x-momentum flux at the surface [kg/m2/s]
    real(RP) :: YMFLX (IA,JA) ! y-momentum flux at the surface [kg/m2/s]
    real(RP) :: ZMFLX (IA,JA) ! z-momentum flux at the surface [kg/m2/s]
    real(RP) :: SWUFLX(IA,JA) ! upward shortwave flux at the surface [W/m2]
    real(RP) :: LWUFLX(IA,JA) ! upward longwave flux at the surface [W/m2]
    real(RP) :: SHFLX (IA,JA) ! sensible heat flux at the surface [W/m2]
    real(RP) :: LHFLX (IA,JA) ! latent heat flux at the surface [W/m2]
    real(RP) :: WHFLX (IA,JA) ! water heat flux at the surface [W/m2]

    real(RP) :: DZ    (IA,JA) ! height from the surface to the lowest atmospheric layer [m]

    real(RP) :: DENS  (IA,JA) ! air density at the lowest atmospheric layer [kg/m3]
    real(RP) :: MOMX  (IA,JA) ! momentum x at the lowest atmospheric layer [kg/m2/s]
    real(RP) :: MOMY  (IA,JA) ! momentum y at the lowest atmospheric layer [kg/m2/s]
    real(RP) :: MOMZ  (IA,JA) ! momentum z at the lowest atmospheric layer [kg/m2/s]
    real(RP) :: RHOS  (IA,JA) ! air density at the sruface [kg/m3]
    real(RP) :: PRES  (IA,JA) ! pressure at the surface [Pa]
    real(RP) :: ATMP  (IA,JA) ! air temperature at the surface [K]
    real(RP) :: QV    (IA,JA) ! ratio of water vapor mass to total mass at the lowest atmospheric layer [kg/kg]
    real(RP) :: PREC  (IA,JA) ! precipitaton flux at the surface [kg/m2/s]
    real(RP) :: SWD   (IA,JA) ! downward short-wave radiation flux at the surface (upward positive) [W/m2]
    real(RP) :: LWD   (IA,JA) ! downward long-wave radiation flux at the surface (upward positive) [W/m2]

    real(RP) :: TW    (IA,JA) ! water temperature [K]
    real(RP) :: ALB   (IA,JA) ! surface albedo in short-wave radiation [no unit]
    real(RP) :: DZW   (IA,JA) ! water depth [m]
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*) '*** Coupler: Atmos-Ocean'

    call CPL_AtmOcn_getAtm2CPL( &
      DENS, MOMX, MOMY, MOMZ, & ! (out)
      RHOS, PRES, ATMP, QV,   & ! (out)
      PREC, SWD, LWD          ) ! (out)

    call CPL_AtmOcn_getOcn2CPL( &
      TW, ALB, DZW ) ! (out)

    DZ(:,:) = CZ(KS,:,:) - FZ(KS-1,:,:)

    call CPL_AtmOcn( &
      SST,                                 & ! (inout)
      update_flag,                         & ! (in)
      XMFLX, YMFLX, ZMFLX,                 & ! (out)
      SWUFLX, LWUFLX, SHFLX, LHFLX, WHFLX, & ! (out)
      DZ, DENS, MOMX, MOMY, MOMZ,          & ! (in)
      RHOS, PRES, ATMP, QV, SWD, LWD,      & ! (in)
      TW, ALB, DZW                         ) ! (in)

    call CPL_AtmOcn_putCPL( &
      XMFLX, YMFLX, ZMFLX,  &
      SWUFLX, LWUFLX,       &
      SHFLX, LHFLX, WHFLX,  &
      PREC                  )

    return
  end subroutine CPL_AtmOcn_driver

end module mod_cpl_atmos_ocean_driver
