!-------------------------------------------------------------------------------
!> module COUPLER / Atmosphere-Land Driver
!!
!! @par Description
!!          Coupler driver: atmosphere-land
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2014-02-25 (T.Yamaura)  [new]
!<
!-------------------------------------------------------------------------------
module mod_cpl_atmos_land_driver
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
  public :: CPL_AtmLnd_driver_setup
  public :: CPL_AtmLnd_driver

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

  subroutine CPL_AtmLnd_driver_setup
    use mod_atmos_phy_sf_driver, only: &
       ATMOS_PHY_SF_driver_final
    use mod_land_phy_bucket, only: &
       LAND_PHY_driver_final
    use mod_cpl_vars, only: &
       CPL_TYPE_AtmLnd,    &
       CPL_flushAtm,       &
       CPL_flushLnd,       &
       CPL_AtmLnd_flushCPL
    use mod_cpl_atmos_land, only: &
       CPL_AtmLnd_setup
    implicit none
    !---------------------------------------------------------------------------

    call CPL_flushAtm
    call CPL_flushLnd
    call CPL_AtmLnd_flushCPL

    call ATMOS_PHY_SF_driver_final
    call LAND_PHY_driver_final

    call CPL_AtmLnd_setup( CPL_TYPE_AtmLnd )
    call CPL_AtmLnd_driver( .false. )

    return
  end subroutine CPL_AtmLnd_driver_setup

  subroutine CPL_AtmLnd_driver( update_flag )
    use mod_grid_real, only: &
       CZ => REAL_CZ, &
       FZ => REAL_FZ
    use mod_cpl_vars, only: &
       CPL_AtmLnd_putCPL,     &
       CPL_AtmLnd_getAtm2CPL, &
       CPL_AtmLnd_getLnd2CPL, &
       LST
    use mod_cpl_atmos_land, only: &
       CPL_AtmLnd
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
    real(RP) :: GHFLX (IA,JA) ! ground heat flux at the surface [W/m2]

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

    real(RP) :: TG    (IA,JA) ! soil temperature [K]
    real(RP) :: QVEF  (IA,JA) ! efficiency of evaporation [no unit]
    real(RP) :: EMIT  (IA,JA) ! emissivity in long-wave radiation [no unit]
    real(RP) :: ALB   (IA,JA) ! surface albedo in short-wave radiation [no unit]
    real(RP) :: TCS   (IA,JA) ! thermal conductivity for soil [W/m/K]
    real(RP) :: DZG   (IA,JA) ! soil depth [m]
    real(RP) :: Z0M   (IA,JA) ! roughness length for momemtum [m]
    real(RP) :: Z0H   (IA,JA) ! roughness length for heat [m]
    real(RP) :: Z0E   (IA,JA) ! roughness length for vapor [m]
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*) '*** Coupler: Atmos-Land'

    call CPL_AtmLnd_getAtm2CPL( &
      DENS, MOMX, MOMY, MOMZ, & ! (out)
      RHOS, PRES, ATMP, QV,   & ! (out)
      PREC, SWD, LWD          ) ! (out)

    call CPL_AtmLnd_getLnd2CPL( &
      TG, QVEF, EMIT, & ! (out)
      ALB, TCS, DZG,  & ! (out)
      Z0M, Z0H, Z0E   ) ! (out)

    DZ(:,:) = CZ(KS,:,:) - FZ(KS-1,:,:)

    call CPL_AtmLnd( &
      LST,                                 & ! (inout)
      update_flag,                         & ! (in)
      XMFLX, YMFLX, ZMFLX,                 & ! (out)
      SWUFLX, LWUFLX, SHFLX, LHFLX, GHFLX, & ! (out)
      DZ, DENS, MOMX, MOMY, MOMZ,          & ! (in)
      RHOS, PRES, ATMP, QV, SWD, LWD,      & ! (in)
      TG, QVEF, EMIT, ALB,                 & ! (in)
      TCS, DZG, Z0M, Z0H, Z0E              ) ! (in)

    call CPL_AtmLnd_putCPL( &
      XMFLX, YMFLX, ZMFLX,  &
      SWUFLX, LWUFLX,       &
      SHFLX, LHFLX, GHFLX,  &
      PREC                  )

    return
  end subroutine CPL_AtmLnd_driver

end module mod_cpl_atmos_land_driver
