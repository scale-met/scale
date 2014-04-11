!-------------------------------------------------------------------------------
!> module COUPLER / Atmosphere-Ocean Surface fluxes
!!
!! @par Description
!!          Surface flux from the atmosphere-ocean coupler
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2014-02-26 (T.Yamaura)  [new]
!<
!-------------------------------------------------------------------------------
module scale_cpl_atmos_ocean
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_stdio
  use scale_grid_index
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: CPL_AtmOcn_setup

  abstract interface
     subroutine cao( &
         SST,                                 & ! (inout)
         XMFLX, YMFLX, ZMFLX,                 & ! (out)
         SWUFLX, LWUFLX, SHFLX, LHFLX, WHFLX, & ! (out)
         SST_UPDATE,                          & ! (in)
         DZ, DENS, MOMX, MOMY, MOMZ,          & ! (in)
         RHOS, PRES, ATMP, QV, SWD, LWD,      & ! (in)
         TW, ALB_SW, ALB_LW,                  & ! (in)
         Z0M, Z0H, Z0E                        ) ! (in)
       use scale_precision
       use scale_grid_index
       implicit none

       real(RP), intent(inout) :: SST (IA,JA) ! sea surface temperature [K]

       real(RP), intent(out) :: XMFLX (IA,JA) ! x-momentum flux at the surface [kg/m2/s]
       real(RP), intent(out) :: YMFLX (IA,JA) ! y-momentum flux at the surface [kg/m2/s]
       real(RP), intent(out) :: ZMFLX (IA,JA) ! z-momentum flux at the surface [kg/m2/s]
       real(RP), intent(out) :: SWUFLX(IA,JA) ! upward shortwave flux at the surface [W/m2]
       real(RP), intent(out) :: LWUFLX(IA,JA) ! upward longwave flux at the surface [W/m2]
       real(RP), intent(out) :: SHFLX (IA,JA) ! sensible heat flux at the surface [W/m2]
       real(RP), intent(out) :: LHFLX (IA,JA) ! latent heat flux at the surface [W/m2]
       real(RP), intent(out) :: WHFLX (IA,JA) ! water heat flux at the surface [W/m2]

       logical,  intent(in) :: SST_UPDATE  ! is sea surface temperature updated?

       real(RP), intent(in) :: DZ  (IA,JA) ! height from the surface to the lowest atmospheric layer [m]
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

       real(RP), intent(in) :: TW    (IA,JA) ! water temperature [K]
       real(RP), intent(in) :: ALB_SW(IA,JA) ! surface albedo for SW [0-1]
       real(RP), intent(in) :: ALB_LW(IA,JA) ! surface albedo for LW [0-1]
       real(RP), intent(in) :: Z0M   (IA,JA) ! roughness length for momentum [m]
       real(RP), intent(in) :: Z0H   (IA,JA) ! roughness length for heat [m]
       real(RP), intent(in) :: Z0E   (IA,JA) ! roughness length for vapor [m]
     end subroutine cao
  end interface
  procedure(cao), pointer :: CPL_AtmOcn => NULL()
  public :: CPL_AtmOcn

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

  subroutine CPL_AtmOcn_setup( CPL_TYPE_AtmOcn )
    use scale_process, only: &
       PRC_MPIstop
#define EXTM(pre, name, post) pre ## name ## post
#define NAME(pre, name, post) EXTM(pre, name, post)
#ifdef CAL
    use NAME(scale_cpl_atmos_ocean_, CAO,), only: &
       NAME(CPL_AtmOcn_, CAO, _setup), &
       NAME(CPL_AtmOcn_, CAO,)
#else
    use scale_cpl_atmos_ocean_const, only: &
       CPL_AtmOcn_const_setup, &
       CPL_AtmOcn_const
    use scale_cpl_atmos_ocean_bulk, only: &
       CPL_AtmOcn_bulk_setup, &
       CPL_AtmOcn_bulk
#endif
    implicit none

    character(len=H_SHORT), intent(in) :: CPL_TYPE_AtmOcn
    !---------------------------------------------------------------------------

    select case( CPL_TYPE_AtmOcn )
    case ( 'CONST' )
       call CPL_AtmOcn_const_setup( CPL_TYPE_AtmOcn )
       CPL_AtmOcn => CPL_AtmOcn_const
    case ( 'BULK' )
       call CPL_AtmOcn_bulk_setup( CPL_TYPE_AtmOcn )
       CPL_AtmOcn => CPL_AtmOcn_bulk
    end select

  end subroutine CPL_AtmOcn_setup

end module scale_cpl_atmos_ocean
