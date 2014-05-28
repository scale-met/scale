!-------------------------------------------------------------------------------
!> module COUPLER / Atmosphere-Land Surface fluxes
!!
!! @par Description
!!          Surface flux from the atmosphere-land coupler
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2014-02-25 (T.Yamaura)  [new]
!<
!-------------------------------------------------------------------------------
module scale_cpl_atmos_land
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
  public :: CPL_AtmLnd_setup

  abstract interface
     subroutine cal( &
         LST,        & ! (inout)
         XMFLX,      & ! (out)
         YMFLX,      & ! (out)
         ZMFLX,      & ! (out)
         SHFLX,      & ! (out)
         LHFLX,      & ! (out)
         GHFLX,      & ! (out)
         U10,        & ! (out)
         V10,        & ! (out)
         T2,         & ! (out)
         Q2,         & ! (out)
         LST_UPDATE, & ! (in)
         RHOA,       & ! (in)
         UA,         & ! (in)
         VA,         & ! (in)
         WA,         & ! (in)
         TMPA,       & ! (in)
         PRSA,       & ! (in)
         QVA,        & ! (in)
         PRSS,       & ! (in)
         SWD,        & ! (in)
         LWD,        & ! (in)
         TG,         & ! (in)
         QVEF,       & ! (in)
         ALB_SW,     & ! (in)
         ALB_LW,     & ! (in)
         TCS,        & ! (in)
         DZG,        & ! (in)
         Z0M,        & ! (in)
         Z0H,        & ! (in)
         Z0E         ) ! (in)
       use scale_precision
       use scale_grid_index
       implicit none

       real(RP), intent(inout) :: LST(IA,JA) ! land surface temperature [K]

       real(RP), intent(out) :: XMFLX(IA,JA) ! x-momentum flux at the surface [kg/m2/s]
       real(RP), intent(out) :: YMFLX(IA,JA) ! y-momentum flux at the surface [kg/m2/s]
       real(RP), intent(out) :: ZMFLX(IA,JA) ! z-momentum flux at the surface [kg/m2/s]
       real(RP), intent(out) :: SHFLX(IA,JA) ! sensible heat flux at the surface [W/m2]
       real(RP), intent(out) :: LHFLX(IA,JA) ! latent heat flux at the surface [W/m2]
       real(RP), intent(out) :: GHFLX(IA,JA) ! ground heat flux at the surface [W/m2]
       real(RP), intent(out) :: U10  (IA,JA) ! velocity u at 10m [m/s]
       real(RP), intent(out) :: V10  (IA,JA) ! velocity v at 10m [m/s]
       real(RP), intent(out) :: T2   (IA,JA) ! temperature at 2m [K]
       real(RP), intent(out) :: Q2   (IA,JA) ! water vapor at 2m [kg/kg]

       logical,  intent(in) :: LST_UPDATE  ! is land surface temperature updated?

       real(RP), intent(in) :: RHOA(IA,JA) ! density at the lowest atmospheric layer [kg/m3]
       real(RP), intent(in) :: UA  (IA,JA) ! velocity u at the lowest atmospheric layer [m/s]
       real(RP), intent(in) :: VA  (IA,JA) ! velocity v at the lowest atmospheric layer [m/s]
       real(RP), intent(in) :: WA  (IA,JA) ! velocity w at the lowest atmospheric layer [m/s]
       real(RP), intent(in) :: TMPA(IA,JA) ! temperature at the lowest atmospheric layer [K]
       real(RP), intent(in) :: PRSA(IA,JA) ! pressure at the lowest atmospheric layer [Pa]
       real(RP), intent(in) :: QVA (IA,JA) ! ratio of water vapor mass to total mass at the lowest atmospheric layer [kg/kg]
       real(RP), intent(in) :: PRSS(IA,JA) ! pressure at the surface [Pa]
       real(RP), intent(in) :: SWD (IA,JA) ! downward short-wave radiation flux at the surface (upward positive) [W/m2]
       real(RP), intent(in) :: LWD (IA,JA) ! downward long-wave radiation flux at the surface (upward positive) [W/m2]

       real(RP), intent(in) :: TG    (IA,JA) ! soil temperature [K]
       real(RP), intent(in) :: QVEF  (IA,JA) ! efficiency of evaporation [0-1]
       real(RP), intent(in) :: ALB_SW(IA,JA) ! surface albedo for SW [0-1]
       real(RP), intent(in) :: ALB_LW(IA,JA) ! surface albedo for LW [0-1]
       real(RP), intent(in) :: TCS   (IA,JA) ! thermal conductivity for soil [W/m/K]
       real(RP), intent(in) :: DZG   (IA,JA) ! soil depth [m]
       real(RP), intent(in) :: Z0M   (IA,JA) ! roughness length for momemtum [m]
       real(RP), intent(in) :: Z0H   (IA,JA) ! roughness length for heat [m]
       real(RP), intent(in) :: Z0E   (IA,JA) ! roughness length for vapor [m]
     end subroutine cal
  end interface
  procedure(cal), pointer :: CPL_AtmLnd => NULL()
  public :: CPL_AtmLnd

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

  subroutine CPL_AtmLnd_setup( CPL_TYPE_AtmLnd )
    use scale_process, only: &
       PRC_MPIstop
#define EXTM(pre, name, post) pre ## name ## post
#define NAME(pre, name, post) EXTM(pre, name, post)
#ifdef CAL
    use NAME(scale_cpl_atmos_land_, CAL,), only: &
       NAME(CPL_AtmLnd_, CAL, _setup), &
       NAME(CPL_AtmLnd_, CAL,)
#else
    use scale_cpl_atmos_land_const, only: &
       CPL_AtmLnd_const_setup, &
       CPL_AtmLnd_const
    use scale_cpl_atmos_land_bulk, only: &
       CPL_AtmLnd_bulk_setup, &
       CPL_AtmLnd_bulk
#endif
    implicit none

    character(len=*), intent(in) :: CPL_TYPE_AtmLnd
    !---------------------------------------------------------------------------

    select case( CPL_TYPE_AtmLnd )
    case ( 'CONST' )
       call CPL_AtmLnd_const_setup( CPL_TYPE_AtmLnd )
       CPL_AtmLnd => CPL_AtmLnd_const
    case ( 'BULK' )
       call CPL_AtmLnd_bulk_setup( CPL_TYPE_AtmLnd )
       CPL_AtmLnd => CPL_AtmLnd_bulk
    end select

  end subroutine CPL_AtmLnd_setup

end module scale_cpl_atmos_land
