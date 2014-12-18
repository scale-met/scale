!-------------------------------------------------------------------------------
!> module LAND / Surface fluxes
!!
!! @par Description
!!          Land surface flux
!!
!! @author Team SCALE
!!
!<
!-------------------------------------------------------------------------------
module scale_land_sfc
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
  public :: LAND_SFC_setup

  abstract interface
     subroutine lndsfc( &
           LST_t,  &
           ZMFLX,  &
           XMFLX,  &
           YMFLX,  &
           SHFLX,  &
           LHFLX,  &
           GHFLX,  &
           U10,    &
           V10,    &
           T2,     &
           Q2,     &
           TMPA,   &
           PRSA,   &
           WA,     &
           UA,     &
           VA,     &
           RHOA,   &
           QVA,    &
           Z1,     &
           PBL,    &
           PRSS,   &
           LWD,    &
           SWD,    &
           TG,     &
           LST,    &
           QVEF,   &
           ALB_LW, &
           ALB_SW, &
           DZG,    &
           TCS,    &
           Z0M,    &
           Z0H,    &
           Z0E,    &
           dt      )
       use scale_precision
       use scale_grid_index
       implicit none

       real(RP), intent(out) :: LST_t(IA,JA) ! tendency of LST
       real(RP), intent(out) :: ZMFLX(IA,JA) ! z-momentum flux at the surface [kg/m2/s]
       real(RP), intent(out) :: XMFLX(IA,JA) ! x-momentum flux at the surface [kg/m2/s]
       real(RP), intent(out) :: YMFLX(IA,JA) ! y-momentum flux at the surface [kg/m2/s]
       real(RP), intent(out) :: SHFLX(IA,JA) ! sensible heat flux at the surface [J/m2/s]
       real(RP), intent(out) :: LHFLX(IA,JA) ! latent heat flux at the surface [J/m2/s]
       real(RP), intent(out) :: GHFLX(IA,JA) ! ground heat flux at the surface [J/m2/s]
       real(RP), intent(out) :: U10  (IA,JA) ! velocity u at 10m [m/s]
       real(RP), intent(out) :: V10  (IA,JA) ! velocity v at 10m [m/s]
       real(RP), intent(out) :: T2   (IA,JA) ! temperature at 2m [K]
       real(RP), intent(out) :: Q2   (IA,JA) ! water vapor at 2m [kg/kg]

       real(RP), intent(in) :: TMPA(IA,JA) ! temperature at the lowest atmospheric layer [K]
       real(RP), intent(in) :: PRSA(IA,JA) ! pressure at the lowest atmospheric layer [Pa]
       real(RP), intent(in) :: WA  (IA,JA) ! velocity w at the lowest atmospheric layer [m/s]
       real(RP), intent(in) :: UA  (IA,JA) ! velocity u at the lowest atmospheric layer [m/s]
       real(RP), intent(in) :: VA  (IA,JA) ! velocity v at the lowest atmospheric layer [m/s]
       real(RP), intent(in) :: RHOA(IA,JA) ! density at the lowest atmospheric layer [kg/m3]
       real(RP), intent(in) :: QVA (IA,JA) ! ratio of water vapor mass to total mass at the lowest atmospheric layer [kg/kg]
       real(RP), intent(in) :: Z1  (IA,JA) ! cell center height at the lowest atmospheric layer [m]
       real(RP), intent(in) :: PBL (IA,JA) ! the top of atmospheric mixing layer [m]
       real(RP), intent(in) :: PRSS(IA,JA) ! pressure at the surface [Pa]
       real(RP), intent(in) :: LWD (IA,JA) ! downward long-wave radiation flux at the surface (upward positive) [J/m2/s]
       real(RP), intent(in) :: SWD (IA,JA) ! downward short-wave radiation flux at the surface (upward positive) [J/m2/s]

       real(RP), intent(in) :: TG    (IA,JA) ! soil temperature [K]
       real(RP), intent(in) :: LST   (IA,JA) ! land surface temperature [K]
       real(RP), intent(in) :: QVEF  (IA,JA) ! efficiency of evaporation [0-1]
       real(RP), intent(in) :: ALB_LW(IA,JA) ! surface albedo for LW [0-1]
       real(RP), intent(in) :: ALB_SW(IA,JA) ! surface albedo for SW [0-1]
       real(RP), intent(in) :: DZG   (IA,JA) ! soil depth [m]
       real(RP), intent(in) :: TCS   (IA,JA) ! thermal conductivity for soil [J/m/K/s]
       real(RP), intent(in) :: Z0M   (IA,JA) ! roughness length for momemtum [m]
       real(RP), intent(in) :: Z0H   (IA,JA) ! roughness length for heat [m]
       real(RP), intent(in) :: Z0E   (IA,JA) ! roughness length for vapor [m]
       real(RP), intent(in) :: dt            ! delta time
     end subroutine lndsfc
  end interface
  procedure(lndsfc), pointer :: LAND_SFC => NULL()
  public :: LAND_SFC

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

  subroutine LAND_SFC_setup( LAND_TYPE )
    use scale_process, only: &
       PRC_MPIstop
    use scale_land_sfc_slab, only: &
       LAND_SFC_SLAB_setup, &
       LAND_SFC_SLAB
    implicit none

    character(len=*), intent(in) :: LAND_TYPE
    !---------------------------------------------------------------------------

    select case( LAND_TYPE )
    case ( 'CONST' )
       call LAND_SFC_SLAB_setup( LAND_TYPE )
       LAND_SFC => LAND_SFC_SLAB
    case ( 'BULK' )
       call LAND_SFC_SLAB_setup( LAND_TYPE )
       LAND_SFC => LAND_SFC_SLAB
    end select

  end subroutine LAND_SFC_setup

end module scale_land_sfc
