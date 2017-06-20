!-------------------------------------------------------------------------------
!> module OCEAN / Surface fluxes
!!
!! @par Description
!!          Ocean surface flux
!!
!! @author Team SCALE
!!
!<
!-------------------------------------------------------------------------------
module scale_ocean_sfc
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
  public :: OCEAN_SFC_setup

  abstract interface
     subroutine ocnsfc( &
           SST_t,  &
           ZMFLX,  &
           XMFLX,  &
           YMFLX,  &
           SHFLX,  &
           LHFLX,  &
           WHFLX,  &
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
           TW,     &
           SST,    &
           ALB_LW, &
           ALB_SW, &
           Z0M,    &
           Z0H,    &
           Z0E,    &
           dt      )
       use scale_precision
       use scale_grid_index
       implicit none

       real(RP), intent(out) :: SST_t(IA,JA) ! tendency of sea surface temperature
       real(RP), intent(out) :: ZMFLX(IA,JA) ! z-momentum flux at the surface [kg/m2/s]
       real(RP), intent(out) :: XMFLX(IA,JA) ! x-momentum flux at the surface [kg/m2/s]
       real(RP), intent(out) :: YMFLX(IA,JA) ! y-momentum flux at the surface [kg/m2/s]
       real(RP), intent(out) :: SHFLX(IA,JA) ! sensible heat flux at the surface [W/m2]
       real(RP), intent(out) :: LHFLX(IA,JA) ! latent heat flux at the surface [W/m2]
       real(RP), intent(out) :: WHFLX(IA,JA) ! water heat flux at the surface [W/m2]
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
       real(RP), intent(in) :: LWD (IA,JA) ! downward long-wave radiation flux at the surface (upward positive) [W/m2]
       real(RP), intent(in) :: SWD (IA,JA) ! downward short-wave radiation flux at the surface (upward positive) [W/m2]

       real(RP), intent(in) :: TW    (IA,JA) ! water temperature [K]
       real(RP), intent(in) :: SST   (IA,JA) ! sea surface temperature [K]
       real(RP), intent(in) :: ALB_LW(IA,JA) ! surface albedo for LW (0-1)
       real(RP), intent(in) :: ALB_SW(IA,JA) ! surface albedo for SW (0-1)
       real(RP), intent(in) :: Z0M   (IA,JA) ! roughness length for momentum [m]
       real(RP), intent(in) :: Z0H   (IA,JA) ! roughness length for heat [m]
       real(RP), intent(in) :: Z0E   (IA,JA) ! roughness length for vapor [m]
       real(DP), intent(in) :: dt            ! delta time
     end subroutine ocnsfc
  end interface

  procedure(ocnsfc), pointer :: OCEAN_SFC => NULL()

  public :: OCEAN_SFC
  public :: OCEAN_SFC_SimpleAlbedo

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

  subroutine OCEAN_SFC_setup( OCEAN_TYPE )
    use scale_process, only: &
       PRC_MPIstop
    use scale_ocean_sfc_const, only: &
       OCEAN_SFC_CONST_setup, &
       OCEAN_SFC_CONST
    use scale_ocean_sfc_slab, only: &
       OCEAN_SFC_SLAB_setup, &
       OCEAN_SFC_SLAB
    implicit none

    character(len=*), intent(in) :: OCEAN_TYPE
    !---------------------------------------------------------------------------

    select case( OCEAN_TYPE )
    case( 'CONST' )
       call OCEAN_SFC_CONST_setup( OCEAN_TYPE )
       OCEAN_SFC => OCEAN_SFC_CONST
    case( 'SLAB', 'FILE' )
       call OCEAN_SFC_SLAB_setup( OCEAN_TYPE )
       OCEAN_SFC => OCEAN_SFC_SLAB
    end select

    return
  end subroutine OCEAN_SFC_setup

  !-----------------------------------------------------------------------------
  subroutine OCEAN_SFC_SimpleAlbedo( &
       SFC_albedo_t, &
       SFC_albedo,   &
       cosSZA,       &
       dt            )
    use scale_grid_index
    use scale_const, only: &
       I_SW  => CONST_I_SW, &
       I_LW  => CONST_I_LW
    implicit none

    ! arguments
    real(RP), intent(out) :: SFC_albedo_t(IA,JA,2) ! tendency of sea surface albedo (0-1)
    real(RP), intent(in)  :: SFC_albedo  (IA,JA,2) ! sea surface                    (0-1)
    real(RP), intent(in)  :: cosSZA      (IA,JA)   ! cos(solar zenith angle)        (0-1)
    real(DP), intent(in)  :: dt                    ! delta time

    ! works
    real(RP) :: SFC_albedo1(IA,JA,2)
    real(RP) :: am1

    real(RP), parameter :: c_ocean_albedo(3) = (/ -0.747900_RP, &
                                                  -4.677039_RP, &
                                                   1.583171_RP  /)

    integer :: i, j
    !---------------------------------------------------------------------------

    do j = JS, JE
    do i = IS, IE
       am1 = max( min( cosSZA(i,j), 0.961_RP ), 0.0349_RP )

       ! SFC_albedo1(i,j,I_LW) = 0.5_RP do nothing
       SFC_albedo1(i,j,I_SW) = exp( ( c_ocean_albedo(3)*am1 + c_ocean_albedo(2) )*am1 + c_ocean_albedo(1) )
    enddo
    enddo

    ! calculate tendency
    do j = JS, JE
    do i = IS, IE
       SFC_albedo_t(i,j,I_LW) = 0.0_RP
       SFC_albedo_t(i,j,I_SW) = ( SFC_albedo1(i,j,I_SW) - SFC_albedo(i,j,I_SW) ) / dt
    enddo
    enddo

    return
  end subroutine OCEAN_SFC_SimpleAlbedo

end module scale_ocean_sfc
